"""Univariate spectrophotometric methods.

Most published univariate methods for mixtures are a *processing pipeline*
followed by a *measurement* and a linear calibration:

======================================  =========================================
Method                                  Pipeline → measurement
======================================  =========================================
Direct (λmax)                           – → amplitude(λ)
Zero-crossing derivative (D1…D4)        derivative → amplitude(λ zero-crossing)
Dual wavelength                         – → difference(λ1, λ2)
Induced dual wavelength                 – → weighted difference(λ1, λ2, F)
Area under curve                        – → area(λ1–λ2)
Derivative ratio (DD1)                  divide → derivative → amplitude(λ)
Ratio difference (RD)                   divide → difference(λ1, λ2)
Mean centering of ratio spectra         divide → mean-center → amplitude(λ)
Successive derivative ratio             divide → derivative → divide → derivative
Double divisor ratio derivative         divide(sum divisor) → derivative → amplitude
Ratio subtraction / constant mult.      ratio_subtraction / constant_multiplication
                                        → amplitude(λmax)
Spectrum subtraction / factorized       spectrum_subtraction / factorized_recovery
                                        → amplitude(λmax)
Constant center                         constant_center → amplitude(λmax)
Peak-to-peak derivative                 derivative → peak_to_peak(λ1, λ2)
======================================  =========================================

Methods that need several signals at once (Q-analysis, absorbance
subtraction, amplitude modulation, H-point standard addition) are provided as
dedicated functions below. Multi-wavelength equation methods (Vierordt,
bivariate, AUC-equations) are in :mod:`spectro.core.multicomponent`.
"""

from __future__ import annotations

import copy
from dataclasses import dataclass, field
from typing import Any

import numpy as np

from spectro.core.operations import Param, Resolver, _on_grid, apply_pipeline
from spectro.core.spectrum import Spectrum
from spectro.core.validation import Regression, linear_regression


# --------------------------------------------------------------------------- #
# Measurements
# --------------------------------------------------------------------------- #
@dataclass
class Measurement:
    key: str
    label: str
    params: list[Param]
    description: str = ""


def _w(name, label, default=None, **kw):
    return Param(name, label, "wavelength", default, **kw)


MEASUREMENTS: dict[str, Measurement] = {m.key: m for m in [
    Measurement("amplitude", "Amplitude at λ", [_w("w1", "λ (nm)", 260.0)],
                "Signal value at one wavelength (λmax, zero-crossing point…)."),
    Measurement("difference", "Difference P(λ1) − P(λ2)",
                [_w("w1", "λ1 (nm)", 250.0), _w("w2", "λ2 (nm)", 280.0)],
                "Dual wavelength / ratio difference / derivative difference."),
    Measurement("weighted_difference", "Induced dual wavelength P(λ1) − F·P(λ2)",
                [_w("w1", "λ1 (nm)", 250.0), _w("w2", "λ2 (nm)", 280.0),
                 Param("factor", "Equality factor F (blank = from reference)", "float",
                       None, optional=True),
                 Param("reference", "Interferent spectrum for F", "spectrum", None,
                       optional=True)],
                "F = P_interferent(λ1)/P_interferent(λ2) cancels the interferent."),
    Measurement("ratio", "Ratio P(λ1) / P(λ2)",
                [_w("w1", "λ1 (nm)", 250.0), _w("w2", "λ2 (nm)", 280.0)]),
    Measurement("peak_to_peak", "Peak-to-trough in range (max − min)",
                [_w("w1", "From (nm)", 240.0), _w("w2", "To (nm)", 300.0)]),
    Measurement("maximum", "Maximum in range", [_w("w1", "From (nm)", 240.0),
                                                _w("w2", "To (nm)", 300.0)]),
    Measurement("minimum", "Minimum in range", [_w("w1", "From (nm)", 240.0),
                                                _w("w2", "To (nm)", 300.0)]),
    Measurement("mean", "Mean in range (plateau)", [_w("w1", "From (nm)", 300.0),
                                                    _w("w2", "To (nm)", 320.0)]),
    Measurement("area", "Area under curve (AUC)",
                [_w("w1", "From (nm)", 240.0), _w("w2", "To (nm)", 300.0),
                 Param("baseline", "Baseline", "choice", "zero", choices=("zero", "linear"))],
                "Trapezoidal integral between λ1 and λ2."),
]}


def measure(s: Spectrum, spec: dict, resolve: Resolver | None = None) -> float:
    kind = spec["kind"]
    p = spec.get("params", {})
    w1, w2 = p.get("w1"), p.get("w2")
    if kind == "amplitude":
        return s.value_at(w1)
    if kind == "difference":
        return s.value_at(w1) - s.value_at(w2)
    if kind == "weighted_difference":
        f = p.get("factor")
        if f is None:
            ref = p.get("reference")
            if ref is None:
                raise ValueError("give an equality factor or an interferent spectrum")
            ref = ref if isinstance(ref, Spectrum) else resolve(ref)
            f = ref.value_at(w1) / ref.value_at(w2)
        return s.value_at(w1) - float(f) * s.value_at(w2)
    if kind == "ratio":
        return s.value_at(w1) / s.value_at(w2)
    _, y = s.region(w1, w2)
    if kind == "peak_to_peak":
        return float(y.max() - y.min())
    if kind == "maximum":
        return float(y.max())
    if kind == "minimum":
        return float(y.min())
    if kind == "mean":
        return float(y.mean())
    if kind == "area":
        x, y = s.region(w1, w2)
        area = float(np.trapezoid(y, x))
        if p.get("baseline") == "linear":
            area -= 0.5 * (y[0] + y[-1]) * (x[-1] - x[0])
        return area
    raise KeyError(f"unknown measurement {kind!r}")


# --------------------------------------------------------------------------- #
# Generic univariate calibration
# --------------------------------------------------------------------------- #
@dataclass
class UnivariateMethod:
    """A pipeline + measurement + linear calibration for one compound."""

    name: str
    compound: str
    steps: list[dict] = field(default_factory=list)
    measurement: dict = field(default_factory=lambda: {"kind": "amplitude",
                                                       "params": {"w1": 260.0}})
    through_origin: bool = False
    regression: Regression | None = None
    # concentration value method: the signal (e.g. the plateau of a ratio
    # spectrum with a unit-concentration divisor) IS the concentration; the
    # calibration line is still computed and reported, but not used to predict
    direct: bool = False

    def signal(self, s: Spectrum, resolve: Resolver | None = None) -> float:
        return measure(apply_pipeline(s, self.steps, resolve), self._measurement(resolve),
                       resolve)

    def _measurement(self, resolve: Resolver | None) -> dict:
        """For P(λ1) − F·P(λ2) with an interferent reference, F is computed on
        the interferent *after the same processing* (e.g. its ratio spectrum),
        which is what the dual amplitude difference method requires."""
        m = self.measurement
        p = m.get("params", {})
        if m["kind"] != "weighted_difference" or p.get("factor") is not None \
                or p.get("reference") is None:
            return m
        ref = p["reference"]
        ref = ref if isinstance(ref, Spectrum) else resolve(ref)
        proc = apply_pipeline(ref, self.steps, resolve)
        f = proc.value_at(p["w1"]) / proc.value_at(p["w2"])
        return {"kind": m["kind"], "params": {**p, "factor": float(f)}}

    def calibrate(self, spectra: list[Spectrum], resolve: Resolver | None = None,
                  concentrations: list[float] | None = None) -> Regression:
        if concentrations is None:
            missing = [s.name for s in spectra if self.compound not in s.concentrations]
            if missing:
                raise ValueError(f"no {self.compound} concentration for: {', '.join(missing)}")
            concentrations = [s.concentrations[self.compound] for s in spectra]
        signals = [self.signal(s, resolve) for s in spectra]
        self.regression = linear_regression(concentrations, signals, self.through_origin)
        return self.regression

    def predict(self, s: Spectrum, resolve: Resolver | None = None) -> float:
        if self.direct:
            return float(self.signal(s, resolve))
        if self.regression is None:
            raise RuntimeError("method is not calibrated")
        return float(self.regression.predict_x(self.signal(s, resolve)))

    def concentration(self, signal: float) -> float:
        """Concentration from a signal (calibration line, or the signal itself
        for the concentration value method)."""
        if self.direct:
            return float(signal)
        if self.regression is None:
            raise RuntimeError("method is not calibrated")
        return float(self.regression.predict_x(signal))

    def to_dict(self) -> dict:
        return {"type": "univariate", "name": self.name, "compound": self.compound,
                "steps": self.steps, "measurement": self.measurement,
                "through_origin": self.through_origin, "direct": self.direct,
                "regression": self.regression.to_dict() if self.regression else None}

    @classmethod
    def from_dict(cls, d: dict) -> "UnivariateMethod":
        reg = d.get("regression")
        return cls(d["name"], d["compound"], d.get("steps", []), d["measurement"],
                   d.get("through_origin", False),
                   Regression.from_dict(reg) if reg else None, d.get("direct", False))


# --------------------------------------------------------------------------- #
# Method templates (pre-filled pipelines for the UI wizard)
# --------------------------------------------------------------------------- #
_D1 = {"op": "derivative", "params": {"order": 1, "delta_lambda": 4.0, "scaling": 10.0}}

TEMPLATES: dict[str, dict[str, Any]] = {
    "Direct (zero order, λmax)": {
        "steps": [], "measurement": {"kind": "amplitude", "params": {"w1": 260.0}}},
    "Zero-crossing derivative (D1)": {
        "steps": [{"op": "derivative", "params": {"order": 1, "delta_lambda": 4.0}}],
        "measurement": {"kind": "amplitude", "params": {"w1": 260.0}}},
    "Zero-crossing derivative (D2)": {
        "steps": [{"op": "derivative", "params": {"order": 2, "delta_lambda": 4.0}}],
        "measurement": {"kind": "amplitude", "params": {"w1": 260.0}}},
    "Derivative peak-to-peak": {
        "steps": [{"op": "derivative", "params": {"order": 1, "delta_lambda": 4.0}}],
        "measurement": {"kind": "peak_to_peak", "params": {"w1": 240.0, "w2": 300.0}}},
    "Dual wavelength": {
        "steps": [], "measurement": {"kind": "difference",
                                     "params": {"w1": 250.0, "w2": 280.0}}},
    "Induced dual wavelength": {
        "steps": [], "measurement": {"kind": "weighted_difference",
                                     "params": {"w1": 250.0, "w2": 280.0}}},
    "Dual wavelength in derivative mode": {
        "steps": [{"op": "derivative", "params": {"order": 1, "delta_lambda": 4.0}}],
        "measurement": {"kind": "difference", "params": {"w1": 250.0, "w2": 280.0}}},
    "Area under curve (single component)": {
        "steps": [], "measurement": {"kind": "area",
                                     "params": {"w1": 240.0, "w2": 280.0, "baseline": "zero"}}},
    "Derivative ratio (DD1)": {
        "steps": [{"op": "divide", "params": {"reference": None}},
                  {"op": "derivative", "params": {"order": 1, "delta_lambda": 4.0}}],
        "measurement": {"kind": "amplitude", "params": {"w1": 260.0}}},
    "Ratio difference (RD)": {
        "steps": [{"op": "divide", "params": {"reference": None}}],
        "measurement": {"kind": "difference", "params": {"w1": 250.0, "w2": 280.0}}},
    "Mean centering of ratio spectra (MCR)": {
        "steps": [{"op": "divide", "params": {"reference": None}},
                  {"op": "mean_center", "params": {"start": 220.0, "end": 300.0}}],
        "measurement": {"kind": "amplitude", "params": {"w1": 260.0}}},
    "Successive derivative ratio (ternary)": {
        "steps": [{"op": "divide", "params": {"reference": None}},
                  {"op": "derivative", "params": {"order": 1, "delta_lambda": 4.0}},
                  {"op": "divide", "params": {"reference": None}},
                  {"op": "derivative", "params": {"order": 1, "delta_lambda": 4.0}}],
        "measurement": {"kind": "amplitude", "params": {"w1": 260.0}}},
    "Double divisor ratio derivative (ternary)": {
        "steps": [{"op": "divide_sum", "params": {"reference": None, "reference2": None}},
                  {"op": "derivative", "params": {"order": 1, "delta_lambda": 4.0}}],
        "measurement": {"kind": "amplitude", "params": {"w1": 260.0}}},
    "Dual amplitude difference (ternary)": {
        "steps": [{"op": "divide", "params": {"reference": None}}],
        "measurement": {"kind": "weighted_difference",
                        "params": {"w1": 240.0, "w2": 260.0, "reference": None}}},
    "Extended ratio subtraction (Y)": {
        "steps": [{"op": "extended_ratio_subtraction",
                   "params": {"divisor": None, "reference": None, "start": 300.0,
                              "end": 320.0}}],
        "measurement": {"kind": "amplitude", "params": {"w1": 275.0}}},
    "Successive spectrum subtraction (ternary)": {
        "steps": [{"op": "spectrum_subtraction",
                   "params": {"reference": None, "wavelength": 340.0, "derivative_order": 0}},
                  {"op": "ratio_subtraction",
                   "params": {"divisor": None, "start": 290.0, "end": 310.0}}],
        "measurement": {"kind": "amplitude", "params": {"w1": 245.0}}},
    "Ratio subtraction (X)": {
        "steps": [{"op": "ratio_subtraction", "params": {"divisor": None, "start": 300.0,
                                                          "end": 320.0}}],
        "measurement": {"kind": "amplitude", "params": {"w1": 260.0}}},
    "Constant multiplication (Y)": {
        "steps": [{"op": "constant_multiplication",
                   "params": {"divisor": None, "start": 300.0, "end": 320.0}}],
        "measurement": {"kind": "amplitude", "params": {"w1": 260.0}}},
    "Ratio plateau / constant value (Y)": {
        "steps": [{"op": "divide", "params": {"reference": None}}],
        "measurement": {"kind": "mean", "params": {"w1": 300.0, "w2": 320.0}}},
    "Spectrum subtraction": {
        "steps": [{"op": "spectrum_subtraction",
                   "params": {"reference": None, "wavelength": 300.0,
                              "derivative_order": 0}}],
        "measurement": {"kind": "amplitude", "params": {"w1": 260.0}}},
    "Factorized zero-order (FZM)": {
        "steps": [{"op": "factorized_recovery",
                   "params": {"reference": None, "wavelength": 260.0,
                              "derivative_order": 1, "delta_lambda": 4.0}}],
        "measurement": {"kind": "amplitude", "params": {"w1": 260.0}}},
    "Constant center": {
        "steps": [{"op": "constant_center",
                   "params": {"divisor": None, "reference": None, "w1": 230.0,
                              "w2": 260.0, "target": "X"}}],
        "measurement": {"kind": "amplitude", "params": {"w1": 260.0}}},
    "Spectrum subtraction–constant multiplication (SS-CM, Y)": {
        "steps": [{"op": "constant_multiplication",
                   "params": {"divisor": None, "start": 300.0, "end": 320.0}}],
        "measurement": {"kind": "amplitude", "params": {"w1": 290.0}}},
    "Derivative subtraction (DS, X in D1)": {
        "steps": [_D1, {"op": "ratio_subtraction",
                        "params": {"divisor": None, "start": 300.0, "end": 320.0,
                                   "divisor_derivative": 1, "divisor_delta_lambda": 4.0}}],
        "measurement": {"kind": "amplitude", "params": {"w1": 245.0}}},
    "Derivative subtraction–constant multiplication (DS-CM, Y in D1)": {
        "steps": [_D1, {"op": "constant_multiplication",
                        "params": {"divisor": None, "start": 300.0, "end": 320.0,
                                   "divisor_derivative": 1, "divisor_delta_lambda": 4.0}}],
        "measurement": {"kind": "amplitude", "params": {"w1": 320.0}}},
    "Derivative ratio of D1 spectra (D1 DR)": {
        "steps": [_D1, {"op": "divide", "params": {"reference": None, "divisor_derivative": 1,
                                                   "divisor_delta_lambda": 4.0}},
                  _D1],
        "measurement": {"kind": "amplitude", "params": {"w1": 270.0}}},
    "Successive ratio subtraction (ternary SRS, Z)": {
        "steps": [{"op": "ratio_subtraction",
                   "params": {"divisor": None, "start": 315.0, "end": 365.0}},
                  {"op": "ratio_subtraction",
                   "params": {"divisor": None, "start": 275.0, "end": 305.0}}],
        "measurement": {"kind": "amplitude", "params": {"w1": 240.0}}},
    "Mean centering of ratio spectra (ternary)": {
        "steps": [{"op": "divide", "params": {"reference": None}},
                  {"op": "mean_center", "params": {"start": 215.0, "end": 275.0}},
                  {"op": "divide_centered_ratio",
                   "params": {"reference": None, "divisor": None, "start": 215.0,
                              "end": 275.0}},
                  {"op": "mean_center", "params": {"start": 215.0, "end": 275.0}}],
        "measurement": {"kind": "amplitude", "params": {"w1": 250.0}}},
    "Ratio difference with normalized divisor": {
        "steps": [{"op": "divide", "params": {"reference": None}}],
        "measurement": {"kind": "difference", "params": {"w1": 264.0, "w2": 225.0}}},
    "Derivative transformation (DT, recover zero order)": {
        "steps": [{"op": "factorized_recovery",
                   "params": {"reference": None, "wavelength": 306.0, "wavelength_end": 313.0,
                              "derivative_order": 1, "delta_lambda": 4.0}}],
        "measurement": {"kind": "amplitude", "params": {"w1": 294.0}}},
    "Concentration value (unit divisor plateau, no regression)": {
        "steps": [{"op": "divide", "params": {"reference": None}}],
        "measurement": {"kind": "mean", "params": {"w1": 330.0, "w2": 333.0}},
        "direct": True},
}


# --------------------------------------------------------------------------- #
# Finders
# --------------------------------------------------------------------------- #
def zero_crossings(s: Spectrum, start: float | None = None,
                   end: float | None = None) -> list[float]:
    """Wavelengths where the signal changes sign (linear interpolation)."""
    x, y = (s.region(start, end) if start is not None and end is not None
            else (s.wavelengths, s.values))
    out = []
    for i in range(len(y) - 1):
        if y[i] == 0:
            out.append(float(x[i]))
        elif y[i] * y[i + 1] < 0:
            out.append(float(x[i] - y[i] * (x[i + 1] - x[i]) / (y[i + 1] - y[i])))
    return out


def extrema(s: Spectrum, prominence: float | None = None) -> dict[str, list[tuple[float, float]]]:
    from scipy.signal import find_peaks

    prom = prominence if prominence is not None else 0.01 * np.ptp(s.values)
    mx, _ = find_peaks(s.values, prominence=prom)
    mn, _ = find_peaks(-s.values, prominence=prom)
    return {"maxima": [(float(s.x[i]), float(s.y[i])) for i in mx],
            "minima": [(float(s.x[i]), float(s.y[i])) for i in mn]}


def isoabsorptive_points(a: Spectrum, b: Spectrum, conc_a: float = 1.0,
                         conc_b: float = 1.0) -> list[float]:
    """Wavelengths where two compounds have equal absorptivity.

    Spectra are normalised by their concentrations first (pass equal
    concentrations, or the actual ones)."""
    diff = a.values / conc_a - _on_grid(b, a) / conc_b
    return zero_crossings(a.with_values(diff))


def find_plateaus(s: Spectrum, min_width: float = 10.0, rel_tol: float = 0.02) -> list[tuple[float, float, float]]:
    """Flat regions of a ratio spectrum: (start, end, mean value)."""
    x, y = s.wavelengths, s.values
    scale = np.ptp(y) or 1.0
    out = []
    i = 0
    while i < len(y):
        j = i
        while j + 1 < len(y) and np.ptp(y[i:j + 2]) <= rel_tol * scale:
            j += 1
        if x[j] - x[i] >= min_width:
            out.append((float(x[i]), float(x[j]), float(y[i:j + 1].mean())))
            i = j + 1
        else:
            i += 1
    return out


# --------------------------------------------------------------------------- #
# Special binary methods
# --------------------------------------------------------------------------- #
def absorptivity(standards: list[Spectrum], compound: str, wavelength: float) -> float:
    """Slope of A(λ) vs concentration through the origin (absorptivity in the
    user's concentration unit, 1 cm path)."""
    c = [s.concentrations[compound] for s in standards]
    a = [s.value_at(wavelength) for s in standards]
    return linear_regression(c, a, through_origin=True).slope


def q_analysis(mixture: Spectrum, iso: float, w2: float, ax_iso: float, ax_2: float,
               ay_2: float) -> dict:
    """Absorbance ratio (Q-analysis) for a binary mixture.

    ``iso`` is the isoabsorptive point (ax_iso = ay_iso), ``w2`` the second
    wavelength (usually λmax of X). a-values are absorptivities.
    """
    a1, a2 = mixture.value_at(iso), mixture.value_at(w2)
    qm, qx, qy = a2 / a1, ax_2 / ax_iso, ay_2 / ax_iso
    cx = (qm - qy) / (qx - qy) * a1 / ax_iso
    cy = (qm - qx) / (qy - qx) * a1 / ax_iso
    return {"X": float(cx), "Y": float(cy), "Qm": qm, "Qx": qx, "Qy": qy}


@dataclass
class AbsorbanceSubtraction:
    """Absorbance subtraction (AS) for binary X+Y with an isoabsorptive point.

    X absorbs alone at ``w2``. The amplitude factor AF = A_X(iso)/A_X(w2) is
    obtained from pure X standards. Total concentration is read from the
    isoabsorptive calibration (either compound's standards)."""

    iso: float
    w2: float
    amplitude_factor: float = float("nan")
    iso_regression: Regression | None = None

    def fit(self, x_standards: list[Spectrum], iso_standards: list[Spectrum]) -> None:
        self.amplitude_factor = float(np.mean(
            [s.value_at(self.iso) / s.value_at(self.w2) for s in x_standards]))
        conc = [sum(s.concentrations.values()) for s in iso_standards]
        self.iso_regression = linear_regression(conc, [s.value_at(self.iso)
                                                       for s in iso_standards])

    def predict(self, mixture: Spectrum) -> dict:
        total_a = mixture.value_at(self.iso)
        x_a = self.amplitude_factor * mixture.value_at(self.w2)
        cx = float(self.iso_regression.predict_x(x_a))
        ctot = float(self.iso_regression.predict_x(total_a))
        return {"X": cx, "Y": ctot - cx, "total": ctot}


@dataclass
class AmplitudeModulation:
    """Amplitude modulation (AM): ratio spectra with divisor Y′; the amplitude
    at the isoabsorptive point gives the total; the plateau (where only Y
    absorbs) gives Y; X = total − Y."""

    iso: float
    plateau: tuple[float, float]
    divisor_conc: float
    total_regression: Regression | None = None

    def _ratio(self, s: Spectrum, divisor: Spectrum) -> Spectrum:
        return s.with_values(s.values / _on_grid(divisor, s))

    def fit(self, standards: list[Spectrum], divisor: Spectrum) -> None:
        conc = [sum(s.concentrations.values()) for s in standards]
        amps = [self._ratio(s, divisor).value_at(self.iso) for s in standards]
        self.total_regression = linear_regression(conc, amps)

    def predict(self, mixture: Spectrum, divisor: Spectrum) -> dict:
        r = self._ratio(mixture, divisor)
        total = float(self.total_regression.predict_x(r.value_at(self.iso)))
        _, y = r.region(*self.plateau)
        cy = float(np.mean(y)) * self.divisor_conc
        return {"X": total - cy, "Y": cy, "total": total}


def hpsam(added: list[float], mixtures: list[Spectrum], w1: float, w2: float) -> dict:
    """H-point standard addition method.

    Standard additions of analyte X to the sample; at w1/w2 the interferent Y
    has equal absorbance. The two addition lines intersect at H(−C_H, A_H):
    C_H is the analyte concentration, A_H the interferent's absorbance."""
    a1 = [m.value_at(w1) for m in mixtures]
    a2 = [m.value_at(w2) for m in mixtures]
    r1, r2 = linear_regression(added, a1), linear_regression(added, a2)
    if abs(r1.slope - r2.slope) < 1e-12:
        raise ValueError("the two addition lines are parallel")
    ch = (r2.intercept - r1.intercept) / (r1.slope - r2.slope)
    ah = r1.intercept + r1.slope * ch
    return {"X": float(-ch), "A_H": float(ah), "line1": r1.to_dict(), "line2": r2.to_dict()}


@dataclass
class InducedAmplitudeModulation:
    """Induced amplitude modulation (IAM) — amplitude modulation without an
    isoabsorptive point.

    The mixture is divided by the *unit-concentration* spectrum of Y (Y′ / C).
    In the ratio spectrum P(λ) = (aX/aY)(λ)·C_X + C_Y. The plateau where only Y
    absorbs gives C_Y; at the chosen λ the factor r = aX/aY (from pure X
    standards) turns the remaining amplitude into C_X:
    C_X = (P(λ) − C_Y) / r. With λ at an isoabsorptive point r = 1 (classic AM).
    """

    wavelength: float
    plateau: tuple[float, float]
    factor: float = float("nan")

    @staticmethod
    def unit_divisor(y_standards: list[Spectrum], compound: str) -> Spectrum:
        """Average unit-concentration spectrum of Y from its standards."""
        ref = y_standards[0]
        rows = [np.interp(ref.wavelengths, s.wavelengths, s.values) / s.concentrations[compound]
                for s in y_standards]
        return ref.copy(values=np.mean(rows, axis=0), name=f"{compound} (unit)")

    def fit(self, x_standards: list[Spectrum], x_name: str, divisor: Spectrum) -> None:
        r = []
        for s in x_standards:
            ratio = s.values / _on_grid(divisor, s)
            r.append(float(np.interp(self.wavelength, s.wavelengths, ratio))
                     / s.concentrations[x_name])
        self.factor = float(np.mean(r))

    def predict(self, mixture: Spectrum, divisor: Spectrum) -> dict:
        ratio = mixture.with_values(mixture.values / _on_grid(divisor, mixture))
        cy = float(np.mean(ratio.region(*self.plateau)[1]))
        cx = (ratio.value_at(self.wavelength) - cy) / self.factor
        return {"X": float(cx), "Y": cy}


@dataclass
class AdvancedAbsorbanceSubtraction:
    """Advanced absorbance subtraction (AAS).

    λ1 and λ2 are chosen where the interferent Y has equal absorbance, so
    ΔA = A(λ2) − A(λ1) depends on X only (dual wavelength) → C_X from the X
    calibration of ΔA. The amplitude factor AF = A_X(λ2)/ΔA_X (pure X) gives
    X's absorbance at λ2 in the mixture; the remainder at λ2 belongs to Y and
    is read from the Y calibration at λ2."""

    w1: float
    w2: float
    amplitude_factor: float = float("nan")
    x_regression: Regression | None = None
    y_regression: Regression | None = None

    def fit(self, x_standards: list[Spectrum], x_name: str, y_standards: list[Spectrum],
            y_name: str) -> None:
        dx = [s.value_at(self.w2) - s.value_at(self.w1) for s in x_standards]
        self.x_regression = linear_regression([s.concentrations[x_name] for s in x_standards], dx)
        self.amplitude_factor = float(np.mean([s.value_at(self.w2) / d
                                               for s, d in zip(x_standards, dx)]))
        self.y_regression = linear_regression([s.concentrations[y_name] for s in y_standards],
                                              [s.value_at(self.w2) for s in y_standards])

    def predict(self, mixture: Spectrum) -> dict:
        d = mixture.value_at(self.w2) - mixture.value_at(self.w1)
        cx = float(self.x_regression.predict_x(d))
        ay = mixture.value_at(self.w2) - self.amplitude_factor * d
        return {"X": cx, "Y": float(self.y_regression.predict_x(ay))}


# --------------------------------------------------------------------------- #
# Progressive (successive) resolution of binary, ternary … mixtures
# --------------------------------------------------------------------------- #
def _ratio(s: Spectrum, divisor: Spectrum) -> Spectrum:
    from spectro.core.operations import op_divide
    return s.with_values(op_divide(s, divisor))


def _conc(s: Spectrum, compound: str) -> float:
    if compound not in s.concentrations:
        raise ValueError(f"no {compound} concentration for '{s.name}'")
    return float(s.concentrations[compound])


def pure_standards(spectra: list[Spectrum], compounds: list[str]) -> dict[str, list[Spectrum]]:
    """Group pure standards by their single non-zero concentration."""
    out: dict[str, list[Spectrum]] = {c: [] for c in compounds}
    for s in spectra:
        nz = [c for c, v in s.concentrations.items() if v]
        if len(nz) == 1 and nz[0] in out:
            out[nz[0]].append(s)
    return out


class _Progressive:
    """Common interface of the progressive methods so they can be saved,
    applied to new spectra, exported and ranked by the optimizer like any other
    method: ``fit_spectra`` / ``predict_spectra`` take plain spectrum lists and
    a resolver for the stored divisor reference."""

    compounds: list[str]

    @property
    def is_fitted(self) -> bool:
        return bool(self.regressions)

    def _divisor(self, resolve: Resolver | None):
        return None

    def fit_spectra(self, spectra: list[Spectrum], resolve: Resolver | None = None) -> dict:
        return self._fit(pure_standards(spectra, self.compounds), self._divisor(resolve))

    def predict_spectra(self, spectra: list[Spectrum], resolve: Resolver | None = None
                        ) -> np.ndarray:
        div = self._divisor(resolve)
        rows = [self._predict(s, div) for s in spectra]
        return np.array([[r[c] for c in self.compounds] for r in rows], float).reshape(
            len(rows), len(self.compounds))


@dataclass
class AmplitudeCentering(_Progressive):
    """Progressive resolution on ONE ratio spectrum at ONE wavelength λc.

    The mixture is divided by one divisor (Z′, usually the normalised spectrum
    of one component). In the ratio spectrum, each component's amplitude at λc
    is obtained in turn:

    * ``plateau`` (optional): the flat region where only the divisor compound
      remains → its constant (amplitude calculation). The constant is
      subtracted before the differences below are measured.
    * ``differences[c] = {"w1", "w2", "factor_from"}``: P(λ1) − F·P(λ2) of the
      ratio spectrum depends on c only (the other components have equal
      amplitudes, F = 1, or F = P(λ1)/P(λ2) of the compound named in
      ``factor_from`` — the amplitude/equality factor). A regression of this
      difference against c's own amplitude at λc (pure c standards) gives
      c's *postulated* amplitude at λc (amplitude difference).
    * ``subtract``: the remaining compound's amplitude = recorded amplitude at
      λc − all the others (amplitude subtraction).

    Each amplitude at λc is converted to concentration with that compound's
    regression at λc (or one *unified* regression when λc is an isoabsorptive
    point). This one class covers advanced amplitude centering (AAC, partial
    and complete overlap), the modified amplitude center method (MACM),
    ratio difference–isoabsorptive (RIDSS) and constant value via amplitude
    difference (CV-AD), and amplitude modulation as the binary special case.
    """

    wavelength: float
    compounds: list[str]
    subtract: str | None = None
    divisor_compound: str | None = None
    plateau: tuple[float, float] | None = None
    differences: dict[str, dict] = field(default_factory=dict)
    unified: bool = False
    factors: dict[str, float] = field(default_factory=dict)
    diff_regressions: dict[str, Regression] = field(default_factory=dict)
    regressions: dict[str, Regression] = field(default_factory=dict)
    divisor: Any = None        # stored reference of the divisor spectrum (saved methods)
    name: str = ""

    def _divisor(self, resolve: Resolver | None) -> Spectrum:
        d = self.divisor
        if d is None:
            raise ValueError("the method has no divisor spectrum")
        if isinstance(d, Spectrum):
            return d
        if resolve is None:
            raise ValueError("no resolver for the divisor spectrum")
        return resolve(d)

    def _fit(self, standards, divisor):
        return self.fit(standards, divisor)

    def _predict(self, s, divisor):
        return self.predict(s, divisor)

    def _check(self) -> None:
        known = set(self.differences) | ({self.subtract} if self.subtract else set())
        if self.plateau is not None:
            if not self.divisor_compound:
                raise ValueError("a plateau needs the divisor compound")
            known.add(self.divisor_compound)
        missing = [c for c in self.compounds if c not in known]
        if missing:
            raise ValueError(f"no way to resolve: {', '.join(missing)} — give a λ pair, the "
                             "plateau or amplitude subtraction for every compound")
        if self.subtract and (self.subtract in self.differences or
                              (self.plateau is not None and self.subtract == self.divisor_compound)):
            raise ValueError(f"{self.subtract} is resolved twice")

    def _difference(self, r: Spectrum, c: str) -> float:
        d = self.differences[c]
        return r.value_at(d["w1"]) - self.factors.get(c, 1.0) * r.value_at(d["w2"])

    def _centred(self, r: Spectrum) -> tuple[Spectrum, float | None]:
        if self.plateau is None:
            return r, None
        const = float(np.mean(r.region(*self.plateau)[1]))
        return r.with_values(r.values - const), const

    def fit(self, standards: dict[str, list[Spectrum]], divisor: Spectrum) -> dict:
        """``standards``: pure standards of each compound (compound → spectra)."""
        self._check()
        for c in self.compounds:
            if len(standards.get(c, [])) < 3:
                raise ValueError(f"check at least three pure {c} standards")
        ratios = {c: [_ratio(s, divisor) for s in standards[c]] for c in self.compounds}
        for c, d in self.differences.items():
            src = d.get("factor_from")
            if src:
                if src not in ratios:
                    raise ValueError(f"factor compound {src} has no standards")
                self.factors[c] = float(np.mean([r.value_at(d["w1"]) / r.value_at(d["w2"])
                                                 for r in ratios[src]]))
            self.diff_regressions[c] = linear_regression(
                [r.value_at(self.wavelength) for r in ratios[c]],
                [self._difference(self._centred(r)[0], c) for r in ratios[c]])
        if self.unified:
            xs, ys = [], []
            for c in self.compounds:
                xs += [_conc(s, c) for s in standards[c]]
                ys += [r.value_at(self.wavelength) for r in ratios[c]]
            reg = linear_regression(xs, ys)
            self.regressions = {c: reg for c in self.compounds}
        else:
            self.regressions = {c: linear_regression([_conc(s, c) for s in standards[c]],
                                                     [r.value_at(self.wavelength)
                                                      for r in ratios[c]])
                                for c in self.compounds}
        return {"factors": dict(self.factors),
                "difference_r": {c: r.r for c, r in self.diff_regressions.items()},
                "r": {c: r.r for c, r in self.regressions.items()}}

    def amplitudes(self, mixture: Spectrum, divisor: Spectrum) -> dict[str, float]:
        r = _ratio(mixture, divisor)
        centred, const = self._centred(r)
        amp: dict[str, float] = {}
        if const is not None:
            amp[self.divisor_compound] = const
        for c in self.differences:
            amp[c] = float(self.diff_regressions[c].predict_x(self._difference(centred, c)))
        if self.subtract:
            amp[self.subtract] = r.value_at(self.wavelength) - sum(amp.values())
        return amp

    def predict(self, mixture: Spectrum, divisor: Spectrum) -> dict[str, float]:
        if not self.regressions:
            raise RuntimeError("not fitted")
        amp = self.amplitudes(mixture, divisor)
        return {c: float(self.regressions[c].predict_x(amp[c])) for c in self.compounds}

    def to_dict(self) -> dict:
        return {"type": "amplitude_centering", "wavelength": self.wavelength,
                "compounds": self.compounds, "subtract": self.subtract,
                "divisor_compound": self.divisor_compound,
                "plateau": list(self.plateau) if self.plateau else None,
                "differences": self.differences, "unified": self.unified,
                "factors": self.factors,
                "diff_regressions": {c: r.to_dict() for c, r in self.diff_regressions.items()},
                "regressions": {c: r.to_dict() for c, r in self.regressions.items()},
                "divisor": None if isinstance(self.divisor, Spectrum) else self.divisor,
                "name": self.name}

    @classmethod
    def from_dict(cls, d: dict) -> "AmplitudeCentering":
        return cls(d["wavelength"], d["compounds"], d.get("subtract"), d.get("divisor_compound"),
                   tuple(d["plateau"]) if d.get("plateau") else None,
                   copy.deepcopy(d.get("differences", {})),
                   d.get("unified", False), dict(d.get("factors", {})),
                   {c: Regression.from_dict(r) for c, r in d.get("diff_regressions", {}).items()},
                   {c: Regression.from_dict(r) for c, r in d.get("regressions", {}).items()},
                   d.get("divisor"), d.get("name", ""))

    def describe(self) -> str:
        parts = [f"÷ {self.divisor_compound or 'divisor'}, λc {self.wavelength:.1f} nm"]
        if self.plateau:
            parts.append(f"{self.divisor_compound} from plateau "
                         f"{self.plateau[0]:.0f}–{self.plateau[1]:.0f} nm")
        for c, d in self.differences.items():
            f = f" − F({d['factor_from']})·" if d.get("factor_from") else " − "
            parts.append(f"{c}: P{d['w1']:.1f}{f}P{d['w2']:.1f}")
        if self.subtract:
            parts.append(f"{self.subtract} by subtraction")
        return "; ".join(parts)


@dataclass
class AbsorptionFactorMethod(_Progressive):
    """Successive absorption (amplitude) factor method — MAFM for ternary
    mixtures and the absorption factor method for binary ones.

    ``order`` lists ``(compound, λ)`` so that the first compound absorbs alone
    at its λ, the second is overlapped only by the first at its λ, and so on.
    For each compound k, F_k(λ) = A_k(λ)/A_k(λ_k) comes from its pure
    standards. In the mixture, R_i = A(λ_i) − Σ_{k<i} F_k(λ_i)·R_k is the
    absorbance of compound i alone at λ_i; its postulated absorbance at the
    quantitation wavelength (``quant[i]``, default λ_i) is F_i(λq)·R_i, read
    from compound i's calibration at λq.
    """

    order: list[tuple[str, float]]
    quant: dict[str, float] = field(default_factory=dict)
    factors: dict[str, dict[str, float]] = field(default_factory=dict)
    regressions: dict[str, Regression] = field(default_factory=dict)
    name: str = ""

    @property
    def compounds(self) -> list[str]:
        return [c for c, _ in self.order]

    def _fit(self, standards, divisor):
        return self.fit(standards)

    def _predict(self, s, divisor):
        return self.predict(s)

    def describe(self) -> str:
        return " → ".join(f"{c} at {w:.1f} nm" + (f" (read at {self.quant[c]:.1f})"
                                                   if self.quant.get(c) else "")
                          for c, w in self.order)

    def _lam(self, c: str) -> float:
        return self.quant.get(c) or dict(self.order)[c]

    def fit(self, standards: dict[str, list[Spectrum]]) -> dict:
        names = [c for c, _ in self.order]
        if len(set(names)) != len(names):
            raise ValueError("each compound may appear once")
        lams = [w for _, w in self.order]
        for i, (c, w) in enumerate(self.order):
            st = standards.get(c, [])
            if len(st) < 3:
                raise ValueError(f"check at least three pure {c} standards")
            need = set(lams[i + 1:]) | {self._lam(c)}
            self.factors[c] = {f"{lam:g}": float(np.mean([s.value_at(lam) / s.value_at(w)
                                                          for s in st])) for lam in need}
            self.regressions[c] = linear_regression([_conc(s, c) for s in st],
                                                    [s.value_at(self._lam(c)) for s in st])
        return {"factors": self.factors, "r": {c: r.r for c, r in self.regressions.items()}}

    def predict(self, mixture: Spectrum) -> dict[str, float]:
        if not self.regressions:
            raise RuntimeError("not fitted")
        resid: dict[str, float] = {}
        out: dict[str, float] = {}
        for c, w in self.order:
            a = mixture.value_at(w) - sum(self.factors[k][f"{w:g}"] * r for k, r in resid.items())
            resid[c] = a
            out[c] = float(self.regressions[c].predict_x(self.factors[c][f"{self._lam(c):g}"] * a))
        return out

    def to_dict(self) -> dict:
        return {"type": "absorption_factor", "order": [list(o) for o in self.order],
                "compounds": self.compounds, "quant": self.quant, "factors": self.factors,
                "regressions": {c: r.to_dict() for c, r in self.regressions.items()},
                "name": self.name}

    @classmethod
    def from_dict(cls, d: dict) -> "AbsorptionFactorMethod":
        return cls([(c, float(w)) for c, w in d["order"]], dict(d.get("quant", {})),
                   copy.deepcopy(d.get("factors", {})),
                   {c: Regression.from_dict(r) for c, r in d.get("regressions", {}).items()},
                   d.get("name", ""))


PROGRESSIVE_TYPES = {"amplitude_centering": AmplitudeCentering,
                     "absorption_factor": AbsorptionFactorMethod}


def progressive_from_dict(d: dict):
    """A saved amplitude-centering or absorption-factor method."""
    cls = PROGRESSIVE_TYPES.get(d.get("type"))
    if cls is None:
        raise ValueError(f"not a progressive method: {d.get('type')}")
    return cls.from_dict(d)


def enrichment_correction(found: float, added: float, claimed: float | None = None) -> dict:
    """Sample enrichment (spiking or spectrum addition) of a minor component:
    the amount added is subtracted from the amount found; % of the claimed
    amount if given."""
    net = float(found) - float(added)
    out = {"found_total": float(found), "added": float(added), "found": net}
    if claimed:
        out["percent_of_claimed"] = 100 * net / claimed
    return out


def equal_amplitude_wavelengths(s: Spectrum, wavelength: float, start: float | None = None,
                                end: float | None = None) -> list[float]:
    """Wavelengths where ``s`` has the same value as at ``wavelength`` (for
    choosing λ pairs in dual wavelength / dual amplitude difference / HPSAM)."""
    target = s.value_at(wavelength)
    out = zero_crossings(s.with_values(s.values - target), start, end)
    return [w for w in out if abs(w - wavelength) > 2 * s.step]


# --------------------------------------------------------------------------- #
# Standard addition and robustness
# --------------------------------------------------------------------------- #
def standard_addition_recovery(predict, spectra: list[Spectrum], added: list[float]) -> dict:
    """Standard addition technique.

    ``predict`` maps a spectrum to a concentration (any calibrated method).
    Spectra with added = 0 are the unspiked sample; for each spiked spectrum the
    amount recovered = found − mean(found of unspiked). Also returns the
    classic extrapolation (x-intercept of found vs added)."""
    found = np.array([predict(s) for s in spectra], float)
    added = np.array(added, float)
    base = found[added == 0]
    if base.size == 0:
        raise ValueError("include at least one unspiked sample (added = 0)")
    b = float(base.mean())
    rows = []
    for s, f, a in zip(spectra, found, added):
        rec = 100 * (f - b) / a if a else None
        rows.append({"name": s.name, "added": float(a), "found": float(f),
                     "recovered": float(f - b) if a else None, "recovery": rec})
    recs = [r["recovery"] for r in rows if r["recovery"] is not None]
    out = {"sample": b, "rows": rows,
           "mean_recovery": float(np.mean(recs)) if recs else float("nan"),
           "sd_recovery": float(np.std(recs, ddof=1)) if len(recs) > 1 else float("nan")}
    if np.unique(added).size >= 2:
        reg = linear_regression(added, found)
        out["extrapolated"] = reg.intercept / reg.slope if reg.slope else float("nan")
        out["regression"] = reg.to_dict()
    return out


# derivative orders are a choice of method, not a parameter to vary by ±1
_NOT_VARIED = {"order", "derivative_order", "divisor_derivative"}


def numeric_parameters(method: UnivariateMethod) -> list[dict]:
    """Parameters of a method that can be varied in a robustness study."""
    from spectro.core.operations import REGISTRY

    out = []
    for i, st in enumerate(method.steps):
        op = REGISTRY[st["op"]]
        for p in op.params:
            v = st.get("params", {}).get(p.name, p.default)
            if p.name in _NOT_VARIED:
                continue
            if p.kind in ("wavelength", "float", "int") and v is not None:
                delta = 1.0 if p.kind in ("wavelength", "int") else abs(v) * 0.05 or 0.1
                out.append({"path": ("steps", i, p.name), "label": f"{op.label}: {p.label}",
                            "value": v, "delta": delta, "kind": p.kind})
    for k, v in method.measurement.get("params", {}).items():
        if k in ("w1", "w2") and v is not None:
            out.append({"path": ("measurement", k), "label": f"Measurement {k.replace('w', 'λ')}",
                        "value": v, "delta": 1.0, "kind": "wavelength"})
    return out


def _with(method: UnivariateMethod, path: tuple, value) -> UnivariateMethod:
    m = UnivariateMethod.from_dict(copy.deepcopy(method.to_dict()))
    m.regression = None
    if path[0] == "steps":
        m.steps[path[1]]["params"][path[2]] = value
    else:
        m.measurement["params"][path[1]] = value
    return m


def robustness_study(method: UnivariateMethod, calibration: list[Spectrum],
                     test: list[Spectrum], variations: list[dict],
                     resolve: Resolver | None = None,
                     progress=None) -> dict:
    """Re-calibrate and re-assay with each parameter moved by ±delta.

    Returns per-variant slope, r and mean recovery / found values, and the
    %RSD of the results across all variants (robust if small)."""
    def run(m: UnivariateMethod) -> dict:
        reg = m.calibrate(calibration, resolve)
        found = [m.predict(s, resolve) for s in test]
        rec = [100 * f / s.concentrations[m.compound] for f, s in zip(found, test)
               if s.concentrations.get(m.compound)]
        return {"slope": reg.slope, "r": reg.r, "found": found,
                "mean_found": float(np.mean(found)),
                "mean_recovery": float(np.mean(rec)) if rec else float("nan"),
                "rsd_recovery": float(np.std(rec, ddof=1) / np.mean(rec) * 100)
                if len(rec) > 1 else float("nan")}

    variants = [{"label": "Nominal", "change": "", **run(method)}]
    total = 2 * len(variations)
    for k, v in enumerate(variations):
        for sign in (-1, 1):
            new = v["value"] + sign * v["delta"]
            if v["kind"] == "int":
                new = int(round(new))
            m = _with(method, tuple(v["path"]), new)
            try:
                res = run(m)
            except Exception as exc:  # reported, not fatal
                res = {"error": str(exc)}
            variants.append({"label": v["label"], "change": f"{v['value']:g} → {new:g}", **res})
            if progress and progress(2 * k + (sign > 0) + 1, total) is False:
                raise InterruptedError("cancelled")
    ok = [v for v in variants if "error" not in v]
    means = [v["mean_found"] for v in ok]
    nominal = variants[0]["mean_found"]
    for v in ok:
        v["deviation"] = 100 * (v["mean_found"] - nominal) / nominal if nominal else float("nan")
    return {"variants": variants,
            "rsd_across_variants": float(np.std(means, ddof=1) / np.mean(means) * 100)
            if len(means) > 1 else float("nan"),
            "max_abs_deviation": float(max(abs(v["deviation"]) for v in ok))}
