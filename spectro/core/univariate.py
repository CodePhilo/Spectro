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

    def signal(self, s: Spectrum, resolve: Resolver | None = None) -> float:
        return measure(apply_pipeline(s, self.steps, resolve), self.measurement, resolve)

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
        if self.regression is None:
            raise RuntimeError("method is not calibrated")
        return float(self.regression.predict_x(self.signal(s, resolve)))

    def to_dict(self) -> dict:
        return {"type": "univariate", "name": self.name, "compound": self.compound,
                "steps": self.steps, "measurement": self.measurement,
                "through_origin": self.through_origin,
                "regression": self.regression.to_dict() if self.regression else None}

    @classmethod
    def from_dict(cls, d: dict) -> "UnivariateMethod":
        reg = d.get("regression")
        return cls(d["name"], d["compound"], d.get("steps", []), d["measurement"],
                   d.get("through_origin", False),
                   Regression.from_dict(reg) if reg else None)


# --------------------------------------------------------------------------- #
# Method templates (pre-filled pipelines for the UI wizard)
# --------------------------------------------------------------------------- #
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
        "steps": [{"op": "divide", "params": {"reference": None}},
                  {"op": "derivative", "params": {"order": 1, "delta_lambda": 4.0}}],
        "measurement": {"kind": "amplitude", "params": {"w1": 260.0}}},
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
