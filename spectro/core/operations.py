"""Spectral processing operations and pipelines.

Every operation is registered with a declarative parameter list so the UI can
build forms automatically and the audit trail can record the exact parameters.
A *pipeline* is a JSON-serialisable list of steps::

    [{"op": "smooth_sg", "params": {"window": 11, "polyorder": 2}},
     {"op": "divide", "params": {"reference": 17}},
     {"op": "derivative", "params": {"order": 1, "delta_lambda": 4}}]

Parameters of kind ``spectrum`` hold a reference (e.g. a database id) that is
turned into a :class:`Spectrum` by a *resolver* callable.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from functools import lru_cache
from typing import Any, Callable

import numpy as np
from scipy import sparse

from scipy.sparse.linalg import spsolve

from spectro.core.spectrum import Spectrum

Resolver = Callable[[Any], Spectrum]


@dataclass
class Param:
    name: str
    label: str
    kind: str  # float | int | bool | choice | spectrum | wavelength
    default: Any = None
    choices: tuple = ()
    minimum: float | None = None
    maximum: float | None = None
    help: str = ""
    optional: bool = False


@dataclass
class Operation:
    key: str
    label: str
    category: str
    func: Callable[..., np.ndarray | Spectrum]
    params: list[Param] = field(default_factory=list)
    description: str = ""
    changes_grid: bool = False

    def defaults(self) -> dict[str, Any]:
        return {p.name: p.default for p in self.params}


REGISTRY: dict[str, Operation] = {}


def register(key: str, label: str, category: str, params: list[Param] | None = None,
             description: str = "", changes_grid: bool = False):
    def deco(func):
        REGISTRY[key] = Operation(key, label, category, func, params or [],
                                  description or (func.__doc__ or "").strip(), changes_grid)
        return func
    return deco


def _wl(name: str, label: str, default: float | None = None, **kw) -> Param:
    return Param(name, label, "wavelength", default, **kw)


def _ref(name: str = "reference", label: str = "Reference spectrum", **kw) -> Param:
    return Param(name, label, "spectrum", None, **kw)


def _on_grid(ref: Spectrum, target: Spectrum) -> np.ndarray:
    if ref.same_grid(target):
        return ref.values
    lo, hi = ref.wavelengths[0], ref.wavelengths[-1]
    if target.wavelengths[0] < lo - 1e-6 or target.wavelengths[-1] > hi + 1e-6:
        raise ValueError(
            f"reference '{ref.name}' ({lo:g}–{hi:g} nm) does not cover "
            f"{target.wavelengths[0]:g}–{target.wavelengths[-1]:g} nm; crop first")
    return np.interp(target.wavelengths, ref.wavelengths, ref.values)


def _region_mean(s: Spectrum, start: float, end: float) -> float:
    return float(np.mean(s.region(start, end)[1]))


@lru_cache(maxsize=256)
def _savgol_matrices(window: int, polyorder: int, deriv: int):
    """Convolution coefficients and the edge projection matrices of a
    Savitzky–Golay filter (unit spacing). The edges reproduce SciPy's
    ``mode='interp'`` exactly: a polynomial fitted to the first / last
    ``window`` points, evaluated (or differentiated) at the edge points —
    linear in the data, so it is one matrix product."""
    from math import factorial

    from scipy.signal import savgol_coeffs

    coeffs = savgol_coeffs(window, polyorder, deriv=deriv, use="conv")
    half = window // 2
    t = np.arange(window, dtype=float) - (window - 1) / 2   # centred: well conditioned
    vander = np.vander(t, polyorder + 1, increasing=True)
    pinv = np.linalg.pinv(vander)

    def evaluate(points):
        d = np.zeros((len(points), polyorder + 1))
        for k in range(deriv, polyorder + 1):
            d[:, k] = factorial(k) / factorial(k - deriv) * points ** (k - deriv)
        return d @ pinv
    left = evaluate(t[:half])
    right = evaluate(t[window - half:])
    for a in (coeffs, left, right):
        a.setflags(write=False)
    return coeffs, left, right


def savgol(y: np.ndarray, window: int, polyorder: int, deriv: int = 0,
           delta: float = 1.0) -> np.ndarray:
    """Savitzky–Golay filter, identical to ``scipy.signal.savgol_filter(...,
    mode='interp')`` but with cached coefficients (≈ 10× faster for the many
    short spectra of the optimizer)."""
    y = np.asarray(y, dtype=float)
    coeffs, left, right = _savgol_matrices(int(window), int(polyorder), int(deriv))
    out = np.convolve(y, coeffs, mode="same")
    half = window // 2
    if half:
        out[:half] = left @ y[:window]
        out[-half:] = right @ y[-window:]
    if deriv:
        out /= delta ** deriv
    return out


def savgol_derivative(s: Spectrum, order: int, window: int, polyorder: int) -> np.ndarray:
    window = int(window) | 1
    if polyorder > 10:
        raise ValueError("Savitzky–Golay polynomial orders above 10 are numerically unstable "
                         "(2–4 is usual for UV spectra)")
    if window <= polyorder:
        raise ValueError("window must be larger than the polynomial order")
    if polyorder < order:
        raise ValueError("polynomial order must be ≥ derivative order")
    if window > s.values.size:
        raise ValueError(f"window ({window} points) is longer than the spectrum "
                         f"({s.values.size} points)")
    x = s.wavelengths
    d = np.diff(x)
    if np.ptp(d) <= 1e-6 * np.median(d):
        return savgol(s.values, window, polyorder, deriv=order, delta=float(np.median(d)))
    # Savitzky–Golay assumes equal spacing: work on a regular grid at the finest
    # step, then return to the original wavelengths.
    step = float(np.min(d))
    reg = np.arange(x[0], x[-1] + step / 2, step)
    y = savgol(np.interp(reg, x, s.values), window, polyorder, deriv=order, delta=step)
    return np.interp(x, reg, y)


def difference_derivative(s: Spectrum, order: int, delta_lambda: float) -> np.ndarray:
    """Symmetric finite-difference derivative with interval Δλ (nm), as in
    instrument software: dA/dλ ≈ [A(λ+Δλ/2) − A(λ−Δλ/2)] / Δλ, applied n times."""
    x, y = s.wavelengths, s.values.copy()
    h = float(delta_lambda) / 2.0
    for _ in range(order):
        y = (np.interp(x + h, x, y) - np.interp(x - h, x, y)) / (2 * h)
        edge = (x - h < x[0]) | (x + h > x[-1])
        y[edge] = np.nan
    good = np.isfinite(y)
    if good.sum() >= 2:
        y[~good] = np.interp(x[~good], x[good], y[good])
    return y


# --------------------------------------------------------------------------- #
# Range / grid
# --------------------------------------------------------------------------- #
@register("crop", "Crop wavelength range", "Range",
          [_wl("start", "From (nm)", 200.0), _wl("end", "To (nm)", 400.0)],
          changes_grid=True)
def op_crop(s: Spectrum, start: float, end: float) -> Spectrum:
    """Keep only data between two wavelengths."""
    x, y = s.region(start, end)
    return s.copy(wavelengths=x, values=y)


@register("resample", "Resample to regular interval", "Range",
          [Param("step", "Interval (nm)", "float", 1.0, minimum=0.001),
           _wl("start", "From (nm)", None, optional=True),
           _wl("end", "To (nm)", None, optional=True)], changes_grid=True)
def op_resample(s: Spectrum, step: float, start: float | None = None,
                end: float | None = None) -> Spectrum:
    """Linear interpolation onto a regular grid."""
    lo = s.wavelengths[0] if start is None else max(start, s.wavelengths[0])
    hi = s.wavelengths[-1] if end is None else min(end, s.wavelengths[-1])
    if step <= 0:
        raise ValueError("the resampling interval must be positive")
    n = int(np.floor((hi - lo) / step + 1e-9)) + 1
    if hi <= lo or n < 2:
        raise ValueError(f"resampling range {lo:g}–{hi:g} nm with interval {step:g} nm "
                         "gives fewer than two points")
    return s.resampled(lo + step * np.arange(n))


# --------------------------------------------------------------------------- #
# Smoothing
# --------------------------------------------------------------------------- #
@register("smooth_sg", "Savitzky–Golay smoothing", "Smoothing",
          [Param("window", "Window (points, odd)", "int", 11, minimum=3),
           Param("polyorder", "Polynomial order", "int", 2, minimum=0, maximum=10)])
def op_smooth_sg(s, window: int, polyorder: int):
    """Least-squares polynomial smoothing."""
    return savgol_derivative(s, 0, window, polyorder)


@register("smooth_ma", "Moving-average smoothing", "Smoothing",
          [Param("window", "Window (points)", "int", 5, minimum=1)])
def op_smooth_ma(s, window: int):
    """Centred moving average (edges use a shrinking window)."""
    w = max(1, int(window))
    kernel = np.ones(w)
    num = np.convolve(s.values, kernel, mode="same")
    den = np.convolve(np.ones_like(s.values), kernel, mode="same")
    return num / den


@register("smooth_whittaker", "Whittaker smoothing", "Smoothing",
          [Param("lam", "Smoothness λ", "float", 10.0, minimum=0.0)])
def op_smooth_whittaker(s, lam: float):
    """Penalised least squares smoother (2nd-order differences)."""
    n = s.values.size
    d = sparse.diags([1.0, -2.0, 1.0], [0, 1, 2], shape=(n - 2, n))
    a = sparse.eye(n) + lam * (d.T @ d)
    return spsolve(a.tocsc(), s.values)


# --------------------------------------------------------------------------- #
# Baseline
# --------------------------------------------------------------------------- #
@register("baseline_offset", "Baseline: subtract offset", "Baseline",
          [_wl("start", "Reference region from (nm)", 380.0),
           _wl("end", "to (nm)", 400.0)])
def op_baseline_offset(s, start: float, end: float):
    """Subtract the mean value of a non-absorbing region."""
    return s.values - _region_mean(s, start, end)


@register("baseline_linear", "Baseline: two-point linear", "Baseline",
          [_wl("w1", "Wavelength 1 (nm)", 220.0), _wl("w2", "Wavelength 2 (nm)", 380.0)])
def op_baseline_linear(s, w1: float, w2: float):
    """Subtract a straight line through two baseline points."""
    if abs(w2 - w1) < 1e-9:
        raise ValueError("choose two different baseline wavelengths")
    a1, a2 = s.value_at(w1), s.value_at(w2)
    return s.values - (a1 + (a2 - a1) * (s.wavelengths - w1) / (w2 - w1))


@register("baseline_poly", "Baseline: iterative polynomial", "Baseline",
          [Param("order", "Polynomial order", "int", 2, minimum=0, maximum=8),
           Param("iterations", "Iterations", "int", 100, minimum=1)])
def op_baseline_poly(s, order: int, iterations: int):
    """Modified polynomial fit (Lieber & Mahadevan-Jansen) of the baseline."""
    x = (s.wavelengths - s.wavelengths.mean()) / np.ptp(s.wavelengths)
    y = s.values.copy()
    for _ in range(int(iterations)):
        base = np.polyval(np.polyfit(x, y, int(order)), x)
        new = np.minimum(y, base)
        if np.allclose(new, y):
            break
        y = new
    return s.values - np.polyval(np.polyfit(x, y, int(order)), x)


@register("baseline_als", "Baseline: asymmetric least squares", "Baseline",
          [Param("lam", "Smoothness λ", "float", 1e5, minimum=0.0),
           Param("p", "Asymmetry p", "float", 0.01, minimum=0.0, maximum=1.0),
           Param("iterations", "Iterations", "int", 10, minimum=1)])
def op_baseline_als(s, lam: float, p: float, iterations: int):
    """Eilers–Boelens asymmetric least squares baseline removal."""
    if not 0.0 < p < 1.0:
        raise ValueError("ALS baseline: the asymmetry p must lie strictly between 0 and 1")
    y = s.values
    n = y.size
    d = sparse.diags([1.0, -2.0, 1.0], [0, 1, 2], shape=(n - 2, n))
    dd = lam * (d.T @ d)
    w = np.ones(n)
    z = y
    for _ in range(int(iterations)):
        wm = sparse.diags(w)
        z = spsolve((wm + dd).tocsc(), w * y)
        w = np.where(y > z, p, 1 - p)    # ties keep a weight (no singular system)
    return y - z


# --------------------------------------------------------------------------- #
# Arithmetic
# --------------------------------------------------------------------------- #
@register("subtract", "Subtract spectrum (blank / component)", "Arithmetic",
          [_ref(), Param("factor", "Multiply reference by", "float", 1.0)])
def op_subtract(s, reference: Spectrum, factor: float = 1.0):
    """A − k·B, e.g. blank or solvent correction."""
    return s.values - factor * _on_grid(reference, s)


@register("add", "Add spectrum", "Arithmetic",
          [_ref(), Param("factor", "Multiply reference by", "float", 1.0)])
def op_add(s, reference: Spectrum, factor: float = 1.0):
    """A + k·B (e.g. build a double-divisor / sum spectrum)."""
    return s.values + factor * _on_grid(reference, s)


_DIVISOR_DERIVATIVE = [
    Param("divisor_derivative", "Differentiate divisor first (order, 0 = no)", "int", 0,
          minimum=0, maximum=4,
          help="Derivative-mode ratio methods (D1 DR, derivative subtraction, DS-CM): "
               "the divisor is differentiated (difference method, Δλ below) before use. "
               "Differentiate the sample with an earlier Derivative step."),
    Param("divisor_delta_lambda", "Divisor Δλ (nm)", "float", 4.0, minimum=0.001),
]


def _divisor_values(s: Spectrum, divisor: Spectrum, order: int = 0,
                    delta_lambda: float = 4.0) -> np.ndarray:
    """Divisor on the grid of ``s``, optionally as its derivative of ``order``."""
    if int(order) > 0:
        divisor = op_derivative(divisor, int(order), "difference", delta_lambda)
    return _on_grid(divisor, s)


@register("divide", "Divide by spectrum (ratio spectrum)", "Ratio spectra",
          [_ref("reference", "Divisor spectrum"),
           Param("threshold", "Ignore divisor below", "float", 1e-4, minimum=0.0,
                 help="Points where |divisor| is below this are interpolated."),
           *_DIVISOR_DERIVATIVE])
def op_divide(s, reference: Spectrum, threshold: float = 1e-4, divisor_derivative: int = 0,
              divisor_delta_lambda: float = 4.0):
    """Ratio spectrum: mixture ÷ divisor (or ÷ derivative of the divisor)."""
    div = _divisor_values(s, reference, divisor_derivative, divisor_delta_lambda)
    bad = (np.abs(div) < threshold) | (div == 0)
    if bad.all():
        raise ValueError("divisor is zero over the whole range")
    out = np.empty_like(s.values)
    out[~bad] = s.values[~bad] / div[~bad]
    if bad.any():
        out[bad] = np.interp(s.wavelengths[bad], s.wavelengths[~bad], out[~bad])
    return out


@register("divide_centered_ratio", "Divide by mean-centred ratio (ternary MCR)", "Ratio spectra",
          [_ref("reference", "Second component Y′"), _ref("divisor", "First divisor Z′"),
           _wl("start", "Mean-centring from (nm)", 220.0), _wl("end", "to (nm)", 300.0),
           Param("threshold", "Ignore divisor below", "float", 1e-4, minimum=0.0)],
          changes_grid=True)
def op_divide_centered_ratio(s, reference: Spectrum, divisor: Spectrum, start: float, end: float,
                             threshold: float = 1e-4) -> Spectrum:
    """Mean centering of ratio spectra for ternary mixtures (Afkhami): the
    mean-centred ratio spectrum MC(mixture/Z′) is divided by MC(Y′/Z′) over the
    same range; a following *Mean centering* step removes the Y constant and
    leaves a signal proportional to X only."""
    s = op_crop(s, start, end)
    y = op_crop(reference, start, end)
    ratio = op_divide(y, divisor)
    centred = ratio - ratio.mean()
    if np.max(np.abs(centred)) <= 1e-9 * max(np.max(np.abs(ratio)), 1e-300):
        raise ValueError("Y′/Z′ is constant (proportional spectra) — choose a different "
                         "second component or divisor")
    return s.with_values(op_divide(s, y.with_values(centred), threshold))


@register("divide_sum", "Divide by sum of two spectra (double divisor)", "Ratio spectra",
          [_ref("reference", "Divisor 1 (e.g. pure Y)"), _ref("reference2", "Divisor 2 (e.g. pure Z)"),
           Param("factor2", "Multiply divisor 2 by", "float", 1.0),
           Param("threshold", "Ignore divisor below", "float", 1e-4, minimum=0.0)])
def op_divide_sum(s, reference: Spectrum, reference2: Spectrum, factor2: float = 1.0,
                  threshold: float = 1e-4):
    """Double divisor ratio spectrum: mixture ÷ (Y′ + k·Z′)."""
    div = s.with_values(_on_grid(reference, s) + factor2 * _on_grid(reference2, s))
    return op_divide(s, div, threshold)


@register("multiply", "Multiply by spectrum", "Arithmetic", [_ref()])
def op_multiply(s, reference: Spectrum):
    """A × B (e.g. restore a ratio spectrum to absorbance)."""
    return s.values * _on_grid(reference, s)


@register("scale", "Multiply by constant", "Arithmetic",
          [Param("factor", "Factor", "float", 1.0)])
def op_scale(s, factor: float):
    """k·A (e.g. dilution factor, derivative scaling factor)."""
    return s.values * factor


@register("offset", "Add constant", "Arithmetic", [Param("value", "Constant", "float", 0.0)])
def op_offset(s, value: float):
    """A + c."""
    return s.values + value


@register("mean_center", "Mean centering (of this vector)", "Ratio spectra",
          [_wl("start", "From (nm)", None, optional=True),
           _wl("end", "To (nm)", None, optional=True)], changes_grid=True)
def op_mean_center(s, start: float | None = None, end: float | None = None) -> Spectrum:
    """Subtract the mean of the vector over the chosen range (MCR-spectra method)."""
    if start is not None and end is not None:
        s = op_crop(s, start, end)
    return s.with_values(s.values - s.values.mean())


@register("subtract_plateau", "Subtract plateau constant", "Ratio spectra",
          [_wl("start", "Plateau from (nm)", 300.0), _wl("end", "to (nm)", 320.0)])
def op_subtract_plateau(s, start: float, end: float):
    """Subtract the mean amplitude of a flat (plateau) region of a ratio spectrum."""
    return s.values - _region_mean(s, start, end)


# --------------------------------------------------------------------------- #
# Normalisation / conversion
# --------------------------------------------------------------------------- #
@register("normalize", "Normalize", "Normalization",
          [Param("mode", "Mode", "choice", "max",
                 choices=("max", "area", "vector", "wavelength", "range", "concentration")),
           _wl("wavelength", "At wavelength (nm)", None, optional=True)])
def op_normalize(s, mode: str, wavelength: float | None = None):
    """Scale to unit maximum, area, vector norm, value at λ, min–max range, or
    unit concentration ('concentration': ÷ the spectrum's concentration — the
    *normalized spectrum* used as divisor in AM, CV and concentration value)."""
    y = s.values
    if mode == "concentration":
        c = [v for v in s.concentrations.values() if v]
        if len(c) != 1:
            raise ValueError("unit-concentration normalisation needs a pure standard with "
                             "exactly one non-zero concentration")
        return y / c[0]
    if mode == "max":
        d = np.max(np.abs(y))
    elif mode == "area":
        d = abs(np.trapezoid(y, s.wavelengths))
    elif mode == "vector":
        d = np.linalg.norm(y)
    elif mode == "wavelength":
        if wavelength is None:
            raise ValueError("choose a wavelength")
        d = s.value_at(wavelength)
    elif mode == "range":
        if np.ptp(y) == 0:
            raise ValueError("cannot normalize a constant spectrum")
        return (y - y.min()) / np.ptp(y)
    else:
        raise ValueError(f"unknown mode {mode}")
    if d == 0:
        raise ValueError("cannot normalize by zero")
    return y / d


@register("snv", "Standard normal variate (SNV)", "Normalization")
def op_snv(s):
    """(A − mean) / SD for each spectrum."""
    sd = s.values.std(ddof=1)
    if sd == 0:
        raise ValueError("SNV: the spectrum is constant (zero standard deviation)")
    return (s.values - s.values.mean()) / sd


@register("t_to_a", "Transmittance → absorbance", "Normalization",
          [Param("percent", "Input is %T", "bool", True)])
def op_t_to_a(s, percent: bool = True):
    """A = −log10(T)."""
    t = s.values / (100.0 if percent else 1.0)
    return -np.log10(np.clip(t, 1e-12, None))


@register("a_to_t", "Absorbance → transmittance", "Normalization",
          [Param("percent", "Output %T", "bool", True)])
def op_a_to_t(s, percent: bool = True):
    """T = 10^−A."""
    return 10.0 ** (-s.values) * (100.0 if percent else 1.0)


# --------------------------------------------------------------------------- #
# Derivatives
# --------------------------------------------------------------------------- #
@register("derivative", "Derivative (D1–D4)", "Derivative",
          [Param("order", "Order", "int", 1, minimum=1, maximum=4),
           Param("method", "Method", "choice", "difference",
                 choices=("difference", "savitzky_golay")),
           Param("delta_lambda", "Δλ (nm) [difference]", "float", 4.0, minimum=0.001),
           Param("window", "Window (points) [S-G]", "int", 11, minimum=3),
           Param("polyorder", "Polynomial order [S-G]", "int", 3, minimum=1, maximum=10),
           Param("scaling", "Scaling factor", "float", 1.0,
                 help="Multiply the result (instrument 'scaling factor').")],
          changes_grid=True)
def op_derivative(s, order: int, method: str = "difference", delta_lambda: float = 4.0,
                  window: int = 11, polyorder: int = 3, scaling: float = 1.0):
    """n-th derivative dⁿA/dλⁿ by finite differences (Δλ) or Savitzky–Golay."""
    if method == "savitzky_golay":
        d = savgol_derivative(s, int(order), int(window), max(int(polyorder), int(order)))
        return s.with_values(d * scaling)
    # a Δλ difference cannot be computed within order·Δλ/2 of the range ends:
    # return only the computable part instead of fabricating edge values
    d = difference_derivative(s, int(order), delta_lambda)
    half = int(order) * float(delta_lambda) / 2
    keep = (s.wavelengths >= s.wavelengths[0] + half - 1e-9) & \
           (s.wavelengths <= s.wavelengths[-1] - half + 1e-9)
    if keep.sum() < 2:
        raise ValueError(f"the spectrum is too short for a D{order} derivative with "
                         f"Δλ = {delta_lambda:g} nm")
    return s.copy(wavelengths=s.wavelengths[keep], values=d[keep] * scaling)


# --------------------------------------------------------------------------- #
# Spectrum resolution (recover pure-component spectra from mixtures)
# --------------------------------------------------------------------------- #
@register("ratio_subtraction", "Ratio subtraction (recover X)", "Spectrum resolution",
          [_ref("divisor", "Divisor Y′ (pure extended component)"),
           _wl("start", "Plateau from (nm)", 300.0), _wl("end", "to (nm)", 320.0),
           *_DIVISOR_DERIVATIVE])
def op_ratio_subtraction(s, divisor: Spectrum, start: float, end: float,
                         divisor_derivative: int = 0, divisor_delta_lambda: float = 4.0):
    """[(X+Y)/Y′ − constant] × Y′ → spectrum of X (= spectrum subtraction of
    the constant-multiplied Y).

    The constant is the plateau of the ratio spectrum where only Y absorbs.
    With a differentiated divisor (and a differentiated sample) this is
    derivative subtraction (DS) and gives the D1…D4 spectrum of X.
    """
    d = _divisor_values(s, divisor, divisor_derivative, divisor_delta_lambda)
    ratio = op_divide(s, s.with_values(d))
    plateau = _region_mean(s.with_values(ratio), start, end)
    # [(X+Y)/Y′ − k]·Y′ written as (X+Y) − k·Y′: identical algebra, but it
    # keeps the spectrum where Y′ ≈ 0 (there the ratio is undefined)
    return s.values - plateau * d


@register("constant_multiplication", "Constant multiplication (recover Y)",
          "Spectrum resolution",
          [_ref("divisor", "Divisor Y′ (pure extended component)"),
           _wl("start", "Plateau from (nm)", 300.0), _wl("end", "to (nm)", 320.0),
           *_DIVISOR_DERIVATIVE])
def op_constant_multiplication(s, divisor: Spectrum, start: float, end: float,
                               divisor_derivative: int = 0, divisor_delta_lambda: float = 4.0):
    """constant × Y′ → spectrum of Y (constant = ratio plateau). With a
    differentiated divisor and sample this is DS-CM and gives Y's derivative."""
    d = _divisor_values(s, divisor, divisor_derivative, divisor_delta_lambda)
    ratio = op_divide(s, s.with_values(d))
    plateau = _region_mean(s.with_values(ratio), start, end)
    return plateau * d


@register("extended_ratio_subtraction", "Extended ratio subtraction (recover Y)",
          "Spectrum resolution",
          [_ref("divisor", "Divisor Y′ (pure extended component)"),
           _ref("reference", "Pure X spectrum X′"),
           _wl("start", "Plateau from (nm)", 300.0), _wl("end", "to (nm)", 320.0)])
def op_extended_ratio_subtraction(s, divisor: Spectrum, reference: Spectrum, start: float,
                                  end: float):
    """Ratio subtraction recovers X; the recovered X is matched to the pure X′
    spectrum (least-squares factor k) and Y = mixture − k·X′."""
    x_rec = op_ratio_subtraction(s, divisor, start, end)
    xp = _on_grid(reference, s)
    k = float(np.dot(x_rec, xp) / np.dot(xp, xp))
    return s.values - k * xp


def _transform(s: Spectrum, derivative_order: int, delta_lambda: float) -> Spectrum:
    if derivative_order <= 0:
        return s
    return op_derivative(s, derivative_order, "difference", delta_lambda)


_FACTOR_PARAMS = [
    _ref("reference", "Pure component spectrum"),
    _wl("wavelength", "Wavelength where only this component responds (nm)", 300.0),
    Param("derivative_order", "Signal: derivative order (0 = absorbance)", "int", 0,
          minimum=0, maximum=4),
    Param("delta_lambda", "Δλ (nm)", "float", 4.0, minimum=0.001),
    _wl("wavelength_end", "…to (nm) — plateau range (blank = single λ)", None, optional=True,
        help="Give a range to take k as the mean of signal(mixture)/signal(pure) over the "
             "plateau (constant value / derivative transformation)."),
]


def _factor(s: Spectrum, ref: Spectrum, wavelength: float, derivative_order: int,
            delta_lambda: float, wavelength_end: float | None) -> float:
    ts = _transform(s, derivative_order, delta_lambda)
    tr = _transform(ref, derivative_order, delta_lambda)
    if wavelength_end is None:
        return ts.value_at(wavelength) / tr.value_at(wavelength)
    lo, hi = sorted((wavelength, wavelength_end))
    xs, ys = ts.region(lo, hi)
    yr = np.interp(xs, tr.wavelengths, tr.values)
    if np.any(np.abs(yr) < 1e-12):
        raise ValueError("the pure-component signal is zero inside the plateau range")
    return float(np.mean(ys / yr))


@register("factorized_recovery", "Factorized spectrum / derivative transformation",
          "Spectrum resolution", _FACTOR_PARAMS)
def op_factorized_recovery(s, reference: Spectrum, wavelength: float,
                           derivative_order: int = 0, delta_lambda: float = 4.0,
                           wavelength_end: float | None = None):
    """Recover a component's zero-order spectrum: k × pure spectrum, where
    k = signal(mixture, λ) / signal(pure, λ). With derivative order 0 this is
    spectrum scaling / constant multiplication; with order ≥ 1 it is the
    factorized zero-order method or, with a plateau range, the derivative
    transformation (DT) method (k = plateau of D(mixture)/D(normalised pure))."""
    ref = s.with_values(_on_grid(reference, s))
    return _factor(s, ref, wavelength, derivative_order, delta_lambda,
                   wavelength_end) * ref.values


@register("spectrum_subtraction", "Spectrum subtraction", "Spectrum resolution",
          _FACTOR_PARAMS)
def op_spectrum_subtraction(s, reference: Spectrum, wavelength: float,
                            derivative_order: int = 0, delta_lambda: float = 4.0,
                            wavelength_end: float | None = None):
    """Mixture − k × pure spectrum (k as in factorized recovery) → the other
    component's spectrum (CM-SS, DT-SS)."""
    return s.values - op_factorized_recovery(s, reference, wavelength,
                                             derivative_order, delta_lambda, wavelength_end)


@register("constant_center", "Constant center (recover X or Y)", "Spectrum resolution",
          [_ref("divisor", "Divisor Y′"),
           _ref("reference", "Pure X spectrum (for amplitude ratio)"),
           _wl("w1", "λ1 (nm)", 230.0), _wl("w2", "λ2 (nm)", 260.0),
           Param("target", "Recover", "choice", "X", choices=("X", "Y"))])
def op_constant_center(s, divisor: Spectrum, reference: Spectrum, w1: float,
                       w2: float, target: str = "X"):
    """Mixture/Y′ = X/Y′ + k. The ratio r = P_X(λ1)/P_X(λ2) is constant for X,
    so k = P2 − (P1 − P2)/(r − 1). X = (P − k)·Y′, Y = k·Y′."""
    d = _on_grid(divisor, s)
    p = s.with_values(op_divide(s, divisor))
    px = s.with_values(op_divide(s.with_values(_on_grid(reference, s)), divisor))
    r = px.value_at(w1) / px.value_at(w2)
    if abs(r - 1) < 1e-9:
        raise ValueError("λ1 and λ2 give equal X amplitudes; choose other wavelengths")
    p1, p2 = p.value_at(w1), p.value_at(w2)
    k = p2 - (p1 - p2) / (r - 1)
    return s.values - k * d if target == "X" else k * d


# --------------------------------------------------------------------------- #
# Pipeline execution
# --------------------------------------------------------------------------- #
def apply_step(s: Spectrum, step: dict, resolve: Resolver | None = None) -> Spectrum:
    op = REGISTRY.get(step["op"])
    if op is None:
        raise KeyError(f"unknown operation {step['op']!r}")
    kwargs: dict[str, Any] = {}
    given = step.get("params", {})
    for p in op.params:
        val = given.get(p.name, p.default)
        if p.kind == "spectrum":
            if val is None:
                raise ValueError(f"{op.label}: choose a {p.label.lower()}")
            if resolve is None:
                raise ValueError("no resolver for spectrum references")
            val = val if isinstance(val, Spectrum) else resolve(val)
        elif val is None:
            if not p.optional:
                raise ValueError(f"{op.label}: '{p.label}' is required")
        elif p.kind == "int":
            val = int(val)
        elif p.kind in ("float", "wavelength"):
            val = float(val)
        elif p.kind == "bool":
            val = bool(val)
        kwargs[p.name] = val
    out = op.func(s, **kwargs)
    out = out if isinstance(out, Spectrum) else s.with_values(np.asarray(out, dtype=float))
    if not np.all(np.isfinite(out.values)):
        raise ValueError(f"{op.label} produced invalid values (NaN/∞) — check its parameters "
                         "(e.g. a wavelength too close to the range end, or division by zero)")
    return out


def apply_pipeline(s: Spectrum, steps: list[dict], resolve: Resolver | None = None) -> Spectrum:
    for step in steps:
        s = apply_step(s, step, resolve)
    return s


def describe_step(step: dict, name_of: Callable[[Any], str] | None = None) -> str:
    op = REGISTRY.get(step["op"])
    if op is None:
        return step["op"]
    parts = []
    for p in op.params:
        v = step.get("params", {}).get(p.name, p.default)
        if v is None:
            continue
        if p.kind == "spectrum":
            v = name_of(v) if name_of else v
        elif isinstance(v, float):
            v = f"{v:g}"
        parts.append(f"{p.name}={v}")
    return f"{op.label}" + (f" ({', '.join(parts)})" if parts else "")


def categories() -> dict[str, list[Operation]]:
    out: dict[str, list[Operation]] = {}
    for op in REGISTRY.values():
        out.setdefault(op.category, []).append(op)
    return out
