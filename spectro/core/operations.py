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
from typing import Any, Callable

import numpy as np
from scipy import sparse
from scipy.signal import savgol_filter
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


def savgol_derivative(s: Spectrum, order: int, window: int, polyorder: int) -> np.ndarray:
    window = int(window) | 1
    if window <= polyorder:
        raise ValueError("window must be larger than the polynomial order")
    if polyorder < order:
        raise ValueError("polynomial order must be ≥ derivative order")
    return savgol_filter(s.values, window, polyorder, deriv=order, delta=s.step)


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
    n = int(np.floor((hi - lo) / step + 1e-9)) + 1
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
    y = s.values
    n = y.size
    d = sparse.diags([1.0, -2.0, 1.0], [0, 1, 2], shape=(n - 2, n))
    dd = lam * (d.T @ d)
    w = np.ones(n)
    z = y
    for _ in range(int(iterations)):
        wm = sparse.diags(w)
        z = spsolve((wm + dd).tocsc(), w * y)
        w = p * (y > z) + (1 - p) * (y < z)
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


@register("divide", "Divide by spectrum (ratio spectrum)", "Ratio spectra",
          [_ref("reference", "Divisor spectrum"),
           Param("threshold", "Ignore divisor below", "float", 1e-4, minimum=0.0,
                 help="Points where |divisor| is below this are interpolated.")])
def op_divide(s, reference: Spectrum, threshold: float = 1e-4):
    """Ratio spectrum: mixture ÷ divisor."""
    div = _on_grid(reference, s)
    bad = np.abs(div) < threshold
    if bad.all():
        raise ValueError("divisor is zero over the whole range")
    out = np.empty_like(s.values)
    out[~bad] = s.values[~bad] / div[~bad]
    if bad.any():
        out[bad] = np.interp(s.wavelengths[bad], s.wavelengths[~bad], out[~bad])
    return out


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
                 choices=("max", "area", "vector", "wavelength", "range")),
           _wl("wavelength", "At wavelength (nm)", None, optional=True)])
def op_normalize(s, mode: str, wavelength: float | None = None):
    """Scale to unit maximum, area, vector norm, value at λ, or min–max range."""
    y = s.values
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
        return (y - y.min()) / np.ptp(y)
    else:
        raise ValueError(f"unknown mode {mode}")
    if d == 0:
        raise ValueError("cannot normalize by zero")
    return y / d


@register("snv", "Standard normal variate (SNV)", "Normalization")
def op_snv(s):
    """(A − mean) / SD for each spectrum."""
    return (s.values - s.values.mean()) / s.values.std(ddof=1)


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
           Param("polyorder", "Polynomial order [S-G]", "int", 3, minimum=1),
           Param("scaling", "Scaling factor", "float", 1.0,
                 help="Multiply the result (instrument 'scaling factor').")])
def op_derivative(s, order: int, method: str = "difference", delta_lambda: float = 4.0,
                  window: int = 11, polyorder: int = 3, scaling: float = 1.0):
    """n-th derivative dⁿA/dλⁿ by finite differences (Δλ) or Savitzky–Golay."""
    if method == "savitzky_golay":
        d = savgol_derivative(s, int(order), int(window), max(int(polyorder), int(order)))
    else:
        d = difference_derivative(s, int(order), delta_lambda)
    return d * scaling


# --------------------------------------------------------------------------- #
# Spectrum resolution (recover pure-component spectra from mixtures)
# --------------------------------------------------------------------------- #
@register("ratio_subtraction", "Ratio subtraction (recover X)", "Spectrum resolution",
          [_ref("divisor", "Divisor Y′ (pure extended component)"),
           _wl("start", "Plateau from (nm)", 300.0), _wl("end", "to (nm)", 320.0)])
def op_ratio_subtraction(s, divisor: Spectrum, start: float, end: float):
    """[(X+Y)/Y′ − constant] × Y′ → zero-order spectrum of X.

    The constant is the plateau of the ratio spectrum where only Y absorbs.
    """
    d = _on_grid(divisor, s)
    ratio = op_divide(s, divisor)
    plateau = _region_mean(s.with_values(ratio), start, end)
    return (ratio - plateau) * d


@register("constant_multiplication", "Constant multiplication (recover Y)",
          "Spectrum resolution",
          [_ref("divisor", "Divisor Y′ (pure extended component)"),
           _wl("start", "Plateau from (nm)", 300.0), _wl("end", "to (nm)", 320.0)])
def op_constant_multiplication(s, divisor: Spectrum, start: float, end: float):
    """constant × Y′ → zero-order spectrum of Y (constant = ratio plateau)."""
    ratio = op_divide(s, divisor)
    plateau = _region_mean(s.with_values(ratio), start, end)
    return plateau * _on_grid(divisor, s)


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
    return s.with_values(difference_derivative(s, derivative_order, delta_lambda))


_FACTOR_PARAMS = [
    _ref("reference", "Pure component spectrum"),
    _wl("wavelength", "Wavelength where only this component responds (nm)", 300.0),
    Param("derivative_order", "Signal: derivative order (0 = absorbance)", "int", 0,
          minimum=0, maximum=4),
    Param("delta_lambda", "Δλ (nm)", "float", 4.0, minimum=0.001),
]


@register("factorized_recovery", "Factorized spectrum / spectrum scaling",
          "Spectrum resolution", _FACTOR_PARAMS)
def op_factorized_recovery(s, reference: Spectrum, wavelength: float,
                           derivative_order: int = 0, delta_lambda: float = 4.0):
    """Recover a component's zero-order spectrum: k × pure spectrum, where
    k = signal(mixture, λ) / signal(pure, λ). With derivative order 0 this is
    spectrum scaling; with order ≥ 1 it is the factorized zero-order method
    (λ = zero-crossing of the other component)."""
    ref = s.with_values(_on_grid(reference, s))
    k = (_transform(s, derivative_order, delta_lambda).value_at(wavelength)
         / _transform(ref, derivative_order, delta_lambda).value_at(wavelength))
    return k * ref.values


@register("spectrum_subtraction", "Spectrum subtraction", "Spectrum resolution",
          _FACTOR_PARAMS)
def op_spectrum_subtraction(s, reference: Spectrum, wavelength: float,
                            derivative_order: int = 0, delta_lambda: float = 4.0):
    """Mixture − k × pure spectrum (k as in factorized recovery) → the other
    component's spectrum."""
    return s.values - op_factorized_recovery(s, reference, wavelength,
                                             derivative_order, delta_lambda)


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
    return (p.values - k) * d if target == "X" else k * d


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
    if isinstance(out, Spectrum):
        return out
    return s.with_values(np.asarray(out, dtype=float))


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
