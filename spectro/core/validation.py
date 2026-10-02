"""Regression statistics and method validation (ICH Q2(R2)) helpers."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field

import numpy as np
from scipy import stats


@dataclass
class Regression:
    slope: float
    intercept: float
    r: float
    r2: float
    n: int
    sd_slope: float
    sd_intercept: float
    sy_x: float  # residual standard deviation
    ci_slope: tuple[float, float]
    ci_intercept: tuple[float, float]
    x: list[float] = field(default_factory=list)
    y: list[float] = field(default_factory=list)
    residuals: list[float] = field(default_factory=list)
    lod: float = float("nan")
    loq: float = float("nan")
    range: tuple[float, float] = (float("nan"), float("nan"))
    through_origin: bool = False

    def predict_x(self, y: float | np.ndarray) -> float | np.ndarray:
        """Concentration from a response."""
        return (np.asarray(y) - self.intercept) / self.slope

    def predict_y(self, x: float | np.ndarray) -> float | np.ndarray:
        return self.slope * np.asarray(x) + self.intercept

    def ci_x(self, y0: float, m: int = 1, alpha: float = 0.05) -> float:
        """Half-width of the (1 − alpha) confidence interval of a concentration
        read from the line, for the mean of ``m`` responses ``y0``:
        s_x0 = (Sy/x / b)·√(1/m + 1/n + (y0 − ȳ)² / (b²·Σ(x − x̄)²)), times t(n − 2)
        (Miller & Miller; through the origin: 1/n and the means drop out, n − 1 df)."""
        x = np.asarray(self.x, dtype=float)
        if x.size < 3 or not self.slope or not np.isfinite(self.sy_x):
            return float("nan")
        b = self.slope
        if self.through_origin:
            s = (self.sy_x / abs(b)) * np.sqrt(1 / m + y0 ** 2 / (b * b * np.sum(x * x)))
            dof = x.size - 1
        else:
            ybar = float(np.mean(self.y))
            sxx = float(np.sum((x - x.mean()) ** 2))
            s = (self.sy_x / abs(b)) * np.sqrt(1 / m + 1 / x.size
                                               + (y0 - ybar) ** 2 / (b * b * sxx))
            dof = x.size - 2
        return float(stats.t.ppf(1 - alpha / 2, dof) * s)

    def ci_at_x(self, x0: float, m: int = 1, alpha: float = 0.05) -> float:
        """As :meth:`ci_x`, for the concentration ``x0`` found."""
        return self.ci_x(float(self.predict_y(x0)), m, alpha)

    def to_dict(self) -> dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, d: dict) -> "Regression":
        d = dict(d)
        for k in ("ci_slope", "ci_intercept", "range"):
            d[k] = tuple(d[k])
        return cls(**d)


def linear_regression(x, y, through_origin: bool = False, alpha: float = 0.05,
                      lod_basis: str = "intercept") -> Regression:
    """Ordinary least squares y = b·x + a with full statistics.

    ``lod_basis``: 'intercept' uses σ = SD of the intercept, 'residual' uses
    σ = Sy/x (both accepted by ICH Q2); LOD = 3.3σ/S, LOQ = 10σ/S.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if x.size != y.size:
        raise ValueError("x and y differ in length")
    n = x.size
    p = 1 if through_origin else 2
    if n < p + 1:
        raise ValueError(f"at least {p + 1} calibration points are needed")
    if not (np.all(np.isfinite(x)) and np.all(np.isfinite(y))):
        raise ValueError("calibration data contain missing or non-numeric values")
    if np.ptp(x) == 0:
        raise ValueError("all calibration concentrations are equal — a calibration line "
                         "needs at least two different concentrations")
    if through_origin:
        b = float(np.sum(x * y) / np.sum(x * x))
        a = 0.0
    else:
        b, a = (float(v) for v in np.polyfit(x, y, 1))
    fit = b * x + a
    res = y - fit
    dof = n - p
    sy_x = float(np.sqrt(np.sum(res ** 2) / dof))
    sxx = float(np.sum((x - x.mean()) ** 2))
    if through_origin:
        sd_b = sy_x / np.sqrt(np.sum(x * x))
        sd_a = 0.0
    else:
        sd_b = sy_x / np.sqrt(sxx)
        sd_a = sy_x * np.sqrt(np.sum(x * x) / (n * sxx))
    t = stats.t.ppf(1 - alpha / 2, dof)
    r = float(np.corrcoef(x, y)[0, 1]) if np.ptp(x) > 0 and np.ptp(y) > 0 else float("nan")
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    r2 = 1 - float(np.sum(res ** 2)) / ss_tot if ss_tot > 0 else float("nan")
    sigma = sd_a if (lod_basis == "intercept" and not through_origin) else sy_x
    return Regression(
        slope=b, intercept=a, r=r, r2=r2, n=n, sd_slope=float(sd_b),
        sd_intercept=float(sd_a), sy_x=sy_x,
        ci_slope=(b - t * sd_b, b + t * sd_b),
        ci_intercept=(a - t * sd_a, a + t * sd_a),
        x=x.tolist(), y=y.tolist(), residuals=res.tolist(),
        lod=abs(3.3 * sigma / b) if b else float("nan"),
        loq=abs(10 * sigma / b) if b else float("nan"),
        range=(float(x.min()), float(x.max())), through_origin=through_origin)


def _need(values, n: int, what: str) -> np.ndarray:
    v = np.asarray(values, float).ravel()
    if v.size < n:
        raise ValueError(f"{what}: at least {n} value{'s' if n > 1 else ''} needed "
                         f"(got {v.size})")
    if not np.all(np.isfinite(v)):
        raise ValueError(f"{what}: contains missing or non-numeric values")
    return v


def lack_of_fit(x, y) -> dict:
    """Lack-of-fit F-test (requires replicate x levels)."""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    levels = np.unique(x)
    n, k = x.size, levels.size
    if k == n or k < 3:
        return {"available": False, "reason": "needs replicates at ≥3 levels"}
    b, a = np.polyfit(x, y, 1)
    ss_pe = sum(np.sum((y[x == lv] - y[x == lv].mean()) ** 2) for lv in levels)
    ss_res = np.sum((y - (b * x + a)) ** 2)
    ss_lof = ss_res - ss_pe
    df_lof, df_pe = k - 2, n - k
    f = (ss_lof / df_lof) / (ss_pe / df_pe) if ss_pe > 0 else float("inf")
    return {"available": True, "F": float(f), "df": (df_lof, df_pe),
            "p": float(stats.f.sf(f, df_lof, df_pe))}


def describe(values) -> dict:
    v = _need(values, 1, "statistics")
    mean = float(v.mean())
    sd = float(v.std(ddof=1)) if v.size > 1 else float("nan")
    return {"n": int(v.size), "mean": mean, "sd": sd,
            "rsd": 100 * sd / mean if mean else float("nan"),
            "sem": sd / np.sqrt(v.size) if v.size > 1 else float("nan"),
            "min": float(v.min()), "max": float(v.max())}


def recovery(found, taken) -> dict:
    """% recovery statistics (accuracy)."""
    found = _need(found, 1, "found values")
    taken = _need(taken, 1, "taken values")
    if found.size != taken.size:
        raise ValueError("found and taken values differ in number")
    if np.any(taken == 0):
        raise ValueError("a taken (added) amount is zero — recovery is undefined")
    rec = 100 * found / taken
    d = describe(rec)
    d["recoveries"] = rec.tolist()
    return d


def precision(groups: dict[str, list[float]]) -> dict:
    """%RSD per group (e.g. per concentration level / day)."""
    return {k: describe(v) for k, v in groups.items()}


def t_test(a, b, equal_var: bool = True, alpha: float = 0.05) -> dict:
    """Student's t-test comparing two methods' results (unpaired)."""
    a, b = _need(a, 2, "first method"), _need(b, 2, "second method")
    t, p = stats.ttest_ind(a, b, equal_var=equal_var)
    dof = a.size + b.size - 2
    return {"t": float(abs(t)), "p": float(p), "t_crit": float(stats.t.ppf(1 - alpha / 2, dof)),
            "dof": dof, "significant": bool(p < alpha)}


def f_test(a, b, alpha: float = 0.05, tails: int = 1) -> dict:
    """Variance-ratio F-test, larger variance on top.

    ``tails=1`` (default) compares with F(1−α; ν1, ν2) — the tabulated value
    used in pharmaceutical method-comparison papers (e.g. 5.05 for n = 6 and 6).
    ``tails=2`` uses F(1−α/2) (Miller & Miller two-sided test, 7.15 for 6/6).
    Both critical values are returned."""
    a, b = _need(a, 2, "first method"), _need(b, 2, "second method")
    return f_test_summary(float(a.std(ddof=1)), a.size, float(b.std(ddof=1)), b.size,
                          alpha, tails)


def _check_summary(sd: float, n: int, what: str) -> None:
    if not (np.isfinite(sd) and sd >= 0):
        raise ValueError(f"{what}: the standard deviation must be a non-negative number")
    if int(n) != n or n < 2:
        raise ValueError(f"{what}: n must be a whole number ≥ 2")


def t_test_summary(mean1: float, sd1: float, n1: int, mean2: float, sd2: float, n2: int,
                   alpha: float = 0.05, equal_var: bool = True) -> dict:
    """Student's t-test from summary statistics (mean, SD, n of each method) —
    e.g. to check a published comparison with a reported/official method.
    Pooled variance (equal_var) or Welch's test."""
    _check_summary(sd1, n1, "first method")
    _check_summary(sd2, n2, "second method")
    v1, v2 = sd1 ** 2, sd2 ** 2
    if equal_var:
        dof = n1 + n2 - 2
        sp2 = ((n1 - 1) * v1 + (n2 - 1) * v2) / dof
        se = np.sqrt(sp2 * (1 / n1 + 1 / n2))
    else:
        se = np.sqrt(v1 / n1 + v2 / n2)
        dof = (v1 / n1 + v2 / n2) ** 2 / ((v1 / n1) ** 2 / (n1 - 1) + (v2 / n2) ** 2 / (n2 - 1))
    if se == 0:
        raise ValueError("both standard deviations are zero — the t-test is undefined")
    t = abs(mean1 - mean2) / se
    p = 2 * stats.t.sf(t, dof)
    return {"t": float(t), "p": float(p), "t_crit": float(stats.t.ppf(1 - alpha / 2, dof)),
            "dof": float(dof) if not equal_var else int(dof), "significant": bool(p < alpha)}


def f_test_summary(sd1: float, n1: int, sd2: float, n2: int, alpha: float = 0.05,
                   tails: int = 1) -> dict:
    """F-test from two standard deviations (larger variance on top)."""
    _check_summary(sd1, n1, "first method")
    _check_summary(sd2, n2, "second method")
    if tails not in (1, 2):
        raise ValueError("tails must be 1 or 2")
    v1, v2 = sd1 ** 2, sd2 ** 2
    if v1 == 0 and v2 == 0:
        raise ValueError("both standard deviations are zero — the F-test is undefined")
    if v1 >= v2:
        f, d1, d2 = (v1 / v2 if v2 else float("inf")), int(n1) - 1, int(n2) - 1
    else:
        f, d1, d2 = v2 / v1, int(n2) - 1, int(n1) - 1
    one = float(stats.f.ppf(1 - alpha, d1, d2))
    two = float(stats.f.ppf(1 - alpha / 2, d1, d2))
    p = float(stats.f.sf(f, d1, d2)) * tails
    return {"F": float(f), "p": min(1.0, p), "F_crit": one if tails == 1 else two,
            "F_crit_one_tailed": one, "F_crit_two_tailed": two, "tails": tails,
            "dof": (d1, d2), "significant": bool(f > (one if tails == 1 else two))}


def anova_summary(groups: dict[str, tuple[float, float, int]], alpha: float = 0.05) -> dict:
    """One-way ANOVA from (mean, SD, n) of each group."""
    if len(groups) < 2:
        raise ValueError("ANOVA needs at least two groups")
    for k, (_, sd, n) in groups.items():
        _check_summary(sd, n, f"group '{k}'")
    m = np.array([g[0] for g in groups.values()], float)
    sd = np.array([g[1] for g in groups.values()], float)
    n = np.array([g[2] for g in groups.values()], float)
    k, total = m.size, n.sum()
    grand = float(np.sum(n * m) / total)
    ss_b = float(np.sum(n * (m - grand) ** 2))
    ss_w = float(np.sum((n - 1) * sd ** 2))
    df_b, df_w = k - 1, int(total - k)
    ms_b, ms_w = ss_b / df_b, ss_w / df_w
    if ms_w == 0:
        raise ValueError("all within-group standard deviations are zero")
    f = ms_b / ms_w
    p = float(stats.f.sf(f, df_b, df_w))
    return {"F": float(f), "p": p, "F_crit": float(stats.f.ppf(1 - alpha, df_b, df_w)),
            "ss_between": ss_b, "ss_within": ss_w, "df_between": df_b, "df_within": df_w,
            "ms_between": float(ms_b), "ms_within": float(ms_w), "grand_mean": grand,
            "significant": bool(p < alpha)}


def anova_oneway(groups: dict[str, list[float]], alpha: float = 0.05) -> dict:
    data = [_need(v, 2, f"group '{k}'") for k, v in groups.items()]
    if len(data) < 2:
        raise ValueError("ANOVA needs at least two groups")
    f, p = stats.f_oneway(*data)
    k = len(data)
    n = sum(d.size for d in data)
    grand = np.concatenate(data).mean()
    ss_b = sum(d.size * (d.mean() - grand) ** 2 for d in data)
    ss_w = sum(np.sum((d - d.mean()) ** 2) for d in data)
    return {"F": float(f), "p": float(p), "F_crit": float(stats.f.ppf(1 - alpha, k - 1, n - k)),
            "ss_between": float(ss_b), "ss_within": float(ss_w),
            "df_between": k - 1, "df_within": n - k,
            "ms_between": float(ss_b / (k - 1)), "ms_within": float(ss_w / (n - k)),
            "significant": bool(p < alpha)}


def interval_hypothesis(test, reference, theta: float = 0.02, alpha: float = 0.05) -> dict:
    """Interval hypothesis test (Hartmann et al.): accept if the 1−2α CI of
    the ratio of means lies within [1−θ, 1+θ]. Returns the lower/upper limits
    of the ratio (as fractions)."""
    a, b = _need(test, 2, "test method"), _need(reference, 2, "reference method")
    na, nb = a.size, b.size
    ma, mb = a.mean(), b.mean()
    sp2 = ((na - 1) * a.var(ddof=1) + (nb - 1) * b.var(ddof=1)) / (na + nb - 2)
    t = stats.t.ppf(1 - alpha, na + nb - 2)
    # Fieller-type solution of the quadratic for the ratio θ = μa/μb
    qa = mb ** 2 - t ** 2 * sp2 / nb
    qb = -2 * ma * mb
    qc = ma ** 2 - t ** 2 * sp2 / na
    disc = qb ** 2 - 4 * qa * qc
    if qa <= 0 or disc < 0:
        return {"lower": float("nan"), "upper": float("nan"), "accepted": False}
    lo = (-qb - np.sqrt(disc)) / (2 * qa)
    hi = (-qb + np.sqrt(disc)) / (2 * qa)
    return {"lower": float(lo), "upper": float(hi),
            "accepted": bool(lo >= 1 - theta and hi <= 1 + theta)}


def standard_addition(added, response) -> dict:
    """Standard addition: concentration in the sample = intercept/slope."""
    reg = linear_regression(added, response)
    return {"regression": reg.to_dict(), "found": reg.intercept / reg.slope}


# Critical values of Dixon's Q (two-sided, Rorabacher 1991), n = 3…10
_DIXON = {0.10: [0.941, 0.765, 0.642, 0.560, 0.507, 0.468, 0.437, 0.412],
          0.05: [0.970, 0.829, 0.710, 0.625, 0.568, 0.526, 0.493, 0.466],
          0.01: [0.994, 0.926, 0.821, 0.740, 0.680, 0.634, 0.598, 0.568]}


def grubbs(values, alpha: float = 0.05) -> dict:
    """Two-sided Grubbs test for one outlier (ISO 5725-2).

    G = max|x − x̄| / s, compared with
    G_crit = (n − 1)/√n · √(t² / (n − 2 + t²)), t = t(alpha / (2n), n − 2)."""
    v = _need(values, 3, "Grubbs test")
    n = v.size
    sd = float(np.std(v, ddof=1))
    if sd == 0:
        return {"test": "Grubbs", "n": n, "G": 0.0, "G_crit": float("nan"), "index": None,
                "value": None, "outlier": False}
    i = int(np.argmax(np.abs(v - v.mean())))
    g = float(abs(v[i] - v.mean()) / sd)
    t = stats.t.ppf(1 - alpha / (2 * n), n - 2)
    gc = float((n - 1) / np.sqrt(n) * np.sqrt(t * t / (n - 2 + t * t)))
    return {"test": "Grubbs", "n": n, "G": g, "G_crit": gc, "index": i, "value": float(v[i]),
            "outlier": g > gc, "alpha": alpha}


def dixon(values, alpha: float = 0.05) -> dict:
    """Dixon's Q test (r10 = gap / range) for 3 ≤ n ≤ 10."""
    v = _need(values, 3, "Dixon test")
    n = v.size
    if n > 10:
        raise ValueError("Dixon's Q test is tabulated for 3–10 values; use Grubbs")
    if alpha not in _DIXON:
        raise ValueError("alpha must be 0.10, 0.05 or 0.01")
    order = np.argsort(v)
    s = v[order]
    rng = float(s[-1] - s[0])
    qc = _DIXON[alpha][n - 3]
    if rng == 0:
        return {"test": "Dixon", "n": n, "Q": 0.0, "Q_crit": qc, "index": None, "value": None,
                "outlier": False, "alpha": alpha}
    q_low, q_high = (s[1] - s[0]) / rng, (s[-1] - s[-2]) / rng
    low = q_low >= q_high
    q = float(q_low if low else q_high)
    i = int(order[0] if low else order[-1])
    return {"test": "Dixon", "n": n, "Q": q, "Q_crit": qc, "index": i, "value": float(v[i]),
            "outlier": q > qc, "alpha": alpha}


def prediction_error(predicted, actual) -> dict:
    p, a = _need(predicted, 1, "predicted values"), _need(actual, 1, "actual values")
    if p.size != a.size:
        raise ValueError("predicted and actual values differ in number")
    err = p - a
    return {"RMSEP": float(np.sqrt(np.mean(err ** 2))),
            "bias": float(err.mean()),
            "SEP": float(err.std(ddof=1)) if err.size > 1 else float("nan"),
            "R2": float(1 - np.sum(err ** 2) / np.sum((a - a.mean()) ** 2))
            if np.ptp(a) > 0 else float("nan")}


# --------------------------------------------------------------------------- #
# Validation plan (ICH Q2(R2)) for a chosen method
# --------------------------------------------------------------------------- #
def _sig(v: float, digits: int = 2) -> float:
    """Round to ``digits`` significant figures (pipette-friendly levels)."""
    if v == 0 or not np.isfinite(v):
        return float(v)
    return float(f"{v:.{digits}g}")


def _levels(low: float, high: float, n: int) -> list[float]:
    """About ``n`` levels from low to high: the ends plus round multiples of
    a 1-2-2.5-5 step in between (easy to prepare from a stock solution)."""
    raw = (high - low) / (n - 1)
    mag = 10 ** np.floor(np.log10(raw))
    step = min(k * mag for k in (1, 2, 2.5, 5, 10) if k * mag >= raw * 0.999)
    inner = [k * step for k in range(int(np.ceil(low / step)), int(high / step) + 1)]
    inner = [v for v in inner if v - low >= 0.5 * step and high - v >= 0.5 * step]
    return [_sig(v, 3) for v in [low] + inner + [high]]


def validation_plan(compound: str, low: float, high: float, unit: str = "µg/mL",
                    label: str = "", steps: list[dict] | None = None,
                    measurement: dict | None = None, wl_uncertainty: float = 1.0,
                    expected_error: float = float("nan"), lod: float = float("nan"),
                    multivariate: bool = False) -> list[list[str]]:
    """A ready-to-run validation protocol: [parameter, design, acceptance].

    Levels are spread over the calibration range of the standards; the
    robustness factors are the method's own settings (wavelengths, Δλ,
    smoothing window, plateau…)."""
    if not (np.isfinite(low) and np.isfinite(high)) or high <= low:
        raise ValueError(f"no concentration range for {compound}: give standards with at "
                         "least two different concentrations")
    lin = _levels(low, high, 6 if high / max(low, 1e-12) > 4 else 5)
    target = min(lin, key=lambda v: abs(v - (low + high) / 2))
    acc = [_sig(target * f) for f in (0.8, 1.0, 1.2)]
    if acc[0] < low or acc[-1] > high:
        acc = [_sig(low + (high - low) * f) for f in (0.2, 0.5, 0.8)]
    f = ", ".join

    def lv(vals):
        return f(f"{v:g}" for v in vals) + f" {unit}"
    rows = [
        ["Method", label or compound, ""],
        ["Specificity",
         "laboratory-prepared mixtures covering the expected ratios of the compounds (at "
         "least 5), the placebo / excipients and the solvent blank; compare the processed "
         "spectra of mixture and pure standard",
         "recovery of each mixture 98–102 %; blank and placebo give no signal at the "
         "measured wavelength(s)"],
        ["Linearity", f"{len(lin)} levels: {lv(lin)}, each in triplicate (independent "
                      "dilutions)" + (" — calibration set of mixtures (e.g. Brereton 5-level "
                                      "design, Tools → Calibration design)"
                                      if multivariate else ""),
         "r ≥ 0.999; residuals random (no trend); 95 % CI of the intercept includes 0 or "
         "|intercept| ≤ 2 % of the response at the target"
         if not multivariate else "RMSECV / RMSEP ≤ 2 % of the mean; no trend in residuals"],
        ["Range", f"{lin[0]:g}–{lin[-1]:g} {unit} (target {target:g} {unit})",
         "linearity, accuracy and precision acceptable over the whole range"],
        ["Accuracy", f"3 levels × 3 replicates: {lv(acc)} (laboratory mixtures, or standard "
                     "addition to the dosage form)",
         "mean recovery 98–102 %; RSD ≤ 2 %"],
        ["Repeatability", f"6 determinations at {target:g} {unit} (or the 3 × 3 accuracy "
                          "set), same day, analyst, instrument",
         "RSD ≤ 2 %"],
        ["Intermediate precision", f"the 3 accuracy levels on 3 different days "
                                   f"({lv(acc)})",
         "RSD ≤ 2 %; one-way ANOVA between days not significant (p > 0.05)"],
        ["LOD / LOQ", "from the calibration: 3.3σ/S and 10σ/S (σ = SD of the intercept or "
                      "Sy/x); confirm the LOQ with 6 replicates"
         + (f" (expected LOD ≈ {lod:.3g} {unit})" if np.isfinite(lod) else ""),
         f"LOQ ≤ {lin[0]:g} {unit} (the lowest level); RSD at the LOQ ≤ 10 %"],
    ]
    rob = []
    p = (measurement or {}).get("params", {})
    for key in ("w1", "w2", "w3"):
        if p.get(key) is not None:
            rob.append(f"wavelength {p[key]:g} ± {wl_uncertainty:g} nm")
    for key in ("p1", "p2"):
        if p.get(key) is not None:
            rob.append(f"plateau limit {p[key]:g} ± {wl_uncertainty:g} nm")
    for st in steps or []:
        sp = st.get("params", {})
        if st.get("op") == "derivative" and sp.get("method", "difference") == "difference":
            rob.append(f"Δλ {sp.get('delta_lambda', 4):g} ± 1 nm")
        if st.get("op") == "derivative" and sp.get("method") == "savgol":
            rob.append(f"S-G window {sp.get('window')} ± 2 points")
        if st.get("op") in ("smooth_sg", "smooth_ma"):
            rob.append(f"smoothing window {sp.get('window')} ± 2 points")
    rob.append("a second batch of solvent; ±10 % in the dilution volume of the stock")
    rows.append(["Robustness", "vary one factor at a time, 3 determinations each at the "
                               f"target: {'; '.join(rob)} (Methods → Univariate → "
                               "Robustness study does this automatically)",
                 "recovery of each variation within ±2 % of the nominal; RSD ≤ 2 %"])
    rows.append(["Comparison", "the dosage form by this method and by the official or a "
                               "reported method, 5–6 determinations each",
                 "Student's t and F below the tabulated values (Tools → Validation "
                 "statistics → Compare methods / Compare (mean, SD, n))"])
    if np.isfinite(expected_error):
        rows.append(["Expected", f"the optimizer predicts an error of about "
                                 f"{expected_error:.2f} % for this method",
                     "recoveries much worse than this point to a preparation or "
                     "instrument problem rather than the method"])
    return rows
