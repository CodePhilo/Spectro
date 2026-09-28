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
    v = np.asarray(values, float)
    mean = float(v.mean())
    sd = float(v.std(ddof=1)) if v.size > 1 else float("nan")
    return {"n": int(v.size), "mean": mean, "sd": sd,
            "rsd": 100 * sd / mean if mean else float("nan"),
            "sem": sd / np.sqrt(v.size) if v.size > 1 else float("nan"),
            "min": float(v.min()), "max": float(v.max())}


def recovery(found, taken) -> dict:
    """% recovery statistics (accuracy)."""
    found = np.asarray(found, float)
    taken = np.asarray(taken, float)
    rec = 100 * found / taken
    d = describe(rec)
    d["recoveries"] = rec.tolist()
    return d


def precision(groups: dict[str, list[float]]) -> dict:
    """%RSD per group (e.g. per concentration level / day)."""
    return {k: describe(v) for k, v in groups.items()}


def t_test(a, b, equal_var: bool = True, alpha: float = 0.05) -> dict:
    """Student's t-test comparing two methods' results (unpaired)."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    t, p = stats.ttest_ind(a, b, equal_var=equal_var)
    dof = a.size + b.size - 2
    return {"t": float(abs(t)), "p": float(p), "t_crit": float(stats.t.ppf(1 - alpha / 2, dof)),
            "dof": dof, "significant": bool(p < alpha)}


def f_test(a, b, alpha: float = 0.05) -> dict:
    """Two-tailed variance-ratio F-test (larger variance on top)."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    va, vb = a.var(ddof=1), b.var(ddof=1)
    if va >= vb:
        f, d1, d2 = va / vb, a.size - 1, b.size - 1
    else:
        f, d1, d2 = vb / va, b.size - 1, a.size - 1
    p = min(1.0, 2 * stats.f.sf(f, d1, d2))
    return {"F": float(f), "p": float(p), "F_crit": float(stats.f.ppf(1 - alpha / 2, d1, d2)),
            "dof": (d1, d2), "significant": bool(p < alpha)}


def anova_oneway(groups: dict[str, list[float]], alpha: float = 0.05) -> dict:
    data = [np.asarray(v, float) for v in groups.values()]
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
    a, b = np.asarray(test, float), np.asarray(reference, float)
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


def prediction_error(predicted, actual) -> dict:
    p, a = np.asarray(predicted, float), np.asarray(actual, float)
    err = p - a
    return {"RMSEP": float(np.sqrt(np.mean(err ** 2))),
            "bias": float(err.mean()),
            "SEP": float(err.std(ddof=1)) if err.size > 1 else float("nan"),
            "R2": float(1 - np.sum(err ** 2) / np.sum((a - a.mean()) ** 2))
            if np.ptp(a) > 0 else float("nan")}
