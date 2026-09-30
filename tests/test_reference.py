"""Every calculation compared with an independent reference implementation
(SciPy, scikit-learn, closed-form results)."""

import numpy as np
import pytest
from scipy import stats

from spectro.core import greenness as gr
from spectro.core import univariate as uv
from spectro.core import validation as val
from spectro.core.multicomponent import SignalEquations, SpectralModel, design_concentrations
from spectro.core.operations import apply_pipeline
from spectro.core.spectrum import Spectrum
from tests.conftest import GRID, mixture

RNG = np.random.default_rng(42)


# --------------------------------------------------------------------------- #
# Statistics
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("seed", range(5))
def test_regression_matches_scipy_linregress(seed):
    rng = np.random.default_rng(seed)
    x = np.sort(rng.uniform(1, 30, 8))
    y = 0.047 * x + 0.012 + rng.normal(0, 0.004, x.size)
    r = val.linear_regression(x, y)
    ref = stats.linregress(x, y)
    assert r.slope == pytest.approx(ref.slope, rel=1e-12)
    assert r.intercept == pytest.approx(ref.intercept, abs=1e-12)
    assert r.r == pytest.approx(ref.rvalue, rel=1e-12)
    assert r.sd_slope == pytest.approx(ref.stderr, rel=1e-10)
    assert r.sd_intercept == pytest.approx(ref.intercept_stderr, rel=1e-10)
    res = y - (ref.slope * x + ref.intercept)
    assert r.sy_x == pytest.approx(np.sqrt(np.sum(res ** 2) / (x.size - 2)), rel=1e-12)
    t = stats.t.ppf(0.975, x.size - 2)
    assert r.ci_slope[0] == pytest.approx(ref.slope - t * ref.stderr, rel=1e-10)
    assert r.lod == pytest.approx(3.3 * ref.intercept_stderr / ref.slope, rel=1e-10)
    assert r.loq == pytest.approx(10 * ref.intercept_stderr / ref.slope, rel=1e-10)
    r2 = val.linear_regression(x, y, lod_basis="residual")
    assert r2.lod == pytest.approx(3.3 * r.sy_x / ref.slope, rel=1e-10)
    # inverse prediction
    assert float(r.predict_x(r.predict_y(7.3))) == pytest.approx(7.3)


def test_regression_through_origin_matches_closed_form():
    x = np.array([2, 4, 6, 8, 10.0])
    y = np.array([0.098, 0.203, 0.296, 0.405, 0.498])
    r = val.linear_regression(x, y, through_origin=True)
    b = np.sum(x * y) / np.sum(x * x)
    sy = np.sqrt(np.sum((y - b * x) ** 2) / (x.size - 1))
    assert r.slope == pytest.approx(b) and r.intercept == 0
    assert r.sd_slope == pytest.approx(sy / np.sqrt(np.sum(x * x)))
    assert r.sy_x == pytest.approx(sy)


def test_regression_rejects_degenerate_input():
    with pytest.raises(ValueError):
        val.linear_regression([1, 2], [1, 2])          # too few points
    with pytest.raises(ValueError):
        val.linear_regression([5, 5, 5, 5], [1, 2, 3, 4])  # no concentration spread
    with pytest.raises(ValueError):
        val.linear_regression([1, 2, 3], [1, 2])        # length mismatch


def test_lack_of_fit_matches_manual_anova():
    x = np.repeat([2, 4, 6, 8, 10.0], 3)
    y = 0.05 * x + 0.002 * x ** 2 + RNG.normal(0, 0.002, x.size)
    lof = val.lack_of_fit(x, y)
    b, a = np.polyfit(x, y, 1)
    ss_res = np.sum((y - (b * x + a)) ** 2)
    ss_pe = sum(np.sum((y[x == v] - y[x == v].mean()) ** 2) for v in np.unique(x))
    f = ((ss_res - ss_pe) / 3) / (ss_pe / 10)
    assert lof["F"] == pytest.approx(f) and lof["p"] == pytest.approx(stats.f.sf(f, 3, 10))
    assert lof["p"] < 0.05  # the curvature is detected


def test_tests_match_scipy():
    a = RNG.normal(100, 0.8, 6)
    b = RNG.normal(100.5, 1.1, 6)
    t = val.t_test(a, b)
    ref = stats.ttest_ind(a, b)
    assert t["t"] == pytest.approx(abs(ref.statistic)) and t["p"] == pytest.approx(ref.pvalue)
    f = val.f_test(a, b)
    va, vb = a.var(ddof=1), b.var(ddof=1)
    F = max(va, vb) / min(va, vb)
    assert f["F"] == pytest.approx(F)
    assert f["p"] == pytest.approx(min(1, 2 * stats.f.sf(F, 5, 5)))
    groups = {"d1": RNG.normal(10, 0.1, 5), "d2": RNG.normal(10.05, 0.1, 5),
              "d3": RNG.normal(10.2, 0.1, 5)}
    an = val.anova_oneway(groups)
    ref = stats.f_oneway(*groups.values())
    assert an["F"] == pytest.approx(ref.statistic) and an["p"] == pytest.approx(ref.pvalue)
    assert an["ss_between"] + an["ss_within"] == pytest.approx(
        np.sum((np.concatenate(list(groups.values())) -
                np.concatenate(list(groups.values())).mean()) ** 2))


def test_interval_hypothesis_behaviour():
    a = np.array([99.8, 100.2, 100.1, 99.9, 100.0, 100.3])
    ok = val.interval_hypothesis(a, a.copy(), theta=0.02)
    assert ok["accepted"] and ok["lower"] < 1 < ok["upper"]
    biased = val.interval_hypothesis(a * 1.05, a, theta=0.02)
    assert not biased["accepted"] and biased["lower"] > 1.02


def test_descriptive_recovery_and_prediction_error():
    d = val.describe([98, 100, 102])
    assert d["mean"] == 100 and d["sd"] == pytest.approx(2) and d["rsd"] == pytest.approx(2)
    rec = val.recovery([9.9, 10.1], [10, 10])
    assert rec["recoveries"] == pytest.approx([99, 101])
    pe = val.prediction_error([1.1, 1.9, 3.0], [1, 2, 3])
    assert pe["RMSEP"] == pytest.approx(np.sqrt((0.01 + 0.01 + 0) / 3))
    assert pe["bias"] == pytest.approx(0)
    sa = val.standard_addition([0, 2, 4, 6], [0.25, 0.35, 0.45, 0.55])
    assert sa["found"] == pytest.approx(5.0)


# --------------------------------------------------------------------------- #
# Derivatives and smoothing
# --------------------------------------------------------------------------- #
def _poly(x, c):
    return sum(ci * x ** i for i, ci in enumerate(c))


@pytest.mark.parametrize("method", ["difference", "savitzky_golay"])
def test_derivatives_of_polynomials_are_exact(method):
    x = np.arange(200.0, 300.0, 0.5)
    t = (x - 250) / 50
    s = Spectrum(x, 0.3 + 0.2 * t - 0.4 * t ** 2)       # quadratic in λ
    d1 = apply_pipeline(s, [{"op": "derivative", "params": {
        "order": 1, "method": method, "delta_lambda": 2, "window": 9, "polyorder": 3}}])
    d2 = apply_pipeline(s, [{"op": "derivative", "params": {
        "order": 2, "method": method, "delta_lambda": 2, "window": 9, "polyorder": 3}}])
    t1 = (d1.wavelengths - 250) / 50
    assert np.allclose(d1.values[5:-5], (0.2 - 0.8 * t1[5:-5]) / 50, atol=1e-10)
    assert np.allclose(d2.values[5:-5], -0.8 / 2500, atol=1e-10)


def test_savitzky_golay_on_non_uniform_grid_is_correct_or_refused():
    x = np.sort(np.concatenate([np.arange(200, 250, 1.0), np.arange(250, 300, 0.5)]))
    s = Spectrum(x, 0.01 * x)
    try:
        d = apply_pipeline(s, [{"op": "derivative", "params": {
            "order": 1, "method": "savitzky_golay", "window": 7, "polyorder": 2}}])
    except ValueError:
        return  # refusing is acceptable
    assert np.allclose(d.values[10:-10], 0.01, rtol=1e-6), "wrong derivative on uneven grid"


def test_derivative_edges_are_not_silently_wrong():
    """Near the ends of the range a Δλ-difference derivative cannot be computed;
    the value there must not be a fake constant copied from the interior."""
    x = np.arange(200.0, 300.0, 1.0)
    s = Spectrum(x, (x - 200) ** 2 / 1e4)                 # derivative = 2(x−200)/1e4
    d = apply_pipeline(s, [{"op": "derivative", "params": {"order": 1, "delta_lambda": 8}}])
    true = 2 * (d.wavelengths - 200) / 1e4
    assert np.allclose(d.values, true, atol=1e-12), "edge derivative fabricated"
    assert d.wavelengths[0] >= 204 and d.wavelengths[-1] <= 295   # only computable range
    with pytest.raises(ValueError):
        d.value_at(201.0)


def test_moving_average_and_whittaker_preserve_linear_signal():
    x = np.arange(200.0, 300.0, 1.0)
    s = Spectrum(x, 0.002 * x)
    ma = apply_pipeline(s, [{"op": "smooth_ma", "params": {"window": 5}}])
    wh = apply_pipeline(s, [{"op": "smooth_whittaker", "params": {"lam": 50}}])
    assert np.allclose(ma.values[3:-3], s.values[3:-3])
    assert np.allclose(wh.values, s.values, atol=1e-8)


# --------------------------------------------------------------------------- #
# Chemometrics vs scikit-learn / linear algebra
# --------------------------------------------------------------------------- #
def _design(noise=0.0005, seed=0):
    conc = design_concentrations({"X": 10, "Y": 10, "Z": 10}, {"X": 4, "Y": 4, "Z": 4})
    return [mixture(c, f"D{i}", noise, seed=seed + i) for i, c in enumerate(conc)]


def test_pls2_matches_sklearn():
    from sklearn.cross_decomposition import PLSRegression
    cal, test = _design(), _design(seed=100)[:6]
    m = SpectralModel("PLS2", ["X", "Y", "Z"], n_components=3).fit(cal)
    X = np.vstack([np.interp(m.grid, s.x, s.y) for s in cal])
    Y = np.array([[s.concentrations[c] for c in "XYZ"] for s in cal])
    ref = PLSRegression(n_components=3, scale=False).fit(X, Y)
    Xt = np.vstack([np.interp(m.grid, s.x, s.y) for s in test])
    assert np.allclose(m.predict(test), ref.predict(Xt), atol=1e-8)


def test_pcr_and_ils_match_linear_algebra():
    cal, test = _design(), _design(seed=200)[:5]
    m = SpectralModel("PCR", ["X", "Y", "Z"], n_components=3).fit(cal)
    X = np.vstack([np.interp(m.grid, s.x, s.y) for s in cal])
    Y = np.array([[s.concentrations[c] for c in "XYZ"] for s in cal])
    xm, ym = X.mean(0), Y.mean(0)
    u, sv, vt = np.linalg.svd(X - xm, full_matrices=False)
    P = vt[:3].T
    B = P @ np.linalg.lstsq((X - xm) @ P, Y - ym, rcond=None)[0]
    Xt = np.vstack([np.interp(m.grid, s.x, s.y) for s in test])
    assert np.allclose(m.predict(test), (Xt - xm) @ B + ym, atol=1e-9)
    wl = [225.0, 245.0, 260.0, 275.0, 330.0]
    ils = SpectralModel("ILS", ["X", "Y", "Z"], wavelengths=wl).fit(cal)
    idx = [int(np.argmin(np.abs(ils.grid - w))) for w in wl]
    Xs = np.vstack([np.interp(ils.grid, s.x, s.y) for s in cal])
    coef = np.linalg.lstsq(Xs - Xs.mean(0), Y - ym, rcond=None)[0]
    Xts = np.vstack([np.interp(ils.grid, s.x, s.y) for s in test])
    assert len(idx) == ils.grid.size
    assert np.allclose(ils.predict(test), (Xts - Xs.mean(0)) @ coef + ym, atol=1e-9)


@pytest.mark.parametrize("mtype", ["CLS", "MCR-ALS"])
def test_noiseless_models_are_exact(mtype):
    cal, test = _design(0.0), _design(0.0, seed=300)[:4]
    test = [mixture({"X": 7, "Y": 13, "Z": 9}), mixture({"X": 12.5, "Y": 6, "Z": 11})]
    m = SpectralModel(mtype, ["X", "Y", "Z"]).fit(cal)
    truth = np.array([[s.concentrations[c] for c in "XYZ"] for s in test])
    tol = 1e-8 if mtype == "CLS" else 0.02
    assert np.allclose(m.predict(test), truth, rtol=tol, atol=tol)


def test_frozen_ann_and_svr_equal_sklearn():
    cal, test = _design(), _design(seed=400)[:5]
    ann = SpectralModel("ANN", ["X", "Y", "Z"], n_components=3).fit(cal)
    Xt = np.vstack([np.interp(ann.grid, s.x, s.y) for s in test])
    st = ann.state
    T = ((Xt - st["x_mean"]) / st["x_sd"]) @ st["P"] / st["ts"]
    ref = np.mean([n.predict(T).reshape(len(test), -1) for n in st["nets"]], axis=0) \
        * st["ysd"] + st["y_mean"]
    assert len(st["nets"]) == 5
    assert np.allclose(ann.predict(test), ref, atol=1e-10)
    for kernel in ("linear", "rbf", "poly"):
        svr = SpectralModel("SVR", ["X", "Y", "Z"], options={"kernel": kernel, "C": 100.0}).fit(cal)
        Xp = (np.vstack([np.interp(svr.grid, s.x, s.y) for s in test]) - svr.state["x_mean"]) \
            / svr.state["x_sd"]
        ref = np.column_stack([m.predict(Xp) for m in svr.state["models"]]) * svr.state["ysd"] \
            + svr.state["y_mean"]
        assert np.allclose(svr.predict(test), ref, atol=1e-8), kernel


def test_vip_and_hotelling_identities():
    cal = _design()
    m = SpectralModel("PLS1", ["X", "Y", "Z"], n_components=3).fit(cal)
    vip = m.vip()
    assert np.mean(vip ** 2) == pytest.approx(1.0, rel=1e-6)   # Σ VIP² = p
    d = m.diagnostics(cal)
    n, k = len(cal), m.state["pca"]["k"]
    assert np.sum(d["T2"]) == pytest.approx(k * (n - 1), rel=1e-8)


def test_cls_loo_cross_validation_matches_manual_loop():
    cal = _design()
    m = SpectralModel("CLS", ["X", "Y", "Z"])
    cv = m.cross_validate(cal, method="loo")
    pred = []
    for i in range(len(cal)):
        mm = SpectralModel("CLS", ["X", "Y", "Z"]).fit(cal[:i] + cal[i + 1:])
        pred.append(mm.predict([cal[i]])[0])
    assert np.allclose(cv["predicted"], pred, atol=1e-9)


# --------------------------------------------------------------------------- #
# Special binary methods cross-checked against each other
# --------------------------------------------------------------------------- #
def test_q_analysis_equals_simultaneous_equations(pure):
    xs = [mixture({"X": c, "Y": 0}) for c in (4, 8, 12, 16, 20)]
    ys = [mixture({"X": 0, "Y": c}) for c in (4, 8, 12, 16, 20)]
    iso = [w for w in uv.isoabsorptive_points(pure["X"], pure["Y"]) if 250 < w < 300][0]
    ax_iso, ax2 = uv.absorptivity(xs, "X", iso), uv.absorptivity(xs, "X", 245)
    ay2 = uv.absorptivity(ys, "Y", 245)
    eq = SignalEquations(["X", "Y"], [{"kind": "amplitude", "params": {"w1": iso}},
                                      {"kind": "amplitude", "params": {"w1": 245.0}}])
    eq.fit(xs + ys)
    for a, b in [(6, 10), (15, 3), (2, 18)]:
        m = mixture({"X": a, "Y": b})
        q = uv.q_analysis(m, iso, 245.0, ax_iso, ax2, ay2)
        e = eq.predict([m])[0]
        assert q["X"] == pytest.approx(e[0], rel=1e-4) and q["Y"] == pytest.approx(e[1], rel=1e-4)


def test_hpsam_returns_interferent_signal(pure):
    y = pure["Y"]
    t = y.value_at(260.0)
    w2 = [w for w in GRID if w > 290 and abs(y.value_at(w) - t) < 2e-4][0]
    mixes = [mixture({"X": 6 + a, "Y": 9}) for a in (0, 4, 8, 12)]
    r = uv.hpsam([0, 4, 8, 12], mixes, 260.0, float(w2))
    assert r["A_H"] == pytest.approx(9 * y.value_at(260.0), rel=0.03)


# --------------------------------------------------------------------------- #
# Greenness
# --------------------------------------------------------------------------- #
def test_greenness_formulas():
    s = [1, 0.5, 0.25, 1, 1, 1, 0, 1, 1, 1, 1, 1]
    w = [4, 1, 1, 1, 1, 1, 4, 1, 1, 1, 1, 1]
    assert gr.agree(s, w)["score"] == pytest.approx(round(np.dot(s, w) / sum(w), 2))
    assert gr.eco_scale({"a": 6, "b": 4})["score"] == 90
    assert gr.eco_scale({"a": 60})["rating"] == "inadequate"
    assert gr.reagent_penalty(5, 2, "danger") == 4 and gr.reagent_penalty(150, 1) == 3
    with pytest.raises(ValueError):
        gr.agree([1.2] + [1] * 11)


@pytest.mark.parametrize("call", [
    lambda: val.describe([]), lambda: val.t_test([1], [1, 2]), lambda: val.f_test([], [1, 2]),
    lambda: val.anova_oneway({"a": [1, 2]}), lambda: val.anova_oneway({"a": [1], "b": [2, 3]}),
    lambda: val.recovery([1, 2], [1]), lambda: val.recovery([1], [0]),
    lambda: val.prediction_error([], []), lambda: val.interval_hypothesis([1], [1, 2]),
    lambda: val.describe([1, float("nan")])])
def test_statistics_refuse_insufficient_or_invalid_input(call):
    with pytest.raises(ValueError):
        call()
