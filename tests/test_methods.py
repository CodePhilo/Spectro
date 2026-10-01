"""Every method must recover known concentrations of synthetic mixtures."""

import numpy as np
import pytest

from spectro.core import univariate as uv
from spectro.core.multicomponent import (SignalEquations, SpectralModel, brereton_design,
                                         design_concentrations, ga_select, ipls,
                                         kaiser_selection)
from spectro.core.operations import REGISTRY, apply_pipeline, describe_step
from spectro.core.spectrum import Spectrum
from spectro.core.validation import linear_regression
from tests.conftest import GRID, mixture

TOL = 0.02  # 2 % relative error


def close(found, expected, tol=TOL):
    return abs(found - expected) <= tol * max(abs(expected), 1.0)


def resolver(pure):
    lib = {f"pure:{k}": v for k, v in pure.items()}
    return lambda ref: lib[ref]


# --------------------------------------------------------------------------- #
# univariate pipelines
# --------------------------------------------------------------------------- #
def test_direct_measurement_at_plateau_region(binary_set, pure):
    xs, ys, mixes = binary_set
    # Y absorbs alone at 350 nm
    m = uv.UnivariateMethod("Y direct", "Y", [], {"kind": "amplitude", "params": {"w1": 350.0}})
    reg = m.calibrate(ys)
    assert reg.r2 > 0.9999
    for mix in mixes:
        assert close(m.predict(mix), mix.concentrations["Y"])


def test_zero_crossing_derivative(binary_set, pure):
    xs, ys, mixes = binary_set
    d1 = apply_pipeline(pure["Y"], [{"op": "derivative", "params": {"order": 1,
                                                                    "delta_lambda": 2}}])
    zc = uv.zero_crossings(d1, 240, 300)
    assert zc, "Y has a D1 zero crossing near its maximum"
    steps = [{"op": "derivative", "params": {"order": 1, "delta_lambda": 2}}]
    m = uv.UnivariateMethod("X D1", "X", steps, {"kind": "amplitude", "params": {"w1": zc[0]}})
    m.calibrate(xs)
    for mix in mixes:
        assert close(m.predict(mix), mix.concentrations["X"])


def test_ratio_difference_and_derivative_ratio(binary_set, pure):
    xs, ys, mixes = binary_set
    res = resolver(pure)
    div = [{"op": "divide", "params": {"reference": "pure:Y"}}]
    rd = uv.UnivariateMethod("RD", "X", div, {"kind": "difference",
                                              "params": {"w1": 240.0, "w2": 260.0}})
    rd.calibrate(xs, res)
    dd = uv.UnivariateMethod("DD1", "X", div + [{"op": "derivative",
                                                  "params": {"order": 1, "delta_lambda": 2}}],
                             {"kind": "amplitude", "params": {"w1": 235.0}})
    dd.calibrate(xs, res)
    mc = uv.UnivariateMethod("MCR", "X", div + [{"op": "mean_center",
                                                  "params": {"start": 220.0, "end": 300.0}}],
                             {"kind": "amplitude", "params": {"w1": 245.0}})
    mc.calibrate(xs, res)
    for mix in mixes:
        for m in (rd, dd, mc):
            assert close(m.predict(mix, res), mix.concentrations["X"]), m.name


def test_dual_and_induced_dual_wavelength(binary_set, pure):
    xs, ys, mixes = binary_set
    res = resolver(pure)
    # Induced dual wavelength: cancel Y using its own ratio F
    idw = uv.UnivariateMethod("IDW", "X", [], {"kind": "weighted_difference",
                                               "params": {"w1": 245.0, "w2": 290.0,
                                                          "reference": "pure:Y"}})
    idw.calibrate(xs, res)
    for mix in mixes:
        assert close(idw.predict(mix, res), mix.concentrations["X"])


def test_ratio_subtraction_and_constant_multiplication(binary_set, pure):
    xs, ys, mixes = binary_set
    res = resolver(pure)
    rs = uv.UnivariateMethod("RS", "X", [{"op": "ratio_subtraction", "params": {
        "divisor": "pure:Y", "start": 340.0, "end": 360.0}}],
        {"kind": "amplitude", "params": {"w1": 245.0}})
    rs.calibrate(xs, res)
    cm = uv.UnivariateMethod("CM", "Y", [{"op": "constant_multiplication", "params": {
        "divisor": "pure:Y", "start": 340.0, "end": 360.0}}],
        {"kind": "amplitude", "params": {"w1": 275.0}})
    cm.calibrate(ys, res)
    for mix in mixes:
        assert close(rs.predict(mix, res), mix.concentrations["X"])
        assert close(cm.predict(mix, res), mix.concentrations["Y"])


def test_spectrum_subtraction_factorized_constant_center(binary_set, pure):
    xs, ys, mixes = binary_set
    res = resolver(pure)
    ss = uv.UnivariateMethod("SS", "X", [{"op": "spectrum_subtraction", "params": {
        "reference": "pure:Y", "wavelength": 350.0}}],
        {"kind": "amplitude", "params": {"w1": 245.0}})
    ss.calibrate(xs, res)
    cc = uv.UnivariateMethod("CC", "X", [{"op": "constant_center", "params": {
        "divisor": "pure:Y", "reference": "pure:X", "w1": 240.0, "w2": 260.0,
        "target": "X"}}], {"kind": "amplitude", "params": {"w1": 245.0}})
    cc.calibrate(xs, res)
    for mix in mixes:
        assert close(ss.predict(mix, res), mix.concentrations["X"])
        assert close(cc.predict(mix, res), mix.concentrations["X"])


def test_area_under_curve_equations(binary_set):
    xs, ys, mixes = binary_set
    sig = [{"kind": "area", "params": {"w1": 230.0, "w2": 255.0}},
           {"kind": "area", "params": {"w1": 320.0, "w2": 345.0}}]
    eq = SignalEquations(["X", "Y"], sig)
    eq.fit(xs + ys)
    pred = eq.predict(mixes)
    for p, mix in zip(pred, mixes):
        assert close(p[0], mix.concentrations["X"]) and close(p[1], mix.concentrations["Y"])


def test_vierordt_and_bivariate(binary_set, pure):
    xs, ys, mixes = binary_set
    sig = [{"kind": "amplitude", "params": {"w1": 245.0}},
           {"kind": "amplitude", "params": {"w1": 330.0}}]
    for intercept in (False, True):
        eq = SignalEquations(["X", "Y"], sig, intercept=intercept)
        eq.fit(xs + ys)
        for p, mix in zip(eq.predict(mixes), mixes):
            assert close(p[0], mix.concentrations["X"])
            assert close(p[1], mix.concentrations["Y"])
    best = kaiser_selection({"X": pure["X"], "Y": pure["Y"]}, np.arange(220, 360, 5.0))
    assert best[0]["det"] > 0


def test_q_analysis_absorbance_subtraction_amplitude_modulation(binary_set, pure):
    xs, ys, mixes = binary_set
    iso = [w for w in uv.isoabsorptive_points(pure["X"], pure["Y"]) if 250 < w < 300][0]
    ax_iso = uv.absorptivity(xs, "X", iso)
    ax2 = uv.absorptivity(xs, "X", 245.0)
    ay2 = uv.absorptivity(ys, "Y", 245.0)
    for mix in mixes:
        q = uv.q_analysis(mix, iso, 245.0, ax_iso, ax2, ay2)
        assert close(q["X"], mix.concentrations["X"]) and close(q["Y"], mix.concentrations["Y"])

    # Y extends alone at 350 nm → Y plays the role of "X" in AS
    ab = uv.AbsorbanceSubtraction(iso, 350.0)
    ab.fit(ys, xs + ys)
    for mix in mixes:
        r = ab.predict(mix)
        assert close(r["X"], mix.concentrations["Y"]) and close(r["Y"], mix.concentrations["X"])

    divisor = pure["Y"].copy(values=pure["Y"].values * 10)
    am = uv.AmplitudeModulation(iso, (340.0, 360.0), divisor_conc=10.0)
    am.fit(xs + ys, divisor)
    for mix in mixes:
        r = am.predict(mix, divisor)
        assert close(r["Y"], mix.concentrations["Y"]) and close(r["X"], mix.concentrations["X"])


def test_hpsam(pure):
    # Sample: X=6, Y=9. Additions of X. Choose w1/w2 where Y has equal absorbance.
    y = pure["Y"]
    target = y.value_at(260.0)
    cands = [w for w in GRID if w > 290 and abs(y.value_at(w) - target) < 2e-4]
    w2 = cands[0]
    added = [0, 4, 8, 12]
    mixes = [mixture({"X": 6 + a, "Y": 9}) for a in added]
    r = uv.hpsam(added, mixes, 260.0, float(w2))
    assert abs(r["X"] - 6) < 0.3


def test_ternary_successive_and_double_divisor(pure):
    cal = [mixture({"X": c, "Y": 0, "Z": 0}) for c in (4, 8, 12, 16, 20)]
    mixes = [mixture({"X": a, "Y": b, "Z": c}) for a, b, c in [(6, 10, 8), (12, 6, 14), (16, 12, 5)]]
    lib = dict(pure)
    lib["sum"] = pure["Y"].copy(values=pure["Y"].values + pure["Z"].values)
    res = lambda k: lib[k]  # noqa: E731
    # Double divisor: divide by (Y + Z) sum spectrum is not exact for unequal
    # ratios; use successive ratio instead: /Y → D1 → /d(Z/Y) → D1
    zy = apply_pipeline(pure["Z"], [{"op": "divide", "params": {"reference": "Y"}},
                                    {"op": "derivative", "params": {"order": 1, "delta_lambda": 1}}],
                        res)
    lib["dZY"] = zy
    steps = [{"op": "crop", "params": {"start": 215.0, "end": 300.0}},
             {"op": "divide", "params": {"reference": "Y"}},
             {"op": "derivative", "params": {"order": 1, "delta_lambda": 1}},
             {"op": "divide", "params": {"reference": "dZY", "threshold": 1e-3}},
             {"op": "derivative", "params": {"order": 1, "delta_lambda": 1}}]
    m = uv.UnivariateMethod("SDR", "X", steps, {"kind": "amplitude", "params": {"w1": 240.0}})
    reg = m.calibrate(cal, res)
    assert reg.r2 > 0.999
    for mix in mixes:
        assert close(m.predict(mix, res), mix.concentrations["X"], 0.03)


# --------------------------------------------------------------------------- #
# chemometrics
# --------------------------------------------------------------------------- #
def _design_set(noise=0.0005):
    conc = design_concentrations({"X": 10, "Y": 10, "Z": 10}, {"X": 4, "Y": 4, "Z": 4})
    return [mixture(c, f"D{i}", noise, seed=i) for i, c in enumerate(conc)]


@pytest.mark.parametrize("mtype,k", [("CLS", 3), ("PCR", 3), ("PLS1", 3), ("PLS2", 3),
                                     ("MCR-ALS", 3), ("SVR", 3), ("ANN", 3)])
def test_multivariate_models(mtype, k):
    cal = _design_set()
    test = [mixture({"X": 8, "Y": 12, "Z": 9}), mixture({"X": 13, "Y": 7, "Z": 11})]
    m = SpectralModel(mtype, ["X", "Y", "Z"], ranges=[(210, 370)], n_components=k)
    m.fit(cal)
    pred = m.predict(test)
    tol = 0.06 if mtype == "ANN" else 0.03
    for p, s in zip(pred, test):
        for j, c in enumerate(["X", "Y", "Z"]):
            assert close(p[j], s.concentrations[c], tol), (mtype, c, p[j])


def test_ils_and_cross_validation():
    cal = _design_set()
    m = SpectralModel("ILS", ["X", "Y", "Z"], wavelengths=[225, 245, 260, 275, 330])
    m.fit(cal)
    p = m.predict([mixture({"X": 8, "Y": 12, "Z": 9})])[0]
    assert close(p[0], 8) and close(p[1], 12) and close(p[2], 9)
    pls = SpectralModel("PLS2", ["X", "Y", "Z"], n_components=3)
    cv = pls.cross_validate(cal, method="venetian", folds=5, max_components=6)
    assert cv["suggested_components"] == 3
    assert pls.vip() is not None
    d = pls.diagnostics(cal)
    assert len(d["T2"]) == len(cal)


def test_variable_selection_runs():
    cal = _design_set()
    m = SpectralModel("PLS2", ["X", "Y", "Z"], n_components=3)
    out = ipls(m, cal, intervals=4)
    assert len(out) == 4
    ga = ga_select(m, cal, intervals=6, population=6, generations=2)
    assert ga["ranges"]


def test_brereton_design_orthogonal():
    d = brereton_design(5)
    assert d.shape == (25, 5)
    corr = np.corrcoef(d.T)
    assert np.allclose(corr - np.eye(5), 0, atol=1e-12)
    for j in range(5):
        assert sorted(np.unique(d[:, j], return_counts=True)[1]) == [5] * 5


# --------------------------------------------------------------------------- #
# operations & regression
# --------------------------------------------------------------------------- #
def test_every_operation_runs(pure):
    s = mixture({"X": 10, "Y": 10})
    res = resolver(pure)
    for key, op in REGISTRY.items():
        params = op.defaults()
        for p in op.params:
            if p.kind == "spectrum":
                params[p.name] = "pure:X" if (key, p.name) in (
                    ("constant_center", "reference"),
                    ("divide_centered_ratio", "divisor")) else "pure:Y"
        if key == "normalize":
            params["mode"] = "max"
        if key == "t_to_a":
            s2 = s.copy(values=10 ** (-s.values) * 100)
            out = apply_pipeline(s2, [{"op": key, "params": params}], res)
            assert np.allclose(out.values, s.values)
            continue
        out = apply_pipeline(s, [{"op": key, "params": params}], res)
        assert np.all(np.isfinite(out.values)), key
        assert describe_step({"op": key, "params": params})


def test_derivative_of_line_is_slope():
    x = np.arange(200, 300, 1.0)
    s = Spectrum(x, 0.01 * x)
    for method in ("difference", "savitzky_golay"):
        d = apply_pipeline(s, [{"op": "derivative", "params": {"order": 1, "method": method,
                                                               "delta_lambda": 4, "window": 7}}])
        assert np.allclose(d.values[10:-10], 0.01)


def test_regression_statistics():
    x = [2, 4, 6, 8, 10]
    y = [0.101, 0.199, 0.302, 0.398, 0.501]
    r = linear_regression(x, y)
    assert r.r > 0.999 and abs(r.slope - 0.05) < 1e-3
    assert r.lod < r.loq
    assert abs(float(r.predict_x(0.25)) - 5) < 0.05


# --------------------------------------------------------------------------- #
# added methods: DAD, ERS, IAM, AAS, standard addition, robustness
# --------------------------------------------------------------------------- #
def test_dual_amplitude_difference_ternary(pure):
    lib = dict(pure)
    res = lambda k: lib[k]  # noqa: E731
    # ratio spectra by Z; choose λ pair where Y/Z has equal amplitude
    yz = apply_pipeline(pure["Y"], [{"op": "crop", "params": {"start": 215.0, "end": 300.0}},
                                    {"op": "divide", "params": {"reference": "Z"}}], res)
    pairs = uv.equal_amplitude_wavelengths(yz, 240.0, 220, 290)
    assert pairs
    steps = [{"op": "crop", "params": {"start": 215.0, "end": 300.0}},
             {"op": "divide", "params": {"reference": "Z"}}]
    m = uv.UnivariateMethod("DAD", "X", steps, {"kind": "weighted_difference", "params": {
        "w1": 240.0, "w2": pairs[0], "reference": "Y"}})
    m.calibrate([mixture({"X": c}) for c in (4, 8, 12, 16, 20)], res)
    for a, b, c in [(6, 10, 8), (12, 6, 14), (16, 12, 5)]:
        assert close(m.predict(mixture({"X": a, "Y": b, "Z": c}), res), a)


def test_extended_ratio_subtraction_and_double_divisor(binary_set, pure):
    xs, ys, mixes = binary_set
    res = resolver(pure)
    ers = uv.UnivariateMethod("ERS", "Y", [{"op": "extended_ratio_subtraction", "params": {
        "divisor": "pure:Y", "reference": "pure:X", "start": 340.0, "end": 360.0}}],
        {"kind": "amplitude", "params": {"w1": 275.0}})
    ers.calibrate(ys, res)
    for mix in mixes:
        assert close(ers.predict(mix, res), mix.concentrations["Y"])
    out = apply_pipeline(mixes[0], [{"op": "divide_sum", "params": {
        "reference": "pure:Y", "reference2": "pure:Z"}}], res)
    assert np.all(np.isfinite(out.values))


def test_induced_amplitude_modulation_and_aas(binary_set, pure):
    xs, ys, mixes = binary_set
    div = uv.InducedAmplitudeModulation.unit_divisor(ys, "Y")
    iam = uv.InducedAmplitudeModulation(250.0, (340.0, 360.0))
    iam.fit(xs, "X", div)
    for mix in mixes:
        r = iam.predict(mix, div)
        assert close(r["X"], mix.concentrations["X"]) and close(r["Y"], mix.concentrations["Y"])
    w1 = uv.equal_amplitude_wavelengths(pure["Y"], 245.0, 200, 400)[0]
    aas = uv.AdvancedAbsorbanceSubtraction(w1, 245.0)
    aas.fit(xs, "X", ys, "Y")
    for mix in mixes:
        r = aas.predict(mix)
        assert close(r["X"], mix.concentrations["X"]) and close(r["Y"], mix.concentrations["Y"])


def test_standard_addition_and_robustness(binary_set, pure):
    xs, ys, mixes = binary_set
    res = resolver(pure)
    m = uv.UnivariateMethod("RD", "X", [{"op": "divide", "params": {"reference": "pure:Y"}}],
                            {"kind": "difference", "params": {"w1": 240.0, "w2": 260.0}})
    m.calibrate(xs, res)
    added = [0, 0, 2, 4, 6]
    spiked = [mixture({"X": 5 + a, "Y": 8}) for a in added]
    sa = uv.standard_addition_recovery(lambda s: m.predict(s, res), spiked, added)
    assert abs(sa["sample"] - 5) < 0.05 and abs(sa["mean_recovery"] - 100) < 1
    assert abs(sa["extrapolated"] - 5) < 0.1
    params = uv.numeric_parameters(m)
    assert [p["label"] for p in params][-2:] == ["Measurement λ1", "Measurement λ2"]
    rob = uv.robustness_study(m, xs, mixes, params[-2:], res)
    assert len(rob["variants"]) == 5 and rob["max_abs_deviation"] < 2


# --------------------------------------------------------------------------- #
# chemometrics extras: frozen models, limits, outliers, progress
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("mtype,opts", [("CLS", {}), ("ILS", {}), ("PCR", {}), ("PLS1", {}),
                                        ("PLS2", {}), ("MCR-ALS", {}), ("ANN", {}),
                                        ("SVR", {"kernel": "linear"}),
                                        ("SVR", {"kernel": "rbf", "C": 1000.0})])
def test_model_roundtrip_predicts_identically(mtype, opts):
    import json
    cal = _design_set()
    test = [mixture({"X": 8, "Y": 12, "Z": 9}), mixture({"X": 13, "Y": 7, "Z": 11})]
    m = SpectralModel(mtype, ["X", "Y", "Z"], ranges=[(210, 370)], n_components=3,
                      wavelengths=[225, 245, 260, 275, 330], options=opts)
    m.fit(cal)
    before = m.predict(test)
    d = json.loads(json.dumps(m.to_dict(include_fit=True)))
    m2 = SpectralModel.from_dict(d)
    assert m2.is_fitted
    assert np.allclose(m2.predict(test), before, rtol=1e-9, atol=1e-9)


def test_limits_and_outlier_flags():
    cal = _design_set()
    m = SpectralModel("PLS2", ["X", "Y", "Z"], n_components=3).fit(cal)
    d = m.diagnostics(cal)
    lim = d["limits"]
    assert lim["t2_lim99"] > lim["t2_lim95"] > 0 and lim["q_lim99"] > lim["q_lim95"] > 0
    # an unexpected fourth absorber is flagged by Q; a wrong reference value by residual
    bad = mixture({"X": 10, "Y": 10, "Z": 10})
    bad = bad.copy(values=bad.values + 0.3 * np.exp(-0.5 * ((GRID - 300) / 6) ** 2))
    wrong = mixture({"X": 10, "Y": 10, "Z": 10})
    wrong.concentrations["X"] = 14.0
    flags = m.diagnostics([bad, wrong])["flags"]
    assert "Q" in flags[0] and "X residual" in flags[1]
    assert m.vip() is not None


def test_progress_and_cancel():
    cal = _design_set()
    m = SpectralModel("PLS2", ["X", "Y", "Z"], n_components=3)
    calls = []
    m.cross_validate(cal, method="venetian", folds=5, max_components=3,
                     progress=lambda d, t: calls.append((d, t)))
    assert calls[-1] == (15, 15)
    with pytest.raises(InterruptedError):
        m.cross_validate(cal, method="venetian", folds=5, max_components=3,
                         progress=lambda d, t: d < 4)
    with pytest.raises(InterruptedError):
        ipls(m, cal, intervals=4, progress=lambda d, t: False)


def test_q_limit_flags_few_calibration_samples():
    """At 95 % about 5 % of normal calibration samples may exceed the limits."""
    cal = [mixture(c, f"D{i}", 0.001, seed=100 + i) for i, c in enumerate(
        design_concentrations({"X": 10, "Y": 10, "Z": 10}, {"X": 4, "Y": 4, "Z": 4}))]
    m = SpectralModel("PLS2", ["X", "Y", "Z"], n_components=3).fit(cal)
    d = m.diagnostics(cal)
    q_out = np.mean(np.array(d["Q"]) > d["limits"]["q_lim95"])
    assert q_out <= 0.15
    assert d["limits"]["q_lim95"] > np.mean(d["Q"])
