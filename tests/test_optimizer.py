import math

import numpy as np

from spectro.core.optimizer import OptimizerInput, materialize, optimize
from spectro.core.spectrum import Spectrum
from tests.conftest import GRID, PURE, mixture


def _std(c, comps, levels=(4, 8, 12, 16, 20), noise=0.0005):
    return [mixture({k: (v if k == c else 0.0) for k, v in zip(comps, [lv] * len(comps))},
                    f"{c} {lv}", noise, seed=hash((c, lv)) % 1000) for lv in levels]


def test_binary_finds_accurate_methods_and_rejects_bad_ones():
    comps = ["X", "Y"]
    std = {c: _std(c, comps) for c in comps}
    div = {c: std[c][2] for c in comps}
    mixes = [mixture({"X": a, "Y": b}, f"m{a}{b}", 0.0005, seed=a * 10 + b)
             for a, b in [(6, 10), (10, 6), (14, 14), (8, 16), (18, 5)]]
    res = optimize(OptimizerInput(comps, std, div, mixes))
    for c in comps:
        ranked = res["ranked"][c]
        best = ranked[0]
        assert best.score < 2.0, (c, best.label, best.score)
        assert abs(best.real_mean_recovery - 100) < 2 or math.isnan(best.real_mean_recovery)
        families = {x.family for x in ranked}
        assert {"ratio_difference", "derivative", "vierordt", "cls"} <= families
        # plain zero order suffers from overlap and must rank below the best
        zero = next(x for x in ranked if x.family == "zero")
        assert zero.score > best.score
    # a univariate winner can be turned into a concrete method definition
    uni = next(x for x in res["ranked"]["X"] if x.measurement)
    d = materialize(uni, {"X": 1, "Y": 2})
    assert all(not str(v).startswith("div:") for st in d["steps"] for v in st["params"].values())


def test_ternary_runs_with_ternary_families():
    comps = ["X", "Y", "Z"]
    std = {c: _std(c, comps) for c in comps}
    div = {c: std[c][2] for c in comps}
    res = optimize(OptimizerInput(comps, std, div, [], families={
        "zero", "derivative", "dual_amplitude", "double_divisor", "vierordt", "cls"}))
    for c in comps:
        ranked = res["ranked"][c]
        assert {"dual_amplitude", "cls"} <= {x.family for x in ranked}
        assert ranked[0].score < 3.0, (c, ranked[0].label, ranked[0].score)


def test_unit_spectra_recover_absorptivity():
    from spectro.core.optimizer import unit_spectra
    comps = ["X", "Y"]
    inp = OptimizerInput(comps, {c: _std(c, comps, noise=0.0) for c in comps},
                         wl_range=(210, 390))
    grid, units, _ = unit_spectra(inp)
    import numpy as np
    ref = np.interp(grid, GRID, PURE["X"])
    assert np.allclose(units["X"], ref, atol=1e-6)
    assert isinstance(Spectrum(grid, units["Y"]), Spectrum)


# --------------------------------------------------------------------------- #
# Progressive methods in the optimizer
# --------------------------------------------------------------------------- #
def _lit_input(pure_mk, comps, mixes, noise=0.0003):
    std = {c: [pure_mk({c: v}, f"{c} {v}", noise=noise, seed=int(v) * 7 + i)
               for i, v in enumerate((4, 8, 12, 16, 20, 24))] for c in comps}
    div = {c: std[c][3] for c in comps}
    mx = [pure_mk(m, f"m{i}", noise=noise, seed=100 + i) for i, m in enumerate(mixes)]
    return OptimizerInput(comps, std, div, mx)


def test_optimizer_ranks_amplitude_centering_for_extended_divisor():
    """Ternary with an extended component (AAC / RIDSS situation): the
    optimizer builds amplitude-centering methods by itself, verifies them on
    the laboratory mixtures, and they recover all three compounds."""
    from tests.test_literature import ter
    inp = _lit_input(ter, ["X", "Y", "Z"],
                     [{"X": 10, "Y": 10, "Z": 10}, {"X": 20, "Y": 10, "Z": 10},
                      {"X": 6, "Y": 6, "Z": 18}, {"X": 12, "Y": 4, "Z": 4}])
    inp.families = {"amplitude_centering"}
    res = optimize(inp)
    for c in inp.compounds:
        prog = [x for x in res["ranked"][c] if x.family == "amplitude_centering"
                and not x.error]
        assert prog, c
        best = prog[0]
        assert best.score < 10, (c, best.label, best.score)
        assert abs(best.real_mean_recovery - 100) < 5
        assert best.model["type"] == "amplitude_centering"
        assert best.model["divisor"].startswith("div:")
    # both structures were generated: plateau (AAC-partial) and λ pairs only
    labels = {x.label for c in inp.compounds for x in res["ranked"][c]}
    assert any("plateau" in lab for lab in labels)
    assert any("plateau" not in lab and "by subtraction" in lab for lab in labels)


def test_optimizer_ranks_absorption_factor_and_materializes_it():
    from spectro.core.univariate import progressive_from_dict
    from tests.test_literature import af
    inp = _lit_input(af, ["MET", "MEB", "DLX"],
                     [{"MET": 10, "MEB": 10, "DLX": 10}, {"MET": 20, "MEB": 4, "DLX": 2},
                      {"MET": 3, "MEB": 20, "DLX": 15}])
    inp.families = {"absorption_factor", "zero"}
    res = optimize(inp)
    best = next(x for x in res["ranked"]["MET"] if x.family == "absorption_factor")
    assert best.score < 2
    order = [c for c, _ in best.model["order"]]
    assert order[0] == "MET"        # the only compound absorbing alone (above 320 nm)
    d = materialize(best, {c: 1 for c in inp.compounds})
    m = progressive_from_dict(d["model"])
    m.fit_spectra([s for c in inp.compounds for s in inp.standards[c]])
    out = m.predict_spectra(inp.mixtures)
    true = [[x.concentrations[c] for c in m.compounds] for x in inp.mixtures]
    assert np.allclose(out, true, rtol=0.05)


def test_materialize_replaces_the_progressive_divisor():
    from spectro.core.optimizer import Candidate
    from spectro.core.univariate import AmplitudeCentering
    m = AmplitudeCentering(260.0, ["X", "Z"], subtract="Z", divisor_compound="Z",
                           differences={"X": {"w1": 260.0, "w2": 240.0}}, divisor="div:Z")
    d = materialize(Candidate("X", "amplitude_centering", "", model=m.to_dict()),
                    {"X": 3, "Z": 7})
    assert d["model"]["divisor"] == 7


# --------------------------------------------------------------------------- #
# Smoothing options
# --------------------------------------------------------------------------- #
def test_sg_window_converts_nm_to_odd_points():
    from spectro.core.optimizer import sg_window
    assert sg_window(6, 0.5) == 13
    assert sg_window(10, 1.0) == 11
    assert sg_window(1, 1.0) == 5          # never below 5 points
    assert all(sg_window(w, s) % 2 == 1 for w in (2, 3, 7.3) for s in (0.1, 0.5, 2))


def _noisy_binary(noise=0.003):
    comps = ["X", "Y"]
    std = {c: [mixture({k: (v if k == c else 0.0) for k in comps}, f"{c} {v}", noise,
                       seed=int(v) * 3 + len(c)) for v in (4, 8, 12, 16, 20)] for c in comps}
    mixes = [mixture({"X": a, "Y": b}, f"m{a}{b}", noise, seed=a * 10 + b)
             for a, b in [(6, 10), (10, 6), (14, 14), (8, 16), (18, 5)]]
    return comps, std, mixes


def test_smoothing_variants_are_screened_and_help_noisy_derivatives():
    """With noisy spectra (0.003 AU) smoothing must lower the error of the
    best derivative methods on the LABORATORY mixtures (independent noise),
    not only in the prediction."""
    from spectro.core.univariate import UnivariateMethod
    comps, std, mixes = _noisy_binary()
    fam = {"derivative", "derivative_ratio", "ratio_difference"}
    plain = optimize(OptimizerInput(comps, std, {c: std[c][2] for c in comps}, mixes,
                                    families=fam))
    smooth = optimize(OptimizerInput(comps, std, {c: std[c][2] for c in comps}, mixes,
                                     families=fam, smoothing=(3.0, 6.0, 10.0)))
    assert not any(st["op"] == "smooth_sg" for c in comps for x in plain["ranked"][c]
                   for st in x.steps)
    for c in comps:
        for f in fam:
            best_plain = next(x for x in plain["ranked"][c] if x.family == f and not x.error)
            best_smooth = next(x for x in smooth["ranked"][c] if x.family == f and not x.error)
            assert best_smooth.real_rmsep <= best_plain.real_rmsep + 1e-9, (c, f)
        ranked = smooth["ranked"][c]
        assert any(x.steps and x.steps[0]["op"] == "smooth_sg" for x in ranked)
        sgd = [x for x in ranked if any(st["op"] == "derivative" and
                                        st["params"].get("method") == "savitzky_golay"
                                        for st in x.steps)]
        assert sgd and all("S-G" in x.label for x in sgd)
    # the derivative winner improves a lot here
    for c in comps:
        bp = next(x for x in plain["ranked"][c] if x.family == "derivative" and not x.error)
        bs = next(x for x in smooth["ranked"][c] if x.family == "derivative" and not x.error)
        assert bs.real_rmsep < 0.6 * bp.real_rmsep
    # a smoothed winner becomes a working method
    best = next(x for x in smooth["ranked"]["X"] if x.measurement and x.steps
                and x.steps[0]["op"] == "smooth_sg")
    d = materialize(best, {"X": 1, "Y": 2})
    lib = {1: std["X"][2], 2: std["Y"][2]}
    m = UnivariateMethod("smoothed", "X", d["steps"], d["measurement"])
    m.calibrate(std["X"], lib.__getitem__)
    rec = [100 * m.predict(s, lib.__getitem__) / s.concentrations["X"] for s in mixes]
    assert all(abs(r - 100) < 5 for r in rec)


def test_smoothing_also_applies_to_multicomponent_models():
    comps, std, mixes = _noisy_binary()
    res = optimize(OptimizerInput(comps, std, {c: std[c][2] for c in comps}, mixes,
                                  families={"vierordt", "cls"}, smoothing=(6.0,)))
    labels = [x.label for x in res["ranked"]["X"]]
    assert any(lab.startswith("Classical least squares") and "SG smoothing 6 nm" in lab
               for lab in labels)
    assert any(lab.startswith("Vierordt") and "SG smoothing 6 nm" in lab for lab in labels)
    smoothed = next(x for x in res["ranked"]["X"] if "SG smoothing" in x.label)
    assert smoothed.model["steps"][0]["op"] == "smooth_sg"


def test_smoothing_reaches_progressive_methods():
    from tests.test_literature import af
    comps = ["MET", "MEB", "DLX"]
    std = {c: [af({c: v}, f"{c}{v}", noise=0.002, seed=int(v) * 7 + i)
               for i, v in enumerate((4, 8, 12, 16, 20, 24))] for c in comps}
    mixes = [af(m, f"m{i}", noise=0.002, seed=100 + i) for i, m in enumerate(
        [{"MET": 10, "MEB": 10, "DLX": 10}, {"MET": 20, "MEB": 4, "DLX": 2},
         {"MET": 3, "MEB": 20, "DLX": 15}])]
    res = optimize(OptimizerInput(comps, std, {c: std[c][3] for c in comps}, mixes,
                                  families={"absorption_factor"}, smoothing=(6.0, 10.0)))
    cands = [x for x in res["ranked"]["DLX"] if x.family == "absorption_factor"]
    sm = [x for x in cands if x.model.get("steps")]
    plain = [x for x in cands if not x.model.get("steps")]
    assert sm and plain
    assert sm[0].model["steps"][0]["op"] == "smooth_sg" and "SG smoothing" in sm[0].label
    assert min(x.real_rmsep for x in sm) < min(x.real_rmsep for x in plain)


def test_parallel_optimizer_gives_the_serial_result():
    comps, std, mixes = _noisy_binary()
    inp = OptimizerInput(comps, std, {c: std[c][2] for c in comps}, mixes,
                         families={"derivative", "ratio_difference"}, smoothing=(6.0,))
    a = optimize(inp)
    b = optimize(inp, workers=2)
    for c in comps:
        assert [(x.label, x.score) for x in a["ranked"][c]] == \
            [(x.label, x.score) for x in b["ranked"][c]]
