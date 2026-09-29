import math

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
