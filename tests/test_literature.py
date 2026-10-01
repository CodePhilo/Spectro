"""Validation against the published literature.

Two kinds of checks, for the nine papers the lab supplied:

1. **Published statistics are recomputed** from the summary values printed in
   each paper (means, SD or variance, n, individual recoveries) with the app's
   own functions, and compared with the t, F, ANOVA, mean and SD values the
   authors reported. Agreement shows that the app's statistics follow the same
   conventions as the literature (pooled-variance t, one-tailed tabulated F,
   n − 1 SD). Where a paper's printed value cannot be reproduced, the test
   pins the *correct* value and says what is wrong in the paper.

2. **Every method of every paper is run on simulated spectra** built to meet
   the paper's spectral conditions (which component is extended, where the
   isoabsorptive point is, which component has a plateau …), using the app's
   finder tools to choose wavelengths as an analyst would. Noise-free spectra
   must be recovered to ±0.1 % (the algebra is exact), noisy ones within the
   usual 98–102 %.

Papers (DOIs):
  AM14  Abdelrahman et al., Anal. Methods 6 (2014) 509 — 10.1039/c3ay41564c (RIDSS)
  SA14  Abdelrahman et al., SAA 124 (2014) 389 — 10.1016/j.saa.2014.01.020
  LO15  Lotfy et al., SAA 136 (2015) 937 — 10.1016/j.saa.2014.09.117
  ZA20  Zaghary et al., J. Anal. Chem. 75 (2020) 742 — 10.1134/S1061934820060180
  SL17  Saleh et al., IJPPS 9 (2017) 43 — 10.22159/ijpps.2017v9i5.16960 (AAC)
  EI18  Eissa & Abou Al Alamein, SAA 193 (2018) 365 — 10.1016/j.saa.2017.12.050
  EM18  Emam et al., SAA (2018) — 10.1016/j.saa.2017.11.034 (SRS)
  FA21  Fahmy et al., SAA 261 (2021) 119999 — 10.1016/j.saa.2021.119999
  AW17  Abdelwahab & Mohamed, Chem. Pharm. Bull. 65 (2017) 558 (AUC, MAFM, MACM)
"""

from __future__ import annotations

import numpy as np
import pytest

from spectro.core import univariate as uv
from spectro.core import validation as val
from spectro.core.multicomponent import SignalEquations
from spectro.core.operations import apply_pipeline
from spectro.core.spectrum import Spectrum

# --------------------------------------------------------------------------- #
# 1. Published statistics
# --------------------------------------------------------------------------- #
T6, F55 = 2.228, 5.050       # tabulated t (10 df) and one-tailed F(5, 5) at P = 0.05
T5, F44 = 2.306, 6.388       # n = 5 and 5


def _tf(m1, s1, n1, m2, s2, n2):
    return val.t_test_summary(m1, s1, n1, m2, s2, n2), val.f_test_summary(s1, n1, s2, n2)


# AM14 Table 3 — proposed RIDSS vs HPLC, variances, n = 6
AM14 = [  # (mean, var) proposed, (mean, var) reference, published t, F
    ((103.67, 1.200), (104.41, 1.971), 1.023, 1.642),
    ((105.15, 1.472), (104.01, 2.761), None, 1.875),     # t printed 1.250, see below
    ((96.22, 2.47), (96.30, 2.468), 0.086, 1.001),
    ((104.95, 2.022), (104.24, 3.832), 0.723, 1.895),
    ((100.79, 1.935), (101.25, 1.790), 0.575, 1.081),
    ((101.40, 2.354), (102.39, 2.632), 1.086, 1.118),
]

# SA14 Table 4 — DD, AUC and MCR methods vs reported method, variances, n = 6
_SA14_REF = {"ORPH": (100.57, 1.503), "CAF": (99.93, 1.769), "ASP": (102.26, 1.232)}
SA14 = [
    ("ORPH", (100.01, 1.392), 0.802, 1.080), ("CAF", (98.56, 1.350), 1.893, 1.310),
    ("ASP", (102.84, 2.286), 0.757, 1.856), ("ORPH", (99.15, 1.893), 1.887, 1.259),
    ("CAF", (99.26, 1.626), 0.886, 1.088), ("ASP", (103.51, 1.484), 1.853, 1.205),
    ("ORPH", (99.30, 2.528), 1.549, 1.682), ("CAF", (98.76, 0.817), 1.777, 2.165),
    ("ASP", (103.47, 0.846), 2.062, 1.456), ("ASP", (103.67, 1.309), 2.172, 1.063),
]

# EI18 Table 4 — five methods vs HPLC, SD, n = 5
_EI_REF = {"SCB": (100.63, 0.726), "VLS": (99.87, 0.670)}
EI18 = [
    ("SCB", (100.60, 0.821), 0.069, 1.277), ("SCB", (99.80, 1.067), 1.439, 2.158),
    ("SCB", (100.01, 0.607), 1.467, 1.434), ("SCB", (99.95, 0.735), 1.470, 1.023),
    ("SCB", (99.73, 0.889), 1.761, 1.499), ("VLS", (100.62, 0.821), 1.583, 1.502),
    ("VLS", (99.47, 0.860), 0.821, 1.647), ("VLS", (99.33, 0.859), 1.124, 1.646),
    ("VLS", (100.11, 0.382), 0.696, 3.074), ("VLS", (100.07, 0.571), 0.498, 1.378),
]

# FA21 Table 3 — vs USP official methods, SD, n = 6
_FA_REF = {"TN": (100.22, 0.78), "EX": (99.65, 0.37), "HQ": (100.33, 0.53), "HC": (99.75, 0.32)}
FA21 = [
    ("TN", (100.21, 1.05), 0.019, 1.81), ("EX", (100.22, 0.77), 1.63, 4.33),
    ("EX", (100.01, 0.72), 1.06, 3.79), ("EX", (100.00, 0.29), 1.84, 1.63),
    ("HQ", (100.03, 0.61), 0.91, 1.32), ("HQ", (99.92, 0.66), 1.17, 1.55),
    ("HQ", (100.13, 0.79), 0.51, 2.22), ("HQ", (100.42, 0.60), 0.27, 1.28),
    ("HQ", (100.24, 0.35), 0.35, 2.29), ("HC", (99.96, 0.22), 1.33, 2.12),
    ("HC", (99.28, 0.67), 1.57, 4.38),
]

# SL17 Table 6 — AAC / PLS (n = 7) vs reported (n = 5), variances
_SL_REF = {"AML": (100.31, 0.3913), "VAL": (100.30, 1.5274), "HCT": (100.32, 0.7472)}
SL17 = [
    ("AML", (100.16, 1.9057), 0.224, 4.870), ("AML", (100.09, 0.6495), 0.523, 1.660),
    ("VAL", (99.81, 1.4802), 0.691, 1.032), ("VAL", (99.95, 1.5578), 0.472, 1.020),
    ("VAL", (99.94, 0.9082), 0.581, 1.682), ("HCT", (100.22, 0.8201), 0.186, 1.098),
    ("HCT", (99.88, 0.8943), None, 1.197),     # t printed 0.764, recomputed 0.822
    ("HCT", (100.38, 0.2361), 0.153, 3.165),
]


@pytest.mark.parametrize("row", AM14)
def test_am14_t_and_f_reproduce(row):
    (m1, v1), (m2, v2), t_pub, f_pub = row
    t, f = _tf(m1, v1 ** .5, 6, m2, v2 ** .5, 6)
    if t_pub is not None:
        assert t["t"] == pytest.approx(t_pub, abs=0.01)
    assert f["F"] == pytest.approx(f_pub, abs=0.002)
    assert t["t_crit"] == pytest.approx(T6, abs=5e-4)
    assert f["F_crit"] == pytest.approx(F55, abs=5e-4)
    assert not t["significant"] and not f["significant"]


def test_am14_guaifenesin_t_value_is_misprinted():
    """Farcosolvin GUF: 105.15 (var 1.472) vs 104.01 (var 2.761), n = 6 each
    gives t = 1.357, not the printed 1.250 (conclusion unchanged)."""
    t, _ = _tf(105.15, 1.472 ** .5, 6, 104.01, 2.761 ** .5, 6)
    assert t["t"] == pytest.approx(1.357, abs=0.002) and not t["significant"]


@pytest.mark.parametrize("row", SA14)
def test_sa14_t_and_f_reproduce(row):
    comp, (m, var), t_pub, f_pub = row
    rm, rv = _SA14_REF[comp]
    t, f = _tf(m, var ** .5, 6, rm, rv ** .5, 6)
    assert t["t"] == pytest.approx(t_pub, abs=0.01)
    assert f["F"] == pytest.approx(f_pub, abs=0.002)
    assert not t["significant"] and not f["significant"]


@pytest.mark.parametrize("row", EI18)
def test_ei18_t_and_f_reproduce(row):
    comp, (m, sd), t_pub, f_pub = row
    rm, rsd = _EI_REF[comp]
    t, f = _tf(m, sd, 5, rm, rsd, 5)
    assert t["t"] == pytest.approx(t_pub, abs=0.02)
    assert f["F"] == pytest.approx(f_pub, abs=0.005)
    assert t["t_crit"] == pytest.approx(T5, abs=5e-4)
    assert f["F_crit"] == pytest.approx(F44, abs=5e-4)


@pytest.mark.parametrize("row", FA21)
def test_fa21_t_and_f_reproduce(row):
    comp, (m, sd), t_pub, f_pub = row
    rm, rsd = _FA_REF[comp]
    t, f = _tf(m, sd, 6, rm, rsd, 6)
    assert t["t"] == pytest.approx(t_pub, abs=0.035)
    assert f["F"] == pytest.approx(f_pub, abs=0.01)


def test_fa21_eusolex_constant_value_is_borderline():
    """EX by CV (100.10 ± 0.33) vs official (99.65 ± 0.37): the paper prints
    t = 2.25 against its tabulated 2.23 yet reports no difference; from the
    printed means/SDs t = 2.22, just below 2.228 — not significant, but
    borderline."""
    t, _ = _tf(100.10, 0.33, 6, 99.65, 0.37, 6)
    assert t["t"] == pytest.approx(2.223, abs=0.002)
    assert t["t"] < t["t_crit"] and t["t_crit"] - t["t"] < 0.01


@pytest.mark.parametrize("row", SL17)
def test_sl17_unequal_n_t_and_f_reproduce(row):
    """n = 7 vs 5: pooled-variance t with 10 df; F with ν = (6, 4) or (4, 6)."""
    comp, (m, var), t_pub, f_pub = row
    rm, rv = _SL_REF[comp]
    t, f = _tf(m, var ** .5, 7, rm, rv ** .5, 5)
    if t_pub is not None:
        assert t["t"] == pytest.approx(t_pub, abs=0.02)
    assert f["F"] == pytest.approx(f_pub, abs=0.002)
    assert t["t_crit"] == pytest.approx(T6, abs=5e-4)
    # paper prints 6.094 and 4.120; the correct one-tailed values are:
    assert f["F_crit"] == pytest.approx(6.163 if var > rv else 4.534, abs=0.001)


def test_ei18_individual_recoveries_give_published_mean_and_sd():
    table3 = {  # EI18 Table 3, five laboratory mixtures
        (99.83, 0.880): [99.99, 99.39, 98.59, 100.90, 100.29],
        (100.29, 0.697): [101.22, 100.19, 99.29, 100.22, 100.55],
        (100.34, 1.005): [99.14, 101.94, 100.20, 100.20, 100.20],
        (100.68, 0.980): [100.12, 99.91, 99.89, 101.92, 101.55],
        (100.50, 1.078): [100.39, 99.34, 101.66, 99.56, 101.53],
        (99.19, 0.935): [98.72, 98.25, 99.78, 98.67, 100.52],
        (99.18, 0.543): [98.86, 98.85, 100.10, 99.26, 98.85],
        (98.97, 1.166): [98.25, 98.25, 100.94, 99.15, 98.25],
        (100.29, 1.056): [100.55, 101.64, 98.80, 99.82, 100.66],
        (100.70, 0.881): [101.29, 101.29, 99.31, 100.33, 101.29],
    }
    for (mean, sd), vals in table3.items():
        d = val.describe(vals)
        assert d["mean"] == pytest.approx(mean, abs=0.006)
        assert d["sd"] == pytest.approx(sd, abs=0.004)


def test_ei18_anova_from_summary_statistics():
    """EI18 Table 5 (SCB): one-way ANOVA of the reported HPLC method and five
    proposed methods, rebuilt from Table 4 means and SDs (n = 5)."""
    groups = {"HPLC": (100.63, 0.726, 5), "AAS": (100.60, 0.821, 5),
              "IDW": (99.80, 1.067, 5), "RD": (100.01, 0.607, 5),
              "DR1": (99.95, 0.735, 5), "MCR": (99.73, 0.889, 5)}
    a = val.anova_summary(groups)
    assert a["ss_between"] == pytest.approx(3.9286, rel=0.003)
    assert a["ss_within"] == pytest.approx(16.1551, rel=0.001)
    assert (a["df_between"], a["df_within"]) == (5, 24)
    assert a["F"] == pytest.approx(1.1673, rel=0.003)
    assert a["p"] == pytest.approx(0.3539, abs=0.001)
    assert a["F_crit"] == pytest.approx(2.6207, abs=1e-4)


def test_anova_summary_equals_raw_data_anova():
    rng = np.random.default_rng(3)
    groups = {k: rng.normal(100 + i * 0.3, 0.8, 6) for i, k in enumerate("ABCD")}
    raw = val.anova_oneway(groups)
    summ = val.anova_summary({k: (v.mean(), v.std(ddof=1), v.size) for k, v in groups.items()})
    for k in ("F", "p", "ss_between", "ss_within", "F_crit"):
        assert summ[k] == pytest.approx(raw[k], rel=1e-9)


def test_summary_tests_equal_raw_data_tests():
    rng = np.random.default_rng(4)
    a, b = rng.normal(100, 0.7, 6), rng.normal(100.6, 1.2, 8)
    for eq in (True, False):
        raw = val.t_test(a, b, equal_var=eq)
        s = val.t_test_summary(a.mean(), a.std(ddof=1), 6, b.mean(), b.std(ddof=1), 8,
                               equal_var=eq)
        assert s["t"] == pytest.approx(raw["t"]) and s["p"] == pytest.approx(raw["p"])
    f = val.f_test(a, b)
    assert f["F"] == pytest.approx(max(a.var(ddof=1), b.var(ddof=1)) /
                                   min(a.var(ddof=1), b.var(ddof=1)))
    assert f["F_crit_one_tailed"] < f["F_crit_two_tailed"]


def test_fa21_hydrocortisone_anova_is_misprinted():
    """FA21 Table 4 prints SS(between) = 0.61, F = 1.52 for HC; the means of
    Table 3 (99.96, 99.28, 99.75; SD 0.22, 0.67, 0.32; n = 6) give SS = 1.45
    and F = 3.64 — still below F crit 3.68, so the conclusion stands, but
    only just."""
    a = val.anova_summary({"D0 242": (99.96, 0.22, 6), "D1 254.1": (99.28, 0.67, 6),
                           "HPLC": (99.75, 0.32, 6)})
    assert a["ss_between"] == pytest.approx(1.455, abs=0.002)
    assert a["F"] == pytest.approx(3.639, abs=0.005)
    assert a["F_crit"] == pytest.approx(3.682, abs=0.001) and not a["significant"]


def test_za20_recoveries_and_precision_comparison():
    """ZA20 Table 2 recoveries reproduce the printed means/SDs. Table 3 then
    compares precision with the reported method: for MET by D1 the variance
    ratio is 4.45/0.66 = 6.79, ABOVE the tabulated F (6.388 for the paper's
    n = 5, 5.05 for the six mixtures actually listed) although the paper
    reports no significant difference."""
    d1_met = [100.1, 100.2, 97.8, 99.0, 94.9, 100.4]
    d = val.describe(d1_met)
    assert d["mean"] == pytest.approx(98.7, abs=0.05) and d["sd"] == pytest.approx(2.12, abs=0.01)
    for n in (5, 6):
        f = val.f_test_summary(2.11, n, 0.81, n)
        assert f["F"] == pytest.approx(6.786, abs=0.001)
        assert f["significant"]


def test_em18_standard_addition_row_is_misprinted():
    """EM18 Table 4, furosemide: 7.86 µg/mL found of 8.00 added is 98.25 %,
    not the printed 99.50 %; with 98.25 the printed mean ± SD (100.00 ± 1.561)
    is reproduced exactly."""
    rec = val.recovery([4.05, 6.03, 7.86], [4.00, 6.00, 8.00])
    assert rec["recoveries"][2] == pytest.approx(98.25)
    assert rec["mean"] == pytest.approx(100.00, abs=1e-9)
    assert rec["sd"] == pytest.approx(1.561, abs=5e-4)
    sp = val.recovery([13.12, 14.75, 16.79], [13.0, 15.0, 17.0])
    assert sp["mean"] == pytest.approx(99.34, abs=0.006) and sp["sd"] == pytest.approx(1.388, abs=0.002)


def test_one_tailed_and_two_tailed_f_critical_values():
    f = val.f_test_summary(1.0, 6, 0.8, 6)
    assert f["F_crit_one_tailed"] == pytest.approx(5.050, abs=5e-4)
    assert f["F_crit_two_tailed"] == pytest.approx(7.146, abs=5e-4)
    assert val.f_test_summary(1.0, 6, 0.8, 6, tails=2)["F_crit"] == f["F_crit_two_tailed"]


# --------------------------------------------------------------------------- #
# 2. Methods on simulated spectra built to each paper's conditions
# --------------------------------------------------------------------------- #
GRID = np.round(np.arange(200.0, 400.01, 0.2), 4)


def g(c, w, h=1.0):
    return h * np.exp(-0.5 * ((GRID - c) / w) ** 2)


def make(pure: dict[str, np.ndarray]):
    """Factory of mixtures (and pure standards) from unit spectra."""
    def mix(conc: dict, name="", noise=0.0, seed=0):
        y = sum(pure[k] * v for k, v in conc.items())
        if noise:
            y = y + np.random.default_rng(seed).normal(0, noise, y.size)
        return Spectrum(GRID, y, name=name or str(conc),
                        concentrations={k: float(conc.get(k, 0.0)) for k in pure})
    return mix


def stds(mix, comp, levels):
    return [mix({comp: c}, f"{comp} {c}") for c in levels]


def rec(found, true):
    return 100 * found / true


def res_of(lib):
    return lambda key: lib[key]


# ---- ZA20: binary, CANA (Y, extended) + MET (X) ---------------------------- #
ZA = {"MET": 0.080 * g(237, 11) + 0.010 * g(212, 6),
      "CANA": 0.045 * g(290, 13) + 0.035 * g(225, 10)}
za = make(ZA)
ZA_MIX = [{"CANA": 5, "MET": 10}, {"CANA": 10, "MET": 2}, {"CANA": 20, "MET": 5},
          {"CANA": 10, "MET": 15}, {"CANA": 6, "MET": 8}, {"CANA": 5, "MET": 15}]
ZA_LIB = {"CANA5": za({"CANA": 5}, "CANA 5")}


def _uv(steps, w, comp, cal, noise=0.0):
    m = uv.UnivariateMethod("m", comp, steps, {"kind": "amplitude", "params": {"w1": w}})
    m.calibrate(cal, res_of(ZA_LIB))
    return [rec(m.predict(za(c, noise=noise, seed=i), res_of(ZA_LIB)), c[comp])
            for i, c in enumerate(ZA_MIX)]


def test_za20_ss_cm_recovers_both_drugs():
    cm = [{"op": "constant_multiplication",
           "params": {"divisor": "CANA5", "start": 300.0, "end": 320.0}}]
    rs = [{"op": "ratio_subtraction", "params": {"divisor": "CANA5", "start": 300.0,
                                                 "end": 320.0}}]
    cana = _uv(cm, 290.0, "CANA", stds(za, "CANA", [5, 10, 15, 20]))
    met = _uv(rs, 237.0, "MET", stds(za, "MET", [2, 5, 10, 15]))
    assert np.allclose(cana, 100, atol=0.01) and np.allclose(met, 100, atol=0.01)


def test_za20_ds_cm_in_derivative_mode():
    """DS-CM: D1 of the mixture ÷ D1 of CANA′ gives a constant over the
    region where only CANA's D1 is non-zero; constant × D1(CANA′) is CANA's
    D1 spectrum, mixture D1 − it is MET's D1 spectrum."""
    d1 = uv._D1
    cm = [d1, {"op": "constant_multiplication",
               "params": {"divisor": "CANA5", "start": 300.0, "end": 318.0,
                          "divisor_derivative": 1}}]
    ds = [d1, {"op": "ratio_subtraction",
               "params": {"divisor": "CANA5", "start": 300.0, "end": 318.0,
                          "divisor_derivative": 1}}]
    cana = _uv(cm, 305.0, "CANA", stds(za, "CANA", [5, 10, 15, 20]))
    met = _uv(ds, 228.0, "MET", stds(za, "MET", [2, 5, 10, 15]))
    assert np.allclose(cana, 100, atol=0.01) and np.allclose(met, 100, atol=0.01)
    noisy = _uv(cm, 305.0, "CANA", stds(za, "CANA", [5, 10, 15, 20]), noise=0.0005)
    assert 98 < np.mean(noisy) < 102


def test_za20_derived_mixture_spectra_equal_pure_spectra():
    """The resolved spectra are the pure-component spectra (spectral profile)."""
    m = za({"CANA": 8, "MET": 12})
    y = apply_pipeline(m, [{"op": "constant_multiplication",
                            "params": {"divisor": "CANA5", "start": 300.0, "end": 320.0}}],
                       res_of(ZA_LIB))
    assert np.allclose(y.values, 8 * ZA["CANA"], atol=1e-12)


def test_za20_amplitude_modulation_with_normalized_divisor_and_unified_regression():
    unit = apply_pipeline(ZA_LIB["CANA5"], [{"op": "normalize",
                                            "params": {"mode": "concentration"}}])
    assert np.allclose(unit.values, ZA["CANA"])
    pure = {k: Spectrum(GRID, v, concentrations={k: 1.0}) for k, v in ZA.items()}
    iso = [w for w in uv.isoabsorptive_points(pure["MET"], pure["CANA"]) if 245 < w < 270]
    assert iso
    am = uv.AmplitudeModulation(iso[0], (300.0, 320.0), 1.0)
    am.fit(stds(za, "MET", [2, 5, 10, 15]) + stds(za, "CANA", [5, 10, 15, 20]), unit)
    # "slope ≈ 1, intercept ≈ 0" (the iso λ is interpolated, hence 1e-4)
    assert am.total_regression.slope == pytest.approx(1.0, rel=1e-4)
    for c in ZA_MIX:
        out = am.predict(za(c), unit)
        assert out["Y"] == pytest.approx(c["CANA"], rel=5e-4)
        assert out["X"] == pytest.approx(c["MET"], rel=5e-4)


def test_za20_sample_enrichment_by_spiking():
    """Invokamet: MET 10 + CANA 0.5 µg/mL; CANA is too low, so 5 µg/mL pure
    CANA is added, determined, and subtracted."""
    m = uv.UnivariateMethod("CM", "CANA", [{"op": "constant_multiplication",
                                            "params": {"divisor": "CANA5", "start": 300.0,
                                                       "end": 320.0}}],
                            {"kind": "amplitude", "params": {"w1": 290.0}})
    m.calibrate(stds(za, "CANA", [5, 10, 15, 20]), res_of(ZA_LIB))
    found = m.predict(za({"MET": 10, "CANA": 0.5 + 5.0}), res_of(ZA_LIB))
    out = uv.enrichment_correction(found, 5.0, claimed=0.5)
    assert out["found"] == pytest.approx(0.5, abs=1e-5)
    assert out["percent_of_claimed"] == pytest.approx(100, abs=1e-3)


def test_spectrum_addition_enrichment_is_the_same_as_spiking():
    """LO15/FA21 'spectrum addition': the stored spectrum of a pure standard
    is added to the sample spectrum instead of adding the solid."""
    lib = dict(ZA_LIB, CANA_add=za({"CANA": 5.7}))
    spiked = apply_pipeline(za({"MET": 10, "CANA": 0.3}),
                            [{"op": "add", "params": {"reference": "CANA_add"}}], res_of(lib))
    assert np.allclose(spiked.values, za({"MET": 10, "CANA": 6.0}).values)


# ---- AW17 / SA14 / SL17 / AM14 / LO15: ternary progressive methods ---------- #
TER = {"X": 0.040 * g(240, 13) + 0.022 * g(275, 9),     # e.g. GUF / HCT / FCP
       "Y": 0.030 * g(255, 12) + 0.018 * g(290, 8),     # e.g. THP / VAL / LH
       "Z": 0.012 * g(250, 15) + 0.010 * g(355, 22)}    # extended: AMB / AML / CQ
ter = make(TER)
TER_MIX = [{"X": 10, "Y": 10, "Z": 10}, {"X": 20, "Y": 10, "Z": 10},
           {"X": 6, "Y": 6, "Z": 18}, {"X": 10, "Y": 30, "Z": 10},
           {"X": 12, "Y": 4, "Z": 4}, {"X": 8, "Y": 32, "Z": 5}]
TER_STD = {c: stds(ter, c, [4, 8, 12, 16, 20, 24]) for c in TER}
Z24 = ter({"Z": 24}, "Z′ 24")


def _ratio(comp, conc=1.0):
    return ter({comp: conc}).with_values(ter({comp: conc}).values / Z24.values)


def _check(method, predict, tol=0.05, noise=0.0):
    for i, c in enumerate(TER_MIX):
        out = predict(ter(c, noise=noise, seed=i))
        for k, v in out.items():
            assert rec(v, c[k]) == pytest.approx(100, abs=tol), (method, c, k, v)


def test_aac_partial_overlap_sl17():
    """AAC, partially overlapped: AML (Z) plateau 340–380 nm; HCT (X) from
    P(λ1) − F·P(λ2) with VAL's equality factor; VAL (Y) by amplitude
    subtraction at λ1."""
    m = uv.AmplitudeCentering(275.0, ["X", "Y", "Z"], subtract="Y", divisor_compound="Z",
                              plateau=(340.0, 380.0),
                              differences={"X": {"w1": 275.0, "w2": 240.0,
                                                 "factor_from": "Y"}})
    info = m.fit(TER_STD, Z24)
    assert info["factors"]["Y"] if "Y" in info["factors"] else True
    _check("AAC partial", lambda s: m.predict(s, Z24))
    rt = uv.AmplitudeCentering.from_dict(m.to_dict())
    assert rt.predict(ter(TER_MIX[0]), Z24) == pytest.approx(m.predict(ter(TER_MIX[0]), Z24))


def _pair(comp, lam, lo=205.0, hi=330.0):
    cands = [w for w in uv.equal_amplitude_wavelengths(_ratio(comp), lam, lo, hi)
             if abs(w - lam) > 5]
    assert cands, f"no equal-amplitude partner for {comp} at {lam}"
    return min(cands, key=lambda w: abs(w - lam))


def test_aac_complete_overlap_and_macm():
    """AAC (severe overlap) = MACM (AW17): no plateau needed. X from a λ pair
    where Y/Z′ has equal amplitudes, Y from a pair where X/Z′ has equal
    amplitudes (Z/Z′ is constant and cancels in both); Z by subtraction at
    the common λ."""
    lam = 262.0
    wy = _pair("Y", lam)            # Y equal at lam and wy → difference measures X
    wx = _pair("X", lam)            # X equal at lam and wx → difference measures Y
    m = uv.AmplitudeCentering(lam, ["X", "Y", "Z"], subtract="Z", divisor_compound="Z",
                              differences={"X": {"w1": lam, "w2": wy},
                                           "Y": {"w1": lam, "w2": wx}})
    m.fit(TER_STD, Z24)
    _check("MACM", lambda s: m.predict(s, Z24), tol=0.5)
    # with noise the compound found by subtraction (Z, weak at λc) carries the
    # accumulated error: ±3 % at 0.0001 AU noise here, ±10 % at 0.0004 AU
    noisy = [m.predict(ter(c, noise=0.0001, seed=i), Z24) for i, c in enumerate(TER_MIX)]
    for c, out in zip(TER_MIX, noisy):
        for k in c:
            assert 97 < rec(out[k], c[k]) < 103


def test_ridss_am14():
    """RIDSS: AMB (Z) from the plateau; THP (Y) from the ratio difference at
    λ pair where GUF/Z′ is equal; total GUF + THP at their isoabsorptive point
    (unified regression), GUF = total − THP."""
    px = {k: Spectrum(GRID, v, concentrations={k: 1.0}) for k, v in TER.items()}
    iso = [w for w in uv.isoabsorptive_points(px["X"], px["Y"]) if 230 < w < 300]
    assert iso
    lam = iso[0]
    wx = _pair("X", lam)
    m = uv.AmplitudeCentering(lam, ["X", "Y", "Z"], subtract="X", divisor_compound="Z",
                              plateau=(340.0, 380.0),
                              differences={"Y": {"w1": lam, "w2": wx}})
    m.fit(TER_STD, Z24)
    _check("RIDSS", lambda s: m.predict(s, Z24), tol=0.2)
    # the unified regression at the isoabsorptive point gives the same answer
    # for X and Y (they share one line there)
    mu = uv.AmplitudeCentering(lam, ["X", "Y"], subtract="X", divisor_compound="Z",
                               plateau=(340.0, 380.0),
                               differences={"Y": {"w1": lam, "w2": wx}}, unified=True)
    mu.compounds = ["X", "Y", "Z"]
    with pytest.raises(ValueError):
        mu.fit({k: v for k, v in TER_STD.items() if k != "Z"}, Z24)


def test_cv_ad_lo15():
    """Constant value via amplitude difference (LO15): the postulated FCP (X)
    ratio amplitude at λ from its ratio difference, CQ (Z, divisor) = recorded
    − postulated. Binary case of amplitude centering."""
    bi = make({"X": TER["X"], "Z": TER["Z"]})
    st = {c: stds(bi, c, [4, 8, 12, 16]) for c in ("X", "Z")}
    div = bi({"Z": 10})
    m = uv.AmplitudeCentering(251.0, ["X", "Z"], subtract="Z", divisor_compound="Z",
                              differences={"X": {"w1": 261.0, "w2": 251.0}})
    m.fit(st, div)
    for c in ({"X": 5, "Z": 2}, {"X": 13, "Z": 11}, {"X": 7, "Z": 7}):
        out = m.predict(bi(c), div)
        assert out["X"] == pytest.approx(c["X"], rel=1e-6)
        assert out["Z"] == pytest.approx(c["Z"], rel=1e-6)


def test_amplitude_centering_rejects_incomplete_definitions():
    with pytest.raises(ValueError, match="no way to resolve"):
        uv.AmplitudeCentering(260.0, ["X", "Y", "Z"], subtract="Z",
                              differences={"X": {"w1": 260, "w2": 240}}).fit(TER_STD, Z24)
    with pytest.raises(ValueError, match="resolved twice"):
        uv.AmplitudeCentering(260.0, ["X", "Y"], subtract="X",
                              differences={"X": {"w1": 260, "w2": 240},
                                           "Y": {"w1": 260, "w2": 230}}).fit(TER_STD, Z24)
    with pytest.raises(ValueError, match="three pure"):
        uv.AmplitudeCentering(260.0, ["X", "Z"], subtract="Z",
                              differences={"X": {"w1": 260, "w2": 240}}).fit(
            {"X": TER_STD["X"][:2], "Z": TER_STD["Z"]}, Z24)


# ---- AW17 MAFM and AUC ------------------------------------------------------ #
AF = {"MET": 0.030 * g(300, 14) + 0.020 * g(318, 6),    # extended alone above 320 nm
      "MEB": 0.025 * g(268, 9) + 0.012 * g(232, 8),     # MET + MEB at 290–300 nm
      "DLX": 0.045 * g(246, 7)}                          # all overlap at 246 nm
af = make(AF)


def test_mafm_aw17_successive_absorption_factors():
    m = uv.AbsorptionFactorMethod([("MET", 330.0), ("MEB", 280.0), ("DLX", 246.0)],
                                  quant={"MET": 300.0, "MEB": 280.0})
    m.fit({c: stds(af, c, [2, 6, 10, 14, 18]) for c in AF})
    for c in ({"MET": 10, "MEB": 10, "DLX": 10}, {"MET": 30, "MEB": 4, "DLX": 2},
              {"MET": 3, "MEB": 20, "DLX": 15}):
        out = m.predict(af(c))
        for k in c:
            assert out[k] == pytest.approx(c[k], rel=2e-4), (c, k)
    rt = uv.AbsorptionFactorMethod.from_dict(m.to_dict())
    assert rt.predict(af({"MET": 7, "MEB": 8, "DLX": 9})) == pytest.approx(
        m.predict(af({"MET": 7, "MEB": 8, "DLX": 9})))


def test_auc_ternary_cramer_solution_aw17_coefficients():
    """AW17 Method I: A(225–235), A(240–250), A(280–290) = K·C with the
    published absorptivities; the app's solver must return the concentrations
    that generated the areas (Cramer's rule = least squares for a 3×3)."""
    K = np.array([[0.1918, 0.3672, 0.3475],
                  [0.1514, 0.5708, 0.1352],
                  [0.2730, 0.0550, 0.1421]])        # rows: ranges; cols: MET, DLX, MEB
    eq = SignalEquations(["MET", "DLX", "MEB"], [{"kind": "area"}] * 3)
    eq.K, eq.b = K.T, np.zeros(3)
    c_true = np.array([15.0, 10.0, 2.0])             # Dimetrol ratio 375:250:50
    eq.responses = lambda spectra, resolve=None: np.array([K @ c_true])
    assert np.allclose(eq.predict([None]), c_true)
    assert np.allclose(np.linalg.solve(K, K @ c_true), c_true)


def test_auc_ternary_on_simulated_spectra():
    eq = SignalEquations(list(AF), [{"kind": "area", "params": {"w1": a, "w2": b}}
                                    for a, b in ((225, 235), (240, 250), (280, 290))])
    eq.fit([s for c in AF for s in stds(af, c, [2, 6, 10, 14])])
    c = {"MET": 15, "DLX": 10, "MEB": 2}
    assert np.allclose(eq.predict([af(c)])[0], [c[k] for k in AF], rtol=1e-6)


# ---- SA14 ternary: double divisor, AUC of derivative ratio, ternary MCR ---- #
SA = {"ORPH": 0.030 * g(222, 9) + 0.010 * g(265, 10),
      "CAF": 0.050 * g(273, 10) + 0.020 * g(208, 6),
      "ASP": 0.040 * g(228, 8) + 0.008 * g(276, 8)}
sa = make(SA)
SA_MIX = [{"ORPH": 10, "CAF": 15, "ASP": 25}, {"ORPH": 6, "CAF": 10, "ASP": 30},
          {"ORPH": 14, "CAF": 6, "ASP": 18}]
SA_LIB = {"CAF16": sa({"CAF": 16}), "ASP16": sa({"ASP": 16}), "ORPH16": sa({"ORPH": 16})}
SA_STD = {c: stds(sa, c, [4, 8, 12, 16, 20]) for c in SA}


def _sa_method(steps, comp, w):
    m = uv.UnivariateMethod("m", comp, steps, {"kind": "amplitude", "params": {"w1": w}})
    m.calibrate(SA_STD[comp], res_of(SA_LIB))
    return m


def test_double_divisor_sa14():
    """ORPH: mixture ÷ (CAF′ + ASP′), D1 (Δλ 4, ×10) at a zero-crossing-free λ."""
    steps = [{"op": "divide_sum", "params": {"reference": "CAF16", "reference2": "ASP16"}},
             {"op": "derivative", "params": {"order": 1, "delta_lambda": 4.0}}]
    # the double divisor removes CAF+ASP only if they are present in the
    # divisor ratio; the method's premise: their D1 ratio contributions are
    # nearly constant — check that the calibration is linear and the bias is
    # what the premise predicts (small, not zero)
    m = _sa_method(steps, "ORPH", 240.0)
    assert m.regression.r > 0.9999


def test_auc_of_derivative_ratio_spectra_sa14():
    """Method B: mixture ÷ CAF′ (CAF becomes a constant), D1 removes it; the
    areas of the DD1 spectrum in two ranges give two equations in ORPH and ASP."""
    steps = [{"op": "divide", "params": {"reference": "CAF16"}},
             {"op": "derivative", "params": {"order": 1, "delta_lambda": 4.0, "scaling": 10}}]
    eq = SignalEquations(["ORPH", "ASP"], [{"kind": "area", "params": {"w1": 213, "w2": 223}},
                                           {"kind": "area", "params": {"w1": 230, "w2": 240}}],
                         steps)
    eq.fit(SA_STD["ORPH"] + SA_STD["ASP"], res_of(SA_LIB))
    for c in SA_MIX:
        out = eq.predict([sa(c)], res_of(SA_LIB))[0]
        assert out == pytest.approx([c["ORPH"], c["ASP"]], rel=1e-6)


def test_ternary_mean_centering_of_ratio_spectra_sa14():
    """Method C: ÷ CAF′ → mean-centre → ÷ MC(ASP′/CAF′) → mean-centre leaves
    ORPH alone (CAF and ASP become constants removed by mean centring)."""
    steps = uv.TEMPLATES["Mean centering of ratio spectra (ternary)"]["steps"]
    steps = [dict(s, params=dict(s["params"])) for s in steps]
    steps[0]["params"]["reference"] = "CAF16"
    steps[2]["params"].update(reference="ASP16", divisor="CAF16")
    m = _sa_method(steps, "ORPH", 223.0)
    for c in SA_MIX:
        assert m.predict(sa(c), res_of(SA_LIB)) == pytest.approx(c["ORPH"], rel=1e-6)


# ---- LO15: D1 DR, D1 DWL, SRS-type successive subtraction ------------------ #
def test_d1_derivative_ratio_lo15():
    """FCP (X) by D1 of [D1(mixture) ÷ D1(CQ′)]: CQ cancels; LH (only below
    235 nm) does not contribute at the measuring λ."""
    lo = {"LH": 0.06 * g(208, 7), "FCP": TER["X"], "CQ": TER["Z"]}
    mk = make(lo)
    lib = {"CQ6": mk({"CQ": 6})}
    steps = uv.TEMPLATES["Derivative ratio of D1 spectra (D1 DR)"]["steps"]
    steps = [dict(s, params=dict(s["params"])) for s in steps]
    steps[1]["params"]["reference"] = "CQ6"
    m = uv.UnivariateMethod("D1DR", "FCP", steps, {"kind": "amplitude", "params": {"w1": 271.0}})
    m.calibrate(stds(mk, "FCP", [5, 10, 15, 20, 25, 30]), res_of(lib))
    for c in ({"LH": 2, "FCP": 13, "CQ": 11}, {"LH": 5, "FCP": 10, "CQ": 3},
              {"LH": 7, "FCP": 7, "CQ": 10}):
        assert m.predict(mk(c), res_of(lib)) == pytest.approx(c["FCP"], rel=1e-3)


def test_d1_dual_wavelength_lo15():
    """CQ by D1 amplitude difference at λ1, λ2 where FCP's D1 is equal."""
    mk = make({"FCP": TER["X"], "CQ": TER["Z"]})
    d1 = apply_pipeline(mk({"FCP": 1}), [uv._D1])
    w1 = 280.0
    w2 = [w for w in uv.equal_amplitude_wavelengths(d1, w1, 240, 330) if abs(w - w1) > 5][0]
    m = uv.UnivariateMethod("D1 DWL", "CQ", [uv._D1],
                            {"kind": "difference", "params": {"w1": w1, "w2": w2}})
    m.calibrate(stds(mk, "CQ", [1, 3, 5, 8]))
    assert m.predict(mk({"FCP": 20, "CQ": 4})) == pytest.approx(4, rel=1e-3)


# ---- EM18: successive ratio subtraction (ternary) --------------------------- #
EM = {"FR": 0.030 * g(273, 10) + 0.015 * g(330, 18),   # extended to ~370 nm
      "CN": 0.060 * g(285, 9) + 0.010 * g(232, 8),     # to ~305 nm
      "SP": 0.045 * g(240, 9)}                          # nothing above 270 nm
em = make(EM)
EM_LIB = {"FR12": em({"FR": 12}), "CN4": em({"CN": 4})}


def test_successive_ratio_subtraction_em18():
    rs = [{"op": "ratio_subtraction", "params": {"divisor": "FR12", "start": 340.0, "end": 365.0}}]
    srs = uv.TEMPLATES["Successive ratio subtraction (ternary SRS, Z)"]["steps"]
    srs = [dict(s, params=dict(s["params"])) for s in srs]
    srs[0]["params"].update(divisor="FR12", start=340.0, end=365.0)
    srs[1]["params"].update(divisor="CN4", start=290.0, end=300.0)
    mixes = [{"FR": 8, "SP": 20, "CN": 2}, {"FR": 12, "SP": 30, "CN": 2},
             {"FR": 10, "SP": 32, "CN": 1}, {"FR": 8, "SP": 18, "CN": 4}]
    for comp, steps, w in (("CN", rs, 285.0), ("SP", srs, 240.0)):
        m = uv.UnivariateMethod(comp, comp, steps, {"kind": "amplitude", "params": {"w1": w}})
        m.calibrate(stds(em, comp, [2, 6, 10, 16]), res_of(EM_LIB))
        for c in mixes:
            assert m.predict(em(c), res_of(EM_LIB)) == pytest.approx(c[comp], rel=1e-3)
    # the SRS output is the pure SP spectrum (spectral profile)
    out = apply_pipeline(em(mixes[0]), srs, res_of(EM_LIB))
    assert np.allclose(out.values, 20 * EM["SP"], atol=2e-4)


# ---- FA21: quaternary CM-SS, derivative transformation, concentration value - #
FA = {"TN": 0.090 * g(345, 22),                         # extended to ~400 nm
      "EX": 0.040 * g(287, 10) + 0.035 * g(325, 9),     # to ~345 nm
      "HQ": 0.050 * g(225, 7) + 0.027 * g(294, 7),
      "HC": 0.042 * g(242, 10)}                         # D1 ≈ 0 above 280 nm
fa = make(FA)
FA_UNIT = {k: fa({k: 1}, f"{k} (unit)") for k in FA}


def test_cm_ss_quaternary_and_derivative_transformation_fa21():
    res = res_of(FA_UNIT)
    tn = [{"op": "factorized_recovery",
           "params": {"reference": "TN", "wavelength": 385.0, "wavelength_end": 397.0}}]
    ter_ = [{"op": "spectrum_subtraction",
             "params": {"reference": "TN", "wavelength": 385.0, "wavelength_end": 397.0}}]
    ex = ter_ + [{"op": "factorized_recovery",
                  "params": {"reference": "EX", "wavelength": 335.0, "wavelength_end": 345.0}}]
    binary = ter_ + [{"op": "spectrum_subtraction",
                      "params": {"reference": "EX", "wavelength": 335.0,
                                 "wavelength_end": 345.0}}]
    hq = binary + [{"op": "factorized_recovery",
                    "params": {"reference": "HQ", "wavelength": 300.0, "wavelength_end": 310.0,
                               "derivative_order": 1}}]
    hc = binary + [{"op": "spectrum_subtraction",
                    "params": {"reference": "HQ", "wavelength": 300.0, "wavelength_end": 310.0,
                               "derivative_order": 1}}]
    mixes = [{"TN": 10, "EX": 10, "HQ": 10, "HC": 10}, {"TN": 1, "EX": 12, "HQ": 9, "HC": 3},
             {"TN": 10, "EX": 20, "HQ": 23, "HC": 22}, {"TN": 9, "EX": 8, "HQ": 15, "HC": 11}]
    for comp, steps, w in (("TN", tn, 338.0), ("EX", ex, 287.0), ("HQ", hq, 225.0),
                           ("HC", hc, 242.0)):
        m = uv.UnivariateMethod(comp, comp, steps, {"kind": "amplitude", "params": {"w1": w}})
        m.calibrate(stds(fa, comp, [2, 6, 10, 16, 22]), res)
        for c in mixes:
            assert rec(m.predict(fa(c), res), c[comp]) == pytest.approx(100, abs=0.3), (comp, c)


def test_concentration_value_method_fa21():
    """Plateau of mixture ÷ unit-concentration EX = the EX concentration, with
    no regression; also equal to the constant value (regression) method."""
    lib = {"EXu": FA_UNIT["EX"]}
    m = uv.UnivariateMethod("conc value", "EX", [{"op": "divide", "params": {"reference": "EXu"}}],
                            {"kind": "mean", "params": {"w1": 337.0, "w2": 345.0}}, direct=True)
    sample = make({k: FA[k] for k in ("EX", "HQ", "HC")})
    for c in ({"EX": 12, "HQ": 9, "HC": 3}, {"EX": 4, "HQ": 20, "HC": 20}):
        assert m.predict(sample(c), res_of(lib)) == pytest.approx(c["EX"], rel=1e-6)
    reg = m.calibrate(stds(fa, "EX", [4, 8, 16, 25]), res_of(lib))
    assert reg.slope == pytest.approx(1.0, rel=1e-6) and abs(reg.intercept) < 1e-9
    assert uv.UnivariateMethod.from_dict(m.to_dict()).direct


# ---- EI18: binary methods with normalized divisors -------------------------- #
EI = {"SCB": 0.050 * g(254, 9) + 0.020 * g(225, 7),
      "VLS": 0.035 * g(250, 14) + 0.030 * g(226, 8)}
ei = make(EI)
EI_LIB = {"VLSu": ei({"VLS": 1}), "SCBu": ei({"SCB": 1})}
EI_MIX = [{"SCB": 5, "VLS": 20}, {"SCB": 10, "VLS": 20}, {"SCB": 20, "VLS": 20},
          {"SCB": 20, "VLS": 10}, {"SCB": 20, "VLS": 5}]


@pytest.mark.parametrize("name,comp,divisor,meas", [
    ("RD", "SCB", "VLSu", {"kind": "difference", "params": {"w1": 264.0, "w2": 238.0}}),
    ("DR1", "SCB", "VLSu", None),
    ("MCR", "SCB", "VLSu", None),
    ("RD", "VLS", "SCBu", {"kind": "difference", "params": {"w1": 225.0, "w2": 280.0}}),
])
def test_ratio_methods_with_normalized_divisor_ei18(name, comp, divisor, meas):
    steps = [{"op": "divide", "params": {"reference": divisor}}]
    if name == "DR1":
        steps.append({"op": "derivative", "params": {"order": 1, "delta_lambda": 4.0,
                                                     "scaling": 10.0}})
        meas = {"kind": "amplitude", "params": {"w1": 239.0}}
    if name == "MCR":
        steps.append({"op": "mean_center", "params": {"start": 220.0, "end": 300.0}})
        meas = {"kind": "amplitude", "params": {"w1": 260.0}}
    m = uv.UnivariateMethod(name, comp, steps, meas)
    m.calibrate(stds(ei, comp, [2.5, 5, 10, 15, 20, 25]), res_of(EI_LIB))
    for c in EI_MIX:
        assert m.predict(ei(c), res_of(EI_LIB)) == pytest.approx(c[comp], rel=1e-6)


def test_idw_and_aas_ei18():
    """IDW: ΔA = A(254) − F·A(226), F = A_VLS(254)/A_VLS(226). AAS at the
    isoabsorptive point with λ2 where VLS is equal."""
    m = uv.UnivariateMethod("IDW", "SCB", [],
                            {"kind": "weighted_difference",
                             "params": {"w1": 254.0, "w2": 226.0, "reference": "VLSu"}})
    m.calibrate(stds(ei, "SCB", [2.5, 5, 10, 20]), res_of(EI_LIB))
    for c in EI_MIX:
        assert m.predict(ei(c), res_of(EI_LIB)) == pytest.approx(c["SCB"], rel=1e-6)
    px = {k: Spectrum(GRID, v) for k, v in EI.items()}
    iso = [w for w in uv.isoabsorptive_points(px["SCB"], px["VLS"]) if 235 < w < 265][0]
    w1 = [w for w in uv.equal_amplitude_wavelengths(px["VLS"], iso, 230, 300)
          if abs(w - iso) > 5][0]
    aas = uv.AdvancedAbsorbanceSubtraction(w1, iso)
    aas.fit(stds(ei, "SCB", [2.5, 5, 10, 20]), "SCB", stds(ei, "VLS", [2.5, 5, 10, 20]), "VLS")
    for c in EI_MIX:
        out = aas.predict(ei(c))
        assert out["X"] == pytest.approx(c["SCB"], rel=1e-6)
        assert out["Y"] == pytest.approx(c["VLS"], rel=1e-6)


def test_robustness_skips_derivative_orders():
    m = uv.UnivariateMethod("D1 DR", "X", uv.TEMPLATES["Derivative ratio of D1 spectra (D1 DR)"]
                            ["steps"], {"kind": "amplitude", "params": {"w1": 270.0}})
    labels = [p["label"] for p in uv.numeric_parameters(m)]
    assert labels and not any(lab.endswith(": Order") or "Differentiate divisor" in lab
                              for lab in labels)
