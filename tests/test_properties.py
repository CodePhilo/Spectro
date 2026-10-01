"""Property-based and metamorphic tests: rules that must hold for *any* input."""

import numpy as np
import pytest
from hypothesis import HealthCheck, assume, given, settings
from hypothesis import strategies as st

from spectro.core import univariate as uv
from spectro.core.io import export_csv, export_excel, load_file, write_jcamp, write_spc
from spectro.core.multicomponent import MODEL_TYPES, SpectralModel, design_concentrations
from spectro.core.operations import REGISTRY, apply_pipeline
from spectro.core.spectrum import Spectrum
from tests.conftest import GRID, PURE, band, mixture

SETTINGS = settings(max_examples=40, deadline=None,
                    suppress_health_check=[HealthCheck.function_scoped_fixture,
                                           HealthCheck.too_slow])


def smooth_spectrum(seed: int, x=GRID) -> Spectrum:
    rng = np.random.default_rng(seed)
    y = sum(rng.uniform(0.01, 0.06) * np.exp(-0.5 * ((x - rng.uniform(210, 390))
                                                     / rng.uniform(6, 30)) ** 2)
            for _ in range(3))
    return Spectrum(x, y + 0.01)


# --------------------------------------------------------------------------- #
# 1. Linearity — the method optimizer relies on it
# --------------------------------------------------------------------------- #
DIV = Spectrum(GRID, PURE["Y"] + 0.002, name="div")
REF = Spectrum(GRID, PURE["X"] + 0.001, name="ref")
LIB = {"div": DIV, "ref": REF}

LINEAR_STEPS = {
    "crop": {"start": 220.0, "end": 380.0},
    "resample": {"step": 2.0},
    "smooth_sg": {"window": 11, "polyorder": 2},
    "smooth_ma": {"window": 5},
    "smooth_whittaker": {"lam": 20.0},
    "baseline_offset": {"start": 380.0, "end": 400.0},
    "baseline_linear": {"w1": 210.0, "w2": 390.0},
    "divide": {"reference": "div"},
    "divide_sum": {"reference": "div", "reference2": "ref"},
    "multiply": {"reference": "ref"},
    "scale": {"factor": 2.5},
    "mean_center": {"start": 220.0, "end": 300.0},
    "subtract_plateau": {"start": 340.0, "end": 360.0},
    "derivative": {"order": 2, "delta_lambda": 4.0},
    "ratio_subtraction": {"divisor": "div", "start": 340.0, "end": 360.0},
    "constant_multiplication": {"divisor": "div", "start": 340.0, "end": 360.0},
    "extended_ratio_subtraction": {"divisor": "div", "reference": "ref", "start": 340.0,
                                   "end": 360.0},
    "factorized_recovery": {"reference": "div", "wavelength": 350.0, "derivative_order": 0},
    "spectrum_subtraction": {"reference": "div", "wavelength": 350.0, "derivative_order": 1},
    "constant_center": {"divisor": "div", "reference": "ref", "w1": 240.0, "w2": 260.0},
    "divide_centered_ratio": {"reference": "ref", "divisor": "div", "start": 220.0,
                              "end": 300.0},
}


@pytest.mark.parametrize("op", sorted(LINEAR_STEPS))
@SETTINGS
@given(seed1=st.integers(0, 10 ** 6), seed2=st.integers(0, 10 ** 6),
       a=st.floats(-3, 3), b=st.floats(-3, 3))
def test_processing_operations_are_linear(op, seed1, seed2, a, b):
    s1, s2 = smooth_spectrum(seed1), smooth_spectrum(seed2)
    steps = [{"op": op, "params": LINEAR_STEPS[op]}]
    res = LIB.__getitem__
    lhs = apply_pipeline(s1.with_values(a * s1.values + b * s2.values), steps, res)
    p1, p2 = apply_pipeline(s1, steps, res), apply_pipeline(s2, steps, res)
    scale = max(1e-12, np.max(np.abs(p1.values)) + np.max(np.abs(p2.values)))
    assert np.allclose(lhs.values, a * p1.values + b * p2.values, atol=1e-9 * scale * 10)


def test_every_registered_operation_is_classified():
    """New operations must be added to LINEAR_STEPS or declared non-linear here,
    so the optimizer can never silently use a non-linear step."""
    nonlinear = {"offset", "subtract", "add", "baseline_poly", "baseline_als", "normalize",
                 "snv", "t_to_a", "a_to_t"}
    assert set(REGISTRY) == set(LINEAR_STEPS) | nonlinear


# --------------------------------------------------------------------------- #
# 2. Invariances of calibration methods
# --------------------------------------------------------------------------- #
def _design(noise=0.0005, seed=0, factor=1.0):
    conc = design_concentrations({"X": 10, "Y": 10, "Z": 10}, {"X": 4, "Y": 4, "Z": 4})
    out = []
    for i, c in enumerate(conc):
        m = mixture(c, f"D{i}", noise, seed=seed + i)
        m.concentrations = {k: v * factor for k, v in m.concentrations.items()}
        out.append(m)
    return out


@pytest.mark.parametrize("mtype", MODEL_TYPES)
def test_models_do_not_depend_on_concentration_units(mtype):
    """µg/mL → ng/mL (× 1000) must scale predictions by exactly 1000."""
    kw = dict(n_components=3, wavelengths=[225, 245, 260, 275, 330], ranges=[(210, 370)])
    test = _design(seed=500)[:5]
    a = SpectralModel(mtype, ["X", "Y", "Z"], **kw).fit(_design()).predict(test)
    b = SpectralModel(mtype, ["X", "Y", "Z"], **kw).fit(_design(factor=1000)).predict(test)
    # a neural network is reproducible only to its training precision (~0.5 %)
    tol = 1e-2 if mtype == "ANN" else 1e-3
    assert np.allclose(b, 1000 * a, rtol=tol), mtype


@pytest.mark.parametrize("mtype", MODEL_TYPES)
def test_models_do_not_depend_on_sample_order(mtype):
    kw = dict(n_components=3, wavelengths=[225, 245, 260, 275, 330], ranges=[(210, 370)])
    cal = _design()
    test = _design(seed=600)[:5]
    a = SpectralModel(mtype, ["X", "Y", "Z"], **kw).fit(cal).predict(test)
    perm = np.random.default_rng(1).permutation(len(cal))
    b = SpectralModel(mtype, ["X", "Y", "Z"], **kw).fit([cal[i] for i in perm]).predict(test)
    tol = 1e-2 if mtype == "ANN" else 1e-3
    assert np.allclose(a, b, rtol=tol, atol=tol), mtype


def test_univariate_units_and_grid_invariance():
    xs = [mixture({"X": c, "Y": 0}, f"X{c}", 0.0003, seed=c) for c in (4, 8, 12, 16, 20)]
    mix = [mixture({"X": 7, "Y": 11}, "m", 0.0003, seed=99)]
    lib = {"Y": Spectrum(GRID, PURE["Y"] * 10)}
    steps = [{"op": "divide", "params": {"reference": "Y"}}]
    meas = {"kind": "difference", "params": {"w1": 240.0, "w2": 260.0}}
    m = uv.UnivariateMethod("rd", "X", steps, meas)
    m.calibrate(xs, lib.__getitem__)
    base = m.predict(mix[0], lib.__getitem__)
    # ×1000 concentration units
    xs_k = [s.copy(concentrations={"X": s.concentrations["X"] * 1000}) for s in xs]
    m.calibrate(xs_k, lib.__getitem__)
    assert m.predict(mix[0], lib.__getitem__) == pytest.approx(1000 * base, rel=1e-9)
    # sample measured on a coarser (1 nm) grid and in descending order
    s = mix[0]
    coarse = Spectrum(s.x[::2][::-1], s.y[::2][::-1])
    m.calibrate(xs, lib.__getitem__)
    assert m.predict(coarse, lib.__getitem__) == pytest.approx(base, rel=2e-3)


# --------------------------------------------------------------------------- #
# 3. Fuzzing: random pipelines never fail with an unexpected exception
# --------------------------------------------------------------------------- #
def _random_params(op, draw):
    params = {}
    for p in REGISTRY[op].params:
        if p.kind == "spectrum":
            params[p.name] = draw(st.sampled_from(["div", "ref"]))
        elif p.kind == "choice":
            params[p.name] = draw(st.sampled_from(list(p.choices)))
        elif p.kind == "bool":
            params[p.name] = draw(st.booleans())
        elif p.kind == "int":
            lo = int(p.minimum) if p.minimum is not None else -5
            hi = int(p.maximum) if p.maximum is not None else 60
            params[p.name] = draw(st.integers(lo, max(lo, hi)))
        elif p.kind == "wavelength":
            params[p.name] = draw(st.one_of(st.none(), st.floats(150, 450))) if p.optional \
                else draw(st.floats(150, 450))
        else:
            lo = p.minimum if p.minimum is not None else -1e3
            hi = p.maximum if p.maximum is not None else 1e3
            params[p.name] = draw(st.floats(lo, hi, allow_nan=False))
    return params


@settings(max_examples=300, deadline=None)
@given(data=st.data(), seed=st.integers(0, 10 ** 6))
def test_random_pipelines_fail_only_with_clear_errors(data, seed):
    ops = sorted(REGISTRY)
    steps = [{"op": op, "params": _random_params(op, data.draw)}
             for op in data.draw(st.lists(st.sampled_from(ops), min_size=1, max_size=3))]
    try:
        out = apply_pipeline(smooth_spectrum(seed), steps, LIB.__getitem__)
    except ValueError:
        return          # clear, user-facing message
    assert out.wavelengths.size >= 2
    assert np.all(np.diff(out.wavelengths) > 0)
    # non-finite output must never leak silently
    assert np.all(np.isfinite(out.values)), steps


@settings(max_examples=150, deadline=None)
@given(data=st.data(), seed=st.integers(0, 10 ** 6))
def test_measurements_fail_only_with_clear_errors(data, seed):
    kind = data.draw(st.sampled_from(sorted(uv.MEASUREMENTS)))
    params = {"w1": data.draw(st.floats(150, 450)), "w2": data.draw(st.floats(150, 450))}
    if kind == "weighted_difference":
        params["factor"] = data.draw(st.floats(-5, 5))
    try:
        v = uv.measure(smooth_spectrum(seed), {"kind": kind, "params": params})
    except ValueError:
        return
    assert np.isfinite(v)


# --------------------------------------------------------------------------- #
# 4. File round trips for random spectra and formats
# --------------------------------------------------------------------------- #
@settings(max_examples=40, deadline=None, suppress_health_check=[HealthCheck.function_scoped_fixture])
@given(n=st.integers(1, 6), step=st.sampled_from([0.1, 0.5, 1.0, 2.0]),
       start=st.floats(190, 260), npts=st.integers(20, 300),
       fmt=st.sampled_from(["csv", "xlsx", "jdx", "spc", "semicolon", "tab_rows", "xy_pairs"]),
       descending=st.booleans(), seed=st.integers(0, 1000))
def test_file_round_trips(tmp_path, n, step, start, npts, fmt, descending, seed):
    rng = np.random.default_rng(seed)
    x = np.round(start + step * np.arange(npts), 4)
    spectra = [Spectrum(x, np.round(rng.uniform(-0.05, 2.5, npts), 4), name=f"S {i + 1}")
               for i in range(n)]
    f = tmp_path / f"f_{fmt}_{seed}_{n}"
    if fmt == "csv":
        f = f.with_suffix(".csv")
        export_csv(spectra, f)
    elif fmt == "xlsx":
        f = f.with_suffix(".xlsx")
        export_excel(spectra, f)
    elif fmt == "jdx":
        assume(n == 1)
        f = f.with_suffix(".jdx")
        write_jcamp(spectra[0], f)
    elif fmt == "spc":
        f = f.with_suffix(".spc")
        write_spc(spectra, f)
    else:
        xs = x[::-1] if descending else x
        order = slice(None, None, -1) if descending else slice(None)
        if fmt == "semicolon":
            lines = ["Wavelength;" + ";".join(s.name for s in spectra)]
            lines += [f"{w:g}".replace(".", ",") + ";" +
                      ";".join(f"{s.values[order][i]:.4f}".replace(".", ",") for s in spectra)
                      for i, w in enumerate(xs)]
        elif fmt == "tab_rows":
            assume(npts >= 5)
            lines = ["Sample\t" + "\t".join(f"{w:g}" for w in xs)]
            lines += [f"{s.name}\t" + "\t".join(f"{v:.4f}" for v in s.values[order])
                      for s in spectra]
        else:
            assume(n >= 2)
            lines = [",".join(f"{s.name}," for s in spectra),
                     ",".join("nm,Abs" for _ in spectra)]
            lines += [",".join(f"{w:g},{s.values[order][i]:.4f}" for s in spectra)
                      for i, w in enumerate(xs)]
        f = f.with_suffix(".csv" if fmt != "tab_rows" else ".tsv")
        f.write_text("\n".join(lines), encoding="utf-8")
    r = load_file(f)
    assert len(r.spectra) == n, (fmt, r.layout)
    tol = 1e-4 if fmt == "spc" else 1e-9  # SPC stores float32
    for a, b in zip(spectra, r.spectra):
        assert np.allclose(b.wavelengths, a.wavelengths, atol=tol * 10 + 1e-6), fmt
        assert np.allclose(b.values, a.values, atol=tol * 10 + 1e-6), fmt
    if fmt not in ("spc", "jdx"):
        assert [s.name for s in r.spectra] == [s.name for s in spectra]


# --------------------------------------------------------------------------- #
# 5. Statistical calibration of the outlier limits
# --------------------------------------------------------------------------- #
def _noisy(conc, rng, noise=0.001, offset=0.002):
    """Realistic measurement: white noise + a whole-spectrum baseline offset."""
    m = mixture(conc, "n", noise, seed=int(rng.integers(1 << 30)))
    return m.copy(values=m.values + rng.normal(0, offset))


@pytest.mark.parametrize("model", ["PLS2", "PCR", "CLS"])
def test_limits_flag_about_five_and_one_percent_of_new_normal_samples(model):
    rng = np.random.default_rng(3)
    counts = {"t2_lim95": 0, "t2_lim99": 0, "q_lim95": 0, "q_lim99": 0}
    total = 0
    for rep in range(8):
        conc = design_concentrations({"X": 10, "Y": 10, "Z": 10}, {"X": 4, "Y": 4, "Z": 4})
        cal = [_noisy(c, rng).copy(concentrations=c) for c in conc]
        m = SpectralModel(model, ["X", "Y", "Z"], n_components=3).fit(cal)
        new = [_noisy({c: rng.uniform(6, 14) for c in "XYZ"}, rng) for _ in range(60)]
        d = m.diagnostics(new)
        for key in counts:
            vals = d["T2"] if key.startswith("t2") else d["Q"]
            counts[key] += int(np.sum(np.array(vals) > d["limits"][key]))
        total += len(new)
    rate = {k: v / total for k, v in counts.items()}
    assert 0.01 <= rate["q_lim95"] <= 0.10, rate
    assert rate["q_lim99"] <= 0.03, rate
    assert rate["t2_lim95"] <= 0.10 and rate["t2_lim99"] <= 0.03, rate


def test_limits_are_conservative_not_alarming_with_pure_white_noise():
    rng = np.random.default_rng(4)
    cal = _design(noise=0.001, seed=77)
    m = SpectralModel("PLS2", ["X", "Y", "Z"], n_components=3).fit(cal)
    new = [mixture({c: rng.uniform(6, 14) for c in "XYZ"}, "n", 0.001, seed=700 + i)
           for i in range(100)]
    d = m.diagnostics(new)
    assert np.mean(np.array(d["Q"]) > d["limits"]["q_lim95"]) <= 0.10


# --------------------------------------------------------------------------- #
# 6. Optimizer finds known-best strategies
# --------------------------------------------------------------------------- #
def test_optimizer_prefers_zero_order_when_a_free_band_exists():
    from spectro.core.optimizer import OptimizerInput, optimize
    x_only = 0.03 * band(340, 8)          # X has a band where Y is transparent
    y_unit = 0.04 * band(250, 15)

    def spec(cx, cy, seed):
        rng = np.random.default_rng(seed)
        return Spectrum(GRID, cx * x_only + cy * y_unit + rng.normal(0, 3e-4, GRID.size),
                        concentrations={"X": cx, "Y": cy})
    std = {"X": [spec(c, 0, c) for c in (4, 8, 12, 16, 20)],
           "Y": [spec(0, c, 100 + c) for c in (4, 8, 12, 16, 20)]}
    res = optimize(OptimizerInput(["X", "Y"], std, {"X": std["X"][2], "Y": std["Y"][2]},
                                  [spec(10, 10, 7), spec(6, 16, 8)],
                                  families={"zero", "derivative", "dual", "cls"}))
    zero = next(c for c in res["ranked"]["X"] if c.family == "zero")
    assert 330 <= zero.measurement["params"]["w1"] <= 350
    assert zero.score < 0.5
