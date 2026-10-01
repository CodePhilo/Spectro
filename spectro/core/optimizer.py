"""Method optimizer: find the best processing to determine each compound of a
mixture simultaneously.

Principle
---------
Every processing operation used for mixture analysis (derivatives, division by
a divisor spectrum, mean centering, ratio subtraction, constant multiplication,
spectrum subtraction…) is *linear in the sample spectrum*, and so are the
measurements (amplitude, difference, weighted difference). Therefore the signal
of a mixture is the sum of the signals of its components::

    S(mixture) = Σ_j c_j · S(unit spectrum of j)

With the unit-concentration spectra estimated from the pure standards, the
error that a candidate method would make for the target compound can be
*predicted* for any mixture composition without preparing it:

* interference (bias)  = Σ_{j≠t} c_j·S_j / (c_t·S_t)
* noise                = σ(S) / (c_t·|S_t|)  (σ from the measured noise, pushed
  through the same processing)
* wavelength robustness: the same evaluated with the wavelengths shifted by the
  instrument's wavelength uncertainty (e.g. ±0.5 nm); knife-edge zero-crossings
  are penalised.

1. **Screening** scans every processing family over all wavelengths / λ pairs
   and keeps the best settings of each family (fast, vectorised).
2. **Verification** builds the real method for each short-listed candidate,
   calibrates it on the real standards and predicts (a) the real laboratory
   mixtures, if given, and (b) simulated mixtures built from the unit spectra
   plus realistic noise — end-to-end, including the calibration regression.
3. Candidates are ranked by the larger of the observed error and the
   predicted error including wavelength uncertainty (conservative).

Progressive methods (amplitude centering, successive absorption factor)
resolve all compounds at once, so they are screened as whole configurations:
generated from the unit spectra for every divisor compound / compound order,
run on noise-free standards, scenario mixtures, noise and ±λ-shifted spectra,
and the best of each structure verified end-to-end like the others.
"""

from __future__ import annotations

import itertools
import math
from dataclasses import dataclass, field
from typing import Any, Callable

import numpy as np

from spectro.core.multicomponent import SignalEquations, SpectralModel, kaiser_selection
from spectro.core.operations import apply_pipeline
from spectro.core.spectrum import Spectrum, align
from spectro.core.univariate import (AbsorptionFactorMethod, AmplitudeCentering, UnivariateMethod,
                                    equal_amplitude_wavelengths, measure, progressive_from_dict)
from spectro.core.validation import linear_regression

DIV = "div:"  # placeholder prefix for divisor references, e.g. "div:CAF"

FAMILIES = {
    "zero": "Zero order (direct)",
    "derivative": "Derivative zero-crossing",
    "dual": "Dual wavelength",
    "induced_dual": "Induced dual wavelength",
    "ratio_difference": "Ratio difference",
    "derivative_ratio": "Derivative ratio",
    "mean_centering": "Mean centering of ratio spectra",
    "ratio_subtraction": "Ratio subtraction / constant multiplication",
    "dual_amplitude": "Dual amplitude difference (ternary)",
    "double_divisor": "Double divisor derivative ratio (ternary)",
    "amplitude_centering": "Amplitude centering (AAC / MACM / RIDSS / CV-AD)",
    "absorption_factor": "Successive absorption factor (MAFM)",
    "vierordt": "Vierordt simultaneous equations",
    "bivariate": "Bivariate (Kaiser)",
    "cls": "Classical least squares (full spectrum)",
    "pls": "PLS2 (mixture training set)",
}


@dataclass
class Candidate:
    compound: str
    family: str
    label: str
    steps: list[dict] = field(default_factory=list)
    measurement: dict | None = None
    model: dict | None = None            # multivariate definition (equations / spectral)
    predicted_error: float = math.inf    # % (screening, nominal λ)
    robust_error: float = math.inf       # % (screening, worst case over ±λ uncertainty)
    interference: float = math.nan       # % (bias part)
    noise: float = math.nan              # % (noise part)
    sensitivity: float = math.nan        # signal per unit concentration
    # verification
    sim_rmsep: float = math.nan          # % relative RMSEP on simulated mixtures
    real_rmsep: float = math.nan         # % relative RMSEP on laboratory mixtures
    real_mean_recovery: float = math.nan
    real_rsd: float = math.nan
    r: float = math.nan
    slope: float = math.nan
    lod: float = math.nan
    error: str = ""

    @property
    def score(self) -> float:
        """Ranking key, % (lower is better): the largest of the error observed on
        laboratory mixtures, on simulated mixtures and the predicted error
        including wavelength uncertainty — deliberately conservative."""
        vals = [v for v in (self.real_rmsep, self.sim_rmsep, self.robust_error)
                if not math.isnan(v)]
        return max(vals) if vals else math.inf

    def to_dict(self) -> dict:
        d = {k: getattr(self, k) for k in self.__dataclass_fields__}
        d["score"] = self.score
        return d


@dataclass
class OptimizerInput:
    compounds: list[str]
    standards: dict[str, list[Spectrum]]      # pure standards per compound
    divisors: dict[str, Spectrum] = field(default_factory=dict)
    mixtures: list[Spectrum] = field(default_factory=list)       # laboratory mixtures
    training: list[Spectrum] = field(default_factory=list)       # mixture training set (PLS)
    scenarios: list[dict[str, float]] | None = None
    wl_range: tuple[float, float] | None = None
    wl_uncertainty: float = 0.5
    noise_sigma: float | None = None
    families: set[str] | None = None
    pair_step: float = 2.0
    keep_per_family: int = 1
    seed: int = 0


# --------------------------------------------------------------------------- #
# Unit spectra, noise, scenarios
# --------------------------------------------------------------------------- #
def unit_spectra(inp: OptimizerInput) -> tuple[np.ndarray, dict[str, np.ndarray], np.ndarray]:
    """Grid, unit-concentration spectra (slope of A vs C per λ) and the
    per-wavelength noise σ estimated from the calibration residuals."""
    grid, units, sigma, _ = noise_model(inp)
    return grid, units, sigma


def noise_model(inp: OptimizerInput):
    """Like :func:`unit_spectra` plus the σ of whole-spectrum baseline offsets
    (from the mean residual of each standard). Without a user range the
    low-UV region where the noise exceeds 3× its median is excluded."""
    every = [s for c in inp.compounds for s in inp.standards[c]]
    every += list(inp.divisors.values())
    grid = align(every)[0].wavelengths
    if inp.wl_range:
        lo, hi = sorted(inp.wl_range)
        grid = grid[(grid >= lo - 1e-9) & (grid <= hi + 1e-9)]
    units, resid = {}, []
    for c in inp.compounds:
        std = inp.standards[c]
        if not std:
            raise ValueError(f"no pure standards for {c}")
        conc = np.array([s.concentrations[c] for s in std], float)
        X = np.vstack([np.interp(grid, s.wavelengths, s.values) for s in std])
        if len(std) >= 3 and np.ptp(conc) > 0:
            A = np.column_stack([conc, np.ones_like(conc)])
            coef, *_ = np.linalg.lstsq(A, X, rcond=None)
            units[c] = coef[0]
            resid.append(X - A @ coef)
        else:
            units[c] = (X / conc[:, None]).mean(0)
    offset = 0.0
    if resid:
        r = np.vstack(resid)
        dof = max(1, r.shape[0] - 2 * len(resid))
        offset = float(np.sqrt(np.sum(r.mean(axis=1) ** 2) / dof))
        r = r - r.mean(axis=1, keepdims=True)
        sd = np.sqrt(np.sum(r ** 2, axis=0) / dof)
        pad = np.pad(sd, 4, mode="edge")  # running median, 9 points
        sigma = np.array([np.median(pad[i:i + 9]) for i in range(sd.size)])
    else:
        x = np.vstack([np.interp(grid, s.wavelengths, s.values) for s in every])
        sigma = np.full(grid.size, float(np.median(np.abs(np.diff(x, axis=1))) / 1.35))
    if inp.noise_sigma is not None:
        sigma = np.full(grid.size, float(inp.noise_sigma))
    sigma = np.maximum(sigma, 1e-6)
    if not inp.wl_range and sigma.size > 20:
        ok = sigma <= 3 * np.median(sigma)
        best, start = None, None
        for i, v in enumerate(np.append(ok, False)):
            if v and start is None:
                start = i
            elif not v and start is not None:
                if best is None or i - start > best[1] - best[0]:
                    best = (start, i)
                start = None
        if best and best[1] - best[0] < sigma.size:
            sl = slice(*best)
            grid, sigma = grid[sl], sigma[sl]
            units = {c: u[sl] for c, u in units.items()}
    return grid, units, sigma, offset


def default_scenarios(inp: OptimizerInput) -> list[dict[str, float]]:
    """Mixture compositions used to judge interference: the laboratory
    mixtures if given, otherwise the corners and centre of the standards'
    concentration ranges."""
    if inp.scenarios:
        return inp.scenarios
    mixes = [m.concentrations for m in inp.mixtures
             if all(m.concentrations.get(c) is not None for c in inp.compounds)]
    if mixes:
        return [{c: float(m[c]) for c in inp.compounds} for m in mixes]
    ranges = {c: (min(s.concentrations[c] for s in inp.standards[c]),
                  max(s.concentrations[c] for s in inp.standards[c])) for c in inp.compounds}
    out = [dict(zip(inp.compounds, corner))
           for corner in itertools.product(*[ranges[c] for c in inp.compounds])]
    out.append({c: sum(ranges[c]) / 2 for c in inp.compounds})
    return out


def _band_max(grid, u) -> float:
    """λ of the strongest absorption band maximum (not a rising range edge)."""
    from scipy.signal import find_peaks

    idx, _ = find_peaks(u, prominence=0.05 * np.ptp(u))
    idx = [i for i in idx if grid[i] - grid[0] > 5 and grid[-1] - grid[i] > 5]
    return float(grid[max(idx, key=lambda i: u[i])]) if idx else float(grid[int(np.argmax(u))])


def _plateau(grid, target_u, other_u, min_width=8.0):
    """Region where only ``other`` absorbs (target < 0.5 %, other > 3 % of max)."""
    ok = (np.abs(target_u) < 0.005 * np.max(np.abs(target_u))) & \
         (other_u > 0.03 * np.max(other_u))
    best, start = None, None
    for i, v in enumerate(np.append(ok, False)):
        if v and start is None:
            start = i
        elif not v and start is not None:
            if grid[i - 1] - grid[start] >= min_width and (
                    best is None or grid[i - 1] - grid[start] > best[1] - best[0]):
                best = (float(grid[start]), float(grid[i - 1]))
            start = None
    return best


def _valid_region(grid, div_u, frac=0.05):
    """Largest contiguous region where the divisor is ≥ frac of its maximum."""
    ok = div_u >= frac * np.max(div_u)
    best, start = None, None
    for i, v in enumerate(np.append(ok, False)):
        if v and start is None:
            start = i
        elif not v and start is not None:
            if best is None or i - start > best[1] - best[0]:
                best = (start, i)
            start = None
    if best is None:
        return None
    return float(grid[best[0]]), float(grid[best[1] - 1])


def _valid_regions(grid, div_u, frac=0.05, min_width=20.0):
    """All contiguous regions (≥ min_width nm) where the divisor is ≥ frac of
    its maximum — a divisor with two bands has two usable regions."""
    ok = div_u >= frac * np.max(div_u)
    out, start = [], None
    for i, v in enumerate(np.append(ok, False)):
        if v and start is None:
            start = i
        elif not v and start is not None:
            if grid[i - 1] - grid[start] >= min_width:
                out.append((float(grid[start]), float(grid[i - 1])))
            start = None
    return out


# --------------------------------------------------------------------------- #
# Screening
# --------------------------------------------------------------------------- #
class _Screen:
    def __init__(self, inp: OptimizerInput):
        self.inp = inp
        self.grid, self.units, self.sigma, self.offset = noise_model(inp)
        self.scen = default_scenarios(inp)
        rng = np.random.default_rng(inp.seed)
        self.noise = [Spectrum(self.grid, rng.normal(0, 1, self.grid.size) * self.sigma
                               + rng.normal(0, self.offset)) for _ in range(40)]
        self.unit_sp = {c: Spectrum(self.grid, u, name=f"{c} (unit)")
                        for c, u in self.units.items()}
        self.resolve = self._resolver()

    def _resolver(self):
        divs = dict(self.inp.divisors)
        for c, u in self.unit_sp.items():
            divs.setdefault(c, u.with_values(u.values * 10))

        def resolve(ref):
            if isinstance(ref, str) and ref.startswith(DIV):
                return divs[ref[len(DIV):]]
            raise KeyError(ref)
        return resolve

    # processed unit spectra + noise, sampled at x (and x ± u)
    def processed(self, steps):
        P = {c: apply_pipeline(s, steps, self.resolve) for c, s in self.unit_sp.items()}
        N = [apply_pipeline(n, steps, self.resolve) for n in self.noise]
        return P, N

    def _errors(self, target, S, Snoise):
        """S: {compound: array of signals}, Snoise: (K, n) → error arrays (%)."""
        St = S[target]
        sd = Snoise.std(axis=0, ddof=1)
        bias2, noise2, tot = 0.0, 0.0, 0.0
        for m in self.scen:
            ct = m[target]
            if ct <= 0:
                continue
            with np.errstate(divide="ignore", invalid="ignore"):
                b = sum(m[j] * S[j] for j in S if j != target) / (ct * St)
                nz = sd / (ct * np.abs(St))
            bias2 = bias2 + b ** 2
            noise2 = noise2 + nz ** 2
            tot += 1
        bias = 100 * np.sqrt(bias2 / max(tot, 1))
        noise = 100 * np.sqrt(noise2 / max(tot, 1))
        err = np.sqrt(bias ** 2 + noise ** 2)
        err[~np.isfinite(err)] = np.inf
        return err, bias, noise

    def single(self, target, steps, measure_kind="amplitude"):
        """Best λ for amplitude of a processed spectrum."""
        P, N = self.processed(steps)
        x = P[target].wavelengths
        u = self.inp.wl_uncertainty
        lo, hi = x[0] + u, x[-1] - u
        cand = x[(x >= lo) & (x <= hi)]
        if cand.size == 0:
            return []

        def at(xx):
            S = {c: np.interp(xx, p.wavelengths, p.values) for c, p in P.items()}
            Sn = np.vstack([np.interp(xx, n.wavelengths, n.values) for n in N])
            return S, Sn
        S0, N0 = at(cand)
        err, bias, noise = self._errors(target, S0, N0)
        worst = err.copy()
        if u > 0:
            for d in (-u, u):
                Sd, _ = at(cand + d)
                # the method keeps the nominal λ; the instrument reads λ + d
                e_d, _, _ = self._errors(target, Sd, N0)
                worst = np.maximum(worst, e_d)
        tmax = np.max(np.abs(S0[target]))
        weak = np.abs(S0[target]) < 0.02 * tmax
        worst[weak] = np.inf
        return self._pick(cand[:, None], worst, err, bias, noise, S0[target])

    def pairs(self, target, steps, weight_ref: str | None = None):
        """Best (λ1, λ2) for P(λ1) − F·P(λ2); F from the processed reference
        (induced / dual amplitude) or 1 (dual wavelength / ratio difference)."""
        P, N = self.processed(steps)
        x = P[target].wavelengths
        u = self.inp.wl_uncertainty
        step = max(self.inp.pair_step, float(np.median(np.diff(x))))
        g = np.arange(x[0] + u, x[-1] - u + 1e-9, step)
        if g.size < 2:
            return []
        i1, i2 = np.triu_indices(g.size, k=1)
        l1, l2 = g[i1], g[i2]

        def vals(p, xx):
            return np.interp(xx, p.wavelengths, p.values)
        if weight_ref is not None:
            r1, r2 = vals(P[weight_ref], l1), vals(P[weight_ref], l2)
            with np.errstate(divide="ignore", invalid="ignore"):
                F = np.where(np.abs(r2) > 1e-12, r1 / r2, np.nan)
        else:
            F = np.ones_like(l1)

        def signals(d1, d2):
            S = {c: vals(p, l1 + d1) - F * vals(p, l2 + d2) for c, p in P.items()}
            Sn = np.vstack([vals(n, l1) - F * vals(n, l2) for n in N])
            return S, Sn
        S0, N0 = signals(0.0, 0.0)
        err, bias, noise = self._errors(target, S0, N0)
        worst = err.copy()
        if u > 0:
            for d1, d2 in ((-u, -u), (u, u)):  # instrument offset shifts both λ
                Sd, _ = signals(d1, d2)
                worst = np.maximum(worst, self._errors(target, Sd, N0)[0])
        tmax = np.max(np.abs(S0[target]))
        weak = (np.abs(S0[target]) < 0.02 * tmax) | ~np.isfinite(F) | (np.abs(F) > 20)
        worst[weak] = np.inf
        extra = {"F": F} if weight_ref is not None else {}
        return self._pick(np.column_stack([l1, l2]), worst, err, bias, noise, S0[target],
                          extra)

    def _pick(self, lam, worst, err, bias, noise, st, extra=None):
        order = np.argsort(worst)
        out, used = [], []
        for i in order:
            if not np.isfinite(worst[i]) or len(out) >= self.inp.keep_per_family:
                break
            if any(np.all(np.abs(lam[i] - lu) < 5) for lu in used):
                continue
            used.append(lam[i])
            out.append({"lam": [float(v) for v in lam[i]], "robust": float(worst[i]),
                        "err": float(err[i]), "bias": float(bias[i]),
                        "noise": float(noise[i]), "sens": float(abs(st[i])),
                        **{k: float(v[i]) for k, v in (extra or {}).items()}})
        return out


def _fmt_wl(v: float) -> str:
    return f"{v:.1f}".rstrip("0").rstrip(".")


def screen(inp: OptimizerInput, progress: Callable[[int, int], Any] | None = None
           ) -> list[Candidate]:
    """Fast screening of univariate processing families for every compound."""
    sc = _Screen(inp)
    fam = inp.families or set(FAMILIES)
    comps = inp.compounds
    jobs: list[tuple] = []
    for t in comps:
        others = [c for c in comps if c != t]
        if "zero" in fam:
            jobs.append((t, "zero", [], "single", None, "λ {0} nm"))
        if "derivative" in fam:
            for n in (1, 2, 3, 4):
                for dl in (2.0, 4.0, 8.0):
                    jobs.append((t, "derivative",
                                 [{"op": "derivative", "params": {"order": n, "delta_lambda": dl}}],
                                 "single", None, f"D{n}, Δλ {dl:g} nm, at {{0}} nm"))
        if len(others) == 1:
            o = others[0]
            if "dual" in fam:
                jobs.append((t, "dual", [], "pair", None, "ΔA {0} − {1} nm"))
            if "induced_dual" in fam:
                jobs.append((t, "induced_dual", [], "pair", o,
                             f"{{0}} − F·{{1}} nm (F from {o})"))
            region = _valid_region(sc.grid, sc.units[o])
            crop = ([{"op": "crop", "params": {"start": region[0], "end": region[1]}}]
                    if region else [])
            div = [{"op": "divide", "params": {"reference": DIV + o}}]
            if "ratio_difference" in fam:
                jobs.append((t, "ratio_difference", crop + div, "pair", None,
                             f"÷ {o}, ΔP {{0}} − {{1}} nm"))
            if "derivative_ratio" in fam:
                for dl in (2.0, 4.0, 8.0):
                    jobs.append((t, "derivative_ratio", crop + div + [
                        {"op": "derivative", "params": {"order": 1, "delta_lambda": dl}}],
                        "single", None, f"÷ {o}, D1 Δλ {dl:g} nm, at {{0}} nm"))
            if "mean_centering" in fam and region:
                jobs.append((t, "mean_centering", crop + div + [
                    {"op": "mean_center", "params": {"start": region[0], "end": region[1]}}],
                    "single", None, f"÷ {o}, mean-centred {region[0]:g}–{region[1]:g} nm, "
                                    "at {0} nm"))
            if "ratio_subtraction" in fam:
                pl = _plateau(sc.grid, sc.units[t], sc.units[o])
                if pl:   # o extends alone → ratio subtraction recovers t
                    jobs.append((t, "ratio_subtraction", [{"op": "ratio_subtraction", "params": {
                        "divisor": DIV + o, "start": pl[0], "end": pl[1]}}], "single", None,
                        f"ratio subtraction ÷ {o} (plateau {pl[0]:g}–{pl[1]:g} nm), at {{0}} nm"))
                pl2 = _plateau(sc.grid, sc.units[o], sc.units[t])
                if pl2:  # t extends alone → constant multiplication recovers t
                    jobs.append((t, "ratio_subtraction", [{"op": "constant_multiplication",
                                                            "params": {"divisor": DIV + t,
                                                                       "start": pl2[0],
                                                                       "end": pl2[1]}}],
                                 "single", None,
                                 f"constant multiplication ÷ {t} (plateau {pl2[0]:g}–"
                                 f"{pl2[1]:g} nm), at {{0}} nm"))
        if len(others) == 2:
            for z, y in (others, others[::-1]):
                region = _valid_region(sc.grid, sc.units[z])
                crop = ([{"op": "crop", "params": {"start": region[0], "end": region[1]}}]
                        if region else [])
                if "dual_amplitude" in fam:
                    jobs.append((t, "dual_amplitude",
                                 crop + [{"op": "divide", "params": {"reference": DIV + z}}],
                                 "pair", y, f"÷ {z}, P{{0}} − F·P{{1}} nm (cancels {y})"))
            if "double_divisor" in fam:
                y, z = others
                region = _valid_region(sc.grid, sc.units[y] + sc.units[z])
                crop = ([{"op": "crop", "params": {"start": region[0], "end": region[1]}}]
                        if region else [])
                for dl in (2.0, 4.0):
                    jobs.append((t, "double_divisor", crop + [
                        {"op": "divide_sum", "params": {"reference": DIV + y,
                                                        "reference2": DIV + z}},
                        {"op": "derivative", "params": {"order": 1, "delta_lambda": dl}}],
                        "single", None, f"÷ ({y}+{z}), D1 Δλ {dl:g} nm, at {{0}} nm"))
    out: list[Candidate] = []
    for k, (t, family, steps, kind, ref, label) in enumerate(jobs):
        if progress is not None and progress(k, len(jobs)) is False:
            raise InterruptedError("cancelled")
        try:
            picks = sc.single(t, steps) if kind == "single" else sc.pairs(t, steps, ref)
        except Exception as exc:  # a family may not apply to these spectra
            out.append(Candidate(t, family, FAMILIES[family], steps, error=str(exc)))
            continue
        for p in picks:
            lam = p["lam"]
            if kind == "single":
                meas = {"kind": "amplitude", "params": {"w1": lam[0]}}
            elif ref is None:
                meas = {"kind": "difference", "params": {"w1": lam[0], "w2": lam[1]}}
            else:
                meas = {"kind": "weighted_difference",
                        "params": {"w1": lam[0], "w2": lam[1], "reference": DIV + ref}}
            out.append(Candidate(
                t, family, f"{FAMILIES[family]}: " + label.format(*[_fmt_wl(v) for v in lam]),
                steps, meas, predicted_error=p["err"], robust_error=p["robust"],
                interference=p["bias"], noise=p["noise"], sensitivity=p["sens"]))
    return out


# --------------------------------------------------------------------------- #
# Verification
# --------------------------------------------------------------------------- #
def simulated_mixtures(inp: OptimizerInput, n_rep: int = 3) -> list[Spectrum]:
    """Mixtures built from the unit spectra + realistic noise, for every scenario."""
    grid, units, sigma, offset = noise_model(inp)
    rng = np.random.default_rng(inp.seed + 1)
    out = []
    for i, m in enumerate(default_scenarios(inp)):
        for r in range(n_rep):
            y = sum(units[c] * m[c] for c in inp.compounds)
            y = y + rng.normal(0, 1, grid.size) * sigma + rng.normal(0, offset)
            out.append(Spectrum(grid, y, name=f"sim {i + 1}.{r + 1}",
                                concentrations=dict(m)))
    return out


def _rel_rmsep(pred, true) -> float:
    pred, true = np.asarray(pred, float), np.asarray(true, float)
    ok = true > 0
    if not ok.any():
        return math.nan
    return float(100 * np.sqrt(np.mean(((pred[ok] - true[ok]) / true[ok]) ** 2)))


def verify(cands: list[Candidate], inp: OptimizerInput, resolve,
           progress: Callable[[int, int], Any] | None = None) -> list[Candidate]:
    """Calibrate each candidate on the real standards and predict the real and
    simulated mixtures (end-to-end)."""
    sims = simulated_mixtures(inp)
    real = [m for m in inp.mixtures]
    for k, c in enumerate(cands):
        if progress is not None and progress(k, len(cands)) is False:
            raise InterruptedError("cancelled")
        if c.measurement is None:
            continue
        m = UnivariateMethod(c.label, c.compound, c.steps, c.measurement)
        try:
            reg = m.calibrate(inp.standards[c.compound], resolve)
            c.r, c.slope, c.lod = reg.r, reg.slope, reg.lod
            sp = [m.predict(s, resolve) for s in sims]
            c.sim_rmsep = _rel_rmsep(sp, [s.concentrations[c.compound] for s in sims])
            rr = [s for s in real if s.concentrations.get(c.compound)]
            if rr:
                f = np.array([m.predict(s, resolve) for s in rr])
                t = np.array([s.concentrations[c.compound] for s in rr])
                rec = 100 * f / t
                c.real_rmsep = _rel_rmsep(f, t)
                c.real_mean_recovery = float(rec.mean())
                c.real_rsd = float(rec.std(ddof=1)) if rec.size > 1 else math.nan
        except Exception as exc:
            c.error = str(exc)
            c.sim_rmsep = math.inf
    return cands


def multivariate(inp: OptimizerInput, resolve) -> list[Candidate]:
    """Vierordt, bivariate (binary), CLS and PLS candidates for all compounds."""
    comps = inp.compounds
    fam = inp.families or set(FAMILIES)
    pure = [s for c in comps for s in inp.standards[c]]
    for s in pure:  # absent compounds are 0 in pure standards
        for c in comps:
            s.concentrations.setdefault(c, 0.0)
    sims = simulated_mixtures(inp)
    grid, units, _ = unit_spectra(inp)
    models: list[tuple[str, str, Any]] = []
    if "vierordt" in fam:
        lmax = [_band_max(grid, units[c]) for c in comps]
        if len(set(lmax)) == len(lmax):
            models.append(("vierordt", "λmax " + " / ".join(_fmt_wl(v) for v in lmax) + " nm",
                           SignalEquations(comps, [{"kind": "amplitude", "params": {"w1": v}}
                                                   for v in lmax])))
    if "bivariate" in fam and len(comps) == 2:
        pure_units = {c: Spectrum(grid, units[c]) for c in comps}
        best = kaiser_selection(pure_units, np.arange(grid[0], grid[-1], 2.0), top=1)[0]
        models.append(("bivariate", f"{_fmt_wl(best['w1'])} / {_fmt_wl(best['w2'])} nm",
                       SignalEquations(comps, [{"kind": "amplitude", "params": {"w1": best["w1"]}},
                                               {"kind": "amplitude", "params": {"w1": best["w2"]}}],
                                       intercept=True)))
    if "cls" in fam:
        rng = (float(grid[0]), float(grid[-1]))
        models.append(("cls", f"{_fmt_wl(rng[0])}–{_fmt_wl(rng[1])} nm",
                       SpectralModel("CLS", comps, ranges=[rng])))
    train = [s for s in inp.training if all(s.concentrations.get(c) is not None for c in comps)]
    if "pls" in fam and len(train) >= len(comps) + 3:
        models.append(("pls", "", SpectralModel("PLS2", comps, ranges=[(float(grid[0]),
                                                                        float(grid[-1]))])))
    out = []
    for family, detail, model in models:
        try:
            if family == "pls":
                cv = model.cross_validate(train, resolve, "venetian", 5,
                                          max_components=min(10, len(train) - 2))
                model.n_components = cv["suggested_components"]
                model.fit(train, resolve)
                detail = f"{model.n_components} LVs, {_fmt_wl(grid[0])}–{_fmt_wl(grid[-1])} nm"
            elif isinstance(model, SignalEquations):
                model.fit(pure, resolve)
            else:
                model.fit(pure, resolve)
            ps = model.predict(sims, resolve)
            real = [m for m in inp.mixtures if all(m.concentrations.get(c) for c in comps)]
            pr = model.predict(real, resolve) if real else None
        except Exception as exc:
            for c in comps:
                out.append(Candidate(c, family, FAMILIES[family], error=str(exc)))
            continue
        for j, c in enumerate(comps):
            cand = Candidate(c, family, f"{FAMILIES[family]}: {detail}",
                             model=model.to_dict(include_fit=True)
                             if isinstance(model, SpectralModel) else model.to_dict())
            cand.sim_rmsep = _rel_rmsep(ps[:, j], [s.concentrations[c] for s in sims])
            cand.robust_error = cand.predicted_error = cand.sim_rmsep
            if pr is not None:
                t = np.array([m.concentrations[c] for m in real])
                cand.real_rmsep = _rel_rmsep(pr[:, j], t)
                rec = 100 * pr[:, j] / t
                cand.real_mean_recovery = float(rec.mean())
                cand.real_rsd = float(rec.std(ddof=1)) if rec.size > 1 else math.nan
            out.append(cand)
    return out


# --------------------------------------------------------------------------- #
# Progressive methods (one divisor / a chain of absorption factors)
# --------------------------------------------------------------------------- #
def _shifted(s: Spectrum, d: float) -> Spectrum:
    """The spectrum as an instrument with a wavelength error d would read it."""
    return s.with_values(np.interp(s.wavelengths - d, s.wavelengths, s.values))


class _ProgressiveScreen:
    """Builds progressive-method configurations from the unit spectra and
    predicts their error per compound. Everything is evaluated by running the
    real method: fitted on noise-free standards made from the unit spectra,
    then applied to the scenario mixtures (interference, which should be ~0
    by construction), to noise spectra (noise) and to the scenario mixtures
    read with a ±λ error (robustness)."""

    def __init__(self, sc: _Screen, n_noise: int = 20):
        self.sc, self.inp = sc, sc.inp
        self.comps = self.inp.compounds
        self.grid = sc.grid
        self.units = sc.units
        levels = {c: sorted({s.concentrations[c] for s in self.inp.standards[c]})
                  for c in self.comps}
        self.std = {c: [Spectrum(self.grid, self.units[c] * v, concentrations={c: v})
                        for v in (levels[c] if len(levels[c]) >= 3
                                  else np.linspace(min(levels[c]), max(levels[c]) * 2, 3))]
                    for c in self.comps}
        self.mixes = [Spectrum(self.grid, sum(self.units[c] * m[c] for c in self.comps),
                               concentrations=dict(m)) for m in sc.scen]
        self.noise = sc.noise[:n_noise]
        self.zero = Spectrum(self.grid, np.zeros(self.grid.size))

    def evaluate(self, method, divisor: Spectrum | None) -> dict[str, dict] | None:
        try:
            if divisor is not None:
                method.divisor = divisor
            method.fit_spectra([x for c in self.comps for x in self.std[c]],
                               lambda ref: divisor)

            def pred(spectra):
                return method.predict_spectra(spectra, lambda ref: divisor)
            p0 = pred(self.mixes)
            pn = pred(self.noise) - pred([self.zero])
            shifts = [pred([_shifted(m, d) for m in self.mixes])
                      for d in (-self.inp.wl_uncertainty, self.inp.wl_uncertainty)
                      if self.inp.wl_uncertainty > 0]
        except Exception:
            return None
        out = {}
        sd = pn.std(axis=0, ddof=1)
        for j, c in enumerate(self.comps):
            true = np.array([m.concentrations[c] for m in self.mixes])
            ok = true > 0
            if not ok.any():
                return None
            bias = 100 * np.sqrt(np.mean(((p0[ok, j] - true[ok]) / true[ok]) ** 2))
            noise = 100 * np.sqrt(np.mean((sd[j] / true[ok]) ** 2))
            err = math.hypot(bias, noise)
            robust = err
            for ps in shifts:
                b = 100 * np.sqrt(np.mean(((ps[ok, j] - true[ok]) / true[ok]) ** 2))
                robust = max(robust, math.hypot(b, noise))
            if not np.isfinite(robust):
                return None
            out[c] = {"err": float(err), "robust": float(robust), "bias": float(bias),
                      "noise": float(noise)}
        return out

    # ---------------------------------------------------------- generation
    def _lambdas(self, lo, hi, step=4.0):
        u = self.inp.wl_uncertainty + 1.0
        return np.arange(lo + u, hi - u + 1e-9, max(step, self.inp.pair_step))

    def amplitude_centering(self) -> list[tuple[AmplitudeCentering, str]]:
        """Configurations for every choice of divisor compound: complete
        overlap (MACM / AAC-complete / CV-AD: λ pairs, divisor compound by
        subtraction) and, when the divisor compound has a plateau, AAC-partial /
        RIDSS / AM (plateau, equality-factor difference, subtraction)."""
        out = []
        for z in self.comps:
            others = [c for c in self.comps if c != z]
            if len(others) not in (1, 2):
                continue
            div = self.sc.resolve(DIV + z)
            dv = np.interp(self.grid, div.wavelengths, div.values)
            pl = _plateau(self.grid, sum(np.abs(self.units[c]) for c in others), self.units[z])
            segments = _valid_regions(self.grid, dv, 0.05)
            if segments:
                out += self._centering(z, others, dv, segments, pl)
        return out

    def _centering(self, z, others, dv, segments, pl):
        """Every wavelength used (λc, partners, plateau) must lie where the
        divisor is usable; equal-amplitude partners are searched within the
        segment of λc (the ratio is meaningless across a divisor gap)."""
        out = []
        R = {c: Spectrum(self.grid, np.where(np.abs(dv) > 1e-12,
                                             self.units[c] / np.where(dv == 0, 1, dv), 0.0))
             for c in others}
        seg_of = {}
        lams = []
        for seg in segments:
            for lc in self._lambdas(*seg):
                seg_of[float(lc)] = seg
                lams.append(float(lc))
        valid = [float(w) for seg in segments for w in self._lambdas(*seg, step=2.0)]
        big = {c: 0.05 * max(abs(R[c].value_at(w)) for w in valid) for c in others}

        def strong(c, lam):          # c's ratio amplitude at λ is usable
            return abs(R[c].value_at(lam)) >= big[c]
        # complete overlap (MACM / AAC-complete / CV-AD)
        for lc in lams:
            lo, hi = seg_of[lc]
            options = []            # per compound: up to two partner wavelengths
            for x in others:
                rest = [y for y in others if y != x]
                if not strong(x, lc):
                    break
                if rest:
                    r = R[rest[0]]
                    seg_spec = Spectrum(*r.region(lo, hi))
                    cands = [w for w in equal_amplitude_wavelengths(seg_spec, lc, lo + 2, hi - 2)
                             if abs(w - lc) >= 4]
                else:
                    cands = [w for w in self._lambdas(lo, hi, step=2.0) if abs(w - lc) >= 4]
                if not cands:
                    break
                cands.sort(key=lambda w: -abs(R[x].value_at(lc) - R[x].value_at(w)))
                options.append([(x, float(w)) for w in cands[:2]])
            else:
                for combo in itertools.product(*options):
                    diffs = {x: {"w1": lc, "w2": w} for x, w in combo}
                    out.append((AmplitudeCentering(lc, list(self.comps), subtract=z,
                                                   divisor_compound=z, differences=diffs), z))
        # plateau of the divisor compound (AAC-partial / RIDSS / AM)
        if pl is None or not any(a <= pl[0] and pl[1] <= b for a, b in segments):
            return out
        for lc in lams:
            if lc >= pl[0] - 2:
                continue
            if len(others) == 1:
                if strong(others[0], lc):
                    out.append((AmplitudeCentering(lc, list(self.comps), subtract=others[0],
                                                   divisor_compound=z, plateau=pl), z))
                continue
            for a, b in (others, others[::-1]):
                if not (strong(a, lc) and strong(b, lc)):
                    continue
                cands = []
                for w in valid:
                    if abs(w - lc) < 4 or w >= pl[0] - 2 or not strong(b, w):
                        continue
                    f = R[b].value_at(lc) / R[b].value_at(w)
                    if abs(f) <= 20:
                        cands.append((abs(R[a].value_at(lc) - f * R[a].value_at(w)), w))
                for _, w2 in sorted(cands, reverse=True)[:2]:
                    out.append((AmplitudeCentering(
                        lc, list(self.comps), subtract=b, divisor_compound=z, plateau=pl,
                        differences={a: {"w1": lc, "w2": float(w2), "factor_from": b}}), z))
        return out

    def absorption_factor(self) -> list[AbsorptionFactorMethod]:
        """Every order of the compounds for which each one has a wavelength
        where the compounds after it in the order do not absorb."""
        cmax = {c: max(m[c] for m in self.sc.scen) for c in self.comps}
        cmin = {c: max(min(m[c] for m in self.sc.scen), 1e-9) for c in self.comps}
        u = self.inp.wl_uncertainty + 1.0
        inner = (self.grid >= self.grid[0] + u) & (self.grid <= self.grid[-1] - u)
        out = []
        for order in itertools.permutations(self.comps):
            lams = []
            for i, c in enumerate(order):
                later = order[i + 1:]
                uc = self.units[c]
                ok = inner & (uc >= 0.05 * uc.max())
                if later:
                    inter = sum(np.abs(self.units[k]) * cmax[k] for k in later)
                    ok &= inter <= 0.002 * uc * cmin[c]
                if not ok.any():
                    break
                idx = np.flatnonzero(ok)
                lams.append(float(self.grid[idx[np.argmax(uc[idx])]]))
            else:
                out.append(AbsorptionFactorMethod(list(zip(order, lams))))
        return out


def progressive(inp: OptimizerInput, sc: _Screen, resolve,
                progress: Callable[[int, int], Any] | None = None) -> list[Candidate]:
    """Screen and verify amplitude-centering and absorption-factor methods.

    Each configuration resolves all compounds at once; the best few by their
    worst compound are verified end-to-end on the real standards, simulated
    and laboratory mixtures and returned as one candidate per compound."""
    fam = inp.families or set(FAMILIES)
    ps = _ProgressiveScreen(sc)
    configs: list[tuple[str, Any, str | None]] = []
    if "amplitude_centering" in fam:
        configs += [("amplitude_centering", m, z) for m, z in ps.amplitude_centering()]
    if "absorption_factor" in fam:
        configs += [("absorption_factor", m, None) for m in ps.absorption_factor()]
    scored = []
    for k, (family, m, z) in enumerate(configs):
        if progress is not None and progress(k, len(configs)) is False:
            raise InterruptedError("cancelled")
        div = sc.resolve(DIV + z) if z else None
        ev = ps.evaluate(m, div)
        if ev is not None:
            scored.append((max(e["robust"] for e in ev.values()), family, m, z, ev))
    # best configurations of each kind (family, divisor, structure)
    best: dict[tuple, list] = {}
    for item in sorted(scored, key=lambda t: t[0]):
        _, family, m, z, _ = item
        key = (family, z, bool(getattr(m, "plateau", None)),
               tuple(sorted(getattr(m, "differences", {}))))
        if len(best.setdefault(key, [])) < inp.keep_per_family:
            best[key].append(item)
    sims = simulated_mixtures(inp)
    real = [x for x in inp.mixtures if all(x.concentrations.get(c) for c in inp.compounds)]
    pure = [x for c in inp.compounds for x in inp.standards[c]]
    out = []
    for items in best.values():
        for _, family, m, z, ev in items:
            m = progressive_from_dict(m.to_dict())        # fresh copy, no unit-spectrum fit
            m.divisor = DIV + z if z else None
            label = f"{FAMILIES[family]}: {m.describe()}"
            try:
                m.fit_spectra(pure, resolve)
                psim = m.predict_spectra(sims, resolve)
                preal = m.predict_spectra(real, resolve) if real else None
            except Exception as exc:
                for c in inp.compounds:
                    out.append(Candidate(c, family, label, model=m.to_dict(), error=str(exc)))
                continue
            for j, c in enumerate(inp.compounds):
                cand = Candidate(c, family, label, model=m.to_dict(),
                                 predicted_error=ev[c]["err"], robust_error=ev[c]["robust"],
                                 interference=ev[c]["bias"], noise=ev[c]["noise"])
                cand.r = m.regressions[c].r
                cand.slope = m.regressions[c].slope
                cand.sim_rmsep = _rel_rmsep(psim[:, j], [x.concentrations[c] for x in sims])
                if preal is not None:
                    t = np.array([x.concentrations[c] for x in real])
                    cand.real_rmsep = _rel_rmsep(preal[:, j], t)
                    rec = 100 * preal[:, j] / t
                    cand.real_mean_recovery = float(rec.mean())
                    cand.real_rsd = float(rec.std(ddof=1)) if rec.size > 1 else math.nan
                out.append(cand)
    return out


def optimize(inp: OptimizerInput, progress: Callable[[int, int], Any] | None = None) -> dict:
    """Screen + verify. Returns ranked candidates per compound and a summary."""
    def sub(offset, span):
        if progress is None:
            return None
        return lambda d, t: progress(offset + int(span * d / max(t, 1)), 100)
    for c in inp.compounds:
        for s in inp.standards[c]:
            s.concentrations = {**{k: 0.0 for k in inp.compounds}, **s.concentrations}
    cands = screen(inp, sub(0, 50))
    sc = _Screen(inp)
    ok = [c for c in cands if not c.error]
    verify(ok, inp, sc.resolve, sub(50, 25))
    prog = progressive(inp, sc, sc.resolve, sub(75, 15))
    multi = multivariate(inp, sc.resolve)
    if progress is not None:
        progress(100, 100)
    everything = ok + [c for c in prog if not c.error] + multi + [c for c in cands if c.error] \
        + [c for c in prog if c.error]
    ranked = {c: sorted([x for x in everything if x.compound == c], key=lambda x: x.score)
              for c in inp.compounds}
    return {"ranked": ranked,
            "best": {c: (r[0] if r else None) for c, r in ranked.items()},
            "scenarios": default_scenarios(inp),
            "grid": sc.grid, "units": sc.units, "noise": sc.sigma}


def materialize(c: Candidate, divisor_ids: dict[str, int]) -> dict:
    """Method definition with divisor placeholders replaced by spectrum ids."""
    def fix(v):
        if isinstance(v, str) and v.startswith(DIV):
            return divisor_ids[v[len(DIV):]]
        return v
    steps = [{"op": s["op"], "params": {k: fix(v) for k, v in s["params"].items()}}
             for s in c.steps]
    meas = None
    if c.measurement:
        meas = {"kind": c.measurement["kind"],
                "params": {k: fix(v) for k, v in c.measurement["params"].items()}}
    out = {"steps": steps, "measurement": meas}
    if c.model is not None:
        model = dict(c.model)
        if "divisor" in model:
            model["divisor"] = fix(model["divisor"])
        out["model"] = model
    return out


def explain(c: Candidate, inp: OptimizerInput, resolve) -> dict[str, Spectrum]:
    """Processed unit spectra of all compounds for plotting why a candidate works."""
    grid, units, _ = unit_spectra(inp)
    return {k: apply_pipeline(Spectrum(grid, u, name=f"{k} (1 unit)"), c.steps, resolve)
            for k, u in units.items()}


__all__ = ["Candidate", "OptimizerInput", "FAMILIES", "optimize", "screen", "verify",
           "multivariate", "progressive", "materialize", "explain", "simulated_mixtures", "unit_spectra",
           "measure", "linear_regression"]
