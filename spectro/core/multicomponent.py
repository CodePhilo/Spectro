"""Multicomponent determination: equation methods and chemometrics.

* :class:`SignalEquations` — Vierordt simultaneous equations, multi-wavelength
  least squares, bivariate calibration (with intercepts), AUC equations.
* :class:`SpectralModel` — full-spectrum chemometrics: CLS, ILS/MLR, PCR,
  PLS-1, PLS-2, MCR-ALS, ANN (MLP) and SVR, with cross-validation.
* Variable selection (iPLS, GA-PLS, VIP), diagnostics (T², Q), and the
  Brereton multilevel multifactor calibration design.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np
from scipy.optimize import nnls

from spectro.core.operations import Resolver, apply_pipeline
from spectro.core.spectrum import Spectrum, align
from spectro.core.univariate import measure
from spectro.core.validation import prediction_error


def concentration_matrix(spectra: list[Spectrum], compounds: list[str]) -> np.ndarray:
    missing = [(s.name, c) for s in spectra for c in compounds if c not in s.concentrations]
    if missing:
        name, comp = missing[0]
        raise ValueError(f"'{name}' has no concentration for {comp} "
                         f"(use 0 for absent compounds)")
    return np.array([[s.concentrations[c] for c in compounds] for s in spectra], float)


# --------------------------------------------------------------------------- #
# Equation methods on selected signals
# --------------------------------------------------------------------------- #
@dataclass
class SignalEquations:
    """R = C·K (+ b): responses at chosen signals are linear in concentrations.

    * signals = amplitudes at n ≥ k wavelengths → Vierordt / multi-wavelength
      least squares (``intercept=False``), bivariate (``intercept=True``).
    * signals = areas → area-under-curve equations.
    * ``steps`` apply a pipeline first (e.g. derivative) if wanted.
    """

    compounds: list[str]
    signals: list[dict]
    steps: list[dict] = field(default_factory=list)
    intercept: bool = False
    K: np.ndarray | None = None
    b: np.ndarray | None = None

    def responses(self, spectra: list[Spectrum], resolve: Resolver | None = None) -> np.ndarray:
        out = []
        for s in spectra:
            p = apply_pipeline(s, self.steps, resolve)
            out.append([measure(p, sig, resolve) for sig in self.signals])
        return np.array(out, float)

    def fit(self, spectra: list[Spectrum], resolve: Resolver | None = None) -> dict:
        if len(self.signals) < len(self.compounds):
            raise ValueError("need at least as many signals as compounds")
        C = concentration_matrix(spectra, self.compounds)
        R = self.responses(spectra, resolve)
        A = np.hstack([C, np.ones((C.shape[0], 1))]) if self.intercept else C
        coef, *_ = np.linalg.lstsq(A, R, rcond=None)
        self.K = coef[: len(self.compounds)]
        self.b = coef[-1] if self.intercept else np.zeros(R.shape[1])
        cond = float(np.linalg.cond(self.K)) if self.K.shape[0] <= self.K.shape[1] else float("inf")
        return {"K": self.K.tolist(), "b": self.b.tolist(), "condition_number": cond}

    def predict(self, spectra: list[Spectrum], resolve: Resolver | None = None) -> np.ndarray:
        if self.K is None:
            raise RuntimeError("not fitted")
        R = self.responses(spectra, resolve) - self.b
        c, *_ = np.linalg.lstsq(self.K.T, R.T, rcond=None)
        return c.T

    def to_dict(self) -> dict:
        return {"type": "equations", "compounds": self.compounds, "signals": self.signals,
                "steps": self.steps, "intercept": self.intercept,
                "K": None if self.K is None else self.K.tolist(),
                "b": None if self.b is None else self.b.tolist()}

    @classmethod
    def from_dict(cls, d: dict) -> "SignalEquations":
        m = cls(d["compounds"], d["signals"], d.get("steps", []), d.get("intercept", False))
        if d.get("K") is not None:
            m.K, m.b = np.array(d["K"]), np.array(d["b"])
        return m


def kaiser_selection(pure: dict[str, Spectrum], candidates: np.ndarray | None = None,
                     top: int = 10) -> list[dict]:
    """Bivariate wavelength-pair selection (Kaiser): maximise |det| of the
    sensitivity matrix of two pure components (unit concentration)."""
    names = list(pure)
    if len(names) != 2:
        raise ValueError("Kaiser selection needs exactly two pure spectra")
    a, b = align([pure[names[0]], pure[names[1]]])
    x = a.wavelengths if candidates is None else np.asarray(candidates)
    sa, sb = np.interp(x, a.x, a.y), np.interp(x, b.x, b.y)
    out = []
    for i in range(len(x)):
        for j in range(i + 1, len(x)):
            det = sa[i] * sb[j] - sa[j] * sb[i]
            out.append({"w1": float(x[i]), "w2": float(x[j]), "det": float(abs(det))})
    out.sort(key=lambda d: -d["det"])
    return out[:top]


# --------------------------------------------------------------------------- #
# Full-spectrum chemometrics
# --------------------------------------------------------------------------- #
MODEL_TYPES = ("CLS", "ILS", "PCR", "PLS1", "PLS2", "MCR-ALS", "ANN", "SVR")


@dataclass
class SpectralModel:
    model_type: str
    compounds: list[str]
    steps: list[dict] = field(default_factory=list)
    ranges: list[tuple[float, float]] = field(default_factory=list)
    wavelengths: list[float] = field(default_factory=list)  # ILS: selected points
    n_components: int = 2
    preprocessing: str = "mean_center"  # none | mean_center | autoscale
    options: dict[str, Any] = field(default_factory=dict)
    # fitted state
    grid: np.ndarray | None = None
    state: dict[str, Any] = field(default_factory=dict)

    # ---------------------------------------------------------------- data
    def matrix(self, spectra: list[Spectrum], resolve: Resolver | None = None) -> np.ndarray:
        proc = [apply_pipeline(s, self.steps, resolve) for s in spectra]
        if self.grid is None:
            aligned = align(proc)
            grid = aligned[0].wavelengths
            mask = self._mask(grid)
            self.grid = grid[mask]
        return np.vstack([np.interp(self.grid, p.wavelengths, p.values) for p in proc])

    def _mask(self, grid: np.ndarray) -> np.ndarray:
        if self.model_type == "ILS" and self.wavelengths:
            idx = [int(np.argmin(np.abs(grid - w))) for w in self.wavelengths]
            mask = np.zeros(grid.size, bool)
            mask[idx] = True
            return mask
        if not self.ranges:
            return np.ones(grid.size, bool)
        mask = np.zeros(grid.size, bool)
        for lo, hi in self.ranges:
            mask |= (grid >= min(lo, hi) - 1e-9) & (grid <= max(lo, hi) + 1e-9)
        if not mask.any():
            raise ValueError("the selected wavelength ranges contain no data")
        return mask

    def _prep_fit(self, X: np.ndarray) -> np.ndarray:
        mu = X.mean(0) if self.preprocessing in ("mean_center", "autoscale") else np.zeros(X.shape[1])
        sd = X.std(0, ddof=1) if self.preprocessing == "autoscale" else np.ones(X.shape[1])
        sd[sd == 0] = 1.0
        self.state["x_mean"], self.state["x_sd"] = mu, sd
        return (X - mu) / sd

    def _prep(self, X: np.ndarray) -> np.ndarray:
        return (X - self.state["x_mean"]) / self.state["x_sd"]

    # ---------------------------------------------------------------- fit
    def fit(self, spectra: list[Spectrum], resolve: Resolver | None = None) -> "SpectralModel":
        self.grid = None
        self.state = {}
        X = self.matrix(spectra, resolve)
        Y = concentration_matrix(spectra, self.compounds)
        self._fit_xy(X, Y)
        self._finish_fit(X, Y)
        return self

    def _finish_fit(self, X: np.ndarray, Y: np.ndarray) -> None:
        self.state["pca"] = self._fit_pca_limits(X)
        dof = max(1, X.shape[0] - (1 if self.model_type == "CLS" else
                                   min(int(self.n_components), X.shape[0] - 1)) - 1)
        res = self._predict_x(X) - Y
        self.state["rmsec"] = np.sqrt(np.sum(res ** 2, axis=0) / dof)

    def _fit_xy(self, X: np.ndarray, Y: np.ndarray) -> None:
        self._train_xy(X, Y)
        self.state["frozen"] = self._freeze(X.shape[1])

    def _train_xy(self, X: np.ndarray, Y: np.ndarray) -> None:
        t = self.model_type
        k = int(self.n_components)
        if t == "CLS":
            # X = C·K (+ mean terms). Centering both keeps an intercept.
            center = self.preprocessing != "none"
            xm = X.mean(0) if center else np.zeros(X.shape[1])
            ym = Y.mean(0) if center else np.zeros(Y.shape[1])
            K, *_ = np.linalg.lstsq(Y - ym, X - xm, rcond=None)
            self.state.update(K=K, xm=xm, ym=ym)
            return
        Xp = self._prep_fit(X)
        ym = Y.mean(0)
        self.state["y_mean"] = ym
        if t == "ILS":
            if Xp.shape[1] >= Xp.shape[0]:
                raise ValueError("ILS needs fewer wavelengths than calibration samples; "
                                 "select wavelengths")
            B, *_ = np.linalg.lstsq(Xp, Y - ym, rcond=None)
            self.state["B"] = B
        elif t == "PCR":
            u, s, vt = np.linalg.svd(Xp, full_matrices=False)
            k = min(k, vt.shape[0])
            P = vt[:k].T
            T = Xp @ P
            Q, *_ = np.linalg.lstsq(T, Y - ym, rcond=None)
            self.state.update(P=P, B=P @ Q, explained=(s ** 2 / np.sum(s ** 2))[:k])
        elif t in ("PLS1", "PLS2"):
            from sklearn.cross_decomposition import PLSRegression

            k = min(k, min(Xp.shape) - 1 if min(Xp.shape) > 1 else 1)
            if t == "PLS2":
                pls = PLSRegression(n_components=k, scale=False).fit(Xp, Y - ym)
                self.state["models"] = [pls]
            else:
                self.state["models"] = [PLSRegression(n_components=k, scale=False)
                                        .fit(Xp, (Y - ym)[:, j:j + 1]) for j in range(Y.shape[1])]
        elif t == "MCR-ALS":
            self._fit_mcr(X, Y)
        elif t == "ANN":
            from sklearn.neural_network import MLPRegressor

            u, s, vt = np.linalg.svd(Xp, full_matrices=False)
            k = min(k, vt.shape[0])
            P = vt[:k].T
            T = Xp @ P
            ts = T.std(0, ddof=1)
            ts[ts == 0] = 1
            ysd = Y.std(0, ddof=1)
            ysd[ysd == 0] = 1
            hidden = tuple(int(h) for h in str(self.options.get("hidden", "10")).split(",") if h.strip())
            # an ensemble of networks from different random starts is averaged:
            # a single network lands in slightly different minima depending on
            # sample order / rounding, the average is stable and more accurate
            nets = []
            for e in range(max(1, int(self.options.get("ensemble", 5)))):
                net = MLPRegressor(hidden_layer_sizes=hidden or (10,),
                                   activation=self.options.get("activation", "tanh"),
                                   solver="lbfgs",
                                   max_iter=int(self.options.get("max_iter", 5000)),
                                   alpha=float(self.options.get("alpha", 1e-4)), tol=1e-9,
                                   random_state=int(self.options.get("seed", 0)) + e)
                nets.append(net.fit(T / ts, (Y - ym) / ysd))
            self.state.update(P=P, ts=ts, ysd=ysd, net=nets[0], nets=nets)
        elif t == "SVR":
            from sklearn.svm import SVR

            # targets are standardised so that C and ε do not depend on the
            # concentration unit (ε is a fraction of each compound's SD)
            ysd = Y.std(0, ddof=1)
            ysd[ysd == 0] = 1
            models = []
            for j in range(Y.shape[1]):
                m = SVR(kernel=self.options.get("kernel", "linear"),
                        C=float(self.options.get("C", 100.0)),
                        epsilon=float(self.options.get("epsilon", 0.01)),
                        gamma=self.options.get("gamma", "scale"), tol=1e-6)
                models.append(m.fit(Xp, (Y[:, j] - ym[j]) / ysd[j]))
            self.state.update(models=models, ysd=ysd)
        else:
            raise ValueError(f"unknown model type {t}")

    def _fit_mcr(self, X: np.ndarray, Y: np.ndarray) -> None:
        """MCR-ALS with non-negativity; initial spectra from calibration
        (CLS estimate), concentrations regressed onto the reference values
        (correlation constraint)."""
        K, *_ = np.linalg.lstsq(Y, X, rcond=None)
        S = np.clip(K, 0, None)  # k × p
        C = Y.copy()
        max_iter = int(self.options.get("max_iter", 200))
        tol = float(self.options.get("tol", 1e-8))
        last = np.inf
        for _ in range(max_iter):
            C = np.array([nnls(S.T, x)[0] for x in X])
            # correlation constraint: rescale each profile to the reference
            for j in range(C.shape[1]):
                b = np.polyfit(C[:, j], Y[:, j], 1) if np.ptp(C[:, j]) > 0 else (1.0, 0.0)
                C[:, j] = np.clip(np.polyval(b, C[:, j]), 0, None)
            S = np.array([nnls(C, X[:, i])[0] for i in range(X.shape[1])]).T
            lof = np.linalg.norm(X - C @ S) / np.linalg.norm(X)
            if abs(last - lof) < tol:
                break
            last = lof
        self.state.update(S=S, lof=float(100 * last))
        raw = np.array([nnls(S.T, x)[0] for x in X])
        self.state["calib"] = [np.polyfit(raw[:, j], Y[:, j], 1) for j in range(Y.shape[1])]

    # ---------------------------------------------------------------- freeze
    def _freeze(self, p: int) -> dict:
        """Reduce the fitted model to plain arrays (JSON-serialisable).

        Prediction always runs from this representation, so a saved or
        exported model predicts exactly as it did when it was fitted."""
        t, st = self.model_type, self.state
        if t == "CLS":
            return {"kind": "cls", "K": st["K"], "xm": st["xm"], "ym": st["ym"]}
        if t == "MCR-ALS":
            return {"kind": "mcr", "S": st["S"], "calib": np.array(st["calib"]),
                    "lof": st["lof"]}
        base = {"x_mean": st["x_mean"], "x_sd": st["x_sd"], "y_mean": st["y_mean"]}
        if t in ("ILS", "PCR"):
            return {"kind": "linear", "B": st["B"], "b0": np.zeros(st["B"].shape[1]), **base}
        if t in ("PLS1", "PLS2") or (t == "SVR" and self.options.get("kernel", "linear") == "linear"):
            # these predictors are affine in the preprocessed spectrum: recover
            # the regression vector exactly from p + 1 evaluations
            f = self._raw_predictor()
            b0 = f(np.zeros((1, p)))[0]
            B = f(np.eye(p)) - b0
            fr = {"kind": "linear", "B": B, "b0": b0, **base}
            if t in ("PLS1", "PLS2"):
                fr["vip"] = self._vip_from_models()
            return fr
        if t == "ANN":
            nets = st.get("nets") or [st["net"]]
            return {"kind": "ann", "P": st["P"], "ts": st["ts"], "ysd": st["ysd"],
                    "activation": nets[0].activation,
                    "ensemble": [{"coefs": list(n.coefs_), "intercepts": list(n.intercepts_)}
                                 for n in nets], **base}
        if t == "SVR":
            return {"kind": "svr", **base, "ysd": st["ysd"], "models": [
                {"kernel": m.kernel, "sv": m.support_vectors_, "dual": m.dual_coef_.ravel(),
                 "b": float(m.intercept_[0]), "gamma": float(m._gamma),
                 "degree": int(m.degree), "coef0": float(m.coef0)}
                for m in st["models"]]}
        raise ValueError(t)

    def _raw_predictor(self):
        st = self.state
        if self.model_type == "PLS2":
            return lambda Z: st["models"][0].predict(Z).reshape(len(Z), -1)
        scale = st.get("ysd", 1.0) if self.model_type == "SVR" else 1.0
        return lambda Z: np.column_stack([m.predict(Z).ravel() for m in st["models"]]) * scale

    def _vip_from_models(self) -> np.ndarray:
        out = []
        for pls in self.state["models"]:
            t, w, q = pls.x_scores_, pls.x_weights_, pls.y_loadings_
            p = w.shape[0]
            ss = np.sum(t ** 2, axis=0) * np.sum(q ** 2, axis=0)
            wn = w / np.linalg.norm(w, axis=0)
            out.append(np.sqrt(p * (wn ** 2 @ ss) / ss.sum()))
        # PLS1: one model per compound → combine as root-mean-square so that the
        # defining identity mean(VIP²) = 1 still holds
        return np.sqrt(np.mean(np.square(out), axis=0))

    @property
    def is_fitted(self) -> bool:
        return "frozen" in self.state and self.grid is not None

    # ---------------------------------------------------------------- predict
    def predict(self, spectra: list[Spectrum], resolve: Resolver | None = None) -> np.ndarray:
        if not self.is_fitted:
            raise RuntimeError("model is not fitted")
        return self._predict_x(self.matrix(spectra, resolve))

    def _predict_x(self, X: np.ndarray) -> np.ndarray:
        fr = self.state["frozen"]
        kind = fr["kind"]
        if kind == "cls":
            c, *_ = np.linalg.lstsq(np.asarray(fr["K"]).T, (X - fr["xm"]).T, rcond=None)
            return c.T + fr["ym"]
        if kind == "mcr":
            S = np.asarray(fr["S"])
            raw = np.array([nnls(S.T, x)[0] for x in X])
            return np.column_stack([np.polyval(fr["calib"][j], raw[:, j])
                                    for j in range(raw.shape[1])])
        Xp = (X - fr["x_mean"]) / fr["x_sd"]
        if kind == "linear":
            return Xp @ fr["B"] + fr["b0"] + fr["y_mean"]
        if kind == "ann":
            t0 = (Xp @ fr["P"]) / fr["ts"]
            act = _ACTIVATIONS[fr["activation"]]
            members = fr.get("ensemble") or [{"coefs": fr["coefs"],
                                              "intercepts": fr["intercepts"]}]
            outs = []
            for net in members:
                h = t0
                n = len(net["coefs"])
                for i, (W, b) in enumerate(zip(net["coefs"], net["intercepts"])):
                    h = h @ np.asarray(W) + np.asarray(b)
                    if i < n - 1:
                        h = act(h)
                outs.append(h.reshape(len(X), -1))
            return np.mean(outs, axis=0) * fr["ysd"] + fr["y_mean"]
        if kind == "svr":
            cols = []
            for m in fr["models"]:
                sv = np.asarray(m["sv"])
                if m["kernel"] == "rbf":
                    d2 = (np.sum(Xp ** 2, 1)[:, None] + np.sum(sv ** 2, 1)[None, :]
                          - 2 * Xp @ sv.T)
                    K = np.exp(-m["gamma"] * d2)
                elif m["kernel"] == "poly":
                    K = (m["gamma"] * Xp @ sv.T + m["coef0"]) ** m["degree"]
                elif m["kernel"] == "sigmoid":
                    K = np.tanh(m["gamma"] * Xp @ sv.T + m["coef0"])
                else:
                    K = Xp @ sv.T
                cols.append(K @ np.asarray(m["dual"]) + m["b"])
            return np.column_stack(cols) * fr.get("ysd", 1.0) + fr["y_mean"]
        raise ValueError(kind)

    # ---------------------------------------------------------------- diagnostics
    def _fit_pca_limits(self, X: np.ndarray) -> dict:
        """PCA of the calibration spectra with 95 % / 99 % limits for
        Hotelling T² (F distribution) and Q residuals (Box's χ² approximation)."""
        from scipy import stats

        n = X.shape[0]
        mu = X.mean(0)
        Xc = X - mu
        k = max(1, min(int(self.n_components), min(Xc.shape) - 1))
        u, s, vt = np.linalg.svd(Xc, full_matrices=False)
        eig = s ** 2 / max(1, n - 1)
        out = {"mean": mu, "P": vt[:k].T, "lam": eig[:k], "k": k, "n": n,
               "explained": s ** 2 / np.sum(s ** 2)}
        # T²: limit for a new observation, k(n²−1)/(n(n−k))·F(k, n−k).
        # Q: Box's g·χ²(h). Mean = noise part (per-wavelength residual variance,
        #    dof-corrected) + subspace-estimation part (from leave-one-out Q,
        #    rescaled from n−1 to n samples); spread from the residual
        #    covariance, including correlated (baseline) noise. Validated to flag
        #    ≈ 5 % / 1 % of normal new samples with realistic noise; with purely
        #    white noise it is conservative.
        E = Xc - (Xc @ out["P"]) @ out["P"].T
        dof = max(1, n - 1 - k)
        sig2 = np.sum(E ** 2, axis=0) / dof
        t1 = float(np.sum(sig2))
        mean_new = t1
        if 4 <= n <= 300:
            q_loo = []
            for i in range(n):
                keep = np.arange(n) != i
                mi = X[keep].mean(0)
                Pi = np.linalg.svd(X[keep] - mi, full_matrices=False)[2][:k].T
                ei = (X[i] - mi) - ((X[i] - mi) @ Pi) @ Pi.T
                q_loo.append(float(ei @ ei))
            mean_new = t1 + max(0.0, float(np.mean(q_loo)) - t1) * (n - 1) / n
        ev = s ** 2 / dof
        ev = ev[k:]
        t2 = max(float(np.sum(sig2 ** 2)), float(np.sum(ev ** 2)) - t1 ** 2 / dof)
        if t1 > 0 and t2 > 0:
            var = 2 * t2 * (mean_new / t1) ** 2
            g, h = var / (2 * mean_new), 2 * mean_new ** 2 / var
        else:
            g, h = np.inf, 1.0
        for a in (0.95, 0.99):
            tag = str(int(a * 100))
            out[f"t2_lim{tag}"] = (k * (n * n - 1) / (n * (n - k)) * stats.f.ppf(a, k, n - k)
                                   if n > k else float("inf"))
            out[f"q_lim{tag}"] = float(g * stats.chi2.ppf(a, h))
        return out

    def _t2_q(self, X: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        pca = self.state["pca"]
        Xc = X - pca["mean"]
        T = Xc @ pca["P"]
        E = Xc - T @ np.asarray(pca["P"]).T
        return np.sum(T ** 2 / pca["lam"], axis=1), np.sum(E ** 2, axis=1), T

    def diagnostics(self, spectra: list[Spectrum], resolve: Resolver | None = None,
                    alpha: str = "95") -> dict:
        """Hotelling T² and Q residuals of ``spectra`` against the calibration
        PCA, with limits and outlier flags (and concentration residuals when
        reference values are known)."""
        if "pca" not in self.state:
            raise RuntimeError("fit the model first")
        pca = self.state["pca"]
        X = self.matrix(spectra, resolve)
        t2, q, T = self._t2_q(X)
        t2_lim, q_lim = pca[f"t2_lim{alpha}"], pca[f"q_lim{alpha}"]
        flags = []
        pred = self._predict_x(X)
        rmsec = self.state.get("rmsec")
        for i, s in enumerate(spectra):
            f = []
            if t2[i] > t2_lim:
                f.append("T²")
            if q[i] > q_lim:
                f.append("Q")
            if rmsec is not None:
                for j, c in enumerate(self.compounds):
                    ref = s.concentrations.get(c)
                    if ref is not None and rmsec[j] > 0 and abs(pred[i, j] - ref) > 3 * rmsec[j]:
                        f.append(f"{c} residual")
            flags.append(f)
        return {"scores": T.tolist(), "loadings": np.asarray(pca["P"]).tolist(),
                "T2": t2.tolist(), "Q": q.tolist(), "explained": np.asarray(pca["explained"]).tolist(),
                "grid": self.grid.tolist(), "flags": flags,
                "limits": {k: float(v) for k, v in pca.items() if k.startswith(("t2_lim", "q_lim"))},
                "predicted": pred.tolist()}

    def vip(self) -> np.ndarray | None:
        """Variable importance in projection (PLS models)."""
        fr = self.state.get("frozen", {})
        return None if fr.get("vip") is None else np.asarray(fr["vip"])

    # ---------------------------------------------------------------- CV
    def cross_validate(self, spectra: list[Spectrum], resolve: Resolver | None = None,
                       method: str = "loo", folds: int = 5,
                       max_components: int | None = None, progress=None) -> dict:
        """Cross-validation; for latent-variable models also RMSECV vs number
        of components (and the Haaland–Thomas choice, F-test α = 0.25)."""
        self.grid = None
        X = self.matrix(spectra, resolve)
        Y = concentration_matrix(spectra, self.compounds)
        n = X.shape[0]
        splits = _cv_splits(n, method, folds)
        latent = self.model_type in ("PCR", "PLS1", "PLS2", "ANN")
        ks = (range(1, (max_components or min(10, n - 2)) + 1) if latent
              else [self.n_components])
        curve = []
        best_pred = None
        saved = self.n_components
        total, done = len(ks) * len(splits), 0
        for k in ks:
            self.n_components = k
            pred = np.zeros_like(Y)
            for test in splits:
                train = np.setdiff1d(np.arange(n), test)
                self.state = {}
                self._fit_xy(X[train], Y[train])
                pred[test] = self._predict_x(X[test])
                done += 1
                if progress is not None and progress(done, total) is False:
                    raise InterruptedError("cancelled")
            press = np.sum((pred - Y) ** 2, axis=0)
            curve.append({"k": k, "PRESS": press.tolist(),
                          "RMSECV": np.sqrt(press / n).tolist()})
            if k == saved or best_pred is None:
                best_pred = pred if k == saved else best_pred
        self.n_components = saved
        if best_pred is None:
            best_pred = pred
        suggested = saved
        if latent and curve:
            from scipy import stats

            tot = np.array([np.sum(c["PRESS"]) for c in curve])
            kmin = int(np.argmin(tot))
            suggested = curve[kmin]["k"]
            for i in range(kmin + 1):
                f = tot[i] / tot[kmin]
                if stats.f.sf(f, n, n) > 0.25:
                    suggested = curve[i]["k"]
                    break
        self.state = {}
        self._fit_xy(X, Y)
        self._finish_fit(X, Y)
        return {"curve": curve, "suggested_components": suggested,
                "predicted": best_pred.tolist(), "actual": Y.tolist(),
                "stats": {c: prediction_error(best_pred[:, j], Y[:, j])
                          for j, c in enumerate(self.compounds)}}

    # ---------------------------------------------------------------- persist
    def to_dict(self, include_fit: bool = False) -> dict:
        d = {"type": "spectral", "model_type": self.model_type,
             "compounds": self.compounds, "steps": self.steps,
             "ranges": [list(r) for r in self.ranges], "wavelengths": self.wavelengths,
             "n_components": self.n_components, "preprocessing": self.preprocessing,
             "options": self.options}
        if include_fit and self.is_fitted:
            d["fitted"] = {"grid": _plain(self.grid), "frozen": _plain(self.state["frozen"]),
                           "pca": _plain(self.state.get("pca")),
                           "rmsec": _plain(self.state.get("rmsec"))}
        return d

    @classmethod
    def from_dict(cls, d: dict) -> "SpectralModel":
        m = cls(d["model_type"], d["compounds"], d.get("steps", []),
                [tuple(r) for r in d.get("ranges", [])], d.get("wavelengths", []),
                d.get("n_components", 2), d.get("preprocessing", "mean_center"),
                d.get("options", {}))
        fit = d.get("fitted")
        if fit:
            m.grid = np.asarray(fit["grid"], float)
            m.state = {"frozen": _arrays(fit["frozen"])}
            if fit.get("pca"):
                m.state["pca"] = _arrays(fit["pca"])
            if fit.get("rmsec") is not None:
                m.state["rmsec"] = np.asarray(fit["rmsec"], float)
        return m


_ACTIVATIONS = {"identity": lambda x: x, "logistic": lambda x: 1 / (1 + np.exp(-x)),
                "tanh": np.tanh, "relu": lambda x: np.maximum(x, 0)}

def _plain(obj):
    """numpy → nested lists (JSON)."""
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, (np.floating, np.integer)):
        return obj.item()
    if isinstance(obj, dict):
        return {k: _plain(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_plain(v) for v in obj]
    return obj


def _arrays(obj, key: str = ""):
    """Nested lists → numpy arrays (inverse of :func:`_plain`)."""
    if isinstance(obj, dict):
        return {k: _arrays(v, k) for k, v in obj.items()}
    if isinstance(obj, list):
        if key in ("coefs", "intercepts"):
            return [np.asarray(v, float) for v in obj]
        if key in ("models", "ensemble"):
            return [_arrays(v) for v in obj]
        try:
            return np.asarray(obj, float)
        except (ValueError, TypeError):
            return [_arrays(v) for v in obj]
    return obj


def _cv_splits(n: int, method: str, folds: int) -> list[np.ndarray]:
    idx = np.arange(n)
    if method == "loo" or folds >= n:
        return [np.array([i]) for i in idx]
    if method == "venetian":
        return [idx[i::folds] for i in range(folds)]
    if method == "contiguous":
        return [a for a in np.array_split(idx, folds)]
    rng = np.random.default_rng(0)
    return [a for a in np.array_split(rng.permutation(idx), folds)]


# --------------------------------------------------------------------------- #
# Variable selection
# --------------------------------------------------------------------------- #
def ipls(model: SpectralModel, spectra: list[Spectrum], resolve: Resolver | None = None,
         intervals: int = 10, progress=None) -> list[dict]:
    """Interval PLS: RMSECV of the model on each of ``intervals`` equal regions."""
    base = SpectralModel.from_dict(model.to_dict())
    base.ranges = []
    base.grid = None
    X = base.matrix(spectra, resolve)
    grid = base.grid
    edges = np.array_split(np.arange(grid.size), intervals)
    out = []
    for n_done, e in enumerate(edges, 1):
        if progress is not None and progress(n_done - 1, len(edges)) is False:
            raise InterruptedError("cancelled")
        if e.size < 2:
            continue
        m = SpectralModel.from_dict(model.to_dict())
        m.ranges = [(float(grid[e[0]]), float(grid[e[-1]]))]
        m.n_components = min(m.n_components, e.size - 1)
        try:
            cv = m.cross_validate(spectra, resolve, "venetian", 5, max_components=m.n_components)
            rm = float(np.mean(cv["curve"][-1]["RMSECV"]))
        except Exception as exc:  # pragma: no cover - reported to user
            rm = float("nan")
            out.append({"range": m.ranges[0], "RMSECV": rm, "error": str(exc)})
            continue
        out.append({"range": m.ranges[0], "RMSECV": rm})
    del X
    return out


def ga_select(model: SpectralModel, spectra: list[Spectrum], resolve: Resolver | None = None,
              intervals: int = 20, population: int = 20, generations: int = 20,
              seed: int = 0, progress=None) -> dict:
    """Genetic-algorithm interval selection (GA-PLS style) minimising RMSECV."""
    rng = np.random.default_rng(seed)
    base = SpectralModel.from_dict(model.to_dict())
    base.ranges = []
    base.matrix(spectra, resolve)
    grid = base.grid
    blocks = [b for b in np.array_split(np.arange(grid.size), intervals) if b.size]

    def fitness(mask: np.ndarray) -> float:
        if not mask.any():
            return np.inf
        m = SpectralModel.from_dict(model.to_dict())
        m.ranges = [(float(grid[b[0]]), float(grid[b[-1]])) for b, on in zip(blocks, mask) if on]
        try:
            cv = m.cross_validate(spectra, resolve, "venetian", 5,
                                  max_components=m.n_components)
            return float(np.mean(cv["curve"][-1]["RMSECV"]))
        except Exception:
            return np.inf

    pop = rng.random((population, len(blocks))) < 0.5
    scores = np.array([fitness(p) for p in pop])
    for gen in range(generations):
        if progress is not None and progress(gen, generations) is False:
            raise InterruptedError("cancelled")
        order = np.argsort(scores)
        parents = pop[order[: population // 2]]
        children = []
        while len(children) < population - len(parents):
            a, b = parents[rng.integers(len(parents), size=2)]
            cut = rng.integers(1, len(blocks))
            child = np.concatenate([a[:cut], b[cut:]])
            flip = rng.random(len(blocks)) < 1.0 / len(blocks)
            children.append(child ^ flip)
        pop = np.vstack([parents, children])
        scores = np.concatenate([scores[order[: population // 2]],
                                 [fitness(c) for c in children]])
    best = pop[int(np.argmin(scores))]
    return {"ranges": [(float(grid[b[0]]), float(grid[b[-1]]))
                       for b, on in zip(blocks, best) if on],
            "RMSECV": float(np.min(scores))}


# --------------------------------------------------------------------------- #
# Calibration design
# --------------------------------------------------------------------------- #
def _gf5_msequence() -> list[int]:
    """Maximal-length sequence over GF(5) (period 24) from x² = x + 3."""
    seq = [1, 0]
    while len(seq) < 24:
        seq.append((seq[-1] * 1 + seq[-2] * 3) % 5)
    return seq


def brereton_design(factors: int) -> np.ndarray:
    """Brereton multilevel (5-level) multifactor calibration design.

    Returns a 25 × ``factors`` matrix of coded levels in {−2, −1, 0, 1, 2}
    (the first row is the centre point). Up to 6 factors are supported with
    orthogonal (uncorrelated) columns.
    """
    if not 1 <= factors <= 6:
        raise ValueError("the 25-sample design supports 1–6 compounds")
    seq = _gf5_msequence()
    level = {0: 0, 1: -2, 2: -1, 3: 1, 4: 2}
    coded = [level[v] for v in seq]
    rows = [[0] * factors]
    for i in range(24):
        rows.append([coded[(i + j) % 24] for j in range(factors)])
    return np.array(rows, float)


def design_concentrations(centres: dict[str, float], steps: dict[str, float]) -> list[dict[str, float]]:
    """Map the coded design to concentrations: centre + level × step."""
    names = list(centres)
    coded = brereton_design(len(names))
    return [{n: round(centres[n] + row[j] * steps[n], 6) for j, n in enumerate(names)}
            for row in coded]
