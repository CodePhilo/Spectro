"""Chemometric models (CLS, ILS, PCR, PLS, MCR-ALS, ANN, SVR)."""

from __future__ import annotations

import numpy as np
import pyqtgraph as pg
from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QComboBox,
    QDoubleSpinBox,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QInputDialog,
    QLabel,
    QLineEdit,
    QPlainTextEdit,
    QPushButton,
    QSpinBox,
    QSplitter,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)

from spectro.core.multicomponent import MODEL_TYPES, SpectralModel, ga_select, ipls
from spectro.ui.dialogs_data import Base, name_of, spectrum_choices
from spectro.ui.methods.common import (
    TRAIN_ROLES,
    CompoundChecks,
    check_ids,
    checklists,
    list_box,
    scatter,
    scroll,
    start_editing,
    store_method,
    summary_text,
    xy_plot,
)
from spectro.ui.widgets import SERIES, PasteTable, PipelineEditor, SpectrumPlot, error, fill_table


# --------------------------------------------------------------------------- #
# Chemometrics
# --------------------------------------------------------------------------- #
class ChemometricsDialog(Base):
    def __init__(self, win):
        super().__init__(win, "Chemometrics")
        self.resize(1450, 900)
        self.model: SpectralModel | None = None
        self.predictions = None
        comps = self.project.compound_names()
        root = QSplitter()
        left = QWidget()
        ll = QVBoxLayout(left)
        f = QFormLayout()
        self.mtype = QComboBox()
        self.mtype.addItems(MODEL_TYPES)
        self.mtype.setCurrentText("PLS2")
        self.mtype.currentTextChanged.connect(self._type_changed)
        self.comps = CompoundChecks(comps)
        self.prep = QComboBox()
        self.prep.addItems(["mean_center", "autoscale", "none"])
        self.ncomp = QSpinBox()
        self.ncomp.setRange(1, 30)
        self.ncomp.setValue(max(1, len(comps)))
        self.ranges = QLineEdit()
        self.ranges.setPlaceholderText("e.g. 210-370; 380-395 (blank = all)")
        self.wls = QLineEdit()
        self.wls.setPlaceholderText("ILS: e.g. 225, 245, 260, 275")
        self.hidden = QLineEdit("10")
        self.act = QComboBox()
        self.act.addItems(["tanh", "logistic", "relu", "identity"])
        self.kernel = QComboBox()
        self.kernel.addItems(["linear", "rbf", "poly"])
        self.svr_c = QDoubleSpinBox()
        self.svr_c.setRange(1e-4, 1e6)
        self.svr_c.setValue(100)
        self.svr_eps = QDoubleSpinBox()
        self.svr_eps.setDecimals(4)
        self.svr_eps.setRange(0, 100)
        self.svr_eps.setValue(0.01)
        self.cv = QComboBox()
        self.cv.addItems(["loo", "venetian", "kfold", "contiguous"])
        self.folds = QSpinBox()
        self.folds.setRange(2, 50)
        self.folds.setValue(5)
        for label, w in (("Model", self.mtype), ("Compounds", self.comps),
                         ("Preprocessing", self.prep), ("Components / LVs", self.ncomp),
                         ("Wavelength ranges", self.ranges), ("ILS wavelengths", self.wls),
                         ("ANN hidden layers", self.hidden), ("ANN activation", self.act),
                         ("SVR kernel", self.kernel), ("SVR C", self.svr_c),
                         ("SVR ε", self.svr_eps), ("Cross-validation", self.cv),
                         ("Folds", self.folds)):
            f.addRow(label, w)
        ll.addLayout(f)
        g = QGroupBox("Optional processing before modelling")
        gl = QVBoxLayout(g)
        self.pipe = PipelineEditor(spectrum_choices(self.project), name_of(self.project))
        gl.addWidget(self.pipe)
        ll.addWidget(g)
        self.cal, self.test = checklists(self.project, cal_roles=TRAIN_ROLES)
        ll.addWidget(list_box("Calibration (training) set", self.cal))
        ll.addWidget(list_box("Validation / samples to predict", self.test))
        btns = QHBoxLayout()
        for t, fn in (("Cross-validate", self._cv), ("Fit + predict", self._fit),
                      ("iPLS", self._ipls), ("GA selection", self._ga)):
            b = QPushButton(t)
            b.clicked.connect(fn)
            btns.addWidget(b)
        ll.addLayout(btns)
        btns = QHBoxLayout()
        for t, fn in (("Save method", self._save_method), ("Save results", self._save_results)):
            b = QPushButton(t)
            b.clicked.connect(fn)
            btns.addWidget(b)
        ll.addLayout(btns)
        root.addWidget(scroll(left))
        self.tabs = QTabWidget()
        self.cv_plot = xy_plot("Number of components", "RMSECV")
        self.tabs.addTab(self.cv_plot, "RMSECV")
        self.pred_plot = xy_plot("Actual concentration", "Predicted concentration")
        self.tabs.addTab(self.pred_plot, "Predicted vs actual")
        rw = QWidget()
        rl = QVBoxLayout(rw)
        self.results = PasteTable()
        self.summary = QLabel()
        self.summary.setWordWrap(True)
        rl.addWidget(self.results)
        rl.addWidget(self.summary)
        self.tabs.addTab(rw, "Results")
        ow = QWidget()
        ol = QVBoxLayout(ow)
        orow = QHBoxLayout()
        self.conf = QComboBox()
        self.conf.addItems(["95", "99"])
        self.conf.currentTextChanged.connect(lambda _: self._show_diagnostics())
        excl = QPushButton("Exclude flagged calibration spectra and refit")
        excl.clicked.connect(self._exclude_outliers)
        orow.addWidget(QLabel("Confidence limit (%)"))
        orow.addWidget(self.conf)
        orow.addStretch(1)
        orow.addWidget(excl)
        ol.addLayout(orow)
        self.diag_plot = xy_plot("Hotelling T²", "Q residuals")
        ol.addWidget(self.diag_plot, 2)
        self.diag_table = PasteTable()
        self.diag_table.export_title = "Outlier diagnostics"
        ol.addWidget(self.diag_table, 1)
        self.tabs.addTab(ow, "Outliers (T² / Q)")
        self.aux_plot = SpectrumPlot(y_label="Value")
        self.tabs.addTab(self.aux_plot, "Loadings / VIP / pure spectra")
        self.info = QPlainTextEdit()
        self.info.setReadOnly(True)
        self.tabs.addTab(self.info, "Model info")
        root.addWidget(self.tabs)
        root.setSizes([520, 930])
        QVBoxLayout(self).addWidget(root)
        self._type_changed(self.mtype.currentText())

    def _type_changed(self, t):
        self.wls.setEnabled(t == "ILS")
        for w in (self.hidden, self.act):
            w.setEnabled(t == "ANN")
        for w in (self.kernel, self.svr_c, self.svr_eps):
            w.setEnabled(t == "SVR")
        self.ncomp.setEnabled(t in ("PCR", "PLS1", "PLS2", "ANN"))

    def _build(self) -> SpectralModel:
        ranges = []
        for part in self.ranges.text().replace(",", ";").split(";"):
            part = part.strip()
            if part:
                a, b = part.replace("–", "-").split("-")
                ranges.append((float(a), float(b)))
        wls = [float(w) for w in self.wls.text().replace(";", ",").split(",") if w.strip()]
        opts = {"hidden": self.hidden.text(), "activation": self.act.currentText(),
                "kernel": self.kernel.currentText(), "C": self.svr_c.value(),
                "epsilon": self.svr_eps.value()}
        comps = self.comps.checked()
        if not comps:
            raise ValueError("check at least one compound")
        return SpectralModel(self.mtype.currentText(), comps, self.pipe.get_steps(), ranges, wls,
                             self.ncomp.value(), self.prep.currentText(), opts)

    def _cal(self):
        ids = self.cal.checked_ids()
        if len(ids) < 3:
            raise ValueError("check at least three calibration spectra")
        return ids, self.project.spectra(ids)

    def _work(self, label, fn, **kwargs):
        """Run fn in a worker thread with a cancellable progress dialog."""
        from spectro.ui.widgets import run_with_progress
        return run_with_progress(self, label, fn, **kwargs)

    def _cv(self):
        from spectro.ui.widgets import snapshot_resolver
        try:
            m = self._build()
            ids, spectra = self._cal()
            latent = m.model_type in ("PCR", "PLS1", "PLS2", "ANN")
            cv = self._work(f"Cross-validating {m.model_type}…", m.cross_validate,
                            spectra=spectra, resolve=snapshot_resolver(self.project, m.steps),
                            method=self.cv.currentText(), folds=self.folds.value(),
                            max_components=min(15, len(ids) - 2) if latent else None)
        except InterruptedError:
            return
        except Exception as exc:
            error(self, exc)
            return
        self.cv_plot.clear()
        ks = [c["k"] for c in cv["curve"]]
        for j, comp in enumerate(m.compounds):
            self.cv_plot.plot(ks, [c["RMSECV"][j] for c in cv["curve"]],
                              pen=pg.mkPen(SERIES[j % 8], width=2), symbol="o", symbolSize=8,
                              symbolBrush=SERIES[j % 8], name=comp)
        if len(ks) > 1:
            self.ncomp.setValue(cv["suggested_components"])
        stats = "\n".join(f"{c}: RMSECV {s['RMSEP']:.4g}, R² {s['R2']:.4f}, bias {s['bias']:.3g}"
                          for c, s in cv["stats"].items())
        self.info.setPlainText(f"Cross-validation ({self.cv.currentText()}):\n{stats}\n\n"
                               f"Suggested number of components (Haaland–Thomas, α = 0.25): "
                               f"{cv['suggested_components']}")
        self.tabs.setCurrentIndex(0)
        self.project.log("CALCULATE", f"Cross-validated {m.model_type}",
                         {"inputs": ids, "model": m.to_dict(),
                          "curve": cv["curve"], "suggested": cv["suggested_components"]})

    def _fit(self):
        from spectro.ui.widgets import snapshot_resolver
        try:
            m = self._build()
            ids, spectra = self._cal()
            test_ids = self.test.checked_ids()
            test = self.project.spectra(test_ids)
            resolve = snapshot_resolver(self.project, m.steps)

            def job(progress):
                progress(0, 3)
                m.fit(spectra, resolve)
                progress(1, 3)
                pred = m.predict(test, resolve) if test else np.zeros((0, len(m.compounds)))
                progress(2, 3)
                d_cal = m.diagnostics(spectra, resolve)
                d_test = m.diagnostics(test, resolve) if test else None
                return pred, d_cal, d_test
            pred, diag, diag_test = self._work(f"Fitting {m.model_type}…", job)
        except InterruptedError:
            return
        except Exception as exc:
            error(self, exc)
            return
        self.model, self.cal_ids = m, ids
        self._diag = (spectra, test, diag, diag_test)
        # predicted vs actual
        self.pred_plot.clear()
        rows, lines = [], []
        lo, hi = np.inf, -np.inf
        from spectro.core.validation import prediction_error
        for j, c in enumerate(m.compounds):
            pts = [(s.concentrations[c], p[j], s.name) for s, p in zip(test, pred)
                   if s.concentrations.get(c) is not None]
            if pts:
                a, b, n = zip(*pts)
                scatter(self.pred_plot, a, b, SERIES[j % 8], c, labels=list(n))
                lo, hi = min(lo, *a, *b), max(hi, *a, *b)
                recs = [100 * y / x for x, y, _ in pts if x]
                pe = prediction_error(b, a)
                lines.append(f"{c}: RMSEP {pe['RMSEP']:.4g}, R² {pe['R2']:.4f}. "
                             + summary_text(recs))
        if np.isfinite(lo):
            self.pred_plot.plot([lo, hi], [lo, hi], pen=pg.mkPen("#8a8984", width=1,
                                                                 style=Qt.DashLine))
        for s, p in zip(test, pred):
            row = [s.name]
            for j, c in enumerate(m.compounds):
                t = s.concentrations.get(c)
                row += [float(p[j]), "" if t is None else t, "" if not t else 100 * p[j] / t]
            rows.append(row)
        fill_table(self.results, ["Spectrum"] + [h for c in m.compounds for h in
                                                 (f"{c} found", f"{c} taken", f"{c} rec %")], rows)
        self.results.export_title = f"{m.model_type} predictions"
        self.summary.setText("\n".join(lines))
        self.results.export_notes = lines
        self._show_diagnostics()
        # auxiliary spectra
        from spectro.core.spectrum import Spectrum
        grid = np.array(diag["grid"])
        if m.model_type == "MCR-ALS":
            aux = [Spectrum(grid, row, name=f"MCR-ALS {c}") for c, row in
                   zip(m.compounds, m.state["S"])]
        elif m.model_type == "CLS":
            aux = [Spectrum(grid, row, name=f"K {c}") for c, row in zip(m.compounds, m.state["K"])]
        elif m.vip() is not None:
            aux = [Spectrum(grid, m.vip(), name="VIP")]
        else:
            L = np.array(diag["loadings"]).T
            aux = [Spectrum(grid, L[i], name=f"PC{i + 1} loading") for i in range(min(4, len(L)))]
        self.aux_plot.plot_spectra(aux)
        info = [f"Model: {m.model_type}; compounds: {', '.join(m.compounds)}; "
                f"{len(ids)} calibration spectra; {grid.size} variables "
                f"({grid[0]:g}–{grid[-1]:g} nm); preprocessing: {m.preprocessing}"]
        if m.model_type in ("PCR", "PLS1", "PLS2", "ANN"):
            info.append(f"Components: {m.n_components}")
        info.append("RMSEC: " + ", ".join(f"{c} {v:.4g}" for c, v in
                                          zip(m.compounds, m.state["rmsec"])))
        info.append("Explained variance (PCA of X): " + ", ".join(
            f"{100 * v:.2f}%" for v in diag["explained"][:6]))
        lim = diag["limits"]
        info.append(f"Limits: T² 95 % {lim['t2_lim95']:.3g}, 99 % {lim['t2_lim99']:.3g}; "
                    f"Q 95 % {lim['q_lim95']:.3g}, 99 % {lim['q_lim99']:.3g}")
        if m.model_type == "MCR-ALS":
            info.append(f"MCR-ALS lack of fit: {m.state['lof']:.4f} %")
        self.info.setPlainText("\n".join(info))
        self.predictions = {"ids": test_ids, "compounds": m.compounds, "found": pred.tolist()}
        self.tabs.setCurrentIndex(1 if rows else 3)
        self.project.log("CALCULATE", f"Fitted {m.model_type} and predicted {len(test_ids)} spectra",
                         {"inputs": ids + test_ids, "model": m.to_dict(),
                          "found": pred.tolist()})

    def _show_diagnostics(self):
        if not getattr(self, "_diag", None):
            return
        spectra, test, _, _ = self._diag
        a = self.conf.currentText()
        from spectro.ui.widgets import snapshot_resolver
        resolve = snapshot_resolver(self.project, self.model.steps)
        diag = self.model.diagnostics(spectra, resolve, alpha=a)
        diag_test = self.model.diagnostics(test, resolve, alpha=a) if test else None
        self.diag_plot.clear()
        scatter(self.diag_plot, diag["T2"], diag["Q"], SERIES[0], "Calibration",
                labels=[s.name for s in spectra])
        if diag_test:
            scatter(self.diag_plot, diag_test["T2"], diag_test["Q"], SERIES[1], "Predicted",
                    labels=[s.name for s in test])
        lim = diag["limits"]
        for tag, style in (("95", Qt.DashLine), ("99", Qt.DotLine)):
            pen = pg.mkPen("#8a8984", width=1, style=style)
            self.diag_plot.addItem(pg.InfiniteLine(lim[f"t2_lim{tag}"], angle=90, pen=pen,
                                                   label=f"T² {tag}%",
                                                   labelOpts={"position": 0.95}))
            self.diag_plot.addItem(pg.InfiniteLine(lim[f"q_lim{tag}"], angle=0, pen=pen,
                                                   label=f"Q {tag}%",
                                                   labelOpts={"position": 0.95}))
        rows = []
        self._flagged = []
        for group, sp, d in (("calibration", spectra, diag), ("predicted", test, diag_test)):
            if not d:
                continue
            for s, t2, q, f in zip(sp, d["T2"], d["Q"], d["flags"]):
                rows.append([s.name, group, t2, q, ", ".join(f) or "—"])
                if f and group == "calibration":
                    self._flagged.append(s.metadata.get("id"))
        fill_table(self.diag_table, ["Spectrum", "Set", "Hotelling T²", "Q residual",
                                     f"Flags ({a} % limits, |residual| > 3·RMSEC)"], rows)

    def _exclude_outliers(self):
        flagged = [i for i in getattr(self, "_flagged", []) if i is not None]
        if not flagged:
            error(self, "No calibration spectra are flagged at the chosen limit.")
            return
        names = ", ".join(self.project.record(i).name for i in flagged)
        from PySide6.QtWidgets import QMessageBox
        if QMessageBox.question(self, "Exclude outliers",
                                f"Uncheck {len(flagged)} calibration spectra and refit?\n\n"
                                f"{names}") != QMessageBox.Yes:
            return
        for i in range(self.cal.count()):
            if self.cal.item(i).data(Qt.UserRole) in flagged:
                self.cal.item(i).setCheckState(Qt.Unchecked)
        self.project.log("CALCULATE", f"Excluded {len(flagged)} outlying calibration spectra "
                                      f"({self.model.model_type})",
                         {"inputs": flagged, "limit": self.conf.currentText() + " %"})
        self._fit()

    def _ipls(self):
        from spectro.ui.widgets import snapshot_resolver
        try:
            m = self._build()
            ids, spectra = self._cal()
            n, ok = QInputDialog.getInt(self, "iPLS", "Number of intervals:", 10, 2, 50)
            if not ok:
                return
            out = self._work("Interval PLS…", ipls, model=m, spectra=spectra,
                             resolve=snapshot_resolver(self.project, m.steps), intervals=n)
        except InterruptedError:
            return
        except Exception as exc:
            error(self, exc)
            return
        fill_table(self.results, ["Interval (nm)", "RMSECV"],
                   [[f"{o['range'][0]:g}–{o['range'][1]:g}", o["RMSECV"]] for o in out])
        best = min(out, key=lambda o: o["RMSECV"] if np.isfinite(o["RMSECV"]) else np.inf)
        self.summary.setText(f"Best interval: {best['range'][0]:g}–{best['range'][1]:g} nm. "
                             "Enter it in 'Wavelength ranges' to use it.")
        self.tabs.setCurrentIndex(2)
        self.project.log("CALCULATE", "iPLS interval selection", {"inputs": ids, "result": out})

    def _ga(self):
        from spectro.ui.widgets import snapshot_resolver
        try:
            m = self._build()
            ids, spectra = self._cal()
            out = self._work("Genetic algorithm variable selection…", ga_select, model=m,
                             spectra=spectra, resolve=snapshot_resolver(self.project, m.steps))
        except InterruptedError:
            return
        except Exception as exc:
            error(self, exc)
            return
        text = "; ".join(f"{a:g}-{b:g}" for a, b in out["ranges"])
        self.ranges.setText(text)
        self.summary.setText(f"GA selected {len(out['ranges'])} intervals (RMSECV "
                             f"{out['RMSECV']:.4g}); ranges were filled in.")
        self.project.log("CALCULATE", "GA variable selection", {"inputs": ids, "result": out})

    def _save_method(self):
        if self.model is None:
            error(self, "Fit first.")
            return
        ed = getattr(self, "editing", None)
        name, ok = QInputDialog.getText(self, "Save method", "Name:",
                                        text=ed["name"] if ed else
                                        f"{self.model.model_type} model")
        if ok:
            d = self.model.to_dict(include_fit=True)
            d["calibration_ids"] = self.cal_ids
            if store_method(self, name.strip() or "Model", d) is not None:
                self.win.statusBar().showMessage("Fitted model saved (applies without "
                                                 "refitting).")

    def load_method(self, md: dict) -> None:
        d = md["definition"]
        self.mtype.setCurrentText(d["model_type"])
        for i in range(self.comps.count()):
            it = self.comps.item(i)
            it.setCheckState(Qt.Checked if it.text() in d["compounds"] else Qt.Unchecked)
        self.pipe.set_steps(d.get("steps", []))
        self.ranges.setText("; ".join(f"{a:g}-{b:g}" for a, b in d.get("ranges", [])))
        self.wls.setText(", ".join(f"{w:g}" for w in d.get("wavelengths", [])))
        self.ncomp.setValue(int(d.get("n_components", 2)))
        self.prep.setCurrentText(d.get("preprocessing", "mean_center"))
        o = d.get("options", {})
        if "hidden" in o:
            self.hidden.setText(str(o["hidden"]))
        if "activation" in o:
            self.act.setCurrentText(o["activation"])
        if "kernel" in o:
            self.kernel.setCurrentText(o["kernel"])
        if "C" in o:
            self.svr_c.setValue(float(o["C"]))
        if "epsilon" in o:
            self.svr_eps.setValue(float(o["epsilon"]))
        if d.get("calibration_ids"):
            check_ids(self.cal, d["calibration_ids"])
        start_editing(self, md)

    def _save_results(self):
        if not self.predictions:
            error(self, "Fit + predict first.")
            return
        self.project.save_result(f"{self.model.model_type} results", "chemometrics",
                                 {"model": self.model.to_dict(), **self.predictions,
                                  "calibration_ids": list(getattr(self, "cal_ids", []))},
                                 self.win.current_trial(), getattr(self, "method_id", None),
                                 self.predictions["ids"])
