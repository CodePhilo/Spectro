"""Univariate calibration dialog and robustness study."""

from __future__ import annotations

import numpy as np
import pyqtgraph as pg
from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QDialog,
    QDoubleSpinBox,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QPushButton,
    QSplitter,
    QTableWidget,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)

from spectro.core import univariate as uv
from spectro.core.operations import apply_pipeline
from spectro.ui.dialogs_data import Base, name_of, spectrum_choices
from spectro.ui.methods.common import (
    CAL_ROLES,
    add_grouped,
    check_ids,
    checklists,
    list_box,
    results_rows,
    scatter,
    scroll,
    start_editing,
    store_method,
    summary_text,
    xy_plot,
)
from spectro.ui.widgets import (
    SERIES,
    MeasurementEditor,
    PasteTable,
    PipelineEditor,
    SpectrumPlot,
    error,
    fill_table,
)


# --------------------------------------------------------------------------- #
# Univariate
# --------------------------------------------------------------------------- #
class UnivariateDialog(Base):
    def __init__(self, win):
        super().__init__(win, "Univariate calibration")
        self.resize(1400, 860)
        self.method: uv.UnivariateMethod | None = None
        self.predictions: list = []
        comps = self.project.compound_names()
        root = QSplitter()
        left = QWidget()
        ll = QVBoxLayout(left)
        f = QFormLayout()
        self.compound = QComboBox()
        self.compound.addItems(comps)
        self.compound.setEditable(True)
        self.template = QComboBox()
        self.template.addItem("(choose a method template)")
        add_grouped(self.template, uv.templates_by_category())
        self.template.currentTextChanged.connect(self._template)
        self.name = QLineEdit()
        self.origin = QCheckBox("Force through origin")
        self.direct = QCheckBox("Concentration value: signal = concentration (no regression)")
        self.direct.setToolTip("For a plateau of a ratio spectrum made with a unit-concentration "
                               "(normalized) divisor. The calibration line is still reported.")
        self.spike = QDoubleSpinBox()
        self.spike.setRange(0, 1e6)
        self.spike.setDecimals(4)
        self.spike.setToolTip("Sample enrichment (spiking / spectrum addition): this amount is "
                              "subtracted from each found value (Found − added column).")
        f.addRow("Compound", self.compound)
        f.addRow("Method", self.template)
        f.addRow("Name", self.name)
        f.addRow("", self.origin)
        f.addRow("", self.direct)
        f.addRow("Enrichment added", self.spike)
        ll.addLayout(f)
        g = QGroupBox("1. Processing steps")
        gl = QVBoxLayout(g)
        self.pipe = PipelineEditor(spectrum_choices(self.project), name_of(self.project))
        self.pipe.changed.connect(self._preview)
        gl.addWidget(self.pipe)
        ll.addWidget(g)
        g = QGroupBox("2. Measurement (drag the λ lines on the plot)")
        gl = QVBoxLayout(g)
        self.meas = MeasurementEditor(spectrum_choices(self.project))
        self.meas.changed.connect(self._markers)
        gl.addWidget(self.meas)
        ll.addWidget(g)
        self.cal, self.test = checklists(self.project, comps[0] if comps else None)
        self.compound.currentTextChanged.connect(self._reselect)
        ll.addWidget(list_box("3. Calibration spectra (with known concentrations)", self.cal))
        ll.addWidget(list_box("4. Spectra to determine", self.test))
        self.cal.itemChanged.connect(lambda _: self._preview())
        btns = QHBoxLayout()
        for t, fn in (("Calibrate", self._calibrate), ("Determine", self._predict),
                      ("Save method", self._save_method), ("Save results", self._save_results)):
            b = QPushButton(t)
            b.clicked.connect(fn)
            btns.addWidget(b)
        ll.addLayout(btns)
        btns = QHBoxLayout()
        for t, fn, tip in (
                ("Standard addition…", self._std_addition,
                 "Spectra to determine = unspiked sample(s) + spiked samples; the compound "
                 "concentration of each spectrum is the amount ADDED (0 = unspiked)."),
                ("Robustness…", self._robustness,
                 "Vary method parameters (λ, Δλ, plateau, window…) by ± a small amount and "
                 "compare the results")):
            b = QPushButton(t)
            b.setToolTip(tip)
            b.clicked.connect(fn)
            btns.addWidget(b)
        ll.addLayout(btns)
        root.addWidget(scroll(left))
        self.tabs = QTabWidget()
        self.plot = SpectrumPlot(y_label="Processed signal")
        self.plot.line_moved.connect(lambda k, v: self.meas.form.set_value(k, round(v, 2)))
        self.tabs.addTab(self.plot, "Processed spectra")
        cw = QWidget()
        cl = QVBoxLayout(cw)
        self.cal_plot = xy_plot("Concentration", "Signal")
        self.res_plot = xy_plot("Concentration", "Residual")
        self.res_plot.setMaximumHeight(180)
        self.stats = QTableWidget()
        self.stats.setMaximumHeight(260)
        cl.addWidget(self.cal_plot, 2)
        cl.addWidget(self.res_plot, 1)
        cl.addWidget(self.stats)
        self.tabs.addTab(cw, "Calibration")
        rw = QWidget()
        rl = QVBoxLayout(rw)
        self.results = PasteTable()
        self.summary = QLabel()
        rl.addWidget(self.results)
        rl.addWidget(self.summary)
        self.tabs.addTab(rw, "Results")
        root.addWidget(self.tabs)
        root.setSizes([520, 880])
        QVBoxLayout(self).addWidget(root)
        self.msg = QLabel()
        self.layout().addWidget(self.msg)
        self._preview()

    def _reselect(self, comp):
        recs = self.project.records()
        ids = {r.id for r in recs if r.kind == "raw" and r.role in CAL_ROLES
               and r.concentrations.get(comp, 0) > 0}
        self.cal.populate(recs, ids)
        self._preview()

    def _template(self, name):
        t = uv.TEMPLATES.get(name)
        if not t:
            return
        self.pipe.set_steps(t["steps"])
        self.meas.set(t["measurement"])
        self.direct.setChecked(bool(t.get("direct")))
        self.name.setText(f"{name} — {self.compound.currentText()}")
        needs = [s for s in t["steps"] for k, v in s["params"].items() if v is None]
        if needs:
            self.msg.setText("Choose the divisor / reference spectra in the highlighted steps.")

    def _current(self) -> uv.UnivariateMethod:
        return uv.UnivariateMethod(self.name.text().strip() or "Univariate method",
                                   self.compound.currentText().strip(), self.pipe.get_steps(),
                                   self.meas.get(), self.origin.isChecked(),
                                   direct=self.direct.isChecked())

    def _preview(self):
        ids = self.cal.checked_ids()[:30] or self.test.checked_ids()[:30]
        try:
            spectra = [apply_pipeline(s, self.pipe.get_steps(), self.project.resolver())
                       for s in self.project.spectra(ids)]
            self.plot.plot_spectra(spectra, keep_range=False)
            self.msg.setText("")
        except Exception as exc:
            self.msg.setText(f"⚠ {exc}")
        self._markers()

    def _markers(self):
        wl = dict(self.meas.wavelengths())
        for k in ("w1", "w2"):
            self.plot.set_marker(k, wl.get(k))

    def _calibrate(self):
        m = self._current()
        ids = self.cal.checked_ids()
        if len(ids) < 3:
            error(self, "Check at least three calibration spectra.")
            return
        try:
            reg = m.calibrate(self.project.spectra(ids), self.project.resolver())
        except Exception as exc:
            error(self, exc)
            return
        self.method = m
        self.cal_ids = ids
        self.cal_plot.clear()
        self.res_plot.clear()
        names = [r.name for r in self.project.records(ids=ids)]
        xs = np.linspace(min(reg.x), max(reg.x), 2)
        scatter(self.cal_plot, reg.x, reg.y, SERIES[0], "Calibration", (xs, reg.predict_y(xs)),
                labels=names)
        scatter(self.res_plot, reg.x, reg.residuals, SERIES[0], labels=names)
        self.res_plot.addLine(y=0, pen=pg.mkPen("#8a8984", width=1))
        rows = [["Slope", reg.slope], ["SD of slope", reg.sd_slope],
                ["95% CI of slope", f"{reg.ci_slope[0]:.5g} – {reg.ci_slope[1]:.5g}"],
                ["Intercept", reg.intercept], ["SD of intercept", reg.sd_intercept],
                ["95% CI of intercept", f"{reg.ci_intercept[0]:.5g} – {reg.ci_intercept[1]:.5g}"],
                ["Correlation coefficient (r)", reg.r], ["r²", reg.r2],
                ["Sy/x (residual SD)", reg.sy_x], ["n", str(reg.n)],
                ["Range", f"{reg.range[0]:g} – {reg.range[1]:g}"],
                ["LOD (3.3σ/S)", reg.lod], ["LOQ (10σ/S)", reg.loq]]
        fill_table(self.stats, ["Parameter", "Value"], rows)
        self.tabs.setCurrentIndex(1)
        self.project.log("CALCULATE", f"Calibrated '{m.name}' for {m.compound}: "
                                      f"r = {reg.r:.5f}", {"method": m.to_dict(), "inputs": ids})

    def _predict(self):
        if self.method is None:
            self._calibrate()
            if self.method is None:
                return
        ids = self.test.checked_ids()
        if not ids:
            error(self, "Check the spectra to determine.")
            return
        comp = self.method.compound
        found, names, taken, signals = [], [], [], []
        try:
            for s in self.project.spectra(ids):
                sig = self.method.signal(s, self.project.resolver())
                signals.append(sig)
                found.append(self.method.concentration(sig))
                names.append(s.name)
                taken.append(s.concentrations.get(comp))
        except Exception as exc:
            error(self, exc)
            return
        rows, recs = results_rows(names, found, taken)
        reg = None if self.method.direct else self.method.regression
        for r, sig in zip(rows, signals):
            r.insert(1, float(sig))
            if reg is not None:
                r.insert(3, reg.ci_x(float(sig)))
        headers = ["Spectrum", "Signal", f"Found ({comp})"] + \
            (["± 95 % CI"] if reg is not None else []) + ["Taken", "Recovery %"]
        spike = self.spike.value()
        if spike:
            headers.append(f"Found − added ({spike:g})")
            for r, f in zip(rows, found):
                r.append(uv.enrichment_correction(f, spike)["found"])
        fill_table(self.results, headers, rows)
        self.summary.setText(summary_text(recs))
        self.predictions = {"ids": ids, "names": names, "found": found, "taken": taken,
                            "signals": signals, "recovery": recs}
        if spike:
            self.predictions["enrichment_added"] = spike
            self.predictions["found_minus_added"] = [f - spike for f in found]
        self.last = (self.method.name, "univariate",
                     {"method": self.method.to_dict(), "compounds": [comp],
                      "found": [[f] for f in found], **{k: v for k, v in self.predictions.items()
                                                         if k != "found"}}, ids)
        self.tabs.setCurrentIndex(2)
        self.project.log("CALCULATE", f"Determined {comp} in {len(ids)} spectra with "
                                      f"'{self.method.name}'",
                         {"inputs": ids, "method": self.method.to_dict(), "found": found})

    def _save_method(self):
        if self.method is None:
            error(self, "Calibrate first.")
            return
        d = self.method.to_dict()
        d["calibration_ids"] = self.cal_ids
        store_method(self, self.method.name, d)

    def load_method(self, md: dict) -> None:
        """Fill the dialog from a saved method so it can be changed and re-saved."""
        d = md["definition"]
        self.compound.setCurrentText(d["compound"])
        self.name.setText(md["name"])
        self.origin.setChecked(bool(d.get("through_origin")))
        self.direct.setChecked(bool(d.get("direct")))
        self.pipe.set_steps(d.get("steps", []))
        self.meas.set(d["measurement"])
        if d.get("calibration_ids"):
            check_ids(self.cal, d["calibration_ids"])
        start_editing(self, md)
        self._preview()

    def _save_results(self):
        if not getattr(self, "last", None):
            error(self, "Determine, run standard addition or a robustness study first.")
            return
        name, kind, data, ids = self.last
        if getattr(self, "cal_ids", None) and kind == "univariate":
            data = {**data, "calibration_ids": list(self.cal_ids)}
        self.project.save_result(name, kind, data, self.win.current_trial(),
                                 getattr(self, "method_id", None), ids)
        self.win.statusBar().showMessage("Results saved.")

    def _ready_method(self) -> bool:
        if self.method is None:
            self._calibrate()
        return self.method is not None

    def _std_addition(self):
        if not self._ready_method():
            return
        ids = self.test.checked_ids()
        comp = self.method.compound
        spectra = self.project.spectra(ids)
        added = [s.concentrations.get(comp, 0.0) for s in spectra]
        if len(ids) < 2 or not any(added) or all(added):
            error(self, "Check the unspiked sample(s) (concentration of "
                        f"{comp} = 0 or blank) and the spiked samples (concentration of {comp} "
                        "= amount added).")
            return
        res = self.project.resolver()
        try:
            sa = uv.standard_addition_recovery(lambda s: self.method.predict(s, res),
                                               spectra, added)
        except Exception as exc:
            error(self, exc)
            return
        rows = [[r["name"], r["added"], r["found"],
                 "" if r["recovered"] is None else r["recovered"],
                 "" if r["recovery"] is None else r["recovery"]] for r in sa["rows"]]
        fill_table(self.results, ["Spectrum", f"{comp} added", "Found (total)",
                                  "Recovered (found − sample)", "Recovery %"], rows)
        self.results.export_title = f"Standard addition {comp}"
        lines = [f"{comp} in the sample (mean of unspiked): {sa['sample']:.4f}",
                 f"Recovery of added {comp}: mean {sa['mean_recovery']:.2f} %"
                 + (f", SD {sa['sd_recovery']:.3f}" if sa["sd_recovery"] == sa["sd_recovery"]
                    else "")]
        if "extrapolated" in sa:
            lines.append(f"Extrapolated (x-intercept of found vs added): {sa['extrapolated']:.4f}")
        self.summary.setText("\n".join(lines))
        self.results.export_notes = lines
        self.tabs.setCurrentIndex(2)
        self.last = (f"{self.method.name} — standard addition", "standard_addition",
                     {"method": self.method.to_dict(), "ids": ids, **sa}, ids)
        self.project.log("CALCULATE", f"Standard addition for {comp} with '{self.method.name}'",
                         {"inputs": ids, "result": sa})

    def _robustness(self):
        if not self._ready_method():
            return
        test_ids = self.test.checked_ids()
        if not test_ids:
            error(self, "Check the spectra to determine (e.g. lab mixtures) first.")
            return
        dlg = RobustnessDialog(self, self.method, self.project.spectra(self.cal_ids),
                               self.project.spectra(test_ids))
        if dlg.exec() and dlg.result_data:
            self.last = (f"{self.method.name} — robustness", "robustness",
                         {"method": self.method.to_dict(), "ids": test_ids,
                          **dlg.result_data}, test_ids)
            self.project.log("CALCULATE", f"Robustness study of '{self.method.name}'",
                             {"inputs": self.cal_ids + test_ids,
                              "rsd_across_variants": dlg.result_data["rsd_across_variants"],
                              "max_abs_deviation": dlg.result_data["max_abs_deviation"]})


class RobustnessDialog(QDialog):
    """Choose parameters and ± deltas; run the study in a worker thread."""

    def __init__(self, parent, method: uv.UnivariateMethod, cal, test):
        super().__init__(parent)
        self.setWindowTitle(f"Robustness — {method.name}")
        self.resize(1000, 640)
        self.method, self.cal, self.test = method, cal, test
        self.project = parent.project
        self.result_data = None
        self.params = uv.numeric_parameters(method)
        lay = QVBoxLayout(self)
        lay.addWidget(QLabel("Each checked parameter is moved to nominal − Δ and nominal + Δ; "
                             "the method is re-calibrated and the spectra are re-assayed. "
                             "Small deviations (e.g. < 2 %) indicate a robust method."))
        lay.itemAt(0).widget().setWordWrap(True)
        self.table = PasteTable()
        fill_table(self.table, ["Use", "Parameter", "Nominal", "Δ (±)"],
                   [["", p["label"], p["value"], f"{p['delta']:g}"] for p in self.params],
                   editable=True)
        for i, p in enumerate(self.params):
            it = self.table.item(i, 0)
            it.setFlags(Qt.ItemIsUserCheckable | Qt.ItemIsEnabled)
            it.setCheckState(Qt.Checked if p["kind"] == "wavelength" else Qt.Unchecked)
            for j in (1, 2):
                self.table.item(i, j).setFlags(Qt.ItemIsEnabled | Qt.ItemIsSelectable)
        self.table.setMaximumHeight(220)
        lay.addWidget(self.table)
        run = QPushButton("Run robustness study")
        run.clicked.connect(self._run)
        lay.addWidget(run)
        self.out = PasteTable()
        self.out.export_title = "Robustness"
        lay.addWidget(self.out, 1)
        self.summary = QLabel()
        self.summary.setWordWrap(True)
        lay.addWidget(self.summary)
        from PySide6.QtWidgets import QDialogButtonBox
        bb = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        bb.button(QDialogButtonBox.Ok).setText("Keep result")
        bb.accepted.connect(self.accept)
        bb.rejected.connect(self.reject)
        lay.addWidget(bb)

    def _run(self):
        from spectro.ui.widgets import run_with_progress, snapshot_resolver

        variations = []
        try:
            for i, p in enumerate(self.params):
                if self.table.item(i, 0).checkState() == Qt.Checked:
                    d = float(self.table.item(i, 3).text().replace(",", "."))
                    variations.append({**p, "delta": d})
        except ValueError as exc:
            error(self, f"Invalid Δ: {exc}")
            return
        if not variations:
            error(self, "Check at least one parameter.")
            return
        resolve = snapshot_resolver(self.project, self.method.steps, self.method.measurement)
        try:
            res = run_with_progress(self, "Running robustness study…", uv.robustness_study,
                                    method=uv.UnivariateMethod.from_dict(self.method.to_dict()),
                                    calibration=self.cal, test=self.test,
                                    variations=variations, resolve=resolve)
        except InterruptedError:
            return
        except Exception as exc:
            error(self, exc)
            return
        rows = []
        for v in res["variants"]:
            if "error" in v:
                rows.append([v["label"], v["change"], "", "", "", "", "", v["error"]])
            else:
                rows.append([v["label"], v["change"], v["slope"], v["r"], v["mean_found"],
                             v["mean_recovery"], v["rsd_recovery"], v["deviation"]])
        fill_table(self.out, ["Parameter", "Change", "Slope", "r", "Mean found",
                              "Mean recovery %", "RSD %", "Deviation from nominal %"], rows)
        verdict = "robust" if res["max_abs_deviation"] < 2 else "sensitive — review"
        self.summary.setText(f"RSD of results across variants: {res['rsd_across_variants']:.3f} % · "
                             f"largest deviation from nominal: {res['max_abs_deviation']:.3f} % "
                             f"→ {verdict}")
        self.out.export_notes = [self.summary.text()]
        self.result_data = res
