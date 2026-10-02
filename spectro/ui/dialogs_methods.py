"""Quantitative method dialogs."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pyqtgraph as pg
from PySide6.QtCore import Qt
from PySide6.QtWidgets import (QCheckBox, QComboBox, QDialog, QDialogButtonBox,
                               QDoubleSpinBox, QFileDialog,
                               QFormLayout,
                               QGroupBox,
                               QHBoxLayout, QInputDialog, QLabel, QLineEdit, QListWidget,
                               QListWidgetItem, QPlainTextEdit, QPushButton, QScrollArea,
                               QSpinBox, QSplitter, QTableWidget, QTableWidgetItem, QTabWidget,
                               QVBoxLayout, QWidget)

from spectro.core import univariate as uv
from spectro.core.multicomponent import (MODEL_TYPES, SignalEquations, SpectralModel, ga_select,
                                         ipls, kaiser_selection)
from spectro.core.operations import apply_pipeline
from spectro.core.validation import describe
from spectro.ui.dialogs_data import Base, name_of, spectrum_choices
from spectro.ui.widgets import (SERIES, TEXT_SECONDARY, MeasurementEditor, PasteTable,
                                PipelineEditor, SpectrumChecklist, SpectrumPlot, error,
                                fill_table)

CAL_ROLES = {"standard"}                    # univariate / equation methods
TRAIN_ROLES = {"calibration", "standard"}   # chemometric training set
TEST_ROLES = {"mixture", "sample", "validation", "unknown"}


def scroll(widget: QWidget) -> QScrollArea:
    sa = QScrollArea()
    sa.setWidgetResizable(True)
    sa.setWidget(widget)
    return sa


def xy_plot(x_label: str, y_label: str) -> pg.PlotWidget:
    p = pg.PlotWidget()
    p.setLabel("bottom", x_label)
    p.setLabel("left", y_label)
    p.showGrid(x=True, y=True, alpha=0.12)
    p.addLegend(offset=(10, 10), labelTextColor=TEXT_SECONDARY)
    from spectro.ui.figure import install_figure_export
    install_figure_export(p)
    return p


def scatter(p: pg.PlotWidget, x, y, colour: str, name: str | None = None, line=None,
            labels: list[str] | None = None) -> None:
    """Markers ≥ 8 px with a surface-coloured ring; hover shows the label."""
    kw = dict(x=list(x), y=list(y), size=9, pen=pg.mkPen(SURFACE_RING, width=2),
              brush=pg.mkBrush(colour), name=name, hoverable=True)
    if labels:
        kw.update(data=list(labels), tip=lambda x, y, data: f"{data}\n({x:.4g}, {y:.4g})")
    p.addItem(pg.ScatterPlotItem(**kw))
    if line is not None:
        lx, ly = line
        p.plot(lx, ly, pen=pg.mkPen(colour, width=2))


SURFACE_RING = "#fcfcfb"


def checklists(project, compound: str | None = None, cal_roles=CAL_ROLES, test_roles=TEST_ROLES,
               trial: int | None = None):
    recs = project.records(trial_id=trial)
    cal, test = SpectrumChecklist(), SpectrumChecklist()
    # defaults use raw spectra only: methods apply their own processing, and
    # derived spectra inherit the role and concentrations of their parent
    raw = [r for r in recs if r.kind == "raw"]
    cal_ids = {r.id for r in raw if r.role in cal_roles and
               (compound is None or r.concentrations.get(compound, 0) > 0)}
    if cal_roles is TRAIN_ROLES and any(r.role == "calibration" for r in raw):
        cal_ids = {r.id for r in raw if r.role == "calibration"}
    test_ids = {r.id for r in raw if r.role in test_roles}
    cal.populate(recs, cal_ids)
    test.populate(recs, test_ids)
    return cal, test


def add_grouped(combo: QComboBox, groups: dict[str, list[str]]) -> None:
    """Items under bold, non-selectable family headings."""
    from PySide6.QtGui import QFont

    from spectro.core.catalog import CATEGORIES
    model = combo.model()
    for key, names in groups.items():
        combo.addItem(f"— {CATEGORIES.get(key, key)} —")
        head = model.item(combo.count() - 1)
        head.setEnabled(False)
        f = QFont(head.font())
        f.setBold(True)
        head.setFont(f)
        for n in names:
            combo.addItem(n)


def check_ids(lst: SpectrumChecklist, ids) -> None:
    """Check exactly the spectra in ``ids`` (unknown ids are ignored)."""
    ids = {int(i) for i in ids or []}
    for i in range(lst.count()):
        it = lst.item(i)
        it.setCheckState(Qt.Checked if it.data(Qt.UserRole) in ids else Qt.Unchecked)


def store_method(dlg, name: str, definition: dict) -> int | None:
    """Save a new method, or — when the dialog was opened with *Edit…* from
    Saved methods — a new version of the edited one (old version archived,
    reason required)."""
    from spectro.ui.widgets import ask_reason
    ed = getattr(dlg, "editing", None)
    if ed:
        reason = ask_reason(dlg, f"Save the changes to method #{ed['id']} '{ed['name']}' as a "
                                 "new version? The current version is kept in the archive.",
                            required=True)
        if not reason:
            return None
        mid = dlg.project.revise_method(ed["id"], name, definition, reason)
        dlg.editing = dlg.project.method(mid)
        dlg.setWindowTitle(dlg.windowTitle().split(" — editing")[0]
                           + f" — editing method #{mid} '{name}'")
    else:
        mid = dlg.project.save_method(name, definition, dlg.win.current_trial())
    dlg.method_id = mid
    dlg.win.statusBar().showMessage(f"Method #{mid} saved.")
    return mid


def start_editing(dlg, md: dict) -> None:
    dlg.editing = md
    dlg.setWindowTitle(dlg.windowTitle() + f" — editing method #{md['id']} '{md['name']}'")


def list_box(title: str, lst: SpectrumChecklist) -> QGroupBox:
    g = QGroupBox(title)
    v = QVBoxLayout(g)
    row = QHBoxLayout()
    for t, state in (("All", True), ("None", False)):
        b = QPushButton(t)
        b.clicked.connect(lambda _=False, s=state: lst.set_all(s))
        row.addWidget(b)
    row.addStretch(1)
    v.addLayout(row)
    v.addWidget(lst)
    lst.setMinimumHeight(130)
    return g


def results_rows(names, found, taken):
    rows, recs = [], []
    for n, f, t in zip(names, found, taken):
        rec = 100 * f / t if t not in (None, 0) else None
        if rec is not None:
            recs.append(rec)
        rows.append([n, float(f), "" if t is None else float(t),
                     "" if rec is None else float(rec)])
    return rows, recs


def summary_text(recs: list[float]) -> str:
    if len(recs) < 2:
        return ""
    d = describe(recs)
    return f"Recovery: mean {d['mean']:.2f} %, SD {d['sd']:.3f}, RSD {d['rsd']:.2f} % (n = {d['n']})"


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
        for r, sig in zip(rows, signals):
            r.insert(1, float(sig))
        headers = ["Spectrum", "Signal", f"Found ({comp})", "Taken", "Recovery %"]
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


# --------------------------------------------------------------------------- #
# Equation methods
# --------------------------------------------------------------------------- #
class CompoundChecks(QListWidget):
    def __init__(self, names: list[str]):
        super().__init__()
        for n in names:
            it = QListWidgetItem(n)
            it.setFlags(it.flags() | Qt.ItemIsUserCheckable)
            it.setCheckState(Qt.Checked)
            self.addItem(it)
        self.setMaximumHeight(90)

    def checked(self) -> list[str]:
        return [self.item(i).text() for i in range(self.count())
                if self.item(i).checkState() == Qt.Checked]


class EquationsDialog(Base):
    def __init__(self, win):
        super().__init__(win, "Equation methods — Vierordt / multi-wavelength / bivariate / AUC")
        self.resize(1300, 820)
        self.model: SignalEquations | None = None
        comps = self.project.compound_names()
        root = QSplitter()
        left = QWidget()
        ll = QVBoxLayout(left)
        ll.addWidget(QLabel("Responses at the chosen signals are modelled as R = C·K (+ b). "
                            "Use n ≥ number of compounds signals. Amplitude signals = "
                            "Vierordt / multi-λ least squares; with intercept = bivariate; "
                            "area signals = AUC method."))
        ll.itemAt(0).widget().setWordWrap(True)
        self.comps = CompoundChecks(comps)
        ll.addWidget(QLabel("Compounds"))
        ll.addWidget(self.comps)
        self.signals = PasteTable(3, 3)
        self.signals.setHorizontalHeaderLabels(["Kind (amplitude/area)", "λ1 (nm)", "λ2 (nm, area)"])
        fill_table(self.signals, ["Kind (amplitude/area)", "λ1 (nm)", "λ2 (nm, area only)"],
                   [["amplitude", "", ""]] * max(2, len(comps)), editable=True)
        ll.addWidget(QLabel("Signals"))
        ll.addWidget(self.signals)
        row = QHBoxLayout()
        add = QPushButton("Add signal")
        add.clicked.connect(lambda: self.signals.setRowCount(self.signals.rowCount() + 1))
        kaiser = QPushButton("Suggest λ pair (Kaiser)…")
        kaiser.clicked.connect(self._kaiser)
        self.intercept = QCheckBox("Intercept (bivariate)")
        row.addWidget(add)
        row.addWidget(kaiser)
        row.addWidget(self.intercept)
        ll.addLayout(row)
        g = QGroupBox("Optional processing before measuring (e.g. derivative)")
        gl = QVBoxLayout(g)
        self.pipe = PipelineEditor(spectrum_choices(self.project), name_of(self.project))
        gl.addWidget(self.pipe)
        ll.addWidget(g)
        self.cal, self.test = checklists(self.project)
        ll.addWidget(list_box("Calibration spectra (pure standards and/or mixtures)", self.cal))
        ll.addWidget(list_box("Spectra to determine", self.test))
        btns = QHBoxLayout()
        for t, fn in (("Fit", self._fit), ("Determine", self._predict),
                      ("Save method", self._save_method), ("Save results", self._save_results)):
            b = QPushButton(t)
            b.clicked.connect(fn)
            btns.addWidget(b)
        ll.addLayout(btns)
        root.addWidget(scroll(left))
        right = QWidget()
        rl = QVBoxLayout(right)
        self.info = QTableWidget()
        self.info.setMaximumHeight(220)
        rl.addWidget(QLabel("<b>Coefficients K (rows = compounds, columns = signals)</b>"))
        rl.addWidget(self.info)
        self.results = PasteTable()
        rl.addWidget(QLabel("<b>Results</b>"))
        rl.addWidget(self.results, 1)
        self.summary = QLabel()
        self.summary.setWordWrap(True)
        rl.addWidget(self.summary)
        root.addWidget(right)
        root.setSizes([560, 740])
        QVBoxLayout(self).addWidget(root)
        self.predictions = None

    def _signals(self) -> list[dict]:
        out = []
        for i in range(self.signals.rowCount()):
            items = [self.signals.item(i, j) for j in range(3)]
            kind = (items[0].text().strip().lower() if items[0] else "") or "amplitude"
            w1 = items[1].text().strip() if items[1] else ""
            w2 = items[2].text().strip() if items[2] else ""
            if not w1:
                continue
            if kind.startswith("area"):
                out.append({"kind": "area", "params": {"w1": float(w1), "w2": float(w2)}})
            else:
                out.append({"kind": "amplitude", "params": {"w1": float(w1)}})
        return out

    def _kaiser(self):
        comps = self.comps.checked()
        if len(comps) != 2:
            error(self, "Kaiser selection needs exactly two compounds.")
            return
        pure = {}
        for c in comps:
            recs = [r for r in self.project.records()
                    if r.concentrations.get(c, 0) > 0 and
                    all(v == 0 for k, v in r.concentrations.items() if k != c)]
            if not recs:
                error(self, f"No pure standard of {c} (concentration > 0, others 0 or unset).")
                return
            r = recs[-1]
            s = self.project.spectrum(r.id)
            pure[c] = s.copy(values=s.values / r.concentrations[c])
        best = kaiser_selection(pure, top=5)
        fill_table(self.signals, ["Kind (amplitude/area)", "λ1 (nm)", "λ2 (nm, area only)"],
                   [["amplitude", f"{best[0]['w1']:g}", ""], ["amplitude", f"{best[0]['w2']:g}", ""]],
                   editable=True)
        self.intercept.setChecked(True)
        self.summary.setText("Kaiser top pairs (|det|): " + "; ".join(
            f"{b['w1']:g}/{b['w2']:g} nm ({b['det']:.3g})" for b in best))

    def _fit(self):
        try:
            self.model = SignalEquations(self.comps.checked(), self._signals(),
                                         self.pipe.get_steps(), self.intercept.isChecked())
            ids = self.cal.checked_ids()
            info = self.model.fit(self.project.spectra(ids), self.project.resolver())
        except Exception as exc:
            error(self, exc)
            self.model = None
            return
        self.cal_ids = ids
        rows = [[c] + [float(v) for v in self.model.K[i]] for i, c in enumerate(self.model.compounds)]
        if self.model.intercept:
            rows.append(["intercept b"] + [float(v) for v in self.model.b])
        labels = [f"{s['kind']} {s['params']['w1']:g}" +
                  (f"–{s['params']['w2']:g}" if s["kind"] == "area" else "") for s in self.model.signals]
        fill_table(self.info, [""] + labels, rows)
        self.summary.setText(f"Condition number of K: {info['condition_number']:.3g} "
                             "(lower is better; > 100 means poorly resolved signals).")
        self.project.log("CALCULATE", "Fitted equation method", {"inputs": ids,
                                                                  "model": self.model.to_dict()})

    def _predict(self):
        if self.model is None:
            self._fit()
            if self.model is None:
                return
        ids = self.test.checked_ids()
        try:
            spectra = self.project.spectra(ids)
            pred = self.model.predict(spectra, self.project.resolver())
        except Exception as exc:
            error(self, exc)
            return
        comps = self.model.compounds
        rows, lines = [], []
        for s, p in zip(spectra, pred):
            row = [s.name]
            for c, v in zip(comps, p):
                t = s.concentrations.get(c)
                row += [float(v), "" if t is None else float(t),
                        "" if not t else float(100 * v / t)]
            rows.append(row)
        for j, c in enumerate(comps):
            recs = [100 * p[j] / s.concentrations[c] for s, p in zip(spectra, pred)
                    if s.concentrations.get(c)]
            if len(recs) > 1:
                lines.append(f"{c}: " + summary_text(recs))
        head = ["Spectrum"] + [h for c in comps for h in (f"{c} found", f"{c} taken", f"{c} rec %")]
        fill_table(self.results, head, rows)
        self.summary.setText("\n".join(lines))
        self.predictions = {"ids": ids, "compounds": comps, "found": pred.tolist()}
        self.project.log("CALCULATE", f"Determined {', '.join(comps)} in {len(ids)} spectra "
                                      "(equation method)", {"inputs": ids, **self.predictions})

    def _save_method(self):
        if self.model is None:
            error(self, "Fit first.")
            return
        ed = getattr(self, "editing", None)
        name, ok = QInputDialog.getText(self, "Save method", "Name:",
                                        text=ed["name"] if ed else "Equation method")
        if ok:
            d = self.model.to_dict()
            d["calibration_ids"] = self.cal_ids
            store_method(self, name.strip() or "Equation method", d)

    def load_method(self, md: dict) -> None:
        d = md["definition"]
        for i in range(self.comps.count()):
            it = self.comps.item(i)
            it.setCheckState(Qt.Checked if it.text() in d["compounds"] else Qt.Unchecked)
        rows = [[s["kind"], f"{s['params']['w1']:g}",
                 f"{s['params']['w2']:g}" if s["kind"] == "area" else ""] for s in d["signals"]]
        fill_table(self.signals, ["Kind (amplitude/area)", "λ1 (nm)", "λ2 (nm, area only)"],
                   rows, editable=True)
        self.pipe.set_steps(d.get("steps", []))
        self.intercept.setChecked(bool(d.get("intercept")))
        if d.get("calibration_ids"):
            check_ids(self.cal, d["calibration_ids"])
        start_editing(self, md)

    def _save_results(self):
        if not self.predictions:
            error(self, "Determine first.")
            return
        self.project.save_result("Equation method results", "equations",
                                 {"model": self.model.to_dict(), **self.predictions,
                                  "calibration_ids": list(getattr(self, "cal_ids", []))},
                                 self.win.current_trial(), getattr(self, "method_id", None),
                                 self.predictions["ids"])


# --------------------------------------------------------------------------- #
# Special binary methods
# --------------------------------------------------------------------------- #
def wl_spin(v: float) -> QDoubleSpinBox:
    s = QDoubleSpinBox()
    s.setRange(0, 100000)
    s.setDecimals(2)
    s.setSuffix(" nm")
    s.setValue(v)
    return s


class SpecialDialog(Base):
    def __init__(self, win):
        super().__init__(win, "Binary two-signal methods")
        self.resize(1250, 800)
        self.comps = self.project.compound_names()
        tabs = self.tabs = QTabWidget()
        # zero-order spectra
        tabs.addTab(self._q_tab(), "Zero order: Q-analysis")
        tabs.addTab(self._as_tab(), "Zero order: absorbance subtraction")
        tabs.addTab(self._aas_tab(), "Zero order: advanced absorbance subtraction")
        # ratio spectra (amplitudes)
        tabs.addTab(self._am_tab(), "Ratio: amplitude modulation")
        tabs.addTab(self._iam_tab(), "Ratio: induced amplitude modulation")
        # standard addition
        tabs.addTab(self._h_tab(), "Standard addition: H-point (HPSAM)")
        lay = QVBoxLayout(self)
        lay.addWidget(tabs, 1)
        self.results = PasteTable()
        self.summary = QLabel()
        self.summary.setWordWrap(True)
        lay.addWidget(QLabel("<b>Results</b>"))
        lay.addWidget(self.results, 1)
        lay.addWidget(self.summary)
        row = QHBoxLayout()
        save = QPushButton("Save results")
        save.clicked.connect(self._save)
        row.addStretch(1)
        row.addWidget(save)
        lay.addLayout(row)
        self.last = None

    def _combo(self, idx=0) -> QComboBox:
        c = QComboBox()
        c.addItems(self.comps)
        c.setCurrentIndex(min(idx, max(0, len(self.comps) - 1)))
        return c

    def _lists(self, *titles):
        recs = self.project.records()
        out = []
        for _ in titles:
            lst = SpectrumChecklist()
            lst.populate(recs)
            out.append(lst)
        return out

    def _find_iso(self, x_list, y_list, cx, cy, target: QDoubleSpinBox):
        try:
            a = self.project.spectrum(x_list.checked_ids()[0])
            b = self.project.spectrum(y_list.checked_ids()[0])
            pts = uv.isoabsorptive_points(a, b, a.concentrations[cx.currentText()],
                                          b.concentrations[cy.currentText()])
        except (IndexError, KeyError):
            error(self, "Check at least one standard of each compound (with concentrations).")
            return
        if not pts:
            error(self, "No isoabsorptive point found.")
            return
        choice, ok = QInputDialog.getItem(self, "Isoabsorptive points", "Use:",
                                          [f"{p:.2f}" for p in pts], 0, False)
        if ok:
            target.setValue(float(choice))

    def _layout(self, form_rows, lists):
        w = QWidget()
        h = QHBoxLayout(w)
        f = QFormLayout()
        for label, widget in form_rows:
            f.addRow(label, widget)
        fw = QWidget()
        fw.setLayout(f)
        h.addWidget(fw)
        for title, lst in lists:
            h.addWidget(list_box(title, lst))
        return w

    def _q_tab(self):
        self.q_x, self.q_y = self._combo(0), self._combo(1)
        self.q_iso, self.q_w2 = wl_spin(260), wl_spin(245)
        self.q_xs, self.q_ys, self.q_mix = self._lists(1, 2, 3)
        iso_btn = QPushButton("Find isoabsorptive point")
        iso_btn.clicked.connect(lambda: self._find_iso(self.q_xs, self.q_ys, self.q_x, self.q_y,
                                                       self.q_iso))
        run = QPushButton("Calculate")
        run.clicked.connect(self._run_q)
        return self._layout([("Compound X", self.q_x), ("Compound Y", self.q_y),
                             ("Isoabsorptive λ", self.q_iso), ("", iso_btn),
                             ("Second λ (λmax of X)", self.q_w2), ("", run)],
                            [("X standards", self.q_xs), ("Y standards", self.q_ys),
                             ("Mixtures", self.q_mix)])

    def _run_q(self):
        x, y = self.q_x.currentText(), self.q_y.currentText()
        try:
            xs = self.project.spectra(self.q_xs.checked_ids())
            ys = self.project.spectra(self.q_ys.checked_ids())
            iso, w2 = self.q_iso.value(), self.q_w2.value()
            ax_iso = uv.absorptivity(xs, x, iso)
            ax2 = uv.absorptivity(xs, x, w2)
            ay2 = uv.absorptivity(ys, y, w2)
            out = [(m, uv.q_analysis(m, iso, w2, ax_iso, ax2, ay2))
                   for m in self.project.spectra(self.q_mix.checked_ids())]
        except Exception as exc:
            error(self, exc)
            return
        self._show("Q-analysis", x, y, out, {"iso": iso, "w2": w2, "ax_iso": ax_iso,
                                              "ax2": ax2, "ay2": ay2},
                   f"aX(iso) = {ax_iso:.5g}, aX(λ2) = {ax2:.5g}, aY(λ2) = {ay2:.5g}")

    def _as_tab(self):
        self.as_x, self.as_y = self._combo(0), self._combo(1)
        self.as_iso, self.as_w2 = wl_spin(260), wl_spin(320)
        self.as_xs, self.as_all, self.as_mix = self._lists(1, 2, 3)
        run = QPushButton("Calculate")
        run.clicked.connect(self._run_as)
        return self._layout([("X (absorbs alone at λ2)", self.as_x), ("Y", self.as_y),
                             ("Isoabsorptive λ", self.as_iso), ("λ2", self.as_w2), ("", run)],
                            [("Pure X standards (amplitude factor)", self.as_xs),
                             ("Standards for iso calibration", self.as_all),
                             ("Mixtures", self.as_mix)])

    def _run_as(self):
        try:
            m = uv.AbsorbanceSubtraction(self.as_iso.value(), self.as_w2.value())
            m.fit(self.project.spectra(self.as_xs.checked_ids()),
                  self.project.spectra(self.as_all.checked_ids()))
            out = [(s, m.predict(s)) for s in self.project.spectra(self.as_mix.checked_ids())]
        except Exception as exc:
            error(self, exc)
            return
        self._show("Absorbance subtraction", self.as_x.currentText(), self.as_y.currentText(),
                   out, {"iso": m.iso, "w2": m.w2, "AF": m.amplitude_factor},
                   f"Amplitude factor = {m.amplitude_factor:.5g}; iso calibration r = "
                   f"{m.iso_regression.r:.5f}")

    def _am_tab(self):
        self.am_x, self.am_y = self._combo(0), self._combo(1)
        self.am_iso = wl_spin(260)
        self.am_p1, self.am_p2 = wl_spin(330), wl_spin(350)
        self.am_div = QComboBox()
        for sid, label in spectrum_choices(self.project)():
            self.am_div.addItem(label, sid)
        self.am_divc = QDoubleSpinBox()
        self.am_divc.setRange(0, 1e6)
        self.am_divc.setDecimals(4)
        self.am_divc.setValue(10)
        self.am_std, self.am_mix = self._lists(1, 2)
        run = QPushButton("Calculate")
        run.clicked.connect(self._run_am)
        return self._layout([("X", self.am_x), ("Y (extended, divisor)", self.am_y),
                             ("Isoabsorptive λ", self.am_iso), ("Plateau from", self.am_p1),
                             ("Plateau to", self.am_p2), ("Divisor spectrum (Y)", self.am_div),
                             ("Divisor concentration", self.am_divc), ("", run)],
                            [("Standards (total calibration)", self.am_std),
                             ("Mixtures", self.am_mix)])

    def _run_am(self):
        try:
            div = self.project.spectrum(self.am_div.currentData())
            m = uv.AmplitudeModulation(self.am_iso.value(), (self.am_p1.value(), self.am_p2.value()),
                                       self.am_divc.value())
            m.fit(self.project.spectra(self.am_std.checked_ids()), div)
            out = [(s, m.predict(s, div)) for s in self.project.spectra(self.am_mix.checked_ids())]
        except Exception as exc:
            error(self, exc)
            return
        self._show("Amplitude modulation", self.am_x.currentText(), self.am_y.currentText(), out,
                   {"iso": m.iso, "plateau": m.plateau, "divisor": self.am_div.currentData()},
                   f"Total calibration at iso: r = {m.total_regression.r:.5f}")

    def _iam_tab(self):
        self.iam_x, self.iam_y = self._combo(0), self._combo(1)
        self.iam_w = wl_spin(250)
        self.iam_p1, self.iam_p2 = wl_spin(330), wl_spin(350)
        self.iam_xs, self.iam_ys, self.iam_mix = self._lists(1, 2, 3)
        run = QPushButton("Calculate")
        run.clicked.connect(self._run_iam)
        info = QLabel("Mixture ÷ unit-concentration Y′ (averaged from the Y standards).\n"
                      "Plateau (Y only) → C_Y; P(λ) − C_Y divided by the factor\n"
                      "aX/aY at λ (from X standards) → C_X. No isoabsorptive point needed.")
        return self._layout([("X", self.iam_x), ("Y (extended)", self.iam_y),
                             ("λ for X", self.iam_w), ("Plateau from", self.iam_p1),
                             ("Plateau to", self.iam_p2), (info, QLabel("")), ("", run)],
                            [("X standards", self.iam_xs), ("Y standards (divisor)", self.iam_ys),
                             ("Mixtures", self.iam_mix)])

    def _run_iam(self):
        x, y = self.iam_x.currentText(), self.iam_y.currentText()
        try:
            div = uv.InducedAmplitudeModulation.unit_divisor(
                self.project.spectra(self.iam_ys.checked_ids()), y)
            m = uv.InducedAmplitudeModulation(self.iam_w.value(),
                                              (self.iam_p1.value(), self.iam_p2.value()))
            m.fit(self.project.spectra(self.iam_xs.checked_ids()), x, div)
            out = [(s, m.predict(s, div)) for s in self.project.spectra(self.iam_mix.checked_ids())]
        except Exception as exc:
            error(self, exc)
            return
        self._show("Induced amplitude modulation", x, y, out,
                   {"wavelength": m.wavelength, "plateau": m.plateau, "factor": m.factor},
                   f"Factor aX/aY at {m.wavelength:g} nm = {m.factor:.5g}")

    def _aas_tab(self):
        self.aas_x, self.aas_y = self._combo(0), self._combo(1)
        self.aas_w1, self.aas_w2 = wl_spin(260), wl_spin(245)
        self.aas_xs, self.aas_ys, self.aas_mix = self._lists(1, 2, 3)
        find = QPushButton("Find λ1 where Y equals its value at λ2")
        find.clicked.connect(self._aas_find)
        run = QPushButton("Calculate")
        run.clicked.connect(self._run_aas)
        return self._layout([("X", self.aas_x), ("Y (interferent)", self.aas_y),
                             ("λ2 (e.g. λmax of X)", self.aas_w2), ("λ1", self.aas_w1),
                             ("", find), ("", run)],
                            [("X standards", self.aas_xs), ("Y standards", self.aas_ys),
                             ("Mixtures", self.aas_mix)])

    def _aas_find(self):
        ids = self.aas_ys.checked_ids()
        if not ids:
            error(self, "Check at least one Y standard.")
            return
        cands = uv.equal_amplitude_wavelengths(self.project.spectrum(ids[-1]), self.aas_w2.value())
        if not cands:
            error(self, "Y does not reach the same absorbance at any other wavelength.")
            return
        choice, ok = QInputDialog.getItem(self, "Equal absorbance of Y", "λ1:",
                                          [f"{c:.2f}" for c in cands], 0, False)
        if ok:
            self.aas_w1.setValue(float(choice))

    def _run_aas(self):
        x, y = self.aas_x.currentText(), self.aas_y.currentText()
        try:
            m = uv.AdvancedAbsorbanceSubtraction(self.aas_w1.value(), self.aas_w2.value())
            m.fit(self.project.spectra(self.aas_xs.checked_ids()), x,
                  self.project.spectra(self.aas_ys.checked_ids()), y)
            out = [(s, m.predict(s)) for s in self.project.spectra(self.aas_mix.checked_ids())]
        except Exception as exc:
            error(self, exc)
            return
        self._show("Advanced absorbance subtraction", x, y, out,
                   {"w1": m.w1, "w2": m.w2, "AF": m.amplitude_factor},
                   f"Amplitude factor = {m.amplitude_factor:.5g}; X (ΔA) r = "
                   f"{m.x_regression.r:.5f}; Y at λ2 r = {m.y_regression.r:.5f}")

    def _h_tab(self):
        self.h_x = self._combo(0)
        self.h_w1, self.h_w2 = wl_spin(260), wl_spin(300)
        (self.h_set,) = self._lists(1)
        run = QPushButton("Calculate")
        run.clicked.connect(self._run_h)
        return self._layout([("Analyte X", self.h_x),
                             ("λ1 (interferent equal absorbance)", self.h_w1), ("λ2", self.h_w2),
                             (QLabel("The X concentration of each spectrum\nis taken as the "
                                     "amount added."), QLabel("")), ("", run)],
                            [("Sample + standard additions", self.h_set)])

    def _run_h(self):
        x = self.h_x.currentText()
        try:
            spectra = self.project.spectra(self.h_set.checked_ids())
            added = [s.concentrations.get(x, 0.0) for s in spectra]
            r = uv.hpsam(added, spectra, self.h_w1.value(), self.h_w2.value())
        except Exception as exc:
            error(self, exc)
            return
        fill_table(self.results, ["Quantity", "Value"],
                   [[f"{x} in sample (−C_H)", r["X"]], ["A_H (interferent signal)", r["A_H"]],
                    ["Line λ1: slope", r["line1"]["slope"]], ["Line λ1: intercept", r["line1"]["intercept"]],
                    ["Line λ2: slope", r["line2"]["slope"]], ["Line λ2: intercept", r["line2"]["intercept"]]])
        self.summary.setText("")
        self.last = ("HPSAM", {"added": added, "ids": self.h_set.checked_ids(), **r})
        self.project.log("CALCULATE", f"HPSAM: {x} = {r['X']:.4g}", {"inputs": self.h_set.checked_ids(), **r})

    def _show(self, method, x, y, out, params, info):
        rows, rx, ry = [], [], []
        for s, r in out:
            tx, ty = s.concentrations.get(x), s.concentrations.get(y)
            rows.append([s.name, r["X"], "" if tx is None else tx,
                         "" if not tx else 100 * r["X"] / tx,
                         r["Y"], "" if ty is None else ty, "" if not ty else 100 * r["Y"] / ty])
            if tx:
                rx.append(100 * r["X"] / tx)
            if ty:
                ry.append(100 * r["Y"] / ty)
        fill_table(self.results, ["Spectrum", f"{x} found", f"{x} taken", f"{x} rec %",
                                  f"{y} found", f"{y} taken", f"{y} rec %"], rows)
        self.summary.setText(info + "\n" + (f"{x}: " + summary_text(rx) if len(rx) > 1 else "")
                             + ("\n" + f"{y}: " + summary_text(ry) if len(ry) > 1 else ""))
        ids = [s.metadata.get("id") for s, _ in out]
        self.last = (method, {"ids": ids, "params": params, "X": x, "Y": y,
                              "found": [r for _, r in out]})
        self.project.log("CALCULATE", f"{method}: determined {x} and {y} in {len(out)} mixtures",
                         {"inputs": ids, "params": params})

    def _save(self):
        if not self.last:
            error(self, "Calculate first.")
            return
        method, data = self.last
        self.project.save_result(method, "binary", data, self.win.current_trial(),
                                 inputs=data.get("ids"))


# --------------------------------------------------------------------------- #
# Progressive resolution (amplitude centering, absorption factor)
# --------------------------------------------------------------------------- #
pure_standards = uv.pure_standards


def _cell(t: QTableWidget, i: int, j: int) -> str:
    it = t.item(i, j)
    return it.text().strip() if it is not None else ""


def _num_cell(t: QTableWidget, i: int, j: int) -> float | None:
    v = _cell(t, i, j).replace(",", ".")
    if not v:
        return None
    try:
        return float(v)
    except ValueError:
        raise ValueError(f"row {i + 1}: '{v}' is not a number") from None


class ProgressiveDialog(Base):
    """Ternary (and binary) progressive methods that resolve several compounds
    from one ratio spectrum or a chain of absorbance factors."""

    def __init__(self, win):
        super().__init__(win, "Progressive resolution (AAC, MACM, RIDSS, CV-AD, MAFM)")
        self.resize(1300, 860)
        self.comps = self.project.compound_names()
        self.last = None
        self.fitted = None          # (method, calibration ids) of the last calculation
        tabs = self.tabs = QTabWidget()
        tabs.addTab(self._ac_tab(), "Ratio: amplitude centering (AAC, MACM, RIDSS, CV-AD)")
        tabs.addTab(self._af_tab(), "Zero order: successive absorption factor (MAFM)")
        lay = QVBoxLayout(self)
        pre = QHBoxLayout()
        self.sm_win = QSpinBox()
        self.sm_win.setRange(0, 501)
        self.sm_win.setSingleStep(2)
        self.sm_win.setSpecialValueText("off")
        self.sm_win.setSuffix(" points")
        self.sm_win.setToolTip("Savitzky–Golay smoothing applied to the standards, the divisor "
                               "and the samples before the method (odd window; 0 = off).")
        self.sm_ord = QSpinBox()
        self.sm_ord.setRange(1, 6)
        self.sm_ord.setValue(2)
        for w in (QLabel("Pre-processing — Savitzky–Golay smoothing window:"), self.sm_win,
                  QLabel("order"), self.sm_ord):
            pre.addWidget(w)
        pre.addStretch(1)
        lay.addLayout(pre)
        lay.addWidget(tabs, 3)
        self.results = PasteTable()
        self.summary = QLabel()
        self.summary.setWordWrap(True)
        lay.addWidget(QLabel("<b>Results</b>"))
        lay.addWidget(self.results, 1)
        lay.addWidget(self.summary)
        row = QHBoxLayout()
        row.addStretch(1)
        for text, fn, tip in (
                ("Save method", self._save_method,
                 "Save the calibrated method with its divisor so it can be applied to new "
                 "samples from Methods → Saved methods, exported as a .spmodel file and "
                 "reported"),
                ("Save results", self._save, "")):
            b = QPushButton(text)
            b.setToolTip(tip)
            b.clicked.connect(fn)
            row.addWidget(b)
        lay.addLayout(row)

    def _steps(self) -> list[dict]:
        w = self.sm_win.value()
        if not w:
            return []
        if w % 2 == 0:
            w += 1
        if self.sm_ord.value() >= w:
            raise ValueError("the smoothing window must be larger than the polynomial order")
        return [{"op": "smooth_sg", "params": {"window": w, "polyorder": self.sm_ord.value()}}]

    def _set_steps(self, steps: list[dict]) -> None:
        sg = next((st for st in steps if st["op"] == "smooth_sg"), None)
        self.sm_win.setValue(int(sg["params"]["window"]) if sg else 0)
        if sg:
            self.sm_ord.setValue(int(sg["params"]["polyorder"]))

    def _lists(self):
        recs = self.project.records()
        std, mix = SpectrumChecklist(), SpectrumChecklist()
        std.populate(recs, {r.id for r in recs if r.kind == "raw" and r.role in CAL_ROLES})
        mix.populate(recs, {r.id for r in recs if r.kind == "raw" and r.role in TEST_ROLES})
        return std, mix

    def _optional_combo(self, none_label="(none)") -> QComboBox:
        c = QComboBox()
        c.addItem(none_label)
        c.addItems(self.comps)
        return c

    def _ac_tab(self):
        w = QWidget()
        h = QHBoxLayout(w)
        f = QFormLayout()
        self.ac_comps = CompoundChecks(self.comps)
        self.ac_div = QComboBox()
        for sid, label in spectrum_choices(self.project)():
            self.ac_div.addItem(label, sid)
        self.ac_divc = self._optional_combo()
        self.ac_w = wl_spin(275.5)
        self.ac_plateau = QCheckBox("Plateau of the divisor compound")
        self.ac_p1, self.ac_p2 = wl_spin(350), wl_spin(380)
        self.ac_diff = PasteTable(4, 4)
        self.ac_diff.setHorizontalHeaderLabels(["Compound", "λ1", "λ2", "Factor from"])
        self.ac_diff.setToolTip("One row per compound found by amplitude difference: "
                                "P(λ1) − F·P(λ2) of the ratio spectrum, where the others cancel. "
                                "'Factor from' = the interferent whose equality factor "
                                "F = P(λ1)/P(λ2) is used (blank: F = 1).")
        self.ac_diff.setMinimumHeight(130)
        self.ac_diff.setMaximumHeight(150)
        self.ac_sub = self._optional_combo()
        self.ac_unified = QCheckBox("Unified regression at λc (isoabsorptive point)")
        info = QLabel("Ratio spectrum = mixture ÷ divisor. Plateau → divisor compound; "
                      "λ pairs → postulated amplitudes at λc; the remaining compound by "
                      "amplitude subtraction at λc. Pure standards are grouped by their "
                      "single non-zero concentration.")
        info.setWordWrap(True)
        run = QPushButton("Calculate")
        run.clicked.connect(self._run_ac)
        for label, wid in (("Compounds", self.ac_comps), ("Divisor spectrum", self.ac_div),
                           ("Divisor compound", self.ac_divc), ("Common λc", self.ac_w),
                           ("", self.ac_plateau), ("Plateau from", self.ac_p1),
                           ("Plateau to", self.ac_p2), ("Amplitude differences", self.ac_diff),
                           ("Amplitude subtraction for", self.ac_sub), ("", self.ac_unified),
                           ("", info), ("", run)):
            f.addRow(label, wid)
        fw = QWidget()
        fw.setLayout(f)
        h.addWidget(scroll(fw), 2)
        self.ac_std, self.ac_mix = self._lists()
        h.addWidget(list_box("Pure standards", self.ac_std), 1)
        h.addWidget(list_box("Mixtures / samples", self.ac_mix), 1)
        return w

    def _ac_method(self) -> uv.AmplitudeCentering:
        comps = self.ac_comps.checked()
        diffs = {}
        for i in range(self.ac_diff.rowCount()):
            c = _cell(self.ac_diff, i, 0)
            if not c:
                continue
            if c not in comps:
                raise ValueError(f"row {i + 1}: '{c}' is not one of the checked compounds")
            w1, w2 = _num_cell(self.ac_diff, i, 1), _num_cell(self.ac_diff, i, 2)
            if w1 is None or w2 is None:
                raise ValueError(f"row {i + 1}: give λ1 and λ2")
            src = _cell(self.ac_diff, i, 3) or None
            diffs[c] = {"w1": w1, "w2": w2, "factor_from": src}
        divc = self.ac_divc.currentText() if self.ac_divc.currentIndex() > 0 else None
        sub = self.ac_sub.currentText() if self.ac_sub.currentIndex() > 0 else None
        return uv.AmplitudeCentering(
            self.ac_w.value(), comps, sub, divc,
            (self.ac_p1.value(), self.ac_p2.value()) if self.ac_plateau.isChecked() else None,
            diffs, self.ac_unified.isChecked())

    def _run_ac(self):
        try:
            m = self._ac_method()
            if self.ac_div.currentData() is None:
                raise ValueError("choose the divisor spectrum")
            m.divisor = self.ac_div.currentData()
            m.steps = self._steps()
            res = self.project.resolver()
            cal = self.ac_std.checked_ids()
            info = m.fit_spectra(self.project.spectra(cal), res)
            mixes = self.project.spectra(self.ac_mix.checked_ids())
            out = [(s, m.predict_one(s, res)) for s in mixes]
        except Exception as exc:
            error(self, exc)
            return
        lines = [f"{c}: calibration at λc r = {r:.5f}" for c, r in info["r"].items()]
        lines += [f"{c}: amplitude-difference line r = {r:.5f}"
                  for c, r in info["difference_r"].items()]
        lines += [f"Equality factor F of {m.differences[c]['factor_from']} (used for {c}): "
                  f"{f:.5g}" for c, f in info["factors"].items()]
        self.fitted = (m, cal)
        self._show("Amplitude centering", m.compounds, out, m.to_dict(), lines)

    def _af_tab(self):
        w = QWidget()
        h = QHBoxLayout(w)
        f = QFormLayout()
        self.af_table = PasteTable(4, 3)
        self.af_table.setHorizontalHeaderLabels(["Compound (in order)", "λ where it is added",
                                                 "Quantitation λ (optional)"])
        self.af_table.setToolTip("Row 1: compound absorbing ALONE at its λ; row 2: compound "
                                 "overlapped only by row 1 at its λ; row 3: all overlap.")
        for i, c in enumerate(self.comps[:3]):
            self.af_table.setItem(i, 0, QTableWidgetItem(c))
        info = QLabel("Modified absorption factor method: absorption factors F = A(λ)/A(λk) "
                      "of each pure compound remove its contribution at the following "
                      "wavelengths; each compound is read from its own calibration.")
        info.setWordWrap(True)
        run = QPushButton("Calculate")
        run.clicked.connect(self._run_af)
        f.addRow("Order and wavelengths", self.af_table)
        f.addRow("", info)
        f.addRow("", run)
        fw = QWidget()
        fw.setLayout(f)
        h.addWidget(fw, 2)
        self.af_std, self.af_mix = self._lists()
        h.addWidget(list_box("Pure standards", self.af_std), 1)
        h.addWidget(list_box("Mixtures / samples", self.af_mix), 1)
        return w

    def _run_af(self):
        try:
            order, quant = [], {}
            for i in range(self.af_table.rowCount()):
                c = _cell(self.af_table, i, 0)
                if not c:
                    continue
                lam = _num_cell(self.af_table, i, 1)
                if lam is None:
                    raise ValueError(f"row {i + 1}: give the wavelength")
                order.append((c, lam))
                q = _num_cell(self.af_table, i, 2)
                if q is not None:
                    quant[c] = q
            m = uv.AbsorptionFactorMethod(order, quant, steps=self._steps())
            comps = [c for c, _ in order]
            cal = self.af_std.checked_ids()
            res = self.project.resolver()
            info = m.fit_spectra(self.project.spectra(cal), res)
            out = [(s, m.predict_one(s, res))
                   for s in self.project.spectra(self.af_mix.checked_ids())]
        except Exception as exc:
            error(self, exc)
            return
        lines = [f"{c}: calibration r = {r:.5f}" for c, r in info["r"].items()]
        lines += [f"Factors of {c}: " + ", ".join(f"F({k} nm) = {v:.4f}" for k, v in fs.items())
                  for c, fs in info["factors"].items()]
        self.fitted = (m, cal)
        self._show("Absorption factor method", comps, out, m.to_dict(), lines)

    def _show(self, method: str, comps: list[str], out, params: dict, lines: list[str]):
        rows, recs = [], {c: [] for c in comps}
        for s, r in out:
            row = [s.name]
            for c in comps:
                t = s.concentrations.get(c)
                rec = 100 * r[c] / t if t else None
                if rec is not None:
                    recs[c].append(rec)
                row += [r[c], "" if t is None else t, "" if rec is None else rec]
            rows.append(row)
        fill_table(self.results, ["Spectrum"] + [h for c in comps for h in
                                                 (f"{c} found", f"{c} taken", f"{c} rec %")],
                   rows)
        lines = lines + [f"{c}: " + summary_text(v) for c, v in recs.items() if len(v) > 1]
        self.summary.setText("\n".join(lines))
        self.results.export_title = method
        self.results.export_notes = lines
        ids = [s.metadata.get("id") for s, _ in out]
        self.last = (method, {"ids": ids, "compounds": comps, "params": params,
                              "found": [[r[c] for c in comps] for _, r in out]})
        self.project.log("CALCULATE", f"{method}: determined {', '.join(comps)} in "
                                      f"{len(out)} spectra", {"inputs": ids, "params": params})

    def _save(self):
        if not self.last:
            error(self, "Calculate first.")
            return
        method, data = self.last
        if self.fitted:
            data = {**data, "model": self.fitted[0].to_dict(),
                    "calibration_ids": list(self.fitted[1])}
        self.project.save_result(method, "progressive", data, self.win.current_trial(),
                                 getattr(self, "method_id", None), data.get("ids"))
        self.win.statusBar().showMessage("Results saved.")

    def _save_method(self):
        if not self.fitted:
            error(self, "Calculate first.")
            return
        m, cal = self.fitted
        kind = ("Amplitude centering" if isinstance(m, uv.AmplitudeCentering)
                else "Absorption factor")
        ed = getattr(self, "editing", None)
        name, ok = QInputDialog.getText(self, "Save method", "Method name:",
                                        text=ed["name"] if ed else f"{kind}: {m.describe()}")
        if not ok:
            return
        m.name = name.strip() or kind
        d = m.to_dict()
        d["calibration_ids"] = list(cal)
        store_method(self, m.name, d)

    def load_method(self, md: dict) -> None:
        d = md["definition"]
        self._set_steps(d.get("steps", []))
        cal = d.get("calibration_ids") or []
        if d["type"] == "amplitude_centering":
            self.tabs.setCurrentIndex(0)
            for i in range(self.ac_comps.count()):
                it = self.ac_comps.item(i)
                it.setCheckState(Qt.Checked if it.text() in d["compounds"] else Qt.Unchecked)
            if d.get("divisor") is not None:
                k = self.ac_div.findData(int(d["divisor"]))
                if k >= 0:
                    self.ac_div.setCurrentIndex(k)
            self.ac_divc.setCurrentText(d.get("divisor_compound") or "(none)")
            self.ac_w.setValue(float(d["wavelength"]))
            self.ac_plateau.setChecked(bool(d.get("plateau")))
            if d.get("plateau"):
                self.ac_p1.setValue(float(d["plateau"][0]))
                self.ac_p2.setValue(float(d["plateau"][1]))
            self.ac_diff.clearContents()
            for i, (c, df) in enumerate(d.get("differences", {}).items()):
                for j, v in enumerate((c, f"{df['w1']:g}", f"{df['w2']:g}",
                                       df.get("factor_from") or "")):
                    self.ac_diff.setItem(i, j, QTableWidgetItem(v))
            self.ac_sub.setCurrentText(d.get("subtract") or "(none)")
            self.ac_unified.setChecked(bool(d.get("unified")))
            if cal:
                check_ids(self.ac_std, cal)
        else:
            self.tabs.setCurrentIndex(1)
            self.af_table.clearContents()
            for i, (c, w) in enumerate(d["order"]):
                q = d.get("quant", {}).get(c)
                for j, v in enumerate((c, f"{w:g}", "" if q is None else f"{q:g}")):
                    self.af_table.setItem(i, j, QTableWidgetItem(v))
            if cal:
                check_ids(self.af_std, cal)
        start_editing(self, md)


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


# --------------------------------------------------------------------------- #
# Saved methods & results
# --------------------------------------------------------------------------- #
def determine(project, d: dict, ids: list[int], refit: bool = False):
    """Apply a method/model definition to spectra ``ids``.

    Returns (compounds, found rows, method). The method is calibrated on its
    ``calibration_ids`` when it carries no calibration — or always with
    ``refit`` (e.g. after calibration concentrations were corrected)."""
    res = project.resolver()
    m = method_from_definition(d)
    spectra = project.spectra(ids)
    cal = d.get("calibration_ids") or []
    if refit and not cal:
        raise ValueError("this method does not record its calibration spectra, so it cannot "
                         "be refitted — untick 'Refit'")
    if d["type"] == "univariate":
        if m.regression is None or refit:
            m.calibrate(project.spectra(cal), res)
        return [m.compound], [[m.predict(s, res)] for s in spectra], m
    if d["type"] == "equations":
        if m.K is None or refit:
            m.fit(project.spectra(cal), res)
        return m.compounds, m.predict(spectra, res).tolist(), m
    if d["type"] in uv.PROGRESSIVE_TYPES:
        if not m.is_fitted or refit:
            m.fit_spectra(project.spectra(cal), res)
        return m.compounds, m.predict_spectra(spectra, res).tolist(), m
    if not m.is_fitted or refit:
        m.fit(project.spectra(cal), res)
    return m.compounds, m.predict(spectra, res).tolist(), m


def method_from_definition(d: dict):
    t = d.get("type")
    if t in uv.PROGRESSIVE_TYPES:
        return uv.progressive_from_dict(d)
    if t == "univariate":
        return uv.UnivariateMethod.from_dict(d)
    if t == "equations":
        return SignalEquations.from_dict(d)
    if t == "spectral":
        return SpectralModel.from_dict(d)
    raise ValueError(f"unknown method type {t}")


class SavedDialog(Base):
    def __init__(self, win):
        super().__init__(win, "Saved methods & results")
        self.resize(1100, 700)
        lay = QVBoxLayout(self)
        tabs = QTabWidget()
        mw = QWidget()
        ml = QVBoxLayout(mw)
        self.mtable = QTableWidget()
        self.mtable.itemSelectionChanged.connect(self._show_method)
        ml.addWidget(self.mtable)
        row = QHBoxLayout()
        apply_btn = QPushButton("Apply to selected spectra in the project tree")
        apply_btn.clicked.connect(self._apply)
        edit = QPushButton("Edit…")
        edit.setToolTip("Open the method in its dialog with all its settings; change anything, "
                        "recalibrate and Save method. The change is stored as a new version "
                        "(the old one is archived, the reason is recorded).")
        edit.clicked.connect(self._edit)
        ren = QPushButton("Rename…")
        ren.clicked.connect(self._rename)
        arch = QPushButton("Archive method…")
        arch.clicked.connect(self._archive)
        exp = QPushButton("Export model file…")
        exp.setToolTip("Save the calibrated method/model to a .spmodel file to use in "
                       "another project or on another PC")
        exp.clicked.connect(self._export_model)
        imp = QPushButton("Import model file…")
        imp.clicked.connect(self._import_model)
        for b in (apply_btn, edit, ren, exp, imp, arch):
            row.addWidget(b)
        row.addStretch(1)
        ml.addLayout(row)
        tabs.addTab(mw, "Methods")
        rw = QWidget()
        rl = QVBoxLayout(rw)
        self.rtable = QTableWidget()
        self.rtable.itemSelectionChanged.connect(self._show_result)
        rl.addWidget(self.rtable)
        xl = QPushButton("Export all results to Excel…")
        xl.clicked.connect(self.win.export_results)
        redit = QPushButton("Edit…")
        redit.setToolTip("Rename, add notes, exclude spectra from the statistics, or "
                         "recalculate with the current method and data. Stored as a new "
                         "version (the old one is archived, the reason is recorded).")
        redit.clicked.connect(self._edit_result)
        rarch = QPushButton("Archive result…")
        rarch.clicked.connect(self._archive_result)
        rrow = QHBoxLayout()
        rrow.addWidget(redit)
        rrow.addWidget(rarch)
        rrow.addWidget(xl)
        rrow.addStretch(1)
        rl.addLayout(rrow)
        tabs.addTab(rw, "Results")
        lay.addWidget(tabs, 1)
        self.detail = QPlainTextEdit()
        self.detail.setReadOnly(True)
        lay.addWidget(self.detail, 1)
        self.out = PasteTable()
        lay.addWidget(self.out, 1)
        self._load()

    def _load(self):
        self.methods = self.project.methods()
        fill_table(self.mtable, ["#", "Name", "Type", "Version", "Created (UTC)"],
                   [[str(m["id"]), m["name"], m["type"],
                     str(m["definition"].get("version", 1))
                     + (f" (replaces #{m['definition']['revision_of']})"
                        if m["definition"].get("revision_of") else ""),
                     m["created_utc"][:19]] for m in self.methods])
        self.res = self.project.results()
        fill_table(self.rtable, ["#", "Name", "Kind", "Method", "Version", "Created (UTC)"],
                   [[str(r["id"]), r["name"], r["kind"], str(r["method_id"] or ""),
                     str(r["data"].get("version", 1))
                     + (f" (replaces #{r['data']['revision_of']})"
                        if r["data"].get("revision_of") else ""),
                     r["created_utc"][:19]] for r in self.res])

    def _show_method(self):
        r = self.mtable.currentRow()
        if 0 <= r < len(self.methods):
            d = dict(self.methods[r]["definition"])
            if d.get("regression"):
                d["regression"] = {k: v for k, v in d["regression"].items()
                                   if k not in ("x", "y", "residuals")}
            if d.get("fitted"):
                d["fitted"] = (f"yes — frozen model, {len(d['fitted']['grid'])} variables "
                               "(applies without refitting)")
            self.detail.setPlainText(json.dumps(d, indent=2, ensure_ascii=False))

    def _show_result(self):
        from spectro.storage.tables import result_table

        r = self.rtable.currentRow()
        if not 0 <= r < len(self.res):
            return
        res = self.res[r]
        headers, rows, summary = result_table(self.project, res["data"])
        fill_table(self.out, headers, rows)
        self.out.export_title = res["name"]
        self.out.export_notes = summary
        self.detail.setPlainText("\n".join([f"{res['name']}  ({res['kind']}, result #{res['id']})"]
                                           + summary))

    def _selected_result(self) -> dict | None:
        r = self.rtable.currentRow()
        if not 0 <= r < len(self.res):
            error(self, "Select a result.")
            return None
        return self.res[r]

    def _edit_result(self):
        res = self._selected_result()
        if res is None:
            return
        dlg = ResultEditDialog(self.win, res["id"])
        self.editor = dlg
        dlg.exec()
        self._load()

    def _archive_result(self):
        from spectro.ui.widgets import ask_reason
        res = self._selected_result()
        if res is None:
            return
        reason = ask_reason(self, f"Archive result #{res['id']} '{res['name']}'?", required=True)
        if reason:
            self.project.archive_result(res["id"], reason)
            self._load()

    EDITORS = {"univariate": "UnivariateDialog", "equations": "EquationsDialog",
               "spectral": "ChemometricsDialog", "amplitude_centering": "ProgressiveDialog",
               "absorption_factor": "ProgressiveDialog"}

    def _selected_method(self) -> dict | None:
        r = self.mtable.currentRow()
        if not 0 <= r < len(self.methods):
            error(self, "Select a method.")
            return None
        return self.methods[r]

    def _edit(self):
        md = self._selected_method()
        if md is None:
            return
        cls = self.EDITORS.get(md["type"])
        if cls is None:
            error(self, f"Methods of type '{md['type']}' cannot be edited.")
            return
        dlg = globals()[cls](self.win)
        try:
            dlg.load_method(md)
        except Exception as exc:
            error(self, f"Could not load the method: {exc}")
            return
        self.editor = dlg
        dlg.exec()
        self._load()

    def _rename(self):
        from spectro.ui.widgets import ask_reason
        md = self._selected_method()
        if md is None:
            return
        name, ok = QInputDialog.getText(self, "Rename method", "New name:", text=md["name"])
        if not ok or not name.strip() or name.strip() == md["name"]:
            return
        reason = ask_reason(self, f"Rename method #{md['id']} to '{name.strip()}'?",
                            required=True)
        if reason:
            self.project.revise_method(md["id"], name.strip(), md["definition"], reason)
            self._load()

    def _archive(self):
        from spectro.ui.widgets import ask_reason
        r = self.mtable.currentRow()
        if r < 0:
            return
        reason = ask_reason(self, "Archive this method?", required=True)
        if reason:
            self.project.archive_method(self.methods[r]["id"], reason)
            self._load()

    def _export_model(self):
        from spectro.storage.modelfile import export_model
        r = self.mtable.currentRow()
        if r < 0:
            error(self, "Select a method.")
            return
        md = self.methods[r]
        path, _ = QFileDialog.getSaveFileName(self, "Export model", f"{md['name']}.spmodel",
                                              "Spectro model (*.spmodel)")
        if not path:
            return
        try:
            doc = export_model(self.project, md)
            Path(path).write_text(json.dumps(doc, indent=1, ensure_ascii=False), encoding="utf-8")
        except Exception as exc:
            error(self, exc)
            return
        self.project.log("EXPORT", f"Exported model '{md['name']}'",
                         {"file": path, "method_id": md["id"], "sha256": doc["sha256"],
                          "embedded_spectra": list(doc["embedded_spectra"])})

    def _import_model(self):
        from spectro.storage.modelfile import import_model
        path, _ = QFileDialog.getOpenFileName(self, "Import model", "",
                                              "Spectro model (*.spmodel);;All files (*.*)")
        if not path:
            return
        try:
            doc = json.loads(Path(path).read_text(encoding="utf-8"))
            import_model(self.project, doc, self.win.current_trial(), path)
        except Exception as exc:
            error(self, exc)
            return
        self._load()

    def _apply(self):
        r = self.mtable.currentRow()
        ids = self.win.selected_ids()
        if r < 0 or not ids:
            error(self, "Select a method here and spectra in the project tree.")
            return
        md = self.methods[r]
        d = md["definition"]
        try:
            spectra = self.project.spectra(ids)
            comps, found, _ = determine(self.project, d, ids)
        except Exception as exc:
            error(self, exc)
            return
        rows = []
        for s, f in zip(spectra, found):
            row = [s.name]
            for c, v in zip(comps, f):
                t = s.concentrations.get(c)
                row += [float(v), "" if t is None else t, "" if not t else 100 * v / t]
            rows.append(row)
        fill_table(self.out, ["Spectrum"] + [h for c in comps for h in
                                             (f"{c} found", f"{c} taken", f"{c} rec %")], rows)
        self.project.save_result(f"{md['name']} applied", "routine",
                                 {"ids": ids, "compounds": comps, "found": found},
                                 self.win.current_trial(), md["id"], ids)
        self._load()


class ResultEditDialog(Base):
    """Edit a saved result. Nothing is overwritten: saving stores a new version
    of the result (the old one is archived and the reason is audited).

    Found values can only change by recalculation — with the result's own
    model, or any saved method — never by typing them in."""

    OWN = -1

    def __init__(self, win, rid: int):
        super().__init__(win, f"Edit result #{rid}")
        self.resize(1000, 680)
        self.res = self.project.result(rid)
        self.data = json.loads(json.dumps(self.res["data"]))
        self.method_id = self.res["method_id"]
        lay = QVBoxLayout(self)
        f = QFormLayout()
        self.name = QLineEdit(self.res["name"])
        f.addRow("Name", self.name)
        self.notes = QPlainTextEdit(self.data.get("notes", ""))
        self.notes.setPlaceholderText("Comments shown with the result and in the Excel export")
        self.notes.setMaximumHeight(70)
        f.addRow("Notes", self.notes)
        lay.addLayout(f)
        self.has_rows = bool(self.data.get("ids")) and isinstance(self.data.get("found"), list)
        g = QGroupBox("Spectra — untick to exclude from the statistics (kept in the table, "
                      "marked 'excluded')")
        gl = QVBoxLayout(g)
        self.table = QTableWidget()
        gl.addWidget(self.table)
        self.excl_note = QLineEdit(self.data.get("exclusion_note", ""))
        self.excl_note.setPlaceholderText("Why are spectra excluded? (e.g. Grubbs outlier, "
                                          "bubble in cuvette)")
        gl.addWidget(self.excl_note)
        self.summary = QLabel()
        self.summary.setWordWrap(True)
        gl.addWidget(self.summary)
        lay.addWidget(g, 1)
        rg = QGroupBox("Recalculate found values with the current data")
        rl = QFormLayout(rg)
        self.source = QComboBox()
        own = self._own_definition()
        if own is not None:
            self.source.addItem("This result's own model", self.OWN)
        current = (self.project.current_method_version(self.method_id)
                   if self.method_id else None)
        for m in self.project.methods():
            if m["type"] in ("univariate", "equations", "spectral") or \
                    m["type"] in uv.PROGRESSIVE_TYPES:
                tag = "  (the method of this result, current version)" \
                    if current and m["id"] == current["id"] else ""
                self.source.addItem(f"Saved method #{m['id']} {m['name']}{tag}", m["id"])
        if current is not None:
            self.source.setCurrentIndex(self.source.findData(current["id"]))
        self.refit = QCheckBox("Refit the calibration on the current calibration spectra "
                               "(use after correcting their concentrations or processing)")
        self.refit.setChecked(bool((own or {}).get("calibration_ids")
                                   or self.data.get("calibration_ids")))
        recalc = QPushButton("Recalculate")
        recalc.clicked.connect(self._recalculate)
        rl.addRow("Method", self.source)
        rl.addRow(self.refit)
        rl.addRow(recalc)
        self.recalc_msg = QLabel()
        self.recalc_msg.setWordWrap(True)
        rl.addRow(self.recalc_msg)
        if not self.has_rows or self.source.count() == 0:
            rg.setEnabled(False)
            self.recalc_msg.setText("This result has no per-spectrum concentrations to "
                                    "recalculate." if not self.has_rows else
                                    "No saved method or stored model to recalculate with.")
        lay.addWidget(rg)
        bb = QDialogButtonBox(QDialogButtonBox.Save | QDialogButtonBox.Cancel)
        bb.accepted.connect(self._save)
        bb.rejected.connect(self.reject)
        lay.addWidget(bb)
        self._fill()
        self.table.itemChanged.connect(self._fill)

    def _own_definition(self) -> dict | None:
        for key in ("model", "method"):
            d = self.data.get(key)
            if isinstance(d, dict) and (d.get("type") in ("univariate", "equations", "spectral")
                                        or d.get("type") in uv.PROGRESSIVE_TYPES):
                d = dict(d)
                if self.data.get("calibration_ids") and not d.get("calibration_ids"):
                    d["calibration_ids"] = self.data["calibration_ids"]
                return d
        return None

    def _excluded(self) -> list[int]:
        out = []
        for i in range(self.table.rowCount()):
            it = self.table.item(i, 0)
            if it is not None and it.checkState() != Qt.Checked:
                out.append(int(it.data(Qt.UserRole)))
        return out

    def _fill(self):
        from spectro.storage.tables import result_table

        if self.table.rowCount():
            self.data["excluded_ids"] = self._excluded()
        headers, rows, summary = result_table(self.project, self.data)
        self.table.blockSignals(True)
        fill_table(self.table, headers, rows)
        self.table.blockSignals(False)
        if self.has_rows and headers and headers[0] == "Spectrum":
            excluded = {int(i) for i in self.data.get("excluded_ids") or []}
            valid = []
            for sid in self.data["ids"]:
                try:
                    self.project.record(int(sid))
                    valid.append(int(sid))
                except (KeyError, TypeError, ValueError):
                    pass
            self.table.blockSignals(True)
            for i, sid in enumerate(valid[:self.table.rowCount()]):
                it = self.table.item(i, 0)
                it.setText(it.text().replace(" (excluded)", ""))
                it.setFlags(it.flags() | Qt.ItemIsUserCheckable)
                it.setData(Qt.UserRole, sid)
                it.setCheckState(Qt.Unchecked if sid in excluded else Qt.Checked)
            self.table.blockSignals(False)
        self.summary.setText("\n".join(summary))

    def _definition(self) -> tuple[dict, int | None]:
        key = self.source.currentData()
        if key == self.OWN:
            return self._own_definition(), self.method_id
        md = self.project.method(int(key))
        return md["definition"], md["id"]

    def _recalculate(self):
        try:
            d, mid = self._definition()
            ids = [int(i) for i in self.data["ids"]]
            comps, found, m = determine(self.project, d, ids, refit=self.refit.isChecked())
        except Exception as exc:
            error(self, exc)
            return
        excluded = self._excluded()
        for k in ("X", "Y", "taken", "recovery", "signals", "names", "found_minus_added"):
            self.data.pop(k, None)
        self.data.update({"ids": ids, "compounds": comps, "found": found,
                          "excluded_ids": excluded})
        if self.data.get("enrichment_added") and len(comps) == 1:
            self.data["found_minus_added"] = [f[0] - self.data["enrichment_added"] for f in found]
        key = "method" if isinstance(self.data.get("method"), dict) else "model"
        self.data[key] = m.to_dict()
        if self.refit.isChecked() and d.get("calibration_ids"):
            self.data["calibration_ids"] = list(d["calibration_ids"])
        self.data["recalculated"] = {"source": self.source.currentText(), "method_id": mid,
                                     "refit": self.refit.isChecked()}
        self.method_id = mid
        self._fill()
        self.recalc_msg.setText(f"Recalculated {len(ids)} spectra with "
                                f"{self.source.currentText()}. Save to keep it.")

    def _save(self):
        from spectro.ui.widgets import ask_reason

        name = self.name.text().strip() or self.res["name"]
        data = self.data
        if self.table.rowCount() and self.has_rows:
            data["excluded_ids"] = self._excluded()
        if not data.get("excluded_ids"):
            data.pop("excluded_ids", None)
        note = self.excl_note.text().strip()
        if note:
            data["exclusion_note"] = note
        else:
            data.pop("exclusion_note", None)
        if self.notes.toPlainText().strip():
            data["notes"] = self.notes.toPlainText().strip()
        else:
            data.pop("notes", None)
        clean = {k: v for k, v in self.res["data"].items() if k not in ("revision_of", "version")}
        if name == self.res["name"] and data == clean:
            self.reject()
            return
        reason = ask_reason(self, f"Save a new version of result #{self.res['id']}?",
                            required=True)
        if reason is None:
            return
        try:
            self.new_id = self.project.revise_result(self.res["id"], name, data, reason,
                                                     self.method_id)
        except Exception as exc:
            error(self, exc)
            return
        self.accept()
