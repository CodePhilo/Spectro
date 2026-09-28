"""Quantitative method dialogs."""

from __future__ import annotations

import json

import numpy as np
import pyqtgraph as pg
from PySide6.QtCore import Qt
from PySide6.QtWidgets import (QCheckBox, QComboBox, QDoubleSpinBox, QFormLayout, QGroupBox,
                               QHBoxLayout, QInputDialog, QLabel, QLineEdit, QListWidget,
                               QListWidgetItem, QPlainTextEdit, QPushButton, QScrollArea,
                               QSpinBox, QSplitter, QTableWidget, QTabWidget, QVBoxLayout,
                               QWidget)

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
    cal_ids = {r.id for r in recs if r.role in cal_roles and
               (compound is None or r.concentrations.get(compound, 0) > 0)}
    if cal_roles is TRAIN_ROLES and any(r.role == "calibration" for r in recs):
        cal_ids = {r.id for r in recs if r.role == "calibration"}
    test_ids = {r.id for r in recs if r.role in test_roles}
    cal.populate(recs, cal_ids)
    test.populate(recs, test_ids)
    return cal, test


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
        self.template.addItems(list(uv.TEMPLATES))
        self.template.currentTextChanged.connect(self._template)
        self.name = QLineEdit()
        self.origin = QCheckBox("Force through origin")
        f.addRow("Compound", self.compound)
        f.addRow("Method", self.template)
        f.addRow("Name", self.name)
        f.addRow("", self.origin)
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
        ids = {r.id for r in recs if r.role in CAL_ROLES and r.concentrations.get(comp, 0) > 0}
        self.cal.populate(recs, ids)
        self._preview()

    def _template(self, name):
        t = uv.TEMPLATES.get(name)
        if not t:
            return
        self.pipe.set_steps(t["steps"])
        self.meas.set(t["measurement"])
        self.name.setText(f"{name} — {self.compound.currentText()}")
        needs = [s for s in t["steps"] for k, v in s["params"].items() if v is None]
        if needs:
            self.msg.setText("Choose the divisor / reference spectra in the highlighted steps.")

    def _current(self) -> uv.UnivariateMethod:
        return uv.UnivariateMethod(self.name.text().strip() or "Univariate method",
                                   self.compound.currentText().strip(), self.pipe.get_steps(),
                                   self.meas.get(), self.origin.isChecked())

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
                found.append(float(self.method.regression.predict_x(sig)))
                names.append(s.name)
                taken.append(s.concentrations.get(comp))
        except Exception as exc:
            error(self, exc)
            return
        rows, recs = results_rows(names, found, taken)
        for r, sig in zip(rows, signals):
            r.insert(1, float(sig))
        fill_table(self.results, ["Spectrum", "Signal", f"Found ({comp})", "Taken", "Recovery %"],
                   rows)
        self.summary.setText(summary_text(recs))
        self.predictions = {"ids": ids, "names": names, "found": found, "taken": taken,
                            "signals": signals, "recovery": recs}
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
        self.method_id = self.project.save_method(self.method.name, d, self.win.current_trial())
        self.win.statusBar().showMessage("Method saved.")

    def _save_results(self):
        if not self.predictions:
            error(self, "Determine first.")
            return
        self.project.save_result(f"{self.method.name}", "univariate", {
            "method": self.method.to_dict(), **self.predictions},
            self.win.current_trial(), getattr(self, "method_id", None), self.predictions["ids"])
        self.win.statusBar().showMessage("Results saved.")


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
        name, ok = QInputDialog.getText(self, "Save method", "Name:", text="Equation method")
        if ok:
            d = self.model.to_dict()
            d["calibration_ids"] = self.cal_ids
            self.method_id = self.project.save_method(name, d, self.win.current_trial())

    def _save_results(self):
        if not self.predictions:
            error(self, "Determine first.")
            return
        self.project.save_result("Equation method results", "equations",
                                 {"model": self.model.to_dict(), **self.predictions},
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
        super().__init__(win, "Binary special methods")
        self.resize(1250, 800)
        self.comps = self.project.compound_names()
        tabs = QTabWidget()
        tabs.addTab(self._q_tab(), "Absorbance ratio (Q-analysis)")
        tabs.addTab(self._as_tab(), "Absorbance subtraction")
        tabs.addTab(self._am_tab(), "Amplitude modulation")
        tabs.addTab(self._h_tab(), "H-point standard addition")
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
        self.diag_plot = xy_plot("Hotelling T²", "Q residuals")
        self.tabs.addTab(self.diag_plot, "Outliers (T² vs Q)")
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

    def _cv(self):
        try:
            m = self._build()
            ids, spectra = self._cal()
            cv = m.cross_validate(spectra, self.project.resolver(), self.cv.currentText(),
                                  self.folds.value(),
                                  max_components=min(15, len(ids) - 2) if m.model_type in
                                  ("PCR", "PLS1", "PLS2", "ANN") else None)
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
        try:
            m = self._build()
            ids, spectra = self._cal()
            m.fit(spectra, self.project.resolver())
            test_ids = self.test.checked_ids()
            test = self.project.spectra(test_ids)
            pred = m.predict(test, self.project.resolver()) if test else np.zeros((0, len(m.compounds)))
            diag = m.diagnostics(spectra, self.project.resolver())
        except Exception as exc:
            error(self, exc)
            return
        self.model, self.cal_ids = m, ids
        # predicted vs actual
        self.pred_plot.clear()
        rows, lines = [], []
        lo, hi = np.inf, -np.inf
        for j, c in enumerate(m.compounds):
            pts = [(s.concentrations[c], p[j], s.name) for s, p in zip(test, pred)
                   if s.concentrations.get(c) is not None]
            if pts:
                a, b, n = zip(*pts)
                scatter(self.pred_plot, a, b, SERIES[j % 8], c, labels=list(n))
                lo, hi = min(lo, *a, *b), max(hi, *a, *b)
                recs = [100 * y / x for x, y, _ in pts if x]
                from spectro.core.validation import prediction_error
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
        self.summary.setText("\n".join(lines))
        # diagnostics
        self.diag_plot.clear()
        names = [s.name for s in spectra]
        scatter(self.diag_plot, diag["T2"], diag["Q"], SERIES[0], "Calibration", labels=names)
        # auxiliary spectra
        from spectro.core.spectrum import Spectrum
        grid = np.array(diag["grid"])
        aux = []
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
        info.append("Explained variance (PCA of X): " + ", ".join(
            f"{100 * v:.2f}%" for v in diag["explained"][:6]))
        if m.model_type == "MCR-ALS":
            info.append(f"MCR-ALS lack of fit: {m.state['lof']:.4f} %")
        self.info.setPlainText("\n".join(info))
        self.predictions = {"ids": test_ids, "compounds": m.compounds, "found": pred.tolist()}
        self.tabs.setCurrentIndex(1 if rows else 5)
        self.project.log("CALCULATE", f"Fitted {m.model_type} and predicted {len(test_ids)} spectra",
                         {"inputs": ids + test_ids, "model": m.to_dict(),
                          "found": pred.tolist()})

    def _ipls(self):
        try:
            m = self._build()
            ids, spectra = self._cal()
            n, ok = QInputDialog.getInt(self, "iPLS", "Number of intervals:", 10, 2, 50)
            if not ok:
                return
            out = ipls(m, spectra, self.project.resolver(), n)
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
        try:
            m = self._build()
            ids, spectra = self._cal()
            self.summary.setText("Running genetic algorithm…")
            self.repaint()
            out = ga_select(m, spectra, self.project.resolver())
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
        name, ok = QInputDialog.getText(self, "Save method", "Name:", text=f"{self.model.model_type} model")
        if ok:
            d = self.model.to_dict()
            d["calibration_ids"] = self.cal_ids
            self.method_id = self.project.save_method(name, d, self.win.current_trial())

    def _save_results(self):
        if not self.predictions:
            error(self, "Fit + predict first.")
            return
        self.project.save_result(f"{self.model.model_type} results", "chemometrics",
                                 {"model": self.model.to_dict(), **self.predictions},
                                 self.win.current_trial(), getattr(self, "method_id", None),
                                 self.predictions["ids"])


# --------------------------------------------------------------------------- #
# Saved methods & results
# --------------------------------------------------------------------------- #
def method_from_definition(d: dict):
    t = d.get("type")
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
        arch = QPushButton("Archive method…")
        arch.clicked.connect(self._archive)
        row.addWidget(apply_btn)
        row.addWidget(arch)
        row.addStretch(1)
        ml.addLayout(row)
        tabs.addTab(mw, "Methods")
        self.rtable = QTableWidget()
        self.rtable.itemSelectionChanged.connect(self._show_result)
        tabs.addTab(self.rtable, "Results")
        lay.addWidget(tabs, 1)
        self.detail = QPlainTextEdit()
        self.detail.setReadOnly(True)
        lay.addWidget(self.detail, 1)
        self.out = PasteTable()
        lay.addWidget(self.out, 1)
        self._load()

    def _load(self):
        self.methods = self.project.methods()
        fill_table(self.mtable, ["#", "Name", "Type", "Created (UTC)"],
                   [[str(m["id"]), m["name"], m["type"], m["created_utc"][:19]] for m in self.methods])
        self.res = self.project.results()
        fill_table(self.rtable, ["#", "Name", "Kind", "Method", "Created (UTC)"],
                   [[str(r["id"]), r["name"], r["kind"], str(r["method_id"] or ""),
                     r["created_utc"][:19]] for r in self.res])

    def _show_method(self):
        r = self.mtable.currentRow()
        if 0 <= r < len(self.methods):
            d = dict(self.methods[r]["definition"])
            if d.get("regression"):
                d["regression"] = {k: v for k, v in d["regression"].items()
                                   if k not in ("x", "y", "residuals")}
            self.detail.setPlainText(json.dumps(d, indent=2, ensure_ascii=False))

    def _show_result(self):
        r = self.rtable.currentRow()
        if not 0 <= r < len(self.res):
            return
        data = self.res[r]["data"]
        ids, found, comps = data.get("ids"), data.get("found"), data.get("compounds")
        if not (ids and comps and isinstance(found, list) and found
                and isinstance(found[0], list)):
            self.out.setRowCount(0)
            self.detail.setPlainText(json.dumps(data, indent=2, ensure_ascii=False))
            return
        rows, recs = [], {c: [] for c in comps}
        for sid, f in zip(ids, found):
            try:
                rec = self.project.record(int(sid))
            except KeyError:
                continue
            row = [rec.name]
            for c, v in zip(comps, f):
                t = rec.concentrations.get(c)
                rv = 100 * v / t if t else None
                if rv is not None:
                    recs[c].append(rv)
                row += [float(v), "" if t is None else float(t), "" if rv is None else rv]
            rows.append(row)
        fill_table(self.out, ["Spectrum"] + [h for c in comps for h in
                                             (f"{c} found", f"{c} taken", f"{c} rec %")], rows)
        lines = [f"{self.res[r]['name']}  ({self.res[r]['kind']}, {len(rows)} spectra)"]
        lines += [f"{c}: " + summary_text(v) for c, v in recs.items() if len(v) > 1]
        self.detail.setPlainText("\n".join(lines))

    def _archive(self):
        from spectro.ui.widgets import ask_reason
        r = self.mtable.currentRow()
        if r < 0:
            return
        reason = ask_reason(self, "Archive this method?", required=True)
        if reason:
            self.project.archive_method(self.methods[r]["id"], reason)
            self._load()

    def _apply(self):
        r = self.mtable.currentRow()
        ids = self.win.selected_ids()
        if r < 0 or not ids:
            error(self, "Select a method here and spectra in the project tree.")
            return
        md = self.methods[r]
        d = md["definition"]
        res = self.project.resolver()
        try:
            m = method_from_definition(d)
            spectra = self.project.spectra(ids)
            if d["type"] == "univariate":
                if m.regression is None:
                    m.calibrate(self.project.spectra(d["calibration_ids"]), res)
                found = [[m.predict(s, res)] for s in spectra]
                comps = [m.compound]
            elif d["type"] == "equations":
                if m.K is None:
                    m.fit(self.project.spectra(d["calibration_ids"]), res)
                found = m.predict(spectra, res).tolist()
                comps = m.compounds
            else:
                m.fit(self.project.spectra(d["calibration_ids"]), res)
                found = m.predict(spectra, res).tolist()
                comps = m.compounds
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
