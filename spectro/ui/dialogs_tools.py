"""Tools: spectral finder, validation statistics, calibration design, greenness."""

from __future__ import annotations

import math

import numpy as np
from PySide6.QtCore import QPointF, QRectF, Qt
from PySide6.QtGui import QBrush, QColor, QFont, QPainter, QPen
from PySide6.QtWidgets import (QCheckBox, QComboBox, QDoubleSpinBox, QFileDialog, QFormLayout,
                               QGridLayout, QHBoxLayout, QLabel, QPushButton, QSpinBox,
                               QSplitter, QTableWidgetItem, QTabWidget, QVBoxLayout, QWidget)

from spectro.core import greenness as gr
from spectro.core import univariate as uv
from spectro.core import validation as val
from spectro.core.multicomponent import design_concentrations
from spectro.ui.dialogs_data import Base
from spectro.ui.dialogs_methods import SERIES, scatter, xy_plot
from spectro.ui.widgets import PasteTable, SpectrumPlot, error, fill_table


def _spin(v=0.0, lo=-1e9, hi=1e9, dec=4, suffix=""):
    s = QDoubleSpinBox()
    s.setRange(lo, hi)
    s.setDecimals(dec)
    s.setValue(v)
    s.setSuffix(suffix)
    return s


# --------------------------------------------------------------------------- #
# Finder
# --------------------------------------------------------------------------- #
class FinderDialog(Base):
    def __init__(self, win):
        super().__init__(win, "Spectral finder")
        self.resize(1200, 720)
        root = QSplitter()
        left = QWidget()
        f = QFormLayout(left)
        ids = win.selected_ids()
        recs = self.project.records()
        self.a, self.b = QComboBox(), QComboBox()
        for r in recs:
            self.a.addItem(f"#{r.id} {r.name}", r.id)
            self.b.addItem(f"#{r.id} {r.name}", r.id)
        if ids:
            self.a.setCurrentIndex(max(0, self.a.findData(ids[0])))
            if len(ids) > 1:
                self.b.setCurrentIndex(max(0, self.b.findData(ids[1])))
        self.norm = QCheckBox("Divide by total concentration (for isoabsorptive points)")
        self.norm.setChecked(True)
        self.lo, self.hi = _spin(200, 0, 1e5, 2, " nm"), _spin(400, 0, 1e5, 2, " nm")
        self.prom = _spin(0.0, 0, 1e6, 5)
        self.width = _spin(10, 0.1, 1e4, 1, " nm")
        self.tol = _spin(2, 0.01, 100, 2, " %")
        f.addRow("Spectrum A", self.a)
        f.addRow("Spectrum B", self.b)
        f.addRow("From", self.lo)
        f.addRow("To", self.hi)
        f.addRow(self.norm)
        f.addRow("Min. prominence (0 = auto)", self.prom)
        f.addRow("Plateau min. width", self.width)
        f.addRow("Plateau flatness tolerance", self.tol)
        for text, fn in (("Zero-crossing points of A", self._zero),
                         ("Isoabsorptive points of A and B", self._iso),
                         ("Maxima and minima of A", self._extrema),
                         ("Plateau regions of A (ratio spectra)", self._plateau),
                         ("λ where A equals its value at 'From' (λ pairs)", self._equal)):
            b = QPushButton(text)
            b.clicked.connect(fn)
            f.addRow(b)
        self.table = PasteTable()
        f.addRow(self.table)
        root.addWidget(left)
        self.plot = SpectrumPlot()
        root.addWidget(self.plot)
        root.setSizes([420, 780])
        QVBoxLayout(self).addWidget(root)
        self.a.currentIndexChanged.connect(self._show)
        self.b.currentIndexChanged.connect(self._show)
        self._show()

    def _spectra(self):
        if self.a.currentData() is None:
            raise ValueError("the project has no spectra yet")
        return self.project.spectrum(self.a.currentData()), self.project.spectrum(self.b.currentData())

    def _show(self, marks: list[float] | None = None):
        if self.a.count() == 0:
            return
        a, b = self._spectra()
        self.plot.plot_spectra([a, b] if a.metadata["id"] != b.metadata["id"] else [a])
        self.plot.clear_markers()
        for i, x in enumerate((marks or [])[:12]):
            self.plot.set_marker(f"{x:.1f}", x, movable=False)

    def _report(self, kind, headers, rows, marks):
        fill_table(self.table, headers, rows)
        self._show(marks)
        self.project.log("CALCULATE", f"Finder: {kind} ({len(rows)} found)",
                         {"inputs": [self.a.currentData(), self.b.currentData()],
                          "range": [self.lo.value(), self.hi.value()], "found": rows})

    def _zero(self):
        try:
            a, _ = self._spectra()
        except ValueError as exc:
            error(self, exc)
            return
        zc = uv.zero_crossings(a, self.lo.value(), self.hi.value())
        self._report("zero-crossings", ["λ (nm)"], [[z] for z in zc], zc)

    def _iso(self):
        try:
            a, b = self._spectra()
        except ValueError as exc:
            error(self, exc)
            return
        ca = sum(a.concentrations.values()) if self.norm.isChecked() else 1.0
        cb = sum(b.concentrations.values()) if self.norm.isChecked() else 1.0
        if not ca or not cb:
            error(self, "Both spectra need concentrations to normalise (or untick the option).")
            return
        pts = [p for p in uv.isoabsorptive_points(a, b, ca, cb)
               if self.lo.value() <= p <= self.hi.value()]
        self._report("isoabsorptive points", ["λ (nm)", "A of A at λ", "A of B at λ"],
                     [[p, a.value_at(p), b.value_at(p)] for p in pts], pts)

    def _extrema(self):
        try:
            a, _ = self._spectra()
        except ValueError as exc:
            error(self, exc)
            return
        x, y = a.region(self.lo.value(), self.hi.value())
        sub = a.copy(wavelengths=x, values=y)
        e = uv.extrema(sub, self.prom.value() or None)
        rows = [["max", w, v] for w, v in e["maxima"]] + [["min", w, v] for w, v in e["minima"]]
        rows.sort(key=lambda r: r[1])
        self._report("extrema", ["Type", "λ (nm)", "Value"], rows, [r[1] for r in rows])

    def _equal(self):
        try:
            a, _ = self._spectra()
        except ValueError as exc:
            error(self, exc)
            return
        ref = self.lo.value()
        try:
            pts = uv.equal_amplitude_wavelengths(a, ref)
        except ValueError as exc:
            error(self, exc)
            return
        self._report(f"equal amplitude to {ref:g} nm", ["λ (nm)", f"Value (= value at {ref:g} nm)"],
                     [[p, a.value_at(p)] for p in pts], [ref] + pts)

    def _plateau(self):
        try:
            a, _ = self._spectra()
        except ValueError as exc:
            error(self, exc)
            return
        x, y = a.region(self.lo.value(), self.hi.value())
        sub = a.copy(wavelengths=x, values=y)
        pl = uv.find_plateaus(sub, self.width.value(), self.tol.value() / 100)
        self._report("plateaus", ["From (nm)", "To (nm)", "Mean value"],
                     [list(p) for p in pl], [v for p in pl for v in p[:2]])


# --------------------------------------------------------------------------- #
# Validation
# --------------------------------------------------------------------------- #
class ValidationDialog(Base):
    def __init__(self, win):
        super().__init__(win, "Validation statistics (ICH Q2(R2))")
        self.resize(1250, 780)
        self.last = None
        lay = QVBoxLayout(self)
        lay.addWidget(QLabel("Paste data from Excel (Ctrl+V) into the input table. Results can "
                             "be copied (select + Ctrl+C) and saved to the project."))
        self.tabs = QTabWidget()
        self.inputs: dict[str, PasteTable] = {}
        specs = {
            "Linearity": ["Concentration", "Response"],
            "Accuracy": ["Taken", "Found"],
            "Precision": ["Level / day 1", "Level / day 2", "Level / day 3"],
            "Compare methods": ["Proposed method", "Reference method"],
            "Standard addition": ["Added", "Response"],
        }
        for name, cols in specs.items():
            w = QWidget()
            h = QHBoxLayout(w)
            t = PasteTable(25, len(cols))
            t.setHorizontalHeaderLabels(cols)
            self.inputs[name] = t
            v = QVBoxLayout()
            v.addWidget(t)
            row = QHBoxLayout()
            if name == "Precision":
                add = QPushButton("Add column")
                add.clicked.connect(lambda _=False, t=t: self._add_col(t))
                row.addWidget(add)
            if name == "Linearity":
                self.origin = QCheckBox("Through origin")
                self.lod = QComboBox()
                self.lod.addItems(["σ = SD of intercept", "σ = residual SD (Sy/x)"])
                row.addWidget(self.origin)
                row.addWidget(self.lod)
            if name == "Compare methods":
                self.theta = _spin(2.0, 0.1, 50, 1, " %")
                row.addWidget(QLabel("Interval θ"))
                row.addWidget(self.theta)
            calc = QPushButton("Calculate")
            calc.clicked.connect(lambda _=False, n=name: self._calc(n))
            row.addWidget(calc)
            row.addStretch(1)
            v.addLayout(row)
            h.addLayout(v, 1)
            self.tabs.addTab(w, name)
        lay.addWidget(self.tabs, 1)
        split = QSplitter()
        self.out = PasteTable()
        self.plot = xy_plot("x", "y")
        split.addWidget(self.out)
        split.addWidget(self.plot)
        lay.addWidget(split, 1)
        row = QHBoxLayout()
        save = QPushButton("Save to project")
        save.clicked.connect(self._save)
        save.setEnabled(win.project is not None)
        row.addStretch(1)
        row.addWidget(save)
        lay.addLayout(row)

    def _add_col(self, t: PasteTable):
        t.setColumnCount(t.columnCount() + 1)
        t.setHorizontalHeaderItem(t.columnCount() - 1,
                                  QTableWidgetItem(f"Level / day {t.columnCount()}"))

    def _pairs(self, t: PasteTable):
        a, b = [], []
        for i in range(t.rowCount()):
            x = t.item(i, 0).text().strip().replace(",", ".") if t.item(i, 0) else ""
            y = t.item(i, 1).text().strip().replace(",", ".") if t.item(i, 1) else ""
            if x and y:
                a.append(float(x))
                b.append(float(y))
        return a, b

    def _calc(self, name):
        t = self.inputs[name]
        self.plot.clear()
        try:
            if name == "Linearity":
                x, y = self._pairs(t)
                r = val.linear_regression(x, y, self.origin.isChecked(),
                                          lod_basis="intercept" if self.lod.currentIndex() == 0
                                          else "residual")
                lof = val.lack_of_fit(x, y)
                rows = [["Slope", r.slope], ["SD slope", r.sd_slope], ["Intercept", r.intercept],
                        ["SD intercept", r.sd_intercept], ["r", r.r], ["r²", r.r2],
                        ["Sy/x", r.sy_x], ["n", str(r.n)],
                        ["95% CI slope", f"{r.ci_slope[0]:.5g} – {r.ci_slope[1]:.5g}"],
                        ["95% CI intercept", f"{r.ci_intercept[0]:.5g} – {r.ci_intercept[1]:.5g}"],
                        ["LOD", r.lod], ["LOQ", r.loq]]
                if lof["available"]:
                    rows += [["Lack-of-fit F", lof["F"]], ["Lack-of-fit p", lof["p"]]]
                xs = np.linspace(min(x), max(x), 2)
                scatter(self.plot, x, y, SERIES[0], "data", (xs, r.predict_y(xs)))
                data = r.to_dict()
            elif name == "Accuracy":
                taken, found = self._pairs(t)
                d = val.recovery(found, taken)
                rows = [[f"{a:g} → {b:g}", rec] for a, b, rec in zip(taken, found, d["recoveries"])]
                rows += [["Mean recovery %", d["mean"]], ["SD", d["sd"]], ["RSD %", d["rsd"]],
                         ["n", str(d["n"])]]
                data = d
            elif name == "Precision":
                groups = {}
                for j in range(t.columnCount()):
                    vals = t.column_floats(j)
                    if vals:
                        head = t.horizontalHeaderItem(j)
                        groups[head.text() if head else f"Group {j + 1}"] = vals
                d = val.precision(groups)
                rows = [[g, s["n"], s["mean"], s["sd"], s["rsd"]] for g, s in d.items()]
                allv = [v for g in groups.values() for v in g]
                tot = val.describe(allv)
                rows.append(["All", tot["n"], tot["mean"], tot["sd"], tot["rsd"]])
                headers = ["Group", "n", "Mean", "SD", "RSD %"]
                data = {"groups": d, "overall": tot}
                if len(groups) > 1:
                    an = val.anova_oneway(groups)
                    data["anova"] = an
                    rows.append(["ANOVA F (p)", "", an["F"], an["p"],
                                 "significant" if an["significant"] else "not significant"])
                fill_table(self.out, headers, rows)
                self._done(name, data)
                return
            elif name == "Compare methods":
                a = t.column_floats(0)
                b = t.column_floats(1)
                tt, ff = val.t_test(a, b), val.f_test(a, b)
                ih = val.interval_hypothesis(a, b, self.theta.value() / 100)
                da, db = val.describe(a), val.describe(b)
                rows = [["Proposed: mean ± SD", f"{da['mean']:.4f} ± {da['sd']:.4f} (n={da['n']})"],
                        ["Reference: mean ± SD", f"{db['mean']:.4f} ± {db['sd']:.4f} (n={db['n']})"],
                        ["Student's t (tabulated)", f"{tt['t']:.4f} ({tt['t_crit']:.4f})"],
                        ["t-test p", tt["p"]],
                        ["F (tabulated)", f"{ff['F']:.4f} ({ff['F_crit']:.4f})"],
                        ["F-test p", ff["p"]],
                        ["Interval hypothesis θL – θU",
                         f"{100 * ih['lower']:.2f} % – {100 * ih['upper']:.2f} %"],
                        ["Conclusion", "no significant difference" if not tt["significant"] and
                         not ff["significant"] else "significant difference"]]
                data = {"t": tt, "F": ff, "interval": ih}
            else:
                x, y = self._pairs(t)
                d = val.standard_addition(x, y)
                reg = d["regression"]
                rows = [["Concentration in sample", d["found"]], ["Slope", reg["slope"]],
                        ["Intercept", reg["intercept"]], ["r", reg["r"]]]
                xs = np.array([-d["found"], max(x)])
                scatter(self.plot, x, y, SERIES[0], "data",
                        (xs, reg["slope"] * xs + reg["intercept"]))
                data = d
        except Exception as exc:
            error(self, exc)
            return
        fill_table(self.out, ["Parameter", "Value"], rows)
        self._done(name, data)

    def _done(self, name, data):
        self.last = (name, data)
        if self.project:
            self.project.log("CALCULATE", f"Validation: {name}", {"result": data})

    def _save(self):
        if not self.last or not self.project:
            error(self, "Calculate first.")
            return
        self.project.save_result(f"Validation — {self.last[0]}", "validation", self.last[1],
                                 self.win.current_trial())


# --------------------------------------------------------------------------- #
# Design
# --------------------------------------------------------------------------- #
class DesignDialog(Base):
    def __init__(self, win):
        super().__init__(win, "Calibration design (Brereton multilevel multifactor, 5 levels)")
        self.resize(1000, 720)
        lay = QVBoxLayout(self)
        lay.addWidget(QLabel("25 mixtures with 5 concentration levels per compound and "
                             "uncorrelated (orthogonal) compound concentrations — the standard "
                             "training set for CLS/PCR/PLS/ANN. Level = centre + k·step, k ∈ "
                             "{−2, −1, 0, 1, 2}. Up to 6 compounds."))
        lay.itemAt(0).widget().setWordWrap(True)
        comps = self.project.compound_names() if self.project else ["A", "B"]
        self.inp = PasteTable()
        fill_table(self.inp, ["Compound", "Centre conc.", "Step", "Stock conc."],
                   [[c, "10", "2", "100"] for c in comps[:6]], editable=True)
        self.inp.setMaximumHeight(200)
        lay.addWidget(self.inp)
        row = QHBoxLayout()
        self.vol = _spin(10, 0.1, 1e4, 2, " mL")
        add = QPushButton("Add compound")
        add.clicked.connect(lambda: self.inp.setRowCount(self.inp.rowCount() + 1))
        gen = QPushButton("Generate")
        gen.clicked.connect(self._gen)
        exp = QPushButton("Export to Excel…")
        exp.clicked.connect(self._export)
        row.addWidget(add)
        row.addWidget(QLabel("Final volume"))
        row.addWidget(self.vol)
        row.addWidget(gen)
        row.addWidget(exp)
        row.addStretch(1)
        lay.addLayout(row)
        self.out = PasteTable()
        lay.addWidget(self.out, 1)
        self.rows = None

    def _gen(self):
        try:
            names, centres, steps, stocks = [], {}, {}, {}
            for i in range(self.inp.rowCount()):
                it = self.inp.item(i, 0)
                if not it or not it.text().strip():
                    continue
                n = it.text().strip()
                names.append(n)
                centres[n] = float(self.inp.item(i, 1).text())
                steps[n] = float(self.inp.item(i, 2).text())
                st = self.inp.item(i, 3)
                stocks[n] = float(st.text()) if st and st.text().strip() else None
            design = design_concentrations(centres, steps)
        except Exception as exc:
            error(self, exc)
            return
        if any(v < 0 for d in design for v in d.values()):
            error(self, "Some levels are negative — increase the centre or reduce the step.")
        headers = ["Mixture"] + [f"{n} conc." for n in names] + \
                  [f"{n} stock (mL)" for n in names if stocks[n]]
        rows = []
        for i, d in enumerate(design):
            row = [f"M{i + 1}"] + [d[n] for n in names]
            row += [d[n] * self.vol.value() / stocks[n] for n in names if stocks[n]]
            rows.append(row)
        fill_table(self.out, headers, rows)
        self.rows = (headers, rows)
        if self.project:
            self.project.log("CALCULATE", f"Generated calibration design for {', '.join(names)}",
                             {"centres": centres, "steps": steps, "design": design})

    def _export(self):
        if not self.rows:
            self._gen()
        if not self.rows:
            return
        path, _ = QFileDialog.getSaveFileName(self, "Export design", "design.xlsx", "Excel (*.xlsx)")
        if not path:
            return
        import openpyxl
        wb = openpyxl.Workbook()
        ws = wb.active
        ws.append(self.rows[0])
        for r in self.rows[1]:
            ws.append(r)
        wb.save(path)
        if self.project:
            self.project.log("EXPORT", "Exported calibration design", {"file": path})


# --------------------------------------------------------------------------- #
# Greenness
# --------------------------------------------------------------------------- #
class AgreeChart(QWidget):
    """AGREE-style clock: 12 segments coloured by score, overall score centre."""

    def __init__(self):
        super().__init__()
        self.scores = [1.0] * 12
        self.weights = [2.0] * 12
        self.total = 1.0
        self.setMinimumSize(320, 320)

    def set(self, scores, weights, total):
        self.scores, self.weights, self.total = scores, weights, total
        self.update()

    def paintEvent(self, e):
        p = QPainter(self)
        p.setRenderHint(QPainter.Antialiasing)
        side = min(self.width(), self.height()) - 20
        cx, cy = self.width() / 2, self.height() / 2
        r_out = side / 2
        for i, (s, w) in enumerate(zip(self.scores, self.weights)):
            r = r_out * (0.62 + 0.38 * w / 4)
            rect = QRectF(cx - r, cy - r, 2 * r, 2 * r)
            p.setBrush(QBrush(QColor(gr.score_colour(s))))
            p.setPen(QPen(QColor("#fcfcfb"), 2))
            start = int((90 - i * 30) * 16)
            p.drawPie(rect, start, -30 * 16)
            ang = math.radians(90 - i * 30 - 15)
            p.setPen(QColor("#0b0b0b"))
            p.drawText(QPointF(cx + 0.8 * r * math.cos(ang) - 6, cy - 0.8 * r * math.sin(ang) + 5),
                       str(i + 1))
        rc = r_out * 0.45
        p.setBrush(QBrush(QColor(gr.score_colour(self.total))))
        p.setPen(QPen(QColor("#fcfcfb"), 2))
        p.drawEllipse(QPointF(cx, cy), rc, rc)
        p.setPen(QColor("#0b0b0b"))
        f = QFont()
        f.setPointSize(22)
        f.setBold(True)
        p.setFont(f)
        p.drawText(QRectF(cx - rc, cy - rc, 2 * rc, 2 * rc), Qt.AlignCenter, f"{self.total:.2f}")


class GapiChart(QWidget):
    def __init__(self):
        super().__init__()
        self.colours = ["green"] * 15
        self.setMinimumSize(360, 300)

    def set(self, colours):
        self.colours = colours
        self.update()

    def paintEvent(self, e):
        p = QPainter(self)
        p.setRenderHint(QPainter.Antialiasing)
        # 15 fields in the five GAPI groups (sample prep 1–4 & 6–9, method 5, reagents 10–12,
        # instrumentation 13–15) drawn as labelled tiles
        groups = [("Sampling", [0, 1, 2, 3]), ("Method", [4]), ("Preparation", [5, 6, 7]),
                  ("Reagents", [8, 9, 10]), ("Instrument", [11, 12, 13, 14])]
        y = 10
        for title, idx in groups:
            p.setPen(QColor("#52514e"))
            p.drawText(10, y + 18, title)
            x = 110
            for i in idx:
                p.setBrush(QColor(gr.GAPI_COLOURS[self.colours[i]]))
                p.setPen(QPen(QColor("#fcfcfb"), 2))
                p.drawRoundedRect(QRectF(x, y, 44, 26), 4, 4)
                p.setPen(QColor("#0b0b0b"))
                p.drawText(QRectF(x, y, 44, 26), Qt.AlignCenter, str(i + 1))
                x += 50
            y += 36


class GreennessDialog(Base):
    def __init__(self, win):
        super().__init__(win, "Greenness assessment")
        self.resize(1100, 720)
        self.last = None
        tabs = QTabWidget()
        # AGREE
        w = QWidget()
        h = QHBoxLayout(w)
        g = QGridLayout()
        g.addWidget(QLabel("<b>Principle</b>"), 0, 0)
        g.addWidget(QLabel("<b>Score 0–1</b>"), 0, 1)
        g.addWidget(QLabel("<b>Weight 1–4</b>"), 0, 2)
        self.a_scores, self.a_weights = [], []
        for i, text in enumerate(gr.AGREE_PRINCIPLES):
            g.addWidget(QLabel(text), i + 1, 0)
            s = _spin(1.0, 0, 1, 2)
            s.setSingleStep(0.05)
            wt = QSpinBox()
            wt.setRange(1, 4)
            wt.setValue(2)
            s.valueChanged.connect(self._agree)
            wt.valueChanged.connect(self._agree)
            g.addWidget(s, i + 1, 1)
            g.addWidget(wt, i + 1, 2)
            self.a_scores.append(s)
            self.a_weights.append(wt)
        g.setRowStretch(len(gr.AGREE_PRINCIPLES) + 1, 1)
        gw = QWidget()
        gw.setLayout(g)
        h.addWidget(gw)
        self.agree_chart = AgreeChart()
        h.addWidget(self.agree_chart, 1)
        tabs.addTab(w, "AGREE")
        # Eco-scale
        w = QWidget()
        v = QVBoxLayout(w)
        v.addWidget(QLabel("Penalty points: reagents = amount points (<10 mL: 1, 10–100: 2, "
                           ">100: 3) × hazard (pictograms × 1 for 'warning', × 2 for 'danger'); "
                           "energy ≤0.1 kWh: 0, ≤1.5: 1, >1.5: 2; occupational hazard 0/3; "
                           "waste 1–5 + treatment."))
        v.itemAt(0).widget().setWordWrap(True)
        self.eco = PasteTable()
        fill_table(self.eco, ["Item", "Penalty points"],
                   [["Methanol (reagent)", "6"], ["Instrument energy", "0"],
                    ["Occupational hazard", "0"], ["Waste", "3"]], editable=True)
        v.addWidget(self.eco)
        row = QHBoxLayout()
        add = QPushButton("Add row")
        add.clicked.connect(lambda: self.eco.setRowCount(self.eco.rowCount() + 1))
        calc = QPushButton("Calculate")
        calc.clicked.connect(self._eco)
        self.eco_out = QLabel()
        row.addWidget(add)
        row.addWidget(calc)
        row.addWidget(self.eco_out, 1)
        v.addLayout(row)
        tabs.addTab(w, "Analytical Eco-Scale")
        # GAPI
        w = QWidget()
        h = QHBoxLayout(w)
        f = QFormLayout()
        self.gapi = []
        for i, n in enumerate(gr.GAPI_FIELDS):
            c = QComboBox()
            c.addItems(["green", "yellow", "red"])
            c.currentIndexChanged.connect(self._gapi)
            f.addRow(f"{i + 1}. {n}", c)
            self.gapi.append(c)
        fw = QWidget()
        fw.setLayout(f)
        h.addWidget(fw)
        self.gapi_chart = GapiChart()
        h.addWidget(self.gapi_chart, 1)
        tabs.addTab(w, "GAPI")
        lay = QVBoxLayout(self)
        lay.addWidget(tabs)
        row = QHBoxLayout()
        save = QPushButton("Save to project")
        save.clicked.connect(self._save)
        save.setEnabled(win.project is not None)
        row.addStretch(1)
        row.addWidget(save)
        lay.addLayout(row)
        self._agree()

    def _agree(self):
        scores = [s.value() for s in self.a_scores]
        weights = [float(w.value()) for w in self.a_weights]
        r = gr.agree(scores, weights)
        self.agree_chart.set(scores, weights, r["score"])
        self.last = ("AGREE", r)

    def _eco(self):
        pens = {}
        try:
            for i in range(self.eco.rowCount()):
                a, b = self.eco.item(i, 0), self.eco.item(i, 1)
                if a and b and a.text().strip() and b.text().strip():
                    pens[a.text().strip()] = float(b.text().replace(",", "."))
        except ValueError as exc:
            error(self, exc)
            return
        r = gr.eco_scale(pens)
        self.eco_out.setText(f"<b>Eco-Scale score = {r['score']:g}</b> ({r['rating']})")
        self.last = ("Eco-Scale", r)

    def _gapi(self):
        cols = [c.currentText() for c in self.gapi]
        self.gapi_chart.set(cols)
        self.last = ("GAPI", dict(zip(gr.GAPI_FIELDS, cols)))

    def _save(self):
        if self.last and self.project:
            self.project.save_result(f"Greenness — {self.last[0]}", "greenness", self.last[1],
                                     self.win.current_trial())
