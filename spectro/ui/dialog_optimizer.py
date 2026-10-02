"""Tools → Method optimizer: rank processing strategies for simultaneous
determination and save the chosen ones as methods."""

from __future__ import annotations

import math

import numpy as np

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (QCheckBox, QComboBox, QDoubleSpinBox, QGroupBox, QHBoxLayout, QLabel,
                               QLineEdit, QListWidget, QListWidgetItem, QPushButton, QSplitter,
                               QTableWidget, QTableWidgetItem, QVBoxLayout, QWidget)

from spectro.core.optimizer import FAMILIES, OptimizerInput, explain, materialize, optimize
from spectro.core.spectrum import Spectrum
from spectro.ui.dialogs_data import Base
from spectro.ui.dialogs_methods import CompoundChecks, list_box, scroll
from spectro.ui.widgets import (PasteTable, SpectrumChecklist, SpectrumPlot, error, fill_table,
                                run_with_progress)

def _pure(records, comp, comps):
    return [r for r in records if r.kind == "raw" and r.concentrations.get(comp, 0) > 0
            and all(not r.concentrations.get(o) for o in comps if o != comp)
            and r.role not in ("sample", "unknown", "mixture", "validation")]


class OptimizerDialog(Base):
    def __init__(self, win):
        super().__init__(win, "Method optimizer — best processing for simultaneous determination")
        self.resize(1500, 900)
        self.result = None
        self.inp = None
        self.records = self.project.records()
        comps = self.project.compound_names()
        root = QSplitter()
        left = QWidget()
        ll = QVBoxLayout(left)
        intro = QLabel(
            "Uses the <b>pure standards</b> of each compound to predict how well every "
            "processing strategy separates it from the others (interference, noise and "
            "±λ uncertainty), then checks the best settings of each strategy end-to-end on "
            "simulated mixtures and on your <b>laboratory mixtures</b>.")
        intro.setWordWrap(True)
        ll.addWidget(intro)
        ll.addWidget(QLabel("Compounds in the mixture"))
        self.comps = CompoundChecks(comps)
        self.comps.itemChanged.connect(lambda _: self._refresh_inputs())
        ll.addWidget(self.comps)
        g = QGroupBox("Pure standards and divisor spectra")
        gl = QVBoxLayout(g)
        self.std_table = QTableWidget()
        self.std_table.setMaximumHeight(150)
        gl.addWidget(self.std_table)
        ll.addWidget(g)
        self.mix = SpectrumChecklist()
        self.train = SpectrumChecklist()
        ll.addWidget(list_box("Laboratory mixtures (known concentrations) for verification",
                              self.mix))
        ll.addWidget(list_box("Mixture training set (optional, enables PLS)", self.train))
        g = QGroupBox("Options")
        gl = QVBoxLayout(g)
        row = QHBoxLayout()
        self.range = QLineEdit()
        self.range.setPlaceholderText("auto (excludes the noisy low-UV edge)")
        self.unc = QDoubleSpinBox()
        self.unc.setRange(0, 5)
        self.unc.setSingleStep(0.1)
        self.unc.setValue(0.5)
        self.unc.setSuffix(" nm")
        row.addWidget(QLabel("λ range"))
        row.addWidget(self.range, 1)
        row.addWidget(QLabel("λ uncertainty ±"))
        row.addWidget(self.unc)
        gl.addLayout(row)
        row = QHBoxLayout()
        self.smooth = QCheckBox("Also try Savitzky–Golay smoothing, widths (nm):")
        self.smooth.setChecked(True)
        self.smooth.setToolTip(
            "Every strategy is screened again with Savitzky–Golay smoothing in front (and "
            "derivatives as Savitzky–Golay derivatives). The ranking shows whether smoothing "
            "removes more noise than it distorts the bands. Widths are in nm and converted to "
            "points for your data interval. Screening takes longer.")
        self.smooth_w = QLineEdit("3, 6, 10")
        self.smooth_w.setMaximumWidth(110)
        self.smooth_p = QComboBox()
        self.smooth_p.addItems(["order 2", "order 3"])
        self.smooth.toggled.connect(self.smooth_w.setEnabled)
        self.smooth.toggled.connect(self.smooth_p.setEnabled)
        for w in (self.smooth, self.smooth_w, self.smooth_p):
            row.addWidget(w)
        row.addStretch(1)
        gl.addLayout(row)
        self.fam = QListWidget()
        for key, label in FAMILIES.items():
            it = QListWidgetItem(label)
            it.setData(Qt.UserRole, key)
            it.setFlags(it.flags() | Qt.ItemIsUserCheckable)
            it.setCheckState(Qt.Checked)
            self.fam.addItem(it)
        self.fam.setMaximumHeight(150)
        gl.addWidget(QLabel("Strategies to try"))
        gl.addWidget(self.fam)
        ll.addWidget(g)
        run = QPushButton("Find the best methods")
        run.setStyleSheet("font-weight: bold; padding: 6px;")
        run.clicked.connect(self._run)
        ll.addWidget(run)
        root.addWidget(scroll(left))

        right = QWidget()
        rl = QVBoxLayout(right)
        self.summary = QLabel("Choose the compounds and press <b>Find the best methods</b>.")
        self.summary.setWordWrap(True)
        self.summary.setTextFormat(Qt.RichText)
        rl.addWidget(self.summary)
        row = QHBoxLayout()
        self.show_comp = QComboBox()
        self.show_comp.currentTextChanged.connect(self._fill)
        row.addWidget(QLabel("Ranking for"))
        row.addWidget(self.show_comp)
        row.addStretch(1)
        save = QPushButton("Save selected as method")
        save.clicked.connect(self._save)
        row.addWidget(save)
        rl.addLayout(row)
        split = QSplitter(Qt.Vertical)
        self.table = PasteTable()
        self.table.setSelectionBehavior(QTableWidget.SelectRows)
        self.table.itemSelectionChanged.connect(self._explain)
        self.table.export_title = "Method optimizer ranking"
        split.addWidget(self.table)
        self.plot = SpectrumPlot(y_label="Processed signal of 1 unit")
        split.addWidget(self.plot)
        split.setSizes([430, 380])
        rl.addWidget(split, 1)
        self.expl = QLabel()
        self.expl.setWordWrap(True)
        rl.addWidget(self.expl)
        root.addWidget(right)
        root.setSizes([520, 980])
        QVBoxLayout(self).addWidget(root)
        self._refresh_inputs()

    # ------------------------------------------------------------ inputs
    def _refresh_inputs(self):
        comps = self.comps.checked()
        rows, self.div_boxes, self.std_ids = [], {}, {}
        self.std_table.clear()
        self.std_table.setColumnCount(3)
        self.std_table.setHorizontalHeaderLabels(["Compound", "Pure standards", "Divisor spectrum"])
        self.std_table.setRowCount(len(comps))
        for i, c in enumerate(comps):
            pure = _pure(self.records, c, comps)
            self.std_ids[c] = [r.id for r in pure]
            self.std_table.setItem(i, 0, QTableWidgetItem(c))
            conc = sorted(r.concentrations[c] for r in pure)
            self.std_table.setItem(i, 1, QTableWidgetItem(
                f"{len(pure)} ({conc[0]:g}–{conc[-1]:g})" if pure else "none found"))
            cb = QComboBox()
            for r in sorted(pure, key=lambda r: r.concentrations[c]):
                cb.addItem(f"#{r.id} {r.name}", r.id)
            if pure:
                cb.setCurrentIndex(len(pure) // 2)
            self.std_table.setCellWidget(i, 2, cb)
            self.div_boxes[c] = cb
            rows.append(c)
        self.std_table.resizeColumnsToContents()
        raw = [r for r in self.records if r.kind == "raw"]
        others = [c for c in self.project.compound_names() if c not in comps]
        # verification / training spectra must contain only the selected compounds
        known = [r for r in raw if comps
                 and all(r.concentrations.get(c) is not None for c in comps)
                 and not any(r.concentrations.get(o) for o in others)]
        self.mix.populate([r for r in known if r.role in ("mixture", "validation")],
                          {r.id for r in known if r.role in ("mixture", "validation")})
        self.train.populate([r for r in known if r.role == "calibration"],
                            {r.id for r in known if r.role == "calibration"})

    def _build_input(self) -> OptimizerInput:
        comps = self.comps.checked()
        if len(comps) < 2:
            raise ValueError("choose at least two compounds")
        missing = [c for c in comps if len(self.std_ids.get(c, [])) < 2]
        if missing:
            raise ValueError("need at least two pure standards (only that compound's "
                             f"concentration > 0) for: {', '.join(missing)}")
        rng = None
        t = self.range.text().strip()
        if t:
            a, b = t.replace("–", "-").split("-")
            rng = (float(a), float(b))
        fam = {self.fam.item(i).data(Qt.UserRole) for i in range(self.fam.count())
               if self.fam.item(i).checkState() == Qt.Checked}
        widths: tuple[float, ...] = ()
        if self.smooth.isChecked():
            try:
                widths = tuple(float(w) for w in
                               self.smooth_w.text().replace(";", ",").split(",") if w.strip())
            except ValueError:
                raise ValueError("smoothing widths: give numbers in nm, e.g. 3, 6, 10") from None
            if any(w <= 0 for w in widths):
                raise ValueError("smoothing widths must be positive")
        return OptimizerInput(
            comps, {c: self.project.spectra(self.std_ids[c]) for c in comps},
            {c: self.project.spectrum(self.div_boxes[c].currentData()) for c in comps},
            self.project.spectra(self.mix.checked_ids()),
            self.project.spectra(self.train.checked_ids()),
            wl_range=rng, wl_uncertainty=self.unc.value(), families=fam, smoothing=widths,
            smoothing_order=2 + self.smooth_p.currentIndex())

    # ------------------------------------------------------------ run
    def _run(self):
        try:
            self.inp = self._build_input()
            self.result = run_with_progress(self, "Screening processing strategies…",
                                            optimize, inp=self.inp)
        except InterruptedError:
            return
        except Exception as exc:
            error(self, exc)
            return
        comps = self.inp.compounds
        self.show_comp.blockSignals(True)
        self.show_comp.clear()
        self.show_comp.addItems(comps)
        self.show_comp.blockSignals(False)
        lines = []
        for c in comps:
            ranked = [x for x in self.result["ranked"][c] if not x.error]
            best = ranked[0] if ranked else None
            best_uni = next((x for x in ranked if x.measurement), None)
            if best is None:
                continue
            txt = f"<b>{c}</b>: {best.label} — expected error ≈ {best.score:.2f} %"
            if best_uni is not None and best_uni is not best:
                txt += (f"; best single-signal method: {best_uni.label} "
                        f"(≈ {best_uni.score:.2f} %)")
            lines.append(txt)
        self.summary.setText("<br>".join(["<b>Recommendations</b> (lower is better; "
                                          "conservative estimate)"] + lines))
        self._fill(comps[0])
        self.project.log("CALCULATE", f"Method optimizer for {', '.join(comps)}", {
            "inputs": [i for c in comps for i in self.std_ids[c]] + self.mix.checked_ids(),
            "divisors": {c: self.div_boxes[c].currentData() for c in comps},
            "wl_uncertainty": self.inp.wl_uncertainty,
            "top": {c: [{"label": x.label, "score": x.score} for x in
                        self.result["ranked"][c][:5]] for c in comps}})

    def _fill(self, comp):
        if not self.result or comp not in self.result["ranked"]:
            return
        self.cands = self.result["ranked"][comp]

        def v(x):
            return "" if x is None or (isinstance(x, float) and math.isnan(x)) else (
                "∞" if isinstance(x, float) and math.isinf(x) else x)
        rows = []
        for i, c in enumerate(self.cands, 1):
            rows.append([i, c.label, v(c.score), v(c.robust_error), v(c.sim_rmsep),
                         v(c.real_rmsep), v(c.real_mean_recovery), v(c.real_rsd),
                         v(c.interference), v(c.noise), v(c.r), c.error or ""])
        fill_table(self.table, ["#", "Method and settings", "Score %", "Predicted % (±λ)",
                                "Simulated RMSEP %", "Lab RMSEP %", "Lab recovery %",
                                "Lab RSD %", "Interference %", "Noise %", "r", "Note"], rows)
        self.table.setColumnWidth(1, 470)
        if rows:
            for j in range(self.table.columnCount()):
                it = self.table.item(0, j)
                if it:
                    it.setForeground(Qt.darkGreen)
            self.table.selectRow(0)

    def _explain(self):
        r = self.table.currentRow()
        if not self.result or not 0 <= r < len(self.cands):
            return
        c = self.cands[r]
        self.plot.clear_markers()
        if c.family == "amplitude_centering" and c.model:
            self._explain_centering(c)
            return
        if not c.measurement:
            grid, units = self.result["grid"], self.result["units"]
            self.plot.plot_spectra([Spectrum(grid, u, name=f"{k} (1 unit)")
                                    for k, u in units.items()])
            self.expl.setText("Multicomponent model: uses the whole spectrum / several "
                              "signals of all compounds at once.")
            return
        from spectro.core.optimizer import _Screen
        try:
            proc = explain(c, self.inp, _Screen(self.inp).resolve)
        except Exception as exc:
            self.expl.setText(str(exc))
            return
        order = [c.compound] + [k for k in proc if k != c.compound]
        self.plot.plot_spectra([proc[k] for k in order])
        p = c.measurement["params"]
        for key in ("w1", "w2"):
            if p.get(key) is not None:
                self.plot.set_marker("λ1" if key == "w1" else "λ2", p[key], movable=False)
        self.expl.setText(
            f"Processed spectra of one concentration unit of each compound. At the marked "
            f"wavelength(s) the interferents contribute ≈ {c.interference:.2f} % and noise "
            f"≈ {c.noise:.2f} % of the {c.compound} signal (averaged over the mixtures); "
            f"with a ±{self.inp.wl_uncertainty:g} nm wavelength error the predicted error "
            f"is ≤ {c.robust_error:.2f} %.")

    def _explain_centering(self, c):
        from spectro.core.optimizer import DIV, _Screen
        m = c.model
        grid, units = self.result["grid"], self.result["units"]
        try:
            div = _Screen(self.inp).resolve(DIV + m["divisor_compound"])
            dv = np.interp(grid, div.wavelengths, div.values)
            ok = np.abs(dv) >= 0.05 * np.max(np.abs(dv))
            self.plot.plot_spectra([Spectrum(grid[ok], u[ok] / dv[ok], name=f"{k} ÷ "
                                             f"{m['divisor_compound']}′ (1 unit)")
                                    for k, u in units.items()])
        except Exception as exc:
            self.expl.setText(str(exc))
            return
        self.plot.set_marker("λc", m["wavelength"], movable=False)
        for k, d in m["differences"].items():
            self.plot.set_marker(f"λ2 {k}", d["w2"], movable=False)
        self.expl.setText(
            "Ratio spectra of one unit of each compound. Every compound is read at the common "
            "λc: from the plateau, from an amplitude difference at a λ pair where the other "
            "compounds cancel, or by subtracting the others from the recorded amplitude. "
            f"Predicted error with ±{self.inp.wl_uncertainty:g} nm: {c.robust_error:.2f} % "
            f"(noise {c.noise:.2f} %).")

    # ------------------------------------------------------------ save
    def _save(self):
        r = self.table.currentRow()
        if not self.result or not 0 <= r < len(self.cands):
            error(self, "Select a method in the ranking.")
            return
        c = self.cands[r]
        if c.error:
            error(self, c.error)
            return
        comps = self.inp.compounds
        tid = self.win.current_trial()
        try:
            if c.measurement:
                from spectro.core.univariate import UnivariateMethod
                d = materialize(c, {k: self.div_boxes[k].currentData() for k in comps})
                m = UnivariateMethod(c.label, c.compound, d["steps"], d["measurement"])
                ids = self.std_ids[c.compound]
                m.calibrate(self.project.spectra(ids), self.project.resolver())
                definition = m.to_dict()
            elif c.family in ("amplitude_centering", "absorption_factor"):
                from spectro.core.univariate import progressive_from_dict
                d = materialize(c, {k: self.div_boxes[k].currentData() for k in comps})
                m = progressive_from_dict(d["model"])
                m.name = c.label
                ids = [i for k in comps for i in self.std_ids[k]]
                m.fit_spectra(self.project.spectra(ids), self.project.resolver())
                definition = m.to_dict()
            else:
                definition = dict(c.model)
                ids = ([i for k in comps for i in self.std_ids[k]] if c.family != "pls"
                       else self.train.checked_ids())
            definition["calibration_ids"] = ids
            definition["optimizer"] = {"score": c.score, "robust_error": c.robust_error,
                                       "sim_rmsep": c.sim_rmsep, "real_rmsep": c.real_rmsep}
            mid = self.project.save_method(c.label, definition, tid)
        except Exception as exc:
            error(self, exc)
            return
        self.win.statusBar().showMessage(f"Saved method #{mid}: {c.label}")
        self.summary.setText(self.summary.text() + f"<br>✔ Saved method #{mid}: {c.label}")


__all__ = ["OptimizerDialog"]
