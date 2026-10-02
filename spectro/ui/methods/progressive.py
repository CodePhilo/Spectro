"""Progressive resolution (amplitude centering, absorption factor)."""

from __future__ import annotations

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QFormLayout,
    QHBoxLayout,
    QInputDialog,
    QLabel,
    QPushButton,
    QSpinBox,
    QTableWidget,
    QTableWidgetItem,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)

from spectro.core import univariate as uv
from spectro.ui.dialogs_data import Base, spectrum_choices
from spectro.ui.methods.common import (
    CAL_ROLES,
    TEST_ROLES,
    CompoundChecks,
    check_ids,
    list_box,
    scroll,
    start_editing,
    store_method,
    summary_text,
    wl_spin,
)
from spectro.ui.widgets import PasteTable, SpectrumChecklist, error, fill_table

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
