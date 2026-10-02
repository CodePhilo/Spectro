"""Binary two-signal methods (Q-analysis, AS, AM, IAM, AAS, HPSAM)."""

from __future__ import annotations

from PySide6.QtWidgets import (
    QComboBox,
    QDoubleSpinBox,
    QFormLayout,
    QHBoxLayout,
    QInputDialog,
    QLabel,
    QPushButton,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)

from spectro.core import univariate as uv
from spectro.ui.dialogs_data import Base, spectrum_choices
from spectro.ui.methods.common import (
    list_box,
    summary_text,
    wl_spin,
)
from spectro.ui.widgets import PasteTable, SpectrumChecklist, error, fill_table

# --------------------------------------------------------------------------- #
# Special binary methods
# --------------------------------------------------------------------------- #


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
