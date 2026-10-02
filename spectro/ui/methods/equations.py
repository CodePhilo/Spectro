"""Simultaneous-equation methods (Vierordt, bivariate, AUC)."""

from __future__ import annotations

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QCheckBox,
    QGroupBox,
    QHBoxLayout,
    QInputDialog,
    QLabel,
    QPushButton,
    QSplitter,
    QTableWidget,
    QVBoxLayout,
    QWidget,
)

from spectro.core.multicomponent import SignalEquations, kaiser_selection
from spectro.ui.dialogs_data import Base, name_of, spectrum_choices
from spectro.ui.methods.common import (
    CompoundChecks,
    check_ids,
    checklists,
    list_box,
    scroll,
    start_editing,
    store_method,
    summary_text,
)
from spectro.ui.widgets import PasteTable, PipelineEditor, error, fill_table


# --------------------------------------------------------------------------- #
# Equation methods
# --------------------------------------------------------------------------- #
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
