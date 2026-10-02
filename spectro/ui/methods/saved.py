"""Saved methods & results, result editor."""

from __future__ import annotations

import json
from pathlib import Path

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QDialogButtonBox,
    QFileDialog,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QInputDialog,
    QLabel,
    QLineEdit,
    QPlainTextEdit,
    QPushButton,
    QTableWidget,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)

from spectro.core import univariate as uv
from spectro.storage.recalc import determine, method_from_definition  # noqa: F401
from spectro.ui.dialogs_data import Base
from spectro.ui.methods.chemometrics import ChemometricsDialog  # noqa: F401
from spectro.ui.methods.equations import EquationsDialog  # noqa: F401
from spectro.ui.methods.progressive import ProgressiveDialog  # noqa: F401
from spectro.ui.methods.univariate import UnivariateDialog  # noqa: F401
from spectro.ui.widgets import PasteTable, error, fill_table


# --------------------------------------------------------------------------- #
# Saved methods & results
# --------------------------------------------------------------------------- #
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
        hist = QPushButton("History…")
        hist.setToolTip("All versions of the method, the differences between them, and "
                        "restore an older version")
        hist.clicked.connect(lambda: self._history("method"))
        for b in (apply_btn, edit, ren, hist, exp, imp, arch):
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
        rhist = QPushButton("History…")
        rhist.clicked.connect(lambda: self._history("result"))
        rrow.addWidget(redit)
        rrow.addWidget(rhist)
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

    def _history(self, kind: str):
        from spectro.ui.dialogs_data import HistoryDialog
        item = self._selected_method() if kind == "method" else self._selected_result()
        if item is None:
            return
        dlg = HistoryDialog(self.win, kind, item["id"])
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
        orow = QHBoxLayout()
        self.out_test = QComboBox()
        self.out_test.addItems(["Grubbs (ISO 5725)", "Dixon Q (n ≤ 10)"])
        self.out_test.currentIndexChanged.connect(self._outliers)
        self.out_msg = QLabel()
        self.out_msg.setWordWrap(True)
        self.out_btn = QPushButton("Exclude flagged")
        self.out_btn.setToolTip("Untick the spectra flagged as outliers and fill in the note "
                                "(you decide — nothing is excluded automatically)")
        self.out_btn.clicked.connect(self._exclude_flagged)
        orow.addWidget(QLabel("Outlier test on recoveries (α = 0.05):"))
        orow.addWidget(self.out_test)
        orow.addWidget(self.out_msg, 1)
        orow.addWidget(self.out_btn)
        gl.addLayout(orow)
        lay.addWidget(g, 1)
        from spectro.storage.recalc import stale_taken
        stale = stale_taken(self.project, self.data)
        if stale:
            warn = QLabel("⚠ Concentrations changed after this result was saved: "
                          + "; ".join(f"{self.project.record(s_).name} {c} {a:g} → "
                                      f"{'not set' if b is None else f'{b:g}'}"
                                      for s_, c, a, b in stale[:6] if a is not None)
                          + (" …" if len(stale) > 6 else "")
                          + ". Press Recalculate (with Refit if a standard changed).")
            warn.setWordWrap(True)
            warn.setStyleSheet("color: #b3261e;")
            lay.insertWidget(1, warn)
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
        self._outliers()

    def _recoveries(self) -> dict[str, list[tuple[int, float]]]:
        """Recovery % per compound for the included spectra: {c: [(id, %)]}."""
        from spectro.storage.recalc import compounds_of, found_rows

        comps, found = compounds_of(self.data), found_rows(self.data)
        if not (comps and found and self.data.get("ids")):
            return {}
        snap = self.data.get("taken_snapshot") or {}
        excluded = set(self._excluded()) if self.table.rowCount() else set()
        out: dict[str, list[tuple[int, float]]] = {c: [] for c in comps}
        for sid, row in zip(self.data["ids"], found):
            sid = int(sid)
            if sid in excluded:
                continue
            try:
                conc = snap.get(str(sid)) or self.project.record(sid).concentrations
            except KeyError:
                continue
            for c, v in zip(comps, row):
                t = conc.get(c)
                if t and v is not None:
                    out[c].append((sid, 100 * float(v) / t))
        return out

    def _outliers(self):
        from spectro.core.validation import dixon, grubbs

        self.flagged = []
        lines = []
        for c, pairs in self._recoveries().items():
            vals = [v for _, v in pairs]
            if len(vals) < 3:
                continue
            try:
                r = (grubbs(vals) if self.out_test.currentIndex() == 0 else dixon(vals))
            except ValueError as exc:
                lines.append(f"{c}: {exc}")
                continue
            stat, crit = ("G", "G_crit") if r["test"] == "Grubbs" else ("Q", "Q_crit")
            if r["outlier"]:
                sid = pairs[r["index"]][0]
                self.flagged.append((sid, c, r))
                lines.append(f"{c}: ⚠ {self.project.record(sid).name} ({r['value']:.2f} %) is an "
                             f"outlier, {stat} = {r[stat]:.3f} > {r[crit]:.3f}")
            else:
                lines.append(f"{c}: no outlier ({stat} = {r[stat]:.3f} ≤ {r[crit]:.3f}, "
                             f"n = {r['n']})")
        self.out_msg.setText("\n".join(lines) or "Needs at least 3 recoveries.")
        self.out_btn.setEnabled(bool(self.flagged))

    def _exclude_flagged(self):
        ids = {sid for sid, _, _ in self.flagged}
        self.table.blockSignals(True)
        for i in range(self.table.rowCount()):
            it = self.table.item(i, 0)
            if it is not None and it.data(Qt.UserRole) in ids:
                it.setCheckState(Qt.Unchecked)
        self.table.blockSignals(False)
        note = "; ".join(f"{self.project.record(sid).name}: {r['test']} outlier for {c} "
                         f"({'G' if r['test'] == 'Grubbs' else 'Q'} = "
                         f"{r['G' if r['test'] == 'Grubbs' else 'Q']:.3f})"
                         for sid, c, r in self.flagged)
        cur = self.excl_note.text().strip()
        self.excl_note.setText(f"{cur}; {note}" if cur else note)
        self._fill()

    def _definition(self) -> tuple[dict, int | None]:
        key = self.source.currentData()
        if key == self.OWN:
            return self._own_definition(), self.method_id
        md = self.project.method(int(key))
        return md["definition"], md["id"]

    def _recalculate(self):
        from spectro.storage.recalc import recalculated_data
        try:
            d, mid = self._definition()
            self.data["excluded_ids"] = self._excluded()
            new, _ = recalculated_data(self.project, self.data, d, self.refit.isChecked())
        except Exception as exc:
            error(self, exc)
            return
        new["recalculated"] = {"source": self.source.currentText(), "method_id": mid,
                               "refit": self.refit.isChecked()}
        self.data = new
        self.method_id = mid
        self._fill()
        self.recalc_msg.setText(f"Recalculated {len(new['ids'])} spectra with "
                                f"{self.source.currentText()}, using the current "
                                "concentrations. Save to keep it.")

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
