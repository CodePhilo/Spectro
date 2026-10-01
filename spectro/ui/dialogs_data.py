"""Dialogs for data management: import, trials, compounds, concentrations,
processing pipelines, audit trail and export."""

from __future__ import annotations

import csv
import json
import re
from pathlib import Path

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (QAbstractItemView, QComboBox, QDialog,
                               QDialogButtonBox, QFileDialog, QFormLayout, QHBoxLayout,
                               QInputDialog, QLabel, QLineEdit, QPlainTextEdit,
                               QPushButton, QSplitter, QTableWidget, QTableWidgetItem,
                               QVBoxLayout, QWidget)

from spectro.core.io import (IMPORT_FILTER, ParseOptions, export_csv, export_excel, load_file,
                             write_jcamp, write_spc)
from spectro.core.naming import concentrations_from_name
from spectro.core.operations import apply_pipeline
from spectro.storage.project import ROLES
from spectro.ui.widgets import (PasteTable, PipelineEditor, SpectrumPlot, ask_reason, error,
                                fill_table)


def spectrum_choices(project):
    return lambda: [(r.id, f"#{r.id}  {r.name}") for r in project.records()]


def name_of(project):
    def f(ref):
        try:
            return project.record(int(ref)).name
        except Exception:
            return str(ref)
    return f


class Base(QDialog):
    def __init__(self, win, title: str):
        super().__init__(win)
        self.win = win
        self.project = win.project
        self.setWindowTitle(title)
        self.setWindowFlag(Qt.WindowMaximizeButtonHint, True)


# --------------------------------------------------------------------------- #
# Import
# --------------------------------------------------------------------------- #
class ImportDialog(Base):
    def __init__(self, win):
        super().__init__(win, "Import spectra")
        self.resize(1150, 700)
        self.files: list[str] = []
        self.results: dict[str, object] = {}
        lay = QVBoxLayout(self)
        top = QHBoxLayout()
        add = QPushButton("Add files…")
        add.clicked.connect(self._add_files)
        clear = QPushButton("Clear")
        clear.clicked.connect(self._clear)
        self.layout_cb = QComboBox()
        self.layout_cb.addItems(["auto", "columns", "xy_pairs", "rows"])
        self.layout_cb.setToolTip("columns: λ + one column per spectrum\n"
                                  "xy_pairs: (λ, A) column pairs\nrows: λ in first row, "
                                  "one spectrum per row")
        self.dec_cb = QComboBox()
        self.dec_cb.addItems(["auto", "point (.)", "comma (,)"])
        for w in (self.layout_cb, self.dec_cb):
            w.currentIndexChanged.connect(self._reparse)
        self.trial_cb = QComboBox()
        self._fill_trials(win.current_trial())
        new_trial = QPushButton("New trial…")
        new_trial.clicked.connect(self._new_trial)
        self.role_cb = QComboBox()
        self.role_cb.addItems(ROLES)
        self.role_cb.setEditable(True)
        self.role_cb.currentTextChanged.connect(self._default_role)
        for w in (add, clear, QLabel("Layout"), self.layout_cb, QLabel("Decimal"), self.dec_cb,
                  QLabel("Trial"), self.trial_cb, new_trial, QLabel("Role"), self.role_cb):
            top.addWidget(w)
        top.addStretch(1)
        lay.addLayout(top)
        info = QLabel("Files from any spectrophotometer: CSV/TXT/ASC/PRN/DAT (any delimiter, "
                      "decimal comma, header lines), Excel (.xlsx/.xls), JCAMP-DX, GRAMS SPC, "
                      "and OLE-container .spc. Layout is detected automatically. You can type "
                      "or paste concentrations below, or fill them from the names.")
        info.setWordWrap(True)
        lay.addWidget(info)
        split = QSplitter(Qt.Vertical)
        self.table = PasteTable()
        self.table.itemSelectionChanged.connect(self._preview)
        split.addWidget(self.table)
        self.plot = SpectrumPlot()
        split.addWidget(self.plot)
        split.setSizes([380, 260])
        lay.addWidget(split, 1)
        row = QHBoxLayout()
        parse = QPushButton("Fill concentrations from names")
        parse.setToolTip("'PAR 10 CAF 5' → by compound name; otherwise numbers in the name "
                         "are assigned to the compounds in order.")
        parse.clicked.connect(self._parse_names)
        row.addWidget(parse)
        self.status = QLabel()
        row.addWidget(self.status, 1)
        lay.addLayout(row)
        bb = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        bb.button(QDialogButtonBox.Ok).setText("Import")
        bb.accepted.connect(self._import)
        bb.rejected.connect(self.reject)
        lay.addWidget(bb)
        self.compounds = self.project.compound_names()

    def _fill_trials(self, select=None):
        self.trial_cb.clear()
        for t in self.project.trials():
            self.trial_cb.addItem(t["name"], t["id"])
        self.trial_cb.addItem("(unassigned)", None)
        if select is not None:
            i = self.trial_cb.findData(select)
            if i >= 0:
                self.trial_cb.setCurrentIndex(i)

    def _new_trial(self):
        name, ok = QInputDialog.getText(self, "New trial", "Trial name:")
        if ok and name.strip():
            tid = self.project.add_trial(name.strip())
            self._fill_trials(tid)

    def _options(self) -> ParseOptions:
        dec = {0: None, 1: False, 2: True}[self.dec_cb.currentIndex()]
        return ParseOptions(layout=self.layout_cb.currentText(), decimal_comma=dec)

    def _add_files(self):
        files, _ = QFileDialog.getOpenFileNames(self, "Import spectra", "", IMPORT_FILTER)
        for f in files:
            if f not in self.files:
                self.files.append(f)
        self._reparse()

    def _clear(self):
        self.files, self.results = [], {}
        self._reparse()

    def _reparse(self):
        self.results = {}
        problems = []
        for f in self.files:
            try:
                self.results[f] = load_file(f, self._options())
            except Exception as exc:
                problems.append(f"{Path(f).name}: {exc}")
        from spectro.core.io.excel import role_from_sheet
        headers = ["Import", "File", "Name", "Points", "Range (nm)", "Role"] + self.compounds
        self.table.clear()
        self.table.setColumnCount(len(headers))
        self.table.setHorizontalHeaderLabels(headers)
        rows = [(f, i, s) for f, r in self.results.items() for i, s in enumerate(r.spectra)]
        self.table.setRowCount(len(rows))
        self.rows = rows
        for k, (f, i, s) in enumerate(rows):
            chk = QTableWidgetItem()
            chk.setFlags(Qt.ItemIsUserCheckable | Qt.ItemIsEnabled | Qt.ItemIsSelectable)
            chk.setCheckState(Qt.Checked)
            self.table.setItem(k, 0, chk)
            for j, text in ((1, Path(f).name), (3, str(s.x.size)),
                            (4, f"{s.x[0]:g}–{s.x[-1]:g}")):
                it = QTableWidgetItem(text)
                it.setFlags(it.flags() & ~Qt.ItemIsEditable)
                self.table.setItem(k, j, it)
            self.table.setItem(k, 2, QTableWidgetItem(s.name))
            sheet_role = role_from_sheet(s.metadata.get("sheet"))
            it = QTableWidgetItem(sheet_role or self.role_cb.currentText())
            it.setData(Qt.UserRole, bool(sheet_role))     # fixed by the sheet name
            it.setToolTip("Role of this spectrum. Excel sheets named Standards, Mixtures, "
                          "Divisor, Samples… set it automatically; otherwise the Role chosen "
                          "above is used. Editable.")
            self.table.setItem(k, 5, it)
            for j, c in enumerate(self.compounds):
                self.table.setItem(k, 6 + j, QTableWidgetItem(""))
        self.table.resizeColumnsToContents()
        layouts = {r.layout for r in self.results.values()}
        msg = f"{len(rows)} spectra in {len(self.results)} files (layout: {', '.join(layouts) or '–'})."
        if problems:
            msg += "  Problems: " + "; ".join(problems)
        self.status.setText(msg)
        if rows:
            self.table.selectRow(0)

    def _default_role(self, role):
        """The Role above applies to every spectrum whose role is not set by its sheet."""
        for k in range(self.table.rowCount()):
            it = self.table.item(k, 5)
            if it is not None and not it.data(Qt.UserRole):
                it.setText(role)

    def _preview(self):
        rows = sorted({i.row() for i in self.table.selectedIndexes()}) or [0]
        spectra = [self.rows[r][2] for r in rows[:50] if r < len(self.rows)]
        self.plot.plot_spectra(spectra)

    def _parse_names(self):
        if not self.compounds:
            error(self, "Define compounds first (Edit → Compounds…).")
            return
        for k in range(self.table.rowCount()):
            found = concentrations_from_name(self.table.item(k, 2).text(), self.compounds)
            for j, c in enumerate(self.compounds):
                if c in found:
                    self.table.item(k, 6 + j).setText(f"{found[c]:g}")

    def _import(self):
        if not self.results:
            self.reject()
            return
        per_file: dict[str, dict] = {}
        try:
            for k, (f, i, s) in enumerate(self.rows):
                d = per_file.setdefault(f, {"only": [], "conc": {}, "names": {}, "roles": {}})
                if self.table.item(k, 0).checkState() != Qt.Checked:
                    continue
                d["only"].append(i)
                d["names"][i] = self.table.item(k, 2).text().strip() or s.name
                d["roles"][i] = self.table.item(k, 5).text().strip()
                conc = {}
                for j, c in enumerate(self.compounds):
                    t = self.table.item(k, 6 + j).text().strip().replace(",", ".")
                    if t:
                        conc[c] = float(t)
                if conc:
                    d["conc"][i] = conc
        except ValueError as exc:
            error(self, f"Invalid concentration: {exc}")
            return
        total = 0
        for f, d in per_file.items():
            if not d["only"]:
                continue
            try:
                total += len(self.project.import_file(
                    f, self.trial_cb.currentData(), self._options(),
                    role=self.role_cb.currentText().strip(), only=d["only"],
                    concentrations=d["conc"], names=d["names"], roles=d["roles"]))
            except Exception as exc:
                error(self, f"{Path(f).name}: {exc}")
        self.win.statusBar().showMessage(f"Imported {total} spectra.")
        self.accept()


# --------------------------------------------------------------------------- #
# Trials / compounds / concentrations
# --------------------------------------------------------------------------- #
class TrialsDialog(Base):
    def __init__(self, win):
        super().__init__(win, "Trials")
        self.resize(760, 420)
        lay = QVBoxLayout(self)
        self.table = QTableWidget()
        self.table.setSelectionBehavior(QAbstractItemView.SelectRows)
        lay.addWidget(self.table)
        row = QHBoxLayout()
        for text, fn in (("Add…", self._add), ("Save edits…", self._save),
                         ("Archive…", self._archive)):
            b = QPushButton(text)
            b.clicked.connect(fn)
            row.addWidget(b)
        row.addStretch(1)
        close = QPushButton("Close")
        close.clicked.connect(self.accept)
        row.addWidget(close)
        lay.addLayout(row)
        self._load()

    def _load(self):
        self.trials = self.project.trials()
        fill_table(self.table, ["Name", "Description", "Status", "Created (UTC)"],
                   [[t["name"], t["description"], t["status"], t["created_utc"][:19]]
                    for t in self.trials], editable=True)

    def _add(self):
        name, ok = QInputDialog.getText(self, "New trial", "Name (e.g. 'Binary mixture PAR + CAF'):")
        if ok and name.strip():
            self.project.add_trial(name.strip())
            self._load()

    def _save(self):
        edits = []
        for i, t in enumerate(self.trials):
            vals = {k: self.table.item(i, j).text() for j, k in
                    enumerate(("name", "description", "status"))}
            if any(vals[k] != t[k] for k in vals):
                edits.append((t["id"], vals))
        if not edits:
            return
        reason = ask_reason(self, f"Save changes to {len(edits)} trials?")
        if reason is None:
            return
        for tid, vals in edits:
            self.project.update_trial(tid, reason, **vals)
        self._load()

    def _archive(self):
        r = self.table.currentRow()
        if r < 0:
            return
        reason = ask_reason(self, f"Archive trial '{self.trials[r]['name']}'?", required=True)
        if reason:
            self.project.archive_trial(self.trials[r]["id"], reason)
            self._load()


class CompoundsDialog(Base):
    def __init__(self, win):
        super().__init__(win, "Compounds")
        self.resize(640, 380)
        lay = QVBoxLayout(self)
        lay.addWidget(QLabel("Analytes used for concentrations and methods. Units are free "
                             "text (µg/mL, mg/L, M…)."))
        self.table = QTableWidget()
        lay.addWidget(self.table)
        row = QHBoxLayout()
        for text, fn in (("Add…", self._add), ("Save edits…", self._save)):
            b = QPushButton(text)
            b.clicked.connect(fn)
            row.addWidget(b)
        row.addStretch(1)
        close = QPushButton("Close")
        close.clicked.connect(self.accept)
        row.addWidget(close)
        lay.addLayout(row)
        self._load()

    def _load(self):
        self.items = self.project.compounds()
        fill_table(self.table, ["Name", "Unit", "MW (g/mol)", "Notes"],
                   [[c["name"], c["unit"] or "", "" if c["mw"] is None else f"{c['mw']:g}",
                     c["notes"] or ""] for c in self.items], editable=True)

    def _add(self):
        name, ok = QInputDialog.getText(self, "Add compound", "Name:")
        if not ok or not name.strip():
            return
        unit, ok = QInputDialog.getText(self, "Add compound", "Concentration unit:", text="µg/mL")
        try:
            self.project.add_compound(name.strip(), unit or "µg/mL")
        except Exception as exc:
            error(self, exc)
        self._load()

    def _save(self):
        edits = []
        for i, c in enumerate(self.items):
            mw = self.table.item(i, 2).text().strip()
            vals = {"name": self.table.item(i, 0).text().strip(),
                    "unit": self.table.item(i, 1).text().strip(),
                    "mw": float(mw) if mw else None,
                    "notes": self.table.item(i, 3).text()}
            if any(vals[k] != c[k] for k in vals):
                edits.append((c["id"], vals))
        if not edits:
            return
        reason = ask_reason(self, "Save compound changes? (Renaming does not rename "
                                  "concentrations already stored on spectra.)")
        if reason is None:
            return
        for cid, vals in edits:
            self.project.update_compound(cid, reason, **vals)
        self._load()


class ConcentrationsDialog(Base):
    """Spreadsheet-style editor of roles and concentrations (paste from Excel)."""

    def __init__(self, win, ids: list[int]):
        super().__init__(win, "Concentration table")
        self.resize(900, 520)
        self.ids = ids
        self.compounds = self.project.compound_names()
        lay = QVBoxLayout(self)
        lay.addWidget(QLabel("Edit or paste (Ctrl+V) from Excel. Leave blank for 'not set'; "
                             "use 0 when a compound is absent from a mixture."))
        if not self.compounds:
            lay.addWidget(QLabel("<b>No compounds defined — use Edit → Compounds… first.</b>"))
        self.table = PasteTable()
        lay.addWidget(self.table)
        recs = self.project.records(ids=ids, include_archived=True)
        self.recs = recs
        rows = []
        for r in recs:
            rows.append([f"#{r.id}", r.name, r.role] +
                        [("" if r.concentrations.get(c) is None else f"{r.concentrations[c]:g}")
                         for c in self.compounds])
        fill_table(self.table, ["ID", "Name", "Role"] + self.compounds, rows, editable=True)
        for i in range(len(rows)):
            it = self.table.item(i, 0)
            it.setFlags(it.flags() & ~Qt.ItemIsEditable)
        bb = QDialogButtonBox(QDialogButtonBox.Save | QDialogButtonBox.Cancel)
        bb.accepted.connect(self._save)
        bb.rejected.connect(self.reject)
        lay.addWidget(bb)

    def _save(self):
        conc_all, other = {}, []
        try:
            for i, r in enumerate(self.recs):
                conc = {k: v for k, v in r.concentrations.items() if k not in self.compounds}
                for j, c in enumerate(self.compounds):
                    t = self.table.item(i, 3 + j).text().strip().replace(",", ".")
                    if t:
                        conc[c] = float(t)
                conc_all[r.id] = conc
                name = self.table.item(i, 1).text().strip()
                role = self.table.item(i, 2).text().strip()
                if name != r.name or role != r.role:
                    other.append((r.id, name, role))
        except ValueError as exc:
            error(self, f"Invalid number: {exc}")
            return
        reason = ask_reason(self, "Save concentration table?")
        if reason is None:
            return
        self.project.set_concentrations_bulk(conc_all, reason)
        for sid, name, role in other:
            self.project.update_spectrum(sid, reason, name=name, role=role)
        self.accept()


# --------------------------------------------------------------------------- #
# Processing
# --------------------------------------------------------------------------- #
class ProcessDialog(Base):
    def __init__(self, win, ids: list[int]):
        super().__init__(win, f"Processing pipeline — {len(ids)} spectra")
        self.resize(1200, 720)
        self.ids = ids
        self.spectra = self.project.spectra(ids)
        lay = QHBoxLayout(self)
        left = QWidget()
        ll = QVBoxLayout(left)
        self.editor = PipelineEditor(spectrum_choices(self.project), name_of(self.project))
        self.editor.changed.connect(self._preview)
        ll.addWidget(QLabel("<b>Steps</b> (applied in order; raw data is never changed — "
                            "results are stored as new derived spectra)"))
        ll.addWidget(self.editor, 1)
        f = QFormLayout()
        self.suffix = QLineEdit()
        self.suffix.setPlaceholderText("automatic (e.g. 'Divide → Derivative')")
        f.addRow("Result label", self.suffix)
        ll.addLayout(f)
        row = QHBoxLayout()
        load = QPushButton("Load pipeline from spectrum…")
        load.setToolTip("Reuse the processing recorded on an existing derived spectrum")
        load.clicked.connect(self._load_from)
        row.addWidget(load)
        row.addStretch(1)
        ll.addLayout(row)
        bb = QDialogButtonBox(QDialogButtonBox.Apply | QDialogButtonBox.Close)
        bb.button(QDialogButtonBox.Apply).clicked.connect(self._apply)
        bb.rejected.connect(self.reject)
        ll.addWidget(bb)
        split = QSplitter()
        split.addWidget(left)
        right = QWidget()
        rl = QVBoxLayout(right)
        self.plot = SpectrumPlot(y_label="Signal")
        self.msg = QLabel()
        self.msg.setWordWrap(True)
        rl.addWidget(QLabel("<b>Live preview</b>"))
        rl.addWidget(self.plot, 1)
        rl.addWidget(self.msg)
        split.addWidget(right)
        split.setSizes([460, 740])
        lay.addWidget(split)
        self._preview()

    def _preview(self):
        steps = self.editor.get_steps()
        try:
            out = [apply_pipeline(s, steps, self.project.resolver()) for s in self.spectra[:30]]
            self.plot.plot_spectra(out, keep_range=False)
            self.msg.setText("" if steps else "Add a step to see its effect.")
            self.msg.setStyleSheet("")
        except Exception as exc:
            self.plot.plot_spectra(self.spectra[:30])
            self.msg.setText(f"⚠ {exc}")
            self.msg.setStyleSheet("color: #b3261e;")

    def _load_from(self):
        recs = [r for r in self.project.records() if r.kind == "derived"]
        if not recs:
            error(self, "No derived spectra yet.")
            return
        labels = [f"#{r.id} {r.name}" for r in recs]
        choice, ok = QInputDialog.getItem(self, "Load pipeline", "Spectrum:", labels, 0, False)
        if ok:
            chain = self.project.lineage(recs[labels.index(choice)].id)
            steps = [st for rec in chain for st in rec.pipeline]
            self.editor.set_steps(steps)

    def _apply(self):
        steps = self.editor.get_steps()
        if not steps:
            error(self, "Add at least one step.")
            return
        try:
            out = self.project.process(self.ids, steps, self.suffix.text().strip() or None)
        except Exception as exc:
            error(self, exc)
            return
        self.win.statusBar().showMessage(f"Created {len(out)} derived spectra.")
        self.accept()


# --------------------------------------------------------------------------- #
# Audit
# --------------------------------------------------------------------------- #
class AuditDialog(Base):
    def __init__(self, win):
        super().__init__(win, "Audit trail")
        self.resize(1250, 720)
        lay = QVBoxLayout(self)
        row = QHBoxLayout()
        self.action = QComboBox()
        self.action.addItems(["(all)", "CREATE", "IMPORT", "PROCESS", "EDIT", "ARCHIVE",
                              "RESTORE", "CALCULATE", "EXPORT", "VERIFY", "OPEN", "CLOSE"])
        self.text = QLineEdit()
        self.text.setPlaceholderText("search summary / details / reason")
        self.action.currentIndexChanged.connect(self._load)
        self.text.returnPressed.connect(self._load)
        row.addWidget(QLabel("Action"))
        row.addWidget(self.action)
        row.addWidget(self.text, 1)
        for t, fn in (("Search", self._load), ("Verify integrity", self._verify),
                      ("Export CSV…", self._export)):
            b = QPushButton(t)
            b.clicked.connect(fn)
            row.addWidget(b)
        lay.addLayout(row)
        split = QSplitter(Qt.Vertical)
        self.table = QTableWidget()
        self.table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.table.itemSelectionChanged.connect(self._details)
        split.addWidget(self.table)
        self.detail = QPlainTextEdit()
        self.detail.setReadOnly(True)
        split.addWidget(self.detail)
        split.setSizes([480, 200])
        lay.addWidget(split)
        self._load()

    def _load(self):
        act = self.action.currentText()
        self.entries = self.project.audit_entries(action=None if act == "(all)" else act,
                                                  text=self.text.text().strip() or None)
        fill_table(self.table, ["#", "Time (UTC)", "User", "Host", "Version", "Action", "Entity",
                                "Summary", "Reason", "Hash"],
                   [[str(e.seq), e.ts_utc.replace("T", " "), e.user, e.host, e.app_version,
                     e.action, f"{e.entity}{'' if e.entity_id is None else ' #' + str(e.entity_id)}",
                     e.summary, e.reason, e.hash[:12] + "…"] for e in self.entries])

    def _details(self):
        r = self.table.currentRow()
        if 0 <= r < len(self.entries):
            e = self.entries[r]
            self.detail.setPlainText(
                f"Entry #{e.seq}   {e.ts_utc}   {e.user}@{e.host}   v{e.app_version}\n"
                f"Previous hash: {e.prev_hash}\nThis hash:     {e.hash}\n\n"
                + json.dumps(e.details, indent=2, ensure_ascii=False))

    def _verify(self):
        self.win.verify()
        self._load()

    def _export(self):
        path, _ = QFileDialog.getSaveFileName(self, "Export audit trail", "audit_trail.csv",
                                              "CSV (*.csv)")
        if not path:
            return
        entries = self.project.audit_entries()
        with open(path, "w", newline="", encoding="utf-8-sig") as fh:
            w = csv.writer(fh)
            w.writerow(["seq", "ts_utc", "user", "host", "app_version", "action", "entity",
                        "entity_id", "summary", "reason", "details", "prev_hash", "hash"])
            for e in reversed(entries):
                w.writerow([e.seq, e.ts_utc, e.user, e.host, e.app_version, e.action, e.entity,
                            e.entity_id, e.summary, e.reason,
                            json.dumps(e.details, ensure_ascii=False), e.prev_hash, e.hash])
        self.project.log("EXPORT", f"Exported audit trail ({len(entries)} entries)",
                         {"file": path})
        self._load()


# --------------------------------------------------------------------------- #
# Export
# --------------------------------------------------------------------------- #
def export_spectra(win, ids: list[int]) -> None:
    project = win.project
    path, flt = QFileDialog.getSaveFileName(
        win, "Export spectra", "spectra.xlsx",
        "Excel (*.xlsx);;CSV (*.csv);;JCAMP-DX, one file per spectrum (*.jdx);;GRAMS SPC (*.spc)")
    if not path:
        return
    spectra = project.spectra(ids)
    try:
        if path.lower().endswith(".csv"):
            export_csv(spectra, path)
            files = [path]
        elif path.lower().endswith(".jdx"):
            base = Path(path)
            files = []
            for s in spectra:
                safe = re.sub(r"[^\w\-. ]", "_", s.name)[:80]
                f = base.with_name(f"{base.stem}_{safe}.jdx") if len(spectra) > 1 else base
                write_jcamp(s, f)
                files.append(str(f))
        elif path.lower().endswith(".spc"):
            write_spc(spectra, path)
            files = [path]
        else:
            if not path.lower().endswith(".xlsx"):
                path += ".xlsx"
            export_excel(spectra, path)
            files = [path]
    except Exception as exc:
        error(win, exc)
        return
    project.log("EXPORT", f"Exported {len(ids)} spectra", {"inputs": ids, "files": files})
    win.statusBar().showMessage(f"Exported {len(ids)} spectra.")
