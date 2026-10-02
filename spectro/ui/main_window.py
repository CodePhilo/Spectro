"""Main application window."""

from __future__ import annotations

import sys
from pathlib import Path

from PySide6.QtCore import QSettings, Qt
from PySide6.QtGui import QAction, QKeySequence
from PySide6.QtWidgets import (QAbstractItemView, QApplication, QComboBox, QDockWidget,
                               QDoubleSpinBox,
                               QFileDialog, QFormLayout, QHBoxLayout, QInputDialog, QLabel,
                               QLineEdit, QListWidget, QMainWindow, QMessageBox, QPushButton,
                               QTabWidget, QTableWidget, QTreeWidget, QTreeWidgetItem,
                               QVBoxLayout, QWidget, QPlainTextEdit)

from spectro import __version__
from spectro.core.operations import describe_step
from spectro.storage.project import ROLES, Project
from spectro.ui import widgets
from spectro.ui.widgets import (VIEW_MODES, PasteTable, SpectrumPlot, ask_reason, error,
                                fill_table)

APP = "Spectro"


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.project: Project | None = None
        self.settings = QSettings("Spectro", "Spectro")
        self.setWindowTitle(APP)
        self.resize(1400, 850)
        self._build_ui()
        self._build_menus()
        self._update_enabled()

    # ------------------------------------------------------------ layout
    def _build_ui(self) -> None:
        self.plot = SpectrumPlot()
        self.plot.hovered.connect(lambda t: self.statusBar().showMessage(t))
        central = QWidget()
        cv = QVBoxLayout(central)
        cv.setContentsMargins(0, 0, 0, 0)
        bar = QHBoxLayout()
        bar.setContentsMargins(6, 4, 6, 0)
        self.view_cb = QComboBox()
        self.view_cb.addItems(list(VIEW_MODES))
        self.view_cb.setToolTip("Overlay · Stacked (offset) · Difference (each minus the first "
                                "selected) · Normalized (max = 1)")
        self.offset_sb = QDoubleSpinBox()
        self.offset_sb.setDecimals(4)
        self.offset_sb.setRange(0, 1e6)
        self.offset_sb.setSingleStep(0.05)
        self.offset_sb.setSpecialValueText("auto")
        self.offset_sb.setToolTip("Offset between stacked spectra (0 = automatic)")
        self.offset_sb.setEnabled(False)
        self.view_cb.currentTextChanged.connect(self._view_changed)
        self.offset_sb.valueChanged.connect(self._view_changed)
        fig_btn = QPushButton("Export figure…")
        fig_btn.clicked.connect(self.export_figure)
        bar.addWidget(QLabel("View"))
        bar.addWidget(self.view_cb)
        bar.addWidget(QLabel("Offset"))
        bar.addWidget(self.offset_sb)
        bar.addStretch(1)
        bar.addWidget(fig_btn)
        cv.addLayout(bar)
        cv.addWidget(self.plot, 1)
        self.setCentralWidget(central)

        self.tree = QTreeWidget()
        self.tree.setHeaderLabels(["Name", "Role", "Pts", "Concentrations"])
        self.tree.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.tree.itemSelectionChanged.connect(self._selection_changed)
        self.tree.setColumnWidth(0, 260)
        self.tree.setContextMenuPolicy(Qt.ActionsContextMenu)
        self.filter_edit = QLineEdit()
        self.filter_edit.setPlaceholderText("Filter spectra…")
        self.filter_edit.textChanged.connect(self._apply_filter)
        self.show_archived = QPushButton("Show archived")
        self.show_archived.setCheckable(True)
        self.show_archived.toggled.connect(self.refresh_tree)
        left = QWidget()
        ll = QVBoxLayout(left)
        ll.setContentsMargins(4, 4, 4, 4)
        top = QHBoxLayout()
        top.addWidget(self.filter_edit, 1)
        top.addWidget(self.show_archived)
        ll.addLayout(top)
        ll.addWidget(self.tree)
        dock = QDockWidget("Project", self)
        dock.setObjectName("project_dock")
        dock.setWidget(left)
        self.addDockWidget(Qt.LeftDockWidgetArea, dock)

        # properties
        self.props = QTabWidget()
        w = QWidget()
        f = QFormLayout(w)
        self.p_name = QLineEdit()
        self.p_role = QComboBox()
        self.p_role.addItems(ROLES)
        self.p_role.setEditable(True)
        self.p_info = QLabel()
        self.p_info.setWordWrap(True)
        self.p_info.setTextInteractionFlags(Qt.TextSelectableByMouse)
        self.p_notes = QPlainTextEdit()
        self.p_notes.setMaximumHeight(70)
        self.p_conc = PasteTable(0, 2)
        self.p_conc.setHorizontalHeaderLabels(["Compound", "Concentration"])
        self.p_conc.setMaximumHeight(150)
        save = QPushButton("Save changes…")
        save.clicked.connect(self._save_properties)
        f.addRow("Name", self.p_name)
        f.addRow("Role", self.p_role)
        f.addRow("Concentrations", self.p_conc)
        f.addRow("Notes", self.p_notes)
        f.addRow(save)
        f.addRow(self.p_info)
        self.props.addTab(w, "Properties")
        self.p_meta = QTableWidget()
        self.props.addTab(self.p_meta, "Metadata")
        lw = QWidget()
        lv = QVBoxLayout(lw)
        self.p_lineage = QListWidget()
        lv.addWidget(QLabel("Lineage (raw → this spectrum):"))
        lv.addWidget(self.p_lineage)
        rb = QPushButton("Verify: recompute from raw data")
        rb.clicked.connect(self._replay)
        lv.addWidget(rb)
        eb = QPushButton("Edit processing…")
        eb.setToolTip("Change the steps that made this spectrum (a new version is stored)")
        eb.clicked.connect(self._edit_prop_processing)
        lv.addWidget(eb)
        self.props.addTab(lw, "Lineage")
        self.p_history = QTableWidget()
        self.props.addTab(self.p_history, "History")
        pdock = QDockWidget("Selected spectrum", self)
        pdock.setObjectName("props_dock")
        pdock.setWidget(self.props)
        self.addDockWidget(Qt.RightDockWidgetArea, pdock)

        # audit
        self.audit_table = QTableWidget()
        self.audit_table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.audit_table.setSelectionBehavior(QAbstractItemView.SelectRows)
        adock = QDockWidget("Audit trail (latest)", self)
        adock.setObjectName("audit_dock")
        adock.setWidget(self.audit_table)
        self.addDockWidget(Qt.BottomDockWidgetArea, adock)
        self.docks = [dock, pdock, adock]
        self.statusBar().showMessage(
            "Create or open a project to start — or try Help → Open demo project.")

    def _build_menus(self) -> None:
        mb = self.menuBar()
        self.project_actions: list[QAction] = []

        def act(menu, text, slot, shortcut=None, needs_project=True, tip=""):
            a = QAction(text, self)
            if shortcut:
                a.setShortcut(QKeySequence(shortcut))
            a.setStatusTip(tip)
            a.triggered.connect(slot)
            menu.addAction(a)
            if needs_project:
                self.project_actions.append(a)
            return a

        m = mb.addMenu("&File")
        act(m, "&New project…", self.new_project, "Ctrl+N", False)
        act(m, "&Open project…", self.open_project, "Ctrl+O", False)
        self.recent_menu = m.addMenu("Open &recent")
        self._fill_recent()
        act(m, "&Close project", self.close_project)
        m.addSeparator()
        act(m, "&Import spectra (any instrument / Excel / CSV)…", self.import_spectra, "Ctrl+I")
        act(m, "&Export selected spectra…", self.export_spectra, "Ctrl+E")
        act(m, "Export &figure (publication quality)…", self.export_figure, "Ctrl+Shift+E",
            needs_project=False)
        act(m, "Export results to E&xcel…", self.export_results)
        act(m, "Trial &report (PDF/HTML)…", self.report)
        m.addSeparator()
        act(m, "E&xit", self.close, "Ctrl+Q", False)

        m = mb.addMenu("&Edit")
        act(m, "&Undo a change…", self.undo_dialog, "Ctrl+Z",
            tip="Restore the previous version of a recent edit (kept in the audit trail)")
        m.addSeparator()
        act(m, "&Trials…", self.edit_trials)
        act(m, "&Compounds…", self.edit_compounds)
        act(m, "Concentration &table for selection…", self.edit_concentrations, "Ctrl+T")
        act(m, "Set &role of selection…", self.set_role)
        act(m, "&Move selection to trial…", self.move_to_trial)
        m.addSeparator()
        act(m, "&Archive selection…", self.archive, "Del")
        act(m, "&Restore selection", self.restore)
        for a in self.project_actions[-6:]:
            self.tree.addAction(a)

        m = mb.addMenu("&Process")
        act(m, "&Processing pipeline…", self.process, "Ctrl+P",
            tip="Smoothing, baseline, derivatives, ratio spectra, resolution…")
        act(m, "&Edit processing of selected spectra…", self.edit_processing, "Ctrl+Shift+P",
            tip="Change the steps of processed spectra; everything built on them is rebuilt")

        m = mb.addMenu("&Methods")
        m.addSection("One compound per method")
        act(m, "&Univariate calibration — zero order, derivative, ratio, spectrum "
               "resolution (29 templates)…", self.univariate, "Ctrl+U")
        m.addSection("All compounds at once")
        act(m, "&Equation methods (Vierordt, bivariate, AUC, multi-λ)…", self.equations)
        act(m, "&Binary two-signal methods (Q, AS, AAS, AM, IAM, HPSAM)…", self.special)
        act(m, "&Progressive resolution (AAC, MACM, RIDSS, CV-AD, MAFM)…", self.progressive,
            tip="Ternary methods that resolve every compound from one ratio spectrum or a "
                "chain of absorption factors")
        act(m, "&Chemometrics (CLS, ILS, PCR, PLS, MCR-ALS, ANN, SVR)…", self.chemometrics,
            "Ctrl+M")
        m.addSeparator()
        act(m, "Saved methods && results…", self.saved)
        act(m, "Method &guide (how to apply each method)…", self.method_guide, needs_project=False,
            tip="Open the illustrated guide to every method")

        m = mb.addMenu("&Tools")
        act(m, "&Method optimizer (find the best processing)…", self.optimizer, "Ctrl+Shift+O",
            tip="Rank processing strategies for simultaneous determination")
        act(m, "&Spectral finder (zero-crossing, isoabsorptive, extrema, plateau)…",
            self.finder, "Ctrl+F")
        act(m, "&Validation statistics (ICH Q2)…", self.validation, needs_project=False)
        act(m, "Calibration &design (multilevel multifactor)…", self.design, needs_project=False)
        act(m, "&Greenness assessment (AGREE, Eco-Scale, GAPI)…", self.greenness,
            needs_project=False)

        m = mb.addMenu("&Audit")
        act(m, "&Audit trail…", self.audit_viewer)
        act(m, "&Verify data integrity", self.verify)

        m = mb.addMenu("&View")
        for d in self.docks:
            m.addAction(d.toggleViewAction())

        m = mb.addMenu("&Help")
        act(m, "Open &demo project…", self.open_demo, needs_project=False,
            tip="Paracetamol / caffeine / aspirin demo with methods and results")
        act(m, "&Copy demo data files…", self.copy_demo, needs_project=False,
            tip="Instrument-style files for trying File → Import spectra")
        m.addSeparator()
        act(m, "&About", self.about, needs_project=False)

    def _update_enabled(self) -> None:
        for a in self.project_actions:
            a.setEnabled(self.project is not None)
        self.setWindowTitle(f"{APP} — {self.project.path.name}" if self.project else APP)

    # ------------------------------------------------------------ project
    def _fill_recent(self) -> None:
        self.recent_menu.clear()
        for p in self.settings.value("recent", [], list) or []:
            a = self.recent_menu.addAction(p)
            a.triggered.connect(lambda _=False, path=p: self._open(path))

    def _remember(self, path: Path) -> None:
        rec = [str(path)] + [p for p in (self.settings.value("recent", [], list) or [])
                             if p != str(path)]
        self.settings.setValue("recent", rec[:8])
        self._fill_recent()

    def new_project(self) -> None:
        path, _ = QFileDialog.getSaveFileName(self, "New project", "", "Spectro project (*.spectro)")
        if not path:
            return
        if not path.endswith(".spectro"):
            path += ".spectro"
        if Path(path).exists():
            error(self, "That file already exists; choose Open project instead.")
            return
        self.close_project()
        try:
            self._attach(Project.create(path))
        except Exception as exc:
            error(self, exc)
            return
        self.project.add_trial("Trial 1")

    def open_project(self) -> None:
        path, _ = QFileDialog.getOpenFileName(self, "Open project", "", "Spectro project (*.spectro)")
        if path:
            self._open(path)

    def _open(self, path: str) -> None:
        self.close_project()
        try:
            proj = Project.open(path)
        except Exception as exc:
            error(self, exc)
            return
        self._attach(proj)
        report = proj.verify()
        if not report["ok"]:
            QMessageBox.critical(self, "Integrity check failed",
                                 "This project failed its integrity check:\n\n"
                                 + "\n".join(report["problems"][:20]))

    def _attach(self, proj: Project) -> None:
        self.project = proj
        widgets.EXPORT_HOOK = lambda summary, details: self.project.log("EXPORT", summary,
                                                                       details)
        proj.subscribe(self._changed)
        self._remember(proj.path)
        self._update_enabled()
        self.refresh_tree()
        self.refresh_audit()

    def close_project(self) -> None:
        widgets.EXPORT_HOOK = None
        if self.project:
            self.project.close()
            self.project = None
        self.tree.clear()
        self.plot.plot_spectra([])
        self.audit_table.setRowCount(0)
        self._update_enabled()

    def closeEvent(self, e) -> None:
        self.close_project()
        super().closeEvent(e)

    def _changed(self, entity: str, entity_id) -> None:
        if entity in ("spectrum", "trial", "compound"):
            self.refresh_tree()
        self.refresh_audit()

    # ------------------------------------------------------------ tree
    def refresh_tree(self) -> None:
        if not self.project:
            return
        sel = set(self.selected_ids())
        expanded = {self.tree.topLevelItem(i).data(0, Qt.UserRole + 1)
                    for i in range(self.tree.topLevelItemCount())
                    if self.tree.topLevelItem(i).isExpanded()}
        self.tree.blockSignals(True)
        self.tree.clear()
        arch = self.show_archived.isChecked()
        trials = self.project.trials(include_archived=arch)
        recs = self.project.records(include_archived=arch)
        by_id = {r.id: r for r in recs}
        groups = [(t["id"], t["name"] + ("  [archived]" if t["archived"] else "")) for t in trials]
        if any(r.trial_id is None for r in recs):
            groups.append((None, "Unassigned"))
        first = not expanded
        for tid, name in groups:
            top = QTreeWidgetItem([name])
            top.setData(0, Qt.UserRole + 1, tid)
            f = top.font(0)
            f.setBold(True)
            top.setFont(0, f)
            top.setFlags(top.flags() & ~Qt.ItemIsSelectable)
            self.tree.addTopLevelItem(top)
            items: dict[int, QTreeWidgetItem] = {}
            for r in [r for r in recs if r.trial_id == tid]:
                conc = ", ".join(f"{k}={v:g}" for k, v in r.concentrations.items())
                it = QTreeWidgetItem([r.name, r.role, str(r.npoints), conc])
                it.setData(0, Qt.UserRole, r.id)
                it.setToolTip(0, f"#{r.id} {r.kind}  {r.name}")
                if r.kind == "derived":
                    it.setForeground(0, Qt.darkBlue)
                if r.archived:
                    it.setForeground(0, Qt.gray)
                parent = items.get(r.parent_ids[0]) if r.parent_ids and r.parent_ids[0] in by_id else None
                (parent or top).addChild(it)
                items[r.id] = it
                if r.id in sel:
                    it.setSelected(True)
            top.setExpanded(first or tid in expanded)
        self.tree.blockSignals(False)
        self._apply_filter(self.filter_edit.text())
        self._selection_changed()

    def _apply_filter(self, text: str) -> None:
        text = text.lower().strip()

        def visit(it: QTreeWidgetItem) -> bool:
            child_vis = any([visit(it.child(i)) for i in range(it.childCount())])
            own = not text or text in " ".join(it.text(c) for c in range(4)).lower()
            vis = own or child_vis
            it.setHidden(not vis)
            return vis
        for i in range(self.tree.topLevelItemCount()):
            top = self.tree.topLevelItem(i)
            for j in range(top.childCount()):
                visit(top.child(j))

    def selected_ids(self) -> list[int]:
        return [it.data(0, Qt.UserRole) for it in self.tree.selectedItems()
                if it.data(0, Qt.UserRole) is not None]

    def current_trial(self) -> int | None:
        it = self.tree.currentItem()
        while it is not None and it.parent() is not None:
            it = it.parent()
        if it is not None:
            return it.data(0, Qt.UserRole + 1)
        trials = self.project.trials() if self.project else []
        return trials[0]["id"] if trials else None

    def _selection_changed(self) -> None:
        if not self.project:
            return
        ids = self.selected_ids()
        try:
            spectra = self.project.spectra(ids[:200])
        except KeyError:
            spectra = []
        self.plot.plot_spectra(spectra)
        self._show_properties(ids[0] if len(ids) == 1 else None)
        if len(ids) > 1:
            self.statusBar().showMessage(f"{len(ids)} spectra selected")

    # ------------------------------------------------------------ properties
    def _show_properties(self, sid: int | None) -> None:
        self._prop_id = sid
        enabled = sid is not None
        self.props.setEnabled(enabled)
        if not enabled:
            self.p_name.clear()
            self.p_info.clear()
            self.p_conc.setRowCount(0)
            self.p_lineage.clear()
            self.p_history.setRowCount(0)
            return
        r = self.project.record(sid)
        self.p_name.setText(r.name)
        self.p_role.setCurrentText(r.role)
        self.p_notes.setPlainText(r.notes)
        comps = self.project.compound_names()
        names = comps + [k for k in r.concentrations if k not in comps]
        self.p_conc.setRowCount(len(names))
        from PySide6.QtWidgets import QTableWidgetItem
        for i, n in enumerate(names):
            a = QTableWidgetItem(n)
            a.setFlags(a.flags() & ~Qt.ItemIsEditable)
            self.p_conc.setItem(i, 0, a)
            v = r.concentrations.get(n)
            self.p_conc.setItem(i, 1, QTableWidgetItem("" if v is None else f"{v:g}"))
        s = self.project.spectrum(sid)
        self.p_info.setText(
            f"#{r.id} · {r.kind} · {r.npoints} points · {s.wavelengths[0]:g}–{s.wavelengths[-1]:g} nm"
            f" (step {s.step:g})\nCreated {r.created_utc}\n"
            + (f"Source: {r.source_file}\nFile SHA-256: {r.source_sha256[:16]}…\n" if r.source_file else "")
            + f"Data SHA-256: {r.data_sha256[:16]}…")
        fill_table(self.p_meta, ["Key", "Value"], [[k, str(v)] for k, v in r.metadata.items()])
        self.p_lineage.clear()
        name_of = lambda ref: self.project.record(int(ref)).name  # noqa: E731
        for rec in self.project.lineage(sid):
            self.p_lineage.addItem(f"#{rec.id}  {rec.name}  ({rec.kind})")
            for st in rec.pipeline:
                self.p_lineage.addItem(f"      ↳ {describe_step(st, name_of)}")
        hist = self.project.spectrum_history(sid)
        fill_table(self.p_history, ["#", "Time (UTC)", "User", "Action", "Summary", "Reason"],
                   [[str(e.seq), e.ts_utc[:19].replace("T", " "), e.user, e.action, e.summary,
                     e.reason] for e in hist])

    def _save_properties(self) -> None:
        sid = getattr(self, "_prop_id", None)
        if sid is None:
            return
        conc = {}
        try:
            for i in range(self.p_conc.rowCount()):
                n = self.p_conc.item(i, 0).text()
                it = self.p_conc.item(i, 1)
                t = it.text().strip().replace(",", ".") if it else ""
                if t:
                    conc[n] = float(t)
        except ValueError as exc:
            error(self, f"Invalid concentration: {exc}")
            return
        reason = ask_reason(self, "Save changes to this spectrum's description?")
        if reason is None:
            return
        before = self.project.record(sid).concentrations
        try:
            self.project.update_spectrum(sid, reason, name=self.p_name.text().strip(),
                                         role=self.p_role.currentText().strip(),
                                         notes=self.p_notes.toPlainText(), concentrations=conc)
        except Exception as exc:
            error(self, exc)
            return
        if conc != before:
            from spectro.ui.dialogs_data import notify_impact
            notify_impact(self, self.project, [sid] + self.project.descendants(sid),
                          "Concentrations saved (processed versions of this sample were "
                          "updated too).")

    def _edit_prop_processing(self) -> None:
        sid = getattr(self, "_prop_id", None)
        if sid is None:
            return
        from spectro.ui.dialogs_data import ProcessDialog
        try:
            dlg = ProcessDialog(self, [sid], edit=True)
        except Exception as exc:
            error(self, exc)
            return
        dlg.exec()

    def _replay(self) -> None:
        sid = getattr(self, "_prop_id", None)
        if sid is None:
            return
        try:
            res = self.project.replay(sid)
        except Exception as exc:
            error(self, exc)
            return
        (QMessageBox.information if res["ok"] else QMessageBox.critical)(
            self, "Verification", res["message"])

    def refresh_audit(self) -> None:
        if not self.project:
            return
        entries = self.project.audit_entries(limit=200)
        fill_table(self.audit_table, ["#", "Time (UTC)", "User", "Host", "Action", "Reason", "Summary"],
                   [[str(e.seq), e.ts_utc[:19].replace("T", " "), e.user, e.host, e.action,
                     e.reason, e.summary] for e in entries])

    # ------------------------------------------------------------ edit actions
    def _need_selection(self) -> list[int] | None:
        ids = self.selected_ids()
        if not ids:
            error(self, "Select one or more spectra in the project tree first.")
            return None
        return ids

    def set_role(self) -> None:
        ids = self._need_selection()
        if not ids:
            return
        role, ok = QInputDialog.getItem(self, "Set role", "Role:", list(ROLES), 0, True)
        if not ok:
            return
        reason = ask_reason(self, f"Set role '{role}' on {len(ids)} spectra?")
        if reason is None:
            return
        for sid in ids:
            self.project.update_spectrum(sid, reason, role=role)

    def move_to_trial(self) -> None:
        ids = self._need_selection()
        if not ids:
            return
        trials = self.project.trials()
        names = [t["name"] for t in trials]
        name, ok = QInputDialog.getItem(self, "Move to trial", "Trial:", names, 0, False)
        if not ok:
            return
        tid = trials[names.index(name)]["id"]
        reason = ask_reason(self, f"Move {len(ids)} spectra to '{name}'?")
        if reason is None:
            return
        for sid in ids:
            self.project.update_spectrum(sid, reason, trial_id=tid)

    def archive(self) -> None:
        ids = self._need_selection()
        if not ids:
            return
        reason = ask_reason(self, f"Archive {len(ids)} spectra? (Nothing is deleted; "
                                  "archived spectra can be restored.)", required=True)
        if reason:
            self.project.archive_spectra(ids, reason)

    def restore(self) -> None:
        ids = self._need_selection()
        if not ids:
            return
        reason = ask_reason(self, f"Restore {len(ids)} spectra?")
        if reason is not None:
            self.project.archive_spectra(ids, reason, archived=False)

    # ------------------------------------------------------------ dialogs
    def _dialog(self, module: str, cls: str, *args):
        from spectro.ui import dialog_optimizer, dialogs_data, dialogs_methods, dialogs_tools

        mod = {"dialogs_data": dialogs_data, "dialogs_methods": dialogs_methods,
               "dialogs_tools": dialogs_tools, "dialog_optimizer": dialog_optimizer}[module]
        dlg = getattr(mod, cls)(self, *args)
        dlg.exec()
        return dlg

    def import_spectra(self):
        self._dialog("dialogs_data", "ImportDialog")

    def export_spectra(self):
        from spectro.ui.dialogs_data import export_spectra
        ids = self._need_selection()
        if ids:
            export_spectra(self, ids)

    def edit_trials(self):
        self._dialog("dialogs_data", "TrialsDialog")

    def edit_compounds(self):
        self._dialog("dialogs_data", "CompoundsDialog")

    def edit_concentrations(self):
        ids = self._need_selection()
        if ids:
            self._dialog("dialogs_data", "ConcentrationsDialog", ids)

    def process(self):
        ids = self._need_selection()
        if ids:
            self._dialog("dialogs_data", "ProcessDialog", ids)

    def undo_dialog(self):
        self._dialog("dialogs_data", "UndoDialog")

    def edit_processing(self):
        ids = self._need_selection()
        if not ids:
            return
        from spectro.ui.dialogs_data import ProcessDialog
        try:
            dlg = ProcessDialog(self, ids, edit=True)
        except Exception as exc:
            error(self, exc)
            return
        dlg.exec()

    def audit_viewer(self):
        self._dialog("dialogs_data", "AuditDialog")

    def univariate(self):
        self._dialog("dialogs_methods", "UnivariateDialog")

    def equations(self):
        self._dialog("dialogs_methods", "EquationsDialog")

    def special(self):
        self._dialog("dialogs_methods", "SpecialDialog")

    def method_guide(self):
        """Open docs/methods/index.html (source checkout or bundled data)."""
        from PySide6.QtCore import QUrl
        from PySide6.QtGui import QDesktopServices
        base = Path(getattr(sys, "_MEIPASS", Path(__file__).resolve().parents[2]))
        here = base / "docs" / "methods" / "index.html"
        if not here.exists():
            QMessageBox.information(self, "Method guide", "The guide is in docs/methods of the "
                                    "Spectro source (README.md / index.html).")
            return
        QDesktopServices.openUrl(QUrl.fromLocalFile(str(here)))

    def progressive(self):
        self._dialog("dialogs_methods", "ProgressiveDialog")

    def chemometrics(self):
        self._dialog("dialogs_methods", "ChemometricsDialog")

    def saved(self):
        self._dialog("dialogs_methods", "SavedDialog")

    def optimizer(self):
        self._dialog("dialog_optimizer", "OptimizerDialog")

    def finder(self):
        self._dialog("dialogs_tools", "FinderDialog")

    def validation(self):
        self._dialog("dialogs_tools", "ValidationDialog")

    def design(self):
        self._dialog("dialogs_tools", "DesignDialog")

    def greenness(self):
        self._dialog("dialogs_tools", "GreennessDialog")

    def report(self):
        from spectro.ui.report import export_report
        export_report(self)

    def verify(self):
        rep = self.project.verify()
        self.project.log("VERIFY", "Integrity check: " + ("passed" if rep["ok"] else "FAILED"), rep)
        if rep["ok"]:
            QMessageBox.information(self, "Integrity check",
                                    f"Passed.\n\n{rep['audit_entries']} audit entries in an "
                                    f"unbroken hash chain; {rep['spectra']} spectra match their "
                                    "data hashes.")
        else:
            QMessageBox.critical(self, "Integrity check FAILED", "\n".join(rep["problems"][:30]))

    def _view_changed(self, *_):
        mode = self.view_cb.currentText()
        self.offset_sb.setEnabled(mode == "Stacked")
        self.plot.set_view(mode, self.offset_sb.value())

    def export_figure(self):
        from spectro.ui.figure import FigureDialog

        if not self.plot.curves:
            error(self, "Select spectra to plot first.")
            return
        FigureDialog(self.plot, self).exec()

    def export_results(self):
        from spectro.storage.tables import export_results_excel

        trials = self.project.trials()
        labels = ["(all trials)"] + [t["name"] for t in trials]
        choice, ok = QInputDialog.getItem(self, "Export results", "Trial:", labels, 0, False)
        if not ok:
            return
        tid = None if choice == labels[0] else trials[labels.index(choice) - 1]["id"]
        path, _ = QFileDialog.getSaveFileName(self, "Export results", "results.xlsx",
                                              "Excel (*.xlsx)")
        if not path:
            return
        if not path.lower().endswith(".xlsx"):
            path += ".xlsx"
        try:
            n = export_results_excel(self.project, path, tid)
        except Exception as exc:
            error(self, exc)
            return
        self.project.log("EXPORT", f"Exported {n} results to Excel", {"file": path,
                                                                      "trial_id": tid})
        self.statusBar().showMessage(f"Exported {n} results to {path}")

    def open_demo(self):
        from spectro.demo import build_demo_project

        default = Path.home() / "Documents" / "Spectro demo.spectro"
        if not default.parent.exists():
            default = Path.home() / "Spectro demo.spectro"
        path, _ = QFileDialog.getSaveFileName(self, "Create demo project", str(default),
                                              "Spectro project (*.spectro)")
        if not path:
            return
        if not path.endswith(".spectro"):
            path += ".spectro"
        if Path(path).exists():
            if QMessageBox.question(self, "Demo project",
                                    f"{Path(path).name} already exists. Open it?") \
                    == QMessageBox.Yes:
                self._open(path)
            return
        self.close_project()
        QApplication.setOverrideCursor(Qt.WaitCursor)
        try:
            proj = build_demo_project(path, lambda msg: (self.statusBar().showMessage(msg),
                                                         QApplication.processEvents()))
        except Exception as exc:
            QApplication.restoreOverrideCursor()
            error(self, exc)
            return
        QApplication.restoreOverrideCursor()
        self._attach(proj)
        QMessageBox.information(
            self, "Demo project",
            "Demo project created.\n\n"
            "• Two trials: binary paracetamol + caffeine (Panadol Extra) and ternary "
            "aspirin + paracetamol + caffeine (Excedrin).\n"
            "• Select spectra in the tree to plot them; derived spectra (ratio, D1) are "
            "nested under their parents.\n"
            "• Methods → Saved methods & results shows the ratio difference, Vierordt and "
            "PLS2 methods and their recoveries.\n"
            "• Audit → Audit trail shows every step used to build this project.\n\n"
            "The spectra are simulated from published absorptivities with realistic "
            "noise and errors — not measured data.")

    def copy_demo(self):
        from spectro.demo import copy_demo_files

        folder = QFileDialog.getExistingDirectory(self, "Copy demo data files to…",
                                                  str(Path.home()))
        if not folder:
            return
        try:
            target = copy_demo_files(folder)
        except Exception as exc:
            error(self, exc)
            return
        QMessageBox.information(self, "Demo data",
                                f"Copied to:\n{target}\n\nUse File → Import spectra… to try "
                                "each instrument format (Shimadzu-style TXT, Cary-style CSV, "
                                "Excel, European CSV, JCAMP-DX, row-wise TSV, GRAMS SPC).")

    def about(self):
        QMessageBox.about(self, "About Spectro",
                          f"<b>Spectro {__version__}</b><br>UV/Vis spectral trial organizer for "
                          "pharmaceutical analysis.<br><br>Raw data is immutable; every "
                          "operation is recorded in a hash-chained audit trail.")


def main(argv: list[str] | None = None) -> int:
    argv = argv if argv is not None else sys.argv
    if "--selfcheck" in argv:
        from spectro.ui._selfcheck import run
        print(run())
        return 0
    app = QApplication(argv)
    app.setApplicationName(APP)
    app.setOrganizationName("Spectro")
    app.setStyle("Fusion")
    win = MainWindow()
    win.show()
    args = app.arguments()[1:]
    if args and args[0].endswith(".spectro") and Path(args[0]).exists():
        win._open(args[0])
    return app.exec()
