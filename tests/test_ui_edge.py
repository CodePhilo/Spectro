"""Every dialog and action on empty / incomplete projects must show a message,
never crash."""

import pytest

pytest.importorskip("PySide6")

from PySide6.QtCore import Qt  # noqa: E402
from PySide6.QtWidgets import QApplication  # noqa: E402

from spectro.storage.project import Project  # noqa: E402
from tests.conftest import mixture  # noqa: E402


@pytest.fixture(scope="module")
def app():
    return QApplication.instance() or QApplication([])


@pytest.fixture
def messages(monkeypatch):
    """Collect user-facing error messages instead of showing message boxes."""
    import spectro.ui.dialog_optimizer as do
    import spectro.ui.dialogs_data as dd
    import spectro.ui.dialogs_methods as dm
    import spectro.ui.dialogs_tools as dt
    import spectro.ui.main_window as mw
    import spectro.ui.widgets as w
    got = []

    def err(parent, exc, title="Spectro"):
        got.append(str(exc))
    for mod in (w, dd, dm, dt, mw, do):
        monkeypatch.setattr(mod, "error", err, raising=False)
        monkeypatch.setattr(mod, "ask_reason", lambda *a, **k: "t", raising=False)
    from PySide6.QtWidgets import QInputDialog, QMessageBox
    for name in ("information", "critical", "warning"):
        monkeypatch.setattr(QMessageBox, name, lambda *a, **k: None)
    monkeypatch.setattr(QMessageBox, "question", lambda *a, **k: QMessageBox.Yes)
    monkeypatch.setattr(QInputDialog, "getText", lambda *a, **k: ("x", False))
    monkeypatch.setattr(QInputDialog, "getItem", lambda *a, **k: ("", False))
    monkeypatch.setattr(QInputDialog, "getInt", lambda *a, **k: (4, True))
    return got


def _window(app, tmp_path, with_data=False, compounds=True):
    import spectro.ui.main_window as mw
    win = mw.MainWindow()
    p = Project.create(tmp_path / "e.spectro")
    if compounds:
        for c in ("X", "Y"):
            p.add_compound(c)
    if with_data:
        tid = p.add_trial("T")
        for i, conc in enumerate([{"X": 5}, {"Y": 5}, {"X": 3, "Y": 4}]):
            p.add_spectrum(mixture(conc, f"s{i}"), tid, role="mixture")
    win._attach(p)
    return win


ACTIONS = {
    "dialogs_methods.UnivariateDialog": ["_calibrate", "_predict", "_save_method",
                                         "_save_results", "_std_addition", "_robustness"],
    "dialogs_methods.EquationsDialog": ["_fit", "_predict", "_save_method", "_save_results",
                                        "_kaiser"],
    "dialogs_methods.SpecialDialog": ["_run_q", "_run_as", "_run_am", "_run_iam", "_run_aas",
                                      "_aas_find", "_run_h", "_save"],
    "dialogs_methods.ProgressiveDialog": ["_run_ac", "_run_af", "_save", "_save_method"],
    "dialogs_methods.ChemometricsDialog": ["_cv", "_fit", "_ipls", "_ga", "_save_method",
                                           "_save_results", "_exclude_outliers"],
    "dialogs_methods.SavedDialog": ["_apply", "_export_model", "_archive"],
    "dialogs_tools.FinderDialog": ["_zero", "_iso", "_extrema", "_plateau", "_equal"],
    "dialogs_tools.DesignDialog": ["_gen"],
    "dialogs_tools.GreennessDialog": ["_eco", "_gapi", "_save"],
    "dialog_optimizer.OptimizerDialog": ["_run", "_save"],
    "dialogs_data.TrialsDialog": ["_save", "_archive"],
    "dialogs_data.CompoundsDialog": ["_save"],
    "dialogs_data.AuditDialog": ["_load"],
}


@pytest.mark.parametrize("with_data,compounds", [(False, False), (False, True), (True, True)])
def test_dialogs_never_crash(app, tmp_path, messages, with_data, compounds):
    import importlib
    win = _window(app, tmp_path, with_data, compounds)
    try:
        for target, actions in ACTIONS.items():
            mod, cls = target.split(".")
            dlg = getattr(importlib.import_module(f"spectro.ui.{mod}"), cls)(win)
            for a in actions:
                try:
                    getattr(dlg, a)()
                except Exception as exc:  # pragma: no cover - reported as failure
                    pytest.fail(f"{cls}.{a} crashed ({with_data=}, {compounds=}): "
                                f"{type(exc).__name__}: {exc}")
            dlg.close()
        from spectro.ui.dialogs_tools import ValidationDialog
        v = ValidationDialog(win)
        for name in v.inputs:
            v._calc(name)            # empty tables → message, not a crash
        v._save()
    finally:
        win.close_project()


def test_main_window_actions_without_project_or_selection(app, tmp_path, messages):
    import spectro.ui.main_window as mw
    win = mw.MainWindow()
    assert all(not a.isEnabled() for a in win.project_actions)
    win.export_figure()                     # no curves → message
    win._attach(Project.create(tmp_path / "n.spectro"))
    assert all(a.isEnabled() for a in win.project_actions)
    for fn in (win.export_spectra, win.edit_concentrations, win.process, win.set_role,
               win.move_to_trial, win.archive, win.restore, win._save_properties, win._replay):
        fn()                                 # nothing selected → message
    for mode in ("Stacked", "Difference", "Normalized", "Overlay"):
        win.view_cb.setCurrentText(mode)     # empty plot in every view
    win.close_project()
    assert all(not a.isEnabled() for a in win.project_actions)
    assert len(messages) >= 5


def test_progressive_dialog_end_to_end(app, tmp_path, messages):
    """Amplitude centering and absorption factor from the dialog on a ternary
    project: recoveries 100 %, results saved and logged."""
    from PySide6.QtWidgets import QTableWidgetItem

    import spectro.ui.main_window as mw
    from spectro.ui.dialogs_methods import ProgressiveDialog
    from tests.test_literature import AF, TER, af, ter

    win = mw.MainWindow()
    p = Project.create(tmp_path / "t.spectro")
    for c in TER:
        p.add_compound(c)
    tid = p.add_trial("Ternary")
    for c in TER:
        for v in (4, 8, 12, 16, 20, 24):
            p.add_spectrum(ter({c: v}, f"{c} {v}"), tid, role="standard")
    div = p.add_spectrum(ter({"Z": 24}, "Z′ 24"), tid, role="divisor")
    for conc in ({"X": 10, "Y": 10, "Z": 10}, {"X": 20, "Y": 10, "Z": 10}):
        p.add_spectrum(ter(conc, "mix"), tid, role="mixture")
    win._attach(p)
    try:
        d = ProgressiveDialog(win)
        d.ac_div.setCurrentIndex(d.ac_div.findData(div))
        d.ac_divc.setCurrentText("Z")
        d.ac_w.setValue(275.0)
        d.ac_plateau.setChecked(True)
        d.ac_p1.setValue(340.0)
        d.ac_p2.setValue(380.0)
        for j, v in enumerate(["X", "275", "240", "Y"]):
            d.ac_diff.setItem(0, j, QTableWidgetItem(v))
        d.ac_sub.setCurrentText("Y")
        d.ac_std.set_all(False)
        d.ac_std.populate(p.records(), {r.id for r in p.records() if r.role == "standard"})
        d._run_ac()
        assert not messages, messages
        _, data = d.last
        assert data["compounds"] == ["X", "Y", "Z"]
        for f, conc in zip(data["found"], ({"X": 10, "Y": 10, "Z": 10}, {"X": 20, "Y": 10, "Z": 10})):
            assert f == pytest.approx([conc["X"], conc["Y"], conc["Z"]], rel=5e-4)
        d._save()
        assert p.results()[-1]["kind"] == "progressive"

        # save the calibrated method, then apply it from Saved methods
        from PySide6.QtWidgets import QInputDialog

        from spectro.ui.dialogs_methods import SavedDialog
        orig = QInputDialog.getText
        QInputDialog.getText = staticmethod(lambda *a, **k: ("AAC test", True))
        try:
            d._save_method()
        finally:
            QInputDialog.getText = orig
        saved = p.methods()[-1]
        assert saved["name"] == "AAC test" and saved["type"] == "amplitude_centering"
        assert saved["definition"]["divisor"] == div and saved["definition"]["calibration_ids"]
        mix_ids = [r.id for r in p.records() if r.role == "mixture"]
        sd = SavedDialog(win)
        sd.mtable.selectRow(len(sd.methods) - 1)
        win.selected_ids = lambda: mix_ids
        sd._apply()
        assert not messages, messages
        routine = p.results()[-1]
        assert routine["kind"] == "routine" and routine["data"]["compounds"] == ["X", "Y", "Z"]
        assert routine["data"]["found"][1] == pytest.approx([20, 10, 10], rel=5e-4)
        del win.selected_ids

        # absorption factor on a second compound set
        q = Project.create(tmp_path / "a.spectro")
        for c in AF:
            q.add_compound(c)
        t2 = q.add_trial("MAFM")
        for c in AF:
            for v in (2, 6, 10, 14, 18):
                q.add_spectrum(af({c: v}, f"{c} {v}"), t2, role="standard")
        q.add_spectrum(af({"MET": 10, "MEB": 10, "DLX": 10}, "mix"), t2, role="mixture")
        win.close_project()
        win._attach(q)
        d2 = ProgressiveDialog(win)
        for i, (c, lam, qlam) in enumerate([("MET", "330", "300"), ("MEB", "280", ""),
                                            ("DLX", "246", "")]):
            for j, v in enumerate((c, lam, qlam)):
                d2.af_table.setItem(i, j, QTableWidgetItem(v))
        d2._run_af()
        assert not messages, messages
        assert d2.last[1]["found"][0] == pytest.approx([10, 10, 10], rel=5e-4)
    finally:
        win.close_project()


def test_validation_compare_from_summary_values(app, tmp_path, messages):
    """Eissa 2018 Table 4/5 typed into the dialog: t/F per method and ANOVA."""
    from PySide6.QtWidgets import QTableWidgetItem

    from spectro.ui.dialogs_tools import ValidationDialog
    win = _window(app, tmp_path)
    try:
        v = ValidationDialog(win)
        t = v.inputs["Compare (mean, SD, n)"]
        rows = [("HPLC", 100.63, 0.726, 5), ("AAS", 100.60, 0.821, 5), ("IDW", 99.80, 1.067, 5),
                ("RD", 100.01, 0.607, 5), ("DR1", 99.95, 0.735, 5), ("MCR", 99.73, 0.889, 5)]
        for i, r in enumerate(rows):
            for j, x in enumerate(r):
                t.setItem(i, j, QTableWidgetItem(str(x)))
        v._calc("Compare (mean, SD, n)")
        assert not messages, messages
        name, data = v.last
        assert data["comparisons"]["IDW"]["t"]["t"] == pytest.approx(1.438, abs=0.002)
        assert data["comparisons"]["IDW"]["F"]["F_crit"] == pytest.approx(6.388, abs=0.001)
        assert data["anova"]["F"] == pytest.approx(1.1673, rel=0.003)
        t.setItem(0, 2, QTableWidgetItem("abc"))
        v._calc("Compare (mean, SD, n)")
        assert messages and "number" in messages[-1]
    finally:
        win.close_project()


def test_optimizer_saves_a_progressive_candidate(app, tmp_path, messages):
    """A ternary project with an extended compound: the optimizer's best
    amplitude-centering candidate is saved as a calibrated method that uses
    the project's divisor spectrum."""
    import spectro.ui.main_window as mw
    from spectro.ui.dialog_optimizer import OptimizerDialog
    from tests.test_literature import TER, ter

    win = mw.MainWindow()
    p = Project.create(tmp_path / "o.spectro")
    for c in TER:
        p.add_compound(c)
    tid = p.add_trial("T")
    for c in TER:
        for v in (4, 8, 12, 16, 20, 24):
            p.add_spectrum(ter({c: v}, f"{c} {v}", noise=0.0003, seed=int(v)), tid,
                           role="standard")
    p.add_spectrum(ter({"X": 10, "Y": 10, "Z": 10}, "mix", noise=0.0003, seed=99), tid,
                   role="mixture")
    win._attach(p)
    try:
        d = OptimizerDialog(win)
        for i in range(d.fam.count()):
            it = d.fam.item(i)
            it.setCheckState(Qt.Checked if it.data(Qt.UserRole) == "amplitude_centering"
                             else Qt.Unchecked)
        d._run()
        assert not messages, messages
        row = next(i for i, c in enumerate(d.cands) if c.family == "amplitude_centering"
                   and not c.error)
        d.table.selectRow(row)
        d._explain()
        d._save()
        assert not messages, messages
        saved = p.methods()[-1]["definition"]
        assert saved["type"] == "amplitude_centering" and saved["regressions"]
        div_ids = {d.div_boxes[c].currentData() for c in TER}
        assert saved["divisor"] in div_ids
    finally:
        win.close_project()


def test_saved_methods_can_be_edited_and_renamed(app, tmp_path, messages, monkeypatch):
    """Edit… reopens each method type in its dialog with all settings; saving
    stores a new version and archives the old one."""
    from PySide6.QtWidgets import QInputDialog, QTableWidgetItem

    import spectro.ui.main_window as mw
    from spectro.ui.dialogs_methods import (ChemometricsDialog, EquationsDialog,
                                            ProgressiveDialog, SavedDialog, UnivariateDialog)
    from tests.test_literature import TER, ter
    monkeypatch.setattr(QInputDialog, "getText", lambda *a, **k: (k.get("text") or "m", True))
    win = mw.MainWindow()
    p = Project.create(tmp_path / "e.spectro")
    for c in TER:
        p.add_compound(c)
    tid = p.add_trial("T")
    for c in TER:
        for v in (4, 8, 12, 16, 20, 24):
            p.add_spectrum(ter({c: v}, f"{c} {v}"), tid, role="standard")
    div = p.add_spectrum(ter({"Z": 24}, "Z′"), tid, role="divisor")
    for i, (a, b, z) in enumerate([(10, 10, 10), (20, 10, 10), (6, 6, 18), (12, 4, 4),
                                   (8, 14, 6), (16, 12, 8), (4, 18, 12), (14, 6, 16)]):
        p.add_spectrum(ter({"X": a, "Y": b, "Z": z}, f"mix {i}"), tid,
                       role="calibration" if i < 6 else "mixture")
    win._attach(p)
    try:
        # univariate: save, edit λ, save again
        d = UnivariateDialog(win)
        d.compound.setCurrentText("Z")
        d.template.setCurrentText("Direct (zero order, λmax)")
        d.meas.set({"kind": "amplitude", "params": {"w1": 355.0}})
        d._calibrate()
        d._save_method()
        m1 = p.methods()[-1]
        e = UnivariateDialog(win)
        e.load_method(m1)
        assert e.name.text() == m1["name"] and e.meas.get()["params"]["w1"] == 355.0
        assert sorted(e.cal.checked_ids()) == sorted(m1["definition"]["calibration_ids"])
        e.meas.set({"kind": "amplitude", "params": {"w1": 360.0}})
        e._calibrate()
        e._save_method()
        live = p.methods()
        assert m1["id"] not in [m["id"] for m in live]
        new = live[-1]["definition"]
        assert new["measurement"]["params"]["w1"] == 360.0 and new["revision_of"] == m1["id"]
        assert p.method(m1["id"])["archived"]

        # equations
        q = EquationsDialog(win)
        q.load_method({"id": 0, "name": "x", "definition": {
            "type": "equations", "compounds": ["X", "Y", "Z"], "steps": [], "intercept": False,
            "signals": [{"kind": "amplitude", "params": {"w1": 240.0}},
                        {"kind": "area", "params": {"w1": 285.0, "w2": 295.0}},
                        {"kind": "amplitude", "params": {"w1": 360.0}}]}})
        assert [s_["kind"] for s_ in q._signals()] == ["amplitude", "area", "amplitude"]

        # progressive (with smoothing): save, edit, re-save
        g = ProgressiveDialog(win)
        g.ac_div.setCurrentIndex(g.ac_div.findData(div))
        g.ac_divc.setCurrentText("Z")
        g.ac_w.setValue(275.0)
        g.ac_plateau.setChecked(True)
        g.ac_p1.setValue(340.0)
        g.ac_p2.setValue(380.0)
        for j, v in enumerate(["X", "275", "240", "Y"]):
            g.ac_diff.setItem(0, j, QTableWidgetItem(v))
        g.ac_sub.setCurrentText("Y")
        g.sm_win.setValue(15)
        g._run_ac()
        g._save_method()
        pm = p.methods()[-1]
        assert pm["definition"]["steps"][0]["params"]["window"] == 15
        h = ProgressiveDialog(win)
        h.load_method(pm)
        assert h.sm_win.value() == 15 and h.ac_w.value() == 275.0
        assert h.ac_diff.item(0, 3).text() == "Y" and h.ac_plateau.isChecked()
        h.sm_win.setValue(21)
        h._run_ac()
        h._save_method()
        assert p.methods()[-1]["definition"]["steps"][0]["params"]["window"] == 21
        assert p.methods()[-1]["definition"]["version"] == 2

        # chemometrics
        c = ChemometricsDialog(win)
        c.mtype.setCurrentText("CLS")
        c._fit()
        c._save_method()
        cm = p.methods()[-1]
        c2 = ChemometricsDialog(win)
        c2.load_method(cm)
        assert c2.mtype.currentText() == "CLS"
        assert sorted(c2.cal.checked_ids()) == sorted(cm["definition"]["calibration_ids"])

        # Saved methods: rename and the Edit… entry point
        sd = SavedDialog(win)
        sd.mtable.selectRow(0)
        before = sd.methods[0]
        monkeypatch.setattr(QInputDialog, "getText", lambda *a, **k: ("renamed", True))
        sd._rename()
        assert any(m["name"] == "renamed" for m in p.methods())
        assert before["id"] not in [m["id"] for m in p.methods()]
        monkeypatch.setattr(UnivariateDialog, "exec", lambda self: 0)
        monkeypatch.setattr(ChemometricsDialog, "exec", lambda self: 0)
        monkeypatch.setattr(ProgressiveDialog, "exec", lambda self: 0)
        for row in range(len(sd.methods)):
            sd.mtable.selectRow(row)
            sd._edit()
            assert "editing method" in sd.editor.windowTitle()
        assert not messages, messages
    finally:
        win.close_project()
