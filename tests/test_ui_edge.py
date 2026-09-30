"""Every dialog and action on empty / incomplete projects must show a message,
never crash."""

import pytest

pytest.importorskip("PySide6")

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
