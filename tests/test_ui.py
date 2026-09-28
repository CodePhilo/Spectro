"""Headless smoke tests of the desktop UI (QT_QPA_PLATFORM=offscreen)."""

import pytest

pytest.importorskip("PySide6")

from PySide6.QtCore import Qt  # noqa: E402
from PySide6.QtWidgets import QApplication  # noqa: E402

from spectro.core.io import export_excel  # noqa: E402
from spectro.core.multicomponent import design_concentrations  # noqa: E402
from spectro.storage.project import Project  # noqa: E402
from tests.conftest import mixture  # noqa: E402


@pytest.fixture(scope="module")
def app():
    return QApplication.instance() or QApplication([])


@pytest.fixture
def win(app, tmp_path, monkeypatch):
    import spectro.ui.dialogs_data as dd
    import spectro.ui.main_window as mw
    import spectro.ui.widgets as w

    for mod in (w, dd, mw):
        monkeypatch.setattr(mod, "ask_reason", lambda *a, **k: "test", raising=False)
    monkeypatch.setattr(mw.QMessageBox, "information", lambda *a, **k: None)
    monkeypatch.setattr(mw.QMessageBox, "critical", lambda *a, **k: None)
    monkeypatch.setattr(w.QMessageBox, "warning",
                        lambda *a, **k: (_ for _ in ()).throw(AssertionError(a[2])))
    win = mw.MainWindow()
    proj = Project.create(tmp_path / "ui.spectro")
    for c in ("X", "Y", "Z"):
        proj.add_compound(c)
    tid = proj.add_trial("Ternary")
    stds = [mixture({"X": c, "Y": 0, "Z": 0}, f"X {c}") for c in (4, 8, 12, 16, 20)]
    stds += [mixture({"X": 0, "Y": c, "Z": 0}, f"Y {c}") for c in (4, 8, 12, 16, 20)]
    design = [mixture(c, f"D{i}", 0.0005, i) for i, c in enumerate(
        design_concentrations({"X": 10, "Y": 10, "Z": 10}, {"X": 4, "Y": 4, "Z": 4}))]
    mixes = [mixture({"X": a, "Y": b, "Z": 0}, f"mix {a}+{b}") for a, b in [(6, 10), (12, 8)]]
    for group, role in ((stds, "standard"), (design, "calibration"), (mixes, "mixture")):
        f = tmp_path / f"{role}.xlsx"
        export_excel(group, f)
        ids = proj.import_file(f, tid, role=role)
        proj.set_concentrations_bulk({i: s.concentrations for i, s in zip(ids, group)})
    win._attach(proj)
    yield win
    win.close_project()


def select(win, ids):
    win.tree.clearSelection()
    it = win.tree.invisibleRootItem()
    stack = [it]
    while stack:
        node = stack.pop()
        for i in range(node.childCount()):
            ch = node.child(i)
            if ch.data(0, Qt.UserRole) in ids:
                ch.setSelected(True)
            stack.append(ch)


def test_tree_properties_and_process(win):
    p = win.project
    recs = p.records()
    assert win.tree.topLevelItemCount() == 1
    select(win, [recs[0].id])
    assert win.p_name.text() == recs[0].name
    from spectro.ui.dialogs_data import ProcessDialog
    ydiv = next(r for r in recs if r.name == "Y 20")
    select(win, [r.id for r in recs[:3]])
    dlg = ProcessDialog(win, [r.id for r in recs[:3]])
    dlg.editor.set_steps([{"op": "divide", "params": {"reference": ydiv.id}},
                          {"op": "derivative", "params": {"order": 1}}])
    dlg._preview()
    assert "⚠" not in dlg.msg.text()
    dlg._apply()
    assert len([r for r in p.records() if r.kind == "derived"]) == 3


def test_univariate_dialog(win):
    from spectro.ui.dialogs_methods import UnivariateDialog
    d = UnivariateDialog(win)
    d.compound.setCurrentText("Y")
    d.template.setCurrentText("Direct (zero order, λmax)")
    d.meas.form.set_value("w1", 350.0)
    d._calibrate()
    assert d.method.regression.r > 0.999
    d.test.set_all(False)
    for i in range(d.test.count()):
        if "mix" in d.test.item(i).text():
            d.test.item(i).setCheckState(Qt.Checked)
    d._predict()
    assert d.predictions and all(abs(r - 100) < 2 for r in d.predictions["recovery"])
    d._save_method()
    d._save_results()


def test_equations_and_chemometrics(win):
    from spectro.ui.dialogs_methods import ChemometricsDialog, EquationsDialog
    e = EquationsDialog(win)
    e.comps.item(2).setCheckState(Qt.Unchecked)  # X, Y only
    from spectro.ui.widgets import fill_table
    fill_table(e.signals, ["Kind", "λ1", "λ2"], [["amplitude", "245", ""],
                                               ["amplitude", "330", ""]], editable=True)
    e.cal.set_all(False)
    for i in range(e.cal.count()):
        if "[standard]" in e.cal.item(i).text():
            e.cal.item(i).setCheckState(Qt.Checked)
    e._fit()
    e._predict()
    assert e.predictions
    c = ChemometricsDialog(win)
    c.mtype.setCurrentText("PLS2")
    c.ncomp.setValue(3)
    c.cv.setCurrentText("venetian")
    c._cv()
    c._fit()
    assert c.model is not None


def test_tools_and_report(win, tmp_path):
    from spectro.ui.dialogs_methods import SpecialDialog
    from spectro.ui.dialogs_tools import (DesignDialog, FinderDialog, GreennessDialog,
                                          ValidationDialog)
    from spectro.ui.report import build_html
    f = FinderDialog(win)
    f._zero()
    f._extrema()
    DesignDialog(win)._gen()
    g = GreennessDialog(win)
    g._eco()
    g._gapi()
    v = ValidationDialog(win)
    t = v.inputs["Linearity"]
    from PySide6.QtWidgets import QTableWidgetItem
    for i, (x, y) in enumerate([(2, 0.1), (4, 0.21), (6, 0.3), (8, 0.41), (10, 0.5)]):
        t.setItem(i, 0, QTableWidgetItem(str(x)))
        t.setItem(i, 1, QTableWidgetItem(str(y)))
    v._calc("Linearity")
    SpecialDialog(win)
    doc = build_html(win.project, win.project.trials()[0]["id"], include_plot=True)
    assert "Audit trail" in doc and "PASSED" in doc
