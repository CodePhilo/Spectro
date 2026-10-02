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


def test_view_modes_and_figure_export(win, tmp_path, monkeypatch):
    from spectro.ui.figure import FigureDialog, collect_series, render
    recs = win.project.records()
    select(win, [r.id for r in recs if r.role == "standard"][:4])
    for mode in ("Stacked", "Difference", "Normalized", "Overlay"):
        win.view_cb.setCurrentText(mode)
        assert len(win.plot.curves) >= 3
    win.view_cb.setCurrentText("Stacked")
    ys = [c.values.max() for _, c in win.plot.curves]
    assert ys == sorted(ys)
    data = collect_series(win.plot)
    assert len(data["series"]) == 4 and data["xlabel"].startswith("Wavelength")
    dlg = FigureDialog(win.plot, win)
    opts = dlg.options()
    for ext in ("png", "tif", "svg", "pdf"):
        fig = render(data, {**opts, "black_white": ext == "pdf"}, dpi=300)
        fig.savefig(tmp_path / f"f.{ext}", dpi=300)
        assert (tmp_path / f"f.{ext}").stat().st_size > 1000
    from PIL import Image
    with Image.open(tmp_path / "f.tif") as im:
        assert abs(im.size[0] - round(85 / 25.4 * 300)) <= 2


def test_table_export_and_results_excel(win, tmp_path, monkeypatch):
    import openpyxl
    from spectro.ui import widgets
    from spectro.ui.widgets import PasteTable, fill_table
    t = PasteTable()
    fill_table(t, ["A", "B"], [["x", 1.23456789]])
    monkeypatch.setattr(widgets.QFileDialog, "getSaveFileName",
                        lambda *a, **k: (str(tmp_path / "t.xlsx"), ""))
    t.export_excel()
    assert openpyxl.load_workbook(tmp_path / "t.xlsx").active["B2"].value == 1.23456789
    assert win.project.audit_entries(action="EXPORT")


def test_standard_addition_and_robustness_ui(win, monkeypatch):
    from spectro.ui.dialogs_methods import RobustnessDialog, UnivariateDialog
    p = win.project
    tid = p.trials()[0]["id"]
    spikes = []
    for i, a in enumerate([0, 0, 2, 4, 6]):
        s = mixture({"X": 5 + a, "Y": 8}, f"spiked {i}")
        s.concentrations = {"X": float(a)}
        spikes.append(p.add_spectrum(s, tid, role="sample"))
    d = UnivariateDialog(win)
    d.compound.setCurrentText("X")
    d.template.setCurrentText("Ratio difference (RD)")
    ydiv = next(r for r in p.records() if r.name == "Y 20")
    st = d.pipe.get_steps()
    st[0]["params"]["reference"] = ydiv.id
    d.pipe.set_steps(st)
    d.meas.form.set_value("w1", 240.0)
    d.meas.form.set_value("w2", 260.0)
    d._calibrate()
    d.test.set_all(False)
    for i in range(d.test.count()):
        if d.test.item(i).data(Qt.UserRole) in spikes:
            d.test.item(i).setCheckState(Qt.Checked)
    d._std_addition()
    kind, data = d.last[1], d.last[2]
    assert kind == "standard_addition" and abs(data["sample"] - 5) < 0.1
    assert abs(data["mean_recovery"] - 100) < 2
    d._save_results()
    rob = RobustnessDialog(d, d.method, p.spectra(d.cal_ids), p.spectra(spikes[2:]))
    rob._run()
    assert rob.result_data and len(rob.result_data["variants"]) == 5


def test_chemometrics_limits_outliers_and_models(win, monkeypatch, tmp_path):
    from spectro.ui import dialogs_methods as dm
    c = dm.ChemometricsDialog(win)
    c.mtype.setCurrentText("PLS2")
    c.ncomp.setValue(3)
    c._fit()
    assert c.model.is_fitted and c.diag_table.rowCount() > 0
    c.conf.setCurrentText("99")
    monkeypatch.setattr(dm.QInputDialog, "getText", lambda *a, **k: ("PLS model", True))
    c._save_method()
    saved = [m for m in win.project.methods() if m["name"] == "PLS model"][0]
    assert saved["definition"]["fitted"]
    c._flagged = [c.cal.item(0).data(Qt.UserRole)]
    from PySide6.QtWidgets import QMessageBox
    monkeypatch.setattr(QMessageBox, "question", lambda *a, **k: QMessageBox.Yes)
    c._exclude_outliers()
    assert c.cal.item(0).checkState() == Qt.Unchecked
    c._cv()
    s = dm.SavedDialog(win)
    s.rtable.selectRow(0)
    s._show_result()
    path = tmp_path / "m.spmodel"
    monkeypatch.setattr(dm.QFileDialog, "getSaveFileName", lambda *a, **k: (str(path), ""))
    monkeypatch.setattr(dm.QFileDialog, "getOpenFileName", lambda *a, **k: (str(path), ""))
    s.mtable.selectRow(len(s.methods) - 1)
    s._export_model()
    n = len(win.project.methods())
    s._import_model()
    assert len(win.project.methods()) == n + 1


def test_new_binary_tabs(win):
    from spectro.ui.dialogs_methods import SpecialDialog
    d = SpecialDialog(win)
    recs = win.project.records()

    def check(lst, pred):
        for i in range(lst.count()):
            rid = lst.item(i).data(Qt.UserRole)
            r = next(x for x in recs if x.id == rid)
            lst.item(i).setCheckState(Qt.Checked if pred(r) else Qt.Unchecked)
    x_std = lambda r: r.name.startswith("X ")  # noqa: E731
    y_std = lambda r: r.name.startswith("Y ")  # noqa: E731
    mix = lambda r: r.name.startswith("mix")  # noqa: E731
    d.iam_y.setCurrentText("Y")
    d.iam_w.setValue(250)
    d.iam_p1.setValue(340)
    d.iam_p2.setValue(360)
    check(d.iam_xs, x_std)
    check(d.iam_ys, y_std)
    check(d.iam_mix, mix)
    d._run_iam()
    found = d.last[1]["found"]
    assert abs(found[0]["X"] - 6) < 0.15 and abs(found[0]["Y"] - 10) < 0.15
    d.aas_y.setCurrentText("Y")
    d.aas_w2.setValue(245)
    check(d.aas_xs, x_std)
    check(d.aas_ys, y_std)
    check(d.aas_mix, mix)
    from spectro.core import univariate as uv
    d.aas_w1.setValue(uv.equal_amplitude_wavelengths(win.project.spectrum(
        next(r.id for r in recs if r.name == "Y 20")), 245.0)[0])
    d._run_aas()
    found = d.last[1]["found"]
    assert abs(found[0]["X"] - 6) < 0.15 and abs(found[0]["Y"] - 10) < 0.15


def test_default_selection_excludes_derived(win):
    from spectro.ui.dialogs_methods import checklists
    p = win.project
    std = [r.id for r in p.records() if r.role == "standard"]
    derived = p.process(std, [{"op": "derivative", "params": {"order": 1}}])
    cal, test = checklists(p, "X")
    assert not set(cal.checked_ids()) & set(derived)
    assert set(cal.checked_ids()) <= set(std)


def test_method_optimizer_dialog(win):
    from spectro.ui.dialog_optimizer import OptimizerDialog
    d = OptimizerDialog(win)
    for i in range(d.comps.count()):  # X + Y only
        d.comps.item(i).setCheckState(Qt.Checked if d.comps.item(i).text() in ("X", "Y")
                                      else Qt.Unchecked)
    assert all(len(d.std_ids[c]) == 5 for c in ("X", "Y"))
    d._run()
    assert d.result and d.table.rowCount() > 5
    assert "Recommendations" in d.summary.text()
    best = d.cands[0]
    assert best.score < 3
    assert d.smooth.isChecked() and d.inp.smoothing == (3.0, 6.0, 10.0)
    assert any("SG smoothing" in c.label for c in d.cands)
    d.smooth_w.setText("6, x")
    try:
        d._build_input()
        raise AssertionError("bad widths accepted")
    except ValueError as exc:
        assert "smoothing widths" in str(exc)
    d.smooth_w.setText("3, 6, 10")
    d.table.selectRow(0)
    d._explain()
    n = len(win.project.methods())
    d._save()
    assert len(win.project.methods()) == n + 1
    uni = next(i for i, c in enumerate(d.cands) if c.measurement)
    d.table.selectRow(uni)
    d._save()
    saved = win.project.methods()[-1]["definition"]
    assert saved["type"] == "univariate" and saved["regression"]["r"] > 0.999
    from PySide6.QtWidgets import QDialog
    orig = QDialog.exec
    QDialog.exec = lambda self: 0
    try:
        d._plan()
    finally:
        QDialog.exec = orig
    plan = {r[0]: r for r in d.plan_rows}
    assert "levels" in plan["Linearity"][1] and "Robustness" in plan
    assert "± 1 nm" in plan["Robustness"][1] or "± 2 points" in plan["Robustness"][1] \
        or "nm" in plan["Robustness"][1]


def test_optimizer_ignores_mixtures_with_unselected_compounds(win):
    from spectro.core.spectrum import Spectrum
    from spectro.ui.dialog_optimizer import OptimizerDialog
    p = win.project
    tern = mixture({"X": 6, "Y": 8, "Z": 5}, "ternary mix")
    tid = p.trials()[0]["id"]
    sid = p.add_spectrum(Spectrum(tern.x, tern.y, name="ternary mix",
                                  concentrations=tern.concentrations), tid, role="mixture")
    d = OptimizerDialog(win)
    for i in range(d.comps.count()):
        d.comps.item(i).setCheckState(Qt.Checked if d.comps.item(i).text() in ("X", "Y")
                                      else Qt.Unchecked)
    ids = [d.mix.item(i).data(Qt.UserRole) for i in range(d.mix.count())]
    assert sid not in ids and ids
