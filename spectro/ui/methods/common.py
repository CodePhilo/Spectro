"""Helpers shared by the method dialogs."""

from __future__ import annotations

import pyqtgraph as pg
from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QComboBox,
    QDoubleSpinBox,
    QGroupBox,
    QHBoxLayout,
    QListWidget,
    QListWidgetItem,
    QPushButton,
    QScrollArea,
    QVBoxLayout,
    QWidget,
)

from spectro.core.validation import describe
from spectro.ui.widgets import TEXT_SECONDARY, SpectrumChecklist

CAL_ROLES = {"standard"}                    # univariate / equation methods
TRAIN_ROLES = {"calibration", "standard"}   # chemometric training set
TEST_ROLES = {"mixture", "sample", "validation", "unknown"}


def scroll(widget: QWidget) -> QScrollArea:
    sa = QScrollArea()
    sa.setWidgetResizable(True)
    sa.setWidget(widget)
    return sa


def xy_plot(x_label: str, y_label: str) -> pg.PlotWidget:
    p = pg.PlotWidget()
    p.setLabel("bottom", x_label)
    p.setLabel("left", y_label)
    p.showGrid(x=True, y=True, alpha=0.12)
    p.addLegend(offset=(10, 10), labelTextColor=TEXT_SECONDARY)
    from spectro.ui.figure import install_figure_export
    install_figure_export(p)
    return p


def scatter(p: pg.PlotWidget, x, y, colour: str, name: str | None = None, line=None,
            labels: list[str] | None = None) -> None:
    """Markers ≥ 8 px with a surface-coloured ring; hover shows the label."""
    kw = dict(x=list(x), y=list(y), size=9, pen=pg.mkPen(SURFACE_RING, width=2),
              brush=pg.mkBrush(colour), name=name, hoverable=True)
    if labels:
        kw.update(data=list(labels), tip=lambda x, y, data: f"{data}\n({x:.4g}, {y:.4g})")
    p.addItem(pg.ScatterPlotItem(**kw))
    if line is not None:
        lx, ly = line
        p.plot(lx, ly, pen=pg.mkPen(colour, width=2))


SURFACE_RING = "#fcfcfb"


def checklists(project, compound: str | None = None, cal_roles=CAL_ROLES, test_roles=TEST_ROLES,
               trial: int | None = None):
    recs = project.records(trial_id=trial)
    cal, test = SpectrumChecklist(), SpectrumChecklist()
    # defaults use raw spectra only: methods apply their own processing, and
    # derived spectra inherit the role and concentrations of their parent
    raw = [r for r in recs if r.kind == "raw"]
    cal_ids = {r.id for r in raw if r.role in cal_roles and
               (compound is None or r.concentrations.get(compound, 0) > 0)}
    if cal_roles is TRAIN_ROLES and any(r.role == "calibration" for r in raw):
        cal_ids = {r.id for r in raw if r.role == "calibration"}
    test_ids = {r.id for r in raw if r.role in test_roles}
    cal.populate(recs, cal_ids)
    test.populate(recs, test_ids)
    return cal, test


def add_grouped(combo: QComboBox, groups: dict[str, list[str]]) -> None:
    """Items under bold, non-selectable family headings."""
    from PySide6.QtGui import QFont

    from spectro.core.catalog import CATEGORIES
    model = combo.model()
    for key, names in groups.items():
        combo.addItem(f"— {CATEGORIES.get(key, key)} —")
        head = model.item(combo.count() - 1)
        head.setEnabled(False)
        f = QFont(head.font())
        f.setBold(True)
        head.setFont(f)
        for n in names:
            combo.addItem(n)


def check_ids(lst: SpectrumChecklist, ids) -> None:
    """Check exactly the spectra in ``ids`` (unknown ids are ignored)."""
    ids = {int(i) for i in ids or []}
    for i in range(lst.count()):
        it = lst.item(i)
        it.setCheckState(Qt.Checked if it.data(Qt.UserRole) in ids else Qt.Unchecked)


def store_method(dlg, name: str, definition: dict) -> int | None:
    """Save a new method, or — when the dialog was opened with *Edit…* from
    Saved methods — a new version of the edited one (old version archived,
    reason required)."""
    from spectro.ui.widgets import ask_reason
    ed = getattr(dlg, "editing", None)
    if ed:
        reason = ask_reason(dlg, f"Save the changes to method #{ed['id']} '{ed['name']}' as a "
                                 "new version? The current version is kept in the archive.",
                            required=True)
        if not reason:
            return None
        mid = dlg.project.revise_method(ed["id"], name, definition, reason)
        dlg.editing = dlg.project.method(mid)
        dlg.setWindowTitle(dlg.windowTitle().split(" — editing")[0]
                           + f" — editing method #{mid} '{name}'")
    else:
        mid = dlg.project.save_method(name, definition, dlg.win.current_trial())
    dlg.method_id = mid
    dlg.win.statusBar().showMessage(f"Method #{mid} saved.")
    return mid


def start_editing(dlg, md: dict) -> None:
    dlg.editing = md
    dlg.setWindowTitle(dlg.windowTitle() + f" — editing method #{md['id']} '{md['name']}'")


def list_box(title: str, lst: SpectrumChecklist) -> QGroupBox:
    g = QGroupBox(title)
    v = QVBoxLayout(g)
    row = QHBoxLayout()
    for t, state in (("All", True), ("None", False)):
        b = QPushButton(t)
        b.clicked.connect(lambda _=False, s=state: lst.set_all(s))
        row.addWidget(b)
    row.addStretch(1)
    v.addLayout(row)
    v.addWidget(lst)
    lst.setMinimumHeight(130)
    return g


def results_rows(names, found, taken):
    rows, recs = [], []
    for n, f, t in zip(names, found, taken):
        rec = 100 * f / t if t not in (None, 0) else None
        if rec is not None:
            recs.append(rec)
        rows.append([n, float(f), "" if t is None else float(t),
                     "" if rec is None else float(rec)])
    return rows, recs


def summary_text(recs: list[float]) -> str:
    if len(recs) < 2:
        return ""
    d = describe(recs)
    return f"Recovery: mean {d['mean']:.2f} %, SD {d['sd']:.3f}, RSD {d['rsd']:.2f} % (n = {d['n']})"




class CompoundChecks(QListWidget):
    def __init__(self, names: list[str]):
        super().__init__()
        for n in names:
            it = QListWidgetItem(n)
            it.setFlags(it.flags() | Qt.ItemIsUserCheckable)
            it.setCheckState(Qt.Checked)
            self.addItem(it)
        self.setMaximumHeight(90)

    def checked(self) -> list[str]:
        return [self.item(i).text() for i in range(self.count())
                if self.item(i).checkState() == Qt.Checked]


def wl_spin(v: float) -> QDoubleSpinBox:
    s = QDoubleSpinBox()
    s.setRange(0, 100000)
    s.setDecimals(2)
    s.setSuffix(" nm")
    s.setValue(v)
    return s
