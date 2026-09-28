"""Reusable Qt widgets: parameter forms, pipeline editor, spectrum plot, tables."""

from __future__ import annotations

import math
from typing import Any, Callable

import numpy as np
import pyqtgraph as pg
from PySide6.QtCore import Qt, Signal
from PySide6.QtGui import QGuiApplication, QKeySequence
from PySide6.QtWidgets import (QAbstractItemView, QCheckBox, QComboBox, QDoubleSpinBox,
                               QFormLayout, QHBoxLayout, QHeaderView, QLabel, QLineEdit,
                               QListWidget, QListWidgetItem, QMessageBox, QPushButton,
                               QSpinBox, QTableWidget, QTableWidgetItem, QVBoxLayout, QWidget,
                               QInputDialog)

from spectro.core.operations import REGISTRY, Param, categories, describe_step
from spectro.core.spectrum import Spectrum
from spectro.core.univariate import MEASUREMENTS

# Validated categorical order (identity), light surface.
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"]
# Single-hue ramp used when more than eight spectra are overlaid.
RAMP = ["#86b6ef", "#6da7ec", "#5598e7", "#3987e5", "#2a78d6", "#256abf", "#1c5cab",
        "#184f95", "#104281", "#0d366b"]
SURFACE = "#fcfcfb"
TEXT_SECONDARY = "#52514e"
GRID_ALPHA = 0.12

pg.setConfigOptions(antialias=True, background=SURFACE, foreground=TEXT_SECONDARY)


def colours_for(n: int) -> list[str]:
    if n <= len(SERIES):
        return SERIES[:n]
    idx = np.linspace(0, len(RAMP) - 1, n).round().astype(int)
    return [RAMP[i] for i in idx]


def fmt(v: Any, digits: int = 4) -> str:
    if v is None:
        return ""
    if isinstance(v, float):
        if math.isnan(v):
            return "–"
        if v != 0 and (abs(v) < 1e-3 or abs(v) >= 1e5):
            return f"{v:.{digits}e}"
        return f"{v:.{digits}f}"
    return str(v)


def error(parent: QWidget, exc: Exception | str, title: str = "Spectro") -> None:
    QMessageBox.warning(parent, title, str(exc))


def ask_reason(parent: QWidget, what: str, required: bool = False) -> str | None:
    """Prompt for a reason for change (recorded in the audit trail)."""
    while True:
        text, ok = QInputDialog.getText(parent, "Reason for change",
                                        f"{what}\n\nReason (recorded in the audit trail):")
        if not ok:
            return None
        if text.strip() or not required:
            return text.strip()


# --------------------------------------------------------------------------- #
# Parameter form (built from Param declarations)
# --------------------------------------------------------------------------- #
class ParamForm(QWidget):
    changed = Signal()

    def __init__(self, params: list[Param], values: dict | None = None,
                 spectrum_choices: Callable[[], list[tuple[Any, str]]] | None = None,
                 parent: QWidget | None = None):
        super().__init__(parent)
        self.params = params
        self.spectrum_choices = spectrum_choices or (lambda: [])
        self.editors: dict[str, QWidget] = {}
        lay = QFormLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        values = values or {}
        for p in params:
            ed = self._editor(p, values.get(p.name, p.default))
            if p.help:
                ed.setToolTip(p.help)
            self.editors[p.name] = ed
            lay.addRow(p.label, ed)
        if not params:
            lay.addRow(QLabel("No parameters."))

    def _editor(self, p: Param, value: Any) -> QWidget:
        if p.kind == "spectrum":
            cb = QComboBox()
            cb.addItem("(choose)" if not p.optional else "(none)", None)
            for key, label in self.spectrum_choices():
                cb.addItem(label, key)
            if value is not None:
                i = cb.findData(value)
                if i >= 0:
                    cb.setCurrentIndex(i)
            cb.currentIndexChanged.connect(self.changed)
            return cb
        if p.kind == "choice":
            cb = QComboBox()
            cb.addItems([str(c) for c in p.choices])
            if value is not None:
                cb.setCurrentText(str(value))
            cb.currentIndexChanged.connect(self.changed)
            return cb
        if p.kind == "bool":
            ch = QCheckBox()
            ch.setChecked(bool(value))
            ch.toggled.connect(self.changed)
            return ch
        if p.optional:
            le = QLineEdit("" if value is None else f"{value:g}")
            le.setPlaceholderText("auto / blank")
            le.editingFinished.connect(self.changed)
            return le
        if p.kind == "int":
            sb = QSpinBox()
            sb.setRange(int(p.minimum) if p.minimum is not None else -10 ** 6,
                        int(p.maximum) if p.maximum is not None else 10 ** 6)
            sb.setValue(int(value if value is not None else 0))
            sb.valueChanged.connect(self.changed)
            return sb
        ds = QDoubleSpinBox()
        ds.setDecimals(6 if p.kind == "float" else 2)
        ds.setRange(p.minimum if p.minimum is not None else -1e12,
                    p.maximum if p.maximum is not None else 1e12)
        ds.setSingleStep(1.0 if p.kind == "wavelength" else 0.1)
        if p.kind == "wavelength":
            ds.setSuffix(" nm")
        ds.setValue(float(value if value is not None else 0.0))
        ds.valueChanged.connect(self.changed)
        return ds

    def values(self) -> dict[str, Any]:
        out: dict[str, Any] = {}
        for p in self.params:
            ed = self.editors[p.name]
            if isinstance(ed, QComboBox) and p.kind == "spectrum":
                out[p.name] = ed.currentData()
            elif isinstance(ed, QComboBox):
                out[p.name] = ed.currentText()
            elif isinstance(ed, QCheckBox):
                out[p.name] = ed.isChecked()
            elif isinstance(ed, QLineEdit):
                t = ed.text().strip().replace(",", ".")
                out[p.name] = None if not t else (int(float(t)) if p.kind == "int" else float(t))
            elif isinstance(ed, QSpinBox):
                out[p.name] = ed.value()
            else:
                out[p.name] = ed.value()
        return out

    def set_value(self, name: str, value: Any) -> None:
        ed = self.editors.get(name)
        if ed is None:
            return
        ed.blockSignals(True)
        if isinstance(ed, QDoubleSpinBox):
            ed.setValue(float(value))
        elif isinstance(ed, QSpinBox):
            ed.setValue(int(value))
        elif isinstance(ed, QLineEdit):
            ed.setText(f"{value:g}")
        ed.blockSignals(False)
        self.changed.emit()


# --------------------------------------------------------------------------- #
# Pipeline editor
# --------------------------------------------------------------------------- #
class PipelineEditor(QWidget):
    """Ordered list of processing steps with a parameter form per step."""

    changed = Signal()

    def __init__(self, spectrum_choices: Callable[[], list[tuple[Any, str]]],
                 name_of: Callable[[Any], str] | None = None, parent: QWidget | None = None):
        super().__init__(parent)
        self.spectrum_choices = spectrum_choices
        self.name_of = name_of
        self.steps: list[dict] = []
        lay = QVBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        add_row = QHBoxLayout()
        self.op_combo = QComboBox()
        for cat, ops in categories().items():
            self.op_combo.addItem(f"— {cat} —", None)
            idx = self.op_combo.count() - 1
            self.op_combo.model().item(idx).setEnabled(False)
            for op in ops:
                self.op_combo.addItem(f"   {op.label}", op.key)
        self.op_combo.setCurrentIndex(1)
        add_btn = QPushButton("Add step")
        add_btn.clicked.connect(self._add)
        add_row.addWidget(self.op_combo, 1)
        add_row.addWidget(add_btn)
        lay.addLayout(add_row)
        self.list = QListWidget()
        self.list.setMaximumHeight(130)
        self.list.currentRowChanged.connect(self._select)
        lay.addWidget(self.list)
        btns = QHBoxLayout()
        for text, fn in (("▲", lambda: self._move(-1)), ("▼", lambda: self._move(1)),
                         ("Remove", self._remove)):
            b = QPushButton(text)
            b.clicked.connect(fn)
            btns.addWidget(b)
        btns.addStretch(1)
        lay.addLayout(btns)
        self.desc = QLabel()
        self.desc.setWordWrap(True)
        self.desc.setStyleSheet(f"color: {TEXT_SECONDARY};")
        lay.addWidget(self.desc)
        self.form_holder = QVBoxLayout()
        lay.addLayout(self.form_holder)
        lay.addStretch(1)
        self.form: ParamForm | None = None

    def set_steps(self, steps: list[dict]) -> None:
        self.steps = [{"op": s["op"], "params": dict(s.get("params", {}))} for s in steps]
        self._refresh(len(self.steps) - 1)
        self.changed.emit()

    def get_steps(self) -> list[dict]:
        self._store()
        return [{"op": s["op"], "params": dict(s["params"])} for s in self.steps]

    def _add(self) -> None:
        key = self.op_combo.currentData()
        if key is None:
            return
        self._store()
        self.steps.append({"op": key, "params": REGISTRY[key].defaults()})
        self._refresh(len(self.steps) - 1)
        self.changed.emit()

    def _remove(self) -> None:
        r = self.list.currentRow()
        if r < 0:
            return
        self.steps.pop(r)
        self._drop_form()
        self._refresh(min(r, len(self.steps) - 1))
        self.changed.emit()

    def _move(self, d: int) -> None:
        r = self.list.currentRow()
        j = r + d
        if r < 0 or not 0 <= j < len(self.steps):
            return
        self._store()
        self.steps[r], self.steps[j] = self.steps[j], self.steps[r]
        self._drop_form()
        self._refresh(j)
        self.changed.emit()

    def _refresh(self, select: int) -> None:
        self.list.blockSignals(True)
        self.list.clear()
        for i, s in enumerate(self.steps):
            self.list.addItem(f"{i + 1}. {describe_step(s, self.name_of)}")
        self.list.blockSignals(False)
        self._drop_form()
        if select >= 0:
            self.list.setCurrentRow(select)
        self._select(select)

    def _store(self) -> None:
        r = self.list.currentRow()
        if self.form is not None and 0 <= r < len(self.steps):
            self.steps[r]["params"] = self.form.values()
            item = self.list.item(r)
            if item:
                item.setText(f"{r + 1}. {describe_step(self.steps[r], self.name_of)}")

    def _drop_form(self) -> None:
        if self.form is not None:
            self.form.setParent(None)
            self.form.deleteLater()
            self.form = None

    def _select(self, row: int) -> None:
        self._drop_form()
        if not 0 <= row < len(self.steps):
            self.desc.setText("")
            return
        step = self.steps[row]
        op = REGISTRY[step["op"]]
        self.desc.setText(op.description)
        self.form = ParamForm(op.params, step["params"], self.spectrum_choices)
        self.form.changed.connect(self._form_changed)
        self.form_holder.addWidget(self.form)

    def _form_changed(self) -> None:
        self._store()
        self.changed.emit()


class MeasurementEditor(QWidget):
    changed = Signal()

    def __init__(self, spectrum_choices, parent: QWidget | None = None):
        super().__init__(parent)
        self.spectrum_choices = spectrum_choices
        lay = QVBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        self.kind = QComboBox()
        for m in MEASUREMENTS.values():
            self.kind.addItem(m.label, m.key)
        self.kind.currentIndexChanged.connect(self._rebuild)
        lay.addWidget(self.kind)
        self.desc = QLabel()
        self.desc.setWordWrap(True)
        self.desc.setStyleSheet(f"color: {TEXT_SECONDARY};")
        lay.addWidget(self.desc)
        self.holder = QVBoxLayout()
        lay.addLayout(self.holder)
        self.form: ParamForm | None = None
        self._values: dict = {}
        self._rebuild()

    def _rebuild(self) -> None:
        if self.form is not None:
            self._values.update({k: v for k, v in self.form.values().items() if v is not None})
            self.form.setParent(None)
            self.form.deleteLater()
        m = MEASUREMENTS[self.kind.currentData()]
        self.desc.setText(m.description)
        self.form = ParamForm(m.params, {p.name: self._values.get(p.name, p.default)
                                         for p in m.params}, self.spectrum_choices)
        self.form.changed.connect(self.changed)
        self.holder.addWidget(self.form)
        self.changed.emit()

    def get(self) -> dict:
        params = {k: v for k, v in self.form.values().items() if v is not None}
        return {"kind": self.kind.currentData(), "params": params}

    def set(self, spec: dict) -> None:
        self._values = dict(spec.get("params", {}))
        i = self.kind.findData(spec["kind"])
        self.kind.blockSignals(True)
        self.kind.setCurrentIndex(max(i, 0))
        self.kind.blockSignals(False)
        if self.form is not None:
            self.form.setParent(None)
            self.form.deleteLater()
            self.form = None
        self._rebuild()
        self._values = dict(spec.get("params", {}))

    def wavelengths(self) -> list[tuple[str, float]]:
        v = self.form.values() if self.form else {}
        return [(k, v[k]) for k in ("w1", "w2") if v.get(k) is not None]


# --------------------------------------------------------------------------- #
# Spectrum plot
# --------------------------------------------------------------------------- #
class SpectrumPlot(pg.PlotWidget):
    """Overlay plot with crosshair readout and hover identification."""

    hovered = Signal(str)
    line_moved = Signal(str, float)

    def __init__(self, parent: QWidget | None = None, y_label: str = "Absorbance"):
        super().__init__(parent)
        self.setLabel("bottom", "Wavelength (nm)")
        self.setLabel("left", y_label)
        self.showGrid(x=True, y=True, alpha=GRID_ALPHA)
        self.legend = self.addLegend(offset=(-10, 10), labelTextColor=TEXT_SECONDARY)
        self.curves: list[tuple[pg.PlotDataItem, Spectrum]] = []
        self.vlines: dict[str, pg.InfiniteLine] = {}
        pen = pg.mkPen("#8a8984", width=1, style=Qt.DashLine)
        self._vx = pg.InfiniteLine(angle=90, pen=pen)
        self._hy = pg.InfiniteLine(angle=0, pen=pen)
        for ln in (self._vx, self._hy):
            ln.setZValue(-5)
            self.addItem(ln, ignoreBounds=True)
        self._proxy = pg.SignalProxy(self.scene().sigMouseMoved, rateLimit=30,
                                     slot=self._mouse)

    def plot_spectra(self, spectra: list[Spectrum], keep_range: bool = False) -> None:
        vr = self.viewRange() if keep_range else None
        for item, _ in self.curves:
            self.removeItem(item)
        self.legend.clear()
        self.curves = []
        cols = colours_for(len(spectra))
        show_legend = len(spectra) <= 12
        for s, c in zip(spectra, cols):
            item = self.plot(s.wavelengths, s.values, pen=pg.mkPen(c, width=2),
                             name=s.name if show_legend else None)
            item.curve.setClickable(True, width=8)
            self.curves.append((item, s))
        self.legend.setVisible(show_legend and bool(spectra))
        if vr:
            self.setRange(xRange=vr[0], yRange=vr[1], padding=0)
        else:
            self.enableAutoRange()

    def set_marker(self, key: str, x: float | None, movable: bool = True,
                   colour: str = "#52514e") -> None:
        if x is None:
            if key in self.vlines:
                self.removeItem(self.vlines.pop(key))
            return
        ln = self.vlines.get(key)
        if ln is None:
            ln = pg.InfiniteLine(pos=x, angle=90, movable=movable,
                                 pen=pg.mkPen(colour, width=2),
                                 label=key, labelOpts={"position": 0.95, "color": colour})
            ln.sigPositionChangeFinished.connect(lambda l, k=key: self.line_moved.emit(k, l.value()))
            self.addItem(ln)
            self.vlines[key] = ln
        else:
            ln.blockSignals(True)
            ln.setValue(x)
            ln.blockSignals(False)

    def clear_markers(self) -> None:
        for ln in self.vlines.values():
            self.removeItem(ln)
        self.vlines = {}

    def _mouse(self, evt) -> None:
        pos = evt[0]
        if not self.sceneBoundingRect().contains(pos):
            return
        p = self.plotItem.vb.mapSceneToView(pos)
        self._vx.setPos(p.x())
        self._hy.setPos(p.y())
        # identify the nearest curve at this wavelength
        best, name = math.inf, ""
        for _, s in self.curves:
            if s.wavelengths[0] <= p.x() <= s.wavelengths[-1]:
                d = abs(np.interp(p.x(), s.wavelengths, s.values) - p.y())
                if d < best:
                    best, name, val = d, s.name, np.interp(p.x(), s.wavelengths, s.values)
        yr = self.viewRange()[1]
        near = name and best < 0.03 * (yr[1] - yr[0])
        text = f"λ = {p.x():.2f} nm   y = {p.y():.5f}"
        if near:
            text += f"   |   {name}: {val:.5f}"
        self.hovered.emit(text)


# --------------------------------------------------------------------------- #
# Table with copy/paste (Excel-compatible)
# --------------------------------------------------------------------------- #
class PasteTable(QTableWidget):
    def __init__(self, rows: int = 0, cols: int = 0, parent: QWidget | None = None):
        super().__init__(rows, cols, parent)
        self.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.horizontalHeader().setSectionResizeMode(QHeaderView.Interactive)
        self.horizontalHeader().setStretchLastSection(True)

    def keyPressEvent(self, e) -> None:
        if e.matches(QKeySequence.Copy):
            self.copy_selection()
        elif e.matches(QKeySequence.Paste):
            self.paste()
        elif e.key() in (Qt.Key_Delete, Qt.Key_Backspace):
            for it in self.selectedItems():
                if it.flags() & Qt.ItemIsEditable:
                    it.setText("")
        else:
            super().keyPressEvent(e)

    def copy_selection(self) -> None:
        rngs = self.selectedRanges()
        if not rngs:
            return
        r = rngs[0]
        lines = []
        for i in range(r.topRow(), r.bottomRow() + 1):
            lines.append("\t".join((self.item(i, j).text() if self.item(i, j) else "")
                                   for j in range(r.leftColumn(), r.rightColumn() + 1)))
        QGuiApplication.clipboard().setText("\n".join(lines))

    def copy_all(self) -> None:
        head = "\t".join(self.horizontalHeaderItem(j).text() if self.horizontalHeaderItem(j) else ""
                         for j in range(self.columnCount()))
        rows = ["\t".join((self.item(i, j).text() if self.item(i, j) else "")
                          for j in range(self.columnCount())) for i in range(self.rowCount())]
        QGuiApplication.clipboard().setText("\n".join([head] + rows))

    def paste(self) -> None:
        text = QGuiApplication.clipboard().text()
        if not text:
            return
        r0, c0 = max(self.currentRow(), 0), max(self.currentColumn(), 0)
        for i, line in enumerate(text.rstrip("\n").split("\n")):
            if r0 + i >= self.rowCount():
                self.setRowCount(r0 + i + 1)
            for j, cell in enumerate(line.split("\t")):
                if c0 + j >= self.columnCount():
                    continue
                it = self.item(r0 + i, c0 + j)
                if it is None:
                    it = QTableWidgetItem()
                    self.setItem(r0 + i, c0 + j, it)
                if it.flags() & Qt.ItemIsEditable:
                    it.setText(cell.strip())

    def column_floats(self, col: int) -> list[float]:
        out = []
        for i in range(self.rowCount()):
            it = self.item(i, col)
            t = it.text().strip().replace(",", ".") if it else ""
            if t:
                out.append(float(t))
        return out


def fill_table(table: QTableWidget, headers: list[str], rows: list[list[Any]],
               editable: bool = False) -> None:
    table.clear()
    table.setColumnCount(len(headers))
    table.setRowCount(len(rows))
    table.setHorizontalHeaderLabels(headers)
    for i, row in enumerate(rows):
        for j, v in enumerate(row):
            it = QTableWidgetItem(fmt(v) if not isinstance(v, str) else v)
            if isinstance(v, (int, float)) and not isinstance(v, bool):
                it.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)
            if not editable:
                it.setFlags(it.flags() & ~Qt.ItemIsEditable)
            table.setItem(i, j, it)
    table.resizeColumnsToContents()
    table.horizontalHeader().setStretchLastSection(True)


class SpectrumChecklist(QListWidget):
    """Checkable list of project spectra (id stored in UserRole)."""

    def __init__(self, parent: QWidget | None = None):
        super().__init__(parent)

    def populate(self, records, checked: set[int] | None = None) -> None:
        self.clear()
        for r in records:
            conc = ", ".join(f"{k}={v:g}" for k, v in r.concentrations.items())
            label = f"#{r.id}  {r.name}" + (f"  [{r.role}]" if r.role else "") + (
                f"  ({conc})" if conc else "")
            it = QListWidgetItem(label)
            it.setData(Qt.UserRole, r.id)
            it.setFlags(it.flags() | Qt.ItemIsUserCheckable)
            it.setCheckState(Qt.Checked if checked and r.id in checked else Qt.Unchecked)
            self.addItem(it)

    def checked_ids(self) -> list[int]:
        return [self.item(i).data(Qt.UserRole) for i in range(self.count())
                if self.item(i).checkState() == Qt.Checked]

    def set_all(self, state: bool) -> None:
        for i in range(self.count()):
            self.item(i).setCheckState(Qt.Checked if state else Qt.Unchecked)
