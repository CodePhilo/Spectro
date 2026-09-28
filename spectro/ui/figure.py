"""Publication-quality figure export for any plot in the application.

The data shown in a pyqtgraph plot (curves, markers, reference lines, axis
labels and view range) is collected and redrawn with matplotlib at journal
settings: size in mm, 300–1200 dpi raster (PNG/TIFF) or vector (SVG/PDF/EPS),
font family and size, line width, legend placement and a black-and-white mode.
"""

from __future__ import annotations

import io
from typing import Any

import numpy as np
import pyqtgraph as pg
from PySide6.QtCore import Qt
from PySide6.QtGui import QAction, QPixmap
from PySide6.QtWidgets import (QCheckBox, QComboBox, QDialog, QDialogButtonBox,
                               QDoubleSpinBox, QFileDialog, QFormLayout, QHBoxLayout, QLabel,
                               QLineEdit, QMessageBox, QSpinBox, QVBoxLayout, QWidget)

FORMATS = {"PNG": "png", "TIFF": "tif", "SVG": "svg", "PDF": "pdf", "EPS": "eps"}
PRESETS = {"Single column (85 mm)": (85, 65), "1.5 column (120 mm)": (120, 85),
           "Double column (170 mm)": (170, 105), "Custom": None}
LINESTYLES = ["-", "--", "-.", ":", (0, (5, 1)), (0, (3, 1, 1, 1, 1, 1))]
MARKERS = ["o", "s", "^", "D", "v", "P", "X", "*"]


def _colour(obj: Any) -> str | None:
    if obj is None:
        return None
    try:
        return pg.mkColor(obj).name()
    except Exception:
        try:
            return pg.mkPen(obj).color().name()
        except Exception:
            return None


def collect_series(plot: pg.PlotWidget) -> dict:
    """Everything needed to redraw ``plot``: series, reference lines, labels, range."""
    pi = plot.getPlotItem()
    series = []
    for item in pi.listDataItems():
        if isinstance(item, pg.PlotDataItem):
            x, y = item.getData()
            if x is None or len(x) == 0:
                continue
            pen = item.opts.get("pen")
            has_line = pen is not None and pg.mkPen(pen).style() != Qt.NoPen
            symbol = item.opts.get("symbol")
            series.append({"x": np.asarray(x), "y": np.asarray(y), "label": item.name(),
                           "color": _colour(pg.mkPen(pen).color()) if has_line
                           else _colour(item.opts.get("symbolBrush")),
                           "line": has_line, "marker": symbol is not None})
        elif isinstance(item, pg.ScatterPlotItem):
            x, y = item.getData()
            if x is None or len(x) == 0:
                continue
            series.append({"x": np.asarray(x), "y": np.asarray(y), "label": item.name(),
                           "color": _colour(item.opts["brush"].color()
                                            if hasattr(item.opts["brush"], "color")
                                            else item.opts["brush"]),
                           "line": False, "marker": True})
    refs = []
    for item in pi.items:
        if isinstance(item, pg.InfiniteLine) and item.zValue() >= 0 and item.isVisible():
            label = item.label.format if getattr(item, "label", None) is not None else ""
            refs.append({"angle": item.angle, "pos": float(item.value()), "label": label,
                         "color": _colour(item.pen.color()),
                         "style": {Qt.DashLine: "--", Qt.DotLine: ":"}.get(item.pen.style(), "-")})
    (x0, x1), (y0, y1) = pi.viewRange()
    xs = [s["x"] for s in series if len(s["x"])]
    if xs and all(s["line"] for s in series):
        # spectra: axis exactly spanning the data (no view padding), unless zoomed in
        dmin = float(min(np.nanmin(x) for x in xs))
        dmax = float(max(np.nanmax(x) for x in xs))
        x0, x1 = max(x0, dmin), min(x1, dmax)
    return {"series": series, "refs": refs,
            "xlabel": pi.getAxis("bottom").labelText or "",
            "ylabel": pi.getAxis("left").labelText or "",
            "xlim": (x0, x1), "ylim": (y0, y1)}


def render(data: dict, opts: dict, dpi: float):
    """Draw with matplotlib; returns the Figure."""
    from matplotlib.figure import Figure

    w, h = opts["width_mm"] / 25.4, opts["height_mm"] / 25.4
    import logging
    logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)
    rc = {"font.family": [opts["font"], "Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
          "font.size": opts["font_size"],
          "axes.linewidth": 0.8, "xtick.direction": "in", "ytick.direction": "in",
          "xtick.major.width": 0.8, "ytick.major.width": 0.8,
          "svg.fonttype": "none", "pdf.fonttype": 42, "ps.fonttype": 42}
    import matplotlib

    with matplotlib.rc_context(rc):
        fig = Figure(figsize=(w, h), dpi=dpi, layout="constrained")
        ax = fig.add_subplot()
        bw = opts["black_white"]
        for i, s in enumerate(data["series"]):
            col = "black" if bw else (s["color"] or "black")
            label = s["label"] if (s["label"] and opts["legend"] != "None") else None
            if s["line"]:
                ax.plot(s["x"], s["y"], color=col, lw=opts["line_width"], label=label,
                        ls=LINESTYLES[i % len(LINESTYLES)] if bw else "-",
                        marker=(MARKERS[i % len(MARKERS)] if bw else "o") if s["marker"] else None,
                        ms=opts["marker_size"])
            else:
                ax.plot(s["x"], s["y"], ls="none", marker=MARKERS[i % len(MARKERS)] if bw else "o",
                        mfc="white" if bw else col, mec="black" if bw else "white",
                        mew=0.8, ms=opts["marker_size"], label=label)
        if opts["show_refs"]:
            for r in data["refs"]:
                col = "0.3" if bw else (r["color"] or "0.3")
                if r["angle"] == 90:
                    ax.axvline(r["pos"], color=col, lw=0.8, ls=r["style"])
                    if r["label"]:
                        ax.annotate(r["label"], (r["pos"], 1), xycoords=("data", "axes fraction"),
                                    xytext=(2, -2), textcoords="offset points", va="top",
                                    fontsize=opts["font_size"] * 0.8, color=col)
                else:
                    ax.axhline(r["pos"], color=col, lw=0.8, ls=r["style"])
        ax.set_xlabel(opts["xlabel"])
        ax.set_ylabel(opts["ylabel"])
        if opts["title"]:
            ax.set_title(opts["title"])
        ax.set_xlim(*opts["xlim"])
        ax.set_ylim(*opts["ylim"])
        if opts["grid"]:
            ax.grid(True, lw=0.4, color="0.85")
        if not opts["box"]:
            ax.spines[["top", "right"]].set_visible(False)
        else:
            ax.tick_params(top=True, right=True)
        if opts["legend"] != "None" and any(s["label"] for s in data["series"]):
            if opts["legend"] == "Outside right":
                ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1), frameon=False,
                          fontsize=opts["font_size"] * 0.85)
            else:
                ax.legend(loc="best", frameon=False, fontsize=opts["font_size"] * 0.85)
    return fig


class FigureDialog(QDialog):
    def __init__(self, plot: pg.PlotWidget, parent: QWidget | None = None,
                 title: str = ""):
        super().__init__(parent)
        self.setWindowTitle("Export publication figure")
        self.data = collect_series(plot)
        self.resize(1050, 620)
        lay = QHBoxLayout(self)
        f = QFormLayout()
        self.preset = QComboBox()
        self.preset.addItems(list(PRESETS))
        self.w, self.h = QDoubleSpinBox(), QDoubleSpinBox()
        for sb, v in ((self.w, 85), (self.h, 65)):
            sb.setRange(20, 500)
            sb.setSuffix(" mm")
            sb.setValue(v)
        self.preset.currentTextChanged.connect(self._preset)
        self.fmt = QComboBox()
        self.fmt.addItems(list(FORMATS))
        self.fmt.setCurrentText("TIFF")
        self.dpi = QSpinBox()
        self.dpi.setRange(72, 2400)
        self.dpi.setValue(600)
        self.font = QComboBox()
        self.font.addItems(["Arial", "Helvetica", "Times New Roman", "Calibri", "DejaVu Sans"])
        self.font.setEditable(True)
        self.fsize = QDoubleSpinBox()
        self.fsize.setRange(5, 20)
        self.fsize.setValue(8)
        self.fsize.setSuffix(" pt")
        self.lw = QDoubleSpinBox()
        self.lw.setRange(0.2, 5)
        self.lw.setSingleStep(0.25)
        self.lw.setValue(1.0)
        self.ms = QDoubleSpinBox()
        self.ms.setRange(1, 15)
        self.ms.setValue(4)
        self.legend = QComboBox()
        self.legend.addItems(["Best", "Outside right", "None"])
        self.bw = QCheckBox("Black and white (line styles / markers)")
        self.box = QCheckBox("Closed box, ticks on all sides")
        self.box.setChecked(True)
        self.grid = QCheckBox("Grid")
        self.refs = QCheckBox("Show marker / reference lines")
        self.refs.setChecked(True)
        self.title = QLineEdit(title)
        self.xlabel = QLineEdit(self.data["xlabel"])
        self.ylabel = QLineEdit(self.data["ylabel"])
        (x0, x1), (y0, y1) = self.data["xlim"], self.data["ylim"]
        self.lims = []
        for v in (x0, x1, y0, y1):
            sb = QDoubleSpinBox()
            sb.setRange(-1e9, 1e9)
            sb.setDecimals(4)
            sb.setValue(v)
            self.lims.append(sb)
        xr, yr = QHBoxLayout(), QHBoxLayout()
        for w_ in self.lims[:2]:
            xr.addWidget(w_)
        for w_ in self.lims[2:]:
            yr.addWidget(w_)
        for label, w_ in (("Size", self.preset), ("Width", self.w), ("Height", self.h),
                          ("Format", self.fmt), ("Resolution (dpi)", self.dpi),
                          ("Font", self.font), ("Font size", self.fsize),
                          ("Line width", self.lw), ("Marker size", self.ms),
                          ("Legend", self.legend), ("Title", self.title),
                          ("X label", self.xlabel), ("Y label", self.ylabel)):
            f.addRow(label, w_)
        f.addRow("X range", xr)
        f.addRow("Y range", yr)
        for w_ in (self.bw, self.box, self.grid, self.refs):
            f.addRow(w_)
        left = QWidget()
        lv = QVBoxLayout(left)
        lv.addLayout(f)
        bb = QDialogButtonBox(QDialogButtonBox.Save | QDialogButtonBox.Close)
        bb.accepted.connect(self._save)
        bb.rejected.connect(self.reject)
        lv.addStretch(1)
        lv.addWidget(bb)
        lay.addWidget(left)
        self.preview = QLabel()
        self.preview.setAlignment(Qt.AlignCenter)
        self.preview.setMinimumWidth(560)
        self.preview.setStyleSheet("background: #f0efec;")
        lay.addWidget(self.preview, 1)
        for w_ in (self.w, self.h, self.fsize, self.lw, self.ms, *self.lims):
            w_.valueChanged.connect(self._update)
        for w_ in (self.font, self.legend):
            w_.currentTextChanged.connect(self._update)
        for w_ in (self.bw, self.box, self.grid, self.refs):
            w_.toggled.connect(self._update)
        for w_ in (self.title, self.xlabel, self.ylabel):
            w_.editingFinished.connect(self._update)
        self._update()

    def _preset(self, name):
        size = PRESETS.get(name)
        if size:
            self.w.setValue(size[0])
            self.h.setValue(size[1])

    def options(self) -> dict:
        return {"width_mm": self.w.value(), "height_mm": self.h.value(),
                "font": self.font.currentText(), "font_size": self.fsize.value(),
                "line_width": self.lw.value(), "marker_size": self.ms.value(),
                "legend": self.legend.currentText(), "black_white": self.bw.isChecked(),
                "box": self.box.isChecked(), "grid": self.grid.isChecked(),
                "show_refs": self.refs.isChecked(), "title": self.title.text(),
                "xlabel": self.xlabel.text(), "ylabel": self.ylabel.text(),
                "xlim": (self.lims[0].value(), self.lims[1].value()),
                "ylim": (self.lims[2].value(), self.lims[3].value())}

    def _update(self):
        try:
            fig = render(self.data, self.options(), dpi=110)
            buf = io.BytesIO()
            fig.savefig(buf, format="png", dpi=110)
            pm = QPixmap()
            pm.loadFromData(buf.getvalue())
            self.preview.setPixmap(pm.scaled(self.preview.size() * 0.95, Qt.KeepAspectRatio,
                                             Qt.SmoothTransformation)
                                   if pm.width() > self.preview.width() else pm)
        except Exception as exc:  # show problems in the preview area
            self.preview.setText(f"Preview failed: {exc}")

    def _save(self):
        fmt = self.fmt.currentText()
        ext = FORMATS[fmt]
        path, _ = QFileDialog.getSaveFileName(self, "Save figure", f"figure.{ext}",
                                              f"{fmt} (*.{ext})")
        if not path:
            return
        if not path.lower().endswith("." + ext):
            path += "." + ext
        opts = self.options()
        try:
            fig = render(self.data, opts, dpi=self.dpi.value())
            kw = {"pil_kwargs": {"compression": "tiff_lzw"}} if ext == "tif" else {}
            fig.savefig(path, dpi=self.dpi.value(), **kw)
        except Exception as exc:
            QMessageBox.warning(self, "Export figure", str(exc))
            return
        from spectro.ui.widgets import log_export
        log_export(f"Exported figure ({fmt}, {self.dpi.value()} dpi)",
                   {"file": path, "format": fmt, "dpi": self.dpi.value(),
                    "series": [s["label"] for s in self.data["series"]],
                    "size_mm": [opts["width_mm"], opts["height_mm"]]})
        QMessageBox.information(self, "Export figure", f"Saved {path}")


def install_figure_export(plot: pg.PlotWidget) -> None:
    """Add 'Export publication figure…' to the plot's right-click menu."""
    menu = plot.getPlotItem().vb.menu
    act = QAction("Export publication figure…", menu)
    act.triggered.connect(lambda: FigureDialog(plot, plot.window()).exec())
    menu.addSeparator()
    menu.addAction(act)
