"""Layout-driven (instrument-independent) parsing of tabular spectral data.

Text files and spreadsheets from any spectrophotometer are reduced to a grid of
cells and then interpreted by *layout*, not by vendor:

``columns``   one wavelength column followed by one column per spectrum
``xy_pairs``  repeating (wavelength, value) column pairs, one pair per spectrum
``rows``      one wavelength row, then one row per spectrum (optionally with a
              sample-name column on the left)

Header/metadata lines above the numeric block are harvested as metadata and
used for spectrum names.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass, field
from typing import Any, Iterable, Sequence

import numpy as np

from spectro.core.spectrum import Spectrum

LAYOUTS = ("auto", "columns", "xy_pairs", "rows")

# Header words that describe an axis rather than name a sample.
_GENERIC = {
    "", "wavelength", "wavelength (nm)", "wavelength nm", "wavelength nm.", "nm",
    "wl", "lambda", "x", "y", "abs", "abs.", "absorbance", "absorbance (au)", "a",
    "au", "%t", "t", "transmittance", "data", "value", "values", "intensity",
    "wavelength[nm]", "abs[au]", "absorbance[au]", "(nm)", "nm.",
}

_NUM_RE = re.compile(r"^[+-]?(\d+([.]\d*)?|[.]\d+)([eE][+-]?\d+)?$")


def to_float(cell: Any, decimal_comma: bool = False) -> float:
    """Convert a cell to float; return NaN for anything non-numeric."""
    if cell is None:
        return math.nan
    if isinstance(cell, bool):
        return math.nan
    if isinstance(cell, (int, float, np.floating, np.integer)):
        return float(cell)
    s = str(cell).strip().strip('"').strip("'").strip()
    if not s:
        return math.nan
    s = s.replace("−", "-").replace(" ", "").replace(" ", "")
    if decimal_comma:
        s = s.replace(",", ".")
    if not _NUM_RE.match(s):
        return math.nan
    return float(s)


def _cell_text(cell: Any) -> str:
    if cell is None:
        return ""
    return str(cell).strip().strip('"').strip()


@dataclass
class ParseOptions:
    layout: str = "auto"
    decimal_comma: bool | None = None  # None = auto-detect (text only)
    name_prefix: str = ""


@dataclass
class ParseResult:
    spectra: list[Spectrum]
    layout: str
    metadata: dict[str, str] = field(default_factory=dict)
    warnings: list[str] = field(default_factory=list)


# --------------------------------------------------------------------------- #
# numeric helpers
# --------------------------------------------------------------------------- #
def _longest_finite_run(v: np.ndarray) -> tuple[int, int]:
    best = (0, 0)
    start = None
    for i, ok in enumerate(np.append(np.isfinite(v), False)):
        if ok and start is None:
            start = i
        elif not ok and start is not None:
            if i - start > best[1] - best[0]:
                best = (start, i)
            start = None
    return best


def _is_axis_like(v: np.ndarray) -> bool:
    """True if ``v`` looks like a wavelength axis (monotonic, roughly regular,
    values in a plausible wavelength / wavenumber range)."""
    if v.size < 3 or not np.all(np.isfinite(v)):
        return False
    d = np.diff(v)
    if not (np.all(d > 0) or np.all(d < 0)):
        return False
    ad = np.abs(d)
    regular = np.median(ad) > 0 and (np.max(ad) / np.median(ad) < 5.0)
    return bool(regular and np.min(np.abs(v)) >= 50.0)


def _meaningful(text: str) -> bool:
    return text.lower() not in _GENERIC and not re.fullmatch(r"[\s\-_=:#]*", text)


# --------------------------------------------------------------------------- #
# main entry
# --------------------------------------------------------------------------- #
def parse_grid(rows: Sequence[Sequence[Any]], options: ParseOptions | None = None,
               default_name: str = "Spectrum") -> ParseResult:
    """Interpret a 2-D grid of cells as one or more spectra."""
    options = options or ParseOptions()
    if options.layout not in LAYOUTS:
        raise ValueError(f"unknown layout {options.layout!r}")
    width = max((len(r) for r in rows), default=0)
    if width == 0:
        raise ValueError("the file contains no data")
    text = [[_cell_text(c) for c in r] + [""] * (width - len(r)) for r in rows]
    dc = bool(options.decimal_comma)
    num = np.array([[to_float(c, dc) for c in r] + [math.nan] * (width - len(r))
                    for r in rows], dtype=float)

    candidates: list[ParseResult] = []
    errors: list[str] = []
    wanted = [options.layout] if options.layout != "auto" else ["xy_pairs", "columns", "rows"]
    for layout in wanted:
        try:
            res = {"columns": _parse_columns, "xy_pairs": _parse_xy_pairs,
                   "rows": _parse_rows}[layout](num, text, default_name)
        except ValueError as exc:
            errors.append(f"{layout}: {exc}")
            continue
        if res.spectra:
            candidates.append(res)
    if not candidates:
        raise ValueError("could not find spectral data. " + "; ".join(errors))

    # Auto: a successful xy_pairs parse is a strong signal (every other column
    # is a wavelength axis) so it wins; otherwise keep the one with most points.
    pairs = [r for r in candidates if r.layout == "xy_pairs"]
    best = pairs[0] if pairs else max(
        candidates, key=lambda r: sum(s.wavelengths.size for s in r.spectra))
    best.metadata.update(_harvest_metadata(text, num))
    if options.name_prefix:
        for s in best.spectra:
            s.name = f"{options.name_prefix}{s.name}"
    return best


def _name_above(text: list[list[str]], row: int, cols: Iterable[int]) -> str:
    for r in range(row - 1, max(-1, row - 6), -1):
        first = next((t for t in text[r] if t), "")
        if first.startswith(("#", "//", ";")):
            continue
        for c in cols:
            if c < len(text[r]) and _meaningful(text[r][c]):
                return text[r][c]
    return ""


def _parse_columns(num, text, default_name) -> ParseResult:
    ncol = num.shape[1]
    for xc in range(ncol):
        a, b = _longest_finite_run(num[:, xc])
        if b - a < 3 or not _is_axis_like(num[a:b, xc]):
            continue
        x = num[a:b, xc]
        spectra: list[Spectrum] = []
        for c in range(ncol):
            if c == xc:
                continue
            y = num[a:b, c]
            ok = np.isfinite(y)
            if ok.sum() < 3 or ok.sum() < 0.8 * (b - a):
                continue
            name = _name_above(text, a, [c]) or f"{default_name} #{len(spectra) + 1}"
            spectra.append(Spectrum(x[ok], y[ok], name=name))
        if spectra:
            if len(spectra) == 1:
                c = next(i for i in range(ncol) if i != xc and np.isfinite(num[a:b, i]).sum() >= 3)
                spectra[0].name = _name_above(text, a, [c, xc]) or default_name
            return ParseResult(spectra, "columns")
    raise ValueError("no wavelength column found")


def _parse_xy_pairs(num, text, default_name) -> ParseResult:
    ncol = num.shape[1]
    if ncol < 4:
        raise ValueError("fewer than two column pairs")
    spectra: list[Spectrum] = []
    first_start = None
    for c in range(0, ncol - 1, 2):
        a, b = _longest_finite_run(num[:, c])
        if b - a < 3:
            # allow trailing empty pairs (e.g. Cary's trailing comma)
            if np.isfinite(num[:, c]).sum() == 0 and np.isfinite(num[:, c + 1]).sum() == 0:
                continue
            raise ValueError(f"column {c + 1} is not a wavelength column")
        x = num[a:b, c]
        if not _is_axis_like(x):
            raise ValueError(f"column {c + 1} is not a wavelength column")
        y = num[a:b, c + 1]
        ok = np.isfinite(y)
        if ok.sum() < 3:
            raise ValueError(f"column {c + 2} has no data")
        first_start = a if first_start is None else first_start
        name = _name_above(text, a, [c, c + 1]) or f"{default_name} #{len(spectra) + 1}"
        spectra.append(Spectrum(x[ok], y[ok], name=name))
    if len(spectra) < 2:
        raise ValueError("fewer than two column pairs")
    return ParseResult(spectra, "xy_pairs")


def _parse_rows(num, text, default_name) -> ParseResult:
    nrow = num.shape[0]
    for r in range(nrow):
        a, b = _longest_finite_run(num[r, :])
        if b - a < 5 or not _is_axis_like(num[r, a:b]):
            continue
        x = num[r, a:b]
        spectra: list[Spectrum] = []
        for rr in range(r + 1, nrow):
            y = num[rr, a:b]
            ok = np.isfinite(y)
            if ok.sum() < 0.8 * (b - a):
                continue
            label = next((text[rr][c] for c in range(a) if _meaningful(text[rr][c])), "")
            if not label and a > 0 and np.isfinite(num[rr, 0]):
                label = f"{default_name} {text[rr][0]}"
            spectra.append(Spectrum(x[ok], y[ok],
                                    name=label or f"{default_name} #{len(spectra) + 1}"))
        if spectra:
            return ParseResult(spectra, "rows")
    raise ValueError("no wavelength row found")


def _harvest_metadata(text: list[list[str]], num: np.ndarray) -> dict[str, str]:
    """Collect ``key: value`` / ``key<TAB>value`` header lines."""
    meta: dict[str, str] = {}
    for r, row in enumerate(text[:80]):
        if np.isfinite(num[r]).sum() >= 2:
            break
        cells = [c for c in row if c]
        if not cells:
            continue
        if len(cells) == 1 and ":" in cells[0]:
            k, v = cells[0].split(":", 1)
        elif len(cells) == 1 and "=" in cells[0]:
            k, v = cells[0].split("=", 1)
        elif len(cells) == 2:
            k, v = cells
        else:
            continue
        k, v = k.strip().strip("#").strip(), v.strip()
        if k and v and len(meta) < 60:
            meta[k] = v
    return meta
