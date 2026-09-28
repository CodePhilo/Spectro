"""File import / export, independent of instrument make or model."""

from __future__ import annotations

import hashlib
from pathlib import Path

from spectro.core.io.excel import load_excel
from spectro.core.io.ole import load_ole
from spectro.core.io.jcamp import load_jcamp, write_jcamp
from spectro.core.io.spc import is_spc, load_spc, write_spc
from spectro.core.io.table import LAYOUTS, ParseOptions, ParseResult
from spectro.core.io.text import load_text
from spectro.core.spectrum import Spectrum, to_matrix

EXCEL_EXT = {".xlsx", ".xlsm", ".xls"}
JCAMP_EXT = {".jdx", ".dx", ".jcm", ".jcamp"}
SPC_EXT = {".spc"}

IMPORT_FILTER = (
    "All supported (*.csv *.tsv *.txt *.asc *.prn *.dat *.xlsx *.xlsm *.xls "
    "*.jdx *.dx *.jcm *.spc *.*);;"
    "Text / CSV (*.csv *.tsv *.txt *.asc *.prn *.dat);;"
    "Excel (*.xlsx *.xlsm *.xls);;JCAMP-DX (*.jdx *.dx *.jcm);;SPC (*.spc);;"
    "All files (*.*)"
)

__all__ = ["load_file", "export_csv", "export_excel", "write_jcamp", "write_spc",
           "ParseOptions", "ParseResult", "LAYOUTS", "file_sha256"]


def file_sha256(path: str | Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_file(path: str | Path, options: ParseOptions | None = None) -> ParseResult:
    """Load spectra from any supported file, detecting its type by content."""
    path = Path(path)
    ext = path.suffix.lower()
    with open(path, "rb") as fh:
        head = fh.read(512)
    if ext in EXCEL_EXT or head.startswith(b"PK\x03\x04") or head.startswith(b"\xd0\xcf\x11\xe0"):
        if head.startswith(b"\xd0\xcf\x11\xe0") and ext not in EXCEL_EXT:
            result = load_ole(path)  # e.g. Shimadzu UVProbe .spc
        else:
            result = load_excel(path, options)
    elif ext in SPC_EXT or is_spc(head):
        result = load_spc(path)
    elif ext in JCAMP_EXT or head.lstrip().startswith(b"##"):
        result = load_jcamp(path)
    else:
        result = load_text(path, options)
    for s in result.spectra:
        s.metadata.setdefault("source_file", str(path))
        s.metadata.setdefault("import_layout", result.layout)
    return result


def export_csv(spectra: list[Spectrum], path: str | Path, delimiter: str = ",") -> None:
    grid, mat = to_matrix(spectra)
    with open(path, "w", encoding="utf-8", newline="") as fh:
        fh.write(delimiter.join(["Wavelength (nm)"] + [_q(s.name, delimiter) for s in spectra]) + "\n")
        for i, wl in enumerate(grid):
            fh.write(delimiter.join([f"{wl:.6g}"] + [f"{v:.8g}" for v in mat[:, i]]) + "\n")


def export_excel(spectra: list[Spectrum], path: str | Path) -> None:
    import openpyxl

    grid, mat = to_matrix(spectra)
    wb = openpyxl.Workbook()
    ws = wb.active
    ws.title = "Spectra"
    ws.append(["Wavelength (nm)"] + [s.name for s in spectra])
    for i, wl in enumerate(grid):
        ws.append([float(wl)] + [float(v) for v in mat[:, i]])
    wb.save(str(path))


def _q(text: str, delim: str) -> str:
    if delim in text or '"' in text:
        return '"' + text.replace('"', '""') + '"'
    return text
