"""Excel workbooks (.xlsx/.xlsm via openpyxl, legacy .xls via xlrd)."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from spectro.core.io.table import ParseOptions, ParseResult, parse_grid


def _sheets(path: Path) -> list[tuple[str, list[list[Any]]]]:
    if path.suffix.lower() == ".xls":
        try:
            import xlrd
        except ImportError as exc:  # pragma: no cover
            raise ValueError("reading .xls files requires the 'xlrd' package") from exc
        book = xlrd.open_workbook(str(path))
        return [(sh.name, [sh.row_values(i) for i in range(sh.nrows)])
                for sh in book.sheets()]
    import openpyxl

    book = openpyxl.load_workbook(str(path), read_only=True, data_only=True)
    try:
        return [(ws.title, [list(r) for r in ws.iter_rows(values_only=True)])
                for ws in book.worksheets]
    finally:
        book.close()


def load_excel(path: str | Path, options: ParseOptions | None = None,
               sheet: str | None = None) -> ParseResult:
    """Parse every sheet that contains spectra (or only ``sheet``)."""
    path = Path(path)
    options = options or ParseOptions()
    sheets = [(n, rows) for n, rows in _sheets(path) if sheet is None or n == sheet]
    if not sheets:
        raise ValueError(f"sheet {sheet!r} not found")
    combined: ParseResult | None = None
    problems: list[str] = []
    multi = len(sheets) > 1
    for name, rows in sheets:
        rows = [r for r in rows if r is not None]
        if not any(c not in (None, "") for r in rows for c in r):
            continue
        try:
            res = parse_grid(rows, options, default_name=name if multi else path.stem)
        except ValueError as exc:
            problems.append(f"sheet '{name}': {exc}")
            continue
        for s in res.spectra:
            s.metadata["sheet"] = name
            if multi and name not in s.name:
                s.name = f"{name} / {s.name}"
        if combined is None:
            combined = res
        else:
            combined.spectra.extend(res.spectra)
            for k, v in res.metadata.items():
                combined.metadata.setdefault(f"{name}: {k}", v)
    if combined is None:
        raise ValueError("no spectra found in workbook. " + "; ".join(problems))
    combined.warnings.extend(problems)
    return combined
