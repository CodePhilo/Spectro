"""Delimited text files (CSV, TSV, TXT, ASC, PRN, DAT...) from any instrument."""

from __future__ import annotations

import csv
import io
import re
from pathlib import Path

from spectro.core.io.table import ParseOptions, ParseResult, parse_grid, to_float

_DELIMS = ["\t", ";", "|", ",", "ws"]
_DECIMAL_COMMA_RE = re.compile(r"(?<![\d,])[+-]?\d+,\d+(?![\d,])")


def read_text_file(path: str | Path) -> str:
    raw = Path(path).read_bytes()
    if raw.startswith((b"\xff\xfe", b"\xfe\xff")):
        return raw.decode("utf-16")
    if raw[:200].count(b"\x00") > 20:  # UTF-16 without BOM
        try:
            return raw.decode("utf-16-le")
        except UnicodeDecodeError:
            pass
    for enc in ("utf-8-sig", "cp1252", "latin-1"):
        try:
            return raw.decode(enc)
        except UnicodeDecodeError:
            continue
    return raw.decode("latin-1", errors="replace")


def split_lines(content: str, delimiter: str) -> list[list[str]]:
    lines = content.splitlines()
    if delimiter == "ws":
        return [line.split() for line in lines]
    return list(csv.reader(io.StringIO("\n".join(lines)), delimiter=delimiter))


def detect_format(content: str) -> tuple[str, bool]:
    """Return ``(delimiter, decimal_comma)`` that best explains the numbers."""
    sample = "\n".join(content.splitlines()[:400])
    best: tuple[int, str, bool] = (-1, "ws", False)
    for delim in _DELIMS:
        rows = split_lines(sample, delim)
        for dc in (False, True):
            if dc and delim == ",":
                continue
            score = 0
            for r in rows:
                n = sum(1 for c in r if to_float(c, dc) == to_float(c, dc))  # not NaN
                if n >= 2:
                    score += n
            if score > best[0]:
                best = (score, delim, dc)
    _, delim, dc = best
    # Decimal comma only if comma-decimals actually occur.
    if dc and not _DECIMAL_COMMA_RE.search(sample):
        dc = False
    return delim, dc


def parse_text(content: str, options: ParseOptions | None = None,
               default_name: str = "Spectrum",
               delimiter: str | None = None) -> ParseResult:
    options = options or ParseOptions()
    auto_delim, auto_dc = detect_format(content)
    delim = delimiter or auto_delim
    dc = auto_dc if options.decimal_comma is None else options.decimal_comma
    rows = split_lines(content, delim)
    opts = ParseOptions(layout=options.layout, decimal_comma=dc,
                        name_prefix=options.name_prefix)
    result = parse_grid(rows, opts, default_name=default_name)
    result.metadata.setdefault("_delimiter", {"\t": "tab", "ws": "whitespace"}.get(delim, delim))
    result.metadata.setdefault("_decimal", "comma" if dc else "point")
    return result


def load_text(path: str | Path, options: ParseOptions | None = None) -> ParseResult:
    path = Path(path)
    return parse_text(read_text_file(path), options, default_name=path.stem)
