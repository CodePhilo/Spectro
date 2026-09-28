"""JCAMP-DX (.jdx/.dx/.jcm) reader and writer, including ASDF compression."""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np

from spectro.core.io.table import ParseResult
from spectro.core.io.text import read_text_file
from spectro.core.spectrum import Spectrum

_SQZ = {"@": 0, **{c: i + 1 for i, c in enumerate("ABCDEFGHI")},
        **{c: -(i + 1) for i, c in enumerate("abcdefghi")}}
_DIF = {"%": 0, **{c: i + 1 for i, c in enumerate("JKLMNOPQR")},
        **{c: -(i + 1) for i, c in enumerate("jklmnopqr")}}
_DUP = {**{c: i + 1 for i, c in enumerate("STUVWXYZ")}, "s": 9}


def decode_asdf_line(line: str) -> list[float]:
    """Decode one XYDATA line (AFFN, SQZ, DIF, DUP forms) into numbers.

    The first returned number is the abscissa.
    """
    tokens: list[tuple[str, str]] = []  # (kind, digits)
    i, n = 0, len(line)
    while i < n:
        ch = line[i]
        if ch in " \t,":
            i += 1
            continue
        if ch in "+-.0123456789" and not (ch == "." and False):
            j = i + 1
            while j < n and (line[j].isdigit() or line[j] in ".eE" or
                             (line[j] in "+-" and line[j - 1] in "eE")):
                j += 1
            tokens.append(("abs", line[i:j]))
            i = j
            continue
        kind = "sqz" if ch in _SQZ else "dif" if ch in _DIF else "dup" if ch in _DUP else None
        if kind is None:
            if ch == "?":  # missing value
                tokens.append(("abs", "nan"))
            i += 1
            continue
        j = i + 1
        while j < n and (line[j].isdigit() or line[j] == "."):
            j += 1
        lead = {"sqz": _SQZ, "dif": _DIF, "dup": _DUP}[kind][ch]
        digits = line[i + 1:j]
        if kind == "dup":
            val = f"{lead}{digits}"
        else:
            sign = "-" if lead < 0 else ""
            val = f"{sign}{abs(lead)}{digits}"
        tokens.append((kind, val))
        i = j

    values: list[float] = []
    last_kind = "abs"
    last_diff = 0.0
    for kind, val in tokens:
        if kind in ("abs", "sqz"):
            values.append(float(val))
            last_kind = "abs"
        elif kind == "dif":
            if not values:
                raise ValueError("DIF value without preceding value")
            last_diff = float(val)
            values.append(values[-1] + last_diff)
            last_kind = "dif"
        else:  # dup
            count = int(float(val))
            for _ in range(count - 1):
                values.append(values[-1] + last_diff if last_kind == "dif" else values[-1])
    return values


def _parse_blocks(content: str) -> list[dict]:
    blocks: list[dict] = []
    cur: dict | None = None
    key = None
    for raw in content.splitlines():
        line = raw.split("$$", 1)[0].rstrip()
        m = re.match(r"^##\s*([^=]+?)\s*=\s*(.*)$", line)
        if m:
            key = re.sub(r"[\s\-_/]", "", m.group(1)).upper()
            val = m.group(2).strip()
            if key == "TITLE":
                cur = {"_data": [], "_datakey": None}
                blocks.append(cur)
            if cur is None:
                cur = {"_data": [], "_datakey": None}
                blocks.append(cur)
            if key in ("XYDATA", "XYPOINTS", "PEAKTABLE", "DATATABLE"):
                cur["_datakey"] = key
                cur["_format"] = val
            elif key == "END":
                cur = None if cur and cur.get("_datakey") else cur
            else:
                cur[key] = val
            continue
        if cur is not None and cur.get("_datakey") and line.strip():
            cur["_data"].append(line)
    return [b for b in blocks if b.get("_datakey")]


def _f(block: dict, key: str, default: float | None = None) -> float | None:
    try:
        return float(block[key])
    except (KeyError, ValueError):
        return default


def _block_to_spectrum(block: dict, fallback_name: str) -> Spectrum:
    xf = _f(block, "XFACTOR", 1.0)
    yf = _f(block, "YFACTOR", 1.0)
    fmt = block.get("_format", "").replace(" ", "").upper()
    if "X++" in fmt:
        firstx, lastx = _f(block, "FIRSTX"), _f(block, "LASTX")
        npts = _f(block, "NPOINTS")
        delta = _f(block, "DELTAX")
        if delta is None and None not in (firstx, lastx, npts) and npts > 1:
            delta = (lastx - firstx) / (npts - 1)
        xs: list[float] = []
        ys: list[float] = []
        for line in block["_data"]:
            vals = decode_asdf_line(line)
            if len(vals) < 2:
                continue
            x0 = vals[0] * xf
            line_ys = vals[1:]
            if delta is None:
                raise ValueError("XYDATA without DELTAX or FIRSTX/LASTX/NPOINTS")
            line_xs = [x0 + k * delta for k in range(len(line_ys))]
            # Y-check value: first point repeats the last point of previous line
            if xs and abs(line_xs[0] - xs[-1]) < abs(delta) / 2:
                line_xs, line_ys = line_xs[1:], line_ys[1:]
            xs.extend(line_xs)
            ys.extend(line_ys)
        x = np.array(xs)
        if firstx is not None and len(xs) > 1 and npts and len(xs) == int(npts):
            x = np.linspace(firstx, lastx if lastx is not None else xs[-1], len(xs))
        y = np.array(ys) * yf
    else:
        nums = [float(t) for t in re.split(r"[\s,;]+", " ".join(block["_data"]))
                if t and re.match(r"^[+-]?[\d.]", t)]
        arr = np.array(nums).reshape(-1, 2)
        x, y = arr[:, 0] * xf, arr[:, 1] * yf
    ok = np.isfinite(y)
    meta = {k: v for k, v in block.items() if not k.startswith("_")}
    return Spectrum(x[ok], y[ok], name=block.get("TITLE") or fallback_name, metadata=meta)


def parse_jcamp(content: str, default_name: str = "JCAMP") -> ParseResult:
    blocks = _parse_blocks(content)
    if not blocks:
        raise ValueError("no ##XYDATA or ##XYPOINTS block found")
    spectra = [_block_to_spectrum(b, f"{default_name} #{i + 1}") for i, b in enumerate(blocks)]
    return ParseResult(spectra, "jcamp", metadata={"_format": "JCAMP-DX"})


def load_jcamp(path: str | Path) -> ParseResult:
    path = Path(path)
    return parse_jcamp(read_text_file(path), default_name=path.stem)


def write_jcamp(spectrum: Spectrum, path: str | Path, yunits: str = "ABSORBANCE") -> None:
    x, y = spectrum.wavelengths, spectrum.values
    lines = [
        f"##TITLE={spectrum.name}",
        "##JCAMP-DX=4.24",
        "##DATA TYPE=UV/VIS SPECTRUM",
        "##ORIGIN=Spectro",
        "##OWNER=",
        "##XUNITS=NANOMETERS",
        f"##YUNITS={yunits}",
        "##XFACTOR=1",
        "##YFACTOR=1",
        f"##FIRSTX={x[0]:.6g}",
        f"##LASTX={x[-1]:.6g}",
        f"##NPOINTS={x.size}",
        f"##FIRSTY={y[0]:.8g}",
        "##XYPOINTS=(XY..XY)",
    ]
    lines += [f"{a:.6g}, {b:.8g}" for a, b in zip(x, y)]
    lines.append("##END=")
    Path(path).write_text("\n".join(lines) + "\n", encoding="utf-8")
