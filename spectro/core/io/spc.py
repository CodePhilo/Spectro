"""Thermo/Galactic SPC binary files (new LSB 0x4B, MSB 0x4C, and old 0x4D).

SPC is the open exchange format written by most UV/Vis instrument software
(GRAMS-compatible export), independent of instrument vendor.
"""

from __future__ import annotations

import struct
from pathlib import Path

import numpy as np

from spectro.core.io.table import ParseResult
from spectro.core.spectrum import Spectrum

TSPREC, TCGRAM, TMULTI, TRANDM, TORDRD, TALABS, TXYXYS, TXVALS = (
    0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80)

_XTYPES = {0: "Arbitrary", 1: "Wavenumber (cm-1)", 2: "Micrometers", 3: "Nanometers",
           4: "Seconds", 5: "Minutes", 6: "Hertz", 13: "Raman shift"}
_YTYPES = {0: "Arbitrary", 1: "Interferogram", 2: "Absorbance", 3: "Kubelka-Munk",
           4: "Counts", 11: "Log(1/R)", 128: "Transmission", 129: "Reflectance"}


def is_spc(data: bytes) -> bool:
    return len(data) >= 256 and data[1] in (0x4B, 0x4C, 0x4D)


def _y_values(raw: bytes, off: int, npts: int, exp: int, prec16: bool,
              endian: str, old: bool = False) -> tuple[np.ndarray, int]:
    if exp == -128 and not old:
        y = np.frombuffer(raw, dtype=f"{endian}f4", count=npts, offset=off).astype(float)
        return y, off + 4 * npts
    if prec16:
        ints = np.frombuffer(raw, dtype=f"{endian}i2", count=npts, offset=off).astype(float)
        return ints * (2.0 ** exp) / 2.0 ** 16, off + 2 * npts
    if old:
        # old format stores 32-bit integers with the 16-bit words swapped
        words = np.frombuffer(raw, dtype="<u2", count=2 * npts, offset=off).reshape(-1, 2)
        ints = ((words[:, 0].astype(np.uint32) << 16) | words[:, 1]).view(np.int32)
        return ints.astype(float) * (2.0 ** exp) / 2.0 ** 32, off + 4 * npts
    ints = np.frombuffer(raw, dtype=f"{endian}i4", count=npts, offset=off).astype(float)
    return ints * (2.0 ** exp) / 2.0 ** 32, off + 4 * npts


def parse_spc(raw: bytes, default_name: str = "SPC") -> ParseResult:
    if not is_spc(raw):
        raise ValueError("not an SPC file")
    flags, version = raw[0], raw[1]
    if version == 0x4D:
        return _parse_old(raw, default_name)
    e = "<" if version == 0x4B else ">"
    (fexp,) = struct.unpack_from("b", raw, 3)
    (fnpts,) = struct.unpack_from(f"{e}i", raw, 4)
    ffirst, flast = struct.unpack_from(f"{e}dd", raw, 8)
    (fnsub,) = struct.unpack_from(f"{e}i", raw, 24)
    fxtype, fytype = raw[28], raw[29]
    comment = raw[88:218].split(b"\x00")[0].decode("latin-1", "replace").strip()
    fnsub = max(fnsub, 1) if flags & TMULTI else 1

    off = 512
    common_x = None
    if flags & TXVALS and not flags & TXYXYS:
        common_x = np.frombuffer(raw, dtype=f"{e}f4", count=fnpts, offset=off).astype(float)
        off += 4 * fnpts
    elif not flags & TXYXYS:
        common_x = np.linspace(ffirst, flast, fnpts)

    spectra = []
    for k in range(fnsub):
        if off + 32 > len(raw):
            break
        subexp = struct.unpack_from("b", raw, off + 1)[0]
        (subnpts,) = struct.unpack_from(f"{e}i", raw, off + 16)
        off += 32
        exp = subexp if flags & TMULTI else fexp
        if flags & TXYXYS:
            npts = subnpts
            x = np.frombuffer(raw, dtype=f"{e}f4", count=npts, offset=off).astype(float)
            off += 4 * npts
        else:
            npts, x = fnpts, common_x
        y, off = _y_values(raw, off, npts, exp, bool(flags & TSPREC), e)
        name = default_name if fnsub == 1 else f"{default_name} #{k + 1}"
        spectra.append(Spectrum(x, y, name=name, metadata={
            "x_units": _XTYPES.get(fxtype, str(fxtype)),
            "y_units": _YTYPES.get(fytype, str(fytype)),
            "comment": comment}))
    return ParseResult(spectra, "spc", metadata={"_format": "SPC", "comment": comment})


def _parse_old(raw: bytes, default_name: str) -> ParseResult:
    flags = raw[0]
    (oexp,) = struct.unpack_from("<h", raw, 2)
    onpts, ofirst, olast = struct.unpack_from("<fff", raw, 4)
    npts = int(onpts)
    x = np.linspace(ofirst, olast, npts)
    off = 256
    spectra = []
    k = 0
    while off + 32 <= len(raw):
        off += 32
        y, off = _y_values(raw, off, npts, oexp, bool(flags & TSPREC), "<", old=True)
        k += 1
        spectra.append(Spectrum(x, y, name=f"{default_name} #{k}"))
        if not flags & TMULTI:
            break
    if len(spectra) == 1:
        spectra[0].name = default_name
    return ParseResult(spectra, "spc", metadata={"_format": "SPC (old)"})


def load_spc(path: str | Path) -> ParseResult:
    path = Path(path)
    return parse_spc(path.read_bytes(), default_name=path.stem)


def write_spc(spectra: list[Spectrum], path: str | Path) -> None:
    """Write spectra to a new-format SPC file (float Y, shared or XY X values)."""
    if not spectra:
        raise ValueError("nothing to write")
    same = all(s.same_grid(spectra[0]) for s in spectra)
    flags = TMULTI if len(spectra) > 1 else 0
    flags |= TXVALS if same else (TXVALS | TXYXYS)
    x0 = spectra[0].wavelengths
    head = bytearray(512)
    struct.pack_into("<BBbb", head, 0, flags, 0x4B, 0, -128)
    struct.pack_into("<i", head, 4, x0.size if same else 0)
    struct.pack_into("<dd", head, 8, float(x0[0]), float(x0[-1]))
    struct.pack_into("<i", head, 24, len(spectra))
    head[28], head[29] = 3, 2
    body = bytearray()
    if same:
        body += x0.astype("<f4").tobytes()
    for i, s in enumerate(spectra):
        sub = bytearray(32)
        struct.pack_into("<Bbh", sub, 0, 0, -128, i)
        struct.pack_into("<i", sub, 16, 0 if same else s.wavelengths.size)
        body += sub
        if not same:
            body += s.wavelengths.astype("<f4").tobytes()
        body += s.values.astype("<f4").tobytes()
    Path(path).write_bytes(bytes(head) + bytes(body))
