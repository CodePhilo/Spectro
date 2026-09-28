"""Best-effort reader for OLE-container spectra (e.g. Shimadzu UVProbe .spc).

Such files store each data set as a pair of streams whose names end in
``X Data.1`` / ``Y Data.1`` (or similar) containing little-endian doubles.
The layout is not publicly documented, so the reader validates what it finds:
an X stream must decode to a regular, monotonic wavelength axis.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from spectro.core.io.table import ParseResult, _is_axis_like
from spectro.core.spectrum import Spectrum


def _decode(data: bytes) -> np.ndarray | None:
    """Try offsets 0/4/8/16 and float64/float32 until the data looks numeric."""
    for dtype, size in (("<f8", 8), ("<f4", 4)):
        for off in (0, 4, 8, 16):
            n = (len(data) - off) // size
            if n < 3:
                continue
            arr = np.frombuffer(data, dtype=dtype, count=n, offset=off).astype(float)
            if np.all(np.isfinite(arr)) and np.all(np.abs(arr) < 1e7):
                yield arr


def load_ole(path: str | Path) -> ParseResult:
    try:
        import olefile
    except ImportError as exc:  # pragma: no cover
        raise ValueError("reading OLE spectra requires the 'olefile' package") from exc
    path = Path(path)
    ole = olefile.OleFileIO(str(path))
    try:
        streams = ["/".join(s) for s in ole.listdir(streams=True)]
        spectra: list[Spectrum] = []
        for xs in streams:
            leaf = xs.rsplit("/", 1)[-1].lower()
            if "x data" not in leaf:
                continue
            ys = xs[: len(xs) - len(leaf)] + xs.rsplit("/", 1)[-1].replace("X", "Y", 1)
            if ys not in streams:
                continue
            xraw, yraw = ole.openstream(xs).read(), ole.openstream(ys).read()
            found = None
            for x in _decode(xraw):
                if not _is_axis_like(x):
                    continue
                for y in _decode(yraw):
                    if y.size == x.size:
                        found = (x, y)
                        break
                if found:
                    break
            if found:
                parent = xs.rsplit("/", 1)[0]
                name = path.stem if not spectra else f"{path.stem} #{len(spectra) + 1}"
                spectra.append(Spectrum(found[0], found[1], name=name,
                                        metadata={"ole_stream": parent}))
    finally:
        ole.close()
    if not spectra:
        raise ValueError(
            "this is a proprietary OLE container whose data could not be located. "
            "Please export it from the instrument software as CSV, TXT, JCAMP-DX "
            "or GRAMS SPC.")
    return ParseResult(spectra, "ole", metadata={"_format": "OLE container"})
