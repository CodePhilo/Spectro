"""Shared helper: store UI screenshots as compact 256-colour PNGs."""

from pathlib import Path


def compact_png(path: str | Path) -> None:
    """Re-encode a screenshot with a 256-colour palette (≈ 3× smaller, no
    visible loss for UI captures). Leaves the file unchanged without Pillow."""
    try:
        from PIL import Image
    except ImportError:
        return
    im = Image.open(path).convert("RGB")
    im.quantize(colors=256, method=Image.Quantize.MEDIANCUT,
                dither=Image.Dither.NONE).save(path, optimize=True)
