import os

import numpy as np
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from spectro.core.spectrum import Spectrum

GRID = np.arange(200.0, 400.01, 0.5)


def band(center, width, height=1.0):
    return height * np.exp(-0.5 * ((GRID - center) / width) ** 2)


# Unit-concentration (1 µg/mL) spectra of three model compounds.
# X absorbs 200–300 nm, Y extends further (alone above ~310 nm), Z in between.
PURE = {
    "X": 0.040 * band(245, 12) + 0.015 * band(270, 8),
    "Y": 0.030 * band(275, 22) + 0.020 * band(330, 18),
    "Z": 0.035 * band(260, 10) + 0.010 * band(225, 8),
}


def mixture(conc: dict, name: str = "", noise: float = 0.0, seed: int = 0) -> Spectrum:
    y = sum(PURE[k] * v for k, v in conc.items())
    if noise:
        y = y + np.random.default_rng(seed).normal(0, noise, y.size)
    return Spectrum(GRID, y, name=name or str(conc),
                    concentrations={k: float(v) for k, v in conc.items()})


@pytest.fixture
def pure():
    return {k: Spectrum(GRID, v, name=k, concentrations={k: 1.0}) for k, v in PURE.items()}


@pytest.fixture
def binary_set():
    """Calibration (pure standards) and lab-prepared mixtures of X + Y."""
    xs = [mixture({"X": c, "Y": 0.0}, f"X {c}") for c in (4, 8, 12, 16, 20)]
    ys = [mixture({"X": 0.0, "Y": c}, f"Y {c}") for c in (4, 8, 12, 16, 20)]
    mixes = [mixture({"X": a, "Y": b}, f"mix {a}+{b}") for a, b in
             [(6, 10), (10, 6), (14, 14), (8, 16), (18, 5)]]
    return xs, ys, mixes
