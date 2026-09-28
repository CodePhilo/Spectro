"""The Spectrum value object used throughout the core library."""

from __future__ import annotations

import hashlib
from dataclasses import dataclass, field
from typing import Any

import numpy as np


@dataclass
class Spectrum:
    """A single UV/Vis spectrum: absorbance (or any ordinate) versus wavelength.

    Wavelengths are always stored in ascending order.
    """

    wavelengths: np.ndarray
    values: np.ndarray
    name: str = ""
    metadata: dict[str, Any] = field(default_factory=dict)
    # Known concentrations keyed by compound name (e.g. {"Paracetamol": 10.0}).
    concentrations: dict[str, float] = field(default_factory=dict)

    def __post_init__(self) -> None:
        x = np.asarray(self.wavelengths, dtype=float).ravel()
        y = np.asarray(self.values, dtype=float).ravel()
        if x.shape != y.shape:
            raise ValueError(
                f"wavelengths ({x.size}) and values ({y.size}) differ in length"
            )
        if x.size < 2:
            raise ValueError("a spectrum needs at least two points")
        order = np.argsort(x, kind="stable")
        self.wavelengths = x[order]
        self.values = y[order]

    # ------------------------------------------------------------------ helpers
    @property
    def x(self) -> np.ndarray:
        return self.wavelengths

    @property
    def y(self) -> np.ndarray:
        return self.values

    @property
    def step(self) -> float:
        """Median wavelength interval (nm)."""
        return float(np.median(np.diff(self.wavelengths)))

    def copy(self, **changes: Any) -> "Spectrum":
        data = {
            "wavelengths": self.wavelengths.copy(),
            "values": self.values.copy(),
            "name": self.name,
            "metadata": dict(self.metadata),
            "concentrations": dict(self.concentrations),
        }
        data.update(changes)
        return Spectrum(**data)

    def with_values(self, values: np.ndarray, name: str | None = None) -> "Spectrum":
        return self.copy(values=np.asarray(values, dtype=float),
                         name=self.name if name is None else name)

    def value_at(self, wavelength: float) -> float:
        """Linearly interpolated ordinate at ``wavelength``."""
        lo, hi = self.wavelengths[0], self.wavelengths[-1]
        if not lo - 1e-9 <= wavelength <= hi + 1e-9:
            raise ValueError(
                f"{wavelength} nm is outside the spectrum range {lo:g}–{hi:g} nm"
            )
        return float(np.interp(wavelength, self.wavelengths, self.values))

    def region(self, start: float, end: float) -> tuple[np.ndarray, np.ndarray]:
        lo, hi = min(start, end), max(start, end)
        mask = (self.wavelengths >= lo - 1e-9) & (self.wavelengths <= hi + 1e-9)
        if mask.sum() == 0:
            raise ValueError(f"no data points between {lo:g} and {hi:g} nm")
        return self.wavelengths[mask], self.values[mask]

    def resampled(self, grid: np.ndarray) -> "Spectrum":
        grid = np.asarray(grid, dtype=float)
        if grid[0] < self.wavelengths[0] - 1e-9 or grid[-1] > self.wavelengths[-1] + 1e-9:
            raise ValueError("resampling grid extends beyond the spectrum range")
        return self.copy(wavelengths=grid,
                         values=np.interp(grid, self.wavelengths, self.values))

    def data_hash(self) -> str:
        """SHA-256 of the numeric content (used for audit integrity)."""
        h = hashlib.sha256()
        h.update(np.ascontiguousarray(self.wavelengths, dtype="<f8").tobytes())
        h.update(np.ascontiguousarray(self.values, dtype="<f8").tobytes())
        return h.hexdigest()

    def same_grid(self, other: "Spectrum", tol: float = 1e-6) -> bool:
        return (self.wavelengths.shape == other.wavelengths.shape
                and bool(np.allclose(self.wavelengths, other.wavelengths, atol=tol)))


def common_grid(spectra: list[Spectrum]) -> np.ndarray:
    """Return a wavelength grid shared by all spectra.

    If every spectrum already has the same grid it is returned unchanged;
    otherwise the overlapping range is sampled at the coarsest step.
    """
    if not spectra:
        raise ValueError("no spectra given")
    first = spectra[0]
    if all(first.same_grid(s) for s in spectra[1:]):
        return first.wavelengths.copy()
    lo = max(s.wavelengths[0] for s in spectra)
    hi = min(s.wavelengths[-1] for s in spectra)
    if hi <= lo:
        raise ValueError("the spectra do not share a common wavelength range")
    step = max(s.step for s in spectra)
    n = int(np.floor((hi - lo) / step + 1e-9)) + 1
    return lo + step * np.arange(n)


def align(spectra: list[Spectrum]) -> list[Spectrum]:
    """Resample all spectra onto a common grid (no-op if already aligned)."""
    grid = common_grid(spectra)
    return [s if s.wavelengths.shape == grid.shape and np.allclose(s.wavelengths, grid)
            else s.resampled(grid) for s in spectra]


def to_matrix(spectra: list[Spectrum]) -> tuple[np.ndarray, np.ndarray]:
    """Stack spectra as rows of a matrix on a common grid -> (grid, matrix)."""
    aligned = align(spectra)
    return aligned[0].wavelengths, np.vstack([s.values for s in aligned])
