"""Generate the demo data set shipped in ``spectro/demo_data``.

The spectra are *simulated* from published band positions and absorptivities
(methanol, 1 cm cell) and then degraded like real measurements:

* Paracetamol (PAR) λmax ≈ 243 nm, A(1 %, 1 cm) ≈ 668
* Caffeine (CAF)    λmax ≈ 273 nm, A(1 %, 1 cm) ≈ 502
* Aspirin (ASA)     λmax ≈ 229 nm, A(1 %, 1 cm) ≈ 480, weak band ≈ 277 nm

Realism added per spectrum: volumetric (preparation) error ≈ 0.3 %, wavelength
jitter ≈ 0.05 nm, baseline offset and tilt, photometric noise that increases
strongly below 215 nm (low lamp energy / solvent cut-off) plus ≈ 0.2 %
proportional noise, and 4-decimal rounding as in instrument exports.

Files are written in the export styles of different instrument software so the
importer's layout detection can be tried on each:

binary/ (Panadol Extra: paracetamol 500 mg + caffeine 65 mg per tablet)
  shimadzu_txt/PAR_*.txt      single spectrum, quoted header, descending λ, 0.5 nm
  CAF_standards_xy_pairs.csv  Cary-WinUV style (λ, Abs) column pairs, 1 nm
  lab_mixtures.xlsx           Excel, λ + one column per mixture, 1 nm
  panadol_extra_tablets.csv   semicolon + decimal comma, metadata header, 1 nm
  methanol_blank.txt          whitespace-delimited with '#' comments
ternary/ (Excedrin: aspirin 250 mg + paracetamol 250 mg + caffeine 65 mg)
  training_set_brereton.jdx   25 JCAMP-DX blocks (5-level multifactor design)
  validation_mixtures.tsv     one spectrum per row (λ in first row)
  excedrin_tablets.spc        GRAMS SPC multi-file

Run:  python tools/make_demo_data.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from spectro.core.io import export_excel, write_spc
from spectro.core.multicomponent import design_concentrations
from spectro.core.spectrum import Spectrum

OUT = Path(__file__).resolve().parents[1] / "spectro" / "demo_data"
RNG = np.random.default_rng(20260928)

# (centre nm, absorptivity per µg/mL at the band maximum, σ nm)
BANDS = {
    "PAR": [(243.0, 0.0668, 14.0), (201.0, 0.085, 8.5), (285.0, 0.004, 10.0)],
    "CAF": [(273.0, 0.0502, 12.0), (206.0, 0.105, 9.0), (228.0, 0.008, 8.0)],
    "ASA": [(229.0, 0.0480, 10.0), (202.0, 0.090, 7.0), (277.0, 0.0056, 9.0)],
}
TRUE_LAMBDA_MAX = {"PAR": 243.0, "CAF": 273.0, "ASA": 229.0}


def absorptivity(compound: str, wl: np.ndarray) -> np.ndarray:
    return sum(h * np.exp(-0.5 * ((wl - c) / s) ** 2) for c, h, s in BANDS[compound])


def measure(conc: dict[str, float], wl: np.ndarray, noise_scale: float = 1.0) -> np.ndarray:
    """One simulated measurement of a solution against a solvent blank."""
    shift = RNG.normal(0, 0.05)
    prep = {k: v * RNG.normal(1.0, 0.003) for k, v in conc.items()}
    y = sum(absorptivity(k, wl + shift) * v for k, v in prep.items())
    y = y + RNG.normal(0, 0.0012) + RNG.normal(0, 4e-6) * (wl - 300)
    sigma = noise_scale * (0.0004 * (1 + 12 * np.exp(-(wl - 200) / 7)) + 0.002 * np.abs(y))
    y = y + RNG.normal(0, 1, wl.size) * sigma
    return np.round(y, 4)


def name_of(conc: dict[str, float]) -> str:
    return " ".join(f"{k} {v:g}" for k, v in conc.items())


def write_shimadzu_txt(path: Path, title: str, wl: np.ndarray, y: np.ndarray) -> None:
    lines = [f'"{title}"', '"Wavelength nm.","Abs."']
    lines += [f"{a:.1f},{b:.4f}" for a, b in zip(wl[::-1], y[::-1])]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_xy_pairs(path: Path, spectra: list[tuple[str, np.ndarray, np.ndarray]]) -> None:
    head1 = ",".join(f"{n}," for n, _, _ in spectra)
    head2 = ",".join("Wavelength (nm),Abs" for _ in spectra) + ","
    rows = [",".join(f"{x[i]:.1f},{y[i]:.4f}" for _, x, y in spectra) + ","
            for i in range(spectra[0][1].size)]
    path.write_text("\n".join([head1, head2] + rows) + "\n", encoding="utf-8")


def write_jcamp_blocks(path: Path, blocks: list[tuple[str, np.ndarray, np.ndarray]]) -> None:
    out = []
    for title, x, y in blocks:
        out += [f"##TITLE={title}", "##JCAMP-DX=4.24", "##DATA TYPE=UV/VIS SPECTRUM",
                "##ORIGIN=Spectro demo data (simulated)", "##OWNER=public domain",
                "##XUNITS=NANOMETERS", "##YUNITS=ABSORBANCE", "##XFACTOR=1", "##YFACTOR=1",
                f"##FIRSTX={x[0]:g}", f"##LASTX={x[-1]:g}", f"##NPOINTS={x.size}",
                "##XYPOINTS=(XY..XY)"]
        out += [f"{a:g}, {b:.4f}" for a, b in zip(x, y)]
        out.append("##END=")
    path.write_text("\n".join(out) + "\n", encoding="utf-8")


def main() -> None:
    binary = OUT / "binary"
    ternary = OUT / "ternary"
    (binary / "shimadzu_txt").mkdir(parents=True, exist_ok=True)
    ternary.mkdir(parents=True, exist_ok=True)
    wl05 = np.round(np.arange(200.0, 400.01, 0.5), 1)
    wl1 = np.arange(200.0, 401.0, 1.0)

    # ---------------- binary: paracetamol + caffeine ----------------
    for c in (2, 4, 6, 8, 10, 12, 15, 20):
        write_shimadzu_txt(binary / "shimadzu_txt" / f"PAR_{c:02d}.txt", f"PAR {c} ug/mL",
                           wl05, measure({"PAR": c}, wl05))
    caf = [(f"CAF {c}", wl1[::-1], measure({"CAF": c}, wl1[::-1]))
           for c in (2, 4, 6, 8, 10, 12, 15, 20)]
    write_xy_pairs(binary / "CAF_standards_xy_pairs.csv", caf)

    mixes = [{"PAR": 10, "CAF": 1.3}, {"PAR": 15, "CAF": 2}, {"PAR": 8, "CAF": 8},
             {"PAR": 12, "CAF": 4}, {"PAR": 6, "CAF": 10}, {"PAR": 4, "CAF": 12},
             {"PAR": 16, "CAF": 6}, {"PAR": 10, "CAF": 5}]
    export_excel([Spectrum(wl1, measure(m, wl1), name=name_of(m)) for m in mixes],
                 binary / "lab_mixtures.xlsx")

    # Tablets: diluted so CAF is in range: 15.38 µg/mL PAR + 2.0 µg/mL CAF
    # (500:65 ratio); true content 99.2 / 100.6 / 99.8 % of label claim
    content = [0.992, 1.006, 0.998]
    tabs = [measure({"PAR": 15.3846 * f, "CAF": 2.0 * f}, wl1) for f in content]
    lines = ["Sample;Panadol Extra tablets (500 mg PAR / 65 mg CAF)", "Solvent;Methanol",
             "Dilution;nominal 15.38 ug/mL PAR + 2.0 ug/mL CAF", "Instrument;double-beam UV-Vis, 1 cm quartz",
             "", "Wavelength;Tablet 1;Tablet 2;Tablet 3"]
    lines += [f"{w:.1f};" + ";".join(f"{t[i]:.4f}" for t in tabs) for i, w in enumerate(wl1)]
    (binary / "panadol_extra_tablets.csv").write_text(
        "\n".join(lines).replace(".", ",") + "\n", encoding="utf-8")

    blank = measure({}, wl1)
    (binary / "methanol_blank.txt").write_text(
        "# Methanol vs methanol (baseline check)\n# nm  Abs\n"
        + "\n".join(f"{w:.1f}  {b:.4f}" for w, b in zip(wl1, blank)) + "\n", encoding="utf-8")

    # ---------------- ternary: aspirin + paracetamol + caffeine ----------------
    design = design_concentrations({"PAR": 10, "CAF": 6, "ASA": 12},
                                   {"PAR": 2, "CAF": 1.5, "ASA": 3})
    write_jcamp_blocks(ternary / "training_set_brereton.jdx",
                       [(f"D{i + 1:02d} {name_of(c)}", wl1, measure(c, wl1))
                        for i, c in enumerate(design)])

    val = [{"PAR": 9, "CAF": 5, "ASA": 14}, {"PAR": 12, "CAF": 7, "ASA": 10},
           {"PAR": 7, "CAF": 4, "ASA": 16}, {"PAR": 13, "CAF": 8, "ASA": 8},
           {"PAR": 11, "CAF": 3.5, "ASA": 11}]
    rows = ["Sample\t" + "\t".join(f"{w:g}" for w in wl1)]
    rows += [f"V{i + 1} {name_of(c)}\t" + "\t".join(f"{v:.4f}" for v in measure(c, wl1))
             for i, c in enumerate(val)]
    (ternary / "validation_mixtures.tsv").write_text("\n".join(rows) + "\n", encoding="utf-8")

    # Excedrin: 250 mg ASA + 250 mg PAR + 65 mg CAF → 12 / 12 / 3.12 µg/mL
    exc = [Spectrum(wl1, measure({"ASA": 12 * f, "PAR": 12 * f, "CAF": 3.12 * f}, wl1),
                    name=f"Excedrin tablet {i + 1}") for i, f in enumerate([1.004, 0.995, 1.001])]
    write_spc(exc, ternary / "excedrin_tablets.spc")
    print(f"demo data written to {OUT}")


if __name__ == "__main__":
    main()
