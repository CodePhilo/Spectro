# Spectro

A Windows desktop application for organizing UV/Vis spectral trials in a
pharmaceutical analytical chemistry lab. It supports every mainstream
spectrophotometric method for **simultaneous determination of compounds in
mixtures**, and it **logs every operation** in a tamper-evident audit trail.

- **Instrument-independent import.** Text/CSV from any spectrophotometer (any
  delimiter, decimal comma, header lines; column, XY-pair or row layouts),
  Excel (`.xlsx`, `.xls`), JCAMP-DX, GRAMS SPC, and OLE-container `.spc`
  (e.g. Shimadzu UVProbe).
- **Non-destructive processing.** Raw data can never be changed. Each
  processing result is stored as a new *derived* spectrum that records its
  parent and its exact pipeline, so it can be recomputed and verified at any
  time.
- **Audit trail.** Every import, edit, processing step, calculation, export
  and archive is logged with the user, workstation, UTC time, parameters,
  before/after values, reason for change and data hashes. Entries are
  hash-chained and the database rejects edits or deletions, so any tampering
  is detected by *Audit → Verify data integrity*.

See [`docs/PLAN.md`](docs/PLAN.md) for the full plan and the review checklist.

## Install and run (from source)

```bash
python -m pip install -r requirements.txt
python run_spectro.py            # or: python -m spectro
```

Python 3.10+ is required. To build a standalone Windows program
(`dist\Spectro\Spectro.exe`), run `build_windows.bat`. The GitHub Actions
workflow `.github/workflows/windows.yml` runs the tests, builds the `.exe`
and uploads it as the `Spectro-windows` artifact on every push.

## Demo data

**Help → Open demo project…** builds a complete project in about a second:

- **Binary trial (Panadol Extra):** paracetamol and caffeine standards
  (2–20 µg/mL), 8 laboratory mixtures, 3 tablet solutions and a solvent blank.
  It includes ratio and D1 spectra, ratio difference methods for both drugs and
  Vierordt equations.
- **Ternary trial (Excedrin):** a 25-mixture Brereton training set, 5
  validation mixtures and 3 tablet solutions, with a cross-validated PLS2
  model.

Recoveries are 98.6–102 % with RSD ≈ 0.5–1 %. The project's audit trail shows
every step used to build it. **Help → Copy demo data files…** copies the
underlying files so you can try the import dialog on seven different
instrument-export styles (see [`spectro/demo_data/README.md`](spectro/demo_data/README.md)).
The same project can be built from the command line with
`python -m spectro.demo "Spectro demo.spectro"`.

> The demo spectra are **simulated** from published band positions and
> absorptivities, with realistic noise, baseline and preparation errors. They
> are not measured data. Regenerate them with `python tools/make_demo_data.py`.

## Workflow gallery

[`docs/workflows/`](docs/workflows/README.md) shows 23 annotated screenshots of
the main workflows on the demo project, each with the benefit it brings. To
regenerate them (for example after changing the UI), run:

```bash
python tools/workflow_screenshots.py            # → docs/workflows/ (index.html + README.md)
python tools/workflow_screenshots.py --out shots --only optimizer,univariate
```

## Typical workflow

1. **File → New project…** creates one `.spectro` file (a SQLite database) that
   holds everything.
2. **Edit → Compounds…** adds the analytes (e.g. PAR, CAF) and their units.
3. **File → Import spectra…** adds files. Check the preview, type or paste
   concentrations (or *Fill concentrations from names*), then choose the trial
   and role (`standard`, `calibration`, `mixture`, `sample`, `blank`,
   `divisor`…).
4. **Process → Processing pipeline…** chains smoothing, baseline correction,
   derivatives, ratio spectra, subtraction and so on, with a live preview.
5. **Methods** calibrates and determines concentrations. Results show
   % recovery, mean, SD and RSD, and can be saved.
6. **Tools → Validation statistics** covers linearity, LOD/LOQ, accuracy,
   precision/ANOVA, t/F tests, the interval hypothesis test and standard
   addition.
7. **File → Trial report** exports a PDF or HTML report with spectra,
   lineage, methods, results and the audit trail.

## Methods

| Family | Methods |
|---|---|
| Pre-processing | crop, resample, Savitzky–Golay, moving average, Whittaker; baseline offset, two-point, iterative polynomial, ALS; blank subtraction and spectral arithmetic; normalisation (max, area, vector, at λ, range), SNV, %T ↔ A |
| Zero order | direct λmax; Vierordt simultaneous equations and multi-wavelength least squares; absorbance ratio (Q-analysis); dual wavelength; induced dual wavelength; absorbance subtraction; advanced absorbance subtraction; amplitude modulation; induced amplitude modulation; area under curve (single and equations); bivariate (with Kaiser λ selection); H-point standard addition |
| Derivative | D1–D4 (finite difference Δλ or Savitzky–Golay, scaling factor); zero-crossing; peak-to-peak; dual wavelength in derivative mode; factorized zero-order |
| Ratio spectra | derivative ratio (DD1); ratio difference; mean centering of ratio spectra; successive derivative ratio, double divisor (sum divisor) and dual amplitude difference (ternary); ratio subtraction and extended ratio subtraction; successive spectrum subtraction (ternary); constant multiplication; constant center; spectrum subtraction; plateau/constant value |
| Chemometrics | CLS, ILS/MLR, PCR, PLS-1, PLS-2, MCR-ALS (non-negativity + correlation constraint), ANN (MLP on PCA scores), SVR; LOO / venetian / k-fold / contiguous cross-validation with Haaland–Thomas selection of components; iPLS and GA interval selection; VIP; Hotelling T² vs Q with 95/99 % limits, outlier flags and one-click exclusion; fitted models saved in the project and exportable as `.spmodel` files (with embedded divisor spectra and checksums); progress bars with cancel for CV, iPLS, GA and fitting; Brereton 5-level multifactor calibration design |
| Method optimizer | *Tools → Method optimizer*: from the pure standards it predicts, for every compound, the interference, noise and ±λ-robustness of every processing strategy (zero order, derivative zero-crossing, dual / induced dual wavelength, ratio difference, derivative ratio, mean centering, ratio subtraction / constant multiplication, dual amplitude difference, double divisor, Vierordt, bivariate, CLS, PLS) over all wavelengths and λ pairs; the best settings of each are verified end-to-end on simulated and laboratory mixtures and ranked; any row can be saved as a ready, calibrated method |
| Finder tools | zero-crossing points, isoabsorptive points, maxima/minima, ratio-spectrum plateaus, equal-amplitude λ pairs |
| Method studies | standard addition inside univariate methods (recovery of each spike + extrapolation); robustness study (each parameter ± Δ, re-calibrated and re-assayed, % deviation) |
| Output | results to Excel (one sheet per result + index + methods), any table to Excel (right-click), publication figures from any plot (right-click / *Export figure*: size in mm, 300–1200 dpi TIFF/PNG or vector SVG/PDF/EPS, fonts, black-and-white mode); plot views: overlay, stacked, difference, normalized |
| Validation (ICH Q2(R2)) | regression with SD and 95 % CI of slope and intercept, LOD/LOQ (σ of intercept or Sy/x), lack-of-fit, recovery, %RSD, one-way ANOVA, Student's t, F-test, interval hypothesis, standard addition, RMSEP/bias/SEP |
| Greenness | AGREE (weighted 12 principles with pictogram), Analytical Eco-Scale, GAPI |

Each method is tested in `tests/test_methods.py` against synthetic binary and
ternary mixtures with known concentrations.

## Development

```bash
python -m pip install -r requirements.txt pytest hypothesis
QT_QPA_PLATFORM=offscreen python -m pytest -q
```

The suite checks every calculation against independent references (SciPy,
scikit-learn, closed-form results), runs property-based and fuzz tests,
malformed-input, data-integrity and end-to-end tests, and opens every dialog on
empty projects. See [`docs/TESTING.md`](docs/TESTING.md), which also lists the
problems these tests found and how each was fixed.

Layout: `spectro/core` (algorithms and file formats, no UI),
`spectro/storage` (project database and audit trail), `spectro/ui` (PySide6
windows and dialogs), `tests/`.
