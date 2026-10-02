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

[`docs/workflows/`](docs/workflows/README.md) shows 26 annotated screenshots of
the main workflows on the demo project, each with the benefit it brings, and
[`docs/methods/`](docs/methods/README.md) is the method guide for beginners
(44 methods, each with the idea in plain words, principle, what to measure,
step-by-step use, checks and common mistakes, a "why it works" figure,
screenshots of every step, and a sample Excel workbook in
[`docs/methods/data/`](docs/methods/data) to import and reproduce the
example; the workbooks also show the calibration as Excel formulas). Both are generated from the running app; regenerate them after
changing the UI:

```bash
python tools/workflow_screenshots.py            # → docs/workflows/ (index.html + README.md)
python tools/workflow_screenshots.py --out shots --only optimizer,univariate
python tools/method_guide.py                    # → docs/methods/ (README.md + index.html)
python tools/method_guide.py --only rd,am       # a few entries
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

Methods are grouped by **which spectrum is manipulated** (zero order,
derivative, ratio) and **what they yield** (a signal at chosen wavelengths,
or the recovered spectrum of each component), with simultaneous-equation,
chemometric and standard-addition methods kept apart. The same families are
used in the app (template list, dialogs) and in the
**[method guide](docs/methods/README.md)**, which explains how to apply every
method with a worked example and a screenshot (*Methods → Method guide*).

| Family | Methods | Where in Spectro |
|---|---|---|
| Zero-order spectra (absorbance) | direct λmax; dual wavelength (DW); induced dual wavelength (IDW); area under curve (single component); absorbance ratio (Q-analysis); absorbance subtraction (AS); advanced absorbance subtraction (AAS); absorption factor / successive absorption factor (MAFM) | Univariate templates; *Binary two-signal methods*; *Progressive resolution* (MAFM) |
| Simultaneous equations | Vierordt and multi-wavelength least squares; bivariate (Kaiser λ selection); area under curve (binary / ternary, Cramer's rule) | *Equation methods* |
| Derivative spectra | D1–D4 (finite difference Δλ or Savitzky–Golay, scaling factor); zero-crossing; peak-to-peak; dual wavelength in derivative mode (D1 DWL) | Univariate templates |
| Ratio spectra — derivative and mean centering | derivative ratio (DD1 / DR1); derivative ratio of D1 spectra (D1 DR); mean centering of ratio spectra (binary and ternary); successive derivative ratio; double divisor ratio derivative | Univariate templates |
| Ratio spectra — amplitude methods | ratio difference (RD); dual amplitude difference; constant value (plateau); concentration value; amplitude modulation (AM); induced amplitude modulation (IAM); constant value via amplitude difference (CV-AD); advanced amplitude centering (AAC, partial / complete overlap); modified amplitude center (MACM); ratio difference–isoabsorptive (RIDSS) | Univariate templates; *Binary two-signal methods* (AM, IAM); *Progressive resolution* (CV-AD, AAC, MACM, RIDSS) |
| Spectrum resolution (recover each spectrum) | ratio subtraction (RS); constant multiplication and SS-CM; extended ratio subtraction (ERS); successive ratio subtraction (SRS); successive spectrum subtraction; spectrum subtraction (SS); factorized zero-order (FZM); constant center (CC); derivative transformation (DT, DT-SS); derivative subtraction (DS) and DS-CM | Univariate templates (and *Process* to store the recovered spectra) |
| Multivariate (chemometrics) | CLS, ILS/MLR, PCR, PLS-1, PLS-2, MCR-ALS (non-negativity + correlation constraint), ANN (MLP on PCA scores), SVR; LOO / venetian / k-fold / contiguous cross-validation with Haaland–Thomas selection of components; iPLS and GA interval selection; VIP; Hotelling T² vs Q with 95/99 % limits, outlier flags and one-click exclusion; Brereton 5-level multifactor calibration design | *Chemometrics*; *Tools → Calibration design* |
| Standard addition and enrichment | standard addition inside any univariate method (recovery of each spike + extrapolation); H-point standard addition (HPSAM); sample enrichment by spiking (found − added) or spectrum addition | Univariate (*Standard addition…*, *Enrichment added*); *Binary two-signal methods* (HPSAM); *Process → Add spectrum* |

| Tools and services | |
|---|---|
| Pre-processing | crop, resample, Savitzky–Golay, moving average, Whittaker; baseline offset, two-point, iterative polynomial, ALS; blank subtraction and spectral arithmetic; normalisation (max, area, vector, at λ, range, unit concentration), SNV, %T ↔ A |
| Method optimizer | *Tools → Method optimizer*: from the pure standards it predicts, for every compound, the interference, noise and ±λ-robustness of every processing strategy (zero order, derivative zero-crossing, dual / induced dual wavelength, ratio difference, derivative ratio, mean centering, ratio subtraction / constant multiplication, dual amplitude difference, double divisor, amplitude centering with every compound as divisor, with and without a plateau, successive absorption factor in every compound order, Vierordt, bivariate, CLS, PLS), optionally each again with Savitzky–Golay smoothing of chosen widths (nm) and with Savitzky–Golay derivatives; the best settings of each are verified end-to-end on simulated and laboratory mixtures and ranked; any row can be saved as a ready, calibrated method |
| Finder tools | zero-crossing points, isoabsorptive points, maxima/minima, ratio-spectrum plateaus, equal-amplitude λ pairs |
| Saved methods | every method type is saved calibrated, applied to new samples, and exported as a `.spmodel` file with its divisor/reference spectra embedded and checksummed |
| Method studies | robustness study (each parameter ± Δ, re-calibrated and re-assayed, % deviation) |
| Output | results to Excel (one sheet per result + index + methods), any table to Excel (right-click), publication figures from any plot (right-click / *Export figure*: size in mm, 300–1200 dpi TIFF/PNG or vector SVG/PDF/EPS, fonts, black-and-white mode); plot views: overlay, stacked, difference, normalized |
| Validation (ICH Q2(R2)) | regression with SD and 95 % CI of slope and intercept, LOD/LOQ (σ of intercept or Sy/x), lack-of-fit, recovery, %RSD, one-way ANOVA, Student's t, F-test (one- or two-tailed critical value), interval hypothesis, standard addition, RMSEP/bias/SEP; t, F and ANOVA also from published summary values (mean, SD, n) to compare with a reported or official method |
| Greenness | AGREE (weighted 12 principles with pictogram), Analytical Eco-Scale, GAPI |

Each method is tested in `tests/test_methods.py` against synthetic binary and
ternary mixtures with known concentrations. `tests/test_literature.py`
validates the app against nine published papers: it recomputes their t, F
and ANOVA values and runs each paper's method on spectra built to the paper's
conditions. See [`docs/LITERATURE.md`](docs/LITERATURE.md) for method coverage
and the problems found (in the app and in the papers).

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
