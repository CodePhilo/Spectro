# Testing Spectro

```bash
python -m pip install -r requirements.txt pytest hypothesis
QT_QPA_PLATFORM=offscreen python -m pytest -q          # whole suite (~30 s)
QT_QPA_PLATFORM=offscreen python -m pytest -q -W error::RuntimeWarning   # strict
```

The suite is organised by *what kind of error it can catch*.

| File | Catches | How |
|---|---|---|
| `test_reference.py` | **wrong numbers** | every statistic and algorithm compared with an independent implementation: regression vs `scipy.stats.linregress` (slope, intercept, r, SDs, CIs, LOD/LOQ), lack-of-fit, t/F/ANOVA vs SciPy, derivatives vs analytic polynomials, PLS vs scikit-learn, PCR/ILS vs linear algebra, CLS/MCR-ALS exact on noiseless data, frozen ANN/SVR vs scikit-learn, VIP and Hotelling identities, LOO CV vs a manual loop, Q-analysis vs simultaneous equations, greenness formulas; statistics refuse insufficient input |
| `test_properties.py` | **hidden assumptions** | Hypothesis property tests: every processing step used by the optimizer is linear; every registered operation is classified; results do not depend on the concentration unit, sample order, wavelength spacing or direction; 300 random pipelines and 150 random measurements fail only with clear `ValueError`s and never leak NaN; random spectra round-trip through every file format; T²/Q limits flag ≈ 5 % / 1 % of normal new samples; the optimizer finds a known best method |
| `test_integrity.py` | **bad input, data integrity, end-to-end** | unusable files, BOM/CRLF/unicode minus/footers, duplicate wavelengths, missing Excel cells, a 120-spectrum file, unicode names; every mutation writes exactly one audit entry and no-op edits none; failed operations leave no trace; tampering with data or the log is detected; two app instances on one file; foreign files refused; chained derivations replay after archiving references; stored demo results equal a from-scratch NumPy calculation; Excel export equals stored numbers; PDF/HTML report |
| `test_literature.py` | **disagreement with the published literature** | 88 published t/F values, ANOVAs and recovery statistics from nine papers recomputed from their summary values; every method of those papers (AAC, MACM, RIDSS, CV-AD, MAFM, ternary AUC and MCR, double divisor, SRS, DS-CM, D1 DR, DT, CM-SS, concentration value, AM with unit divisor, enrichment) run on spectra built to each paper's conditions; see [`LITERATURE.md`](LITERATURE.md) |
| `test_docs.py` | **docs drifting from the app** | every template belongs to a method family consistent with its processing; the template list shows family headings; the method guide explains every template and family, all its images exist and its examples recover 95–105 %; both screenshot generators run; every method has a sample workbook that imports with the right roles and concentrations, and its Excel-formula calculation reproduces the app's results |
| `test_ui_edge.py` | **crashes** | every dialog action on an empty project, a project without compounds and a project with data; main-window actions without a project or selection |
| `test_methods.py`, `test_optimizer.py`, `test_io.py`, `test_project.py`, `test_tables.py`, `test_demo.py`, `test_ui.py` | feature tests | every method recovers known concentrations; feature flows |

## Problems found by this suite (and fixed)

| # | Area | Problem | Fix |
|---|---|---|---|
| 1 | Derivatives | Savitzky–Golay derivative assumed equal spacing: on a grid mixing 1 nm and 0.5 nm steps the derivative was **off by 2×** in part of the range | resample to a regular grid internally |
| 2 | Derivatives | Δλ-difference derivatives **fabricated values** within Δλ/2 of the range ends (copied from the interior) | the derivative spectrum now covers only the computable range; measuring outside it is an error |
| 3 | Regression | equal calibration concentrations returned NaN/∞ slope and LOD **silently** | clear error |
| 4 | SVR | results **depended on the concentration unit** (ng/mL gave nonsense vs µg/mL) because ε and C were absolute | targets standardised; ε is relative |
| 5 | ANN | predictions varied ≈ 0.3 % with sample order / units (single network, loose convergence) | ensemble of 5 networks, tight tolerance |
| 6 | Q limit | the 95 % Q limit flagged **0 %** of normal new samples (too lenient to detect real problems) | new limit (noise + LOO subspace term), validated at 4–5 % / 1 % |
| 7 | VIP | PLS1 VIP combined by averaging broke mean(VIP²) = 1 | root-mean-square combination |
| 8 | Pipelines | some parameter combinations returned **all-NaN spectra silently** | every step checks its output; clear message |
| 9 | Resample | end before start raised an internal `IndexError` | clear message |
| 10 | Spectrum | NaN values were accepted into spectra | refused at construction |
| 11 | Import | a repeated wavelength (overlap point) made the whole file unreadable | duplicates averaged |
| 12 | Import | an empty extra Excel sheet renamed all spectra "Sheet / …" | only sheets with data count |
| 13 | Statistics | empty / one-value inputs produced NaN tables | minimum-size checks with messages |
| 14 | Baseline / SNV | equal baseline wavelengths or a constant spectrum divided by zero | clear messages |
| 15 | UI | the spectral finder crashed on a project without spectra | message instead |
| 16 | Ratio subtraction | [(mixture ÷ Y′) − k] × Y′ was **wrong wherever Y′ ≈ 0** (ratio interpolated), e.g. 0.30 AU error at 229 nm in successive ratio subtraction; also constant center | computed as mixture − k·Y′ (same algebra, exact everywhere) |
| 17 | F-test | critical value was two-tailed (7.15 for 6/6) while the literature uses the one-tailed table value (5.05) | one-tailed by default, both reported, selectable |
| 18 | Divide | an exact-zero divisor point raised a divide-by-zero warning with threshold 0 | treated as below threshold |
| 19 | ALS baseline | asymmetry p = 0 or 1 gave a singular matrix | refused with a message; ties keep a weight |
| 20 | Robustness | derivative orders were varied by ±1 | not varied |
| 21 | Savitzky–Golay | polynomial orders up to the window length were accepted; order 18 gave an ill-conditioned fit (NumPy RankWarning) | orders above 10 refused (2–4 is usual) |
| 22 | Amplitude centering | the *unified regression* pooled the divisor compound with the two compounds sharing the isoabsorptive point → RIDSS recoveries ≈ 124 % (found while generating the method guide; the old test did not exercise it) | the divisor compound keeps its own line; test now checks recoveries |
| 23 | Optimizer | amplitude-centering candidates were generated only in the largest region where the divisor is usable; a divisor with two bands lost the region where the other compounds absorb, and plateau configurations vanished | every usable region is searched, λc must give each compound a usable ratio amplitude, two partner λ are screened per λc; ternary test errors fell from 3.5 / 6.6 / 3.0 % to 1.6 / 2.6 / 0.4 % |
| 24 | Concentration edits | correcting a sample's concentrations did not reach its processed (derivative, ratio…) versions, so methods calibrated on processed spectra and their recoveries still used the old values | the correction is copied to every processed version in the same audited transaction, and the saved methods/results that use the spectra are listed so they can be refitted / recalculated |

## Known limitations (documented, not bugs)

* **ANN** is reproducible only to its training precision (≈ 0.5 %); the invariance
  tests use 1 % for it.
* **Q limit with purely white noise** is conservative (flags < 1 % instead of 5 %).
  Real spectra always contain correlated (baseline) noise, for which it is
  calibrated.
* The method optimizer assumes additivity of absorbances (Beer–Lambert, no
  interactions); verify with laboratory mixtures.
