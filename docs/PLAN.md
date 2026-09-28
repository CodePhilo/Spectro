# Spectro — Project Plan (Draft for Review)

A Windows desktop application for organizing UV/Vis spectral trials in a
pharmaceutical analytical chemistry lab: import spectra, process them with the
full range of univariate and multivariate methods used for **simultaneous
determination of compounds in mixtures**, validate methods, and keep a
**complete, tamper‑evident log of every operation** performed on every spectrum.

> Status: **Approved with changes (review round 1).**
> Decisions from review:
> - **No user roles / login / e‑signatures.** The audit trail still records the
>   operating‑system user and workstation for every entry.
> - **Instrument‑agnostic import.** No parser is tied to a specific vendor or
>   model: a generic text/CSV parser (any delimiter, decimal comma, header
>   lines, column/row/XY‑pair layouts), Excel (.xlsx/.xlsm/.xls), and the
>   open exchange formats JCAMP‑DX and Galactic SPC that most instrument
>   software can export.

---

## 1. Goals and non‑goals

### Goals
1. Organize work as **Projects → Studies/Trials → Samples → Spectra**, with
   metadata (analyst, instrument, solvent, cell path length, slit width, scan
   speed, date, batch/lot, dilution factor).
2. Import spectra from common UV/Vis instruments and generic formats.
3. Provide every mainstream spectrophotometric method for binary, ternary and
   quaternary mixtures (zero‑order, derivative, ratio‑spectra family,
   chemometrics).
4. **Non‑destructive processing**: raw data is never altered; every result is
   produced by a recorded pipeline that can be replayed.
5. **Audit trail** of every action (import, edit, process, delete, sign) —
   who, what, when, why, before/after values.
6. Method validation per **ICH Q2(R2)** and statistical comparison with a
   reference method.
7. Report generation (PDF / Excel) with the processing history attached.

### Non‑goals (v1)
- Direct instrument control (driving the spectrophotometer). **[DECIDE]** —
  see §11.
- Multi‑site server/database. v1 is single workstation (optionally a shared
  network folder). See §11.
- Non‑UV techniques (IR, Raman, fluorescence) — the architecture will allow
  them later, but they are not in scope.

---

## 2. Recommended technology stack **[REC]**

| Concern | Choice | Why |
|---|---|---|
| Language | **Python 3.12** | NumPy/SciPy/scikit‑learn ecosystem is the standard for spectral math; easy to validate numerically. |
| GUI | **PySide6 (Qt 6)** | Native Windows look, docking panels, tables, LGPL license. |
| Plotting | **pyqtgraph** (interactive) + **matplotlib** (report figures) | pyqtgraph handles hundreds of overlaid spectra smoothly; matplotlib gives publication‑quality exports. |
| Numerics | NumPy, SciPy, scikit‑learn, (optional) `pymcr` for MCR‑ALS | Proven, well‑tested implementations. |
| Storage | **SQLite** (one `.spectro` project file) + optional raw‑file vault | Single portable file, ACID transactions, easy backup. |
| Reports | ReportLab (PDF), openpyxl (Excel) | |
| Packaging | PyInstaller → signed Windows installer (Inno Setup) | Runs on lab PCs without a Python install. |
| Testing | pytest + reference datasets with known answers | Needed for software validation (CSV/GAMP 5). |

Alternative considered: C#/.NET WPF. Better native feel, but the chemometrics
libraries would have to be re‑implemented or bridged — higher risk for
numerical correctness. Python is recommended.

---

## 3. Architecture

```
┌──────────────────────────── UI (PySide6) ─────────────────────────────┐
│ Project tree │ Spectrum viewer (overlay/derivative/ratio) │ Properties │
│ Method wizards │ Calibration & prediction tables │ Audit log viewer    │
└──────────────┬────────────────────────────────────────────────────────┘
               │ commands (every user action is a Command object)
┌──────────────▼──────────────┐     ┌──────────────────────────────────┐
│ Application services        │────▶│ Audit service (append‑only,      │
│ (projects, pipelines,       │     │ hash‑chained, reason‑for‑change) │
│  methods, validation, users)│     └──────────────────────────────────┘
└──────────────┬──────────────┘
┌──────────────▼──────────────┐     ┌──────────────────────────────────┐
│ Core (pure Python, no UI)   │     │ Persistence (SQLite repository)  │
│ io/ preprocessing/          │     │ raw spectra immutable (BLOB +    │
│ univariate/ ratio/          │     │ SHA‑256), derived data cached,   │
│ chemometrics/ validation/   │     │ pipelines stored as JSON         │
└─────────────────────────────┘     └──────────────────────────────────┘
```

Key design rules:
- **Core has zero UI dependencies** → every algorithm is unit‑testable and can
  also be used from a script/Jupyter for cross‑checking.
- **Every operation is a `Command`** (name, parameters, inputs, outputs,
  version). The command bus executes it *and* writes the audit record in the
  same DB transaction — it is impossible to process without logging.
- **Processing is a pipeline (DAG)**: `raw → smooth(SG, 11, 2) → D1 → ratio(÷ divisor X) → D1 → read @ 262.4 nm`.
  Changing a parameter creates a new pipeline version; old results remain
  reproducible.
- **Plugin registry** for methods: each method declares its parameters, input
  requirements (e.g. "needs pure component spectra") and outputs, so the UI
  wizard is generated automatically and new methods are cheap to add.

Proposed repo layout:
```
spectro/
  core/        io/, preprocessing/, univariate/, ratio/, chemometrics/, validation/
  app/         commands, services, audit, users/auth
  storage/     sqlite schema, migrations, repositories
  ui/          main window, viewers, wizards, dialogs
  reports/
tests/         unit + reference datasets (known concentrations)
docs/
```

---

## 4. Data model (main entities)

- **Project** – name, description, owner, created.
- **Compound** – name, MW, λmax, reference absorptivity (a / ε), CAS.
- **Trial / Study** – objective (e.g. "Binary mixture PAR + CAF"), compounds,
  solvent, instrument, cell (cm), status (draft → in review → approved).
- **Sample** – type (pure standard, lab‑prepared mixture, pharmaceutical
  formulation, blank, spiked/standard addition), known concentrations
  (µg/mL), dilution factor, label claim, lot.
- **Spectrum** – λ array, A array, acquisition metadata, source file + SHA‑256,
  role (calibration / validation / unknown / divisor / blank).
- **Pipeline** – ordered processing steps + parameters (JSON), version.
- **MethodModel** – calibration result: equation / regression / PLS model,
  selected wavelengths, figures of merit.
- **Result** – predicted concentrations, % recovery, % of label claim.
- **AuditEntry** – see §6.
- **User / Role / Signature**.

---

## 5. Functional scope

### 5.1 Import / export
- CSV/TXT (λ, A columns; multi‑column), Excel.
- **JCAMP‑DX** (.jdx/.dx), **Thermo Galactic SPC** (.spc).
- **Any instrument model**: parsers are layout‑driven, not vendor‑driven.
  Auto‑detects delimiter, decimal separator, header/metadata lines, and
  layout (λ + N columns, XY column pairs, or spectra in rows).
- Batch import with filename → metadata mapping rules (e.g.
  `PAR_10_CAF_5_rep2.csv` → PAR 10 µg/mL, CAF 5 µg/mL, replicate 2).
- Resampling/interpolation to a common λ grid, with a check for mismatched
  ranges/intervals.
- Export processed spectra (CSV, JCAMP‑DX) and results (Excel).

### 5.2 Viewing
- Overlay, stacked, and difference views; zoom/pan; crosshair readout of λ/A.
- Live preview of any pipeline step (e.g. drag the derivative order or
  Δλ and see the curve update).
- Auto‑detection of **zero‑crossing points**, **isoabsorptive points**,
  maxima/minima, plateau regions of ratio spectra.
- Overlay of pure‑component spectra with the mixture for method design.

### 5.3 Pre‑processing
- Blank/solvent subtraction; baseline correction (offset, linear, polynomial,
  ALS/asymmetric least squares).
- Smoothing: Savitzky–Golay, moving average, Whittaker.
- Wavelength range cropping, interpolation/resampling.
- Normalization (max, area, vector), SNV, MSC, mean‑centering, autoscaling.
- Spectral arithmetic: add, subtract, multiply/divide by constant or by
  another spectrum.

### 5.4 Univariate methods — zero order
- Direct measurement at λmax (single component / non‑overlapping).
- **Simultaneous equations (Vierordt's)** — binary and multi‑component
  (least‑squares with n wavelengths).
- **Absorbance ratio (Q‑analysis)** using isoabsorptive point.
- **Dual wavelength** (and **induced dual wavelength**).
- **Absorbance subtraction** / **advanced absorbance subtraction**.
- **Amplitude factor**, **amplitude modulation**, **induced amplitude
  modulation**.
- **Area under the curve (AUC)**.
- **Factorized zero‑order** / factorized response methods.
- **H‑point standard addition method (HPSAM)**.
- **Bivariate calibration**.

### 5.5 Derivative methods
- Derivatives order 1–4 (Savitzky–Golay or finite difference with Δλ and
  scaling factor), peak‑to‑baseline and peak‑to‑peak amplitudes.
- **Zero‑crossing technique** (D1, D2, …).
- **Dual wavelength in derivative mode**.

### 5.6 Ratio‑spectra family
- **Derivative ratio (DD1 / DR)**.
- **Successive derivative ratio** (ternary mixtures).
- **Double divisor ratio spectra derivative** (ternary).
- **Ratio difference (RD)**.
- **Ratio subtraction** and **extended ratio subtraction**.
- **Constant multiplication** and **constant center** (± spectrum subtraction).
- **Mean centering of ratio spectra (MCR‑spectra)**.
- **Spectrum subtraction** / **successive spectrum subtraction**.
- **Dual amplitude difference**, **amplitude center**.
- Divisor selection helper: try several divisor concentrations and show the
  effect on linearity/noise (common practical pain point).

### 5.7 Multivariate / chemometrics
- **CLS** (classical least squares), **ILS/MLR**.
- **PCR**, **PLS‑1 / PLS‑2**.
- **MCR‑ALS** (with non‑negativity, closure, unimodality constraints).
- **ANN** (MLP regression) — optional **[DECIDE]**.
- Variable selection: interval PLS (iPLS), GA‑PLS, VIP / selectivity ratio.
- Cross‑validation (LOO, k‑fold, venetian blinds), RMSECV / RMSEP, choice of
  latent variables with plots.
- **Calibration‑set design**: multilevel multifactor design (Brereton) — the
  app generates the list of mixtures to prepare.
- Scores/loadings/residual plots, outlier detection (Hotelling T², Q).
- Model export/import for routine use.

### 5.8 Validation & statistics (ICH Q2(R2))
- Linearity: regression, r, r², slope/intercept with SD and CI, residual plot,
  lack‑of‑fit test.
- LOD / LOQ (3.3σ/S, 10σ/S; σ from intercept SD or residual SD).
- Accuracy (% recovery ± SD), standard addition.
- Precision: repeatability and intermediate precision (%RSD).
- Specificity: lab‑prepared mixtures at different ratios.
- Robustness (optional, parameter variation table).
- Application to dosage forms: % label claim.
- Comparison with reference/official method: Student's t, F‑test, one‑way
  ANOVA, and interval hypothesis tests.
- (Optional) greenness scores (Eco‑Scale, GAPI, AGREE) commonly required by
  journals **[DECIDE]**.

### 5.9 Reporting
- Trial report (PDF): metadata, spectra figures, pipeline steps, calibration
  tables, validation tables, results, **audit trail excerpt**, signatures.
- Excel export of all tables.
- Figure export (PNG/SVG/TIFF, 300–600 dpi) for publications.

---

## 6. Operation logging / audit trail

Designed to meet the intent of **FDA 21 CFR Part 11**, **EU GMP Annex 11** and
**ALCOA+** data‑integrity principles.

Each `AuditEntry` records:
- sequence number, UTC timestamp (+ local time zone), user id, workstation,
  app version;
- action type (IMPORT, CREATE, EDIT, PROCESS, DELETE/ARCHIVE, SIGN, LOGIN,
  EXPORT, CONFIG_CHANGE…);
- target entity and id;
- operation name + full parameters (JSON);
- before/after values for edits;
- **reason for change** (mandatory for edits/deletes once a trial is out of
  draft);
- hash of input data and output data;
- `prev_hash` + `entry_hash` (SHA‑256 hash chain) → any tampering with the log
  is detectable by a built‑in **integrity check**.

Rules:
- Audit table is **append‑only** (no UPDATE/DELETE path in the code; DB
  triggers reject them).
- Nothing is physically deleted — records are archived/voided with a reason.
- Audit viewer with filters (user, date, trial, action) and export to PDF/CSV.
- Per‑spectrum "history" tab showing its full lineage (raw file → every
  processing step → results that used it).
- "Replay" button: re‑executes a pipeline from raw data and confirms the
  result matches the stored one (reproducibility check).

---

## 7. Users & security — *removed from scope (review round 1)*

No roles, login or e‑signatures. Every audit entry records the OS user name
and workstation automatically.

## 8. Quality & verification of the software

- Unit tests for every algorithm against **published datasets / hand‑computed
  examples** (known concentrations must be recovered within tolerance).
- Cross‑check PLS/PCR against scikit‑learn and, where possible, published
  reference results.
- Golden‑file tests for importers (sample files from each instrument).
- Test that every UI command produces exactly one audit entry.
- Deliverables for lab CSV (computer system validation): requirements
  spec (URS), traceability matrix, test report — **[DECIDE]** whether you
  need these formally (GAMP 5 Category 4/5).

---

## 9. Milestones

| # | Milestone | Content | Exit criteria |
|---|---|---|---|
| M0 | Foundation | Repo scaffold, CI, core data model, SQLite schema, audit service with hash chain, command bus | Every command writes an audit entry; integrity check passes/fails correctly |
| M1 | Import & view | CSV/TXT/Excel/JCAMP/SPC importers, project tree, spectrum viewer, metadata editor | Import 100 spectra and overlay smoothly |
| M2 | Pre‑processing & pipelines | Smoothing, baseline, derivatives, arithmetic, pipeline editor with live preview, lineage view | Pipeline replay reproduces results bit‑for‑bit |
| M3 | Univariate methods | Zero‑order + derivative + ratio‑spectra family (§5.4–5.6), zero‑crossing / isoabsorptive detection | Reference datasets recovered within ±2 % |
| M4 | Validation & reports | ICH Q2 module, statistics, PDF/Excel reports | Report for a full binary‑mixture trial |
| M5 | Chemometrics | CLS, PCR, PLS, MCR‑ALS, CV, calibration design, variable selection | Matches scikit‑learn / published results |
| M7 | Packaging | Windows installer, user manual, sample project | Installs on a clean Windows 10/11 PC |

Suggested order puts the audit trail **first** (M0) — retro‑fitting it later
is the most common failure in lab software.

---

## 10. Risks

| Risk | Mitigation |
|---|---|
| Vendor file formats are undocumented | Get sample exports early; fall back to CSV export from the instrument software. |
| Numerical differences vs. published papers (derivative scaling, Δλ) | Expose all parameters; document conventions; test against papers. |
| Regulatory expectations vary | Build Part 11‑style controls from the start; confirm scope with QA. |
| Scope creep (many ratio methods) | Plugin architecture; implement methods in priority order you provide. |

---

## 11. Open questions **[DECIDE]**

1. Which spectrophotometer(s) and software (Shimadzu UV‑1800/UVProbe, Agilent
   Cary, PerkinElmer Lambda, Jasco V‑series…)? Can you share a few sample
   export files?
2. Single PC, or several PCs sharing projects (network folder vs. a small
   server DB such as PostgreSQL)?
3. Is the app used in a **GMP/QC** setting (needs Part 11 controls and CSV
   documentation) or **R&D/academic** (lighter controls acceptable)?
4. Which methods are highest priority? (e.g. your current trials are
   binary/ternary, which methods do you publish most?)
5. Is ANN / advanced chemometrics needed in v1, or can it wait for v2?
6. Report language/branding (lab logo, SOP numbers, template)?
7. Greenness assessment tools (AGREE/GAPI/Eco‑Scale) — needed?
8. Is Python + PySide6 acceptable, or is there an IT requirement for .NET?

---

## 12. Review checklist

Please tick what you accept, strike what you reject, and add comments.

**Architecture & stack**
- [ ] Python 3.12 + PySide6 + pyqtgraph + NumPy/SciPy/scikit‑learn
- [ ] Single‑file SQLite project database
- [ ] Core library separated from UI (scriptable, testable)
- [ ] Windows installer via PyInstaller + Inno Setup

**Data integrity (recommended)**
- [ ] Raw spectra immutable, stored with SHA‑256 of the source file
- [ ] Non‑destructive, versioned processing pipelines with replay
- [ ] Append‑only, hash‑chained audit trail with integrity check
- [ ] Mandatory reason‑for‑change after a trial leaves Draft
- [ ] No hard deletes (archive/void only)
- [x] ~~User roles, login, e‑signatures~~ (removed)
- [ ] Audit trail included in PDF reports

**Scientific scope**
- [ ] Pre‑processing set (§5.3)
- [ ] Zero‑order methods (§5.4)
- [ ] Derivative methods (§5.5)
- [ ] Ratio‑spectra family (§5.6)
- [ ] Chemometrics CLS/PCR/PLS/MCR‑ALS (§5.7)
- [ ] ANN (optional)
- [ ] Calibration‑set experimental design generator
- [ ] ICH Q2(R2) validation module + statistical comparison (§5.8)
- [ ] Greenness metrics (optional)

**Process**
- [ ] Milestone order M0 → M7 (audit trail first)
- [ ] Algorithm verification against reference datasets
- [ ] Formal CSV documents (URS, traceability matrix, test report)

**Your answers to §11 open questions**
- [ ] Instruments / sample files provided
- [ ] Deployment (single PC / shared)
- [ ] GMP vs R&D setting
- [ ] Method priority list
