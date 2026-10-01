# Methods from the literature: coverage and validation

The lab supplied nine papers on resolving binary, ternary and quaternary
mixtures. This page lists every method they use, says where it is in
Spectro, and records what the validation tests (`tests/test_literature.py`)
found. The PDFs themselves are not in the repository.

| Code | Paper |
|---|---|
| AM14 | Abdelrahman et al., *Anal. Methods* 6 (2014) 509 — 10.1039/c3ay41564c — ambroxol / guaifenesin / theophylline |
| SA14 | Abdelrahman et al., *Spectrochim. Acta A* 124 (2014) 389 — 10.1016/j.saa.2014.01.020 — orphenadrine / caffeine / aspirin |
| LO15 | Lotfy et al., *Spectrochim. Acta A* 136 (2015) 937 — 10.1016/j.saa.2014.09.117 — lidocaine / fluocortolone / cinchocaine |
| AW17 | Abdelwahab & Mohamed, *Chem. Pharm. Bull.* 65 (2017) 558 — "Three new methods for resolving ternary mixture…" — metronidazole / diloxanide / mebeverine |
| SL17 | Saleh et al., *IJPPS* 9 (2017) 43 — 10.22159/ijpps.2017v9i5.16960 — amlodipine / valsartan / HCT (advanced amplitude centering) |
| EI18 | Eissa & Abou Al Alamein, *Spectrochim. Acta A* 193 (2018) 365 — 10.1016/j.saa.2017.12.050 — sacubitril / valsartan |
| EM18 | Emam et al., *Spectrochim. Acta A* (2018) — 10.1016/j.saa.2017.11.034 — furosemide / spironolactone / canrenone |
| ZA20 | Zaghary et al., *J. Anal. Chem.* 75 (2020) 742 — 10.1134/S1061934820060180 — canagliflozin / metformin |
| FA21 | Fahmy et al., *Spectrochim. Acta A* 261 (2021) 119999 — 10.1016/j.saa.2021.119999 — tretinoin / oxybenzone / hydroquinone / hydrocortisone acetate |

## Method coverage

✅ = already in Spectro before this review · 🆕 = **added in this review** (was missing) · 🔧 = was present, improved

| Method (papers) | Where in Spectro | Status |
|---|---|---|
| Zero order at λmax, D1 zero-crossing, D1 peak-to-peak, Δλ = 4 nm, scaling 10 (all) | *Methods → Univariate*, templates *Direct*, *Zero-crossing derivative*, *Derivative peak-to-peak* (Derivative step has Δλ and scaling factor) | ✅ |
| Dual wavelength (EI18), induced dual wavelength (EI18) | templates *Dual wavelength*, *Induced dual wavelength* | ✅ |
| D1 dual wavelength, D1 DWL (LO15) | template *Dual wavelength in derivative mode* | ✅ |
| Advanced absorbance subtraction (EI18) | *Binary special methods → Advanced absorbance subtraction* | ✅ |
| Ratio difference, derivative ratio DR1, mean centering of ratio spectra (SA14, LO15, EI18) | templates *Ratio difference*, *Derivative ratio (DD1)*, *Mean centering of ratio spectra* | ✅ |
| Ratio methods with a **normalized (unit-concentration) divisor** (EI18, ZA20, LO15, FA21) | *Process → Normalize → mode "concentration"* turns a pure standard into its unit spectrum; template *Ratio difference with normalized divisor* | 🆕 |
| Double divisor ratio derivative (SA14) | template *Double divisor ratio derivative (ternary)* | ✅ |
| Area under curve, binary and **ternary (3 × 3, Cramer's rule)** (AW17) | *Methods → Equations* with area signals (any number of compounds) | ✅ |
| AUC of derivative ratio spectra (SA14 method B) | *Methods → Equations*: processing (divide + D1) then area signals | ✅ |
| **Ternary mean centering of ratio spectra** — second division by the mean-centred ratio of the other component (SA14 method C) | new step *Divide by mean-centred ratio (ternary MCR)*; template *Mean centering of ratio spectra (ternary)* | 🆕 |
| Ratio subtraction, constant multiplication, spectrum subtraction, SS-CM, CM-SS (ZA20, LO15, FA21) | steps *Ratio subtraction*, *Constant multiplication*, *Spectrum subtraction*; template *SS-CM* | ✅ 🔧 |
| Extended ratio subtraction (LO15, EM18) | step *Extended ratio subtraction* | ✅ |
| **Successive ratio subtraction, SRS** (EM18, LO15) | template *Successive ratio subtraction (ternary SRS, Z)* (two ratio-subtraction steps with different divisors) | 🆕 template |
| **Derivative subtraction (DS) and DS-CM** (LO15, ZA20) | new option *Differentiate divisor first* on ratio subtraction / constant multiplication; templates *DS* and *DS-CM* | 🆕 |
| **D1 derivative ratio, D1 DR** — D1 of [D1 mixture ÷ D1 divisor] (LO15) | same option on *Divide*; template *Derivative ratio of D1 spectra (D1 DR)* | 🆕 |
| Constant value, plateau (LO15, FA21) | template *Ratio plateau / constant value* | ✅ |
| **Concentration value** — plateau with unit divisor *is* the concentration, no regression (FA21) | checkbox *Concentration value* in Univariate; template *Concentration value* | 🆕 |
| Constant center, CC-SS (LO15) | step *Constant center* | ✅ 🔧 |
| **Derivative transformation, DT and DT-SS** (FA21) | *Factorized spectrum / derivative transformation* and *Spectrum subtraction* now take a plateau range (k = mean of D(mixture)/D(pure) over it); template *Derivative transformation* | 🆕 |
| Amplitude modulation with normalized divisor and unified regression (ZA20, LO15) | *Binary special methods → Amplitude modulation* (unit divisor, total regression) | ✅ |
| **Advanced amplitude centering, AAC** — partial and complete overlap (SL17) | *Methods → Progressive resolution → Amplitude centering* | 🆕 |
| **Modified amplitude center method, MACM** (AW17) | same (complete-overlap form) | 🆕 |
| **Ratio difference–isoabsorptive, RIDSS** (AM14) | same (plateau + ratio difference + subtraction at the isoabsorptive point, unified regression) | 🆕 |
| **Constant value via amplitude difference, CV-AD** (LO15) | same (binary case, divisor compound by subtraction) | 🆕 |
| **Modified absorption factor method, MAFM** (AW17) | *Methods → Progressive resolution → Absorption factor (successive)* | 🆕 |
| **Sample enrichment** by spiking or spectrum addition (ZA20, LO15, FA21, SL17) | spectrum addition: *Add spectrum* step; spiking: *Enrichment added* in Univariate (Found − added column) | 🆕 (spiking) ✅ (addition) |
| Standard addition (all) | *Standard addition…* in Univariate | ✅ |
| PLS with 5-level calibration design, LOO CV (SL17) | *Methods → Chemometrics*, *Tools → Calibration design* | ✅ |
| Student's t, F, one-way ANOVA vs official / reported method (all) | *Tools → Validation statistics*: raw data (*Compare methods*, *Precision*) and **from published mean, SD, n** (*Compare (mean, SD, n)*) | 🆕 (summary values) 🔧 (F convention) |

**Nothing used in these nine papers is now missing.** The progressive
methods (amplitude centering, absorption factor) are saved like any other
method (*Save method* in the dialog; *Saved methods → Apply*, `.spmodel`
export with the divisor embedded), and the method optimizer builds, screens
and verifies them automatically (families *Amplitude centering* and
*Successive absorption factor*).

## What the validation tests check

1. **Published statistics are recomputed** from each paper's summary values
   with the app's functions. 88 published t and F values in AM14, SA14, EI18,
   FA21 and SL17 reproduce to within rounding (t ± 0.01–0.03, F ± 0.002–0.01),
   as do EI18's ANOVA and the means and SDs of the individual recoveries in
   EI18 and ZA20. So Spectro follows the same conventions as the
   literature: pooled-variance t, n − 1 standard deviations, F with the
   larger variance on top compared with the **one-tailed** tabulated value
   (5.05 for 6 vs 6).
2. **Every method is run on simulated spectra** built to meet the paper's
   conditions (extended component, isoabsorptive point, plateau, …), with
   wavelengths chosen by the app's finder tools as an analyst would.
   Noise-free mixtures are recovered to ± 0.1 % (± 0.3 % for the
   four-step quaternary chain) and noisy ones within 97–103 %. The resolved spectra are checked against the pure
   spectra (spectral profile).

## Problems found

### In Spectro (fixed)

| # | Problem | Fix |
|---|---|---|
| L1 | **Ratio subtraction, SRS and constant center lost the spectrum wherever the divisor was ≈ 0.** Computed as [(mixture ÷ Y′) − k] × Y′, the ratio is undefined where Y′ ≈ 0 and was interpolated, so the recovered spectrum came out wrong (the SRS test showed a 0.30 AU error at 229 nm) | computed as mixture − k·Y′. The algebra is identical, but the result is now exact at every wavelength |
| L2 | The F-test compared with the **two-tailed** critical value F(0.975) = 7.15 for 6 vs 6, while every paper (and most pharmacopoeial comparisons) use the one-tailed table value 5.05. Spectro was more lenient than the papers | the default is now one-tailed. The two-tailed value is always shown too, and either can be chosen in the dialog |
| L3 | An exact-zero divisor point raised a divide-by-zero warning when the threshold was 0 | such points are treated as below the threshold |
| L4 | ALS baseline with asymmetry p = 0 or 1 gave a singular system | p must lie strictly between 0 and 1; ties keep a weight |
| L5 | The robustness study varied derivative orders by ±1 | derivative orders are not varied |

### In the papers (worth knowing before reusing their numbers)

| Paper | Item | Printed | From the paper's own data | Effect |
|---|---|---|---|---|
| ZA20 | Table 3, MET by D1 vs reported: F | "less than tabulated" | F = 4.45 / 0.66 = **6.79 > 6.388** (n = 5) and > 5.05 (n = 6, the six mixtures actually listed) | **significant difference in precision** that the paper does not report. AM MET (F = 5.22) also exceeds 5.05 if n = 6 |
| ZA20 | Table 3 | n = 5, t(crit) = 2.571 | Table 2 lists 6 mixtures; for n = 5 + 5 the critical t is 2.306 | inconsistent n and table values |
| LO15 | Table 5 | t(crit) = 2.57 | n = 6 + 6 → 10 df → **2.228** | the critical t is wrong. Several t/F values could not be reproduced from the table's means and SDs as extracted (e.g. FCP D1: t = 3.39). Check the printed PDF before reusing them |
| FA21 | Table 4, HC ANOVA | SS between 0.61, F 1.52 | SS 1.455, **F 3.64** (F crit 3.68) | still not significant, but only just |
| FA21 | Table 3, EX constant value | t = 2.25 (crit 2.23), "no difference" | t = 2.22 < 2.228 | borderline; the printed value contradicts the paper's own conclusion |
| EM18 | Table 4, furosemide standard addition | 7.86 of 8.00 → 99.50 % | **98.25 %**. With it, the printed mean ± SD (100.00 ± 1.561) is exact | typo in one row |
| AM14 | Table 3, Farcosolvin GUF | t = 1.250 | t = 1.357 | conclusion unchanged |
| SL17 | Table 6 | F crit 6.094 / 4.120 | F(0.95; 6, 4) = **6.163**, F(0.95; 4, 6) = **4.534**; HCT (severe) t = 0.822 (printed 0.764) | conclusions unchanged |
| AW17 | AUC equations | A(240–250) MEB coefficient 0.1352 in the text, 0.1421 in Table 1 (also used for A(280–290)) | — | one of the two is a typo |
