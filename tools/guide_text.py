"""Beginner-level text of the method guide (used by tools/method_guide.py).

For every method: ``plain`` (the idea in plain words), ``need`` (which
spectra to measure) and ``check`` (how to know it worked / common mistakes).
"""

BEGINNER: dict[str, dict[str, str]] = {
    # ------------------------------------------------------------- zero order
    "direct": dict(
        plain="If only one compound absorbs light at a certain wavelength, the absorbance there "
              "tells you its concentration directly — exactly as for a pure solution.",
        need="5–6 standards of the compound (e.g. 4–20 µg/mL) and the mixtures.",
        check="Overlay the pure spectra of all compounds first: at the chosen λ every other "
              "compound must read ≈ 0. If it does not, the result is biased high."),
    "dw": dict(
        plain="Pick two wavelengths where the interfering compound absorbs the SAME amount. "
              "Subtracting the two readings removes the interferent completely, and what is "
              "left belongs to your analyte only.",
        need="Standards of the analyte (for the calibration), one standard of the interferent "
             "(to find the two wavelengths), the mixtures.",
        check="Check with a pure interferent standard: its ΔA must be ≈ 0. Pick a pair where "
              "the analyte's two absorbances differ a lot (large ΔA = better sensitivity)."),
    "idw": dict(
        plain="Like dual wavelength, but the two readings of the interferent need not be equal: "
              "the second reading is multiplied by a factor F (= ratio of the interferent's "
              "two absorbances) so that the interferent still cancels.",
        need="Analyte standards, one interferent standard (to compute F), mixtures.",
        check="F must come from a pure interferent spectrum measured the same way. A pure "
              "interferent standard processed like a sample must give ≈ 0."),
    "auc1": dict(
        plain="Instead of one absorbance value, use the area under the spectrum between two "
              "wavelengths — it is proportional to concentration and less sensitive to noise.",
        need="Standards of the compound, mixtures; a band where only this compound absorbs.",
        check="Use the same integration limits for standards and samples; check that other "
              "compounds have no absorbance inside the range."),
    "q": dict(
        plain="At the isoabsorptive point both compounds absorb equally, so the absorbance "
              "there measures their SUM. A second wavelength tells how that sum is split.",
        need="Standards of both compounds (to find the isoabsorptive point and the "
             "absorptivities), mixtures.",
        check="The isoabsorptive point must be found from spectra of EQUAL concentrations "
              "(the finder divides by concentration for you)."),
    "as": dict(
        plain="One compound (X) absorbs alone at a long wavelength. From X's reading there and "
              "X's own spectrum shape, Spectro works out how much X absorbs at the "
              "isoabsorptive point; the rest of the absorbance there belongs to the other one.",
        need="Standards of X (for the amplitude factor), standards of both (for the iso "
             "calibration), mixtures.",
        check="X must be the extended compound (absorbing alone at λ2)."),
    "aas": dict(
        plain="Two wavelengths where the interferent Y reads the same (one of them the "
              "isoabsorptive point): their difference gives X; X's share at the iso point is "
              "then subtracted to give Y.",
        need="Standards of X, standards of Y, mixtures.",
        check="Use the finder button to get λ1 exactly; a wrong λ1 leaves part of Y in ΔA."),
    "mafm": dict(
        plain="Peel the mixture like an onion: the compound that absorbs alone at the longest "
              "wavelength is measured first; its contribution at the next wavelength is "
              "calculated from its own spectrum shape and removed; then the next compound, "
              "and so on.",
        need="Standards of every compound, mixtures.",
        check="Order matters: the first compound must truly absorb alone at its wavelength, "
              "the second must be overlapped only by the first, and so on. Errors carry "
              "forward to the later compounds."),
    # ------------------------------------------------------------- equations
    "vierordt": dict(
        plain="Measure at two (or more) wavelengths. Each reading is the sum of both compounds' "
              "contributions, so you get two equations with two unknowns and solve them.",
        need="Standards of each compound (zero concentration for the other), mixtures.",
        check="The condition number shown after fitting should be small (< 10 is very good; "
              "> 100 means the wavelengths do not distinguish the compounds)."),
    "bivariate": dict(
        plain="The same idea as simultaneous equations, but Spectro chooses the two "
              "wavelengths that separate the compounds best (Kaiser method).",
        need="Standards of both compounds, mixtures.",
        check="Look at the top pairs listed: choose one inside a region with good absorbance, "
              "not at the noisy ends of the spectrum."),
    "auc_eq": dict(
        plain="Simultaneous equations with areas instead of single absorbances: one area range "
              "per compound, solved together.",
        need="Standards of every compound (with zeros for the others), mixtures.",
        check="Each range should be dominated by a different compound."),
    # ------------------------------------------------------------- derivative
    "zc": dict(
        plain="Differentiating a spectrum gives its slope. The slope of the interferent is "
              "zero at the top of its band (zero-crossing): read your analyte's derivative "
              "exactly there and the interferent contributes nothing.",
        need="Analyte standards, one interferent standard (to find its zero-crossings), "
             "mixtures.",
        check="Keep Δλ and the scaling factor identical for standards and samples. Avoid "
              "zero-crossings where the analyte's derivative is small."),
    "p2p": dict(
        plain="Measure the height from the lowest to the highest point of a derivative band — "
              "a bigger signal than a single point and immune to baseline shifts.",
        need="Standards and mixtures; a derivative band free from the other compounds.",
        check="Make the range wide enough to contain both the maximum and the minimum."),
    "d1dwl": dict(
        plain="Dual wavelength applied to derivative spectra: two wavelengths where the "
              "interferent's derivative is equal.",
        need="Analyte standards, one interferent standard, mixtures.",
        check="Find the pair on the derivative of the interferent (process it first)."),
    # ------------------------------------------------------------- ratio derivative
    "dd1": dict(
        plain="Divide the mixture spectrum by the interferent's spectrum. The interferent "
              "becomes a flat line (a constant) and the derivative of a constant is zero, so "
              "the derivative of the ratio spectrum depends only on your analyte.",
        need="Analyte standards, a divisor (one standard of the interferent), mixtures.",
        check="Use the same divisor for standards and samples. Avoid wavelengths where the "
              "divisor absorbs almost nothing (the ratio becomes noisy there)."),
    "d1dr": dict(
        plain="First take the derivative of the mixture, then divide by the derivative of the "
              "interferent, then differentiate again: the interferent cancels.",
        need="Analyte standards, a divisor, mixtures.",
        check="Spikes appear where the divisor's derivative crosses zero — read the signal "
              "away from them."),
    "mcr": dict(
        plain="Divide by the interferent (it becomes a constant), then subtract the average of "
              "the ratio spectrum — a constant disappears when you subtract its average.",
        need="Analyte standards, a divisor, mixtures.",
        check="The mean-centring range must be the same for every spectrum."),
    "mcr3": dict(
        plain="Two divisions for three compounds: the first divisor turns Z into a constant "
              "(removed by mean centring), the second turns Y into a constant (removed again).",
        need="Standards of X, divisors of Y and Z, mixtures.",
        check="Use the same range in every mean-centring step."),
    "sdr": dict(
        plain="Ratio and derivative applied twice: first Y is removed (÷ Y′, derivative), then Z "
              "(÷ derivative of Z′/Y′, derivative).",
        need="Standards of X, divisors Y′ and Z′, mixtures; the second divisor is made once in "
             "Process.",
        check="Read the signal away from spikes caused by zeros of the second divisor."),
    "dd": dict(
        plain="Divide by the sum of the other two standards at once, then differentiate.",
        need="Standards of X, divisors Y′ and Z′, mixtures.",
        check="Exact only when Y : Z in the sample equals Y′ : Z′; test mixtures of other ratios "
              "before trusting it."),
    # ------------------------------------------------------------- ratio amplitude
    "rd": dict(
        plain="Divide by the interferent (it becomes a constant). The difference between two "
              "points of the ratio spectrum removes the constant — no derivative needed.",
        need="Analyte standards, a divisor, mixtures.",
        check="Choose λ1 at a peak and λ2 at a trough of the analyte's ratio spectrum."),
    "dad": dict(
        plain="For three compounds: divide by Z (constant) and choose two points where Y's "
              "ratio spectrum is equal: their difference cancels both Y and Z.",
        need="Standards of X, divisors Z′ and Y′, mixtures.",
        check="The two wavelengths must give equal Y/Z′ amplitudes (finder)."),
    "cv": dict(
        plain="Where only the extended compound Y absorbs, mixture ÷ Y′ is a flat line whose "
              "height is the concentration ratio C_Y / C_Y′.",
        need="Standards of Y, the divisor Y′, mixtures.",
        check="The plateau must be flat in every mixture; if it slopes, another compound "
              "absorbs there."),
    "concval": dict(
        plain="Constant value with a divisor of exactly 1 µg/mL: the flat line's height IS "
              "the concentration — no calibration line needed.",
        need="A divisor normalised to 1 µg/mL (made in Spectro), mixtures.",
        check="Still run the calibration once: slope ≈ 1 and intercept ≈ 0 confirm it."),
    "am": dict(
        plain="With a 1 µg/mL divisor of Y, the ratio spectrum at the isoabsorptive point "
              "equals C_X + C_Y and its plateau equals C_Y; X = total − Y.",
        need="Standards of X and Y, a unit-concentration Y divisor, mixtures.",
        check="Needs an isoabsorptive point AND a plateau region of Y."),
    "iam": dict(
        plain="Amplitude modulation without an isoabsorptive point: the X contribution at any "
              "wavelength is converted with the ratio of absorptivities.",
        need="Standards of X and Y, mixtures.",
        check="Y must have a plateau region."),
    "cvad": dict(
        plain="Divide by Y′: X's ratio difference between two wavelengths tells how big X's "
              "share is at one of them; what remains there is Y.",
        need="Standards of X and Y, a divisor Y′, mixtures.",
        check="The amplitude-difference line (shown after calculating) should have r ≈ 1."),
    "aac": dict(
        plain="One divisor, one wavelength, three compounds: the plateau gives Z, an "
              "amplitude difference gives X, and Y is what is left.",
        need="Standards of X, Y and Z, a divisor Z′, mixtures.",
        check="The plateau must be in the region where only Z absorbs."),
    "macm": dict(
        plain="Like AAC but without a plateau: two different amplitude differences give X and "
              "Y, and Z is what is left at the common wavelength.",
        need="Standards of all three, one divisor, mixtures.",
        check="Let the method optimizer propose λc and the λ pairs: hand-picked pairs on steep "
              "slopes are sensitive to small wavelength errors."),
    "ridss": dict(
        plain="Z from the plateau, Y from a ratio difference, and X + Y together at their "
              "isoabsorptive point — so X = total − Y.",
        need="Standards of X, Y, Z, a divisor Z′, mixtures.",
        check="Tick 'unified regression' only when λc is the isoabsorptive point of X and Y."),
    # ------------------------------------------------------------- resolution
    "rs": dict(
        plain="Remove the extended compound Y from the mixture: its amount is read from the "
              "plateau and that much of Y′ is subtracted. What remains looks exactly like the "
              "spectrum of pure X, which you then read at its λmax.",
        need="Standards of X, a divisor Y′, mixtures.",
        check="The recovered spectra should overlay the pure X spectra (look at the plot)."),
    "cm": dict(
        plain="The plateau height × Y′ rebuilds Y's own spectrum in the mixture; read Y at its "
              "λmax.",
        need="Standards of Y, a divisor Y′, mixtures.",
        check="The rebuilt spectra are copies of Y′ scaled — check they match the mixture in "
              "the plateau region."),
    "ers": dict(
        plain="After ratio subtraction gives X, X is matched to its pure spectrum and "
              "subtracted from the mixture: Y's full spectrum remains.",
        need="Standards of Y, a divisor Y′, a pure X spectrum, mixtures.",
        check="Useful when Y must be read at a band overlapped by X."),
    "srs": dict(
        plain="Remove the compounds one by one: first the most extended (U) with its plateau, "
              "then the next (V) with its own plateau; W's pure spectrum remains.",
        need="Standards of W, divisors U′ and V′, mixtures.",
        check="Each plateau must be where only that compound remains."),
    "sss": dict(
        plain="Same as SRS but the first removal uses a single wavelength (spectrum "
              "subtraction) instead of a ratio plateau.",
        need="Standards of W, references U′ and V′, mixtures.",
        check="Use a wavelength (or range) where only U absorbs."),
    "ss": dict(
        plain="Scale the pure Y spectrum to match the mixture where only Y absorbs, then "
              "subtract it.",
        need="Standards of X, a reference Y′, mixtures.",
        check="A range (λ … to) averages noise better than a single wavelength."),
    "fzm": dict(
        plain="Use a derivative zero-crossing to find how much X is in the mixture, then scale "
              "X's pure spectrum by that amount: you get X's whole zero-order spectrum.",
        need="Standards of X, a pure X reference, mixtures.",
        check="Gives the same number as zero-crossing D1, plus the recovered spectrum."),
    "cc": dict(
        plain="In mixture ÷ Y′, X's ratio spectrum has a known shape; comparing two points "
              "reveals the hidden constant of Y without needing a plateau.",
        need="Standards of X, a divisor Y′, a pure X reference, mixtures.",
        check="λ1 and λ2 must give clearly different X ratio amplitudes."),
    "dt": dict(
        plain="In derivative spectra the interferent is flat (zero) in a region where Y still "
              "changes; the ratio there gives Y's amount, which rebuilds Y's zero-order "
              "spectrum.",
        need="Standards of Y, a reference Y′, mixtures.",
        check="Choose the region where the interferent's derivative is zero."),
    "ds": dict(
        plain="Ratio subtraction and constant multiplication done on derivative spectra — "
              "useful when the zero-order plateau is too small or noisy.",
        need="Standards, a divisor, mixtures.",
        check="The divisor is differentiated automatically ('Differentiate divisor first')."),
    # ------------------------------------------------------------- multivariate
    "chemo": dict(
        plain="Instead of choosing wavelengths, a statistical model learns from many mixtures "
              "of known composition how the whole spectrum depends on each concentration.",
        need="A training set of 15–25 mixtures (a calibration design), validation mixtures.",
        check="Cross-validate to choose the number of components; check RMSEP and recoveries "
              "of mixtures NOT used for training; look at the T²/Q outlier plot."),
    # ------------------------------------------------------------- standard addition
    "stdadd": dict(
        plain="Add known amounts of the pure analyte to the real sample and check you find "
              "exactly that much more — this proves excipients do not interfere.",
        need="The sample alone and the sample with 2–3 additions; a calibrated method.",
        check="Recoveries of the added amounts should be 98–102 %."),
    "hpsam": dict(
        plain="Plot the response vs amount added at two wavelengths where the interferent is "
              "equal: the two lines cross at a point whose position gives the analyte "
              "concentration, free of the interferent.",
        need="The sample with 3–4 standard additions of the analyte.",
        check="The two lines must have clearly different slopes."),
    "enrich": dict(
        plain="When one component is far below its linear range, add a known amount of it, "
              "measure, and subtract the amount added.",
        need="Standards, the sample spiked with a known amount of the minor component.",
        check="Report the content after subtracting the added amount."),
}

PRIMER = """
### 1. What these methods are for

A pharmaceutical mixture (e.g. a tablet with two or three drugs) gives one UV
spectrum that is the **sum** of the spectra of its components (Beer–Lambert law:
A = a·b·C for each compound, and absorbances add up). When the bands overlap, a
single wavelength cannot tell the drugs apart. Every method in this guide is a
way to separate that sum:

* **choose wavelengths** where the other compounds cancel (zero order, derivative);
* **divide by an interferent's spectrum** so it becomes a constant that can be
  removed (ratio spectra);
* **rebuild each component's own spectrum** (spectrum resolution);
* **solve equations** with several signals, or **train a model** on many mixtures
  (chemometrics).

In the method descriptions **X** is the compound you determine and **Y**, **Z**
the interfering ones; a prime (Y′) marks a standard spectrum used as divisor or
reference. The worked examples use the compounds of the guide systems
(A, B; X, Y, Z; U, V, W — see section 6).

### 2. Words you will meet

| Term | Meaning |
|---|---|
| Standard | Solution of ONE pure compound at a known concentration; 5–6 levels make a calibration line |
| Calibration line | Signal vs concentration of the standards; its slope and intercept convert a sample's signal into a concentration |
| Laboratory-prepared mixture | Mixture made by you from standards, so the true concentrations are known; used to test a method |
| Recovery % | Found ÷ taken × 100; 98–102 % is the usual acceptance range |
| RSD % | Relative standard deviation of repeated results; ≤ 2 % is usual |
| λmax | Wavelength of maximum absorbance of a band |
| Isoabsorptive (isosbestic) point | Wavelength where two compounds of EQUAL concentration have equal absorbance |
| Extended component | The compound whose spectrum continues alone at longer wavelengths |
| Divisor (Y′) | A standard spectrum of an interfering compound that the mixture is divided by |
| Ratio spectrum | Mixture ÷ divisor; the divisor compound becomes a flat line (constant) |
| Plateau | The flat part of a ratio spectrum where only the divisor compound absorbs |
| Unit (normalised) divisor | Divisor rescaled to 1 µg/mL (*Process → Normalize → concentration*) |
| Derivative D1, D2 | Slope (and slope of the slope) of a spectrum; Δλ = interval used, scaling factor = multiplier for readability |
| Zero-crossing | Wavelength where a derivative is zero (top of a band) |
| Amplitude | Value of a (processed) spectrum at a wavelength |
| Equality factor F | Ratio of an interferent's amplitudes at two wavelengths, used to cancel it |
| Mean centering | Subtracting the average of a spectrum; removes constants |
| LOD / LOQ | Smallest detectable / quantifiable concentration |

### 3. Which method should I try?

1. **Let Spectro choose:** *Tools → Method optimizer* ranks the methods for your
   standards and mixtures. Read its first rows, then look them up here.
2. **One compound absorbs alone somewhere?** → direct measurement for it, and
   ratio subtraction / constant multiplication / constant value for the rest.
3. **There is an isoabsorptive point?** → Q-analysis, absorbance subtraction,
   amplitude modulation, RIDSS.
4. **Complete overlap, two compounds?** → dual wavelength, derivative ratio,
   ratio difference, mean centering.
5. **Three compounds?** → amplitude centering (AAC/MACM), successive ratio
   subtraction, double divisor, ternary mean centering, or chemometrics.
6. **A minor component below its range?** → sample enrichment.

### 4. Using the sample Excel files

Every method below has a sample workbook (*Sample data*) with the spectra of
its example, ready to import:

* sheets **Standards**, **Divisor**, **Mixtures** (or **Samples**, **Calibration**)
  hold the spectra — first column wavelength, one column per spectrum; the
  sheet name becomes the role of its spectra on import;
* the column names contain the concentrations (e.g. `Std: A 8, B 0`), so
  *Fill concentrations from names* fills the table for you;
* sheets starting with `_` are notes and are not imported: **_Read me**
  (what to do), **_Calculation** (the calibration done with Excel formulas, so
  you can follow every number) and **_Expected results** (what Spectro gives).

Steps:

1. *File → New project…*, then *Edit → Compounds…* and add the compounds named in
   the workbook (e.g. A and B).
2. *File → Import spectra… → Add files…* and choose the workbook. Press
   *Fill concentrations from names*. The *Role* column is already filled from
   the sheet names. Press *Import*.
3. Open the method (see *How to apply* below) and compare your results with the
   *_Expected results* sheet.

### 5. Reading the results

* **Calibration:** r ≥ 0.999 and a small intercept; residuals scattered randomly.
* **Mixtures:** mean recovery 98–102 %, RSD ≤ 2 %.
* If a method fails for one compound only, that compound is not cancelled at
  the chosen wavelength(s) — recheck them with the finder.
"""

__all__ = ["BEGINNER", "PRIMER"]
