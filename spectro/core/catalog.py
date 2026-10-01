"""Method families used throughout the app and the documentation.

UV/Vis methods for mixtures are classified by **which spectrum is
manipulated** (zero order, derivative, ratio) and by **what the method
yields** (a signal at chosen wavelengths, or the recovered spectrum of each
component), with the multicomponent, chemometric and standard-addition
families kept apart:

==================  ===========================================================
zero_order          absorbance spectra: λmax, two-wavelength and factor methods
equations           several signals solved together (Vierordt, bivariate, AUC)
derivative          D1–D4 of the absorbance spectra
ratio_derivative    ratio spectra, then derivative or mean centring
ratio_amplitude     amplitudes of ratio spectra (difference, plateau, modulation,
                    amplitude centring)
resolution          spectrum resolution: the zero-order (or derivative) spectrum
                    of each component is recovered, then read at its λmax
multivariate        full-spectrum chemometric models
standard_addition   standard addition, H-point, sample enrichment
==================  ===========================================================
"""

from __future__ import annotations

CATEGORIES: dict[str, str] = {
    "zero_order": "Zero-order spectra (absorbance)",
    "equations": "Simultaneous equations (several signals)",
    "derivative": "Derivative spectra",
    "ratio_derivative": "Ratio spectra — derivative and mean centering",
    "ratio_amplitude": "Ratio spectra — amplitude methods",
    "resolution": "Spectrum resolution (recover each component's spectrum)",
    "multivariate": "Multivariate calibration (chemometrics)",
    "standard_addition": "Standard addition and sample enrichment",
}

DESCRIPTIONS: dict[str, str] = {
    "zero_order": "Work directly on the absorbance spectra: one wavelength where only the "
                  "analyte absorbs, or two wavelengths chosen so that the interferent cancels.",
    "equations": "Measure as many signals as there are compounds (or more) and solve the "
                 "linear system — all compounds at once.",
    "derivative": "Differentiate the spectra: broad overlapping bands are narrowed and a "
                  "compound can be read where the other one's derivative is zero.",
    "ratio_derivative": "Divide by the spectrum of an interferent (the divisor) — the "
                        "interferent becomes a constant — then remove that constant by "
                        "differentiation or mean centering.",
    "ratio_amplitude": "Divide by a divisor and use amplitudes of the ratio spectrum: "
                       "differences that cancel the constant, the plateau that gives the "
                       "divisor compound, or amplitudes at an isoabsorptive point / common λ.",
    "resolution": "Recover the zero-order (or derivative) spectrum of each component from the "
                  "mixture, then measure each at its own λmax as if it were pure.",
    "multivariate": "Calibrate on whole spectra of mixtures (or pure spectra) and predict all "
                    "compounds at once.",
    "standard_addition": "Add known amounts of analyte to the sample (physically or by "
                         "adding a stored spectrum) to check or correct for the matrix, or to "
                         "lift a minor component into its linear range.",
}

__all__ = ["CATEGORIES", "DESCRIPTIONS"]
