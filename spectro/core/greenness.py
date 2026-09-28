"""Greenness assessment of analytical methods (AGREE, Analytical Eco-Scale, GAPI)."""

from __future__ import annotations

AGREE_PRINCIPLES = [
    "1. Direct analysis (avoid sample treatment)",
    "2. Minimal sample size and number of samples",
    "3. In-situ / on-line measurement",
    "4. Integration of processes, fewer steps",
    "5. Automation and miniaturisation",
    "6. Avoid derivatization",
    "7. Minimise waste volume",
    "8. Multi-analyte methods, high throughput",
    "9. Minimise energy consumption",
    "10. Reagents from renewable sources",
    "11. Eliminate toxic reagents",
    "12. Operator safety",
]


def agree(scores: list[float], weights: list[float] | None = None) -> dict:
    """AGREE overall score = Σ wᵢsᵢ / Σ wᵢ with each sᵢ in [0, 1]."""
    if len(scores) != 12:
        raise ValueError("AGREE needs 12 principle scores")
    weights = weights or [2.0] * 12
    if any(not 0 <= s <= 1 for s in scores):
        raise ValueError("scores must be between 0 and 1")
    total = sum(w * s for w, s in zip(weights, scores)) / sum(weights)
    return {"score": round(total, 2), "scores": scores, "weights": weights,
            "colour": score_colour(total)}


def score_colour(value: float) -> str:
    """Red (0) → yellow (0.5) → green (1), as in the AGREE pictogram."""
    v = min(max(value, 0.0), 1.0)
    if v < 0.5:
        r, g = 255, int(510 * v)
    else:
        r, g = int(510 * (1 - v)), 200 + int(55 * (1 - v))
    return f"#{r:02x}{g:02x}40"


def eco_scale(penalties: dict[str, float]) -> dict:
    """Analytical Eco-Scale: 100 − Σ penalty points.

    >75 excellent, >50 acceptable, ≤50 inadequate green analysis."""
    total = sum(penalties.values())
    score = 100 - total
    rating = "excellent" if score > 75 else "acceptable" if score > 50 else "inadequate"
    return {"score": score, "penalties": penalties, "total_penalty": total, "rating": rating}


# Reagent penalty = amount points × hazard points (Gałuszka et al. 2012)
def reagent_penalty(volume_ml: float, hazard_pictograms: int, signal_word: str = "warning") -> float:
    amount = 1 if volume_ml < 10 else 2 if volume_ml <= 100 else 3
    hazard = hazard_pictograms * (2 if signal_word.lower() == "danger" else 1)
    return float(amount * hazard)


def energy_penalty(kwh_per_sample: float) -> float:
    return 0.0 if kwh_per_sample <= 0.1 else 1.0 if kwh_per_sample <= 1.5 else 2.0


GAPI_FIELDS = [
    "Collection", "Preservation", "Transport", "Storage", "Type of method",
    "Scale of extraction", "Solvents/reagents used", "Additional treatments",
    "Amount (reagents/solvents)", "Health hazard", "Safety hazard",
    "Energy (per sample)", "Occupational hazard", "Waste", "Waste treatment",
]
GAPI_COLOURS = {"green": "#2ca02c", "yellow": "#f2c80f", "red": "#d62728"}
