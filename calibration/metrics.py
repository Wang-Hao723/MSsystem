"""Held-out probability diagnostics kept separate from target-decoy q-values."""

from __future__ import annotations

import math
from collections.abc import Sequence


def evaluate_probability_calibration(
    probabilities: Sequence[float],
    labels: Sequence[int | bool],
    *,
    evaluation_set_independent: bool,
    n_bins: int = 10,
) -> dict:
    """Compute Brier score and ECE only for declared independent labels."""
    if not evaluation_set_independent:
        raise ValueError("Brier/ECE require an independent labeled evaluation set")
    if len(probabilities) != len(labels) or not probabilities:
        raise ValueError("probabilities and labels must be non-empty and aligned")
    if n_bins < 2:
        raise ValueError("n_bins must be at least 2")
    probs = [float(value) for value in probabilities]
    ys = [int(value) for value in labels]
    if any(not math.isfinite(value) or not 0.0 <= value <= 1.0 for value in probs):
        raise ValueError("probabilities must be finite values in [0, 1]")
    if any(value not in {0, 1} for value in ys):
        raise ValueError("labels must be binary")
    brier = sum((prob - label) ** 2 for prob, label in zip(probs, ys)) / len(probs)
    bins = []
    ece = 0.0
    for index in range(n_bins):
        lower = index / n_bins
        upper = (index + 1) / n_bins
        members = [
            i
            for i, prob in enumerate(probs)
            if (
                lower <= prob <= upper
                if index == n_bins - 1
                else lower <= prob < upper
            )
        ]
        if not members:
            continue
        mean_probability = sum(probs[i] for i in members) / len(members)
        observed_fraction = sum(ys[i] for i in members) / len(members)
        weight = len(members) / len(probs)
        ece += weight * abs(mean_probability - observed_fraction)
        bins.append(
            {
                "lower": lower,
                "upper": upper,
                "count": len(members),
                "mean_probability": mean_probability,
                "observed_fraction": observed_fraction,
            }
        )
    return {
        "n": len(probs),
        "brier_score": brier,
        "expected_calibration_error": ece,
        "n_bins": n_bins,
        "bins": bins,
        "evaluation_set_independent": True,
        "qvalues_used_as_probabilities": False,
    }
