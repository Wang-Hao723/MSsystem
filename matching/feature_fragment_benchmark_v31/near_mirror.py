"""Confidence-bounded validation gate for near-mirror links."""

from __future__ import annotations

import math
from typing import Any

from scipy.stats import beta


def minimum_zero_event_trials(maximum_rate: float, confidence: float) -> int:
    if not 0 < maximum_rate < 1 or not 0 < confidence < 1:
        raise ValueError("maximum_rate and confidence must be between zero and one")
    return int(math.ceil(math.log(1.0 - confidence) / math.log(1.0 - maximum_rate)))


def clopper_pearson_upper(events: int, trials: int, confidence: float) -> float:
    if not 0 <= events <= trials or trials <= 0:
        raise ValueError("require 0 <= events <= trials and trials > 0")
    if events == trials:
        return 1.0
    return float(beta.ppf(confidence, events + 1, trials - events))


def clopper_pearson_lower(events: int, trials: int, confidence: float) -> float:
    if not 0 <= events <= trials or trials <= 0:
        raise ValueError("require 0 <= events <= trials and trials > 0")
    if events == 0:
        return 0.0
    return float(beta.ppf(1.0 - confidence, events, trials - events + 1))


def calibrate_near_mirror_gate(
    *,
    false_links: int,
    independent_negative_pairs: int,
    true_links: int,
    independent_positive_pairs: int,
    maximum_false_link_probability: float = 0.01,
    minimum_true_link_recall: float = 0.9,
    confidence: float = 0.95,
    holdout_isolated: bool,
) -> dict[str, Any]:
    false_upper = clopper_pearson_upper(
        false_links, independent_negative_pairs, confidence
    )
    recall_lower = clopper_pearson_lower(
        true_links, independent_positive_pairs, confidence
    )
    failures = []
    if not holdout_isolated:
        failures.append("calibration_and_validation_pairs_not_isolated")
    if false_upper > maximum_false_link_probability:
        failures.append("false_link_probability_upper_bound_exceeds_limit")
    if recall_lower < minimum_true_link_recall:
        failures.append("true_link_recall_lower_bound_below_limit")
    return {
        "status": "pass" if not failures else "fail_closed",
        "confidence_one_sided": confidence,
        "false_links": false_links,
        "independent_negative_pairs": independent_negative_pairs,
        "false_link_probability_upper": false_upper,
        "maximum_false_link_probability": maximum_false_link_probability,
        "true_links": true_links,
        "independent_positive_pairs": independent_positive_pairs,
        "true_link_recall_lower": recall_lower,
        "minimum_true_link_recall": minimum_true_link_recall,
        "holdout_isolated": holdout_isolated,
        "failures": failures,
        "near_mirror_edges_enabled": not failures,
    }
