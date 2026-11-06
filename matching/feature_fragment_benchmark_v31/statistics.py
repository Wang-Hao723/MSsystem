"""Prospective per-stratum sample-size and interval-precision calculations."""

from __future__ import annotations

import math
from typing import Any

from scipy.stats import norm


Z_975 = 1.959963984540054


def expected_haldane_log_lr_interval(
    *, target_total: int, background_total: int, target_frequency: float, background_frequency: float
) -> dict[str, float]:
    target_events = target_total * target_frequency
    background_events = background_total * background_frequency
    a = target_events + 0.5
    nt = target_total + 1.0
    c = background_events + 0.5
    nb = background_total + 1.0
    point = (a / nt) / (c / nb)
    variance = max(0.0, 1.0 / a - 1.0 / nt + 1.0 / c - 1.0 / nb)
    half = Z_975 * math.sqrt(variance)
    return {
        "point": point,
        "lower_95": math.exp(math.log(point) - half),
        "upper_95": math.exp(math.log(point) + half),
        "log_half_width_95": half,
    }


def two_proportion_required_per_arm(
    *, p_background: float, p_target: float, alpha: float, power: float
) -> int:
    if not 0 < p_background < p_target < 1:
        raise ValueError("require 0 < p_background < p_target < 1")
    if not 0 < alpha < 1 or not 0 < power < 1:
        raise ValueError("alpha and power must be between zero and one")
    pooled = (p_background + p_target) / 2.0
    numerator = (
        norm.ppf(1.0 - alpha / 2.0) * math.sqrt(2.0 * pooled * (1.0 - pooled))
        + norm.ppf(power)
        * math.sqrt(
            p_background * (1.0 - p_background)
            + p_target * (1.0 - p_target)
        )
    ) ** 2
    return int(math.ceil(numerator / ((p_target - p_background) ** 2)))


def interval_precision_required_per_arm(
    *,
    p_background: float,
    p_target: float,
    maximum_log_half_width: float,
    minimum_n: int = 2,
    maximum_n: int = 100000,
) -> int:
    for n in range(minimum_n, maximum_n + 1):
        interval = expected_haldane_log_lr_interval(
            target_total=n,
            background_total=n,
            target_frequency=p_target,
            background_frequency=p_background,
        )
        if interval["lower_95"] > 1.0 and interval["log_half_width_95"] <= maximum_log_half_width:
            return n
    raise ValueError("interval planning search exceeded maximum_n")


def prospective_planning_table(
    *,
    family_wise_alpha: float,
    family_sizes: tuple[int, ...] = (10, 50, 100),
    background_frequencies: tuple[float, ...] = (0.05, 0.1, 0.2),
    likelihood_ratios: tuple[float, ...] = (2.0, 3.0, 5.0),
    power: float = 0.8,
    maximum_log_half_width: float = math.log(2.0),
    minimum_n: int = 30,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for family_size in family_sizes:
        alpha = family_wise_alpha / family_size
        for p_background in background_frequencies:
            for lr_positive in likelihood_ratios:
                p_target = p_background * lr_positive
                if p_target >= 1.0:
                    continue
                power_n = two_proportion_required_per_arm(
                    p_background=p_background,
                    p_target=p_target,
                    alpha=alpha,
                    power=power,
                )
                precision_n = interval_precision_required_per_arm(
                    p_background=p_background,
                    p_target=p_target,
                    maximum_log_half_width=maximum_log_half_width,
                    minimum_n=minimum_n,
                )
                rows.append(
                    {
                        "holm_family_size": family_size,
                        "per_test_two_sided_alpha": alpha,
                        "power": power,
                        "background_frequency": p_background,
                        "target_frequency": p_target,
                        "lr_positive": lr_positive,
                        "two_proportion_power_required_per_arm": power_n,
                        "interval_precision_required_per_arm": precision_n,
                        "planning_required_per_arm": max(minimum_n, power_n, precision_n),
                    }
                )
    return rows
