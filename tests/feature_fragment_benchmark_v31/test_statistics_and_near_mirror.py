from __future__ import annotations

from matching.feature_fragment_benchmark_v31.near_mirror import (
    calibrate_near_mirror_gate,
    minimum_zero_event_trials,
)
from matching.feature_fragment_benchmark_v31.statistics import (
    prospective_planning_table,
)


def test_near_mirror_zero_false_link_requirement_is_confidence_bounded() -> None:
    assert minimum_zero_event_trials(0.01, 0.95) == 299
    passing = calibrate_near_mirror_gate(
        false_links=0,
        independent_negative_pairs=299,
        true_links=29,
        independent_positive_pairs=29,
        holdout_isolated=True,
    )
    assert passing["status"] == "pass"
    assert passing["near_mirror_edges_enabled"] is True
    failing = calibrate_near_mirror_gate(
        false_links=0,
        independent_negative_pairs=298,
        true_links=29,
        independent_positive_pairs=29,
        holdout_isolated=True,
    )
    assert failing["status"] == "fail_closed"
    assert failing["near_mirror_edges_enabled"] is False


def test_frozen_power_scenario_requires_135_identities_per_arm() -> None:
    rows = prospective_planning_table(family_wise_alpha=0.05)
    selected = next(
        row
        for row in rows
        if row["holm_family_size"] == 50
        and row["background_frequency"] == 0.1
        and row["lr_positive"] == 3.0
    )
    assert selected["two_proportion_power_required_per_arm"] == 135
    assert selected["interval_precision_required_per_arm"] == 86
    assert selected["planning_required_per_arm"] == 135
