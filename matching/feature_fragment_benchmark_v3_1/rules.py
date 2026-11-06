"""Frozen V3 rule procedures bound to the external V3.1 spectral-QC policy."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any

from matching.feature_fragment_benchmark_v3.rules import (
    calibrate_rules,
    discover_rules,
    evaluate_sealed_rules,
)


def _validate_frozen_qc_policy(qc_policy: Mapping[str, Any]) -> None:
    peak = qc_policy["peak_presence"]
    tolerance = peak["mass_tolerance"]
    expected = (
        float(peak["minimum_relative_intensity"]),
        str(tolerance["rule"]),
        float(tolerance["absolute_da"]),
        float(tolerance["ppm"]),
    )
    if expected != (0.01, "maximum_of_absolute_da_and_ppm", 0.01, 20.0):
        raise ValueError("v3_rule_engine_qc_policy_mismatch_requires_new_implementation")


def discover_rules_with_policy(
    observations: Sequence[Mapping[str, Any]],
    *,
    qc_policy: Mapping[str, Any],
    minimum_target_identities: int = 20,
    minimum_background_identities: int = 20,
    maximum_q_value: float = 0.05,
    minimum_lr_positive: float = 2.0,
) -> list[dict[str, Any]]:
    """Delegate only after the external policy matches the frozen V3 engine."""
    _validate_frozen_qc_policy(qc_policy)
    rules = discover_rules(
        observations,
        minimum_target_identities=minimum_target_identities,
        minimum_background_identities=minimum_background_identities,
        maximum_q_value=maximum_q_value,
        minimum_lr_positive=minimum_lr_positive,
    )
    for rule in rules:
        rule["schema_version"] = "3.1.0"
        rule["spectral_qc_policy_id"] = qc_policy["policy_id"]
        rule["mass_tolerance"]["policy"] = "maximum_of_absolute_da_and_ppm"
    return rules


__all__ = [
    "calibrate_rules",
    "discover_rules_with_policy",
    "evaluate_sealed_rules",
]
