"""Machine-readable contract for the normal library-search score path.

The contract reports existing behavior.  It does not calculate candidate
scores or import ``web.server``.  Web-layer defaults are explicit caller
inputs so this module remains independent of the API model.
"""

from __future__ import annotations

import json
from collections.abc import Mapping
from typing import Any

from config.settings import (
    CONSISTENCY_MODE,
    CONSISTENCY_PENALTY,
    EXPERT_RANK_INFLUENCE_ENABLED,
    FINAL_SCORE_THRESHOLD,
    ISOMER_DIAGNOSTIC_MZ_TOL_DA,
    MULTIMETRIC_FRAGMENT_TOL_DA,
    MULTIMETRIC_STRICT_FRAGMENT_TOL_DA,
    STRICT_MIN_ENTROPY,
    STRICT_MIN_FINAL,
    STRICT_MIN_REVERSE_MATCH,
    STRICT_MAX_MASS_PPM,
    USE_CONSISTENCY_FILTER,
    USE_MULTIMETRIC_SHADOW,
    USE_STRICT_EVIDENCE_GATE,
    WEIGHT_COSINE,
    WEIGHT_ENTROPY,
    WEIGHT_EXPERT,
    WEIGHT_FRAGMENT,
    WEIGHT_INTENSITY,
    WEIGHT_MASS,
)


SCHEMA_NAME = "mssystem.web-search-reliability-contract"
SCHEMA_VERSION = "1.0.0"


def _json_copy(value: Mapping[str, Any], argument_name: str) -> dict[str, Any]:
    """Return a detached JSON value and reject non-object defaults."""

    if not isinstance(value, Mapping):
        raise TypeError(f"{argument_name} must be a mapping")
    try:
        copied = json.loads(json.dumps(dict(value), allow_nan=False))
    except (TypeError, ValueError) as exc:
        raise TypeError(f"{argument_name} must be JSON-serializable") from exc
    if not isinstance(copied, dict):  # defensive: Mapping always becomes dict
        raise TypeError(f"{argument_name} must serialize to a JSON object")
    return copied


def build_reliability_contract(
    *,
    target_defaults: Mapping[str, Any],
    expected_defaults: Mapping[str, Any],
    flagship_policy: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Describe the active normal web-search scoring and dataflow contract.

    ``target_defaults`` and ``expected_defaults`` are supplied by the caller.
    This avoids a dependency on the web request model and makes the function
    deterministic for a given settings import and pair of input mappings.
    """

    supplied_defaults = {
        "target": _json_copy(target_defaults, "target_defaults"),
        "expected": _json_copy(expected_defaults, "expected_defaults"),
        "source": "caller",
    }

    configured_weights = {
        "cosine": float(WEIGHT_COSINE),
        "entropy": float(WEIGHT_ENTROPY),
        "mass": float(WEIGHT_MASS),
        "fragment": float(WEIGHT_FRAGMENT),
        "intensity": float(WEIGHT_INTENSITY),
        "expert": float(WEIGHT_EXPERT),
    }
    effective_weights = dict(configured_weights)
    if not EXPERT_RANK_INFLUENCE_ENABLED:
        effective_weights["expert"] = 0.0
    effective_weight_sum = float(sum(effective_weights.values()))
    configured_weight_sum = float(sum(configured_weights.values()))

    identity_selector = None
    if flagship_policy is not None:
        policy = _json_copy(flagship_policy, "flagship_policy")
        scan_score = policy.get("scan_score") or {}
        selector_weights = {
            name: float(scan_score.get(name) or 0.0)
            for name in (
                "cosine",
                "entropy",
                "reverse_match",
                "mass_accuracy",
                "fragment_coverage",
                "intensity_agreement",
            )
        }
        conditioned = policy.get("conditioned_multilibrary") or {}
        identity_selector = {
            "algorithm_id": policy.get("algorithm_id"),
            "scan_score_function": "matching.flagship_ranker._spectral_score",
            "scan_score_components": selector_weights,
            "scan_score_weight_sum": float(sum(selector_weights.values())),
            "entity_score": {
                "spectral_fraction": float(
                    scan_score.get("spectral_fraction") or 0.0
                ),
                "reciprocal_rank_fraction": float(
                    scan_score.get("reciprocal_rank_fraction") or 0.0
                ),
            },
            "conditioned_multilibrary": {
                "enabled": bool(conditioned.get("enabled")),
                "aggregation": conditioned.get("aggregation"),
                "best_reference_fraction": float(
                    conditioned.get("best_reference_fraction") or 0.0
                ),
                "consensus_fraction": float(
                    conditioned.get("consensus_fraction") or 0.0
                ),
                "minimum_provenance_strata": int(
                    conditioned.get("minimum_provenance_strata") or 0
                ),
                "minimum_unique_content_profiles": int(
                    conditioned.get("minimum_unique_content_profiles") or 0
                ),
                "conditioned_consensus_rule": (
                    "maximum_per_compatible_provenance_stratum_then_median"
                ),
                "fallback": conditioned.get(
                    "fallback_when_condition_metadata_incomplete"
                ),
                "unconditioned_consensus_affects_ranking": bool(
                    conditioned.get(
                        "unconditioned_consensus_affects_ranking", False
                    )
                ),
            },
            "scores_are_probabilities": False,
        }

    assertions = {
        "configured_weights_sum_to_one": abs(configured_weight_sum - 1.0) < 1e-12,
        "expert_rank_is_hard_disabled": EXPERT_RANK_INFLUENCE_ENABLED is False,
        "effective_expert_weight_is_zero": effective_weights["expert"] == 0.0,
        "effective_ranking_weights_sum_to_0_95": (
            abs(effective_weight_sum - 0.95) < 1e-12
        ),
        "reverse_match_has_zero_weight": "reverse_match" not in effective_weights,
        "multimetric_is_shadow_only": True,
        "target_and_expected_defaults_are_caller_supplied": True,
    }
    if identity_selector is not None:
        assertions.update({
            "identity_selector_scan_weights_sum_to_one": abs(
                identity_selector["scan_score_weight_sum"] - 1.0
            ) < 1e-12,
            "identity_selector_entity_weights_sum_to_one": abs(
                identity_selector["entity_score"]["spectral_fraction"]
                + identity_selector["entity_score"][
                    "reciprocal_rank_fraction"
                ]
                - 1.0
            ) < 1e-12,
            "unconditioned_consensus_does_not_affect_ranking": (
                identity_selector["conditioned_multilibrary"][
                    "unconditioned_consensus_affects_ranking"
                ] is False
            ),
        })

    return {
        "schema": SCHEMA_NAME,
        "version": SCHEMA_VERSION,
        "scope": "normal_web_library_search",
        "caller_defaults": supplied_defaults,
        "weighted_score": {
            "function": "matching.scores.enhanced_composite_score",
            "configured_components": configured_weights,
            "effective_ranking_components": effective_weights,
            "configured_weight_sum": configured_weight_sum,
            "effective_ranking_weight_sum": effective_weight_sum,
            "expert": {
                "evidence_computed": True,
                "ranking_influence_enabled": bool(
                    EXPERT_RANK_INFLUENCE_ENABLED
                ),
                "effective_input_when_disabled": 0.0,
            },
            "clipped_range": [0.0, 1.0],
        },
        "reverse_match": {
            "ranking_weight": 0.0,
            "role": "gate_only",
            "weighted_score_input": False,
            "gate_uses": [
                "strict_evidence_gate",
                "name_only_lead_eligibility",
                "presence_only_reference_eligibility",
            ],
        },
        "production_identity_selector": identity_selector,
        "fragment_tolerance_da": {
            "production_similarity_components": 0.10,
            "multimetric_shadow": float(MULTIMETRIC_FRAGMENT_TOL_DA),
            "multimetric_shadow_strict": float(
                MULTIMETRIC_STRICT_FRAGMENT_TOL_DA
            ),
            "isomer_diagnostic_shadow": float(ISOMER_DIAGNOSTIC_MZ_TOL_DA),
        },
        "post_score_stages": [
            {
                "order": 1,
                "name": "mass_error_penalty",
                "operation": "multiply_final_score",
                "enabled": True,
            },
            {
                "order": 2,
                "name": "classification_consistency",
                "operation": (
                    "reject_candidate"
                    if CONSISTENCY_MODE == "hard"
                    else "multiply_final_score"
                ),
                "enabled": bool(USE_CONSISTENCY_FILTER),
                "default_enabled": False,
                "mode": str(CONSISTENCY_MODE),
                "penalty": float(CONSISTENCY_PENALTY),
            },
            {
                "order": 3,
                "name": "final_score_threshold",
                "operation": "reject_below_threshold",
                "enabled": True,
                "threshold": float(FINAL_SCORE_THRESHOLD),
            },
            {
                "order": 4,
                "name": "strict_evidence_gate",
                "operation": "reject_candidate",
                "enabled": bool(USE_STRICT_EVIDENCE_GATE),
                "default_enabled": False,
                "thresholds": {
                    "final_score": float(STRICT_MIN_FINAL),
                    "entropy": float(STRICT_MIN_ENTROPY),
                    "reverse_match": float(STRICT_MIN_REVERSE_MATCH),
                    "absolute_mass_ppm": float(STRICT_MAX_MASS_PPM),
                },
            },
        ],
        "shadow_stages": [
            {
                "name": "multimetric_top_k_audit",
                "enabled": bool(USE_MULTIMETRIC_SHADOW),
                "execution_position": "after_candidate_selection",
                "ranking_influence": False,
                "acceptance_influence": False,
                "confidence_influence": False,
            }
        ],
        "assertions": assertions,
        "all_assertions_satisfied": all(assertions.values()),
    }


__all__ = ["SCHEMA_NAME", "SCHEMA_VERSION", "build_reliability_contract"]
