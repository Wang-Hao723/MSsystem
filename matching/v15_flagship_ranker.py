"""V15 production entry point with an immutable V14-compatible audit lane."""

from __future__ import annotations

import json
from collections.abc import Iterable, Mapping
from pathlib import Path
from typing import Any

from matching.flagship_ranker import build_flagship_result


MODULE_ID = "v15_conditioned_multilibrary_flagship"


def default_v15_flagship_policy_path() -> Path:
    return (
        Path(__file__).resolve().parent.parent
        / "config"
        / "v15_flagship_policy.json"
    )


def load_v15_flagship_policy(path: Path | None = None) -> dict[str, Any]:
    selected = path or default_v15_flagship_policy_path()
    payload = json.loads(selected.read_text(encoding="utf-8"))
    if payload.get("module_id") != MODULE_ID:
        raise ValueError("V15 flagship module id is not locked")
    if payload.get("deployment_state") != "production":
        raise ValueError("V15 flagship policy must declare production deployment")
    conditioned = payload.get("conditioned_multilibrary") or {}
    if conditioned.get("enabled") is not True:
        raise ValueError("V15 conditioned multi-library aggregation is required")
    if conditioned.get("cross_instrument_ce_conversion_allowed") is not False:
        raise ValueError("V15 cannot infer a cross-instrument CE conversion")
    if conditioned.get("cross_unit_ce_conversion_allowed") is not False:
        raise ValueError("V15 cannot convert unlike collision-energy units")
    if conditioned.get("ce_comparability_requires_activation_and_fixed_mode") is not True:
        raise ValueError("V15 CE comparability requires activation and fixed mode")
    if conditioned.get("stepped_ramp_or_composite_ce_reduced_to_midpoint") is not False:
        raise ValueError("V15 cannot reduce composite CE spectra to a midpoint")
    if conditioned.get("cross_vendor_qtof_ce_equivalence_assumed") is not False:
        raise ValueError("V15 cannot assume cross-vendor QTOF CE equivalence")
    if conditioned.get("missing_content_fingerprint_can_corroborate") is not False:
        raise ValueError("References without peak fingerprints cannot corroborate")
    if conditioned.get("unconditioned_consensus_affects_ranking") is not False:
        raise ValueError("All-condition consensus must remain display-only")
    if conditioned.get("fallback_when_condition_metadata_incomplete") != (
        "legacy_best_reference_with_display_only_all_condition_median"
    ):
        raise ValueError("V15 condition-metadata fallback is not fail-closed")
    if conditioned.get("score_is_probability") is not False:
        raise ValueError("V15 reference utility cannot be labelled a probability")
    ce_gate = payload.get("source_paired_ce_gate") or {}
    if ce_gate.get("enabled") is not True:
        raise ValueError("V15 source-paired CE gate is required")
    if ce_gate.get("expanded_candidates_can_be_selected") is not False:
        raise ValueError("V15 CE gate cannot promote an expanded candidate")
    if ce_gate.get("cross_instrument_conversion_allowed") is not False:
        raise ValueError("V15 CE gate cannot infer a cross-instrument conversion")
    if ce_gate.get("cross_unit_ce_conversion_allowed") is not False:
        raise ValueError("V15 CE gate cannot convert unlike CE units")
    if ce_gate.get("score_is_probability") is not False:
        raise ValueError("V15 CE consensus utility cannot be a probability")
    chemistry = payload.get("chemistry_gate") or {}
    if chemistry.get(
        "lcms_adduct_hard_veto_requires_validated_network_or_group_evidence"
    ) is not True:
        raise ValueError(
            "V15 must ignore unvalidated LC-MS adduct hypotheses for hard vetoes"
        )
    if chemistry.get(
        "mass_consistent_library_neutral_loss_can_mark_precursor_transformation_ambiguity"
    ) is not True:
        raise ValueError(
            "V15 must preserve mass-consistent precursor neutral-loss ambiguity"
        )
    if chemistry.get("precursor_transformation_ambiguity_can_confirm_identity") is not False:
        raise ValueError(
            "V15 precursor-transformation ambiguity cannot confirm identity"
        )
    if chemistry.get(
        "precursor_transformation_ambiguity_can_receive_2d_connectivity_decision"
    ) is not False:
        raise ValueError(
            "V15 precursor-transformation ambiguity cannot receive a 2D decision"
        )
    abstention = payload.get("cross_platform_abstention_gate") or {}
    if abstention.get("enabled") is not True:
        raise ValueError("V15 cross-platform abstention gate is required")
    if abstention.get("top1_can_change") is not False:
        raise ValueError("V15 abstention gate cannot change Top1")
    if abstention.get("top5_can_change") is not False:
        raise ValueError("V15 abstention gate cannot change Top5")
    if abstention.get("entity_score_can_change") is not False:
        raise ValueError("V15 abstention gate cannot change entity scores")
    if abstention.get("external_independent_validation_complete") is not False:
        raise ValueError(
            "V15 abstention calibration must not claim independent validation"
        )
    resolution = payload.get("identity_resolution") or {}
    if resolution.get("enabled") is not True:
        raise ValueError("V15 identity-scope resolution is required")
    if resolution.get("candidate_names_affect_resolution") is not False:
        raise ValueError("V15 identity resolution cannot use candidate names")
    if resolution.get("ground_truth_affects_resolution") is not False:
        raise ValueError("V15 identity resolution cannot use ground truth")
    if resolution.get("competition_set_rule") != (
        "all_non_vetoed_v15_ranked_shortlist_candidates"
    ):
        raise ValueError("V15 identity resolution must inspect the full shortlist")
    if int(resolution.get("maximum_ranked_shortlist_size") or 0) != int(
        (payload.get("candidate_pool") or {}).get("shortlist_size") or 0
    ):
        raise ValueError(
            "V15 identity-resolution bound must equal the ranking shortlist size"
        )
    for field in (
        "fixed_score_threshold_used",
        "score_ratio_threshold_used",
        "collision_energy_threshold_used_for_candidate_exclusion",
    ):
        if resolution.get(field) is not False:
            raise ValueError(f"V15 identity resolution cannot enable {field}")
    for field in ("top1_can_change", "top5_can_change", "entity_score_can_change"):
        if resolution.get(field) is not False:
            raise ValueError(f"V15 identity resolution cannot change {field}")
    forbidden = set(payload.get("forbidden_ranking_inputs") or [])
    if "feature_fragment_evidence" not in forbidden:
        raise ValueError("V15 feature-fragment evidence must be excluded from ranking")
    return payload


def build_v15_flagship_result(
    task_id: str,
    rows: Iterable[Mapping[str, Any]],
    *,
    policy: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    selected = dict(policy or load_v15_flagship_policy())
    return build_flagship_result(task_id, rows, policy=selected)


__all__ = [
    "MODULE_ID",
    "build_v15_flagship_result",
    "default_v15_flagship_policy_path",
    "load_v15_flagship_policy",
]
