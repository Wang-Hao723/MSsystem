"""Unified, non-ranking decision audit for each V15 precursor entity."""

from __future__ import annotations

import copy
import hashlib
import json
from collections import Counter
from collections.abc import Mapping
from pathlib import Path
from typing import Any


MODULE_ID = "v15_decision_assurance"
ALGORITHM_ID = "mssystem-v15-decision-assurance-1.0.0"


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _project_root() -> Path:
    return Path(__file__).resolve().parent.parent


def load_project_calibration_state(
    *,
    registry_path: Path | None = None,
    method_path: Path | None = None,
) -> dict[str, Any]:
    """Return fail-closed project calibration state without assigning q-values."""

    root = _project_root()
    registry_file = registry_path or root / "config/calibration_registry.json"
    method_file = method_path or root / "config/calibration_method.json"
    registry = json.loads(registry_file.read_text(encoding="utf-8"))
    method = json.loads(method_file.read_text(encoding="utf-8"))
    production_id = registry.get("production_calibration_id")
    return {
        "status": (
            "registered_project_calibration_requires_entity_assignment"
            if production_id
            else "unavailable"
        ),
        "production_calibration_id": production_id,
        "method_id": method.get("method_id"),
        "method_version": method.get("method_version"),
        "required_score_field": method.get("score_field"),
        "q_value": None,
        "q_connectivity": None,
        "q_isomer_group": None,
        "q_formula": None,
        "q_class": None,
        "reason": (
            "entity_level_calibration_assignment_not_bound_to_flagship_artifact"
            if production_id
            else "no_production_project_target_decoy_calibration_registered"
        ),
        "uncalibrated_entity_exposes_qvalue": False,
        "raw_score_is_probability": False,
        "qvalue_is_probability": False,
        "calibration_changes_ranking": False,
        "registry_sha256": _sha256(registry_file),
        "method_sha256": _sha256(method_file),
    }


def attach_v15_decision_assurance(
    payload: Mapping[str, Any],
    *,
    calibration_state: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Attach one consolidated audit while preserving rank, score and scope."""

    output = copy.deepcopy(dict(payload))
    state = dict(calibration_state or load_project_calibration_state())
    before = [
        (
            str(entity.get("entity_id") or ""),
            str(((entity.get("winner") or {}).get("identity_key") or "")),
            copy.deepcopy(entity.get("decision") or {}),
            [
                (
                    candidate.get("identity_key"),
                    candidate.get("entity_rank"),
                    candidate.get("entity_score"),
                )
                for candidate in entity.get("candidates") or []
            ],
        )
        for entity in output.get("entities") or []
    ]
    assurance_statuses: Counter[str] = Counter()
    for entity in output.get("entities") or []:
        decision = entity.get("decision") or {}
        recall = entity.get("candidate_recall_support") or {}
        expanded = recall.get("expanded_candidates") or []
        query_quality = decision.get("query_spectrum_release_gate") or {}
        abstention = decision.get("cross_platform_abstention_gate") or {}
        status = (
            "candidate_only_or_abstained"
            if decision.get("status") != "selected"
            else "released_without_project_qvalue"
            if state.get("status") == "unavailable"
            else "released_with_calibration_assignment_pending"
        )
        assurance_statuses[status] += 1
        entity["decision_assurance"] = {
            "module_id": MODULE_ID,
            "algorithm_id": ALGORITHM_ID,
            "status": status,
            "decision_status": decision.get("status"),
            "decision_scope": decision.get("scope"),
            "candidate_recall": {
                "status": recall.get("status")
                or ("evaluated" if recall else "unavailable"),
                "expanded_candidate_count": len(expanded),
                "target_total_connectivities": recall.get(
                    "target_total_connectivities"
                ),
                "expanded_candidates_can_replace_winner": False,
                "expanded_candidates_can_receive_2d_connectivity_decision": False,
            },
            "query_spectrum_quality": {
                "status": query_quality.get("status") or "not_applicable",
                "deployment_state": query_quality.get("deployment_state"),
                "affects_release_scope": bool(
                    query_quality.get("affects_release_scope")
                ),
                "failed_gates": list(
                    query_quality.get("failed_gates") or []
                ),
                "quality_qualified_ce_bins": query_quality.get(
                    "observed_quality_qualified_ce_bins"
                ),
            },
            "cross_platform_abstention": {
                "status": abstention.get("status") or "not_applicable",
                "failed_gates": list(abstention.get("failed_gates") or []),
                "external_independent_validation_complete": bool(
                    abstention.get("external_independent_validation_complete")
                ),
            },
            "project_target_decoy_calibration": copy.deepcopy(state),
            "formal_claim_compiler_is_separate": True,
            "winner_rank_changed": False,
            "winner_score_changed": False,
            "decision_scope_changed": False,
            "score_is_probability": False,
            "ground_truth_used_at_inference": False,
        }

    after = [
        (
            str(entity.get("entity_id") or ""),
            str(((entity.get("winner") or {}).get("identity_key") or "")),
            copy.deepcopy(entity.get("decision") or {}),
            [
                (
                    candidate.get("identity_key"),
                    candidate.get("entity_rank"),
                    candidate.get("entity_score"),
                )
                for candidate in entity.get("candidates") or []
            ],
        )
        for entity in output.get("entities") or []
    ]
    if after != before:
        raise RuntimeError("Decision-assurance attachment changed V15 output")
    output.setdefault("metadata", {}).update(
        {
            "decision_assurance_status": "available",
            "decision_assurance_module_id": MODULE_ID,
            "decision_assurance_algorithm_id": ALGORITHM_ID,
            "decision_assurance_status_counts": dict(
                sorted(assurance_statuses.items())
            ),
            "production_project_calibration_available": bool(
                state.get("production_calibration_id")
            ),
            "entity_qvalues_exposed": 0,
        }
    )
    output.setdefault("invariants", {}).update(
        {
            "decision_assurance_changes_top1": False,
            "decision_assurance_changes_top5": False,
            "decision_assurance_changes_entity_score": False,
            "decision_assurance_changes_decision_scope": False,
            "uncalibrated_flagship_entities_expose_qvalues": False,
        }
    )
    return output


__all__ = [
    "ALGORITHM_ID",
    "MODULE_ID",
    "attach_v15_decision_assurance",
    "load_project_calibration_state",
]
