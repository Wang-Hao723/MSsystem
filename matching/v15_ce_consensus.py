"""Production-gated multi-CE, source-paired evidence for V15.

The V15 base ranker first scores the candidate pool without this module.  This
module then examines only close, same-formula candidates already present in the
V15 shortlist.  A switch requires agreement across query collision-energy
conditions, two symmetrically processed spectral views, mirror-collapsed
library provenance profiles, and a query/reference analyzer-family anchor.

Reference collision energy is never imputed.  When a reference record lacks CE
metadata it can support the analyzer-family anchor, but it cannot be reported
as an exact or near-CE reference match.  Scores are ranking utilities, not
probabilities, and a selected connectivity remains an MS/MS library candidate
rather than an authentic-standard confirmation.
"""

from __future__ import annotations

import copy
import hashlib
import json
import math
from collections import Counter
from collections.abc import Mapping, Sequence
from typing import Any

from matching.library_quality import normalize_formula_text


MODULE_ID = "v15_source_paired_platform_ce_gate"
ALGORITHM_ID = "mssystem-v15-source-paired-platform-ce-gate-1.0.0"


def _canonical(value: Any) -> bytes:
    return json.dumps(
        value,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
        default=str,
    ).encode("utf-8")


def _sha256(value: Any) -> str:
    return hashlib.sha256(_canonical(value)).hexdigest()


def _finite(value: Any) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def _identity(candidate: Mapping[str, Any] | None) -> str:
    return str((candidate or {}).get("identity_key") or "").strip()


def _formula(candidate: Mapping[str, Any] | None) -> str:
    display = (candidate or {}).get("candidate") or candidate or {}
    return normalize_formula_text(display.get("formula"))


def _winner(entity: Mapping[str, Any]) -> Mapping[str, Any] | None:
    winner = entity.get("winner")
    return winner if isinstance(winner, Mapping) else None


def _eligibility(entity: Mapping[str, Any], policy: Mapping[str, Any]) -> dict[str, Any]:
    config = policy.get("source_paired_ce_gate") or {}
    winner = _winner(entity)
    candidates = [
        candidate
        for candidate in (entity.get("candidates") or [])
        if isinstance(candidate, Mapping) and not candidate.get("vetoed")
    ]
    baseline_identity = _identity(winner)
    baseline_formula = _formula(winner)
    same_formula = [
        candidate
        for candidate in candidates
        if baseline_formula and _formula(candidate) == baseline_formula
    ]
    failures: list[str] = []
    if winner is None or not baseline_identity:
        failures.append("no_v15_selected_baseline")
    if not baseline_formula:
        failures.append("baseline_formula_unavailable")
    minimum = int(config.get("minimum_same_formula_connectivities", 2))
    maximum = int(config.get("maximum_same_formula_connectivities", 12))
    if len(same_formula) < minimum:
        failures.append("insufficient_same_formula_connectivities")
    if len(same_formula) > maximum:
        failures.append("same_formula_candidate_set_exceeds_bound")
    if int(entity.get("distinct_collision_energy_bins") or 0) < int(
        config.get("minimum_query_collision_energy_conditions", 2)
    ):
        failures.append("insufficient_query_collision_energy_conditions")
    scores = sorted(
        (
            float(candidate.get("entity_score") or 0.0),
            _identity(candidate),
        )
        for candidate in candidates
    )
    baseline_score = _finite((winner or {}).get("entity_score"))
    runner_scores = [score for score, identity in scores if identity != baseline_identity]
    margin = (
        baseline_score - max(runner_scores)
        if baseline_score is not None and runner_scores
        else None
    )
    if margin is None or margin < 0.0:
        failures.append("invalid_v15_top1_margin")
    elif margin > float(config.get("maximum_v15_top1_margin", 0.10)):
        failures.append("v15_margin_above_gate")
    return {
        "eligible": not failures,
        "failures": failures,
        "baseline_identity_key": baseline_identity or None,
        "formula": baseline_formula or None,
        "same_formula_connectivity_count": len(same_formula),
        "v15_top1_margin": round(margin, 8) if margin is not None else None,
    }


def _selected_comparison(consensus: Mapping[str, Any]) -> Mapping[str, Any] | None:
    selected = str(consensus.get("selected_identity_key") or "")
    return next(
        (
            comparison
            for comparison in (consensus.get("comparisons") or [])
            if isinstance(comparison, Mapping)
            and str(comparison.get("challenger_identity_key") or "") == selected
        ),
        None,
    )


def _production_gate(
    entity: Mapping[str, Any],
    eligibility: Mapping[str, Any],
    consensus: Mapping[str, Any],
    policy: Mapping[str, Any],
) -> dict[str, Any]:
    config = policy.get("source_paired_ce_gate") or {}
    platform_config = config.get("reference_platform_anchor") or {}
    failures = list(eligibility.get("failures") or [])
    selected_identity = str(consensus.get("selected_identity_key") or "")
    baseline_identity = str(eligibility.get("baseline_identity_key") or "")
    comparison = _selected_comparison(consensus)
    candidates = [
        candidate
        for candidate in (entity.get("candidates") or [])
        if isinstance(candidate, Mapping)
    ]
    candidate_by_identity = {_identity(candidate): candidate for candidate in candidates}
    selected_candidate = candidate_by_identity.get(selected_identity)

    if consensus.get("status") != "qualified_switch" or not comparison:
        failures.append("broad_source_paired_consensus_not_qualified")
    if str(consensus.get("baseline_identity_key") or "") != baseline_identity:
        failures.append("source_paired_baseline_differs_from_v15_baseline")
    if selected_candidate is None:
        failures.append("selected_candidate_not_in_v15_shortlist")
    elif selected_candidate.get("vetoed"):
        failures.append("selected_candidate_failed_chemistry_gate")
    elif bool((selected_candidate.get("candidate") or {}).get("potential_contaminant")):
        failures.append("selected_candidate_flagged_as_contaminant")
    if selected_candidate is not None and _formula(selected_candidate) != str(
        eligibility.get("formula") or ""
    ):
        failures.append("selected_candidate_is_not_same_formula")
    selected_rank = next(
        (
            index
            for index, candidate in enumerate(candidates, start=1)
            if _identity(candidate) == selected_identity
        ),
        None,
    )
    if selected_rank is None or selected_rank > int(
        config.get("maximum_selected_v15_rank", 5)
    ):
        failures.append("selected_candidate_rank_exceeds_bound")

    platform = (
        comparison.get("reference_platform_audit")
        if isinstance(comparison, Mapping)
        else {}
    ) or {}
    checks = (
        (
            "platform_paired_observations_below_gate",
            int(platform.get("paired_observation_count") or 0),
            int(platform_config.get("minimum_paired_observations", 6)),
        ),
        (
            "platform_provenance_strata_below_gate",
            int(platform.get("common_provenance_stratum_count") or 0),
            int(platform_config.get("minimum_common_provenance_strata", 1)),
        ),
        (
            "platform_query_ce_conditions_below_gate",
            int(platform.get("distinct_query_collision_energy_condition_count") or 0),
            int(platform_config.get("minimum_query_collision_energy_conditions", 2)),
        ),
        (
            "platform_spectral_views_below_gate",
            int(platform.get("spectral_view_count") or 0),
            int(platform_config.get("minimum_spectral_views", 2)),
        ),
        (
            "platform_unique_profiles_below_gate",
            int(platform.get("unique_provenance_evidence_profile_count") or 0),
            int(platform_config.get("minimum_unique_evidence_profiles", 1)),
        ),
    )
    for reason, observed, required in checks:
        if observed < required:
            failures.append(reason)
    float_checks = (
        (
            "platform_median_delta_below_gate",
            _finite(platform.get("median_score_delta")),
            float(platform_config.get("minimum_median_score_delta", 0.01)),
        ),
        (
            "platform_observation_win_fraction_below_gate",
            _finite(platform.get("observation_win_fraction")),
            float(platform_config.get("minimum_observation_win_fraction", 2 / 3)),
        ),
        (
            "platform_profile_win_fraction_below_gate",
            _finite(platform.get("unique_provenance_profile_win_fraction")),
            float(platform_config.get("minimum_unique_profile_win_fraction", 1.0)),
        ),
    )
    for reason, observed, required in float_checks:
        if observed is None or observed + 1e-12 < required:
            failures.append(reason)

    return {
        "module_id": MODULE_ID,
        "algorithm_id": ALGORITHM_ID,
        "status": "qualified_switch" if not failures else "baseline_retained",
        "reason": (
            "qualified_source_paired_platform_anchored_multi_ce_switch"
            if not failures
            else failures[0]
        ),
        "failed_gates": failures,
        "baseline_identity_key": baseline_identity or None,
        "selected_identity_key": selected_identity or baseline_identity or None,
        "selected_candidate_base_rank": selected_rank,
        "same_formula_only": True,
        "expanded_candidate_selection_allowed": False,
        "reference_platform_audit": copy.deepcopy(dict(platform)),
        "reference_exact_ce_audit": copy.deepcopy(
            dict((comparison or {}).get("reference_condition_audit") or {})
        ),
        "missing_reference_ce_treated_as_exact_match": False,
        "cross_unit_ce_conversion_applied": False,
        "score_is_probability": False,
        "ground_truth_used_at_inference": False,
        "production_effect": not failures,
    }


def build_v15_ce_evidence_bundle(
    v15_payload: Mapping[str, Any],
    source_rows: Sequence[Mapping[str, Any]],
    *,
    retriever: Any | None = None,
    candidate_recall_bundle: Mapping[str, Any] | None = None,
    detection_method_id: str = "unspecified",
) -> dict[str, Any]:
    """Retrieve expensive source-paired evidence only for eligible V15 entities."""

    policy = v15_payload.get("policy") or {}
    config = policy.get("source_paired_ce_gate") or {}
    if config.get("enabled") is not True:
        raise ValueError("V15 source-paired CE gate is disabled by policy")
    metadata = v15_payload.get("metadata") or {}
    expected_rows = int(metadata.get("source_row_count") or -1)
    if len(source_rows) != expected_rows:
        raise ValueError("Source row count does not match the V15 base artifact")
    recall_by_entity: dict[str, Mapping[str, Any]] | None = None
    recall_metadata: Mapping[str, Any] = {}
    if candidate_recall_bundle is not None:
        recall_metadata = candidate_recall_bundle.get("metadata") or {}
        if recall_metadata.get("base_source_semantic_sha256") != metadata.get(
            "source_semantic_sha256"
        ):
            raise ValueError("Candidate recall is bound to another source artifact")
        if recall_metadata.get("base_artifact_sha256") != _sha256(v15_payload):
            raise ValueError("Candidate recall is bound to another V15 base ranking")
        if str(recall_metadata.get("detection_method_id") or "") != str(
            detection_method_id or "unspecified"
        ):
            raise ValueError("Candidate recall uses another detection method")
        recall_by_entity = {
            str(record.get("entity_id") or ""): record
            for record in candidate_recall_bundle.get("entities") or []
            if isinstance(record, Mapping)
        }
    records: list[dict[str, Any]] = []
    retrieval_count = 0
    for entity in v15_payload.get("entities") or []:
        if not isinstance(entity, Mapping):
            continue
        eligibility = _eligibility(entity, policy)
        record: dict[str, Any] = {
            "entity_id": str(entity.get("entity_id") or ""),
            "eligibility": eligibility,
            "status": "not_applicable",
            "production_gate": {
                "module_id": MODULE_ID,
                "algorithm_id": ALGORITHM_ID,
                "status": "baseline_retained",
                "reason": (eligibility.get("failures") or ["not_eligible"])[0],
                "failed_gates": list(eligibility.get("failures") or []),
                "production_effect": False,
            },
        }
        if eligibility.get("eligible"):
            entity_id = str(entity.get("entity_id") or "")
            if recall_by_entity is not None:
                recall = recall_by_entity.get(entity_id)
                if recall is None:
                    raise ValueError("Candidate recall is missing a V15 entity")
            else:
                if retriever is None:
                    raise ValueError("V15 CE evidence requires a retriever or recall bundle")
                indices = [int(value) for value in entity.get("member_result_indices") or []]
                spectrum_ids = [str(value) for value in entity.get("member_spectrum_ids") or []]
                if len(indices) != len(spectrum_ids):
                    raise ValueError("V15 entity member index/id arrays differ in length")
                scans = []
                for index, expected_id in zip(indices, spectrum_ids):
                    if index < 0 or index >= len(source_rows):
                        raise ValueError("V15 entity member index is out of range")
                    row = source_rows[index]
                    if str(row.get("spectrum_id") or "") != expected_id:
                        raise ValueError("V15 entity member spectrum binding failed")
                    scans.append(row)
                recall = retriever.retrieve_entity(
                    entity_id,
                    scans,
                    existing_candidates=[
                        candidate
                        for candidate in (entity.get("candidates") or [])
                        if isinstance(candidate, Mapping)
                    ],
                    detection_method_id=str(detection_method_id or "unspecified"),
                )
            consensus = recall.get("source_paired_consensus") or {}
            gate = _production_gate(entity, eligibility, consensus, policy)
            record.update(
                {
                    "status": "evaluated",
                    "source_paired_consensus": consensus,
                    "production_gate": gate,
                    "scan_audit": recall.get("scan_audit") or [],
                    "expanded_candidate_count": int(
                        recall.get("expanded_candidate_count") or 0
                    ),
                    "expanded_candidates_used_for_production": False,
                }
            )
            retrieval_count += 1
        records.append(record)

    status_counts = Counter(record["production_gate"]["status"] for record in records)
    return {
        "metadata": {
            "schema_version": "mssystem-v15-ce-evidence/1.0",
            "module_id": MODULE_ID,
            "algorithm_id": ALGORITHM_ID,
            "base_task_id": metadata.get("task_id"),
            "base_source_semantic_sha256": metadata.get("source_semantic_sha256"),
            "base_artifact_sha256": _sha256(v15_payload),
            "library_snapshot_id": (
                recall_metadata.get("library_snapshot_id")
                if recall_by_entity is not None
                else getattr(retriever, "library_snapshot_id", None)
            ),
            "retrieval_policy_sha256": (
                recall_metadata.get("policy_sha256")
                if recall_by_entity is not None
                else getattr(retriever, "policy_sha256", None)
            ),
            "candidate_recall_bundle_reused": recall_by_entity is not None,
            "detection_method_id": str(detection_method_id or "unspecified"),
            "entity_count": len(records),
            "entities_retrieved": retrieval_count,
            "qualified_switch_count": status_counts.get("qualified_switch", 0),
            "production_effect": bool(status_counts.get("qualified_switch", 0)),
            "score_is_probability": False,
        },
        "invariants": {
            "ground_truth_used_at_inference": False,
            "candidate_names_affect_scoring": False,
            "expanded_candidates_can_be_selected": False,
            "same_formula_only": True,
            "mirrored_reference_content_profiles_create_extra_votes": False,
            "technical_repeats_increase_condition_count": False,
            "missing_reference_ce_treated_as_exact_match": False,
            "cross_unit_ce_conversion_applied": False,
            "score_is_probability": False,
        },
        "entities": records,
    }


def apply_v15_ce_evidence(
    v15_payload: Mapping[str, Any], evidence_bundle: Mapping[str, Any]
) -> dict[str, Any]:
    """Apply only fully qualified switches while retaining base scores for audit."""

    payload = copy.deepcopy(dict(v15_payload))
    metadata = payload.get("metadata") or {}
    evidence_metadata = evidence_bundle.get("metadata") or {}
    if evidence_metadata.get("base_source_semantic_sha256") != metadata.get(
        "source_semantic_sha256"
    ):
        raise ValueError("V15 CE evidence is bound to another source artifact")
    if evidence_metadata.get("base_artifact_sha256") != _sha256(v15_payload):
        raise ValueError("V15 CE evidence is bound to another base ranking")
    evidence_by_entity = {
        str(record.get("entity_id") or ""): record
        for record in (evidence_bundle.get("entities") or [])
        if isinstance(record, Mapping)
    }
    changed = 0
    for entity in payload.get("entities") or []:
        entity_id = str(entity.get("entity_id") or "")
        evidence = evidence_by_entity.get(entity_id)
        gate = (evidence or {}).get("production_gate") or {}
        entity["source_paired_ce_evidence"] = copy.deepcopy(evidence) if evidence else {
            "status": "unavailable",
            "production_effect": False,
        }
        entity["ranking_changed_by_source_paired_ce_gate"] = False
        if gate.get("status") != "qualified_switch":
            continue
        selected_identity = str(gate.get("selected_identity_key") or "")
        candidates = list(entity.get("candidates") or [])
        selected = next(
            (
                candidate
                for candidate in candidates
                if _identity(candidate) == selected_identity
            ),
            None,
        )
        if selected is None:
            raise ValueError("Qualified V15 CE candidate is absent from the shortlist")
        baseline = _winner(entity)
        baseline_identity = _identity(baseline)
        if baseline_identity != str(gate.get("baseline_identity_key") or ""):
            raise ValueError("Qualified V15 CE gate baseline changed before application")
        reordered = [selected, *(candidate for candidate in candidates if candidate is not selected)]
        for final_rank, candidate in enumerate(reordered, start=1):
            candidate["base_conditioned_multilibrary_rank"] = candidate.get("entity_rank")
            candidate["entity_rank"] = final_rank
            candidate["source_paired_ce_selection_qualified"] = candidate is selected
        selected["selection_basis"] = (
            "qualified_source_paired_platform_anchored_multi_ce_consensus"
        )
        selected["entity_score_retained_from_conditioned_multilibrary_ranker"] = True
        entity["candidates"] = reordered
        entity["winner"] = selected
        entity["leading_candidate_when_abstained"] = None
        entity["decision"] = {
            **dict(entity.get("decision") or {}),
            "status": "selected_candidate_only",
            "scope": "source_paired_same_formula_connectivity_candidate",
            "selection_basis": selected["selection_basis"],
            "base_conditioned_multilibrary_winner_identity_key": baseline_identity,
            "source_paired_median_score_delta": (
                (gate.get("reference_platform_audit") or {}).get("median_score_delta")
            ),
            "blockers": ["same_formula_connectivity_requires_orthogonal_confirmation"],
            "stereochemistry_claimed": False,
            "score_is_probability": False,
        }
        entity["ranking_changed_by_source_paired_ce_gate"] = (
            selected_identity != baseline_identity
        )
        changed += int(selected_identity != baseline_identity)

    metadata["entities_reranked_by_source_paired_ce_gate"] = changed
    metadata["source_paired_ce_gate_evaluated_entities"] = int(
        evidence_metadata.get("entities_retrieved") or 0
    )
    metadata["source_paired_ce_gate_library_snapshot_id"] = evidence_metadata.get(
        "library_snapshot_id"
    )
    metadata["source_paired_ce_gate_algorithm_id"] = ALGORITHM_ID
    metadata["source_paired_ce_gate_production_effect"] = bool(changed)
    payload["metadata"] = metadata
    payload.setdefault("invariants", {}).update(
        {
            "source_paired_ce_gate_ground_truth_used_at_inference": False,
            "source_paired_ce_gate_expanded_candidates_can_win": False,
            "source_paired_ce_gate_same_formula_only": True,
            "missing_reference_ce_treated_as_exact_match": False,
            "source_paired_ce_gate_score_is_probability": False,
        }
    )
    return payload


__all__ = [
    "ALGORITHM_ID",
    "MODULE_ID",
    "apply_v15_ce_evidence",
    "build_v15_ce_evidence_bundle",
]
