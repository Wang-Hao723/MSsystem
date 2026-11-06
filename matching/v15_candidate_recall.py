"""Bounded V15 candidate-hypothesis support selected by locked transfer tests."""

from __future__ import annotations

import copy
import hashlib
import json
from collections.abc import Mapping
from pathlib import Path
from typing import Any

from matching.evirank_retrieval import load_evirecall_policy


MODULE_ID = "v15_candidate_recall_support"
ALGORITHM_ID = "mssystem-v15-candidate-recall-support-1.0.0"


def _canonical(value: Any) -> bytes:
    return json.dumps(
        value,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
        default=str,
    ).encode("utf-8")


def semantic_sha256(value: Any) -> str:
    return hashlib.sha256(_canonical(value)).hexdigest()


def default_v15_candidate_recall_policy_path() -> Path:
    return (
        Path(__file__).resolve().parent.parent
        / "config"
        / "v15_candidate_recall_policy.json"
    )


def load_v15_candidate_recall_policy(
    path: Path | None = None,
) -> dict[str, Any]:
    selected = path or default_v15_candidate_recall_policy_path()
    overlay = json.loads(selected.read_text(encoding="utf-8"))
    if overlay.get("module_id") != MODULE_ID:
        raise ValueError("V15 candidate-recall module id is not locked")
    if overlay.get("algorithm_id") != ALGORITHM_ID:
        raise ValueError("V15 candidate-recall algorithm id is not locked")
    base_path = selected.parent.parent / "config" / "v12_candidate_recall_policy.json"
    expected_hash = str(overlay.get("base_policy_file_sha256") or "")
    observed_hash = hashlib.sha256(base_path.read_bytes()).hexdigest()
    if observed_hash != expected_hash:
        raise ValueError("V15 candidate-recall base policy hash changed")
    policy = load_evirecall_policy(base_path)
    overrides = overlay.get("candidate_pool_overrides") or {}
    policy["algorithm_id"] = ALGORITHM_ID
    policy["candidate_pool"]["target_total_connectivities"] = int(
        overrides.get("target_total_connectivities") or 0
    )
    policy["candidate_pool"]["maximum_additions_per_entity"] = int(
        overrides.get("maximum_additions_per_entity") or 0
    )
    policy["v15_integration"] = copy.deepcopy(overlay)
    if policy["candidate_pool"]["target_total_connectivities"] != 40:
        raise ValueError("V15 candidate recall must retain the selected Top40 bound")
    if policy["candidate_pool"]["maximum_additions_per_entity"] != 40:
        raise ValueError("V15 candidate recall addition bound must remain 40")
    if [view.get("name") for view in policy.get("spectral_views") or []] != [
        "library_match",
        "sensitive",
    ]:
        raise ValueError("V15 candidate recall must retain the tested two views")
    for field in (
        "expanded_candidates_can_replace_winner",
        "expanded_candidates_can_receive_2d_connectivity_decision",
        "expanded_candidates_enter_residual_model_scoring",
        "score_is_probability",
        "ground_truth_used_at_inference",
    ):
        if overlay.get(field) is not False:
            raise ValueError(f"V15 candidate-recall invariant changed: {field}")
    return policy


def attach_v15_candidate_recall(
    payload: Mapping[str, Any],
    recall_bundle: Mapping[str, Any],
    *,
    base_payload: Mapping[str, Any],
) -> dict[str, Any]:
    """Attach expanded candidates without changing rank, score, or decision."""

    base_metadata = base_payload.get("metadata") or {}
    recall_metadata = recall_bundle.get("metadata") or {}
    if recall_metadata.get("base_source_semantic_sha256") != base_metadata.get(
        "source_semantic_sha256"
    ):
        raise ValueError("V15 candidate recall is bound to another source artifact")
    if recall_metadata.get("base_artifact_sha256") != semantic_sha256(base_payload):
        raise ValueError("V15 candidate recall is bound to another base ranking")
    selected_policy = load_v15_candidate_recall_policy()
    expected_policy_hash = semantic_sha256(selected_policy)
    if recall_metadata.get("policy_sha256") != expected_policy_hash:
        raise ValueError("V15 candidate-recall policy binding changed")

    output = copy.deepcopy(dict(payload))
    before = [
        (
            str(entity.get("entity_id") or ""),
            str(((entity.get("winner") or {}).get("identity_key") or "")),
            copy.deepcopy(entity.get("decision") or {}),
            [
                (
                    str(candidate.get("identity_key") or ""),
                    candidate.get("entity_score"),
                    candidate.get("entity_rank"),
                )
                for candidate in (entity.get("candidates") or [])
            ],
        )
        for entity in output.get("entities") or []
    ]
    records = {
        str(record.get("entity_id") or ""): record
        for record in recall_bundle.get("entities") or []
        if isinstance(record, Mapping)
    }
    expected_ids = {
        str(entity.get("entity_id") or "")
        for entity in output.get("entities") or []
    }
    if set(records) != expected_ids:
        raise ValueError("V15 candidate-recall entity membership changed")

    expanded_total = 0
    entities_with_expansion = 0
    for entity in output.get("entities") or []:
        entity_id = str(entity.get("entity_id") or "")
        record = copy.deepcopy(dict(records[entity_id]))
        expanded = []
        for candidate in record.get("expanded_candidates") or []:
            item = copy.deepcopy(dict(candidate))
            item.update(
                {
                    "claim_scope": "candidate_hypothesis_support_only",
                    "can_replace_winner": False,
                    "can_receive_2d_connectivity_decision": False,
                    "requires_primary_ranking_or_orthogonal_evidence": True,
                }
            )
            expanded.append(item)
        record["expanded_candidates"] = expanded
        record.update(
            {
                "integration_module_id": MODULE_ID,
                "integration_algorithm_id": ALGORITHM_ID,
                "claim_scope": "candidate_hypothesis_support_only",
                "expanded_candidates_can_replace_winner": False,
                "expanded_candidates_can_receive_2d_connectivity_decision": False,
                "production_identity_rank_unchanged": True,
            }
        )
        entity["candidate_recall_support"] = record
        expanded_total += len(expanded)
        entities_with_expansion += int(bool(expanded))

    after = [
        (
            str(entity.get("entity_id") or ""),
            str(((entity.get("winner") or {}).get("identity_key") or "")),
            copy.deepcopy(entity.get("decision") or {}),
            [
                (
                    str(candidate.get("identity_key") or ""),
                    candidate.get("entity_score"),
                    candidate.get("entity_rank"),
                )
                for candidate in (entity.get("candidates") or [])
            ],
        )
        for entity in output.get("entities") or []
    ]
    if after != before:
        raise RuntimeError("Candidate-recall attachment changed the primary ranking")

    output.setdefault("metadata", {}).update(
        {
            "candidate_recall_support_status": "evaluated",
            "candidate_recall_support_module_id": MODULE_ID,
            "candidate_recall_support_algorithm_id": ALGORITHM_ID,
            "candidate_recall_support_policy_sha256": expected_policy_hash,
            "candidate_recall_support_library_snapshot_id": recall_metadata.get(
                "library_snapshot_id"
            ),
            "candidate_recall_target_total_connectivities": 40,
            "entities_with_candidate_recall_expansion": entities_with_expansion,
            "expanded_candidate_hypothesis_count": expanded_total,
            "candidate_recall_support_affects_primary_ranking": False,
        }
    )
    output.setdefault("invariants", {}).update(
        {
            "candidate_recall_support_ground_truth_used_at_inference": False,
            "candidate_recall_support_affects_primary_ranking": False,
            "candidate_recall_support_affects_primary_score": False,
            "expanded_candidates_can_replace_winner": False,
            "expanded_candidates_can_receive_2d_connectivity_decision": False,
            "expanded_candidates_enter_residual_model_scoring": False,
            "candidate_recall_score_is_probability": False,
        }
    )
    return output


__all__ = [
    "ALGORITHM_ID",
    "MODULE_ID",
    "attach_v15_candidate_recall",
    "default_v15_candidate_recall_policy_path",
    "load_v15_candidate_recall_policy",
    "semantic_sha256",
]
