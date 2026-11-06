"""Expose observed-only replicate support without promoting a regressing ranker."""
from __future__ import annotations

import copy

from pipeline.spectral_tie_ranker import build_flagship_result

RULE = "observed_candidates_only_v16"


def attach_replicate_missingness_audit(payload, source_rows):
    config = payload.get("policy", {}).get("replicate_missingness_audit") or {}
    if config.get("enabled") is not True:
        return payload
    if (config.get("observed_support_rule") != RULE
            or config.get("affects_primary_ranking") is not False
            or config.get("missing_is_measured_zero") is not False):
        raise ValueError("Replicate audit cannot modify primary ranking or impute missing scores")
    if payload["policy"]["candidate_pool"].get("same_energy_replicate_support") != "legacy_zero_fill":
        raise ValueError("Primary replicate rule differs from the retained compatibility baseline")
    output = copy.deepcopy(payload)
    repeated = {entity["entity_id"] for entity in payload.get("entities") or []
                if any(bin_.get("technical_replicates", 0) > 1 for bin_ in entity.get("collision_energy_bins") or [])}
    affected = changes = 0
    if repeated:
        policy = copy.deepcopy(payload["policy"])
        policy["candidate_pool"]["same_energy_replicate_support"] = RULE
        alternative = build_flagship_result(payload["metadata"]["task_id"], source_rows, policy=policy)
        if alternative["metadata"]["source_semantic_sha256"] != payload["metadata"]["source_semantic_sha256"]:
            raise ValueError("Observed-only diagnostic is bound to different source scan evidence")
        by_id = {entity["entity_id"]: entity for entity in alternative["entities"]}
        if set(by_id) != {entity["entity_id"] for entity in payload["entities"]}:
            raise ValueError("Observed-only diagnostic changed precursor entity grouping")
        for entity in output["entities"]:
            if entity["entity_id"] not in repeated:
                continue
            observed = by_id[entity["entity_id"]]
            candidates = []
            for candidate in observed["candidates"]:
                evidence = candidate.get("collision_energy_evidence") or []
                candidates.append({"identity_key": candidate["identity_key"],
                    "observed_only_entity_score": candidate["entity_score"],
                    "collision_energy_evidence": [{key: value.get(key) for key in (
                        "collision_energy_bin", "support", "technical_replicates",
                        "replicates_with_candidate", "replicates_missing_from_truncated_topk",
                        "truncated_topk_absence_imputed_as_zero")} for value in evidence]})
            missing = any((bin_.get("replicates_missing_from_truncated_topk") or 0) > 0
                          for candidate in candidates for bin_ in candidate["collision_energy_evidence"])
            if not missing:
                continue
            old_id = (entity.get("winner") or {}).get("identity_key")
            new_id = (observed.get("winner") or {}).get("identity_key")
            affected += 1
            changes += old_id != new_id
            entity["replicate_missingness_audit"] = {
                "status": "incomplete_truncated_topk_vectors",
                "primary_support_rule": "legacy_zero_fill_compatibility_baseline",
                "observed_support_rule": RULE, "observed_only_candidates": candidates,
                "primary_winner_identity_key": old_id, "observed_only_winner_identity_key": new_id,
                "counterfactual_winner_differs": old_id != new_id,
                "primary_ranking_changed": False, "score_is_probability": False,
                "next_computation": "symmetric_candidate_by_query_rescoring_before_rank_promotion",
            }
    output["metadata"].update({"replicate_missingness_audit_status": "evaluated",
        "replicate_missingness_tie_algorithm_id": (
            (payload["policy"].get("spectral_tie_handling") or {}).get("algorithm_id")
        ),
        "entities_with_incomplete_replicate_candidates": affected,
        "observed_only_counterfactual_winner_changes": changes,
        "replicate_missingness_affects_primary_ranking": False})
    output["invariants"]["replicate_missingness_audit_affects_primary_ranking"] = False
    return output
