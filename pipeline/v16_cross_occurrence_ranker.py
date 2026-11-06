"""Controlled V16 cross-occurrence consensus for an explicit LC-MS cohort.

This result-postprocessing module can only reorder an entity's existing,
non-vetoed, same-formula
Top-5 candidates.  Alignment uses measured m/z, RT, polarity and task identity;
candidate identities and scores do not participate in edge construction.
Candidate names, sample filenames and reference-standard labels are never read.
"""

from __future__ import annotations

import copy
import hashlib
import json
import math
import re
import statistics
from collections import Counter, defaultdict
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from matching.conditioned_multilibrary import reference_content_fingerprint
from matching.evidence_governance import _cluster_features
from matching.library_quality import normalize_formula_text
from matching.scores import (
    batch_neutral_loss_cosine_scores,
)
from pipeline.v16_complete_candidate_panel import project_query
from utils.spectrum import process_peaks, resolve_peak_processing_config


MODULE_ID = "v16_cross_occurrence_consensus_ranker"
ALGORITHM_ID = "mssystem-v16-cross-occurrence-consensus-top5-1.0.0"


def _canonical(value: Any) -> bytes:
    return json.dumps(
        value,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
        default=str,
    ).encode("utf-8")


def semantic_sha256(value: Any) -> str:
    return hashlib.sha256(_canonical(value)).hexdigest()


def default_policy_path() -> Path:
    return Path(__file__).resolve().parents[1] / "config" / "v16_cross_occurrence_policy.json"


def load_policy(path: Path | None = None) -> dict[str, Any]:
    selected = path or default_policy_path()
    policy = json.loads(selected.read_text(encoding="utf-8"))
    if policy.get("module_id") != MODULE_ID or policy.get("algorithm_id") != ALGORITHM_ID:
        raise ValueError("V16 cross-occurrence policy identity changed")
    if policy.get("deployment_state") != "controlled_opt_in":
        raise ValueError("Cross-occurrence ranking must remain controlled opt-in")
    activation = policy.get("activation") or {}
    if activation.get("explicit_completed_task_ids_required") is not True:
        raise ValueError("Explicit completed task IDs are required")
    if activation.get("same_lc_method_confirmation_required") is not True:
        raise ValueError("Same-method confirmation is required")
    if activation.get("automatic_filename_grouping_allowed") is not False:
        raise ValueError("Filename grouping must remain disabled")
    if int(activation.get("minimum_task_count") or 0) != 3:
        raise ValueError("At least three independently submitted tasks are required")
    maximum = int(activation.get("maximum_task_count") or 0)
    if maximum < 3 or maximum > 100:
        raise ValueError("Cross-occurrence task bound is invalid")
    alignment = policy.get("alignment") or {}
    if float(alignment.get("mz_tolerance_ppm") or 0.0) != 15.0:
        raise ValueError("The validated precursor tolerance changed")
    if float(alignment.get("rt_window_min") or 0.0) != 0.15:
        raise ValueError("The validated retention-time window changed")
    if alignment.get("linkage") != "complete_link":
        raise ValueError("Complete-link alignment is required")
    scoring = policy.get("candidate_scoring") or {}
    if scoring.get("candidate_scope") != "existing_non_vetoed_native_ranked_shortlist":
        raise ValueError("Only the existing native shortlist may be reranked")
    if int(scoring.get("maximum_native_candidates") or 0) != 5:
        raise ValueError("The native candidate bound changed")
    if scoring.get("formula_lock") != "current_ranking_lead_formula":
        raise ValueError("The current ranking-lead formula lock is required")
    if list(scoring.get("views") or []) != ["library_match", "sensitive"]:
        raise ValueError("The two validated spectral views changed")
    if float(scoring.get("direct_weight") or 0.0) != 0.5:
        raise ValueError("Direct-spectrum weight changed")
    if float(scoring.get("neutral_loss_weight") or 0.0) != 0.5:
        raise ValueError("Neutral-loss weight changed")
    if scoring.get("expanded_candidates_can_replace_winner") is not False:
        raise ValueError("Expanded candidates may not replace the winner")
    peer = policy.get("peer_rule") or {}
    if int(peer.get("minimum_peer_tasks") or 0) != 2:
        raise ValueError("The peer-task threshold changed")
    if not math.isclose(float(peer.get("minimum_winner_vote_fraction") or 0.0), 2.0 / 3.0):
        raise ValueError("The peer vote threshold changed")
    if peer.get("target_direct_top1_must_equal_peer_winner") is not True:
        raise ValueError("Target direct agreement must remain required")
    boundary = policy.get("claim_boundary") or {}
    required_false = (
        "decision_scope_may_be_upgraded",
        "authentic_standard_confirmation_created",
        "candidate_score_changed",
        "compound_count_may_increase",
        "candidate_names_used_for_alignment_or_scoring",
        "ground_truth_used_at_inference",
    )
    if any(boundary.get(field) is not False for field in required_false):
        raise ValueError("A cross-occurrence claim boundary changed")
    return policy


def _finite(value: Any) -> float | None:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def connectivity_key(item: Mapping[str, Any] | None) -> str:
    if not isinstance(item, Mapping):
        return ""
    display = item.get("candidate") if isinstance(item.get("candidate"), Mapping) else {}
    value = str(
        item.get("connectivity_identity")
        or item.get("identity_key")
        or display.get("standardized_connectivity_inchikey")
        or display.get("inchikey")
        or ""
    ).strip()
    value = value.split(":")[-1].split("-")[0].upper()
    return value if re.fullmatch(r"[A-Z]{14}", value) else ""


def candidate_formula(item: Mapping[str, Any]) -> str:
    display = item.get("candidate") if isinstance(item.get("candidate"), Mapping) else item
    return normalize_formula_text(display.get("formula") or display.get("molecular_formula"))


def project_entity(task_id: str, entity: Mapping[str, Any]) -> dict[str, Any] | None:
    """Return the identity-minimal event used by the cohort ranker."""

    selected = entity.get("winner") or entity.get("leading_candidate_when_abstained")
    baseline = connectivity_key(selected if isinstance(selected, Mapping) else None)
    mz = _finite(
        entity.get("representative_search_precursor_mz")
        or entity.get("representative_precursor_mz")
    )
    rt = _finite(entity.get("representative_retention_time"))
    modes = {
        str(value or "").strip().upper()
        for value in (entity.get("ion_modes") or [])
        if str(value or "").strip().upper() in {"POS", "NEG"}
    }
    if not baseline or mz is None or mz <= 0.0 or rt is None or len(modes) != 1:
        return None
    maximum = 5
    candidates: list[dict[str, Any]] = []
    seen: set[str] = set()
    for fallback_rank, item in enumerate(entity.get("candidates") or [], 1):
        if not isinstance(item, Mapping) or bool(item.get("vetoed")):
            continue
        identity = connectivity_key(item)
        if not identity or identity in seen:
            continue
        seen.add(identity)
        candidates.append(
            {
                "connectivity": identity,
                "formula": candidate_formula(item),
                "rank": int(item.get("entity_rank") or fallback_rank),
            }
        )
    candidates.sort(key=lambda item: (item["rank"], item["connectivity"]))
    if len(candidates) > maximum:
        raise ValueError("A V16 entity exceeds the locked native Top-5 candidate bound")
    baseline_row = next((item for item in candidates if item["connectivity"] == baseline), None)
    if baseline_row is None or not baseline_row["formula"]:
        return None
    formula_candidates = [
        item for item in candidates if item["formula"] == baseline_row["formula"]
    ]
    if not formula_candidates:
        return None
    return {
        "event_id": f"{task_id}::{entity.get('entity_id')}",
        "task_id": str(task_id),
        "entity_id": str(entity.get("entity_id") or ""),
        "polarity": next(iter(modes)),
        "measured_precursor_mz": mz,
        "retention_time_min": rt,
        "baseline_formula": baseline_row["formula"],
        "baseline_top1_connectivity": baseline,
        "native_same_formula_candidates": formula_candidates,
        "member_result_indices": [int(value) for value in entity.get("member_result_indices") or []],
    }


def cluster_events(
    events: Sequence[Mapping[str, Any]], policy: Mapping[str, Any]
) -> list[list[dict[str, Any]]]:
    features: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for event in events:
        features[str(event["task_id"])].append(
            {
                "task_id": str(event["task_id"]),
                "feature_id": str(event["event_id"]),
                "polarity": str(event["polarity"]),
                "precursor_mz": float(event["measured_precursor_mz"]),
                "retention_time_min": float(event["retention_time_min"]),
            }
        )
    alignment = policy["alignment"]
    return _cluster_features(
        dict(features),
        mz_ppm=float(alignment["mz_tolerance_ppm"]),
        rt_window=float(alignment["rt_window_min"]),
    )


def score_event_direct_neutral(
    event: Mapping[str, Any],
    source_rows: Sequence[Mapping[str, Any]],
    retriever: Any,
    policy: Mapping[str, Any],
) -> dict[str, Any]:
    """Score the frozen same-formula native panel without reading names."""

    candidates = list(event.get("native_same_formula_candidates") or [])
    baseline = str(event.get("baseline_top1_connectivity") or "")
    queries = [project_query(row) for row in source_rows]
    if not queries:
        raise ValueError("An aligned entity has no saved query spectra")
    candidate_ids = {str(item["connectivity"]) for item in candidates}
    ranks = {str(item["connectivity"]): int(item["rank"]) for item in candidates}
    query_fingerprints = {
        fingerprint
        for fingerprint in (reference_content_fingerprint(query["peaks"]) for query in queries)
        if fingerprint
    }
    view_policy = {
        str(view.get("name") or ""): view
        for view in retriever.policy.get("spectral_views") or []
    }
    expected_views = list((policy.get("candidate_scoring") or {}).get("views") or [])
    if list(view_policy) != expected_views:
        raise ValueError("Runtime library scorer does not expose the locked two views")
    grouped: dict[str, dict[tuple[str, str, str], list[float]]] = {
        identity: defaultdict(list) for identity in candidate_ids
    }
    scored_pairs = 0
    for query in queries:
        scan = retriever.retrieve_scan(query)
        records, _audit = retriever._fetch_window(query["precursor_mz"], query["polarity"])
        references = {str(record.get("spectrum_id") or ""): record for record in records}
        for identity in sorted(candidate_ids):
            record = (scan.get("candidates") or {}).get(identity)
            if not isinstance(record, Mapping):
                continue
            for view_name in expected_views:
                direct = _finite((record.get("view_scores") or {}).get(view_name))
                reference_stub = (record.get("reference_records") or {}).get(view_name) or {}
                spectrum_id = str(reference_stub.get("spectrum_id") or "")
                reference = references.get(spectrum_id)
                fingerprint = str(reference_stub.get("content_fingerprint") or "")
                if direct is None or reference is None or (fingerprint and fingerprint in query_fingerprints):
                    continue
                view = view_policy[view_name]
                config = resolve_peak_processing_config(view["peak_processing_preset"])
                tolerance = float(view["fragment_tolerance_da"])
                query_peaks, _ = process_peaks(
                    query["peaks"],
                    query["precursor_mz"],
                    config=config,
                    precursor_charge=query.get("charge"),
                )
                reference_peaks, _ = process_peaks(
                    reference.get("peaks") or [],
                    reference.get("precursor_mz"),
                    config=config,
                    precursor_charge=reference.get("charge"),
                )
                neutral = float(
                    batch_neutral_loss_cosine_scores(
                        query_peaks,
                        [reference_peaks],
                        query["precursor_mz"],
                        [reference.get("precursor_mz")],
                        mz_tol=tolerance,
                    )[0]
                )
                condition_id = str(scan.get("condition_id") or "")
                profile = fingerprint or spectrum_id
                grouped[identity][(condition_id, view_name, profile)].append(
                    0.5 * float(direct) + 0.5 * neutral
                )
                scored_pairs += 1
    scores: dict[str, float] = {}
    for identity, groups in grouped.items():
        collapsed = [statistics.median(values) for values in groups.values() if values]
        if collapsed:
            scores[identity] = round(float(statistics.median(collapsed)), 8)
    if baseline not in scores:
        direct_top1 = baseline
    else:
        direct_top1 = min(
            scores,
            key=lambda identity: (-scores[identity], ranks.get(identity, 10**9), identity),
        )
    return {
        "direct_top1_connectivity": direct_top1,
        "candidate_scores": dict(sorted(scores.items())),
        "scored_pair_count": scored_pairs,
        "baseline_score_available": baseline in scores,
        "single_candidate_shortcut": False,
    }


def build_consensus_plan(
    events: Sequence[Mapping[str, Any]], policy: Mapping[str, Any]
) -> dict[str, Any]:
    """Create deterministic proposals from already scored, truth-free events."""

    by_id = {str(event["event_id"]): dict(event) for event in events}
    clusters = cluster_events(events, policy)
    peer_rule = policy["peer_rule"]
    minimum_peers = int(peer_rule["minimum_peer_tasks"])
    fraction = float(peer_rule["minimum_winner_vote_fraction"])
    proposals: list[dict[str, Any]] = []
    cluster_audits: list[dict[str, Any]] = []
    for number, cluster in enumerate(clusters, 1):
        members = [by_id[str(feature["feature_id"])] for feature in cluster]
        cluster_id = f"XC{number:06d}"
        cluster_audits.append(
            {
                "cluster_id": cluster_id,
                "event_ids": sorted(str(member["event_id"]) for member in members),
                "task_count": len({str(member["task_id"]) for member in members}),
            }
        )
        for target in members:
            peers = [
                peer
                for peer in members
                if peer["event_id"] != target["event_id"]
                and peer["baseline_formula"] == target["baseline_formula"]
            ]
            scored_peers = [
                peer
                for peer in peers
                if peer.get("baseline_score_available") is True
            ]
            votes = Counter(
                str(peer.get("direct_top1_connectivity") or "")
                for peer in scored_peers
                if str(peer.get("direct_top1_connectivity") or "")
            )
            ordered = sorted(votes.items(), key=lambda item: (-item[1], item[0]))
            failures: list[str] = []
            peer_winner = ""
            vote_count = sum(votes.values())
            winner_fraction = None
            if target.get("baseline_score_available") is not True:
                failures.append("target_baseline_direct_neutral_score_unavailable")
            if vote_count < minimum_peers:
                failures.append("insufficient_same_formula_peer_tasks")
            elif not ordered:
                failures.append("no_peer_vote")
            elif len(ordered) > 1 and ordered[0][1] == ordered[1][1]:
                failures.append("peer_vote_tie")
            else:
                peer_winner = ordered[0][0]
                winner_fraction = ordered[0][1] / vote_count
                if winner_fraction + 1e-12 < fraction:
                    failures.append("peer_vote_fraction_below_two_thirds")
            native = {
                str(item["connectivity"])
                for item in target.get("native_same_formula_candidates") or []
            }
            if peer_winner and peer_winner not in native:
                failures.append("peer_winner_absent_from_target_native_same_formula_top5")
            direct = str(target.get("direct_top1_connectivity") or "")
            if peer_winner and direct != peer_winner:
                failures.append("target_direct_top1_disagrees_with_peer_winner")
            baseline = str(target["baseline_top1_connectivity"])
            if peer_winner == baseline:
                failures.append("consensus_retains_current_top1")
            applied = bool(peer_winner and not failures)
            proposals.append(
                {
                    "event_id": str(target["event_id"]),
                    "task_id": str(target["task_id"]),
                    "entity_id": str(target["entity_id"]),
                    "cluster_id": cluster_id,
                    "baseline_top1_connectivity": baseline,
                    "proposed_top1_connectivity": peer_winner if applied else baseline,
                    "direct_top1_connectivity": direct,
                    "baseline_formula": str(target["baseline_formula"]),
                    "aligned_peer_count": len(members) - 1,
                    "same_formula_peer_count": len(peers),
                    "scored_same_formula_peer_count": len(scored_peers),
                    "vote_count": vote_count,
                    "votes": [
                        {"connectivity": identity, "count": count}
                        for identity, count in ordered
                    ],
                    "winner_vote_fraction": winner_fraction,
                    "candidate_scores": dict(target.get("candidate_scores") or {}),
                    "applied": applied,
                    "failure_reasons": failures,
                }
            )
    return {
        "schema_version": "mssystem.v16.cross_occurrence_consensus.plan/1.0.0",
        "module_id": MODULE_ID,
        "algorithm_id": ALGORITHM_ID,
        "event_count": len(events),
        "cluster_count": len(clusters),
        "multi_task_cluster_count": sum(row["task_count"] > 1 for row in cluster_audits),
        "applied_count": sum(bool(row["applied"]) for row in proposals),
        "clusters": cluster_audits,
        "proposals": sorted(proposals, key=lambda row: row["event_id"]),
        "candidate_names_used": False,
        "ground_truth_used": False,
        "score_is_probability": False,
    }


def entity_semantic_sha256(payload: Mapping[str, Any]) -> str:
    digest = hashlib.sha256()
    for entity in payload.get("entities") or []:
        digest.update(_canonical(entity))
        digest.update(b"\n")
    return digest.hexdigest()


def apply_consensus_to_payload(
    payload: Mapping[str, Any],
    proposals: Sequence[Mapping[str, Any]],
    *,
    cohort_id: str,
    policy: Mapping[str, Any],
) -> tuple[dict[str, Any], int]:
    """Apply verified proposals while preserving candidate membership/scores."""

    output = copy.deepcopy(dict(payload))
    by_entity = {
        str(row["entity_id"]): row
        for row in proposals
        if row.get("applied") and str(row.get("task_id") or "") == str((payload.get("metadata") or {}).get("task_id") or "")
    }
    changed = 0
    for entity in output.get("entities") or []:
        if not isinstance(entity, dict):
            continue
        proposal = by_entity.get(str(entity.get("entity_id") or ""))
        if proposal is None:
            continue
        candidates = [item for item in entity.get("candidates") or [] if isinstance(item, dict)]
        before_pairs = sorted(
            (connectivity_key(item), item.get("entity_score")) for item in candidates
        )
        baseline = str(proposal["baseline_top1_connectivity"])
        promoted = str(proposal["proposed_top1_connectivity"])
        current_selected = entity.get("winner") or entity.get("leading_candidate_when_abstained")
        if connectivity_key(current_selected) != baseline:
            raise ValueError("Flagship ranking changed after the consensus plan was built")
        selected = next((item for item in candidates if connectivity_key(item) == promoted), None)
        if selected is None or bool(selected.get("vetoed")):
            raise ValueError("The proposed candidate is absent or vetoed")
        if candidate_formula(selected) != str(proposal["baseline_formula"]):
            raise ValueError("Cross-occurrence proposal crossed a molecular-formula boundary")
        previous_lead = next(
            (item for item in candidates if connectivity_key(item) == baseline),
            None,
        )
        prior_semantics = (
            (current_selected or {}).get("reporting_semantics")
            if isinstance(current_selected, Mapping)
            else None
        )
        # reporting_semantics marks the displayed ranking lead, not an intrinsic
        # molecular property.  Move it with the lead so the reordered candidate
        # list and the winner projection cannot disagree.
        if previous_lead is not None and previous_lead is not selected:
            previous_lead.pop("reporting_semantics", None)
        if prior_semantics:
            selected["reporting_semantics"] = prior_semantics
        reordered = [selected] + [item for item in candidates if item is not selected]
        for rank, item in enumerate(reordered, 1):
            item["entity_rank"] = rank
        if sorted((connectivity_key(item), item.get("entity_score")) for item in reordered) != before_pairs:
            raise RuntimeError("Cross-occurrence ranking changed candidate scores or membership")
        previous_margin = (entity.get("decision") or {}).get("winner_margin")
        replacement = copy.deepcopy(selected)
        entity["candidates"] = reordered
        if entity.get("winner") is not None:
            entity["winner"] = replacement
        else:
            entity["leading_candidate_when_abstained"] = replacement
        decision = dict(entity.get("decision") or {})
        decision.update(
            {
                "ranking_lead_identity_key": replacement.get("identity_key"),
                "winner_margin": None,
                "spectral_winner_margin_before_cross_occurrence": previous_margin,
                "winner_margin_not_applicable_after_cross_occurrence": True,
                "ranking_basis": ALGORITHM_ID,
                "cross_occurrence_consensus_applied": True,
            }
        )
        entity["decision"] = decision
        resolution = dict(entity.get("identity_resolution") or {})
        if resolution:
            resolution["ranking_lead_identity_key"] = replacement.get("identity_key")
            resolution["competitive_candidate_identity_keys"] = [
                item.get("identity_key") for item in reordered if not item.get("vetoed")
            ]
            resolution["ranking_changed"] = False
            entity["identity_resolution"] = resolution
        entity["ranking_changed_from_scan_top1"] = bool(
            replacement.get("identity_key")
            != (entity.get("scan_top1_audit_baseline") or {}).get("identity_key")
        )
        entity["ranking_changed_from_v14"] = bool(
            replacement.get("identity_key")
            != (entity.get("v14_audit_baseline") or {}).get("identity_key")
        )
        entity["ranking_changed_by_cross_occurrence_consensus"] = True
        entity["cross_occurrence_consensus"] = {
            "module_id": MODULE_ID,
            "algorithm_id": ALGORITHM_ID,
            "cohort_id": cohort_id,
            "cluster_id": proposal["cluster_id"],
            "baseline_top1_connectivity": baseline,
            "promoted_top1_connectivity": promoted,
            "same_formula": True,
            "same_formula_peer_count": proposal["same_formula_peer_count"],
            "votes": copy.deepcopy(proposal["votes"]),
            "winner_vote_fraction": proposal["winner_vote_fraction"],
            "target_direct_top1_connectivity": proposal["direct_top1_connectivity"],
            "candidate_scores": copy.deepcopy(proposal.get("candidate_scores") or {}),
            "candidate_score_changed": False,
            "decision_scope_changed": False,
            "authentic_standard_confirmation_created": False,
            "score_is_probability": False,
        }
        assurance = dict(entity.get("decision_assurance") or {})
        if assurance:
            assurance["cross_occurrence_consensus"] = {
                "applied": True,
                "cohort_id": cohort_id,
                "decision_scope_changed": False,
                "standard_confirmed": False,
            }
            entity["decision_assurance"] = assurance
        changed += 1
    metadata = output.setdefault("metadata", {})
    metadata.update(
        {
            "cross_occurrence_consensus_status": "activated",
            "cross_occurrence_consensus_module_id": MODULE_ID,
            "cross_occurrence_consensus_algorithm_id": ALGORITHM_ID,
            "cross_occurrence_consensus_policy_sha256": semantic_sha256(policy),
            "cross_occurrence_consensus_cohort_id": cohort_id,
            "entities_reranked_by_cross_occurrence_consensus": changed,
            "entities_reranked_from_first_scan": sum(
                bool(entity.get("ranking_changed_from_scan_top1"))
                for entity in output.get("entities") or []
            ),
            "entities_reranked_from_v14": sum(
                bool(entity.get("ranking_changed_from_v14"))
                for entity in output.get("entities") or []
            ),
        }
    )
    output.setdefault("invariants", {}).update(
        {
            "cross_occurrence_candidate_names_used": False,
            "cross_occurrence_ground_truth_used": False,
            "cross_occurrence_expanded_candidates_can_replace_winner": False,
            "cross_occurrence_candidate_scores_changed": False,
            "cross_occurrence_decision_scope_changed": False,
            "cross_occurrence_authentic_standard_confirmation_created": False,
        }
    )
    metadata["entity_semantic_sha256"] = entity_semantic_sha256(output)
    return output, changed


__all__ = [
    "ALGORITHM_ID",
    "MODULE_ID",
    "apply_consensus_to_payload",
    "build_consensus_plan",
    "candidate_formula",
    "cluster_events",
    "connectivity_key",
    "entity_semantic_sha256",
    "load_policy",
    "project_entity",
    "score_event_direct_neutral",
    "semantic_sha256",
]
