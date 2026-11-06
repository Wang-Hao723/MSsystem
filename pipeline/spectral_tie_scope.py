"""Do not report a complete identity group from a truncated tied shortlist."""
from __future__ import annotations

import copy
import hashlib
from collections import Counter
from pathlib import Path

ALGORITHM_ID = "mssystem-v16-incomplete-spectral-tie-scope-1.0.0"
LOADED_IMPLEMENTATION_SHA256 = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()


def implementation_sha256():
    if hashlib.sha256(Path(__file__).read_bytes()).hexdigest() != LOADED_IMPLEMENTATION_SHA256:
        raise ValueError("Spectral tie scope code changed; reload the reviewed implementation")
    return LOADED_IMPLEMENTATION_SHA256


def _ranking_snapshot(payload):
    return [(entity.get("entity_id"), (entity.get("winner") or {}).get("identity_key"),
             [(c.get("identity_key"), c.get("entity_rank"), c.get("entity_score"))
              for c in entity.get("candidates") or []])
            for entity in payload.get("entities") or []]


def apply_spectral_tie_scope_guard(payload):
    if (payload.get("policy", {}).get("spectral_tie_handling") or {}).get("enabled") is not True:
        return payload
    digest = implementation_sha256()
    before = _ranking_snapshot(payload)
    output = copy.deepcopy(payload)
    incomplete = 0
    for entity in output.get("entities") or []:
        decision = entity.get("decision") or {}
        evidence = decision.get("score_tie_handling") or {}
        outside = evidence.get("leading_entity_tie_outside_display") or []
        truncated = evidence.get("truncated_top_spectral_ties") or []
        if not (outside or truncated):
            continue
        incomplete += 1
        previous_guard = entity.get("spectral_tie_release_guard") or {}
        if previous_guard.get("algorithm_id") != ALGORITHM_ID:
            entity["spectral_tie_release_guard"] = {
                "algorithm_id": ALGORITHM_ID, "status": "incomplete_spectral_tie_group",
                "upstream_decision": copy.deepcopy(decision),
                "upstream_identity_resolution": copy.deepcopy(entity.get("identity_resolution")),
                "tie_evidence": copy.deepcopy(evidence),
                "displayed_candidate_set_complete": False,
                "ranking_changed": False, "scores_changed": False,
                "standard_confirmed": False,
            }
        decision.update({
            "status": "unresolved_formula_event", "scope": "unresolved_precursor_event",
            "reportable_identity_key": None, "reportable_formula": None,
            "ranking_lead_is_reportable_single_connectivity": False,
            "blockers": list(dict.fromkeys([*(decision.get("blockers") or []),
                                             "incomplete_spectral_tie_group"])),
        })
        entity["decision"] = decision
        resolution = entity.setdefault("identity_resolution", {})
        resolution.update({
            "status": "unresolved_incomplete_spectral_tie_group",
            "reportable_scope": "unresolved_precursor_event", "reportable_formula": None,
            "displayed_candidate_set_complete": False,
            "spectral_tie_scope_guard_algorithm_id": ALGORITHM_ID,
        })
        selected = entity.get("winner") or entity.get("leading_candidate_when_abstained")
        if selected:
            selected["reporting_semantics"] = "ranking_lead_not_reportable_compound_identity"
    assert _ranking_snapshot(output) == before, "Spectral tie reporting guard changed ranking"
    metadata = output.setdefault("metadata", {})
    metadata.update({
        "spectral_tie_scope_guard_algorithm_id": ALGORITHM_ID,
        "spectral_tie_scope_guard_source_sha256": digest,
        "entities_with_incomplete_spectral_tie_groups": incomplete,
        "decision_status_counts": dict(Counter((e.get("decision") or {}).get("status")
                                                for e in output.get("entities") or [])),
        "decision_scope_counts": dict(Counter((e.get("decision") or {}).get("scope")
                                               for e in output.get("entities") or [])),
        "identity_resolution_status_counts": dict(Counter(
            (e.get("identity_resolution") or {}).get("status")
            for e in output.get("entities") or [])),
        "unresolved_same_formula_isomer_entity_count": sum(
            (e.get("identity_resolution") or {}).get("status") == "unresolved_same_formula_isomer_group"
            for e in output.get("entities") or []),
        "unresolved_cross_formula_entity_count": sum(
            (e.get("identity_resolution") or {}).get("status") == "unresolved_cross_formula_event"
            for e in output.get("entities") or []),
        "unresolved_incomplete_formula_entity_count": sum(
            (e.get("identity_resolution") or {}).get("status") == "unresolved_incomplete_formula_event"
            for e in output.get("entities") or []),
    })
    output.setdefault("invariants", {})["truncated_tied_candidate_set_can_claim_unique_or_group_high"] = False
    implementation_sha256()
    return output
