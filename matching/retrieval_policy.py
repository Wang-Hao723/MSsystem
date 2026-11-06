"""Auditable conditional retrieval and library-tier policy."""

from __future__ import annotations

from typing import Any


RETRIEVAL_POLICY_ID = "conditional-retrieval-v1"
SEARCHABLE_TIERS = (
    "curated",
    "experimental",
    "massbank_supplement",
    "merged",
    "legacy_unresolved",
)
EXCLUDED_TIERS = (
    "predicted",
    "in_silico",
    "literature_presence_only",
    "user_unverified",
    "metadata_conflict",
    "quarantined",
    "unknown",
)


def build_retrieval_audit(
    *,
    direct_candidates: int,
    neutral_candidates: int,
    candidates_after_source_gate: int,
    candidates_after_polarity_gate: int,
    candidates_scored: int,
    reliable_specific_adduct: bool,
    query_polarity: str | None,
    candidate_budget: int,
    mmap_index_used: bool,
) -> dict[str, Any]:
    return {
        "schema_version": "1.0.0",
        "policy_id": RETRIEVAL_POLICY_ID,
        "retrieval_mode": (
            "direct_precursor_plus_neutral_mass"
            if reliable_specific_adduct and neutral_candidates
            else "direct_precursor_only"
        ),
        "direct_precursor_candidates": int(direct_candidates),
        "neutral_mass_candidates": int(neutral_candidates),
        "candidates_after_source_gate": int(candidates_after_source_gate),
        "candidates_after_polarity_gate": int(candidates_after_polarity_gate),
        "candidates_scored": int(candidates_scored),
        "candidate_budget": int(candidate_budget),
        "dense_window_prescreen_applied": candidates_after_polarity_gate > candidate_budget,
        "source_gate_before_candidate_budget": True,
        "polarity_gate_before_candidate_budget": True,
        "reliable_specific_adduct_required_for_neutral_mass": True,
        "query_specific_adduct_reliable": bool(reliable_specific_adduct),
        "query_polarity": query_polarity,
        "searchable_source_tiers": SEARCHABLE_TIERS,
        "excluded_source_tiers": EXCLUDED_TIERS,
        "mmap_mass_index_used": bool(mmap_index_used),
        "ranking_changed": False,
    }
