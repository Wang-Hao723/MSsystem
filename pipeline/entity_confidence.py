"""Preserve measured single-energy High evidence at an unchanged identity scope.

This layer cannot change a candidate, score, shortlist, chemistry veto or
standard-confirmation decision. Its retrospective assessment and row-level
error ledger are in research/20260913_single_energy_high_validation.
"""
from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any

from matching.flagship_ranker import connectivity_identity
from pipeline.compound_role_projection import spectral_contaminant_veto

ALGORITHM_ID = "single-energy-scope-preserving-existing-high-20260914-v2-report-role-separation"
_ALLOWED_BLOCKERS = {
    "insufficient_distinct_collision_energies",
    "competitive_same_formula_connectivity",
    "competitive_same_formula_connectivities",
}


def _high_scope(match: Mapping[str, Any]) -> str | None:
    scope = match.get("identification_scope")
    isomer = scope == "isomer_group" or (not scope and match.get("isomer_ambiguous"))
    if isomer and match.get("isomer_group_confidence") == "High":
        return "isomer_group"
    if not isomer and (match.get("two_dimensional_connectivity_confidence") or match.get("confidence")) == "High":
        return "2d_connectivity"
    return None


def _connectivity_key(identity: str) -> str:
    namespace, _, value = str(identity or "").partition(":")
    return value if namespace in {"std2d", "ik14"} and len(value) == 14 else identity


def single_energy_high_evidence(
    entity: Mapping[str, Any], rows: Sequence[Mapping[str, Any]],
) -> dict[str, Any] | None:
    """Return auditable support; absent evidence always leaves the old tier intact."""
    decision = entity.get("decision") or {}
    winner = entity.get("winner") or {}
    if (
        not winner
        or winner.get("vetoed")
        or (winner.get("chemistry_gate") or {}).get("precursor_transformation_ambiguity")
        or entity.get("distinct_collision_energy_bins") != 1
        or decision.get("scope") not in {"single_collision_energy_candidate", "molecular_formula_isomer_group"}
        or set(decision.get("blockers") or []) - _ALLOWED_BLOCKERS
    ):
        return None
    # Existing multi-energy High rules are preserved separately by the caller.
    if (entity.get("identity_resolution") or {}).get("upstream_decision_status") == "selected":
        return None
    ids = {c["identity_key"] for c in entity.get("candidates") or [] if not c.get("vetoed")}
    if winner.get("identity_key") not in ids:
        return None
    scope = "isomer_group" if decision.get("scope") == "molecular_formula_isomer_group" else "2d_connectivity"
    supports = []
    for index in entity.get("member_result_indices") or []:
        if not isinstance(index, int) or not 0 <= index < len(rows):
            raise ValueError("Entity High evidence refers to a missing source result row")
        row = rows[index]
        match = row.get("match") or {}
        raw_scope = _high_scope(match)
        if not raw_scope or spectral_contaminant_veto(match) or any(match.get(k) for k in (
            "name_only_lead", "legacy_unvalidated_expert_influence_detected",
        )):
            continue
        identity = connectivity_identity(match)[0]
        if identity != winner["identity_key"]:
            continue
        if scope == "2d_connectivity" and (len(ids) != 1 or raw_scope != "2d_connectivity"):
            continue
        if raw_scope == "isomer_group":
            raw_ids = {identity}
            raw_ids.update(connectivity_identity(c)[0] for c in match.get("isomer_set") or [] if c.get("spectrally_competitive"))
            if scope != "isomer_group" or not {_connectivity_key(v) for v in raw_ids}.issubset({_connectivity_key(v) for v in ids}):
                continue
        supports.append({
            "result_index": index, "spectrum_id": row.get("spectrum_id"),
            "identity_key": identity, "raw_high_scope": raw_scope,
            "lib_spectrum_id": match.get("lib_spectrum_id"),
            "dataset_source": match.get("dataset_source"),
            "cosine": match.get("cosine"), "mass_ppm": match.get("mass_ppm"),
            "direct_matched_product_ions": match.get("direct_matched_product_ions"),
        })
        if match.get("potential_contaminant") and not spectral_contaminant_veto(match):
            supports[-1]["report_role_excluded_from_plant_table_only"] = True
    if not supports:
        return None
    return {
        "algorithm_id": ALGORITHM_ID, "tier": "High", "scope": scope,
        "entity_id": entity.get("entity_id"), "supporting_spectra": supports,
        "source": "existing_measured_ms2_high_for_same_selected_connectivity",
        "candidate_set_narrowed": False, "ranking_changed": False,
        "score_changed": False, "standard_confirmed": False,
        "confidence_is_probability": False, "confidence_is_project_qvalue": False,
    }
