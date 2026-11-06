"""EviTwin-MS: authentic-standard candidate ranking for LC-QTOF data.

EviTwin-MS is an additive ranking lane built for the common situation where a
public-library candidate has already been recalled but the wrong isomer is
listed first.  It searches only an explicitly installed authentic-standard
library, compares collision-energy matched spectra, and keeps the formal claim
scope separate from the display ranking.

The implementation does not edit the primary V15 result or legacy audit
artifacts.  It emits a source-bound sidecar that a report or UI can place ahead
of public-library candidates.  Truth labels, expected compound names, and
formula lists are not accepted by the inference API.
"""

from __future__ import annotations

import hashlib
import json
import math
import os
import statistics
import tempfile
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Mapping, Sequence

from matching.authentic_anchor import (
    _mass_filtered_reference_ids,
    evaluate_authentic_anchor,
    evaluate_authentic_anchor_actual_ce_shadow,
    validate_anchor_library,
)
from matching.evitwin_coelution import (
    EVITWIN_COELUTION_ALGORITHM_ID,
    evaluate_coeluting_reference_pairs,
    validate_coelution_result,
)
from matching.library_quality import normalize_formula_text
from matching.rule_quality import normalize_adduct


EVITWIN_SCHEMA_VERSION = "evitwin-ms/1.4"
EVITWIN_MODULE_ID = "evitwin_ms"
EVITWIN_ALGORITHM_ID = "mssystem-v15-evitwin-ms-authentic-ranker-1.4.0"


@dataclass(frozen=True)
class EviTwinPolicy:
    """Frozen operational thresholds for the authentic-reference lane."""

    precursor_tolerance_ppm: float = 10.0
    minimum_collision_energy_count: int = 2
    minimum_collision_energy_votes: int = 2
    minimum_priority_similarity: float = 0.95
    minimum_singleton_similarity: float = 0.98
    minimum_connectivity_margin: float = 0.05
    review_similarity: float = 0.50
    maximum_candidates_per_feature: int = 10
    minimum_reference_rt_order_separation_min: float = 0.05
    minimum_observed_rt_order_separation_min: float = 0.15
    allow_display_rank_change: bool = True
    allow_formal_level1_release: bool = False

    def validated(self) -> "EviTwinPolicy":
        if not 0 < float(self.precursor_tolerance_ppm) <= 100:
            raise ValueError("precursor_tolerance_ppm must be in (0, 100]")
        if int(self.minimum_collision_energy_count) < 2:
            raise ValueError("at least two collision-energy conditions are required")
        if int(self.minimum_collision_energy_votes) < 2:
            raise ValueError("at least two collision-energy votes are required")
        if not 0 <= float(self.review_similarity) <= float(
            self.minimum_priority_similarity
        ) <= 1:
            raise ValueError("review/priority similarity thresholds are invalid")
        if not float(self.minimum_priority_similarity) <= float(
            self.minimum_singleton_similarity
        ) <= 1:
            raise ValueError("singleton similarity must be at least the priority threshold")
        if not 0 <= float(self.minimum_connectivity_margin) <= 1:
            raise ValueError("minimum_connectivity_margin must be in [0, 1]")
        if int(self.maximum_candidates_per_feature) < 2:
            raise ValueError("maximum_candidates_per_feature must be at least two")
        if not 0 < float(self.minimum_reference_rt_order_separation_min) <= 2:
            raise ValueError("reference RT-order separation must be in (0, 2] min")
        if not 0 < float(self.minimum_observed_rt_order_separation_min) <= 2:
            raise ValueError("observed RT-order separation must be in (0, 2] min")
        if not self.allow_display_rank_change:
            raise ValueError("EviTwin-MS must be allowed to change its own display ranking")
        if self.allow_formal_level1_release:
            raise ValueError(
                "formal Level-1 release requires a separate same-method RT contract"
            )
        return self


def _finite(value: Any) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def _full_inchikey(value: Any) -> str:
    token = str(value or "").strip().upper().replace("INCHIKEY=", "")
    parts = token.split("-")
    if (
        len(parts) == 3
        and len(parts[0]) == 14
        and len(parts[1]) == 10
        and len(parts[2]) == 1
        and all(part.isalpha() for part in parts)
    ):
        return token
    return ""


def _connectivity(value: Any) -> str:
    full = _full_inchikey(value)
    if full:
        return full[:14]
    token = str(value or "").strip().upper().split("-", 1)[0]
    return token if len(token) == 14 and token.isalpha() else ""


def _candidate_payload(value: Any) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        return {}
    nested = value.get("candidate")
    return nested if isinstance(nested, Mapping) else value


def _candidate_connectivity(value: Any) -> str:
    """Return only a structure-derived 2D key from a V15 candidate record."""

    if not isinstance(value, Mapping):
        return ""
    candidate = _candidate_payload(value)
    for raw in (
        candidate.get("standardized_connectivity_inchikey"),
        candidate.get("inchikey"),
        value.get("connectivity_identity"),
    ):
        if connectivity := _connectivity(raw):
            return connectivity
    identity = str(value.get("identity_key") or "").strip().upper()
    if identity.startswith("STD2D:"):
        return _connectivity(identity.split(":", 1)[1])
    return ""


def _candidate_formula(value: Any) -> str:
    return normalize_formula_text(_candidate_payload(value).get("formula"))


def _eligible_candidate_pool(entity: Mapping[str, Any]) -> dict[str, Any]:
    """Bind EviTwin to one formula inside the V15 Top5/Top40 support pool.

    Expanded Top40 records are accepted only when the V15 attachment records
    every demotion-only invariant.  Missing formula or identity-scope evidence
    is an abstention, never an invitation to search the authentic library.
    """

    native = [
        candidate
        for candidate in (entity.get("candidates") or ())
        if isinstance(candidate, Mapping) and not bool(candidate.get("vetoed"))
    ]
    resolution = entity.get("identity_resolution")
    resolution = resolution if isinstance(resolution, Mapping) else {}
    resolution_status = str(resolution.get("status") or "").strip()
    if resolution_status in {
        "unresolved_cross_formula_event",
        "unresolved_incomplete_formula_event",
    }:
        return {
            "status": "abstained",
            "reason_codes": ["v15_identity_scope_is_unresolved_precursor_event"],
            "formula_lock": None,
            "native_candidate_connectivities": [],
            "expanded_candidate_connectivities": [],
            "eligible_candidate_connectivities": [],
            "expanded_support_contract_verified": False,
        }

    native_formulae = {
        formula for candidate in native if (formula := _candidate_formula(candidate))
    }
    missing_native_formula_count = sum(
        not bool(_candidate_formula(candidate)) for candidate in native
    )
    locked_formula = normalize_formula_text(resolution.get("reportable_formula"))
    if not locked_formula:
        locked_formula = next(iter(native_formulae)) if len(native_formulae) == 1 else ""
    if not locked_formula:
        return {
            "status": "abstained",
            "reason_codes": ["single_v15_molecular_formula_not_resolved"],
            "formula_lock": None,
            "native_candidate_connectivities": [],
            "expanded_candidate_connectivities": [],
            "eligible_candidate_connectivities": [],
            "expanded_support_contract_verified": False,
        }
    if missing_native_formula_count:
        return {
            "status": "abstained",
            "reason_codes": ["native_v15_candidate_formula_incomplete"],
            "formula_lock": None,
            "native_candidate_connectivities": [],
            "expanded_candidate_connectivities": [],
            "eligible_candidate_connectivities": [],
            "expanded_support_contract_verified": False,
        }
    if any(formula and formula != locked_formula for formula in native_formulae):
        return {
            "status": "abstained",
            "reason_codes": ["native_v15_candidate_formula_conflict"],
            "formula_lock": None,
            "native_candidate_connectivities": [],
            "expanded_candidate_connectivities": [],
            "eligible_candidate_connectivities": [],
            "expanded_support_contract_verified": False,
        }

    native_connectivities = sorted({
        connectivity
        for candidate in native
        if _candidate_formula(candidate) == locked_formula
        and (connectivity := _candidate_connectivity(candidate))
    })
    support = entity.get("candidate_recall_support")
    support = support if isinstance(support, Mapping) else {}
    expanded_contract = bool(
        support.get("integration_module_id") == "v15_candidate_recall_support"
        and support.get("expanded_candidates_can_replace_winner") is False
        and support.get("expanded_candidates_can_receive_2d_connectivity_decision")
        is False
        and support.get("production_identity_rank_unchanged") is True
    )
    expanded_connectivities: list[str] = []
    if expanded_contract:
        expanded_connectivities = sorted({
            connectivity
            for candidate in (support.get("expanded_candidates") or ())
            if isinstance(candidate, Mapping)
            and candidate.get("can_replace_winner") is False
            and candidate.get("can_receive_2d_connectivity_decision") is False
            and candidate.get("requires_primary_ranking_or_orthogonal_evidence") is True
            and _candidate_formula(candidate) == locked_formula
            and (connectivity := _candidate_connectivity(candidate))
        })
    eligible = sorted(set(native_connectivities) | set(expanded_connectivities))
    return {
        "status": "eligible" if eligible else "abstained",
        "reason_codes": (
            ["same_formula_v15_candidate_pool_bound"]
            if eligible
            else ["no_structured_candidate_in_locked_formula_pool"]
        ),
        "formula_lock": locked_formula,
        "native_candidate_connectivities": native_connectivities,
        "expanded_candidate_connectivities": expanded_connectivities,
        "eligible_candidate_connectivities": eligible,
        "expanded_support_contract_verified": expanded_contract,
    }


def _canonical_adduct(value: Any) -> str:
    normalized = normalize_adduct(value)
    aliases = {
        "[M+FA-H]-": "[M+HCOO]-",
        "[M+HAC-H]-": "[M+CH3COO]-",
    }
    normalized = aliases.get(normalized.upper(), normalized)
    allowed = {
        "[M-H]-",
        "[M+Cl]-",
        "[M+HCOO]-",
        "[M+CH3COO]-",
        "[M-2H]2-",
        "[M+H]+",
        "[M+Na]+",
        "[M+K]+",
        "[M+NH4]+",
        "[M+2H]2+",
    }
    return normalized if normalized in allowed else ""


def _candidate_independent_entity_adduct(
    entity: Mapping[str, Any],
    polarity: str,
    locked_formula: str,
) -> tuple[str, dict[str, Any]]:
    selected = entity.get("winner") or entity.get("leading_candidate_when_abstained")
    selected = selected if isinstance(selected, Mapping) else {}
    gate = selected.get("chemistry_gate")
    gate = gate if isinstance(gate, Mapping) else {}
    source = str(gate.get("adduct_source") or "").strip()
    adduct = _canonical_adduct(gate.get("adduct"))
    formula = normalize_formula_text(gate.get("formula"))
    expected_suffix = "+" if polarity == "POS" else "-"
    status_counts = gate.get("status_counts")
    status_counts = status_counts if isinstance(status_counts, Mapping) else {}
    evaluated_bins = int(gate.get("evaluated_collision_energy_bins") or 0)
    all_bins_consistent = bool(
        evaluated_bins >= 1
        and int(status_counts.get("consistent") or 0) == evaluated_bins
        and sum(int(value or 0) for value in status_counts.values()) == evaluated_bins
    )
    source_is_reliable = bool(
        (
            source == "ion_identity_network_adduct"
            and gate.get("ion_identity_network_evidence_used") is True
            and gate.get("ion_identity_network_evidence_available") is True
        )
        or (
            source == "lcms_feature_adduct"
            and gate.get("lcms_validated_adduct_group_evidence_available") is True
        )
        or source == "reliable_query_adduct"
    )
    independent_checks = [
        check
        for check in (gate.get("candidate_independent_checks") or ())
        if isinstance(check, Mapping)
        and check.get("consistent") is True
        and _canonical_adduct(check.get("adduct")) == adduct
        and str(check.get("source") or "").strip() == source
    ]
    resolved = bool(
        gate.get("status") == "consistent"
        and gate.get("candidate_conditioned") is False
        and gate.get("veto") is False
        and gate.get("precursor_transformation_ambiguity") is False
        and gate.get(
            "consistent_metadata_overrides_isolated_conflicting_library_metadata"
        )
        is False
        and source_is_reliable
        and independent_checks
        and all_bins_consistent
        and formula
        and formula == normalize_formula_text(locked_formula)
        and adduct
        and adduct.endswith(expected_suffix)
    )
    reasons = []
    if resolved:
        reasons.append("v15_candidate_independent_entity_adduct_passed_all_gates")
    else:
        if not gate:
            reasons.append("v15_entity_chemistry_gate_missing")
        if gate and gate.get("candidate_conditioned") is not False:
            reasons.append("v15_entity_adduct_is_candidate_conditioned")
        if gate and not source_is_reliable:
            reasons.append("v15_entity_adduct_source_not_independently_validated")
        if gate and not all_bins_consistent:
            reasons.append("v15_entity_adduct_not_consistent_across_evaluated_bins")
        if gate and formula != normalize_formula_text(locked_formula):
            reasons.append("v15_entity_adduct_formula_conflicts_with_formula_lock")
        if gate and (
            gate.get("veto") is not False
            or gate.get("precursor_transformation_ambiguity") is not False
        ):
            reasons.append("v15_entity_adduct_gate_is_vetoed_or_ambiguous")
        if gate and not independent_checks:
            reasons.append("v15_entity_adduct_lacks_matching_independent_check")
        if gate and (not adduct or not adduct.endswith(expected_suffix)):
            reasons.append("v15_entity_adduct_missing_or_polarity_mismatched")
    return adduct if resolved else "", {
        "status": "resolved" if resolved else "abstained",
        "resolved_adduct": adduct if resolved else None,
        "adduct_source": source or None,
        "formula": formula or None,
        "locked_formula": normalize_formula_text(locked_formula) or None,
        "evaluated_collision_energy_bins": evaluated_bins,
        "status_counts": dict(status_counts),
        "candidate_conditioned": gate.get("candidate_conditioned"),
        "ion_identity_network_evidence_used": gate.get(
            "ion_identity_network_evidence_used"
        ),
        "lcms_validated_adduct_group_evidence_available": gate.get(
            "lcms_validated_adduct_group_evidence_available"
        ),
        "reason_codes": reasons or ["v15_entity_adduct_gate_not_resolved"],
        "adduct_inferred_from_polarity": False,
        "compound_name_used": False,
        "candidate_identity_used": False,
    }


def _resolved_query_adduct(
    queries: Sequence[Mapping[str, Any]],
    polarity: str,
    *,
    entity: Mapping[str, Any] | None = None,
    locked_formula: str = "",
) -> tuple[str, dict[str, Any]]:
    observed: list[str] = []
    missing_or_unreliable = 0
    for row in queries:
        match = row.get("match")
        match = match if isinstance(match, Mapping) else {}
        reliable = bool(
            row.get("query_adduct_reliable") is True
            or match.get("query_adduct_reliable") is True
        )
        adduct = _canonical_adduct(
            row.get("query_adduct") or match.get("query_adduct")
        )
        if not reliable or not adduct:
            missing_or_unreliable += 1
            continue
        observed.append(adduct)
    unique = sorted(set(observed))
    expected_suffix = "+" if polarity == "POS" else "-"
    explicit_resolved = (
        unique[0]
        if len(unique) == 1
        and missing_or_unreliable == 0
        and unique[0].endswith(expected_suffix)
        else ""
    )
    entity_adduct, entity_audit = _candidate_independent_entity_adduct(
        entity or {},
        polarity,
        locked_formula,
    )
    if explicit_resolved and entity_adduct and explicit_resolved != entity_adduct:
        resolved = ""
        source = None
        reasons = ["explicit_and_entity_candidate_independent_adducts_conflict"]
    elif explicit_resolved:
        resolved = explicit_resolved
        source = "query_spectra_explicit_metadata"
        reasons = ["all_query_spectra_share_one_reliable_adduct"]
    elif entity_adduct and len(unique) <= 1 and (
        not unique or unique[0] == entity_adduct
    ):
        resolved = entity_adduct
        source = "v15_candidate_independent_chemistry_gate"
        reasons = ["candidate_independent_entity_adduct_used_for_missing_query_metadata"]
    elif len(unique) > 1:
        resolved = ""
        source = None
        reasons = ["query_spectra_have_conflicting_reliable_adducts"]
    elif unique and entity_adduct and unique[0] != entity_adduct:
        resolved = ""
        source = None
        reasons = ["partial_explicit_and_entity_candidate_independent_adducts_conflict"]
    elif missing_or_unreliable:
        resolved = ""
        source = None
        reasons = ["query_adduct_missing_or_unreliable"]
    else:
        resolved = ""
        source = None
        reasons = ["query_adduct_polarity_mismatch"]
    return resolved, {
        "status": "resolved" if resolved else "abstained",
        "resolved_adduct": resolved or None,
        "resolution_source": source,
        "observed_reliable_adducts": unique,
        "missing_or_unreliable_query_spectrum_count": missing_or_unreliable,
        "reason_codes": reasons,
        "entity_candidate_independent_adduct_audit": entity_audit,
        "adduct_inferred_from_polarity": False,
    }


def _reference_adduct(entry: Mapping[str, Any]) -> str:
    qc = entry.get("reference_qc")
    qc = qc if isinstance(qc, Mapping) else {}
    selected = _canonical_adduct(qc.get("selected_adduct"))
    spectra = {
        adduct
        for spectrum in (entry.get("spectra") or ())
        if isinstance(spectrum, Mapping)
        and (adduct := _canonical_adduct(spectrum.get("adduct")))
    }
    if selected and (not spectra or spectra == {selected}):
        return selected
    if not selected and len(spectra) == 1:
        return next(iter(spectra))
    return ""


def _winner_summary(entity: Mapping[str, Any] | None) -> dict[str, Any] | None:
    if not isinstance(entity, Mapping):
        return None
    winner = entity.get("winner") or entity.get("leading_candidate_when_abstained")
    if not isinstance(winner, Mapping):
        return None
    candidate = _candidate_payload(winner)
    return {
        "compound_name": candidate.get("compound_name"),
        "compound_name_cn": candidate.get("compound_name_cn"),
        "formula": candidate.get("formula"),
        "inchikey": candidate.get("inchikey"),
        "connectivity": _connectivity(
            candidate.get("standardized_connectivity_inchikey")
            or candidate.get("inchikey")
        ),
        "score": _finite(
            winner.get("evirank_mv_score")
            if winner.get("evirank_mv_score") is not None
            else winner.get("entity_score")
        ),
    }


def _canonical_sha256(value: Any) -> str:
    return hashlib.sha256(
        json.dumps(
            value,
            ensure_ascii=False,
            sort_keys=True,
            separators=(",", ":"),
        ).encode("utf-8")
    ).hexdigest()


def _query_rows(
    entity: Mapping[str, Any], result_rows: Sequence[Mapping[str, Any]]
) -> list[Mapping[str, Any]]:
    selected = []
    seen: set[int] = set()
    for raw_index in entity.get("member_result_indices") or ():
        try:
            index = int(raw_index)
        except (TypeError, ValueError):
            continue
        if index in seen or not 0 <= index < len(result_rows):
            continue
        seen.add(index)
        row = result_rows[index]
        if not isinstance(row, Mapping) or row.get("blank_excluded"):
            continue
        peaks = (
            row.get("query_peaks")
            or row.get("raw_peaks")
            or row.get("peaks")
            or ()
        )
        if len(peaks) < 3:
            continue
        selected.append(row)
    return selected


def _compatible_reference_ids(
    library: Mapping[str, Any],
    *,
    precursor_mz: float,
    polarity: str,
    excluded_source_task_id: str,
    tolerance_ppm: float,
    locked_formula: str,
    eligible_connectivities: Sequence[str],
    query_adduct: str,
) -> tuple[list[str], dict[str, Any]]:
    reference_ids, mass_audit = _mass_filtered_reference_ids(
        library,
        precursor_mz=precursor_mz,
        polarity=polarity,
        tolerance_ppm=tolerance_ppm,
    )
    by_id = {
        str(entry.get("reference_id") or ""): entry
        for entry in library.get("entries") or ()
        if isinstance(entry, Mapping)
    }
    mass_by_id = {
        str(row.get("reference_id") or ""): row
        for row in mass_audit
        if isinstance(row, Mapping)
    }
    allowed = set(eligible_connectivities)
    selected = []
    exclusion_counts = {
        "self_reference": 0,
        "candidate_pool": 0,
        "formula": 0,
        "reference_adduct": 0,
        "mass_hypothesis_adduct": 0,
    }
    for reference_id in reference_ids:
        entry = by_id.get(reference_id) or {}
        if str(entry.get("source_task_id") or "") == excluded_source_task_id:
            exclusion_counts["self_reference"] += 1
            continue
        if _connectivity(entry.get("inchikey")) not in allowed:
            exclusion_counts["candidate_pool"] += 1
            continue
        if normalize_formula_text(entry.get("formula")) != locked_formula:
            exclusion_counts["formula"] += 1
            continue
        if _reference_adduct(entry) != query_adduct:
            exclusion_counts["reference_adduct"] += 1
            continue
        hypothesis = _canonical_adduct(
            (mass_by_id.get(reference_id) or {}).get("adduct_hypothesis")
        )
        if hypothesis != query_adduct:
            exclusion_counts["mass_hypothesis_adduct"] += 1
            continue
        selected.append(reference_id)
    return selected, {
        "mass_compatible_reference_count": len(reference_ids),
        "eligible_reference_count": len(selected),
        "exclusion_counts": exclusion_counts,
        "formula_lock": locked_formula,
        "query_adduct": query_adduct,
        "precursor_tolerance_ppm": float(tolerance_ppm),
    }


def _same_2d_candidate_set(scores: Sequence[Mapping[str, Any]]) -> bool:
    full = {
        key
        for score in scores
        if (key := _full_inchikey(score.get("inchikey")))
    }
    connectivities = {
        key
        for score in scores
        if (key := _connectivity(score.get("inchikey")))
    }
    return len(full) >= 2 and len(connectivities) < len(full)


def _feature_decision(
    anchor: Mapping[str, Any],
    *,
    policy: EviTwinPolicy,
    compatible_anchor_count: int,
) -> dict[str, Any]:
    scores = [
        score
        for score in (anchor.get("candidate_scores") or ())
        if isinstance(score, Mapping)
    ]
    winner = scores[0] if scores else None
    similarity = _finite((winner or {}).get("median_similarity")) or 0.0
    ce_count = int((winner or {}).get("collision_energy_count") or 0)
    votes = int(anchor.get("winner_collision_energy_votes") or 0)
    margin = _finite(anchor.get("winner_margin")) or 0.0
    priority_threshold = (
        policy.minimum_singleton_similarity
        if compatible_anchor_count == 1
        else policy.minimum_priority_similarity
    )
    priority_evidence = bool(
        winner
        and similarity >= priority_threshold
        and ce_count >= policy.minimum_collision_energy_count
        and votes >= policy.minimum_collision_energy_votes
    )
    margin_qualified = bool(
        compatible_anchor_count == 1
        or margin >= policy.minimum_connectivity_margin
    )
    same_2d = _same_2d_candidate_set(scores)
    if priority_evidence and margin_qualified and not same_2d:
        status = "qualified_authentic_rank"
        display_tier = "priority"
        scope = "authentic_reference_2d_candidate"
        reasons = ["multi_ce_authentic_spectrum_priority_gate_passed"]
    elif priority_evidence:
        status = "unresolved_authentic_candidate_group"
        display_tier = "review"
        scope = "isomer_group"
        reasons = [
            "same_2d_stereoisomer_set"
            if same_2d
            else "connectivity_margin_below_frozen_threshold"
        ]
    elif winner and similarity >= policy.review_similarity:
        status = "review_candidate"
        display_tier = "review"
        scope = "candidate_only"
        reasons = ["authentic_similarity_below_priority_gate"]
    else:
        status = "abstained"
        display_tier = "hidden"
        scope = "none"
        reasons = ["authentic_evidence_below_review_gate"]
    return {
        "status": status,
        "display_tier": display_tier,
        "reportable_scope": scope,
        "reason_codes": reasons,
        "priority_evidence": priority_evidence,
        "margin_qualified": margin_qualified,
        "same_2d_stereoisomer_set": same_2d,
        "formal_level1_released": False,
    }


def _reference_retention_time(entry: Mapping[str, Any]) -> float | None:
    values = [
        value
        for spectrum in (entry.get("spectra") or ())
        if isinstance(spectrum, Mapping)
        and (value := _finite(spectrum.get("retention_time_min"))) is not None
    ]
    return statistics.median(values) if values else None


def _canonical_method_id(value: Any) -> str:
    token = " ".join(str(value or "").strip().casefold().split())
    return "" if token in {"", "unspecified", "unknown", "not_available"} else token


def _resolved_query_method(
    declared_method_id: Any,
    queries: Sequence[Mapping[str, Any]],
) -> tuple[str, dict[str, Any]]:
    declared = _canonical_method_id(declared_method_id)
    row_methods = [
        method
        for row in queries
        if (method := _canonical_method_id(row.get("rt_method_id")))
    ]
    unique = sorted(set(row_methods))
    missing_rows = len(queries) - len(row_methods)
    if declared and all(method == declared for method in unique):
        resolved = declared
        reasons = ["task_lc_method_explicit_and_query_rows_nonconflicting"]
    elif not declared and len(unique) == 1 and missing_rows == 0:
        resolved = unique[0]
        reasons = ["all_query_rows_share_one_explicit_lc_method"]
    elif declared and any(method != declared for method in unique):
        resolved = ""
        reasons = ["task_and_query_row_lc_methods_conflict"]
    elif len(unique) > 1:
        resolved = ""
        reasons = ["query_rows_have_conflicting_lc_methods"]
    else:
        resolved = ""
        reasons = ["query_lc_method_missing"]
    return resolved, {
        "status": "resolved" if resolved else "not_applicable",
        "declared_query_method_id": declared or None,
        "observed_query_row_method_ids": unique,
        "query_rows_missing_method_id": missing_rows,
        "resolved_query_method_id": resolved or None,
        "reason_codes": reasons,
    }


def _winner_from_score(score: Mapping[str, Any], library: Mapping[str, Any]) -> dict[str, Any]:
    entry = next(
        (
            item
            for item in library.get("entries") or ()
            if item.get("reference_id") == score.get("reference_id")
        ),
        {},
    )
    return {
        "reference_id": score.get("reference_id"),
        "compound_code": score.get("compound_code"),
        "compound_name": score.get("compound_name"),
        "formula": entry.get("formula"),
        "inchikey": score.get("inchikey"),
        "connectivity": _connectivity(score.get("inchikey")),
        "median_similarity": _finite(score.get("median_similarity")),
        "collision_energy_count": int(score.get("collision_energy_count") or 0),
        "collision_energy_votes": int(score.get("collision_energy_votes") or 0),
        "winner_margin": _finite(score.get("winner_margin")),
    }


def _apply_retention_order_assignment(
    features: list[dict[str, Any]],
    library: Mapping[str, Any],
    policy: EviTwinPolicy,
) -> int:
    """Apply a one-to-one monotonic RT assignment inside an isomer family.

    Absolute RTs are deliberately not compared across the single-standard and
    mixture gradients.  Only the order is used, and only when every competing
    authentic candidate and the same number of separated query features are
    observed.  This prevents two closely related peaks from both being named
    as whichever authentic spectrum wins by a tiny cosine margin.
    """

    for feature in features:
        query_method_id = str(
            (feature.get("query_method_audit") or {}).get(
                "resolved_query_method_id"
            )
            or ""
        )
        feature["retention_order_assignment"] = {
            "status": "not_applied",
            "reason_codes": [
                "same_lc_method_not_proven"
                if not query_method_id
                else "retention_order_gates_not_all_passed"
            ],
            "query_method_id": query_method_id or None,
            "absolute_rt_used_for_cross_method_scoring": False,
            "formal_rt_confirmation": False,
        }

    entries = {
        str(entry.get("reference_id") or ""): entry
        for entry in library.get("entries") or ()
        if isinstance(entry, Mapping)
    }
    grouped: dict[tuple[str, tuple[str, ...]], list[dict[str, Any]]] = {}
    for feature in features:
        query_method_id = str(
            (feature.get("query_method_audit") or {}).get(
                "resolved_query_method_id"
            )
            or ""
        )
        if not query_method_id:
            continue
        key = tuple(sorted(str(value) for value in feature.get(
            "mass_compatible_reference_ids", []
        ) if str(value)))
        if len(key) < 2:
            continue
        grouped.setdefault((query_method_id, key), []).append(feature)
    applied = 0
    for (query_method_id, reference_ids), group in grouped.items():
        if len(reference_ids) > policy.maximum_candidates_per_feature:
            continue
        reference_methods = {
            _canonical_method_id((entries.get(reference_id) or {}).get("method_id"))
            for reference_id in reference_ids
        }
        if reference_methods != {query_method_id}:
            for feature in group:
                feature["retention_order_assignment"]["reason_codes"] = [
                    "reference_lc_method_missing_or_mismatched"
                ]
                feature["retention_order_assignment"][
                    "reference_method_ids"
                ] = sorted(method for method in reference_methods if method)
            continue
        reference_rows = []
        for reference_id in reference_ids:
            rt = _reference_retention_time(entries.get(reference_id) or {})
            if rt is None:
                reference_rows = []
                break
            reference_rows.append((rt, reference_id))
        if len(reference_rows) != len(reference_ids):
            continue
        reference_rows.sort()
        if any(
            reference_rows[index + 1][0] - reference_rows[index][0]
            < policy.minimum_reference_rt_order_separation_min
            for index in range(len(reference_rows) - 1)
        ):
            continue
        # Select the strongest distinct chromatographic features before using
        # their order.  Repeated DDA entities within one RT region cannot
        # masquerade as multiple separated compounds.
        candidates = sorted(
            (
                feature
                for feature in group
                if _finite(feature.get("representative_retention_time_min"))
                is not None
                and max(
                    (
                        float(score.get("median_similarity") or 0.0)
                        for score in feature.get("candidate_scores") or ()
                    ),
                    default=0.0,
                )
                >= policy.review_similarity
            ),
            key=lambda feature: (
                -max(
                    float(score.get("median_similarity") or 0.0)
                    for score in feature.get("candidate_scores") or ()
                ),
                float(feature.get("representative_retention_time_min") or 0.0),
            ),
        )
        selected: list[dict[str, Any]] = []
        for feature in candidates:
            rt = float(feature["representative_retention_time_min"])
            if any(
                abs(rt - float(item["representative_retention_time_min"]))
                < policy.minimum_observed_rt_order_separation_min
                for item in selected
            ):
                continue
            selected.append(feature)
            if len(selected) == len(reference_rows):
                break
        if len(selected) != len(reference_rows):
            continue
        selected.sort(key=lambda feature: float(feature[
            "representative_retention_time_min"
        ]))
        if any(
            float(selected[index + 1]["representative_retention_time_min"])
            - float(selected[index]["representative_retention_time_min"])
            < policy.minimum_observed_rt_order_separation_min
            for index in range(len(selected) - 1)
        ):
            continue
        assignments = []
        valid = True
        for order_index, (feature, (_, reference_id)) in enumerate(
            zip(selected, reference_rows), start=1
        ):
            score = next(
                (
                    item
                    for item in feature.get("candidate_scores") or ()
                    if str(item.get("reference_id") or "") == reference_id
                ),
                None,
            )
            if (
                score is None
                or float(score.get("median_similarity") or 0.0)
                < policy.minimum_priority_similarity
                or int(score.get("collision_energy_count") or 0)
                < policy.minimum_collision_energy_count
                or int(score.get("collision_energy_votes") or 0)
                < policy.minimum_collision_energy_votes
            ):
                valid = False
                break
            assignments.append((feature, score, order_index))
        if not valid:
            continue
        assignment_id = hashlib.sha256(
            "|".join(
                f"{feature.get('feature_id')}={score.get('reference_id')}"
                for feature, score, _ in assignments
            ).encode("utf-8")
        ).hexdigest()[:16]
        for feature, score, order_index in assignments:
            feature["evitwin_winner"] = _winner_from_score(score, library)
            feature["decision"] = {
                "status": "qualified_retention_order_rank",
                "display_tier": "priority",
                "reportable_scope": "isomer_group",
                "reason_codes": [
                    "multi_ce_authentic_spectrum_and_monotonic_rt_order_passed"
                ],
                "priority_evidence": True,
                "margin_qualified": False,
                "same_2d_stereoisomer_set": _same_2d_candidate_set(
                    feature.get("candidate_scores") or ()
                ),
                "formal_level1_released": False,
            }
            feature["retention_order_assignment"] = {
                "status": "applied",
                "assignment_id": assignment_id,
                "assignment_size": len(assignments),
                "observed_order_index": order_index,
                "reference_order_index": order_index,
                "query_method_id": query_method_id,
                "reference_method_id": query_method_id,
                "same_lc_method_verified": True,
                "absolute_rt_used_for_cross_method_scoring": False,
                "formal_rt_confirmation": False,
            }
            applied += 1
    return applied


def build_evitwin_result(
    primary_payload: Mapping[str, Any],
    result_rows: Sequence[Mapping[str, Any]],
    library: Mapping[str, Any],
    *,
    task_id: str,
    polarity: str,
    query_instrument_family: str,
    query_method_id: str = "",
    evirank_payload: Mapping[str, Any] | None = None,
    policy: EviTwinPolicy | None = None,
) -> dict[str, Any]:
    """Build a truth-free authentic-reference ranking sidecar for one task."""

    validate_anchor_library(library)
    selected_policy = (policy or EviTwinPolicy()).validated()
    task_token = str(task_id or "").strip()
    if not task_token:
        raise ValueError("task_id is required")
    mode = str(polarity or "").strip().upper()
    if mode not in {"POS", "NEG"}:
        raise ValueError("polarity must be POS or NEG")
    instrument = str(query_instrument_family or "").strip()
    if not instrument:
        raise ValueError("query_instrument_family is required")
    source_entities = [
        entity
        for entity in (primary_payload.get("entities") or ())
        if isinstance(entity, Mapping)
    ]
    challenger_by_id = {
        str(entity.get("entity_id") or ""): entity
        for entity in (evirank_payload or {}).get("entities") or ()
        if isinstance(entity, Mapping)
    }
    features = []
    actual_ce_shadow_features = []
    applicability_audits = []
    mass_compatible_entities = 0
    candidate_pool_eligible_entities = 0
    reliable_adduct_entities = 0
    adduct_resolution_source_counts: dict[str, int] = {}
    reference_eligible_entities = 0
    self_reference_exclusion_count = 0
    for entity in source_entities:
        entity_id = str(entity.get("entity_id") or "")
        precursor = _finite(entity.get("representative_precursor_mz"))
        applicability = {
            "feature_id": entity_id or None,
            "status": "abstained",
            "reason_codes": [],
        }
        if not entity_id or precursor is None or precursor <= 0:
            applicability["reason_codes"] = [
                "entity_id_or_representative_precursor_missing"
            ]
            applicability_audits.append(applicability)
            continue
        all_ids, _ = _mass_filtered_reference_ids(
            library,
            precursor_mz=precursor,
            polarity=mode,
            tolerance_ppm=selected_policy.precursor_tolerance_ppm,
        )
        applicability["mass_compatible_reference_count"] = len(all_ids)
        mass_compatible_entities += int(bool(all_ids))
        queries = _query_rows(entity, result_rows)
        if not queries:
            applicability["reason_codes"] = ["no_usable_query_msms_spectrum"]
            applicability_audits.append(applicability)
            continue
        candidate_pool = _eligible_candidate_pool(entity)
        applicability["candidate_pool_audit"] = candidate_pool
        if candidate_pool["status"] != "eligible":
            applicability["reason_codes"] = list(candidate_pool["reason_codes"])
            applicability_audits.append(applicability)
            continue
        candidate_pool_eligible_entities += 1
        query_adduct, adduct_audit = _resolved_query_adduct(
            queries,
            mode,
            entity=entity,
            locked_formula=str(candidate_pool["formula_lock"]),
        )
        applicability["query_adduct_audit"] = adduct_audit
        if not query_adduct:
            applicability["reason_codes"] = list(adduct_audit["reason_codes"])
            applicability_audits.append(applicability)
            continue
        reliable_adduct_entities += 1
        adduct_resolution_source = str(
            adduct_audit.get("resolution_source") or "unknown"
        )
        adduct_resolution_source_counts[adduct_resolution_source] = (
            adduct_resolution_source_counts.get(adduct_resolution_source, 0) + 1
        )
        resolved_method, method_audit = _resolved_query_method(
            query_method_id,
            queries,
        )
        reference_ids, reference_filter_audit = _compatible_reference_ids(
            library,
            precursor_mz=precursor,
            polarity=mode,
            excluded_source_task_id=task_token,
            tolerance_ppm=selected_policy.precursor_tolerance_ppm,
            locked_formula=str(candidate_pool["formula_lock"]),
            eligible_connectivities=candidate_pool[
                "eligible_candidate_connectivities"
            ],
            query_adduct=query_adduct,
        )
        self_reference_exclusion_count += int(
            reference_filter_audit["exclusion_counts"]["self_reference"]
        )
        applicability["reference_filter_audit"] = reference_filter_audit
        applicability["query_method_audit"] = method_audit
        if not reference_ids:
            applicability["reason_codes"] = [
                "no_reference_passed_candidate_formula_adduct_and_self_exclusion_gates"
            ]
            applicability_audits.append(applicability)
            continue
        reference_eligible_entities += 1
        public_entity = challenger_by_id.get(entity_id) or entity
        public_winner = _winner_summary(public_entity)
        actual_ce_shadow = evaluate_authentic_anchor_actual_ce_shadow(
            queries,
            library,
            candidate_identity_keys=reference_ids,
            query_polarity=mode,
            query_instrument_family=instrument,
        )
        actual_ce_vote_map = dict(
            actual_ce_shadow.get("collision_energy_votes") or {}
        )
        actual_ce_scores = [
            {
                **dict(score),
                "collision_energy_votes": int(
                    actual_ce_vote_map.get(
                        str(score.get("reference_id") or ""), 0
                    )
                ),
                "winner_margin": _finite(actual_ce_shadow.get("winner_margin")),
            }
            for score in (actual_ce_shadow.get("candidate_scores") or ())
            if isinstance(score, Mapping)
        ][: selected_policy.maximum_candidates_per_feature]
        actual_ce_comparable = len(actual_ce_scores) >= 2
        if actual_ce_comparable:
            actual_ce_winner = actual_ce_scores[0]
            actual_ce_public_connectivity = str(
                (public_winner or {}).get("connectivity") or ""
            )
            actual_ce_winner_connectivity = _connectivity(
                actual_ce_winner.get("inchikey")
            )
            actual_ce_shadow_features.append({
                "feature_id": entity_id,
                "representative_precursor_mz": precursor,
                "representative_retention_time_min": _finite(
                    entity.get("representative_retention_time")
                ),
                "query_spectrum_count": len(queries),
                "collision_energy_conditions": sorted({
                    round(value, 3)
                    for row in queries
                    if (
                        value := _finite(
                            row.get("collision_energy_ev")
                            if row.get("collision_energy_ev") is not None
                            else row.get("collision_energy")
                        )
                    )
                    is not None
                }),
                "candidate_pool_audit": candidate_pool,
                "query_adduct_audit": adduct_audit,
                "reference_filter_audit": reference_filter_audit,
                "public_library_winner": public_winner,
                "shadow_winner": _winner_from_score(actual_ce_winner, library),
                "candidate_scores": actual_ce_scores,
                "decision": {
                    "status": actual_ce_shadow.get("decision"),
                    "reason_codes": list(
                        actual_ce_shadow.get("reason_codes") or ()
                    ),
                    "reportable_scope": "unresolved_candidate_group",
                    "display_rank_change_allowed": False,
                    "production_rank_change_allowed": False,
                    "formal_level1_released": False,
                },
                "shadow_lead_differs_from_public_library": bool(
                    actual_ce_winner_connectivity
                    and actual_ce_winner_connectivity
                    != actual_ce_public_connectivity
                ),
                "scalar_similarity_interpolation_used": bool(
                    actual_ce_shadow.get("scalar_similarity_interpolation_used")
                ),
                "outside_reference_range_clamped": False,
                "legacy_result_mutated": False,
            })
        anchor = evaluate_authentic_anchor(
            queries,
            library,
            candidate_identity_keys=reference_ids,
            query_polarity=mode,
            query_instrument_family=instrument,
            baseline_separated=None,
        )
        vote_map = dict(anchor.get("collision_energy_votes") or {})
        scores = [
            {
                **dict(score),
                "collision_energy_votes": int(
                    vote_map.get(str(score.get("reference_id") or ""), 0)
                ),
                "winner_margin": _finite(anchor.get("winner_margin")),
            }
            for score in (anchor.get("candidate_scores") or ())
            if isinstance(score, Mapping)
            and int(score.get("collision_energy_count") or 0) > 0
        ][: selected_policy.maximum_candidates_per_feature]
        if not scores:
            applicability["status"] = (
                "actual_ce_shadow_scored" if actual_ce_comparable else "abstained"
            )
            applicability["reason_codes"] = [
                "exact_ce_ranker_abstained_actual_ce_shadow_available"
                if actual_ce_comparable
                else (
                    "actual_ce_shadow_has_fewer_than_two_ce_range_compatible_candidates"
                    if actual_ce_scores
                    else "no_same_instrument_authentic_reference_score"
                )
            ]
            applicability_audits.append(applicability)
            continue
        anchor = {**anchor, "candidate_scores": scores}
        decision = _feature_decision(
            anchor,
            policy=selected_policy,
            compatible_anchor_count=len(reference_ids),
        )
        coelution_support = evaluate_coeluting_reference_pairs(
            queries,
            library,
            reference_ids,
        )
        winner = scores[0]
        winner_connectivity = _connectivity(winner.get("inchikey"))
        public_connectivity = str((public_winner or {}).get("connectivity") or "")
        display_change = bool(
            decision["display_tier"] == "priority"
            and winner_connectivity
            and winner_connectivity != public_connectivity
        )
        applicability.update({
            "status": "scored",
            "reason_codes": ["all_evitwin_applicability_gates_passed"],
            "scored_reference_count": len(scores),
        })
        applicability_audits.append(applicability)
        features.append(
            {
                "feature_id": entity_id,
                "source_entity_ids": [entity_id],
                "representative_precursor_mz": precursor,
                "representative_retention_time_min": _finite(
                    entity.get("representative_retention_time")
                ),
                "query_spectrum_count": len(queries),
                "collision_energy_conditions": sorted(
                    {
                        round(value, 3)
                        for row in queries
                        if (
                            value := _finite(
                                row.get("collision_energy_ev")
                                if row.get("collision_energy_ev") is not None
                                else row.get("collision_energy")
                            )
                        )
                        is not None
                    }
                ),
                "mass_compatible_anchor_count": len(reference_ids),
                "mass_compatible_reference_ids": list(reference_ids),
                "candidate_pool_audit": candidate_pool,
                "query_adduct_audit": adduct_audit,
                "query_method_audit": {
                    **method_audit,
                    "resolved_query_method_id": resolved_method or None,
                },
                "reference_filter_audit": reference_filter_audit,
                "public_library_winner": public_winner,
                "evitwin_winner": _winner_from_score(winner, library),
                "candidate_scores": scores,
                "decision": decision,
                "coelution_support": coelution_support,
                "display_rank_changed_from_public_library": display_change,
                "legacy_result_mutated": False,
            }
        )

    retention_order_assignment_count = _apply_retention_order_assignment(
        features,
        library,
        selected_policy,
    )
    # The joint RT-order pass can change the selected authentic candidate.
    # Recompute display changes only after that pass has frozen the feature
    # ranking.
    for feature in features:
        winner_connectivity = str(
            (feature.get("evitwin_winner") or {}).get("connectivity") or ""
        )
        public_connectivity = str(
            (feature.get("public_library_winner") or {}).get("connectivity") or ""
        )
        feature["display_rank_changed_from_public_library"] = bool(
            (feature.get("decision") or {}).get("display_tier") == "priority"
            and winner_connectivity
            and winner_connectivity != public_connectivity
        )

    priority = [
        feature
        for feature in features
        if (feature.get("decision") or {}).get("display_tier") == "priority"
    ]
    coeluting_pairs = [
        {
            "feature_id": feature.get("feature_id"),
            **dict(
                ((feature.get("coelution_support") or {}).get("decision") or {}).get(
                    "supported_pair"
                )
                or {}
            ),
            "reportable_scope": "coeluting_pair_presence_support",
            "formal_level1_released": False,
            "exact_top1_changed": False,
        }
        for feature in features
        if ((feature.get("coelution_support") or {}).get("decision") or {}).get(
            "status"
        )
        == "qualified_coeluting_pair_support"
    ]
    # A compound can appear through repeated DDA precursor events.  Keep every
    # feature for audit, while exposing one strongest row per 2D connectivity
    # (or reference id when connectivity is unavailable) to the user.
    best_by_identity: dict[str, dict[str, Any]] = {}
    for feature in priority:
        winner = feature["evitwin_winner"]
        key = str(
            winner.get("connectivity")
            or winner.get("reference_id")
            or feature.get("feature_id")
        )
        previous = best_by_identity.get(key)
        strength = (
            float(winner.get("median_similarity") or 0.0),
            int(winner.get("collision_energy_count") or 0),
            float(winner.get("winner_margin") or 0.0),
            -float(feature.get("representative_retention_time_min") or 0.0),
        )
        if previous is None:
            best_by_identity[key] = feature
            continue
        old = previous["evitwin_winner"]
        old_strength = (
            float(old.get("median_similarity") or 0.0),
            int(old.get("collision_energy_count") or 0),
            float(old.get("winner_margin") or 0.0),
            -float(previous.get("representative_retention_time_min") or 0.0),
        )
        if strength > old_strength:
            best_by_identity[key] = feature
    compounds = sorted(
        best_by_identity.values(),
        key=lambda feature: (
            -float(feature["evitwin_winner"].get("median_similarity") or 0.0),
            str(feature["evitwin_winner"].get("connectivity") or ""),
            str(feature.get("feature_id") or ""),
        ),
    )
    status_counts: dict[str, int] = {}
    for feature in features:
        status = str((feature.get("decision") or {}).get("status") or "unknown")
        status_counts[status] = status_counts.get(status, 0) + 1
    return {
        "metadata": {
            "schema_version": EVITWIN_SCHEMA_VERSION,
            "module_id": EVITWIN_MODULE_ID,
            "algorithm_id": EVITWIN_ALGORITHM_ID,
            "coelution_algorithm_id": EVITWIN_COELUTION_ALGORITHM_ID,
            "algorithm_name": "EviTwin-MS 标品谱联合排序",
            "deployment_state": "active_authentic_reference_display_ranker",
            "task_id": task_token,
            "polarity": mode,
            "query_instrument_family": instrument,
            "declared_query_method_id": _canonical_method_id(query_method_id) or None,
            "source_primary_semantic_sha256": _canonical_sha256(primary_payload),
            "source_primary_module_id": (
                (primary_payload.get("metadata") or {}).get("module_id")
            ),
            "source_primary_algorithm_id": (
                (primary_payload.get("metadata") or {}).get("algorithm_id")
            ),
            "source_primary_is_v15": str(
                (primary_payload.get("metadata") or {}).get("module_id") or ""
            ).startswith("v15"),
            "source_v11_semantic_sha256": (
                _canonical_sha256(primary_payload)
                if str((primary_payload.get("metadata") or {}).get("module_id") or "")
                .startswith("v11")
                else None
            ),
            "source_results_semantic_sha256": _canonical_sha256(list(result_rows)),
            "source_evirank_semantic_sha256": (
                _canonical_sha256(evirank_payload)
                if isinstance(evirank_payload, Mapping)
                else None
            ),
            "authentic_library_content_sha256": library.get("content_sha256"),
            "policy": asdict(selected_policy),
            "truth_available_at_inference": False,
            "expected_name_available_at_inference": False,
            "actual_ce_shadow_scoring_mode": (
                "recorded_query_ce_scalar_similarity_interpolation"
            ),
        },
        "summary": {
            "source_entity_count": len(source_entities),
            "mass_compatible_entity_count": mass_compatible_entities,
            "candidate_pool_eligible_entity_count": candidate_pool_eligible_entities,
            "reliable_adduct_entity_count": reliable_adduct_entities,
            "adduct_resolution_source_counts": dict(
                sorted(adduct_resolution_source_counts.items())
            ),
            "reference_eligible_entity_count": reference_eligible_entities,
            "applicability_abstention_count": sum(
                audit.get("status") != "scored" for audit in applicability_audits
            ),
            "scored_feature_count": len(features),
            "priority_feature_count": len(priority),
            "priority_compound_count": len(compounds),
            "display_rank_change_count": sum(
                bool(feature.get("display_rank_changed_from_public_library"))
                for feature in priority
            ),
            "self_reference_exclusion_count": self_reference_exclusion_count,
            "retention_order_assignment_count": retention_order_assignment_count,
            "coelution_evaluated_feature_count": sum(
                int((feature.get("coelution_support") or {}).get(
                    "evaluated_pair_count"
                ) or 0) > 0
                for feature in features
            ),
            "coelution_supported_feature_count": len(coeluting_pairs),
            "coelution_supported_pair_count": len(coeluting_pairs),
            "decision_status_counts": dict(sorted(status_counts.items())),
        },
        "invariants": {
            "truth_labels_used_for_ranking": False,
            "compound_names_used_for_scoring": False,
            "database_accessions_used_for_scoring": False,
            "database_source_order_used_for_scoring": False,
            "authentic_library_used_for_candidate_generation": False,
            "single_standard_self_match_allowed": False,
            "legacy_v11_artifact_mutated": False,
            "legacy_v12_artifact_mutated": False,
            "source_primary_artifact_mutated": False,
            "formal_level1_release_performed": False,
            "same_2d_stereoisomer_exact_claim_allowed": False,
            "coelution_support_can_change_exact_top1": False,
            "coelution_support_is_separated_peak_assignment": False,
            "coelution_support_can_release_level1": False,
            "candidate_pool_can_expand_beyond_v15_top40": False,
            "cross_formula_rank_change_allowed": False,
            "unreliable_adduct_rank_change_allowed": False,
            "cross_method_rt_order_allowed": False,
            "actual_ce_scalar_interpolation_can_change_display_order": False,
            "actual_ce_scalar_interpolation_can_change_production_rank": False,
            "actual_ce_scalar_interpolation_extrapolation_allowed": False,
            "authentic_reference_lane_can_change_display_order": True,
        },
        "actual_ce_shadow": {
            "runtime_mode": "shadow_only",
            "scoring_mode": (
                "recorded_query_ce_scalar_similarity_interpolation"
            ),
            "candidate_generation_used": False,
            "query_collision_energy_reconstructed": False,
            "reference_spectra_interpolated": False,
            "outside_reference_range_clamped": False,
            "display_rank_change_allowed": False,
            "production_rank_change_allowed": False,
            "formal_level1_release_allowed": False,
            "summary": {
                "scored_feature_count": len(actual_ce_shadow_features),
                "single_ce_feature_count": sum(
                    len(feature.get("collision_energy_conditions") or ()) == 1
                    for feature in actual_ce_shadow_features
                ),
                "multi_ce_feature_count": sum(
                    len(feature.get("collision_energy_conditions") or ()) >= 2
                    for feature in actual_ce_shadow_features
                ),
                "shadow_lead_difference_count": sum(
                    bool(feature.get("shadow_lead_differs_from_public_library"))
                    for feature in actual_ce_shadow_features
                ),
                "display_rank_change_count": 0,
                "production_rank_change_count": 0,
            },
            "features": actual_ce_shadow_features,
        },
        "compounds": compounds,
        "coeluting_pairs": coeluting_pairs,
        "features": features,
        "applicability_audits": applicability_audits,
    }


def validate_evitwin_result(payload: Mapping[str, Any]) -> None:
    metadata = payload.get("metadata") or {}
    if metadata.get("schema_version") != EVITWIN_SCHEMA_VERSION:
        raise ValueError("unsupported EviTwin-MS schema")
    if metadata.get("module_id") != EVITWIN_MODULE_ID:
        raise ValueError("unexpected EviTwin-MS module id")
    if metadata.get("algorithm_id") != EVITWIN_ALGORITHM_ID:
        raise ValueError("unexpected EviTwin-MS algorithm id")
    if metadata.get("truth_available_at_inference") is not False:
        raise ValueError("truth must be unavailable at EviTwin-MS inference")
    if len(str(metadata.get("source_primary_semantic_sha256") or "")) != 64:
        raise ValueError("EviTwin-MS primary source binding is missing")
    policy = EviTwinPolicy(**dict(metadata.get("policy") or {})).validated()
    if asdict(policy) != dict(metadata.get("policy") or {}):
        raise ValueError("EviTwin-MS policy must be complete and canonical")
    invariants = payload.get("invariants") or {}
    required_false = (
        "truth_labels_used_for_ranking",
        "compound_names_used_for_scoring",
        "database_accessions_used_for_scoring",
        "database_source_order_used_for_scoring",
        "authentic_library_used_for_candidate_generation",
        "single_standard_self_match_allowed",
        "legacy_v11_artifact_mutated",
        "legacy_v12_artifact_mutated",
        "source_primary_artifact_mutated",
        "formal_level1_release_performed",
        "same_2d_stereoisomer_exact_claim_allowed",
        "coelution_support_can_change_exact_top1",
        "coelution_support_is_separated_peak_assignment",
        "coelution_support_can_release_level1",
        "candidate_pool_can_expand_beyond_v15_top40",
        "cross_formula_rank_change_allowed",
        "unreliable_adduct_rank_change_allowed",
        "cross_method_rt_order_allowed",
        "actual_ce_scalar_interpolation_can_change_display_order",
        "actual_ce_scalar_interpolation_can_change_production_rank",
        "actual_ce_scalar_interpolation_extrapolation_allowed",
    )
    if any(invariants.get(key) is not False for key in required_false):
        raise ValueError("EviTwin-MS safety invariants are incomplete")
    if invariants.get("authentic_reference_lane_can_change_display_order") is not True:
        raise ValueError("EviTwin-MS display-rank contract is missing")
    for feature in payload.get("features") or ():
        pool = feature.get("candidate_pool_audit") or {}
        if pool.get("status") != "eligible" or not pool.get("formula_lock"):
            raise ValueError("every scored EviTwin-MS feature must be formula locked")
        adduct = feature.get("query_adduct_audit") or {}
        if adduct.get("status") != "resolved" or not adduct.get("resolved_adduct"):
            raise ValueError("every scored EviTwin-MS feature requires one reliable adduct")
        allowed = set(pool.get("eligible_candidate_connectivities") or ())
        scored = {
            _connectivity(score.get("inchikey"))
            for score in (feature.get("candidate_scores") or ())
            if isinstance(score, Mapping)
        }
        if not scored or not scored <= allowed:
            raise ValueError("EviTwin-MS scored a structure outside the bound candidate pool")
        decision = feature.get("decision") or {}
        if (
            feature.get("display_rank_changed_from_public_library") is True
            and decision.get("display_tier") != "priority"
        ):
            raise ValueError("review or abstained evidence cannot change display rank")
        if decision.get("status") == "qualified_retention_order_rank":
            rt = feature.get("retention_order_assignment") or {}
            if (
                rt.get("status") != "applied"
                or rt.get("same_lc_method_verified") is not True
                or rt.get("query_method_id") != rt.get("reference_method_id")
            ):
                raise ValueError("retention-order rank lacks an exact LC-method gate")
        coelution = feature.get("coelution_support")
        if not isinstance(coelution, Mapping):
            raise ValueError("every EviTwin-MS feature requires a coelution audit")
        validate_coelution_result(coelution)
    actual_ce_shadow = payload.get("actual_ce_shadow") or {}
    if (
        actual_ce_shadow.get("runtime_mode") != "shadow_only"
        or actual_ce_shadow.get("display_rank_change_allowed") is not False
        or actual_ce_shadow.get("production_rank_change_allowed") is not False
        or actual_ce_shadow.get("formal_level1_release_allowed") is not False
        or actual_ce_shadow.get("outside_reference_range_clamped") is not False
        or actual_ce_shadow.get("query_collision_energy_reconstructed") is not False
        or actual_ce_shadow.get("reference_spectra_interpolated") is not False
    ):
        raise ValueError("actual-CE shadow safety contract is incomplete")
    for feature in actual_ce_shadow.get("features") or ():
        pool = feature.get("candidate_pool_audit") or {}
        adduct = feature.get("query_adduct_audit") or {}
        decision = feature.get("decision") or {}
        if pool.get("status") != "eligible" or not pool.get("formula_lock"):
            raise ValueError("actual-CE shadow feature is not formula locked")
        if adduct.get("status") != "resolved" or not adduct.get("resolved_adduct"):
            raise ValueError("actual-CE shadow feature lacks a reliable adduct")
        if (
            decision.get("display_rank_change_allowed") is not False
            or decision.get("production_rank_change_allowed") is not False
            or decision.get("formal_level1_released") is not False
            or feature.get("outside_reference_range_clamped") is not False
            or feature.get("legacy_result_mutated") is not False
        ):
            raise ValueError("actual-CE shadow feature can change a released rank")
        allowed = set(pool.get("eligible_candidate_connectivities") or ())
        scored = {
            _connectivity(score.get("inchikey"))
            for score in (feature.get("candidate_scores") or ())
            if isinstance(score, Mapping)
        }
        if not scored or not scored <= allowed:
            raise ValueError(
                "actual-CE shadow scored a structure outside the candidate pool"
            )


def write_evitwin_result(path: str | Path, payload: Mapping[str, Any]) -> Path:
    validate_evitwin_result(payload)
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(
        prefix=f".{destination.name}.", suffix=".tmp", dir=destination.parent
    )
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, ensure_ascii=False, indent=2)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, destination)
    except Exception:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass
        raise
    return destination


def read_evitwin_result(path: str | Path) -> dict[str, Any]:
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    validate_evitwin_result(payload)
    return payload


__all__ = [
    "EVITWIN_ALGORITHM_ID",
    "EVITWIN_MODULE_ID",
    "EVITWIN_SCHEMA_VERSION",
    "EviTwinPolicy",
    "build_evitwin_result",
    "read_evitwin_result",
    "validate_evitwin_result",
    "write_evitwin_result",
]
