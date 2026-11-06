"""Lossless scientific adapter from the V7/V9 legacy result dictionaries."""

from __future__ import annotations

import hashlib
import math
from collections.abc import Mapping
from typing import Any

from .decision import (
    CalibrationStatus,
    CalibrationVector,
    ClaimDecision,
    ClaimScope,
)
from .evidence import (
    EvidenceApplicability,
    EvidenceNode,
    EvidenceProvenance,
    EvidenceStatus,
    ModuleState,
)
from .identification import (
    CandidateIdentity,
    IdentificationRecord,
    ObservationReference,
    ScoreBundle,
)


_STATUS_MAP = {item.value: item for item in EvidenceStatus}
_APPLICABILITY_MAP = {item.value: item for item in EvidenceApplicability}
_SCOPE_MAP = {
    "authentic_standard": ClaimScope.AUTHENTIC_STANDARD,
    "standard_confirmed": ClaimScope.AUTHENTIC_STANDARD,
    "2d_connectivity": ClaimScope.CONNECTIVITY,
    "connectivity": ClaimScope.CONNECTIVITY,
    "connectivity_level_library_candidate": ClaimScope.CONNECTIVITY,
    "molecular_formula_isomer_group": ClaimScope.ISOMER_GROUP,
    "isomer_group": ClaimScope.ISOMER_GROUP,
    "molecular_formula": ClaimScope.MOLECULAR_FORMULA,
    "formula": ClaimScope.MOLECULAR_FORMULA,
    "chemical_class": ClaimScope.CHEMICAL_CLASS,
    "class": ClaimScope.CHEMICAL_CLASS,
    "putative_or_class_level": ClaimScope.CHEMICAL_CLASS,
    "analog": ClaimScope.ANALOG,
    "unresolved": ClaimScope.UNRESOLVED,
    "spectral_candidate_only": ClaimScope.UNRESOLVED,
    "abstain": ClaimScope.ABSTAIN,
}


def _finite(value: Any) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def _string(value: Any) -> str | None:
    text = str(value or "").strip()
    return text or None


def _connectivity_key(inchikey: str | None) -> str | None:
    if not inchikey:
        return None
    first = inchikey.upper().split("-")[0]
    return first if len(first) == 14 and first.isalpha() else None


def _identity_key(match: Mapping[str, Any]) -> str:
    inchikey = _string(match.get("inchikey"))
    if inchikey:
        return f"inchikey:{inchikey.upper()}"
    formula = _string(match.get("formula")) or "unknown_formula"
    name = _string(match.get("compound_name")) or "unknown_compound"
    digest = hashlib.sha256(f"{formula}|{name}".encode("utf-8")).hexdigest()[:20]
    return f"legacy:{digest}"


def _status(value: Any) -> EvidenceStatus:
    return _STATUS_MAP.get(str(value or "").strip().lower(), EvidenceStatus.UNKNOWN)


def _applicability(value: Any) -> EvidenceApplicability:
    text = str(value or "").strip().lower()
    return _APPLICABILITY_MAP.get(text, EvidenceApplicability.NOT_APPLICABLE)


def _scope(value: Any, *, isomer_ambiguous: bool) -> ClaimScope:
    text = str(value or "").strip().lower()
    if text in _SCOPE_MAP:
        return _SCOPE_MAP[text]
    return ClaimScope.ISOMER_GROUP if isomer_ambiguous else ClaimScope.UNRESOLVED


def _legacy_evidence_nodes(
    task_id: str,
    row: Mapping[str, Any],
    match: Mapping[str, Any],
    runtime_snapshot_id: str | None,
) -> list[EvidenceNode]:
    passport = match.get("evidence_governance")
    passport = passport if isinstance(passport, Mapping) else {}
    nodes: list[EvidenceNode] = []
    spectrum_id = _string(row.get("spectrum_id")) or "unknown"
    for index, raw in enumerate(passport.get("nodes") or []):
        if not isinstance(raw, Mapping):
            continue
        node_id = _string(raw.get("key")) or f"legacy_node_{index + 1}"
        used_for = tuple(str(item) for item in raw.get("used_for") or [])
        raw_provenance = raw.get("provenance")
        raw_provenance = raw_provenance if isinstance(raw_provenance, Mapping) else {}
        provenance = EvidenceProvenance(
            source_type="legacy_evigate_node",
            source_id=node_id,
            task_id=task_id,
            spectrum_id=spectrum_id,
            library_spectrum_id=_string(
                raw_provenance.get("library_spectrum_id")
                or match.get("lib_spectrum_id")
            ),
            library_accession=_string(
                raw_provenance.get("library_accession")
                or match.get("lib_source_accession")
            ),
            method_id=_string(match.get("query_rt_method_id")),
            instrument_family=_string(
                row.get("instrument_type") or match.get("query_instrument_type")
            ),
            runtime_snapshot_id=runtime_snapshot_id,
            **{
                str(key): value
                for key, value in raw_provenance.items()
                if key
                not in {
                    "library_spectrum_id",
                    "library_accession",
                    "query_spectrum_id",
                }
            },
        )
        nodes.append(
            EvidenceNode(
                node_id=node_id,
                label=_string(raw.get("label")) or node_id,
                status=_status(raw.get("status")),
                direction=str(raw.get("direction") or "neutral"),
                applicability=_applicability(raw.get("applicability")),
                reason=str(raw.get("reason") or "Legacy result did not store a reason."),
                value=raw.get("value"),
                provenance=provenance,
                independence_group=str(raw.get("independence_group") or node_id),
                independent_measurement=bool(raw.get("independent_measurement")),
                used_for=used_for,
                warnings=tuple(str(item) for item in raw.get("warnings") or []),
                runtime_influence=bool(
                    {"ranking", "confidence"}.intersection(used_for)
                ),
            )
        )

    isopair = match.get("isopair_evidence")
    isopair = isopair if isinstance(isopair, Mapping) else {}
    adduct = isopair.get("adduct_context")
    adduct = adduct if isinstance(adduct, Mapping) else {}
    if bool(adduct.get("formula_conditioned")):
        nodes.append(
            EvidenceNode(
                node_id="recap_formula_conditioned_adduct",
                label="Formula-conditioned adduct hypothesis",
                status=(
                    EvidenceStatus.INFORMATIONAL
                    if adduct.get("status") == "resolved"
                    else EvidenceStatus.INSUFFICIENT
                ),
                direction="research_counterfactual_only",
                applicability=(
                    EvidenceApplicability.APPLICABLE
                    if adduct.get("status") == "resolved"
                    else EvidenceApplicability.NOT_APPLICABLE
                ),
                reason=(
                    "The provisional formula and precursor mass define an auditable, "
                    "connectivity-neutral hypothesis, but real MIX33 validation forbids "
                    "production rank, confidence or claim influence."
                ),
                value={
                    "status": adduct.get("status"),
                    "resolved_adduct": adduct.get("resolved_adduct"),
                    "mass_error_da": adduct.get("mass_error_da"),
                    "mass_error_ppm": adduct.get("mass_error_ppm"),
                    "production_abstention_reason": isopair.get(
                        "production_abstention_reason"
                    ),
                },
                provenance=EvidenceProvenance(
                    source_type="recap_formula_conditioned_adduct",
                    source_id="matching.recap_adduct",
                    task_id=task_id,
                    spectrum_id=spectrum_id,
                    runtime_snapshot_id=runtime_snapshot_id,
                ),
                independence_group="candidate_conditioned_formula_context",
                independent_measurement=False,
                used_for=("shadow_audit",),
                warnings=("not_independent_query_evidence",),
                runtime_influence=False,
                candidate_conditioned=True,
            )
        )
    return nodes


def _module_states(match: Mapping[str, Any]) -> dict[str, ModuleState]:
    states = {"evidence_governance": ModuleState.PRODUCTION}
    isopair = match.get("isopair_evidence")
    if isinstance(isopair, Mapping):
        adduct = isopair.get("adduct_context")
        if isinstance(adduct, Mapping) and bool(adduct.get("formula_conditioned")):
            states["recap_formula_adduct"] = ModuleState.SHADOW
        if bool(isopair.get("evaluated")):
            states["isopair_ms"] = (
                ModuleState.PRODUCTION
                if bool(isopair.get("runtime_influence"))
                else ModuleState.SHADOW
            )
        else:
            states["isopair_ms"] = ModuleState.QUALIFIED
    if isinstance(match.get("contrastive_expert_shadow"), Mapping):
        states["contrastive_expert"] = ModuleState.SHADOW
    if isinstance(match.get("ce_profile_evidence"), Mapping):
        states["collision_energy_profile"] = ModuleState.SHADOW
    return states


def legacy_row_to_identification(
    task_id: str,
    row: Mapping[str, Any],
    *,
    row_index: int = 0,
    runtime_snapshot_id: str | None = None,
) -> IdentificationRecord | None:
    match = row.get("match") if isinstance(row, Mapping) else None
    if not isinstance(match, Mapping):
        return None

    spectrum_id = _string(row.get("spectrum_id")) or f"index:{row_index}"
    identity_key = _identity_key(match)
    record_hash = hashlib.sha256(
        f"{task_id}|{spectrum_id}|{identity_key}".encode("utf-8")
    ).hexdigest()[:24]
    inchikey = _string(match.get("inchikey"))
    passport = match.get("evidence_governance")
    passport = passport if isinstance(passport, Mapping) else {}
    isopair = match.get("isopair_evidence")
    isopair = isopair if isinstance(isopair, Mapping) else {}
    ranking_changed = bool(isopair.get("ranking_changed"))
    calibration = CalibrationVector(
        status=CalibrationStatus.UNAVAILABLE,
        reason="project_level_target_decoy_qvalue_not_yet_available",
    )
    ambiguity_keys = tuple(
        str(item.get("identity_key"))
        for item in passport.get("ambiguity_set") or []
        if isinstance(item, Mapping) and item.get("identity_key")
    )
    decision = ClaimDecision(
        policy_id=str(passport.get("policy_id") or "legacy-uncompiled"),
        status=str(passport.get("claim_status") or "legacy_uncompiled"),
        scope=_scope(
            passport.get("claim_scope") or match.get("identification_scope"),
            isomer_ambiguous=bool(match.get("isomer_ambiguous")),
        ),
        released=False,
        eligible_before_replication=bool(
            passport.get("release_eligible_before_replication")
        ),
        requires_independent_replication=True,
        reason_codes=tuple(str(item) for item in passport.get("release_blockers") or []),
        warnings=tuple(str(item) for item in passport.get("warnings") or []),
        independent_support_groups=tuple(
            str(item) for item in passport.get("independence_groups") or []
        ),
        independent_injection_count=0,
        ambiguity_identity_keys=ambiguity_keys,
        coverage_guarantee=bool(passport.get("coverage_guarantee")),
        ranking_changed=ranking_changed,
        confidence_field_overwritten=bool(
            passport.get("confidence_field_overwritten")
        ),
        calibration=calibration,
    )
    return IdentificationRecord(
        record_id=f"idrec:{record_hash}",
        observation=ObservationReference(
            task_id=task_id,
            spectrum_id=spectrum_id,
            precursor_event_id=_string(row.get("precursor_feature_id")),
            chromatographic_feature_id=_string(row.get("lcms_feature_id")),
            ion_identity_group_id=_string(
                row.get("ion_identity_network_id") or row.get("lcms_adduct_group_id")
            ),
            neutral_analyte_id=None,
            injection_group_id=None,
            polarity=_string(row.get("ion_mode")),
            precursor_mz=_finite(row.get("precursor_mz")),
            retention_time_minutes=_finite(row.get("retention_time")),
        ),
        candidate=CandidateIdentity(
            identity_key=identity_key,
            compound_name=_string(match.get("compound_name")) or "Unresolved compound",
            inchikey=inchikey,
            connectivity_key=_connectivity_key(inchikey),
            formula=_string(match.get("formula")),
            smiles=_string(match.get("smiles")),
            adduct=_string(match.get("adduct")),
            source_tier=_string(match.get("source_tier")),
            library_spectrum_id=_string(match.get("lib_spectrum_id")),
            library_accession=_string(match.get("lib_source_accession")),
        ),
        scores=ScoreBundle(
            selected_rank=1,
            cosine=_finite(match.get("cosine")),
            entropy=_finite(match.get("entropy")),
            reverse_match=_finite(match.get("reverse_match")),
            mass_ppm=_finite(match.get("mass_ppm")),
            final_score=_finite(match.get("final_score")),
            rank_score=_finite(match.get("rank_score")),
            raw_rank_score=_finite(match.get("raw_rank_score")),
            candidate_margin=_finite(match.get("candidate_margin")),
            legacy_confidence=_string(
                match.get("confidence") or match.get("confidence_class")
            ),
            ranking_source=str(
                match.get("candidate_comparison_source")
                or "legacy_production_search"
            ),
            ranking_changed=ranking_changed,
        ),
        evidence=tuple(
            _legacy_evidence_nodes(
                task_id,
                row,
                match,
                runtime_snapshot_id,
            )
        ),
        decision=decision,
        module_states=_module_states(match),
        provenance={
            "source_schema": "legacy_v7_v9_result_row",
            "source_row_index": row_index,
            "runtime_snapshot_id": runtime_snapshot_id,
            "evidence_policy_id": passport.get("policy_id"),
            "expected_compounds_counted": False,
        },
    )
