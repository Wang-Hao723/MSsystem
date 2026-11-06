"""Build candidate-pair ambiguity graphs from persisted identification rows."""

from __future__ import annotations

import hashlib
import json
from collections.abc import Iterable, Mapping
from typing import Any

from schemas.acquisition import (
    AcquisitionNeed,
    AmbiguityCandidate,
    AmbiguityEdge,
)


def _digest(prefix: str, payload: Any, length: int = 16) -> str:
    encoded = json.dumps(
        payload, ensure_ascii=False, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")
    return f"{prefix}_{hashlib.sha256(encoded).hexdigest()[:length]}"


def candidate_identity_key(candidate: Mapping[str, Any]) -> str:
    inchikey = str(candidate.get("inchikey") or "").strip().upper()
    if inchikey:
        return f"ik::{inchikey}"
    smiles = str(candidate.get("smiles") or "").strip()
    if smiles:
        return f"sm::{smiles}"
    return "nf::" + "|".join([
        str(candidate.get("compound_name") or candidate.get("name") or "unknown")
        .strip().casefold(),
        str(candidate.get("formula") or "").strip(),
    ])


def _candidate(candidate: Mapping[str, Any]) -> AmbiguityCandidate | None:
    name = str(
        candidate.get("compound_name") or candidate.get("name") or ""
    ).strip()
    if not name:
        return None
    rank_score = candidate.get("rank_score", candidate.get("final_score"))
    try:
        parsed_rank = float(rank_score) if rank_score is not None else None
    except (TypeError, ValueError):
        parsed_rank = None
    return AmbiguityCandidate(
        identity_key=candidate_identity_key(candidate),
        compound_name=name,
        formula=str(candidate.get("formula") or "").strip() or None,
        inchikey=str(candidate.get("inchikey") or "").strip() or None,
        smiles=str(candidate.get("smiles") or "").strip() or None,
        adduct=str(candidate.get("adduct") or "").strip() or None,
        rank_score=parsed_rank,
        source_tier=str(candidate.get("source_tier") or "").strip() or None,
    )


def _candidate_pool(match: Mapping[str, Any]) -> tuple[AmbiguityCandidate, ...]:
    raw: list[Mapping[str, Any]] = [match]
    raw.extend(
        item for item in match.get("top_candidates") or []
        if isinstance(item, Mapping)
    )
    passport = match.get("evidence_governance") or {}
    raw.extend(
        item for item in passport.get("ambiguity_set") or []
        if isinstance(item, Mapping)
    )
    candidates: list[AmbiguityCandidate] = []
    seen: set[str] = set()
    selected_formula = str(match.get("formula") or "").strip()
    for item in raw:
        parsed = _candidate(item)
        if parsed is None or parsed.identity_key in seen:
            continue
        # The current V10 ambiguity planner targets same-formula/near-isobaric
        # competition.  Formula-missing candidates remain eligible because a
        # missing formula is itself unresolved evidence, not proof of mismatch.
        if selected_formula and parsed.formula and parsed.formula != selected_formula:
            continue
        seen.add(parsed.identity_key)
        candidates.append(parsed)
        if len(candidates) >= 6:
            break
    return tuple(candidates)


def _entity_id(task_id: str, row: Mapping[str, Any], result_index: int) -> str:
    entity_token = (
        row.get("neutral_analyte_id")
        or row.get("lcms_feature_id")
        or row.get("precursor_feature_id")
        or row.get("ion_identity_network_id")
        or row.get("spectrum_id")
        or f"result:{result_index}"
    )
    # Deliberately excludes compound name, candidate identity and score.
    return _digest("entity", [task_id, str(entity_token)])


def _missing_evidence(row: Mapping[str, Any], match: Mapping[str, Any]) -> tuple[str, ...]:
    missing = ["independent_injection_replicate"]
    if not bool(match.get("rt_confirmed")):
        missing.append("same_method_retention_time_or_authentic_standard")
    ce = (
        match.get("ce_profile_evidence")
        or match.get("collision_energy_evidence")
        or {}
    )
    if not ce or str(ce.get("status") or "") in {"", "unavailable", "not_run"}:
        missing.append("energy_resolved_fragmentation")
    deconvolution = match.get("chimeric_deconvolution") or {}
    if str(deconvolution.get("isolation_status") or "") in {"mixed", "severe"}:
        missing.append("narrow_isolation_confirmation")
    isotope = match.get("formula_isotope_evidence") or {}
    if str(isotope.get("status") or "") != "supporting":
        missing.append("candidate_independent_isotope_support")
    if not row.get("lcms_feature_id"):
        missing.append("continuous_ms1_feature_link")
    return tuple(dict.fromkeys(missing))


def build_ambiguity_edges(
    candidates: tuple[AmbiguityCandidate, ...],
) -> tuple[AmbiguityEdge, ...]:
    edges: list[AmbiguityEdge] = []
    for left_index, left in enumerate(candidates):
        for right in candidates[left_index + 1 :]:
            reasons = []
            if left.formula and right.formula and left.formula == right.formula:
                reasons.append("same_formula_connectivity_competition")
            if left.inchikey and right.inchikey:
                if left.inchikey[:14] != right.inchikey[:14]:
                    reasons.append("different_2d_connectivity")
                elif left.inchikey != right.inchikey:
                    reasons.append("stereochemical_or_protonation_variant")
            if left.rank_score is not None and right.rank_score is not None:
                if abs(left.rank_score - right.rank_score) <= 0.10:
                    reasons.append("near_tie_rank_score")
            if not reasons:
                reasons.append("unresolved_candidate_competition")
            edge_id = _digest(
                "edge", sorted([left.identity_key, right.identity_key])
            )
            edges.append(AmbiguityEdge(
                edge_id=edge_id,
                left=left,
                right=right,
                reasons=tuple(reasons),
            ))
    return tuple(edges)


def build_acquisition_needs(
    task_id: str,
    results: Iterable[Mapping[str, Any]],
    *,
    max_needs: int = 25,
) -> list[AcquisitionNeed]:
    """Select one highest-ranked ambiguous observation per physical entity."""

    representatives: dict[str, tuple[float, AcquisitionNeed]] = {}
    for result_index, row in enumerate(results):
        if not isinstance(row, Mapping):
            continue
        match = row.get("match")
        if not isinstance(match, Mapping):
            continue
        candidates = _candidate_pool(match)
        if len(candidates) < 2:
            continue
        edges = build_ambiguity_edges(candidates)
        if not edges:
            continue
        try:
            precursor_mz = float(row.get("precursor_mz"))
        except (TypeError, ValueError):
            continue
        if precursor_mz <= 0:
            continue
        retention_time = row.get("retention_time")
        try:
            rt = float(retention_time) if retention_time is not None else None
        except (TypeError, ValueError):
            rt = None
        entity_id = _entity_id(task_id, row, result_index)
        need_id = _digest("need", [task_id, entity_id, result_index])
        selected = candidates[0]
        need = AcquisitionNeed(
            need_id=need_id,
            task_id=task_id,
            entity_id=entity_id,
            source_result_index=result_index,
            spectrum_id=str(row.get("spectrum_id") or "").strip() or None,
            precursor_mz=precursor_mz,
            retention_time_min=rt,
            polarity=str(row.get("ion_mode") or row.get("polarity") or "").strip() or None,
            candidate_independent_linkage_fields=(
                "polarity", "precursor_mz", "retention_time_min"
            ),
            candidates=candidates,
            ambiguity_edges=edges,
            missing_evidence=_missing_evidence(row, match),
            baseline_snapshot={
                "selected_identity_key": selected.identity_key,
                "selected_compound_name": selected.compound_name,
                "selected_rank_score": selected.rank_score,
                "candidate_margin": match.get("candidate_margin"),
                "legacy_confidence": match.get("confidence"),
                "claim_status": (
                    (match.get("evidence_governance") or {}).get("claim_status")
                ),
            },
        )
        score = float(selected.rank_score or 0)
        previous = representatives.get(entity_id)
        if previous is None or score > previous[0]:
            representatives[entity_id] = (score, need)
    needs = [item[1] for item in representatives.values()]
    needs.sort(key=lambda need: (
        -len(need.ambiguity_edges),
        -float(need.candidates[0].rank_score or 0),
        need.retention_time_min if need.retention_time_min is not None else 1e9,
        need.precursor_mz,
    ))
    return needs[:max_needs]
