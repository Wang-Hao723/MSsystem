"""Fit and apply predeclared, project-local calibration without reranking."""

from __future__ import annotations

import hashlib
import json
from collections import Counter, defaultdict
from collections.abc import Iterable, Sequence
from pathlib import Path
from typing import Any

from schemas.decision import CalibrationStatus, CalibrationVector, ClaimScope
from schemas.identification import IdentificationRecord

from .models import (
    CalibrationAssignment,
    CalibrationMethod,
    CalibrationProject,
    CalibrationReadiness,
    CompetitionPair,
    CompetitionRole,
    ProjectCalibrationArtifact,
    ProjectCalibrationStatus,
    StratumCalibration,
    WinnerCall,
)
from .qvalue import compete_pairs, estimate_qvalue_curve, qvalue_for_score


def _canonical(value: Any) -> bytes:
    return json.dumps(
        value, ensure_ascii=False, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")


def load_calibration_method(path: Path) -> CalibrationMethod:
    return CalibrationMethod.model_validate_json(path.read_text(encoding="utf-8"))


def method_semantic_sha256(method: CalibrationMethod) -> str:
    return hashlib.sha256(_canonical(method.model_dump(mode="json"))).hexdigest()


def _field_value(call: WinnerCall, field: str) -> str:
    value = getattr(call, field)
    return str(value.value if hasattr(value, "value") else value)


def audit_project_readiness(
    project: CalibrationProject,
    method: CalibrationMethod,
    pairs: Sequence[CompetitionPair],
) -> CalibrationReadiness:
    reasons: list[str] = []
    warnings: list[str] = []
    winners: tuple[WinnerCall, ...] = ()
    if project.method_id != method.method_id:
        reasons.append("method_id_mismatch")
    if project.score_field != method.score_field:
        reasons.append("score_field_mismatch")
    if not project.method_frozen_before_data:
        reasons.append("method_not_frozen_before_project_data")
    if not project.independent_from_threshold_tuning:
        reasons.append("project_reused_for_threshold_tuning")
    if not project.complete_target_decoy_competition:
        reasons.append("incomplete_target_decoy_competition")
    if project.decoy_identity_overlap_count:
        reasons.append("target_decoy_identity_overlap")
    pool_ratio = project.target_pool_size / project.decoy_pool_size
    if not (
        1.0 / method.max_target_decoy_pool_ratio
        <= pool_ratio
        <= method.max_target_decoy_pool_ratio
    ):
        reasons.append("target_decoy_pool_ratio_out_of_policy")
    try:
        winners = compete_pairs(pairs)
    except ValueError as exc:
        reasons.append("competition_uniqueness_violation")
        warnings.append(str(exc))
    if winners:
        invalid_units = sorted(
            {
                call.unit_level.value
                for call in winners
                if call.unit_level not in method.allowed_unit_levels
            }
        )
        invalid_scopes = sorted(
            {
                call.claim_scope.value
                for call in winners
                if call.claim_scope not in method.allowed_claim_scopes
            }
        )
        if invalid_units:
            reasons.append("unsupported_calibration_unit")
            warnings.append(f"unsupported unit levels: {invalid_units}")
        if invalid_scopes:
            reasons.append("unsupported_claim_scope")
            warnings.append(f"unsupported claim scopes: {invalid_scopes}")
    counts = Counter(call.role for call in winners)
    if len(winners) < method.min_competitions:
        reasons.append("insufficient_competitions")
    if counts[CompetitionRole.TARGET] < method.min_target_winners:
        reasons.append("insufficient_target_winners")
    if counts[CompetitionRole.DECOY] < method.min_decoy_winners:
        reasons.append("insufficient_decoy_winners")
    if any(
        value.strip().lower() in {"", "unknown", "none"}
        for call in winners
        for value in (call.polarity, call.instrument_family, call.source_tier)
    ):
        warnings.append("one_or_more_stratum_values_are_unknown")
    reasons = list(dict.fromkeys(reasons))
    warnings = list(dict.fromkeys(warnings))
    return CalibrationReadiness(
        eligible=not reasons,
        reason_codes=tuple(reasons),
        warnings=tuple(warnings),
        competition_count=len(winners),
        target_winner_count=counts[CompetitionRole.TARGET],
        decoy_winner_count=counts[CompetitionRole.DECOY],
        decoy_scale=pool_ratio,
    )


def _fit_strata(
    project: CalibrationProject,
    method: CalibrationMethod,
    winners: Sequence[WinnerCall],
    decoy_scale: float,
) -> tuple[StratumCalibration, ...]:
    models: list[StratumCalibration] = []
    for fields in method.stratum_fallback:
        grouped: dict[tuple[str, ...], list[WinnerCall]] = defaultdict(list)
        for call in winners:
            grouped[tuple(_field_value(call, field) for field in fields)].append(call)
        for values, calls in sorted(grouped.items()):
            counts = Counter(call.role for call in calls)
            if (
                len(calls) < method.min_competitions
                or counts[CompetitionRole.TARGET] < method.min_target_winners
                or counts[CompetitionRole.DECOY] < method.min_decoy_winners
            ):
                continue
            curve = estimate_qvalue_curve(
                calls,
                decoy_scale=decoy_scale,
                pseudocount=method.pseudocount,
            )
            signature = {
                "project_id": project.project_id,
                "method_id": method.method_id,
                "fields": fields,
                "values": values,
                "competition_ids": sorted(call.competition_id for call in calls),
            }
            model_id = f"qstratum:{hashlib.sha256(_canonical(signature)).hexdigest()[:24]}"
            models.append(
                StratumCalibration(
                    model_id=model_id,
                    fields=fields,
                    values=values,
                    competition_count=len(calls),
                    target_winner_count=counts[CompetitionRole.TARGET],
                    decoy_winner_count=counts[CompetitionRole.DECOY],
                    decoy_scale=decoy_scale,
                    curve=curve,
                )
            )
    return tuple(models)


def calibrate_project(
    project: CalibrationProject,
    method: CalibrationMethod,
    pairs: Iterable[CompetitionPair],
) -> ProjectCalibrationArtifact:
    materialized = tuple(pairs)
    readiness = audit_project_readiness(project, method, materialized)
    method_hash = method_semantic_sha256(method)
    if not readiness.eligible:
        return ProjectCalibrationArtifact(
            status=ProjectCalibrationStatus.UNAVAILABLE,
            method_id=method.method_id,
            method_semantic_sha256=method_hash,
            project=project,
            readiness=readiness,
        )
    winners = compete_pairs(materialized)
    models = _fit_strata(
        project,
        method,
        winners,
        readiness.decoy_scale or 1.0,
    )
    model_lookup = {
        (model.fields, model.values): model
        for model in models
    }
    assignments: list[CalibrationAssignment] = []
    target_winners = [
        call for call in winners if call.role == CompetitionRole.TARGET
    ]
    for call in target_winners:
        selected = None
        for fields in method.stratum_fallback:
            values = tuple(_field_value(call, field) for field in fields)
            selected = model_lookup.get((fields, values))
            if selected is not None:
                break
        if selected is None:
            continue
        assignments.append(
            CalibrationAssignment(
                competition_id=call.competition_id,
                entity_id=call.entity_id,
                target_record_id=call.target_record_id,
                unit_level=call.unit_level,
                claim_scope=call.claim_scope,
                score=call.score,
                q_value=qvalue_for_score(selected.curve, call.score),
                score_field=method.score_field,
                method_id=method.method_id,
                project_id=project.project_id,
                stratum_model_id=selected.model_id,
                stratum_fields=selected.fields,
                stratum_values=selected.values,
            )
        )
    unassigned = len(target_winners) - len(assignments)
    status = (
        ProjectCalibrationStatus.CALIBRATED
        if assignments and unassigned == 0
        else ProjectCalibrationStatus.PARTIAL
        if assignments
        else ProjectCalibrationStatus.UNAVAILABLE
    )
    if status == ProjectCalibrationStatus.UNAVAILABLE:
        readiness = readiness.model_copy(
            update={
                "eligible": False,
                "reason_codes": tuple(
                    [*readiness.reason_codes, "no_predeclared_stratum_met_minimums"]
                ),
            }
        )
    return ProjectCalibrationArtifact(
        status=status,
        method_id=method.method_id,
        method_semantic_sha256=method_hash,
        project=project,
        readiness=readiness,
        strata=models,
        assignments=tuple(assignments),
        target_winners_unassigned=unassigned,
    )


def _calibration_vector(assignment: CalibrationAssignment) -> CalibrationVector:
    q_fields: dict[str, float] = {}
    if assignment.claim_scope in {ClaimScope.AUTHENTIC_STANDARD, ClaimScope.CONNECTIVITY}:
        q_fields["q_connectivity"] = assignment.q_value
    elif assignment.claim_scope == ClaimScope.ISOMER_GROUP:
        q_fields["q_isomer_group"] = assignment.q_value
    elif assignment.claim_scope == ClaimScope.MOLECULAR_FORMULA:
        q_fields["q_formula"] = assignment.q_value
    elif assignment.claim_scope == ClaimScope.CHEMICAL_CLASS:
        q_fields["q_class"] = assignment.q_value
    else:
        raise ValueError(f"Unsupported calibrated claim scope: {assignment.claim_scope}")
    return CalibrationVector(
        status=CalibrationStatus.CALIBRATED,
        method_id=assignment.method_id,
        project_id=assignment.project_id,
        stratum_model_id=assignment.stratum_model_id,
        unit_level=assignment.unit_level.value,
        score_field=assignment.score_field,
        transfer_scope=assignment.transfer_scope,
        **q_fields,
    )


def apply_project_calibration(
    records: Iterable[IdentificationRecord],
    artifact: ProjectCalibrationArtifact,
) -> tuple[IdentificationRecord, ...]:
    """Attach earned q-values while preserving identity, scores and legacy tier."""
    assignments = {
        assignment.target_record_id: assignment
        for assignment in artifact.assignments
    }
    output = []
    for record in records:
        assignment = assignments.get(record.record_id)
        if assignment is None:
            output.append(record)
            continue
        if record.decision.scope != assignment.claim_scope:
            raise ValueError(
                f"Calibration claim scope mismatch for {record.record_id}: "
                f"{record.decision.scope.value} != {assignment.claim_scope.value}"
            )
        record_score = getattr(record.scores, assignment.score_field, None)
        if record_score is None or abs(float(record_score) - assignment.score) > 1e-12:
            raise ValueError(
                f"Calibration score mismatch for {record.record_id} "
                f"on {assignment.score_field}"
            )
        record_runtime = record.provenance.get("runtime_snapshot_id")
        if (
            record_runtime
            and record_runtime != artifact.project.runtime_snapshot_id
        ):
            raise ValueError(
                f"Calibration runtime mismatch for {record.record_id}: "
                f"{record_runtime} != {artifact.project.runtime_snapshot_id}"
            )
        calibration = _calibration_vector(assignment)
        decision = record.decision.model_copy(update={"calibration": calibration})
        provenance = dict(record.provenance)
        provenance["project_calibration"] = {
            "project_id": assignment.project_id,
            "method_id": assignment.method_id,
            "stratum_model_id": assignment.stratum_model_id,
            "qvalue_is_probability": False,
            "production_ranking_changed": False,
            "legacy_confidence_overwritten": False,
        }
        output.append(
            record.model_copy(
                update={"decision": decision, "provenance": provenance}
            )
        )
    return tuple(output)
