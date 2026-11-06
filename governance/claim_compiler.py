"""EviGate 2.0: the single fail-closed compiler for final chemical claims."""

from __future__ import annotations

from collections import Counter
from collections.abc import Iterable
from enum import Enum

from pydantic import BaseModel, ConfigDict, Field, model_validator

from schemas.decision import (
    CalibrationStatus,
    CalibrationVector,
    ClaimDecision,
    ClaimScope,
)
from schemas.identification import IdentificationRecord


EVIGATE2_POLICY_ID = "evigate-ms-v2.0.0-20260803"
DEFAULT_QVALUE_THRESHOLD = 0.05
DEFAULT_MIN_INDEPENDENT_INJECTIONS = 2
DEFAULT_MIN_IDENTITY_AGREEMENT = 2.0 / 3.0


class SampleRole(str, Enum):
    SAMPLE = "sample"
    BLANK = "blank"
    QC = "qc"
    AUTHENTIC_STANDARD = "authentic_standard"


class ClaimCompilationContext(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=True)

    cohort_id: str
    sample_role: SampleRole = SampleRole.SAMPLE
    independent_injection_groups: tuple[str, ...] = ()
    identity_agreement_fraction: float = Field(default=0.0, ge=0.0, le=1.0)
    blank_recurrent: bool = False
    same_method_rt_confirmed: bool = False
    authentic_standard_reference_id: str | None = None
    qvalue_threshold: float = Field(
        default=DEFAULT_QVALUE_THRESHOLD, gt=0.0, le=1.0
    )
    min_independent_injections: int = Field(
        default=DEFAULT_MIN_INDEPENDENT_INJECTIONS, ge=2
    )
    min_identity_agreement: float = Field(
        default=DEFAULT_MIN_IDENTITY_AGREEMENT, gt=0.5, le=1.0
    )

    @model_validator(mode="after")
    def validate_standard_confirmation(self) -> "ClaimCompilationContext":
        if self.same_method_rt_confirmed and not self.authentic_standard_reference_id:
            raise ValueError(
                "Same-method RT confirmation requires an authentic standard reference"
            )
        return self


class ClaimCompilationSummary(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=True)

    policy_id: str = EVIGATE2_POLICY_ID
    record_count: int = Field(ge=0)
    released_count: int = Field(ge=0)
    abstained_count: int = Field(ge=0)
    status_counts: dict[str, int]
    scope_counts: dict[str, int]
    calibrated_count: int = Field(ge=0)
    production_ranking_changed_records: int = Field(ge=0)
    legacy_confidence_overwritten_records: int = Field(ge=0)
    qvalues_treated_as_probabilities: bool = False


def default_task_claim_context(
    task_id: str,
    *,
    sample_role: SampleRole = SampleRole.SAMPLE,
) -> ClaimCompilationContext:
    """One task is one injection group; it cannot self-confirm a claim."""
    return ClaimCompilationContext(
        cohort_id=f"single-task:{task_id}",
        sample_role=sample_role,
        independent_injection_groups=(f"injection:{task_id}",),
        identity_agreement_fraction=1.0,
    )


def _scope_qvalue(calibration: CalibrationVector, scope: ClaimScope) -> float | None:
    if scope in {ClaimScope.AUTHENTIC_STANDARD, ClaimScope.CONNECTIVITY}:
        return calibration.q_connectivity
    if scope == ClaimScope.ISOMER_GROUP:
        return calibration.q_isomer_group
    if scope == ClaimScope.MOLECULAR_FORMULA:
        return calibration.q_formula
    if scope == ClaimScope.CHEMICAL_CLASS:
        return calibration.q_class
    return None


def _effective_scope(
    record: IdentificationRecord,
    context: ClaimCompilationContext,
) -> tuple[ClaimScope, list[str]]:
    reasons: list[str] = []
    if context.same_method_rt_confirmed:
        return ClaimScope.AUTHENTIC_STANDARD, reasons
    scope = record.decision.scope
    ambiguity = tuple(dict.fromkeys(record.decision.ambiguity_identity_keys))
    if scope == ClaimScope.CONNECTIVITY and len(ambiguity) > 1:
        scope = ClaimScope.ISOMER_GROUP
        reasons.append("scope_downgraded_due_to_unresolved_connectivity_ambiguity")
    if not record.candidate.connectivity_key and scope in {
        ClaimScope.CONNECTIVITY,
        ClaimScope.AUTHENTIC_STANDARD,
    }:
        scope = ClaimScope.UNRESOLVED
        reasons.append("canonical_connectivity_key_unavailable")
    return scope, reasons


def compile_claim_decision(
    record: IdentificationRecord,
    context: ClaimCompilationContext,
) -> ClaimDecision:
    """Compile one final decision; retrieval scores and legacy tier are inputs only."""
    scope, scope_reasons = _effective_scope(record, context)
    reason_codes = list(dict.fromkeys(scope_reasons))
    warnings = list(record.decision.warnings)
    if "qvalue_is_not_posterior_probability" not in warnings:
        warnings.append("qvalue_is_not_posterior_probability")
    groups = tuple(
        sorted(
            {
                str(group).strip()
                for group in context.independent_injection_groups
                if str(group).strip()
            }
        )
    )
    calibration = record.decision.calibration
    released = False
    status = "blocked_unresolved_claim_scope"
    requires_replication = True

    if context.sample_role == SampleRole.BLANK:
        scope = ClaimScope.ABSTAIN
        status = "abstain_blank_background_context"
        reason_codes.append("sample_role_is_blank")
    elif context.blank_recurrent:
        scope = ClaimScope.ABSTAIN
        status = "abstain_blank_recurrent_candidate"
        reason_codes.append("candidate_recurrent_in_declared_blank")
    elif not record.decision.eligible_before_replication:
        status = "blocked_pre_release_evidence_gate"
        reason_codes.extend(record.decision.reason_codes or ("legacy_pre_gate_failed",))
    elif scope in {ClaimScope.UNRESOLVED, ClaimScope.ANALOG, ClaimScope.ABSTAIN}:
        status = "blocked_unresolved_claim_scope"
        reason_codes.append("claim_scope_not_release_eligible")
    elif scope == ClaimScope.AUTHENTIC_STANDARD:
        requires_replication = False
        if (
            context.same_method_rt_confirmed
            and context.authentic_standard_reference_id
            and context.sample_role in {SampleRole.SAMPLE, SampleRole.AUTHENTIC_STANDARD, SampleRole.QC}
        ):
            released = True
            status = "released_authentic_standard_same_method_rt"
        else:
            status = "blocked_authentic_standard_evidence_incomplete"
            reason_codes.append("same_method_authentic_standard_rt_required")
    elif calibration.status != CalibrationStatus.CALIBRATED:
        status = "blocked_project_calibration_unavailable"
        reason_codes.append("project_level_qvalue_unavailable")
    else:
        qvalue = _scope_qvalue(calibration, scope)
        if qvalue is None:
            status = "blocked_scope_specific_qvalue_unavailable"
            reason_codes.append(f"qvalue_unavailable_for_{scope.value}")
        elif qvalue > context.qvalue_threshold:
            status = "blocked_project_qvalue_threshold"
            reason_codes.append("project_qvalue_above_release_threshold")
        elif len(groups) < context.min_independent_injections:
            status = "blocked_independent_injection_replication"
            reason_codes.append("insufficient_independent_injection_groups")
        elif context.identity_agreement_fraction < context.min_identity_agreement:
            status = "blocked_cross_injection_identity_conflict"
            reason_codes.append("identity_agreement_below_release_threshold")
        else:
            released = True
            status = "released_project_calibrated_replicate_stable"

    return ClaimDecision(
        policy_id=EVIGATE2_POLICY_ID,
        status=status,
        scope=scope,
        released=released,
        eligible_before_replication=record.decision.eligible_before_replication,
        requires_independent_replication=requires_replication,
        reason_codes=tuple(dict.fromkeys(reason_codes)),
        warnings=tuple(dict.fromkeys(warnings)),
        independent_support_groups=groups,
        independent_injection_count=len(groups),
        ambiguity_identity_keys=record.decision.ambiguity_identity_keys,
        coverage_guarantee=False,
        ranking_changed=record.scores.ranking_changed,
        confidence_field_overwritten=False,
        calibration=calibration,
    )


def compile_identification_record(
    record: IdentificationRecord,
    context: ClaimCompilationContext,
) -> IdentificationRecord:
    decision = compile_claim_decision(record, context)
    provenance = dict(record.provenance)
    provenance["claim_compilation"] = {
        "policy_id": EVIGATE2_POLICY_ID,
        "cohort_id": context.cohort_id,
        "sample_role": context.sample_role.value,
        "qvalue_threshold": context.qvalue_threshold,
        "ranking_changed": False,
        "legacy_confidence_overwritten": False,
    }
    return record.model_copy(
        update={"decision": decision, "provenance": provenance}
    )


def compile_identification_records(
    records: Iterable[IdentificationRecord],
    context: ClaimCompilationContext,
) -> tuple[tuple[IdentificationRecord, ...], ClaimCompilationSummary]:
    compiled = tuple(
        compile_identification_record(record, context) for record in records
    )
    status_counts = Counter(record.decision.status for record in compiled)
    scope_counts = Counter(record.decision.scope.value for record in compiled)
    return compiled, ClaimCompilationSummary(
        record_count=len(compiled),
        released_count=sum(record.decision.released for record in compiled),
        abstained_count=sum(
            record.decision.scope == ClaimScope.ABSTAIN for record in compiled
        ),
        status_counts=dict(sorted(status_counts.items())),
        scope_counts=dict(sorted(scope_counts.items())),
        calibrated_count=sum(
            record.decision.calibration.status == CalibrationStatus.CALIBRATED
            for record in compiled
        ),
        production_ranking_changed_records=sum(
            record.scores.ranking_changed for record in compiled
        ),
        legacy_confidence_overwritten_records=sum(
            record.decision.confidence_field_overwritten for record in compiled
        ),
    )
