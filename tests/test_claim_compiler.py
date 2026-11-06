from __future__ import annotations

from schemas.decision import CalibrationStatus, CalibrationVector, ClaimScope
from schemas.legacy import legacy_row_to_identification
from governance.claim_compiler import (
    ClaimCompilationContext,
    SampleRole,
    compile_claim_decision,
    compile_identification_record,
)


def _record(*, ambiguity: bool = False):
    row = {
        "spectrum_id": "scan=1",
        "precursor_mz": 301.0,
        "ion_mode": "NEG",
        "match": {
            "compound_name": "Candidate",
            "inchikey": "AAAAAAAAAAAAAA-UHFFFAOYSA-N",
            "rank_score": 0.91,
            "final_score": 0.91,
            "confidence": "High",
            "evidence_governance": {
                "policy_id": "evigate-v1",
                "claim_status": "provisional_needs_independent_injection",
                "claim_scope": "connectivity_level_library_candidate",
                "release_eligible_before_replication": True,
                "ambiguity_set": (
                    [
                        {"identity_key": "ik14:AAAAAAAAAAAAAA"},
                        {"identity_key": "ik14:BBBBBBBBBBBBBB"},
                    ]
                    if ambiguity
                    else [{"identity_key": "ik14:AAAAAAAAAAAAAA"}]
                ),
                "nodes": [],
            },
        },
    }
    record = legacy_row_to_identification("task-1", row)
    assert record is not None
    return record


def _context(**updates):
    payload = {
        "cohort_id": "cohort-1",
        "sample_role": SampleRole.SAMPLE,
        "independent_injection_groups": ("inj-1", "inj-2"),
        "identity_agreement_fraction": 1.0,
    }
    payload.update(updates)
    return ClaimCompilationContext(**payload)


def _with_qvalue(record, *, scope=ClaimScope.CONNECTIVITY, qvalue=0.01):
    fields = {
        ClaimScope.CONNECTIVITY: {"q_connectivity": qvalue},
        ClaimScope.ISOMER_GROUP: {"q_isomer_group": qvalue},
        ClaimScope.MOLECULAR_FORMULA: {"q_formula": qvalue},
        ClaimScope.CHEMICAL_CLASS: {"q_class": qvalue},
    }[scope]
    calibration = CalibrationVector(
        status=CalibrationStatus.CALIBRATED,
        method_id="method-1",
        project_id="project-1",
        stratum_model_id="stratum-1",
        unit_level="precursor_event",
        score_field="rank_score",
        transfer_scope="project_local",
        **fields,
    )
    decision = record.decision.model_copy(update={"calibration": calibration})
    return record.model_copy(update={"decision": decision})


def test_uncalibrated_legacy_high_is_not_released():
    decision = compile_claim_decision(_record(), _context())
    assert decision.released is False
    assert decision.status == "blocked_project_calibration_unavailable"
    assert "project_level_qvalue_unavailable" in decision.reason_codes


def test_qvalue_alone_cannot_replace_independent_injections():
    record = _with_qvalue(_record())
    decision = compile_claim_decision(
        record,
        _context(independent_injection_groups=("technical-group-1",)),
    )
    assert decision.released is False
    assert decision.status == "blocked_independent_injection_replication"
    assert decision.independent_injection_count == 1


def test_calibrated_replicate_stable_connectivity_can_release_without_reranking():
    record = _with_qvalue(_record())
    compiled = compile_identification_record(record, _context())
    assert compiled.decision.released is True
    assert compiled.decision.scope == ClaimScope.CONNECTIVITY
    assert compiled.scores == record.scores
    assert compiled.candidate == record.candidate
    assert compiled.decision.confidence_field_overwritten is False


def test_connectivity_ambiguity_downgrades_and_requires_isomer_qvalue():
    record = _with_qvalue(_record(ambiguity=True), scope=ClaimScope.CONNECTIVITY)
    blocked = compile_claim_decision(record, _context())
    assert blocked.scope == ClaimScope.ISOMER_GROUP
    assert blocked.status == "blocked_scope_specific_qvalue_unavailable"

    record = _with_qvalue(_record(ambiguity=True), scope=ClaimScope.ISOMER_GROUP)
    released = compile_claim_decision(record, _context())
    assert released.released is True
    assert released.scope == ClaimScope.ISOMER_GROUP


def test_blank_recurrence_forces_abstention_even_with_qvalue():
    record = _with_qvalue(_record())
    decision = compile_claim_decision(record, _context(blank_recurrent=True))
    assert decision.released is False
    assert decision.scope == ClaimScope.ABSTAIN
    assert decision.status == "abstain_blank_recurrent_candidate"


def test_same_method_authentic_standard_can_release_without_fake_qvalue():
    decision = compile_claim_decision(
        _record(),
        _context(
            sample_role=SampleRole.AUTHENTIC_STANDARD,
            same_method_rt_confirmed=True,
            authentic_standard_reference_id="std-001",
            independent_injection_groups=("inj-1",),
        ),
    )
    assert decision.released is True
    assert decision.scope == ClaimScope.AUTHENTIC_STANDARD
    assert decision.calibration.status == CalibrationStatus.UNAVAILABLE
    assert decision.requires_independent_replication is False


def test_declared_blank_task_is_qc_context_not_a_sample_claim():
    record = _with_qvalue(_record())
    decision = compile_claim_decision(
        record, _context(sample_role=SampleRole.BLANK)
    )
    assert decision.scope == ClaimScope.ABSTAIN
    assert decision.status == "abstain_blank_background_context"
