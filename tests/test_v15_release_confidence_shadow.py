from __future__ import annotations

from pathlib import Path

import pytest

from calibration.release_confidence import (
    IsolationRecord,
    PairedReleaseOutcome,
    PartitionRole,
    QualificationEvidence,
    ShadowCandidate,
    audit_crossfit_assignments,
    audit_partition_isolation,
    build_crossfit_assignments,
    evaluate_runtime_review_qualification,
    evaluate_shadow_candidate,
    exact_one_sided_mcnemar_pvalue,
    exact_two_sided_mcnemar_pvalue,
    file_sha256,
    load_shadow_protocol,
    minimum_zero_event_sample_size,
    paired_cluster_bootstrap_difference,
    shadow_rule_from_protocol,
    zero_event_one_sided_upper_bound,
)
from eval.v15_release_confidence_shadow import build_receipt


ROOT = Path(__file__).resolve().parents[1]
PROTOCOL = ROOT / "config/v15_release_confidence_shadow_protocol.json"


def _isolation_record(
    index: int,
    partition: PartitionRole,
    **updates,
) -> IsolationRecord:
    payload = {
        "record_id": f"record-{index}",
        "partition": partition,
        "connectivity_id": f"connectivity-{index}",
        "molecular_formula": f"C{index + 10}H{index + 12}O2",
        "formula_family_id": f"formula-{index}",
        "source_dataset_id": f"dataset-{index}",
        "source_record_ids": (f"source-record-{index}",),
        "source_content_fingerprints": (f"fingerprint-{index}",),
        "laboratory_id": f"laboratory-{index}",
        "instrument_instance_id": f"instrument-{index}",
        "instrument_domain_id": "agilent-qtof-6546-cid",
        "acquisition_batch_id": f"batch-{index}",
        "polarity": "NEG",
        "adduct": "[M-H]-",
    }
    payload.update(updates)
    return IsolationRecord(**payload)


def _eligible_shadow_candidate(**updates) -> ShadowCandidate:
    payload = {
        "record_id": "candidate-1",
        "base_formal_2d_release": True,
        "current_formal_2d_release": False,
        "current_gate_failed_checks": (
            "winner_margin_below_cross_platform_release_gate",
        ),
        "winner_entity_score": 0.91,
        "winner_margin": 0.08,
        "winner_ce_support": 3,
        "query_quality_status": "passed",
        "query_quality_failed_checks": (),
        "same_formula_connectivity_count": 2,
        "unique_compatible_reference_profiles": 2,
        "compatible_reference_provenance_strata": 2,
        "independent_acquisition_support_count": 2,
        "independent_acquisition_top1_consistent": True,
        "development_domain_case_count": 5,
        "formula_conflict": False,
        "adduct_conflict": False,
        "precursor_conflict": False,
        "sample_source_compatibility_conflict": False,
        "historical_error_sentinel": False,
        "eligibility_metadata_complete": True,
    }
    payload.update(updates)
    return ShadowCandidate(**payload)


def test_partition_audit_accepts_new_labs_and_instruments_in_same_domain():
    records = [
        _isolation_record(1, PartitionRole.DEVELOPMENT),
        _isolation_record(2, PartitionRole.TUNING),
        _isolation_record(3, PartitionRole.SEALED_VALIDATION),
    ]
    audit = audit_partition_isolation(records)

    assert audit.eligible is True
    assert audit.reason_codes == ()
    assert audit.transfer_scope == (
        "same_domain_new_laboratory_and_instrument_instance"
    )
    assert audit.coarse_instrument_domain_overlap == ("agilent-qtof-6546-cid",)


def test_partition_audit_fails_on_formula_mirror_lab_and_instrument_leakage():
    development = _isolation_record(1, PartitionRole.DEVELOPMENT)
    tuning = _isolation_record(2, PartitionRole.TUNING)
    sealed = _isolation_record(
        3,
        PartitionRole.SEALED_VALIDATION,
        formula_family_id=development.formula_family_id,
        source_content_fingerprints=development.source_content_fingerprints,
        laboratory_id=development.laboratory_id,
        instrument_instance_id=development.instrument_instance_id,
    )
    audit = audit_partition_isolation([development, tuning, sealed])

    assert audit.eligible is False
    assert "cross_partition_formula_family_id_overlap" in audit.reason_codes
    assert "cross_partition_source_content_fingerprints_overlap" in audit.reason_codes
    assert "cross_partition_laboratory_id_overlap" in audit.reason_codes
    assert "cross_partition_instrument_instance_id_overlap" in audit.reason_codes


def test_grouped_crossfit_keeps_formula_and_mirror_components_together():
    records = [
        _isolation_record(
            1,
            PartitionRole.DEVELOPMENT,
            formula_family_id="shared-formula",
        ),
        _isolation_record(
            2,
            PartitionRole.DEVELOPMENT,
            formula_family_id="shared-formula",
        ),
        _isolation_record(
            3,
            PartitionRole.DEVELOPMENT,
            source_content_fingerprints=("shared-mirror",),
        ),
        _isolation_record(
            4,
            PartitionRole.DEVELOPMENT,
            source_content_fingerprints=("shared-mirror",),
        ),
        _isolation_record(5, PartitionRole.DEVELOPMENT),
        _isolation_record(6, PartitionRole.DEVELOPMENT),
    ]
    assignments = build_crossfit_assignments(records, n_folds=3, seed=20260827)
    audit = audit_crossfit_assignments(records, assignments)
    fold_by_record = {item.record_id: item.fold for item in assignments}

    assert audit.eligible is True
    assert fold_by_record["record-1"] == fold_by_record["record-2"]
    assert fold_by_record["record-3"] == fold_by_record["record-4"]
    assert len(audit.fold_counts) == 3


def test_fixed_margin_only_rule_is_label_free_and_never_changes_runtime():
    protocol = load_shadow_protocol(PROTOCOL)
    rule = shadow_rule_from_protocol(protocol)
    decision = evaluate_shadow_candidate(_eligible_shadow_candidate(), rule)

    assert "top1_correct" not in ShadowCandidate.model_fields
    assert decision.eligible is True
    assert decision.would_restore_2d_connectivity_release is True
    assert decision.runtime_influence is False
    assert decision.production_release_changed is False
    assert decision.top1_changed is False
    assert decision.entity_score_changed is False
    assert decision.q_value is None
    assert decision.probability_correct is None
    assert decision.score_is_probability is False


def test_historical_error_sentinel_can_never_pass_the_shadow_rule():
    rule = shadow_rule_from_protocol(load_shadow_protocol(PROTOCOL))
    decision = evaluate_shadow_candidate(
        _eligible_shadow_candidate(
            record_id="A0509:NEG",
            historical_error_sentinel=False,
        ),
        rule,
    )

    assert decision.eligible is False
    assert decision.would_restore_2d_connectivity_release is False
    assert "historical_error_sentinel_locked" in decision.reason_codes


def test_non_margin_failure_and_missing_metadata_fail_closed():
    rule = shadow_rule_from_protocol(load_shadow_protocol(PROTOCOL))
    non_margin = evaluate_shadow_candidate(
        _eligible_shadow_candidate(
            current_gate_failed_checks=(
                "entity_score_below_cross_platform_release_gate",
                "winner_margin_below_cross_platform_release_gate",
            )
        ),
        rule,
    )
    missing = evaluate_shadow_candidate(
        _eligible_shadow_candidate(
            eligibility_metadata_complete=False,
            unique_compatible_reference_profiles=None,
            formula_conflict=None,
        ),
        rule,
    )

    assert "not_a_margin_only_current_gate_failure" in non_margin.reason_codes
    assert "eligibility_metadata_incomplete" in missing.reason_codes
    assert "missing_unique_compatible_reference_profiles" in missing.reason_codes
    assert "missing_formula_conflict_audit" in missing.reason_codes


def test_zero_event_bound_makes_119_conditional_on_the_selected_margin():
    assert minimum_zero_event_sample_size(0.025, alpha=0.05) == 119
    assert zero_event_one_sided_upper_bound(118) > 0.025
    assert zero_event_one_sided_upper_bound(119) <= 0.025


def test_runtime_review_gate_requires_sealed_safety_and_recovery():
    base = {
        "split_isolation_passed": True,
        "crossfit_audit_passed": True,
        "rule_frozen_before_sealed_labels": True,
        "external_validation_complete": True,
        "multiplicity_plan_followed": True,
        "historical_error_sentinels_rescued": 0,
        "sealed_evaluable_units": 119,
        "sealed_new_wrong_releases": 0,
        "sealed_correct_releases_recovered": 5,
        "sealed_correct_releases_lost": 0,
        "project_target_decoy_ready": False,
    }
    qualified = evaluate_runtime_review_qualification(QualificationEvidence(**base))
    too_small = evaluate_runtime_review_qualification(
        QualificationEvidence(**{**base, "sealed_evaluable_units": 118})
    )
    regression = evaluate_runtime_review_qualification(
        QualificationEvidence(
            **{
                **base,
                "historical_error_sentinels_rescued": 1,
            }
        )
    )
    descriptive_only = evaluate_runtime_review_qualification(
        QualificationEvidence(
            **{
                **base,
                "sealed_correct_releases_recovered": 1,
            }
        )
    )
    target_decoy_claimed_ready = evaluate_runtime_review_qualification(
        QualificationEvidence(
            **{
                **base,
                "project_target_decoy_ready": True,
            }
        )
    )

    assert qualified.eligible_for_runtime_review is True
    assert qualified.sealed_recovery_one_sided_mcnemar_pvalue == pytest.approx(0.03125)
    assert qualified.minimum_correct_gains_when_no_losses == 5
    assert qualified.production_activation_automatic is False
    assert qualified.qvalues_must_remain_null is True
    assert too_small.eligible_for_runtime_review is False
    assert "sealed_safety_interval_too_wide" in too_small.reason_codes
    assert regression.eligible_for_runtime_review is False
    assert "historical_error_sentinel_regression" in regression.reason_codes
    assert descriptive_only.eligible_for_runtime_review is False
    assert "sealed_recovery_mcnemar_not_significant" in descriptive_only.reason_codes
    assert target_decoy_claimed_ready.eligible_for_runtime_review is True
    assert target_decoy_claimed_ready.qvalues_must_remain_null is True


def test_paired_statistics_keep_records_clustered():
    outcomes = [
        PairedReleaseOutcome(
            record_id=f"record-{index}",
            connectivity_cluster_id=f"identity-{index // 2}",
            formula_family_id=f"formula-{index // 2}",
            current_correct_release=False,
            shadow_correct_release=index < 2,
            current_wrong_release=False,
            shadow_wrong_release=False,
        )
        for index in range(6)
    ]
    interval = paired_cluster_bootstrap_difference(
        outcomes,
        endpoint="correct_release",
        cluster_field="formula_family_id",
        resamples=200,
        seed=7,
    )

    assert interval.estimate == pytest.approx(2 / 6)
    assert interval.lower <= interval.estimate <= interval.upper
    assert exact_one_sided_mcnemar_pvalue(5, 0) == pytest.approx(0.03125)
    assert exact_one_sided_mcnemar_pvalue(4, 0) == pytest.approx(0.0625)
    assert exact_two_sided_mcnemar_pvalue(6, 0) == pytest.approx(0.03125)


def test_historical_receipt_rebuilds_counts_and_keeps_qvalues_null():
    receipt = build_receipt(PROTOCOL)

    assert receipt["protocol"]["file_sha256"] == file_sha256(PROTOCOL)
    assert receipt["rebuilt_counts"]["combined"] == {
        "eligible_records": 138,
        "base_releases": 115,
        "base_correct_releases": 94,
        "base_wrong_releases": 21,
        "current_releases": 93,
        "current_correct_releases": 88,
        "current_wrong_releases": 5,
        "demoted_records": 22,
        "demoted_correct_records": 6,
        "demoted_wrong_records": 16,
    }
    assert receipt["threshold_provenance_audit"]["mix33_stored_pretransfer_0_75"] == {
        "releases": 47,
        "correct_releases": 45,
        "wrong_releases": 2,
    }
    assert receipt["threshold_provenance_audit"]["mix33_recomputed_current_0_70"] == {
        "releases": 49,
        "correct_releases": 47,
        "wrong_releases": 2,
    }
    shadow = receipt["historical_shadow_evaluation"]
    assert shadow["demoted_records_evaluated"] == 22
    assert shadow["eligible_shadow_records"] == 0
    assert shadow["historical_error_sentinels_restored"] == 0
    assert receipt["project_target_decoy"]["q_value"] is None
    assert receipt["project_target_decoy"]["qvalues_must_remain_null"] is True
    assert receipt["qualification"]["minimum_correct_gains_when_no_losses"] == 5
    assert receipt["sealed_validation_feasibility"] == {
        "historical_unique_connectivity_universe": 86,
        "historical_universe_is_eligible_as_sealed_validation": False,
        "currently_qualified_sealed_units": 0,
        "minimum_sealed_units_if_zero_new_wrong_releases": 119,
        "confirmatory_safety_endpoint_reachable_now": False,
        "external_data_ingestion_performed": False,
        "new_authentic_standard_purchase_or_request_performed": False,
        "fallback_claim": (
            "internal_stability_supported_research_hypothesis_not_production_qualification"
        ),
    }
    assert receipt["execution"]["production_release_modified"] is False
    assert receipt["decision"]["runtime_integration_allowed"] is False
    for item in receipt["implementation"]["files"]:
        path = Path(item["path"])
        assert path.exists()
        assert item["sha256"] == file_sha256(path)
