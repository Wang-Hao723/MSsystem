from __future__ import annotations

from copy import deepcopy

from matching.feature_fragment_benchmark_v3_1.policy import (
    apply_v31_policy,
    composite_program_weight_status,
)


def test_fixed_energy_keeps_lr_positive_but_unknown_purity_blocks_lr_negative(
    v31_record_factory, qc_policy
) -> None:
    row = apply_v31_policy(v31_record_factory(), qc_policy=qc_policy)
    assert row["quantitative_eligibility"]["positive_lr_eligible"] is True
    assert row["quantitative_eligibility"]["negative_lr_eligible"] is False
    assert "v31_precursor_purity_unavailable_for_lr_negative" in row[
        "quantitative_eligibility"
    ]["ineligible_reason_codes"]
    assert row["v3_1_amendment"]["nce_to_ev_conversion_applied"] is False
    assert row["v3_1_amendment"]["nce_to_ev_derived_value"] is None


def test_frozen_purity_threshold_can_qualify_fixed_precursor_selected_spectrum(
    v31_record_factory, qc_policy
) -> None:
    policy = deepcopy(qc_policy)
    policy["precursor_purity"]["threshold_status"] = "calibrated_and_frozen"
    policy["precursor_purity"]["minimum_fraction"] = 0.8
    record = v31_record_factory()
    record["acquisition"]["precursor_purity"] = {
        "status": "verified",
        "fraction": 0.91,
    }
    row = apply_v31_policy(record, qc_policy=policy)
    assert row["quantitative_eligibility"]["positive_lr_eligible"] is True
    assert row["quantitative_eligibility"]["negative_lr_eligible"] is True


def test_missing_composite_weights_blocks_quantitative_likelihood_ratios(
    v31_record_factory, qc_policy
) -> None:
    record = v31_record_factory()
    record["acquisition"]["collision_energy"].update(
        {"mode": "stepped", "values": [20.0, 30.0, 40.0], "weights": [None] * 3}
    )
    row = apply_v31_policy(record, qc_policy=qc_policy)
    assert row["quantitative_eligibility"]["positive_lr_eligible"] is False
    assert row["quantitative_eligibility"]["negative_lr_eligible"] is False
    assert row["v3_1_amendment"]["observed_evidence_grade_ceiling"] == "G1"
    assert "v31_composite_numeric_weights_incomplete" in row[
        "quantitative_eligibility"
    ]["ineligible_reason_codes"]


def test_complete_stepped_weights_allow_lr_positive_but_never_lr_negative(
    v31_record_factory, qc_policy
) -> None:
    record = v31_record_factory()
    record["acquisition"]["collision_energy"].update(
        {"mode": "stepped", "values": [20.0, 30.0, 40.0], "weights": [1.0, 1.0, 1.0]}
    )
    row = apply_v31_policy(record, qc_policy=qc_policy)
    assert row["quantitative_eligibility"]["positive_lr_eligible"] is True
    assert row["quantitative_eligibility"]["negative_lr_eligible"] is False
    assert row["quantitative_eligibility"]["positive_only_composite_or_in_source"] is True


def test_ramp_endpoints_and_weights_do_not_replace_program_validation(
    v31_record_factory, qc_policy
) -> None:
    record = v31_record_factory()
    record["acquisition"]["collision_energy"].update(
        {"mode": "ramp", "values": [10.0, 40.0], "weights": [0.5, 0.5]}
    )
    status = composite_program_weight_status(record["acquisition"]["collision_energy"])
    assert status["numeric_weights_complete"] is True
    assert status["same_model_program_detection_floor_validated"] is False
    row = apply_v31_policy(record, qc_policy=qc_policy)
    assert row["quantitative_eligibility"]["positive_lr_eligible"] is False
    assert row["quantitative_eligibility"]["negative_lr_eligible"] is False


def test_v3_1_never_promotes_a_v3_ineligible_record(v31_record_factory, qc_policy) -> None:
    record = v31_record_factory()
    record["quantitative_eligibility"]["positive_lr_eligible"] = False
    record["quantitative_eligibility"]["negative_lr_eligible"] = False
    row = apply_v31_policy(record, qc_policy=qc_policy)
    assert row["quantitative_eligibility"]["positive_lr_eligible"] is False
    assert row["quantitative_eligibility"]["negative_lr_eligible"] is False


def test_nce_stays_native_and_has_no_derived_ev(v31_record_factory, qc_policy) -> None:
    record = v31_record_factory()
    record["acquisition"]["collision_energy"].update(
        {"unit": "NCE", "values": [30.0], "weights": [1.0]}
    )
    row = apply_v31_policy(record, qc_policy=qc_policy)
    assert row["acquisition"]["collision_energy"]["unit"] == "NCE"
    assert row["v3_1_amendment"]["nce_to_ev_conversion_applied"] is False
    assert row["v3_1_amendment"]["nce_to_ev_derived_value"] is None
    assert row["v3_1_amendment"]["spectral_qc"]["minimum_peak_gate_active"] is False
