from matching.library_quality import formula_exact_mass
import pytest

from matching.validation_anchor import (
    audit_validation_anchor,
    audit_validation_anchor_batch,
    require_metric_eligible_anchors,
)


def test_valid_positive_anchor_is_metric_eligible():
    mass = formula_exact_mass("C15H20O2")
    truth = {
        "formula": "C15H20O2",
        "mode": "POS",
        "adduct": "[M+H]+",
        "mz": mass + 1.007276466621,
    }
    result = audit_validation_anchor(
        truth,
        [{"entity_id": "e1", "representative_precursor_mz": truth["mz"], "ion_mode": "POS"}],
    )
    assert result["status"] == "valid_observed_anchor"
    assert result["metric_eligible"] is True


def test_opposite_polarity_truth_anchor_exposes_same_mode_correction_candidate():
    mass = formula_exact_mass("C15H20O2")
    truth = {
        "formula": "C15H20O2",
        "mode": "POS",
        "mz": mass - 1.007276466621,
    }
    result = audit_validation_anchor(
        truth,
        [{"entity_id": "positive", "precursor_mz": mass + 1.007276466621, "ion_mode": "POS"}],
    )
    assert result["status"] == "opposite_polarity_formula_anchor"
    assert result["metric_eligible"] is False
    assert result["correction_candidates"][0]["adduct"] == "[M+H]+"
    assert result["automatic_truth_rewrite_allowed"] is False


def test_formula_consistent_but_unobserved_anchor_is_not_search_failure():
    mass = formula_exact_mass("C29H50O")
    truth = {
        "formula": "C29H50O",
        "mode": "POS",
        "adduct": "[M+H]+",
        "mz": mass + 1.007276466621,
    }
    result = audit_validation_anchor(truth, [])
    assert result["status"] == "formula_anchor_not_observed"
    assert result["metric_eligible"] is False
    assert result["correction_candidates"] == []


def test_declared_anchor_missing_but_alternate_adduct_observed():
    mass = formula_exact_mass("C10H12O3")
    truth = {
        "formula": "C10H12O3",
        "mode": "POS",
        "mz": mass + 22.989218,
    }
    result = audit_validation_anchor(
        truth,
        [{"entity_id": "protonated", "precursor_mz": mass + 1.007276466621, "ion_mode": "POS"}],
    )
    assert result["status"] == "declared_anchor_missing_alternate_observed"
    assert result["correction_candidates"][0]["adduct"] == "[M+H]+"


def test_formula_neutral_mass_conflict_fails_closed():
    result = audit_validation_anchor(
        {
            "formula": "C10H12O3",
            "neutral_mass": 999.0,
            "mode": "POS",
            "adduct": "[M+H]+",
        },
        [{"precursor_mz": 1000.007276, "ion_mode": "POS"}],
    )
    assert result["status"] == "formula_neutral_mass_conflict"
    assert result["metric_eligible"] is False


def test_doubly_charged_hypothesis_is_supported():
    mass = formula_exact_mass("C40H80N8O10")
    mz = (mass + 2 * 1.007276466621) / 2
    result = audit_validation_anchor(
        {"formula": "C40H80N8O10", "mode": "POS", "adduct": "[M+2H]2+", "mz": mz},
        [{"precursor_mz": mz, "ion_mode": "POS", "precursor_charge": 2}],
    )
    assert result["status"] == "valid_observed_anchor"
    assert result["declared_formula_adduct_matches"][0]["charge_state"] == 2


def test_permanent_positive_ion_is_a_general_adduct_class():
    mass = formula_exact_mass("C20H24NO4")
    result = audit_validation_anchor(
        {"formula": "C20H24NO4+", "mode": "POS", "adduct": "[M]+", "mz": mass},
        [{"precursor_mz": mass, "ion_mode": "POS"}],
    )
    assert result["status"] == "valid_observed_anchor"
    assert result["intrinsic_formula_charge"] == 1


def test_charged_formula_with_extra_adduct_fails_closed_for_review():
    mass = formula_exact_mass("C20H24NO4")
    result = audit_validation_anchor(
        {"formula": "C20H24NO4+", "neutral_mass": mass, "mode": "POS", "adduct": "[M+2H]2+"},
        [{"precursor_mz": (mass + 2 * 1.007276466621) / 2, "ion_mode": "POS"}],
    )
    assert result["status"] == "intrinsic_formula_charge_adduct_review_required"
    assert result["metric_eligible"] is False


def test_nominal_mz_is_not_treated_as_exact_hrms_anchor():
    mass = formula_exact_mass("C15H12O5")
    exact = mass - 1.007276466621
    result = audit_validation_anchor(
        {"formula": "C15H12O5", "mode": "NEG", "mz": 271.1},
        [{"precursor_mz": exact, "ion_mode": "NEG"}],
    )
    assert result["status"] == "nominal_mz_not_exact_hrms_anchor"
    assert result["inferred_declared_mz_resolution_da"] == 0.1
    assert result["metric_eligible"] is False


def test_batch_gate_excludes_invalid_anchor_and_never_rewrites_truth():
    mass = formula_exact_mass("C15H20O2")
    batch = audit_validation_anchor_batch(
        [
            {
                "record_id": "valid",
                "formula": "C15H20O2",
                "mode": "POS",
                "adduct": "[M+H]+",
                "mz": mass + 1.007276466621,
            },
            {
                "record_id": "wrong-polarity",
                "formula": "C15H20O2",
                "mode": "POS",
                "mz": mass - 1.007276466621,
            },
        ],
        {
            "valid": [
                {"precursor_mz": mass + 1.007276466621, "ion_mode": "POS"}
            ],
            "wrong-polarity": [
                {"precursor_mz": mass + 1.007276466621, "ion_mode": "POS"}
            ],
        },
    )

    assert batch["metric_eligible_record_ids"] == ["valid"]
    assert batch["summary"]["excluded_from_accuracy_metrics"] == 1
    assert batch["automatic_truth_rewrite_allowed"] is False
    with pytest.raises(ValueError, match="anchor gate failed"):
        require_metric_eligible_anchors(batch)
