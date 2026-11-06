from __future__ import annotations

from pathlib import Path

from matching.feature_fragment_benchmark_v31.chemistry import (
    grade_neutral_loss_balance,
    grade_product_formula,
)
from matching.feature_fragment_benchmark_v31.common import read_json
from matching.feature_fragment_benchmark_v31.rule_contract import validate_rule_record


ROOT = Path(__file__).resolve().parents[2]


def _valid_rule() -> dict:
    return {
        "rule_id": "test-rule",
        "rule_kind": "single_fragment",
        "source_evidence": [
            {
                "source_id": "source-a",
                "source_type": "public_experimental_spectrum",
                "stable_locator": "https://example.invalid/record-a",
                "supported_claim": "fragment event observed in the frozen stratum",
            }
        ],
        "chemical_scope": {"ontology": "ChemOnt", "leaf_classes": ["Flavonoids"]},
        "adducts": ["[M-H]-"],
        "polarities": ["negative"],
        "instrument_scope": {
            "instrument_domain_ids": ["domain-a"],
            "activation": ["HCD"],
            "collision_cell_regions": ["thermo_hcd_cell"],
        },
        "collision_energy_scope": {
            "native_unit": "NCE",
            "acquisition_mode": "fixed",
            "minimum": 20.0,
            "maximum": 25.0,
            "program_values": [20.0],
            "program_weights": [1.0],
            "cross_unit_conversion_applied": False,
        },
        "mass_tolerance": {
            "mz_tolerance_da": 0.01,
            "mz_tolerance_ppm": 10.0,
            "calibration_basis": "sealed calibration identities",
        },
        "target_frequency": {
            "event_count": 60,
            "identity_denominator": 100,
            "frequency": 0.6,
        },
        "background_frequency": {
            "event_count": 10,
            "identity_denominator": 100,
            "frequency": 0.1,
        },
        "fragment_pathway_group_id": "fragpath:a",
        "positive_evidence_limit": {"claim": "class_support_only"},
        "negative_evidence_limit": {"verified_isolation_required": True},
        "likelihood_ratio_positive": {"point": 6.0, "lower_95": 2.0},
        "likelihood_ratio_negative": {"point": 0.44, "upper_95": 0.8},
        "evidence_grade": "G4",
        "validation": {
            "identity_isolated": True,
            "formula_isolated": True,
            "source_isolated": True,
            "laboratory_isolated": True,
            "mirror_isolated": True,
            "raw_acquisition_isolated": True,
        },
        "identity_upgrade_allowed": False,
        "formula_upgrade_allowed": False,
        "isomer_resolution_allowed": False,
        "standard_confirmation_allowed": False,
    }


def test_formula_compatibility_is_graded_and_score_neutral() -> None:
    compatible = grade_product_formula(
        "C6H12O6", "C6H10O5", ion_state_complete=True
    )
    assert compatible["status"] == "compatible"
    assert compatible["score_contribution"] == 0.0
    ambiguous = grade_product_formula(
        "C6H12O6", "C7H10O5", ion_state_complete=False
    )
    assert ambiguous["status"] == "ambiguous"
    incompatible = grade_product_formula(
        "C6H12O6", "C7H10O5", ion_state_complete=True
    )
    assert incompatible["status"] == "incompatible"
    balance = grade_neutral_loss_balance(
        "C6H12O6", "C6H10O5", "H2O", ion_state_complete=True
    )
    assert balance["status"] == "compatible"


def test_rule_contract_requires_auditable_fields_and_forbids_fixed_points() -> None:
    schema = read_json(
        ROOT / "config/feature_fragment_benchmark_v31/rule_schema.json"
    )
    valid = validate_rule_record(_valid_rule(), schema)
    assert valid == {"status": "valid", "errors": [], "runtime_export_allowed": True}
    invalid_rule = _valid_rule()
    invalid_rule["one_fragment_points"] = 5
    invalid_rule["source_evidence"] = []
    invalid = validate_rule_record(invalid_rule, schema)
    assert invalid["status"] == "invalid"
    assert invalid["runtime_export_allowed"] is False
    assert "forbidden_fixed_score_field:one_fragment_points" in invalid["errors"]
    assert "source_evidence_missing" in invalid["errors"]
