from __future__ import annotations

from copy import deepcopy
from pathlib import Path

from matching.feature_fragment_benchmark_v31.common import iter_gzip_jsonl, read_json
from matching.feature_fragment_benchmark_v31.metadata import (
    apply_static_metadata_derivations,
    recompute_v31_eligibility,
)
from matching.feature_fragment_benchmark_v31.screening import screen_source_descriptor


ROOT = Path(__file__).resolve().parents[2]


def test_vendor_manual_derivation_changes_only_whitelisted_static_field() -> None:
    record = next(
        iter(
            iter_gzip_jsonl(
                ROOT
                / "docs/audits/feature_fragment_benchmark_v3/record_manifest.jsonl.gz"
            )
        )
    )
    contract = read_json(
        ROOT
        / "config/feature_fragment_benchmark_v31/metadata_derivation_contract.json"
    )
    before = deepcopy(record)
    derived, events = apply_static_metadata_derivations(record, contract)
    assert events[0]["status"] == "derived"
    assert derived["acquisition"]["collision_cell_region"] == (
        "agilent_axial_acceleration_hexapole_collision_cell"
    )
    assert derived["acquisition"]["activation"] == before["acquisition"]["activation"]
    assert derived["acquisition"]["collision_energy"] == before["acquisition"][
        "collision_energy"
    ]
    assert derived["acquisition"]["precursor_isolation"] == before["acquisition"][
        "precursor_isolation"
    ]
    assert derived["quantitative_eligibility"] == before["quantitative_eligibility"]
    assert recompute_v31_eligibility(derived)["positive_lr_eligible"] is False


def test_gnps_aggregate_is_audit_only_and_nist_commercial_is_excluded() -> None:
    gnps = screen_source_descriptor(
        {
            "repository": "GNPS",
            "stable_per_spectrum_accession": True,
            "stable_record_url": True,
            "license_permits_analysis": True,
        }
    )
    assert gnps["status"] == "audit_only"
    assert "dataset_specific_origin" in gnps["missing_requirement_codes"]
    nist = screen_source_descriptor(
        {"repository": "NIST Tandem Mass Spectral Library"}
    )
    assert nist["status"] == "excluded"
    assert nist["quantitative_use_allowed"] is False


def test_complete_dataset_descriptor_has_independent_negative_gate() -> None:
    descriptor = {
        "repository": "MassBank",
        "origin_dataset_id": "dataset-a",
        "stable_per_spectrum_accession": True,
        "stable_record_url": True,
        "license_permits_analysis": True,
        "reference_identity_complete": True,
        "source_origin_id": "study-a",
        "laboratory_id": "lab-a",
        "exact_instrument_model": True,
        "activation_literal": True,
        "native_collision_energy_program": True,
        "polarity_adduct_and_charge": True,
        "scan_range_contains_precursor": True,
        "precursor_selection_verified": True,
        "raw_acquisition_mapping": True,
        "verified_single_compound_material": False,
        "verified_precursor_isolation_window": False,
        "not_composite_or_in_source": True,
    }
    result = screen_source_descriptor(descriptor)
    assert result["status"] == "quantitative_candidate"
    assert result["quantitative_use_allowed"] is True
    assert result["negative_evidence_allowed"] is False
    assert result["negative_evidence_missing_requirement_codes"] == [
        "verified_precursor_isolation_window"
    ]
