from __future__ import annotations

from copy import deepcopy

from matching.feature_fragment_benchmark_v31.readiness import (
    build_eligibility_gap_matrix,
)


def _record(identity: str, formula: str, leaf: str, source: str, lab: str) -> dict:
    return {
        "spectrum_id": identity,
        "compound": {
            "connectivity": identity,
            "formula": formula,
            "chemical_class": {"leaf_class": leaf, "path": [leaf]},
        },
        "source": {
            "source_id": source,
            "laboratory_id": lab,
            "license": {
                "tier": "primary",
                "eligible_for_primary_analysis": True,
            },
        },
        "acquisition": {
            "instrument_domain_id": "domain-a",
            "instrument_domain_resolution": {
                "status": "resolved",
                "exact_condition_metadata_complete": True,
            },
            "instrument_vendor": "Vendor",
            "instrument_model": "Model",
            "analyzer": "QTOF",
            "activation": "CID",
            "activation_location": "collision_cell",
            "collision_cell_region": "cell",
            "polarity": "positive",
            "adduct": "[M+H]+",
            "charge": 1,
            "precursor_mz": 100.0,
            "precursor_isolation": {"status": "verified"},
            "collision_energy": {
                "status": "resolved",
                "conversion_applied": False,
                "unit": "eV",
                "mode": "fixed",
                "center": 20.0,
            },
        },
        "quantitative_eligibility": {
            "range_contains_precursor": True,
            "precursor_selection_verified": True,
            "positive_lr_eligible": True,
            "negative_lr_eligible": True,
        },
        "independence": {
            "source_origin_id": source,
            "laboratory_id": lab,
            "acquisition_replication_group_id": f"acq:{identity}",
            "mirror_cluster_id": f"mirror:{identity}",
        },
        "metadata_completeness": {"missing_fields": []},
    }


def test_gap_matrix_uses_other_classes_in_the_identical_stratum_as_background() -> None:
    thresholds = {
        "minimum_target_identities": 1,
        "minimum_background_identities": 1,
        "minimum_target_formula_families": 1,
        "minimum_background_formula_families": 1,
        "minimum_target_sources": 1,
        "minimum_background_sources": 1,
        "minimum_target_laboratories": 1,
        "minimum_background_laboratories": 1,
    }
    rows, summary = build_eligibility_gap_matrix(
        [
            _record("ABCDEFGHIJKLMN", "C6H12O6", "Class A", "source-a", "lab-a"),
            _record("NOPQRSTUVWXYZA", "C5H10O5", "Class B", "source-b", "lab-b"),
        ],
        thresholds=thresholds,
        planning_required_per_arm=1,
    )
    assert len(rows) == 2
    for row in rows:
        lane = row["evidence_lanes"]["v31_positive"]
        assert lane["target"]["identity_count"] == 1
        assert lane["background"]["identity_count"] == 1
        assert lane["gap"]["passes_count_and_precision_planning_gate"] is True
    assert summary["cells_passing_count_and_precision_planning_gate"][
        "v31_positive"
    ] == 2


def test_identity_with_conflicting_leaf_classes_is_excluded() -> None:
    first = _record("ABCDEFGHIJKLMN", "C6H12O6", "Class A", "source-a", "lab-a")
    second = deepcopy(first)
    second["spectrum_id"] = "second"
    second["compound"]["chemical_class"]["leaf_class"] = "Class B"
    rows, summary = build_eligibility_gap_matrix(
        [first, second],
        thresholds={
            "minimum_target_identities": 1,
            "minimum_background_identities": 1,
            "minimum_target_formula_families": 1,
            "minimum_background_formula_families": 1,
            "minimum_target_sources": 1,
            "minimum_background_sources": 1,
            "minimum_target_laboratories": 1,
            "minimum_background_laboratories": 1,
        },
        planning_required_per_arm=1,
    )
    assert rows == []
    assert summary["identity_class_conflict_count"] == 1
    assert summary["records_excluded_for_identity_class_conflict"] == 2
