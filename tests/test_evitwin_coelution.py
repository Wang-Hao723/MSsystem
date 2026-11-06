from __future__ import annotations

from dataclasses import asdict

import pytest

from matching.authentic_anchor import AuthenticAnchorPolicy, finalize_anchor_library
from matching.evitwin_coelution import (
    EviTwinCoelutionPolicy,
    evaluate_coeluting_reference_pairs,
    validate_coelution_result,
)
from matching.evitwin_ms import build_evitwin_result, validate_evitwin_result


def _library(*, identical: bool = False) -> dict:
    release = AuthenticAnchorPolicy().validated()
    left_10 = [[70.0, 100.0], [90.0, 80.0], [110.0, 60.0], [130.0, 40.0]]
    left_20 = [[71.0, 80.0], [91.0, 100.0], [111.0, 50.0], [131.0, 30.0]]
    right_10 = (
        left_10
        if identical
        else [[75.0, 90.0], [95.0, 70.0], [115.0, 55.0], [135.0, 35.0]]
    )
    right_20 = (
        left_20
        if identical
        else [[76.0, 70.0], [96.0, 90.0], [116.0, 45.0], [136.0, 25.0]]
    )
    return finalize_anchor_library(
        {
            "library_id": "coelution-unit-test",
            "runtime_mode": "shadow_only",
            "construction_policy": {
                "explicit_task_id_only": True,
                "candidate_generation_allowed": False,
                "authentic_library_shadow_search_allowed": True,
            },
            "release_policy": asdict(release),
            "entries": [
                {
                    "reference_id": "LEFT:NEG",
                    "compound_code": "LEFT",
                    "compound_name": "Left isomer",
                    "formula": "C10H10O4",
                    "inchikey": "AAAAAAAAAAAAAA-BBBBBBBBBB-C",
                    "neutral_mass": 194.05791,
                    "polarity": "NEG",
                    "instrument_family": "agilent-qtof-6546",
                    "method_id": "c18-test-method",
                    "authentic_standard": True,
                    "source_task_id": "left-standard",
                    "reference_qc": {"selected_adduct": "[M-H]-"},
                    "spectra": [
                        {
                            "collision_energy_ev": 10.0,
                            "adduct": "[M-H]-",
                            "peaks": left_10,
                        },
                        {
                            "collision_energy_ev": 20.0,
                            "adduct": "[M-H]-",
                            "peaks": left_20,
                        },
                    ],
                },
                {
                    "reference_id": "RIGHT:NEG",
                    "compound_code": "RIGHT",
                    "compound_name": "Right isomer",
                    "formula": "C10H10O4",
                    "inchikey": "CCCCCCCCCCCCCC-DDDDDDDDDD-E",
                    "neutral_mass": 194.05791,
                    "polarity": "NEG",
                    "instrument_family": "agilent-qtof-6546",
                    "method_id": "c18-test-method",
                    "authentic_standard": True,
                    "source_task_id": "right-standard",
                    "reference_qc": {"selected_adduct": "[M-H]-"},
                    "spectra": [
                        {
                            "collision_energy_ev": 10.0,
                            "adduct": "[M-H]-",
                            "peaks": right_10,
                        },
                        {
                            "collision_energy_ev": 20.0,
                            "adduct": "[M-H]-",
                            "peaks": right_20,
                        },
                    ],
                },
            ],
        }
    )


def _mixed_rows() -> list[dict]:
    return [
        {
            "collision_energy": 10.0,
            "query_adduct": "[M-H]-",
            "query_adduct_reliable": True,
            "rt_method_id": "c18-test-method",
            "query_peaks": [
                [70.0, 100.0],
                [90.0, 80.0],
                [110.0, 60.0],
                [130.0, 40.0],
                [75.0, 90.0],
                [95.0, 70.0],
                [115.0, 55.0],
                [135.0, 35.0],
            ],
        },
        {
            "collision_energy": 20.0,
            "query_adduct": "[M-H]-",
            "query_adduct_reliable": True,
            "rt_method_id": "c18-test-method",
            "query_peaks": [
                [71.0, 80.0],
                [91.0, 100.0],
                [111.0, 50.0],
                [131.0, 30.0],
                [76.0, 70.0],
                [96.0, 90.0],
                [116.0, 45.0],
                [136.0, 25.0],
            ],
        },
    ]


def test_multice_dual_authentic_fit_supports_presence_but_not_exact_top1():
    result = evaluate_coeluting_reference_pairs(
        _mixed_rows(),
        _library(),
        ["LEFT:NEG", "RIGHT:NEG"],
    )

    validate_coelution_result(result)
    assert result["decision"]["status"] == "qualified_coeluting_pair_support"
    assert result["decision"]["reportable_scope"] == (
        "coeluting_pair_presence_support"
    )
    pair = result["decision"]["supported_pair"]
    assert pair["supporting_collision_energy_count"] == 2
    assert pair["aggregate_left_specific_peak_count"] >= 2
    assert pair["aggregate_right_specific_peak_count"] >= 2
    assert result["invariants"]["exact_top1_changed"] is False
    assert result["invariants"]["formal_level1_released"] is False


def test_one_collision_energy_cannot_support_a_coeluting_pair():
    result = evaluate_coeluting_reference_pairs(
        _mixed_rows()[:1],
        _library(),
        ["LEFT:NEG", "RIGHT:NEG"],
    )

    assert result["decision"]["status"] == "abstained"
    assert result["pair_audit"][0]["supporting_collision_energy_count"] == 1
    assert (
        result["pair_audit"][0]["aggregate_gates"][
            "multi_collision_energy_support"
        ]
        is False
    )


def test_indistinguishable_reference_spectra_force_abstention():
    library = _library(identical=True)
    rows = [
        {
            "collision_energy": spectrum["collision_energy_ev"],
            "query_peaks": spectrum["peaks"],
        }
        for spectrum in library["entries"][0]["spectra"]
    ]
    result = evaluate_coeluting_reference_pairs(
        rows,
        library,
        ["LEFT:NEG", "RIGHT:NEG"],
    )

    assert result["decision"]["status"] == "abstained"
    assert result["pair_audit"][0]["supporting_collision_energy_count"] == 0
    assert all(
        row["support_gates"]["left_sparse_fit"] is False
        for row in result["pair_audit"][0]["collision_energy_audit"]
    )


def test_policy_forbids_exact_top1_and_level1_release():
    with pytest.raises(ValueError, match="exact Top-1"):
        EviTwinCoelutionPolicy(allow_exact_top1_change=True).validated()
    with pytest.raises(ValueError, match="formal Level 1"):
        EviTwinCoelutionPolicy(allow_formal_level1_release=True).validated()


def test_evitwin_sidecar_exposes_coelution_separately_from_exact_ranking():
    rows = _mixed_rows()
    result = build_evitwin_result(
        {
            "entities": [
                {
                    "entity_id": "coeluting-entity",
                    "representative_precursor_mz": 193.0506335,
                    "representative_retention_time": 4.2,
                    "member_result_indices": [0, 1],
                    "identity_resolution": {
                        "status": "unresolved_same_formula_isomer_group",
                        "reportable_formula": "C10H10O4",
                    },
                    "candidates": [
                        {
                            "candidate": {
                                "formula": "C10H10O4",
                                "inchikey": "AAAAAAAAAAAAAA-BBBBBBBBBB-C",
                            }
                        },
                        {
                            "candidate": {
                                "formula": "C10H10O4",
                                "inchikey": "CCCCCCCCCCCCCC-DDDDDDDDDD-E",
                            }
                        },
                    ],
                }
            ]
        },
        rows,
        _library(),
        task_id="mixture-query",
        polarity="NEG",
        query_instrument_family="agilent-qtof-6546",
    )

    validate_evitwin_result(result)
    assert result["summary"]["coelution_supported_feature_count"] == 1
    assert len(result["coeluting_pairs"]) == 1
    assert result["coeluting_pairs"][0]["exact_top1_changed"] is False
    assert result["coeluting_pairs"][0]["formal_level1_released"] is False
    assert result["invariants"]["coelution_support_can_change_exact_top1"] is False
