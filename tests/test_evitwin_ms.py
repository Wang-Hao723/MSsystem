from __future__ import annotations

import copy
from dataclasses import asdict

import pytest

from matching.authentic_anchor import AuthenticAnchorPolicy, finalize_anchor_library
from matching.evitwin_ms import (
    EVITWIN_ALGORITHM_ID,
    EviTwinPolicy,
    build_evitwin_result,
    validate_evitwin_result,
)


def _library(*, source_task_id: str = "reference-task") -> dict:
    policy = AuthenticAnchorPolicy().validated()
    return finalize_anchor_library(
        {
            "library_id": "unit-test",
            "runtime_mode": "shadow_only",
            "construction_policy": {
                "explicit_task_id_only": True,
                "candidate_generation_allowed": False,
                "authentic_library_shadow_search_allowed": True,
            },
            "release_policy": asdict(policy),
            "entries": [
                {
                    "reference_id": "RIGHT:NEG",
                    "compound_code": "RIGHT",
                    "compound_name": "Right candidate",
                    "formula": "C10H10O4",
                    "inchikey": "AAAAAAAAAAAAAA-BBBBBBBBBB-C",
                    "neutral_mass": 194.05791,
                    "polarity": "NEG",
                    "instrument_family": "agilent-qtof-6546",
                    "method_id": "c18-test-method",
                    "authentic_standard": True,
                    "source_task_id": source_task_id,
                    "reference_qc": {"selected_adduct": "[M-H]-"},
                    "spectra": [
                        {
                            "collision_energy_ev": 10.0,
                            "adduct": "[M-H]-",
                            "peaks": [[80.0, 100.0], [100.0, 60.0], [120.0, 40.0]],
                        },
                        {
                            "collision_energy_ev": 20.0,
                            "adduct": "[M-H]-",
                            "peaks": [[80.0, 80.0], [100.0, 100.0], [140.0, 30.0]],
                        },
                    ],
                },
                {
                    "reference_id": "WRONG:NEG",
                    "compound_code": "WRONG",
                    "compound_name": "Wrong candidate",
                    "formula": "C10H10O4",
                    "inchikey": "CCCCCCCCCCCCCC-DDDDDDDDDD-E",
                    "neutral_mass": 194.05791,
                    "polarity": "NEG",
                    "instrument_family": "agilent-qtof-6546",
                    "method_id": "c18-test-method",
                    "authentic_standard": True,
                    "source_task_id": "other-reference-task",
                    "reference_qc": {"selected_adduct": "[M-H]-"},
                    "spectra": [
                        {
                            "collision_energy_ev": 10.0,
                            "adduct": "[M-H]-",
                            "peaks": [[81.0, 100.0], [101.0, 60.0], [121.0, 40.0]],
                        },
                        {
                            "collision_energy_ev": 20.0,
                            "adduct": "[M-H]-",
                            "peaks": [[81.0, 80.0], [101.0, 100.0], [141.0, 30.0]],
                        },
                    ],
                },
            ],
        }
    )


def _v11() -> dict:
    return {
        "metadata": {
            "task_id": "query-task",
            "module_id": "v15_conditioned_multilibrary_flagship",
            "algorithm_id": "mssystem-v15-test",
        },
        "entities": [
            {
                "entity_id": "entity-1",
                "representative_precursor_mz": 193.0506335,
                "representative_retention_time": 4.2,
                "member_result_indices": [0, 1],
                "identity_resolution": {
                    "status": "unresolved_same_formula_isomer_group",
                    "reportable_formula": "C10H10O4",
                },
                "candidates": [
                    {
                        "identity_key": "std2d:CCCCCCCCCCCCCC",
                        "vetoed": False,
                        "candidate": {
                            "formula": "C10H10O4",
                            "inchikey": "CCCCCCCCCCCCCC-DDDDDDDDDD-E",
                        },
                    },
                    {
                        "identity_key": "std2d:AAAAAAAAAAAAAA",
                        "vetoed": False,
                        "candidate": {
                            "formula": "C10H10O4",
                            "inchikey": "AAAAAAAAAAAAAA-BBBBBBBBBB-C",
                        },
                    },
                ],
                "winner": {
                    "entity_score": 0.82,
                    "candidate": {
                        "compound_name": "Wrong public winner",
                        "formula": "C10H10O4",
                        "inchikey": "CCCCCCCCCCCCCC-DDDDDDDDDD-E",
                    },
                },
            }
        ],
    }


def _rows() -> list[dict]:
    return [
        {
            "collision_energy": 10.0,
            "query_adduct": "[M-H]-",
            "query_adduct_reliable": True,
            "rt_method_id": "c18-test-method",
            "query_peaks": [[80.0, 100.0], [100.0, 60.0], [120.0, 40.0]],
        },
        {
            "collision_energy": 20.0,
            "query_adduct": "[M-H]-",
            "query_adduct_reliable": True,
            "rt_method_id": "c18-test-method",
            "query_peaks": [[80.0, 80.0], [100.0, 100.0], [140.0, 30.0]],
        },
    ]


def test_evitwin_promotes_multi_ce_authentic_candidate_without_mutating_v11():
    source = _v11()
    result = build_evitwin_result(
        source,
        _rows(),
        _library(),
        task_id="query-task",
        polarity="NEG",
        query_instrument_family="agilent-qtof-6546",
    )

    validate_evitwin_result(result)
    assert result["metadata"]["algorithm_id"] == EVITWIN_ALGORITHM_ID
    assert result["summary"]["priority_compound_count"] == 1
    assert result["summary"]["display_rank_change_count"] == 1
    assert result["metadata"]["source_primary_is_v15"] is True
    feature = result["compounds"][0]
    assert feature["evitwin_winner"]["compound_code"] == "RIGHT"
    assert feature["decision"]["status"] == "qualified_authentic_rank"
    assert feature["decision"]["formal_level1_released"] is False
    assert source["entities"][0]["winner"]["candidate"]["compound_name"] == (
        "Wrong public winner"
    )


def test_evitwin_accepts_raw_peaks_from_v15_result_contract():
    rows = []
    for row in _rows():
        rows.append(
            {
                "collision_energy": row["collision_energy"],
                "raw_peaks": row["query_peaks"],
                "query_adduct": row["query_adduct"],
                "query_adduct_reliable": row["query_adduct_reliable"],
                "rt_method_id": row["rt_method_id"],
            }
        )
    result = build_evitwin_result(
        _v11(),
        rows,
        _library(),
        task_id="query-task",
        polarity="NEG",
        query_instrument_family="agilent-qtof-6546",
    )

    assert result["summary"]["priority_compound_count"] == 1
    assert result["compounds"][0]["evitwin_winner"]["compound_code"] == "RIGHT"


def test_evitwin_excludes_reference_built_from_the_query_task():
    result = build_evitwin_result(
        _v11(),
        _rows(),
        _library(source_task_id="query-task"),
        task_id="query-task",
        polarity="NEG",
        query_instrument_family="agilent-qtof-6546",
    )

    assert result["summary"]["self_reference_exclusion_count"] == 1
    assert result["compounds"] == []
    assert result["invariants"]["single_standard_self_match_allowed"] is False


def test_evitwin_never_imports_an_authentic_structure_outside_v15_candidate_pool():
    source = copy.deepcopy(_v11())
    source["entities"][0]["candidates"] = [
        source["entities"][0]["candidates"][0]
    ]
    source["entities"][0]["identity_resolution"] = {
        "status": "single_competitive_connectivity",
        "reportable_formula": "C10H10O4",
    }

    result = build_evitwin_result(
        source,
        _rows(),
        _library(),
        task_id="query-task",
        polarity="NEG",
        query_instrument_family="agilent-qtof-6546",
    )

    assert result["summary"]["display_rank_change_count"] == 0
    assert {
        score["connectivity"] for score in result["features"][0]["candidate_scores"]
    } == {"CCCCCCCCCCCCCC"}
    assert result["features"][0]["candidate_pool_audit"][
        "eligible_candidate_connectivities"
    ] == ["CCCCCCCCCCCCCC"]


def test_evitwin_can_use_verified_top40_support_without_changing_primary_v15():
    source = copy.deepcopy(_v11())
    source["entities"][0]["candidates"] = [
        source["entities"][0]["candidates"][0]
    ]
    source["entities"][0]["identity_resolution"] = {
        "status": "single_competitive_connectivity",
        "reportable_formula": "C10H10O4",
    }
    source["entities"][0]["candidate_recall_support"] = {
        "integration_module_id": "v15_candidate_recall_support",
        "expanded_candidates_can_replace_winner": False,
        "expanded_candidates_can_receive_2d_connectivity_decision": False,
        "production_identity_rank_unchanged": True,
        "expanded_candidates": [
            {
                "identity_key": "std2d:AAAAAAAAAAAAAA",
                "connectivity_identity": "AAAAAAAAAAAAAA",
                "candidate": {
                    "formula": "C10H10O4",
                    "inchikey": "AAAAAAAAAAAAAA-BBBBBBBBBB-C",
                },
                "can_replace_winner": False,
                "can_receive_2d_connectivity_decision": False,
                "requires_primary_ranking_or_orthogonal_evidence": True,
            }
        ],
    }
    original = copy.deepcopy(source)

    result = build_evitwin_result(
        source,
        _rows(),
        _library(),
        task_id="query-task",
        polarity="NEG",
        query_instrument_family="agilent-qtof-6546",
    )

    assert result["summary"]["display_rank_change_count"] == 1
    assert result["compounds"][0]["evitwin_winner"]["connectivity"] == (
        "AAAAAAAAAAAAAA"
    )
    assert source == original


def test_evitwin_abstains_for_cross_formula_or_unreliable_adduct():
    cross_formula = copy.deepcopy(_v11())
    cross_formula["entities"][0]["identity_resolution"] = {
        "status": "unresolved_cross_formula_event",
        "reportable_formula": None,
    }
    cross_formula_result = build_evitwin_result(
        cross_formula,
        _rows(),
        _library(),
        task_id="query-task",
        polarity="NEG",
        query_instrument_family="agilent-qtof-6546",
    )
    unreliable_rows = [
        {**row, "query_adduct_reliable": False} for row in _rows()
    ]
    unreliable_result = build_evitwin_result(
        _v11(),
        unreliable_rows,
        _library(),
        task_id="query-task",
        polarity="NEG",
        query_instrument_family="agilent-qtof-6546",
    )

    assert cross_formula_result["features"] == []
    assert cross_formula_result["applicability_audits"][0]["reason_codes"] == [
        "v15_identity_scope_is_unresolved_precursor_event"
    ]
    assert unreliable_result["features"] == []
    assert unreliable_result["summary"]["reliable_adduct_entity_count"] == 0
    assert unreliable_result["applicability_audits"][0]["reason_codes"] == [
        "query_adduct_missing_or_unreliable"
    ]


def test_evitwin_reuses_candidate_independent_entity_adduct_when_rows_lack_it():
    source = copy.deepcopy(_v11())
    source["entities"][0]["winner"]["chemistry_gate"] = {
        "status": "consistent",
        "formula": "C10H10O4",
        "adduct": "[M-H]-",
        "adduct_source": "ion_identity_network_adduct",
        "candidate_conditioned": False,
        "veto": False,
        "ion_identity_network_evidence_used": True,
        "ion_identity_network_evidence_available": True,
        "lcms_validated_adduct_group_evidence_available": False,
        "precursor_transformation_ambiguity": False,
        "consistent_metadata_overrides_isolated_conflicting_library_metadata": False,
        "status_counts": {"consistent": 2},
        "evaluated_collision_energy_bins": 2,
        "candidate_independent_checks": [
            {
                "source": "ion_identity_network_adduct",
                "adduct": "[M-H]-",
                "consistent": True,
            }
        ],
    }
    rows = [
        {**row, "query_adduct": "", "query_adduct_reliable": False}
        for row in _rows()
    ]
    result = build_evitwin_result(
        source,
        rows,
        _library(),
        task_id="query-task",
        polarity="NEG",
        query_instrument_family="agilent-qtof-6546",
    )

    validate_evitwin_result(result)
    assert result["summary"]["reliable_adduct_entity_count"] == 1
    assert result["summary"]["adduct_resolution_source_counts"] == {
        "v15_candidate_independent_chemistry_gate": 1
    }
    audit = result["features"][0]["query_adduct_audit"]
    assert audit["resolved_adduct"] == "[M-H]-"
    assert audit["resolution_source"] == (
        "v15_candidate_independent_chemistry_gate"
    )


def test_evitwin_rejects_candidate_conditioned_entity_adduct():
    source = copy.deepcopy(_v11())
    source["entities"][0]["winner"]["chemistry_gate"] = {
        "status": "consistent",
        "formula": "C10H10O4",
        "adduct": "[M-H]-",
        "adduct_source": "candidate_formula_mass_unique_inference",
        "candidate_conditioned": True,
        "veto": False,
        "precursor_transformation_ambiguity": False,
        "consistent_metadata_overrides_isolated_conflicting_library_metadata": False,
        "status_counts": {"consistent": 2},
        "evaluated_collision_energy_bins": 2,
        "candidate_independent_checks": [],
    }
    rows = [
        {**row, "query_adduct": "", "query_adduct_reliable": False}
        for row in _rows()
    ]
    result = build_evitwin_result(
        source,
        rows,
        _library(),
        task_id="query-task",
        polarity="NEG",
        query_instrument_family="agilent-qtof-6546",
    )

    assert result["features"] == []
    assert result["actual_ce_shadow"]["features"] == []
    assert result["summary"]["reliable_adduct_entity_count"] == 0


def test_evitwin_exposes_actual_ce_interpolation_as_non_ranking_shadow():
    source = copy.deepcopy(_v11())
    source["entities"][0]["member_result_indices"] = [0]
    row = {
        **_rows()[0],
        "collision_energy": 15.0,
    }
    result = build_evitwin_result(
        source,
        [row],
        _library(),
        task_id="query-task",
        polarity="NEG",
        query_instrument_family="agilent-qtof-6546",
    )

    validate_evitwin_result(result)
    shadow = result["actual_ce_shadow"]
    assert result["compounds"] == []
    assert result["features"] == []
    assert shadow["summary"]["scored_feature_count"] == 1
    assert shadow["summary"]["single_ce_feature_count"] == 1
    assert shadow["summary"]["display_rank_change_count"] == 0
    assert shadow["features"][0]["scalar_similarity_interpolation_used"] is True
    assert shadow["features"][0]["decision"]["display_rank_change_allowed"] is False


def test_exact_spectral_ties_cannot_be_broken_by_name_accession_or_library_order():
    tied = copy.deepcopy(_library())
    tied["entries"][1]["spectra"] = copy.deepcopy(tied["entries"][0]["spectra"])
    tied["entries"][0]["compound_name"] = "ZZZ name"
    tied["entries"][1]["compound_name"] = "AAA name"
    tied["entries"] = list(reversed(tied["entries"]))
    tied = finalize_anchor_library(tied)

    result = build_evitwin_result(
        _v11(),
        _rows(),
        tied,
        task_id="query-task",
        polarity="NEG",
        query_instrument_family="agilent-qtof-6546",
    )

    assert result["summary"]["priority_feature_count"] == 0
    assert result["summary"]["display_rank_change_count"] == 0
    assert result["features"][0]["decision"]["display_tier"] == "review"
    assert all(
        score["collision_energy_votes"] == 0
        for score in result["features"][0]["candidate_scores"]
    )
    assert result["invariants"]["compound_names_used_for_scoring"] is False
    assert result["invariants"]["database_accessions_used_for_scoring"] is False
    assert result["invariants"]["database_source_order_used_for_scoring"] is False


def test_evitwin_policy_forbids_unvalidated_level1_release():
    with pytest.raises(ValueError, match="same-method RT"):
        EviTwinPolicy(allow_formal_level1_release=True).validated()


def test_evitwin_assigns_separated_isomer_peaks_one_to_one_by_rt_order():
    base_policy = AuthenticAnchorPolicy().validated()
    entries = []
    for reference_id, code, connectivity, retention_time in (
        ("EARLY:NEG", "EARLY", "EEEEEEEEEEEEEE", 2.0),
        ("LATE:NEG", "LATE", "LLLLLLLLLLLLLL", 2.6),
    ):
        entries.append(
            {
                "reference_id": reference_id,
                "compound_code": code,
                "compound_name": f"{code.title()} isomer",
                "formula": "C10H10O4",
                "inchikey": f"{connectivity}-MMMMMMMMMM-N",
                "neutral_mass": 194.05791,
                "polarity": "NEG",
                "instrument_family": "agilent-qtof-6546",
                "method_id": "c18-test-method",
                "authentic_standard": True,
                "source_task_id": f"standard-{code.lower()}",
                "reference_qc": {"selected_adduct": "[M-H]-"},
                "spectra": [
                    {
                        "collision_energy_ev": 10.0,
                        "adduct": "[M-H]-",
                        "retention_time_min": retention_time,
                        "peaks": [
                            [80.0 + (code == "LATE"), 100.0],
                            [100.0 + (code == "LATE"), 60.0],
                            [120.0 + (code == "LATE"), 40.0],
                        ],
                    },
                    {
                        "collision_energy_ev": 20.0,
                        "adduct": "[M-H]-",
                        "retention_time_min": retention_time,
                        "peaks": [
                            [80.0 + (code == "LATE"), 80.0],
                            [100.0 + (code == "LATE"), 100.0],
                            [140.0 + (code == "LATE"), 30.0],
                        ],
                    },
                ],
            }
        )
    library = finalize_anchor_library(
        {
            "library_id": "rt-order-test",
            "runtime_mode": "shadow_only",
            "construction_policy": {
                "explicit_task_id_only": True,
                "candidate_generation_allowed": False,
                "authentic_library_shadow_search_allowed": True,
            },
            "release_policy": asdict(base_policy),
            "entries": entries,
        }
    )
    v11 = {
        "entities": [
            {
                "entity_id": "early-feature",
                "representative_precursor_mz": 193.0506335,
                "representative_retention_time": 4.0,
                "member_result_indices": [0, 1],
                "identity_resolution": {
                    "status": "unresolved_same_formula_isomer_group",
                    "reportable_formula": "C10H10O4",
                },
                "candidates": [
                    {
                        "candidate": {
                            "formula": "C10H10O4",
                            "inchikey": "EEEEEEEEEEEEEE-MMMMMMMMMM-N",
                        }
                    },
                    {
                        "candidate": {
                            "formula": "C10H10O4",
                            "inchikey": "LLLLLLLLLLLLLL-MMMMMMMMMM-N",
                        }
                    },
                ],
            },
            {
                "entity_id": "late-feature",
                "representative_precursor_mz": 193.0506335,
                "representative_retention_time": 5.0,
                "member_result_indices": [2, 3],
                "identity_resolution": {
                    "status": "unresolved_same_formula_isomer_group",
                    "reportable_formula": "C10H10O4",
                },
                "candidates": [
                    {
                        "candidate": {
                            "formula": "C10H10O4",
                            "inchikey": "EEEEEEEEEEEEEE-MMMMMMMMMM-N",
                        }
                    },
                    {
                        "candidate": {
                            "formula": "C10H10O4",
                            "inchikey": "LLLLLLLLLLLLLL-MMMMMMMMMM-N",
                        }
                    },
                ],
            },
        ]
    }
    rows = _rows() + [
        {
            **row,
            "query_peaks": [
                [float(peak[0]) + 1.0, peak[1]] for peak in row["query_peaks"]
            ],
        }
        for row in _rows()
    ]

    result = build_evitwin_result(
        v11,
        rows,
        library,
        task_id="mixture-query",
        polarity="NEG",
        query_instrument_family="agilent-qtof-6546",
        query_method_id="c18-test-method",
    )

    assert result["summary"]["retention_order_assignment_count"] == 2
    by_feature = {item["feature_id"]: item for item in result["features"]}
    assert by_feature["early-feature"]["evitwin_winner"]["compound_code"] == "EARLY"
    assert by_feature["late-feature"]["evitwin_winner"]["compound_code"] == "LATE"
    assert all(
        item["decision"]["status"] == "qualified_retention_order_rank"
        for item in by_feature.values()
    )
    assert all(
        item["retention_order_assignment"][
            "absolute_rt_used_for_cross_method_scoring"
        ]
        is False
        for item in by_feature.values()
    )

    mismatched = build_evitwin_result(
        v11,
        rows,
        library,
        task_id="mixture-query",
        polarity="NEG",
        query_instrument_family="agilent-qtof-6546",
        query_method_id="different-lc-method",
    )
    assert mismatched["summary"]["retention_order_assignment_count"] == 0
    assert all(
        item["retention_order_assignment"]["status"] == "not_applied"
        for item in mismatched["features"]
    )
    assert all(
        item["query_method_audit"]["reason_codes"]
        == ["task_and_query_row_lc_methods_conflict"]
        for item in mismatched["features"]
    )
