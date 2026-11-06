from __future__ import annotations

import copy

import pytest

from matching.authentic_anchor import (
    AuthenticAnchorPolicy,
    evaluate_authentic_anchor,
    evaluate_authentic_anchor_actual_ce_shadow,
    finalize_anchor_library,
    search_authentic_anchor_feature,
    search_authentic_anchor_library,
    validate_anchor_library,
)


def _entry(reference_id, inchikey, peaks_a, peaks_b, *, instrument="agilent-qtof-6546"):
    return {
        "reference_id": reference_id,
        "compound_code": reference_id,
        "compound_name": reference_id,
        "formula": "C10H10O2",
        "inchikey": inchikey,
        "neutral_mass": 162.0681,
        "polarity": "NEG",
        "instrument_family": instrument,
        "method_id": "c18-t3-8min-neg",
        "column_mode": "achiral_c18",
        "authentic_standard": True,
        "source_task_id": f"task-{reference_id}",
        "spectra": [
            {"collision_energy_ev": 10.0, "peaks": peaks_a},
            {"collision_energy_ev": 20.0, "peaks": peaks_b},
        ],
    }


def _library(entries):
    policy = AuthenticAnchorPolicy().validated()
    return finalize_anchor_library({
        "library_id": "test-authentic-anchor",
        "runtime_mode": "shadow_only",
        "construction_policy": {
            "explicit_task_id_only": True,
            "candidate_generation_allowed": False,
            "authentic_library_shadow_search_allowed": True,
            "truth_labels_used_at_query_time": False,
            "single_standard_injections_are_reference_only": True,
        },
        "release_policy": policy.__dict__,
        "entries": entries,
    })


PEAKS_A10 = [[50.0, 100.0], [75.0, 60.0], [90.0, 30.0]]
PEAKS_A20 = [[55.0, 100.0], [80.0, 70.0], [95.0, 25.0]]
PEAKS_B10 = [[51.0, 100.0], [76.0, 60.0], [91.0, 30.0]]
PEAKS_B20 = [[56.0, 100.0], [81.0, 70.0], [96.0, 25.0]]


def _queries():
    return [
        {"collision_energy": 10.0, "query_peaks": [[50.005, 100], [75.005, 60], [90.005, 30]]},
        {"collision_energy": 20.0, "query_peaks": [[55.005, 100], [80.005, 70], [95.005, 25]]},
    ]


def test_multi_ce_anchor_supports_only_existing_distinct_connectivity_candidate():
    library = _library([
        _entry("A0001", "AAAAAAAAAAAAAA-BBBBBBBBBB-C", PEAKS_A10, PEAKS_A20),
        _entry("B0002", "CCCCCCCCCCCCCC-DDDDDDDDDD-E", PEAKS_B10, PEAKS_B20),
    ])
    result = evaluate_authentic_anchor(
        _queries(),
        library,
        candidate_identity_keys=["A0001", "B0002"],
        query_polarity="NEG",
        query_instrument_family="agilent-qtof-6546",
        baseline_separated=True,
    )

    assert result["decision"] == "anchor_connectivity_supported"
    assert result["winner_reference_id"] == "A0001"
    assert result["winner_collision_energy_votes"] == 2
    assert result["candidate_generation_used"] is False
    assert result["primary_rank_changed"] is False


def test_same_2d_stereoisomers_are_forced_to_group_on_achiral_msms():
    library = _library([
        _entry("STEREO-R", "AAAAAAAAAAAAAA-BBBBBBBBBB-C", PEAKS_A10, PEAKS_A20),
        _entry("STEREO-S", "AAAAAAAAAAAAAA-FFFFFFFFFF-G", PEAKS_B10, PEAKS_B20),
    ])
    result = evaluate_authentic_anchor(
        _queries(),
        library,
        candidate_identity_keys=["STEREO-R", "STEREO-S"],
        query_polarity="NEG",
        query_instrument_family="agilent-qtof-6546",
        baseline_separated=True,
    )

    assert result["decision"] == "isomer_group_only"
    assert result["same_2d_stereoisomer_set"] is True
    assert "same_2d_stereoisomers_not_resolved_by_achiral_msms" in result["reason_codes"]


def test_coelution_and_single_collision_energy_fail_closed():
    library = _library([
        _entry("A0001", "AAAAAAAAAAAAAA-BBBBBBBBBB-C", PEAKS_A10, PEAKS_A20),
        _entry("B0002", "CCCCCCCCCCCCCC-DDDDDDDDDD-E", PEAKS_B10, PEAKS_B20),
    ])
    coeluted = evaluate_authentic_anchor(
        _queries(), library,
        candidate_identity_keys=["A0001", "B0002"],
        query_polarity="NEG",
        query_instrument_family="agilent-qtof-6546",
        baseline_separated=False,
    )
    one_ce = evaluate_authentic_anchor(
        _queries()[:1], library,
        candidate_identity_keys=["A0001", "B0002"],
        query_polarity="NEG",
        query_instrument_family="agilent-qtof-6546",
        baseline_separated=True,
    )

    assert coeluted["decision"] == "isomer_group_only"
    assert "chromatography_not_proven_baseline_separated" in coeluted["reason_codes"]
    assert one_ce["decision"] == "abstained"
    assert "insufficient_collision_energy_coverage" in one_ce["reason_codes"]


def test_actual_ce_scalar_interpolation_is_single_vote_shadow_only():
    library = _library([
        _entry("A0001", "AAAAAAAAAAAAAA-BBBBBBBBBB-C", PEAKS_A10, PEAKS_A20),
        _entry("B0002", "CCCCCCCCCCCCCC-DDDDDDDDDD-E", PEAKS_B10, PEAKS_B20),
    ])
    result = evaluate_authentic_anchor_actual_ce_shadow(
        [{"collision_energy": 15.0, "query_peaks": _queries()[0]["query_peaks"]}],
        library,
        candidate_identity_keys=["A0001", "B0002"],
        query_polarity="NEG",
        query_instrument_family="agilent-qtof-6546",
    )

    assert result["decision"] == "single_ce_shadow_candidate_group"
    assert result["query_collision_energy_condition_count"] == 1
    assert result["winner_collision_energy_votes"] == 1
    assert result["scalar_similarity_interpolation_used"] is True
    assert result["outside_reference_range_clamped"] is False
    assert result["evitwin_display_rank_changed"] is False
    assert result["ranking_use_allowed"] is False


def test_actual_ce_shadow_does_not_extrapolate_outside_reference_range():
    library = _library([
        _entry("A0001", "AAAAAAAAAAAAAA-BBBBBBBBBB-C", PEAKS_A10, PEAKS_A20),
        _entry("B0002", "CCCCCCCCCCCCCC-DDDDDDDDDD-E", PEAKS_B10, PEAKS_B20),
    ])
    result = evaluate_authentic_anchor_actual_ce_shadow(
        [{"collision_energy": 25.0, "query_peaks": _queries()[1]["query_peaks"]}],
        library,
        candidate_identity_keys=["A0001", "B0002"],
        query_polarity="NEG",
        query_instrument_family="agilent-qtof-6546",
    )

    assert result["decision"] == "abstained"
    assert result["candidate_scores"] == []
    assert result["outside_reference_range_clamped"] is False


def test_multi_candidate_scoring_requires_shared_query_ce_coverage():
    complete = _entry(
        "A0001", "AAAAAAAAAAAAAA-BBBBBBBBBB-C", PEAKS_A10, PEAKS_A20
    )
    incomplete = _entry(
        "B0002", "CCCCCCCCCCCCCC-DDDDDDDDDD-E", PEAKS_B10, PEAKS_B20
    )
    incomplete["spectra"] = [
        {"collision_energy_ev": 30.0, "peaks": PEAKS_B10},
        {"collision_energy_ev": 40.0, "peaks": PEAKS_B20},
    ]
    result = evaluate_authentic_anchor(
        _queries(),
        _library([complete, incomplete]),
        candidate_identity_keys=["A0001", "B0002"],
        query_polarity="NEG",
        query_instrument_family="agilent-qtof-6546",
        baseline_separated=True,
    )

    assert result["decision"] == "abstained"
    assert result["shared_candidate_collision_energy_conditions"] == []
    assert all(
        score["collision_energy_count"] == 0
        for score in result["candidate_scores"]
    )
    assert sorted(
        bool(score["available_collision_energy_scores"])
        for score in result["candidate_scores"]
    ) == [False, True]


def test_actual_ce_shadow_rejects_a_partially_covered_candidate_pair():
    complete = _entry(
        "A0001", "AAAAAAAAAAAAAA-BBBBBBBBBB-C", PEAKS_A10, PEAKS_A20
    )
    outside = _entry(
        "B0002", "CCCCCCCCCCCCCC-DDDDDDDDDD-E", PEAKS_B10, PEAKS_B20
    )
    outside["spectra"] = [
        {"collision_energy_ev": 20.0, "peaks": PEAKS_B10},
        {"collision_energy_ev": 40.0, "peaks": PEAKS_B20},
    ]
    result = evaluate_authentic_anchor_actual_ce_shadow(
        [{"collision_energy": 15.0, "query_peaks": _queries()[0]["query_peaks"]}],
        _library([complete, outside]),
        candidate_identity_keys=["A0001", "B0002"],
        query_polarity="NEG",
        query_instrument_family="agilent-qtof-6546",
    )

    assert result["decision"] == "abstained"
    assert result["shared_candidate_collision_energy_conditions"] == []
    assert result["candidate_scores"] == []


def test_candidate_set_is_mandatory_and_instrument_mismatch_cannot_release():
    library = _library([
        _entry("A0001", "AAAAAAAAAAAAAA-BBBBBBBBBB-C", PEAKS_A10, PEAKS_A20),
        _entry("B0002", "CCCCCCCCCCCCCC-DDDDDDDDDD-E", PEAKS_B10, PEAKS_B20),
    ])
    with pytest.raises(ValueError, match="existing candidate set"):
        evaluate_authentic_anchor(
            _queries(), library,
            candidate_identity_keys=[],
            query_polarity="NEG",
            query_instrument_family="agilent-qtof-6546",
            baseline_separated=True,
        )
    mismatch = evaluate_authentic_anchor(
        _queries(), library,
        candidate_identity_keys=["A0001", "B0002"],
        query_polarity="NEG",
        query_instrument_family="orbitrap",
        baseline_separated=True,
    )
    assert mismatch["decision"] == "abstained"
    assert mismatch["instrument_mismatch_count"] == 2


def test_content_hash_detects_tampering():
    library = _library([
        _entry("A0001", "AAAAAAAAAAAAAA-BBBBBBBBBB-C", PEAKS_A10, PEAKS_A20),
    ])
    tampered = copy.deepcopy(library)
    tampered["entries"][0]["compound_name"] = "changed after freeze"

    with pytest.raises(ValueError, match="content hash mismatch"):
        validate_anchor_library(tampered)


def test_mass_filtered_anchor_search_can_recall_without_changing_primary_pool():
    library = _library([
        _entry("A0001", "AAAAAAAAAAAAAA-BBBBBBBBBB-C", PEAKS_A10, PEAKS_A20),
        _entry("B0002", "CCCCCCCCCCCCCC-DDDDDDDDDD-E", PEAKS_B10, PEAKS_B20),
    ])
    result = search_authentic_anchor_library(
        _queries(),
        library,
        precursor_mz=161.060823533379,
        query_polarity="NEG",
        query_instrument_family="agilent-qtof-6546",
        baseline_separated=True,
    )

    assert result["decision"] == "anchor_connectivity_supported"
    assert result["winner_reference_id"] == "A0001"
    assert result["mass_compatible_anchor_count"] == 2
    assert result["candidate_generation_used"] is True
    assert result["primary_candidate_pool_changed"] is False


def test_feature_search_never_averages_different_adduct_precursors():
    library = _library([
        _entry("A0001", "AAAAAAAAAAAAAA-BBBBBBBBBB-C", PEAKS_A10, PEAKS_A20),
        _entry("B0002", "CCCCCCCCCCCCCC-DDDDDDDDDD-E", PEAKS_B10, PEAKS_B20),
    ])
    queries = [
        {**_queries()[0], "search_precursor_mz": 161.060823533379},
        {**_queries()[1], "search_precursor_mz": 207.066303533379},
    ]
    # Give both precursor/adduct groups full CE coverage.  A naive median of
    # 161 and 207 would be a fictitious ion at m/z 184 and find no anchor.
    queries.extend([
        {**_queries()[1], "search_precursor_mz": 161.060823533379},
        {**_queries()[0], "search_precursor_mz": 207.066303533379},
    ])
    result = search_authentic_anchor_feature(
        queries,
        library,
        query_polarity="NEG",
        query_instrument_family="agilent-qtof-6546",
        baseline_separated=True,
    )

    assert result["decision"] == "anchor_connectivity_supported"
    assert result["winner_reference_id"] == "A0001"
    assert result["precursor_group_count"] == 2
    assert result["selected_precursor_group_mz"] in {
        161.060823533379, 207.066303533379,
    }
