"""V7 regression tests for shadow metrics and explicit isomer-group evidence."""

from __future__ import annotations

import copy
import json

import pytest

from matching.honest_calibration import annotate
from matching.multimetric_shadow import (
    build_multimetric_shadow,
    build_peak_view_shadow,
)
from matching.search import _apply_isomer_fragment_rerank


IK_A = "AAAAAAAAAAAAAA-BBBBBBBBBB-C"
IK_B = "CCCCCCCCCCCCCC-DDDDDDDDDD-E"


def _isomer_candidate(
    name: str,
    inchikey: str,
    peaks: list[list[float]],
    *,
    adduct: str = "[M+H]+",
    collision_energy: float = 20.0,
) -> dict:
    return {
        "compound_name": name,
        "inchikey": inchikey,
        "formula": "C10H12O2",
        "smiles": "CCOc1ccccc1",
        "has_structure": True,
        "classification": "flavonoid",
        "classification_source": "curated",
        "source_tier": "experimental",
        "cosine": 0.95,
        "final_score": 0.70,
        "rank_score": 0.70,
        "reverse_match": 0.90,
        "mass_ppm": 2.0,
        "mass_metadata_consistent": True,
        "adduct": adduct,
        "query_adduct": "[M+H]+",
        "query_adduct_reliable": True,
        "lib_precursor_mz": 200.0,
        "lib_collision_energy": collision_energy,
        "lib_collision_energy_unit": "eV",
        "lib_instrument_type": "QTOF",
        "lib_peaks": peaks,
    }


def test_multimetric_shadow_is_json_safe_and_cannot_change_production_state():
    query = [[50.0, 20.0], [100.0, 100.0], [150.0, 40.0]]
    library = [[50.0, 10.0], [100.0, 50.0], [150.0, 20.0]]
    before_query = copy.deepcopy(query)
    before_library = copy.deepcopy(library)

    evidence = build_multimetric_shadow(
        query,
        library,
        query_precursor_mz=200.0,
        library_precursor_mz=200.0,
        production_modified_cosine=0.987,
        spectral_entropy=0.975,
        reverse_coverage=1.0,
    )

    assert evidence["status"] == "available"
    assert evidence["runtime_influence"] is False
    assert evidence["score_contribution"] == 0.0
    assert evidence["ranking_changed"] is False
    assert evidence["direct_cosine_hungarian"]["score"] == pytest.approx(1.0)
    assert evidence["modified_cosine_hungarian"]["score"] == pytest.approx(1.0)
    assert evidence["neutral_loss_cosine"]["score"] == pytest.approx(1.0)
    assert query == before_query
    assert library == before_library
    json.dumps(evidence)


def test_peak_view_shadow_exposes_low_abundance_support_without_influence():
    query = [[50.0 + i, 100.0 - i] for i in range(15)] + [[120.0, 4.0]]
    library = [[120.0, 100.0], [130.0, 20.0], [140.0, 10.0]]

    evidence = build_peak_view_shadow(
        query,
        library,
        query_precursor_mz=200.0,
        library_precursor_mz=200.0,
    )

    production = evidence["views"]["library_match"]
    sensitive = evidence["views"]["sensitive"]
    assert evidence["status"] == "available"
    assert evidence["runtime_influence"] is False
    assert evidence["score_contribution"] == 0.0
    assert evidence["ranking_changed"] is False
    assert production["query_peak_count"] == 15
    assert sensitive["query_peak_count"] == 16
    assert (
        sensitive["modified_cosine_hungarian"]["score"]
        > production["modified_cosine_hungarian"]["score"]
    )
    assert evidence["sensitive_gain_over_production_view"] > 0
    json.dumps(evidence)


def test_candidate_specific_product_ions_and_neutral_losses_are_explicit():
    candidates = [
        _isomer_candidate(
            "isomer A",
            IK_A,
            [[80.0, 0.2], [100.0, 1.0], [150.0, 0.5]],
        ),
        _isomer_candidate(
            "isomer B",
            IK_B,
            [[90.0, 0.2], [110.0, 1.0], [150.0, 0.5]],
        ),
    ]
    original_ranks = [candidate["rank_score"] for candidate in candidates]

    _apply_isomer_fragment_rerank(
        candidates,
        [[100.0, 1.0], [150.0, 0.5]],
        query_pmz=200.0,
        query_collision_energy={"value": 20.0, "unit": "eV"},
    )

    detail = candidates[0]["isomer_fragment_detail"]
    observed_products = [
        item for item in detail["candidate_specific_product_ions"]
        if item["observed"]
    ]
    observed_losses = [
        item for item in detail["candidate_specific_neutral_losses"]
        if item["observed"]
    ]
    assert detail["schema_version"] == 3
    assert detail["diagnostic_status"] == (
        "single_candidate_specific_feature_observed"
    )
    assert any(item["library_mz"] == 100.0 for item in observed_products)
    assert any(item["library_loss_da"] == 100.0 for item in observed_losses)
    assert detail["eligible_for_exact_identity_promotion"] is False
    assert detail["metadata_comparability"] == "incomplete_or_mismatched"
    assert detail["resolution_status"] == "single_selected_feature_insufficient"
    assert detail["opposing_supported_features"] == 0
    assert detail["diagnostic_mz_tolerance_da"] == pytest.approx(0.02)
    assert [candidate["rank_score"] for candidate in candidates] == original_ranks
    assert candidates[0]["isomer_fragment_rank_applied"] is False


def test_known_adduct_mismatch_blocks_candidate_specific_contrast():
    candidates = [
        _isomer_candidate("protonated", IK_A, [[100.0, 1.0]]),
        _isomer_candidate(
            "sodiated",
            IK_B,
            [[110.0, 1.0]],
            adduct="[M+Na]+",
        ),
    ]

    _apply_isomer_fragment_rerank(
        candidates,
        [[100.0, 1.0]],
        query_pmz=200.0,
        query_collision_energy={"value": 20.0, "unit": "eV"},
    )

    for candidate in candidates:
        detail = candidate["isomer_fragment_detail"]
        assert detail["diagnostic_status"] == (
            "candidate_metadata_not_comparable"
        )
        assert detail["comparable_candidate_count"] == 0
        assert detail["unique_total"] == 0
        assert candidate["isomer_fragment_bonus"] == 0.0


def test_isomer_group_high_is_distinct_from_exact_and_class_high():
    best = _isomer_candidate(
        "isomer A",
        IK_A,
        [[100.0, 1.0], [150.0, 0.5]],
    )
    competitor = _isomer_candidate(
        "isomer B",
        IK_B,
        [[110.0, 1.0], [150.0, 0.5]],
    )
    competitor.update({"cosine": 0.935, "final_score": 0.69, "rank_score": 0.69})
    best["isomer_fragment_detail"] = {
        "diagnostic_status": "multiple_candidate_specific_features_observed",
        "unique_supported": 2,
        "neutral_loss_supported": 2,
        "eligible_for_exact_identity_promotion": False,
    }
    best["top_candidates"] = [best.copy(), competitor]

    result = annotate(best)

    assert result["isomer_ambiguous"] is True
    assert result["exact_structure_confidence"] == "Medium"
    assert result["compound_confidence"] == "Medium"
    assert result["isomer_group_confidence"] == "High"
    assert result["class_confidence"] == "High"
    assert result["identification_scope"] == "isomer_group"
    assert result["isomer_group_formula"] == "C10H12O2"
    assert result["isomer_group_size"] == 2
    assert result["evidence_tier"] == "Isomer-group High"
    assert result["isomer_group_evidence"][
        "diagnostic_evidence_runtime_influence"
    ] is False


def test_nonambiguous_group_confidence_does_not_exceed_exact_confidence():
    best = _isomer_candidate(
        "single candidate",
        IK_A,
        [[100.0, 1.0], [150.0, 0.5]],
    )
    best["classification"] = "other"
    best["top_candidates"] = [best.copy()]

    result = annotate(best)

    assert result["isomer_ambiguous"] is False
    assert result["isomer_group_confidence"] == result[
        "exact_structure_confidence"
    ]


def test_fastapi_serialization_preserves_all_three_scopes_and_shadow():
    from web.server import _serialize_results

    best = _isomer_candidate(
        "isomer A",
        IK_A,
        [[100.0, 1.0], [150.0, 0.5]],
    )
    competitor = _isomer_candidate(
        "isomer B",
        IK_B,
        [[110.0, 1.0], [150.0, 0.5]],
    )
    competitor.update({"cosine": 0.94, "final_score": 0.69})
    # Let the production serializer derive the neutral formula from SMILES,
    # avoiding a deliberately inconsistent synthetic metadata record.
    best["formula"] = ""
    competitor["formula"] = ""
    shadow = build_multimetric_shadow(
        [[100.0, 1.0], [150.0, 0.5]],
        best["lib_peaks"],
        query_precursor_mz=200.0,
        library_precursor_mz=200.0,
        production_modified_cosine=best["cosine"],
        spectral_entropy=0.90,
        reverse_coverage=best["reverse_match"],
    )
    best["multimetric_shadow"] = shadow
    best["top_candidates"] = [
        {
            **best,
            "multimetric_shadow": shadow,
            "top_candidates": [],
        },
        competitor,
    ]
    serialized = _serialize_results([{
        "spectrum_id": "scan=1",
        "precursor_mz": 200.0,
        "query_peaks": [[100.0, 1.0], [150.0, 0.5]],
        "raw_peaks": [[100.0, 1.0], [150.0, 0.5]],
        "match": best,
    }])[0]["match"]

    assert serialized["exact_structure_confidence"] == "Medium"
    assert serialized["isomer_group_confidence"] == "High"
    assert serialized["confidence_class"] == "High"
    assert serialized["identification_scope"] == "isomer_group"
    assert serialized["isomer_ambiguity_type"] == (
        "constitutional_or_positional_isomer_group"
    )
    assert serialized["isomer_resolution_status"] == "unresolved"
    assert serialized["multimetric_shadow"]["runtime_influence"] is False
    assert serialized["top_candidates"][0][
        "multimetric_shadow"
    ]["ranking_changed"] is False
    json.dumps(serialized)


def test_exports_and_ui_name_the_scopes_and_shadow_metrics_explicitly():
    from web.report_common import REPORT_COLUMN_NAMES

    assert "Exact_Structure_Confidence" in REPORT_COLUMN_NAMES
    assert "Isomer_Group_Confidence" in REPORT_COLUMN_NAMES
    assert "Identification_Scope" in REPORT_COLUMN_NAMES
    assert "Isomer_Ambiguity_Type" in REPORT_COLUMN_NAMES
    assert "Isomer_Resolution_Status" in REPORT_COLUMN_NAMES
    assert "Isomer_Spectrally_Competitive_Candidates" in REPORT_COLUMN_NAMES
    assert "Isomer_Opposing_Features_Supported" in REPORT_COLUMN_NAMES
    assert "Isomer_All_Same_Formula_Opposing_Features" in REPORT_COLUMN_NAMES
    assert "Isomer_Reference_Replication" in REPORT_COLUMN_NAMES
    assert "Isomer_Diagnostic_Tolerance_Da" in REPORT_COLUMN_NAMES
    assert "Shadow_Modified_Cosine_Hungarian" in REPORT_COLUMN_NAMES
    with open("web/static/js/results.js", encoding="utf-8") as handle:
        script = handle.read()
    assert "异构体组 High" in script
    assert "仅类别 High" in script
    assert "candidate_specific_product_ions" in script
    assert "opposing_supported_features" in script
    assert "reference_replication" in script
    assert "achiral LC-MS/MS" in script
    assert "hasExplicitIsomerState" in script
    assert "no_spectrally_competitive_isomer" in script
    assert "neutral_loss_cosine" in script
