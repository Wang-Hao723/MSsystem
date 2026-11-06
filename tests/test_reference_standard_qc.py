from __future__ import annotations

from copy import deepcopy
from pathlib import Path

import pytest

from matching.reference_standard_qc import (
    PROTON_MASS,
    _truth_rows,
    build_reference_standard_qc,
    list_reference_standard_sets,
    read_reference_standard_qc_artifact,
    validate_reference_standard_set_id,
    write_reference_standard_qc_artifact,
)
from web.reference_standard_qc_export import (
    reference_standard_qc_csv_text,
    reference_standard_qc_summary_rows,
    render_reference_standard_qc_html,
)
from web.server import AnalysisParams


def test_truth_sets_are_absent_from_ordinary_upload_and_analysis_contract():
    root = Path(__file__).resolve().parents[1]
    html = (root / "web/static/index.html").read_text(encoding="utf-8")
    upload_js = (root / "web/static/js/upload.js").read_text(encoding="utf-8")

    assert "referenceStandardSetId" not in html
    assert "混标／标准品质控真值集" not in html
    assert "reference_standard_set_id" not in upload_js
    assert "loadReferenceStandardSets" not in upload_js
    assert (
        "reference_standard_set_id"
        not in AnalysisParams.model_json_schema()["properties"]
    )


def test_eight_reference_standard_sets_match_curated_truth_counts():
    sets = list_reference_standard_sets()

    assert [(item["id"], item["truth_identity_count"]) for item in sets] == [
        ("F001", 20),
        ("F002", 9),
        ("F003", 8),
        ("F004", 9),
        ("F005", 11),
        ("F006", 10),
        ("F007", 9),
        ("F008", 12),
    ]
    assert all(item["selection_required"] for item in sets)
    assert all(item["used_for_scoring"] is False for item in sets)


def test_reference_standard_set_is_explicit_and_never_inferred_from_filename():
    assert validate_reference_standard_set_id(None) == ""
    assert validate_reference_standard_set_id("") == ""
    assert validate_reference_standard_set_id("f001") == "F001"
    with pytest.raises(ValueError, match="Unknown reference-standard QC set"):
        validate_reference_standard_set_id("F999")


def test_posthoc_qc_separates_declared_observable_retrieved_and_high_scope():
    truth = _truth_rows("F001")[0]
    precursor_mz = truth["monoisotopic_mass"] - PROTON_MASS
    row = {
        "spectrum_id": "scan-1",
        "precursor_mz": precursor_mz,
        "search_precursor_mz": precursor_mz,
        "search_precursor_mz_source": "measured_continuous_ms1_feature_centroid",
        "lcms_feature_used_for_search": True,
        "precursor_charge": 1,
        "ion_mode": "NEG",
        "match": {
            "compound_name": truth["name_en"],
            "inchikey": truth["inchikey"],
            "formula": truth["formula"],
            "confidence": "Medium",
            "exact_structure_confidence": "Medium",
            "isomer_group_confidence": "Medium",
            "isomer_ambiguous": False,
            "cosine": 0.88,
            "final_score": 0.70,
            "mass_ppm": 0.0,
            "top_candidates": [],
        },
    }
    original = deepcopy(row)

    payload = build_reference_standard_qc(
        "synthetic",
        [row],
        set_id="F001",
        polarity="NEG",
        filename="F001_NEG_1.mzML",
    )

    assert row == original
    summary = payload["summary"]
    assert summary["declared_truth_identities"] == 20
    assert summary["single_injection_compatible_precursor_identities"] == 1
    assert summary["observable_saved_topk_connectivity"] == 1
    assert summary["observable_top1_connectivity"] == 1
    assert summary["observable_reported_scope_supported"] == 1
    assert summary["observable_high_reportable_scope_supported"] == 0
    assert payload["used_for_candidate_generation"] is False
    assert payload["used_for_scoring"] is False
    assert payload["used_for_ranking"] is False
    assert payload["production_results_mutated"] is False


def test_reference_standard_qc_artifact_roundtrip_preserves_no_influence(tmp_path):
    payload = {
        "reference_standard_set": {"id": "F001"},
        "summary": {"declared_truth_identities": 20},
        "used_for_scoring": False,
        "used_for_ranking": False,
    }
    path = tmp_path / "task_reference_standard_qc.json"

    summary = write_reference_standard_qc_artifact(path, payload)
    restored = read_reference_standard_qc_artifact(path)

    assert summary["set_id"] == "F001"
    assert summary["used_for_scoring"] is False
    assert restored == payload


def test_candidate_on_incompatible_precursor_is_not_counted_as_recovered():
    truth = _truth_rows("F001")[0]
    row = {
        "spectrum_id": "unrelated-precursor",
        "precursor_mz": 123.456,
        "search_precursor_mz": 123.456,
        "precursor_charge": 1,
        "ion_mode": "NEG",
        "match": {
            "compound_name": truth["name_en"],
            "inchikey": truth["inchikey"],
            "formula": truth["formula"],
            "confidence": "High",
            "exact_structure_confidence": "High",
            "top_candidates": [],
        },
    }

    payload = build_reference_standard_qc(
        "synthetic-incompatible",
        [row],
        set_id="F001",
        polarity="NEG",
    )

    standard = next(
        item
        for item in payload["per_standard"]
        if item["inchikey"] == truth["inchikey"]
    )
    assert standard["observable_in_injection"] is False
    assert standard["saved_topk_connectivity"] is False
    assert standard["top1_connectivity_any_scan"] is False
    assert standard["high_reportable_scope_supported_any_scan"] is False


def test_report_projections_show_denominators_and_mark_unobservable_as_na():
    truth = _truth_rows("F001")[0]
    precursor_mz = truth["monoisotopic_mass"] - PROTON_MASS
    payload = build_reference_standard_qc(
        "synthetic-export",
        [{
            "spectrum_id": "scan-1",
            "precursor_mz": precursor_mz,
            "ion_mode": "NEG",
            "match": {
                "compound_name": truth["name_en"],
                "inchikey": truth["inchikey"],
                "formula": truth["formula"],
                "confidence": "Medium",
                "exact_structure_confidence": "Medium",
                "top_candidates": [],
            },
        }],
        set_id="F001",
        polarity="NEG",
    )

    summary_rows = reference_standard_qc_summary_rows(payload, language="zh")
    html = render_reference_standard_qc_html(payload, language="zh")
    csv_text = reference_standard_qc_csv_text(payload, language="zh")

    assert summary_rows[1]["结果"] == "1/20 (5.0%)"
    assert summary_rows[2]["结果"] == "1/1 (100.0%)"
    assert "1/20 (5.0%)" in html
    assert "不适用" in html
    assert "1/1 (100.0%)" in csv_text
    assert "不适用" in csv_text
    assert "不参与检索或排序" in csv_text
