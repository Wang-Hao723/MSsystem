from __future__ import annotations

from pathlib import Path

import pytest

from web.app import _gen_rich
from web.raw_spectrum_locator import (
    format_rt_clock,
    raw_spectrum_locator,
    raw_spectrum_locator_text,
    retention_time_minutes,
)
from web.report import generate_html
from web.report_en import generate_english_html


def _match() -> dict:
    return {
        "compound_name": "Quercetin",
        "formula": "C15H10O7",
        "smiles": "O=c1cc(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12",
        "inchikey": "REFJWTPEDVJJIY-UHFFFAOYSA-N",
        "classification": "flavonoid",
        "dataset_source": "MassBank",
        "source_tier": "experimental",
        "cosine": 0.9132,
        "entropy": 0.8811,
        "reverse_match": 0.82,
        "mass_ppm": 2.1,
        "mass_score": 0.95,
        "fragment_score": 0.84,
        "intensity_score": 0.80,
        "final_score": 0.7712,
        "rank_score": 0.773,
        "confidence": "High",
        "confidence_class": "High",
        "evidence_tier": "High",
        "adduct": "[M-H]-",
    }


def _row() -> dict:
    return {
        "spectrum_id": "controllerType=0 controllerNumber=1 scan=123",
        "precursor_scan_ref": "controllerType=0 controllerNumber=1 scan=122",
        "precursor_mz": 301.05561,
        "retention_time": 315.0,
        "retention_time_unit": "seconds",
        "collision_energy": 35.0,
        "collision_energy_unit": "eV",
        "peaks_count": 2,
        "query_peaks": [[151.0037, 100.0], [179.0, 44.0]],
        "lib_peaks": [[151.0040, 88.0], [180.0, 30.0]],
        "match": _match(),
    }


def test_locator_normalizes_explicit_rt_units_without_guessing_scan_identity():
    row = _row()

    assert retention_time_minutes(row) == pytest.approx(5.25)
    assert format_rt_clock(5.25) == "05:15.0"
    locator = raw_spectrum_locator(row, "/private/input/blind_sample.mzML")
    assert locator["source_file"] == "blind_sample.mzML"
    assert locator["retention_time_min"] == pytest.approx(5.25)
    assert locator["spectrum_id"].endswith("scan=123")
    assert locator["precursor_scan_ref"].endswith("scan=122")
    text = raw_spectrum_locator_text(
        row, "/private/input/blind_sample.mzML", language="zh"
    )
    assert "RT 5.250 min (05:15.0)" in text
    assert "MS2 谱图编号/native ID" in text
    assert "scan=123" in text


def test_chinese_standard_and_detailed_html_place_rt_next_to_ms2_fragments(
    tmp_path: Path,
):
    row = _row()
    report_row = {**row, "matches": [row["match"]]}
    output = tmp_path / "standard.html"

    generate_html(
        [report_row],
        output,
        "blind_sample.mzML",
        high_confidence_scope_groups=8,
        high_confidence_exact_structures=0,
        high_confidence_isomer_groups=8,
    )
    standard = output.read_text(encoding="utf-8")
    detailed = _gen_rich(
        [report_row], source_name="blind_sample.mzML"
    )

    assert "原始 MS2 人工复核" in standard
    assert "High 级鉴定组" in standard
    assert "单一二维连接 0 · 未解异构体组 8" in standard
    assert "RT 5.250 min · MS2" in standard
    assert "scan=123" in standard
    assert "原始 MS2 定位信息" in detailed
    assert "RT 5.250 min (05:15.0)" in detailed
    assert "blind_sample.mzML" in detailed
    assert "scan=123" in detailed


def test_english_report_uses_the_same_observed_ms2_locator():
    rendered = generate_english_html(
        [_row()], "blind_sample.mzML", detailed=True
    )

    assert "Raw MS2 locator" in rendered
    assert "RT 5.250 min" in rendered
    assert "05:15.0" in rendered
    assert "spectrum/native ID" in rendered
    assert "scan=123" in rendered


def test_locator_separates_isolation_target_from_measured_ms1_search_mass():
    row = {
        **_row(),
        "search_precursor_mz": 301.05701,
        "search_precursor_mz_source": (
            "measured_continuous_ms1_feature_centroid"
        ),
        "search_precursor_delta_ppm": 4.65,
        "lcms_feature_used_for_search": True,
    }

    locator = raw_spectrum_locator(row, "blind_sample.mzML")
    rendered = raw_spectrum_locator_text(
        row, "blind_sample.mzML", language="zh"
    )

    assert locator["selected_ion_mz"] == pytest.approx(301.05561)
    assert locator["search_precursor_mz"] == pytest.approx(301.05701)
    assert "MS2 隔离目标 m/z 301.05561" in rendered
    assert "检索用实测 MS1 m/z 301.05701" in rendered

    detailed = _gen_rich(
        [{**row, "matches": [row["match"]]}],
        source_name="blind_sample.mzML",
    )
    assert "MS2 隔离目标 301.05561" in detailed
    assert "检索用实测 MS1 301.05701" in detailed
    assert "检索/计算质量（实测 MS1 特征质心）" in detailed


def test_integrated_report_exposes_locator_and_separate_review_index_download():
    root = Path(__file__).resolve().parents[1]
    script = (root / "web/static/js/report-generator.js").read_text(
        encoding="utf-8"
    )
    overview = (root / "web/static/js/sample-overview.js").read_text(
        encoding="utf-8"
    )
    index = (root / "web/static/index.html").read_text(encoding="utf-8")

    assert "如何回查原始 mzML" in script
    assert "最佳 MS2 原始定位" in script
    assert "retention_time_unit" in script
    assert "MS2_spectrum_native_ID" in script
    assert "downloadReviewIndex" in script
    assert "代表 MS2 原始定位" in overview
    assert "下载 MS2 人工复核索引（CSV）" in index
