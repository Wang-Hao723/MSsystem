from __future__ import annotations

import csv
import hashlib
import html
import io
import zipfile

from fastapi.testclient import TestClient
from openpyxl import load_workbook

from matching.expected_compounds import write_expected_compounds_artifact
from web.expected_compounds_export import (
    EXPECTED_FEATURE_COLUMNS,
    append_expected_csv_section,
    expected_feature_groups,
    expected_feature_rows,
    expected_raw_feature_rows,
    expected_summary_row,
    inject_expected_html,
    render_expected_html_report,
)
from web.report_common import (
    MS2_READER_COLUMN_NAMES,
    RESULTS_COLUMN_NAMES,
    row_to_dict,
)


def _payload() -> dict:
    return {
        "schema_version": 1,
        "task_id": "expected-export",
        "method": "Expected Compounds / continuous-MS1 targeted XIC",
        "status": "ok",
        "requested": True,
        "catalog": {
            "target_rows": 12,
            "expected_ion_hypotheses": 24,
            "dataset_sources": ["literature:test"],
            "quality_audit": {
                "active_target_rows": 10,
                "quarantined_target_rows": 2,
                "quarantine_reason_counts": {
                    "invalid_nonempty_molecular_formula": 2,
                },
            },
        },
        "summary": {
            "status": "ok",
            "target_catalog_rows": 12,
            "active_target_rows": 10,
            "quarantined_target_rows": 2,
            "expected_ion_hypotheses": 24,
            "detected_feature_count": 1,
            "supported_features": 1,
            "mass_only_features": 0,
            "ambiguous_feature_count": 0,
            "candidate_target_rows": 1,
            "candidate_target_names": 1,
            "ms1_scan_count": 100,
            "usable_ms1_scan_count": 99,
            "evidence_counts": {
                "MS1_EXACT_MASS_CANDIDATE": 0,
                "MS1_ISOTOPE_SUPPORTED": 1,
                "MS1_MS2_DIAGNOSTIC_PROBABLE": 0,
                "EXPERIMENTAL_LIBRARY_CONFIRMED": 0,
            },
            "counted_in_standard_spectral_library_matches": False,
            "standard_spectral_library_match_increment": 0,
        },
        "features": [{
            "feature_id": "ECI_1_P001",
            "candidate_name": "Glucose",
            "candidate_names": ["Glucose"],
            "formula": "C6H12O6",
            "adduct": "[M+H]+",
            "polarity": "POS",
            "target_mz": 181.07066457,
            "observed_apex_mz": 181.07070000,
            "mass_error_ppm": 0.1957,
            "apex_rt_minutes": 1.25,
            "boundary_start_rt_minutes": 1.15,
            "boundary_end_rt_minutes": 1.35,
            "fwhm_seconds": 6.0,
            "apex_intensity": 6000.0,
            "xic_area_intensity_seconds": 25000.0,
            "signal_to_noise": 12.0,
            "scan_points": 7,
            "xic_nonzero_points": 7,
            "evidence_level": "MS1_ISOTOPE_SUPPORTED",
            "cd_pattern_accepted": True,
            "isomer_ambiguous": False,
            "identity_count": 1,
            "identification_scope": "expected_compound_candidate",
            "formula_isotope_evidence": {"status": "supporting"},
            "associated_ms2": [],
            "ms2_diagnostic_support": [],
            "experimental_library_support": [],
            "target_members": [{
                "target_id": "t1",
                "inchikey": "WQZGKKKJIJFFOK-GASJEMHNSA-N",
            }],
            "dataset_sources": ["literature:test"],
            "generation_sources": ["formula_exact_mass"],
            "evidence_interpretation": (
                "Exact-mass feature with a supporting isotope envelope."
            ),
            "ordinary_library_match_counted": False,
        }],
        "scientific_boundaries": {
            "ordinary_library_matched_count_unchanged": True,
        },
    }


def test_expected_export_rows_keep_independent_ms1_contract():
    payload = _payload()
    summary = expected_summary_row(payload)
    rows = expected_feature_rows(payload)

    assert summary["Detected_Features"] == 1
    assert summary["Supported_Features"] == 1
    assert summary["Active_Target_Rows"] == 10
    assert summary["Quarantined_Target_Rows"] == 2
    assert summary["Ordinary_MS2_Match_Increment"] == 0
    assert len(rows) == 1
    assert list(rows[0]) == EXPECTED_FEATURE_COLUMNS
    assert rows[0]["Candidate_Name"] == "Glucose"
    assert rows[0]["Candidate_Count"] == 1
    assert rows[0]["Target_Database"] == "test"
    assert rows[0]["Evidence"] == "MS1 + 同位素支持"
    assert rows[0]["Background_Status"] == ""
    assert not any("Isomer" in column for column in EXPECTED_FEATURE_COLUMNS)
    assert not any(isinstance(value, bool) for value in rows[0].values())
    groups = expected_feature_groups(payload)
    assert groups["high"] == []
    assert groups["isotope"][0]["Candidate_Name"] == "Glucose"
    assert groups["mass_only"] == []


def test_standard_targeted_row_compacts_long_candidate_groups():
    payload = _payload()
    payload["features"][0]["candidate_names"] = [
        "Glucose",
        "D-Glucose",
        "L-Glucose",
    ]

    row = expected_feature_rows(payload)[0]

    assert row["Candidate_Name"] == "Glucose"
    assert row["Candidate_Count"] == 1


def test_standard_targeted_table_groups_peaks_but_audit_keeps_each_segment():
    payload = _payload()
    second = dict(payload["features"][0])
    second.update({
        "feature_id": "ECI_1_P002",
        "adduct": "[M+Na]+",
        "target_mz": 203.052609,
        "observed_apex_mz": 203.052620,
        "apex_rt_minutes": 1.27,
        "xic_area_intensity_seconds": 5000.0,
        "baseline_corrected_area_intensity_seconds": 4500.0,
        "apex_intensity": 3000.0,
        "evidence_level": "MS1_EXACT_MASS_CANDIDATE",
        "peak_quality_score": 7.5,
    })
    payload["features"].append(second)
    payload["summary"]["detected_feature_count"] = 2

    standard = expected_feature_rows(payload)
    audit = expected_raw_feature_rows(payload)

    assert len(standard) == 1
    assert len(audit) == 2
    assert standard[0]["Raw_Peak_Count"] == 2
    assert standard[0]["Related_Adducts"] == "[M+H]+ | [M+Na]+"
    assert standard[0]["Peak_Area"] == 25000.0
    assert {row["Feature_ID"] for row in audit} == {
        "ECI_1_P001", "ECI_1_P002"
    }


def test_resolved_support_does_not_split_the_parent_ambiguity_group():
    payload = _payload()
    members = [
        {
            "target_id": "t1",
            "identity_key": "ik::AAAAAAAAAAAAAA-BBBBBBBBBB-C",
            "compound_name": "Candidate A",
        },
        {
            "target_id": "t2",
            "identity_key": "ik::CCCCCCCCCCCCCC-DDDDDDDDDD-E",
            "compound_name": "Candidate B",
        },
    ]
    payload["features"][0].update({
        "candidate_name": "Candidate A / Candidate B",
        "candidate_names": ["Candidate A", "Candidate B"],
        "target_members": members,
        "identity_count": 2,
        "evidence_level": "MS1_EXACT_MASS_CANDIDATE",
        "cd_pattern_accepted": False,
        "formula_isotope_evidence": {},
    })
    supported = dict(payload["features"][0])
    supported.update({
        "feature_id": "ECI_1_P002",
        "apex_rt_minutes": 2.5,
        "evidence_level": "MS1_MS2_DIAGNOSTIC_PROBABLE",
        "cd_pattern_accepted": True,
        "ms2_diagnostic_support": [{
            "spectrum_id": "scan-2",
            "target_id": "t1",
            "matched_diagnostic_ions": [101.0],
        }],
    })
    payload["features"].append(supported)
    payload["summary"]["detected_feature_count"] = 2

    rows = expected_feature_rows(payload)
    summary = expected_summary_row(payload)

    assert len(rows) == 1
    assert rows[0]["Candidate_Count"] == 2
    assert rows[0]["Raw_Peak_Count"] == 2
    assert rows[0]["Candidate_Name"] == "Candidate A"
    assert summary["Compound_Summary_Rows"] == 1
    assert summary["Supported_Compounds"] == 1


def test_standard_ms2_row_hides_library_isomer_metadata_and_alias_lists():
    raw_name = (
        "MassBank:PR101019 Kaempferol-3-Glucoside"
        "|Kaempferol alias (known isomers: 0; isobaric peaks: 0)"
    )
    result = {
        "spectrum_id": "scan-1",
        "match": {"compound_name": raw_name},
    }

    assert row_to_dict(result, report_mode=False)["Compound_Name"] == (
        "Kaempferol-3-Glucoside"
    )
    assert row_to_dict(result, report_mode=True)["Compound_Name"] == raw_name
    result["match"]["compound_name"] = (
        "3,4-diaminobenzoic acid "
        "(known isomers: 0; isobaric peaks: 0)"
    )
    assert row_to_dict(result, report_mode=False)["Compound_Name"] == (
        "3,4-diaminobenzoic acid"
    )


def test_standard_ms2_row_restores_complete_acquisition_fragment_and_class_fields():
    result = {
        "spectrum_id": "334243",
        "precursor_mz": 433.112091064453,
        "retention_time": 5.570583333333,
        "precursor_charge": 1,
        "collision_energy": 35.0,
        "collision_energy_unit": "eV",
        "instrument_type": "qtof",
        "instrument": "Agilent 6545XT",
        "ion_mode": "NEG",
        "peaks_count": 1305,
        "query_peaks": [
            [119.04992902463847, 6161.962890625],
            [151.00279386177098, 17701.708984375],
            [271.0600820817677, 14616.34765625],
        ],
        "raw_peaks": [
            [151.00279386177098, 17701.708984375],
            [271.0600820817677, 14616.34765625],
            [119.04992902463847, 6161.962890625],
        ],
        "match": {
            "compound_name": "Prunin, Naringenin-7-O-glucoside",
            "source_compound_name": (
                '"Prunin, Naringenin-7-O-glucoside '
                'CollisionEnergy:102040" M-H'
            ),
            "formula": "C21H22O10",
            "adduct": "",
            "inchikey": "DLIKSSGEMUFQOK-UHFFFAOYSA-N",
            "smiles": (
                "O=C1CC(c2ccc(O)cc2)Oc2cc("
                "OC3OC(CO)C(O)C(O)C3O)cc(O)c21"
            ),
            "confidence": "Medium",
            "compound_confidence": "Medium",
            "evidence_tier": "Isomer-group High",
            "final_score": 0.9156,
            "cosine": 0.9912,
            "mass_ppm": 4.45,
            "fragment_count": 3,
            "ms_ms_fragments": "119.0499; 151.0028; 271.0601",
            "classification": "flavonoid",
            "classification_source": "structure_rule",
            "classification_confidence": 0.97,
            "dataset_source": "Flavonoid",
            "source_tier": "experimental",
            "lib_spectrum_id": "Flavonoid_2127269",
            "lib_collision_energy": "",
            "lib_collision_energy_raw": "",
            "lib_collision_energy_unit": "",
            "lib_instrument_type": "",
            "lib_instrument": "",
            "query_instrument_type": "qtof",
            "query_instrument": "Agilent 6545XT",
        },
    }

    row = row_to_dict(result, report_mode=False)

    assert list(row) == RESULTS_COLUMN_NAMES
    assert row["Compound_Name"] == "Naringenin 7-O-glucoside"
    assert row["Adduct"] == "[M-H]-"
    assert row["Library_Reported_Adduct"] == ""
    assert row["Adduct_Source"] == "candidate_formula_precursor_unique_10ppm"
    assert 0 < row["Adduct_Inference_Error_PPM"] < 10
    assert row["Ion_Mode"] == "NEG"
    assert abs(row["Search_Precursor_mz"] - result["precursor_mz"]) < 1e-4
    assert row["MS1_Feature_Used_For_Search"] is False
    assert row["Query_Collision_Energy"] == 35
    assert row["Query_Collision_Energy_Unit"] == "eV"
    assert row["Query_Instrument_Type"] == "qtof"
    assert row["Query_Instrument"] == "Agilent 6545XT"
    assert row["Fragment_Count"] == 3
    assert "151.00279:17701.7" in row["Matched_Fragment_Details"]
    assert "151.00279:17702" in row["MS_MS_Spectrum"]
    assert row["Classification"] == "flavonoid"
    assert row["Classification_Source"] == "structure_rule"
    assert row["Classification_Confidence"] == 0.97
    assert row["Library_Spectrum"] == "Flavonoid_2127269"
    assert row["Library_Collision_Energy"] is None
    assert row["Library_CE_Metadata_Status"] == ""
    assert row["Final_Claim_Released"] == ""  # no claim contract on this fixture
    assert isinstance(row["Classification_Reliable"], bool)
    assert isinstance(row["Candidate_Pathway_Eligible_For_Fusion"], bool)
    assert isinstance(row["Candidate_Pathway_Runtime_Influence"], bool)
    assert "Isomer_Best_Candidate" in RESULTS_COLUMN_NAMES
    assert "Isomer_Competing_Candidates" in RESULTS_COLUMN_NAMES
    assert "Isomer_Candidate_Details_JSON" in RESULTS_COLUMN_NAMES
    assert not any(
        token in column.lower()
        for column in RESULTS_COLUMN_NAMES
        for token in ("upload_time", "analysis_time", "blank_")
    )


def test_ms2_csv_complete_schema_quotes_cells_and_contains_no_boolean_dump():
    from web.server import _ms2_csv_text

    result = {
        "spectrum_id": "scan-1",
        "precursor_mz": 154.0504,
        "retention_time": 1.25,
        "precursor_charge": 1,
        "collision_energy": 20,
        "collision_energy_unit": "NCE",
        "instrument_type": "Orbitrap",
        "ion_mode": "POS",
        "query_peaks": [[80.0, 100.0], [136.0398, 500.0]],
        "raw_peaks": [[136.0398, 500.0], [80.0, 100.0]],
        "peaks_count": 2,
        "match": {
            "compound_name": "3,4-diaminobenzoic acid",
            "formula": "C7H8N2O2",
            "adduct": "[M+H]+",
            "smiles": "N=C(N)c1ccc(C(=O)O)cc1",
            "inchikey": "ABCDEFGHIJKLMN-ABCDEFGHIJ-A",
            "confidence": "Medium",
            "compound_confidence": "Medium",
            "evidence_tier": "Isomer-group High",
            "final_score": 0.8,
            "cosine": 0.9,
            "mass_ppm": 2.0,
            "fragment_count": 1,
            "ms_ms_fragments": "136.0398",
            "classification": "benzenoid",
            "classification_source": "library",
            "classification_confidence": 0.9,
            "dataset_source": "MassBank",
            "source_tier": "experimental",
            "lib_spectrum_id": "MB1",
            "lib_collision_energy": 18,
            "lib_collision_energy_unit": "eV",
            "lib_instrument_type": "QTOF",
            "lib_instrument": "Model X",
        },
    }

    text = _ms2_csv_text([result])
    rows = list(csv.DictReader(io.StringIO(text)))

    assert len(rows) == 1
    assert len(rows[0]) == len(MS2_READER_COLUMN_NAMES)
    assert rows[0]["鉴定组"] == "G001"
    assert rows[0]["候选角色"] == "单一候选"
    assert rows[0]["组内候选数"] == "1"
    assert rows[0]["候选化合物"] == "3,4-diaminobenzoic acid"
    assert rows[0]["MS2 隔离目标 m/z"] == "154.0504"
    assert rows[0]["检索用实测/回退 m/z"] == "154.0504"
    assert rows[0]["样本 CE"] == "20"
    assert rows[0]["样本 CE 单位"] == "NCE"
    assert rows[0]["样品离子形式依据"] == "谱库候选离子形式（非样品独立测定）"
    assert rows[0]["代表馆藏 CE"] == "18"
    assert rows[0]["代表馆藏 CE 单位"] == "eV"
    assert rows[0]["代表谱匹配碎片 (m/z:强度)"] == "136.03980:500"
    assert rows[0]["原始 MS2 碎片离子（已保存输入峰，m/z:强度）"] == (
        "136.03980:500 80.00000:100"
    )
    assert rows[0]["候选化学类别"] == "benzenoid"
    assert "鉴定组 ID" not in rows[0]
    assert "候选记录状态" not in rows[0]
    assert "True" not in text
    assert "False" not in text
    assert "Isomer-group" not in text


def test_ms2_csv_preserves_search_and_raw_peaks_with_explicit_cell_limit():
    """CSV previews name a complete peak export that retains all saved values."""
    import gzip
    import hashlib
    import json
    from web.lossless_peak_export import canonical_peak_json, iter_raw_peak_jsonl_gzip
    from web.server import _ms2_csv_text

    raw_peaks = [
        [50.0 + index * 0.01, 1000.0 + index]
        for index in range(4885)
    ]
    query_peaks = [
        [100.0 + index, 5000.0 - index * 100]
        for index in range(15)
    ]
    result = {
        "spectrum_id": "profile-scan",
        "precursor_mz": 289.0715,
        "retention_time": 7.5396,
        "precursor_charge": 1,
        "collision_energy": 35,
        "collision_energy_unit": "eV",
        "instrument_type": "qtof",
        "ion_mode": "NEG",
        "query_peaks": query_peaks,
        "raw_peaks": raw_peaks,
        "peaks_count": len(raw_peaks),
        "match": {
            "compound_name": "EC (C)",
            "formula": "C15H14O6",
            "adduct": "[M-H]-",
            "smiles": "Oc1cc(O)c2C(C(O)COc2c1)c1ccc(O)c(O)c1",
            "inchikey": "AAAAAAAAAAAAAA-BBBBBBBBBB-C",
            "confidence": "Medium",
            "final_score": 0.8,
            "cosine": 0.9,
            "mass_ppm": 1.0,
            "fragment_count": 1,
            "ms_ms_fragments": "100.0000",
            "classification": "flavonoid",
            "classification_source": "structure_rule",
            "classification_confidence": 0.9,
            "compound_role": "unclassified",
            "potential_contaminant": False,
            "natural_product_candidate": False,
            "dataset_source": "MassBank",
            "source_tier": "experimental",
            "lib_spectrum_id": "MB-profile",
        },
    }

    text = _ms2_csv_text([result], task_id="profile-export", results=[result])
    parsed = list(csv.reader(io.StringIO(text)))
    rows = list(csv.DictReader(io.StringIO(text)))

    assert len(parsed) == 2
    assert all(len(row) == len(MS2_READER_COLUMN_NAMES) for row in parsed)
    assert rows[0]["实际检索 MS2 峰数"] == "15"
    assert len(rows[0]["样本 MS2 碎片离子（实际检索谱，m/z:强度）"].split()) == 15
    assert rows[0]["原始 MS2 输入峰数"] == "4885"
    assert "50.00000:1000" in rows[0][
        "原始 MS2 碎片离子（已保存输入峰，m/z:强度）"
    ]
    assert "此处为已保存峰的节选" in rows[0][
        "原始 MS2 碎片离子（已保存输入峰，m/z:强度）"
    ]
    assert "14_原始峰索引" in rows[0]["备注"]
    assert "15_原始峰JSON" in rows[0]["备注"]
    assert "完整原始峰 JSONL.gz" in rows[0]["备注"]
    assert "/api/tasks/profile-export/download/raw-peaks" in rows[0]["备注"]
    assert max(len(cell) for row in parsed for cell in row) <= 32000
    archive = b"".join(iter_raw_peak_jsonl_gzip([result], "profile-export"))
    complete_records = [json.loads(line) for line in gzip.decompress(archive).splitlines()]
    assert len(complete_records) == 1
    saved = complete_records[0]
    assert saved["spectrum_id"] == "profile-scan"
    assert saved["result_index"] == 0
    assert saved["peak_count"] == 4885
    assert saved["raw_peaks"] == raw_peaks
    assert saved["raw_peaks_sha256"] == hashlib.sha256(
        canonical_peak_json(raw_peaks).encode("ascii")
    ).hexdigest()
    assert result["query_peaks"] == query_peaks
    assert result["raw_peaks"] == raw_peaks


def test_expected_csv_and_html_are_visibly_separate_sections():
    payload = _payload()
    buffer = io.StringIO()
    writer = csv.writer(buffer)
    writer.writerow(["Ordinary_MS2_Matches"])
    append_expected_csv_section(writer, payload)
    csv_text = buffer.getvalue()

    assert "独立连续 MS1 靶向层" in csv_text
    assert "不计入普通 MS2 谱库匹配数" in csv_text
    assert "Glucose" in csv_text
    assert "False" not in csv_text
    assert "Isomer_" not in csv_text
    assert "未提供空白" not in csv_text

    html = inject_expected_html(
        "<html><body><main><h1>Ordinary report</h1></main></body></html>",
        payload,
    )
    assert 'id="expected-compounds-report"' in html
    assert "Glucose" in html
    assert html.index("Ordinary report") < html.index(
        "expected-compounds-report"
    )


def test_export_page_discloses_expected_compounds_report_contract():
    from pathlib import Path

    root = Path(__file__).resolve().parents[1]
    page = (root / "web" / "static" / "index.html").read_text(
        encoding="utf-8"
    )
    export_script = (
        root / "web" / "static" / "js" / "export.js"
    ).read_text(encoding="utf-8")
    assert "结果导航 → 靶向高证据 → MS2 谱库鉴定" in page
    assert "靶向 CSV 与 MS2 CSV 分开下载" in page
    assert "不计入普通 MS2 的匹配数、High 数或去重候选数" in page
    assert "targeted_features" in page
    assert "/static/js/export.js?v=20261022-core-evidence-repair-v1" in page
    assert 'id="rawPeaksExportCard"' in page
    assert 'onclick="ExportMod.downloadRawPeaks()"' in page
    assert "完整原始峰 JSONL.gz" in page
    assert "queueSingleExport('raw-peaks')" in export_script
    assert 'id="plantCandidateExportCard"' in page
    assert "'plant_candidates'" in page
    assert "'external_review'" in page
    assert "download_nonce: String(Date.now())" in export_script
    assert "result_layer: resultLayer" in export_script
    assert "ms2_schema: 'reader_v3'" in export_script
    assert "_downloadTargetedSummaryCompat" in export_script
    assert "_downloadMs2Compat" in export_script
    assert "MS2_READER_HEADERS" in export_script
    assert 'id="batchExportType"' in page
    assert "CSV — MS2 鉴定结果（含异构体候选；每样本 1 份）" in page
    assert "CSV — 靶向离子 / 原始峰审计（逐峰）" in page
    assert "CSV 完整包（每样本 3 份）" in page
    assert "仅 CSV — 异构体竞争候选明细（逐候选）" in page
    assert "ExportMod.batchExport('csv')" in page
    assert 'id="batchCsvBundleButton"' in page
    assert 'id="batchExportButton"' in page
    assert "HTML — 仅靶向库匹配" in page
    assert "batchExportSelected" in export_script


def test_targeted_only_html_uses_normal_names_and_no_internal_group_heading():
    payload = _payload()
    payload["features"][0]["candidate_name"] = "Glucose"
    payload["features"][0]["candidate_names"] = ["Glucose", "Dextrose"]
    rendered = render_expected_html_report(
        payload,
        "sample.mzML",
        analysis_warnings=[{
            "code": "deconvolution_isolation_evidence_missing",
            "modules": ["DecoID 类嵌合谱去卷积"],
            "message": "缺少隔离窗证据，去卷积未运行。",
            "action": "保留 precursor isolation window 元数据。",
        }],
    )

    assert "靶向库匹配报告" in rendered
    assert "Glucose" in rendered
    assert "异构体组" not in rendered
    assert "异构组" not in rendered
    assert "Dextrose" not in rendered
    assert "歧义组" not in rendered
    assert "同质量歧义" not in rendered
    assert "同分异构体" not in rendered
    assert 'class="expected-raw-feature-audit"' in rendered
    assert '<details id="ms1-module-warnings"' in rendered
    assert '<details open id="ms1-module-warnings"' not in rendered
    assert "缺少隔离窗证据，去卷积未运行。" in rendered


def test_targeted_html_collapses_very_long_compound_names_without_losing_them():
    payload = _payload()
    long_name = (
        "(2S,3R,4S,5S,6R)-2-[[(3S,5R,6S)-6-hydroxy-3,5-dimethyl-"
        "cyclohexyl]oxy]-6-(hydroxymethyl)oxane-3,4,5-triol"
    )
    payload["features"][0]["candidate_name"] = long_name
    payload["features"][0]["candidate_names"] = [long_name]

    rendered = render_expected_html_report(payload, "sample.mzML")

    assert 'class="expected-full-compound-name"' in rendered
    assert "查看完整名称" in rendered
    assert html.escape(long_name) in rendered
    assert html.escape(long_name[:53].rstrip()) + "…" in rendered


def test_server_exports_expected_compounds_to_all_main_report_formats(
    monkeypatch, tmp_path,
):
    import web.server as server

    task_id = "expected-export"
    monkeypatch.setattr(server, "TASK_DIR", tmp_path)
    server._tasks[task_id] = {
        "id": task_id,
        "filename": "sample.mzML",
        "status": "completed",
        "total_spectra": 0,
        "matched_spectra": 0,
        "match_rate": 0.0,
        "results": [],
        "expected_compounds_summary": _payload()["summary"],
        "analysis_warnings": [{
            "code": "deconvolution_isolation_evidence_missing",
            "severity": "warning",
            "modules": ["DecoID 类嵌合谱去卷积"],
            "message": "缺少隔离窗证据，去卷积未运行。",
            "action": "保留 precursor isolation window 元数据。",
        }],
    }
    write_expected_compounds_artifact(
        tmp_path / f"{task_id}_expected_compounds.json",
        _payload(),
    )

    client = TestClient(server.app)
    try:
        targeted_csv = client.get(
            f"/api/tasks/{task_id}/download/csv?result_layer=targeted"
        )
        assert targeted_csv.status_code == 200
        assert "证据等级,主要候选,候选数,靶向数据库" in targeted_csv.text
        assert "Glucose" in targeted_csv.text
        assert "False" not in targeted_csv.text
        assert "Isomer_" not in targeted_csv.text
        assert "未提供空白" not in targeted_csv.text
        assert targeted_csv.headers["x-mssystem-result-layer"] == "targeted"
        assert "no-store" in targeted_csv.headers["cache-control"]
        assert "targeted_compounds.csv" in targeted_csv.headers[
            "content-disposition"
        ]

        targeted_audit_csv = client.get(
            f"/api/tasks/{task_id}/download/csv"
            "?result_layer=targeted_features"
        )
        assert targeted_audit_csv.status_code == 200
        assert "靶向峰 ID,证据等级,候选名称" in (
            targeted_audit_csv.text
        )
        assert "ECI_1_P001" in targeted_audit_csv.text
        assert targeted_audit_csv.headers[
            "x-mssystem-result-layer"
        ] == "targeted_features"
        assert "targeted_peak_audit.csv" in targeted_audit_csv.headers[
            "content-disposition"
        ]

        ms2_csv = client.get(f"/api/tasks/{task_id}/download/csv")
        assert ms2_csv.status_code == 200
        assert ms2_csv.text.startswith(
            "\ufeff鉴定组,离子特征编号,候选化合物,组内排名,候选角色,组内候选数,分子式,鉴定结论"
        )
        assert "样本 CE,样本 CE 单位,样本仪器类型" in (
            ms2_csv.text
        )
        assert "代表谱匹配碎片 (m/z:强度)" in ms2_csv.text
        assert "原始 MS2 碎片离子（已保存输入峰，m/z:强度）" in ms2_csv.text
        assert (
            "候选化学类别,鉴定组化学类别,分类置信度,分类依据,"
            "分类证据状态,分类可靠性说明"
        ) in ms2_csv.text
        assert "鉴定组 ID" not in ms2_csv.text
        assert "候选记录状态" not in ms2_csv.text
        assert "Glucose" not in ms2_csv.text
        assert ms2_csv.headers["x-mssystem-result-layer"] == "ms2"
        assert ms2_csv.headers["x-mssystem-ms2-schema"] == "reader_v3"
        hashes = {
            hashlib.sha256(response.content).hexdigest()
            for response in (targeted_csv, targeted_audit_csv, ms2_csv)
        }
        assert len(hashes) == 3

        excel_response = client.get(f"/api/tasks/{task_id}/download/excel")
        assert excel_response.status_code == 200
        workbook = load_workbook(
            io.BytesIO(excel_response.content), read_only=False
        )
        assert workbook.sheetnames[:7] == [
            "00_置信度说明",
            "01_结果导航",
            "02_靶向_高证据",
            "03_MS2谱库鉴定",
            "MS2_异构体候选明细",
            "04_靶向_同位素",
            "05_靶向_仅质量",
        ]
        assert "10_靶向峰审计" in workbook.sheetnames
        assert "11_模块提示" in workbook.sheetnames
        confidence_values = "\n".join(
            str(cell.value or "")
            for row in workbook["00_置信度说明"].iter_rows()
            for cell in row
        )
        assert "High 不等于已经确证" in confidence_values
        assert "Final score／综合匹配分" in confidence_values
        feature_sheet = workbook["04_靶向_同位素"]
        headers = [cell.value for cell in feature_sheet[1]]
        assert "主要候选" in headers
        assert feature_sheet["B2"].value == "Glucose"
        assert "Isomer_Ambiguous" not in headers
        assert feature_sheet.row_dimensions[2].height == 21
        assert feature_sheet["B2"].alignment.wrap_text is not True
        assert all(
            not isinstance(cell.value, bool)
            for row in feature_sheet.iter_rows()
            for cell in row
        )
        summary_sheet = workbook["01_结果导航"]
        summary_values = {
            summary_sheet.cell(row, 1).value: summary_sheet.cell(row, 2).value
            for row in range(2, summary_sheet.max_row + 1)
        }
        assert summary_values["靶向检出特征数"] == 1
        assert summary_values["MS1 模块适用性提示"] == 1
        assert "High 级鉴定组数" in summary_values
        assert "High 单一二维连接数" in summary_values
        assert "High 未解异构体组数" in summary_values
        assert "高证据支持化合物数" not in summary_values
        assert "Upload_Time" not in summary_values
        assert "Analysis_Time" not in summary_values
        peak_audit_sheet = workbook["10_靶向峰审计"]
        assert peak_audit_sheet.max_row == 2
        assert peak_audit_sheet["A2"].value == "ECI_1_P001"
        warning_sheet = workbook["11_模块提示"]
        assert warning_sheet["D2"].value == "缺少隔离窗证据，去卷积未运行。"

        html_response = client.get(
            f"/api/tasks/{task_id}/download/html?format_type=standard"
        )
        assert html_response.status_code == 200
        assert 'id="expected-compounds-report"' in html_response.text
        assert "Glucose" in html_response.text
        assert 'class="expected-raw-feature-audit"' in html_response.text
        assert '<details id="ms1-module-warnings"' in html_response.text
        assert 'id="mssystem-confidence-guide"' in html_response.text
        assert "High 不等于已经确证" in html_response.text

        targeted_html = client.get(
            f"/api/tasks/{task_id}/download/html?format_type=targeted"
        )
        assert targeted_html.status_code == 200
        assert "靶向库匹配报告" in targeted_html.text
        assert "Glucose" in targeted_html.text
        assert 'class="expected-raw-feature-audit"' in targeted_html.text
        assert '<details id="ms1-module-warnings"' in targeted_html.text
        assert 'id="mssystem-confidence-guide"' in targeted_html.text
        assert targeted_html.headers["x-mssystem-report-layer"] == "targeted"
        assert "targeted_library_report.html" in targeted_html.headers[
            "content-disposition"
        ]

        pdf_response = client.get(f"/api/tasks/{task_id}/download/pdf")
        assert pdf_response.status_code == 200
        assert pdf_response.content.startswith(b"%PDF")

        batch_response = client.post(
            "/api/export/batch/excel",
            json=[task_id],
        )
        assert batch_response.status_code == 200
        with zipfile.ZipFile(io.BytesIO(batch_response.content)) as archive:
            workbook_names = [
                name for name in archive.namelist()
                if name.endswith(".xlsx")
            ]
            assert len(workbook_names) == 1
            batch_workbook = load_workbook(
                io.BytesIO(archive.read(workbook_names[0])),
                read_only=True,
            )
            assert "04_靶向_同位素" in batch_workbook.sheetnames
            assert batch_workbook["04_靶向_同位素"]["B2"].value == (
                "Glucose"
            )

        combined_csv_response = client.post(
            "/api/export/batch/csv",
            json=[task_id],
        )
        assert combined_csv_response.status_code == 200
        assert combined_csv_response.headers[
            "x-mssystem-batch-format"
        ] == "csv"
        with zipfile.ZipFile(
            io.BytesIO(combined_csv_response.content)
        ) as archive:
            names = archive.namelist()
            assert "00_confidence_and_score_guide.csv" in names
            layer_names = {
                "ms2": next(
                    name for name in names if name.endswith("_ms2.csv")
                ),
                "targeted": next(
                    name
                    for name in names
                    if name.endswith("_targeted_ms1.csv")
                ),
                "targeted_audit": next(
                    name
                    for name in names
                    if name.endswith("_targeted_peak_audit.csv")
                ),
            }
            layer_bytes = {
                key: archive.read(name)
                for key, name in layer_names.items()
            }
            assert len({
                hashlib.sha256(content).hexdigest()
                for content in layer_bytes.values()
            }) == 3
            assert not any(
                name.endswith("_isomer_competing_candidates.csv")
                for name in names
            )
            assert layer_bytes["ms2"].decode("utf-8-sig").startswith(
                "鉴定组,离子特征编号,候选化合物,组内排名,候选角色"
            )
            assert layer_bytes["targeted"].decode(
                "utf-8-sig"
            ).startswith(
                "证据等级,主要候选,候选数,靶向数据库"
            )
            assert layer_bytes["targeted_audit"].decode(
                "utf-8-sig"
            ).startswith(
                "靶向峰 ID,证据等级,候选名称,候选数"
            )

        for batch_format, expected_suffix in (
            ("csv_ms2", "_ms2_identifications.csv"),
            ("csv_isomer_candidates", "_isomer_competing_candidates.csv"),
            ("csv_targeted", "_targeted_compounds.csv"),
            ("csv_targeted_features", "_targeted_peak_audit.csv"),
            ("csv_audit", "_ms2_detailed_audit.csv"),
            ("html_standard", "_standard_report.html"),
            ("html_detailed", "_detailed_report.html"),
            ("html_targeted", "_targeted_library_report.html"),
        ):
            response = client.post(
                f"/api/export/batch/{batch_format}", json=[task_id]
            )
            assert response.status_code == 200, (
                batch_format, response.text
            )
            assert response.headers["x-mssystem-batch-format"] == (
                batch_format
            )
            assert response.headers["x-mssystem-exported-tasks"] == "1"
            assert response.headers["x-mssystem-failed-tasks"] == "0"
            assert "no-store" in response.headers["cache-control"]
            with zipfile.ZipFile(io.BytesIO(response.content)) as archive:
                assert "00_confidence_and_score_guide.csv" in archive.namelist()
                assert any(
                    name.endswith(expected_suffix)
                    for name in archive.namelist()
                ), (batch_format, archive.namelist())
                if batch_format.startswith("html_"):
                    report_name = next(
                        name for name in archive.namelist()
                        if name.endswith(expected_suffix)
                    )
                    assert b'mssystem-confidence-guide' in archive.read(
                        report_name
                    )
    finally:
        client.close()
        server._tasks.pop(task_id, None)


def test_target_free_main_exports_never_load_target_catalogs(monkeypatch, tmp_path):
    """Generic exports retain candidate evidence without selecting target lists."""
    import copy
    import matching.expected_compounds as expected
    import web.server as server

    calls = []

    def forbidden(*args, **kwargs):
        calls.append((args, kwargs))
        raise AssertionError("A generic export accessed target data")

    for module, names in (
        (expected, ("load_expected_compound_catalog", "list_expected_datasets",
                    "resolve_expected_dataset_selection", "run_expected_compounds_workflow")),
        (server, ("list_expected_datasets", "resolve_expected_dataset_selection",
                  "run_expected_compounds_workflow")),
    ):
        for name in names:
            monkeypatch.setattr(module, name, forbidden)

    task_id = "generic-no-targets"
    results = [{
        "spectrum_id": "scan-generic", "precursor_mz": 47.0491,
        "retention_time": 2.0, "ion_mode": "POS", "peaks_count": 1,
        "query_peaks": [[31.0, 100.0]], "raw_peaks": [[31.0, 100.0]],
        "lib_peaks": [[31.0, 100.0]],
        "match": {
            "compound_name": "Ethanol", "formula": "C2H6O",
            "smiles": "CCO", "inchikey": "LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
            "dataset_source": "synthetic-test", "lib_spectrum_id": "reference-1",
            "cosine": 0.65, "final_score": 0.65, "rank_score": 0.65,
            "confidence": "Medium", "mass_ppm": 1.0,
            "lib_peaks": [[31.0, 100.0]], "top_candidates": [],
            "final_claim_decision": {"policy_id": "test-only", "released": False,
                                     "scope": "candidate_only"},
        },
    }]
    task = {"id": task_id, "filename": "RUN-A.mzML", "status": "completed",
            "total_spectra": 1, "matched_spectra": 1, "match_rate": 100.0,
            "params": {"enable_expected_compounds": False}, "results": results}
    monkeypatch.setattr(server, "TASK_DIR", tmp_path)
    server._tasks[task_id] = task
    # Compare at the real presentation boundary. Legacy normalization is an
    # existing input step, not an export-induced score/governance change.
    results = server._get_task_results(task_id)
    original = copy.deepcopy(results)
    client = TestClient(server.app)
    try:
        for suffix in ("csv", "excel", "html?format_type=standard",
                       "html?format_type=detailed", "pdf"):
            response = client.get(f"/api/tasks/{task_id}/download/{suffix}")
            assert response.status_code == 200, (suffix, response.text[:300])
            assert response.content
        for format_name in ("csv", "excel", "html_standard", "html_detailed"):
            response = client.post(f"/api/export/batch/{format_name}", json=[task_id])
            assert response.status_code == 200, format_name
            assert response.headers["x-mssystem-failed-tasks"] == "0"
        assert calls == []
        # HTML assigns only this non-scientific row navigation index.
        assert [{k: v for k, v in row.items() if k != "_result_index"}
                for row in results] == original
        assert not (tmp_path / f"{task_id}_expected_compounds.json").exists()
    finally:
        client.close()
        server._tasks.pop(task_id, None)


def test_batch_detailed_html_uses_deduplicated_ms2_representatives(
    monkeypatch,
):
    import web.server as server

    task_id = "batch-detailed-dedup"
    task = {
        "id": task_id,
        "filename": "project_sample.mzML",
        "status": "completed",
        "matched_spectra": 2,
        "match_rate": 100.0,
        "params": {},
    }

    def result(spectrum_id: str, final_score: float) -> dict:
        match = {
            "compound_name": "Ethanol",
            "inchikey": "LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
            "smiles": "CCO",
            "formula": "C2H6O",
            "classification": "other",
            "dataset_source": "test",
            "lib_spectrum_id": f"lib-{spectrum_id}",
            "cosine": final_score,
            "final_score": final_score,
            "rank_score": final_score,
            "mass_ppm": 1.0,
            "lib_peaks": [[31.0, 100.0]],
            "top_candidates": [],
        }
        return {
            "spectrum_id": spectrum_id,
            "precursor_mz": 47.0491,
            "peaks_count": 1,
            "query_peaks": [[31.0, 100.0]],
            "raw_peaks": [[31.0, 100.0]],
            "lib_peaks": [[31.0, 100.0]],
            "match": match,
        }

    results = [result("scan-1", 0.95), result("scan-2", 0.80)]
    server._tasks[task_id] = task
    monkeypatch.setattr(
        server,
        "_load_task_results_data",
        lambda requested: (task, results) if requested == task_id else (None, []),
    )
    monkeypatch.setattr(
        server,
        "_load_expected_compounds_payload",
        lambda *_args, **_kwargs: {},
    )

    client = TestClient(server.app)
    try:
        response = client.post(
            "/api/export/batch/html_detailed",
            json=[task_id],
        )
        assert response.status_code == 200, response.text
        with zipfile.ZipFile(io.BytesIO(response.content)) as archive:
            report_name = next(
                name
                for name in archive.namelist()
                if name.endswith("_detailed_report.html")
            )
            report = archive.read(report_name).decode("utf-8")
        assert report.count('class="main-row"') == 1
        assert 'data-name="ethanol"' in report
        assert "0.950" in report
        assert "0.800" not in report
    finally:
        client.close()
        server._tasks.pop(task_id, None)
