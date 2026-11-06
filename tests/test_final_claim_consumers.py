from __future__ import annotations

import io
from pathlib import Path

from fastapi.testclient import TestClient
import pdfplumber

from web import server
from web.report import generate_html
from web.report_common import RESULTS_COLUMN_NAMES, row_to_dict
from web.report_en import generate_english_html


def _projected_match() -> dict:
    return {
        "compound_name": "Candidate",
        "formula": "C10H12O2",
        "smiles": "CC1=CC=CC=C1O",
        "inchikey": "AAAAAAAAAAAAAA-UHFFFAOYSA-N",
        "cosine": 0.91,
        "final_score": 0.81,
        "confidence": "High",
        "compound_confidence": "High",
        "confidence_class": "High",
        "honesty_tier": "putative",
        "classification": "phenol",
        "final_claim_decision": {
            "policy_id": "evigate-ms-v2.0.0-20260803",
            "status": "blocked_project_calibration_unavailable",
            "scope": "2d_connectivity",
            "released": False,
            "independent_injection_count": 1,
            "reason_codes": ["project_level_qvalue_unavailable"],
            "calibration": {
                "status": "unavailable",
                "reason": "project_level_target_decoy_qvalue_not_yet_available",
            },
        },
    }


def _row() -> dict:
    return {
        "spectrum_id": "scan=1",
        "precursor_mz": 165.091,
        "ion_mode": "POSITIVE",
        "peaks_count": 3,
        "query_peaks": [[77.0, 100.0], [107.0, 50.0]],
        "match": _projected_match(),
    }


def test_standard_export_separates_legacy_confidence_from_final_claim():
    exported = row_to_dict(_row(), report_mode=False)

    assert "Final_Claim_Status" in RESULTS_COLUMN_NAMES
    assert exported["Confidence"] == "High"
    assert exported["Final_Claim_Released"] is False
    assert exported["Calibration_Status"] == "unavailable"
    assert exported["Project_Q_Value"] is None


def test_chinese_and_english_html_show_final_claim(tmp_path):
    match = _projected_match()
    report_row = {
        "precursor_mz": 165.091,
        "retention_time": 1.2,
        "peaks_count": 2,
        "query_peaks": [[77.0, 100.0]],
        "matches": [match],
    }
    output = tmp_path / "report.html"
    generate_html([report_row], output, "sample.mzML")
    chinese = output.read_text(encoding="utf-8")
    english = generate_english_html([_row()], "sample.mzML")

    assert "鉴定结论" in chinese
    assert "未达到当前结论条件" in chinese
    assert "Identification conclusion" in english
    assert "Criteria not met" in english


def test_standard_html_uses_fixed_name_width_and_reportable_scope_confidence(tmp_path):
    match = _projected_match()
    match.update({
        "confidence": "Medium",
        "compound_confidence": "Medium",
        "exact_structure_confidence": "Medium",
        "two_dimensional_connectivity_confidence": "High",
        "identification_scope": "exact_2d_connectivity",
        "honesty_tier": "connectivity_high_stereo_unresolved",
    })
    output = tmp_path / "fixed-width-report.html"
    generate_html([{
        "precursor_mz": 165.091,
        "retention_time": 1.2,
        "peaks_count": 2,
        "query_peaks": [[77.0, 100.0]],
        "matches": [match],
    }], output, "sample.mzML")
    rendered = output.read_text(encoding="utf-8")

    assert "#matchTable{min-width:2240px;table-layout:fixed}" in rendered
    assert '<th style="width:190px">Compound Name</th>' in rendered
    assert "High（二维连接）" in rendered
    assert 'data-conf="High"' in rendered
    assert "完整身份：Medium" in rendered


def test_pdf_separates_matched_ms2_from_dedup_and_keeps_scope_high():
    match = _projected_match()
    match.update({
        "confidence": "Medium",
        "compound_confidence": "Medium",
        "exact_structure_confidence": "Medium",
        "two_dimensional_connectivity_confidence": "High",
        "identification_scope": "exact_2d_connectivity",
        "honesty_tier": "connectivity_high_stereo_unresolved",
    })
    row_one = _row()
    row_one["match"] = match
    row_two = {**row_one, "spectrum_id": "scan=2", "precursor_mz": 165.092}
    payload = server._build_pdf_report(
        {
            "id": "pdf-scope-test",
            "filename": "scope-test.mzML",
            "compound_conclusion_summary": {
                "high_confidence_scope_group_count": 1,
                "high_confidence_2d_connectivity_count": 1,
                "high_confidence_isomer_group_count": 0,
            },
        },
        [row_one, row_two, {"spectrum_id": "scan=3", "match": None}],
        reference_standard_qc_payload={},
    )
    with pdfplumber.open(io.BytesIO(payload)) as document:
        text = "\n".join(page.extract_text() or "" for page in document.pages)

    assert "2/3" in text
    assert "去重候选: 1" in text
    assert "High（二维连接）" in text
    assert "完整身份 Medium" in text


def test_task_resources_expose_contract_roles_without_local_paths(monkeypatch, tmp_path):
    task_id = "resource-task"
    task_dir = tmp_path / "tasks"
    task_dir.mkdir()
    (task_dir / f"{task_id}_results.json").write_text("[]\n", encoding="utf-8")
    monkeypatch.setattr(server, "TASK_DIR", task_dir)
    monkeypatch.setattr(
        server,
        "_tasks",
        {task_id: {"id": task_id, "filename": "sample.mzML", "status": "completed"}},
    )

    response = TestClient(server.app).get(f"/api/tasks/{task_id}/resources")

    assert response.status_code == 200
    payload = response.json()
    resources = {item["resource"]: item for item in payload["resources"]}
    assert resources["legacy_results"]["status"] == "available"
    assert resources["final_claims"]["contract_role"] == (
        "evigate2_final_claim_contract"
    )
    assert all("path" not in item for item in payload["resources"])


def test_static_results_view_has_distinct_final_claim_column():
    source = (
        Path(__file__).resolve().parents[1] / "web" / "static" / "js" / "results.js"
    ).read_text(encoding="utf-8")
    assert "可作结论 / 已评估" in source
    assert "_finalClaimBadge" in source


def test_summary_projection_keeps_role_provenance_and_claim_reasons():
    row = _row()
    row["match"].update(
        {
            "compound_role": "external_contaminant",
            "potential_contaminant": True,
            "contaminant_category": "phthalate_plasticizer",
            "identity_source": "v9_curated_contaminant_baseline",
            "role_catalog_id": "mssystem-compound-role-v10-2026-08-04",
            "role_catalog_version": "10.1.0",
            "role_annotation_scope": "exact_identity_curation",
            "screening_domain": "pesticide_residue",
            "screening_catalog_id": "screening-test",
        }
    )

    projected = server._result_summary_projection(row, 0)["match"]

    assert projected["compound_role"] == "external_contaminant"
    assert projected["contaminant_category"] == "phthalate_plasticizer"
    assert projected["role_catalog_id"] == (
        "mssystem-compound-role-v10-2026-08-04"
    )
    assert projected["role_annotation_scope"] == "exact_identity_curation"
    assert projected["screening_catalog_id"] == "screening-test"
    assert projected["final_claim_decision"]["reason_codes"] == [
        "project_level_qvalue_unavailable"
    ]
