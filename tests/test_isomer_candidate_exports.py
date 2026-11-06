import csv
import io
import json
import re
import zipfile

import pdfplumber
from fastapi.testclient import TestClient
from openpyxl import load_workbook

from tests.test_english_report import _sample_results
from web.report_common import (
    isomer_candidate_export_rows,
    isomer_candidate_export_rows_from_compounds,
    isomer_candidate_projection,
    ms2_unified_candidate_rows,
    render_isomer_candidate_competition_html,
    row_to_dict,
)


def test_canonical_isomer_projection_exports_lead_and_every_competitor():
    results = _sample_results()
    projection = isomer_candidate_projection(results[0]["match"])

    assert projection["status"] == "unresolved_competition"
    assert projection["candidate_count"] == 2
    assert projection["competitive_candidate_count"] == 1
    assert projection["missing_candidate_record_count"] == 0
    assert [
        candidate["compound_name"]
        for candidate in projection["candidate_members"]
    ] == ["Quercetin", "Morin"]
    assert projection["candidate_members"][0]["is_ranking_lead"] is True
    assert projection["candidate_members"][1]["spectrally_competitive"] is True

    main_row = row_to_dict(results[0], report_mode=False)
    assert main_row["Isomer_Best_Candidate"] == "Quercetin"
    assert "Morin" in main_row["Isomer_Competing_Candidates"]
    payload = json.loads(main_row["Isomer_Candidate_Details_JSON"])
    assert [candidate["compound_name"] for candidate in payload["candidates"]] == [
        "Quercetin",
        "Morin",
    ]

    long_rows = isomer_candidate_export_rows(results, language="zh")
    assert len(long_rows) == 2
    assert [row["Candidate_Rank"] for row in long_rows] == [1, 2]
    assert [row["Compound_Name"] for row in long_rows] == ["Quercetin", "Morin"]
    assert long_rows[0]["Is_Ranking_Lead"] is True
    assert long_rows[1]["Is_Ranking_Lead"] is False


def test_final_compound_candidates_override_broader_search_provenance_names():
    representative = _sample_results()[0]
    compounds = [{
        "identity_scope": "isomer_group",
        "formula": "C15H10O7",
        "candidate_member_count": 2,
        "candidate_names": ["Quercetin", "Morin", "Kaempferol"],
        "candidate_members": [
            {
                "compound_name": "Quercetin",
                "formula": "C15H10O7",
                "inchikey": "REF-QUERCETIN",
                "cosine": 0.98,
            },
            {
                "compound_name": "Morin",
                "formula": "C15H10O7",
                "inchikey": "REF-MORIN",
                "cosine": 0.96,
            },
        ],
        "conclusion": {
            "isomer_group_only": True,
            "isomer_group_confidence": "High",
        },
        "representative": representative,
    }]

    rows = isomer_candidate_export_rows_from_compounds(
        compounds, language="zh"
    )

    assert [row["Compound_Name"] for row in rows] == [
        "Quercetin",
        "Morin",
    ]
    assert {row["Candidate_Count"] for row in rows} == {2}
    assert [row["Candidate_Rank"] for row in rows] == [1, 2]
    assert "Kaempferol" not in {row["Compound_Name"] for row in rows}

    unified = ms2_unified_candidate_rows(
        [], compounds=compounds, language="zh"
    )
    assert [row["Candidate_Compound_Name"] for row in unified] == [
        "Quercetin",
        "Morin",
    ]
    assert {row["Identification_Group_ID"] for row in unified} == {
        "compound-1"
    }
    assert {row["Record_Type"] for row in unified} == {
        "isomer_candidate"
    }


def test_csv_excel_and_batch_zip_carry_complete_isomer_candidate_table(
    monkeypatch, tmp_path,
):
    import web.server as server

    task_id = "isomer-export-contract"
    results = _sample_results()
    task = {
        "id": task_id,
        "filename": "F001_NEG_1.mzML",
        "status": "completed",
        "results": results,
        "total_spectra": len(results),
        "matched_spectra": 1,
    }
    monkeypatch.setattr(server, "TASK_DIR", tmp_path)
    server._tasks[task_id] = task
    client = TestClient(server.app)
    try:
        main_response = client.get(f"/api/tasks/{task_id}/download/csv")
        assert main_response.status_code == 200
        main_rows = list(csv.DictReader(io.StringIO(main_response.text.lstrip("\ufeff"))))
        assert len(main_rows) == 2
        assert [row["候选化合物"] for row in main_rows] == [
            "Quercetin",
            "Morin",
        ]
        assert [row["组内排名"] for row in main_rows] == ["1", "2"]
        assert [row["候选角色"] for row in main_rows] == [
            "最佳候选（异构体未唯一解析）",
            "竞争异构体候选",
        ]
        assert {row["鉴定组"] for row in main_rows} == {"G001"}
        assert {row["组内候选数"] for row in main_rows} == {"2"}
        assert "样本 MS2 碎片离子（实际检索谱，m/z:强度）" in main_rows[0]
        assert "候选化学类别" in main_rows[0]
        assert "鉴定组 ID" not in main_rows[0]
        assert "候选记录状态" not in main_rows[0]
        assert "异构体候选机器可读明细 (JSON)" not in main_rows[0]

        candidate_response = client.get(
            f"/api/tasks/{task_id}/download/csv?result_layer=isomer_candidates"
        )
        assert candidate_response.status_code == 200
        assert candidate_response.headers["x-mssystem-result-layer"] == (
            "isomer_candidates"
        )
        assert "F001_NEG_1_isomer_competing_candidates.csv" in (
            candidate_response.headers["content-disposition"]
        )
        candidate_rows = list(csv.DictReader(io.StringIO(
            candidate_response.text.lstrip("\ufeff")
        )))
        assert [row["候选化合物名称"] for row in candidate_rows] == [
            "Quercetin",
            "Morin",
        ]
        assert [row["候选排名"] for row in candidate_rows] == ["1", "2"]

        excel_response = client.get(f"/api/tasks/{task_id}/download/excel")
        assert excel_response.status_code == 200
        workbook = load_workbook(io.BytesIO(excel_response.content), read_only=True)
        assert "MS2_异构体候选明细" in workbook.sheetnames
        worksheet = workbook["MS2_异构体候选明细"]
        values = list(worksheet.values)
        headers = list(values[0])
        name_index = headers.index("候选化合物名称")
        assert [row[name_index] for row in values[1:]] == ["Quercetin", "Morin"]

        batch_response = client.post(
            "/api/export/batch/csv?language=zh",
            json=[task_id],
        )
        assert batch_response.status_code == 200
        with zipfile.ZipFile(io.BytesIO(batch_response.content)) as archive:
            names = archive.namelist()
            assert not any(
                name for name in names
                if name.endswith("_isomer_competing_candidates.csv")
            )
            main_name = next(name for name in names if name.endswith("_ms2.csv"))
            unified_csv = archive.read(main_name).decode("utf-8-sig")
        assert "Quercetin" in unified_csv
        assert "Morin" in unified_csv
        assert "候选化合物" in unified_csv

        ms2_package_response = client.post(
            "/api/export/batch/csv_ms2?language=zh",
            json=[task_id],
        )
        assert ms2_package_response.status_code == 200
        with zipfile.ZipFile(
            io.BytesIO(ms2_package_response.content)
        ) as archive:
            names = archive.namelist()
            assert any(
                name.endswith("_ms2_identifications.csv")
                for name in names
            )
            assert not any(
                name.endswith("_isomer_competing_candidates.csv")
                for name in names
            )
            main_name = next(
                name for name in names
                if name.endswith("_ms2_identifications.csv")
            )
            unified_csv = archive.read(main_name).decode("utf-8-sig")
        assert "Quercetin" in unified_csv
        assert "Morin" in unified_csv
    finally:
        server._tasks.pop(task_id, None)


def test_offline_html_and_pdf_do_not_require_live_page_to_read_candidates(
    tmp_path,
):
    from web.app import _gen_rich
    from web.report import generate_html
    from web.report_en import generate_english_html
    from web.server import _build_pdf_report

    results = _sample_results()
    report_rows = [
        {**row, "matches": [row["match"]] if row.get("match") else []}
        for row in results
    ]

    standard_path = tmp_path / "standard.html"
    generate_html(report_rows, standard_path, "F001_NEG_1.mzML")
    standard = standard_path.read_text(encoding="utf-8")
    detailed = _gen_rich(report_rows, "F001_NEG_1.mzML")
    english = generate_english_html(results, "F001_NEG_1.mzML", detailed=True)
    for rendered in (standard, detailed):
        assert "未解析异构体竞争候选——离线报告完整清单" in rendered
        assert "Quercetin" in rendered
        assert "Morin" in rendered
        assert '<details class="icx-group">' in rendered
        assert '<details class="icx-candidate">' in rendered
        assert '<details class="icx-group" open' not in rendered
        assert '<details class="icx-candidate" open' not in rendered
        assert "展开组内候选" in rendered
        assert "查看分数与证据" in rendered
        assert ".icx-group:not([open])>.icx-group-body" in rendered
    assert "Unresolved isomer competition — complete offline candidate list" in english
    assert "Quercetin" in english
    assert "Morin" in english
    assert "Show candidates" in english
    assert "Show scores and evidence" in english

    pdf = _build_pdf_report(
        {"id": "pdf-isomer", "filename": "F001_NEG_1.mzML"},
        results,
        language="en",
    )
    with pdfplumber.open(io.BytesIO(pdf)) as document:
        text = "\n".join(page.extract_text() or "" for page in document.pages)
    assert "Unresolved isomer competition" in text
    assert "Quercetin" in text
    assert "Morin" in text


def test_offline_isomer_appendix_uses_two_closed_disclosure_levels():
    results = _sample_results()
    results[0]["match"]["compound_confidence"] = "High"
    results[0]["match"]["isomer_resolution_status"] = (
        "diagnostically_uninformative"
    )
    rendered = render_isomer_candidate_competition_html(
        results, language="zh"
    )

    assert rendered.count('<details class="icx-group">') == 1
    assert rendered.count('<details class="icx-candidate">') == 2
    assert rendered.count('class="icx-candidate-body"') == 2
    assert "MS/MS 相似度" in rendered
    assert "谱熵相似度" in rendered
    assert "反向匹配分" in rendered
    assert "综合匹配分" in rendered
    assert "候选排序分" in rendered
    assert "质量误差（ppm）" in rendered
    assert "解析状态：" in rendered
    assert "当前 MS/MS 碎片不能区分组内异构体" in rendered
    assert "当前排序第 1（异构体未唯一解析）" in rendered
    assert "候选明细已保存在任务结果中" in rendered
    assert "diagnostically_uninformative" not in rendered
    assert "ranking_lead" not in rendered


def test_export_page_offers_direct_and_batch_isomer_candidate_downloads():
    from pathlib import Path

    root = Path(__file__).resolve().parents[1]
    page = (root / "web" / "static" / "index.html").read_text(encoding="utf-8")
    script = (root / "web" / "static" / "js" / "export.js").read_text(
        encoding="utf-8"
    )
    overview = (
        root / "web" / "static" / "js" / "sample-overview.js"
    ).read_text(encoding="utf-8")
    assert "downloadCSV('results','isomer_candidates')" in page
    assert 'value="csv_isomer_candidates|zh"' in page
    assert "CSV_ISOMER_CANDIDATES" in script
    assert "_downloadIsomerCandidatesCompat" in script
    assert "downloadMs2Csv" in script
    assert "_buildCompatCsvBatch" in script
    assert "const outcome = await _buildCompatCsvBatch" not in script
    assert "/export/jobs/${encodeURIComponent(format)}" in script
    assert "/export/jobs/${job.id}/download" in script
    assert "statusResponse" in script
    assert "let _batchExportInFlight = false" in script
    assert "if (_batchExportInFlight)" in script
    assert "URL.revokeObjectURL(url), 60_000" in script
    assert "_fetchAllCompoundConclusions" in script
    assert "_supportsIsomerCandidateCsv === false" in script
    assert "MS2_READER_HEADERS" in script
    assert "样本 MS2 碎片离子（实际检索谱，m/z:强度）" in script
    assert "原始 MS2 碎片离子（已保存输入峰，m/z:强度）" in script
    assert "候选化学类别" in script
    assert "候选化合物" in script
    assert "异构体候选机器可读明细 (JSON)" not in script
    assert "异构体候选已经逐行合并进 MS2 表" in page
    assert "SampleReport.downloadIsomerCandidateCsv()" in page
    assert "/static/js/confidence-guide.js?v=14.0.4" in page
    for script_name in ("report-generator", "sample-overview"):
        assert re.search(r'<script\s+src="/static/js/' + script_name + r'\.js\?v=[^"\s]+"', page)
    assert "const suppliedMemberNames" in overview
    assert "const candidateNames = suppliedMemberNames.length" in overview
    assert "if (!candidateMemberCount) candidateMemberCount = candidateNames.length" in (
        overview
    )
