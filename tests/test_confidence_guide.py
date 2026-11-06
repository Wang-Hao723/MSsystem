from __future__ import annotations

import io
import json
from pathlib import Path

from fastapi.testclient import TestClient
from openpyxl import load_workbook

from web.confidence_guide import (
    CONFIDENCE_GUIDE_HTML_ID,
    CONFIDENCE_GUIDE_POLICY_ID,
    confidence_guide_html,
    confidence_guide_payload,
    confidence_guide_table_rows,
    format_confidence_guide_worksheet,
    inject_confidence_guide_html,
)


ROOT = Path(__file__).resolve().parents[1]


def test_confidence_policy_keeps_identification_and_confirmation_separate():
    payload = confidence_guide_payload("zh")
    levels = {item["key"]: item for item in payload["levels"]}

    assert payload["policy_id"] == CONFIDENCE_GUIDE_POLICY_ID
    assert "High 不等于已经确证" in payload["read_first"]
    assert payload["safeguards"] == {
        "confidence_is_probability": False,
        "confidence_is_project_qvalue": False,
        "high_means_confirmed": False,
        "unmatched_means_absent": False,
        "auxiliary_scores_change_identity_rank_by_definition": False,
    }
    assert "不能单凭 High" in levels["high"]["avoid"]
    assert "候选短名单" in levels["medium"]["use"]
    assert "不要计入已鉴定化合物数量" in levels["low"]["avoid"]
    assert "未匹配不等于该化合物不存在" in levels["unknown"]["avoid"]
    assert "不计作样品支持的主要化合物" in levels["background"]["avoid"]


def test_identity_scopes_and_all_main_score_families_are_explained():
    payload = confidence_guide_payload("zh")
    scopes = "\n".join(item["label"] + item["report"] for item in payload["identity_scopes"])
    score_keys = {item["key"] for item in payload["scores"]}

    assert "单一二维连接结构" in scopes
    assert "异构体组中的某一成员" in scopes
    assert {
        "cosine", "entropy", "reverse_match", "mass_error",
        "fragment_support", "final_score", "rank_score",
        "candidate_margin", "confidence_label", "formal_claim",
        "auxiliary_scores",
    } <= score_keys
    final = next(item for item in payload["scores"] if item["key"] == "final_score")
    assert "不是后验概率" in final["boundary"]
    auxiliary = next(item for item in payload["scores"] if item["key"] == "auxiliary_scores")
    assert "不能替代化合物身份置信度" in auxiliary["boundary"]


def test_html_guide_is_complete_bilingual_and_injection_is_idempotent():
    zh = confidence_guide_html("zh")
    en = confidence_guide_html("en")
    source = "<html><body><header><h1>Report</h1></header><p>Body</p></body></html>"
    injected = inject_confidence_guide_html(source, "zh")

    assert f'id="{CONFIDENCE_GUIDE_HTML_ID}"' in zh
    assert "各项评分和辅助字段怎样读" in zh
    assert "证据等级速查" in zh
    assert "mcg-tier-map" in zh
    assert "<details open" not in zh
    assert "High does not mean" in en
    assert injected.index(CONFIDENCE_GUIDE_HTML_ID) < injected.index("Body")
    assert inject_confidence_guide_html(injected, "zh") == injected


def test_excel_report_starts_with_plain_language_confidence_guide():
    from web.server import _build_analysis_excel_bytes

    workbook_bytes = _build_analysis_excel_bytes(
        "confidence-guide-test",
        {
            "id": "confidence-guide-test",
            "filename": "demo.mzML",
            "analysis_warnings": [],
        },
        [],
        language="zh",
    )
    workbook = load_workbook(io.BytesIO(workbook_bytes), read_only=False)

    assert workbook.sheetnames[0] == "00_置信度说明"
    guide = workbook["00_置信度说明"]
    values = "\n".join(
        str(cell.value or "")
        for row in guide.iter_rows()
        for cell in row
    )
    assert "High 不等于已经确证" in values
    assert "Medium（中等证据候选）" in values
    assert "Final score／综合匹配分" in values
    assert guide.merged_cells.ranges
    assert guide["A1"].value == "怎样理解置信度与各项评分"
    assert guide["A1"].fill.fgColor.rgb.endswith("123B50")
    assert guide["A4"].alignment.wrap_text is True
    assert guide.row_dimensions[4].height == 52
    assert guide.freeze_panes == "A7"
    assert guide.page_setup.orientation == "landscape"


def test_api_and_frontend_use_the_same_canonical_guide():
    from web import server

    response = TestClient(server.app).get(
        "/api/interpretation/confidence-guide?language=zh"
    )
    assert response.status_code == 200
    assert response.json()["policy_id"] == CONFIDENCE_GUIDE_POLICY_ID

    index = (ROOT / "web/static/index.html").read_text(encoding="utf-8")
    overview = (ROOT / "web/static/js/sample-overview.js").read_text(encoding="utf-8")
    report = (ROOT / "web/static/js/report-generator.js").read_text(encoding="utf-8")
    renderer = (ROOT / "web/static/js/confidence-guide.js").read_text(encoding="utf-8")
    static_payloads = json.loads(
        (ROOT / "web/static/policies/confidence-guide.json").read_text(
            encoding="utf-8"
        )
    )
    server_source = (ROOT / "web/server.py").read_text(encoding="utf-8")

    assert "confidence-guide.js?v=14.0.4" in index
    assert static_payloads["zh"] == confidence_guide_payload("zh")
    assert static_payloads["en"] == confidence_guide_payload("en")
    assert "/static/policies/confidence-guide.json?v=1.0.0" in renderer
    assert "ConfidenceGuide.load(language)" in overview
    assert "sample-overview-confidence-guide" in overview
    assert "sample-report-confidence-guide" in report
    assert "/interpretation/confidence-guide" in renderer
    assert "confidence_guide_pdf_story" in server_source
    assert "inject_confidence_guide_html" in server_source
    assert '"00_confidence_and_score_guide.csv"' in server_source


def test_web_guide_uses_compact_progressive_disclosure_layout():
    renderer = (
        ROOT / "web/static/js/confidence-guide.js"
    ).read_text(encoding="utf-8")

    assert "confidence-guide-shell" in renderer
    assert "confidence-guide-report-shell" in renderer
    assert "confidence-guide-report-tier-map" in renderer
    assert "confidence-guide-level-switcher" in renderer
    assert "confidence-guide-score-grid" in renderer
    assert "High</b> ≠" in renderer
    assert "if (reportMode)" in renderer
    assert "<details class=\"confidence-guide-shell\">" in renderer


def test_flat_guide_rows_cover_levels_identity_scopes_and_scores():
    rows = confidence_guide_table_rows("zh")
    categories = {row["说明类型"] for row in rows}

    assert {"必读", "证据等级", "身份范围", "评分／字段"} <= categories
    assert any(row["名称"] == "High（高可信证据）" for row in rows)
    assert any(row["名称"] == "Candidate margin／候选差距" for row in rows)
