from __future__ import annotations

from copy import deepcopy
from io import BytesIO
from pathlib import Path
import subprocess

from fastapi.testclient import TestClient
from openpyxl import load_workbook

from web.governance_blocked_projection import (
    append_governance_blocked_csv_section,
    inject_governance_blocked_html,
    project_governance_blocked_candidates,
    render_governance_blocked_html,
)


def _row(index: int, status: str | None) -> dict:
    passport = None if status is None else {
        "policy_id": "evigate-ms-locked-v1-2026-08-02",
        "claim_status": status,
        "claim_scope": "spectral_candidate_only",
        "release_eligible_before_replication": False,
        "release_blockers": ["reason_one", "reason,two"],
    }
    match = {
        "compound_name": "Sinomenine" if index == 1 else "Ordinary candidate",
        "formula": "C19H23NO4",
        "inchikey": f"IDENTITY{index}",
        "smiles": "COC1=CC=CC=C1",
        "cosine": 0.91,
        "confidence": "High",
        "final_score": 0.72,
        "final_claim_decision": {
            "released": False,
            "status": "blocked_project_calibration_unavailable",
            "reason_codes": ["project_calibration_unavailable"],
        },
    }
    if passport is not None:
        match["evidence_governance"] = passport
    return {
        "spectrum_id": f"scan-{index}",
        "lcms_feature_id": f"LF-{index}",
        "retention_time": 3.25 + index,
        "retention_time_unit": "minute",
        "precursor_mz": 330.1 + index,
        "match": match,
    }


def test_projection_extracts_only_two_raw_blocked_states_without_mutation():
    rows = [
        _row(1, "reported_but_governance_blocked"),
        _row(2, "insufficient_for_release"),
        _row(3, "provisional_needs_independent_injection"),
        _row(4, None),
    ]
    before = deepcopy(rows)

    projected = project_governance_blocked_candidates(rows)

    assert rows == before
    assert [row["source_result_index"] for row in projected] == [0, 1]
    assert [row["claim_status"] for row in projected] == [
        "reported_but_governance_blocked",
        "insufficient_for_release",
    ]
    assert projected[0]["release_blockers"] == ["reason_one", "reason,two"]
    assert projected[0]["compound_name_cn"] == ["青藤碱", "清风藤碱"]
    assert projected[0]["identity_key"] == "IDENTITY1"
    assert rows[0]["match"]["confidence"] == "High"
    assert rows[0]["match"]["final_score"] == 0.72
    assert rows[0]["match"]["final_claim_decision"]["released"] is False


def test_projection_preserves_explicit_source_index_after_pagination():
    row = _row(1, "reported_but_governance_blocked")
    row["_result_index"] = 37

    projected = project_governance_blocked_candidates([row])

    assert projected[0]["source_result_index"] == 37


def test_csv_section_is_append_only_and_empty_projection_is_byte_stable():
    original = "A,B\n1,2\n"
    assert append_governance_blocked_csv_section(original, []) == original

    projected = project_governance_blocked_candidates([
        _row(1, "reported_but_governance_blocked")
    ])
    output = append_governance_blocked_csv_section(original, projected)

    assert output.startswith("A,B\n1,2\n\n")
    assert "治理未放行候选：1 条保存的扫描候选记录" in output
    assert "Identity_Key" in output
    assert '"reason_one; reason,two"' in output


def test_html_uses_closed_independent_tbody_and_empty_search_restores_rows():
    projected = project_governance_blocked_candidates([
        _row(1, "reported_but_governance_blocked"),
        _row(2, "insufficient_for_release"),
    ])
    html = render_governance_blocked_html(projected)

    assert '<details class="governance-blocked-panel"' in html
    assert '<details class="governance-blocked-panel" open' not in html
    assert '<tbody id="governanceBlockedBody">' in html
    assert "#governanceBlockedBody .governance-blocked-row" in html
    assert "#compoundResultsBody" not in html
    assert "#tableBody" not in html
    assert "var show = !query ||" in html
    assert html.count('class="governance-blocked-row"') == 2
    assert 'data-mz="331.1"' in html

    script = html.split("<script>", 1)[1].split("</script>", 1)[0]
    node_program = f"""
const blocked = [
  {{dataset: {{search: 'sinomenine', mz: '331.1'}}, style: {{display: 'none'}}}},
  {{dataset: {{search: 'ordinary', mz: '332.1'}}, style: {{display: 'none'}}}}
];
const mainRow = {{style: {{display: 'main-table-sentinel'}}}};
const counter = {{textContent: ''}};
global.document = {{
  querySelectorAll: selector => {{
    if (selector !== '#governanceBlockedBody .governance-blocked-row')
      throw new Error('selector escaped blocked tbody: ' + selector);
    return blocked;
  }},
  getElementById: id => id === 'governanceBlockedSearchCount' ? counter : null,
}};
{script}
filterGovernanceBlockedRows('absent');
if (!blocked.every(row => row.style.display === 'none')) throw new Error('filter failed');
filterGovernanceBlockedRows('');
if (!blocked.every(row => row.style.display === '')) throw new Error('empty search did not restore');
if (mainRow.style.display !== 'main-table-sentinel') throw new Error('main tbody changed');
"""
    completed = subprocess.run(
        ["node", "-e", node_program], capture_output=True, text=True, check=False
    )
    assert completed.returncode == 0, completed.stderr


def test_html_injection_covers_standard_report_without_changing_empty_report():
    document = "<html><body><main>existing report</main></body></html>"
    assert inject_governance_blocked_html(document, []) == document

    projected = project_governance_blocked_candidates([
        _row(1, "reported_but_governance_blocked")
    ])
    output = inject_governance_blocked_html(document, projected, language="en")

    assert "existing report" in output
    assert "Governance not released: 1 saved scan candidate records" in output
    assert output.index("Governance not released") < output.index("</body>")


def test_results_api_pages_only_the_blocked_rows_on_each_page(monkeypatch, tmp_path):
    from web import server

    task_id = "governance-pagination"
    rows = [
        _row(1, "reported_but_governance_blocked"),
        {
            "spectrum_id": "scan-unannotated",
            "precursor_mz": 250.0,
            "retention_time": 5.0,
            "match": None,
        },
        _row(3, "insufficient_for_release"),
    ]
    monkeypatch.setitem(server._tasks, task_id, {
        "id": task_id,
        "filename": "fixture.mzML",
        "status": "failed",
        "results": rows,
    })
    monkeypatch.setattr(server, "TASK_DIR", tmp_path)
    client = TestClient(server.app)

    first = client.get(
        f"/api/tasks/{task_id}/results?projection=summary&limit=1&offset=0"
    )
    second = client.get(
        f"/api/tasks/{task_id}/results?projection=summary&limit=1&offset=1"
    )
    third = client.get(
        f"/api/tasks/{task_id}/results?projection=summary&limit=1&offset=2"
    )

    assert first.status_code == second.status_code == third.status_code == 200
    assert first.json()["governance_blocked_count"] == 2
    assert [row["source_result_index"] for row in
            first.json()["governance_blocked_candidates"]] == [0]
    assert second.json()["governance_blocked_candidates"] == []
    assert [row["source_result_index"] for row in
            third.json()["governance_blocked_candidates"]] == [2]
    assert first.json()["results"][0]["match"]["final_score"] == 0.72
    assert first.json()["results"][0]["match"]["final_claim_decision"][
        "released"
    ] is False


def test_primary_csv_excel_and_pdf_share_the_blocked_projection(
    monkeypatch, tmp_path
):
    from web import server

    rows = [
        _row(1, "reported_but_governance_blocked"),
        _row(2, "provisional_needs_independent_injection"),
    ]
    csv_text = server._ms2_csv_text(rows, results=rows, language="en")
    assert "Governance not released: 1 saved scan candidate records" in csv_text
    assert "reported_but_governance_blocked" in csv_text
    assert "provisional_needs_independent_injection" not in csv_text.split(
        "=== Governance not released", 1
    )[1]

    task = {
        "id": "governance-export",
        "filename": "fixture.mzML",
        "status": "failed",
        "total_spectra": 2,
        "matched_spectra": 2,
        "match_rate": 100.0,
        "params": {},
    }
    monkeypatch.setattr(server, "TASK_DIR", tmp_path)
    monkeypatch.setattr(
        server,
        "_flagship_artifact_path",
        lambda *_args, **_kwargs: Path(tmp_path) / "absent.json",
    )
    excel = server._build_analysis_excel_bytes(
        task["id"], task, rows, language="en"
    )
    workbook = load_workbook(BytesIO(excel), read_only=True)
    blocked_sheet = workbook["Governance_Blocked"]
    assert blocked_sheet.max_row == 2
    assert blocked_sheet["D2"].value == "IDENTITY1"
    assert blocked_sheet["L2"].value == "reported_but_governance_blocked"

    pdf = server._build_pdf_report(
        {key: value for key, value in task.items() if key != "id"},
        rows,
        language="en",
    )
    assert pdf.startswith(b"%PDF")
    story = server._governance_blocked_pdf_story(
        rows, language="en", regular_font="Helvetica", bold_font="Helvetica-Bold"
    )
    story_text = " ".join(
        str(getattr(item, "text", ""))
        + " "
        + " ".join(
            getattr(cell, "getPlainText", lambda: "")()
            for table_row in getattr(item, "_cellvalues", [])
            for cell in table_row
        )
        for item in story
    )
    assert "Governance not released" in story_text
    assert "Sinomenine" in story_text
    assert "reason_one" in story_text


def test_legacy_excel_keeps_existing_sheets_and_adds_isolated_blocked_sheet():
    from web.app import _export_excel

    rows = [
        _row(1, "reported_but_governance_blocked"),
        _row(2, "provisional_needs_independent_injection"),
        {"spectrum_id": "scan-u", "precursor_mz": 120.0, "match": None},
    ]
    data = _export_excel({"filename": "fixture.mzML"}, rows)
    workbook = load_workbook(BytesIO(data), read_only=True)

    assert {"Matched", "Unmatched", "Summary", "治理未放行候选"}.issubset(
        workbook.sheetnames
    )
    blocked_sheet = workbook["治理未放行候选"]
    assert blocked_sheet.max_row == 2
    assert blocked_sheet["D2"].value == "IDENTITY1"


def test_detailed_html_accepts_preprojected_full_scan_rows():
    from web.app import _gen_rich

    projected = project_governance_blocked_candidates([
        _row(1, "reported_but_governance_blocked")
    ])
    html = _gen_rich([], governance_blocked_candidates=projected)

    assert "治理未放行候选 1 条" in html
    assert 'id="governanceBlockedBody"' in html
    assert "reason_one; reason,two" in html


def test_filtered_csv_keeps_original_saved_scan_index(monkeypatch):
    import web.server as server
    from fastapi.testclient import TestClient
    task_id = "blocked-filter-index"
    rows = [_row(1, "insufficient_for_release"), _row(2, "insufficient_for_release")]
    rows[0]["blank_clean_view_excluded"] = True
    monkeypatch.setitem(server._tasks, task_id, {"id": task_id, "filename": "test.mzML", "status": "completed"})
    monkeypatch.setattr(server, "_get_task_results", lambda _: rows)
    monkeypatch.setattr(server, "_load_expected_compounds_payload", lambda *a, **k: {})
    monkeypatch.setattr(server, "_load_reference_standard_qc_payload", lambda *a, **k: {})
    monkeypatch.setattr(server, "_ms2_csv_text", lambda matched, **kw: str(project_governance_blocked_candidates(kw["results"])[0]["source_result_index"]))
    response = TestClient(server.app).get(f"/api/tasks/{task_id}/download/csv?blank_view=clean")
    assert response.status_code == 200
    assert response.text.lstrip("\ufeff") == "1"
    assert "_result_index" not in rows[1]
