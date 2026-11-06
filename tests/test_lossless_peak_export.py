import copy
import csv
import gzip
import hashlib
import io
import json
import zipfile

import pytest
from fastapi.testclient import TestClient
from openpyxl import Workbook, load_workbook

from web.lossless_peak_export import (
    CHUNK_CHARACTERS, add_raw_peak_sheets, canonical_peak_json,
    iter_raw_peak_jsonl_gzip, peak_sheet_names, prepare_reader_peak_cells,
)


def _rows():
    return [
        {"spectrum_id": "=duplicate", "ion_mode": "NEG", "raw_peaks":
         [[50.1234567890123 + i / 987.0, 1e-8 if i % 3 else 1.2345678912345e9]
          for i in range(3001)], "raw_peaks_complete": True,
         "raw_peaks_input_count": 3001, "raw_peaks_saved_count": 3001},
        {"spectrum_id": "=duplicate", "ion_mode": "POS", "raw_peaks": [],
         "raw_peaks_complete": True, "raw_peaks_input_count": 0},
        {"spectrum_id": "missing", "ion_mode": "NEG", "raw_peaks_complete": False,
         "raw_peaks_input_count": 4000},
        {"spectrum_id": "null", "raw_peaks": None},
        {"spectrum_id": "short", "raw_peaks": [[-0.0, 42], [1e-20, 1.2e20]],
         "raw_peaks_complete": False, "raw_peaks_input_count": 9},
    ]


def _reconstruct(workbook, source, language="zh"):
    _, chunk_name = peak_sheet_names(language)
    parts = {}
    for name in workbook.sheetnames:
        if not name.startswith(chunk_name):
            continue
        for row in list(workbook[name].values)[1:]:
            parts.setdefault(row[1], []).append(row)
            assert len(row[-1]) <= CHUNK_CHARACTERS
    assert set(parts) == set(range(len(source)))
    for index, expected in enumerate(source):
        selected = sorted(parts[index], key=lambda row: row[9])
        assert [r[9] for r in selected] == list(range(1, len(selected) + 1))
        assert {r[10] for r in selected} == {len(selected)}
        text = "".join(row[-1] for row in selected)
        assert json.loads(text) == expected.get("raw_peaks")
        assert text == canonical_peak_json(expected.get("raw_peaks"))
        assert {row[8] for row in selected} == {hashlib.sha256(text.encode("ascii")).hexdigest()}
        assert {row[7] for row in selected} == {"raw_peaks" in expected}
        assert {row[6] for row in selected} == {expected.get("raw_peaks_complete")}


@pytest.mark.parametrize("language", ["zh", "en"])
def test_saved_arrays_roundtrip_and_other_cells_remain_exact(language):
    rows = _rows()
    original = copy.deepcopy(rows)
    book = Workbook()
    book.active.title = "Existing"
    book.active.append(["Header", "Value"])
    book.active.append(["untouched", 1.23456789])
    before = list(book.active.values)
    stats = add_raw_peak_sheets(book, rows, "task", language=language)
    assert stats["saved_peaks"] == 3003
    assert stats["chunks"] > len(rows)
    buf = io.BytesIO()
    book.save(buf)
    loaded = load_workbook(io.BytesIO(buf.getvalue()))
    _reconstruct(loaded, rows, language)
    assert list(loaded["Existing"].values) == before
    assert loaded[peak_sheet_names(language)[0]]["C2"].data_type == "s"
    assert rows == original


def test_page_overflow_splits_sheets_without_losing_chunks(monkeypatch):
    import web.lossless_peak_export as export
    monkeypatch.setattr(export, "EXCEL_MAX_ROWS", 4)
    book = Workbook()
    stats = add_raw_peak_sheets(book, _rows(), "task")
    assert len(stats["index_sheets"]) == 2
    assert len(stats["chunk_sheets"]) >= 2
    _reconstruct(book, _rows())


def test_gzip_preserves_empty_missing_incomplete_and_duplicate_ids():
    rows = _rows()
    compressed = b"".join(iter_raw_peak_jsonl_gzip(iter(rows), "t", "ab" * 32))
    records = [json.loads(line) for line in gzip.decompress(compressed).splitlines()]
    assert len(records) == len(rows)
    for index, (record, expected) in enumerate(zip(records, rows)):
        assert record["result_index"] == index
        assert record["raw_peaks"] == expected.get("raw_peaks")
        assert record["raw_peaks_present"] == ("raw_peaks" in expected)
        assert record["raw_peaks_complete"] == expected.get("raw_peaks_complete")
        assert record["source_results_sha256"] == "ab" * 32
    assert b"".join(iter_raw_peak_jsonl_gzip([], "t")) == b"".join(iter_raw_peak_jsonl_gzip([], "t"))
    assert gzip.decompress(b"".join(iter_raw_peak_jsonl_gzip([], "t"))) == b""


def test_preview_retains_other_fields_and_names_complete_artifact():
    row = {"Spectrum": "=duplicate", "Raw_MS_MS_Spectrum": "123:456;" * 8000,
           "Candidate_Final_Score": 0.981, "Candidate_Confidence": "High"}
    result = prepare_reader_peak_cells([row], task_id="t", results=_rows())[0]
    assert len(result["Raw_MS_MS_Spectrum"]) < CHUNK_CHARACTERS
    assert "0,1" in result["Raw_MS_MS_Spectrum"]
    assert "/api/tasks/t/download/raw-peaks" in result["Raw_MS_MS_Spectrum"]
    assert "15_原始峰JSON" in result["Raw_MS_MS_Spectrum"]
    assert result["Candidate_Final_Score"] == row["Candidate_Final_Score"]
    assert result["Candidate_Confidence"] == "High"
    assert prepare_reader_peak_cells([{"Raw_MS_MS_Spectrum": "12:34"}]) == [{"Raw_MS_MS_Spectrum": "12:34"}]


def test_live_routes_full_workbook_csv_and_saved_archive(monkeypatch, tmp_path):
    import web.server as server
    from tests.test_english_report import _sample_results
    rows = _sample_results()
    rows[0].update(_rows()[0])
    rows[0]["spectrum_id"] = "long-spectrum"
    rows.extend([{**r, "match": None, "precursor_mz": 500.0} for r in _rows()[1:]])
    original = copy.deepcopy(rows)
    task_id = "lossless-export-fixture"
    task = {"id": task_id, "filename": "测试_NEG.mzML", "status": "completed",
            "results": rows, "total_spectra": len(rows), "matched_spectra": 1}
    monkeypatch.setattr(server, "TASK_DIR", tmp_path)
    # This fixture contains synthetic rows, not the production SQL/cache pair.
    import runtime.verified_sql_retrieval as sql
    monkeypatch.setattr(sql, "sql_metadata_repair_preflight", lambda: {"compatible": True})
    monkeypatch.setattr(server, "_compound_conclusions_for_export", lambda *a: None)
    monkeypatch.setattr(server, "_compound_scope_summary_for_export", lambda *a: {})
    server._tasks[task_id] = task
    client = TestClient(server.app)
    try:
        # Exercise production builder without spawning a second full server.
        workbook = load_workbook(io.BytesIO(server._build_analysis_excel_bytes(task_id, task, rows)))
        _reconstruct(workbook, rows)
        response = client.get(f"/api/tasks/{task_id}/download/csv")
        assert response.status_code == 200
        csv_rows = list(csv.DictReader(io.StringIO(response.text.lstrip("\ufeff"))))
        raw_header = "原始 MS2 碎片离子（已保存输入峰，m/z:强度）"
        assert "15_原始峰JSON" in csv_rows[0][raw_header]
        headers = list(next(workbook["03_MS2谱库鉴定"].values))
        raw_col = headers.index(raw_header)
        assert list(workbook["03_MS2谱库鉴定"].values)[1][raw_col] == csv_rows[0][raw_header]
        saved = tmp_path / f"{task_id}_results.json"
        saved.write_text(json.dumps(rows), encoding="utf-8")
        before = saved.read_bytes()
        response = client.get(f"/api/tasks/{task_id}/download/raw-peaks")
        assert response.status_code == 200
        records = [json.loads(line) for line in gzip.decompress(response.content).splitlines()]
        assert [r["raw_peaks"] for r in records] == [r.get("raw_peaks") for r in original]
        assert response.headers["X-MSsystem-Source-Results-SHA256"] == hashlib.sha256(before).hexdigest()
        assert "filename*=UTF-8''" in response.headers["content-disposition"]
        assert saved.read_bytes() == before
        batch = client.post("/api/export/batch/csv_ms2?language=zh", json=[task_id])
        assert batch.status_code == 200
        with zipfile.ZipFile(io.BytesIO(batch.content)) as archive:
            member = next(name for name in archive.namelist() if name.endswith("_raw_peaks.jsonl.gz"))
            batch_records = [json.loads(line) for line in gzip.decompress(archive.read(member)).splitlines()]
        assert [r["raw_peaks"] for r in batch_records] == [r.get("raw_peaks") for r in original]
        task["status"] = "running"
        assert client.get(f"/api/tasks/{task_id}/download/raw-peaks").status_code == 409
    finally:
        server._tasks.pop(task_id, None)


def test_mixed_external_group_keeps_all_members_and_original_rank():
    from web.plant_candidate_export import split_plant_reader_rows
    rows = [
        {"Group_Number": "G001", "Candidate_Compound_Name": "External",
         "Candidate_InChIKey": "FKMHSNTVILORFA-UHFFFAOYSA-N",
         "Candidate_Rank": 1, "Group_Candidate_Count": 2, "Candidate_Confidence": "High"},
        {"Group_Number": "G001", "Candidate_Compound_Name": "Unclassified competitor",
         "Candidate_Rank": 2, "Group_Candidate_Count": 2, "Candidate_Confidence": "Medium"},
        {"Group_Number": "G002", "Candidate_Compound_Name": "Unknown origin retained",
         "Candidate_Rank": 1, "Group_Candidate_Count": 1},
    ]
    plant, review = split_plant_reader_rows(rows)
    assert [r["Group_Number"] for r in plant] == ["G002"]
    assert [r["Candidate_Rank"] for r in review] == [1, 2]
    assert [r["Candidate_Confidence"] for r in review] == ["High", "Medium"]
    assert all(r["Group_Candidate_Count"] == 2 for r in review)


def test_persisted_group_filter_empty_selection_does_not_fall_back_to_all():
    from web.plant_candidate_export import plant_reader_views
    from tests.test_english_report import _sample_results
    results = _sample_results()
    compound = {
        "identity_scope": "isomer_group", "candidate_member_count": 2,
        "candidate_members": [
            {"compound_name": "Quercetin", "formula": "C15H10O7", "inchikey": "REF-QUERCETIN", "cosine": .98},
            {"compound_name": "Morin", "formula": "C15H10O7", "inchikey": "REF-MORIN", "cosine": .96},
        ],
        "conclusion": {"isomer_group_only": True, "isomer_group_confidence": "High", "counted_as_plant_candidate": False},
        "plant_report": {"excluded_from_plant_candidate_table": True},
        "representative": results[0],
    }
    plant, review = plant_reader_views(results, compounds=[compound])
    assert plant == []
    assert len(review) == 2
    assert [row["Candidate_Rank"] for row in review] == [1, 2]
    assert [row["Group_Candidate_Count"] for row in review] == [2, 2]
