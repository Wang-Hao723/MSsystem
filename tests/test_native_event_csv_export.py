"""Native MS/MS CSV exports preserve one reversible row per saved event."""

import csv
import io
import json

from fastapi.testclient import TestClient

from web import server


def _saved_rows():
    return [
        {
            "spectrum_id": "scan=matched",
            "precursor_mz": 301.1234,
            "search_precursor_mz": 300.1161,
            "search_precursor_mz_source": "lcms_feature_mz",
            "retention_time": 5.25,
            "retention_time_unit": "minute",
            "ion_mode": "POS",
            "precursor_charge": 1,
            "collision_energy": 20.0,
            "collision_energy_unit": "eV",
            "collision_energy_mode": "fixed",
            "activation_method": "HCD",
            "instrument_configuration_id": "IC-matched",
            "lcms_feature_id": "F001",
            "lcms_feature_mz": 300.1161,
            "lcms_feature_rt_minutes": 5.24,
            "lcms_feature_used_for_search": True,
            "blank_flagged": False,
            "raw_peaks": [[75.0, 10.0], [120.0, 50.0]],
            "match": {
                "compound_name": "Matched candidate",
                "confidence": "High",
                "evidence_tier": "MS/MS library candidate",
                "cosine": 0.8765,
                "mass_ppm": -1.25,
                "nested_audit": {"source": "saved-match"},
            },
        },
        {
            "spectrum_id": "scan=unmatched",
            "precursor_mz": 402.2,
            "search_precursor_mz": None,
            "retention_time": None,
            "ion_mode": "NEG",
            "precursor_charge": -1,
            "collision_energy": None,
            "lcms_feature_id": None,
            "match": None,
            "raw_peaks": [[88.0, 4.0]],
            "unmatched_audit": {"reason": "no_candidate"},
        },
        {
            "spectrum_id": "scan=blank",
            "precursor_mz": 503.3,
            "search_precursor_mz": 503.3,
            "retention_time": 7.5,
            "ion_mode": "POS",
            "precursor_charge": 1,
            "collision_energy": 35.0,
            "collision_energy_unit": "eV",
            "activation_method": "CID",
            "precursor_feature_id": "PF-blank",
            "blank_flagged": True,
            "blank_clean_view_excluded": True,
            "blank_sample_ratio": 0.4,
            "match": {"compound_name": "Blank candidate", "cosine": 0.55},
            "blank_evidence": {"sample": 40.0, "blank": 100.0},
        },
    ]


def test_complete_v2_exports_every_saved_event_and_reversible_json(monkeypatch):
    rows = _saved_rows()
    monkeypatch.setitem(
        server._tasks, "native-events",
        {"id": "native-events", "filename": "case1.mzML", "status": "completed"},
    )
    monkeypatch.setattr(server, "_get_task_results", lambda task_id: rows)

    response = TestClient(server.app).get(
        "/api/tasks/native-events/download/csv"
        "?format_type=results&result_layer=ms2&ms2_schema=complete_v2"
        "&blank_view=all"
    )
    assert response.status_code == 200
    assert "raw_msms_events_all.csv" in response.headers["content-disposition"]
    assert response.headers["x-mssystem-ms2-schema"] == "complete_v2"
    assert response.headers["x-mssystem-record-unit"] == "raw-msms-event"

    exported = list(csv.DictReader(io.StringIO(response.text.lstrip("\ufeff"))))
    assert len(exported) == len(rows) == 3
    assert [row["result_index"] for row in exported] == ["0", "1", "2"]
    assert [row["spectrum_id"] for row in exported] == [
        "scan=matched", "scan=unmatched", "scan=blank",
    ]
    assert exported[0]["precursor_mz"] == "301.1234"
    assert exported[0]["search_precursor_mz"] == "300.1161"
    assert exported[0]["top1_cosine"] == "0.8765"
    assert exported[0]["top1_mass_ppm"] == "-1.25"
    assert exported[1]["top1_compound_name"] == ""
    assert exported[1]["retention_time"] == ""
    assert exported[1]["collision_energy"] == ""
    assert exported[2]["blank_status"] == "flagged"
    assert [json.loads(row["native_record_json"]) for row in exported] == rows


def test_complete_v2_rejects_report_and_default_reader_is_unchanged(monkeypatch):
    rows = _saved_rows()
    monkeypatch.setitem(
        server._tasks, "native-events-mode",
        {"id": "native-events-mode", "filename": "case1.mzML", "status": "completed"},
    )
    monkeypatch.setattr(server, "_get_task_results", lambda task_id: rows)
    monkeypatch.setattr(server, "_load_expected_compounds_payload", lambda *a, **k: {})
    monkeypatch.setattr(server, "_load_reference_standard_qc_payload", lambda *a, **k: {})
    monkeypatch.setattr(server, "_compound_conclusions_for_export", lambda *a, **k: None)
    monkeypatch.setattr(server, "_ms2_csv_text", lambda *a, **k: "reader-v3-sentinel\n")
    client = TestClient(server.app)

    rejected = client.get(
        "/api/tasks/native-events-mode/download/csv"
        "?format_type=report&result_layer=ms2&ms2_schema=complete_v2"
    )
    assert rejected.status_code == 400
    assert "format_type=results" in rejected.json()["detail"]

    default = client.get("/api/tasks/native-events-mode/download/csv")
    assert default.status_code == 200
    assert default.text == "\ufeffreader-v3-sentinel\n"
    assert default.headers["x-mssystem-ms2-schema"] == "reader_v3"
    assert "_ms2_all.csv" in default.headers["content-disposition"]


def test_complete_v2_reads_saved_artifact_before_display_claim_projection(monkeypatch, tmp_path):
    import copy
    rows = _saved_rows()
    (tmp_path / "native-saved_results.json").write_text(json.dumps(rows))
    projected = copy.deepcopy(rows)
    projected[0]["match"]["final_claim_decision"] = "candidate_only"
    monkeypatch.setattr(server, "TASK_DIR", tmp_path)
    monkeypatch.setitem(server._tasks, "native-saved", {"filename": "saved.mzML", "status": "completed"})
    monkeypatch.setattr(server, "_get_task_results", lambda task_id: projected)
    response = TestClient(server.app).get("/api/tasks/native-saved/download/csv?ms2_schema=complete_v2")
    assert response.status_code == 200
    exported = list(csv.DictReader(io.StringIO(response.text.lstrip("\ufeff"))))
    assert [json.loads(row["native_record_json"]) for row in exported] == rows
