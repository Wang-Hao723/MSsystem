from __future__ import annotations

import json

from fastapi.testclient import TestClient

from web import server


def _legacy_result() -> dict:
    return {
        "spectrum_id": "scan=7",
        "precursor_mz": 301.0712,
        "retention_time": 1.25,
        "ion_mode": "NEG",
        "match": {
            "compound_name": "Candidate",
            "inchikey": "AAAAAAAAAAAAAA-UHFFFAOYSA-N",
            "formula": "C15H10O7",
            "smiles": "O=c1ccoc2ccccc12",
            "adduct": "[M-H]-",
            "source_tier": "experimental",
            "cosine": 0.91,
            "entropy": 0.87,
            "final_score": 0.73,
            "rank_score": 0.73,
            "confidence": "High",
            "isopair_evidence": {
                "evaluated": False,
                "runtime_influence": False,
                "ranking_changed": False,
            },
            "evidence_governance": {
                "policy_id": "evigate-ms-locked-v1-2026-08-02",
                "claim_status": "provisional_needs_independent_injection",
                "claim_scope": "connectivity_level_library_candidate",
                "release_eligible_before_replication": True,
                "ranking_changed": False,
                "nodes": [],
            },
        },
    }


def test_identification_contract_api_materializes_without_mutating_legacy(
    tmp_path, monkeypatch
):
    task_id = "contract-task"
    result_path = tmp_path / f"{task_id}_results.json"
    result_path.write_text(json.dumps([_legacy_result()]), encoding="utf-8")
    before = result_path.read_bytes()
    monkeypatch.setattr(server, "TASK_DIR", tmp_path)
    monkeypatch.setattr(
        server,
        "_tasks",
        {
            task_id: {
                "id": task_id,
                "status": "completed",
                "filename": "sample.mzML",
            }
        },
    )
    monkeypatch.setattr(server, "_save_tasks", lambda: None)
    client = TestClient(server.app)

    response = client.get(
        f"/api/tasks/{task_id}/identification-records?limit=1"
    )

    assert response.status_code == 200
    payload = response.json()
    assert payload["metadata"]["record_count"] == 1
    assert payload["metadata"]["production_ranking_changed_records"] == 0
    assert payload["records"][0]["candidate"]["inchikey"] == (
        "AAAAAAAAAAAAAA-UHFFFAOYSA-N"
    )
    assert payload["records"][0]["decision"]["scope"] == "2d_connectivity"
    assert payload["legacy_results_unchanged"] is True
    assert payload["expected_compounds_counted"] is False
    assert result_path.read_bytes() == before
    assert (tmp_path / f"{task_id}_identification_records.json").is_file()


def test_identification_contract_api_can_refuse_lazy_materialization(
    tmp_path, monkeypatch
):
    task_id = "contract-task-no-cache"
    (tmp_path / f"{task_id}_results.json").write_text(
        json.dumps([_legacy_result()]), encoding="utf-8"
    )
    monkeypatch.setattr(server, "TASK_DIR", tmp_path)
    monkeypatch.setattr(
        server,
        "_tasks",
        {
            task_id: {
                "id": task_id,
                "status": "completed",
                "filename": "sample.mzML",
            }
        },
    )
    client = TestClient(server.app)

    response = client.get(
        f"/api/tasks/{task_id}/identification-records?materialize=false"
    )

    assert response.status_code == 404
    assert "not been materialized" in response.json()["detail"]
