from __future__ import annotations

import json

from fastapi.testclient import TestClient

from governance.claim_compiler import default_task_claim_context
from pipeline.claim_contract import claim_contract_path, write_claim_contract
from pipeline.result_contract import identification_contract_path, write_identification_contract
from web import server


def _row():
    return {
        "spectrum_id": "scan=1",
        "precursor_mz": 301.0,
        "match": {
            "compound_name": "Candidate",
            "inchikey": "AAAAAAAAAAAAAA-UHFFFAOYSA-N",
            "rank_score": 0.9,
            "confidence": "High",
            "evidence_governance": {
                "policy_id": "v1",
                "claim_status": "provisional_needs_independent_injection",
                "claim_scope": "connectivity_level_library_candidate",
                "release_eligible_before_replication": True,
                "nodes": [],
            },
        },
    }


def test_claim_api_materializes_and_keeps_legacy_high_separate(tmp_path, monkeypatch):
    task_id = "claim-task"
    identification = identification_contract_path(tmp_path, task_id)
    write_identification_contract(task_id, [_row()], identification)
    monkeypatch.setattr(server, "TASK_DIR", tmp_path)
    monkeypatch.setattr(
        server,
        "_tasks",
        {task_id: {"id": task_id, "status": "completed", "filename": "sample.mzML"}},
    )
    monkeypatch.setattr(server, "_save_tasks", lambda: None)

    response = TestClient(server.app).get(f"/api/tasks/{task_id}/claims?limit=1")

    assert response.status_code == 200
    payload = response.json()
    assert payload["metadata"]["released_count"] == 0
    assert payload["decisions"][0]["legacy_confidence"] == "High"
    assert payload["decisions"][0]["decision"]["released"] is False
    assert payload["legacy_high_is_project_qvalue_tier"] is False
    assert payload["qvalues_are_probabilities"] is False
    assert claim_contract_path(tmp_path, task_id).is_file()
