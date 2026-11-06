from __future__ import annotations

import json

from fastapi.testclient import TestClient

from pipeline.evograph_shadow import (
    evograph_sidecar_path,
    materialize_evograph_sidecar,
    read_evograph_sidecar,
)
from web import server


def _legacy_row() -> dict:
    return {
        "spectrum_id": "scan=11",
        "precursor_mz": 181.070664,
        "precursor_charge": 1,
        "ion_mode": "POS",
        "raw_peaks": [[163.060100, 100.0], [135.065200, 30.0]],
        "match": {
            "compound_name": "Candidate",
            "formula": "C6H12O6",
            "adduct": "[M+H]+",
            "final_score": 0.91,
            "rank_score": 0.91,
            "confidence": "High",
        },
    }


def test_evograph_sidecar_hashes_both_sides_and_preserves_legacy_result(tmp_path):
    task_id = "evograph-sidecar"
    legacy = tmp_path / f"{task_id}_results.json"
    legacy.write_text(json.dumps([_legacy_row()]), encoding="utf-8")
    before = legacy.read_bytes()
    artifact = evograph_sidecar_path(tmp_path, task_id)

    metadata = materialize_evograph_sidecar(task_id, legacy, artifact)
    payload = read_evograph_sidecar(artifact, limit=1)

    assert metadata["source_results_sha256"] == metadata["source_results_sha256_after"]
    assert metadata["source_results_unchanged"] is True
    assert metadata["production_ranking_changed_records"] == 0
    assert metadata["legacy_confidence_changed_records"] == 0
    assert metadata["high_promotions"] == 0
    assert metadata["safety_invariants_passed"] is True
    assert "match" not in payload["records"][0]
    assert legacy.read_bytes() == before
    assert json.loads(legacy.read_text())[0]["match"]["final_score"] == 0.91
    assert json.loads(legacy.read_text())[0]["match"]["confidence"] == "High"


def test_evograph_api_lazy_materialization_is_paginated_and_shadow_only(
    tmp_path, monkeypatch
):
    task_id = "evograph-api"
    legacy = tmp_path / f"{task_id}_results.json"
    legacy.write_text(json.dumps([_legacy_row()]), encoding="utf-8")
    before = legacy.read_bytes()
    monkeypatch.setattr(server, "TASK_DIR", tmp_path)
    monkeypatch.setattr(
        server,
        "_tasks",
        {
            task_id: {
                "id": task_id,
                "status": "completed",
                "filename": "sample.mgf",
                "result_run_id": "run-1",
            }
        },
    )
    monkeypatch.setattr(server, "_save_tasks", lambda: None)

    response = TestClient(server.app).get(
        f"/api/tasks/{task_id}/evograph-fragmentation?limit=1&result_index=0"
    )

    assert response.status_code == 200
    payload = response.json()
    assert payload["metadata"]["record_count"] == 1
    assert len(payload["records"]) == 1
    assert payload["records"][0]["result_index"] == 0
    assert payload["module_state"] == "shadow"
    assert payload["used_for_scoring"] is False
    assert payload["production_ranking_changed"] is False
    assert payload["legacy_confidence_changed"] is False
    assert payload["heuristic_marginals_are_probabilities"] is False
    assert payload["log_bayes_factor_available"] is False
    assert legacy.read_bytes() == before
    assert evograph_sidecar_path(tmp_path, task_id).is_file()


def test_evograph_api_can_refuse_materialization(tmp_path, monkeypatch):
    task_id = "evograph-api-no-cache"
    (tmp_path / f"{task_id}_results.json").write_text(
        json.dumps([_legacy_row()]), encoding="utf-8"
    )
    monkeypatch.setattr(server, "TASK_DIR", tmp_path)
    monkeypatch.setattr(
        server,
        "_tasks",
        {task_id: {"id": task_id, "status": "completed", "filename": "sample.mgf"}},
    )

    response = TestClient(server.app).get(
        f"/api/tasks/{task_id}/evograph-fragmentation?materialize=false"
    )

    assert response.status_code == 404
    assert "not been materialized" in response.json()["detail"]


def test_detailed_report_joins_evograph_by_immutable_result_index(
    tmp_path, monkeypatch
):
    task_id = "evograph-detailed-report"
    row = _legacy_row()
    row.update(
        {
            "retention_time": 1.25,
            "retention_time_unit": "minutes",
            "query_peaks": list(row["raw_peaks"]),
            "peaks_count": len(row["raw_peaks"]),
            "lib_peaks": [[163.0601, 90.0], [135.0652, 25.0]],
        }
    )
    row["match"].update(
        {
            "smiles": "OCC(O)C(O)C(O)CO",
            "classification": "carbohydrate",
            "cosine": 0.92,
            "mass_ppm": 1.0,
            "source_tier": "experimental",
            "dataset_source": "test",
        }
    )
    legacy = tmp_path / f"{task_id}_results.json"
    legacy.write_text(json.dumps([row]), encoding="utf-8")
    materialize_evograph_sidecar(
        task_id,
        legacy,
        evograph_sidecar_path(tmp_path, task_id),
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
                "results": [row],
                "analysis_warnings": [],
            }
        },
    )
    monkeypatch.setattr(server, "_save_tasks", lambda: None)
    monkeypatch.setattr(server, "_prepare_results_for_display", lambda value: value)
    monkeypatch.setattr(server, "_attach_claim_decisions", lambda *args: None)
    monkeypatch.setattr(server, "_load_expected_compounds_payload", lambda *args, **kwargs: {})

    response = TestClient(server.app).get(
        f"/api/tasks/{task_id}/download/html?format_type=detailed&language=zh"
    )

    assert response.status_code == 200
    assert "观测碎片组成图（EvoGraph-F）" in response.text
    assert "候选条件下的分子式守恒检查" in response.text
    assert "EvoGraph-F 不参与检索、排序、置信度或 High 升档" in response.text
