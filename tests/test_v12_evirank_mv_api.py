from __future__ import annotations

import hashlib
import json

from fastapi.testclient import TestClient

from web import server


def _row(collision_energy: float, strength: float) -> dict:
    candidate = {
        "compound_name": "Synthetic V12 candidate",
        "inchikey": "LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
        "formula": "C2H6O",
        "smiles": "CCO",
        "cosine": strength,
        "entropy": strength,
        "reverse_match": strength,
        "fragment_score": strength,
        "intensity_score": strength,
        "mass_ppm": 0.0,
        "final_score": strength,
        "rank_score": strength,
        "dataset_source": "synthetic-test",
    }
    return {
        "spectrum_id": f"scan-{collision_energy:g}",
        "precursor_feature_id": "feature-v12-1",
        "precursor_mz": 300.0,
        "retention_time": 1.0,
        "retention_time_unit": "minutes",
        "collision_energy": collision_energy,
        "collision_energy_unit": "eV",
        "ion_mode": "POS",
        "tic": 1000.0 + collision_energy,
        "raw_peaks": [[100.0, 100.0], [150.0, 50.0]],
        "query_peaks": [[100.0, 100.0], [150.0, 50.0]],
        "match": {**candidate, "top_candidates": [candidate]},
    }


def _prepare_completed_task(tmp_path, monkeypatch, task_id: str) -> None:
    source = tmp_path / f"{task_id}_results.json"
    source.write_text(
        json.dumps([_row(10.0, 0.90), _row(20.0, 0.88)]),
        encoding="utf-8",
    )
    monkeypatch.setattr(server, "TASK_DIR", tmp_path)
    monkeypatch.setattr(
        server,
        "_tasks",
        {
            task_id: {
                "id": task_id,
                "status": "completed",
                "flagship_summary": {},
                "evirank_mv_summary": {},
            }
        },
    )
    monkeypatch.setattr(server, "_save_tasks", lambda: None)


def test_evirank_endpoint_materializes_shadow_without_mutating_v11(
    tmp_path, monkeypatch
):
    task_id = "evirank-api"
    _prepare_completed_task(tmp_path, monkeypatch, task_id)
    client = TestClient(server.app)

    response = client.get(
        f"/api/tasks/{task_id}/evirank-mv",
        params={"offset": 0, "limit": 1, "detection_method_id": "qtof-pos-cid-v1"},
    )

    assert response.status_code == 200
    payload = response.json()
    assert payload["metadata"]["algorithm_name"] == "V12-Flagship / EviRank-MV"
    assert payload["metadata"]["deployment_state"] == "shadow"
    assert payload["metadata"]["detection_method_id"] == "qtof-pos-cid-v1"
    assert payload["metadata"]["residual_model_state"] == "unavailable_no_frozen_model"
    assert payload["shadow_challenger"] is True
    assert payload["production_identity_selector"] is False
    assert payload["v11_production_result_unchanged"] is True
    assert payload["score_is_probability"] is False
    assert payload["entities_total"] == 1
    assert payload["entities"][0]["ranking_changed_from_v11"] is False

    v11_path = tmp_path / f"{task_id}_v11_flagship.json"
    v12_path = tmp_path / f"{task_id}_v12_evirank_mv.json"
    assert v11_path.is_file()
    assert v12_path.is_file()
    v11_hash = hashlib.sha256(v11_path.read_bytes()).hexdigest()
    v12_hash = hashlib.sha256(v12_path.read_bytes()).hexdigest()
    assert server._tasks[task_id]["flagship_summary"]["production_identity_selector"] is True
    assert server._tasks[task_id]["evirank_mv_summary"]["production_identity_selector"] is False

    repeated = client.get(f"/api/tasks/{task_id}/evirank-mv")
    assert repeated.status_code == 200
    assert hashlib.sha256(v11_path.read_bytes()).hexdigest() == v11_hash
    assert hashlib.sha256(v12_path.read_bytes()).hexdigest() == v12_hash

    resource_response = client.get(f"/api/tasks/{task_id}/resources")
    assert resource_response.status_code == 200
    resource = next(
        item
        for item in resource_response.json()["resources"]
        if item["resource"] == "v12_evirank_mv_shadow"
    )
    assert resource["status"] == "available"
    assert resource["contract_role"].endswith("no_production_effect")


def test_evirank_endpoint_refuses_silent_method_rebinding(tmp_path, monkeypatch):
    task_id = "evirank-method-binding"
    _prepare_completed_task(tmp_path, monkeypatch, task_id)
    client = TestClient(server.app)
    first = client.get(
        f"/api/tasks/{task_id}/evirank-mv",
        params={"detection_method_id": "orbitrap-pos-hcd-v1"},
    )
    assert first.status_code == 200
    artifact = tmp_path / f"{task_id}_v12_evirank_mv.json"
    before = artifact.read_bytes()

    rebound = client.get(
        f"/api/tasks/{task_id}/evirank-mv",
        params={"detection_method_id": "qtof-neg-cid-v1"},
    )

    assert rebound.status_code == 409
    assert "silently rebinding" in rebound.json()["detail"]
    assert artifact.read_bytes() == before


def test_evirank_endpoint_honors_nonmaterializing_and_entity_filters(
    tmp_path, monkeypatch
):
    task_id = "evirank-read-contract"
    _prepare_completed_task(tmp_path, monkeypatch, task_id)
    client = TestClient(server.app)

    absent = client.get(
        f"/api/tasks/{task_id}/evirank-mv",
        params={"materialize": False},
    )
    assert absent.status_code == 404
    assert not (tmp_path / f"{task_id}_v11_flagship.json").exists()
    assert not (tmp_path / f"{task_id}_v12_evirank_mv.json").exists()

    created = client.get(f"/api/tasks/{task_id}/evirank-mv")
    assert created.status_code == 200
    entity_id = created.json()["entities"][0]["entity_id"]
    selected = client.get(
        f"/api/tasks/{task_id}/evirank-mv",
        params={"entity_id": entity_id},
    )
    assert selected.status_code == 200
    assert selected.json()["entities"][0]["entity_id"] == entity_id

    missing = client.get(
        f"/api/tasks/{task_id}/evirank-mv",
        params={"entity_id": "entity:missing"},
    )
    assert missing.status_code == 404

