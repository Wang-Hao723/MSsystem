from __future__ import annotations

import json

from fastapi.testclient import TestClient

from web import server


def _row(collision_energy: float, strength: float) -> dict:
    candidate = {
        "compound_name": "Synthetic flagship candidate",
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
        "precursor_feature_id": "feature-1",
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


def test_flagship_endpoint_materializes_and_pages_primary_decisions(
    tmp_path, monkeypatch
):
    task_id = "flagship-api"
    source = tmp_path / f"{task_id}_results.json"
    source.write_text(
        json.dumps([_row(10.0, 0.90), _row(20.0, 0.88)]),
        encoding="utf-8",
    )
    monkeypatch.setattr(server, "TASK_DIR", tmp_path)
    monkeypatch.setattr(
        server,
        "_tasks",
        {task_id: {"id": task_id, "status": "completed", "flagship_summary": {}}},
    )
    monkeypatch.setattr(server, "_save_tasks", lambda: None)

    response = TestClient(server.app).get(
        f"/api/tasks/{task_id}/flagship?offset=0&limit=1"
    )

    assert response.status_code == 200
    payload = response.json()
    assert payload["primary_identity_selector"] is True
    assert payload["legacy_scan_results_unchanged"] is True
    assert payload["expert_rank_influence"] is False
    assert payload["entities_total"] == 1
    assert payload["has_more"] is False
    entity = payload["entities"][0]
    assert entity["winner"]["identity_key"] == "std2d:LFQSCWFLJHTTHZ"
    assert payload["entities"][0]["distinct_collision_energy_bins"] == 2
    artifact = tmp_path / f"{task_id}_v11_flagship.json"
    assert artifact.is_file()
    assert server._tasks[task_id]["flagship_summary"]["production_identity_selector"] is True

    structure = TestClient(server.app).get(
        f"/api/tasks/{task_id}/flagship/structure",
        params={"entity_id": entity["entity_id"], "format": "svg"},
    )
    assert structure.status_code == 200
    assert structure.headers["content-type"].startswith("image/svg+xml")
    assert b"<svg" in structure.content

    missing = TestClient(server.app).get(
        f"/api/tasks/{task_id}/flagship?entity_id=entity%3Amissing"
    )
    assert missing.status_code == 404
