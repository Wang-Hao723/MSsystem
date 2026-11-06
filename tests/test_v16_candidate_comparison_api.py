from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path

import pytest
from fastapi import FastAPI
from fastapi.testclient import TestClient

from web import v16_candidate_comparison as module

A, B = "A" * 14, "B" * 14


def setup_service(tmp_path, monkeypatch):
    source = tmp_path / "saved_results.json"
    source.write_text(json.dumps([{"spectrum_id": "q1", "precursor_mz": 200., "ion_mode": "NEG",
        "raw_peaks": [[70., 100.], [90., 60.], [110., 20.]], "raw_peaks_complete": True,
        "ground_truth": "must not enter scoring", "compound_name": "must not enter scoring"}]))
    artifact = tmp_path / "saved_v16_flagship.json"
    entity = {"entity_id": "entity:saved:one", "member_result_indices": [0],
        "candidates": [{"identity_key": "std2d:" + A, "entity_rank": 1}],
        "winner": {"identity_key": "std2d:" + A}}
    artifact.write_text(json.dumps({"entities": [entity]}))
    database = tmp_path / "library.db"
    database.write_bytes(b"test-only-database-binding")
    tasks = {"saved": {"status": "completed"}}
    calls = []

    class Retriever:
        db_path = database
        library_snapshot_id = "frozen-test"
        policy = {"test": "fixed"}

        def retrieve_entity(self, entity_id, rows, *, existing_candidates):
            return {"expanded_candidates": [{"connectivity_identity": B}]}

    def scorer(rows, pool, retriever, **kwargs):
        calls.append({"rows": copy.deepcopy(rows), "pool": list(pool), "kwargs": kwargs})
        return {"candidate_count": len(pool), "comparison_evidence": {"primary_ranking_changed": False}}

    monkeypatch.setattr(module, "rescore_complete_panel", scorer)
    service = module.CandidateComparisonService(task_directory=lambda: tmp_path,
        task_lookup=tasks.get, artifact_path=lambda directory, task_id: directory / f"{task_id}_v16_flagship.json",
        artifact_current=lambda path: True,
        read_entity=lambda path, entity_id: json.loads(path.read_text()),
        retriever_factory=Retriever, runtime_snapshot=lambda: "runtime-test")
    app = FastAPI()
    module.register_candidate_comparison_routes(app, service)
    return TestClient(app), service, tasks, source, artifact, database, calls


def request(client, **body):
    return client.post("/api/tasks/saved/candidate-comparison", json={"entity_id": "entity:saved:one", **body})


def test_endpoint_is_read_only_and_projects_labels_out(tmp_path, monkeypatch):
    client, _, tasks, source, artifact, database, calls = setup_service(tmp_path, monkeypatch)
    snapshots = {path: path.read_bytes() for path in (source, artifact, database)}
    before = copy.deepcopy(tasks)
    response = request(client)
    assert response.status_code == 200, response.text
    data = response.json()
    assert data["candidate_display_order"] == [A, B]
    assert data["production_results_modified"] is False
    assert data["primary_ranking_changed"] is False
    assert data["fragment_only_sorting_accepted_for_production"] is False
    assert calls[0]["pool"] == [A, B]
    assert "ground_truth" not in calls[0]["rows"][0]
    assert "compound_name" not in calls[0]["rows"][0]
    assert tasks == before
    assert {path: path.read_bytes() for path in snapshots} == snapshots


def test_completed_comparison_is_cached_and_input_changes_invalidate_it(tmp_path, monkeypatch):
    client, service, _, source, _, database, calls = setup_service(tmp_path, monkeypatch)
    assert request(client).json()["cache_hit"] is False
    assert request(client).json()["cache_hit"] is True
    assert len(calls) == 1
    source.write_text(source.read_text() + "\n")
    assert request(client).json()["cache_hit"] is False
    database.write_bytes(database.read_bytes() + b"changed")
    assert request(client).json()["cache_hit"] is False
    assert len(calls) == 3
    assert len(service.cache) <= 3


def test_running_task_is_not_restarted_or_analyzed(tmp_path, monkeypatch):
    client, _, tasks, _, _, _, calls = setup_service(tmp_path, monkeypatch)
    tasks["saved"]["status"] = "running"
    assert request(client).status_code == 409
    assert not calls


def test_busy_comparison_does_not_spawn_another_worker(tmp_path, monkeypatch):
    client, service, _, _, _, _, calls = setup_service(tmp_path, monkeypatch)
    service.lock.acquire()
    try:
        assert request(client).status_code == 429
        assert not calls
    finally:
        service.lock.release()
    assert request(client).status_code == 200


@pytest.mark.parametrize("body", [{"apply_rank_override": True}, {"ground_truth": A}, {"candidate_pool": [A]}])
def test_request_cannot_inject_truth_change_candidate_pool_or_enable_ranking(tmp_path, monkeypatch, body):
    client, _, _, _, _, _, calls = setup_service(tmp_path, monkeypatch)
    assert request(client, **body).status_code == 422
    assert not calls


def test_missing_or_foreign_entity_is_rejected(tmp_path, monkeypatch):
    client, _, _, _, _, _, calls = setup_service(tmp_path, monkeypatch)
    assert request(client, entity_id="entity:other:one").status_code == 404
    assert client.post("/api/tasks/unknown/candidate-comparison", json={"entity_id": "entity:saved:one"}).status_code == 404
    assert not calls


@pytest.mark.parametrize("indices,status", [([], 422), ([-1], 422), ([True], 422), (["0"], 422), ([1], 409), (list(range(201)), 422)])
def test_invalid_missing_or_excessive_query_indices_are_not_truncated(tmp_path, monkeypatch, indices, status):
    client, _, _, _, artifact, _, calls = setup_service(tmp_path, monkeypatch)
    payload = json.loads(artifact.read_text())
    payload["entities"][0]["member_result_indices"] = indices
    artifact.write_text(json.dumps(payload))
    assert request(client).status_code == status
    assert not calls


def test_stale_flagship_is_not_materialized_by_comparison(tmp_path, monkeypatch):
    client, service, _, _, _, _, calls = setup_service(tmp_path, monkeypatch)
    service.artifact_current = lambda _: False
    assert request(client).status_code == 409
    assert not calls


def test_vetoed_candidate_is_not_reintroduced_by_recall(tmp_path, monkeypatch):
    client, _, _, _, artifact, _, calls = setup_service(tmp_path, monkeypatch)
    payload = json.loads(artifact.read_text())
    payload["entities"][0]["candidates"].append({"identity_key": "std2d:" + B, "vetoed": True})
    artifact.write_text(json.dumps(payload))
    data = request(client).json()
    assert calls[0]["pool"] == [A]
    assert data["excluded_vetoed_connectivities"] == [B]


def test_abstained_entity_does_not_get_a_fabricated_primary_winner(tmp_path, monkeypatch):
    client, _, _, _, artifact, _, calls = setup_service(tmp_path, monkeypatch)
    payload = json.loads(artifact.read_text())
    payload["entities"][0]["winner"] = None
    artifact.write_text(json.dumps(payload))
    assert request(client).status_code == 200
    assert calls[0]["kwargs"]["baseline_connectivity"] is None


def test_changed_inputs_discard_computation_and_release_lock(tmp_path, monkeypatch):
    client, service, _, source, _, _, calls = setup_service(tmp_path, monkeypatch)
    def mutating_scorer(*args, **kwargs):
        source.write_text(source.read_text() + " ")
        return {}
    monkeypatch.setattr(module, "rescore_complete_panel", mutating_scorer)
    assert request(client).status_code == 409
    assert not service.lock.locked()
    assert not service.cache


def test_route_is_registered_in_real_app():
    from web import server
    assert any(getattr(route, "path", "") == "/api/tasks/{task_id}/candidate-comparison" for route in server.app.routes)


def test_export_download_is_bound_to_the_displayed_comparison(tmp_path, monkeypatch):
    client, _, _, _, _, _, calls = setup_service(tmp_path, monkeypatch)
    data = request(client).json()
    exported = client.get("/api/tasks/saved/candidate-comparison/export", params={
        "entity_id": data["entity_id"], "comparison_id": data["comparison_id"]})
    assert exported.status_code == 200
    assert exported.headers["content-disposition"] == 'attachment; filename="saved_v16_candidate_comparison.json"'
    assert exported.headers["cache-control"].startswith("no-store")
    assert exported.json()["panel"] == data["panel"]
    assert len(calls) == 1


def test_export_refuses_a_stale_displayed_result(tmp_path, monkeypatch):
    client, _, _, source, _, _, _ = setup_service(tmp_path, monkeypatch)
    data = request(client).json()
    source.write_text(source.read_text() + " ")
    response = client.get("/api/tasks/saved/candidate-comparison/export", params={
        "entity_id": data["entity_id"], "comparison_id": data["comparison_id"]})
    assert response.status_code == 409


@pytest.mark.parametrize("identifier", ["", "not-a-source-digest"])
def test_export_requires_an_explicit_comparison_digest(tmp_path, monkeypatch, identifier):
    client, _, _, _, _, _, calls = setup_service(tmp_path, monkeypatch)
    assert client.get("/api/tasks/saved/candidate-comparison/export", params={
        "entity_id": "entity:saved:one", "comparison_id": identifier}).status_code == 422
    assert not calls
