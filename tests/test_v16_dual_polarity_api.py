from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path

import pytest
from fastapi.testclient import TestClient
from pydantic import ValidationError

from web import server


_PAIR_ENDPOINT = "/api/workflows/dual-polarity-shadow"
_PRODUCTION_FALSE_INVARIANTS = (
    "production_top1_changed",
    "production_top5_changed",
    "production_score_changed",
    "production_decision_changed",
    "production_confidence_changed",
    "production_compound_count_changed",
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _task(task_id: str, polarity: str, feature_xml: Path) -> dict:
    misleading = "NEG" if polarity == "POS" else "POS"
    return {
        "id": task_id,
        "status": "completed",
        "ion_mode": polarity,
        "ion_mode_counts": {polarity: 4},
        "filename": f"filename-claims-{misleading}.mzML",
        "untargeted_workflow": {"feature_xml": str(feature_xml)},
    }


def _shadow_payload() -> dict:
    return {
        "schema_version": "mssystem.v16.dual-polarity-shadow.sidecar/1.0.0",
        "metadata": {
            "module_state": "shadow",
            "sidecar_semantic_sha256": "a" * 64,
        },
        "summary": {
            "status": "evaluated_no_edge",
            "neutral_entity_count": 0,
            "compound_count_increment": 0,
        },
        "candidate_independent_graph": {
            "semantic_sha256": "b" * 64,
            "anchors": {"POS": [], "NEG": []},
        },
        "neutral_entities": [],
        "invariants": {
            **{field: False for field in _PRODUCTION_FALSE_INVARIANTS},
            "sample_filename_used_for_pairing": False,
        },
    }


def _prepare_pair(tmp_path: Path, monkeypatch):
    pos_feature = tmp_path / "pos.featureXML"
    neg_feature = tmp_path / "neg.featureXML"
    pos_feature.write_bytes(b"<featureMap polarity='POS'/>")
    neg_feature.write_bytes(b"<featureMap polarity='NEG'/>")

    task_ids = {"POS": "explicit-pos", "NEG": "explicit-neg"}
    tasks = {
        task_ids["POS"]: _task(task_ids["POS"], "POS", pos_feature),
        task_ids["NEG"]: _task(task_ids["NEG"], "NEG", neg_feature),
    }
    source_paths: dict[str, dict[str, Path]] = {}
    for polarity, task_id in task_ids.items():
        result_path = tmp_path / f"{task_id}_results.json"
        result_path.write_text(
            json.dumps([{"ion_mode": polarity, "spectrum_id": task_id}]),
            encoding="utf-8",
        )
        flagship_path = server._flagship_artifact_path(tmp_path, task_id)
        flagship_path.write_text(
            json.dumps(
                {
                    "metadata": {
                        "task_id": task_id,
                        "artifact_role": "primary_v15_test_fixture",
                        "entity_semantic_sha256": polarity.lower() * 16,
                    },
                    "entities": [],
                }
            ),
            encoding="utf-8",
        )
        source_paths[task_id] = {
            "results": result_path,
            "flagship": flagship_path,
            "feature_xml": pos_feature if polarity == "POS" else neg_feature,
        }

    calls: list[dict] = []

    def fake_builder(task_a, results_a, task_b, results_b, **kwargs):
        calls.append(
            {
                "task_a": copy.deepcopy(task_a),
                "results_a": copy.deepcopy(results_a),
                "task_b": copy.deepcopy(task_b),
                "results_b": copy.deepcopy(results_b),
                "kwargs": copy.deepcopy(kwargs),
            }
        )
        return _shadow_payload()

    def forbidden(*_args, **_kwargs):
        raise AssertionError("the read-only shadow endpoint used a mutating/materializing helper")

    monkeypatch.setattr(server, "TASK_DIR", tmp_path)
    monkeypatch.setattr(server, "_tasks", tasks)
    monkeypatch.setattr(
        server,
        "_v16_current_flagship_metadata",
        lambda path, *, task_id: json.loads(
            Path(path).read_text(encoding="utf-8")
        )["metadata"],
    )
    monkeypatch.setattr(
        server,
        "_read_flagship_artifact",
        lambda path: json.loads(Path(path).read_text(encoding="utf-8")),
    )
    monkeypatch.setattr(server, "_build_v16_dual_polarity_shadow", fake_builder)
    monkeypatch.setattr(server, "_get_task_results", forbidden)
    monkeypatch.setattr(server, "_save_tasks", forbidden)
    monkeypatch.setattr(server, "_current_runtime_snapshot_id", lambda: "v16-test")
    return task_ids, tasks, source_paths, calls


def _request(task_ids: dict[str, str]) -> dict:
    return {
        "positive_task_id": task_ids["POS"],
        "negative_task_id": task_ids["NEG"],
        "same_sample_confirmed": True,
        "same_lc_method_confirmed": True,
        "pair_label": "explicit engineering pair",
    }


@pytest.mark.parametrize(
    "missing_field",
    ["same_sample_confirmed", "same_lc_method_confirmed"],
)
def test_request_model_requires_both_operator_confirmations(missing_field):
    values = {
        "positive_task_id": "explicit-pos",
        "negative_task_id": "explicit-neg",
        "same_sample_confirmed": True,
        "same_lc_method_confirmed": True,
    }
    values.pop(missing_field)

    with pytest.raises(ValidationError):
        server.DualPolarityShadowRequest.model_validate(values)


def test_post_uses_exact_task_polarity_and_is_a_read_only_sidecar(
    tmp_path, monkeypatch
):
    task_ids, tasks, source_paths, calls = _prepare_pair(tmp_path, monkeypatch)
    tasks_before = copy.deepcopy(tasks)
    source_hashes_before = {
        (task_id, kind): _sha256(path)
        for task_id, paths in source_paths.items()
        for kind, path in paths.items()
        if kind in {"results", "flagship"}
    }

    response = TestClient(server.app).post(_PAIR_ENDPOINT, json=_request(task_ids))

    assert response.status_code == 200, response.text
    envelope = response.json()
    payload = envelope["sidecar"]
    assert envelope["status"] == "available"
    assert envelope["reused"] is False
    assert Path(envelope["artifact"]).parent == tmp_path / "dual_polarity_shadow"
    assert len(calls) == 1
    assert calls[0]["task_a"]["ion_mode"] == "POS"
    assert calls[0]["task_b"]["ion_mode"] == "NEG"
    assert "filename-claims-NEG" in calls[0]["task_a"]["filename"]
    assert "filename-claims-POS" in calls[0]["task_b"]["filename"]
    assert tasks == tasks_before
    assert {
        (task_id, kind): _sha256(path)
        for task_id, paths in source_paths.items()
        for kind, path in paths.items()
        if kind in {"results", "flagship"}
    } == source_hashes_before
    for field in _PRODUCTION_FALSE_INVARIANTS:
        assert payload["invariants"][field] is False
    assert payload["summary"]["compound_count_increment"] == 0

    sidecar_root = tmp_path / "dual_polarity_shadow"
    sidecars = list(sidecar_root.glob("*.json"))
    assert len(sidecars) == 1
    assert sidecars[0].parent == sidecar_root
    assert not list(tmp_path.glob("*results*.tmp*"))
    assert not list(tmp_path.glob("*flagship*.tmp*"))


def test_post_rejects_filename_inference_noncompleted_and_same_task(
    tmp_path, monkeypatch
):
    task_ids, tasks, _source_paths, calls = _prepare_pair(tmp_path, monkeypatch)
    client = TestClient(server.app)

    tasks[task_ids["NEG"]]["ion_mode"] = "POS"
    wrong_polarity = client.post(_PAIR_ENDPOINT, json=_request(task_ids))
    assert wrong_polarity.status_code == 409
    assert not calls

    tasks[task_ids["NEG"]]["ion_mode"] = "NEG"
    tasks[task_ids["NEG"]]["status"] = "running"
    noncompleted = client.post(_PAIR_ENDPOINT, json=_request(task_ids))
    assert noncompleted.status_code == 409
    assert not calls

    same = _request(task_ids)
    same["negative_task_id"] = task_ids["POS"]
    duplicate = client.post(_PAIR_ENDPOINT, json=same)
    assert duplicate.status_code == 422
    assert not calls
    assert not (tmp_path / "dual_polarity_shadow").exists()


@pytest.mark.parametrize("missing_kind", ["results", "flagship", "feature_xml"])
def test_post_fails_closed_when_a_required_source_is_missing(
    tmp_path, monkeypatch, missing_kind
):
    task_ids, _tasks, source_paths, calls = _prepare_pair(tmp_path, monkeypatch)
    source_paths[task_ids["NEG"]][missing_kind].unlink()

    response = TestClient(server.app).post(_PAIR_ENDPOINT, json=_request(task_ids))

    assert response.status_code == 409
    assert not calls
    assert not (tmp_path / "dual_polarity_shadow").exists()


def test_post_fails_closed_for_stale_v15_flagship(tmp_path, monkeypatch):
    task_ids, _tasks, _source_paths, calls = _prepare_pair(tmp_path, monkeypatch)
    def stale(_path, *, task_id):
        raise RuntimeError(
            f"Task {task_id} current V15 flagship artifact is absent or stale"
        )

    monkeypatch.setattr(server, "_v16_current_flagship_metadata", stale)

    response = TestClient(server.app).post(_PAIR_ENDPOINT, json=_request(task_ids))

    assert response.status_code == 409
    assert not calls
    assert not (tmp_path / "dual_polarity_shadow").exists()


@pytest.mark.parametrize(
    "malicious_pair_id",
    ["..%2Ftasks_index.json", "%2E%2E%2Ftasks_index.json", ".hidden"],
)
def test_get_rejects_pair_id_path_traversal(tmp_path, monkeypatch, malicious_pair_id):
    monkeypatch.setattr(server, "TASK_DIR", tmp_path)
    secret = tmp_path / "tasks_index.json"
    secret.write_text('{"secret": true}', encoding="utf-8")

    response = TestClient(server.app).get(
        f"{_PAIR_ENDPOINT}/{malicious_pair_id}"
    )

    assert response.status_code in {400, 404, 422}
    assert "secret" not in response.text
