from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path

import pytest

from scripts.evaluate_v16_iimn_feature_centroid_shadow import (
    _atomic_write_json,
    _explicit_task_ids,
    _stream_candidate_independent_rows,
    evaluate_tasks,
)


def _write_inputs(tmp_path: Path, task_id: str = "task-a"):
    task_root = tmp_path / "tasks"
    task_root.mkdir()
    raw = tmp_path / "source.mzML"
    raw.write_bytes(b"source-ms-data")
    result = task_root / f"{task_id}_results.json"
    result.write_text(
        json.dumps(
            [
                {
                    "spectrum_id": "scan-1",
                    "precursor_mz": 100.0,
                    "retention_time": 1.0,
                    "ion_mode": "NEG",
                    "lcms_feature_id": "LF1",
                    "lcms_feature_mz": 100.002,
                    "lcms_feature_link_status": "linked_openms_feature",
                    "match": {"compound_name": "forbidden", "final_score": 0.9},
                    "candidates": [{"name": "forbidden"}],
                }
            ]
        ),
        encoding="utf-8",
    )
    index = tmp_path / "tasks_index.json"
    index.write_text(
        json.dumps(
            {
                task_id: {
                    "id": task_id,
                    "status": "completed",
                    "ion_mode": "NEG",
                    "filename": "MISLEADING_POS_9.mzML",
                    "analysis_file_path": str(raw),
                }
            }
        ),
        encoding="utf-8",
    )
    return index, task_root, raw, result


def _fake_builder(task, rows, *, source_results_path, source_data_path):
    result_hash = hashlib.sha256(Path(source_results_path).read_bytes()).hexdigest()
    raw_hash = hashlib.sha256(Path(source_data_path).read_bytes()).hexdigest()
    task_id = task["id"]
    return {
        "metadata": {
            "task_id": task_id,
            "polarity": task["ion_mode"],
            "task_metadata_sha256": "task-sha",
            "candidate_independent_projection_sha256": "projection-sha",
            "sidecar_semantic_sha256": f"sidecar-{task_id}",
        },
        "source_binding": {
            "results_file": {"path": str(source_results_path), "sha256": result_hash},
            "source_data_file": {"path": str(source_data_path), "sha256": raw_hash},
            "current_trace_semantic_sha256": "trace-current",
            "feature_centroid_trace_semantic_sha256": "trace-shadow",
            "source_files_unchanged": True,
        },
        "trace_center_audit": {
            "counts": {
                "input_rows": len(rows),
                "rows_using_feature_centroid": 1,
                "rows_using_precursor_fallback": 0,
                "feature_centroid_outside_current_xic_window": 1,
            }
        },
        "comparison": {
            "current_precursor_center": {
                "summary": {
                    "network_count": 1,
                    "networked_features": 2,
                    "edge_count": 1,
                },
                "network_semantic_sha256": "network-current",
            },
            "feature_centroid_center": {
                "summary": {
                    "network_count": 3,
                    "networked_features": 6,
                    "edge_count": 3,
                },
                "network_semantic_sha256": "network-shadow",
            },
            "delta": {
                "network_count": 2,
                "networked_feature_count": 4,
                "edge_count": 2,
                "added_network_fingerprints": ["a", "b"],
                "removed_network_fingerprints": [],
            },
        },
        "invariants": {
            "candidate_identity_used_for_network_construction": False,
            "source_files_unchanged": True,
        },
    }


def test_explicit_task_ids_accept_parameters_or_json_and_reject_duplicates(tmp_path):
    task_file = tmp_path / "task_ids.json"
    task_file.write_text(json.dumps({"task_ids": ["b", "c"]}), encoding="utf-8")
    assert _explicit_task_ids(["a"], task_file) == ["a", "b", "c"]
    with pytest.raises(ValueError, match="Duplicate"):
        _explicit_task_ids(["a", "a"], None)
    with pytest.raises(ValueError, match="explicit task ID"):
        _explicit_task_ids([], None)


def test_stream_projection_drops_candidate_objects(tmp_path):
    index, task_root, _, result = _write_inputs(tmp_path)
    del index, task_root
    rows = _stream_candidate_independent_rows(
        result,
        allowed_fields=[
            "spectrum_id",
            "precursor_mz",
            "retention_time",
            "ion_mode",
            "lcms_feature_id",
            "lcms_feature_mz",
            "lcms_feature_link_status",
        ],
    )
    assert rows[0]["lcms_feature_mz"] == pytest.approx(100.002)
    assert "match" not in rows[0]
    assert "candidates" not in rows[0]


def test_evaluate_uses_task_metadata_polarity_not_misleading_filename(tmp_path):
    index, task_root, raw, result = _write_inputs(tmp_path)
    raw_before = raw.read_bytes()
    result_before = result.read_bytes()
    payload = evaluate_tasks(
        task_ids=["task-a"],
        tasks_index_path=index,
        tasks_dir=task_root,
        builder=_fake_builder,
    )
    assert payload["tasks"][0]["polarity"] == "NEG"
    assert payload["filename_used_for_polarity_or_pairing"] is False
    assert payload["pairing_performed"] is False
    assert payload["aggregate"]["ALL"]["network_count_delta"] == 2
    assert payload["aggregate"]["NEG"]["rows_using_feature_centroid"] == 1
    assert payload["aggregate"]["ALL"]["tasks_with_network_gain"] == 1
    assert raw.read_bytes() == raw_before
    assert result.read_bytes() == result_before


def test_repeat_check_requires_identical_sidecar_semantics(tmp_path):
    index, task_root, _, _ = _write_inputs(tmp_path)
    calls = 0

    def unstable(*args, **kwargs):
        nonlocal calls
        calls += 1
        payload = copy.deepcopy(_fake_builder(*args, **kwargs))
        payload["metadata"]["sidecar_semantic_sha256"] = f"run-{calls}"
        return payload

    with pytest.raises(RuntimeError, match="repeat semantic hash changed"):
        evaluate_tasks(
            task_ids=["task-a"],
            tasks_index_path=index,
            tasks_dir=task_root,
            repeat_check=True,
            builder=unstable,
        )


@pytest.mark.parametrize("failure", ["missing-task", "missing-result", "missing-raw"])
def test_missing_task_or_source_fails_before_output(tmp_path, failure):
    index, task_root, raw, result = _write_inputs(tmp_path)
    task_id = "task-a"
    if failure == "missing-task":
        task_id = "absent"
    elif failure == "missing-result":
        result.unlink()
    else:
        raw.unlink()
    with pytest.raises((KeyError, FileNotFoundError)):
        evaluate_tasks(
            task_ids=[task_id],
            tasks_index_path=index,
            tasks_dir=task_root,
            builder=_fake_builder,
        )


def test_semantic_hash_ignores_runtime_and_atomic_output_is_valid(tmp_path):
    index, task_root, _, _ = _write_inputs(tmp_path)
    first = evaluate_tasks(
        task_ids=["task-a"],
        tasks_index_path=index,
        tasks_dir=task_root,
        builder=_fake_builder,
    )
    second = evaluate_tasks(
        task_ids=["task-a"],
        tasks_index_path=index,
        tasks_dir=task_root,
        builder=_fake_builder,
    )
    assert first["evaluation_semantic_sha256"] == second[
        "evaluation_semantic_sha256"
    ]
    output = _atomic_write_json(tmp_path / "evaluation.json", first)
    assert json.loads(output.read_text(encoding="utf-8"))[
        "evaluation_semantic_sha256"
    ] == first["evaluation_semantic_sha256"]
