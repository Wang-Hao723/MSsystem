from __future__ import annotations

import copy
import json
import threading
import time
from pathlib import Path

import pytest

from scripts import evaluate_v16_dual_polarity_pairs as evaluator


ROOT = Path(__file__).resolve().parents[1]
MANIFEST_PATH = ROOT / "config" / "v16_dual_polarity_evaluation_pairs.json"

EXPECTED_PAIRS = [
    ("F001_replicate_1_POS_NEG", "c5e1b939eb5a", "4bdd291eda7b"),
    ("F001_replicate_2_POS_NEG", "4acc981f7b0e", "9305a02e24c6"),
    ("F001_replicate_3_POS_NEG", "ade106cc9a74", "9810577e800b"),
    ("F002_replicate_1_POS_NEG", "43c927c80aa6", "0f799e3aa94e"),
    ("F002_replicate_2_POS_NEG", "464a49946f22", "f8e539746744"),
    ("F002_replicate_3_POS_NEG", "5108373b5f64", "8573b9ac79d1"),
    ("F003_replicate_1_POS_NEG", "9e5b1cdc3bf5", "702ab2a35faa"),
    ("F003_replicate_2_POS_NEG", "0013e0b02f5b", "7ad306cb1311"),
    ("F003_replicate_3_POS_NEG", "12df693e554a", "55412720c190"),
    ("F004_replicate_1_POS_NEG", "6b225cc5e3d9", "6b2d440b9363"),
    ("F004_replicate_2_POS_NEG", "55e4b275779c", "87f10abca698"),
    ("F004_replicate_3_POS_NEG", "72f03dae27be", "77d248e6b094"),
    ("F005_replicate_1_POS_NEG", "a4c649288d95", "ffb8d0271c93"),
    ("F005_replicate_2_POS_NEG", "47a14bf56475", "f64a735be1f6"),
    ("F005_replicate_3_POS_NEG", "9fd2ed4184cc", "890f0a541998"),
    ("F006_replicate_1_POS_NEG", "3873259acd5f", "373546c131c4"),
    ("F006_replicate_3_POS_NEG", "cdf120e88e24", "79fabfed7634"),
    ("F007_replicate_1_POS_NEG", "203e426750b0", "550305f42173"),
    ("F007_replicate_2_POS_NEG", "4cc7b7861473", "6ed1c7e608e0"),
    ("F007_replicate_3_POS_NEG", "ab45d501a8a9", "fdd5e34c7726"),
    ("F008_replicate_1_POS_NEG", "22c62a0b9ea8", "556e00e48f52"),
    ("F008_replicate_2_POS_NEG", "9cbc1eb931a2", "3acd92a38c56"),
    ("F008_replicate_3_POS_NEG", "278f067846cc", "3a7c081c8566"),
    ("MIX_33_replicate_1_POS_NEG", "d5f1e22a6285", "6ff25e8cffcb"),
    ("MIX_33_replicate_2_POS_NEG", "9cf967274994", "0107e880d3af"),
]


def _manifest(pairs: list[dict]) -> dict:
    return {
        "schema_version": evaluator.MANIFEST_SCHEMA,
        "manifest_id": "test-explicit-pairs",
        "selection_method": evaluator.CONFIRMATION_SOURCE,
        "task_status_required": "completed",
        "positive_task_ion_mode_required": "POS",
        "negative_task_ion_mode_required": "NEG",
        "candidate_ground_truth_included": False,
        "expected_pair_count": len(pairs),
        "cohort_pair_counts": {"test": len(pairs)},
        "pairs": pairs,
    }


def _pair(label: str, positive: str, negative: str) -> dict:
    return {
        "pair_label": label,
        "positive_task_id": positive,
        "negative_task_id": negative,
        "same_sample_confirmed": True,
        "same_lc_method_confirmed": True,
        "confirmation_source": evaluator.CONFIRMATION_SOURCE,
    }


def _sidecar(pair: dict, *, scale: int = 1) -> dict:
    openms = {}
    for polarity, task_field in (
        ("POS", "positive_task_id"),
        ("NEG", "negative_task_id"),
    ):
        multiplier = scale if polarity == "POS" else scale + 1
        openms[pair[task_field]] = {
            "status": "validated",
            "input_feature_count": 10 * multiplier,
            "input_group_count": 5 * multiplier,
            "eligible_group_count": 2 * multiplier,
            "eligible_feature_count": 4 * multiplier,
            "rejected_group_count": 3 * multiplier,
            "projected_result_row_count": 20 * multiplier,
            "bound_result_row_count": 4 * multiplier,
            "result_binding_mismatch_count": 0,
            "manifest_rejection_reasons": [],
            "saved_result_rows_modified": False,
            "production_runtime_influence": False,
        }
    return {
        "metadata": {
            "sidecar_semantic_sha256": f"{scale:x}" * 64,
            "openms_adduct_group_validation": openms,
            "source_file_sha256_unchanged_after_build": True,
            "task_compound_counts_unchanged": True,
        },
        "summary": {
            "status": "evaluated_with_edges",
            "reliable_anchor_counts": {"POS": 3 * scale, "NEG": 2 * scale},
            "candidate_independent_compatible_edge_count": 2 * scale,
            "candidate_independent_accepted_edge_count": scale,
            "ambiguous_edge_count": scale,
            "ambiguous_anchor_count": 2 * scale,
            "not_evaluable_anchor_count": 0,
            "neutral_entity_count": scale,
            "accepted_neutral_entity_count": scale,
            "precursor_event_evidence_record_count": 2 * scale,
            "dual_polarity_supported_connectivity_count": scale,
            "compound_count_increment": 0,
        },
        "invariants": {
            **{
                field: False for field in evaluator.PRODUCTION_UNCHANGED_FIELDS
            },
            **evaluator.SHADOW_GOVERNANCE_VALUES,
        },
    }


def _health() -> dict:
    return {
        "status": "ok",
        "service_generation": "V16",
        "service_version": "16.0.0",
        "runtime_snapshot_id": "v16-test",
        "runtime_manifest": {
            "status": "fresh",
            "fresh": True,
            "release_ready": True,
            "restart_required": False,
        },
        "capabilities": {
            "v16_dual_polarity_shadow_evidence_graph": True,
            "dual_polarity_shadow_endpoint": evaluator.PAIR_ENDPOINT,
        },
        "storage": {"library_write_enabled": False, "ok": True},
    }


def test_checked_in_manifest_is_the_explicit_23_plus_2_pair_list():
    payload = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))
    pairs = evaluator.validate_manifest(payload)

    assert payload["cohort_pair_counts"] == {"F001_F008": 23, "MIX_33": 2}
    assert [
        (pair["pair_label"], pair["positive_task_id"], pair["negative_task_id"])
        for pair in pairs
    ] == EXPECTED_PAIRS
    assert all(pair["same_sample_confirmed"] is True for pair in pairs)
    assert all(pair["same_lc_method_confirmed"] is True for pair in pairs)
    assert all(
        pair["confirmation_source"] == evaluator.CONFIRMATION_SOURCE
        for pair in pairs
    )
    assert "F006_replicate_2_POS_NEG" not in {pair[0] for pair in EXPECTED_PAIRS}


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("same_sample_confirmed", False),
        ("same_lc_method_confirmed", False),
        ("confirmation_source", "derived_from_filename"),
    ],
)
def test_manifest_rejects_missing_explicit_confirmation(field, value):
    pair = _pair("pair-one", "111111111111", "222222222222")
    pair[field] = value

    with pytest.raises(evaluator.EvaluationError):
        evaluator.validate_manifest(_manifest([pair]))


def test_manifest_rejects_reused_tasks_and_script_has_no_discovery_path():
    pairs = [
        _pair("pair-one", "111111111111", "222222222222"),
        _pair("pair-two", "111111111111", "333333333333"),
    ]
    with pytest.raises(evaluator.EvaluationError, match="reuses a task"):
        evaluator.validate_manifest(_manifest(pairs))

    source = Path(evaluator.__file__).read_text(encoding="utf-8")
    assert "tasks_index" not in source
    assert "glob(" not in source
    assert "list_tasks" not in source


def test_sequential_replay_aggregates_metrics_and_uses_service_hash_verification(
    tmp_path, monkeypatch
):
    pairs = [
        _pair("pair-one", "111111111111", "222222222222"),
        _pair("pair-two", "333333333333", "444444444444"),
    ]
    manifest = _manifest(pairs)
    calls = []
    sidecars = {
        "pair-one": _sidecar(pairs[0], scale=1),
        "pair-two": _sidecar(pairs[1], scale=2),
    }
    pair_ids = {"pair-one": "DPS-1111111111111111", "pair-two": "DPS-2222222222222222"}

    def fake_request(*, method, url, payload, timeout):
        calls.append((method, url, copy.deepcopy(payload), timeout))
        if url.endswith(evaluator.HEALTH_ENDPOINT):
            return _health()
        if method == "POST":
            label = payload["pair_label"]
            assert set(payload) == set(evaluator.POST_FIELDS)
            return {
                "status": "available",
                "pair_id": pair_ids[label],
                "reused": label == "pair-one",
                "artifact": f"/audit/{pair_ids[label]}.json",
                "sidecar": sidecars[label],
            }
        pair_id = url.rsplit("/", 1)[-1]
        label = next(label for label, value in pair_ids.items() if value == pair_id)
        return {
            "status": "available",
            "pair_id": pair_id,
            "sidecar": sidecars[label],
            "artifact": f"/audit/{pair_id}.json",
        }

    manifest_path = tmp_path / "manifest.json"
    output_path = tmp_path / "audit.json"
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    real_replace = evaluator.os.replace
    replacements = []

    def recording_replace(source, destination):
        replacements.append((Path(source), Path(destination)))
        real_replace(source, destination)

    monkeypatch.setattr(evaluator.os, "replace", recording_replace)
    audit = evaluator.run(
        manifest_path=manifest_path,
        output_path=output_path,
        base_url="http://v16.test",
        request_json=fake_request,
    )

    assert [call[0] for call in calls] == ["GET", "POST", "GET", "POST", "GET"]
    assert audit["summary"]["succeeded_pair_count"] == 2
    assert audit["summary"]["failed_pair_count"] == 0
    assert audit["summary"]["reused_pair_count"] == 1
    assert audit["summary"]["evaluation_status_counts"] == {
        "evaluated_with_edges": 2
    }
    aggregate = audit["summary"]["aggregate"]
    assert aggregate["anchors"] == {"POS": 9, "NEG": 6, "total": 15}
    assert aggregate["edges"] == {"compatible": 6, "accepted": 3}
    assert aggregate["ambiguity"]["edge_count"] == 3
    assert aggregate["entities"]["accepted_neutral_entity_count"] == 3
    assert aggregate["dual_candidate_support"]["connectivity_count"] == 3
    assert aggregate["openms_validation"]["POS"]["eligible_group_count"] == 6
    assert aggregate["openms_validation"]["NEG"]["eligible_group_count"] == 10
    assert all(
        row["sidecar_semantic_hash_verification"] == {
            "authority": "v16_service_get_endpoint",
            "verified": True,
            "sidecar_semantic_sha256": row[
                "sidecar_semantic_hash_verification"
            ]["sidecar_semantic_sha256"],
            "locally_recomputed": False,
        }
        for row in audit["pairs"]
    )
    assert audit["execution"]["candidate_ground_truth_loaded"] is False
    assert audit["execution"]["task_mutation_endpoint_called"] is False
    assert replacements and replacements[-1][1] == output_path
    assert replacements[-1][0].parent == output_path.parent
    assert json.loads(output_path.read_text(encoding="utf-8")) == audit
    assert not list(tmp_path.glob("*.tmp"))


def test_pair_failure_is_isolated_and_audit_is_still_written(tmp_path):
    pairs = [
        _pair("pair-one", "111111111111", "222222222222"),
        _pair("pair-two", "333333333333", "444444444444"),
        _pair("pair-three", "555555555555", "666666666666"),
    ]
    manifest = _manifest(pairs)
    sidecars = {
        "pair-one": _sidecar(pairs[0], scale=1),
        "pair-three": _sidecar(pairs[2], scale=3),
    }
    pair_ids = {
        "pair-one": "DPS-1111111111111111",
        "pair-three": "DPS-3333333333333333",
    }
    post_labels = []

    def fake_request(*, method, url, payload, timeout):
        if url.endswith(evaluator.HEALTH_ENDPOINT):
            return _health()
        if method == "POST":
            label = payload["pair_label"]
            post_labels.append(label)
            if label == "pair-two":
                raise evaluator.HTTPRequestError("synthetic 409", status_code=409)
            return {
                "status": "available",
                "pair_id": pair_ids[label],
                "reused": False,
                "artifact": f"/audit/{pair_ids[label]}.json",
                "sidecar": sidecars[label],
            }
        pair_id = url.rsplit("/", 1)[-1]
        label = next(label for label, value in pair_ids.items() if value == pair_id)
        return {"status": "available", "pair_id": pair_id, "sidecar": sidecars[label]}

    manifest_path = tmp_path / "manifest.json"
    output_path = tmp_path / "audit.json"
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    audit = evaluator.run(
        manifest_path=manifest_path,
        output_path=output_path,
        base_url="http://v16.test",
        request_json=fake_request,
    )

    assert post_labels == ["pair-one", "pair-two", "pair-three"]
    assert [row["status"] for row in audit["pairs"]] == [
        "succeeded",
        "failed",
        "succeeded",
    ]
    assert audit["pairs"][1]["failure"]["http_status_code"] == 409
    assert audit["summary"]["succeeded_pair_count"] == 2
    assert audit["summary"]["failed_pair_count"] == 1
    assert audit["summary"]["failed_pairs"][0]["pair_label"] == "pair-two"
    assert audit["summary"]["all_manifest_pairs_verified"] is False
    assert audit["v15_production_unchanged_invariants"][
        "all_successful_pairs_verified"
    ] is True
    assert output_path.is_file()


def test_optional_concurrency_two_is_bounded_and_preserves_manifest_order():
    pairs = [
        _pair("pair-one", "111111111111", "222222222222"),
        _pair("pair-two", "333333333333", "444444444444"),
        _pair("pair-three", "555555555555", "666666666666"),
    ]
    sidecars = {
        pair["pair_label"]: _sidecar(pair, scale=index + 1)
        for index, pair in enumerate(pairs)
    }
    pair_ids = {
        pair["pair_label"]: f"DPS-{index + 1:016x}"
        for index, pair in enumerate(pairs)
    }
    state = {"active": 0, "maximum": 0}
    lock = threading.Lock()

    def fake_request(*, method, url, payload, timeout):
        if url.endswith(evaluator.HEALTH_ENDPOINT):
            return _health()
        if method == "POST":
            label = payload["pair_label"]
            with lock:
                state["active"] += 1
                state["maximum"] = max(state["maximum"], state["active"])
            try:
                time.sleep(0.01)
            finally:
                with lock:
                    state["active"] -= 1
            return {
                "status": "available",
                "pair_id": pair_ids[label],
                "reused": False,
                "sidecar": sidecars[label],
            }
        pair_id = url.rsplit("/", 1)[-1]
        label = next(label for label, value in pair_ids.items() if value == pair_id)
        return {"status": "available", "pair_id": pair_id, "sidecar": sidecars[label]}

    audit = evaluator.evaluate_pairs(
        _manifest(pairs),
        base_url="http://v16.test",
        concurrency=2,
        request_json=fake_request,
    )

    assert state["maximum"] == 2
    assert [row["pair_label"] for row in audit["pairs"]] == [
        "pair-one",
        "pair-two",
        "pair-three",
    ]


@pytest.mark.parametrize("concurrency", [0, 3, 10])
def test_concurrency_is_strictly_limited_to_one_or_two(concurrency):
    manifest = _manifest([_pair("pair-one", "111111111111", "222222222222")])
    with pytest.raises(evaluator.EvaluationError, match="concurrency"):
        evaluator.evaluate_pairs(
            manifest,
            base_url="http://v16.test",
            concurrency=concurrency,
            request_json=lambda **_kwargs: _health(),
        )


def test_wrong_service_generation_fails_before_any_post():
    manifest = _manifest([_pair("pair-one", "111111111111", "222222222222")])
    methods = []

    def fake_request(*, method, url, payload, timeout):
        methods.append(method)
        health = _health()
        health["service_generation"] = "V15"
        return health

    with pytest.raises(evaluator.EvaluationError, match="non-V16"):
        evaluator.evaluate_pairs(
            manifest,
            base_url="http://wrong.test",
            request_json=fake_request,
        )
    assert methods == ["GET"]


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        ({"status": "degraded"}, "degraded"),
        ({"runtime_manifest.fresh": False}, "stale"),
        ({"runtime_manifest.release_ready": False}, "stale"),
        ({"runtime_manifest.restart_required": True}, "stale"),
        ({"storage.ok": False}, "storage"),
        ({"storage.library_write_enabled": True}, "read-only"),
    ],
)
def test_unhealthy_or_mutable_runtime_fails_before_any_post(mutation, message):
    manifest = _manifest([_pair("pair-one", "111111111111", "222222222222")])
    methods = []

    def fake_request(*, method, url, payload, timeout):
        methods.append(method)
        health = _health()
        for dotted, value in mutation.items():
            target = health
            parts = dotted.split(".")
            for part in parts[:-1]:
                target = target[part]
            target[parts[-1]] = value
        return health

    with pytest.raises(evaluator.EvaluationError, match=message):
        evaluator.evaluate_pairs(
            manifest,
            base_url="http://unhealthy.test",
            request_json=fake_request,
        )
    assert methods == ["GET"]
