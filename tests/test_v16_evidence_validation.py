from __future__ import annotations

import copy
import json
import random
from pathlib import Path

import pytest

from eval.v16_evidence_validation.evaluate import (
    PROTOCOL,
    base_receipt,
    bind_retrospective_labels,
    derangement,
    edge_lists,
    evaluate_controls,
    explicit_groups,
    file_sha256,
    load_controls_receipt,
    read_bound_json,
    receipt_semantic_sha256,
    score_recovery,
    validate_graph,
    write_new_json,
)
from matching.v16_dual_polarity import (
    load_v16_dual_polarity_policy,
    semantic_sha256,
)


@pytest.fixture
def policy():
    return load_v16_dual_polarity_policy()


@pytest.fixture
def protocol():
    result = json.loads(PROTOCOL.read_text())
    result["controls"]["rt_derangement_rounds"] = 6
    return result


def anchor(task: str, polarity: str, number: int, mass: float, rt: float) -> dict:
    return {
        "anchor_id": f"{task}-{number}", "task_id": task, "polarity": polarity,
        "neutral_mass": mass, "retention_time_min": rt,
        "aligned_retention_time_min": rt, "candidate_independent": True,
        "source_row_indices": [number],
    }


def pair(label: str, group: str, policy: dict, *, single_negative=False) -> dict:
    tasks = {p: f"{label}-{p}" for p in ("POS", "NEG")}
    anchors = {
        p: [anchor(tasks[p], p, i, mass, rt)
            for i, (mass, rt) in enumerate(((300.0, 3.0), (600.0, 8.0)))]
        for p in ("POS", "NEG")
    }
    if single_negative:
        anchors["NEG"] = anchors["NEG"][:1]
    graph = {
        "task_pair": tasks, "policy_sha256": semantic_sha256(policy),
        "anchors": anchors, **edge_lists(anchors["POS"], anchors["NEG"], policy),
        "candidate_fields_present": False, "spectrum_similarity_present": False,
    }
    graph["semantic_sha256"] = semantic_sha256(graph)
    return {"pair_label": label, "study_group": group, "task_pair": tasks,
            "graph": graph, "rt_alignment": {"status": "not_supplied_raw_rt_used"}}


def test_frozen_graph_is_reproduced_by_the_same_edge_algorithm(policy):
    fixture = pair("one", "A", policy)
    validate_graph(fixture["graph"], fixture["task_pair"], policy)
    assert len(fixture["graph"]["accepted_unique_edges"]) == 2


def test_tampered_graph_is_rejected_before_controls(policy):
    fixture = pair("one", "A", policy)
    fixture["graph"]["anchors"]["NEG"][0]["neutral_mass"] += 1
    with pytest.raises(ValueError, match="hash mismatch"):
        validate_graph(fixture["graph"], fixture["task_pair"], policy)


def test_rehashed_inconsistent_graph_cannot_pass_edge_reproduction(policy):
    fixture = pair("one", "A", policy)
    fixture["graph"]["accepted_unique_edges"] = []
    fixture["graph"]["semantic_sha256"] = semantic_sha256(
        {k: v for k, v in fixture["graph"].items() if k != "semantic_sha256"})
    with pytest.raises(ValueError, match="edge reproduction"):
        validate_graph(fixture["graph"], fixture["task_pair"], policy)


def test_controls_preserve_input_and_do_not_claim_fdr_or_identity_accuracy(policy, protocol):
    fixtures = [pair("one", "A", policy), pair("two", "B", policy)]
    saved = copy.deepcopy(fixtures)
    result = evaluate_controls(fixtures, policy, protocol)
    assert fixtures == saved
    assert result["observed"]["summary"]["accepted"] == 4
    assert result["synthetic_zero_controls_passed"] is True
    assert result["confirmed_false_positive_rate"] is None
    assert result["empirical_fdr"] is None
    assert result["true_identity_pair_recall"] is None
    assert result["production_promotion_allowed"] is False


def test_cross_sample_controls_exclude_all_same_study_replicates(policy, protocol):
    fixtures = [pair("one", "A", policy), pair("two", "A", policy),
                pair("three", "B", policy)]
    result = evaluate_controls(fixtures, policy, protocol)["cross_sample"]
    assert result["summary"]["pair_count"] == 4
    assert result["matches_are_not_verified_false_positives"] is True
    forbidden = {("one", "two"), ("two", "one")}
    assert not any((x["positive_pair_label"], x["negative_pair_label"]) in forbidden
                   for x in result["pairs"])


def test_rt_controls_are_repeatable_and_singletons_are_explicitly_excluded(policy, protocol):
    fixtures = [pair("one", "A", policy), pair("singleton", "B", policy, single_negative=True)]
    result = evaluate_controls(fixtures, policy, protocol)
    assert result == evaluate_controls(fixtures, policy, protocol)
    assert result["rt_derangement"]["permutable_pair_count"] == 1
    assert result["rt_derangement"]["excluded_pair_labels"] == ["singleton"]
    assert result["rt_derangement"]["accepted_edge_count_max"] == 0
    assert result["rt_alignment_status_counts"] == {"not_supplied_raw_rt_used": 2}


def test_missing_anchors_remain_not_evaluable_in_control_denominators(policy, protocol):
    fixture = pair("missing", "A", policy)
    fixture["graph"]["anchors"]["NEG"] = []
    result = evaluate_controls([fixture], policy, protocol)
    assert result["observed"]["summary"]["pair_count"] == 1
    assert result["observed"]["summary"]["evaluable_pair_count"] == 0
    assert result["rt_derangement"]["permutable_pair_count"] == 0


@pytest.mark.parametrize("size", [2, 3, 10])
def test_derangement_has_no_fixed_positions_and_preserves_marginals(size):
    shuffled = derangement(size, random.Random(42))
    assert sorted(shuffled) == list(range(size))
    assert all(i != j for i, j in enumerate(shuffled))
    with pytest.raises(ValueError):
        derangement(1, random.Random(42))


def test_output_is_atomic_and_never_overwrites_existing_audit(tmp_path):
    path = tmp_path / "audit.json"
    write_new_json(path, {"kept": True})
    before = path.read_bytes()
    with pytest.raises(FileExistsError):
        write_new_json(path, {"replacement": True})
    assert path.read_bytes() == before
    assert list(tmp_path.iterdir()) == [path]


def test_input_hash_binding_and_duplicate_group_rejection(tmp_path):
    path = tmp_path / "input.json"
    path.write_text('{"n": 1}')
    assert read_bound_json(path, file_sha256(path)) == {"n": 1}
    with pytest.raises(ValueError, match="Frozen input hash mismatch"):
        read_bound_json(path, "0" * 64)
    with pytest.raises(ValueError, match="Duplicate explicit pair"):
        explicit_groups({"explicit_study_groups": {"A": ["x"], "B": ["x"]}})


def test_score_requires_integrity_bound_controls_receipt(tmp_path, protocol):
    payload = {**base_receipt(protocol, "controls"), "results": {}}
    payload["semantic_sha256"] = receipt_semantic_sha256(payload)
    path = tmp_path / "controls.json"
    path.write_text(json.dumps(payload))
    assert load_controls_receipt(path, protocol)["phase"] == "controls"
    payload["results"]["pretend_gain"] = 100
    path.write_text(json.dumps(payload))
    with pytest.raises(ValueError, match="not frozen"):
        load_controls_receipt(path, protocol)


def label(task_id, *, duplicate=False):
    return {"record_id": f"{task_id}:standard", "truth_connectivity": "A" * 14,
            "same_injection_duplicate_selected_entity": duplicate}


def test_recovery_counts_frozen_anchor_links_without_calling_them_true_identity(policy):
    fixture = pair("one", "A", policy)
    labels, bindings = {}, {}
    for polarity, task_id in fixture["task_pair"].items():
        row = label(task_id)
        labels[task_id] = [row]
        bindings[(task_id, row["record_id"])] = {"status": "verified_measured_query",
            "anchor_ids": [fixture["graph"]["anchors"][polarity][0]["anchor_id"]]}
    result = score_recovery([fixture], bindings, labels)
    assert result["accepted_mass_selected_anchor_pair_count"] == 1
    assert result["study_group_by_connectivity_with_any_accepted_link"] == 1
    assert result["true_identity_pair_recall"] is None
    assert result["top1_accuracy_gain"] is None


def test_failed_query_binding_never_scores_an_accepted_standard_link(policy):
    fixture = pair("one", "A", policy)
    labels, bindings = {}, {}
    for polarity, task_id in fixture["task_pair"].items():
        row = label(task_id)
        labels[task_id] = [row]
        bindings[(task_id, row["record_id"])] = {"status": "measured_query_hash_mismatch",
            "anchor_ids": [fixture["graph"]["anchors"][polarity][0]["anchor_id"]]}
    result = score_recovery([fixture], bindings, labels)
    assert result["accepted_mass_selected_anchor_pair_count"] == 0
    assert result["pair_status_counts"] == {"unverified_source_or_measured_query_binding": 1}


def test_ambiguous_reference_selection_is_excluded_not_credited(policy):
    fixture = pair("one", "A", policy)
    labels = {task: [label(task, duplicate=True)] for task in fixture["task_pair"].values()}
    result = score_recovery([fixture], {}, labels)
    assert result["pair_status_counts"] == {"excluded_ambiguous_reference_selection": 1}
    assert result["accepted_mass_selected_anchor_pair_count"] == 0


@pytest.mark.parametrize("query_hash_matches", [True, False])
def test_changed_flagship_requires_exact_measured_query_hash(tmp_path, policy, query_hash_matches):
    from eval.v15_top40_restricted_rerank.evaluate_phase1_mix33_transfer import (
        _query_projection, _stream_rows,
    )

    fixture = pair("one", "A", policy)
    task_id = fixture["task_pair"]["POS"]
    task_root = tmp_path / "output/tasks"
    task_root.mkdir(parents=True)
    raw_path, flagship_path = task_root / "raw.json", task_root / "flagship.json"
    raw = {"spectrum_id": "scan-1", "precursor_mz": 301.0, "ion_mode": "POS",
           "raw_peaks": [[100.0, 1.0]], "raw_peaks_complete": True}
    entity = {"entity_id": "entity-1", "member_result_indices": [0],
              "member_spectrum_ids": ["scan-1"]}
    raw_path.write_text(json.dumps([raw]))
    flagship_path.write_text(json.dumps({"entities": [entity]}))
    # The original label producer uses ijson Decimal values in peak arrays;
    # reproduce that byte-level contract rather than normalizing it away.
    persisted_rows = _stream_rows(raw_path, {0})
    projection, _ = _query_projection("test", task_id, 0, entity, [persisted_rows[0]])
    row = {**label(task_id), "run": "run-1", "mode": "POS",
           "selected_entity_available": True, "selected_entity_id": "entity-1",
           "query_spectrum_sha256": projection["query_spectrum_sha256"] if query_hash_matches else "0" * 64}
    labels_path, sources_path = tmp_path / "labels.json", tmp_path / "sources.json"
    labels_path.write_text(json.dumps({"records": [row]}))
    sources_path.write_text(json.dumps({"source_files": {"run-1": {
        "task_id": task_id, "mode": "POS", "raw_results": {"sha256": file_sha256(raw_path)},
        "v15_artifact": {"sha256": "1" * 64}}}}))
    fixture["source_file_bindings"] = {task_id: {
        "results_path": str(raw_path), "results_file_sha256": file_sha256(raw_path),
        "v15_flagship_path": str(flagship_path), "v15_flagship_file_sha256": file_sha256(flagship_path)}}
    protocol = {"inputs": {
        "retrospective_labels": {"path": "labels.json", "sha256": file_sha256(labels_path)},
        "retrospective_source_manifest": {"path": "sources.json", "sha256": file_sha256(sources_path)}}}
    bindings, _, audit = bind_retrospective_labels(tmp_path, [fixture], protocol)
    bound = bindings[(task_id, row["record_id"])]
    assert audit["tasks"][task_id]["flagship_version_or_content_changed"] is True
    if query_hash_matches:
        assert bound["status"] == "verified_measured_query"
        assert bound["anchor_ids"] == [f"{task_id}-0"]
    else:
        assert bound["status"] == "measured_query_hash_mismatch"
        assert bound["anchor_ids"] == []
