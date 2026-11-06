from __future__ import annotations

import copy
import json

import pytest

from eval.v16_evidence_validation import centroid_challenger as challenger
from eval.v16_evidence_validation import evaluate as audit
from matching.v16_dual_polarity import (
    _classify_edges, _frozen_candidate_independent_graph,
    load_v16_dual_polarity_policy, semantic_sha256,
)
from matching.v16_iimn_feature_centroid_shadow import (
    build_v16_iimn_feature_centroid_shadow,
    load_v16_iimn_feature_centroid_shadow_policy,
)


def anchor(polarity, number, mass, rt=2.0):
    return {
        "anchor_id": f"{polarity}-{number}", "task_id": polarity,
        "polarity": polarity, "source": "reliable_iimn", "group_id": str(number),
        "neutral_mass": mass, "retention_time_min": rt,
        "aligned_retention_time_min": rt, "rt_basis": "raw_task_rt",
        "precursor_mz": mass, "feature_ids": [f"{polarity}-{number}-a", f"{polarity}-{number}-b"],
        "source_row_indices": [2 * number, 2 * number + 1],
        "adducts": ["[M-H]-" if polarity == "NEG" else "[M+H]+"],
        "verified_member_feature_count": 2, "candidate_independent": True,
    }


@pytest.fixture
def policy():
    return load_v16_dual_polarity_policy()


def paired(policy):
    anchors = {p: [anchor(p, 0, 300)] for p in ("POS", "NEG")}
    anchors["POS"].append(anchor("POS", 1, 500))
    tasks = {p: p for p in anchors}
    graph = _frozen_candidate_independent_graph(anchors,
        *_classify_edges(anchors["POS"], anchors["NEG"], policy),
        task_pair=tasks, policy_sha256=semantic_sha256(policy))
    graph["semantic_sha256"] = semantic_sha256(graph)
    return {"pair_label": "test_pair", "study_group": "TEST", "task_pair": tasks,
        "graph": graph, "rt_alignment": {"status": "not_supplied_raw_rt_used"}}


def test_disjoint_centroid_extension_preserves_all_anchors_and_adds_a_pair(policy):
    pair = paired(policy)
    frozen = copy.deepcopy(pair)
    updated, changes = challenger.extend_pair(pair, [anchor("NEG", 1, 500)], policy)
    assert pair == frozen
    assert updated["graph"]["anchors"]["POS"] == pair["graph"]["anchors"]["POS"]
    assert pair["graph"]["anchors"]["NEG"][0] in updated["graph"]["anchors"]["NEG"]
    assert len(changes["added_anchor_ids"]) == 1
    assert len(changes["added_accepted_edges"]) == 1
    assert not changes["lost_accepted_edges"]


@pytest.mark.parametrize("field", ["source_row_indices", "feature_ids"])
def test_overlapping_centroid_groups_are_quarantined_without_overwriting(policy, field):
    pair = paired(policy)
    extra = anchor("NEG", 1, 500)
    extra[field] = pair["graph"]["anchors"]["NEG"][0][field]
    updated, changes = challenger.extend_pair(pair, [extra], policy)
    assert len(changes["quarantined_overlap"]) == 1
    assert not changes["added_anchor_ids"]
    assert updated["graph"] == pair["graph"]


def test_new_ambiguity_is_not_hidden_to_keep_baseline_acceptance(policy):
    pair = paired(policy)
    updated, changes = challenger.extend_pair(pair, [anchor("NEG", 1, 300)], policy)
    assert len(changes["lost_accepted_edges"]) == 1
    assert not updated["graph"]["accepted_unique_edges"]
    assert len(updated["graph"]["ambiguous_no_merge_edges"]) == 2
    assert changes["all_existing_negative_anchors_unchanged"] is True


@pytest.mark.parametrize("mutation", [
    {"polarity": "POS"}, {"task_id": "wrong-task"},
    {"candidate_independent": False}, {"source": "unlinked_default_adduct"},
])
def test_unbound_or_candidate_dependent_anchors_are_rejected(policy, mutation):
    extra = {**anchor("NEG", 1, 500), **mutation}
    with pytest.raises(ValueError, match="bound candidate-independent NEG"):
        challenger.extend_pair(paired(policy), [extra], policy)


def test_extension_is_order_independent_and_does_not_export_candidate_names(policy):
    extras = [anchor("NEG", 1, 500), anchor("NEG", 2, 700)]
    extras[0]["candidate_name"] = "must not reach graph"
    first, _ = challenger.extend_pair(paired(policy), extras, policy)
    second, _ = challenger.extend_pair(paired(policy), list(reversed(extras)), policy)
    assert first == second
    assert "candidate_name" not in json.dumps(first)


def test_mixed_rt_alignment_is_rejected(policy):
    pair = paired(policy)
    pair["rt_alignment"]["status"] = "applied"
    with pytest.raises(ValueError, match="aligned and unaligned"):
        challenger.extend_pair(pair, [], policy)


@pytest.mark.parametrize("label,group,has_anchor,allowed", [
    ("MIX_33_replicate_1_POS_NEG", "MIX_33", False, True),
    ("MIX_33_replicate_2_POS_NEG", "MIX_33", False, True),
    ("MIX_33_replicate_1_POS_NEG", "MIX_33", True, False),
    ("F001_replicate_1_POS_NEG", "F001", False, False),
    ("MIX_33_replicate_3_POS_NEG", "MIX_33", False, False),
])
def test_missing_centroid_source_is_allowed_only_for_two_frozen_empty_mix_pairs(
        label, group, has_anchor, allowed):
    protocol = json.loads(challenger.PROTOCOL.read_text())
    pair = {"pair_label": label, "study_group": group,
            "graph": {"anchors": {"NEG": ["anchor"] if has_anchor else []}}}
    assert challenger.retain_without_centroid_source(pair, protocol) is allowed


def reconstruction_fixture(tmp_path, monkeypatch):
    rows = [
        {"spectrum_id": "s1", "precursor_feature_id": "PF1", "lcms_feature_id": "LF1",
         "lcms_feature_link_status": "linked_openms_feature", "lcms_feature_candidate_count": 1,
         "precursor_mz": 298.992723533379, "lcms_feature_mz": 298.9969,
         "retention_time": 1.0, "retention_time_unit": "minutes", "ion_mode": "NEG",
         "match": {"compound_name": "not permitted in reconstruction", "final_score": .99}},
        {"spectrum_id": "s2", "precursor_feature_id": "PF2", "lcms_feature_id": "LF2",
         "lcms_feature_link_status": "linked_openms_feature", "lcms_feature_candidate_count": 1,
         "precursor_mz": 344.998201, "lcms_feature_mz": 345.0025,
         "retention_time": 1.02, "retention_time_unit": "minutes", "ion_mode": "NEG"},
    ]
    traces = {feature: {"rt": [.97, .98, .99, 1, 1.01, 1.02, 1.03],
                         "intensity": [0, 1, 3, 10, 3, 1, 0]}
              for feature in ("LF1", "LF2")}
    result_path, source_path = tmp_path / "results.json", tmp_path / "source.mzML"
    result_path.write_text(json.dumps(rows), encoding="utf-8")
    source_path.write_bytes(b"test placeholder: traces injected")
    previous = build_v16_iimn_feature_centroid_shadow(
        {"id": "NEG", "status": "completed", "ion_mode": "NEG"}, rows,
        source_results_path=result_path, source_data_path=source_path,
        current_traces=traces, feature_centroid_traces=traces)
    binding = {**previous["source_binding"],
        "candidate_independent_projection_sha256": previous["metadata"]["candidate_independent_projection_sha256"]}
    record = {"task_id": "NEG", "polarity": "NEG", "source_binding": binding,
        "candidate_independent_result_rows": 2,
        "linked_feature_centroid_center": previous["comparison"]["feature_centroid_center"]}
    pair = {"task_pair": {"POS": "POS", "NEG": "NEG"}, "source_file_bindings": {
        "NEG": {"results_path": str(result_path),
                "results_file_sha256": audit.file_sha256(result_path)}}}
    monkeypatch.setattr(challenger, "extract_targeted_ms1_traces",
                        lambda *args, **kwargs: (copy.deepcopy(traces), {"status": "ok"}))
    return record, pair


def test_reconstruction_requires_exact_prior_trace_network_and_source_bindings(tmp_path, monkeypatch, policy):
    record, pair = reconstruction_fixture(tmp_path, monkeypatch)
    snapshots = {}
    anchors, result = challenger.reconstruct_negative(record, pair, policy,
        load_v16_iimn_feature_centroid_shadow_policy(), snapshots)
    assert len(anchors) == 1
    assert anchors[0]["source_row_indices"] == [0, 1]
    assert result["reproduces_prior_trace_and_network"] is True
    assert len(snapshots) == 2
    assert "compound_name" not in json.dumps(anchors)


@pytest.mark.parametrize("field", ["feature_centroid_trace_semantic_sha256", "network_semantic_sha256"])
def test_reconstruction_fails_when_prior_semantics_do_not_reproduce(tmp_path, monkeypatch, policy, field):
    record, pair = reconstruction_fixture(tmp_path, monkeypatch)
    if field.startswith("feature_centroid"):
        record["source_binding"][field] = "incorrect"
    else:
        record["linked_feature_centroid_center"][field] = "incorrect"
    with pytest.raises(ValueError, match="differs from prior evaluation"):
        challenger.reconstruct_negative(record, pair, policy,
            load_v16_iimn_feature_centroid_shadow_policy(), {})


def test_rebinding_compares_graphs_using_the_same_measured_rows(policy):
    pair = paired(policy)
    bindings = {("NEG", "label"): {"status": "verified_measured_query",
        "member_result_indices": [0], "anchor_ids": ["wrong-old-id"]}}
    rebound = challenger.rebind_anchors(bindings, [pair])
    assert rebound[("NEG", "label")]["anchor_ids"] == ["NEG-0"]
    assert bindings[("NEG", "label")]["anchor_ids"] == ["wrong-old-id"]


def test_bad_controls_are_rejected_before_any_labels_are_read(tmp_path, monkeypatch):
    path = tmp_path / "controls.json"
    path.write_text(json.dumps({"phase": "score"}), encoding="utf-8")
    monkeypatch.setattr(audit, "bind_retrospective_labels",
        lambda *args: pytest.fail("labels opened before validating receipt"))
    with pytest.raises(ValueError, match="receipt or evaluator binding"):
        challenger.score({"protocol_id": "test"}, {}, path)


def test_empty_promotion_metric_is_not_manufactured_from_recovery():
    row = {"pair_label": "pair", "truth_connectivity": "structure",
        "both_measured_queries_bound": True, "positive_anchor_ids": ["p"],
        "negative_anchor_ids": [], "status": "missing_reliable_graph_anchor"}
    before = {"details": [row], "accepted_mass_selected_anchor_pair_count": 0,
              "both_measured_queries_bound_count": 1}
    after = {**before, "details": [{**row, "negative_anchor_ids": ["n"],
        "status": "accepted_mass_selected_anchor_link"}], "accepted_mass_selected_anchor_pair_count": 1}
    comparison = challenger.compare_recovery(before, after)
    assert len(comparison["gained_accepted_pairs"]) == 1
    assert comparison["negative_anchor_missing_before"] == 1
    assert comparison["negative_anchor_missing_after"] == 0
    assert comparison["top1_accuracy_gain"] is None
    assert comparison["independent_validation"] is False
