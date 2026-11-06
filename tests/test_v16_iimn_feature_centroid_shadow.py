from __future__ import annotations

import copy
import json

import pytest

from matching.v16_iimn_feature_centroid_shadow import (
    ALGORITHM_ID,
    MODULE_ID,
    build_v16_iimn_feature_centroid_shadow,
    load_v16_iimn_feature_centroid_shadow_policy,
    project_candidate_independent_rows,
)


def _task() -> dict:
    return {
        "id": "neg-task",
        "status": "completed",
        "ion_mode": "NEG",
        "result_run_id": "run-1",
        "untargeted_workflow": {"feature_xml": "/read-only/neg.featureXML"},
    }


def _rows() -> list[dict]:
    return [
        {
            "spectrum_id": "scan-1",
            "precursor_feature_id": "PF1",
            "lcms_feature_id": "LF1",
            "lcms_feature_link_status": "linked_openms_feature",
            "lcms_feature_candidate_count": 1,
            "precursor_mz": 298.992723533379,
            "lcms_feature_mz": 298.9969,
            "retention_time": 1.00,
            "lcms_feature_rt_minutes": 1.01,
            "retention_time_unit": "minutes",
            "ion_mode": "NEG",
            "match": {
                "compound_name": "must not be projected",
                "formula": "C10H20O10",
                "inchikey": "AAAAAAAAAAAAAA-BBBBBBBBBB-C",
                "final_score": 0.99,
            },
        },
        {
            "spectrum_id": "scan-2",
            "precursor_feature_id": "PF2",
            "lcms_feature_id": "LF2",
            "lcms_feature_link_status": "linked_openms_feature",
            "lcms_feature_candidate_count": 2,
            "precursor_mz": 344.998201,
            "lcms_feature_mz": 345.0025,
            "retention_time": 1.02,
            "lcms_feature_rt_minutes": 1.01,
            "retention_time_unit": "minutes",
            "ion_mode": "NEG",
            "candidates": [{"name": "must not be projected"}],
        },
    ]


def _trace(left: list[float], right: list[float]):
    rt = [0.97, 0.98, 0.99, 1.00, 1.01, 1.02, 1.03]
    return {
        "LF1": {"rt": rt, "intensity": left, "target_mz": 298.99},
        "LF2": {"rt": rt, "intensity": right, "target_mz": 345.00},
    }


def test_policy_is_locked_shadow_only():
    policy = load_v16_iimn_feature_centroid_shadow_policy()
    assert policy["module_id"] == MODULE_ID
    assert policy["algorithm_id"] == ALGORITHM_ID
    assert policy["execution"]["mode"] == "shadow_only"
    assert policy["execution"]["production_results_mutated"] is False

    broken = copy.deepcopy(policy)
    broken["execution"]["production_ranking_changed"] = True
    with pytest.raises(ValueError, match="production_ranking_changed"):
        build_v16_iimn_feature_centroid_shadow(
            _task(),
            _rows(),
            current_traces=_trace([0, 1, 3, 10, 3, 1, 0], [0, 1, 3, 10, 3, 1, 0]),
            feature_centroid_traces=_trace(
                [0, 1, 3, 10, 3, 1, 0], [0, 1, 3, 10, 3, 1, 0]
            ),
            policy=broken,
        )


def test_projection_excludes_candidates_and_is_candidate_invariant():
    first = _rows()
    second = copy.deepcopy(first)
    second[0]["match"] = {"compound_name": "different", "final_score": 0.01}
    second[1]["candidates"] = [{"name": "different"}]
    assert project_candidate_independent_rows(first) == (
        project_candidate_independent_rows(second)
    )
    projection_text = json.dumps(
        project_candidate_independent_rows(first), sort_keys=True
    )
    assert "compound_name" not in projection_text
    assert "final_score" not in projection_text
    assert "candidates" not in projection_text


def test_feature_centroid_trace_can_add_network_without_changing_pair_inputs():
    current = _trace(
        [0, 1, 3, 10, 3, 1, 0],
        [0, 8, 1, 9, 1, 8, 0],
    )
    shadow = _trace(
        [0, 1, 3, 10, 3, 1, 0],
        [0, 2, 6, 20, 6, 2, 0],
    )
    original = _rows()
    frozen = copy.deepcopy(original)
    payload = build_v16_iimn_feature_centroid_shadow(
        _task(),
        original,
        current_traces=current,
        feature_centroid_traces=shadow,
    )
    assert original == frozen
    current_graph = payload["comparison"]["current_precursor_center"]
    shadow_graph = payload["comparison"]["feature_centroid_center"]
    assert current_graph["summary"]["network_count"] == 0
    assert shadow_graph["summary"]["network_count"] == 1
    assert payload["comparison"]["delta"]["network_count"] == 1
    assert payload["invariants"]["neutral_mass_hypothesis_input_unchanged"] is True
    assert payload["invariants"]["pairing_rt_input_unchanged"] is True
    assert payload["invariants"]["production_ranking_changed"] is False


def test_row_audit_uses_only_unique_linked_feature_and_records_fallbacks():
    rows = _rows() + [
        {
            "spectrum_id": "scan-3",
            "precursor_feature_id": "PF3",
            "lcms_feature_id": "LF3",
            "lcms_feature_link_status": "ambiguous_best_feature",
            "precursor_mz": 400.0,
            "lcms_feature_mz": 400.01,
            "retention_time": 2.0,
            "ion_mode": "NEG",
        },
        {
            "spectrum_id": "scan-4",
            "precursor_feature_id": "PF4",
            "lcms_feature_link_status": "linked_openms_feature",
            "precursor_mz": 500.0,
            "lcms_feature_mz": 500.01,
            "retention_time": 3.0,
            "ion_mode": "NEG",
        },
    ]
    trace = _trace(
        [0, 1, 3, 10, 3, 1, 0],
        [0, 2, 6, 20, 6, 2, 0],
    )
    payload = build_v16_iimn_feature_centroid_shadow(
        _task(), rows, current_traces=trace, feature_centroid_traces=trace
    )
    audit = payload["trace_center_audit"]
    assert audit["counts"]["rows_using_feature_centroid"] == 2
    assert audit["counts"]["rows_using_precursor_fallback"] == 2
    assert audit["rows"][0]["shadow_trace_center_source"] == "lcms_feature_mz"
    assert audit["rows"][2]["fallback_reason"] == (
        "feature_binding_missing_or_nonunique"
    )
    assert audit["rows"][3]["fallback_reason"] == "linked_feature_id_missing"


def test_candidate_changes_do_not_change_network_semantics():
    traces = _trace(
        [0, 1, 3, 10, 3, 1, 0],
        [0, 2, 6, 20, 6, 2, 0],
    )
    first_rows = _rows()
    second_rows = copy.deepcopy(first_rows)
    second_rows[0]["match"] = {
        "compound_name": "unrelated",
        "formula": "C1H1",
        "inchikey": "ZZZZZZZZZZZZZZ-YYYYYYYYYY-X",
        "final_score": 1.0,
    }
    first = build_v16_iimn_feature_centroid_shadow(
        _task(), first_rows, current_traces=traces, feature_centroid_traces=traces
    )
    second = build_v16_iimn_feature_centroid_shadow(
        _task(), second_rows, current_traces=traces, feature_centroid_traces=traces
    )
    assert first["metadata"]["candidate_independent_projection_sha256"] == (
        second["metadata"]["candidate_independent_projection_sha256"]
    )
    assert first["comparison"]["feature_centroid_center"][
        "network_semantic_sha256"
    ] == second["comparison"]["feature_centroid_center"][
        "network_semantic_sha256"
    ]


def test_source_hashes_are_bound_and_files_are_not_modified(tmp_path):
    results_path = tmp_path / "results.json"
    raw_path = tmp_path / "source.mzML"
    results_path.write_text(json.dumps(_rows()), encoding="utf-8")
    raw_path.write_bytes(b"read-only-source")
    result_before = results_path.read_bytes()
    raw_before = raw_path.read_bytes()
    traces = _trace(
        [0, 1, 3, 10, 3, 1, 0],
        [0, 2, 6, 20, 6, 2, 0],
    )
    payload = build_v16_iimn_feature_centroid_shadow(
        _task(),
        _rows(),
        source_results_path=results_path,
        source_data_path=raw_path,
        current_traces=traces,
        feature_centroid_traces=traces,
    )
    assert results_path.read_bytes() == result_before
    assert raw_path.read_bytes() == raw_before
    assert payload["source_binding"]["results_file"]["sha256"]
    assert payload["source_binding"]["source_data_file"]["sha256"]
    assert payload["source_binding"]["source_files_unchanged"] is True


@pytest.mark.parametrize(
    "task, error",
    [
        ({"id": "x", "status": "running", "ion_mode": "NEG"}, "completed"),
        ({"id": "x", "status": "completed", "ion_mode": "AUTO"}, "POS or NEG"),
    ],
)
def test_invalid_task_state_or_polarity_fails_closed(task, error):
    traces = _trace(
        [0, 1, 3, 10, 3, 1, 0],
        [0, 2, 6, 20, 6, 2, 0],
    )
    with pytest.raises(ValueError, match=error):
        build_v16_iimn_feature_centroid_shadow(
            task,
            _rows(),
            current_traces=traces,
            feature_centroid_traces=traces,
        )
