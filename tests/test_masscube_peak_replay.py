"""Checks for measured binding and equal-budget repeat endpoints."""
import math

from eval.evaluate_masscube_peak_evidence import (
    reciprocal_match,
    repeat_endpoints,
    uniquely_bound_feature_id,
)


def test_parent_binding_accepts_only_existing_unique_attachment_enum():
    assert uniquely_bound_feature_id({"feature_id": "LF1", "link_status": "linked_openms_feature"}) == "LF1"
    for status in (None, "linked", "ambiguous_best_feature", "no_feature_in_mz_rt_window"):
        assert uniquely_bound_feature_id({"feature_id": "LF1", "link_status": status}) is None


def test_reciprocal_match_is_one_to_one_and_respects_mass_and_rt():
    first = [
        {"id": "a", "mz": 300.0, "rt_minutes": 2.0},
        {"id": "b", "mz": 300.0001, "rt_minutes": 2.08},
        {"id": "c", "mz": 400.0, "rt_minutes": 2.0},
    ]
    second = [
        {"id": "x", "mz": 300.0, "rt_minutes": 2.0},
        {"id": "y", "mz": 400.0, "rt_minutes": 2.11},
    ]
    assert reciprocal_match(first, second) == {"a": "x"}


def test_repeat_comparison_keeps_equal_budget_and_full_denominator():
    protocol = {"runs": [{"task_id": str(i), "technical_repeat_group": "same"} for i in range(3)]}
    results = {}
    for run in range(3):
        features = []
        for i in range(5):
            features.append({
                "id": str(i), "mz": 200.0 + i + (0.2 if i == 4 and run else 0),
                "rt_minutes": 2.0, "area": i + 1,
                "quality": {"status": "ok", "gaussian_similarity": 1 - i * 0.1, "noise_score": 0.0, "asymmetry_factor": 1.0},
            })
        results[str(run)] = {"features": features, "isf_relations": []}
    endpoint = repeat_endpoints(protocol, results)
    expected = 3 * math.ceil(5 * .5)
    assert endpoint["pooled"]["baseline_area"]["selected_count"] == expected
    assert endpoint["pooled"]["quality"]["selected_count"] == expected
    assert endpoint["pooled"]["quality"]["matched_both_other_repeats"] == expected
    assert all(len(x["baseline_area"]["selected_ids"]) == len(x["quality"]["selected_ids"]) == 3 for x in endpoint["per_run"])
