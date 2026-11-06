import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest

from matching.masscube_peak_evidence import assess_isf_pair, evaluate_peak_quality


AUTHOR_SOURCE = (
    Path(__file__).resolve().parent
    / "fixtures"
    / "masscube_1_0_22_feature_evaluation.py"
)


def _load_author_module():
    spec = importlib.util.spec_from_file_location("masscube_author_feature_evaluation", AUTHOR_SOURCE)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_valid_peak_metrics_exactly_match_author_and_are_invariant():
    author = _load_author_module()
    rt = np.array([1.00, 1.02, 1.04, 1.06, 1.08, 1.10, 1.12, 1.14, 1.16])
    intensity = np.array([1.0, 6.0, 28.0, 72.0, 100.0, 70.0, 25.0, 5.0, 1.0])
    rt_before = rt.copy()
    intensity_before = intensity.copy()

    result = evaluate_peak_quality(rt, intensity)
    shifted_scaled = evaluate_peak_quality(rt + 8.0, intensity * 37.0)

    assert result["status"] == "ok"
    assert result["gaussian_similarity"] == pytest.approx(
        author.calculate_gaussian_similarity(rt, intensity), abs=1e-15
    )
    assert result["noise_score"] == author.calculate_noise_level(intensity)
    assert result["asymmetry_factor"] == author.calculate_asymmetry_factor(intensity)
    assert shifted_scaled["gaussian_similarity"] == pytest.approx(result["gaussian_similarity"])
    assert shifted_scaled["noise_score"] == result["noise_score"]
    assert shifted_scaled["asymmetry_factor"] == result["asymmetry_factor"]
    assert shifted_scaled["apex_rt_minutes"] == pytest.approx(result["apex_rt_minutes"] + 8.0)
    np.testing.assert_array_equal(rt, rt_before)
    np.testing.assert_array_equal(intensity, intensity_before)


def test_smooth_peak_has_less_noise_than_oscillating_trace():
    rt = np.arange(9, dtype=float)
    smooth = evaluate_peak_quality(rt, [1, 5, 20, 60, 100, 60, 20, 5, 1])
    noisy = evaluate_peak_quality(rt, [10, 35, 20, 80, 100, 55, 70, 25, 10])
    assert smooth["status"] == noisy["status"] == "ok"
    assert smooth["noise_score"] < noisy["noise_score"]


@pytest.mark.parametrize(
    "values",
    [
        [1, 8, 35, 100, 82, 61, 42, 29, 15, 8, 1],
        [1, 8, 15, 29, 42, 61, 82, 100, 35, 8, 1],
        [8, 25, 19, 75, 100, 62, 77, 23, 34, 8, 1],
        [0, 0, 1, 5, 100, 5, 1, 0, 0, 0, 0],
        [10, 10, 25, 50, 100, 50, 25, 10, 10, 8, 2],
        [60, 75, 90, 100, 85, 80, 65, 50, 40, 25, 15],
    ],
)
def test_asymmetric_noisy_and_sparse_trace_metrics_match_author(values):
    author = _load_author_module()
    # Uneven scan intervals check the actual RT-based Gaussian calculation.
    rt = np.cumsum(np.linspace(.005, .02, len(values)))
    intensity = np.asarray(values, dtype=float)
    actual = evaluate_peak_quality(rt, intensity)
    assert actual["status"] == "ok"
    assert actual["gaussian_similarity"] == pytest.approx(author.calculate_gaussian_similarity(rt, intensity), abs=1e-15)
    assert actual["noise_score"] == author.calculate_noise_level(intensity)
    assert actual["asymmetry_factor"] == author.calculate_asymmetry_factor(intensity)


@pytest.mark.parametrize(
    ("rt", "intensity", "flag"),
    [
        ([1, 2, 3, 4], [1, 2, 3, 1], "too_few_points"),
        ([1, 2, 3, 4, 5], [5, 4, 3, 2, 1], "endpoint_apex"),
        ([1, 2, 2, 4, 5], [1, 3, 5, 3, 1], "non_increasing_rt"),
        ([1, 2, 3, 4, 5], [0, 0, 0, 0, 0], "zero_intensity"),
        ([1, 2, 3, 4, 5], [2, 2, 2, 2, 2], "constant_intensity"),
        ([1, 2, 3, 4, 5], [1, 2, np.nan, 2, 1], "nonfinite_values"),
        ([1, 2, 3, 4, 5], [1, 2, -1, 2, 1], "negative_intensity"),
    ],
)
def test_invalid_peak_traces_fail_explicitly_and_json_safely(rt, intensity, flag):
    result = evaluate_peak_quality(rt, intensity)
    assert result["status"] == "invalid_peak"
    assert flag in result["flags"]
    assert result["gaussian_similarity"] is None
    assert result["noise_score"] is None
    assert result["asymmetry_factor"] is None
    json.dumps(result, allow_nan=False)


def _trace(values, *, shift=0.0, scan_ids=True):
    result = {
        "rt": [1.00 + shift, 1.02 + shift, 1.04 + shift, 1.06 + shift, 1.08 + shift, 1.10 + shift, 1.12 + shift],
        "intensity": list(values),
    }
    if scan_ids:
        result["scan_ids"] = [101, 102, 103, 104, 105, 106, 107]
    return result


def test_true_isf_pair_requires_ms2_mass_coelution_and_correlated_common_scans():
    parent = _trace([2, 12, 50, 100, 52, 13, 2])
    child = _trace([1, 6, 25, 51, 27, 7, 1])
    parent_before = json.loads(json.dumps(parent))
    child_before = json.loads(json.dumps(child))

    result = assess_isf_pair(
        301.1234,
        151.0612,
        [(80.0, 20.0), (151.064, 500.0)],
        parent,
        child,
        parent_polarity="positive",
        child_polarity="POSITIVE",
    )

    assert result["supported"] is True
    assert result["status"] == "supported"
    assert result["matched_fragment_mz"] == 151.064
    assert result["mass_error_da"] == pytest.approx(0.0028)
    assert result["common_points"] == 7
    assert result["correlation"] > 0.99
    assert result["apex_rt_delta_min"] == 0.0
    assert result["evidence_level"] == "candidate_isf_relation"
    assert result["automatic_remove"] is False
    assert parent == parent_before and child == child_before
    json.dumps(result, allow_nan=False)


def test_missing_parent_ms2_fragment_fails_despite_perfect_coelution():
    trace = _trace([1, 7, 30, 80, 30, 7, 1])
    result = assess_isf_pair(300.0, 150.0, [(149.98, 1000)], trace, trace)
    assert result["supported"] is False
    assert "child_mass_absent_from_parent_ms2" in result["reasons"]
    assert result["correlation"] == pytest.approx(1.0)


def test_unrelated_peak_shape_fails_correlation_gate():
    parent = _trace([1, 5, 20, 100, 30, 8, 2])
    child = _trace([50, 10, 70, 5, 90, 20, 60])
    result = assess_isf_pair(300, 150, [(150, 100)], parent, child)
    assert result["supported"] is False
    assert "correlation_not_above_threshold" in result["reasons"]


@pytest.mark.parametrize(
    ("parent", "child", "kwargs", "reason"),
    [
        (
            _trace([1, 2, 4, 8, 4, 2, 1]),
            {"rt": [1.00, 1.02, 1.04, 1.06], "intensity": [1, 2, 4, 2]},
            {},
            "child_trace_too_short",
        ),
        (
            _trace([1, 2, 4, 8, 4, 2, 1]),
            _trace([3, 3, 3, 3, 3, 3, 3]),
            {},
            "child_trace_constant",
        ),
        (
            _trace([1, 2, 4, 8, 4, 2, 1]),
            _trace([1, 2, 4, 8, 4, 2, 1], shift=0.001, scan_ids=False),
            {},
            "too_few_common_points",
        ),
        (
            _trace([1, 2, 4, 8, 4, 2, 1]),
            _trace([1, 2, 4, 8, 4, 2, 1]),
            {"parent_polarity": "positive", "child_polarity": "negative"},
            "polarity_conflict",
        ),
    ],
)
def test_isf_pair_fails_closed_for_insufficient_or_conflicting_evidence(parent, child, kwargs, reason):
    result = assess_isf_pair(300, 150, [(150, 100)], parent, child, **kwargs)
    assert result["supported"] is False
    assert reason in result["reasons"]
    assert result["automatic_remove"] is False
    json.dumps(result, allow_nan=False)


def test_rt_keys_are_used_when_scan_ids_are_not_available():
    parent = _trace([1, 4, 15, 50, 16, 4, 1], scan_ids=False)
    child = _trace([2, 8, 30, 100, 32, 8, 2], scan_ids=False)
    result = assess_isf_pair(400, 200, [(200.0, 10)], parent, child)
    assert result["supported"] is True
    assert result["common_points"] == 7


def test_single_coincident_spike_cannot_pass_via_zero_padding():
    trace = {
        "rt": [1 + i * .01 for i in range(21)],
        "intensity": [0.] * 10 + [100.] + [0.] * 10,
        "scan_ids": list(range(21)),
    }
    result = assess_isf_pair(300, 150, [(150, 100)], trace, trace)
    assert result["common_points"] == 21
    assert result["common_positive_points"] == 1
    assert result["correlation"] == pytest.approx(1.)
    assert result["supported"] is False
    assert "too_few_common_positive_points" in result["reasons"]


def test_five_common_measured_points_with_zero_padding_remain_eligible():
    trace = {
        "rt": [1 + i * .01 for i in range(9)],
        "intensity": [0, 0, 2, 12, 50, 12, 2, 0, 0],
        "scan_ids": list(range(9)),
    }
    result = assess_isf_pair(300, 150, [(150, 100)], trace, trace)
    assert result["common_positive_points"] == 5
    assert result["supported"] is True
