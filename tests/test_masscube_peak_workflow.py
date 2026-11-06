import copy
import json

from matching.masscube_peak_workflow import annotate_feature_evidence


def _trace(scale=1.0, *, scan_offset=0):
    return {
        "rt": [1.00, 1.02, 1.04, 1.06, 1.08, 1.10, 1.12],
        "intensity": [scale * value for value in (1, 7, 30, 80, 30, 7, 1)],
        "scan_ids": [scan_offset + value for value in (101, 102, 103, 104, 105, 106, 107)],
    }


def _feature(feature_id, mz, *, rt=1.06, height=100.0, isotope=False):
    return {
        "id": feature_id,
        "mz": mz,
        "rt_minutes": rt,
        "height": height,
        "is_isotope": isotope,
    }


def _by_id(result):
    return {feature["id"]: feature for feature in result["features"]}


def test_known_relation_requires_fragment_and_annotates_quality_without_mutation():
    features = [_feature("parent", 301.1, height=1000), _feature("child", 151.05, height=200)]
    records = [
        {
            "feature_id": "parent",
            "peaks": [(70.0, 4), (151.054, 90)],
            "tic": 1000,
            "spectrum_id": "p1",
        }
    ]
    traces = {"parent": _trace(10), "child": _trace(2)}
    untouched = copy.deepcopy((features, records, traces))

    result = annotate_feature_evidence(features, records, traces, polarity="positive")

    assert result["summary"]["feature_count"] == 2
    assert result["summary"]["quality_ok_count"] == 2
    assert result["summary"]["candidate_prefilter_count"] == 1
    assert result["summary"]["assessed_pair_count"] == 1
    assert result["summary"]["isf_relation_count"] == 1
    relation = result["isf_relations"][0]
    assert (relation["parent_id"], relation["child_id"]) == ("parent", "child")
    assert relation["parent_spectrum_id"] == "p1"
    assert relation["supported"] is True
    assert relation["automatic_remove"] is False
    assert _by_id(result)["child"]["candidate_isf_parent_ids"] == ["parent"]
    assert all(feature["automatic_remove"] is False for feature in result["features"])
    assert (features, records, traces) == untouched
    json.dumps(result, allow_nan=False)


def test_perfectly_coeluting_feature_without_parent_ms2_fragment_is_not_related():
    features = [_feature("p", 300, height=500), _feature("c", 150, height=100)]
    records = [{"feature_id": "p", "peaks": [(149.98, 100)], "tic": 10, "spectrum_id": "s"}]
    result = annotate_feature_evidence(features, records, {"p": _trace(), "c": _trace(0.5)})
    assert result["isf_relations"] == []
    assert _by_id(result)["c"]["candidate_isf_parent_ids"] == []


def test_highest_tic_spectrum_is_fixed_and_spectrum_id_breaks_ties():
    features = [_feature("p", 300, height=500), _feature("c", 150, height=100)]
    traces = {"p": _trace(), "c": _trace(0.5)}
    higher_tic_without_fragment = [
        {"feature_id": "p", "peaks": [(150.0, 100)], "tic": 10, "spectrum_id": "low"},
        {"feature_id": "p", "peaks": [(100.0, 100)], "tic": 20, "spectrum_id": "high"},
    ]
    assert annotate_feature_evidence(features, higher_tic_without_fragment, traces)["isf_relations"] == []

    tied = [
        {"feature_id": "p", "peaks": [(100.0, 100)], "tic": 20, "spectrum_id": "B"},
        {"feature_id": "p", "peaks": [(150.0, 100)], "tic": 20, "spectrum_id": "A"},
    ]
    relation = annotate_feature_evidence(features, tied, traces)["isf_relations"][0]
    assert relation["parent_spectrum_id"] == "A"


def test_missing_ms1_trace_is_explicit_and_prevents_relation():
    features = [_feature("p", 300, height=500), _feature("c", 150, height=100)]
    records = [{"feature_id": "p", "peaks": [(150, 100)], "tic": 10, "spectrum_id": "s"}]
    result = annotate_feature_evidence(features, records, {"p": _trace()})
    assert _by_id(result)["c"]["quality"]["status"] == "missing_trace"
    assert _by_id(result)["c"]["quality"]["flags"] == ["missing_trace"]
    assert result["summary"]["quality_missing_trace_count"] == 1
    assert result["isf_relations"] == []


def test_equal_or_higher_mass_isotope_height_and_rt_filters_exclude_children():
    features = [
        _feature("p", 300, height=100),
        _feature("equal", 300, height=10),
        _feature("higher", 301, height=10),
        _feature("isotope", 150.000, height=10, isotope=True),
        _feature("too_tall", 150.004, height=300),
        _feature("far_rt", 150.008, rt=1.12, height=10),
    ]
    records = [
        {
            "feature_id": "p",
            "peaks": [(150.004, 100), (300, 100), (301, 100)],
            "tic": 1,
            "spectrum_id": "s",
        }
    ]
    traces = {feature["id"]: _trace() for feature in features}
    result = annotate_feature_evidence(features, records, traces)
    assert result["isf_relations"] == []


def test_multiple_supported_parents_are_retained_and_marked_ambiguous():
    features = [
        _feature("p-high-mz", 400, height=500),
        _feature("p-low-mz", 350, height=1000),
        _feature("child", 150, height=100),
    ]
    records = [
        {"feature_id": "p-low-mz", "peaks": [(150, 20)], "tic": 50, "spectrum_id": "lo"},
        {"feature_id": "p-high-mz", "peaks": [(150, 30)], "tic": 60, "spectrum_id": "hi"},
    ]
    traces = {feature["id"]: _trace(index + 1) for index, feature in enumerate(features)}
    result = annotate_feature_evidence(features, records, traces)

    assert [(item["parent_id"], item["child_id"]) for item in result["isf_relations"]] == [
        ("p-high-mz", "child"),
        ("p-low-mz", "child"),
    ]
    assert all(item["ambiguous_parent"] is True for item in result["isf_relations"])
    assert _by_id(result)["child"]["candidate_isf_parent_ids"] == ["p-high-mz", "p-low-mz"]
    assert result["summary"]["ambiguous_child_count"] == 1


def test_extra_identity_or_truth_fields_are_not_consumed_or_propagated():
    feature = {
        **_feature("x", 200),
        "candidate_identity": "must not enter evidence",
        "library_score": 0.999,
        "truth_label": "positive",
    }
    result = annotate_feature_evidence([feature], [], {"x": _trace()})
    emitted = result["features"][0]
    assert "candidate_identity" not in emitted
    assert "library_score" not in emitted
    assert "truth_label" not in emitted
    assert result["summary"]["automatic_remove_count"] == 0
