from __future__ import annotations

from matching.feature_fragment_benchmark_v3.independence import (
    calibrate_near_mirror_threshold,
    cluster_mirrors,
)
from matching.feature_fragment_benchmark_v3.split import build_origin_isolated_split


def _identity(index: int) -> str:
    alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    digits = []
    value = index
    for _ in range(14):
        digits.append(alphabet[value % 26])
        value //= 26
    return "".join(reversed(digits))


def test_exact_mirrors_cluster_but_uncalibrated_near_mirrors_do_not_link(
    v3_record_factory,
) -> None:
    first = v3_record_factory("mirror-a")
    second = v3_record_factory("mirror-b")
    third = v3_record_factory(
        "near-only",
        peaks=[[50.001, 100.0], [100.001, 50.0], [150.001, 24.0]],
    )
    result = cluster_mirrors([first, second, third])

    assert result["exact_duplicate_edge_count"] == 1
    assert result["near_duplicate_edge_count"] == 0
    assert result["assignment"]["mirror-a"] == result["assignment"]["mirror-b"]
    assert result["assignment"]["mirror-a"] != result["assignment"]["near-only"]
    assert result["near_candidate_status"] == "not_evaluated_without_labelled_calibration"


def test_near_mirror_calibration_requires_labelled_positive_and_negative_pairs() -> None:
    insufficient = calibrate_near_mirror_threshold(
        [{"is_duplicate": True, "cosine": 0.99, "coverage_left": 1.0, "coverage_right": 1.0}]
    )
    assert insufficient["near_mirror_links_allowed"] is False

    labelled = []
    for index in range(20):
        labelled.append(
            {
                "is_duplicate": True,
                "cosine": 0.995 - index / 10000,
                "coverage_left": 0.95,
                "coverage_right": 0.95,
            }
        )
        labelled.append(
            {
                "is_duplicate": False,
                "cosine": 0.90 - index / 1000,
                "coverage_left": 0.95,
                "coverage_right": 0.95,
            }
        )
    calibrated = calibrate_near_mirror_threshold(labelled)
    assert calibrated["status"] == "calibrated"
    assert calibrated["near_mirror_links_allowed"] is True
    assert calibrated["false_positive_rate"] <= 0.01


def test_source_lab_identity_formula_mirror_and_acquisition_isolation(
    v3_record_factory,
) -> None:
    records = []
    ordinal = 0
    for lab_index in range(6):
        for source_index in range(2):
            source = f"source-{lab_index}-{source_index}"
            laboratory = f"lab-{lab_index}"
            for _ in range(8):
                identity = _identity(ordinal)
                records.append(
                    v3_record_factory(
                        f"spectrum-{ordinal:03d}",
                        connectivity=identity,
                        formula=f"C{ordinal + 10}H{ordinal + 20}O2",
                        source_id=source,
                        laboratory_id=laboratory,
                        peaks=[
                            [50.0, 100.0],
                            [100.0 + ordinal / 1000.0, 50.0],
                            [150.0 + ordinal / 1000.0, 25.0],
                        ],
                    )
                )
                ordinal += 1
    result = build_origin_isolated_split(records, seed="v3-split-test")

    assert result["split_feasible"] is True
    assert result["isolation_audit"]["zero_forbidden_overlap"] is True
    assert all(result["split_record_counts"][name] > 0 for name in ("discovery", "calibration", "test"))
    overlaps = result["isolation_audit"]["overlap_counts"]
    for dimension in ("identity", "formula", "mirror", "acquisition"):
        assert all(value == 0 for value in overlaps[dimension].values())
    for dimension in ("source", "laboratory"):
        assert overlaps[dimension]["discovery__test"] == 0
        assert overlaps[dimension]["calibration__test"] == 0
