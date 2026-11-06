from __future__ import annotations

from copy import deepcopy

from matching.feature_fragment_benchmark_v3_1.feasibility import (
    build_domain_feasibility_report,
)
from matching.feature_fragment_benchmark_v3_1.independence import (
    canonical_usi_core,
    cluster_mirrors_v31,
    origin_raw_locator,
)


def test_usi_origin_and_peak_hash_are_separate_exact_edges(v31_record_factory) -> None:
    usi_a = v31_record_factory("usi-a")
    usi_b = v31_record_factory("usi-b")
    usi_b["independence"]["exact_peak_hash"] = "b" * 64
    usi_a["source"]["usi"] = "mzspec:PXD000001:run1:scan:7:AAAAAA/2"
    usi_b["source"]["usi"] = "mzspec:PXD000001:run1:scan:7:BBBBBB/2"
    origin_a = v31_record_factory("origin-a")
    origin_b = v31_record_factory("origin-b")
    origin_a["independence"]["exact_peak_hash"] = "c" * 64
    origin_b["independence"]["exact_peak_hash"] = "d" * 64
    for record in (origin_a, origin_b):
        record["source"]["origin_dataset_id"] = "MSV000001"
        record["raw_acquisition"]["raw_file_id"] = "sample.raw"
        record["raw_acquisition"]["scan_id"] = "scan=42"
    hash_a = v31_record_factory("hash-a")
    hash_b = v31_record_factory("hash-b")
    hash_a["independence"]["exact_peak_hash"] = "e" * 64
    hash_b["independence"]["exact_peak_hash"] = "e" * 64
    result = cluster_mirrors_v31(
        [usi_a, usi_b, origin_a, origin_b, hash_a, hash_b]
    )
    assert result["assignment"]["usi-a"] == result["assignment"]["usi-b"]
    assert result["assignment"]["origin-a"] == result["assignment"]["origin-b"]
    assert result["assignment"]["hash-a"] == result["assignment"]["hash-b"]
    assert result["exact_edge_count_by_type"] == {
        "canonical_usi_core": 1,
        "origin_raw_locator": 1,
        "exact_peak_hash": 1,
    }


def test_different_usi_collection_does_not_overmerge_without_other_edge(
    v31_record_factory,
) -> None:
    left = v31_record_factory("left")
    right = v31_record_factory("right")
    left["source"]["usi"] = "mzspec:PXD000001:run1:scan:7"
    right["source"]["usi"] = "mzspec:PXD000002:run1:scan:7"
    left["independence"]["exact_peak_hash"] = "1" * 64
    right["independence"]["exact_peak_hash"] = "2" * 64
    result = cluster_mirrors_v31([left, right])
    assert result["assignment"]["left"] != result["assignment"]["right"]
    assert result["near_duplicate_edge_count"] == 0
    assert result["near_candidate_status"] == "not_evaluated_without_labelled_calibration"


def test_usi_and_origin_helpers_require_complete_locators(v31_record_factory) -> None:
    assert canonical_usi_core("mzspec:PXD1:run:scan:9:PEPTIDE") == "mzspec:PXD1:run:scan:9"
    assert canonical_usi_core("not-a-usi") is None
    record = v31_record_factory()
    assert origin_raw_locator(record) is None
    record["raw_acquisition"].update({"raw_file_id": "a.raw", "scan_id": "10"})
    assert origin_raw_locator(record).startswith("origin:")


def test_feasibility_is_release_bounded_and_partial_pooling_cannot_rescue(
    v31_record_factory,
) -> None:
    ineligible = v31_record_factory("ineligible")
    ineligible["quantitative_eligibility"]["positive_lr_eligible"] = False
    ineligible["quantitative_eligibility"]["ineligible_reason_codes"] = [
        "exact_condition_metadata_incomplete"
    ]
    report = build_domain_feasibility_report(
        [ineligible],
        instrument_domains=[{"domain_id": "domain-a"}],
        thresholds={
            "minimum_target_identities": 1,
            "minimum_background_identities": 1,
            "minimum_target_formula_families": 1,
            "minimum_background_formula_families": 1,
            "minimum_target_sources": 1,
            "minimum_background_sources": 1,
            "minimum_target_laboratories": 1,
            "minimum_background_laboratories": 1,
        },
        frozen_release_id="test-release",
        frozen_corpus_sha256="f" * 64,
    )
    assert report["domains"][0]["status"] == "infeasible_for_frozen_release"
    assert report["permanent_infeasibility_claim_allowed"] is False
    assert report["partial_pooling"]["may_rescue_failed_exact_domain"] is False
    assert "infeasible_forever" not in str(report)


def test_capacity_pass_is_not_misreported_as_rule_feasibility(v31_record_factory) -> None:
    records = []
    for index in range(2):
        row = v31_record_factory(f"eligible-{index}")
        row["compound"]["connectivity"] = f"IDENTITY{index:06d}"
        row["independence"]["two_dimensional_identity_id"] = row["compound"]["connectivity"]
        records.append(row)
    report = build_domain_feasibility_report(
        records,
        instrument_domains=[{"domain_id": "domain-a"}],
        thresholds={
            "minimum_target_identities": 1,
            "minimum_background_identities": 1,
            "minimum_target_formula_families": 1,
            "minimum_background_formula_families": 1,
            "minimum_target_sources": 1,
            "minimum_background_sources": 1,
            "minimum_target_laboratories": 1,
            "minimum_background_laboratories": 1,
        },
        frozen_release_id="test-release",
        frozen_corpus_sha256="f" * 64,
    )
    domain = report["domains"][0]
    assert domain["status"] == "capacity_sufficient_requires_rule_specific_check"
    assert domain["rule_specific_target_background_distribution_assessed"] is False
