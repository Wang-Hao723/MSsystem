from __future__ import annotations

import json
from copy import deepcopy

import pytest

from matching.feature_fragment_benchmark.split import stratum_id, stratum_key
from matching.feature_fragment_benchmark.statistics import (
    _exploratory_absence_eligible,
    _strict_absence_eligible,
    build_identity_observations,
    calibrate_rules,
    discover_rules,
    test_frozen_rules as evaluate_frozen_rules,
)


FEATURE_A = "product_ion:10000:0.01"
FEATURE_B = "product_ion:12000:0.01"


def test_independent_group_count_is_not_mirror_deduplicated_record_count(
    benchmark_record_factory,
) -> None:
    records = [
        benchmark_record_factory(
            f"technical-repeat-{index}",
            identity="shared-identity",
            formula="C10H20O2",
        )
        for index in range(3)
    ]
    for record in records:
        record["independence"]["independent_fragmentation_group_id"] = "experiment-1"
    records[0]["independence"]["mirror_cluster_id"] = "mirror-a"
    records[1]["independence"]["mirror_cluster_id"] = "mirror-a"
    records[2]["independence"]["mirror_cluster_id"] = "mirror-b"

    observations = build_identity_observations(records)

    assert len(observations) == 1
    assert observations[0]["mirror_cluster_count"] == 2
    assert observations[0]["mirror_deduplicated_record_count"] == 2
    assert observations[0]["independent_fragmentation_group_count"] == 1


def _observation(
    benchmark_record_factory,
    ordinal: int,
    *,
    chemical_class: str,
    features: set[str],
    source_id: str = "source-a",
    laboratory_id: str = "lab-a",
) -> dict:
    record = benchmark_record_factory(
        f"observation-{chemical_class}-{ordinal}",
        identity=f"identity-{chemical_class}-{ordinal}",
        formula=f"C{ordinal + 10}H{ordinal + 20}O2",
        chemical_class=chemical_class,
        source_id=source_id,
        laboratory_id=laboratory_id,
    )
    return {
        "observation_id": f"observation-{chemical_class}-{ordinal}",
        "identity": record["compound"]["connectivity"],
        "formula": record["compound"]["formula"],
        "chemical_class": chemical_class,
        "chemical_class_provenance": record["compound"]["chemical_class"],
        "stratum_id": stratum_id(record),
        "stratum_key": list(stratum_key(record)),
        "features": set(features),
        "record_count": 1,
        "mirror_cluster_count": 1,
        "independent_fragmentation_group_count": 1,
        "sources": [source_id],
        "laboratories": [laboratory_id],
        "research_ids": [record["source"]["research_id"]],
        "precursor_mz": record["acquisition"]["precursor_mz"],
        "representative_record": record,
        "independent_records": [record],
    }


@pytest.mark.parametrize(
    "false_quality_field",
    [
        "centroid_peak_list",
        "positive_finite_peaks",
        "at_least_three_peaks",
        "collision_energy_verified",
        "instrument_domain_verified",
        "reported_mass_range_present",
        "reported_mass_range_contains_precursor",
    ],
)
def test_false_quality_gate_cannot_enter_strict_absence_denominator(
    benchmark_record_factory, false_quality_field: str
) -> None:
    record = benchmark_record_factory(
        "quality-false", quality_overrides={false_quality_field: False}
    )
    assert _strict_absence_eligible(record, FEATURE_A) is False


@pytest.mark.parametrize(
    "false_quality_field",
    [
        "centroid_peak_list",
        "positive_finite_peaks",
        "at_least_three_peaks",
        "collision_energy_verified",
        "instrument_domain_verified",
    ],
)
def test_false_quality_gate_cannot_enter_exploratory_absence_denominator(
    benchmark_record_factory, false_quality_field: str
) -> None:
    record = benchmark_record_factory(
        "exploratory-quality-false",
        quality_overrides={false_quality_field: False},
    )
    assert _exploratory_absence_eligible(record, FEATURE_A) is False


def test_missing_reported_mass_range_cannot_enter_strict_absence_denominator(
    benchmark_record_factory,
) -> None:
    record = benchmark_record_factory("no-mass-range")
    record["acquisition"]["reported_mass_range"] = None
    record["quality"]["reported_mass_range_present"] = False
    record["quality"]["reported_mass_range_contains_precursor"] = False

    assert _strict_absence_eligible(record, FEATURE_A) is False


def test_sealed_test_changes_do_not_mutate_or_rediscover_rules(
    benchmark_record_factory,
) -> None:
    discovery = [
        *[
            _observation(
                benchmark_record_factory,
                index,
                chemical_class="Target",
                features={FEATURE_A, FEATURE_B},
            )
            for index in range(6)
        ],
        *[
            _observation(
                benchmark_record_factory,
                index,
                chemical_class="Background",
                features=set(),
            )
            for index in range(6)
        ],
    ]
    calibration = deepcopy(discovery)
    discovered = discover_rules(
        discovery,
        minimum_target_identities=4,
        minimum_background_identities=4,
        minimum_target_present=4,
        maximum_q_value=1.0,
        minimum_lr_positive=1.1,
    )
    assert discovered
    calibrated = calibrate_rules(
        discovered,
        calibration,
        minimum_target_identities=4,
        minimum_background_identities=4,
    )
    frozen_discovery = json.dumps(discovered, sort_keys=True)
    frozen_calibration = json.dumps(calibrated, sort_keys=True)

    test_a = deepcopy(discovery)
    test_b = deepcopy(discovery)
    for row in test_b:
        row["features"] = {FEATURE_B} if row["chemical_class"] == "Background" else set()

    evaluated_a = evaluate_frozen_rules(calibrated, test_a)
    evaluated_b = evaluate_frozen_rules(calibrated, test_b)

    assert json.dumps(discovered, sort_keys=True) == frozen_discovery
    assert json.dumps(calibrated, sort_keys=True) == frozen_calibration
    assert [rule["rule_id"] for rule in evaluated_a] == [
        rule["rule_id"] for rule in evaluated_b
    ]
    assert [rule["frozen_parameter_sha256"] for rule in evaluated_a] == [
        rule["frozen_parameter_sha256"] for rule in evaluated_b
    ]


def test_small_sealed_test_stays_provisional_when_lr_interval_is_imprecise(
    benchmark_record_factory,
) -> None:
    discovery = [
        *[
            _observation(
                benchmark_record_factory,
                index,
                chemical_class="Target",
                features={FEATURE_A},
            )
            for index in range(4)
        ],
        *[
            _observation(
                benchmark_record_factory,
                index,
                chemical_class="Background",
                features=set(),
            )
            for index in range(4)
        ],
    ]
    discovered = discover_rules(
        discovery,
        minimum_target_identities=2,
        minimum_background_identities=2,
        minimum_target_present=2,
        maximum_q_value=1.0,
        minimum_lr_positive=1.1,
    )
    calibrated = calibrate_rules(
        discovered,
        deepcopy(discovery),
        minimum_target_identities=2,
        minimum_background_identities=2,
    )

    tested = evaluate_frozen_rules(calibrated, deepcopy(discovery))

    assert tested
    assert {rule["evidence_grade"] for rule in tested} == {"provisional_shadow"}
    for rule in tested:
        unmet = set(rule["empirical_validation_gate_unmet"])
        assert "test_target_identities_lt_30" in unmet
        assert "test_background_identities_lt_30" in unmet
        assert "lr_positive_interval_too_wide" in unmet


def test_fragment_group_is_frozen_before_test_and_not_selected_by_dynamic_lr(
    benchmark_record_factory,
) -> None:
    discovery = [
        *[
            _observation(
                benchmark_record_factory,
                index,
                chemical_class="Target",
                features={FEATURE_A, FEATURE_B},
            )
            for index in range(6)
        ],
        *[
            _observation(
                benchmark_record_factory,
                index,
                chemical_class="Background",
                features=set(),
            )
            for index in range(6)
        ],
    ]
    discovered = discover_rules(
        discovery,
        minimum_target_identities=4,
        minimum_background_identities=4,
        minimum_target_present=4,
        maximum_q_value=1.0,
        minimum_lr_positive=1.1,
        fragment_group_jaccard=0.85,
    )

    assert len(discovered) >= 2
    group_ids = {rule["independent_fragmentation_group"] for rule in discovered}
    assert len(group_ids) == 1
    assert {
        rule["fragment_group_policy"]["representative_selection"]
        for rule in discovered
    } == {"frozen_before_test_not_dynamic_max_lr"}

    calibrated = calibrate_rules(
        discovered,
        deepcopy(discovery),
        minimum_target_identities=4,
        minimum_background_identities=4,
    )
    test_observations = deepcopy(discovery)
    for row in test_observations:
        row["features"] = (
            {FEATURE_A} if row["chemical_class"] == "Target" else {FEATURE_B}
        )
    tested = evaluate_frozen_rules(calibrated, test_observations)

    assert {rule["independent_fragmentation_group"] for rule in tested} == group_ids
    assert all(
        rule["fragment_group_policy"]["runtime_combination"]
        == "at_most_one_frozen_group_contribution"
        for rule in tested
    )
