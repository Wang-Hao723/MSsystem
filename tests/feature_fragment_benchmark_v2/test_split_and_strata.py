from __future__ import annotations

from matching.feature_fragment_benchmark.split import (
    build_isolated_split,
    stratum_id,
    stratum_key,
)
from matching.feature_fragment_benchmark.statistics import build_identity_observations


def test_split_has_zero_identity_formula_mirror_overlap_and_deterministic_seal(
    benchmark_record_factory,
) -> None:
    records = [
        benchmark_record_factory(
            "bridge-a", identity="identity-a", formula="C10H20O2", mirror_cluster_id="m-a"
        ),
        benchmark_record_factory(
            "bridge-b", identity="identity-b", formula="C10H20O2", mirror_cluster_id="m-b"
        ),
        benchmark_record_factory(
            "bridge-c", identity="identity-c", formula="C11H22O2", mirror_cluster_id="m-b"
        ),
    ]
    records.extend(
        benchmark_record_factory(
            f"independent-{index}",
            identity=f"identity-{index + 10}",
            formula=f"C{index + 12}H{2 * index + 24}O2",
            mirror_cluster_id=f"m-{index + 10}",
        )
        for index in range(12)
    )

    forward = build_isolated_split(records, seed="deterministic-seal")
    reverse = build_isolated_split(list(reversed(records)), seed="deterministic-seal")

    assert forward["split_feasible"] is True
    assert forward["isolation_audit"]["forbidden_overlap_total"] == 0
    assert forward["isolation_audit"]["zero_forbidden_overlap"] is True
    assert forward["split_manifest_sha256"] == reverse["split_manifest_sha256"]
    assert (
        forward["sealed_test_spectrum_ids_sha256"]
        == reverse["sealed_test_spectrum_ids_sha256"]
    )
    assert forward["assignment"] == reverse["assignment"]
    bridge_splits = {
        forward["assignment"][name]["split"]
        for name in ("bridge-a", "bridge-b", "bridge-c")
    }
    assert len(bridge_splits) == 1


def test_ce_unit_mode_adduct_and_polarity_are_distinct_strata(
    benchmark_record_factory,
) -> None:
    records = [
        benchmark_record_factory("ev-fixed", identity="shared-stratum-probe"),
        benchmark_record_factory(
            "nce-fixed", identity="shared-stratum-probe", ce_unit="NCE"
        ),
        benchmark_record_factory(
            "ev-ramp",
            identity="shared-stratum-probe",
            ce_mode="ramp",
            ce_value=None,
            ce_minimum=10.0,
            ce_maximum=30.0,
        ),
        benchmark_record_factory(
            "sodium", identity="shared-stratum-probe", adduct="[M+Na]+"
        ),
        benchmark_record_factory(
            "negative",
            identity="shared-stratum-probe",
            adduct="[M-H]-",
            polarity="negative",
            charge=-1,
        ),
    ]

    keys = [stratum_key(record) for record in records]
    assert len(set(keys)) == len(records)
    assert len({stratum_id(record) for record in records}) == len(records)
    assert len(build_identity_observations(records)) == len(records)
    assert keys[0][3:5] == ("eV", "fixed")
    assert keys[1][3:5] == ("NCE", "fixed")
    assert keys[2][3:5] == ("eV", "ramp")
    assert keys[3][6] == "[M+Na]+"
    assert keys[4][7:9] == ("negative", -1)


def test_technical_repeats_collapse_to_one_identity_observation(
    benchmark_record_factory,
) -> None:
    records = [
        benchmark_record_factory(
            "repeat-a",
            identity="shared-identity",
            mirror_cluster_id="mirror-a",
        ),
        benchmark_record_factory(
            "repeat-b",
            identity="shared-identity",
            mirror_cluster_id="mirror-a",
            peaks=[[200.0, 25.0], [50.0, 100.0], [100.0, 50.0]],
        ),
        benchmark_record_factory(
            "repeat-c",
            identity="shared-identity",
            mirror_cluster_id="mirror-c",
            peaks=[[50.0, 1.0], [100.0, 0.45], [201.0, 0.20]],
        ),
    ]
    for record in records:
        record["independence"]["independent_fragmentation_group_id"] = (
            "technical-series-1"
        )

    observations = build_identity_observations(records)

    assert len(observations) == 1
    assert observations[0]["identity"] == "shared-identity"
    assert observations[0]["record_count"] == 3
    assert observations[0]["mirror_cluster_count"] == 2
