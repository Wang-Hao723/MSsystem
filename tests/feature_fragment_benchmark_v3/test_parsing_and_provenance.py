from __future__ import annotations

from copy import deepcopy

from matching.feature_fragment_benchmark.corpus import parse_scan_range as parse_base_range
from matching.feature_fragment_benchmark_v3.parsing import (
    ce_program_has_complete_vector,
    parse_collision_energy_program,
    parse_scan_range,
)
from matching.feature_fragment_benchmark_v3.provenance import (
    make_acquisition_replication_group,
)


def test_real_massbank_hyphenated_scan_range_does_not_create_negative_bound() -> None:
    assert parse_base_range("55-1254") == (55.0, 1254.0)
    assert parse_scan_range("55–1254")["minimum_mz"] == 55.0
    assert parse_scan_range("55 to 1254")["maximum_mz"] == 1254.0
    assert parse_scan_range("-55-1254")["status"] == "unresolved"


def test_collision_energy_modes_keep_native_units_and_complete_programs() -> None:
    fixed = parse_collision_energy_program("35 eV")
    nominal = parse_collision_energy_program("25 % (nominal)")
    ramp = parse_collision_energy_program("Ramp 10-40 eV")
    stepped = parse_collision_energy_program("Stepped NCE 20, 30, 40 NCE")
    ces = parse_collision_energy_program("CE 35 V; CES ±15 V")
    rolling = parse_collision_energy_program("rolling CE = 0.0575*m/z + 9 V")
    mse = parse_collision_energy_program("MSE ramp 10-40 eV")

    assert (fixed["unit"], fixed["mode"], fixed["values"]) == (
        "eV",
        "fixed",
        [35.0],
    )
    assert nominal["unit"] == "percent_nominal"
    assert ramp["values"] == [10.0, 40.0]
    assert stepped["values"] == [20.0, 30.0, 40.0]
    assert ces["values"] == [20.0, 35.0, 50.0]
    assert rolling["equation"] == "0.0575*m/z + 9 V"
    assert mse["precursor_selected"] is False
    assert all(
        row["conversion_applied"] is False
        for row in (fixed, nominal, ramp, stepped, ces, rolling, mse)
    )
    assert all(
        ce_program_has_complete_vector(row)
        for row in (fixed, nominal, ramp, stepped, ces, rolling, mse)
    )


def test_peak_changes_do_not_create_new_acquisition_replication_group(
    v3_record_factory,
) -> None:
    first = v3_record_factory("scan-a")
    second = deepcopy(first)
    second["peaks"] = [[51.0, 1.0], [101.0, 0.5], [151.0, 0.25]]
    second["independence"]["exact_peak_hash"] = "c" * 64

    assert make_acquisition_replication_group(first) == make_acquisition_replication_group(
        second
    )
    assert first["independence"]["peak_hash_used_in_acquisition_replication_group"] is False


def test_missing_isolation_disables_quantitative_lr(v3_record_factory) -> None:
    record = v3_record_factory("no-isolation", isolation_status="unreported")

    assert record["quantitative_eligibility"]["positive_lr_eligible"] is False
    assert record["quantitative_eligibility"]["negative_lr_eligible"] is False
    reasons = set(record["quantitative_eligibility"]["ineligible_reason_codes"])
    assert "precursor_selection_unverified" in reasons
    assert "precursor_isolation_unverified_for_negative_evidence" in reasons


def test_share_alike_license_is_quarantined(v3_record_factory) -> None:
    record = v3_record_factory("share-alike", license_label="CC BY-SA")

    assert record["source"]["license"]["tier"] == "quarantine"
    assert record["quantitative_eligibility"]["positive_lr_eligible"] is False
