from __future__ import annotations

import csv
from pathlib import Path

import pytest

from paper_v10.validation import (
    REQUIRED_COLUMNS,
    LockboxValidationError,
    evaluate_lockbox,
    load_lockbox_csv,
)


def _row(index: int, **updates):
    truth = f"TRUTHBLOCK{index:04d}"[:14]
    payload = {
        "project_id": "lockbox-1",
        "method_id": "v10-frozen-1",
        "evaluation_role": "primary",
        "evaluation_unit_id": f"unit-{index}",
        "injection_id": f"injection-{index}",
        "block_id": "day-1",
        "entity_id": f"neutral-{index}",
        "entity_level": "neutral_analyte",
        "standard_id": f"STD{index:03d}",
        "family_id": "F01" if index <= 2 else "",
        "is_hard_family": "true" if index <= 2 else "false",
        "truth_inchikey": "",
        "truth_connectivity": truth,
        "truth_formula": "C10H10O2",
        "library_presence": "exact_present",
        "instrument_family": "QE",
        "polarity": "NEG",
        "candidate_inchikey": "",
        "candidate_connectivity": truth,
        "candidate_formula": "C10H10O2",
        "claim_scope": "2d_connectivity",
        "released": "true",
        "abstained": "false",
        "calibration_status": "calibrated",
        "q_value": "0.01",
        "rank_score": "0.9",
        "baseline_candidate_connectivity": "WRONGBLOCK0000",
    }
    payload.update(updates)
    return payload


def _write(path: Path, rows):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=sorted(REQUIRED_COLUMNS))
        writer.writeheader()
        writer.writerows(rows)


def test_lockbox_evaluation_uses_primary_standard_units_not_repeat_scans(tmp_path):
    primary_one = _row(1)
    primary_two = _row(
        2,
        candidate_connectivity="WRONGBLOCK0000",
        baseline_candidate_connectivity="TRUTHBLOCK0002",
        q_value="0.04",
    )
    open_world = _row(
        3,
        library_presence="exact_absent_analog_present",
        candidate_connectivity="",
        candidate_formula="C10H10O2",
        claim_scope="molecular_formula",
        released="true",
        abstained="false",
        calibration_status="calibrated",
        q_value="",
    )
    repeat = _row(
        4,
        evaluation_role="technical_repeat",
        evaluation_unit_id="unit-1-repeat",
        standard_id="STD001",
        candidate_connectivity=primary_one["truth_connectivity"],
        truth_connectivity=primary_one["truth_connectivity"],
        family_id="F01",
        is_hard_family="true",
    )
    path = tmp_path / "calls.csv"
    _write(path, [primary_one, primary_two, open_world, repeat])
    rows = load_lockbox_csv(path)
    result = evaluate_lockbox(rows, bootstrap_iterations=100, bootstrap_seed=7)

    assert result["counts"]["all_rows"] == 4
    assert result["counts"]["primary_units"] == 3
    assert result["closed_set"]["top1_connectivity_accuracy"]["estimate"] == pytest.approx(0.5)
    assert result["open_world"]["connectivity_abstention_sensitivity"]["estimate"] == 1.0
    assert result["paired_baseline_comparison"]["v10_only_correct"] == 1
    assert result["paired_baseline_comparison"]["baseline_only_correct"] == 1
    assert result["technical_repeatability"]["technical_repeat_rows"] == 1


def test_duplicate_primary_unit_fails_closed(tmp_path):
    rows = [_row(1), _row(2, evaluation_unit_id="unit-1")]
    path = tmp_path / "duplicate.csv"
    _write(path, rows)
    with pytest.raises(LockboxValidationError, match="duplicate primary"):
        load_lockbox_csv(path)


def test_released_connectivity_requires_calibrated_qvalue(tmp_path):
    path = tmp_path / "uncalibrated.csv"
    _write(path, [_row(1, calibration_status="unavailable", q_value="")])
    with pytest.raises(LockboxValidationError, match="calibrated q_value"):
        load_lockbox_csv(path)


def test_curated_panel_has_120_unique_standards_and_38_complete_families():
    root = Path(__file__).resolve().parents[1]
    panel = root / "paper_v10" / "templates" / "standard_panel_120.csv"
    with panel.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 120
    assert len({row["standard_id"] for row in rows}) == 120
    families = {row["family_id"] for row in rows if row["family_id"]}
    assert len(families) == 38
    assert sum(row["procurement_phase"] == "A_core_hard64" for row in rows) == 64
    assert all(row["pubchem_cid"] and row["inchikey"] for row in rows)
