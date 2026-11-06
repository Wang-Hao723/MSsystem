from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest

from calibration.metrics import evaluate_probability_calibration
from calibration.models import (
    CalibrationMethod,
    CalibrationProject,
    CalibrationUnitLevel,
    CompetitionPair,
    CompetitionRole,
    ProjectCalibrationStatus,
)
from calibration.project import apply_project_calibration, calibrate_project
from calibration.qvalue import compete_pairs, estimate_qvalue_curve
from schemas.decision import CalibrationStatus, ClaimScope
from schemas.legacy import legacy_row_to_identification


def _method(**updates) -> CalibrationMethod:
    base = CalibrationMethod(
        method_id="test-target-decoy-v1",
        method_version="1",
        score_field="rank_score",
        pseudocount=1.0,
        min_competitions=4,
        min_target_winners=2,
        min_decoy_winners=2,
        allowed_unit_levels=(CalibrationUnitLevel.PRECURSOR_EVENT,),
        allowed_claim_scopes=(ClaimScope.CONNECTIVITY,),
        stratum_fallback=(
            ("unit_level", "claim_scope", "polarity", "instrument_family"),
            ("unit_level", "claim_scope"),
        ),
    )
    return base.model_copy(update=updates)


def _project(**updates) -> CalibrationProject:
    payload = {
        "project_id": "project-1",
        "method_id": "test-target-decoy-v1",
        "runtime_snapshot_id": "runtime-1",
        "library_snapshot_id": "library-1",
        "candidate_space_id": "candidate-space-1",
        "score_field": "rank_score",
        "target_pool_size": 100,
        "decoy_pool_size": 100,
        "decoy_design_id": "mass-matched-real-compound-v1",
        "decoy_generation_seed": 42,
        "method_frozen_before_data": True,
        "independent_from_threshold_tuning": True,
        "complete_target_decoy_competition": True,
        "decoy_identity_overlap_count": 0,
        "instrument_domain": "test-instrument",
    }
    payload.update(updates)
    return CalibrationProject(**payload)


def _pair(index: int, *, target_wins: bool, target_record_id: str | None = None):
    high = 1.0 - index * 0.01
    linked_record = target_record_id is not None and index == 0 and target_wins
    return CompetitionPair(
        competition_id=f"competition-{index}",
        entity_id=f"event-{index}",
        unit_level=CalibrationUnitLevel.PRECURSOR_EVENT,
        claim_scope=ClaimScope.CONNECTIVITY,
        target_record_id=target_record_id or f"record-{index}",
        target_score=0.91 if linked_record else high if target_wins else high - 0.02,
        decoy_score=0.89 if linked_record else high - 0.02 if target_wins else high,
        polarity="NEG",
        instrument_family="QE",
        source_tier="experimental",
    )


def test_exact_score_ties_fail_conservatively_to_decoy():
    pair = _pair(1, target_wins=True).model_copy(
        update={"target_score": 0.8, "decoy_score": 0.8}
    )
    assert compete_pairs([pair])[0].role == CompetitionRole.DECOY


def test_qvalue_curve_is_tie_aware_and_monotone():
    calls = compete_pairs(
        [
            _pair(0, target_wins=True),
            _pair(1, target_wins=True),
            _pair(2, target_wins=False),
            _pair(3, target_wins=True),
            _pair(4, target_wins=False),
        ]
    )
    curve = estimate_qvalue_curve(calls, decoy_scale=1.0, pseudocount=1.0)
    qvalues = [point.q_value for point in curve]
    assert qvalues == sorted(qvalues)
    assert all(0.0 <= value <= 1.0 for value in qvalues)
    assert curve[-1].target_winners == 3
    assert curve[-1].decoy_winners == 2


def test_project_calibration_fails_closed_when_development_data_was_tuned():
    artifact = calibrate_project(
        _project(independent_from_threshold_tuning=False),
        _method(),
        [_pair(i, target_wins=i < 4) for i in range(6)],
    )
    assert artifact.status == ProjectCalibrationStatus.UNAVAILABLE
    assert artifact.assignments == ()
    assert "project_reused_for_threshold_tuning" in artifact.readiness.reason_codes


def test_project_calibration_uses_predeclared_fallback_and_never_reranks():
    pairs = [_pair(i, target_wins=i not in {2, 5}) for i in range(8)]
    artifact = calibrate_project(_project(), _method(), pairs)
    assert artifact.status == ProjectCalibrationStatus.CALIBRATED
    assert len(artifact.assignments) == 6
    assert artifact.target_winners_unassigned == 0
    assert artifact.production_ranking_changed is False
    assert artifact.probability_available is False
    assert all(0.0 <= item.q_value <= 1.0 for item in artifact.assignments)


def test_calibration_attaches_scope_specific_qvalue_without_overwriting_score():
    row = {
        "spectrum_id": "scan=1",
        "precursor_mz": 301.0,
        "ion_mode": "NEG",
        "match": {
            "compound_name": "Candidate",
            "inchikey": "AAAAAAAAAAAAAA-UHFFFAOYSA-N",
            "rank_score": 0.91,
            "final_score": 0.91,
            "confidence": "High",
            "evidence_governance": {
                "policy_id": "test-policy",
                "claim_status": "provisional_needs_independent_injection",
                "claim_scope": "connectivity_level_library_candidate",
                "nodes": [],
            },
        },
    }
    record = legacy_row_to_identification("task-1", row)
    assert record is not None
    pairs = [
        _pair(
            i,
            target_wins=i not in {2, 5},
            target_record_id=record.record_id if i == 0 else None,
        )
        for i in range(8)
    ]
    artifact = calibrate_project(_project(), _method(), pairs)
    calibrated = apply_project_calibration([record], artifact)[0]

    assert calibrated.candidate == record.candidate
    assert calibrated.scores == record.scores
    assert calibrated.scores.legacy_confidence == "High"
    assert calibrated.decision.released is False
    assert calibrated.decision.calibration.status == CalibrationStatus.CALIBRATED
    assert calibrated.decision.calibration.q_connectivity is not None
    assert calibrated.decision.calibration.probability is None


def test_brier_and_ece_require_independent_labels():
    with pytest.raises(ValueError, match="independent"):
        evaluate_probability_calibration(
            [0.8, 0.2], [1, 0], evaluation_set_independent=False
        )
    result = evaluate_probability_calibration(
        [0.8, 0.2], [1, 0], evaluation_set_independent=True, n_bins=2
    )
    assert result["brier_score"] == pytest.approx(0.04)
    assert result["qvalues_used_as_probabilities"] is False


def test_calibration_registry_is_hash_verified_and_has_no_production_entry():
    root = Path(__file__).resolve().parents[1]
    registry = json.loads(
        (root / "config" / "calibration_registry.json").read_text(encoding="utf-8")
    )
    assert registry["production_calibration_id"] is None
    for entry in registry["entries"]:
        path = root / entry["artifact"]
        assert hashlib.sha256(path.read_bytes()).hexdigest() == entry["sha256"]
        assert entry["production_eligible"] is False
