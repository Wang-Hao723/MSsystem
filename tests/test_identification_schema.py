from __future__ import annotations

import json
from pathlib import Path

import pytest
from pydantic import ValidationError

from pipeline.result_contract import (
    read_identification_contract,
    write_identification_contract,
)
from schemas.decision import CalibrationStatus, CalibrationVector
from schemas.evidence import (
    EvidenceApplicability,
    EvidenceNode,
    EvidenceProvenance,
    EvidenceStatus,
)
from schemas.legacy import legacy_row_to_identification


def _legacy_row() -> dict:
    return {
        "spectrum_id": "scan=42",
        "precursor_mz": 301.0712,
        "retention_time": 5.5,
        "ion_mode": "NEG",
        "precursor_feature_id": "PF0001",
        "match": {
            "compound_name": "Example",
            "inchikey": "AAAAAAAAAAAAAA-UHFFFAOYSA-N",
            "formula": "C15H10O7",
            "smiles": "O=c1cc(O)c2cc(O)cc(O)c2o1",
            "adduct": "[M-H]-",
            "source_tier": "experimental",
            "lib_spectrum_id": "LIB1",
            "cosine": 0.91,
            "entropy": 0.88,
            "reverse_match": 0.9,
            "mass_ppm": 2.1,
            "final_score": 0.76,
            "rank_score": 0.76,
            "confidence": "High",
            "isomer_ambiguous": True,
            "isopair_evidence": {
                "evaluated": True,
                "runtime_influence": False,
                "ranking_changed": False,
                "counterfactual_ranking_changed": True,
                "production_abstention_reason": "shadow_only_after_mix33_validation_regressed",
                "adduct_context": {
                    "status": "resolved",
                    "formula_conditioned": True,
                    "resolved_adduct": "[M-H]-",
                    "mass_error_ppm": 2.1,
                },
            },
            "evidence_governance": {
                "policy_id": "evigate-ms-locked-v1-2026-08-02",
                "claim_status": "isomer_group_only",
                "claim_scope": "molecular_formula_isomer_group",
                "release_eligible_before_replication": False,
                "release_blockers": ["competitive_same_formula_isomer"],
                "ranking_changed": False,
                "coverage_guarantee": False,
                "nodes": [
                    {
                        "key": "spectral_similarity",
                        "label": "Library MS/MS similarity",
                        "status": "supporting",
                        "direction": "supports_candidate",
                        "applicability": "applicable",
                        "reason": "Correlated spectral metrics form one family.",
                        "value": {"cosine": 0.91},
                        "provenance": {"library_spectrum_id": "LIB1"},
                        "independence_group": "query_ms2_vs_library_spectrum",
                        "independent_measurement": False,
                        "used_for": ["ranking", "confidence"],
                        "warnings": [],
                    }
                ],
            },
        },
    }


def test_legacy_adapter_preserves_identity_scores_and_unavailable_calibration():
    record = legacy_row_to_identification(
        "task1", _legacy_row(), runtime_snapshot_id="v9-test"
    )
    assert record is not None
    assert record.candidate.inchikey == "AAAAAAAAAAAAAA-UHFFFAOYSA-N"
    assert record.scores.rank_score == 0.76
    assert record.scores.legacy_confidence == "High"
    assert record.decision.scope.value == "isomer_group"
    assert record.decision.calibration.status == CalibrationStatus.UNAVAILABLE
    assert record.decision.calibration.q_connectivity is None
    assert record.module_states["recap_formula_adduct"].value == "shadow"


def test_formula_conditioned_evidence_is_forced_to_shadow_only():
    record = legacy_row_to_identification("task1", _legacy_row())
    recap = next(
        node for node in record.evidence if node.node_id == "recap_formula_conditioned_adduct"
    )
    assert recap.candidate_conditioned
    assert recap.runtime_influence is False
    assert recap.used_for == ("shadow_audit",)

    with pytest.raises(ValidationError):
        EvidenceNode(
            node_id="bad",
            label="Bad circular evidence",
            status=EvidenceStatus.SUPPORTING,
            direction="supports_candidate",
            applicability=EvidenceApplicability.APPLICABLE,
            reason="candidate-conditioned",
            provenance=EvidenceProvenance(),
            independence_group="candidate_context",
            used_for=("ranking",),
            runtime_influence=True,
            candidate_conditioned=True,
        )


def test_uncalibrated_vector_rejects_qvalues():
    with pytest.raises(ValidationError):
        CalibrationVector(status=CalibrationStatus.UNAVAILABLE, q_connectivity=0.01)


def test_contract_roundtrip_is_additive_and_paginated(tmp_path: Path):
    path = tmp_path / "task1_identification_records.json"
    rows = [_legacy_row(), {"spectrum_id": "unmatched", "match": None}, _legacy_row()]
    rows[2] = json.loads(json.dumps(rows[2]))
    rows[2]["spectrum_id"] = "scan=43"
    metadata = write_identification_contract(
        "task1", rows, path, runtime_snapshot_id="v9-test"
    )
    page = read_identification_contract(path, offset=1, limit=1)

    assert metadata["record_count"] == 2
    assert metadata["production_ranking_changed_records"] == 0
    assert metadata["expected_compounds_counted"] is False
    assert page["metadata"]["record_count"] == 2
    assert len(page["records"]) == 1
    assert page["records"][0]["observation"]["spectrum_id"] == "scan=43"
    assert page["has_more"] is False


def test_historical_v9_row_migrates_without_changing_selected_identity():
    root = Path(__file__).resolve().parents[1]
    path = root / "output" / "tasks" / "4515e871f701_results.json"
    if not path.exists():
        pytest.skip("V9 golden task artifact is not available")
    import ijson

    with path.open("rb") as handle:
        row = next(
            item
            for item in ijson.items(handle, "item", use_float=True)
            if isinstance(item.get("match"), dict)
        )
    record = legacy_row_to_identification("4515e871f701", row)
    assert record.candidate.inchikey == row["match"].get("inchikey")
    assert record.scores.rank_score == row["match"].get("rank_score")
    assert record.scores.ranking_changed is False
    assert record.decision.ranking_changed is False
