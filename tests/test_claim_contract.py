from __future__ import annotations

from governance.claim_compiler import default_task_claim_context
from pipeline.claim_contract import read_claim_contract, write_claim_contract
from schemas.legacy import legacy_row_to_identification


def _record():
    record = legacy_row_to_identification(
        "task-1",
        {
            "spectrum_id": "scan=1",
            "precursor_mz": 301.0,
            "match": {
                "compound_name": "Candidate",
                "inchikey": "AAAAAAAAAAAAAA-UHFFFAOYSA-N",
                "rank_score": 0.9,
                "confidence": "High",
                "evidence_governance": {
                    "policy_id": "v1",
                    "claim_status": "provisional_needs_independent_injection",
                    "claim_scope": "connectivity_level_library_candidate",
                    "release_eligible_before_replication": True,
                    "nodes": [],
                },
            },
        },
    )
    assert record is not None
    return record


def test_claim_contract_is_compact_paginated_and_fail_closed(tmp_path):
    path = tmp_path / "task-1_claim_decisions.json"
    metadata = write_claim_contract(
        "task-1",
        [_record(), _record().model_copy(update={"record_id": "record-2"})],
        path,
        context=default_task_claim_context("task-1"),
    )
    page = read_claim_contract(path, offset=1, limit=1)
    assert metadata["record_count"] == 2
    assert metadata["released_count"] == 0
    assert metadata["production_ranking_changed_records"] == 0
    assert metadata["legacy_confidence_overwritten_records"] == 0
    assert page["decisions"][0]["record_id"] == "record-2"
    assert page["decisions"][0]["legacy_confidence"] == "High"
    assert page["decisions"][0]["decision"]["status"] == (
        "blocked_project_calibration_unavailable"
    )
    assert page["has_more"] is False
