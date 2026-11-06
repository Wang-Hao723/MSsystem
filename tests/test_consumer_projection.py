from __future__ import annotations

from governance.claim_compiler import default_task_claim_context
from pipeline.claim_contract import write_claim_contract
from pipeline.consumer_projection import attach_claim_decisions
from schemas.legacy import legacy_row_to_identification


def test_claim_projection_joins_by_spectrum_and_identity_without_disk_mutation(tmp_path):
    row = {
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
    }
    record = legacy_row_to_identification("task-1", row)
    assert record is not None
    path = tmp_path / "claims.json"
    write_claim_contract(
        "task-1", [record], path, context=default_task_claim_context("task-1")
    )
    before = path.read_bytes()

    summary = attach_claim_decisions([row], path)

    assert summary["attached"] == 1
    assert row["match"]["final_claim_decision"]["released"] is False
    assert row["match"]["final_claim_decision"]["status"] == (
        "blocked_project_calibration_unavailable"
    )
    assert row["match"]["confidence"] == "High"
    assert path.read_bytes() == before
