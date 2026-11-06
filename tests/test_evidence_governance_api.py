from fastapi.testclient import TestClient

from web import server


def _candidate():
    return {
        "compound_name": "Candidate A",
        "inchikey": "AAAAAAAAAAAAAA-BBBBBBBBBB-C",
        "smiles": "CC1=CC=CC=C1O",
        "formula": "C10H12O2",
        "adduct": "[M+H]+",
        "cosine": 0.94,
        "entropy": 0.90,
        "reverse_match": 0.92,
        "final_score": 0.74,
        "rank_score": 0.74,
        "expert_score": 0.55,
        "expert_details": {
            "dimension_assessments": {
                "neutral_loss": {"status": "scored"},
            },
        },
        "mass_ppm": 1.5,
        "mass_metadata_consistent": True,
        "source_tier": "experimental",
        "has_structure": True,
        "quality_flags": [],
        "confidence": "High",
        "isomer_ambiguous": False,
        "top_candidates": [],
    }


def _task_rows(rt):
    rows = []
    for scan, energy in enumerate((20.0, 40.0, 60.0), start=1):
        rows.append({
            "spectrum_id": str(scan),
            "precursor_mz": 165.091 + (rt - 5.0) * 0.01,
            "retention_time": rt,
            "ion_mode": "POS",
            "collision_energy": energy,
            "precursor_feature_id": "PF000001",
            "precursor_feature_spectra": 3,
            "precursor_feature_collision_energies": [20.0, 40.0, 60.0],
            "precursor_feature_identity_consensus": 1.0,
            "precursor_feature_identity_consistent": True,
            "blank_status": "not_requested",
            "blank_flagged": False,
            "match": _candidate(),
        })
    return rows


def test_policy_and_passport_endpoints_expose_non_reranking_contract(monkeypatch):
    rows = _task_rows(5.0)
    monkeypatch.setattr(server, "_get_task_results", lambda task_id: rows)
    client = TestClient(server.app)

    policy = client.get("/api/evidence-governance/policy")
    assert policy.status_code == 200
    assert policy.json()["module"] == "EviGate-MS"
    assert policy.json()["invariants"]["ranking_changed"] is False
    assert policy.json()["invariants"][
        "stepped_energy_scans_are_independent"
    ] is False

    passport = client.get("/api/tasks/task-a/evidence-passport/0")
    assert passport.status_code == 200
    payload = passport.json()
    assert payload["candidate"]["compound_name"] == "Candidate A"
    assert payload["passport"]["ranking_changed"] is False
    assert payload["passport"]["claim_status"] == (
        "provisional_needs_independent_injection"
    )

    ledger = client.get(
        "/api/tasks/task-a/evidence-governance?claim_status="
        "provisional_needs_independent_injection&limit=1"
    )
    assert ledger.status_code == 200
    assert ledger.json()["total"] == 3
    assert len(ledger.json()["ledger"]) == 1
    assert ledger.json()["ledger"][0]["passport_url"].endswith(
        "/evidence-passport/0"
    )


def test_cohort_endpoint_validates_tasks_and_returns_paginated_consensus(
    monkeypatch,
):
    task_rows = {"inj-a": _task_rows(5.0), "inj-b": _task_rows(5.02)}
    monkeypatch.setattr(server, "_tasks", {
        "inj-a": {"id": "inj-a", "status": "completed", "filename": "A.mzML"},
        "inj-b": {"id": "inj-b", "status": "completed", "filename": "B.mzML"},
    })
    monkeypatch.setattr(
        server, "_get_task_results", lambda task_id: task_rows[task_id]
    )
    client = TestClient(server.app)

    response = client.post("/api/evidence-governance/cohort", json={
        "task_ids": ["inj-a", "inj-b"],
        "independence_groups": {"inj-a": "run-a", "inj-b": "run-b"},
        "limit": 1,
    })

    assert response.status_code == 200
    payload = response.json()
    assert payload["policy_id"] == server.EVIGATE_POLICY_ID
    assert payload["summary"]["released_unique_identities"] == 1
    assert payload["clusters_total"] == 1
    assert payload["clusters"][0]["action"] == (
        "release_replicate_stable_candidate"
    )
    assert payload["clusters"][0]["ranking_changed"] is False


def test_cohort_endpoint_rejects_unfinished_task(monkeypatch):
    monkeypatch.setattr(server, "_tasks", {
        "done": {"id": "done", "status": "completed", "filename": "done.mzML"},
        "busy": {"id": "busy", "status": "running", "filename": "busy.mzML"},
    })
    client = TestClient(server.app)

    response = client.post("/api/evidence-governance/cohort", json={
        "task_ids": ["done", "busy"],
    })

    assert response.status_code == 409
    assert "requires completed tasks" in response.json()["detail"]
