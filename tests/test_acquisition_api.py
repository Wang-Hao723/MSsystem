from fastapi.testclient import TestClient

from web import server


def _row():
    return {
        "spectrum_id": "scan-1",
        "precursor_mz": 250.1,
        "retention_time": 2.0,
        "ion_mode": "POS",
        "instrument_type": "tof",
        "instrument": "Agilent instrument model",
        "collision_energy": 15.0,
        "collision_energy_unit": "eV",
        "precursor_feature_collision_energies": [15.0],
        "lcms_feature_id": "feature-1",
        "query_peaks": [[100.0, 100.0]],
        "match": {
            "compound_name": "A",
            "inchikey": "AAAAAAAAAAAAAA-BBBBBBBBBB-C",
            "formula": "C12H12O6",
            "rank_score": 0.8,
            "lib_peaks": [[100.0, 100.0], [150.0, 30.0]],
            "top_candidates": [{
                "compound_name": "B",
                "inchikey": "CCCCCCCCCCCCCC-DDDDDDDDDD-E",
                "formula": "C12H12O6",
                "rank_score": 0.77,
                "lib_peaks": [[100.0, 100.0], [175.0, 35.0]],
            }],
        },
    }


def test_acquisition_plan_api_and_advisory_export(monkeypatch, tmp_path):
    task = {
        "id": "api-plan",
        "filename": "sample.mzML",
        "status": "completed",
        "progress": 100,
        "vendor_conversion": {
            "format": "Agilent MassHunter .d to mzML",
            "source": {"acquisition": {"instrument": "Q-TOF"}},
        },
    }
    monkeypatch.setattr(server, "_tasks", {"api-plan": task})
    monkeypatch.setattr(server, "_get_task_results", lambda task_id: [_row()])
    monkeypatch.setattr(server, "TASK_DIR", tmp_path)
    client = TestClient(server.app)

    created = client.post(
        "/api/tasks/api-plan/acquisition-plans",
        json={"max_plans": 5, "max_actions_per_plan": 4, "refresh": False},
    )

    assert created.status_code == 200
    payload = created.json()
    assert payload["metadata"]["plan_count"] == 1
    action_kinds = {
        action["action_kind"]
        for action in payload["plans"][0]["ranked_actions"]
    }
    assert "authentic_standard_rt" not in action_kinds
    assert "ion_mobility_ccs" not in action_kinds
    assert payload["metadata"]["resource_gate"][
        "sample_only_primary_workflow"
    ] is True
    plan_id = payload["plans"][0]["plan_id"]

    listed = client.get("/api/tasks/api-plan/acquisition-plans")
    assert listed.status_code == 200
    assert listed.json()["total"] == 1

    exported = client.get(
        f"/api/tasks/api-plan/acquisition-plans/{plan_id}/export?format=thermo-qe"
    )
    assert exported.status_code == 200
    assert "not an instrument method file" in exported.text
    assert "targeted_ms2_multi_ce" in exported.text
    assert "agilent_qtof_ev_sequential_v1.0.0" in exported.text
    assert "resistant_high_energy_extension=60.0" in exported.text
    assert ",False,False,True,False,True" in exported.text
