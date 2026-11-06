import copy
import json
from pathlib import Path

from acquisition.planner import generate_acquisition_plans, plan_store_path


def _results():
    return [{
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
    }]


def test_planner_materializes_additive_shadow_artifact(tmp_path):
    results = _results()
    original = copy.deepcopy(results)
    policy = json.loads(Path("config/acquisition_policy.json").read_text())
    task = {
        "id": "task-plan",
        "filename": "sample.mzML",
        "source_raw_filename": "sample.d",
        "vendor_conversion": {
            "format": "Agilent MassHunter .d to mzML",
            "source": {
                "acquisition": {
                    "instrument": "Q-TOF",
                    "acquisition_mode": "AutoMS2",
                }
            },
        },
    }

    store = generate_acquisition_plans(
        task,
        results,
        task_dir=tmp_path,
        project_root=Path.cwd(),
        runtime_snapshot_id="runtime-test",
        policy=policy,
        max_plans=10,
        max_actions_per_plan=4,
    )

    assert plan_store_path(tmp_path, "task-plan").is_file()
    assert store["metadata"]["plan_count"] == 1
    assert store["plans"][0]["module_state"] == "shadow"
    assert store["plans"][0]["production_ranking_changed"] is False
    assert store["plans"][0]["predicted_evidence_confirmatory"] is False
    assert results == original
    action_kinds = {
        action["action_kind"]
        for action in store["plans"][0]["ranked_actions"]
    }
    assert "authentic_standard_rt" not in action_kinds
    assert "ion_mobility_ccs" not in action_kinds
    resource_gate = store["metadata"]["resource_gate"]
    assert resource_gate["mode"] == "sample_only_fail_closed"
    assert resource_gate["sample_only_primary_workflow"] is True
    assert resource_gate["authentic_standard_actions_planned"] == 0
    assert resource_gate["ion_mobility_actions_planned"] == 0
    energy_action = next(
        action for action in store["plans"][0]["ranked_actions"]
        if action["action_kind"] == "targeted_ms2_multi_ce"
    )
    assert energy_action["collision_energies"] == [10.0, 20.0, 40.0]
    assert energy_action["collision_energy_design"]["automatic_ceiling"] == 60.0
    assert energy_action["collision_energy_design"][
        "universal_complete_fragmentation_assumed"
    ] is False
    energy_summary = store["metadata"]["collision_energy_design"]
    assert energy_summary["initial_energy_schedules"] == [[10.0, 20.0, 40.0]]
    assert energy_summary["conditional_energy_schedules"][
        "resistant_high_energy_extension"
    ] == [[60.0]]
    assert energy_summary["single_energy_no_response_curve_plan_count"] == 1

    cached = generate_acquisition_plans(
        task,
        results,
        task_dir=tmp_path,
        project_root=Path.cwd(),
        runtime_snapshot_id="runtime-test",
        policy=policy,
    )
    assert cached["cache_hit"] is True
