import copy
import json
from pathlib import Path

import pytest

from acquisition.action_space import candidate_actions
from acquisition.ambiguity_graph import build_acquisition_needs
from acquisition.planner import generate_acquisition_plans, read_plan_store
from acquisition.resources import resolve_acquisition_resources


IDENTITY_A = "AAAAAAAAAAAAAA-BBBBBBBBBB-C"
IDENTITY_B = "CCCCCCCCCCCCCC-DDDDDDDDDD-E"


def _policy():
    return json.loads(Path("config/acquisition_policy.json").read_text())


def _row(*, spectrum_id="scan-1", feature_id="feature-1"):
    return {
        "spectrum_id": spectrum_id,
        "precursor_mz": 250.1,
        "retention_time": 2.0,
        "ion_mode": "POS",
        "instrument_type": "tof",
        "instrument": "Agilent Q-TOF",
        "collision_energy": 15.0,
        "collision_energy_unit": "eV",
        "precursor_feature_collision_energies": [15.0],
        "lcms_feature_id": feature_id,
        "query_peaks": [[100.0, 100.0], [120.0, 40.0]],
        "match": {
            "compound_name": "A",
            "inchikey": IDENTITY_A,
            "formula": "C12H12O6",
            "rank_score": 0.8,
            "lib_peaks": [[100.0, 100.0], [150.0, 30.0]],
            "top_candidates": [{
                "compound_name": "B",
                "inchikey": IDENTITY_B,
                "formula": "C12H12O6",
                "rank_score": 0.77,
                "lib_peaks": [[100.0, 100.0], [175.0, 35.0]],
            }],
        },
    }


def test_default_action_space_omits_unavailable_standard_and_ccs():
    row = _row()
    need = build_acquisition_needs("task", [row])[0]

    kinds = {
        action["action_kind"]
        for action in candidate_actions(need, row, _policy())
    }

    assert "targeted_ms2_multi_ce" in kinds
    assert "authentic_standard_rt" not in kinds
    assert "ion_mobility_ccs" not in kinds


def test_candidate_specific_standard_enters_only_with_inventory_and_budget():
    row = _row()
    need = build_acquisition_needs("task", [row])[0]
    context = {
        "acquisition_resources": {
            "declaration_id": "cabinet-20260804",
            "authentic_standard_identity_keys": [IDENTITY_A],
            "authentic_standard_action_budget": 1,
        }
    }

    actions = candidate_actions(need, row, _policy(), context=context)
    standard = next(
        action for action in actions
        if action["action_kind"] == "authentic_standard_rt"
    )

    assert standard["resource_status"] == "explicitly_available"
    assert standard["resource_identity_keys"] == (f"ik::{IDENTITY_A}",)
    assert standard["resource_declaration_id"] == "cabinet-20260804"
    assert "ion_mobility_ccs" not in {
        action["action_kind"] for action in actions
    }


def test_ccs_requires_instrument_reference_and_positive_budget():
    row = _row()
    need = build_acquisition_needs("task", [row])[0]
    context = {
        "acquisition_resources": {
            "ion_mobility_available": True,
            "ccs_reference_identity_keys": [IDENTITY_A, IDENTITY_B],
            "ion_mobility_action_budget": 1,
        }
    }

    actions = candidate_actions(need, row, _policy(), context=context)
    ccs = next(
        action for action in actions if action["action_kind"] == "ion_mobility_ccs"
    )

    assert ccs["resource_status"] == "explicitly_available"
    assert ccs["required_resources"] == (
        "ion_mobility_instrument", "ccs_reference"
    )

    with pytest.raises(ValueError, match="instrument availability"):
        resolve_acquisition_resources({
            "acquisition_resources": {
                "ccs_reference_identity_keys": [IDENTITY_A],
                "ion_mobility_action_budget": 1,
            }
        }, _policy())


def test_global_standard_budget_caps_number_of_plans(tmp_path):
    rows = [
        _row(spectrum_id="scan-1", feature_id="feature-1"),
        _row(spectrum_id="scan-2", feature_id="feature-2"),
    ]
    rows[1] = copy.deepcopy(rows[1])
    rows[1]["precursor_mz"] = 251.1
    rows[1]["retention_time"] = 2.5
    task = {
        "id": "budgeted-task",
        "filename": "sample.mzML",
        "acquisition_resources": {
            "declaration_id": "cabinet-20260804",
            "authentic_standard_identity_keys": [IDENTITY_A],
            "authentic_standard_action_budget": 1,
        },
    }

    store = generate_acquisition_plans(
        task,
        rows,
        task_dir=tmp_path,
        project_root=Path.cwd(),
        runtime_snapshot_id="runtime-test",
        policy=_policy(),
        max_plans=10,
        max_actions_per_plan=6,
    )

    standard_actions = [
        action
        for plan in store["plans"]
        for action in plan["ranked_actions"]
        if action["action_kind"] == "authentic_standard_rt"
    ]
    assert len(standard_actions) == 1
    assert store["metadata"]["resource_gate"][
        "authentic_standard_action_budget"
    ] == 1
    assert store["metadata"]["resource_gate"][
        "authentic_standard_actions_planned"
    ] == 1


def test_legacy_unverified_resource_actions_are_removed_on_read(tmp_path):
    path = tmp_path / "legacy.json"
    path.write_text(json.dumps({
        "schema_version": "1.1.0",
        "metadata": {"plan_count": 1, "planned_actions": 2},
        "plans": [{
            "plan_id": "legacy-plan",
            "ranked_actions": [
                {"action_kind": "authentic_standard_rt"},
                {"action_kind": "targeted_ms2_multi_ce"},
            ],
        }],
        "followups": [],
    }), encoding="utf-8")

    store = read_plan_store(path)

    assert [
        action["action_kind"]
        for action in store["plans"][0]["ranked_actions"]
    ] == ["targeted_ms2_multi_ce"]
    assert store["metadata"]["planned_actions"] == 1
    assert store["metadata"]["resource_gate"][
        "legacy_unverified_resource_actions_omitted"
    ] == 1
