import json
from pathlib import Path

from acquisition.action_space import candidate_actions
from acquisition.ambiguity_graph import build_acquisition_needs
from acquisition.information_gain import score_action
from acquisition.observable_catalog import build_observable_catalog


def _row():
    return {
        "spectrum_id": "scan-1",
        "precursor_mz": 250.1,
        "retention_time": 2.0,
        "ion_mode": "POS",
        "instrument_type": "tof",
        "instrument": "Agilent Q-TOF",
        "collision_energy": 15.0,
        "collision_energy_unit": "eV",
        "precursor_feature_collision_energies": [15.0],
        "query_peaks": [[100.0, 100.0], [120.0, 40.0]],
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
            "chimeric_deconvolution": {"isolation_status": "mixed"},
        },
    }


def test_maximin_action_utility_is_bounded_and_pair_specific():
    row = _row()
    need = build_acquisition_needs("task", [row])[0]
    policy = json.loads(Path("config/acquisition_policy.json").read_text())
    proposal = next(
        item for item in candidate_actions(need, row, policy)
        if item["action_kind"] == "targeted_ms2_multi_ce"
    )
    action = score_action(
        need, row, proposal, build_observable_catalog(row), policy
    )

    assert len(action.pairwise_separation) == 1
    assert action.minimum_pair_separation > 0.2
    assert action.diagnostic_ions == (150.0, 175.0)
    assert action.collision_energies == (10.0, 20.0, 40.0)
    assert action.collision_energy_design is not None
    assert action.collision_energy_design.evidence_status == (
        "single_energy_only_no_response_curve"
    )
    assert 0 <= action.utility_interval[0] <= action.utility <= action.utility_interval[1] <= 1
    assert "不是身份概率" in "".join(action.rationale)
