from __future__ import annotations

import copy
import csv
import json
import subprocess
from pathlib import Path

from fastapi.testclient import TestClient

from matching.evidiscrim import (
    build_evidiscrim_shadow,
    load_evidiscrim_policy,
    results_sha256,
    stream_select_evidiscrim_rows,
)
from matching.openworld_shadow import (
    FIORA_ALLOWED_ADDUCTS,
    FioraRuntime,
    build_fiora_evidiscrim_shadow,
)
from acquisition.planner import generate_acquisition_plans
from web import server


ROOT = Path(__file__).resolve().parents[1]
POLICY = load_evidiscrim_policy(ROOT / "config" / "evidiscrim_policy.json")


def _candidate(
    name: str,
    inchikey: str,
    smiles: str,
    *,
    formula: str = "C2H6O",
    score: float = 0.8,
) -> dict:
    return {
        "compound_name": name,
        "inchikey": inchikey,
        "smiles": smiles,
        "formula": formula,
        "adduct": "[M+H]+",
        "rank_score": score,
    }


def _two_candidate_result(*, instrument: str = "HCD") -> tuple[list[dict], dict, dict]:
    ethanol = _candidate(
        "ethanol",
        "LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
        "CCO",
        score=0.8,
    )
    ether = _candidate(
        "dimethyl ether",
        "LCGLNKUTAGEVQW-UHFFFAOYSA-N",
        "COC",
        score=0.7,
    )
    ethanol["top_candidates"] = [dict(ethanol), ether]
    results = [{
        "spectrum_id": "scan=1",
        "precursor_mz": 47.0,
        "collision_energy": 20.0,
        "collision_energy_unit": "eV",
        "instrument_type": instrument,
        "raw_peaks": [[31.0, 100.0], [29.0, 10.0]],
        "match": ethanol,
    }]
    return results, ethanol, ether


def _record(
    identity: str,
    energy: float,
    peaks: list[list[float]],
    *,
    replicates: list[list[list[float]]] | None = None,
    qualified: bool = True,
) -> dict:
    record = {
        "row_index": 0,
        "identity_key": f"ik::{identity}",
        "collision_energy": energy,
        "collision_energy_unit": "eV",
        "predicted_peaks": peaks,
        "energy_transfer_qualified": qualified,
        "provider_id": "synthetic-model-a",
        "architecture_family": "synthetic-architecture-a",
    }
    if replicates is not None:
        record["prediction_members"] = [
            {
                "provider_id": f"synthetic-model-{index + 1}",
                "architecture_family": f"synthetic-architecture-{index + 1}",
                "predicted_peaks": member_peaks,
                "energy_transfer_qualified": qualified,
            }
            for index, member_peaks in enumerate(replicates)
        ]
    else:
        record["prediction_members"] = [
            {
                "provider_id": "synthetic-model-a",
                "architecture_family": "synthetic-architecture-a",
                "predicted_peaks": peaks,
                "energy_transfer_qualified": qualified,
            },
            {
                "provider_id": "synthetic-model-b",
                "architecture_family": "synthetic-architecture-b",
                "predicted_peaks": peaks,
                "energy_transfer_qualified": qualified,
            },
        ]
    return record


def test_candidate_diagnostic_fit_resolves_pair_without_mutating_production():
    results, ethanol, ether = _two_candidate_result()
    before = copy.deepcopy(results)
    records = [
        _record(ethanol["inchikey"], 20.0, [[31.0, 1.0], [29.0, 0.1]]),
        _record(ether["inchikey"], 20.0, [[45.0, 1.0], [29.0, 0.1]]),
    ]

    artifact = build_evidiscrim_shadow(
        results,
        task_id="task",
        prediction_records=records,
        policy=POLICY,
    )

    spectrum = artifact["spectra"][0]
    pair = spectrum["pairwise_comparisons"][0]
    assert pair["status"] == "resolved_counterfactual"
    assert pair["counterfactual_winner_identity_key"] == (
        f"ik::{ethanol['inchikey']}"
    )
    assert {ion["mz"] for ion in pair["diagnostic_ions"]} == {31.0, 45.0}
    absent = next(ion for ion in pair["diagnostic_ions"] if ion["mz"] == 45.0)
    assert absent["observed_present"] is False
    assert absent["absence_penalty_weight"] == 0.35
    assert artifact["ranking_changed"] is False
    assert artifact["confidence_changed"] is False
    assert artifact["high_promotions"] == 0
    assert artifact["source_results_sha256"] == artifact["source_results_sha256_after"]
    assert results == before


def test_energy_selection_uses_uncertainty_adjusted_maximin_separation():
    results, ethanol, ether = _two_candidate_result()
    records = [
        _record(ethanol["inchikey"], 10.0, [[31.0, 1.0], [45.0, 0.5]]),
        _record(ether["inchikey"], 10.0, [[31.0, 1.0], [45.0, 0.5]]),
        _record(ethanol["inchikey"], 20.0, [[31.0, 1.0], [45.0, 0.1]]),
        _record(ether["inchikey"], 20.0, [[31.0, 0.1], [45.0, 1.0]]),
        _record(ethanol["inchikey"], 40.0, [[31.0, 1.0], [45.0, 0.5]]),
        _record(ether["inchikey"], 40.0, [[31.0, 0.5], [45.0, 1.0]]),
    ]

    artifact = build_evidiscrim_shadow(
        results,
        task_id="task",
        prediction_records=records,
        policy=POLICY,
    )

    selection = artifact["spectra"][0]["energy_selection"]
    assert selection["status"] == "recommended_shadow"
    assert selection["recommended_collision_energy"] == 20.0
    assert selection["objective"] == (
        "maximize_worst_pair_spectral_distribution_divergence"
    )
    assert selection["planning_only"] is True
    assert selection["cross_instrument_conversion_applied"] is False


def test_prediction_ensemble_variance_reduces_diagnostic_reliability():
    results, ethanol, ether = _two_candidate_result()
    low_variance = [
        _record(
            ethanol["inchikey"], 20.0, [],
            replicates=[
                [[31.0, 1.0], [45.0, 0.3]],
                [[31.0, 1.0], [45.0, 0.3]],
            ],
        ),
        _record(
            ether["inchikey"], 20.0, [],
            replicates=[
                [[31.0, 0.3], [45.0, 1.0]],
                [[31.0, 0.3], [45.0, 1.0]],
            ],
        ),
    ]
    high_variance = [
        _record(
            ethanol["inchikey"], 20.0, [],
            replicates=[
                [[31.0, 1.0], [45.0, 0.1]],
                [[31.0, 1.0], [45.0, 0.5]],
            ],
        ),
        _record(
            ether["inchikey"], 20.0, [],
            replicates=[
                [[31.0, 0.1], [45.0, 1.0]],
                [[31.0, 0.5], [45.0, 1.0]],
            ],
        ),
    ]

    low = build_evidiscrim_shadow(
        results,
        task_id="low",
        prediction_records=low_variance,
        policy=POLICY,
    )["spectra"][0]["pairwise_comparisons"][0]
    high = build_evidiscrim_shadow(
        results,
        task_id="high",
        prediction_records=high_variance,
        policy=POLICY,
    )["spectra"][0]["pairwise_comparisons"][0]

    assert high["total_diagnostic_weight"] < low["total_diagnostic_weight"]
    assert high["energy_separations"][0]["uncertainty_reliability_factor"] < (
        low["energy_separations"][0]["uncertainty_reliability_factor"]
    )


def test_same_connectivity_stereoisomers_are_explicitly_out_of_scope():
    left = _candidate(
        "stereo-left", "", "F[C@H](Cl)Br", formula="CHBrClF"
    )
    right = _candidate(
        "stereo-right", "", "F[C@@H](Cl)Br", formula="CHBrClF", score=0.7
    )
    left["top_candidates"] = [dict(left), right]
    results = [{
        "spectrum_id": "stereo",
        "precursor_mz": 165.0,
        "collision_energy": 20.0,
        "collision_energy_unit": "eV",
        "raw_peaks": [[80.0, 100.0], [100.0, 50.0]],
        "match": left,
    }]
    records = [
        {
            "row_index": 0,
            "identity_key": f"sm::{left['smiles']}",
            "collision_energy": 20.0,
            "collision_energy_unit": "eV",
            "predicted_peaks": [[80.0, 1.0]],
            "energy_transfer_qualified": True,
        },
        {
            "row_index": 0,
            "identity_key": f"sm::{right['smiles']}",
            "collision_energy": 20.0,
            "collision_energy_unit": "eV",
            "predicted_peaks": [[100.0, 1.0]],
            "energy_transfer_qualified": True,
        },
    ]

    artifact = build_evidiscrim_shadow(
        results,
        task_id="stereo",
        prediction_records=records,
        policy=POLICY,
    )

    pair = artifact["spectra"][0]["pairwise_comparisons"][0]
    assert pair["eligible"] is False
    assert pair["status"] == "ineligible"
    assert pair["abstention_reason"] == (
        "same_2d_connectivity_stereoisomer_out_of_scope"
    )
    assert pair["counterfactual_winner_identity_key"] is None


def test_invalid_structure_and_nondiagnostic_predictions_fail_closed():
    results, ethanol, ether = _two_candidate_result()
    results[0]["match"]["top_candidates"][-1]["smiles"] = "not-a-smiles"
    records = [
        _record(ethanol["inchikey"], 20.0, [[31.0, 1.0]]),
        _record(ether["inchikey"], 20.0, [[31.0, 1.0]]),
    ]
    invalid = build_evidiscrim_shadow(
        results,
        task_id="invalid",
        prediction_records=records,
        policy=POLICY,
    )["spectra"][0]["pairwise_comparisons"][0]
    assert invalid["eligible"] is False
    assert invalid["abstention_reason"] == "invalid_candidate_smiles"

    valid_results, ethanol, ether = _two_candidate_result()
    nondiagnostic = build_evidiscrim_shadow(
        valid_results,
        task_id="nondiagnostic",
        prediction_records=[
            _record(ethanol["inchikey"], 20.0, [[31.0, 1.0], [45.0, 0.5]]),
            _record(ether["inchikey"], 20.0, [[31.0, 1.0], [45.0, 0.5]]),
        ],
        policy=POLICY,
    )["spectra"][0]["pairwise_comparisons"][0]
    assert nondiagnostic["status"] == "abstained"
    assert nondiagnostic["abstention_reason"] == (
        "insufficient_candidate_diagnostic_ions"
    )
    nondiagnostic_selection = build_evidiscrim_shadow(
        valid_results,
        task_id="nondiagnostic-energy",
        prediction_records=[
            _record(ethanol["inchikey"], 20.0, [[31.0, 1.0], [45.0, 0.5]]),
            _record(ether["inchikey"], 20.0, [[31.0, 1.0], [45.0, 0.5]]),
        ],
        policy=POLICY,
    )["spectra"][0]["energy_selection"]
    assert nondiagnostic_selection["status"] == (
        "abstained_predicted_worst_pair_separation_below_threshold"
    )
    assert nondiagnostic_selection["recommended_collision_energy"] is None


def test_input_prediction_order_is_deterministic():
    results, ethanol, ether = _two_candidate_result()
    records = [
        _record(ethanol["inchikey"], 20.0, [[31.0, 1.0], [29.0, 0.1]]),
        _record(ether["inchikey"], 20.0, [[45.0, 1.0], [29.0, 0.1]]),
    ]
    forward = build_evidiscrim_shadow(
        results, task_id="task", prediction_records=records, policy=POLICY
    )
    reverse = build_evidiscrim_shadow(
        results,
        task_id="task",
        prediction_records=list(reversed(records)),
        policy=POLICY,
    )
    assert forward["spectra"] == reverse["spectra"]
    assert forward["summary"] == reverse["summary"]


def test_candidate_limit_override_does_not_create_unrequested_missing_pairs():
    first = _candidate("1-propanol", "", "CCCO", formula="C3H8O")
    second = _candidate("2-propanol", "", "CC(O)C", formula="C3H8O", score=0.7)
    third = _candidate("methoxyethane", "", "CCOC", formula="C3H8O", score=0.6)
    first["top_candidates"] = [dict(first), second, third]
    results = [{
        "spectrum_id": "limit",
        "precursor_mz": 61.0,
        "collision_energy": 20.0,
        "collision_energy_unit": "eV",
        "raw_peaks": [[31.0, 100.0], [43.0, 20.0]],
        "match": first,
    }]
    records = [
        {
            "row_index": 0,
            "identity_key": "sm::CCCO",
            "collision_energy": 20.0,
            "collision_energy_unit": "eV",
            "predicted_peaks": [[31.0, 1.0], [45.0, 0.1]],
            "energy_transfer_qualified": True,
            "prediction_members": [
                {
                    "provider_id": "synthetic-a",
                    "architecture_family": "synthetic-a",
                    "predicted_peaks": [[31.0, 1.0], [45.0, 0.1]],
                    "energy_transfer_qualified": True,
                },
                {
                    "provider_id": "synthetic-b",
                    "architecture_family": "synthetic-b",
                    "predicted_peaks": [[31.0, 1.0], [45.0, 0.1]],
                    "energy_transfer_qualified": True,
                },
            ],
        },
        {
            "row_index": 0,
            "identity_key": "sm::CC(O)C",
            "collision_energy": 20.0,
            "collision_energy_unit": "eV",
            "predicted_peaks": [[31.0, 0.1], [45.0, 1.0]],
            "energy_transfer_qualified": True,
            "prediction_members": [
                {
                    "provider_id": "synthetic-a",
                    "architecture_family": "synthetic-a",
                    "predicted_peaks": [[31.0, 0.1], [45.0, 1.0]],
                    "energy_transfer_qualified": True,
                },
                {
                    "provider_id": "synthetic-b",
                    "architecture_family": "synthetic-b",
                    "predicted_peaks": [[31.0, 0.1], [45.0, 1.0]],
                    "energy_transfer_qualified": True,
                },
            ],
        },
    ]

    artifact = build_evidiscrim_shadow(
        results,
        task_id="limit",
        prediction_records=records,
        policy=POLICY,
        max_candidates_per_spectrum=2,
    )

    spectrum = artifact["spectra"][0]
    assert artifact["summary"]["maximum_candidates_per_spectrum"] == 2
    assert spectrum["candidate_count"] == 2
    assert spectrum["eligible_pair_count"] == 1
    assert spectrum["energy_selection"]["status"] == "recommended_shadow"


def test_large_result_stream_digest_matches_json_load_and_preserves_source_index(
    tmp_path,
):
    results, _, _ = _two_candidate_result()
    unmatched = {
        "spectrum_id": "unmatched",
        "precursor_mz": 100.0,
        "collision_energy": 20.0,
        "raw_peaks": [[50.0, 1.0]],
        "match": None,
    }
    repeated = copy.deepcopy(results[0])
    repeated["spectrum_id"] = "repeat-same-candidate-pairs"
    payload = [unmatched, *results, repeated]
    path = tmp_path / "results.json"
    path.write_text(json.dumps(payload, ensure_ascii=False), encoding="utf-8")

    selected, digest, audit = stream_select_evidiscrim_rows(
        path,
        max_spectra=2,
        max_candidates_per_spectrum=2,
        supported_adducts=FIORA_ALLOWED_ADDUCTS,
        policy=POLICY,
    )

    assert digest == results_sha256(payload)
    assert audit["total_source_rows"] == 3
    assert audit["selected_source_result_indices"] == [1]
    assert audit["duplicate_candidate_pair_signatures_skipped"] == 1
    assert selected[0]["_evidiscrim_source_result_index"] == 1


def test_fiora_adapter_generates_candidate_energy_grid(monkeypatch, tmp_path):
    executable = tmp_path / "fiora-predict"
    executable.write_text("stub", encoding="utf-8")
    executable.chmod(0o755)
    model = tmp_path / "model.pt"
    model.write_bytes(b"model")
    results, _, _ = _two_candidate_result(instrument="Orbitrap HCD")
    results[0]["collision_energy_unit"] = "NCE"

    def fake_run(command, **kwargs):
        input_path = Path(command[command.index("-i") + 1])
        output_path = Path(command[command.index("-o") + 1])
        with input_path.open(encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle))
        blocks = []
        for row in rows:
            if row["SMILES"] == "CCO":
                peaks = "31 1.0\n45 0.1"
            else:
                peaks = "31 0.1\n45 1.0"
            blocks.append(
                "BEGIN IONS\n"
                f"TITLE={row['Name']}\n"
                "PEPMASS=47\n"
                f"{peaks}\n"
                "END IONS\n"
            )
        output_path.write_text("".join(blocks), encoding="utf-8")
        return subprocess.CompletedProcess(command, 0, stdout="ok", stderr="")

    monkeypatch.setattr(subprocess, "run", fake_run)
    artifact = build_fiora_evidiscrim_shadow(
        results,
        task_id="task",
        runtime=FioraRuntime(executable=executable, model_path=model),
        policy=POLICY,
        collision_energy_grid=(10.0, 20.0),
    )

    assert artifact["status"] == "completed_shadow"
    assert artifact["summary"]["candidate_energy_requests"] == 4
    assert artifact["summary"]["usable_prediction_records"] == 4
    selection = artifact["spectra"][0]["energy_selection"]
    assert selection["status"] == (
        "unavailable_insufficient_independent_architecture_families"
    )
    assert selection["recommended_collision_energy"] is None
    assert artifact["summary"]["insufficient_architecture_pairs"] > 0
    assert all(
        pair["common_architecture_family_count"] == 1
        for pair in artifact["spectra"][0]["pairwise_comparisons"]
    )
    assert artifact["runtime_influence"] is False


def test_fiora_domain_shift_blocks_actionable_energy_recommendation(
    monkeypatch, tmp_path
):
    executable = tmp_path / "fiora-predict"
    executable.write_text("stub", encoding="utf-8")
    executable.chmod(0o755)
    model = tmp_path / "model.pt"
    model.write_bytes(b"model")
    results, _, _ = _two_candidate_result(instrument="Agilent QTOF")

    def fake_run(command, **kwargs):
        input_path = Path(command[command.index("-i") + 1])
        output_path = Path(command[command.index("-o") + 1])
        with input_path.open(encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle))
        output_path.write_text(
            "".join(
                "BEGIN IONS\n"
                f"TITLE={row['Name']}\n"
                "PEPMASS=47\n"
                + ("31 1\n45 0.1\n" if row["SMILES"] == "CCO" else "31 0.1\n45 1\n")
                + "END IONS\n"
                for row in rows
            ),
            encoding="utf-8",
        )
        return subprocess.CompletedProcess(command, 0, stdout="ok", stderr="")

    monkeypatch.setattr(subprocess, "run", fake_run)
    artifact = build_fiora_evidiscrim_shadow(
        results,
        task_id="task",
        runtime=FioraRuntime(executable=executable, model_path=model),
        policy=POLICY,
        collision_energy_grid=(10.0, 20.0),
    )

    selection = artifact["spectra"][0]["energy_selection"]
    assert selection["status"] == (
        "unavailable_insufficient_independent_architecture_families"
    )
    assert selection["recommended_collision_energy"] is None
    assert artifact["summary"]["collision_energy_recommendations"] == 0
    assert all(
        separation["instrument_domain_qualified"] is False
        for pair in artifact["spectra"][0]["pairwise_comparisons"]
        for separation in pair["energy_separations"]
    )


def test_evidiscrim_api_materializes_and_reads_paginated_sidecar(
    monkeypatch, tmp_path
):
    import matching.openworld_shadow as openworld

    task_id = "evidiscrim-api-task"
    results, ethanol, ether = _two_candidate_result()
    records = [
        _record(ethanol["inchikey"], 20.0, [[31.0, 1.0], [45.0, 0.1]]),
        _record(ether["inchikey"], 20.0, [[31.0, 0.1], [45.0, 1.0]]),
    ]
    artifact = build_evidiscrim_shadow(
        results,
        task_id=task_id,
        prediction_records=records,
        policy=POLICY,
    )
    monkeypatch.setattr(server, "TASK_DIR", tmp_path)
    monkeypatch.setattr(server, "_get_task_results", lambda _: results)
    monkeypatch.setattr(
        openworld,
        "build_fiora_evidiscrim_shadow",
        lambda *args, **kwargs: copy.deepcopy(artifact),
    )
    with server._tasks_lock:
        server._tasks[task_id] = {"id": task_id, "status": "completed"}
    try:
        client = TestClient(server.app)
        created = client.post(
            f"/api/tasks/{task_id}/evidiscrim",
            json={"max_spectra": 1, "max_candidates_per_spectrum": 2},
        )
        assert created.status_code == 200
        assert created.json()["status"] == "completed_shadow"
        assert created.json()["runtime_influence"] is False

        fetched = client.get(
            f"/api/tasks/{task_id}/evidiscrim?limit=1&offset=0"
        )
        assert fetched.status_code == 200
        assert fetched.json()["total"] == 1
        assert len(fetched.json()["spectra"]) == 1
        assert fetched.json()["has_more"] is False
    finally:
        with server._tasks_lock:
            server._tasks.pop(task_id, None)


def test_eviacquire_consumes_only_current_qualified_shadow_context(tmp_path):
    results, ethanol, ether = _two_candidate_result(instrument="Agilent QTOF")
    results[0].update({
        "retention_time": 2.0,
        "ion_mode": "POS",
        "instrument": "Agilent instrument model",
        "precursor_feature_collision_energies": [20.0],
        "lcms_feature_id": "feature-1",
        "query_peaks": results[0]["raw_peaks"],
    })
    predictions = []
    for energy in (10.0, 20.0, 40.0):
        predictions.extend([
            _record(
                ethanol["inchikey"], energy,
                [[31.0, 1.0], [45.0, 0.1]], qualified=True,
            ),
            _record(
                ether["inchikey"], energy,
                [[31.0, 0.1], [45.0, 1.0]], qualified=True,
            ),
        ])
    evidiscrim = build_evidiscrim_shadow(
        results,
        task_id="task-plan",
        prediction_records=predictions,
        policy=POLICY,
    )
    acquisition_policy = json.loads(
        (ROOT / "config" / "acquisition_policy.json").read_text(encoding="utf-8")
    )
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
        project_root=ROOT,
        runtime_snapshot_id="runtime-test",
        policy=acquisition_policy,
        refresh=True,
        evidiscrim_artifact=evidiscrim,
    )

    context = store["metadata"]["evidiscrim_context"]
    assert context["status"] == "accepted_shadow_planning_context"
    assert context["planning_context_used"] is True
    action = next(
        item for item in store["plans"][0]["ranked_actions"]
        if item["action_kind"] == "targeted_ms2_multi_ce"
    )
    expected = max(
        item["spectral_distribution_divergence"]
        for item in evidiscrim["spectra"][0]["pairwise_comparisons"][0][
            "energy_separations"
        ]
    )
    assert max(action["pairwise_separation"].values()) == round(expected, 4)
    assert "EviDiscrim-MS" in "".join(action["rationale"])
    assert store["metadata"]["production_ranking_changed"] is False
    assert store["metadata"]["predicted_evidence_confirmatory"] is False

    projected_results = copy.deepcopy(results)
    projected_results[0]["presentation_only_annotation"] = {
        "status": "added_after_commit"
    }
    projected_store = generate_acquisition_plans(
        task,
        projected_results,
        task_dir=tmp_path / "projected",
        project_root=ROOT,
        runtime_snapshot_id="runtime-test",
        policy=acquisition_policy,
        refresh=True,
        evidiscrim_artifact=evidiscrim,
    )
    assert projected_store["metadata"]["evidiscrim_context"]["status"] == (
        "accepted_source_row_fingerprint_matched_shadow_context"
    )
    assert projected_store["metadata"]["evidiscrim_context"][
        "planning_context_used"
    ] is True

    stale_results = copy.deepcopy(results)
    stale_results[0]["precursor_mz"] = 47.01
    stale_store = generate_acquisition_plans(
        {**task, "id": "task-stale"},
        stale_results,
        task_dir=tmp_path,
        project_root=ROOT,
        runtime_snapshot_id="runtime-test",
        policy=acquisition_policy,
        refresh=True,
        evidiscrim_artifact=evidiscrim,
    )
    assert stale_store["metadata"]["evidiscrim_context"]["status"] in {
        "rejected_task_id_mismatch",
        "rejected_stale_source_results",
    }
    assert stale_store["metadata"]["evidiscrim_context"][
        "planning_context_used"
    ] is False

    legacy = copy.deepcopy(evidiscrim)
    legacy["schema_version"] = "1.0.0"
    legacy_store = generate_acquisition_plans(
        task,
        results,
        task_dir=tmp_path / "legacy",
        project_root=ROOT,
        runtime_snapshot_id="runtime-test",
        policy=acquisition_policy,
        refresh=True,
        evidiscrim_artifact=legacy,
    )
    assert legacy_store["metadata"]["evidiscrim_context"]["status"] == (
        "rejected_legacy_evidiscrim_uncertainty_contract"
    )
    assert legacy_store["metadata"]["evidiscrim_context"][
        "planning_context_used"
    ] is False
