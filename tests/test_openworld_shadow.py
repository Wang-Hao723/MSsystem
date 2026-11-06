from __future__ import annotations

import copy
import csv
import json
import subprocess
from pathlib import Path

from matching.openworld_shadow import (
    FiddleRuntime,
    FioraRuntime,
    build_fiddle_formula_shadow,
    build_fiora_forward_shadow,
)


def _candidate(name: str, inchikey: str, smiles: str, top1: bool = False) -> dict:
    return {
        "compound_name": name,
        "inchikey": inchikey,
        "smiles": smiles,
        "formula": "C2H6O",
        "adduct": "[M+H]+",
        "rank_score": 0.8 if top1 else 0.7,
    }


def test_missing_runtime_fails_closed(tmp_path):
    payload = build_fiora_forward_shadow(
        [],
        task_id="task",
        runtime=FioraRuntime(
            executable=tmp_path / "missing-executable",
            model_path=tmp_path / "missing-model.pt",
        ),
    )

    assert payload["status"] == "unavailable"
    assert payload["runtime_influence"] is False
    assert payload["ranking_changed"] is False


def test_fiora_shadow_can_disagree_but_never_mutates_production(monkeypatch, tmp_path):
    executable = tmp_path / "fiora-predict"
    executable.write_text("stub", encoding="utf-8")
    executable.chmod(0o755)
    model = tmp_path / "model.pt"
    model.write_bytes(b"model")
    top1 = _candidate(
        "production A",
        "AAAAAAAAAAAAAA-BBBBBBBBBB-C",
        "CCO",
        top1=True,
    )
    alternative = _candidate(
        "alternative B",
        "CCCCCCCCCCCCCC-DDDDDDDDDD-E",
        "COC",
    )
    top1["top_candidates"] = [dict(top1), alternative]
    results = [{
        "spectrum_id": "scan=1",
        "precursor_mz": 100.0,
        "collision_energy": 20.0,
        "instrument_type": "QTOF",
        "raw_peaks": [[60.0, 100.0], [80.0, 50.0]],
        "match": top1,
    }]
    before = copy.deepcopy(results)

    def fake_run(command, **kwargs):
        input_path = Path(command[command.index("-i") + 1])
        output_path = Path(command[command.index("-o") + 1])
        with input_path.open(encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle))
        blocks = []
        for row in rows:
            # Candidate 0 intentionally disagrees; candidate 1 matches query.
            peaks = "30 100\n40 50" if row["Name"].endswith("c0") else "60 100\n80 50"
            blocks.append(
                "BEGIN IONS\n"
                f"TITLE={row['Name']}\n"
                "PEPMASS=100\n"
                f"{peaks}\n"
                "END IONS\n"
            )
        output_path.write_text("".join(blocks), encoding="utf-8")
        return subprocess.CompletedProcess(command, 0, stdout="ok", stderr="")

    monkeypatch.setattr(subprocess, "run", fake_run)
    payload = build_fiora_forward_shadow(
        results,
        task_id="task",
        runtime=FioraRuntime(executable=executable, model_path=model),
    )

    assert payload["status"] == "completed_shadow"
    assert payload["summary"]["counterfactual_alternative_winner"] == 1
    assert payload["summary"]["high_promotions"] == 0
    assert payload["ranking_changed"] is False
    assert payload["spectra"][0]["production_top1_preserved"] is False
    assert results == before


def _fiddle_runtime(tmp_path: Path) -> FiddleRuntime:
    python = tmp_path / "python"
    python.write_text("stub", encoding="utf-8")
    python.chmod(0o755)
    worker = tmp_path / "worker.py"
    worker.write_text("stub", encoding="utf-8")
    model = tmp_path / "formula.pt"
    model.write_bytes(b"formula")
    rescore = tmp_path / "rescore.pt"
    rescore.write_bytes(b"rescore")
    config = tmp_path / "config.yml"
    config.write_text("encoding: {}", encoding="utf-8")
    return FiddleRuntime(
        python_executable=python,
        worker_path=worker,
        model_path=model,
        rescore_model_path=rescore,
        config_path=config,
    )


def test_missing_fiddle_runtime_fails_closed(tmp_path):
    payload = build_fiddle_formula_shadow(
        [],
        task_id="task",
        runtime=FiddleRuntime(
            python_executable=tmp_path / "missing-python",
            worker_path=tmp_path / "missing-worker.py",
            model_path=tmp_path / "missing-model.pt",
            rescore_model_path=tmp_path / "missing-rescore.pt",
            config_path=tmp_path / "missing-config.yml",
        ),
    )

    assert payload["status"] == "unavailable"
    assert payload["runtime_influence"] is False
    assert payload["confidence_changed"] is False


def test_fiddle_adds_unmatched_formula_candidates_without_high_promotion(
    monkeypatch, tmp_path
):
    matched = _candidate(
        "production A",
        "AAAAAAAAAAAAAA-BBBBBBBBBB-C",
        "CCO",
        top1=True,
    )
    matched["confidence"] = "High"
    results = [
        {
            "spectrum_id": "scan=1",
            "precursor_mz": 100.0,
            "collision_energy": 20.0,
            "instrument_type": "QTOF",
            "ion_mode": "POS",
            "raw_peaks": [
                [50.0, 10.0],
                [51.0, 20.0],
                [52.0, 30.0],
                [53.0, 40.0],
                [54.0, 50.0],
            ],
            "match": matched,
        },
        {
            "spectrum_id": "scan=2",
            "precursor_mz": 91.0,
            "collision_energy": 30.0,
            "instrument_type": "QTOF",
            "ion_mode": "NEG",
            "raw_peaks": [
                [50.0, 10.0],
                [60.0, 20.0],
                [70.0, 30.0],
                [80.0, 40.0],
                [90.0, 50.0],
            ],
            "match": None,
        },
    ]
    before = copy.deepcopy(results)
    captured_requests = []

    def fake_run(command, **kwargs):
        input_path = Path(command[command.index("--input") + 1])
        output_path = Path(command[command.index("--output") + 1])
        request = json.loads(input_path.read_text(encoding="utf-8"))
        captured_requests.extend(request["spectra"])
        predictions = []
        for spectrum in request["spectra"]:
            formula = "C2H6O" if spectrum["id"] == "f0" else "C3H8O"
            predictions.append({
                "id": spectrum["id"],
                "candidates": [
                    {"formula": formula, "score": 0.8, "mass": 46.04}
                ],
                "metadata": {"instrument_type": "qtof"},
            })
        output_path.write_text(
            json.dumps({"predictions": predictions}), encoding="utf-8"
        )
        return subprocess.CompletedProcess(command, 0, stdout="ok", stderr="")

    monkeypatch.setattr(subprocess, "run", fake_run)
    payload = build_fiddle_formula_shadow(
        results,
        task_id="task",
        runtime=_fiddle_runtime(tmp_path),
        precursor_correction={
            "status": "qualified_shadow",
            "module_id": "candidate_independent_lockmass_recalibration",
            "candidate_independent": True,
            "chemical_identity_used_for_fit": False,
            "correction": {"estimated_error_ppm": 10.0},
        },
    )

    assert payload["status"] == "completed_shadow"
    assert payload["summary"]["unmatched_spectra_with_formula_candidates"] == 1
    assert payload["summary"]["library_formula_top5_recalled"] == 1
    assert payload["summary"]["high_promotions"] == 0
    assert payload["summary"]["precursor_correction_applied_spectra"] == 2
    assert captured_requests[0]["precursor_mz"] < results[0]["precursor_mz"]
    assert payload["ranking_changed"] is False
    assert payload["spectra"][1]["adduct_evidence"]["hypothesis"] is True
    assert results == before
