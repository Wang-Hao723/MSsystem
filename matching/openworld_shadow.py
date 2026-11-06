"""Fail-closed adapters for open-world formula/structure evidence.

Only shadow artifacts are produced here.  Predicted spectra are
candidate-conditioned simulations, not independent measurements, and cannot
change the production library Top-1, confidence tier, or release decision.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

from matching.library_quality import compound_identity_key, normalize_formula_text
from matching.evidiscrim import build_evidiscrim_shadow
from matching.scores import cosine_score, entropy_score, reverse_match_ratio
from schemas.evidiscrim import EviDiscrimArtifact
from utils.mgf_io import parse_mgf_file
from utils.spectrum import process_peaks, resolve_peak_processing_config


SCHEMA_VERSION = "1.0.0"
FIORA_MODULE_ID = "forward_spectrum_fiora"
FIORA_ALLOWED_ADDUCTS = frozenset({"[M+H]+", "[M-H]-", "[M]+", "[M]-"})
FIDDLE_MODULE_ID = "formula_prediction_fiddle"
FIDDLE_SUPPORTED_ADDUCTS = frozenset({"[M+H]+", "[M-H]-"})


def _number(value: Any) -> float | None:
    if isinstance(value, Mapping):
        value = value.get("value")
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def _sha256(path: Path) -> str | None:
    if not path.is_file():
        return None
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _json_sha256(payload: Any) -> str:
    canonical = json.dumps(
        payload,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    return hashlib.sha256(canonical).hexdigest()


def _candidate_rows(match: Mapping[str, Any], limit: int) -> list[Mapping[str, Any]]:
    rows = [match, *(match.get("top_candidates") or [])]
    passport = match.get("evidence_governance") or {}
    rows.extend(passport.get("ambiguity_set") or [])
    output: list[Mapping[str, Any]] = []
    seen: set[str] = set()
    for candidate in rows:
        if not isinstance(candidate, Mapping):
            continue
        identity = compound_identity_key(candidate)
        if not identity or identity in seen:
            continue
        seen.add(identity)
        output.append(candidate)
        if len(output) >= limit:
            break
    return output


def _observed_peaks(row: Mapping[str, Any]) -> list[list[float]]:
    for key in ("raw_peaks", "query_peaks", "lcms_consensus_peaks"):
        peaks = row.get(key)
        if isinstance(peaks, list) and peaks:
            return peaks
    return []


def _instrument_domain(row: Mapping[str, Any], match: Mapping[str, Any]) -> dict[str, Any]:
    query = str(
        row.get("instrument_type")
        or match.get("query_instrument_type")
        or row.get("instrument")
        or match.get("query_instrument")
        or ""
    ).strip()
    normalized = query.lower()
    matched = any(token in normalized for token in ("hcd", "orbitrap", "ftms"))
    return {
        "query_instrument": query or "unknown",
        "model_instrument": "HCD",
        "domain_match": matched,
        "projection_applied": not matched,
        "reason": (
            "query_instrument_compatible_with_model_HCD_domain"
            if matched
            else "FIORA_open_model_accepts_HCD_only_query_projected_for_shadow"
        ),
    }


@dataclass(frozen=True)
class FioraRuntime:
    executable: Path
    model_path: Path
    timeout_seconds: int = 1800
    device: str = "cpu"
    min_probability: float = 0.001

    def availability(self) -> dict[str, Any]:
        executable_ok = self.executable.is_file() and os.access(
            self.executable, os.X_OK
        )
        model_ok = self.model_path.is_file()
        return {
            "available": bool(executable_ok and model_ok),
            "executable": str(self.executable),
            "executable_available": bool(executable_ok),
            "model_path": str(self.model_path),
            "model_available": bool(model_ok),
            "model_sha256": _sha256(self.model_path),
            "device": self.device,
        }


@dataclass(frozen=True)
class FiddleRuntime:
    python_executable: Path
    worker_path: Path
    model_path: Path
    rescore_model_path: Path
    config_path: Path
    timeout_seconds: int = 3600
    device: str = "cpu"
    batch_size: int = 32

    def availability(self) -> dict[str, Any]:
        python_ok = self.python_executable.is_file() and os.access(
            self.python_executable, os.X_OK
        )
        worker_ok = self.worker_path.is_file()
        model_ok = self.model_path.is_file()
        rescore_ok = self.rescore_model_path.is_file()
        config_ok = self.config_path.is_file()
        return {
            "available": bool(
                python_ok and worker_ok and model_ok and rescore_ok and config_ok
            ),
            "python_executable": str(self.python_executable),
            "python_available": bool(python_ok),
            "worker_path": str(self.worker_path),
            "worker_available": bool(worker_ok),
            "worker_sha256": _sha256(self.worker_path),
            "model_path": str(self.model_path),
            "model_available": bool(model_ok),
            "model_sha256": _sha256(self.model_path),
            "rescore_model_path": str(self.rescore_model_path),
            "rescore_model_available": bool(rescore_ok),
            "rescore_model_sha256": _sha256(self.rescore_model_path),
            "config_path": str(self.config_path),
            "config_available": bool(config_ok),
            "config_sha256": _sha256(self.config_path),
            "device": self.device,
            "batch_size": max(1, int(self.batch_size)),
        }


def _base_artifact(
    task_id: str,
    results: Sequence[Mapping[str, Any]],
    runtime: FioraRuntime,
) -> dict[str, Any]:
    return {
        "schema_version": SCHEMA_VERSION,
        "task_id": task_id,
        "module_id": FIORA_MODULE_ID,
        "status": "unavailable",
        "deployment_state": "shadow",
        "runtime_influence": False,
        "ranking_changed": False,
        "confidence_changed": False,
        "claim_scope_changed": False,
        "candidate_conditioned": True,
        "independent_measurement": False,
        "source_results_sha256": _json_sha256(results),
        "runtime": runtime.availability(),
        "model": {
            "id": "fiora-os-v1.0.0-msnlib-v7",
            "role": "candidate_structure_to_predicted_msms",
            "training_domain": "MSnLib_v7_HCD",
            "supported_adducts": sorted(FIORA_ALLOWED_ADDUCTS),
            "score_semantics": "predicted_vs_measured_spectral_similarity_not_probability",
        },
        "spectra": [],
        "summary": {},
        "interpretation": (
            "Forward-model evidence is a candidate-conditioned shadow audit. "
            "It cannot promote High or replace an experimental reference spectrum."
        ),
    }


def _prepare_requests(
    results: Sequence[Mapping[str, Any]],
    *,
    max_candidates_per_spectrum: int,
    max_spectra: int | None,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    requests: list[dict[str, Any]] = []
    skipped: list[dict[str, Any]] = []
    evaluated_spectra = 0
    for local_row_index, row in enumerate(results):
        try:
            row_index = int(
                row.get("_evidiscrim_source_result_index", local_row_index)
            ) if isinstance(row, Mapping) else local_row_index
        except (TypeError, ValueError):
            row_index = local_row_index
        match = row.get("match") if isinstance(row, Mapping) else None
        if not isinstance(match, Mapping):
            continue
        if max_spectra is not None and evaluated_spectra >= max_spectra:
            break
        observed = _observed_peaks(row)
        precursor = _number(row.get("precursor_mz"))
        collision_energy = _number(row.get("collision_energy"))
        if not observed or precursor is None or collision_energy is None:
            skipped.append({
                "row_index": row_index,
                "spectrum_id": row.get("spectrum_id"),
                "reason": "missing_observed_peaks_precursor_or_collision_energy",
            })
            continue
        domain = _instrument_domain(row, match)
        row_requests = 0
        for candidate_index, candidate in enumerate(
            _candidate_rows(match, max_candidates_per_spectrum)
        ):
            smiles = str(candidate.get("smiles") or "").strip()
            adduct = str(candidate.get("adduct") or "").strip()
            if not smiles:
                continue
            if adduct not in FIORA_ALLOWED_ADDUCTS:
                continue
            request_id = f"r{row_index}_c{candidate_index}"
            requests.append({
                "request_id": request_id,
                "row_index": row_index,
                "spectrum_id": str(row.get("spectrum_id") or row_index),
                "observed_peaks": observed,
                "query_precursor_mz": precursor,
                "collision_energy": collision_energy,
                "instrument_domain": domain,
                "candidate_index": candidate_index,
                "identity_key": compound_identity_key(candidate),
                "compound_name": str(candidate.get("compound_name") or ""),
                "formula": str(candidate.get("formula") or ""),
                "smiles": smiles,
                "inchikey": str(candidate.get("inchikey") or ""),
                "adduct": adduct,
                "collision_energy_unit": str(
                    row.get("collision_energy_unit")
                    or "instrument_native_unspecified"
                ),
                "production_top1": candidate_index == 0,
                "production_rank_score": _number(candidate.get("rank_score")),
            })
            row_requests += 1
        if row_requests:
            evaluated_spectra += 1
        else:
            skipped.append({
                "row_index": row_index,
                "spectrum_id": row.get("spectrum_id"),
                "reason": "no_candidate_with_structure_and_supported_explicit_adduct",
            })
    return requests, skipped


def _write_fiora_csv(requests: Iterable[Mapping[str, Any]], path: Path) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=("Name", "SMILES", "Precursor_type", "CE", "Instrument_type"),
        )
        writer.writeheader()
        for request in requests:
            writer.writerow({
                "Name": request["request_id"],
                "SMILES": request["smiles"],
                "Precursor_type": request["adduct"],
                "CE": request["collision_energy"],
                "Instrument_type": "HCD",
            })


def _prediction_map(path: Path) -> dict[str, dict[str, Any]]:
    predictions: dict[str, dict[str, Any]] = {}
    for spectrum in parse_mgf_file(path):
        title = str(spectrum.get("title") or "").strip()
        if title:
            predictions[title] = spectrum
    return predictions


def _score_prediction(
    request: Mapping[str, Any], prediction: Mapping[str, Any]
) -> dict[str, Any]:
    query_precursor = request["query_precursor_mz"]
    predicted_precursor = _number(prediction.get("precursor_mz")) or query_precursor
    config = resolve_peak_processing_config("sensitive")
    query, query_audit = process_peaks(
        request["observed_peaks"], query_precursor, config=config
    )
    predicted, predicted_audit = process_peaks(
        prediction.get("peaks") or [], predicted_precursor, config=config
    )
    cosine = cosine_score(
        query,
        predicted,
        query_pmz=query_precursor,
        lib_pmz=predicted_precursor,
    )
    entropy = entropy_score(query, predicted)
    reverse = reverse_match_ratio(
        query,
        predicted,
        query_pmz=query_precursor,
        lib_pmz=predicted_precursor,
    )
    descriptive = 0.5 * cosine + 0.5 * entropy
    return {
        "status": "available",
        "identity_key": request["identity_key"],
        "compound_name": request["compound_name"],
        "formula": request["formula"],
        "adduct": request["adduct"],
        "production_top1": request["production_top1"],
        "production_rank_score": request["production_rank_score"],
        "forward_similarity": round(descriptive, 6),
        "modified_cosine": round(cosine, 6),
        "spectral_entropy": round(entropy, 6),
        "predicted_peak_coverage": round(reverse, 6),
        "query_peak_count": len(query),
        "predicted_peak_count": len(predicted),
        "query_retained_pct": query_audit.get("retained_pct"),
        "predicted_retained_pct": predicted_audit.get("retained_pct"),
        "runtime_influence": False,
        "score_contribution": 0.0,
    }


def build_fiora_forward_shadow(
    results: Sequence[Mapping[str, Any]],
    *,
    task_id: str,
    runtime: FioraRuntime,
    max_candidates_per_spectrum: int = 5,
    max_spectra: int | None = None,
) -> dict[str, Any]:
    """Run FIORA once for a task and return a separate shadow artifact."""
    artifact = _base_artifact(task_id, results, runtime)
    if not artifact["runtime"]["available"]:
        artifact["summary"] = {"reason": "executable_or_model_unavailable"}
        return artifact
    requests, skipped = _prepare_requests(
        results,
        max_candidates_per_spectrum=max(1, int(max_candidates_per_spectrum)),
        max_spectra=max_spectra,
    )
    artifact["summary"]["skipped_spectra"] = len(skipped)
    artifact["skipped"] = skipped
    if not requests:
        artifact["status"] = "not_applicable"
        artifact["summary"]["reason"] = "no_applicable_candidate_requests"
        return artifact

    with tempfile.TemporaryDirectory(prefix="mssystem-fiora-") as temp_dir:
        temp = Path(temp_dir)
        input_csv = temp / "requests.csv"
        output_mgf = temp / "predictions.mgf"
        _write_fiora_csv(requests, input_csv)
        command = [
            str(runtime.executable),
            "-i", str(input_csv),
            "-o", str(output_mgf),
            "--model", str(runtime.model_path),
            "--dev", runtime.device,
            "--min_prob", str(float(runtime.min_probability)),
            "--no-rt",
            "--no-ccs",
            "--no-annotation",
        ]
        try:
            completed = subprocess.run(
                command,
                check=False,
                capture_output=True,
                text=True,
                timeout=max(1, int(runtime.timeout_seconds)),
            )
        except (OSError, subprocess.SubprocessError) as exc:
            artifact["status"] = "failed"
            artifact["summary"].update({
                "reason": "fiora_process_failed",
                "error": f"{type(exc).__name__}: {exc}"[:500],
                "candidate_requests": len(requests),
            })
            return artifact
        artifact["runtime"]["command"] = command
        artifact["runtime"]["returncode"] = completed.returncode
        artifact["runtime"]["stdout_tail"] = completed.stdout[-1000:]
        artifact["runtime"]["stderr_tail"] = completed.stderr[-1000:]
        if completed.returncode != 0 or not output_mgf.is_file():
            artifact["status"] = "failed"
            artifact["summary"].update({
                "reason": "fiora_nonzero_or_missing_output",
                "candidate_requests": len(requests),
            })
            return artifact
        predictions = _prediction_map(output_mgf)

    grouped: dict[int, list[dict[str, Any]]] = {}
    request_by_row: dict[int, dict[str, Any]] = {}
    prediction_failures = 0
    for request in requests:
        row_index = int(request["row_index"])
        request_by_row.setdefault(row_index, request)
        prediction = predictions.get(request["request_id"])
        if prediction is None:
            prediction_failures += 1
            continue
        grouped.setdefault(row_index, []).append(
            _score_prediction(request, prediction)
        )

    alternative_winners = 0
    stable_top1 = 0
    domain_shifted = 0
    for row_index in sorted(grouped):
        candidate_scores = sorted(
            grouped[row_index],
            key=lambda item: (
                -float(item["forward_similarity"]),
                str(item["identity_key"]),
            ),
        )
        for rank, candidate in enumerate(candidate_scores, start=1):
            candidate["forward_rank"] = rank
        production = next(
            (item for item in candidate_scores if item["production_top1"]),
            None,
        )
        winner = candidate_scores[0]
        top1_preserved = bool(
            production and winner["identity_key"] == production["identity_key"]
        )
        stable_top1 += int(top1_preserved)
        alternative_winners += int(not top1_preserved)
        request = request_by_row[row_index]
        domain = request["instrument_domain"]
        domain_shifted += int(not domain["domain_match"])
        artifact["spectra"].append({
            "row_index": row_index,
            "spectrum_id": request["spectrum_id"],
            "query_precursor_mz": request["query_precursor_mz"],
            "collision_energy": request["collision_energy"],
            "instrument_domain": domain,
            "candidate_scores": candidate_scores,
            "forward_top1_identity_key": winner["identity_key"],
            "forward_top1_compound_name": winner["compound_name"],
            "production_top1_preserved": top1_preserved,
            "production_top1_forward_rank": (
                production["forward_rank"] if production else None
            ),
            "counterfactual_candidate_order_changed": not top1_preserved,
            "runtime_influence": False,
        })

    artifact["status"] = "completed_shadow"
    artifact["summary"].update({
        "candidate_requests": len(requests),
        "predictions_returned": len(predictions),
        "prediction_failures": prediction_failures,
        "spectra_evaluated": len(artifact["spectra"]),
        "production_top1_preserved": stable_top1,
        "counterfactual_alternative_winner": alternative_winners,
        "query_instrument_domain_shifted_spectra": domain_shifted,
        "production_ranking_changed": False,
        "high_promotions": 0,
    })
    return artifact


def _prepare_evidiscrim_requests(
    results: Sequence[Mapping[str, Any]],
    *,
    policy: Mapping[str, Any],
    collision_energy_grid: Sequence[float] | None,
    max_candidates_per_spectrum: int,
    max_spectra: int | None,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Expand one FIORA request per candidate and instrument-native CE value."""

    base_requests, skipped = _prepare_requests(
        results,
        max_candidates_per_spectrum=max_candidates_per_spectrum,
        max_spectra=max_spectra,
    )
    configured_grid = (
        collision_energy_grid
        if collision_energy_grid is not None
        else (policy.get("energy_selection") or {}).get(
            "default_collision_energy_grid", ()
        )
    )
    parsed_grid = sorted({
        float(value)
        for value in configured_grid or ()
        if _number(value) is not None and float(value) > 0
    })
    requests: list[dict[str, Any]] = []
    for base in base_requests:
        observed_energy = float(base["collision_energy"])
        unit = str(
            base.get("collision_energy_unit")
            or "instrument_native_unspecified"
        )
        # Without an explicit native unit, adding a numeric CE grid would imply
        # an unverified conversion.  The observed number may still be audited.
        energies = (
            sorted(set(parsed_grid) | {observed_energy})
            if not unit.endswith("unspecified")
            else [observed_energy]
        )
        for energy_index, energy in enumerate(energies):
            request = dict(base)
            request["request_id"] = (
                f"{base['request_id']}_e{energy_index}_{energy:g}"
            )
            request["collision_energy"] = energy
            request["observed_collision_energy"] = observed_energy
            request["energy_transfer_qualified"] = bool(
                (base.get("instrument_domain") or {}).get("domain_match")
                and not unit.endswith("unspecified")
            )
            requests.append(request)
    return requests, skipped


def _run_fiora_request_batch(
    requests: Sequence[Mapping[str, Any]],
    runtime: FioraRuntime,
) -> dict[str, Any]:
    """Execute one FIORA batch and return predictions plus bounded logs."""

    with tempfile.TemporaryDirectory(prefix="mssystem-evidiscrim-fiora-") as temp_dir:
        temp = Path(temp_dir)
        input_csv = temp / "requests.csv"
        output_mgf = temp / "predictions.mgf"
        _write_fiora_csv(requests, input_csv)
        command = [
            str(runtime.executable),
            "-i", str(input_csv),
            "-o", str(output_mgf),
            "--model", str(runtime.model_path),
            "--dev", runtime.device,
            "--min_prob", str(float(runtime.min_probability)),
            "--no-rt",
            "--no-ccs",
            "--no-annotation",
        ]
        try:
            completed = subprocess.run(
                command,
                check=False,
                capture_output=True,
                text=True,
                timeout=max(1, int(runtime.timeout_seconds)),
            )
        except (OSError, subprocess.SubprocessError) as exc:
            return {
                "status": "failed",
                "reason": "fiora_process_failed",
                "error": f"{type(exc).__name__}: {exc}"[:500],
                "command": command,
                "predictions": {},
            }
        if completed.returncode != 0 or not output_mgf.is_file():
            return {
                "status": "failed",
                "reason": "fiora_nonzero_or_missing_output",
                "returncode": completed.returncode,
                "stdout_tail": completed.stdout[-1000:],
                "stderr_tail": completed.stderr[-1000:],
                "command": command,
                "predictions": {},
            }
        return {
            "status": "completed",
            "returncode": completed.returncode,
            "stdout_tail": completed.stdout[-1000:],
            "stderr_tail": completed.stderr[-1000:],
            "command": command,
            "predictions": _prediction_map(output_mgf),
        }


def _evidiscrim_prediction_record(
    request: Mapping[str, Any], prediction: Mapping[str, Any]
) -> dict[str, Any]:
    query_precursor = float(request["query_precursor_mz"])
    predicted_precursor = _number(prediction.get("precursor_mz")) or query_precursor
    config = resolve_peak_processing_config("sensitive")
    predicted, audit = process_peaks(
        prediction.get("peaks") or [], predicted_precursor, config=config
    )
    return {
        "row_index": int(request["row_index"]),
        "spectrum_id": request["spectrum_id"],
        "identity_key": request["identity_key"],
        "compound_name": request["compound_name"],
        "formula": request["formula"],
        "smiles": request["smiles"],
        "inchikey": request.get("inchikey") or "",
        "adduct": request["adduct"],
        "production_top1": bool(request["production_top1"]),
        "collision_energy": float(request["collision_energy"]),
        "collision_energy_unit": str(request["collision_energy_unit"]),
        "observed_collision_energy": float(request["observed_collision_energy"]),
        "predicted_precursor_mz": predicted_precursor,
        "predicted_peaks": [
            [round(float(mz), 6), round(float(intensity), 8)]
            for mz, intensity in predicted
        ],
        "predicted_peak_count": len(predicted),
        "predicted_retained_pct": audit.get("retained_pct"),
        "provider_id": "fiora-os-v1.0.0-msnlib-v7",
        "model_id": "fiora-os-v1.0.0-msnlib-v7",
        "architecture_family": "fiora_local_neighborhood_fragmentation_gnn",
        "base_architecture_family": (
            "fiora_local_neighborhood_fragmentation_gnn"
        ),
        "checkpoint_id": "msnlib-v7-hcd",
        "instrument_domain": dict(request.get("instrument_domain") or {}),
        "energy_transfer_qualified": bool(
            request.get("energy_transfer_qualified", False)
        ),
        "candidate_conditioned": True,
        "independent_measurement": False,
    }


def build_fiora_evidiscrim_shadow(
    results: Sequence[Mapping[str, Any]],
    *,
    task_id: str,
    runtime: FioraRuntime,
    policy: Mapping[str, Any],
    collision_energy_grid: Sequence[float] | None = None,
    max_candidates_per_spectrum: int = 5,
    max_spectra: int | None = None,
    source_results_sha256_override: str | None = None,
    source_results_immutability_verified: bool = False,
) -> dict[str, Any]:
    """Run FIORA over candidate/CE hypotheses and build EviDiscrim-MS."""

    runtime_state = runtime.availability()
    model = {
        "id": "fiora-os-v1.0.0-msnlib-v7",
        "role": "candidate_structure_collision_energy_to_predicted_msms",
        "training_domain": "MSnLib_v7_HCD",
        "supported_adducts": sorted(FIORA_ALLOWED_ADDUCTS),
        "score_semantics": "candidate_conditioned_prediction_not_probability",
        "uncertainty_semantics": (
            "single_architecture_only_not_qualified_for_pair_resolution_until_"
            "combined_with_an_independent_forward_architecture"
        ),
    }
    if not runtime_state["available"]:
        artifact = build_evidiscrim_shadow(
            results,
            task_id=task_id,
            prediction_records=(),
            policy=policy,
            model=model,
            runtime=runtime_state,
            max_candidates_per_spectrum=max_candidates_per_spectrum,
            source_results_sha256_override=source_results_sha256_override,
            source_results_immutability_verified=(
                source_results_immutability_verified
            ),
        )
        artifact["status"] = "unavailable"
        artifact["summary"]["reason"] = "executable_or_model_unavailable"
        return EviDiscrimArtifact.model_validate(artifact).model_dump(mode="json")

    requests, skipped = _prepare_evidiscrim_requests(
        results,
        policy=policy,
        collision_energy_grid=collision_energy_grid,
        max_candidates_per_spectrum=max(1, int(max_candidates_per_spectrum)),
        max_spectra=max_spectra,
    )
    if not requests:
        artifact = build_evidiscrim_shadow(
            results,
            task_id=task_id,
            prediction_records=(),
            policy=policy,
            model=model,
            runtime=runtime_state,
            skipped=skipped,
            max_candidates_per_spectrum=max_candidates_per_spectrum,
            source_results_sha256_override=source_results_sha256_override,
            source_results_immutability_verified=(
                source_results_immutability_verified
            ),
        )
        artifact["summary"]["reason"] = "no_applicable_candidate_energy_requests"
        return EviDiscrimArtifact.model_validate(artifact).model_dump(mode="json")

    execution = _run_fiora_request_batch(requests, runtime)
    runtime_state.update({
        key: value for key, value in execution.items()
        if key != "predictions"
    })
    if execution["status"] != "completed":
        artifact = build_evidiscrim_shadow(
            results,
            task_id=task_id,
            prediction_records=(),
            policy=policy,
            model=model,
            runtime=runtime_state,
            skipped=skipped,
            max_candidates_per_spectrum=max_candidates_per_spectrum,
            source_results_sha256_override=source_results_sha256_override,
            source_results_immutability_verified=(
                source_results_immutability_verified
            ),
        )
        artifact["status"] = "failed"
        artifact["summary"].update({
            "reason": execution.get("reason"),
            "candidate_energy_requests": len(requests),
        })
        return EviDiscrimArtifact.model_validate(artifact).model_dump(mode="json")

    predictions = execution["predictions"]
    records: list[dict[str, Any]] = []
    prediction_failures = 0
    for request in requests:
        prediction = predictions.get(request["request_id"])
        if prediction is None:
            prediction_failures += 1
            continue
        record = _evidiscrim_prediction_record(request, prediction)
        if record["predicted_peaks"]:
            records.append(record)
        else:
            prediction_failures += 1
    artifact = build_evidiscrim_shadow(
        results,
        task_id=task_id,
        prediction_records=records,
        policy=policy,
        model=model,
        runtime=runtime_state,
        skipped=skipped,
        max_candidates_per_spectrum=max_candidates_per_spectrum,
        source_results_sha256_override=source_results_sha256_override,
        source_results_immutability_verified=(
            source_results_immutability_verified
        ),
    )
    artifact["summary"].update({
        "candidate_energy_requests": len(requests),
        "fiora_predictions_returned": len(predictions),
        "usable_prediction_records": len(records),
        "prediction_failures": prediction_failures,
        "query_instrument_domain_shifted_requests": sum(
            not bool((request.get("instrument_domain") or {}).get("domain_match"))
            for request in requests
        ),
    })
    return EviDiscrimArtifact.model_validate(artifact).model_dump(mode="json")


def _formula_base_artifact(
    task_id: str,
    results: Sequence[Mapping[str, Any]],
    runtime: FiddleRuntime,
) -> dict[str, Any]:
    return {
        "schema_version": SCHEMA_VERSION,
        "task_id": task_id,
        "module_id": FIDDLE_MODULE_ID,
        "status": "unavailable",
        "deployment_state": "shadow",
        "runtime_influence": False,
        "ranking_changed": False,
        "confidence_changed": False,
        "claim_scope_changed": False,
        "candidate_conditioned": True,
        "independent_measurement": False,
        "source_results_sha256": _json_sha256(results),
        "runtime": runtime.availability(),
        "model": {
            "id": "msfiddle-2.1.0-qtof",
            "role": "measured_msms_to_molecular_formula_candidates",
            "training_domain": "published_QTOF_formula_prediction_model",
            "supported_adducts": sorted(FIDDLE_SUPPORTED_ADDUCTS),
            "score_semantics": (
                "model_candidate_score_not_probability_not_project_qvalue"
            ),
            "output_scope": "molecular_formula_only_not_2d_or_stereochemical_identity",
        },
        "spectra": [],
        "summary": {},
        "interpretation": (
            "Formula predictions are an open-world shadow candidate source. "
            "Agreement with a selected library formula is descriptive and is not "
            "ground-truth accuracy. Predictions cannot promote High without "
            "identity-disjoint project validation and empirical calibration."
        ),
    }


def _formula_instrument_domain(row: Mapping[str, Any]) -> dict[str, Any]:
    query = str(row.get("instrument_type") or row.get("instrument") or "").strip()
    normalized = query.lower().replace("-", "")
    matched = "tof" in normalized
    return {
        "query_instrument": query or "unknown",
        "model_instrument": "QTOF",
        "domain_match": matched,
        "reason": (
            "query_instrument_matches_QTOF_domain"
            if matched
            else "query_instrument_not_confirmed_as_QTOF"
        ),
    }


def _formula_adduct_hypothesis(
    row: Mapping[str, Any], match: Mapping[str, Any] | None
) -> dict[str, Any] | None:
    lcms_adduct = str(row.get("lcms_adduct") or "").strip()
    if lcms_adduct in FIDDLE_SUPPORTED_ADDUCTS:
        return {
            "adduct": lcms_adduct,
            "source": "lcms_adduct",
            "hypothesis": False,
            "candidate_conditioned": False,
            "library_formula_recall_eligible": True,
        }
    match_adduct = str((match or {}).get("adduct") or "").strip()
    if match_adduct in FIDDLE_SUPPORTED_ADDUCTS:
        return {
            "adduct": match_adduct,
            "source": "library_match_adduct",
            "hypothesis": False,
            "candidate_conditioned": True,
            "library_formula_recall_eligible": True,
        }
    ion_mode = str(row.get("ion_mode") or row.get("polarity") or "").upper()
    if "NEG" in ion_mode or ion_mode in {"-", "NEGATIVE"}:
        inferred = "[M-H]-"
    elif "POS" in ion_mode or ion_mode in {"+", "POSITIVE"}:
        inferred = "[M+H]+"
    else:
        return None
    return {
        "adduct": inferred,
        "source": "ion_mode_hypothesis",
        "hypothesis": True,
        "candidate_conditioned": False,
        "library_formula_recall_eligible": False,
        "conflicting_explicit_adduct": lcms_adduct or match_adduct or None,
    }


def _filter_formula_peaks(
    peaks: Sequence[Any], precursor_mz: float
) -> tuple[list[float], list[float], dict[str, Any]]:
    upper_bound = min(1500.0, float(precursor_mz) + 1.5)
    filtered: list[tuple[float, float]] = []
    for peak in peaks:
        if not isinstance(peak, Sequence) or isinstance(peak, (str, bytes)):
            continue
        if len(peak) < 2:
            continue
        mz = _number(peak[0])
        intensity = _number(peak[1])
        if mz is None or intensity is None or intensity <= 0:
            continue
        if 50.0 <= mz <= upper_bound:
            filtered.append((mz, intensity))
    if len(filtered) > 500:
        filtered = sorted(filtered, key=lambda item: -item[1])[:500]
    filtered.sort(key=lambda item: item[0])
    mz_array = [item[0] for item in filtered]
    intensity_array = [item[1] for item in filtered]
    return mz_array, intensity_array, {
        "raw_peak_count": len(peaks),
        "retained_peak_count": len(filtered),
        "retained_pct": round(100.0 * len(filtered) / max(1, len(peaks)), 3),
        "minimum_mz": 50.0,
        "maximum_mz": round(upper_bound, 6),
        "maximum_peaks": 500,
        "precursor_region_policy": "retain_product_ions_up_to_precursor_plus_1.5_Da",
    }


def _prepare_formula_requests(
    results: Sequence[Mapping[str, Any]],
    *,
    max_spectra: int | None,
    precursor_correction: Mapping[str, Any] | None = None,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    requests: list[dict[str, Any]] = []
    skipped: list[dict[str, Any]] = []
    for row_index, row in enumerate(results):
        if max_spectra is not None and len(requests) >= max_spectra:
            break
        match = row.get("match") if isinstance(row, Mapping) else None
        match = match if isinstance(match, Mapping) else None
        observed_precursor_mz = _number(row.get("precursor_mz"))
        precursor_mz = observed_precursor_mz
        correction_applied = False
        if observed_precursor_mz is not None and precursor_correction is not None:
            from matching.lockmass_recalibration import corrected_mz

            corrected = corrected_mz(observed_precursor_mz, precursor_correction)
            if corrected is not None:
                precursor_mz = corrected
                correction_applied = True
        collision_energy = _number(row.get("collision_energy"))
        observed = _observed_peaks(row)
        if precursor_mz is None or collision_energy is None or not observed:
            skipped.append({
                "row_index": row_index,
                "spectrum_id": row.get("spectrum_id"),
                "reason": "missing_observed_peaks_precursor_or_collision_energy",
            })
            continue
        charge = _number(row.get("precursor_charge"))
        if charge is not None and abs(charge) != 1:
            skipped.append({
                "row_index": row_index,
                "spectrum_id": row.get("spectrum_id"),
                "reason": "non_singly_charged_precursor_outside_simple_adduct_adapter",
            })
            continue
        adduct_evidence = _formula_adduct_hypothesis(row, match)
        if adduct_evidence is None:
            skipped.append({
                "row_index": row_index,
                "spectrum_id": row.get("spectrum_id"),
                "reason": "missing_supported_adduct_or_polarity_hypothesis",
            })
            continue
        mz_array, intensity_array, peak_audit = _filter_formula_peaks(
            observed, precursor_mz
        )
        if len(mz_array) < 5:
            skipped.append({
                "row_index": row_index,
                "spectrum_id": row.get("spectrum_id"),
                "reason": "fewer_than_five_valid_product_ions",
                "peak_processing": peak_audit,
            })
            continue
        if max(intensity_array) == min(intensity_array):
            skipped.append({
                "row_index": row_index,
                "spectrum_id": row.get("spectrum_id"),
                "reason": "flat_intensity_profile_not_encodable",
                "peak_processing": peak_audit,
            })
            continue
        library_formula = normalize_formula_text((match or {}).get("formula"))
        request_id = f"f{row_index}"
        requests.append({
            "request_id": request_id,
            "row_index": row_index,
            "spectrum_id": str(row.get("spectrum_id") or row_index),
            "matched": match is not None,
            "existing_confidence": str((match or {}).get("confidence") or ""),
            "identity_key": compound_identity_key(match or {}),
            "library_formula": library_formula,
            "library_formula_recall_eligible": bool(
                library_formula and adduct_evidence["library_formula_recall_eligible"]
            ),
            "adduct_evidence": adduct_evidence,
            "instrument_domain": _formula_instrument_domain(row),
            "peak_processing": peak_audit,
            "precursor_mz_observed": observed_precursor_mz,
            "precursor_mz_model_input": precursor_mz,
            "precursor_correction_applied": correction_applied,
            "worker_spectrum": {
                "id": request_id,
                "mz_array": mz_array,
                "intensity_array": intensity_array,
                "precursor_mz": precursor_mz,
                "adduct": adduct_evidence["adduct"],
                "collision_energy": collision_energy,
            },
        })
    return requests, skipped


def _formula_candidates(prediction: Mapping[str, Any]) -> list[dict[str, Any]]:
    candidates: list[dict[str, Any]] = []
    for rank, raw in enumerate(prediction.get("candidates") or [], start=1):
        if not isinstance(raw, Mapping):
            continue
        formula = normalize_formula_text(raw.get("formula"))
        if not formula:
            continue
        candidates.append({
            "rank": rank,
            "formula": formula,
            "score": _number(raw.get("score")),
            "neutral_mass": _number(raw.get("mass")),
            "score_semantics": "model_candidate_score_not_probability",
            "runtime_influence": False,
        })
    return candidates


def _rate(numerator: int, denominator: int) -> float | None:
    if denominator <= 0:
        return None
    return round(numerator / denominator, 6)


def build_fiddle_formula_shadow(
    results: Sequence[Mapping[str, Any]],
    *,
    task_id: str,
    runtime: FiddleRuntime,
    top_k: int = 5,
    max_spectra: int | None = None,
    precursor_correction: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Run msFIDDLE in an isolated runtime and return formula-only evidence."""
    artifact = _formula_base_artifact(task_id, results, runtime)
    if not artifact["runtime"]["available"]:
        artifact["summary"] = {"reason": "python_worker_or_model_assets_unavailable"}
        return artifact
    requests, skipped = _prepare_formula_requests(
        results,
        max_spectra=max_spectra,
        precursor_correction=precursor_correction,
    )
    artifact["precursor_mass_correction"] = {
        "requested": precursor_correction is not None,
        "status": (
            str((precursor_correction or {}).get("status") or "not_requested")
        ),
        "module_id": (precursor_correction or {}).get("module_id"),
        "candidate_independent": bool(
            (precursor_correction or {}).get("candidate_independent")
        ),
        "chemical_identity_used_for_fit": bool(
            (precursor_correction or {}).get("chemical_identity_used_for_fit")
        ),
        "correction": (precursor_correction or {}).get("correction"),
    }
    artifact["skipped"] = skipped
    artifact["summary"]["skipped_spectra"] = len(skipped)
    if not requests:
        artifact["status"] = "not_applicable"
        artifact["summary"]["reason"] = "no_applicable_spectra"
        return artifact

    with tempfile.TemporaryDirectory(prefix="mssystem-fiddle-") as temp_dir:
        temp = Path(temp_dir)
        input_json = temp / "spectra.json"
        output_json = temp / "predictions.json"
        input_json.write_text(
            json.dumps(
                {"spectra": [request["worker_spectrum"] for request in requests]},
                ensure_ascii=False,
            ),
            encoding="utf-8",
        )
        command = [
            str(runtime.python_executable),
            str(runtime.worker_path),
            "--input", str(input_json),
            "--output", str(output_json),
            "--model", str(runtime.model_path),
            "--rescore-model", str(runtime.rescore_model_path),
            "--config", str(runtime.config_path),
            "--instrument-type", "qtof",
            "--device", runtime.device,
            "--batch-size", str(max(1, int(runtime.batch_size))),
            "--top-k", str(max(1, int(top_k))),
        ]
        try:
            completed = subprocess.run(
                command,
                check=False,
                capture_output=True,
                text=True,
                timeout=max(1, int(runtime.timeout_seconds)),
            )
        except (OSError, subprocess.SubprocessError) as exc:
            artifact["status"] = "failed"
            artifact["summary"].update({
                "reason": "fiddle_process_failed",
                "error": f"{type(exc).__name__}: {exc}"[:500],
                "spectrum_requests": len(requests),
            })
            return artifact
        artifact["runtime"]["command"] = command
        artifact["runtime"]["returncode"] = completed.returncode
        artifact["runtime"]["stdout_tail"] = completed.stdout[-1000:]
        artifact["runtime"]["stderr_tail"] = completed.stderr[-1000:]
        if completed.returncode != 0 or not output_json.is_file():
            artifact["status"] = "failed"
            artifact["summary"].update({
                "reason": "fiddle_nonzero_or_missing_output",
                "spectrum_requests": len(requests),
            })
            return artifact
        try:
            raw_output = json.loads(output_json.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            artifact["status"] = "failed"
            artifact["summary"].update({
                "reason": "fiddle_invalid_output",
                "error": f"{type(exc).__name__}: {exc}"[:500],
            })
            return artifact

    predictions = {
        str(item.get("id")): item
        for item in raw_output.get("predictions", [])
        if isinstance(item, Mapping) and item.get("id") is not None
    }
    evaluable = top1_hits = topk_hits = 0
    high_evaluable = high_top1_hits = high_topk_hits = 0
    lcms_evaluable = lcms_topk_hits = 0
    matched_requests = unmatched_requests = 0
    unmatched_with_candidates = 0
    hypothesis_requests = 0
    domain_shifted = 0
    prediction_failures = 0
    for request in requests:
        matched_requests += int(request["matched"])
        unmatched_requests += int(not request["matched"])
        hypothesis_requests += int(request["adduct_evidence"]["hypothesis"])
        domain_shifted += int(not request["instrument_domain"]["domain_match"])
        prediction = predictions.get(request["request_id"])
        if prediction is None:
            prediction_failures += 1
            continue
        candidates = _formula_candidates(prediction)
        predicted_formulas = [item["formula"] for item in candidates]
        library_formula = request["library_formula"]
        top1_agreement = bool(
            library_formula and predicted_formulas and predicted_formulas[0] == library_formula
        )
        topk_agreement = bool(library_formula and library_formula in predicted_formulas)
        recall_eligible = bool(request["library_formula_recall_eligible"])
        if recall_eligible:
            evaluable += 1
            top1_hits += int(top1_agreement)
            topk_hits += int(topk_agreement)
            if request["existing_confidence"] == "High":
                high_evaluable += 1
                high_top1_hits += int(top1_agreement)
                high_topk_hits += int(topk_agreement)
            if request["adduct_evidence"]["source"] == "lcms_adduct":
                lcms_evaluable += 1
                lcms_topk_hits += int(topk_agreement)
        if not request["matched"] and candidates:
            unmatched_with_candidates += 1
        artifact["spectra"].append({
            "row_index": request["row_index"],
            "spectrum_id": request["spectrum_id"],
            "matched": request["matched"],
            "existing_confidence": request["existing_confidence"],
            "identity_key": request["identity_key"],
            "library_formula": library_formula or None,
            "library_formula_recall_eligible": recall_eligible,
            "library_formula_top1_agreement": top1_agreement,
            "library_formula_topk_agreement": topk_agreement,
            "adduct_evidence": request["adduct_evidence"],
            "instrument_domain": request["instrument_domain"],
            "peak_processing": request["peak_processing"],
            "precursor_mz_observed": request["precursor_mz_observed"],
            "precursor_mz_model_input": request["precursor_mz_model_input"],
            "precursor_correction_applied": request["precursor_correction_applied"],
            "formula_candidates": candidates,
            "prediction_metadata": prediction.get("metadata") or {},
            "runtime_influence": False,
        })

    artifact["status"] = "completed_shadow"
    artifact["summary"].update({
        "spectrum_requests": len(requests),
        "predictions_returned": len(predictions),
        "prediction_failures": prediction_failures,
        "spectra_reported": len(artifact["spectra"]),
        "matched_spectrum_requests": matched_requests,
        "unmatched_spectrum_requests": unmatched_requests,
        "unmatched_spectra_with_formula_candidates": unmatched_with_candidates,
        "ion_mode_adduct_hypothesis_requests": hypothesis_requests,
        "query_instrument_domain_shifted_spectra": domain_shifted,
        "precursor_correction_applied_spectra": sum(
            request["precursor_correction_applied"] for request in requests
        ),
        "library_formula_evaluable_explicit_supported_adduct": evaluable,
        "library_formula_top1_recalled": top1_hits,
        "library_formula_top5_recalled": topk_hits,
        "library_formula_top1_recall_rate": _rate(top1_hits, evaluable),
        "library_formula_top5_recall_rate": _rate(topk_hits, evaluable),
        "existing_high_spectra_evaluable": high_evaluable,
        "existing_high_formula_top1_recalled": high_top1_hits,
        "existing_high_formula_top5_recalled": high_topk_hits,
        "existing_high_formula_top5_recall_rate": _rate(
            high_topk_hits, high_evaluable
        ),
        "lcms_adduct_evaluable_spectra": lcms_evaluable,
        "lcms_adduct_formula_top5_recalled": lcms_topk_hits,
        "lcms_adduct_formula_top5_recall_rate": _rate(
            lcms_topk_hits, lcms_evaluable
        ),
        "evaluation_semantics": (
            "selected_library_formula_consistency_not_ground_truth_accuracy"
        ),
        "production_ranking_changed": False,
        "high_promotions": 0,
    })
    return artifact


def write_shadow_artifact(payload: Mapping[str, Any], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp-{os.getpid()}")
    with temporary.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)


__all__ = [
    "FIDDLE_SUPPORTED_ADDUCTS",
    "FIORA_ALLOWED_ADDUCTS",
    "FiddleRuntime",
    "FioraRuntime",
    "build_fiddle_formula_shadow",
    "build_fiora_evidiscrim_shadow",
    "build_fiora_forward_shadow",
    "write_shadow_artifact",
]
