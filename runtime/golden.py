"""Golden scientific-result projections for zero-impact refactoring.

Large task artifacts contain raw peaks and presentation metadata, so bytewise
hashes are too sensitive to harmless serialization changes.  This module hashes
an explicit scientific projection for every spectrum and separately freezes the
task-level units and evidence-governance summaries.
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Iterator, Mapping, Sequence


GOLDEN_SCHEMA_VERSION = "1.0.0"
MISSING = "__MSSYSTEM_MISSING__"
SUMMARY_PATHS: tuple[str, ...] = (
    "filename",
    "status",
    "total_spectra",
    "matched_spectra",
    "unique_compounds",
    "match_rate",
    "ms1_scans",
    "ion_mode",
    "ion_mode_counts",
    "compound_confidence_summary.unique_library_candidates",
    "compound_confidence_summary.high_supported_unique_compounds",
    "compound_confidence_summary.representative_high_unique_compounds",
    "compound_confidence_summary.ranking_changed",
    "identification_funnel.ms2_scans",
    "identification_funnel.ms2_precursor_events_approximate",
    "identification_funnel.matched_ms2_scans",
    "identification_funnel.matched_ms2_precursor_events",
    "identification_funnel.unique_library_candidates",
    "evidence_governance_summary.schema_version",
    "evidence_governance_summary.policy_id",
    "evidence_governance_summary.matched_spectra_audited",
    "evidence_governance_summary.unique_identities_audited",
    "evidence_governance_summary.claim_status_unique_identities",
    "evidence_governance_summary.ranking_changed",
    "blank_summary.requested",
    "blank_summary.feature_matches",
    "blank_summary.excluded_spectra",
    "blank_summary.assigned_blank_count",
    "blank_summary.skipped_blank_count",
    "expected_compounds_summary.status",
    "expected_compounds_summary.detected_feature_count",
    "expected_compounds_summary.candidate_target_rows",
    "expected_compounds_summary.counted_in_standard_spectral_library_matches",
    "params.db_mode",
    "params.cosine_threshold",
    "params.final_score_threshold",
    "params.mz_tolerance_da",
    "params.mz_tolerance_ppm",
    "params.top_n_peaks",
    "params.strict_evidence_gate",
    "params.blank_filter_policy",
)


def _canonical_json(value: Any) -> str:
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"))


def _get_path(record: Mapping[str, Any], dotted_path: str) -> Any:
    value: Any = record
    for part in dotted_path.split("."):
        if not isinstance(value, Mapping) or part not in value:
            return MISSING
        value = value[part]
    return value


def task_summary_projection(task: Mapping[str, Any]) -> dict[str, Any]:
    return {path: _get_path(task, path) for path in SUMMARY_PATHS}


def _match_projection(match: Mapping[str, Any] | None) -> dict[str, Any] | None:
    if not isinstance(match, Mapping):
        return None
    isopair = match.get("isopair_evidence") if isinstance(match.get("isopair_evidence"), Mapping) else {}
    governance = (
        match.get("evidence_governance")
        if isinstance(match.get("evidence_governance"), Mapping)
        else {}
    )
    adduct_context = (
        isopair.get("adduct_context")
        if isinstance(isopair.get("adduct_context"), Mapping)
        else {}
    )
    return {
        "compound_name": match.get("compound_name"),
        "inchikey": match.get("inchikey"),
        "reported_inchikey": match.get("reported_inchikey"),
        "formula": match.get("formula"),
        "adduct": match.get("adduct"),
        "lib_spectrum_id": match.get("lib_spectrum_id"),
        "cosine": match.get("cosine"),
        "entropy": match.get("entropy"),
        "final_score": match.get("final_score"),
        "rank_score": match.get("rank_score"),
        "raw_rank_score": match.get("raw_rank_score"),
        "confidence": match.get("confidence"),
        "confidence_class": match.get("confidence_class"),
        "identification_scope": match.get("identification_scope"),
        "annotation_level": match.get("annotation_level"),
        "isomer_ambiguous": match.get("isomer_ambiguous"),
        "isomer_resolution_status": match.get("isomer_resolution_status"),
        "isopair": {
            "status": isopair.get("status"),
            "reason": isopair.get("reason"),
            "evaluated": isopair.get("evaluated"),
            "runtime_influence": isopair.get("runtime_influence"),
            "ranking_changed": isopair.get("ranking_changed"),
            "counterfactual_ranking_changed": isopair.get("counterfactual_ranking_changed"),
            "production_abstention_reason": isopair.get("production_abstention_reason"),
            "adduct_status": adduct_context.get("status"),
            "adduct_source": adduct_context.get("evidence_source"),
        },
        "claim": {
            "status": governance.get("claim_status"),
            "scope": governance.get("allowed_claim_scope"),
            "release_eligible": governance.get("release_eligible_before_replication"),
            "ranking_changed": governance.get("ranking_changed"),
        },
    }


def spectrum_projection(record: Mapping[str, Any]) -> dict[str, Any]:
    return {
        "spectrum_id": record.get("spectrum_id"),
        "precursor_mz": record.get("precursor_mz"),
        "retention_time": record.get("retention_time"),
        "ion_mode": record.get("ion_mode"),
        "precursor_feature_id": record.get("precursor_feature_id"),
        "lcms_feature_id": record.get("lcms_feature_id"),
        "blank_status": record.get("blank_status"),
        "blank_flagged": record.get("blank_flagged"),
        "match": _match_projection(record.get("match")),
    }


def _iter_result_records(path: Path) -> Iterator[Mapping[str, Any]]:
    try:
        import ijson

        with path.open("rb") as handle:
            for record in ijson.items(handle, "item", use_float=True):
                if isinstance(record, Mapping):
                    yield record
        return
    except ImportError:
        pass
    with path.open("r", encoding="utf-8") as handle:
        payload = json.load(handle)
    if not isinstance(payload, list):
        raise ValueError(f"Expected a result array in {path}")
    for record in payload:
        if isinstance(record, Mapping):
            yield record


def semantic_result_digest(path: Path) -> dict[str, Any]:
    digest = hashlib.sha256()
    counts = {
        "spectra": 0,
        "matched": 0,
        "isopair_evaluated": 0,
        "isopair_runtime_ranking_changed": 0,
        "isopair_counterfactual_ranking_changed": 0,
    }
    for record in _iter_result_records(path):
        projection = spectrum_projection(record)
        digest.update(_canonical_json(projection).encode("utf-8"))
        digest.update(b"\n")
        counts["spectra"] += 1
        match = record.get("match")
        if not isinstance(match, Mapping):
            continue
        counts["matched"] += 1
        isopair = match.get("isopair_evidence")
        if isinstance(isopair, Mapping):
            counts["isopair_evaluated"] += int(bool(isopair.get("evaluated")))
            counts["isopair_runtime_ranking_changed"] += int(
                bool(isopair.get("ranking_changed"))
            )
            counts["isopair_counterfactual_ranking_changed"] += int(
                bool(isopair.get("counterfactual_ranking_changed"))
            )
    return {"semantic_sha256": digest.hexdigest(), **counts}


def _load_json(path: Path) -> Any:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def _task_ids(task_sets: Mapping[str, Any]) -> list[str]:
    ids: list[str] = []
    seen: set[str] = set()
    for cohort in task_sets.get("cohorts", []):
        for task_id in cohort.get("task_ids", []):
            if task_id not in seen:
                ids.append(str(task_id))
                seen.add(str(task_id))
    return ids


def build_golden_baseline(
    tasks_index_path: Path,
    task_sets_path: Path,
    results_root: Path,
    *,
    include_semantic_digests: bool = True,
) -> dict[str, Any]:
    tasks_index = _load_json(tasks_index_path)
    task_sets = _load_json(task_sets_path)
    tasks: dict[str, Any] = {}
    for task_id in _task_ids(task_sets):
        if task_id not in tasks_index:
            raise KeyError(f"Golden task {task_id} is missing from {tasks_index_path}")
        task_record: dict[str, Any] = {
            "summary": task_summary_projection(tasks_index[task_id]),
        }
        result_path = results_root / f"{task_id}_results.json"
        if include_semantic_digests:
            if not result_path.is_file():
                raise FileNotFoundError(result_path)
            task_record["result_semantics"] = semantic_result_digest(result_path)
        tasks[task_id] = task_record
    return {
        "schema_version": GOLDEN_SCHEMA_VERSION,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "task_sets": task_sets,
        "summary_paths": list(SUMMARY_PATHS),
        "tasks": tasks,
        "semantics_note": (
            "semantic_sha256 covers identity, scientific scores, confidence, blank state, "
            "IsoPair influence and EviGate claim fields; raw peak arrays and presentation-only "
            "metadata are intentionally excluded."
        ),
    }


def _diff(expected: Any, actual: Any, prefix: str = "") -> list[dict[str, Any]]:
    if isinstance(expected, Mapping) and isinstance(actual, Mapping):
        differences: list[dict[str, Any]] = []
        for key in sorted(set(expected) | set(actual)):
            path = f"{prefix}.{key}" if prefix else str(key)
            if key not in expected:
                differences.append({"path": path, "expected": MISSING, "actual": actual[key]})
            elif key not in actual:
                differences.append({"path": path, "expected": expected[key], "actual": MISSING})
            else:
                differences.extend(_diff(expected[key], actual[key], path))
        return differences
    if expected != actual:
        return [{"path": prefix, "expected": expected, "actual": actual}]
    return []


def verify_golden_baseline(
    baseline_path: Path,
    tasks_index_path: Path,
    results_root: Path,
    *,
    verify_semantic_digests: bool = True,
) -> dict[str, Any]:
    baseline = _load_json(baseline_path)
    tasks_index = _load_json(tasks_index_path)
    differences: list[dict[str, Any]] = []
    checked_tasks = 0
    for task_id, expected in baseline.get("tasks", {}).items():
        checked_tasks += 1
        if task_id not in tasks_index:
            differences.append(
                {"path": f"tasks.{task_id}", "expected": "present", "actual": MISSING}
            )
            continue
        actual_summary = task_summary_projection(tasks_index[task_id])
        differences.extend(
            _diff(expected.get("summary", {}), actual_summary, f"tasks.{task_id}.summary")
        )
        if verify_semantic_digests and "result_semantics" in expected:
            path = results_root / f"{task_id}_results.json"
            if not path.is_file():
                differences.append(
                    {
                        "path": f"tasks.{task_id}.result_semantics",
                        "expected": expected["result_semantics"],
                        "actual": MISSING,
                    }
                )
            else:
                actual_semantics = semantic_result_digest(path)
                differences.extend(
                    _diff(
                        expected["result_semantics"],
                        actual_semantics,
                        f"tasks.{task_id}.result_semantics",
                    )
                )
    return {
        "ok": not differences,
        "checked_tasks": checked_tasks,
        "semantic_digests_checked": verify_semantic_digests,
        "difference_count": len(differences),
        "differences": differences,
    }
