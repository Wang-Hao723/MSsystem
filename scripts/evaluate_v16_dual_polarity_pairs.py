#!/usr/bin/env python3
"""Replay explicitly declared V16 POS/NEG shadow-evidence pairs.

This command reads only the supplied manifest.  It does not scan task names,
infer pairings, load candidate truth, or call a task-mutating endpoint.  The
V16 service remains responsible for validating sidecar semantic hashes.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import tempfile
import time
import urllib.error
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable, Mapping, Sequence


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MANIFEST = ROOT / "config" / "v16_dual_polarity_evaluation_pairs.json"
DEFAULT_OUTPUT = (
    ROOT
    / "output"
    / "validation_eval"
    / "v16_dual_polarity_evaluation_pairs.json"
)
PAIR_ENDPOINT = "/api/workflows/dual-polarity-shadow"
HEALTH_ENDPOINT = "/api/health"
MANIFEST_SCHEMA = "mssystem.v16.dual-polarity-evaluation-pairs/1.0.0"
AUDIT_SCHEMA = "mssystem.v16.dual-polarity-evaluation-audit/1.0.0"
CONFIRMATION_SOURCE = "explicit_study_design_manifest_not_filename_inference"
PAIR_ID_PATTERN = re.compile(r"^DPS-[0-9a-f]{16}$")
TASK_ID_PATTERN = re.compile(r"^[0-9a-f]{12}$")
SHA256_PATTERN = re.compile(r"^[0-9a-f]{64}$")

PRODUCTION_UNCHANGED_FIELDS = (
    "production_top1_changed",
    "production_top5_changed",
    "production_score_changed",
    "production_decision_changed",
    "production_confidence_changed",
    "production_compound_count_changed",
)
SHADOW_GOVERNANCE_VALUES = {
    "ground_truth_used_at_inference": False,
    "sample_filename_used_for_pairing": False,
    "filenames_used_for_pairing": False,
    "cross_polarity_ms2_similarity_is_auxiliary_only": True,
}
OPENMS_COUNT_FIELDS = (
    "input_feature_count",
    "input_group_count",
    "eligible_group_count",
    "eligible_feature_count",
    "rejected_group_count",
    "projected_result_row_count",
    "bound_result_row_count",
    "result_binding_mismatch_count",
)
POST_FIELDS = (
    "positive_task_id",
    "negative_task_id",
    "same_sample_confirmed",
    "same_lc_method_confirmed",
    "pair_label",
)

RequestJSON = Callable[..., Mapping[str, Any]]


class EvaluationError(RuntimeError):
    """Raised for a malformed manifest or an invalid service response."""


class HTTPRequestError(EvaluationError):
    """HTTP failure with the remote status retained for the audit."""

    def __init__(self, message: str, *, status_code: int | None = None):
        super().__init__(message)
        self.status_code = status_code


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _canonical_sha256(value: Any) -> str:
    encoded = json.dumps(
        value,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _atomic_write_json(path: Path, payload: Mapping[str, Any]) -> None:
    """Write an audit through a same-directory fsync + atomic replace."""

    path = path.expanduser().resolve()
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            json.dump(
                payload,
                handle,
                ensure_ascii=False,
                indent=2,
                sort_keys=True,
                allow_nan=False,
            )
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
        try:
            directory_descriptor = os.open(path.parent, os.O_RDONLY)
        except OSError:
            return
        try:
            os.fsync(directory_descriptor)
        finally:
            os.close(directory_descriptor)
    finally:
        if temporary.exists():
            temporary.unlink()


def _load_json(path: Path) -> Any:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise EvaluationError(f"Could not read JSON manifest {path}: {exc}") from exc


def validate_manifest(payload: Any) -> list[dict[str, Any]]:
    """Validate explicit task IDs and confirmations without inspecting names."""

    if not isinstance(payload, Mapping):
        raise EvaluationError("Pair manifest must be a JSON object")
    if payload.get("schema_version") != MANIFEST_SCHEMA:
        raise EvaluationError("Pair manifest schema_version is unsupported")
    if payload.get("selection_method") != CONFIRMATION_SOURCE:
        raise EvaluationError("Pair selection must be an explicit study-design manifest")
    if payload.get("task_status_required") != "completed":
        raise EvaluationError("Pair manifest must require completed tasks")
    if payload.get("positive_task_ion_mode_required") != "POS":
        raise EvaluationError("Pair manifest must require exact POS tasks")
    if payload.get("negative_task_ion_mode_required") != "NEG":
        raise EvaluationError("Pair manifest must require exact NEG tasks")
    if payload.get("candidate_ground_truth_included") is not False:
        raise EvaluationError("Candidate ground truth must not be included")

    raw_pairs = payload.get("pairs")
    if not isinstance(raw_pairs, list) or not raw_pairs:
        raise EvaluationError("Pair manifest must contain at least one explicit pair")
    expected_count = payload.get("expected_pair_count")
    if not isinstance(expected_count, int) or expected_count != len(raw_pairs):
        raise EvaluationError("expected_pair_count does not match the explicit pair list")

    labels: set[str] = set()
    positive_ids: set[str] = set()
    negative_ids: set[str] = set()
    pairs: list[dict[str, Any]] = []
    required = set(POST_FIELDS) | {"confirmation_source"}
    for index, raw in enumerate(raw_pairs):
        if not isinstance(raw, Mapping):
            raise EvaluationError(f"Pair {index + 1} must be an object")
        missing = sorted(required - set(raw))
        if missing:
            raise EvaluationError(
                f"Pair {index + 1} is missing required fields: {', '.join(missing)}"
            )
        pair = {key: raw[key] for key in required}
        label = str(pair["pair_label"] or "").strip()
        positive = str(pair["positive_task_id"] or "").strip()
        negative = str(pair["negative_task_id"] or "").strip()
        if not label or label in labels:
            raise EvaluationError(f"Pair {index + 1} has an empty or duplicate label")
        if not TASK_ID_PATTERN.fullmatch(positive):
            raise EvaluationError(f"Pair {label} has an invalid positive task ID")
        if not TASK_ID_PATTERN.fullmatch(negative):
            raise EvaluationError(f"Pair {label} has an invalid negative task ID")
        if positive == negative:
            raise EvaluationError(f"Pair {label} must contain two distinct task IDs")
        if positive in positive_ids or negative in negative_ids:
            raise EvaluationError(f"Pair {label} reuses a task within one polarity")
        if pair["same_sample_confirmed"] is not True:
            raise EvaluationError(f"Pair {label} lacks same-sample confirmation")
        if pair["same_lc_method_confirmed"] is not True:
            raise EvaluationError(f"Pair {label} lacks same-LC-method confirmation")
        if pair["confirmation_source"] != CONFIRMATION_SOURCE:
            raise EvaluationError(f"Pair {label} has an unsupported confirmation source")
        labels.add(label)
        positive_ids.add(positive)
        negative_ids.add(negative)
        pair["pair_label"] = label
        pair["positive_task_id"] = positive
        pair["negative_task_id"] = negative
        pairs.append(pair)

    if positive_ids & negative_ids:
        raise EvaluationError("A task ID cannot be used as both POS and NEG")
    return pairs


def _normalise_base_url(value: str) -> str:
    url = str(value or "").strip().rstrip("/")
    parsed = urllib.parse.urlsplit(url)
    if parsed.scheme not in {"http", "https"} or not parsed.netloc:
        raise EvaluationError("base URL must be an absolute HTTP(S) URL")
    if parsed.query or parsed.fragment:
        raise EvaluationError("base URL must not contain a query or fragment")
    return url


def _request_json(
    *,
    method: str,
    url: str,
    payload: Mapping[str, Any] | None,
    timeout: float,
) -> Mapping[str, Any]:
    body = None
    headers = {"Accept": "application/json"}
    if payload is not None:
        body = json.dumps(payload, ensure_ascii=False, allow_nan=False).encode("utf-8")
        headers["Content-Type"] = "application/json"
    request = urllib.request.Request(
        url,
        data=body,
        method=method,
        headers=headers,
    )
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:
            status = int(getattr(response, "status", 200))
            raw = response.read()
    except urllib.error.HTTPError as exc:
        try:
            detail = exc.read().decode("utf-8", errors="replace")
        except Exception:
            detail = ""
        raise HTTPRequestError(
            f"{method} {url} returned HTTP {exc.code}: {detail[:1000]}",
            status_code=int(exc.code),
        ) from exc
    except (urllib.error.URLError, TimeoutError, OSError) as exc:
        raise HTTPRequestError(f"{method} {url} failed: {exc}") from exc
    if status < 200 or status >= 300:
        raise HTTPRequestError(
            f"{method} {url} returned HTTP {status}", status_code=status
        )
    try:
        decoded = json.loads(raw.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise HTTPRequestError(f"{method} {url} returned invalid JSON") from exc
    if not isinstance(decoded, Mapping):
        raise HTTPRequestError(f"{method} {url} did not return a JSON object")
    return decoded


def _service_preflight(
    base_url: str,
    *,
    timeout: float,
    request_json: RequestJSON,
) -> dict[str, Any]:
    health = request_json(
        method="GET",
        url=f"{base_url}{HEALTH_ENDPOINT}",
        payload=None,
        timeout=timeout,
    )
    generation = str(health.get("service_generation") or "").upper()
    capabilities = health.get("capabilities") or {}
    if generation != "V16":
        raise EvaluationError(f"Refusing non-V16 service generation: {generation or 'missing'}")
    runtime = health.get("runtime_manifest") or {}
    if health.get("status") != "ok":
        raise EvaluationError("Refusing a degraded V16 service")
    if not isinstance(runtime, Mapping):
        raise EvaluationError("V16 health lacks runtime-manifest state")
    if (
        runtime.get("status") != "fresh"
        or runtime.get("fresh") is not True
        or runtime.get("release_ready") is not True
        or runtime.get("restart_required") is not False
    ):
        raise EvaluationError(
            "Refusing a stale, restart-required, or unreleased V16 runtime"
        )
    if not isinstance(capabilities, Mapping) or capabilities.get(
        "v16_dual_polarity_shadow_evidence_graph"
    ) is not True:
        raise EvaluationError("V16 dual-polarity shadow capability is unavailable")
    if capabilities.get("dual_polarity_shadow_endpoint") != PAIR_ENDPOINT:
        raise EvaluationError("V16 service advertises an unexpected shadow endpoint")
    storage = health.get("storage") or {}
    if not isinstance(storage, Mapping) or storage.get("ok") is not True:
        raise EvaluationError("V16 runtime storage is not healthy")
    if storage.get("library_write_enabled") is not False:
        raise EvaluationError(
            "Evaluation requires a read-only production library"
        )
    return {
        "status": health.get("status"),
        "service_generation": generation,
        "service_version": health.get("service_version"),
        "runtime_snapshot_id": health.get("runtime_snapshot_id"),
        "runtime_manifest_status": runtime.get("status"),
        "runtime_manifest_fresh": runtime.get("fresh"),
        "release_ready": runtime.get("release_ready"),
        "restart_required": runtime.get("restart_required"),
        "dual_polarity_shadow_endpoint": PAIR_ENDPOINT,
        "library_write_enabled": storage.get("library_write_enabled"),
        "storage_ok": storage.get("ok"),
    }


def _nonnegative_int(value: Any, *, field: str) -> int:
    if value is None or value == "":
        return 0
    if isinstance(value, bool):
        raise EvaluationError(f"{field} must be a non-negative integer")
    try:
        converted = int(value)
    except (TypeError, ValueError) as exc:
        raise EvaluationError(f"{field} must be a non-negative integer") from exc
    if converted < 0:
        raise EvaluationError(f"{field} must be a non-negative integer")
    return converted


def _summary_count(summary: Mapping[str, Any], field: str) -> int:
    return _nonnegative_int(summary.get(field), field=f"summary.{field}")


def _openms_projection(
    metadata: Mapping[str, Any], pair: Mapping[str, Any]
) -> dict[str, Any]:
    audits = metadata.get("openms_adduct_group_validation")
    if not isinstance(audits, Mapping):
        raise EvaluationError("Sidecar lacks OpenMS validation audits")
    projected: dict[str, Any] = {}
    for polarity, task_field in (
        ("POS", "positive_task_id"),
        ("NEG", "negative_task_id"),
    ):
        task_id = str(pair[task_field])
        raw = audits.get(task_id)
        if not isinstance(raw, Mapping):
            raise EvaluationError(f"Sidecar lacks {polarity} OpenMS validation audit")
        values = {
            field: _nonnegative_int(
                raw.get(field), field=f"openms.{polarity}.{field}"
            )
            for field in OPENMS_COUNT_FIELDS
        }
        if raw.get("saved_result_rows_modified") is not False:
            raise EvaluationError(
                f"{polarity} OpenMS validation did not preserve saved result rows"
            )
        if raw.get("production_runtime_influence") is not False:
            raise EvaluationError(
                f"{polarity} OpenMS validation changed production runtime"
            )
        projected[polarity] = {
            "task_id": task_id,
            "status": raw.get("status"),
            **values,
            "manifest_rejection_reasons": list(
                raw.get("manifest_rejection_reasons") or []
            ),
            "saved_result_rows_modified": raw.get("saved_result_rows_modified"),
            "production_runtime_influence": raw.get(
                "production_runtime_influence"
            ),
        }
    projected["total"] = {
        field: sum(projected[polarity][field] for polarity in ("POS", "NEG"))
        for field in OPENMS_COUNT_FIELDS
    }
    return projected


def _metrics_from_sidecar(
    sidecar: Mapping[str, Any], pair: Mapping[str, Any]
) -> dict[str, Any]:
    summary = sidecar.get("summary")
    metadata = sidecar.get("metadata")
    invariants = sidecar.get("invariants")
    if not isinstance(summary, Mapping):
        raise EvaluationError("Sidecar lacks a summary object")
    if not isinstance(metadata, Mapping):
        raise EvaluationError("Sidecar lacks a metadata object")
    if not isinstance(invariants, Mapping):
        raise EvaluationError("Sidecar lacks an invariants object")

    production = {field: invariants.get(field) for field in PRODUCTION_UNCHANGED_FIELDS}
    invalid = [field for field, value in production.items() if value is not False]
    if invalid:
        raise EvaluationError(
            "V15 production invariant is not false: " + ", ".join(invalid)
        )
    invalid_governance = [
        field
        for field, expected in SHADOW_GOVERNANCE_VALUES.items()
        if invariants.get(field) is not expected
    ]
    if invalid_governance:
        raise EvaluationError(
            "Shadow governance invariant has an unexpected value: "
            + ", ".join(invalid_governance)
        )
    if metadata.get("source_file_sha256_unchanged_after_build") is not True:
        raise EvaluationError("Source V15 evidence was not verified unchanged")
    if metadata.get("task_compound_counts_unchanged") is not True:
        raise EvaluationError("Task compound counts were not verified unchanged")
    compound_increment = _summary_count(summary, "compound_count_increment")
    if compound_increment != 0:
        raise EvaluationError("Dual-polarity shadow incremented the compound count")

    anchor_counts = summary.get("reliable_anchor_counts") or {}
    if not isinstance(anchor_counts, Mapping):
        raise EvaluationError("summary.reliable_anchor_counts must be an object")
    anchors = {
        polarity: _nonnegative_int(
            anchor_counts.get(polarity), field=f"summary.reliable_anchor_counts.{polarity}"
        )
        for polarity in ("POS", "NEG")
    }
    anchors["total"] = anchors["POS"] + anchors["NEG"]

    return {
        "evaluation_status": summary.get("status"),
        "anchors": anchors,
        "edges": {
            "compatible": _summary_count(
                summary, "candidate_independent_compatible_edge_count"
            ),
            "accepted": _summary_count(
                summary, "candidate_independent_accepted_edge_count"
            ),
        },
        "ambiguity": {
            "edge_count": _summary_count(summary, "ambiguous_edge_count"),
            "anchor_count": _summary_count(summary, "ambiguous_anchor_count"),
            "not_evaluable_anchor_count": _summary_count(
                summary, "not_evaluable_anchor_count"
            ),
        },
        "entities": {
            "neutral_entity_count": _summary_count(summary, "neutral_entity_count"),
            "accepted_neutral_entity_count": _summary_count(
                summary, "accepted_neutral_entity_count"
            ),
            "precursor_event_evidence_record_count": _summary_count(
                summary, "precursor_event_evidence_record_count"
            ),
        },
        "dual_candidate_support": {
            "connectivity_count": _summary_count(
                summary, "dual_polarity_supported_connectivity_count"
            )
        },
        "openms_validation": _openms_projection(metadata, pair),
        "v15_production_unchanged_invariants": {
            **production,
            "compound_count_increment": compound_increment,
            "source_file_sha256_unchanged_after_build": True,
            "task_compound_counts_unchanged": True,
            "all_verified": True,
        },
    }


def _evaluate_pair(
    pair: Mapping[str, Any],
    *,
    base_url: str,
    timeout: float,
    request_json: RequestJSON,
) -> dict[str, Any]:
    started = time.monotonic()
    manifest_pair = dict(pair)
    post_payload = {field: pair[field] for field in POST_FIELDS}
    try:
        envelope = request_json(
            method="POST",
            url=f"{base_url}{PAIR_ENDPOINT}",
            payload=post_payload,
            timeout=timeout,
        )
        if envelope.get("status") != "available":
            raise EvaluationError("POST did not return an available sidecar")
        pair_id = str(envelope.get("pair_id") or "")
        if not PAIR_ID_PATTERN.fullmatch(pair_id):
            raise EvaluationError("POST returned an invalid pair_id")
        post_sidecar = envelope.get("sidecar")
        if not isinstance(post_sidecar, Mapping):
            raise EvaluationError("POST response lacks the sidecar")
        post_hash = str(
            (post_sidecar.get("metadata") or {}).get("sidecar_semantic_sha256")
            or ""
        )
        if not SHA256_PATTERN.fullmatch(post_hash):
            raise EvaluationError("POST sidecar lacks a declared semantic hash")

        # GET verifies the declared semantic hash inside the V16 service.  This
        # evaluator intentionally does not duplicate the service hash algorithm.
        verified = request_json(
            method="GET",
            url=f"{base_url}{PAIR_ENDPOINT}/{pair_id}",
            payload=None,
            timeout=timeout,
        )
        if verified.get("status") != "available" or verified.get("pair_id") != pair_id:
            raise EvaluationError("Service-side sidecar verification failed")
        verified_sidecar = verified.get("sidecar")
        if not isinstance(verified_sidecar, Mapping):
            raise EvaluationError("Verified GET response lacks the sidecar")
        verified_hash = str(
            (verified_sidecar.get("metadata") or {}).get(
                "sidecar_semantic_sha256"
            )
            or ""
        )
        if verified_hash != post_hash:
            raise EvaluationError("POST and service-verified GET hashes differ")
        metrics = _metrics_from_sidecar(verified_sidecar, pair)
        return {
            **manifest_pair,
            "status": "succeeded",
            "pair_id": pair_id,
            "reused": envelope.get("reused") is True,
            "artifact": envelope.get("artifact"),
            "sidecar_semantic_hash_verification": {
                "authority": "v16_service_get_endpoint",
                "verified": True,
                "sidecar_semantic_sha256": verified_hash,
                "locally_recomputed": False,
            },
            "metrics": metrics,
            "failure": None,
            "duration_seconds": round(time.monotonic() - started, 6),
        }
    except Exception as exc:
        failure: dict[str, Any] = {
            "type": type(exc).__name__,
            "message": str(exc),
        }
        status_code = getattr(exc, "status_code", None)
        if status_code is not None:
            failure["http_status_code"] = status_code
        return {
            **manifest_pair,
            "status": "failed",
            "pair_id": None,
            "reused": None,
            "sidecar_semantic_hash_verification": {
                "authority": "v16_service_get_endpoint",
                "verified": False,
                "sidecar_semantic_sha256": None,
                "locally_recomputed": False,
            },
            "metrics": None,
            "failure": failure,
            "duration_seconds": round(time.monotonic() - started, 6),
        }


def _empty_aggregate() -> dict[str, Any]:
    return {
        "anchors": {"POS": 0, "NEG": 0, "total": 0},
        "edges": {"compatible": 0, "accepted": 0},
        "ambiguity": {
            "edge_count": 0,
            "anchor_count": 0,
            "not_evaluable_anchor_count": 0,
        },
        "entities": {
            "neutral_entity_count": 0,
            "accepted_neutral_entity_count": 0,
            "precursor_event_evidence_record_count": 0,
        },
        "dual_candidate_support": {"connectivity_count": 0},
        "openms_validation": {
            polarity: {field: 0 for field in OPENMS_COUNT_FIELDS}
            for polarity in ("POS", "NEG", "total")
        },
    }


def _aggregate_pair_results(results: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    aggregate = _empty_aggregate()
    succeeded = [result for result in results if result.get("status") == "succeeded"]
    for result in succeeded:
        metrics = result["metrics"]
        for field in ("POS", "NEG", "total"):
            aggregate["anchors"][field] += metrics["anchors"][field]
        for field in ("compatible", "accepted"):
            aggregate["edges"][field] += metrics["edges"][field]
        for field in (
            "edge_count",
            "anchor_count",
            "not_evaluable_anchor_count",
        ):
            aggregate["ambiguity"][field] += metrics["ambiguity"][field]
        for field in (
            "neutral_entity_count",
            "accepted_neutral_entity_count",
            "precursor_event_evidence_record_count",
        ):
            aggregate["entities"][field] += metrics["entities"][field]
        aggregate["dual_candidate_support"]["connectivity_count"] += metrics[
            "dual_candidate_support"
        ]["connectivity_count"]
        for polarity in ("POS", "NEG", "total"):
            for field in OPENMS_COUNT_FIELDS:
                aggregate["openms_validation"][polarity][field] += metrics[
                    "openms_validation"
                ][polarity][field]
    return aggregate


def evaluate_pairs(
    manifest: Mapping[str, Any],
    *,
    base_url: str,
    concurrency: int = 1,
    timeout: float = 7200.0,
    request_json: RequestJSON = _request_json,
) -> dict[str, Any]:
    """Evaluate a validated manifest and return a compact audit payload."""

    if concurrency not in {1, 2}:
        raise EvaluationError("concurrency must be 1 or 2")
    if timeout <= 0:
        raise EvaluationError("timeout must be positive")
    base_url = _normalise_base_url(base_url)
    pairs = validate_manifest(manifest)
    generated_at = _utc_now()
    service = _service_preflight(
        base_url, timeout=min(timeout, 60.0), request_json=request_json
    )

    if concurrency == 1:
        results = [
            _evaluate_pair(
                pair,
                base_url=base_url,
                timeout=timeout,
                request_json=request_json,
            )
            for pair in pairs
        ]
    else:
        ordered: list[dict[str, Any] | None] = [None] * len(pairs)
        with ThreadPoolExecutor(max_workers=concurrency) as executor:
            futures = {
                executor.submit(
                    _evaluate_pair,
                    pair,
                    base_url=base_url,
                    timeout=timeout,
                    request_json=request_json,
                ): index
                for index, pair in enumerate(pairs)
            }
            for future in as_completed(futures):
                ordered[futures[future]] = future.result()
        results = [result for result in ordered if result is not None]

    successful_count = sum(result["status"] == "succeeded" for result in results)
    failed_count = len(results) - successful_count
    reused_count = sum(result.get("reused") is True for result in results)
    production_verified = all(
        result.get("status") == "succeeded"
        and result["metrics"]["v15_production_unchanged_invariants"][
            "all_verified"
        ]
        is True
        for result in results
    )
    evaluation_status_counts: dict[str, int] = {}
    for result in results:
        if result.get("status") != "succeeded":
            continue
        status = str(result["metrics"].get("evaluation_status") or "missing")
        evaluation_status_counts[status] = evaluation_status_counts.get(status, 0) + 1
    failed_pairs = [
        {
            "pair_label": result["pair_label"],
            "positive_task_id": result["positive_task_id"],
            "negative_task_id": result["negative_task_id"],
            "failure": result["failure"],
        }
        for result in results
        if result.get("status") == "failed"
    ]
    return {
        "schema_version": AUDIT_SCHEMA,
        "generated_at": generated_at,
        "completed_at": _utc_now(),
        "manifest": {
            "manifest_id": manifest.get("manifest_id"),
            "schema_version": manifest.get("schema_version"),
            "semantic_sha256": _canonical_sha256(manifest),
            "selection_method": manifest.get("selection_method"),
            "explicit_pair_count": len(pairs),
            "cohort_pair_counts": manifest.get("cohort_pair_counts"),
            "candidate_ground_truth_included": False,
        },
        "service": service,
        "execution": {
            "base_url": base_url,
            "endpoint": PAIR_ENDPOINT,
            "concurrency": concurrency,
            "timeout_seconds": timeout,
            "pair_order_preserved_in_audit": True,
            "failure_isolation": "per_pair_continue",
            "candidate_ground_truth_loaded": False,
            "task_mutation_endpoint_called": False,
            "sidecar_semantic_hash_authority": "v16_service_get_endpoint",
        },
        "summary": {
            "pair_count": len(results),
            "succeeded_pair_count": successful_count,
            "failed_pair_count": failed_count,
            "reused_pair_count": reused_count,
            "new_sidecar_pair_count": successful_count - reused_count,
            "all_manifest_pairs_verified": production_verified,
            "evaluation_status_counts": dict(sorted(evaluation_status_counts.items())),
            "failed_pairs": failed_pairs,
            "aggregate": _aggregate_pair_results(results),
        },
        "v15_production_unchanged_invariants": {
            "required_values": {
                field: False for field in PRODUCTION_UNCHANGED_FIELDS
            },
            "compound_count_increment_required": 0,
            "successful_pairs_checked": successful_count,
            "failed_pairs_not_evaluable": failed_count,
            "all_successful_pairs_verified": all(
                result["metrics"]["v15_production_unchanged_invariants"][
                    "all_verified"
                ]
                is True
                for result in results
                if result.get("status") == "succeeded"
            ),
            "all_manifest_pairs_verified": production_verified,
        },
        "pairs": results,
    }


def run(
    *,
    manifest_path: Path,
    output_path: Path,
    base_url: str,
    concurrency: int = 1,
    timeout: float = 7200.0,
    request_json: RequestJSON = _request_json,
) -> dict[str, Any]:
    manifest = _load_json(manifest_path)
    audit = evaluate_pairs(
        manifest,
        base_url=base_url,
        concurrency=concurrency,
        timeout=timeout,
        request_json=request_json,
    )
    audit["manifest"]["path"] = str(manifest_path.expanduser().resolve())
    audit["output_path"] = str(output_path.expanduser().resolve())
    _atomic_write_json(output_path, audit)
    return audit


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Replay explicitly declared V16 dual-polarity shadow pairs and "
            "write an atomic JSON audit."
        )
    )
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--base-url", default="http://127.0.0.1:8516")
    parser.add_argument(
        "--concurrency",
        type=int,
        choices=(1, 2),
        default=1,
        help="Sequential by default; at most two pairs may run concurrently.",
    )
    parser.add_argument(
        "--timeout-seconds",
        type=float,
        default=7200.0,
        help="Per-pair HTTP timeout.",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        audit = run(
            manifest_path=args.manifest,
            output_path=args.output,
            base_url=args.base_url,
            concurrency=args.concurrency,
            timeout=args.timeout_seconds,
        )
    except EvaluationError as exc:
        print(f"V16 dual-polarity evaluation failed before replay: {exc}")
        return 2
    summary = audit["summary"]
    print(
        "V16 dual-polarity evaluation: "
        f"{summary['succeeded_pair_count']}/{summary['pair_count']} succeeded; "
        f"{summary['failed_pair_count']} failed; "
        f"audit={audit['output_path']}"
    )
    return 1 if summary["failed_pair_count"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
