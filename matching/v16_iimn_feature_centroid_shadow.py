"""Candidate-independent V16 shadow for OpenMS-centred IIMN XICs.

The production IIMN implementation remains unchanged.  This module rebuilds
the same candidate-independent network twice: first with the saved MS2
``precursor_mz`` as the XIC target, then with ``lcms_feature_mz`` only when the
row records a unique OpenMS feature binding.  Missing, ambiguous, or invalid
feature bindings retain ``precursor_mz``.

Only a locked projection of measured precursor, RT, polarity, and OpenMS
feature-link fields is retained.  Library matches, candidate names, formulas,
structures, truth labels, and scores never reach the network builder.
"""

from __future__ import annotations

import copy
import hashlib
import json
import math
from collections import defaultdict
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from matching.ion_identity import (
    _feature_rows,
    annotate_ion_identity_networks,
    extract_targeted_ms1_traces,
)


POLICY_SCHEMA_VERSION = (
    "mssystem.v16.iimn-feature-centroid-shadow.policy/1.0.0"
)
SIDECAR_SCHEMA_VERSION = (
    "mssystem.v16.iimn-feature-centroid-shadow.sidecar/1.0.0"
)
MODULE_ID = "v16_iimn_feature_centroid_shadow"
ALGORITHM_ID = "mssystem-v16-iimn-feature-centroid-shadow-1.0.0"
DEFAULT_POLICY_PATH = (
    Path(__file__).resolve().parent.parent
    / "config"
    / "v16_iimn_feature_centroid_shadow_policy.json"
)

_FALSE_EXECUTION_FIELDS = (
    "runtime_effect",
    "production_results_mutated",
    "production_ranking_changed",
    "production_confidence_changed",
    "production_compound_count_changed",
    "automatic_activation_allowed",
)


def _canonical(value: Any) -> bytes:
    return json.dumps(
        value,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
        default=str,
    ).encode("utf-8")


def semantic_sha256(value: Any) -> str:
    return hashlib.sha256(_canonical(value)).hexdigest()


def file_sha256(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _finite(value: Any) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def _normalise_polarity(value: Any) -> str:
    text = str(value or "").strip().upper()
    return {
        "POSITIVE": "POS",
        "+": "POS",
        "NEGATIVE": "NEG",
        "-": "NEG",
    }.get(text, text if text in {"POS", "NEG"} else "")


def _validate_policy(policy: Mapping[str, Any]) -> dict[str, Any]:
    selected = copy.deepcopy(dict(policy))
    if selected.get("schema_version") != POLICY_SCHEMA_VERSION:
        raise ValueError("V16 IIMN feature-centroid policy schema is not locked")
    if selected.get("module_id") != MODULE_ID:
        raise ValueError("V16 IIMN feature-centroid module id is not locked")
    if selected.get("algorithm_id") != ALGORITHM_ID:
        raise ValueError("V16 IIMN feature-centroid algorithm id is not locked")

    execution = selected.get("execution") or {}
    if execution.get("mode") != "shadow_only":
        raise ValueError("V16 IIMN feature-centroid execution must remain shadow-only")
    for field in _FALSE_EXECUTION_FIELDS:
        if execution.get(field) is not False:
            raise ValueError(f"V16 IIMN execution invariant changed: {field}")

    projection = selected.get("candidate_independent_projection") or {}
    if projection.get("candidate_fields_projected") is not False:
        raise ValueError("Candidate fields cannot enter the IIMN shadow projection")
    for field in (
        "candidate_identity_used_for_trace_selection",
        "candidate_identity_used_for_network_construction",
    ):
        if projection.get(field) is not False:
            raise ValueError(f"Candidate-independent invariant changed: {field}")
    allowed = projection.get("allowed_result_fields") or []
    required = {
        "precursor_mz",
        "retention_time",
        "ion_mode",
        "lcms_feature_id",
        "precursor_feature_id",
        "lcms_feature_mz",
        "lcms_feature_link_status",
    }
    if not required.issubset(set(allowed)):
        raise ValueError("Candidate-independent result projection is incomplete")

    trace = selected.get("trace_center") or {}
    expected = {
        "current_source": "precursor_mz",
        "shadow_source": "lcms_feature_mz",
        "required_link_status": "linked_openms_feature",
        "feature_id_required": True,
        "finite_positive_feature_mz_required": True,
        "fallback_source": "precursor_mz",
        "fallback_when_missing_or_nonunique": True,
    }
    for field, value in expected.items():
        if trace.get(field) != value:
            raise ValueError(f"V16 trace-centre invariant changed: {field}")
    if (_finite(trace.get("mz_tolerance_da_floor")) or 0) <= 0:
        raise ValueError("Trace-centre absolute tolerance must be positive")
    if (_finite(trace.get("mz_tolerance_ppm")) or 0) <= 0:
        raise ValueError("Trace-centre ppm tolerance must be positive")

    network = selected.get("network_construction") or {}
    for field in (
        "neutral_mass_hypothesis_input_unchanged",
        "pairing_rt_input_unchanged",
        "ion_species_table_unchanged",
        "shape_correlation_threshold_unchanged",
        "overlap_threshold_unchanged",
        "apex_rt_threshold_unchanged",
    ):
        if network.get(field) is not True:
            raise ValueError(f"V16 IIMN network invariant changed: {field}")
    return selected


def load_v16_iimn_feature_centroid_shadow_policy(
    path: str | Path | None = None,
) -> dict[str, Any]:
    selected = Path(path) if path is not None else DEFAULT_POLICY_PATH
    return _validate_policy(json.loads(selected.read_text(encoding="utf-8")))


def project_candidate_independent_rows(
    rows: Sequence[Mapping[str, Any]],
    *,
    policy: Mapping[str, Any] | None = None,
) -> list[dict[str, Any]]:
    """Return the locked result projection consumed by this shadow module."""

    selected = _validate_policy(
        policy or load_v16_iimn_feature_centroid_shadow_policy()
    )
    allowed = tuple(
        selected["candidate_independent_projection"]["allowed_result_fields"]
    )
    output = []
    for row in rows:
        if not isinstance(row, Mapping):
            continue
        output.append(
            {
                field: copy.deepcopy(row.get(field))
                for field in allowed
                if row.get(field) is not None
            }
        )
    return output


def _task_projection(task: Mapping[str, Any], polarity: str) -> dict[str, Any]:
    workflow = task.get("untargeted_workflow") or {}
    if not isinstance(workflow, Mapping):
        workflow = {}
    return {
        "task_id": str(task.get("id") or task.get("task_id") or ""),
        "status": str(task.get("status") or ""),
        "polarity": polarity,
        "result_run_id": str(task.get("result_run_id") or ""),
        "feature_xml": str(
            workflow.get("feature_xml") or workflow.get("output_features") or ""
        ),
    }


def _validate_task_and_rows(
    task: Mapping[str, Any], rows: Sequence[Mapping[str, Any]]
) -> tuple[str, str]:
    task_id = str(task.get("id") or task.get("task_id") or "").strip()
    if not task_id:
        raise ValueError("A stable task id is required")
    if str(task.get("status") or "") != "completed":
        raise ValueError("The IIMN feature-centroid shadow requires a completed task")
    polarity = _normalise_polarity(task.get("ion_mode"))
    if polarity not in {"POS", "NEG"}:
        raise ValueError("Task ion_mode must be POS or NEG")
    for row in rows:
        observed = {
            value
            for field in ("ion_mode", "polarity")
            if (value := _normalise_polarity(row.get(field)))
        }
        if any(value != polarity for value in observed):
            raise ValueError("A projected row conflicts with task ion_mode")
    return task_id, polarity


def _trace_center_projection(
    rows: Sequence[Mapping[str, Any]], policy: Mapping[str, Any]
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], dict[str, int]]:
    trace = policy["trace_center"]
    required_status = str(trace["required_link_status"])
    floor = float(trace["mz_tolerance_da_floor"])
    ppm = float(trace["mz_tolerance_ppm"])
    targets = copy.deepcopy(list(rows))
    audit = []
    counts = defaultdict(int)

    for row_index, (source, target) in enumerate(zip(rows, targets)):
        precursor = _finite(source.get("precursor_mz"))
        feature_mz = _finite(source.get("lcms_feature_mz"))
        feature_id = str(source.get("lcms_feature_id") or "").strip()
        link_status = str(source.get("lcms_feature_link_status") or "").strip()

        use_feature = bool(
            link_status == required_status
            and feature_id
            and feature_mz is not None
            and feature_mz > 0
        )
        if use_feature:
            selected = feature_mz
            selected_source = "lcms_feature_mz"
            fallback_reason = None
            target["precursor_mz"] = feature_mz
            counts["rows_using_feature_centroid"] += 1
        else:
            selected = precursor
            selected_source = "precursor_mz"
            if link_status != required_status:
                fallback_reason = "feature_binding_missing_or_nonunique"
            elif not feature_id:
                fallback_reason = "linked_feature_id_missing"
            else:
                fallback_reason = "linked_feature_mz_invalid"
            counts["rows_using_precursor_fallback"] += 1
            counts[f"fallback_{fallback_reason}"] += 1

        delta_da = (
            feature_mz - precursor
            if feature_mz is not None and precursor is not None
            else None
        )
        delta_ppm = (
            delta_da / feature_mz * 1e6
            if delta_da is not None and feature_mz not in {None, 0.0}
            else None
        )
        tolerance = (
            max(floor, abs(precursor) * ppm * 1e-6)
            if precursor is not None
            else None
        )
        outside = bool(
            use_feature
            and delta_da is not None
            and tolerance is not None
            and abs(delta_da) > tolerance
        )
        if outside:
            counts["feature_centroid_outside_current_xic_window"] += 1
        audit.append(
            {
                "row_index": row_index,
                "spectrum_id": str(source.get("spectrum_id") or ""),
                "precursor_feature_id": str(
                    source.get("precursor_feature_id") or ""
                ),
                "lcms_feature_id": feature_id,
                "lcms_feature_link_status": link_status,
                "lcms_feature_candidate_count": source.get(
                    "lcms_feature_candidate_count"
                ),
                "current_trace_center_mz": precursor,
                "feature_centroid_mz": feature_mz,
                "shadow_trace_center_mz": selected,
                "shadow_trace_center_source": selected_source,
                "fallback_reason": fallback_reason,
                "feature_minus_precursor_da": (
                    round(delta_da, 10) if delta_da is not None else None
                ),
                "feature_minus_precursor_ppm": (
                    round(delta_ppm, 6) if delta_ppm is not None else None
                ),
                "current_xic_tolerance_da": (
                    round(tolerance, 10) if tolerance is not None else None
                ),
                "feature_centroid_outside_current_xic_window": outside,
                "candidate_fields_used": False,
            }
        )
    counts["input_rows"] = len(rows)
    return targets, audit, dict(sorted(counts.items()))


def _trace_projection(traces: Mapping[str, Mapping[str, Any]]) -> dict[str, Any]:
    output = {}
    for feature_id, raw_trace in sorted(
        ((str(key), value) for key, value in traces.items()),
        key=lambda item: item[0],
    ):
        trace = raw_trace or {}
        output[feature_id] = {
            "rt": list(trace.get("rt") or []),
            "intensity": list(trace.get("intensity") or []),
            "target_mz": trace.get("target_mz"),
            "points": trace.get("points"),
            "nonzero_points": trace.get("nonzero_points"),
            "apex_intensity": trace.get("apex_intensity"),
            "apex_rt": trace.get("apex_rt"),
            "status": trace.get("status"),
        }
    return output


def _network_projection(
    rows: Sequence[Mapping[str, Any]], summary: Mapping[str, Any]
) -> dict[str, Any]:
    groups: dict[str, dict[str, dict[str, Any]]] = defaultdict(dict)
    for row_index, row in enumerate(rows):
        if str(row.get("ion_identity_status") or "") != "networked":
            continue
        network_id = str(row.get("ion_identity_network_id") or "").strip()
        if not network_id:
            continue
        feature_id = str(
            row.get("lcms_feature_id")
            or row.get("precursor_feature_id")
            or row.get("spectrum_id")
            or f"row-{row_index}"
        )
        groups[network_id][feature_id] = {
            "feature_id": feature_id,
            "ion_identity": str(row.get("ion_identity") or ""),
            "neutral_mass": _finite(row.get("ion_identity_neutral_mass")),
        }

    networks = []
    for network_id in sorted(groups):
        members = [groups[network_id][key] for key in sorted(groups[network_id])]
        material = {
            "member_feature_ids": [item["feature_id"] for item in members],
            "ion_identities": sorted(
                {item["ion_identity"] for item in members if item["ion_identity"]}
            ),
            "neutral_masses": sorted(
                {
                    round(float(item["neutral_mass"]), 8)
                    for item in members
                    if item["neutral_mass"] is not None
                }
            ),
        }
        networks.append(
            {
                "network_id": network_id,
                **material,
                "network_fingerprint_sha256": semantic_sha256(material),
                "candidate_fields_used": False,
            }
        )
    summary_fields = (
        "status",
        "feature_count",
        "usable_trace_count",
        "network_count",
        "networked_features",
        "collapsed_ion_forms",
        "edge_count",
        "adduct_edges",
        "in_source_fragment_edges",
        "runtime_influence",
    )
    projection = {
        "summary": {field: summary.get(field) for field in summary_fields},
        "networks": networks,
        "candidate_fields_used": False,
    }
    projection["network_semantic_sha256"] = semantic_sha256(projection)
    return projection


def _source_path_binding(path: str | Path | None) -> dict[str, Any] | None:
    if path is None:
        return None
    selected = Path(path).expanduser().resolve()
    if not selected.is_file():
        raise FileNotFoundError(f"Shadow source file is unavailable: {selected}")
    return {"path": str(selected), "sha256": file_sha256(selected)}


def build_v16_iimn_feature_centroid_shadow(
    task: Mapping[str, Any],
    results: Sequence[Mapping[str, Any]],
    *,
    source_results_path: str | Path | None = None,
    source_data_path: str | Path | None = None,
    current_traces: Mapping[str, Mapping[str, Any]] | None = None,
    feature_centroid_traces: Mapping[str, Mapping[str, Any]] | None = None,
    policy: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Build a deterministic current-vs-feature-centred IIMN sidecar in memory."""

    selected_policy = _validate_policy(
        policy or load_v16_iimn_feature_centroid_shadow_policy()
    )
    policy_sha = semantic_sha256(selected_policy)
    projected = project_candidate_independent_rows(results, policy=selected_policy)
    projected_before = semantic_sha256(projected)
    task_id, polarity = _validate_task_and_rows(task, projected)
    task_projection = _task_projection(task, polarity)

    result_binding_before = _source_path_binding(source_results_path)
    data_binding_before = _source_path_binding(source_data_path)
    shadow_targets, row_audit, audit_counts = _trace_center_projection(
        projected, selected_policy
    )

    supplied_current = current_traces is not None
    supplied_shadow = feature_centroid_traces is not None
    if supplied_current != supplied_shadow:
        raise ValueError("Current and feature-centroid traces must be supplied together")
    if supplied_current:
        current_trace_map = copy.deepcopy(dict(current_traces or {}))
        shadow_trace_map = copy.deepcopy(dict(feature_centroid_traces or {}))
        current_trace_summary = {"status": "provided_traces"}
        shadow_trace_summary = {"status": "provided_traces"}
    else:
        if source_data_path is None:
            raise ValueError("source_data_path or both trace maps are required")
        current_trace_map, current_trace_summary = extract_targeted_ms1_traces(
            source_data_path,
            _feature_rows(copy.deepcopy(projected)),
        )
        shadow_trace_map, shadow_trace_summary = extract_targeted_ms1_traces(
            source_data_path,
            _feature_rows(copy.deepcopy(shadow_targets)),
        )

    current_rows = copy.deepcopy(projected)
    shadow_rows = copy.deepcopy(projected)
    current_summary = annotate_ion_identity_networks(
        current_rows, traces=current_trace_map
    )
    shadow_summary = annotate_ion_identity_networks(
        shadow_rows, traces=shadow_trace_map
    )
    current_graph = _network_projection(current_rows, current_summary)
    shadow_graph = _network_projection(shadow_rows, shadow_summary)

    current_fingerprints = {
        item["network_fingerprint_sha256"] for item in current_graph["networks"]
    }
    shadow_fingerprints = {
        item["network_fingerprint_sha256"] for item in shadow_graph["networks"]
    }

    projected_after = semantic_sha256(projected)
    result_binding_after = _source_path_binding(source_results_path)
    data_binding_after = _source_path_binding(source_data_path)
    if projected_before != projected_after:
        raise AssertionError("Candidate-independent input projection was mutated")
    if result_binding_before != result_binding_after:
        raise AssertionError("Source result file changed during shadow execution")
    if data_binding_before != data_binding_after:
        raise AssertionError("Source data file changed during shadow execution")

    payload = {
        "schema_version": SIDECAR_SCHEMA_VERSION,
        "metadata": {
            "module_id": MODULE_ID,
            "algorithm_id": ALGORITHM_ID,
            "execution_mode": "shadow_only",
            "task_id": task_id,
            "polarity": polarity,
            "policy_sha256": policy_sha,
            "task_metadata_sha256": semantic_sha256(task_projection),
            "candidate_independent_projection_sha256": projected_before,
        },
        "source_binding": {
            "task": task_projection,
            "results_file": result_binding_before,
            "source_data_file": data_binding_before,
            "current_trace_semantic_sha256": semantic_sha256(
                _trace_projection(current_trace_map)
            ),
            "feature_centroid_trace_semantic_sha256": semantic_sha256(
                _trace_projection(shadow_trace_map)
            ),
            "source_files_unchanged": True,
        },
        "trace_center_audit": {
            "counts": audit_counts,
            "rows": row_audit,
            "current_trace_extraction": current_trace_summary,
            "feature_centroid_trace_extraction": shadow_trace_summary,
        },
        "comparison": {
            "current_precursor_center": current_graph,
            "feature_centroid_center": shadow_graph,
            "delta": {
                "network_count": int(
                    (shadow_summary.get("network_count") or 0)
                    - (current_summary.get("network_count") or 0)
                ),
                "networked_feature_count": int(
                    (shadow_summary.get("networked_features") or 0)
                    - (current_summary.get("networked_features") or 0)
                ),
                "edge_count": int(
                    (shadow_summary.get("edge_count") or 0)
                    - (current_summary.get("edge_count") or 0)
                ),
                "added_network_fingerprints": sorted(
                    shadow_fingerprints - current_fingerprints
                ),
                "removed_network_fingerprints": sorted(
                    current_fingerprints - shadow_fingerprints
                ),
            },
        },
        "invariants": {
            "candidate_independent_projection_only": True,
            "candidate_identity_used_for_trace_selection": False,
            "candidate_identity_used_for_network_construction": False,
            "neutral_mass_hypothesis_input_unchanged": True,
            "pairing_rt_input_unchanged": True,
            "production_results_mutated": False,
            "production_ranking_changed": False,
            "production_confidence_changed": False,
            "production_compound_count_changed": False,
            "automatic_activation_allowed": False,
            "source_files_unchanged": True,
        },
    }
    payload["metadata"]["sidecar_semantic_sha256"] = semantic_sha256(payload)
    return payload
