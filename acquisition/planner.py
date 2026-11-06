"""Deterministic materialization of EviAcquire-MS shadow plans."""

from __future__ import annotations

import hashlib
import json
import os
from collections.abc import Mapping, Sequence
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from acquisition.action_space import candidate_actions
from acquisition.ambiguity_graph import build_acquisition_needs
from acquisition.information_gain import score_action
from acquisition.observable_catalog import build_observable_catalog
from acquisition.resources import (
    RESOURCE_GATED_ACTIONS,
    allocate_resource_actions,
    public_resource_summary,
    resolve_acquisition_resources,
)
from matching.evidiscrim import results_sha256, source_row_sha256
from schemas.acquisition import AcquisitionPlan, EvidenceDelta
from schemas.evidiscrim import EviDiscrimArtifact


STORE_SCHEMA_VERSION = "1.2.0"


def _canonical(value: Any) -> bytes:
    return json.dumps(
        value, ensure_ascii=False, sort_keys=True, separators=(",", ":"),
        default=str,
    ).encode("utf-8")


def _sha256(value: Any) -> str:
    return hashlib.sha256(_canonical(value)).hexdigest()


def plan_store_path(task_dir: Path, task_id: str) -> Path:
    return task_dir / f"{task_id}_acquisition_plans.json"


def _write_atomic(path: Path, payload: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp-{os.getpid()}")
    with temporary.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, ensure_ascii=False, separators=(",", ":"))
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)


def read_plan_store(path: Path) -> dict[str, Any]:
    if not path.is_file():
        return {
            "schema_version": STORE_SCHEMA_VERSION,
            "metadata": {"status": "not_generated", "plan_count": 0},
            "plans": [],
            "followups": [],
        }
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError("Acquisition plan store must be a JSON object")
    payload.setdefault("plans", [])
    payload.setdefault("followups", [])
    payload.setdefault("metadata", {})
    return _sanitize_resource_gated_actions(payload)


def _sanitize_resource_gated_actions(payload: dict[str, Any]) -> dict[str, Any]:
    """Fail closed for legacy plans that never proved resource availability."""

    sanitized_plans: list[dict[str, Any]] = []
    omitted = 0
    for raw_plan in payload.get("plans") or []:
        if not isinstance(raw_plan, Mapping):
            continue
        plan = dict(raw_plan)
        actions = []
        for raw_action in plan.get("ranked_actions") or []:
            if not isinstance(raw_action, Mapping):
                continue
            action = dict(raw_action)
            kind = str(action.get("action_kind") or "")
            if kind in RESOURCE_GATED_ACTIONS and not (
                action.get("resource_status") == "explicitly_available"
                and action.get("required_resources")
                and action.get("resource_identity_keys")
            ):
                omitted += 1
                continue
            actions.append(action)
        if not actions:
            continue
        plan["ranked_actions"] = actions
        sanitized_plans.append(plan)
    payload["plans"] = sanitized_plans
    metadata = dict(payload.get("metadata") or {})
    metadata["plan_count"] = len(sanitized_plans)
    metadata["planned_actions"] = sum(
        len(plan.get("ranked_actions") or []) for plan in sanitized_plans
    )
    resource_gate = dict(metadata.get("resource_gate") or {})
    if not resource_gate:
        resource_gate = {
            "mode": "sample_only_fail_closed",
            "default_fail_closed": True,
            "authentic_standard_inventory_count": 0,
            "authentic_standard_action_budget": 0,
            "authentic_standard_actions_planned": 0,
            "ion_mobility_available": False,
            "ccs_reference_count": 0,
            "ion_mobility_action_budget": 0,
            "ion_mobility_actions_planned": 0,
            "omitted_action_kinds": sorted(RESOURCE_GATED_ACTIONS),
            "sample_only_primary_workflow": True,
        }
    if omitted:
        resource_gate["legacy_unverified_resource_actions_omitted"] = omitted
    metadata["resource_gate"] = resource_gate
    payload["metadata"] = metadata
    return payload


def _library_snapshot_id(project_root: Path) -> str | None:
    try:
        manifest = json.loads(
            (project_root / "config" / "runtime_manifest.json").read_text(
                encoding="utf-8"
            )
        )
    except (OSError, ValueError):
        return None
    return (
        ((manifest.get("library") or {}).get("production_full") or {}).get(
            "snapshot_id"
        )
    )


def _validated_evidiscrim_context(
    results: Sequence[Mapping[str, Any]],
    task_id: str,
    artifact: Mapping[str, Any] | None,
) -> tuple[dict[int, dict[str, Any]], dict[str, Any], str | None]:
    """Fail closed when a sidecar is stale, malformed, or not pure shadow."""

    if not isinstance(artifact, Mapping):
        return {}, {
            "status": "not_provided",
            "planning_context_used": False,
            "spectra_available": 0,
        }, None
    raw_digest = _sha256(artifact)
    try:
        validated = EviDiscrimArtifact.model_validate(artifact)
    except Exception as exc:
        return {}, {
            "status": "rejected_invalid_artifact",
            "planning_context_used": False,
            "spectra_available": 0,
            "reason": f"{type(exc).__name__}: {exc}"[:500],
        }, raw_digest
    if not validated.schema_version.startswith("2."):
        return {}, {
            "status": "rejected_legacy_evidiscrim_uncertainty_contract",
            "planning_context_used": False,
            "spectra_available": len(validated.spectra),
            "reason": (
                "EviAcquire requires schema v2 independent-architecture "
                "uncertainty; singleton-floor artifacts are audit-only"
            ),
        }, raw_digest
    if validated.task_id != task_id:
        return {}, {
            "status": "rejected_task_id_mismatch",
            "planning_context_used": False,
            "spectra_available": 0,
        }, raw_digest
    full_digest_matches = validated.source_results_sha256 == results_sha256(results)
    if not full_digest_matches:
        candidate_limit = int(
            (validated.summary or {}).get("maximum_candidates_per_spectrum") or 5
        )
        row_fingerprints_match = bool(validated.spectra) and all(
            spectrum.row_index < len(results)
            and spectrum.source_row_sha256 == source_row_sha256(
                results[spectrum.row_index],
                max_candidates_per_spectrum=candidate_limit,
            )
            for spectrum in validated.spectra
        )
        if not row_fingerprints_match:
            return {}, {
                "status": "rejected_stale_source_results",
                "planning_context_used": False,
                "spectra_available": 0,
            }, raw_digest
    if validated.status != "completed_shadow":
        return {}, {
            "status": f"not_usable_{validated.status}",
            "planning_context_used": False,
            "spectra_available": len(validated.spectra),
        }, raw_digest
    index = {
        spectrum.row_index: spectrum.model_dump(mode="json")
        for spectrum in validated.spectra
    }
    return index, {
        "status": (
            "accepted_shadow_planning_context"
            if full_digest_matches
            else "accepted_source_row_fingerprint_matched_shadow_context"
        ),
        "planning_context_used": False,
        "spectra_available": len(index),
        "candidate_conditioned": True,
        "independent_measurement": False,
        "production_ranking_changed": False,
    }, raw_digest


def generate_acquisition_plans(
    task: Mapping[str, Any],
    results: Sequence[Mapping[str, Any]],
    *,
    task_dir: Path,
    project_root: Path,
    runtime_snapshot_id: str | None,
    policy: Mapping[str, Any],
    max_plans: int = 25,
    max_actions_per_plan: int = 6,
    refresh: bool = False,
    evidiscrim_artifact: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    task_id = str(task.get("id") or "")
    if not task_id:
        raise ValueError("Task ID is required")
    path = plan_store_path(task_dir, task_id)
    baseline_digest = _sha256([
        {
            "spectrum_id": row.get("spectrum_id"),
            "precursor_mz": row.get("precursor_mz"),
            "retention_time": row.get("retention_time"),
            "selected_identity": (
                (row.get("match") or {}).get("inchikey")
                or (row.get("match") or {}).get("compound_name")
            ),
            "rank_score": (row.get("match") or {}).get("rank_score"),
        }
        for row in results if isinstance(row, Mapping) and row.get("match")
    ])
    resources = resolve_acquisition_resources(task, policy)
    resource_digest = str(resources["resource_context_sha256"])
    evidiscrim_index, evidiscrim_status, evidiscrim_digest = (
        _validated_evidiscrim_context(results, task_id, evidiscrim_artifact)
    )
    existing = read_plan_store(path)
    if (
        not refresh
        and existing.get("plans")
        and (existing.get("metadata") or {}).get("baseline_result_sha256")
        == baseline_digest
        and (existing.get("metadata") or {}).get("policy_id")
        == policy.get("policy_id")
        and (existing.get("metadata") or {}).get("resource_context_sha256")
        == resource_digest
        and (existing.get("metadata") or {}).get("evidiscrim_artifact_sha256")
        == evidiscrim_digest
    ):
        existing["cache_hit"] = True
        return existing

    needs = build_acquisition_needs(
        task_id, results, max_needs=max_plans
    )
    resource_allocations = allocate_resource_actions(needs, resources)
    generated_at = datetime.now(timezone.utc).isoformat()
    library_snapshot_id = _library_snapshot_id(project_root)
    plans: list[AcquisitionPlan] = []
    evidiscrim_plans_used = 0
    for need in needs:
        row = results[need.source_result_index]
        catalog = build_observable_catalog(
            row,
            evidiscrim_spectrum=evidiscrim_index.get(need.source_result_index),
        )
        actions = [
            score_action(need, row, proposal, catalog, policy)
            for proposal in candidate_actions(
                need,
                row,
                policy,
                context=task,
                enabled_resource_actions=resource_allocations.get(
                    need.need_id, frozenset()
                ),
                resolved_resources=resources,
            )
        ]
        actions.sort(key=lambda action: (
            -action.utility,
            -action.minimum_pair_separation,
            action.duty_cycle_cost,
            action.action_kind,
        ))
        actions = actions[:max_actions_per_plan]
        evidiscrim_used = any(
            any("EviDiscrim-MS" in rationale for rationale in action.rationale)
            for action in actions
        )
        evidiscrim_plans_used += int(evidiscrim_used)
        unsupported = tuple(
            edge.edge_id for edge in need.ambiguity_edges
            if max(
                (action.pairwise_separation.get(edge.edge_id, 0.0) for action in actions),
                default=0.0,
            ) < float(policy.get("minimum_actionable_separation", 0.25))
        )
        plan_id = "plan_" + hashlib.sha256(_canonical([
            task_id,
            need.entity_id,
            policy.get("policy_id"),
            baseline_digest,
            resource_digest,
        ])).hexdigest()[:18]
        plans.append(AcquisitionPlan(
            plan_id=plan_id,
            policy_id=str(policy.get("policy_id") or "eviacquire-ms-unversioned"),
            task_id=task_id,
            entity_id=need.entity_id,
            generated_at=generated_at,
            runtime_snapshot_id=runtime_snapshot_id,
            library_snapshot_id=library_snapshot_id,
            baseline_result_sha256=baseline_digest,
            need=need,
            ranked_actions=tuple(actions),
            unsupported_pairs=unsupported,
            provenance={
                "planner": "candidate-pair maximin active evidence",
                "candidate_conditioned_observables_used_for": "planning_only",
                "candidate_independent_followup_linkage": True,
                "source_result_index": need.source_result_index,
                "source_task_filename": task.get("filename"),
                "resource_mode": resources.get("mode"),
                "resource_context_sha256": resource_digest,
                "evidiscrim_planning_context_used": evidiscrim_used,
                "evidiscrim_candidate_conditioned": bool(evidiscrim_used),
                "evidiscrim_independent_measurement": False,
            },
        ))
    plans.sort(key=lambda plan: (
        -plan.ranked_actions[0].utility,
        -len(plan.need.ambiguity_edges),
        plan.need.retention_time_min
        if plan.need.retention_time_min is not None else 1e9,
        plan.need.precursor_mz,
    ))
    energy_actions = [
        action
        for plan in plans
        for action in plan.ranked_actions
        if action.action_kind == "targeted_ms2_multi_ce"
        and action.collision_energy_design is not None
    ]
    conditional_stages: dict[str, set[tuple[float, ...]]] = {}
    for action in energy_actions:
        for stage in action.collision_energy_design.stages[1:]:
            conditional_stages.setdefault(stage.stage_id, set()).add(
                stage.collision_energies
            )
    energy_summary = {
        "profile_ids": sorted({
            action.collision_energy_design.profile_id for action in energy_actions
        }),
        "initial_energy_schedules": [
            list(values) for values in sorted({
                action.collision_energies for action in energy_actions
            })
        ],
        "conditional_energy_schedules": {
            stage_id: [list(values) for values in sorted(schedules)]
            for stage_id, schedules in sorted(conditional_stages.items())
        },
        "single_energy_no_response_curve_plan_count": sum(
            action.collision_energy_design.evidence_status
            == "single_energy_only_no_response_curve"
            for action in energy_actions
        ),
        "universal_complete_fragmentation_assumed": False,
        "cross_instrument_conversion_applied": False,
    }
    standard_actions_planned = sum(
        action.action_kind == "authentic_standard_rt"
        for plan in plans
        for action in plan.ranked_actions
    )
    mobility_actions_planned = sum(
        action.action_kind == "ion_mobility_ccs"
        for plan in plans
        for action in plan.ranked_actions
    )
    payload = {
        "schema_version": STORE_SCHEMA_VERSION,
        "metadata": {
            "status": "available",
            "task_id": task_id,
            "policy_id": policy.get("policy_id"),
            "generated_at": generated_at,
            "plan_count": len(plans),
            "planned_actions": sum(len(plan.ranked_actions) for plan in plans),
            "ambiguity_edge_count": sum(
                len(plan.need.ambiguity_edges) for plan in plans
            ),
            "baseline_result_sha256": baseline_digest,
            "resource_context_sha256": resource_digest,
            "evidiscrim_artifact_sha256": evidiscrim_digest,
            "evidiscrim_context": {
                **evidiscrim_status,
                "planning_context_used": bool(evidiscrim_plans_used),
                "plans_using_context": evidiscrim_plans_used,
                "score_semantics": (
                    "uncertainty_adjusted_spectral_distribution_divergence_"
                    "for_planning_not_probability"
                ),
            },
            "resource_gate": public_resource_summary(
                resources,
                standard_actions_planned=standard_actions_planned,
                mobility_actions_planned=mobility_actions_planned,
            ),
            "production_ranking_changed": False,
            "predicted_evidence_confirmatory": False,
            "legacy_results_mutated": False,
            "module_state": "shadow",
            "collision_energy_design": energy_summary,
        },
        "plans": [plan.model_dump(mode="json") for plan in plans],
        "followups": existing.get("followups") or [],
    }
    _write_atomic(path, payload)
    payload["artifact"] = str(path.resolve())
    payload["cache_hit"] = False
    return payload


def get_plan(path: Path, plan_id: str) -> dict[str, Any]:
    store = read_plan_store(path)
    for plan in store.get("plans") or []:
        if str(plan.get("plan_id") or "") == plan_id:
            return plan
    raise KeyError(plan_id)


def append_followup(path: Path, delta: EvidenceDelta) -> dict[str, Any]:
    store = read_plan_store(path)
    followups = list(store.get("followups") or [])
    payload = delta.model_dump(mode="json")
    followups = [item for item in followups if item.get("delta_id") != delta.delta_id]
    followups.append(payload)
    store["followups"] = followups
    store.setdefault("metadata", {})["followups_linked"] = len(followups)
    _write_atomic(path, store)
    return payload


def plan_store_summary(path: Path) -> dict[str, Any]:
    store = read_plan_store(path)
    metadata = store.get("metadata") or {}
    return {
        "requested": path.is_file(),
        "applicable": bool(metadata.get("plan_count")),
        "plan_count": int(metadata.get("plan_count") or 0),
        "planned_actions": int(metadata.get("planned_actions") or 0),
        "followups_linked": len(store.get("followups") or []),
        "resource_gate": metadata.get("resource_gate") or {},
        "module_state": "shadow",
        "production_ranking_changed": False,
    }
