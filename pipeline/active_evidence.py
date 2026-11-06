"""Pipeline boundary for materializing additive EviAcquire-MS artifacts."""

from __future__ import annotations

import json
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from acquisition.planner import generate_acquisition_plans
from matching.evidiscrim import evidiscrim_artifact_path


def load_acquisition_policy(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict) or not payload.get("policy_id"):
        raise ValueError("Acquisition policy must be a versioned JSON object")
    return payload


def materialize_active_evidence_plans(
    task: Mapping[str, Any],
    results: Sequence[Mapping[str, Any]],
    *,
    task_dir: Path,
    project_root: Path,
    runtime_snapshot_id: str | None,
    max_plans: int,
    max_actions_per_plan: int,
    refresh: bool,
) -> dict[str, Any]:
    policy = load_acquisition_policy(
        project_root / "config" / "acquisition_policy.json"
    )
    evidiscrim_path = evidiscrim_artifact_path(task_dir, str(task.get("id") or ""))
    evidiscrim_artifact: dict[str, Any] | None = None
    if evidiscrim_path.is_file():
        try:
            loaded = json.loads(evidiscrim_path.read_text(encoding="utf-8"))
            if isinstance(loaded, dict):
                evidiscrim_artifact = loaded
        except (OSError, ValueError) as exc:
            # Preserve fail-closed provenance without making ordinary sample-only
            # acquisition planning depend on an optional shadow sidecar.
            evidiscrim_artifact = {
                "status": "invalid_load",
                "load_error": f"{type(exc).__name__}: {exc}"[:500],
            }
    return generate_acquisition_plans(
        task,
        results,
        task_dir=task_dir,
        project_root=project_root,
        runtime_snapshot_id=runtime_snapshot_id,
        policy=policy,
        max_plans=max_plans,
        max_actions_per_plan=max_actions_per_plan,
        refresh=refresh,
        evidiscrim_artifact=evidiscrim_artifact,
    )
