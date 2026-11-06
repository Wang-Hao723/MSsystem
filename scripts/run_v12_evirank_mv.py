#!/usr/bin/env python3
"""Materialize a V12-Flagship / EviRank-MV shadow sidecar."""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from pipeline.evirank_mv import evirank_artifact_path, write_evirank_artifact
from pipeline.result_contract import runtime_snapshot_id_from_manifest


def _task_id(path: Path, explicit: str | None) -> str:
    if explicit:
        value = explicit.strip()
    else:
        value = path.name
        for suffix in ("_v11_flagship.json", ".json"):
            if value.endswith(suffix):
                value = value[: -len(suffix)]
                break
    if not re.fullmatch(r"[A-Za-z0-9_.-]+", value):
        raise ValueError("Task id contains unsupported characters")
    return value


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Build the V12-Flagship / EviRank-MV shadow challenger from V11"
    )
    parser.add_argument("v11_artifact", type=Path, help="Existing *_v11_flagship.json")
    parser.add_argument("--task-id")
    parser.add_argument("--detection-method-id", default="unspecified")
    parser.add_argument("--view-bundle", type=Path)
    parser.add_argument(
        "--candidate-recall-bundle",
        type=Path,
        help="Optional *_v12_candidate_recall.json shadow bundle",
    )
    parser.add_argument("--residual-model", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    source = args.v11_artifact.expanduser().resolve()
    if not source.is_file():
        parser.error(f"V11 artifact does not exist: {source}")
    task_id = _task_id(source, args.task_id)
    output = (
        args.output.expanduser().resolve()
        if args.output
        else evirank_artifact_path(PROJECT_ROOT / "output" / "v12_evirank_mv", task_id)
    )
    metadata = write_evirank_artifact(
        task_id,
        source,
        output,
        detection_method_id=args.detection_method_id,
        view_bundle=args.view_bundle,
        candidate_recall_bundle=args.candidate_recall_bundle,
        residual_model=args.residual_model,
        runtime_snapshot_id=runtime_snapshot_id_from_manifest(PROJECT_ROOT),
    )
    print(json.dumps({
        "input": str(source),
        "output": str(output),
        "algorithm_id": metadata.get("algorithm_id"),
        "algorithm_name": metadata.get("algorithm_name"),
        "deployment_state": metadata.get("deployment_state"),
        "residual_model_state": metadata.get("residual_model_state"),
        "candidate_recall_state": metadata.get("candidate_recall_state"),
        "entities_with_candidate_recall_expansion": metadata.get(
            "entities_with_candidate_recall_expansion"
        ),
        "entities_with_qualified_source_paired_switch": metadata.get(
            "entities_with_qualified_source_paired_switch"
        ),
        "expanded_candidates_entered_shadow_candidate_union": metadata.get(
            "expanded_candidates_entered_shadow_candidate_union"
        ),
        "expanded_candidates_entered_residual_model_scoring": metadata.get(
            "expanded_candidates_entered_residual_model_scoring"
        ),
        "entity_count": metadata.get("entity_count"),
        "entities_reranked_from_v11": metadata.get("entities_reranked_from_v11"),
        "production_ranking_changed": metadata.get("production_ranking_changed"),
        "size_bytes": metadata.get("size_bytes"),
    }, ensure_ascii=False, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
