#!/usr/bin/env python3
"""Materialize a V11 flagship sidecar from an existing scan-result artifact."""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from pipeline.flagship import flagship_artifact_path, write_flagship_artifact
from pipeline.result_contract import iter_legacy_result_artifact
from pipeline.result_contract import runtime_snapshot_id_from_manifest


def _task_id(path: Path, explicit: str | None) -> str:
    if explicit:
        value = explicit.strip()
        if not re.fullmatch(r"[A-Za-z0-9_.-]+", value):
            raise ValueError("--task-id contains unsupported characters")
        return value
    name = path.name
    for suffix in ("_results.json", ".json"):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
            break
    value = re.sub(r"[^A-Za-z0-9_.-]+", "-", name).strip("-.")
    if not value:
        raise ValueError("Could not infer a task id from the input path")
    return value


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Build the primary V11 entity-level result from immutable scan candidates"
        )
    )
    parser.add_argument("input", type=Path, help="Existing *_results.json artifact")
    parser.add_argument("--task-id", help="Task identifier; inferred from the filename by default")
    parser.add_argument(
        "--output",
        type=Path,
        help=(
            "Output sidecar. Defaults to V11/output/v11_flagship/<task>_v11_flagship.json"
        ),
    )
    args = parser.parse_args()

    source = args.input.expanduser().resolve()
    if not source.is_file():
        parser.error(f"input does not exist: {source}")
    task_id = _task_id(source, args.task_id)
    output = (
        args.output.expanduser().resolve()
        if args.output
        else flagship_artifact_path(PROJECT_ROOT / "output" / "v11_flagship", task_id)
    )
    metadata = write_flagship_artifact(
        task_id,
        iter_legacy_result_artifact(source),
        output,
        runtime_snapshot_id=runtime_snapshot_id_from_manifest(PROJECT_ROOT),
    )
    print(
        json.dumps(
            {
                "input": str(source),
                "output": str(output),
                "algorithm_id": metadata.get("algorithm_id"),
                "policy_sha256": metadata.get("policy_sha256"),
                "source_row_count": metadata.get("source_row_count"),
                "entity_count": metadata.get("entity_count"),
                "entities_with_candidate": metadata.get("entities_with_candidate"),
                "entities_abstained": metadata.get("entities_abstained"),
                "entities_unmatched": metadata.get("entities_unmatched"),
                "entities_reranked_from_first_scan": metadata.get(
                    "entities_reranked_from_first_scan"
                ),
                "fragmentation_status_counts": metadata.get(
                    "fragmentation_status_counts"
                ),
                "size_bytes": metadata.get("size_bytes"),
            },
            ensure_ascii=False,
            indent=2,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
