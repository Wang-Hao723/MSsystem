#!/usr/bin/env python3
"""Materialize a hash-bound V12 EviRecall-MV shadow sidecar."""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from pipeline.evirank_retrieval import (  # noqa: E402
    evirecall_artifact_path,
    write_evirecall_artifact,
)
from pipeline.result_contract import runtime_snapshot_id_from_manifest  # noqa: E402


def _task_id(path: Path, explicit: str | None) -> str:
    value = str(explicit or "").strip()
    if not value:
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
        description="Build the read-only V12 EviRecall-MV candidate-expansion sidecar"
    )
    parser.add_argument("v11_artifact", type=Path)
    parser.add_argument("legacy_results", type=Path)
    parser.add_argument("--task-id")
    parser.add_argument("--detection-method-id", default="unspecified")
    parser.add_argument(
        "--database",
        type=Path,
        default=Path("/Users/zizi/Desktop/Mssystem_shared_data/local_library.db"),
    )
    parser.add_argument(
        "--runtime-manifest",
        type=Path,
        default=PROJECT_ROOT / "config" / "runtime_manifest.json",
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    v11 = args.v11_artifact.expanduser().resolve()
    source = args.legacy_results.expanduser().resolve()
    database = args.database.expanduser().resolve()
    manifest_path = args.runtime_manifest.expanduser().resolve()
    for path in (v11, source, database, manifest_path):
        if not path.is_file():
            parser.error(f"Required input does not exist: {path}")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    library = manifest["library"]["production_full"]
    if database != Path(library["path"]).resolve():
        parser.error("Database path is not the frozen runtime-manifest library")
    task_id = _task_id(v11, args.task_id)
    output = (
        args.output.expanduser().resolve()
        if args.output
        else evirecall_artifact_path(
            PROJECT_ROOT / "output" / "v12_candidate_recall",
            task_id,
        )
    )
    metadata = write_evirecall_artifact(
        task_id,
        v11,
        source,
        database,
        output,
        library_snapshot_id=str(library["snapshot_id"]),
        library_sampled_sha256=str(library["fingerprint"]["sampled_sha256"]),
        detection_method_id=args.detection_method_id,
        runtime_snapshot_id=runtime_snapshot_id_from_manifest(PROJECT_ROOT),
    )
    print(
        json.dumps(
            {
                "input_v11": str(v11),
                "input_legacy_results": str(source),
                "output": str(output),
                "algorithm_id": metadata.get("algorithm_id"),
                "algorithm_name": metadata.get("algorithm_name"),
                "deployment_state": metadata.get("deployment_state"),
                "library_snapshot_id": metadata.get("library_snapshot_id"),
                "entity_count": metadata.get("entity_count"),
                "entities_with_expansion": metadata.get("entities_with_expansion"),
                "expanded_candidate_count": metadata.get("expanded_candidate_count"),
                "entities_with_qualified_source_paired_switch": metadata.get(
                    "entities_with_qualified_source_paired_switch"
                ),
                "production_effect": metadata.get("production_effect"),
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
