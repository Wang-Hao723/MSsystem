#!/usr/bin/env python3
"""Generate the V16 runtime manifest and inherited frozen baselines."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from runtime.golden import build_golden_baseline
from runtime.manifest import (
    _runtime_snapshot_id,
    build_runtime_manifest,
    refresh_runtime_manifest_source,
    sha256_file,
    write_json_atomic,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--manifest",
        type=Path,
        default=PROJECT_ROOT / "config" / "runtime_manifest.json",
    )
    parser.add_argument(
        "--archive-manifest",
        type=Path,
        help="Optional second path for an immutable phase checkpoint copy.",
    )
    parser.add_argument(
        "--golden",
        type=Path,
        default=PROJECT_ROOT / "config" / "golden_baseline.json",
    )
    parser.add_argument(
        "--task-sets",
        type=Path,
        default=PROJECT_ROOT / "config" / "golden_task_sets.json",
    )
    parser.add_argument(
        "--tasks-index",
        type=Path,
        default=PROJECT_ROOT / "output" / "tasks_index.json",
    )
    parser.add_argument(
        "--results-root",
        type=Path,
        default=PROJECT_ROOT / "output" / "tasks",
    )
    parser.add_argument("--health-url", default="http://127.0.0.1:8516/api/health")
    parser.add_argument("--restore-archive", type=Path)
    parser.add_argument(
        "--summary-only",
        action="store_true",
        help="Skip per-spectrum semantic digests (faster but weaker).",
    )
    parser.add_argument(
        "--source-only",
        action="store_true",
        help=(
            "Refresh source/environment/catalog identity while preserving the "
            "previously verified read-only library snapshot."
        ),
    )
    parser.add_argument(
        "--rebind-library",
        action="store_true",
        help=(
            "With --source-only, explicitly refresh the verified read-only "
            "SQLite library snapshot after a controlled database repair."
        ),
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.source_only:
        if args.rebind_library:
            try:
                previous = json.loads(args.manifest.read_text(encoding="utf-8"))
            except (OSError, TypeError, ValueError):
                previous = {}
            manifest = build_runtime_manifest(PROJECT_ROOT)
            for key in (
                "golden_baseline_refresh",
                "screening_coverage_refresh",
                "restore_archive",
            ):
                if previous.get(key) is not None:
                    manifest[key] = previous[key]
            golden_path = PROJECT_ROOT / "config/golden_baseline.json"
            golden_payload = (
                json.loads(golden_path.read_text(encoding="utf-8"))
                if golden_path.is_file()
                else {}
            )
            manifest["golden_baseline"] = {
                "path": str(golden_path.resolve()),
                "exists": golden_path.is_file(),
                "sha256": sha256_file(golden_path) if golden_path.is_file() else None,
                "schema_version": golden_payload.get("schema_version"),
                "task_count": len(golden_payload.get("tasks") or {}),
            }
            if isinstance(manifest.get("golden_baseline_refresh"), dict):
                manifest["golden_baseline_refresh"] = {
                    **manifest["golden_baseline_refresh"],
                    "reason": (
                        "V16 output/ does not contain the frozen V9-derived "
                        "golden task artifacts"
                    ),
                }
            if isinstance(manifest.get("screening_coverage_refresh"), dict):
                manifest["screening_coverage_refresh"] = {
                    **manifest["screening_coverage_refresh"],
                    "coverage_artifact": str(
                        (
                            PROJECT_ROOT
                            / "config/screening_catalog_coverage.json"
                        ).resolve()
                    ),
                }
            manifest["refresh_scope"] = (
                "source_environment_catalogs_and_explicit_read_only_library_rebind"
            )
            manifest["library"]["startup_verification"] = (
                "explicit_read_only_snapshot_rebind_after_controlled_repair"
            )
            manifest["runtime_snapshot_id"] = _runtime_snapshot_id(manifest)
            write_json_atomic(args.manifest, manifest)
        else:
            manifest = refresh_runtime_manifest_source(PROJECT_ROOT, args.manifest)
        if args.archive_manifest is not None:
            write_json_atomic(args.archive_manifest, manifest)
        print(f"runtime_snapshot_id={manifest['runtime_snapshot_id']}")
        print(f"source_tree_sha256={manifest['source']['source_tree_sha256']}")
        print(
            "compound_role_catalog="
            f"{manifest['role_catalogs']['compound_roles'].get('catalog_id')} "
            f"external_contaminants="
            f"{manifest['role_catalogs']['compound_roles'].get('external_contaminant_count')}"
        )
        print(
            "screening_catalog="
            f"{manifest['role_catalogs']['screening'].get('catalog_id')} "
            f"identities="
            f"{manifest['role_catalogs']['screening'].get('indexed_full_inchikey_count')}"
        )
        return 0
    baseline = build_golden_baseline(
        args.tasks_index,
        args.task_sets,
        args.results_root,
        include_semantic_digests=not args.summary_only,
    )
    write_json_atomic(args.golden, baseline)
    manifest = build_runtime_manifest(
        PROJECT_ROOT,
        health_url=args.health_url,
        restore_archive=args.restore_archive,
    )
    manifest["golden_baseline"] = {
        "path": str(args.golden.resolve()),
        "schema_version": baseline["schema_version"],
        "task_count": len(baseline["tasks"]),
    }
    write_json_atomic(args.manifest, manifest)
    if args.archive_manifest is not None:
        write_json_atomic(args.archive_manifest, manifest)
    print(f"runtime_snapshot_id={manifest['runtime_snapshot_id']}")
    print(f"source_tree_sha256={manifest['source']['source_tree_sha256']}")
    print(f"library_snapshot_id={manifest['library']['production_full'].get('snapshot_id')}")
    print(f"golden_tasks={len(baseline['tasks'])}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
