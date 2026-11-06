#!/usr/bin/env python3
"""Build a provenance-preserving PubMed artifact for one task result file."""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from matching.isomer_literature import (  # noqa: E402
    enrich_isomer_groups_with_pubmed,
    write_json_atomic,
)


def _default_output(results_path: Path) -> Path:
    stem = results_path.name
    if stem.endswith("_results.json"):
        stem = stem[:-len("_results.json")]
    else:
        stem = results_path.stem
    return results_path.with_name(f"{stem}_isomer_literature.json")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Search PubMed for unresolved isomer groups. Bibliographic leads "
            "are stored separately and never modify match scores."
        )
    )
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    parser.add_argument(
        "--cache-dir",
        type=Path,
        default=PROJECT_ROOT / "output" / "isomer_literature_cache",
    )
    parser.add_argument(
        "--max-groups",
        type=int,
        default=0,
        help="0 processes every unique competitive ambiguity group.",
    )
    parser.add_argument("--max-results", type=int, default=5)
    parser.add_argument("--timeout", type=float, default=12.0)
    args = parser.parse_args()

    results_path = args.results.resolve()
    output_path = (
        args.output.resolve()
        if args.output else _default_output(results_path)
    )
    with results_path.open("r", encoding="utf-8") as handle:
        results = json.load(handle)
    if not isinstance(results, list):
        raise ValueError("Task results must be a JSON array")

    def progress(index, total, group, lookup):
        names = " / ".join(group.get("candidate_names") or [])[:100]
        status = lookup.get("status", "unknown")
        print(
            f"[{index}/{total}] {group.get('formula') or 'no-formula'} "
            f"{names}: {status}",
            flush=True,
        )

    artifact = enrich_isomer_groups_with_pubmed(
        results,
        cache_dir=args.cache_dir.resolve(),
        max_groups=max(0, args.max_groups),
        max_results=max(1, min(args.max_results, 10)),
        timeout=args.timeout,
        progress_callback=progress,
    )
    artifact["source_results_file"] = str(results_path)
    artifact["source_results_size_bytes"] = results_path.stat().st_size
    artifact["source_results_sha256"] = _sha256(results_path)
    write_json_atomic(output_path, artifact)
    print(
        json.dumps({
            "output": str(output_path),
            "total_unique_ambiguous_groups": artifact[
                "total_unique_ambiguous_groups"
            ],
            "groups_processed": artifact["groups_processed"],
            "status_counts": artifact["status_counts"],
            "used_for_scoring": artifact["used_for_scoring"],
        }, ensure_ascii=False),
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
