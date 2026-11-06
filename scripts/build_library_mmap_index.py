#!/usr/bin/env python3
"""Build content-validated precursor/neutral mass mmap indices for a library."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from config.settings import DB_PATH, DB_PATH_SMALL
from matching.offline_mz_index import build_library_mass_indices
from matching.search import load_library


def _snapshot_id() -> str:
    payload = json.loads(
        (PROJECT_ROOT / "config" / "runtime_manifest.json").read_text(encoding="utf-8")
    )
    library = payload.get("library") or {}
    production = library.get("production_full") or {}
    return str(production.get("snapshot_id") or library.get("snapshot_id") or "unknown")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db-mode", choices=("small", "full"), default="full")
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()
    db_path = DB_PATH_SMALL if args.db_mode == "small" else DB_PATH
    output_dir = args.output_dir or (
        PROJECT_ROOT / "output" / "library_indices" / f"{args.db_mode}-mmap-v1"
    )
    records, _ = load_library(db_path)
    manifest = build_library_mass_indices(
        records,
        output_dir,
        library_snapshot_id=_snapshot_id(),
        source_db=db_path,
        complete_library=True,
    )
    print(json.dumps(manifest, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
