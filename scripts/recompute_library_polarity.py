#!/usr/bin/env python3
"""Recompute polarity coverage with the same rule used by runtime search."""
from __future__ import annotations

import argparse
import json
import os
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from config.settings import DB_PATH  # noqa: E402
from matching.library_polarity import (  # noqa: E402
    compound_quality_snapshot,
    recompute_polarity_snapshot,
)


def _atomic_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=str(path.parent)
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, ensure_ascii=False, indent=2)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        Path(temporary).replace(path)
    except Exception:
        Path(temporary).unlink(missing_ok=True)
        raise


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--db", type=Path, default=DB_PATH)
    parser.add_argument(
        "--cache",
        type=Path,
        default=ROOT / "config" / "library_polarity_snapshot.json",
    )
    parser.add_argument(
        "--registry",
        type=Path,
        default=ROOT / "config" / "public_source_registry.json",
    )
    parser.add_argument("--apply", action="store_true")
    args = parser.parse_args()

    snapshot = recompute_polarity_snapshot(args.db)
    structure = compound_quality_snapshot(args.db)
    result = {
        "database_snapshot": snapshot,
        "structure_quality_snapshot": structure,
    }
    if args.apply:
        _atomic_json(args.cache, snapshot)
        registry = json.loads(args.registry.read_text(encoding="utf-8"))
        registry["database_snapshot"] = snapshot
        registry["structure_quality_snapshot"] = structure
        _atomic_json(args.registry, registry)
        result["written"] = [str(args.cache.resolve()), str(args.registry.resolve())]
    else:
        result["mode"] = "dry_run"
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
