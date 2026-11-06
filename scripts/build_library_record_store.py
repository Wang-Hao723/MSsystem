#!/usr/bin/env python3
"""Build a new read-only normalized-record cache without a full RAM preload."""
from __future__ import annotations

import argparse
import json
import resource
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from config.settings import DB_PATH  # noqa: E402
from matching.library_record_store import (  # noqa: E402
    build_record_store, canonical_bytes, processing_identity,
)
from matching.search import iter_normalized_library_records  # noqa: E402


def peak_rss_bytes() -> int:
    value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return int(value if sys.platform == "darwin" else value * 1024)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-db", type=Path, default=DB_PATH)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--maximum-rss-mib", type=int, default=8192)
    parser.add_argument("--limit", type=int, help="Offline partial benchmark only; never eligible for runtime")
    args = parser.parse_args()
    if args.maximum_rss_mib <= 0 or (args.limit is not None and args.limit <= 0):
        parser.error("Memory and optional record limits must be positive")
    if args.output_dir.exists():
        parser.error("Output generation already exists; choose a new directory")
    identity = processing_identity()
    started = time.perf_counter()
    def stream():
        iterator = iter_normalized_library_records(args.source_db)
        try:
            for number, record in enumerate(iterator, 1):
                if number % 1000 == 0 and peak_rss_bytes() > args.maximum_rss_mib * 1024**2:
                    raise MemoryError("Owned cache build exceeded its explicit RSS guard")
                yield record
                if args.limit is not None and number >= args.limit:
                    break
        finally:
            iterator.close()
    def progress(count):
        print(f"{count} normalized records; elapsed={time.perf_counter()-started:.1f}s; "
              f"peak_RSS={peak_rss_bytes()/1024**2:.1f}MiB", flush=True)
    manifest = build_record_store(stream(), args.output_dir, source_db=args.source_db,
        identity=identity, complete_library=args.limit is None, progress=progress)
    if processing_identity() != identity:
        raise RuntimeError("Processing code/options changed during cache construction; cache not activated")
    receipt = {"schema_version": "mssystem.library-record-store-build/1.0.0",
        "cache_dir": str(args.output_dir.resolve()), "source_db": str(args.source_db.resolve()),
        "manifest_semantic_sha256": manifest["semantic_sha256"],
        "record_count": manifest["record_count"], "complete_library": manifest["complete_library"],
        "elapsed_seconds": time.perf_counter() - started, "peak_rss_bytes": peak_rss_bytes(),
        "configured_rss_guard_mib": args.maximum_rss_mib,
        "shared_database_modified": False, "runtime_activated": False,
        "cache_bytes": sum(x["size_bytes"] for x in manifest["files"].values())}
    with (args.output_dir / "build_receipt.json").open("xb") as handle:
        handle.write(canonical_bytes(receipt) + b"\n")
    print(json.dumps(receipt, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
