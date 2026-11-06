#!/usr/bin/env python3
"""Benchmark mmap mass lookup on real Full-DB mass density without cutover."""

from __future__ import annotations

import argparse
import json
import os
import sqlite3
import sys
import time
from pathlib import Path

import numpy as np


PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from config.settings import DB_PATH, MZ_BUCKET_DALTON
from matching.offline_mz_index import (
    build_library_mass_indices,
    load_library_mass_indices,
)
from matching.search import _query_candidates


def _bucket_index(records: list[dict]) -> dict[int, list[int]]:
    index: dict[int, list[int]] = {}
    for record_index, record in enumerate(records):
        mass = float(record["precursor_mz"])
        index.setdefault(int(mass / MZ_BUCKET_DALTON), []).append(record_index)
    for values in index.values():
        values.sort(key=lambda i: records[i]["precursor_mz"])
    return index


def _estimated_bucket_bytes(index: dict[int, list[int]]) -> int:
    # CPython integers are normally 28 bytes; the estimate intentionally
    # includes list capacity and dictionary/key overhead but not record data.
    return (
        sys.getsizeof(index)
        + sum(sys.getsizeof(key) + sys.getsizeof(values) for key, values in index.items())
        + sum(28 * len(values) for values in index.values())
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rows", type=int, default=250_000)
    parser.add_argument("--queries", type=int, default=5_000)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=PROJECT_ROOT / "output" / "benchmarks" / "p4_mmap_real_mass",
    )
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    connection = sqlite3.connect(f"file:{DB_PATH.resolve()}?mode=ro", uri=True)
    rows = connection.execute(
        "SELECT spectrum_id, precursor_mz FROM Spectra "
        "WHERE precursor_mz IS NOT NULL AND precursor_mz > 0 LIMIT ?",
        (args.rows,),
    ).fetchall()
    connection.close()
    records = [
        {"spectrum_id": str(spectrum_id), "precursor_mz": float(mass)}
        for spectrum_id, mass in rows
    ]
    rng = np.random.default_rng(args.seed)
    query_count = min(args.queries, len(records))
    selected = rng.choice(len(records), size=query_count, replace=False)
    queries = [records[int(index)]["precursor_mz"] for index in selected]

    started = time.perf_counter()
    bucket = _bucket_index(records)
    bucket_build_seconds = time.perf_counter() - started

    index_dir = args.output_dir / "index"
    started = time.perf_counter()
    manifest = build_library_mass_indices(
        records,
        index_dir,
        library_snapshot_id="raw-full-db-mass-benchmark-not-production",
        source_db=DB_PATH,
        complete_library=False,
    )
    mmap_build_seconds = time.perf_counter() - started
    started = time.perf_counter()
    mmap_index, _, _ = load_library_mass_indices(
        index_dir,
        expected_record_count=len(records),
        expected_record_alignment_sha256=manifest["record_alignment_sha256"],
        require_complete_library=False,
    )
    mmap_open_seconds = time.perf_counter() - started

    started = time.perf_counter()
    bucket_results = [
        _query_candidates(mass, records, bucket, tol_da=0.02, tol_ppm=10.0)
        for mass in queries
    ]
    bucket_query_seconds = time.perf_counter() - started
    started = time.perf_counter()
    mmap_results = [
        _query_candidates(mass, records, mmap_index, tol_da=0.02, tol_ppm=10.0)
        for mass in queries
    ]
    mmap_query_seconds = time.perf_counter() - started
    mismatches = [
        index
        for index, (left, right) in enumerate(zip(bucket_results, mmap_results))
        if left != right
    ]
    mmap_bytes = mmap_index.mapped_bytes
    bucket_bytes = _estimated_bucket_bytes(bucket)
    report = {
        "schema_version": "1.0.0",
        "source_db": str(DB_PATH.resolve()),
        "source_db_size": DB_PATH.stat().st_size,
        "benchmark_rows": len(records),
        "benchmark_queries": query_count,
        "mass_tolerance_da": 0.02,
        "mass_tolerance_ppm": 10.0,
        "bucket_build_seconds": bucket_build_seconds,
        "mmap_build_seconds": mmap_build_seconds,
        "mmap_cold_open_seconds": mmap_open_seconds,
        "bucket_query_seconds": bucket_query_seconds,
        "mmap_query_seconds": mmap_query_seconds,
        "query_speed_ratio_bucket_over_mmap": (
            bucket_query_seconds / mmap_query_seconds if mmap_query_seconds else None
        ),
        "estimated_bucket_index_bytes": bucket_bytes,
        "mmap_index_file_bytes": mmap_bytes,
        "estimated_index_memory_reduction_fraction": (
            1.0 - mmap_bytes / bucket_bytes if bucket_bytes else None
        ),
        "candidate_order_mismatch_count": len(mismatches),
        "candidate_order_equivalent": not mismatches,
        "mismatch_query_indices_sample": mismatches[:20],
        "score_fields_stored_in_index": False,
        "ranking_fields_stored_in_index": False,
        "production_cutover_eligible": False,
        "cutover_blockers": [
            "benchmark_uses_partial_raw_db_mass_rows_not_complete_filtered_runtime_records",
            "complete_full_library_index_not_built",
            "golden_end_to_end_topk_equivalence_not_yet_executed_with_mmap_runtime",
            "full_runtime_rss_and_cold_start_not_yet_measured",
        ],
    }
    output_path = args.output_dir / "benchmark.json"
    temporary = output_path.with_name(f".{output_path.name}.tmp-{os.getpid()}")
    with temporary.open("w", encoding="utf-8") as handle:
        json.dump(report, handle, ensure_ascii=False, indent=2)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, output_path)
    print(json.dumps(report, ensure_ascii=False, indent=2))
    return 0 if not mismatches else 1


if __name__ == "__main__":
    raise SystemExit(main())
