#!/usr/bin/env python3
"""Profile fixed-position saved queries without reading identities or truth."""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import time


def canonical_sha(value) -> str:
    data = json.dumps(
        value, ensure_ascii=False, sort_keys=True, separators=(",", ":"),
        allow_nan=True, default=str,
    ).encode()
    return hashlib.sha256(data).hexdigest()


def fixed_queries(path: Path, positions: tuple[int, ...]) -> list[dict]:
    import ijson

    wanted = set(positions)
    queries = []
    with path.open("rb") as handle:
        for index, row in enumerate(ijson.items(handle, "item")):
            if index not in wanted:
                if index > positions[-1]:
                    break
                continue
            peaks = row.get("raw_peaks") or row.get("query_peaks") or []
            queries.append({
                "position": index,
                "query_peaks": [[float(p[0]), float(p[1])] for p in peaks],
                "query_pmz": float(row.get("search_precursor_mz") or row["precursor_mz"]),
                "query_polarity": row.get("ion_mode") or row.get("polarity"),
                "query_adduct": row.get("query_adduct") or row.get("adduct"),
                "query_adduct_reliable": row.get("query_adduct_reliable"),
                "query_charge": row.get("precursor_charge"),
                "query_collision_energy": row.get("collision_energy"),
                "query_instrument_type": row.get("instrument_type"),
                "query_instrument": row.get("instrument"),
                "query_activation": row.get("activation_method"),
                "query_collision_energy_mode": row.get("collision_energy_mode"),
                "query_instrument_domain_id": row.get("instrument_domain_id"),
                "query_collision_cell_region": row.get("collision_cell_region"),
                "query_spectrum_representation": row.get("spectrum_representation"),
                "query_raw_peaks_complete": not bool(row.get("limited_fragment_list")),
                "query_rt": row.get("retention_time"),
                "query_rt_unit": row.get("retention_time_unit") or "minutes",
            })
    if [q["position"] for q in queries] != list(positions):
        raise RuntimeError("saved result does not contain every fixed position")
    return queries


class CountingRecords:
    def __init__(self, records):
        self.records = records
        self.requests = 0
        self.uncached_requests = 0

    def __len__(self):
        return len(self.records)

    def __getitem__(self, index):
        self.requests += 1
        return self.records[index]

    def read_uncached(self, index):
        self.uncached_requests += 1
        return self.records.read_uncached(index)

    def __getattr__(self, name):
        return getattr(self.records, name)


def load_proposed(path: Path):
    spec = importlib.util.spec_from_file_location("profile_query_search", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def one_arm(name, search_module, queries, cache_dir, manifest, indices, base, checkpoint=None):
    from matching.library_record_store import MMapRecordStore

    records = CountingRecords(MMapRecordStore(cache_dir, manifest))
    base._lib_records = records
    base._mz_index, base._nm_index = indices
    base._lib_source_families = list(manifest.get("source_families") or [])
    before = records.cache_info()
    outputs = []
    query_rows = []
    started = time.perf_counter()
    for query in queries:
        kwargs = {key: value for key, value in query.items()
                  if key not in {"position", "query_peaks", "query_pmz"}}
        query_start = time.perf_counter()
        result = search_module.search(
            query["query_peaks"], query["query_pmz"],
            records=records, mz_index=indices[0], **kwargs,
        )
        outputs.append(result)
        query_rows.append({
            "position": query["position"],
            "wall_seconds": time.perf_counter() - query_start,
            "result_sha256": canonical_sha(result),
        })
        if checkpoint is not None:
            checkpoint({
                "arm": name,
                "completed_queries": list(query_rows),
                "record_getitem_requests": records.requests,
                "record_uncached_requests": records.uncached_requests,
                "decoded_cache_current": records.cache_info(),
                "elapsed_seconds": time.perf_counter() - started,
            })
    after = records.cache_info()
    row = {
        "arm": name, "wall_seconds": time.perf_counter() - started,
        "record_getitem_requests": records.requests,
        "record_uncached_requests": records.uncached_requests,
        "decoded_cache_before": before, "decoded_cache_after": after,
        "decoded_cache_miss_delta": after["misses"] - before["misses"],
        "decoded_cache_hit_delta": after["hits"] - before["hits"],
        "query_results_sha256": canonical_sha(outputs), "queries": query_rows,
    }
    records.records.close()
    return row


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-root", type=Path, required=True)
    parser.add_argument("--cache-dir", type=Path, required=True)
    parser.add_argument("--saved-results", type=Path, required=True)
    parser.add_argument("--proposed-query-search", type=Path, required=True)
    parser.add_argument("--proposed-record-memo", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--positions", default="0,1048,2096,3144")
    args = parser.parse_args()
    positions = tuple(int(v) for v in args.positions.split(","))
    if tuple(sorted(set(positions))) != positions or not positions:
        raise ValueError("positions must be unique ascending integers")
    sys.path.insert(0, str(args.source_root))
    from matching.library_record_store import MMapRecordStore  # noqa: F401
    from matching.offline_mz_index import MMapMassIndex
    import matching.search as base
    import pipeline.query_search as baseline

    memo_spec = importlib.util.spec_from_file_location(
        "pipeline.query_record_memo", args.proposed_record_memo,
    )
    memo_module = importlib.util.module_from_spec(memo_spec)
    assert memo_spec.loader is not None
    memo_spec.loader.exec_module(memo_module)
    sys.modules["pipeline.query_record_memo"] = memo_module
    proposed = load_proposed(args.proposed_query_search)
    manifest = json.loads((args.cache_dir / "manifest.json").read_text())
    indices = (
        MMapMassIndex(args.cache_dir / "precursor_masses.npy",
                      args.cache_dir / "precursor_record_indices.npy"),
        MMapMassIndex(args.cache_dir / "neutral_masses.npy",
                      args.cache_dir / "neutral_record_indices.npy"),
    )
    queries = fixed_queries(args.saved_results, positions)
    query_input_sha = canonical_sha(queries)
    partial_path = args.output.with_suffix(".partial.json")
    partial_rows = []

    def checkpoint(row):
        current = [item for item in partial_rows if item.get("arm") != row["arm"]]
        current.append(row)
        temporary = partial_path.with_suffix(partial_path.suffix + ".tmp")
        temporary.parent.mkdir(parents=True, exist_ok=True)
        temporary.write_text(json.dumps({
            "schema": "mssystem-query-record-memo-profile-partial/1",
            "positions": positions,
            "query_input_sha256": query_input_sha,
            "arms": current,
        }, indent=2) + "\n")
        temporary.replace(partial_path)
        partial_rows[:] = current

    baseline_row = one_arm(
        "baseline", baseline, queries, args.cache_dir, manifest, indices, base,
        checkpoint=checkpoint,
    )
    proposed_row = one_arm(
        "query_memo", proposed, queries, args.cache_dir, manifest, indices, base,
        checkpoint=checkpoint,
    )
    if baseline_row["query_results_sha256"] != proposed_row["query_results_sha256"]:
        raise RuntimeError("query result bytes differ between profile arms")
    payload = {
        "schema": "mssystem-query-record-memo-profile/1",
        "selection": "fixed saved-result array positions; no identity or truth fields read",
        "positions": positions,
        "query_input_sha256": query_input_sha,
        "source_root": str(args.source_root), "cache_dir": str(args.cache_dir),
        "saved_results": str(args.saved_results),
        "arms": [baseline_row, proposed_row],
        "result_equivalent": True,
        "record_getitem_reduction": baseline_row["record_getitem_requests"] - proposed_row["record_getitem_requests"],
        "decoded_miss_reduction": baseline_row["decoded_cache_miss_delta"] - proposed_row["decoded_cache_miss_delta"],
        "wall_seconds_reduction": baseline_row["wall_seconds"] - proposed_row["wall_seconds"],
    }
    args.output.parent.mkdir(parents=True, exist_ok=False)
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps({
        "result_equivalent": True,
        "baseline_seconds": baseline_row["wall_seconds"],
        "proposed_seconds": proposed_row["wall_seconds"],
        "baseline_requests": baseline_row["record_getitem_requests"],
        "proposed_requests": proposed_row["record_getitem_requests"],
        "baseline_uncached_requests": baseline_row["record_uncached_requests"],
        "proposed_uncached_requests": proposed_row["record_uncached_requests"],
        "baseline_misses": baseline_row["decoded_cache_miss_delta"],
        "proposed_misses": proposed_row["decoded_cache_miss_delta"],
        "output": str(args.output),
    }))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
