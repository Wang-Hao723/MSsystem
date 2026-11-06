#!/usr/bin/env python3
"""Bounded single-query profile for query-local mmap record memoization."""
from __future__ import annotations
import argparse, hashlib, importlib, importlib.util, json, os, signal, sys, time
from pathlib import Path


def sha(value):
    return hashlib.sha256(json.dumps(value, ensure_ascii=False, sort_keys=True,
        separators=(",", ":"), allow_nan=True, default=str).encode()).hexdigest()


def atomic_json(path, value):
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(value, ensure_ascii=False, indent=2, default=str) + "\n")
    tmp.replace(path)


def fixed_query(path, position):
    import ijson
    with path.open("rb") as handle:
        for index, row in enumerate(ijson.items(handle, "item")):
            if index == position:
                peaks = row.get("raw_peaks") or row.get("query_peaks") or []
                return {
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
                }
            if index > position:
                break
    raise RuntimeError(f"saved result lacks fixed position {position}")


class CountingRecords:
    def __init__(self, records): self.records, self.requests = records, 0
    def __len__(self): return len(self.records)
    def __getitem__(self, index):
        self.requests += 1
        return self.records[index]
    def __getattr__(self, name): return getattr(self.records, name)


def load_file(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--arm", choices=("baseline", "query_memo"), required=True)
    ap.add_argument("--source-root", type=Path, required=True)
    ap.add_argument("--cache-dir", type=Path, required=True)
    ap.add_argument("--saved-results", type=Path, required=True)
    ap.add_argument("--proposed-query-search", type=Path, required=True)
    ap.add_argument("--proposed-record-memo", type=Path, required=True)
    ap.add_argument("--output-dir", type=Path, required=True)
    ap.add_argument("--position", type=int, default=1048)
    ap.add_argument("--search-timeout", type=int, default=120)
    args = ap.parse_args()
    if os.environ.get("PYTHONHASHSEED") is None:
        raise RuntimeError("set PYTHONHASHSEED identically for both profile arms")
    args.output_dir.mkdir(parents=True, exist_ok=False)
    checkpoint_path = args.output_dir / "CHECKPOINT.json"
    state = {"schema": "mssystem-query-record-memo-single/1", "arm": args.arm,
             "position": args.position, "events": []}
    def event(stage, **fields):
        state["events"].append({"stage": stage, "monotonic": time.monotonic(), **fields})
        atomic_json(checkpoint_path, state)

    event("start")
    sys.path.insert(0, str(args.source_root))
    from matching.library_record_store import MMapRecordStore
    from matching.offline_mz_index import MMapMassIndex
    base = importlib.import_module("matching.search")
    baseline = importlib.import_module("pipeline.query_search")
    query = fixed_query(args.saved_results, args.position)
    event("query_loaded", query_input_sha256=sha(query), precursor_mz=query["query_pmz"],
          peak_count=len(query["query_peaks"]), polarity=query["query_polarity"])
    manifest = json.loads((args.cache_dir / "manifest.json").read_text())
    records = CountingRecords(MMapRecordStore(args.cache_dir, manifest))
    mz_index = MMapMassIndex(args.cache_dir / "precursor_masses.npy",
                             args.cache_dir / "precursor_record_indices.npy")
    nm_index = MMapMassIndex(args.cache_dir / "neutral_masses.npy",
                             args.cache_dir / "neutral_record_indices.npy")
    base._lib_records, base._mz_index, base._nm_index = records, mz_index, nm_index
    base._lib_source_families = tuple(manifest.get("source_families") or ())
    tol = max(float(base.MZ_TOLERANCE_DA), query["query_pmz"] * float(base.MZ_TOLERANCE_PPM) / 1e6)
    event("mmap_indices_bound", record_count=len(records),
          coarse_precursor_window_count=len(mz_index.range_query(query["query_pmz"], tol)),
          tolerance_da=tol, process_cache_before=records.cache_info())
    if args.arm == "baseline":
        search_module = baseline
    else:
        memo = load_file("pipeline.query_record_memo", args.proposed_record_memo)
        sys.modules["pipeline.query_record_memo"] = memo
        search_module = load_file("profile_query_search", args.proposed_query_search)

    originals = (base._query_candidates, base._prescreen_dense_candidates)
    def query_candidates_probe(*a, **kw):
        started = time.perf_counter(); answer = originals[0](*a, **kw)
        event("query_candidates", wall_seconds=time.perf_counter()-started,
              returned_count=len(answer), record_requests=records.requests)
        return answer
    def prescreen_probe(*a, **kw):
        started = time.perf_counter(); input_count=len(a[0]); answer = originals[1](*a, **kw)
        event("dense_prescreen", wall_seconds=time.perf_counter()-started,
              input_count=input_count, selected_count=len(answer[0]), record_requests=records.requests)
        return answer
    base._query_candidates, base._prescreen_dense_candidates = query_candidates_probe, prescreen_probe
    kwargs = {k: v for k, v in query.items() if k not in {"position", "query_peaks", "query_pmz"}}
    event("search_started")
    started = time.perf_counter()
    def timeout_handler(_signum, _frame):
        raise TimeoutError(f"search exceeded {args.search_timeout}s")
    previous_alarm = signal.signal(signal.SIGALRM, timeout_handler)
    signal.alarm(args.search_timeout)
    try:
        result = search_module.search(query["query_peaks"], query["query_pmz"],
                                      records=records, mz_index=mz_index, **kwargs)
    except TimeoutError as error:
        event("search_timeout", wall_seconds=time.perf_counter()-started,
              error=str(error), record_requests=records.requests,
              process_cache_current=records.cache_info())
        return 124
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, previous_alarm)
        base._query_candidates, base._prescreen_dense_candidates = originals
    wall = time.perf_counter() - started
    result_path = args.output_dir / "RESULT.json"
    atomic_json(result_path, result)
    state.update({"completed": True, "wall_seconds": wall,
                  "result_sha256": sha(result), "record_requests": records.requests,
                  "process_cache_after": records.cache_info(),
                  "result_path": str(result_path)})
    event("search_completed", wall_seconds=wall, result_sha256=state["result_sha256"],
          record_requests=records.requests, process_cache_after=records.cache_info())
    atomic_json(args.output_dir / "PROFILE.json", state)
    records.records.close()
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
