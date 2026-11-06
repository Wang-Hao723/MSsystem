#!/usr/bin/env python3
"""Compare the stored V6 network, hydrated V6, and V7 evidence network.

This is a reproducible diagnostic utility, not a production analysis path. It
loads one saved task, restores full MS2 peaks from the paired mzML by spectrum
ID, and reports topology/peak-coverage differences.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import sys
import time
from copy import deepcopy
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from matching.molecular_network import build_network as build_v7_network
from utils.mzml_io import parse_mzml_file


def _load_legacy_builder(path: Path):
    spec = importlib.util.spec_from_file_location("mssystem_v6_network_audit", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load legacy network module: {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.build_network


def _hydrate(results: list[dict], mzml_path: Path) -> dict:
    spectra = parse_mzml_file(mzml_path)
    by_id = {
        str(spectrum.get("spectrum_id")): spectrum.get("peaks") or []
        for spectrum in spectra
    }
    hydrated = 0
    missing = []
    for row in results:
        spectrum_id = str(row.get("spectrum_id"))
        peaks = by_id.get(spectrum_id)
        if peaks:
            row["raw_peaks"] = [
                [float(peak[0]), float(peak[1])] for peak in peaks
            ]
            hydrated += 1
        else:
            missing.append(spectrum_id)
    return {
        "parsed_ms2": len(spectra),
        "hydrated_results": hydrated,
        "missing_result_ids": missing[:20],
        "complete": hydrated == len(results),
    }


def _run(label: str, builder, results: list[dict], **kwargs) -> dict:
    started = time.perf_counter()
    network = builder(results, **kwargs)
    elapsed = time.perf_counter() - started
    stats = network.get("stats") or {}
    return {
        "label": label,
        "elapsed_sec": round(elapsed, 3),
        "n_nodes": stats.get("n_nodes"),
        "n_edges": stats.get("n_edges"),
        "n_families_reported": stats.get("n_families"),
        "n_molecular_families": stats.get(
            "n_molecular_families", (
                (stats.get("n_families") or 0) - (stats.get("orphan_count") or 0)
            )
        ),
        "singletons": stats.get("singleton_count", stats.get("orphan_count")),
        "connected_nodes": stats.get("connected_node_count"),
        "largest_family": stats.get("largest_family"),
        "threshold_edges": stats.get("n_edges_threshold"),
        "after_topk_edges": stats.get("n_edges_after_topk"),
        "removed_topk_edges": stats.get("n_edges_removed_topk"),
        "raw_peak_coverage_pct": stats.get("raw_peak_coverage_pct"),
        "spectral_edges": stats.get("spectral_edges"),
        "analog_edges": stats.get("analog_edges"),
        "neutral_loss_edges": stats.get("neutral_loss_edges"),
        "ion_identity_edges": stats.get("ion_identity_edges"),
        "topk_policy": stats.get("topk_policy"),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", required=True, type=Path)
    parser.add_argument("--mzml", required=True, type=Path)
    parser.add_argument(
        "--legacy-module",
        type=Path,
        default=Path("/Users/zizi/Desktop/Mssystem_v6/matching/molecular_network.py"),
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    results = json.loads(args.results.read_text(encoding="utf-8"))
    if not isinstance(results, list):
        raise ValueError("results JSON must be a list")
    stored = deepcopy(results)
    hydrated = deepcopy(results)
    hydration = _hydrate(hydrated, args.mzml)
    if not hydration["complete"]:
        raise RuntimeError(
            f"Peak hydration incomplete: {hydration['hydrated_results']}/{len(results)}"
        )

    legacy_builder = _load_legacy_builder(args.legacy_module)
    report = {
        "results_path": str(args.results),
        "mzml_path": str(args.mzml),
        "input_results": len(results),
        "stored_raw_peak_rows": sum(bool(row.get("raw_peaks")) for row in stored),
        "hydration": hydration,
        "runs": [
            _run(
                "V6 stored task (matched rows have full peaks only)",
                legacy_builder,
                stored,
            ),
            _run(
                "V6 algorithm with every row hydrated",
                legacy_builder,
                hydrated,
            ),
            _run(
                "V7 feature-consensus multi-evidence network",
                build_v7_network,
                hydrated,
            ),
        ],
    }
    text = json.dumps(report, ensure_ascii=False, indent=2)
    print(text)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
