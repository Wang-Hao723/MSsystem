#!/usr/bin/env python3
"""Build a frozen AuthenticAnchor-MS artifact from explicit standard tasks."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping


PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from matching.authentic_anchor import AuthenticAnchorPolicy, write_anchor_library


PROTON_MASS = 1.007276466621
SODIUM_ION_MASS = 22.989218
POTASSIUM_ION_MASS = 38.963158
AMMONIUM_ION_MASS = 18.033823
CHLORIDE_ION_MASS = 34.968853
FORMIC_ACID_MASS = 46.005479
ACETIC_ACID_MASS = 60.021129


def _finite(value: Any) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def _ion_hypotheses(neutral_mass: float, polarity: str) -> list[tuple[str, float]]:
    if polarity == "NEG":
        return [
            ("[M-H]-", neutral_mass - PROTON_MASS),
            ("[M+Cl]-", neutral_mass + CHLORIDE_ION_MASS),
            ("[M+HCOO]-", neutral_mass + FORMIC_ACID_MASS - PROTON_MASS),
            ("[M+CH3COO]-", neutral_mass + ACETIC_ACID_MASS - PROTON_MASS),
            ("[M-2H]2-", (neutral_mass - 2 * PROTON_MASS) / 2),
        ]
    return [
        ("[M+H]+", neutral_mass + PROTON_MASS),
        ("[M+Na]+", neutral_mass + SODIUM_ION_MASS),
        ("[M+K]+", neutral_mass + POTASSIUM_ION_MASS),
        ("[M+NH4]+", neutral_mass + AMMONIUM_ION_MASS),
        ("[M+2H]2+", (neutral_mass + 2 * PROTON_MASS) / 2),
    ]


def _target_rows(
    rows: list[dict[str, Any]],
    standard: Mapping[str, Any],
    *,
    tolerance_ppm: float,
) -> list[dict[str, Any]]:
    polarity = str(standard.get("polarity") or "").upper()
    neutral_mass = _finite(standard.get("neutral_mass"))
    if polarity not in {"POS", "NEG"} or neutral_mass is None:
        return []
    hypotheses = _ion_hypotheses(neutral_mass, polarity)
    output = []
    for row in rows:
        precursor = _finite(row.get("search_precursor_mz") or row.get("precursor_mz"))
        peaks = row.get("query_peaks") or row.get("peaks") or []
        if precursor is None or len(peaks) < 3 or row.get("blank_excluded"):
            continue
        ppm, label, theoretical = min(
            (
                abs((precursor - theoretical) / theoretical * 1e6),
                label,
                theoretical,
            )
            for label, theoretical in hypotheses
        )
        if ppm <= tolerance_ppm:
            clean = dict(row)
            clean["anchor_adduct"] = label
            clean["anchor_theoretical_mz"] = theoretical
            clean["anchor_mass_error_abs_ppm"] = ppm
            output.append(clean)
    return output


def _row_signal(row: Mapping[str, Any]) -> float:
    return max(
        _finite(row.get("precursor_ms1_intensity")) or 0.0,
        _finite(row.get("tic")) or 0.0,
        sum(
            _finite(peak[1]) or 0.0
            for peak in (row.get("query_peaks") or row.get("peaks") or ())
            if isinstance(peak, (list, tuple)) and len(peak) >= 2
        ),
    )


def _retention_clusters(
    rows: list[dict[str, Any]], *, tolerance_min: float
) -> list[list[dict[str, Any]]]:
    """Cluster target spectra by adduct and chromatographic retention region."""

    by_adduct: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        rt = _finite(row.get("retention_time"))
        ce = _finite(row.get("collision_energy"))
        if rt is None or ce is None:
            continue
        by_adduct[str(row.get("anchor_adduct") or "unknown")].append(row)
    clusters: list[list[dict[str, Any]]] = []
    for _, group in sorted(by_adduct.items()):
        ordered = sorted(
            group,
            key=lambda row: (
                float(row.get("retention_time")),
                float(row.get("collision_energy")),
                str(row.get("spectrum_id") or ""),
            ),
        )
        for row in ordered:
            rt = float(row["retention_time"])
            if (
                not clusters
                or str(clusters[-1][0].get("anchor_adduct") or "unknown")
                != str(row.get("anchor_adduct") or "unknown")
                or rt - float(clusters[-1][-1]["retention_time"]) > tolerance_min
            ):
                clusters.append([row])
            else:
                clusters[-1].append(row)
    return clusters


def _select_references(
    rows: list[dict[str, Any]],
    max_per_collision_energy: int,
    *,
    rt_cluster_tolerance_min: float,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    clusters = _retention_clusters(
        rows,
        tolerance_min=rt_cluster_tolerance_min,
    )
    if not clusters:
        return [], {
            "status": "no_retention_cluster",
            "candidate_cluster_count": 0,
            "selected_spectrum_count": 0,
        }

    def cluster_strength(cluster: list[dict[str, Any]]) -> tuple[Any, ...]:
        by_ce: dict[float, list[dict[str, Any]]] = defaultdict(list)
        for row in cluster:
            ce = _finite(row.get("collision_energy"))
            if ce is not None:
                by_ce[round(ce, 3)].append(row)
        per_ce_signal = [
            max(_row_signal(row) for row in group)
            for _, group in sorted(by_ce.items())
        ]
        mean_log_signal = (
            sum(math.log1p(value) for value in per_ce_signal) / len(per_ce_signal)
            if per_ce_signal
            else 0.0
        )
        median_rt = sorted(float(row["retention_time"]) for row in cluster)[
            len(cluster) // 2
        ]
        # Multi-energy coherence is a hard first gate.  Within eligible
        # clusters, aggregate log signal identifies the actual chromatographic
        # analyte even when one CE was not triggered; CE count then breaks
        # near-ties.  A single intense background scan can never outrank a
        # coherent multi-energy authentic-standard cluster.
        return (
            len(by_ce) >= 2,
            mean_log_signal,
            len(by_ce),
            min(per_ce_signal, default=0.0),
            sum(len(row.get("query_peaks") or row.get("peaks") or ()) for row in cluster),
            -median_rt,
        )

    selected_cluster = max(clusters, key=cluster_strength)
    grouped: dict[float, list[dict[str, Any]]] = defaultdict(list)
    for row in selected_cluster:
        ce = _finite(row.get("collision_energy"))
        if ce is not None:
            grouped[round(ce, 3)].append(row)
    output = []
    for ce, group in sorted(grouped.items()):
        selected = sorted(
            group,
            key=lambda row: (
                -_row_signal(row),
                _finite(row.get("anchor_mass_error_abs_ppm")) or math.inf,
                str(row.get("spectrum_id") or ""),
            ),
        )[:max_per_collision_energy]
        for row in selected:
            output.append(
                {
                    "source_spectrum_id": row.get("spectrum_id"),
                    "collision_energy_ev": ce,
                    "retention_time_min": _finite(row.get("retention_time")),
                    "precursor_mz": _finite(
                        row.get("search_precursor_mz") or row.get("precursor_mz")
                    ),
                    "adduct": row.get("anchor_adduct"),
                    "mass_error_abs_ppm": round(
                        float(row.get("anchor_mass_error_abs_ppm") or 0.0), 6
                    ),
                    "peaks": row.get("query_peaks") or row.get("peaks") or [],
                }
            )
    selected_rts = [
        float(row["retention_time"])
        for row in selected_cluster
        if _finite(row.get("retention_time")) is not None
    ]
    return output, {
        "status": "coherent_rt_cluster_selected",
        "selection_rule": (
            "require_multi_ce_then_maximize_log_signal_within_one_adduct_and_"
            "rt_cluster"
        ),
        "candidate_cluster_count": len(clusters),
        "selected_adduct": selected_cluster[0].get("anchor_adduct"),
        "selected_rt_cluster_min": round(min(selected_rts), 6),
        "selected_rt_cluster_max": round(max(selected_rts), 6),
        "selected_rt_cluster_center_min": round(
            sum(selected_rts) / len(selected_rts), 6
        ),
        "selected_collision_energies_ev": sorted(grouped),
        "selected_spectrum_count": len(output),
        "discarded_mass_compatible_spectrum_count": len(rows) - len(output),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument(
        "--task-dir", type=Path, default=PROJECT_ROOT / "output" / "tasks"
    )
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--tolerance-ppm", type=float, default=10.0)
    parser.add_argument("--max-per-collision-energy", type=int, default=3)
    parser.add_argument("--rt-cluster-tolerance-min", type=float, default=0.08)
    parser.add_argument("--allow-missing", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    manifest_bytes = args.manifest.read_bytes()
    manifest = json.loads(manifest_bytes)
    if manifest.get("selection_policy") != "explicit_task_id_only":
        raise ValueError("manifest selection_policy must be explicit_task_id_only")
    if manifest.get("candidate_generation_allowed") is not False:
        raise ValueError("manifest must forbid candidate generation")
    standards = manifest.get("standards") or []
    if not standards:
        raise ValueError("manifest contains no standards")
    entries = []
    rejected = []
    for standard in standards:
        task_id = str(standard.get("source_task_id") or "").strip()
        if not task_id:
            rejected.append({"compound_code": standard.get("compound_code"), "reason": "missing_task_id"})
            continue
        result_path = args.task_dir / f"{task_id}_results.json"
        if not result_path.is_file():
            rejected.append({"compound_code": standard.get("compound_code"), "reason": "missing_result"})
            continue
        rows = json.loads(result_path.read_text(encoding="utf-8"))
        targets = _target_rows(rows, standard, tolerance_ppm=args.tolerance_ppm)
        spectra, reference_qc = _select_references(
            targets,
            args.max_per_collision_energy,
            rt_cluster_tolerance_min=args.rt_cluster_tolerance_min,
        )
        if not spectra:
            rejected.append({"compound_code": standard.get("compound_code"), "reason": "no_usable_target_ms2"})
            continue
        entries.append(
            {
                "reference_id": str(
                    standard.get("reference_id")
                    or f"{standard.get('compound_code')}:{standard.get('polarity')}"
                ),
                "compound_code": standard.get("compound_code"),
                "compound_name": standard.get("compound_name"),
                "formula": standard.get("formula"),
                "inchikey": standard.get("inchikey"),
                "neutral_mass": _finite(standard.get("neutral_mass")),
                "polarity": str(standard.get("polarity") or "").upper(),
                "instrument_family": standard.get("instrument_family"),
                "method_id": standard.get("method_id"),
                "column_mode": standard.get("column_mode"),
                "authentic_standard": True,
                "source_task_id": task_id,
                "source_filename": standard.get("source_filename"),
                "source_standard_lot": standard.get("source_standard_lot"),
                "reference_qc": reference_qc,
                "spectra": spectra,
            }
        )
    if rejected and not args.allow_missing:
        raise RuntimeError(
            f"{len(rejected)} standards could not be built; rerun with --allow-missing only for an audited partial artifact"
        )
    policy = AuthenticAnchorPolicy().validated()
    artifact = {
        "library_id": manifest.get("library_id"),
        "created_at": datetime.now(timezone.utc).isoformat(),
        "runtime_mode": "shadow_only",
        "construction_policy": {
            "explicit_task_id_only": True,
            "candidate_generation_allowed": False,
            "authentic_library_shadow_search_allowed": True,
            "truth_labels_used_at_query_time": False,
            "single_standard_injections_are_reference_only": True,
        },
        "release_policy": policy.__dict__,
        "source_manifest_path": str(args.manifest.resolve()),
        "source_manifest_sha256": hashlib.sha256(manifest_bytes).hexdigest(),
        "entry_count": len(entries),
        "rejected_count": len(rejected),
        "build_audit": {
            "requested_standard_count": len(standards),
            "built_standard_count": len(entries),
            "rejected": rejected,
            "precursor_tolerance_ppm": args.tolerance_ppm,
            "max_references_per_collision_energy": args.max_per_collision_energy,
            "rt_cluster_tolerance_min": args.rt_cluster_tolerance_min,
            "reference_selection_rule": (
                "coherent_multi_ce_adduct_and_rt_cluster_before_per_ce_signal_"
                "selection"
            ),
        },
        "entries": entries,
    }
    write_anchor_library(args.output, artifact)
    print(json.dumps({
        "output": str(args.output.resolve()),
        "entry_count": len(entries),
        "rejected_count": len(rejected),
        "runtime_mode": "shadow_only",
    }, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
