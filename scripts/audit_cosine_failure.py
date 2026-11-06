#!/usr/bin/env python3
"""Read-only diagnostic: WHY do matched MIX33 spectra score cosine < 0.8?

The High-confidence funnel audit (``scripts/audit_high_confidence_funnel.py``)
tells us *how many* matched spectra are blocked from High by ``cosine < 0.8``.
It does not tell us *why*.  Lowering the cosine gate to rescue them would be
wrong: most of them are genuinely wrong candidates, wrong adducts, or sparse
references.  Only a subset (reverse coverage and precursor mass both pass, but
cosine still lags) is plausibly salvageable by better preprocessing / collision
energy handling / intensity transforms.

This script consumes existing task JSON only.  It never opens or mutates the
shared spectral-library database, never changes a score, and never relaxes a
threshold.  Each cosine-failing matched spectrum is given:

* the raw features needed to localise the failure (precursor ppm, reverse
  coverage, forward coverage, library / query peak counts, intensity Pearson);
* a set of *all* applicable failure tags (tags intentionally overlap, because
  the underlying causes are not mutually exclusive);
* a single *primary* failure label, chosen by a documented priority, so the
  aggregate histogram is a clean partition.

The eight primary labels are exactly the failure modes named in the recovery
plan::

    wrong_candidate          mass/adduct_conflict   sparse_reference
    collision_energy_mismatch intensity_mismatch     co_isolation
    peak_truncation          library_coverage_failure

Per-spectrum records are dumped alongside the aggregate so the Phase-2
multi-view ablation can reuse the same feature table without re-parsing peaks.
"""

from __future__ import annotations

import argparse
import json
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path
import sys
from typing import Any, Iterable

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

DEFAULT_TASK_DIR = PROJECT_ROOT / "output" / "tasks"
DEFAULT_INDEX = PROJECT_ROOT / "output" / "tasks_index.json"

# ── Diagnosis thresholds (mirror production gates + the recovery plan) ──
# These are read gates for *classifying* failures, never for changing a score.
COS_FAIL_GATE = 0.80          # analyse matched spectra below this cosine
REV_COVERAGE_GATE = 0.60      # reverse (library->query intensity coverage)
PPM_CONSISTENCY_GATE = 20.0   # precursor mass consistency, |ppm|
EFF_LIB_SPARSE = 4            # <= this many effective reference fragments
CLEAN_CEILING = 15            # query Top15/5% preprocessing ceiling
RAW_TRUNCATION_FLOOR = 15     # raw query had more peaks than the ceiling kept
COISOLATION_FORWARD = 0.50    # forward (query->library) coverage below this
INTENSITY_MISMATCH = 0.50     # intensity_score (matched-fragment Pearson)
MZ_MATCH_TOL = 0.10           # peak-pair tolerance used by production cosine


def _f(value: Any) -> float | None:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return None
    return None if out != out else out


def _effective_reference_count(lib_peaks: list, lib_precursor: float | None) -> int:
    """Product-ion count used as the reference richness proxy.

    The stored ``lib_peaks`` is the candidate's full library spectrum.  A raw
    length <= 4 over-counts "sparse" because it includes the precursor / high-m/z
    noise.  We therefore count fragments strictly below the precursor (the
    interpretable product ions), which matches the "effective fragments" notion
    used in the recovery plan.
    """
    if not lib_peaks:
        return 0
    if lib_precursor is None or lib_precursor <= 0:
        # No precursor anchor: fall back to raw count so we never silently
        # inflate richness.
        return len(lib_peaks)
    cut = lib_precursor - 0.5
    return sum(1 for mz, _ in lib_peaks if _f(mz) is not None and float(mz) < cut and float(_) != 0)


def _coverage_ratio(source: list, target: list, tol: float = MZ_MATCH_TOL) -> float:
    """Intensity-weighted fraction of ``source`` peaks present in ``target``.

    This is the forward analogue of ``reverse_match_ratio``.  We deliberately
    use the *raw* query peaks for the forward direction so co-isolation (many
    sample-only peaks) is not masked by the Top15 ceiling.
    """
    if not source or not target:
        return 0.0
    target_mz = np.asarray([p[0] for p in target], dtype=float)
    if target_mz.size == 0:
        return 0.0
    total = 0.0
    matched = 0.0
    for mz, inten in source:
        ia = abs(float(inten))
        if ia <= 0:
            continue
        total += ia
        if float(np.min(np.abs(target_mz - float(mz)))) <= tol:
            matched += ia
    return float(matched / total) if total > 0 else 0.0


# Primary-label priority: most decisive / least recoverable first.  A spectrum
# carries every applicable tag, but the aggregate histogram reports only this
# single primary label so the partition is clean.
PRIMARY_PRIORITY = [
    "wrong_candidate",
    "mass_adduct_conflict",
    "sparse_reference",
    "library_coverage_failure",
    "co_isolation",
    "peak_truncation",
    "collision_energy_mismatch",
    "intensity_mismatch",
]


@dataclass
class SpectrumDiagnosis:
    task_id: str
    spectrum_id: str
    ion_mode: str
    precursor_mz: float | None
    cosine: float | None
    reverse: float | None
    ppm: float | None
    intensity_score: float | None
    mass_metadata_consistent: bool | None
    query_clean_n: int
    query_raw_n: int
    lib_raw_n: int
    lib_effective_n: int
    forward_ratio_clean: float
    forward_ratio_raw: float
    ce_query: Any
    ce_available: bool
    top12_cos_gap: float | None
    recoverable: bool
    tags: list[str] = field(default_factory=list)
    primary: str = "unclassified"


def _classify(d: SpectrumDiagnosis) -> tuple[list[str], str]:
    tags: list[str] = []
    rev_ok = d.reverse is not None and d.reverse >= REV_COVERAGE_GATE
    rev_bad = d.reverse is not None and d.reverse < REV_COVERAGE_GATE
    ppm_ok = d.ppm is None or abs(d.ppm) <= PPM_CONSISTENCY_GATE
    ppm_bad = d.ppm is not None and abs(d.ppm) > PPM_CONSISTENCY_GATE
    inten_low = d.intensity_score is not None and d.intensity_score < INTENSITY_MISMATCH
    sparse = d.lib_effective_n <= EFF_LIB_SPARSE
    trunc = d.query_clean_n >= CLEAN_CEILING and d.query_raw_n > RAW_TRUNCATION_FLOOR
    # Co-isolation / matrix interference: the library explains itself (reverse
    # high) but the query carries substantial intensity the library cannot
    # explain (forward low).  We use the *cleaned* (Top15/5%) forward ratio,
    # because that is the actual view cosine sees; the raw top-100 view is
    # contaminated by low-intensity noise and over-flags this cause.
    coisol = d.forward_ratio_clean < COISOLATION_FORWARD and rev_ok

    # 1. Precursor mass inconsistent -> wrong adduct or metadata conflict.
    if ppm_bad:
        tags.append("mass_adduct_conflict")
    # 2. Library's characteristic fragments absent from the query.
    #    Combined with a bad mass this is an unambiguous wrong candidate;
    #    with a consistent mass but neither side explains the other it is a
    #    structurally different compound sharing only the precursor.
    if rev_bad and ppm_bad:
        tags.append("wrong_candidate")
    if rev_bad and ppm_ok and d.forward_ratio_raw < COISOLATION_FORWARD and not sparse:
        tags.append("wrong_candidate")
    # 3. Reference spectrum too sparse to support a reliable call.
    if sparse:
        tags.append("sparse_reference")
    # 4. Reverse coverage fails but the reference is rich and mass is fine:
    #    the library spectrum simply does not cover what the sample shows
    #    (incomplete library entry / different sub-structure), not a wrong mass.
    if rev_bad and ppm_ok and not sparse:
        tags.append("library_coverage_failure")
    # 5. Query carries many peaks the library cannot explain while the
    #    library's peaks are largely present -> co-fragmentation / co-isolation.
    if coisol:
        tags.append("co_isolation")
    # 6. The Top15/5% ceiling discarded low-abundance diagnostic ions.
    if trunc:
        tags.append("peak_truncation")
    # 7. Recoverable core (reverse ok, mass ok): the remaining cosine drag is
    #    most plausibly a collision-energy / instrument difference.  We can only
    #    name CE directly when a CE was actually recorded; otherwise the
    #    intensity profile is the only visible proxy.
    if rev_ok and ppm_ok:
        if d.ce_available:
            tags.append("collision_energy_mismatch")
        if inten_low:
            tags.append("intensity_mismatch")
        if not d.ce_available and not inten_low and not coisol and not trunc:
            # No other cause but cosine still lags: treat as CE/instrument
            # mismatch by exclusion (CE chain is broken for MIX33, see report).
            tags.append("collision_energy_mismatch")

    primary = next((label for label in PRIMARY_PRIORITY if label in tags), "unclassified")
    # De-duplicate while preserving priority order.
    seen: set[str] = set()
    ordered = [t for t in tags if not (t in seen or seen.add(t))]
    return ordered, primary


def diagnose_row(task_id: str, row: dict[str, Any]) -> SpectrumDiagnosis | None:
    match = row.get("match")
    if not isinstance(match, dict):
        return None
    cos = _f(match.get("cosine"))
    if cos is None or cos >= COS_FAIL_GATE:
        return None

    raw_peaks = row.get("raw_peaks") or []
    query_peaks = row.get("query_peaks") or []
    lib_peaks = row.get("lib_peaks") or []
    lib_precursor = _f(match.get("lib_precursor_mz"))
    ce_q = match.get("collision_energy")
    ce_available = ce_q not in (None, "", 0, "0", "None")

    top_cands = match.get("top_candidates") or []
    top12_gap = None
    cands_cos = [_f(c.get("cosine")) for c in top_cands if _f(c.get("cosine")) is not None]
    if len(cands_cos) >= 2:
        top12_gap = float(sorted(cands_cos, reverse=True)[0] - sorted(cands_cos, reverse=True)[1])

    diag = SpectrumDiagnosis(
        task_id=task_id,
        spectrum_id=str(row.get("spectrum_id") or ""),
        ion_mode=str(row.get("ion_mode") or ""),
        precursor_mz=_f(row.get("precursor_mz")),
        cosine=cos,
        reverse=_f(match.get("reverse_match")),
        ppm=_f(match.get("mass_ppm")),
        intensity_score=_f(match.get("intensity_score")),
        mass_metadata_consistent=match.get("mass_metadata_consistent"),
        query_clean_n=len(query_peaks),
        query_raw_n=len(raw_peaks),
        lib_raw_n=len(lib_peaks),
        lib_effective_n=_effective_reference_count(lib_peaks, lib_precursor),
        forward_ratio_clean=_coverage_ratio(query_peaks, lib_peaks),
        forward_ratio_raw=_coverage_ratio(raw_peaks, lib_peaks),
        ce_query=ce_q,
        ce_available=ce_available,
        top12_cos_gap=top12_gap,
        recoverable=bool(
            (_f(match.get("reverse_match")) or 0.0) >= REV_COVERAGE_GATE
            and (match.get("mass_ppm") is None or abs(_f(match.get("mass_ppm")) or 0.0) <= PPM_CONSISTENCY_GATE)
        ),
    )
    diag.tags, diag.primary = _classify(diag)
    return diag


def audit_task(task_id: str, result_path: Path) -> list[SpectrumDiagnosis]:
    rows = json.loads(result_path.read_text(encoding="utf-8"))
    out: list[SpectrumDiagnosis] = []
    for row in rows:
        if not isinstance(row, dict):
            continue
        diag = diagnose_row(task_id, row)
        if diag is not None:
            out.append(diag)
    return out


def _iter_selected(index: dict, selected: Iterable[str]):
    selected_set = set(selected)
    for task_id, meta in index.items():
        if selected_set and task_id not in selected_set:
            continue
        if meta.get("status") == "completed" and meta.get("has_results_file"):
            yield task_id, meta


def aggregate(diagnoses: list[SpectrumDiagnosis]) -> dict[str, Any]:
    primary = Counter(d.primary for d in diagnoses)
    recoverable = [d for d in diagnoses if d.recoverable]
    recoverable_primary = Counter(d.primary for d in recoverable)

    # Reproduce the recovery-plan baseline subsets for cross-validation.
    rev_bad = sum(1 for d in diagnoses if d.reverse is not None and d.reverse < REV_COVERAGE_GATE)
    ppm_bad = sum(1 for d in diagnoses if d.ppm is not None and abs(d.ppm) > PPM_CONSISTENCY_GATE)
    both_bad = sum(
        1 for d in diagnoses
        if (d.reverse is not None and d.reverse < REV_COVERAGE_GATE)
        and (d.ppm is not None and abs(d.ppm) > PPM_CONSISTENCY_GATE)
    )
    trunc = sum(
        1 for d in diagnoses
        if d.query_clean_n >= CLEAN_CEILING and d.query_raw_n > RAW_TRUNCATION_FLOOR
    )
    sparse = sum(1 for d in diagnoses if d.lib_effective_n <= EFF_LIB_SPARSE)
    ce_available = sum(1 for d in diagnoses if d.ce_available)

    # Secondary-tag co-occurrence (how often each cause fires, regardless of
    # the primary bucket).
    tag_counts = Counter()
    for d in diagnoses:
        tag_counts.update(d.tags)

    return {
        "cosine_failing_spectra": len(diagnoses),
        "recoverable_spectra": len(recoverable),
        "primary_label_counts": dict(primary.most_common()),
        "recoverable_primary_label_counts": dict(recoverable_primary.most_common()),
        "secondary_tag_counts": dict(tag_counts.most_common()),
        "validation_subsets": {
            "reverse_below_0.60": rev_bad,
            "ppm_over_20": ppm_bad,
            "reverse_below_0.60_and_ppm_over_20": both_bad,
            "peak_truncation": trunc,
            "sparse_reference_effective_le_4": sparse,
            "collision_energy_available": ce_available,
        },
    }


def _render_markdown(payload: dict[str, Any]) -> str:
    agg = payload["aggregate"]
    lines = [
        "# V6 cosine-failure diagnostic (read-only)",
        "",
        "Explains *why* matched MIX33 spectra land below cosine 0.80, without",
        "changing any score or threshold.  Counts are spectrum-level (not",
        "deduplicated to compounds) and match the funnel-audit blocker set.",
        "",
        "## Cross-validation against the recovery-plan baseline",
        "",
        "| Subset | Spectra |",
        "|---|---:|",
    ]
    v = agg["validation_subsets"]
    lines.append(f"| cosine < 0.80 (analysed) | {agg['cosine_failing_spectra']} |")
    lines.append(f"| recoverable (reverse >= 0.60 & |ppm| <= 20) | {agg['recoverable_spectra']} |")
    lines.append(f"| reverse < 0.60 | {v['reverse_below_0.60']} |")
    lines.append(f"| ppm > 20 | {v['ppm_over_20']} |")
    lines.append(f"| reverse < 0.60 & ppm > 20 | {v['reverse_below_0.60_and_ppm_over_20']} |")
    lines.append(f"| peak truncation (Top15 ceiling hit) | {v['peak_truncation']} |")
    lines.append(f"| sparse reference (effective fragments <= 4) | {v['sparse_reference_effective_le_4']} |")
    lines.append(f"| collision energy actually recorded | {v['collision_energy_available']} |")

    lines += ["", "## Primary failure label (clean partition)", "",
              "| Primary cause | Spectra |", "|---|---:|"]
    for label, count in agg["primary_label_counts"].items():
        lines.append(f"| {label} | {count} |")

    lines += ["", "## Recoverable subset — what is holding cosine back", "",
              "These pass reverse coverage and precursor mass, so the candidate is",
              "plausibly correct; the cosine drag is a preprocessing / CE / intensity",
              "effect the Phase-2 multi-view ablation targets.", "",
              "| Primary cause | Spectra |", "|---|---:|"]
    for label, count in agg["recoverable_primary_label_counts"].items():
        lines.append(f"| {label} | {count} |")

    lines += ["", "## Secondary tags (overlapping; a spectrum may carry several)", "",
              "| Tag | Spectra |", "|---|---:|"]
    for label, count in agg["secondary_tag_counts"].items():
        lines.append(f"| {label} | {count} |")
    return "\n".join(lines) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("task_ids", nargs="*",
                        help="Task IDs; default audits all completed tasks")
    parser.add_argument("--task-dir", type=Path, default=DEFAULT_TASK_DIR)
    parser.add_argument("--index", type=Path, default=DEFAULT_INDEX)
    parser.add_argument("--json-output", type=Path)
    parser.add_argument("--markdown-output", type=Path)
    args = parser.parse_args()

    index = json.loads(args.index.read_text(encoding="utf-8"))
    per_task: list[dict[str, Any]] = []
    all_diagnoses: list[SpectrumDiagnosis] = []
    for task_id, meta in _iter_selected(index, args.task_ids):
        result_path = args.task_dir / f"{task_id}_results.json"
        if not result_path.exists():
            continue
        diagnoses = audit_task(task_id, result_path)
        if not diagnoses:
            continue
        per_task.append({
            "task_id": task_id,
            "filename": str(meta.get("filename") or ""),
            "diagnoses": [d.__dict__ for d in diagnoses],
        })
        all_diagnoses.extend(diagnoses)

    payload = {
        "schema_version": "1.0",
        "task_count": len(per_task),
        "thresholds": {
            "cos_fail_gate": COS_FAIL_GATE,
            "reverse_coverage_gate": REV_COVERAGE_GATE,
            "ppm_consistency_gate": PPM_CONSISTENCY_GATE,
            "effective_lib_sparse": EFF_LIB_SPARSE,
            "clean_ceiling": CLEAN_CEILING,
            "coisolation_forward": COISOLATION_FORWARD,
            "intensity_mismatch": INTENSITY_MISMATCH,
        },
        "aggregate": aggregate(all_diagnoses),
        "tasks": per_task,
    }
    if args.json_output:
        args.json_output.parent.mkdir(parents=True, exist_ok=True)
        args.json_output.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
    markdown = _render_markdown(payload)
    if args.markdown_output:
        args.markdown_output.parent.mkdir(parents=True, exist_ok=True)
        args.markdown_output.write_text(markdown, encoding="utf-8")
    print(markdown)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
