#!/usr/bin/env python3
"""Audit where V6 spectra and compounds leave the High-confidence funnel.

The audit consumes existing task JSON only.  It never opens or mutates the
shared spectral-library database.  Counts are reported at both spectrum and
unique-compound level so repeated DDA/CE scans are not confused with chemical
identifications.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from dataclasses import dataclass
import json
from pathlib import Path
import re
import sys
from typing import Any, Iterable

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from matching.honest_calibration import annotate
from matching.library_quality import compound_identity_key, resolve_structure_formula
from config.settings import (
    ISOMER_AMBIGUITY_MAX_COSINE_GAP,
    ISOMER_AMBIGUITY_MAX_RANK_GAP,
)


DEFAULT_TASK_DIR = PROJECT_ROOT / "output" / "tasks"
DEFAULT_INDEX = PROJECT_ROOT / "output" / "tasks_index.json"


def _float(value: Any) -> float | None:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return None if result != result else result


def _identity(match: dict[str, Any]) -> str:
    key = compound_identity_key(match)
    if key:
        return key
    name = re.sub(r"[^a-z0-9]+", "", str(match.get("compound_name") or "").lower())
    formula = str(match.get("formula") or "").strip()
    return f"fallback:{name}:{formula}" if name or formula else "fallback:unknown"


def _same_formula_competitors(match: dict[str, Any]) -> list[dict[str, Any]]:
    formula = resolve_structure_formula(
        str(match.get("smiles") or ""),
        str(match.get("formula") or ""),
        str(match.get("adduct") or ""),
    )
    if not formula:
        return []
    best_id = _identity(match)
    competitors: list[dict[str, Any]] = []
    for candidate in match.get("top_candidates") or []:
        candidate_formula = resolve_structure_formula(
            str(candidate.get("smiles") or ""),
            str(candidate.get("formula") or ""),
            str(candidate.get("adduct") or ""),
        )
        if candidate_formula != formula:
            continue
        if _identity(candidate) == best_id:
            continue
        competitors.append(candidate)
    return competitors


def _best_competitor_gaps(match: dict[str, Any]) -> dict[str, float] | None:
    competitors = _same_formula_competitors(match)
    if not competitors:
        return None

    best_rank = _float(match.get("rank_score"))
    best_final = _float(match.get("final_score"))
    best_cos = _float(match.get("cosine"))
    if best_rank is None or best_final is None or best_cos is None:
        return None

    comp_ranks = [
        _float(item.get("rank_score"))
        if _float(item.get("rank_score")) is not None
        else (_float(item.get("final_score")) or 0.0)
        for item in competitors
    ]
    comp_finals = [_float(item.get("final_score")) or 0.0 for item in competitors]
    comp_cosines = [_float(item.get("cosine")) or 0.0 for item in competitors]
    return {
        # Clearance requires beating *every* same-formula alternative.  The
        # strongest rank and strongest cosine can belong to different entries.
        "rank": best_rank - max(comp_ranks),
        "final": best_final - max(comp_finals),
        "cosine": best_cos - max(comp_cosines),
    }


def _gap_bin(value: float) -> str:
    if value < 0.01:
        return "<0.01"
    if value < 0.02:
        return "0.01-0.02"
    if value < 0.05:
        return "0.02-0.05"
    if value < 0.10:
        return "0.05-0.10"
    return ">=0.10"


def _peak_bin(count: int) -> str:
    if count < 3:
        return "0-2"
    if count < 5:
        return "3-4"
    if count < 10:
        return "5-9"
    if count < 20:
        return "10-19"
    return ">=20"


def _high_blockers(match: dict[str, Any], ann: dict[str, Any]) -> list[str]:
    blockers: list[str] = []
    cos = _float(match.get("cosine")) or 0.0
    final = _float(match.get("final_score")) or 0.0
    ppm = _float(match.get("mass_ppm"))
    reverse = _float(match.get("reverse_match"))
    source = str(match.get("source_tier") or "unknown").lower()

    if not match.get("has_structure") or not str(match.get("smiles") or "").strip():
        blockers.append("no_validated_structure")
    if cos < 0.80:
        blockers.append("cosine_below_medium")
    elif cos < 0.90:
        blockers.append("cosine_below_base_high")
    if final < 0.50:
        blockers.append("final_below_high")
    if ppm is not None and abs(ppm) > 100:
        blockers.append("ppm_over_100")
    elif ppm is not None and abs(ppm) > 20:
        blockers.append("ppm_over_20")
    if reverse is not None and reverse < 0.60:
        blockers.append("reverse_below_0.60")
    if ann.get("isomer_ambiguous"):
        blockers.append("same_formula_competitor")
    if (_float(ann.get("class_agreement")) or 0.0) < 0.50:
        blockers.append("class_consensus_below_0.50")
    if source in {"predicted", "quarantined", "literature_presence_only", "user_unverified"}:
        blockers.append(f"source_cap:{source}")
    return blockers


@dataclass
class TaskAudit:
    task_id: str
    filename: str
    metrics: dict[str, Any]


def audit_task(task_id: str, filename: str, result_path: Path) -> TaskAudit:
    with result_path.open("r", encoding="utf-8") as handle:
        rows = json.load(handle)

    counts: Counter[str] = Counter()
    blockers: Counter[str] = Counter()
    peak_bins_matched: Counter[str] = Counter()
    peak_bins_unmatched: Counter[str] = Counter()
    ppm_high: Counter[str] = Counter()
    ambiguity_rank_gaps: Counter[str] = Counter()
    ambiguity_cos_gaps: Counter[str] = Counter()
    identities: dict[str, list[dict[str, Any]]] = defaultdict(list)
    high_identities: set[str] = set()
    class_high_identities: set[str] = set()
    medium_identities: set[str] = set()
    clear_winner_downgrades: set[str] = set()

    for row in rows:
        counts["spectra"] += 1
        match = row.get("match") if isinstance(row, dict) else None
        peak_count = int(row.get("peaks_count") or len(row.get("query_peaks") or []))
        if not isinstance(match, dict):
            counts["unmatched_spectra"] += 1
            peak_bins_unmatched[_peak_bin(peak_count)] += 1
            continue

        counts["matched_spectra"] += 1
        peak_bins_matched[_peak_bin(peak_count)] += 1
        ann = annotate(match, match.get("top_candidates"))
        identity = _identity(match)
        identities[identity].append({"row": row, "match": match, "ann": ann})

        confidence = str(ann.get("compound_confidence") or "Unknown")
        class_confidence = str(ann.get("class_confidence") or "Unknown")
        counts[f"confidence:{confidence}"] += 1
        counts[f"class_confidence:{class_confidence}"] += 1
        if class_confidence == "High":
            class_high_identities.add(identity)
        if confidence == "High":
            high_identities.add(identity)
            ppm = _float(match.get("mass_ppm"))
            if ppm is None:
                ppm_high["missing"] += 1
            elif abs(ppm) <= 10:
                ppm_high["<=10"] += 1
            elif abs(ppm) <= 20:
                ppm_high["10-20"] += 1
            elif abs(ppm) <= 50:
                ppm_high["20-50"] += 1
            elif abs(ppm) <= 100:
                ppm_high["50-100"] += 1
            else:
                ppm_high[">100"] += 1
        elif confidence == "Medium":
            medium_identities.add(identity)

        if ann.get("isomer_ambiguous"):
            counts["ambiguous_spectra"] += 1
            gaps = _best_competitor_gaps(match)
            if gaps:
                ambiguity_rank_gaps[_gap_bin(gaps["rank"])] += 1
                ambiguity_cos_gaps[_gap_bin(gaps["cosine"])] += 1
                # A same-formula entry far below the winner should remain visible
                # in Top-K but should not automatically mean MS/MS is ambiguous.
                if (
                    gaps["rank"] > ISOMER_AMBIGUITY_MAX_RANK_GAP
                    and gaps["cosine"] > ISOMER_AMBIGUITY_MAX_COSINE_GAP
                ):
                    counts["clear_winner_but_marked_ambiguous_spectra"] += 1
                    clear_winner_downgrades.add(identity)

        if confidence != "High":
            blockers.update(_high_blockers(match, ann))

        if match.get("rt_confirmed"):
            counts["rt_confirmed_spectra"] += 1
        isotope = match.get("formula_isotope_evidence")
        if isinstance(isotope, dict):
            counts["formula_isotope_available_spectra"] += 1
        pathway = ((match.get("expert_details") or {}).get("pathway_align") or {})
        if (_float(pathway.get("score")) or 0.0) > 0:
            counts["pathway_evidence_spectra"] += 1

        tc = match.get("top_candidates") or []
        if tc and _identity(tc[0]) != identity:
            counts["top_candidates_first_not_best"] += 1

    counts["unique_compounds"] = len(identities)
    counts["unique_high_compounds"] = len(high_identities)
    counts["unique_class_high_compounds"] = len(class_high_identities)
    counts["unique_medium_compounds"] = len(medium_identities)
    counts["unique_clear_winner_but_marked_ambiguous"] = len(clear_winner_downgrades)

    metrics = {
        "counts": dict(counts),
        "non_high_blockers_spectrum_counts": dict(blockers.most_common()),
        "matched_peak_count_bins": dict(peak_bins_matched),
        "unmatched_peak_count_bins": dict(peak_bins_unmatched),
        "high_mass_ppm_bins": dict(ppm_high),
        "ambiguous_rank_gap_bins": dict(ambiguity_rank_gaps),
        "ambiguous_cosine_gap_bins": dict(ambiguity_cos_gaps),
    }
    return TaskAudit(task_id=task_id, filename=filename, metrics=metrics)


def _iter_selected(index: dict[str, Any], selected: Iterable[str]) -> Iterable[tuple[str, dict[str, Any]]]:
    selected_set = set(selected)
    for task_id, meta in index.items():
        if selected_set and task_id not in selected_set:
            continue
        if meta.get("status") == "completed" and meta.get("has_results_file"):
            yield task_id, meta


def _render_markdown(audits: list[TaskAudit]) -> str:
    lines = [
        "# V6 High-confidence funnel audit",
        "",
        "Counts distinguish matched spectra from unique compound identities.",
        "",
        "| Task | File | MS2 | Matched | Unique | High spectra | Unique High | Class-High spectra | Unique Class-High | Ambiguous | Clear winner falsely marked ambiguous |",
        "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    aggregate = Counter()
    aggregate_blockers = Counter()
    aggregate_ppm = Counter()
    for audit in audits:
        counts = audit.metrics["counts"]
        aggregate.update(counts)
        aggregate_blockers.update(audit.metrics["non_high_blockers_spectrum_counts"])
        aggregate_ppm.update(audit.metrics["high_mass_ppm_bins"])
        lines.append(
            "| {task} | {file} | {spectra} | {matched} | {unique} | {high} | {uhigh} | {chigh} | {uchigh} | {amb} | {clear} |".format(
                task=audit.task_id,
                file=audit.filename.replace("|", "\\|"),
                spectra=counts.get("spectra", 0),
                matched=counts.get("matched_spectra", 0),
                unique=counts.get("unique_compounds", 0),
                high=counts.get("confidence:High", 0),
                uhigh=counts.get("unique_high_compounds", 0),
                chigh=counts.get("class_confidence:High", 0),
                uchigh=counts.get("unique_class_high_compounds", 0),
                amb=counts.get("ambiguous_spectra", 0),
                clear=counts.get("clear_winner_but_marked_ambiguous_spectra", 0),
            )
        )

    lines.extend([
        "",
        "## Aggregate spectrum-level blockers",
        "",
        "| Blocker | Spectra |",
        "|---|---:|",
    ])
    for reason, count in aggregate_blockers.most_common():
        lines.append(f"| {reason} | {count} |")

    lines.extend([
        "",
        "## Mass error among spectra currently labelled High",
        "",
        "| Absolute ppm | Spectra |",
        "|---|---:|",
    ])
    for label in ("<=10", "10-20", "20-50", "50-100", ">100", "missing"):
        lines.append(f"| {label} | {aggregate_ppm.get(label, 0)} |")

    lines.extend([
        "",
        "## Interpretation guard",
        "",
        "A clear-winner flag means every same-formula Top-K competitor is more than "
        f"{ISOMER_AMBIGUITY_MAX_RANK_GAP:.2f} lower in rank score and "
        f"{ISOMER_AMBIGUITY_MAX_COSINE_GAP:.2f} lower in cosine. It is a calibration "
        "candidate, not an automatic identification claim.",
        "",
    ])
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("task_ids", nargs="*", help="Task IDs; default audits all completed tasks")
    parser.add_argument("--task-dir", type=Path, default=DEFAULT_TASK_DIR)
    parser.add_argument("--index", type=Path, default=DEFAULT_INDEX)
    parser.add_argument("--json-output", type=Path)
    parser.add_argument("--markdown-output", type=Path)
    args = parser.parse_args()

    index = json.loads(args.index.read_text(encoding="utf-8"))
    audits: list[TaskAudit] = []
    for task_id, meta in _iter_selected(index, args.task_ids):
        result_path = args.task_dir / f"{task_id}_results.json"
        if not result_path.exists():
            continue
        audits.append(audit_task(task_id, str(meta.get("filename") or ""), result_path))

    payload = {
        "schema_version": "1.0",
        "task_count": len(audits),
        "tasks": [
            {"task_id": audit.task_id, "filename": audit.filename, **audit.metrics}
            for audit in audits
        ],
    }
    if args.json_output:
        args.json_output.parent.mkdir(parents=True, exist_ok=True)
        args.json_output.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
    markdown = _render_markdown(audits)
    if args.markdown_output:
        args.markdown_output.parent.mkdir(parents=True, exist_ok=True)
        args.markdown_output.write_text(markdown, encoding="utf-8")
    print(markdown)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
