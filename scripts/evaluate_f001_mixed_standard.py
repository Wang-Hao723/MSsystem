#!/usr/bin/env python3
"""Reproducibly evaluate one F001 mixed-standard task against curated truth.

This evaluator is deliberately post hoc: the 20-standard manifest is never
passed to candidate generation or ranking.  It evaluates saved Top-1/Top-5/
Top-K identities by full InChIKey and by 2D-connectivity block, and separately
reports whether a polarity-compatible precursor event was actually observed
in this single injection.  Historical pooled, multi-polarity results are
included only as context and never used as this task's denominator.
"""

from __future__ import annotations

import argparse
import json
import math
import statistics
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_TRUTH = (
    PROJECT_ROOT
    / "outputs"
    / "019fc5bf-15dc-76a1-9d64-57597ecfcda1"
    / "inventory_reconcile_corrected"
    / "八个方_库存采购核对.json"
)
DEFAULT_HISTORICAL = (
    PROJECT_ROOT
    / "eval"
    / "v7val"
    / "results"
    / "external_compare_fractions.json"
)

PROTON_MASS = 1.007276466621
CHLORIDE_ION_MASS = 34.969401259
FORMIC_ACID_MASS = 46.005479304
ACETIC_ACID_MASS = 60.021129368
SODIUM_ION_MASS = 22.989218
POTASSIUM_ION_MASS = 38.963158
AMMONIUM_ION_MASS = 18.033823
AMMONIA_MASS = 17.026549101


def _read_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def _finite(value: Any) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def _inchikey(value: Any) -> str:
    text = str(value or "").strip().upper().replace("INCHIKEY=", "")
    return text if len(text) >= 14 else ""


def _connectivity(value: Any) -> str:
    return _inchikey(value)[:14]


def _candidate_identity(candidate: dict[str, Any]) -> tuple[str, str]:
    key = _inchikey(candidate.get("inchikey"))
    fallback = "|".join(
        [
            str(candidate.get("compound_name") or "").strip().casefold(),
            str(candidate.get("formula") or "").strip(),
        ]
    )
    return key, key or fallback


def _saved_candidates(row: dict[str, Any]) -> list[dict[str, Any]]:
    match = row.get("match")
    if not isinstance(match, dict):
        return []
    candidates: list[dict[str, Any]] = [match]
    candidates.extend(
        candidate
        for candidate in (match.get("top_candidates") or [])
        if isinstance(candidate, dict)
    )
    output: list[dict[str, Any]] = []
    seen: set[str] = set()
    for candidate in candidates:
        _, identity = _candidate_identity(candidate)
        if not identity or identity in seen:
            continue
        seen.add(identity)
        output.append(candidate)
    return output


def _truth_rows(path: Path, formula_id: str) -> list[dict[str, Any]]:
    payload = _read_json(path)
    rows = []
    for row in payload.get("procurement") or []:
        formula_ids = {
            token.strip()
            for token in str(row.get("formula_ids") or "").split("；")
            if token.strip()
        }
        if formula_id not in formula_ids:
            continue
        key = _inchikey(row.get("inchikey"))
        exact_mass = _finite(row.get("monoisotopic_mass"))
        if not key or exact_mass is None:
            raise ValueError(
                f"Curated truth is incomplete for {row.get('name_en') or row.get('name_zh')}"
            )
        rows.append(
            {
                "name_zh": str(row.get("name_zh") or ""),
                "name_en": str(row.get("name_en") or ""),
                "formula": str(row.get("molecular_formula") or ""),
                "monoisotopic_mass": exact_mass,
                "inchikey": key,
                "connectivity": key[:14],
                "acceptable_connectivities": [key[:14]],
                "evaluation_identity_scope": "declared_standard_connectivity",
                "standard_codes": str(row.get("inventory_ids") or ""),
                "current_system_issue": str(row.get("current_system_issue") or ""),
                "software_next_action": str(row.get("software_next_action") or ""),
            }
        )
        if str(row.get("name_en") or "").strip().casefold() == (
            "monoammonium glycyrrhizinate"
        ):
            # ESI observes the glycyrrhizinate/acid molecular moiety after
            # ammonium-counterion dissociation.  Parent-moiety recovery is a
            # valid salt-normalized result, but is kept separate from exact
            # salt-form identification.
            rows[-1]["acceptable_connectivities"] = [
                key[:14],
                "LPLVUJXQOOQHMX",
                "QCQGNTXXPCJRFZ",
            ]
            rows[-1]["evaluation_identity_scope"] = (
                "salt_normalized_glycyrrhizinate_parent_moiety"
            )
    rows.sort(key=lambda row: (row["name_en"].casefold(), row["inchikey"]))
    if formula_id == "F001" and len(rows) != 20:
        raise ValueError(
            f"F001 curated truth must contain exactly 20 identities; found {len(rows)}"
        )
    return rows


def _historical_pool(path: Path, formula_id: str) -> dict[str, dict[str, Any]]:
    if not path.exists():
        return {}
    payload = _read_json(path)
    fraction = (payload.get("per_fraction") or {}).get(formula_id) or {}
    return {
        str(row.get("name") or ""): {
            "standard_code": row.get("code"),
            "in_historical_shared_pool": bool(row.get("in_pool")),
            "historical_pooled_msms_scans": int(row.get("n_scans") or 0),
            "historical_v7_pooled_top1": bool(row.get("v7")),
        }
        for row in (fraction.get("per_standard") or [])
    }


def _ion_hypotheses(truth: dict[str, Any], polarity: str) -> list[dict[str, Any]]:
    mass = float(truth["monoisotopic_mass"])
    polarity = polarity.upper()
    if polarity == "NEG":
        hypotheses = [
            ("[M-H]-", mass - PROTON_MASS, 1),
            ("[M+Cl]-", mass + CHLORIDE_ION_MASS, 1),
            ("[M+HCOO]-", mass + FORMIC_ACID_MASS - PROTON_MASS, 1),
            ("[M+CH3COO]-", mass + ACETIC_ACID_MASS - PROTON_MASS, 1),
            ("[M-2H]2-", (mass - 2 * PROTON_MASS) / 2, 2),
        ]
        if truth["name_en"].casefold() == "monoammonium glycyrrhizinate":
            parent_mass = mass - AMMONIA_MASS
            hypotheses.append(
                ("[glycyrrhizic acid-H]- (salt-normalized)", parent_mass - PROTON_MASS, 1)
            )
    else:
        hypotheses = [
            ("[M+H]+", mass + PROTON_MASS, 1),
            ("[M+Na]+", mass + SODIUM_ION_MASS, 1),
            ("[M+NH4]+", mass + AMMONIUM_ION_MASS, 1),
            ("[M+K]+", mass + POTASSIUM_ION_MASS, 1),
            ("[M+2H]2+", (mass + 2 * PROTON_MASS) / 2, 2),
        ]
    return [
        {"label": label, "mz": mz, "charge_state": charge}
        for label, mz, charge in hypotheses
        if mz > 0
    ]


def _row_mass(row: dict[str, Any]) -> float | None:
    return _finite(row.get("search_precursor_mz")) or _finite(
        row.get("precursor_mz")
    )


def _near_events(
    rows: Iterable[dict[str, Any]],
    truth: dict[str, Any],
    polarity: str,
    tolerance_ppm: float,
) -> list[dict[str, Any]]:
    events = []
    for row in rows:
        observed = _row_mass(row)
        if observed is None:
            continue
        row_polarity = str(row.get("ion_mode") or row.get("polarity") or polarity).upper()
        if row_polarity and row_polarity != polarity.upper():
            continue
        row_charge = _finite(row.get("precursor_charge"))
        row_charge_state = int(round(abs(row_charge))) if row_charge not in (None, 0) else None
        for hypothesis in _ion_hypotheses(truth, polarity):
            if row_charge_state is not None and row_charge_state != hypothesis["charge_state"]:
                continue
            ppm = (observed - hypothesis["mz"]) / hypothesis["mz"] * 1e6
            if abs(ppm) > tolerance_ppm:
                continue
            events.append(
                {
                    "spectrum_id": row.get("spectrum_id"),
                    "retention_time_min": _finite(row.get("retention_time")),
                    "observed_search_mz": observed,
                    "selected_ion_mz": _finite(row.get("precursor_mz")),
                    "search_mass_source": row.get("search_precursor_mz_source"),
                    "lcms_feature_used_for_search": bool(
                        row.get("lcms_feature_used_for_search")
                    ),
                    "hypothesis": hypothesis["label"],
                    "theoretical_mz": hypothesis["mz"],
                    "mass_error_ppm": round(ppm, 4),
                    "has_library_match": isinstance(row.get("match"), dict),
                }
            )
            break
    events.sort(key=lambda event: (abs(event["mass_error_ppm"]), event["retention_time_min"] or 0))
    return events


def _identity_rank(
    candidates: list[dict[str, Any]], truth: dict[str, Any], exact: bool
) -> int | None:
    wanted = (
        {truth["inchikey"]}
        if exact
        else set(truth.get("acceptable_connectivities") or [truth["connectivity"]])
    )
    for index, candidate in enumerate(candidates, 1):
        observed = _inchikey(candidate.get("inchikey"))
        observed = observed if exact else observed[:14]
        if observed and observed in wanted:
            return index
    return None


def _truth_containing_isomer_group(
    row: dict[str, Any], truth: dict[str, Any]
) -> bool:
    match = row.get("match")
    if not isinstance(match, dict) or not bool(match.get("isomer_ambiguous")):
        return False
    if str(match.get("isomer_group_confidence") or "").upper() not in {
        "HIGH",
        "MEDIUM",
    }:
        return False
    wanted = set(
        truth.get("acceptable_connectivities") or [truth["connectivity"]]
    )
    candidates = []
    candidates.extend(match.get("isomer_set") or [])
    candidates.extend(match.get("same_formula_alternatives") or [])
    return any(
        _connectivity(candidate.get("inchikey")) in wanted
        for candidate in candidates
        if isinstance(candidate, dict)
    )


def _pct(count: int, total: int) -> float | None:
    return round(100.0 * count / total, 1) if total else None


def evaluate(
    task_id: str,
    formula_id: str,
    truth_path: Path,
    historical_path: Path,
    mass_tolerance_ppm: float,
) -> dict[str, Any]:
    results_path = PROJECT_ROOT / "output" / "tasks" / f"{task_id}_results.json"
    if not results_path.exists():
        raise FileNotFoundError(results_path)
    results = _read_json(results_path)
    if not isinstance(results, list):
        raise ValueError("Task results must be a JSON list")
    task_index = _read_json(PROJECT_ROOT / "output" / "tasks_index.json")
    task = dict(task_index.get(task_id) or {})
    truth_rows = _truth_rows(truth_path, formula_id)
    historical = _historical_pool(historical_path, formula_id)
    polarity = str(task.get("ion_mode") or "NEG").upper()

    per_standard = []
    for truth in truth_rows:
        exact_ranks: list[int] = []
        connectivity_ranks: list[int] = []
        exact_top1_rows = []
        connectivity_top1_rows = []
        truth_containing_group_rows = []
        for row in results:
            candidates = _saved_candidates(row)
            if not candidates:
                continue
            exact_rank = _identity_rank(candidates, truth, exact=True)
            connectivity_rank = _identity_rank(candidates, truth, exact=False)
            if exact_rank is not None:
                exact_ranks.append(exact_rank)
            if connectivity_rank is not None:
                connectivity_ranks.append(connectivity_rank)
            if exact_rank == 1:
                exact_top1_rows.append(row)
            if connectivity_rank == 1:
                connectivity_top1_rows.append(row)
            if _truth_containing_isomer_group(row, truth):
                truth_containing_group_rows.append(row)

        near = _near_events(
            results, truth, polarity=polarity, tolerance_ppm=mass_tolerance_ppm
        )
        historical_row = historical.get(truth["name_zh"], {})
        best_row = (exact_top1_rows or connectivity_top1_rows or [None])[0]
        best_match = (best_row or {}).get("match") or {}
        per_standard.append(
            {
                **truth,
                **historical_row,
                "polarity": polarity,
                "compatible_precursor_event_count": len(near),
                "best_compatible_precursor_event": near[0] if near else None,
                "saved_topk_exact": bool(exact_ranks),
                "saved_topk_connectivity": bool(connectivity_ranks),
                "top5_exact_any_scan": any(rank <= 5 for rank in exact_ranks),
                "top5_connectivity_any_scan": any(
                    rank <= 5 for rank in connectivity_ranks
                ),
                "top1_exact_any_scan": bool(exact_top1_rows),
                "top1_connectivity_any_scan": bool(connectivity_top1_rows),
                "truth_containing_isomer_group_any_scan": bool(
                    truth_containing_group_rows
                ),
                "reported_scope_supported_any_scan": bool(
                    connectivity_top1_rows or truth_containing_group_rows
                ),
                "best_saved_rank_exact": min(exact_ranks) if exact_ranks else None,
                "best_saved_rank_connectivity": (
                    min(connectivity_ranks) if connectivity_ranks else None
                ),
                "top1_exact_scan_count": len(exact_top1_rows),
                "top1_connectivity_scan_count": len(connectivity_top1_rows),
                "truth_containing_isomer_group_scan_count": len(
                    truth_containing_group_rows
                ),
                "best_correct_top1": (
                    {
                        "spectrum_id": best_row.get("spectrum_id"),
                        "retention_time_min": _finite(best_row.get("retention_time")),
                        "selected_ion_mz": _finite(best_row.get("precursor_mz")),
                        "search_precursor_mz": _row_mass(best_row),
                        "search_mass_source": best_row.get(
                            "search_precursor_mz_source"
                        ),
                        "candidate_name": best_match.get("compound_name"),
                        "candidate_inchikey": best_match.get("inchikey"),
                        "confidence": best_match.get("confidence"),
                        "cosine": _finite(best_match.get("cosine")),
                        "final_score": _finite(best_match.get("final_score")),
                        "mass_ppm": _finite(best_match.get("mass_ppm")),
                    }
                    if best_row
                    else None
                ),
            }
        )

    observable = [
        row for row in per_standard if row["compatible_precursor_event_count"] > 0
    ]
    feature_rows = [
        row for row in results if bool(row.get("lcms_feature_used_for_search"))
    ]
    feature_deltas = [
        abs(value)
        for value in (
            _finite(row.get("search_precursor_delta_ppm")) for row in feature_rows
        )
        if value is not None
    ]

    def count(key: str, rows: list[dict[str, Any]] = per_standard) -> int:
        return sum(bool(row.get(key)) for row in rows)

    summary = {
        "declared_truth_identities": len(per_standard),
        "historical_shared_pool_covered": sum(
            bool(row.get("in_historical_shared_pool")) for row in per_standard
        ),
        "single_injection_compatible_precursor_identities": len(observable),
        "saved_topk_exact_recall": count("saved_topk_exact"),
        "saved_topk_connectivity_recall": count("saved_topk_connectivity"),
        "top5_exact_any_scan": count("top5_exact_any_scan"),
        "top5_connectivity_any_scan": count("top5_connectivity_any_scan"),
        "top1_exact_any_scan": count("top1_exact_any_scan"),
        "top1_connectivity_any_scan": count("top1_connectivity_any_scan"),
        "truth_containing_isomer_group_any_scan": count(
            "truth_containing_isomer_group_any_scan"
        ),
        "reported_scope_supported_any_scan": count(
            "reported_scope_supported_any_scan"
        ),
        "observable_top1_exact": count("top1_exact_any_scan", observable),
        "observable_top1_connectivity": count(
            "top1_connectivity_any_scan", observable
        ),
        "observable_reported_scope_supported": count(
            "reported_scope_supported_any_scan", observable
        ),
    }
    summary.update(
        {
            "saved_topk_connectivity_recall_pct_all_declared": _pct(
                summary["saved_topk_connectivity_recall"], len(per_standard)
            ),
            "top1_connectivity_pct_all_declared": _pct(
                summary["top1_connectivity_any_scan"], len(per_standard)
            ),
            "top1_connectivity_pct_observable": _pct(
                summary["observable_top1_connectivity"], len(observable)
            ),
            "reported_scope_supported_pct_observable": _pct(
                summary["observable_reported_scope_supported"], len(observable)
            ),
        }
    )
    params = task.get("params") or {}
    return {
        "schema_version": 1,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "task_id": task_id,
        "filename": task.get("filename"),
        "formula_id": formula_id,
        "evaluation_design": (
            "post_hoc_curated_mixed_standard_truth; truth excluded from inference"
        ),
        "scope": (
            f"single {polarity} injection; not pooled across polarity or replicate files"
        ),
        "truth_source": str(truth_path),
        "historical_context_source": str(historical_path),
        "mass_tolerance_ppm_for_observable_event_audit": mass_tolerance_ppm,
        "expected_compound_selection": {
            "selected_datasets": list(params.get("expected_compound_datasets") or []),
            "selection_source": params.get(
                "expected_compound_dataset_selection_source"
            ),
            "selection_status": params.get(
                "expected_compound_dataset_selection_status"
            ),
            "used_as_truth_for_ranking": False,
        },
        "task_counts": {
            "ms2_spectra": len(results),
            "matched_spectra": sum(
                isinstance(row.get("match"), dict) for row in results
            ),
            "ms1_feature_centroid_search_spectra": len(feature_rows),
            "selected_ion_fallback_search_spectra": len(results) - len(feature_rows),
            "median_abs_selected_to_ms1_delta_ppm": (
                round(statistics.median(feature_deltas), 4)
                if feature_deltas
                else None
            ),
        },
        "summary": summary,
        "denominator_warning": (
            "The 20 declared standards are mixture truth. A single NEG injection is not "
            "expected to yield a polarity-compatible selected MS2 event for every standard. "
            "Metrics using the observable denominator are reported separately."
        ),
        "historical_warning": (
            "Historical V7 F001 18/20 is pooled across a shared eight-fraction retrieval "
            "design and is not a direct comparator for this single NEG task."
        ),
        "per_standard": per_standard,
    }


def _markdown(payload: dict[str, Any]) -> str:
    summary = payload["summary"]
    rows = payload["per_standard"]
    lines = [
        f"# {payload['formula_id']} mixed-standard audit — task `{payload['task_id']}`",
        "",
        f"- File: `{payload.get('filename') or ''}`",
        f"- Scope: {payload['scope']}",
        f"- Truth: {summary['declared_truth_identities']} curated identities; never supplied to ranking",
        f"- Expected-Compounds selection: `{payload['expected_compound_selection']['selection_status']}`; selected datasets = {payload['expected_compound_selection']['selected_datasets']}",
        f"- MS1 measured-centroid search: {payload['task_counts']['ms1_feature_centroid_search_spectra']}/{payload['task_counts']['ms2_spectra']} spectra",
        "",
        "## Results",
        "",
        f"Saved Top-K connectivity recall: **{summary['saved_topk_connectivity_recall']}/{summary['declared_truth_identities']}** ({summary['saved_topk_connectivity_recall_pct_all_declared']}%).",
        f"Any-scan Top-1 connectivity: **{summary['top1_connectivity_any_scan']}/{summary['declared_truth_identities']}** ({summary['top1_connectivity_pct_all_declared']}%).",
        f"Among identities with a compatible precursor event in this NEG injection: **{summary['observable_top1_connectivity']}/{summary['single_injection_compatible_precursor_identities']}** Top-1 ({summary['top1_connectivity_pct_observable']}%).",
        f"At the actually reportable scope (Top-1 parent/connectivity or a truth-containing isomer group): **{summary['observable_reported_scope_supported']}/{summary['single_injection_compatible_precursor_identities']}** ({summary['reported_scope_supported_pct_observable']}%).",
        "",
        f"> {payload['denominator_warning']}",
        "",
        f"> {payload['historical_warning']}",
        "",
        "## Per standard",
        "",
        "| Standard | Formula | NEG-compatible event | Best saved rank (2D/parent) | Top-1 (2D/parent) | Truth-containing isomer group | Reported scope supported | Best candidate / confidence |",
        "|---|---:|---:|---:|---:|---:|---:|---|",
    ]
    for row in rows:
        best = row.get("best_correct_top1") or {}
        candidate = str(best.get("candidate_name") or "—").replace("|", "\\|")
        confidence = best.get("confidence") or "—"
        lines.append(
            "| {name} | {formula} | {event} | {rank} | {top1} | {group} | {scope} | {candidate} / {confidence} |".format(
                name=(row["name_en"] or row["name_zh"]).replace("|", "\\|"),
                formula=row["formula"],
                event=row["compatible_precursor_event_count"],
                rank=row["best_saved_rank_connectivity"] or "—",
                top1="yes" if row["top1_connectivity_any_scan"] else "no",
                group=(
                    "yes" if row["truth_containing_isomer_group_any_scan"] else "no"
                ),
                scope="yes" if row["reported_scope_supported_any_scan"] else "no",
                candidate=candidate,
                confidence=confidence,
            )
        )
    lines.extend(
        [
            "",
            "Exact-stereo and connectivity metrics are deliberately separate because many public library records omit or disagree on stereochemical layers. Connectivity success must not be relabelled as authentic-standard Level 1 confirmation.",
            "",
        ]
    )
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--task-id", required=True)
    parser.add_argument("--formula-id", default="F001")
    parser.add_argument("--truth-json", type=Path, default=DEFAULT_TRUTH)
    parser.add_argument("--historical-json", type=Path, default=DEFAULT_HISTORICAL)
    parser.add_argument("--mass-tolerance-ppm", type=float, default=10.0)
    parser.add_argument("--output-prefix", type=Path)
    args = parser.parse_args()

    payload = evaluate(
        task_id=args.task_id,
        formula_id=args.formula_id,
        truth_path=args.truth_json.resolve(),
        historical_path=args.historical_json.resolve(),
        mass_tolerance_ppm=args.mass_tolerance_ppm,
    )
    prefix = args.output_prefix or (
        PROJECT_ROOT
        / "docs"
        / "audits"
        / f"MSSYSTEM_V12_{args.formula_id}_MIXED_STANDARD_{args.task_id}"
    )
    prefix.parent.mkdir(parents=True, exist_ok=True)
    json_path = prefix.with_suffix(".json")
    markdown_path = prefix.with_suffix(".md")
    json_path.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    markdown_path.write_text(_markdown(payload), encoding="utf-8")
    print(
        json.dumps(
            {
                "json": str(json_path),
                "markdown": str(markdown_path),
                "summary": payload["summary"],
            },
            ensure_ascii=False,
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
