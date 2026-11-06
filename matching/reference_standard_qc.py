"""Post-hoc reference-standard recovery audit.

Reference-standard truth is an evaluation label, never an identification
feature.  This module therefore accepts only already-committed result rows and
cannot call retrieval, scoring, ranking, or confidence calibration code.  A
set must be selected explicitly; filenames are deliberately ignored.
"""

from __future__ import annotations

import json
import math
import os
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Mapping


PROJECT_ROOT = Path(__file__).resolve().parents[1]
REGISTRY_PATH = PROJECT_ROOT / "config" / "reference_standard_sets.json"
ARTIFACT_SUFFIX = "_reference_standard_qc.json"
REFERENCE_STANDARD_QC_SCHEMA_VERSION = 2

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


def _truth_path(registry: Mapping[str, Any]) -> Path:
    configured = Path(str(registry.get("truth_source") or ""))
    return configured if configured.is_absolute() else PROJECT_ROOT / configured


def load_reference_standard_registry() -> dict[str, Any]:
    registry = _read_json(REGISTRY_PATH)
    if registry.get("selection_policy") != "explicit_only_no_filename_inference":
        raise ValueError("Reference-standard registry must require explicit selection")
    if registry.get("used_for_candidate_generation") is not False:
        raise ValueError("Reference-standard truth cannot be used for candidate generation")
    if registry.get("used_for_ranking") is not False:
        raise ValueError("Reference-standard truth cannot be used for ranking")
    sets = registry.get("sets") or []
    ids = [str(item.get("id") or "").strip().upper() for item in sets]
    if not ids or len(ids) != len(set(ids)) or any(not item for item in ids):
        raise ValueError("Reference-standard set identifiers must be non-empty and unique")
    truth_path = _truth_path(registry)
    if not truth_path.is_file():
        raise FileNotFoundError(truth_path)
    registry = dict(registry)
    registry["truth_source_resolved"] = str(truth_path.resolve())
    return registry


def list_reference_standard_sets() -> list[dict[str, Any]]:
    registry = load_reference_standard_registry()
    truth = _read_json(_truth_path(registry))
    actual_counts: dict[str, int] = {}
    for row in truth.get("procurement") or []:
        for set_id in str(row.get("formula_ids") or "").split("；"):
            set_id = set_id.strip().upper()
            if set_id:
                actual_counts[set_id] = actual_counts.get(set_id, 0) + 1
    output = []
    for item in registry.get("sets") or []:
        clean = dict(item)
        set_id = str(clean.get("id") or "").strip().upper()
        declared = int(clean.get("declared_standard_count") or 0)
        actual = int(actual_counts.get(set_id) or 0)
        if declared != actual:
            raise ValueError(
                f"{set_id} registry/truth count mismatch: {declared} != {actual}"
            )
        clean.update(
            {
                "id": set_id,
                "truth_identity_count": actual,
                "selection_required": True,
                "used_for_scoring": False,
                "used_for_ranking": False,
            }
        )
        output.append(clean)
    return output


def validate_reference_standard_set_id(value: Any) -> str:
    set_id = str(value or "").strip().upper()
    if not set_id:
        return ""
    allowed = {item["id"] for item in list_reference_standard_sets()}
    if set_id not in allowed:
        raise ValueError(
            f"Unknown reference-standard QC set {set_id!r}; "
            f"available sets: {', '.join(sorted(allowed))}"
        )
    return set_id


def _truth_rows(set_id: str) -> list[dict[str, Any]]:
    registry = load_reference_standard_registry()
    payload = _read_json(_truth_path(registry))
    rows: list[dict[str, Any]] = []
    for row in payload.get("procurement") or []:
        formula_ids = {
            token.strip().upper()
            for token in str(row.get("formula_ids") or "").split("；")
            if token.strip()
        }
        if set_id not in formula_ids:
            continue
        key = _inchikey(row.get("inchikey"))
        exact_mass = _finite(row.get("monoisotopic_mass"))
        if not key or exact_mass is None:
            raise ValueError(
                "Curated reference-standard truth is incomplete for "
                f"{row.get('name_en') or row.get('name_zh')}"
            )
        truth = {
            "name_zh": str(row.get("name_zh") or ""),
            "name_en": str(row.get("name_en") or ""),
            "formula": str(row.get("molecular_formula") or ""),
            "monoisotopic_mass": exact_mass,
            "inchikey": key,
            "connectivity": key[:14],
            "acceptable_connectivities": [key[:14]],
            "evaluation_identity_scope": "declared_standard_connectivity",
            "standard_codes": str(row.get("inventory_ids") or ""),
        }
        if str(row.get("name_en") or "").strip().casefold() == (
            "monoammonium glycyrrhizinate"
        ):
            # ESI commonly observes the acid/glycyrrhizinate parent after the
            # ammonium counter-ion dissociates.  Keep the salt-normalized scope
            # explicit instead of calling it an exact salt-form identification.
            truth["acceptable_connectivities"] = [
                key[:14],
                "LPLVUJXQOOQHMX",
                "QCQGNTXXPCJRFZ",
            ]
            truth["evaluation_identity_scope"] = (
                "salt_normalized_glycyrrhizinate_parent_moiety"
            )
        rows.append(truth)
    rows.sort(key=lambda item: (item["name_en"].casefold(), item["inchikey"]))
    expected = {
        item["id"]: int(item["declared_standard_count"])
        for item in list_reference_standard_sets()
    }[set_id]
    if len(rows) != expected:
        raise ValueError(f"{set_id} truth must contain {expected} identities; found {len(rows)}")
    return rows


def _candidate_identity(candidate: Mapping[str, Any]) -> tuple[str, str]:
    key = _inchikey(candidate.get("inchikey"))
    fallback = "|".join(
        [
            str(candidate.get("compound_name") or "").strip().casefold(),
            str(candidate.get("formula") or "").strip(),
        ]
    )
    return key, key or fallback


def _saved_candidates(row: Mapping[str, Any]) -> list[dict[str, Any]]:
    match = row.get("match")
    if not isinstance(match, Mapping):
        return []
    candidates: list[Mapping[str, Any]] = [match]
    candidates.extend(
        candidate
        for candidate in (match.get("top_candidates") or [])
        if isinstance(candidate, Mapping)
    )
    output: list[dict[str, Any]] = []
    seen: set[str] = set()
    for candidate in candidates:
        _, identity = _candidate_identity(candidate)
        if not identity or identity in seen:
            continue
        seen.add(identity)
        output.append(dict(candidate))
    return output


def _ion_hypotheses(truth: Mapping[str, Any], polarity: str) -> list[dict[str, Any]]:
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
        if str(truth.get("name_en") or "").casefold() == (
            "monoammonium glycyrrhizinate"
        ):
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


def _row_mass(row: Mapping[str, Any]) -> float | None:
    return _finite(row.get("search_precursor_mz")) or _finite(row.get("precursor_mz"))


def _compatible_precursor_event(
    row: Mapping[str, Any],
    truth: Mapping[str, Any],
    polarity: str,
    tolerance_ppm: float,
) -> dict[str, Any] | None:
    """Return the compatible ion event for one saved MS2 row, if any."""

    observed = _row_mass(row)
    if observed is None:
        return None
    row_polarity = str(
        row.get("ion_mode") or row.get("polarity") or polarity
    ).upper()
    if row_polarity and row_polarity != polarity.upper():
        return None
    charge = _finite(row.get("precursor_charge"))
    charge_state = int(round(abs(charge))) if charge not in (None, 0) else None
    for hypothesis in _ion_hypotheses(truth, polarity):
        if charge_state is not None and charge_state != hypothesis["charge_state"]:
            continue
        ppm = (observed - hypothesis["mz"]) / hypothesis["mz"] * 1e6
        if abs(ppm) > tolerance_ppm:
            continue
        return {
            "spectrum_id": row.get("spectrum_id"),
            "retention_time_min": _finite(row.get("retention_time")),
            "observed_search_mz": observed,
            "search_mass_source": row.get("search_precursor_mz_source"),
            "lcms_feature_used_for_search": bool(
                row.get("lcms_feature_used_for_search")
            ),
            "hypothesis": hypothesis["label"],
            "theoretical_mz": hypothesis["mz"],
            "mass_error_ppm": round(ppm, 4),
            "has_library_match": isinstance(row.get("match"), Mapping),
        }
    return None


def _near_events(
    rows: Iterable[Mapping[str, Any]],
    truth: Mapping[str, Any],
    polarity: str,
    tolerance_ppm: float,
) -> list[dict[str, Any]]:
    events = [
        event
        for row in rows
        if (
            event := _compatible_precursor_event(
                row, truth, polarity, tolerance_ppm
            )
        )
        is not None
    ]
    events.sort(key=lambda item: (abs(item["mass_error_ppm"]), item["retention_time_min"] or 0))
    return events


def _identity_rank(
    candidates: list[dict[str, Any]], truth: Mapping[str, Any], exact: bool
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
    row: Mapping[str, Any], truth: Mapping[str, Any], *, high_only: bool = False
) -> bool:
    match = row.get("match")
    if (
        not isinstance(match, Mapping)
        or str(match.get("identification_scope") or "") != "isomer_group"
    ):
        return False
    group_confidence = str(match.get("isomer_group_confidence") or "").upper()
    if group_confidence not in ({"HIGH"} if high_only else {"HIGH", "MEDIUM"}):
        return False
    wanted = set(truth.get("acceptable_connectivities") or [truth["connectivity"]])
    candidates: list[Mapping[str, Any]] = [match]
    candidates.extend(
        candidate
        for candidate in (match.get("isomer_set") or [])
        if isinstance(candidate, Mapping)
    )
    candidates.extend(
        candidate
        for candidate in (match.get("same_formula_alternatives") or [])
        if isinstance(candidate, Mapping)
    )
    return any(_connectivity(candidate.get("inchikey")) in wanted for candidate in candidates)


def _pct(count: int, total: int) -> float | None:
    return round(100.0 * count / total, 1) if total else None


def build_reference_standard_qc(
    task_id: str,
    rows: Iterable[Mapping[str, Any]],
    *,
    set_id: str,
    polarity: str,
    filename: str = "",
    mass_tolerance_ppm: float = 10.0,
    selection_source: str = "explicit_user_selection",
) -> dict[str, Any]:
    """Evaluate saved results against one explicitly selected standard set."""

    set_id = validate_reference_standard_set_id(set_id)
    if not set_id:
        raise ValueError("Reference-standard QC requires an explicit set_id")
    polarity = str(polarity or "").strip().upper()
    if polarity not in {"POS", "NEG"}:
        raise ValueError("Reference-standard QC requires resolved POS or NEG polarity")
    if mass_tolerance_ppm <= 0 or not math.isfinite(float(mass_tolerance_ppm)):
        raise ValueError("mass_tolerance_ppm must be positive and finite")
    result_rows = [dict(row) for row in rows if isinstance(row, Mapping)]
    truth_rows = _truth_rows(set_id)
    per_standard: list[dict[str, Any]] = []

    for truth in truth_rows:
        near: list[dict[str, Any]] = []
        exact_ranks: list[int] = []
        connectivity_ranks: list[int] = []
        exact_top1_rows: list[dict[str, Any]] = []
        connectivity_top1_rows: list[dict[str, Any]] = []
        group_rows: list[dict[str, Any]] = []
        high_group_rows: list[dict[str, Any]] = []
        high_connectivity_rows: list[dict[str, Any]] = []
        for row in result_rows:
            compatible_event = _compatible_precursor_event(
                row, truth, polarity, float(mass_tolerance_ppm)
            )
            if compatible_event is None:
                continue
            near.append(compatible_event)
            candidates = _saved_candidates(row)
            if candidates:
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
                    match = row.get("match") or {}
                    if (
                        str(
                            match.get("two_dimensional_connectivity_confidence")
                            or match.get("confidence")
                            or ""
                        ).upper() == "HIGH"
                        and str(match.get("identification_scope") or "")
                        != "isomer_group"
                    ):
                        high_connectivity_rows.append(row)
            if _truth_containing_isomer_group(row, truth):
                group_rows.append(row)
            if _truth_containing_isomer_group(row, truth, high_only=True):
                high_group_rows.append(row)

        near.sort(
            key=lambda item: (
                abs(item["mass_error_ppm"]),
                item["retention_time_min"] or 0,
            )
        )
        correct_top1_rows = (
            connectivity_top1_rows or exact_top1_rows or group_rows
        )
        best_row = max(
            correct_top1_rows,
            key=lambda row: (
                str(
                    (row.get("match") or {}).get(
                        "two_dimensional_connectivity_confidence"
                    )
                    or (row.get("match") or {}).get("confidence")
                    or ""
                ).upper() == "HIGH",
                str(
                    (row.get("match") or {}).get("isomer_group_confidence")
                    or ""
                ).upper() == "HIGH",
                _finite((row.get("match") or {}).get("rank_score")) or 0.0,
                _finite((row.get("match") or {}).get("cosine")) or 0.0,
            ),
            default=None,
        )
        best_match = (best_row or {}).get("match") or {}
        reportable = bool(connectivity_top1_rows or group_rows)
        high_reportable = bool(high_connectivity_rows or high_group_rows)
        per_standard.append(
            {
                **truth,
                "polarity": polarity,
                "compatible_precursor_event_count": len(near),
                "best_compatible_precursor_event": near[0] if near else None,
                "observable_in_injection": bool(near),
                "saved_topk_exact": bool(exact_ranks),
                "saved_topk_connectivity": bool(connectivity_ranks),
                "top5_exact_any_scan": any(rank <= 5 for rank in exact_ranks),
                "top5_connectivity_any_scan": any(rank <= 5 for rank in connectivity_ranks),
                "top1_exact_any_scan": bool(exact_top1_rows),
                "top1_connectivity_any_scan": bool(connectivity_top1_rows),
                "truth_containing_isomer_group_any_scan": bool(group_rows),
                "truth_containing_high_isomer_group_any_scan": bool(high_group_rows),
                "reported_scope_supported_any_scan": reportable,
                "high_reportable_scope_supported_any_scan": high_reportable,
                "best_saved_rank_exact": min(exact_ranks) if exact_ranks else None,
                "best_saved_rank_connectivity": min(connectivity_ranks) if connectivity_ranks else None,
                "best_correct_top1": (
                    {
                        "spectrum_id": best_row.get("spectrum_id"),
                        "retention_time_min": _finite(best_row.get("retention_time")),
                        "candidate_name": best_match.get("compound_name"),
                        "candidate_inchikey": best_match.get("inchikey"),
                        "exact_structure_confidence": best_match.get(
                            "exact_structure_confidence", best_match.get("confidence")
                        ),
                        "two_dimensional_connectivity_confidence": best_match.get(
                            "two_dimensional_connectivity_confidence",
                            best_match.get("confidence"),
                        ),
                        "isomer_group_confidence": best_match.get("isomer_group_confidence"),
                        "cosine": _finite(best_match.get("cosine")),
                        "final_score": _finite(best_match.get("final_score")),
                        "mass_ppm": _finite(best_match.get("mass_ppm")),
                    }
                    if best_row
                    else None
                ),
            }
        )

    observable = [row for row in per_standard if row["observable_in_injection"]]

    def count(key: str, source: list[dict[str, Any]] = per_standard) -> int:
        return sum(bool(row.get(key)) for row in source)

    summary = {
        "declared_truth_identities": len(per_standard),
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
        "high_reportable_scope_supported_any_scan": count(
            "high_reportable_scope_supported_any_scan"
        ),
        "observable_saved_topk_connectivity": count("saved_topk_connectivity", observable),
        "observable_top1_exact": count("top1_exact_any_scan", observable),
        "observable_top1_connectivity": count("top1_connectivity_any_scan", observable),
        "observable_reported_scope_supported": count(
            "reported_scope_supported_any_scan", observable
        ),
        "observable_high_reportable_scope_supported": count(
            "high_reportable_scope_supported_any_scan", observable
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
            "observable_saved_topk_connectivity_pct": _pct(
                summary["observable_saved_topk_connectivity"], len(observable)
            ),
            "top1_connectivity_pct_observable": _pct(
                summary["observable_top1_connectivity"], len(observable)
            ),
            "reported_scope_supported_pct_observable": _pct(
                summary["observable_reported_scope_supported"], len(observable)
            ),
            "high_reportable_scope_supported_pct_observable": _pct(
                summary["observable_high_reportable_scope_supported"], len(observable)
            ),
        }
    )
    registry = load_reference_standard_registry()
    set_metadata = next(
        item for item in list_reference_standard_sets() if item["id"] == set_id
    )
    return {
        "schema_version": REFERENCE_STANDARD_QC_SCHEMA_VERSION,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "task_id": str(task_id),
        "filename": str(filename or ""),
        "reference_standard_set": set_metadata,
        "evaluation_design": "post_hoc_curated_truth_excluded_from_inference",
        "selection_source": str(selection_source or "explicit_user_selection"),
        "selection_policy": registry["selection_policy"],
        "scope": f"single {polarity} injection; not pooled across polarity or replicate files",
        "scope_zh": (
            f"单次 {polarity} 进样；不跨正负极性或重复文件合并统计"
        ),
        "truth_source": registry["truth_source_resolved"],
        "mass_tolerance_ppm_for_observable_event_audit": float(mass_tolerance_ppm),
        "used_for_candidate_generation": False,
        "used_for_scoring": False,
        "used_for_ranking": False,
        "production_results_mutated": False,
        "task_counts": {
            "ms2_spectra": len(result_rows),
            "matched_spectra": sum(isinstance(row.get("match"), Mapping) for row in result_rows),
        },
        "summary": summary,
        "denominator_warning": (
            "Declared mixture size, polarity-compatible observed standards, candidate "
            "recovery, Top-1 recovery and blind High groups are different metrics. "
            "A standard with no compatible selected MS2 event is not a search false negative."
        ),
        "denominator_warning_zh": (
            "混标声明数、本针当前极性可观测数、候选层找回数、Top1 找回数和盲检 High 组数"
            "是不同指标。当前进样没有兼容选中 MS2 事件的标准品，不计作检索假阴性。"
        ),
        "per_standard": per_standard,
    }


def reference_standard_qc_artifact_path(task_dir: Path, task_id: str) -> Path:
    return task_dir / f"{task_id}{ARTIFACT_SUFFIX}"


def write_reference_standard_qc_artifact(
    path: Path, payload: Mapping[str, Any]
) -> dict[str, Any]:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp-{os.getpid()}")
    with temporary.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, ensure_ascii=False, separators=(",", ":"), default=str)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)
    return {
        "status": "available",
        "set_id": (payload.get("reference_standard_set") or {}).get("id"),
        "artifact": str(path.resolve()),
        "size_bytes": path.stat().st_size,
        "used_for_scoring": False,
        "used_for_ranking": False,
        **dict(payload.get("summary") or {}),
    }


def read_reference_standard_qc_artifact(path: Path) -> dict[str, Any]:
    payload = _read_json(path)
    if payload.get("used_for_scoring") is not False or payload.get("used_for_ranking") is not False:
        raise ValueError("Invalid reference-standard QC artifact influence declaration")
    return payload


__all__ = [
    "REFERENCE_STANDARD_QC_SCHEMA_VERSION",
    "ARTIFACT_SUFFIX",
    "build_reference_standard_qc",
    "list_reference_standard_sets",
    "load_reference_standard_registry",
    "read_reference_standard_qc_artifact",
    "reference_standard_qc_artifact_path",
    "validate_reference_standard_set_id",
    "write_reference_standard_qc_artifact",
]
