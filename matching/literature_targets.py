"""Independent screening against literature target/diagnostic-ion lists.

This intentionally does not call the spectral-library ranking engine.  A target
hit means that precursor mass, polarity, and reported diagnostic ions were
observed; it does not mean an experimental reference spectrum or compound-level
identification was obtained.
"""

from __future__ import annotations

import bisect
import math
import sqlite3
from collections import defaultdict
from pathlib import Path
from typing import Any

from utils.adduct import adduct_polarity


def _table_exists(conn: sqlite3.Connection, table: str) -> bool:
    return conn.execute(
        "SELECT 1 FROM sqlite_master WHERE type='table' AND name=?", (table,)
    ).fetchone() is not None


def list_target_datasets(db_path: Path) -> list[dict[str, Any]]:
    conn = sqlite3.connect(f"file:{Path(db_path).resolve()}?mode=ro", uri=True)
    try:
        if not _table_exists(conn, "TargetCompounds"):
            return []
        rows = conn.execute(
            "SELECT dataset_source, COUNT(*), COUNT(DISTINCT lower(compound_name) || '|' || "
            "coalesce(molecular_formula,'')) FROM TargetCompounds "
            "WHERE source_tier='literature_presence_only' GROUP BY dataset_source "
            "ORDER BY dataset_source"
        ).fetchall()
        return [
            {"dataset_source": row[0] or "", "targets": int(row[1]), "unique_compounds": int(row[2])}
            for row in rows
        ]
    finally:
        conn.close()


def load_targets(db_path: Path, dataset_source: str) -> list[dict[str, Any]]:
    conn = sqlite3.connect(f"file:{Path(db_path).resolve()}?mode=ro", uri=True)
    conn.row_factory = sqlite3.Row
    try:
        if not (_table_exists(conn, "TargetCompounds") and _table_exists(conn, "DiagnosticIons")):
            return []
        rows = conn.execute(
            "SELECT t.id, t.target_id, t.compound_name, t.molecular_formula, t.precursor_mz, "
            "t.adduct, t.plant_source, t.quality_flags_json, d.fragment_mz "
            "FROM TargetCompounds t LEFT JOIN DiagnosticIons d ON d.target_id=t.id "
            "WHERE t.dataset_source=? AND t.source_tier='literature_presence_only' "
            "ORDER BY t.id, d.ordinal",
            (dataset_source,),
        )
        grouped: dict[int, dict[str, Any]] = {}
        for row in rows:
            target = grouped.setdefault(int(row["id"]), {
                "target_id": row["target_id"],
                "compound_name": row["compound_name"] or "",
                "formula": row["molecular_formula"] or "",
                "precursor_mz": float(row["precursor_mz"] or 0),
                "adduct": row["adduct"] or "",
                "polarity": adduct_polarity(row["adduct"]),
                "plant_source": row["plant_source"] or "",
                "quality_flags_json": row["quality_flags_json"] or "[]",
                "diagnostic_ions": [],
            })
            if row["fragment_mz"] is not None:
                target["diagnostic_ions"].append(float(row["fragment_mz"]))
        return list(grouped.values())
    finally:
        conn.close()


def _normalize_polarity(value: Any) -> str | None:
    text = str(value or "").strip().upper()
    if text.startswith("POS") or text == "+":
        return "POS"
    if text.startswith("NEG") or text == "-":
        return "NEG"
    return None


def _matched_diagnostic_ions(peaks, diagnostic_ions, tolerance: float) -> list[float]:
    observed = sorted(float(p[0]) for p in peaks or [] if p and len(p) >= 2)
    diagnostics = [float(ion) for ion in diagnostic_ions]
    # Maximum-cardinality one-to-one assignment.  A single centroid must never
    # satisfy two nearby literature ions and artificially pass a two-ion gate.
    edges: list[list[int]] = []
    for ion in diagnostics:
        left = bisect.bisect_left(observed, ion - tolerance)
        candidates: list[int] = []
        position = left
        while position < len(observed) and observed[position] <= ion + tolerance:
            candidates.append(position)
            position += 1
        candidates.sort(key=lambda index: abs(observed[index] - ion))
        edges.append(candidates)

    observed_to_diagnostic: dict[int, int] = {}

    def assign(diagnostic_index: int, seen: set[int]) -> bool:
        for observed_index in edges[diagnostic_index]:
            if observed_index in seen:
                continue
            seen.add(observed_index)
            previous = observed_to_diagnostic.get(observed_index)
            if previous is None or assign(previous, seen):
                observed_to_diagnostic[observed_index] = diagnostic_index
                return True
        return False

    assigned: set[int] = set()
    for diagnostic_index in range(len(diagnostics)):
        if assign(diagnostic_index, set()):
            assigned.add(diagnostic_index)
    # Augmenting paths can move an earlier assignment; derive the final set.
    assigned = set(observed_to_diagnostic.values())
    return [ion for index, ion in enumerate(diagnostics) if index in assigned]


def _compound_key(target: dict[str, Any]) -> str:
    return (
        f"{str(target.get('compound_name') or '').strip().casefold()}|"
        f"{str(target.get('formula') or '').strip()}"
    )


def _evidence_signature(target: dict[str, Any]) -> tuple[Any, ...]:
    """Signature for targets that MS/MS cannot distinguish under this table.

    Literature rows for positional isomers often repeat the same precursor and
    diagnostic-ion list under several names.  Counting each name as an
    independently identified compound would overstate the evidence.
    """
    return (
        _normalize_polarity(target.get("polarity")) or "",
        round(float(target.get("precursor_mz") or 0), 3),
        tuple(sorted(round(float(ion), 2) for ion in target.get("diagnostic_ions") or [])),
    )


def _hit_rank(hit: dict[str, Any]) -> tuple[float, float, float]:
    return (
        float(len(hit.get("matched_diagnostic_ions") or [])),
        float(hit.get("diagnostic_coverage") or 0),
        -abs(float(hit.get("precursor_ppm") or 0)),
    )


def screen_results(
    results: list[dict[str, Any]],
    targets: list[dict[str, Any]],
    precursor_ppm: float = 20.0,
    fragment_tolerance_da: float = 0.05,
    min_diagnostic_ions: int = 2,
    min_diagnostic_coverage: float = 0.25,
) -> dict[str, Any]:
    spectra = []
    for row in results or []:
        precursor = row.get("precursor_mz")
        try:
            precursor = float(precursor)
        except (TypeError, ValueError):
            continue
        spectra.append({
            "spectrum_id": row.get("spectrum_id", ""),
            "precursor_mz": precursor,
            "retention_time": row.get("retention_time"),
            "polarity": _normalize_polarity(row.get("ion_mode")),
            # Prefer the same denoised/top-N query spectrum used by the library
            # search.  Matching against every raw centroid in a dense LC-MS/MS
            # scan makes two-ion literature rules pass by chance.
            "peaks": row.get("query_peaks") or row.get("raw_peaks") or [],
        })

    sample_polarities = {
        spectrum["polarity"] for spectrum in spectra if spectrum.get("polarity")
    }
    screenable_targets = [
        target for target in targets
        if len(target.get("diagnostic_ions") or []) >= max(
            min_diagnostic_ions,
            int(math.ceil(
                len(target.get("diagnostic_ions") or []) * min_diagnostic_coverage
            )),
        )
    ]
    names_by_signature: dict[tuple[Any, ...], set[str]] = defaultdict(set)
    plants_by_signature: dict[tuple[Any, ...], set[str]] = defaultdict(set)
    formulas_by_signature: dict[tuple[Any, ...], set[str]] = defaultdict(set)
    adducts_by_signature: dict[tuple[Any, ...], set[str]] = defaultdict(set)
    for target in screenable_targets:
        signature = _evidence_signature(target)
        if str(target.get("compound_name") or "").strip():
            names_by_signature[signature].add(str(target["compound_name"]).strip())
        if str(target.get("plant_source") or "").strip():
            plants_by_signature[signature].add(str(target["plant_source"]).strip())
        if str(target.get("formula") or "").strip():
            formulas_by_signature[signature].add(str(target["formula"]).strip())
        if str(target.get("adduct") or "").strip():
            adducts_by_signature[signature].add(str(target["adduct"]).strip())

    polarity_eligible_keys: set[str] = set()
    precursor_keys: set[str] = set()
    one_ion_keys: set[str] = set()
    polarity_eligible_rows = 0
    precursor_rows = 0
    one_ion_rows = 0
    hits: list[dict[str, Any]] = []
    for target in targets:
        target_mz = float(target.get("precursor_mz") or 0)
        if target_mz <= 0:
            continue
        target_polarity = _normalize_polarity(target.get("polarity"))
        diagnostic = target.get("diagnostic_ions") or []
        required_diagnostic_ions = max(
            min_diagnostic_ions,
            int(math.ceil(len(diagnostic) * min_diagnostic_coverage)),
        )
        if len(diagnostic) < required_diagnostic_ions:
            # Do not silently weaken a two-ion rule to one/zero just because a
            # target row is incomplete.  Such rows need curation, not a hit.
            continue
        compound_key = _compound_key(target)
        signature = _evidence_signature(target)
        polarity_eligible = (
            (not target_polarity and bool(spectra))
            or target_polarity in sample_polarities
        )
        if polarity_eligible:
            polarity_eligible_rows += 1
            polarity_eligible_keys.add(compound_key)
        best = None
        row_has_precursor = False
        row_has_one_ion = False
        for spectrum in spectra:
            if target_polarity and target_polarity != spectrum["polarity"]:
                continue
            ppm = (spectrum["precursor_mz"] - target_mz) / target_mz * 1e6
            if abs(ppm) > precursor_ppm:
                continue
            row_has_precursor = True
            matched_ions = _matched_diagnostic_ions(
                spectrum["peaks"], diagnostic, fragment_tolerance_da
            )
            if matched_ions:
                row_has_one_ion = True
            if len(matched_ions) < required_diagnostic_ions:
                continue
            coverage = len(matched_ions) / len(diagnostic) if diagnostic else 0.0
            candidate = {
                "target_id": target.get("target_id"),
                "compound_name": target.get("compound_name", ""),
                "formula": target.get("formula", ""),
                "adduct": target.get("adduct", ""),
                "plant_source": target.get("plant_source", ""),
                "target_precursor_mz": target_mz,
                "spectrum_id": spectrum["spectrum_id"],
                "sample_precursor_mz": spectrum["precursor_mz"],
                "retention_time": spectrum["retention_time"],
                "precursor_ppm": round(ppm, 3),
                "matched_diagnostic_ions": matched_ions,
                "diagnostic_ion_count": len(diagnostic),
                "required_diagnostic_ions": required_diagnostic_ions,
                "diagnostic_coverage": round(coverage, 3),
                "evidence_level": "literature_target_screening_only",
                "_evidence_signature": signature,
            }
            rank = _hit_rank(candidate)
            if best is None or rank > best[0]:
                best = (rank, candidate)
        if row_has_precursor:
            precursor_rows += 1
            precursor_keys.add(compound_key)
        if row_has_one_ion:
            one_ion_rows += 1
            one_ion_keys.add(compound_key)
        if best:
            hits.append(best[1])

    best_by_compound: dict[str, dict[str, Any]] = {}
    target_compounds = set()
    for target in targets:
        target_compounds.add(_compound_key(target))
    for hit in hits:
        key = _compound_key(hit)
        previous = best_by_compound.get(key)
        rank = _hit_rank(hit)
        if previous is None or rank > _hit_rank(previous):
            best_by_compound[key] = hit

    # Collapse names that share indistinguishable precursor/diagnostic evidence
    # into one displayed feature, while retaining the nominal candidate count.
    best_by_signature: dict[tuple[Any, ...], dict[str, Any]] = {}
    for hit in hits:
        signature = hit["_evidence_signature"]
        previous = best_by_signature.get(signature)
        if previous is None or _hit_rank(hit) > _hit_rank(previous):
            best_by_signature[signature] = hit
    evidence_groups: list[dict[str, Any]] = []
    for signature, hit in best_by_signature.items():
        grouped = {k: v for k, v in hit.items() if k != "_evidence_signature"}
        names = sorted(names_by_signature.get(signature) or {grouped.get("compound_name", "")})
        names = [name for name in names if name]
        plants = sorted(plants_by_signature.get(signature) or set())
        formulas = sorted(formulas_by_signature.get(signature) or set())
        adducts = sorted(adducts_by_signature.get(signature) or set())
        grouped["ambiguous_target_names"] = names
        grouped["nominal_compound_count"] = len(names)
        grouped["isomer_ambiguous"] = len(names) > 1
        grouped["plant_sources"] = plants
        grouped["ambiguous_target_formulas"] = formulas
        grouped["ambiguous_target_adducts"] = adducts
        if len(names) > 1:
            grouped["compound_name"] = " / ".join(names)
        if len(formulas) > 1:
            grouped["formula"] = " / ".join(formulas)
        if len(adducts) > 1:
            grouped["adduct"] = " / ".join(adducts)
        evidence_groups.append(grouped)

    for hit in best_by_compound.values():
        hit.pop("_evidence_signature", None)

    evidence_groups.sort(
        key=lambda item: (
            -len(item["matched_diagnostic_ions"]),
            -item["diagnostic_coverage"],
            abs(item["precursor_ppm"]),
        )
    )

    return {
        "target_rows": len(targets),
        "target_unique_compounds": len(target_compounds),
        "screenable_target_rows": len(screenable_targets),
        "target_distinguishable_groups": len(names_by_signature),
        "matched_target_rows": len(hits),
        "matched_unique_compounds": len(best_by_compound),
        "matched_distinguishable_groups": len(evidence_groups),
        "screening_funnel": {
            "sample_spectra": len(spectra),
            "sample_polarities": sorted(sample_polarities),
            "polarity_eligible_target_rows": polarity_eligible_rows,
            "polarity_eligible_unique_compounds": len(polarity_eligible_keys),
            "precursor_matched_target_rows": precursor_rows,
            "precursor_matched_unique_compounds": len(precursor_keys),
            "one_diagnostic_ion_target_rows": one_ion_rows,
            "one_diagnostic_ion_unique_compounds": len(one_ion_keys),
            "minimum_diagnostic_ions_target_rows": len(hits),
            "minimum_diagnostic_ions_unique_compounds": len(best_by_compound),
        },
        "criteria": {
            "precursor_ppm": precursor_ppm,
            "fragment_tolerance_da": fragment_tolerance_da,
            "min_diagnostic_ions": min_diagnostic_ions,
            "min_diagnostic_coverage": min_diagnostic_coverage,
        },
        "interpretation": (
            "Target-screening evidence only; not an experimental library identification."
        ),
        "matches": evidence_groups,
    }
