"""Compound Discoverer-style expected-compound screening of continuous MS1.

The ordinary MSsystem search starts from acquired MS/MS spectra.  That is the
right input for spectral-library identification, but it cannot recover a
reported compound that was visible in MS1 and never selected for DDA.  This
module implements a separate, explicitly targeted evidence layer:

* read presence-only ``TargetCompounds`` rows without treating them as spectra;
* generate reported and formula/exact-mass ion hypotheses;
* stream every survey MS1 scan once and build sparse XICs;
* detect chromatographic peaks using explicit intensity/S/N/scan-count gates;
* evaluate the candidate formula against the isotope envelope at the apex;
* associate existing MS2 result rows and reported diagnostic ions.

The returned features never enter the ordinary spectral-library ``matched``
count.  Exact mass alone is not an isomer-specific identification, and names
sharing the same formula/adduct/m/z are retained as one ambiguity group.
"""

from __future__ import annotations

import bisect
import hashlib
import json
import math
import os
import re
import sqlite3
import tempfile
from collections import Counter, defaultdict
from functools import lru_cache
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import numpy as np

from matching.library_quality import (
    formula_exact_mass,
    normalize_formula_text,
    normalize_inchikey,
)
from matching.literature_targets import _matched_diagnostic_ions
from utils.adduct import adduct_polarity, ion_mz
from utils.isotope import evaluate_formula_isotope_consistency


DEFAULT_MASS_PPM = 5.0
DEFAULT_MINIMUM_INTENSITY = 1000.0
DEFAULT_SIGNAL_TO_NOISE = 1.5
DEFAULT_MINIMUM_SCAN_POINTS = 5
DEFAULT_MS2_PRECURSOR_PPM = 10.0
DEFAULT_DIAGNOSTIC_TOLERANCE_DA = 0.02
DEFAULT_MINIMUM_DIAGNOSTIC_IONS = 2
DEFAULT_MINIMUM_DIAGNOSTIC_COVERAGE = 0.25
DEFAULT_HIGH_ACCURACY_FRAGMENT_TOLERANCE_DA = 0.0025
DEFAULT_MINIMUM_HIGH_ACCURACY_FRAGMENT_PAIRS = 3
DEFAULT_MINIMUM_HIGH_ACCURACY_LIBRARY_COVERAGE = 0.50
DEFAULT_REPORTED_PRECURSOR_CONFLICT_PPM = 20.0
DEFAULT_MAXIMUM_GAP_RATIO = 0.35
DEFAULT_MAXIMUM_FWHM_SECONDS = 60.0
DEFAULT_MINIMUM_RELATIVE_VALLEY_DEPTH = 0.10
DEFAULT_OPENMS_RT_TOLERANCE_MINUTES = 0.05
ELECTRON_MASS_DA = 0.000548579909065

SAFE_COMMON_ADDUCTS = {
    "POS": ("[M+H]+", "[M+Na]+", "[M+NH4]+"),
    # Without acquisition/mobile-phase metadata, formate and acetate adducts
    # are not safe defaults.  They remain executable when explicitly reported
    # by the source catalogue; the fallback for an unspecified NEG target is
    # limited to the deprotonated ion.
    "NEG": ("[M-H]-",),
}

_HARD_TARGET_QUALITY_FLAGS = frozenset({
    "identity_status:rejected_identity_formula_mismatch",
})

# Each of the eight formula catalogues is sample-specific.  An empty target
# selection must never expand to every catalogue: doing so turns a targeted
# workflow into a cross-prescription mass lookup and produces biologically
# impossible candidate groups.  The filename resolver is intentionally strict
# (token boundaries, not a loose substring) so that, for example, ``XF0012``
# cannot be interpreted as formula F001.
FORMULA_EXPECTED_DATASET_MAP = {
    "F001": "F001_独活寄生汤_targeted_20260817",
    "F002": "F002_四神煎_targeted_20260817",
    "F003": "F003_乌头汤_targeted_20260817",
    "F004": "F004_双合汤_targeted_20260817",
    "F005": "XBT_targeted_20260817",
    "F006": "STZYT_targeted_20260817",
    "F007": "F007_HQGZWW_targeted_20260817",
    "F008": "F008_QHSS_targeted_20260817",
}

_FORMULA_FILENAME_TOKEN_RE = re.compile(
    r"(?:^|[_\-\s])(?P<formula>F00[1-8])(?=$|[_\-\s.])",
    re.IGNORECASE,
)
EXPERIMENTAL_LIBRARY_TIERS = frozenset({
    "experimental", "merged", "curated", "massbank_supplement",
})
EVIDENCE_LEVELS = (
    "MS1_EXACT_MASS_CANDIDATE",
    "MS1_ISOTOPE_SUPPORTED",
    "MS1_MS2_DIAGNOSTIC_PROBABLE",
    "EXPERIMENTAL_LIBRARY_CONFIRMED",
)
# The last two codes are retained for backward-compatible files and API
# consumers.  Their V15 display semantics are screening/library-candidate
# support, never authentic-standard confirmation.
EVIDENCE_DISPLAY_LABELS = {
    "MS1_EXACT_MASS_CANDIDATE": "MS1 exact-mass candidate",
    "MS1_ISOTOPE_SUPPORTED": "MS1 formula/isotope support",
    "MS1_MS2_DIAGNOSTIC_PROBABLE": "MS1 + reported diagnostic-ion screening support",
    "EXPERIMENTAL_LIBRARY_CONFIRMED": "MS/MS experimental-library candidate support",
}

COMPOUND_DISCOVERER_METHOD_REFERENCE = {
    "reference_software": "Thermo Scientific Compound Discoverer 3.5",
    "reference_workflows": [
        "Find Expected Compounds",
        "Detect Features",
        "Search mzCloud",
    ],
    "reference_urls": [
        (
            "https://docs.thermofisher.com/r/"
            "Compound-Discoverer-3.5-User-Guide/"
            "2011680011v2en-US2022809867"
        ),
        (
            "https://docs.thermofisher.com/r/"
            "Compound-Discoverer-3.5-User-Guide/"
            "2245584907v2en-US2248193419"
        ),
        (
            "https://docs.thermofisher.com/r/"
            "Compound-Discoverer-3.5-User-Guide/"
            "2011823115v1en-US2027019787"
        ),
    ],
    "relationship": (
        "independent open implementation benchmarked against published "
        "workflow semantics; not Thermo code and not asserted to reproduce "
        "Compound Discoverer feature-for-feature"
    ),
    "compound_discoverer_reference_defaults": {
        "expected_mass_tolerance_ppm": 5.0,
        "minimum_peak_intensity": 1000.0,
        "minimum_signal_to_noise": 1.5,
        "minimum_total_isotope_peaks": 2,
        "maximum_peak_gap_ratio": 0.35,
        "maximum_fwhm_minutes": 1.0,
        "minimum_relative_valley_depth": 0.10,
        "minimum_scan_points": 5,
        "find_expected_ms2_target_to_feature_tolerance_da": 0.025,
        "mssystem_high_resolution_ms2_association_ppm": 10.0,
        "detect_features_default_minimum_peak_intensity": 10000.0,
        "q_exactive_detect_features_recommended_minimum_peak_intensity_range": [
            500000.0,
            1000000.0,
        ],
        "high_accuracy_fragment_tolerance_da": 0.0025,
    },
}


def _available_dataset_source_names(
    available_sources: Sequence[str | Mapping[str, Any]] | None,
) -> list[str]:
    names: list[str] = []
    for item in available_sources or []:
        value = (
            item.get("dataset_source")
            if isinstance(item, Mapping)
            else item
        )
        text = str(value or "").strip()
        if text and text not in names:
            names.append(text)
    return names


def infer_expected_dataset_sources(
    filename: str | Path,
    available_sources: Sequence[str | Mapping[str, Any]] | None,
) -> list[str]:
    """Infer exactly one formula target catalogue from a sample filename.

    Only the explicit F001--F008 token is used.  No match, or a mapped
    catalogue that is not present in the active database, returns an empty
    list; callers must then fail closed instead of searching all catalogues.
    """

    match = _FORMULA_FILENAME_TOKEN_RE.search(Path(str(filename)).name)
    if not match:
        return []
    dataset = FORMULA_EXPECTED_DATASET_MAP.get(
        match.group("formula").upper(), ""
    )
    available = set(_available_dataset_source_names(available_sources))
    return [dataset] if dataset and dataset in available else []


def resolve_expected_dataset_selection(
    filename: str | Path,
    requested_sources: Sequence[str] | None,
    available_sources: Sequence[str | Mapping[str, Any]] | None,
    *,
    allow_filename_autoselect: bool = False,
) -> dict[str, Any]:
    """Resolve an auditable, fail-closed target-catalogue selection.

    A filename token such as ``F001`` is not a sample-type declaration: it can
    identify a mixed-standard injection, a sequence position, or a formula
    sample.  Production therefore requires an explicit catalogue selection.
    Filename inference remains available only as an opt-in migration helper.
    """

    requested = list(dict.fromkeys(
        str(item).strip()
        for item in (requested_sources or [])
        if str(item).strip()
    ))
    available = _available_dataset_source_names(available_sources)
    available_set = set(available)
    if requested:
        unknown = [item for item in requested if item not in available_set]
        return {
            "requested_dataset_sources": requested,
            "selected_dataset_sources": [
                item for item in requested if item in available_set
            ],
            "selection_source": "user_explicit",
            "selection_status": (
                "invalid_unknown_requested_sources"
                if unknown else "resolved_explicit"
            ),
            "inferred_formula_code": "",
            "unknown_requested_sources": unknown,
        }

    name = Path(str(filename)).name
    match = _FORMULA_FILENAME_TOKEN_RE.search(name)
    formula_code = match.group("formula").upper() if match else ""
    mapped = FORMULA_EXPECTED_DATASET_MAP.get(formula_code, "")
    if not allow_filename_autoselect:
        return {
            "requested_dataset_sources": [],
            "selected_dataset_sources": [],
            "selection_source": "none",
            "selection_status": "unresolved_explicit_selection_required",
            "inferred_formula_code": formula_code,
            "mapped_dataset_source": mapped,
            "filename_inference_applied": False,
            "unknown_requested_sources": [],
        }
    selected = infer_expected_dataset_sources(name, available)
    if selected:
        status = "resolved_from_filename"
    elif mapped:
        status = "unresolved_inferred_dataset_unavailable"
    else:
        status = "unresolved_no_selection"
    return {
        "requested_dataset_sources": [],
        "selected_dataset_sources": selected,
        "selection_source": (
            "filename_autoselect" if mapped else "none"
        ),
        "selection_status": status,
        "inferred_formula_code": formula_code,
        "mapped_dataset_source": mapped,
        "filename_inference_applied": True,
        "unknown_requested_sources": [],
    }


def _number(value: Any) -> float | None:
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        return None
    return parsed if math.isfinite(parsed) else None


def _normalise_polarity(value: Any) -> str | None:
    text = str(value or "").strip().upper()
    if text.startswith("POS") or text == "+":
        return "POS"
    if text.startswith("NEG") or text == "-":
        return "NEG"
    return None


def _normalise_adduct(value: Any) -> str:
    """Return only an exact, executable single-molecule adduct.

    Literature tables frequently contain m/z-prefixed strings, partial
    brackets, neutral-loss notes and dimer notation in an ``adduct`` column.
    Those strings remain available as reported metadata but must never become
    active extraction hypotheses.
    """

    text = re.sub(r"\s+", "", str(value or "")).upper()
    aliases = {
        "[M]+": "[M]+", "M+": "[M]+",
        "[M]-": "[M]-", "M-": "[M]-",
        "[M+H]+": "[M+H]+", "M+H": "[M+H]+",
        "[M+NA]+": "[M+Na]+", "M+NA": "[M+Na]+",
        "[M+K]+": "[M+K]+", "M+K": "[M+K]+",
        "[M+NH4]+": "[M+NH4]+", "M+NH4": "[M+NH4]+",
        "[M+H-H2O]+": "[M+H-H2O]+",
        "[M-H2O+H]+": "[M+H-H2O]+",
        "[M+H-NH3]+": "[M+H-NH3]+",
        "[M-H]-": "[M-H]-", "M-H": "[M-H]-",
        "[M+CL]-": "[M+Cl]-", "M+CL": "[M+Cl]-",
        "[M+FA-H]-": "[M+HCOO]-", "[M+HCOO]-": "[M+HCOO]-",
        "[M+HAC-H]-": "[M+CH3COO]-",
        "[M+AC-H]-": "[M+CH3COO]-",
        "[M+CH3COO]-": "[M+CH3COO]-",
        "[M-H-H2O]-": "[M-H-H2O]-",
    }
    return aliases.get(text, "")


def _table_exists(connection: sqlite3.Connection, table: str) -> bool:
    return connection.execute(
        "SELECT 1 FROM sqlite_master WHERE type='table' AND name=?",
        (table,),
    ).fetchone() is not None


def _json_list(value: Any) -> list:
    if isinstance(value, list):
        return value
    try:
        parsed = json.loads(str(value or "[]"))
    except (TypeError, ValueError, json.JSONDecodeError):
        return []
    if isinstance(parsed, list):
        return parsed
    if isinstance(parsed, dict) and isinstance(parsed.get("flags"), list):
        return parsed["flags"]
    return []


@lru_cache(maxsize=50000)
def _target_structure_evidence(smiles: str) -> dict[str, Any]:
    """Validate one reported target structure without making it authoritative."""

    text = str(smiles or "").strip()
    if not text:
        return {"status": "not_reported", "formula": "", "inchikey": ""}
    try:
        from rdkit import Chem, RDLogger
        from rdkit.Chem import rdMolDescriptors
    except ImportError:
        return {
            "status": "validator_unavailable",
            "formula": "",
            "inchikey": "",
        }
    RDLogger.DisableLog("rdApp.error")
    mol = Chem.MolFromSmiles(text)
    if mol is None:
        return {"status": "invalid", "formula": "", "inchikey": ""}
    try:
        inchikey = normalize_inchikey(Chem.MolToInchiKey(mol))
    except Exception:
        inchikey = ""
    return {
        "status": "valid",
        "formula": normalize_formula_text(re.sub(
            r"[+-](?:\d+)?$", "", rdMolDescriptors.CalcMolFormula(mol)
        )),
        "inchikey": inchikey,
        "formal_charge": int(Chem.GetFormalCharge(mol)),
    }


def _formula_only_intrinsic_charge(target: Mapping[str, Any]) -> int:
    """Retain a checked catalogue charge without selecting one stereoisomer."""
    flags = set(target.get("quality_flags") or [])
    if (
        target.get("source_tier") != "literature_presence_only"
        or "calculated_ion_hypothesis" not in flags
    ):
        return 0
    declarations = [z for z in (-1, 1) if f"formal_charge:{z}" in flags]
    if len(declarations) != 1:
        return 0
    charge = declarations[0]
    if _normalise_adduct(target.get("adduct")) != ("[M]+" if charge > 0 else "[M]-"):
        return 0
    mass = formula_exact_mass(normalize_formula_text(target.get("formula")))
    precursor = _number(target.get("precursor_mz"))
    if mass is None or precursor is None:
        return 0
    if abs(precursor - (mass - charge * ELECTRON_MASS_DA)) > 1e-6:
        return 0
    return charge


def _target_runtime_activation(target: Mapping[str, Any]) -> dict[str, Any]:
    """Return a non-destructive runtime gate for one presence-only row."""

    flags = {
        str(value or "").strip()
        for value in target.get("quality_flags") or []
        if str(value or "").strip()
    }
    raw_formula = str(
        target.get("reported_formula") or target.get("formula") or ""
    ).strip()
    formula = normalize_formula_text(target.get("formula"))
    formula_mass = formula_exact_mass(formula)
    exact_mass = _number(target.get("exact_mass"))
    reported_precursor = (
        _number(target.get("reported_precursor_mz"))
        if _number(target.get("reported_precursor_mz")) is not None
        else _number(target.get("precursor_mz"))
    )
    reasons: list[str] = []
    for flag in sorted(flags & _HARD_TARGET_QUALITY_FLAGS):
        reasons.append(flag)
    if raw_formula and formula_mass is None:
        reasons.append("invalid_nonempty_molecular_formula")
    if formula_mass is None and exact_mass is None and reported_precursor is None:
        reasons.append("no_usable_formula_exact_mass_or_precursor")
    active_adduct = _normalise_adduct(target.get("adduct"))
    if (
        formula_mass is None
        and exact_mass is None
        and reported_precursor is not None
        and (
            not active_adduct
            or adduct_polarity(active_adduct) not in {"POS", "NEG"}
        )
    ):
        # A bare source m/z has no executable polarity/ion identity.  Calling
        # such a row active previously created a false catalogue KPI even
        # though neither POS nor NEG could generate an XIC hypothesis.
        reasons.append("reported_precursor_without_executable_adduct")

    structure_status = str(
        target.get("structure_validation_status") or "not_reported"
    )
    identity_conflict = "identity_status:identity_conflict" in flags
    rejected_identity = any(
        value in flags for value in _HARD_TARGET_QUALITY_FLAGS
    )
    structure_conflict = structure_status in {
        "invalid_reported_structure_quarantined",
        "structure_formula_conflict_quarantined",
        "structure_connectivity_conflict_quarantined",
    }
    active = not reasons
    verified_identity = bool(target.get("structure_identity_verified"))
    exact_identity_upgrade_eligible = bool(
        active
        and verified_identity
        and normalize_inchikey(target.get("inchikey"))
        and not identity_conflict
        and not rejected_identity
        and not structure_conflict
    )
    identity_claim_scope = (
        "runtime_quarantined"
        if not active
        else "formula_or_mass_group_only"
        if identity_conflict or structure_conflict
        else "verified_full_identity"
        if exact_identity_upgrade_eligible
        else "presence_only_candidate"
    )
    return {
        "runtime_active": active,
        "runtime_quarantine_reasons": sorted(set(reasons)),
        "identity_claim_scope": identity_claim_scope,
        "exact_identity_upgrade_eligible": exact_identity_upgrade_eligible,
        "formula_mass_available": formula_mass is not None,
        "reported_exact_mass_available": exact_mass is not None,
        "reported_precursor_available": reported_precursor is not None,
    }


def target_catalog_quality_summary(
    targets: Sequence[Mapping[str, Any]],
) -> dict[str, Any]:
    reason_counts: Counter[str] = Counter()
    scope_counts: Counter[str] = Counter()
    active = 0
    for target in targets:
        is_active = target.get("runtime_active") is not False
        active += int(is_active)
        reason_counts.update(
            str(value)
            for value in target.get("runtime_quarantine_reasons") or []
        )
        scope_counts[str(
            target.get("identity_claim_scope") or "not_evaluated"
        )] += 1
    total = len(targets)
    return {
        "policy": "mssystem-v12-presence-target-runtime-quarantine-1.0.0",
        "source_rows_preserved": total,
        "active_target_rows": active,
        "quarantined_target_rows": total - active,
        "quarantine_reason_counts": dict(sorted(reason_counts.items())),
        "identity_claim_scope_counts": dict(sorted(scope_counts.items())),
        "database_rows_deleted": 0,
        "database_rows_rewritten": 0,
    }


def list_expected_datasets(db_path: str | Path) -> list[dict[str, Any]]:
    """List available presence-only catalogues without loading spectra."""

    path = Path(db_path).expanduser().resolve()
    connection = sqlite3.connect(f"file:{path}?mode=ro", uri=True)
    try:
        if not _table_exists(connection, "TargetCompounds"):
            return []
        rows = connection.execute(
            "SELECT dataset_source, COUNT(*), "
            "COUNT(DISTINCT lower(trim(compound_name)) || '|' || "
            "coalesce(molecular_formula,'')) "
            "FROM TargetCompounds "
            "WHERE source_tier='literature_presence_only' "
            "GROUP BY dataset_source ORDER BY dataset_source"
        )
        return [
            {
                "dataset_source": str(source or ""),
                "target_rows": int(row_count or 0),
                "named_formula_entities": int(entity_count or 0),
            }
            for source, row_count, entity_count in rows
        ]
    finally:
        connection.close()


def load_expected_compound_catalog(
    db_path: str | Path,
    dataset_sources: Sequence[str] | None = None,
) -> list[dict[str, Any]]:
    """Read the presence-only target catalog and its diagnostic-ion lists.

    The SQLite connection is opened in read-only mode.  No target row is
    promoted to a library spectrum by this function.
    """

    sources = sorted({
        str(item).strip() for item in (dataset_sources or []) if str(item).strip()
    })
    # Fail closed.  Historically an empty list omitted the SQL WHERE clause and
    # loaded every prescription catalogue, contaminating F001 with F002--F008
    # candidates.  All callers must now resolve an explicit per-sample list.
    if not sources:
        return []
    path = Path(db_path).expanduser().resolve()
    connection = sqlite3.connect(f"file:{path}?mode=ro", uri=True)
    connection.row_factory = sqlite3.Row
    try:
        if not _table_exists(connection, "TargetCompounds"):
            return []
        target_columns = {
            str(row[1])
            for row in connection.execute(
                "PRAGMA table_info(TargetCompounds)"
            )
        }
        optional_target_select = {
            name: (
                f"t.{name} AS {name}"
                if name in target_columns
                else f"NULL AS {name}"
            )
            for name in (
                "reported_precursor_mz",
                "reported_adduct",
                "reported_retention_time",
                "reported_retention_time_unit",
                "reported_collision_energy_raw",
                "instrument",
                "target_searchability",
            )
        }
        has_diagnostics = _table_exists(connection, "DiagnosticIons")
        where = ["t.source_tier='literature_presence_only'"]
        parameters: list[Any] = []
        marks = ",".join("?" for _ in sources)
        where.append(f"t.dataset_source IN ({marks})")
        parameters.extend(sources)
        diagnostic_select = (
            "d.fragment_mz AS fragment_mz, d.ordinal AS diagnostic_ordinal"
            if has_diagnostics
            else "NULL AS fragment_mz, NULL AS diagnostic_ordinal"
        )
        diagnostic_join = (
            "LEFT JOIN DiagnosticIons d ON d.target_id=t.id"
            if has_diagnostics
            else ""
        )
        query = f"""
            SELECT
                t.id, t.target_id, t.dataset_source, t.compound_name,
                t.molecular_formula, t.smiles, t.inchi, t.inchikey,
                t.reported_inchikey, t.structure_identity_verified,
                t.structure_identity_source, t.exact_mass, t.precursor_mz,
                t.adduct, t.charge, t.plant_source, t.source_file,
                t.source_content_sha256, t.dataset_version,
                t.source_row_index, t.provenance, t.doi,
                t.dataset_accession, t.spectrum_accession,
                t.quality_flags_json,
                {', '.join(optional_target_select.values())},
                {diagnostic_select}
            FROM TargetCompounds t
            {diagnostic_join}
            WHERE {' AND '.join(where)}
            ORDER BY t.id, diagnostic_ordinal
        """
        grouped: dict[int, dict[str, Any]] = {}
        for row in connection.execute(query, parameters):
            identifier = int(row["id"])
            target = grouped.setdefault(
                identifier,
                {
                    "catalog_row_id": identifier,
                    "target_id": str(row["target_id"] or ""),
                    "dataset_source": str(row["dataset_source"] or ""),
                    "compound_name": str(row["compound_name"] or "").strip(),
                    "reported_formula": str(
                        row["molecular_formula"] or ""
                    ).strip(),
                    "formula": normalize_formula_text(
                        row["molecular_formula"]
                    ),
                    "smiles": str(row["smiles"] or "").strip(),
                    "inchi": str(row["inchi"] or "").strip(),
                    "inchikey": normalize_inchikey(row["inchikey"]),
                    "reported_inchikey": normalize_inchikey(
                        row["reported_inchikey"]
                    ),
                    "structure_identity_verified": bool(
                        row["structure_identity_verified"]
                    ),
                    "reported_structure_identity_verified": bool(
                        row["structure_identity_verified"]
                    ),
                    "structure_identity_source": str(
                        row["structure_identity_source"] or ""
                    ),
                    "exact_mass": _number(row["exact_mass"]),
                    "precursor_mz": _number(row["precursor_mz"]),
                    "adduct": _normalise_adduct(row["adduct"]),
                    "raw_adduct": str(row["adduct"] or "").strip(),
                    "reported_precursor_mz": _number(
                        row["reported_precursor_mz"]
                    ),
                    "reported_adduct": str(
                        row["reported_adduct"] or row["adduct"] or ""
                    ).strip(),
                    "reported_retention_time": _number(
                        row["reported_retention_time"]
                    ),
                    "reported_retention_time_unit": str(
                        row["reported_retention_time_unit"] or ""
                    ),
                    "reported_collision_energy_raw": str(
                        row["reported_collision_energy_raw"] or ""
                    ),
                    "instrument": str(row["instrument"] or ""),
                    "target_searchability": str(
                        row["target_searchability"] or ""
                    ),
                    "charge": row["charge"],
                    "plant_source": str(row["plant_source"] or ""),
                    "source_file": str(row["source_file"] or ""),
                    "source_content_sha256": str(
                        row["source_content_sha256"] or ""
                    ),
                    "dataset_version": str(row["dataset_version"] or ""),
                    "source_row_index": row["source_row_index"],
                    "provenance": str(row["provenance"] or ""),
                    "doi": str(row["doi"] or ""),
                    "dataset_accession": str(
                        row["dataset_accession"] or ""
                    ),
                    "spectrum_accession": str(
                        row["spectrum_accession"] or ""
                    ),
                    "quality_flags": _json_list(row["quality_flags_json"]),
                    "source_tier": "literature_presence_only",
                    "diagnostic_ions": [],
                },
            )
            diagnostic = _number(row["fragment_mz"])
            if diagnostic is not None and diagnostic > 0:
                target["diagnostic_ions"].append(diagnostic)
        for target in grouped.values():
            raw_adduct = str(target.get("raw_adduct") or "").strip()
            if raw_adduct and not target.get("adduct"):
                flags = list(target.get("quality_flags") or [])
                if "invalid_reported_adduct_quarantined" not in flags:
                    flags.append("invalid_reported_adduct_quarantined")
                target["quality_flags"] = flags
                target["adduct_validation_status"] = "quarantined"
            elif target.get("adduct"):
                target["adduct_validation_status"] = "canonical_supported"
            else:
                target["adduct_validation_status"] = "not_reported"

            # A source table can mark an identity as verified while carrying a
            # malformed or formula-conflicting SMILES.  Keep the raw value for
            # audit, but do not expose such a structure as an active depiction
            # or structure-level claim.
            reported_smiles = str(target.get("smiles") or "").strip()
            target["reported_smiles"] = reported_smiles
            structure = _target_structure_evidence(reported_smiles)
            status = str(structure.get("status") or "not_reported")
            formula = normalize_formula_text(target.get("formula"))
            inchikey = normalize_inchikey(target.get("inchikey"))
            structure_formula = normalize_formula_text(structure.get("formula"))
            structure_inchikey = normalize_inchikey(structure.get("inchikey"))
            flags = list(target.get("quality_flags") or [])
            if status == "invalid":
                target["smiles"] = ""
                status = "invalid_reported_structure_quarantined"
                flags.append("invalid_reported_structure_quarantined")
            elif status == "valid" and formula and structure_formula != formula:
                target["smiles"] = ""
                status = "structure_formula_conflict_quarantined"
                flags.append("structure_formula_conflict_quarantined")
            elif (
                status == "valid" and inchikey and structure_inchikey
                and structure_inchikey[:14] != inchikey[:14]
            ):
                target["smiles"] = ""
                status = "structure_connectivity_conflict_quarantined"
                flags.append("structure_connectivity_conflict_quarantined")
            elif (
                status == "valid" and inchikey and structure_inchikey
                and structure_inchikey != inchikey
            ):
                status = (
                    "connectivity_valid_stereochemistry_or_"
                    "protonation_underspecified"
                )
                flags.append(
                    "reported_structure_stereochemistry_or_"
                    "protonation_underspecified"
                )
            elif status == "valid":
                status = "internally_consistent"
            target["structure_validation_status"] = status
            target["structure_formula"] = structure_formula
            target["structure_inchikey"] = structure_inchikey
            # Elemental composition and intrinsic ionic charge are separate.
            # Quaternary alkaloids retain their charged structure and [M]+ ion.
            target["formal_charge"] = (
                int(structure.get("formal_charge") or 0)
                if status in {"internally_consistent",
                              "connectivity_valid_stereochemistry_or_protonation_underspecified"}
                else _formula_only_intrinsic_charge(target)
                if status == "not_reported"
                else 0
            )
            if status in {
                "invalid_reported_structure_quarantined",
                "structure_formula_conflict_quarantined",
                "structure_connectivity_conflict_quarantined",
            }:
                target["structure_identity_verified"] = False
            target["quality_flags"] = sorted(set(flags))
            target.update(_target_runtime_activation(target))
        return list(grouped.values())
    finally:
        connection.close()


def _target_identity(member: Mapping[str, Any]) -> str:
    key = normalize_inchikey(member.get("inchikey"))
    if key:
        return f"ik::{key}"
    name = re.sub(
        r"[^a-z0-9]+",
        "",
        str(member.get("compound_name") or "").casefold(),
    )
    formula = normalize_formula_text(member.get("formula"))
    if name or formula:
        return f"nf::{name}|{formula}"
    return f"target::{member.get('target_id') or member.get('catalog_row_id')}"


def _member_projection(target: Mapping[str, Any]) -> dict[str, Any]:
    return {
        "target_id": str(target.get("target_id") or ""),
        "catalog_row_id": target.get("catalog_row_id"),
        "compound_name": str(target.get("compound_name") or ""),
        "formula": normalize_formula_text(target.get("formula")),
        "reported_formula": str(target.get("reported_formula") or ""),
        "smiles": str(target.get("smiles") or ""),
        "reported_smiles": str(target.get("reported_smiles") or ""),
        "inchikey": normalize_inchikey(target.get("inchikey")),
        "reported_inchikey": normalize_inchikey(
            target.get("reported_inchikey")
        ),
        "structure_identity_verified": bool(
            target.get("structure_identity_verified")
        ),
        "reported_structure_identity_verified": bool(
            target.get("reported_structure_identity_verified")
        ),
        "structure_validation_status": str(
            target.get("structure_validation_status") or "not_reported"
        ),
        "structure_formula": str(target.get("structure_formula") or ""),
        "structure_inchikey": normalize_inchikey(
            target.get("structure_inchikey")
        ),
        "formal_charge": int(target.get("formal_charge") or 0),
        "dataset_source": str(target.get("dataset_source") or ""),
        "backing_dataset_source": str(
            target.get("backing_dataset_source") or ""
        ),
        "target_polarity_scope": str(
            target.get("target_polarity_scope") or ""
        ),
        "plant_source": str(target.get("plant_source") or ""),
        "source_file": str(target.get("source_file") or ""),
        "source_content_sha256": str(
            target.get("source_content_sha256") or ""
        ),
        "dataset_version": str(target.get("dataset_version") or ""),
        "source_row_index": target.get("source_row_index"),
        "provenance": str(target.get("provenance") or ""),
        "doi": str(target.get("doi") or ""),
        "dataset_accession": str(target.get("dataset_accession") or ""),
        "spectrum_accession": str(target.get("spectrum_accession") or ""),
        "reported_precursor_mz": _number(
            target.get("reported_precursor_mz")
        ),
        "reported_adduct": str(target.get("reported_adduct") or ""),
        "active_adduct": str(target.get("adduct") or ""),
        "adduct_validation_status": str(
            target.get("adduct_validation_status") or "not_reported"
        ),
        "reported_retention_time": _number(
            target.get("reported_retention_time")
        ),
        "reported_retention_time_unit": str(
            target.get("reported_retention_time_unit") or ""
        ),
        "reported_collision_energy_raw": str(
            target.get("reported_collision_energy_raw") or ""
        ),
        "instrument": str(target.get("instrument") or ""),
        "target_searchability": str(
            target.get("target_searchability") or ""
        ),
        "quality_flags": list(target.get("quality_flags") or []),
        "runtime_active": target.get("runtime_active") is not False,
        "runtime_quarantine_reasons": list(
            target.get("runtime_quarantine_reasons") or []
        ),
        "identity_claim_scope": str(
            target.get("identity_claim_scope") or "not_evaluated"
        ),
        "exact_identity_upgrade_eligible": (
            bool(target.get("exact_identity_upgrade_eligible"))
            if "exact_identity_upgrade_eligible" in target
            else None
        ),
        "diagnostic_ions": sorted({
            round(float(value), 6)
            for value in target.get("diagnostic_ions") or []
            if _number(value) is not None and float(value) > 0
        }),
        "identity_key": _target_identity(target),
    }


def _neutral_mass_evidence(target: Mapping[str, Any]) -> dict[str, Any]:
    formula = normalize_formula_text(target.get("formula"))
    formula_mass = formula_exact_mass(formula)
    reported_mass = _number(target.get("exact_mass"))
    selected = formula_mass if formula_mass is not None else reported_mass
    ppm = None
    if formula_mass and reported_mass:
        ppm = (reported_mass - formula_mass) / formula_mass * 1e6
    return {
        "selected_mass": selected,
        "selected_source": (
            "molecular_formula"
            if formula_mass is not None
            else "reported_exact_mass"
            if reported_mass is not None
            else "unavailable"
        ),
        "formula_mass": formula_mass,
        "reported_exact_mass": reported_mass,
        "formula_exact_mass_ppm_difference": ppm,
    }


def build_expected_ion_hypotheses(
    targets: Sequence[Mapping[str, Any]],
    polarity: str,
    *,
    include_safe_common_adducts: bool = True,
    include_reported_precursor_conflicts: bool = False,
    ambiguity_ppm: float = DEFAULT_MASS_PPM,
    reported_conflict_ppm: float = (
        DEFAULT_REPORTED_PRECURSOR_CONFLICT_PPM
    ),
) -> list[dict[str, Any]]:
    """Create and collapse expected ions for one acquisition polarity."""

    mode = _normalise_polarity(polarity)
    if mode not in {"POS", "NEG"}:
        return []
    claims: list[dict[str, Any]] = []
    for raw_target in targets:
        target = dict(raw_target)
        target_scope = str(
            target.get("target_polarity_scope") or ""
        ).strip().upper()
        if target_scope and target_scope != mode:
            continue
        if target.get("runtime_active") is False:
            continue
        member = _member_projection(target)
        formula = member["formula"]
        neutral_evidence = _neutral_mass_evidence(target)
        neutral_mass = neutral_evidence["selected_mass"]
        formal_charge = int(target.get("formal_charge") or 0)
        calculated_catalog_ion = (
            "calculated_ion_hypothesis" in (target.get("quality_flags") or [])
            and _number(target.get("reported_precursor_mz")) is None
        )
        reported_mz = (
            _number(target.get("reported_precursor_mz"))
            if _number(target.get("reported_precursor_mz")) is not None
            else _number(target.get("precursor_mz"))
        )
        if calculated_catalog_ion:
            reported_mz = None
        reported_adduct = _normalise_adduct(target.get("adduct"))
        raw_reported_adduct = str(
            target.get("reported_adduct") or target.get("raw_adduct") or ""
        ).strip()
        invalid_reported_adduct = bool(
            raw_reported_adduct and not reported_adduct
        )
        reported_mode = adduct_polarity(reported_adduct)
        used_adducts: set[str] = set()

        if formal_charge:
            # A permanent monocation/anion is screened as its molecular ion.
            # Adding a proton or generating the opposite-polarity neutral
            # adduct would change both charge and chemical identity.
            if (
                abs(formal_charge) != 1
                or mode != ("POS" if formal_charge > 0 else "NEG")
                or neutral_mass is None
            ):
                continue
            ionic_mz = neutral_mass - formal_charge * ELECTRON_MASS_DA
            claims.append({
                "target_mz": ionic_mz,
                "formula": formula,
                "adduct": "[M]+" if formal_charge > 0 else "[M]-",
                "polarity": mode,
                "generation_source": "formula_intrinsic_ion_exact_mass",
                "ion_origin": "intrinsic_charge_molecular_ion",
                "reported_precursor_priority": False,
                "formula_theoretical_primary": True,
                "reported_precursor_conflict_candidate": False,
                "requires_orthogonal_evidence_for_upgrade": False,
                "theoretical_mz": ionic_mz,
                "reported_precursor_mz": reported_mz,
                "reported_theoretical_ppm": (
                    (reported_mz - ionic_mz) / ionic_mz * 1e6
                    if reported_mz is not None else None
                ),
                "ion_quality_flags": [],
                "neutral_mass_evidence": {
                    **neutral_evidence,
                    "intrinsic_formal_charge": formal_charge,
                    "electron_mass_correction_da": -formal_charge * ELECTRON_MASS_DA,
                },
                "member": member,
            })
            continue

        reported_usable = (
            reported_mz is not None
            and reported_mz > 0
            and not invalid_reported_adduct
            and (reported_mode is None or reported_mode == mode)
        )
        theoretical = (
            ion_mz(neutral_mass, reported_adduct)
            if neutral_mass is not None
            and reported_adduct
            and reported_mode == mode
            else None
        )
        reported_theoretical_ppm = (
            (reported_mz - theoretical) / theoretical * 1e6
            if reported_usable and theoretical
            else None
        )

        # The chemically reconstructed formula/adduct ion is primary whenever
        # it is available.  The source precursor is never silently erased:
        # a >20 ppm conflict becomes a separate recall-only ion hypothesis.
        if theoretical is not None and theoretical > 0:
            claims.append(
                {
                    "target_mz": theoretical,
                    "formula": formula,
                    "adduct": reported_adduct,
                    "polarity": mode,
                    "generation_source": (
                        "formula_exact_mass"
                        if neutral_evidence["selected_source"]
                        == "molecular_formula"
                        else "reported_exact_mass"
                    ),
                    "ion_origin": (
                        "calculated_catalog_adduct_theoretical_ion"
                        if calculated_catalog_ion else "reported_adduct_theoretical_ion"
                    ),
                    "reported_precursor_priority": False,
                    "formula_theoretical_primary": True,
                    "reported_precursor_conflict_candidate": False,
                    "requires_orthogonal_evidence_for_upgrade": False,
                    "theoretical_mz": theoretical,
                    "reported_precursor_mz": reported_mz,
                    "reported_theoretical_ppm": reported_theoretical_ppm,
                    "ion_quality_flags": (
                        ["reported_precursor_conflicts_with_theoretical_ion"]
                        if reported_theoretical_ppm is not None
                        and abs(reported_theoretical_ppm)
                        > float(reported_conflict_ppm)
                        else []
                    ),
                    "neutral_mass_evidence": neutral_evidence,
                    "member": member,
                }
            )
            used_adducts.add(reported_adduct.upper())
            if (
                include_reported_precursor_conflicts
                and reported_usable
                and reported_theoretical_ppm is not None
                and abs(reported_theoretical_ppm)
                > float(reported_conflict_ppm)
            ):
                claims.append(
                    {
                        "target_mz": reported_mz,
                        "formula": formula,
                        "adduct": reported_adduct,
                        "polarity": mode,
                        "generation_source": (
                            "reported_precursor_conflict_candidate"
                        ),
                        "ion_origin": (
                            "reported_precursor_conflict_candidate"
                        ),
                        "reported_precursor_priority": False,
                        "formula_theoretical_primary": False,
                        "reported_precursor_conflict_candidate": True,
                        "requires_orthogonal_evidence_for_upgrade": True,
                        "theoretical_mz": theoretical,
                        "reported_precursor_mz": reported_mz,
                        "reported_theoretical_ppm": (
                            reported_theoretical_ppm
                        ),
                        "ion_quality_flags": [
                            "reported_precursor_conflicts_with_theoretical_ion",
                            (
                                "reported_precursor_retained_as_"
                                "alternative_hypothesis"
                            ),
                        ],
                        "neutral_mass_evidence": neutral_evidence,
                        "member": member,
                    }
                )
        elif reported_usable:
            label = reported_adduct or f"reported-ion-{mode}"
            claims.append(
                {
                    "target_mz": reported_mz,
                    "formula": formula,
                    "adduct": label,
                    "polarity": mode,
                    "generation_source": "reported_precursor",
                    "ion_origin": "reported_precursor_without_theoretical_ion",
                    "reported_precursor_priority": True,
                    "formula_theoretical_primary": False,
                    "reported_precursor_conflict_candidate": False,
                    "requires_orthogonal_evidence_for_upgrade": False,
                    "theoretical_mz": None,
                    "reported_precursor_mz": reported_mz,
                    "reported_theoretical_ppm": None,
                    "ion_quality_flags": [
                        "reported_precursor_used_without_theoretical_ion"
                    ],
                    "neutral_mass_evidence": neutral_evidence,
                    "member": member,
                }
            )
            used_adducts.add(_normalise_adduct(label).upper())

        if include_safe_common_adducts and neutral_mass is not None:
            for adduct in SAFE_COMMON_ADDUCTS[mode]:
                if _normalise_adduct(adduct).upper() in used_adducts:
                    continue
                calculated_mz = ion_mz(neutral_mass, adduct)
                if calculated_mz is None:
                    continue
                claims.append(
                    {
                        "target_mz": calculated_mz,
                        "formula": formula,
                        "adduct": adduct,
                        "polarity": mode,
                        "generation_source": (
                            "formula_exact_mass"
                            if neutral_evidence["selected_source"]
                            == "molecular_formula"
                            else "reported_exact_mass"
                        ),
                        "ion_origin": "generated_common_adduct",
                        "reported_precursor_priority": False,
                        "formula_theoretical_primary": True,
                        "reported_precursor_conflict_candidate": False,
                        "requires_orthogonal_evidence_for_upgrade": False,
                        "theoretical_mz": calculated_mz,
                        "reported_precursor_mz": reported_mz,
                        "reported_theoretical_ppm": None,
                        "ion_quality_flags": [],
                        "neutral_mass_evidence": neutral_evidence,
                        "member": member,
                    }
                )

    # Cluster only chemically identical formula/adduct hypotheses within the
    # expected-compound mass tolerance.  Blank formulas retain target identity
    # in the key so unrelated rows are never collapsed merely by m/z.
    buckets: dict[
        tuple[str, str, str, bool], list[dict[str, Any]]
    ] = defaultdict(list)
    for claim in claims:
        formula_key = normalize_formula_text(claim["formula"])
        blank_guard = (
            ""
            if formula_key
            else str(claim["member"].get("identity_key") or "")
        )
        key = (
            formula_key,
            _normalise_adduct(claim["adduct"]).upper(),
            blank_guard,
            bool(claim.get("reported_precursor_conflict_candidate")),
        )
        buckets[key].append(claim)

    groups: list[dict[str, Any]] = []
    for (formula, adduct_key, _, _is_conflict), bucket in sorted(
        buckets.items()
    ):
        bucket.sort(key=lambda item: float(item["target_mz"]))
        clusters: list[list[dict[str, Any]]] = []
        for claim in bucket:
            if not clusters:
                clusters.append([claim])
                continue
            current_mzs = [float(item["target_mz"]) for item in clusters[-1]]
            centre = float(np.median(current_mzs))
            tolerance = centre * float(ambiguity_ppm) * 1e-6
            if abs(float(claim["target_mz"]) - centre) <= tolerance:
                clusters[-1].append(claim)
            else:
                clusters.append([claim])

        for cluster in clusters:
            reported = [
                item for item in cluster if item["reported_precursor_priority"]
            ]
            mass_source = reported or cluster
            target_mz = float(np.median([
                float(item["target_mz"]) for item in mass_source
            ]))
            members_by_id: dict[str, dict[str, Any]] = {}
            for item in cluster:
                member = item["member"]
                key = str(
                    member.get("catalog_row_id")
                    or member.get("target_id")
                    or member.get("identity_key")
                )
                members_by_id[key] = member
            members = sorted(
                members_by_id.values(),
                key=lambda item: (
                    str(item.get("compound_name") or "").casefold(),
                    str(item.get("target_id") or ""),
                ),
            )
            names = sorted({
                str(item.get("compound_name") or "").strip()
                for item in members
                if str(item.get("compound_name") or "").strip()
            })
            identities = {
                str(item.get("identity_key") or "") for item in members
            }
            datasets = sorted({
                str(item.get("dataset_source") or "")
                for item in members
                if str(item.get("dataset_source") or "")
            })
            sources = sorted({
                str(item["generation_source"]) for item in cluster
            })
            origins = sorted({
                str(item.get("ion_origin") or "") for item in cluster
                if str(item.get("ion_origin") or "")
            })
            quality_flags = sorted({
                str(flag)
                for item in cluster
                for flag in item.get("ion_quality_flags") or []
            })
            signature = json.dumps(
                {
                    "formula": formula,
                    "adduct": adduct_key,
                    "mz": round(target_mz, 6),
                    "targets": sorted(members_by_id),
                },
                sort_keys=True,
                ensure_ascii=False,
            )
            identifier = "ECI_" + hashlib.sha256(
                signature.encode("utf-8")
            ).hexdigest()[:16]
            groups.append(
                {
                    "expected_ion_id": identifier,
                    "target_mz": round(target_mz, 8),
                    "formula": formula,
                    "adduct": cluster[0]["adduct"],
                    "polarity": mode,
                    "charge": 1,
                    "generation_sources": sources,
                    "ion_origins": origins,
                    "ion_quality_flags": quality_flags,
                    "reported_precursor_priority": bool(reported),
                    "formula_theoretical_primary": any(
                        bool(item.get("formula_theoretical_primary"))
                        for item in cluster
                    ),
                    "reported_precursor_conflict_candidate": any(
                        bool(item.get(
                            "reported_precursor_conflict_candidate"
                        ))
                        for item in cluster
                    ),
                    "requires_orthogonal_evidence_for_upgrade": any(
                        bool(item.get(
                            "requires_orthogonal_evidence_for_upgrade"
                        ))
                        for item in cluster
                    ),
                    "alternative_target_mzs": sorted({
                        round(float(item["target_mz"]), 8)
                        for item in cluster
                    }),
                    "candidate_names": names,
                    "nominal_target_count": len(names),
                    "catalog_row_count": len(members),
                    "identity_count": len(identities),
                    "isomer_ambiguous": len(identities) > 1 or len(names) > 1,
                    "dataset_sources": datasets,
                    "members": members,
                    "neutral_mass_evidence": cluster[0][
                        "neutral_mass_evidence"
                    ],
                    "reported_formula_adduct_mz_checks": [
                        {
                            "target_id": item["member"].get("target_id"),
                            "reported_theoretical_ppm": (
                                round(
                                    float(item["reported_theoretical_ppm"]), 3
                                )
                                if item["reported_theoretical_ppm"] is not None
                                else None
                            ),
                        }
                        for item in cluster
                        if item["reported_theoretical_ppm"] is not None
                    ],
                }
            )
    return sorted(
        groups,
        key=lambda item: (
            float(item["target_mz"]),
            str(item["formula"]),
            str(item["adduct"]),
            str(item["expected_ion_id"]),
        ),
    )


def _scan_rt_minutes(spectrum: Mapping[str, Any]) -> float | None:
    scans = ((spectrum.get("scanList") or {}).get("scan") or [])
    scan = scans[0] if scans else {}
    raw = scan.get("scan start time")
    value = _number(raw)
    if value is None:
        return None
    unit_text = " ".join(
        (
            str(getattr(raw, "unit_info", "") or ""),
            str(getattr(raw, "unit_name", "") or ""),
            str(scan),
        )
    ).lower()
    if "second" in unit_text and "minute" not in unit_text:
        return value / 60.0
    return value


def _peak_near(
    mzs: np.ndarray,
    intensities: np.ndarray,
    target: float,
    tolerance: float,
) -> tuple[float | None, float]:
    start = int(np.searchsorted(mzs, target - tolerance, side="left"))
    end = int(np.searchsorted(mzs, target + tolerance, side="right"))
    if end <= start:
        return None, 0.0
    local = int(np.argmax(intensities[start:end])) + start
    return float(mzs[local]), float(intensities[local])


def _isotope_envelope_at_scan(
    mzs: np.ndarray,
    intensities: np.ndarray,
    observed_mz: float,
    observed_intensity: float,
    *,
    charge: int = 1,
    mass_ppm: float = DEFAULT_MASS_PPM,
) -> dict[str, Any] | None:
    if observed_intensity <= 0:
        return None
    z = max(1, abs(int(charge or 1)))
    spacing = 1.00335483507 / z
    # Keep isotope extraction on the same high-resolution mass scale as the
    # expected-ion XIC.  The former 0.01 Da/10 ppm minimum could turn a nearby
    # co-eluting ion into false isotope support at low m/z.
    tolerance = max(0.0005, observed_mz * float(mass_ppm) * 1e-6)
    previous_mz, previous_intensity = _peak_near(
        mzs,
        intensities,
        observed_mz - spacing,
        tolerance,
    )
    previous_ratio = (
        previous_intensity / observed_intensity
        if previous_mz is not None and observed_intensity > 0
        else None
    )
    pattern_integrity = (
        "possible_subpattern"
        if previous_ratio is not None and previous_ratio >= 0.20
        else "passed"
    )
    peaks = [
        {
            "order": 0,
            "mz": round(observed_mz, 6),
            "intensity": round(observed_intensity, 3),
            "relative_to_m": 1.0,
        }
    ]
    for order in (1, 2, 3):
        found_mz, found_intensity = _peak_near(
            mzs,
            intensities,
            observed_mz + spacing * order,
            tolerance,
        )
        peaks.append(
            {
                "order": order,
                "mz": round(found_mz, 6) if found_mz is not None else None,
                "intensity": (
                    round(found_intensity, 3)
                    if found_mz is not None
                    else 0.0
                ),
                "relative_to_m": (
                    round(found_intensity / observed_intensity, 5)
                    if found_mz is not None
                    else None
                ),
            }
        )
    observed = sum(item["mz"] is not None for item in peaks[1:])
    return {
        "source": "adjacent_ms1_scan",
        "source_detail": "expected_compound_apex_ms1",
        "status": "observed" if observed else "monoisotopic_peak_only",
        "charge": z,
        "spacing_da": round(spacing, 6),
        "mz_tolerance_da": round(tolerance, 6),
        "mz_tolerance_ppm": round(float(mass_ppm), 6),
        "mz_tolerance_source": "expected_compound_mass_ppm_with_centroid_rounding_floor",
        "n_isotopes_observed": observed,
        "peaks": peaks,
        "preceding_isotope_peak": {
            "mz": round(previous_mz, 6) if previous_mz is not None else None,
            "intensity": (
                round(previous_intensity, 3)
                if previous_mz is not None
                else 0.0
            ),
            "relative_to_m": (
                round(previous_ratio, 5)
                if previous_ratio is not None
                else None
            ),
        },
        "pattern_integrity": pattern_integrity,
    }


def extract_sparse_expected_xics(
    spectra: Iterable[Mapping[str, Any]],
    ion_hypotheses: Sequence[Mapping[str, Any]],
    *,
    mass_ppm: float = DEFAULT_MASS_PPM,
) -> dict[str, Any]:
    """Extract all expected-ion XICs in one pass over survey MS1 spectra."""

    ordered = sorted(
        enumerate(ion_hypotheses),
        key=lambda item: float(item[1]["target_mz"]),
    )
    target_mzs = np.asarray(
        [float(item[1]["target_mz"]) for item in ordered],
        dtype=float,
    )
    original_indices = [item[0] for item in ordered]
    tolerances = target_mzs * float(mass_ppm) * 1e-6
    hits: list[dict[int, dict[str, Any]]] = [
        {} for _ in ion_hypotheses
    ]
    scan_rts: list[float] = []
    scan_ids: list[str] = []
    ms1_scan_count = 0
    invalid_ms1_scan_count = 0

    for raw_index, spectrum in enumerate(spectra):
        if spectrum.get("ms level") != 1:
            continue
        ms1_scan_count += 1
        rt = _scan_rt_minutes(spectrum)
        mz_array = spectrum.get("m/z array")
        intensity_array = spectrum.get("intensity array")
        if rt is None or mz_array is None or intensity_array is None:
            invalid_ms1_scan_count += 1
            continue
        mzs = np.asarray(mz_array, dtype=float)
        intensities = np.asarray(intensity_array, dtype=float)
        if mzs.size == 0 or mzs.size != intensities.size:
            invalid_ms1_scan_count += 1
            continue
        finite = np.isfinite(mzs) & np.isfinite(intensities)
        mzs = mzs[finite]
        intensities = intensities[finite]
        if mzs.size == 0:
            invalid_ms1_scan_count += 1
            continue
        if np.any(np.diff(mzs) < 0):
            order = np.argsort(mzs, kind="stable")
            mzs = mzs[order]
            intensities = intensities[order]

        usable_index = len(scan_rts)
        scan_rts.append(float(rt))
        scan_ids.append(str(spectrum.get("id") or f"index={raw_index}"))
        if target_mzs.size == 0:
            continue
        starts = np.searchsorted(
            mzs, target_mzs - tolerances, side="left"
        )
        ends = np.searchsorted(
            mzs, target_mzs + tolerances, side="right"
        )
        for sorted_index in np.flatnonzero(ends > starts):
            start = int(starts[sorted_index])
            end = int(ends[sorted_index])
            local = int(np.argmax(intensities[start:end])) + start
            intensity = float(intensities[local])
            if intensity <= 0:
                continue
            original_index = original_indices[int(sorted_index)]
            observed_mz = float(mzs[local])
            isotope = None
            if str(ion_hypotheses[original_index].get("formula") or ""):
                isotope = _isotope_envelope_at_scan(
                    mzs,
                    intensities,
                    observed_mz,
                    intensity,
                    charge=int(
                        ion_hypotheses[original_index].get("charge") or 1
                    ),
                    mass_ppm=mass_ppm,
                )
            hits[original_index][usable_index] = {
                "intensity": intensity,
                "observed_mz": observed_mz,
                "isotope_envelope": isotope,
            }

    return {
        "status": (
            "ok"
            if ms1_scan_count and scan_rts
            else "unavailable_no_ms1_scans"
            if ms1_scan_count == 0
            else "unavailable_no_usable_ms1_scans"
        ),
        "ms1_scan_count": ms1_scan_count,
        "usable_ms1_scan_count": len(scan_rts),
        "invalid_ms1_scan_count": invalid_ms1_scan_count,
        "scan_rts": scan_rts,
        "scan_ids": scan_ids,
        "hits": hits,
    }


def _trace_noise(intensities: np.ndarray) -> tuple[float, float]:
    if not intensities.size:
        return 0.0, 1.0
    # Estimate chromatographic baseline from the quietest part of the XIC.
    # Using the median of the whole trace makes a broad real peak become its
    # own "noise" when it occupies more than half of a short acquisition.
    ordered = np.sort(np.asarray(intensities, dtype=float))
    quiet_count = max(1, int(math.ceil(ordered.size * 0.30)))
    quiet = ordered[:quiet_count]
    baseline = float(np.median(quiet))
    mad = (
        float(np.median(np.abs(quiet - baseline)))
        if quiet.size
        else 0.0
    )
    noise = max(1.0, 1.4826 * mad, baseline * 0.05)
    return baseline, noise


def _active_segments(
    active_indices: np.ndarray,
    *,
    maximum_missing_scans: int = 1,
) -> list[list[int]]:
    if active_indices.size == 0:
        return []
    segments = [[int(active_indices[0])]]
    allowed_difference = max(1, int(maximum_missing_scans) + 1)
    for raw in active_indices[1:]:
        index = int(raw)
        if index - segments[-1][-1] <= allowed_difference:
            segments[-1].append(index)
        else:
            segments.append([index])
    return segments


def _smoothed_trace(values: np.ndarray) -> np.ndarray:
    """Return a light three-point smoothing used only for peak boundaries."""

    array = np.asarray(values, dtype=float)
    if array.size < 3:
        return array.copy()
    padded = np.pad(array, (1, 1), mode="edge")
    return (
        0.25 * padded[:-2]
        + 0.50 * padded[1:-1]
        + 0.25 * padded[2:]
    )


def _split_segment_by_valleys(
    segment: Sequence[int],
    intensities: np.ndarray,
    *,
    minimum_relative_valley_depth: float,
    minimum_scan_points: int,
) -> list[list[int]]:
    """Split a continuous active region when two apices have a real valley."""

    active = sorted({int(index) for index in segment})
    if len(active) < 2 * int(minimum_scan_points):
        return [active]
    start, end = active[0], active[-1]
    trace = _smoothed_trace(intensities[start : end + 1])
    if trace.size < 3:
        return [active]
    maxima = [
        index
        for index in range(1, trace.size - 1)
        if trace[index] >= trace[index - 1]
        and trace[index] > trace[index + 1]
        and trace[index] > 0
    ]
    if len(maxima) < 2:
        return [active]

    split_positions: list[int] = []
    for left_apex, right_apex in zip(maxima, maxima[1:]):
        if right_apex - left_apex < 2:
            continue
        valley_local = left_apex + int(
            np.argmin(trace[left_apex : right_apex + 1])
        )
        shorter = min(float(trace[left_apex]), float(trace[right_apex]))
        if shorter <= 0:
            continue
        depth = (shorter - float(trace[valley_local])) / shorter
        if depth < float(minimum_relative_valley_depth):
            continue
        split_global = start + valley_local
        left_count = sum(index <= split_global for index in active)
        right_count = len(active) - left_count
        if (
            left_count >= int(minimum_scan_points)
            and right_count >= int(minimum_scan_points)
        ):
            split_positions.append(split_global)

    if not split_positions:
        return [active]
    output: list[list[int]] = []
    lower = -math.inf
    for split in sorted(set(split_positions)):
        piece = [index for index in active if lower < index <= split]
        if len(piece) >= int(minimum_scan_points):
            output.append(piece)
            lower = split
    tail = [index for index in active if index > lower]
    if len(tail) >= int(minimum_scan_points):
        output.append(tail)
    return output or [active]


def _open_peak_quality(
    trace: np.ndarray,
    *,
    apex_offset: int,
    gap_ratio: float,
    signal_to_noise: float,
) -> tuple[float, dict[str, float]]:
    """Transparent peak-shape score; it is not Thermo's proprietary PQF."""

    values = np.maximum(0.0, np.asarray(trace, dtype=float))
    if values.size == 0 or float(np.max(values)) <= 0:
        return 0.0, {
            "gap_ratio": round(float(gap_ratio), 4),
            "jaggedness_index": 1.0,
            "zigzag_index": 1.0,
            "apex_balance": 0.0,
        }
    smoothed = _smoothed_trace(values)
    residual = float(np.sum(np.abs(values - smoothed)))
    jaggedness = min(1.0, residual / max(float(np.sum(smoothed)), 1e-12))
    differences = np.diff(smoothed)
    signs = np.sign(differences[np.abs(differences) > 1e-12])
    zigzags = (
        int(np.sum(signs[1:] != signs[:-1]))
        if signs.size > 1
        else 0
    )
    zigzag_index = min(1.0, zigzags / max(int(values.size) - 2, 1))
    left = max(int(apex_offset), 0)
    right = max(int(values.size) - int(apex_offset) - 1, 0)
    apex_balance = (
        min(left, right) / max(left, right)
        if max(left, right) > 0
        else 0.0
    )
    continuity = max(
        0.0,
        1.0 - float(gap_ratio) / max(DEFAULT_MAXIMUM_GAP_RATIO, 1e-9),
    )
    snr_component = min(
        1.0,
        math.log1p(max(float(signal_to_noise), 0.0)) / math.log(11.0),
    )
    score = 10.0 * (
        0.25 * continuity
        + 0.25 * (1.0 - jaggedness)
        + 0.15 * (1.0 - zigzag_index)
        + 0.15 * apex_balance
        + 0.20 * snr_component
    )
    return round(max(0.0, min(10.0, score)), 3), {
        "gap_ratio": round(float(gap_ratio), 4),
        "jaggedness_index": round(jaggedness, 4),
        "zigzag_index": round(zigzag_index, 4),
        "apex_balance": round(apex_balance, 4),
    }


def detect_sparse_xic_peaks(
    ion_hypothesis: Mapping[str, Any],
    scan_rts: Sequence[float],
    scan_ids: Sequence[str],
    hits: Mapping[int, Mapping[str, Any]],
    *,
    minimum_intensity: float = DEFAULT_MINIMUM_INTENSITY,
    signal_to_noise: float = DEFAULT_SIGNAL_TO_NOISE,
    minimum_scan_points: int = DEFAULT_MINIMUM_SCAN_POINTS,
    maximum_missing_scans: int = 1,
    maximum_gap_ratio: float = DEFAULT_MAXIMUM_GAP_RATIO,
    maximum_fwhm_seconds: float = DEFAULT_MAXIMUM_FWHM_SECONDS,
    minimum_relative_valley_depth: float = (
        DEFAULT_MINIMUM_RELATIVE_VALLEY_DEPTH
    ),
) -> list[dict[str, Any]]:
    """Detect chromatographic peaks in one sparse expected-ion XIC."""

    if not scan_rts or not hits:
        return []
    intensities = np.zeros(len(scan_rts), dtype=float)
    for index, item in hits.items():
        if 0 <= int(index) < intensities.size:
            intensities[int(index)] = max(
                0.0, float(item.get("intensity") or 0.0)
            )
    baseline, noise = _trace_noise(intensities)
    snr_values = np.maximum(0.0, intensities - baseline) / noise
    active = np.flatnonzero(
        (intensities >= float(minimum_intensity))
        & (snr_values >= float(signal_to_noise))
    )
    peaks: list[dict[str, Any]] = []
    target_mz = float(ion_hypothesis["target_mz"])
    coarse_segments = _active_segments(
        active, maximum_missing_scans=maximum_missing_scans
    )
    segments = [
        resolved
        for coarse in coarse_segments
        for resolved in _split_segment_by_valleys(
            coarse,
            intensities,
            minimum_relative_valley_depth=minimum_relative_valley_depth,
            minimum_scan_points=minimum_scan_points,
        )
    ]
    for segment in segments:
        if len(segment) < int(minimum_scan_points):
            continue
        apex_index = max(segment, key=lambda index: intensities[index])
        apex_hit = hits.get(apex_index)
        if not apex_hit:
            continue
        boundary_start = segment[0]
        boundary_end = segment[-1]
        segment_slice = slice(boundary_start, boundary_end + 1)
        segment_intensities = intensities[segment_slice]
        segment_rts = np.asarray(
            scan_rts[boundary_start : boundary_end + 1], dtype=float
        )
        boundary_scan_points = boundary_end - boundary_start + 1
        observed_scan_points = int(sum(
            intensities[index] > 0
            for index in range(boundary_start, boundary_end + 1)
        ))
        gap_ratio = (
            1.0 - observed_scan_points / boundary_scan_points
            if boundary_scan_points
            else 1.0
        )
        if gap_ratio > float(maximum_gap_ratio):
            continue
        if segment_rts.size > 1:
            area = float(
                np.trapezoid(segment_intensities, segment_rts * 60.0)
                if hasattr(np, "trapezoid")
                else np.trapz(segment_intensities, segment_rts * 60.0)
            )
            corrected = np.maximum(
                0.0, segment_intensities - float(baseline)
            )
            baseline_corrected_area = float(
                np.trapezoid(corrected, segment_rts * 60.0)
                if hasattr(np, "trapezoid")
                else np.trapz(corrected, segment_rts * 60.0)
            )
        else:
            area = 0.0
            baseline_corrected_area = 0.0
        observed_mz = float(apex_hit["observed_mz"])
        ppm = (observed_mz - target_mz) / target_mz * 1e6
        formula = str(ion_hypothesis.get("formula") or "")
        envelope = apex_hit.get("isotope_envelope")
        isotope_evidence = (
            evaluate_formula_isotope_consistency(formula, envelope)
            if formula and envelope
            else None
        )
        half_height = intensities[apex_index] / 2.0
        half_indices = [
            index for index in range(boundary_start, boundary_end + 1)
            if intensities[index] >= half_height
        ]
        fwhm_seconds = (
            (float(scan_rts[half_indices[-1]]) - float(scan_rts[half_indices[0]]))
            * 60.0
            if len(half_indices) >= 2
            else 0.0
        )
        if (
            fwhm_seconds > 0
            and fwhm_seconds > float(maximum_fwhm_seconds)
        ):
            continue
        peak_quality, peak_quality_metrics = _open_peak_quality(
            np.maximum(0.0, segment_intensities - float(baseline)),
            apex_offset=apex_index - boundary_start,
            gap_ratio=gap_ratio,
            signal_to_noise=float(snr_values[apex_index]),
        )
        peaks.append(
            {
                "expected_ion_id": ion_hypothesis["expected_ion_id"],
                "target_mz": round(target_mz, 8),
                "observed_apex_mz": round(observed_mz, 8),
                "mass_error_ppm": round(ppm, 4),
                "apex_rt_minutes": round(float(scan_rts[apex_index]), 6),
                "apex_scan_id": str(scan_ids[apex_index]),
                "apex_intensity": round(
                    float(intensities[apex_index]), 3
                ),
                "xic_area_intensity_seconds": round(max(0.0, area), 3),
                "baseline_corrected_area_intensity_seconds": round(
                    max(0.0, baseline_corrected_area), 3
                ),
                "boundary_start_rt_minutes": round(
                    float(scan_rts[boundary_start]), 6
                ),
                "boundary_end_rt_minutes": round(
                    float(scan_rts[boundary_end]), 6
                ),
                "fwhm_seconds": round(max(0.0, fwhm_seconds), 3),
                "scan_points": len(segment),
                "boundary_scan_points": boundary_scan_points,
                "xic_nonzero_points": observed_scan_points,
                "gap_ratio": round(gap_ratio, 4),
                "signal_to_noise": round(
                    float(snr_values[apex_index]), 3
                ),
                "peak_quality_score": peak_quality,
                "peak_quality_status": "accepted",
                "peak_quality_metrics": peak_quality_metrics,
                "noise_baseline": round(baseline, 3),
                "noise_scale": round(noise, 3),
                "isotope_envelope": envelope,
                "formula_isotope_evidence": isotope_evidence,
            }
        )
    return peaks


def _row_rt_minutes(row: Mapping[str, Any]) -> float | None:
    value = _number(row.get("retention_time"))
    if value is None:
        return None
    unit = str(row.get("retention_time_unit") or "minutes").lower()
    return value / 60.0 if "second" in unit and "minute" not in unit else value


def _library_match_is_exact_experimental(
    row: Mapping[str, Any],
    member: Mapping[str, Any],
) -> bool:
    if not _library_match_has_exact_experimental_identity(row, member):
        return False
    return bool(_high_accuracy_fragment_support(row)["supporting"])


def _library_match_has_exact_experimental_identity(
    row: Mapping[str, Any],
    member: Mapping[str, Any],
) -> bool:
    match = row.get("match")
    if not isinstance(match, Mapping):
        return False
    if member.get("exact_identity_upgrade_eligible") is False:
        return False
    target_key = normalize_inchikey(member.get("inchikey"))
    library_key = normalize_inchikey(match.get("inchikey"))
    if not target_key or target_key != library_key:
        return False
    if str(match.get("source_tier") or "").strip().lower() not in (
        EXPERIMENTAL_LIBRARY_TIERS
    ):
        return False
    if str(match.get("confidence") or "").strip() != "High":
        return False
    if bool(match.get("isomer_ambiguous")):
        return False
    return bool(match.get("has_structure", True))


def _high_accuracy_fragment_support(
    row: Mapping[str, Any],
    *,
    tolerance_da: float = DEFAULT_HIGH_ACCURACY_FRAGMENT_TOLERANCE_DA,
    minimum_pairs: int = DEFAULT_MINIMUM_HIGH_ACCURACY_FRAGMENT_PAIRS,
    minimum_library_coverage: float = (
        DEFAULT_MINIMUM_HIGH_ACCURACY_LIBRARY_COVERAGE
    ),
) -> dict[str, Any]:
    """Recheck a targeted library upgrade at high-accuracy fragment scale.

    Production spectral similarity may use a broader window for candidate
    recall across heterogeneous public libraries.  That broad score cannot by
    itself upgrade an Expected Compounds target.  This audit uses direct,
    one-to-one fragment matches at the Compound Discoverer Search mzCloud
    high-accuracy default (2.5 mDa); precursor-shifted analogue pairs are not
    accepted as exact product-ion evidence here.
    """

    match = row.get("match")
    query_raw = row.get("query_peaks") or row.get("raw_peaks") or []
    library_raw = (
        (match.get("lib_peaks") or match.get("library_peaks") or [])
        if isinstance(match, Mapping)
        else []
    )

    def usable(peaks: Sequence[Any]) -> list[tuple[float, float]]:
        output: list[tuple[float, float]] = []
        for peak in peaks or []:
            try:
                mz = float(peak[0])
                intensity = max(0.0, float(peak[1]))
            except (TypeError, ValueError, IndexError):
                continue
            if math.isfinite(mz) and math.isfinite(intensity) and mz > 0:
                output.append((mz, intensity))
        return output

    query = usable(query_raw)
    library = usable(library_raw)
    used_query: set[int] = set()
    pairs: list[dict[str, Any]] = []
    for library_mz, library_intensity in sorted(
        library, key=lambda item: (-item[1], item[0])
    ):
        eligible = [
            (abs(query_mz - library_mz), index, query_mz, query_intensity)
            for index, (query_mz, query_intensity) in enumerate(query)
            if index not in used_query
            and abs(query_mz - library_mz) <= float(tolerance_da)
        ]
        if not eligible:
            continue
        error, index, query_mz, query_intensity = min(eligible)
        used_query.add(index)
        pairs.append({
            "query_mz": round(query_mz, 6),
            "library_mz": round(library_mz, 6),
            "mass_error_da": round(query_mz - library_mz, 6),
            "query_intensity": round(query_intensity, 6),
            "library_intensity": round(library_intensity, 6),
        })

    total_library_intensity = sum(item[1] for item in library)
    matched_library_intensity = sum(
        float(item["library_intensity"]) for item in pairs
    )
    coverage = (
        matched_library_intensity / total_library_intensity
        if total_library_intensity > 0
        else 0.0
    )
    supporting = bool(
        len(pairs) >= int(minimum_pairs)
        and coverage >= float(minimum_library_coverage)
    )
    return {
        "status": "supporting" if supporting else "insufficient",
        "supporting": supporting,
        "tolerance_da": float(tolerance_da),
        "matched_pair_count": len(pairs),
        "minimum_pair_count": int(minimum_pairs),
        "query_peak_count": len(query),
        "library_peak_count": len(library),
        "library_intensity_coverage": round(coverage, 4),
        "minimum_library_intensity_coverage": float(
            minimum_library_coverage
        ),
        "pairs": sorted(pairs, key=lambda item: item["library_mz"]),
        "precursor_shifted_pairs_allowed": False,
        "identity_rank_influence": False,
    }


def _library_match_conflicts_target(
    row: Mapping[str, Any],
    member: Mapping[str, Any],
) -> bool:
    """Let stronger measured-library evidence veto a weak diagnostic claim."""

    match = row.get("match")
    if not isinstance(match, Mapping):
        return False
    if not bool(match.get("has_structure", True)):
        return False
    if str(match.get("source_tier") or "").strip().lower() not in (
        EXPERIMENTAL_LIBRARY_TIERS
    ):
        return False
    target_key = normalize_inchikey(member.get("inchikey"))
    library_key = normalize_inchikey(match.get("inchikey"))
    if not target_key or not library_key or target_key == library_key:
        return False
    target_blocks = target_key.split("-")
    library_blocks = library_key.split("-")
    if target_blocks[0] != library_blocks[0]:
        confidence = str(match.get("confidence") or "").strip()
        final_score = _number(match.get("final_score")) or 0.0
        cosine = _number(match.get("cosine")) or 0.0
        return bool(
            confidence in {"High", "Medium"}
            or (final_score >= 0.50 and cosine >= 0.55)
        )

    # Same connectivity with an undefined stereochemical layer does not
    # contradict a stereospecific literature target.  Two different explicit
    # stereochemical identities do conflict when the library call is High.
    target_stereo = target_blocks[1] if len(target_blocks) > 1 else ""
    library_stereo = library_blocks[1] if len(library_blocks) > 1 else ""
    if "UHFFFAOYSA" in {target_stereo, library_stereo}:
        return False
    return str(match.get("confidence") or "").strip() == "High"


def _associate_ms2(
    features: list[dict[str, Any]],
    ion_hypotheses: Mapping[str, Mapping[str, Any]],
    ms2_rows: Sequence[Mapping[str, Any]],
    *,
    precursor_ppm: float = DEFAULT_MS2_PRECURSOR_PPM,
    diagnostic_tolerance_da: float = DEFAULT_DIAGNOSTIC_TOLERANCE_DA,
    minimum_diagnostic_ions: int = DEFAULT_MINIMUM_DIAGNOSTIC_IONS,
    minimum_diagnostic_coverage: float = (
        DEFAULT_MINIMUM_DIAGNOSTIC_COVERAGE
    ),
) -> None:
    indexed: dict[str, list[tuple[float, Mapping[str, Any]]]] = defaultdict(list)
    for row in ms2_rows or []:
        rt = _row_rt_minutes(row)
        polarity = _normalise_polarity(row.get("ion_mode"))
        precursor = _number(
            row.get("search_precursor_mz") or row.get("precursor_mz")
        )
        if rt is None or polarity is None or precursor is None:
            continue
        indexed[polarity].append((rt, row))
    rt_axes: dict[str, list[float]] = {}
    for polarity, rows in indexed.items():
        rows.sort(key=lambda item: item[0])
        rt_axes[polarity] = [item[0] for item in rows]

    for feature in features:
        ion = ion_hypotheses[feature["expected_ion_id"]]
        polarity = str(ion["polarity"])
        candidates = indexed.get(polarity, [])
        axis = rt_axes.get(polarity, [])
        raw_start_rt = float(feature["boundary_start_rt_minutes"])
        raw_end_rt = float(feature["boundary_end_rt_minutes"])
        apex_rt = float(feature["apex_rt_minutes"])
        fwhm_seconds = max(
            0.0, float(feature.get("fwhm_seconds") or 0.0)
        )
        association_half_width = min(
            0.25,
            max(0.05, fwhm_seconds / 120.0 + 0.03),
        )
        start_rt = max(
            raw_start_rt, apex_rt - association_half_width
        )
        end_rt = min(
            raw_end_rt, apex_rt + association_half_width
        )
        if end_rt < start_rt:
            start_rt = apex_rt - 0.05
            end_rt = apex_rt + 0.05
        feature["ms2_association_start_rt_minutes"] = round(
            start_rt, 6
        )
        feature["ms2_association_end_rt_minutes"] = round(end_rt, 6)
        feature["ms2_association_half_width_minutes"] = round(
            association_half_width, 6
        )
        left = bisect.bisect_left(axis, start_rt)
        right = bisect.bisect_right(axis, end_rt)
        target_mz = float(feature["target_mz"])
        associated: list[Mapping[str, Any]] = []
        for _, row in candidates[left:right]:
            precursor = float(
                row.get("search_precursor_mz") or row["precursor_mz"]
            )
            ppm = (precursor - target_mz) / target_mz * 1e6
            if abs(ppm) <= float(precursor_ppm):
                associated.append(row)

        diagnostic_support: list[dict[str, Any]] = []
        diagnostic_conflicts: list[dict[str, Any]] = []
        library_support: list[dict[str, Any]] = []
        library_rejections: list[dict[str, Any]] = []
        for row in associated:
            peaks = row.get("query_peaks") or row.get("raw_peaks") or []
            for member in ion.get("members") or []:
                diagnostic_ions = list(member.get("diagnostic_ions") or [])
                # A reported product-ion list belongs to its isolated adduct.
                # Common-adduct generation reuses candidate membership, but
                # cannot transfer a protonated/negative-ion fragmentation rule
                # to a different precursor species.
                diagnostic_adduct = _normalise_adduct(
                    member.get("active_adduct") or member.get("reported_adduct")
                )
                if (
                    not diagnostic_adduct
                    or diagnostic_adduct != _normalise_adduct(ion.get("adduct"))
                ):
                    diagnostic_ions = []
                required = max(
                    int(minimum_diagnostic_ions),
                    int(
                        math.ceil(
                            len(diagnostic_ions)
                            * float(minimum_diagnostic_coverage)
                        )
                    ),
                )
                matched = _matched_diagnostic_ions(
                    peaks, diagnostic_ions, diagnostic_tolerance_da
                )
                if diagnostic_ions and len(matched) >= required:
                    if _library_match_conflicts_target(row, member):
                        match = row.get("match") or {}
                        diagnostic_conflicts.append(
                            {
                                "target_id": member.get("target_id"),
                                "compound_name": member.get("compound_name"),
                                "spectrum_id": str(
                                    row.get("spectrum_id") or ""
                                ),
                                "conflicting_library_inchikey": str(
                                    match.get("inchikey") or ""
                                ),
                                "conflicting_library_compound_name": str(
                                    match.get("compound_name") or ""
                                ),
                                "conflicting_library_confidence": str(
                                    match.get("confidence") or ""
                                ),
                                "conflicting_library_final_score": (
                                    match.get("final_score")
                                ),
                                "reason": "orthogonal_library_identity_conflict",
                            }
                        )
                        continue
                    diagnostic_support.append(
                        {
                            "target_id": member.get("target_id"),
                            "compound_name": member.get("compound_name"),
                            "spectrum_id": str(row.get("spectrum_id") or ""),
                            "matched_diagnostic_ions": matched,
                            "diagnostic_ion_count": len(diagnostic_ions),
                            "required_diagnostic_ions": required,
                            "diagnostic_coverage": round(
                                len(matched) / len(diagnostic_ions), 4
                            ),
                        }
                    )
                if _library_match_is_exact_experimental(row, member):
                    match = row["match"]
                    high_accuracy = _high_accuracy_fragment_support(row)
                    library_support.append(
                        {
                            "target_id": member.get("target_id"),
                            "compound_name": member.get("compound_name"),
                            "target_inchikey": member.get("inchikey"),
                            "spectrum_id": str(row.get("spectrum_id") or ""),
                            "library_spectrum_id": str(
                                match.get("lib_spectrum_id") or ""
                            ),
                            "library_dataset": str(
                                match.get("dataset_source") or ""
                            ),
                            "library_source_tier": str(
                                match.get("source_tier") or ""
                            ),
                            "library_match_confidence": str(
                                match.get("confidence") or ""
                            ),
                            "cosine": match.get("cosine"),
                            "final_score": match.get("final_score"),
                            "high_accuracy_fragment_support": high_accuracy,
                        }
                    )
                elif _library_match_has_exact_experimental_identity(
                    row, member
                ):
                    high_accuracy = _high_accuracy_fragment_support(row)
                    library_rejections.append({
                        "target_id": member.get("target_id"),
                        "compound_name": member.get("compound_name"),
                        "spectrum_id": str(row.get("spectrum_id") or ""),
                        "reason": (
                            "insufficient_high_accuracy_fragment_support"
                        ),
                        "high_accuracy_fragment_support": high_accuracy,
                    })

        diagnostic_support.sort(
            key=lambda item: (
                -len(item["matched_diagnostic_ions"]),
                -float(item["diagnostic_coverage"]),
                str(item["compound_name"]),
            )
        )
        library_support.sort(
            key=lambda item: (
                -float(item.get("cosine") or 0),
                str(item["compound_name"]),
            )
        )
        unique_library_support: list[dict[str, Any]] = []
        seen_library_support: set[tuple[str, str, str]] = set()
        for item in library_support:
            key = (
                str(item.get("target_inchikey") or ""),
                str(item.get("spectrum_id") or ""),
                str(item.get("library_spectrum_id") or ""),
            )
            if key in seen_library_support:
                continue
            seen_library_support.add(key)
            unique_library_support.append(item)
        feature["associated_ms2"] = [
            {
                "spectrum_id": str(row.get("spectrum_id") or ""),
                "retention_time_minutes": _row_rt_minutes(row),
                "precursor_mz": row.get("precursor_mz"),
                "has_library_match": isinstance(row.get("match"), Mapping),
                "library_compound_name": (
                    str((row.get("match") or {}).get("compound_name") or "")
                    if isinstance(row.get("match"), Mapping)
                    else ""
                ),
            }
            for row in associated
        ]
        feature["ms2_diagnostic_support"] = diagnostic_support
        feature["ms2_diagnostic_conflicts"] = diagnostic_conflicts
        feature["experimental_library_support"] = unique_library_support
        feature["experimental_library_rejections"] = library_rejections


def _assign_evidence_level(feature: Mapping[str, Any]) -> str:
    if feature.get("experimental_library_support"):
        return "EXPERIMENTAL_LIBRARY_CONFIRMED"
    if feature.get("ms2_diagnostic_support"):
        return "MS1_MS2_DIAGNOSTIC_PROBABLE"
    isotope = feature.get("formula_isotope_evidence") or {}
    envelope = feature.get("isotope_envelope") or {}
    if (
        isotope.get("status") == "supporting"
        and envelope.get("pattern_integrity") != "possible_subpattern"
    ):
        return "MS1_ISOTOPE_SUPPORTED"
    return "MS1_EXACT_MASS_CANDIDATE"


def link_expected_features_to_openms(
    features: Sequence[dict[str, Any]],
    openms_features: Sequence[Mapping[str, Any]],
    *,
    mass_ppm: float = DEFAULT_MASS_PPM,
    rt_tolerance_minutes: float = DEFAULT_OPENMS_RT_TOLERANCE_MINUTES,
) -> dict[str, Any]:
    """Cross-check targeted XIC peaks against independent OpenMS features."""

    usable = [
        feature
        for feature in openms_features or []
        if (_number(feature.get("mz")) or 0) > 0
        and _number(feature.get("rt_minutes")) is not None
    ]
    usable.sort(key=lambda item: float(item["mz"]))
    mz_axis = [float(item["mz"]) for item in usable]
    linked = 0
    ambiguous = 0
    for feature in features or []:
        observed = (
            _number(feature.get("observed_apex_mz"))
            or _number(feature.get("target_mz"))
        )
        rt = _number(feature.get("apex_rt_minutes"))
        if observed is None or rt is None or not usable:
            feature["openms_feature_supported"] = False
            feature["openms_feature_link_status"] = (
                "openms_features_unavailable"
                if not usable
                else "target_feature_missing_mz_or_rt"
            )
            continue
        tolerance = max(
            1e-12, observed * float(mass_ppm) * 1e-6
        )
        left = bisect.bisect_left(mz_axis, observed - tolerance)
        right = bisect.bisect_right(mz_axis, observed + tolerance)
        candidates: list[tuple[float, float, str, Mapping[str, Any]]] = []
        target_start = _number(feature.get("boundary_start_rt_minutes"))
        target_end = _number(feature.get("boundary_end_rt_minutes"))
        for candidate in usable[left:right]:
            candidate_rt = float(candidate["rt_minutes"])
            candidate_start_value = _number(
                candidate.get("rt_start_minutes")
            )
            candidate_end_value = _number(
                candidate.get("rt_end_minutes")
            )
            candidate_start = (
                candidate_rt
                if candidate_start_value is None
                else float(candidate_start_value)
            )
            candidate_end = (
                candidate_rt
                if candidate_end_value is None
                else float(candidate_end_value)
            )
            target_overlaps = (
                target_start is not None
                and target_end is not None
                and candidate_start - rt_tolerance_minutes
                <= target_end
                and candidate_end + rt_tolerance_minutes
                >= target_start
            )
            if not target_overlaps and abs(candidate_rt - rt) > float(
                rt_tolerance_minutes
            ):
                continue
            mass_error = abs(float(candidate["mz"]) - observed) / tolerance
            rt_error = abs(candidate_rt - rt) / max(
                float(rt_tolerance_minutes), 1e-9
            )
            candidates.append(
                (
                    mass_error + 0.25 * rt_error,
                    -float(candidate.get("intensity") or 0.0),
                    str(candidate.get("feature_id") or ""),
                    candidate,
                )
            )
        candidates.sort(key=lambda item: item[:3])
        if not candidates:
            feature["openms_feature_supported"] = False
            feature["openms_feature_link_status"] = (
                "no_openms_feature_in_mz_rt_window"
            )
            feature["openms_feature_candidate_count"] = 0
            continue
        selected = candidates[0][3]
        is_ambiguous = (
            len(candidates) > 1
            and abs(candidates[1][0] - candidates[0][0]) <= 0.05
        )
        if is_ambiguous:
            ambiguous += 1
        linked += 1
        feature.update(
            {
                "openms_feature_supported": True,
                "openms_feature_link_status": (
                    "ambiguous_openms_features"
                    if is_ambiguous
                    else "linked_openms_feature"
                ),
                "openms_feature_candidate_count": len(candidates),
                "openms_feature_id": selected.get("feature_id"),
                "openms_feature_mz": selected.get("mz"),
                "openms_feature_rt_minutes": selected.get("rt_minutes"),
                "openms_feature_rt_start_minutes": selected.get(
                    "rt_start_minutes"
                ),
                "openms_feature_rt_end_minutes": selected.get(
                    "rt_end_minutes"
                ),
                "openms_feature_intensity": selected.get("intensity"),
                "openms_feature_quality": selected.get("quality"),
                "openms_isotope_trace_count": selected.get(
                    "isotope_trace_count"
                ),
                "openms_adduct": selected.get("adduct"),
                "openms_neutral_mass": selected.get("neutral_mass"),
                "openms_adduct_group_id": selected.get("adduct_group_id"),
            }
        )
    return {
        "method": "targeted XIC to OpenMS feature cross-check (m/z + RT)",
        "targeted_feature_count": len(features or []),
        "openms_feature_count": len(usable),
        "linked_targeted_features": linked,
        "ambiguous_targeted_features": ambiguous,
        "unlinked_targeted_features": len(features or []) - linked,
        "mass_tolerance_ppm": float(mass_ppm),
        "rt_tolerance_minutes": float(rt_tolerance_minutes),
        "identity_confidence_increment": 0,
    }


def analyse_expected_compound_spectra(
    spectra: Iterable[Mapping[str, Any]],
    ion_hypotheses: Sequence[Mapping[str, Any]],
    *,
    ms2_rows: Sequence[Mapping[str, Any]] | None = None,
    mass_ppm: float = DEFAULT_MASS_PPM,
    minimum_intensity: float = DEFAULT_MINIMUM_INTENSITY,
    signal_to_noise: float = DEFAULT_SIGNAL_TO_NOISE,
    minimum_scan_points: int = DEFAULT_MINIMUM_SCAN_POINTS,
    ms2_precursor_ppm: float = DEFAULT_MS2_PRECURSOR_PPM,
    diagnostic_tolerance_da: float = DEFAULT_DIAGNOSTIC_TOLERANCE_DA,
    maximum_gap_ratio: float = DEFAULT_MAXIMUM_GAP_RATIO,
    maximum_fwhm_seconds: float = DEFAULT_MAXIMUM_FWHM_SECONDS,
    minimum_relative_valley_depth: float = (
        DEFAULT_MINIMUM_RELATIVE_VALLEY_DEPTH
    ),
) -> dict[str, Any]:
    """Run expected-compound detection against an MS spectrum iterable."""

    extracted = extract_sparse_expected_xics(
        spectra, ion_hypotheses, mass_ppm=mass_ppm
    )
    features: list[dict[str, Any]] = []
    if extracted["status"] == "ok":
        for index, ion in enumerate(ion_hypotheses):
            peaks = detect_sparse_xic_peaks(
                ion,
                extracted["scan_rts"],
                extracted["scan_ids"],
                extracted["hits"][index],
                minimum_intensity=minimum_intensity,
                signal_to_noise=signal_to_noise,
                minimum_scan_points=minimum_scan_points,
                maximum_gap_ratio=maximum_gap_ratio,
                maximum_fwhm_seconds=maximum_fwhm_seconds,
                minimum_relative_valley_depth=(
                    minimum_relative_valley_depth
                ),
            )
            for peak_number, peak in enumerate(peaks, start=1):
                feature = {
                    **peak,
                    "feature_id": (
                        f"{ion['expected_ion_id']}_P{peak_number:03d}"
                    ),
                    "formula": ion.get("formula") or "",
                    "adduct": ion.get("adduct") or "",
                    "polarity": ion.get("polarity") or "",
                    "generation_sources": list(
                        ion.get("generation_sources") or []
                    ),
                    "ion_origins": list(ion.get("ion_origins") or []),
                    "ion_quality_flags": list(
                        ion.get("ion_quality_flags") or []
                    ),
                    "reported_precursor_priority": bool(
                        ion.get("reported_precursor_priority")
                    ),
                    "formula_theoretical_primary": bool(
                        ion.get("formula_theoretical_primary")
                    ),
                    "reported_precursor_conflict_candidate": bool(
                        ion.get(
                            "reported_precursor_conflict_candidate"
                        )
                    ),
                    "requires_orthogonal_evidence_for_upgrade": bool(
                        ion.get(
                            "requires_orthogonal_evidence_for_upgrade"
                        )
                    ),
                    "reported_formula_adduct_mz_checks": list(
                        ion.get(
                            "reported_formula_adduct_mz_checks"
                        ) or []
                    ),
                    "candidate_names": list(
                        ion.get("candidate_names") or []
                    ),
                    "nominal_target_count": int(
                        ion.get("nominal_target_count") or 0
                    ),
                    "catalog_row_count": int(
                        ion.get("catalog_row_count") or 0
                    ),
                    "identity_count": int(
                        ion.get("identity_count") or 0
                    ),
                    "isomer_ambiguous": bool(
                        ion.get("isomer_ambiguous")
                    ),
                    "dataset_sources": list(
                        ion.get("dataset_sources") or []
                    ),
                    "target_members": list(ion.get("members") or []),
                    "counted_in_standard_spectral_library_matches": False,
                }
                features.append(feature)
        ion_map = {
            str(item["expected_ion_id"]): item for item in ion_hypotheses
        }
        _associate_ms2(
            features,
            ion_map,
            ms2_rows or [],
            precursor_ppm=ms2_precursor_ppm,
            diagnostic_tolerance_da=diagnostic_tolerance_da,
        )
        for feature in features:
            feature["evidence_level"] = _assign_evidence_level(feature)
            feature["legacy_evidence_level_code"] = feature["evidence_level"]
            feature["evidence_display_label"] = EVIDENCE_DISPLAY_LABELS[
                feature["evidence_level"]
            ]
            feature["standard_confirmed"] = False
            feature["exact_identity_confirmed"] = False
            feature["feature_fragment_calibration_used"] = False
            feature["diagnostic_rule_validation_status"] = (
                "legacy_fixed_threshold_not_v15_conditionally_calibrated"
                if feature["evidence_level"]
                in {
                    "MS1_MS2_DIAGNOSTIC_PROBABLE",
                    "EXPERIMENTAL_LIBRARY_CONFIRMED",
                }
                else "not_applicable"
            )
            feature["evidence_rank"] = (
                EVIDENCE_LEVELS.index(feature["evidence_level"]) + 1
            )
            feature["cd_pattern_accepted"] = (
                feature["evidence_level"] != "MS1_EXACT_MASS_CANDIDATE"
            )
            feature["ordinary_library_match_counted"] = False
            names = list(feature.get("candidate_names") or [])
            feature["candidate_name"] = (
                names[0]
                if len(names) == 1
                else " / ".join(names[:3])
                + (f" (+{len(names) - 3})" if len(names) > 3 else "")
            )
            feature["identification_scope"] = (
                "ambiguity_group"
                if feature["isomer_ambiguous"]
                else "expected_compound_candidate"
            )
            feature["evidence_interpretation"] = {
                "MS1_EXACT_MASS_CANDIDATE": (
                    "Chromatographic MS1 feature at the expected ion mass; "
                    "not an isomer-specific identification."
                ),
                "MS1_ISOTOPE_SUPPORTED": (
                    "Exact-mass feature with an apex isotope envelope that "
                    "supports the reported molecular formula."
                ),
                "MS1_MS2_DIAGNOSTIC_PROBABLE": (
                    "Exact-mass feature with at least two reported diagnostic "
                    "ions under the legacy fixed threshold. The individual "
                    "ions were not conditionally calibrated by the V15 "
                    "feature-fragment module; this is screening support only."
                ),
                "EXPERIMENTAL_LIBRARY_CONFIRMED": (
                    "Associated measured MS2 has a High exact-InChIKey match "
                    "to an experimental/merged/curated/massbank_supplement "
                    "library spectrum. This is an MS/MS library candidate, "
                    "not an authentic-standard confirmation of the sample."
                ),
            }[feature["evidence_level"]]

    evidence_counts = Counter(
        str(feature["evidence_level"]) for feature in features
    )
    detected_ions = {
        str(feature["expected_ion_id"]) for feature in features
    }
    detected_target_ids = {
        str(member.get("target_id") or "")
        for feature in features
        for member in feature.get("target_members") or []
        if str(member.get("target_id") or "")
    }
    detected_names = {
        str(name)
        for feature in features
        for name in feature.get("candidate_names") or []
        if str(name)
    }
    return {
        "status": extracted["status"],
        "ms1_scan_count": int(extracted["ms1_scan_count"]),
        "usable_ms1_scan_count": int(extracted["usable_ms1_scan_count"]),
        "invalid_ms1_scan_count": int(
            extracted["invalid_ms1_scan_count"]
        ),
        "features": sorted(
            features,
            key=lambda item: (
                float(item["apex_rt_minutes"]),
                float(item["target_mz"]),
                str(item["feature_id"]),
            ),
        ),
        "summary": {
            "status": extracted["status"],
            "expected_ion_hypotheses": len(ion_hypotheses),
            "detected_feature_count": len(features),
            "detected_ion_hypotheses": len(detected_ions),
            "candidate_target_rows": len(detected_target_ids),
            "candidate_target_names": len(detected_names),
            "ambiguous_feature_count": sum(
                bool(feature.get("isomer_ambiguous"))
                for feature in features
            ),
            "supported_features": sum(
                bool(feature.get("cd_pattern_accepted"))
                for feature in features
            ),
            "mass_only_features": sum(
                not bool(feature.get("cd_pattern_accepted"))
                for feature in features
            ),
            "evidence_counts": {
                level: int(evidence_counts.get(level, 0))
                for level in EVIDENCE_LEVELS
            },
            "ms1_scan_count": int(extracted["ms1_scan_count"]),
            "usable_ms1_scan_count": int(
                extracted["usable_ms1_scan_count"]
            ),
            "counted_in_standard_spectral_library_matches": False,
            "standard_spectral_library_match_increment": 0,
        },
    }


def unavailable_payload(
    *,
    status: str,
    input_path: str | Path,
    requested: bool = True,
    warning: str,
    catalog: Mapping[str, Any] | None = None,
    summary_audit: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    audit_fields = dict(summary_audit or {})
    return {
        "schema_version": 1,
        "method": "Expected Compounds / continuous-MS1 targeted XIC",
        "method_reference": COMPOUND_DISCOVERER_METHOD_REFERENCE,
        "status": status,
        "requested": bool(requested),
        "input": str(input_path),
        "summary": {
            "status": status,
            "expected_ion_hypotheses": 0,
            "detected_feature_count": 0,
            "detected_ion_hypotheses": 0,
            "candidate_target_rows": 0,
            "candidate_target_names": 0,
            "ambiguous_feature_count": 0,
            "supported_features": 0,
            "mass_only_features": 0,
            "evidence_counts": {level: 0 for level in EVIDENCE_LEVELS},
            "ms1_scan_count": 0,
            "usable_ms1_scan_count": 0,
            "counted_in_standard_spectral_library_matches": False,
            "standard_spectral_library_match_increment": 0,
            "warning": warning,
            **audit_fields,
        },
        "parameters": {},
        "catalog": dict(catalog or {}),
        "features": [],
        "scientific_boundaries": {
            "presence_only_rows_are_not_reference_spectra": True,
            "exact_mass_is_not_isomer_specific": True,
            "ordinary_library_matched_count_unchanged": True,
        },
    }


def run_expected_compounds_workflow(
    file_path: str | Path,
    *,
    db_path: str | Path,
    polarity: str,
    ms2_rows: Sequence[Mapping[str, Any]] | None = None,
    dataset_sources: Sequence[str] | None = None,
    include_safe_common_adducts: bool = True,
    include_reported_precursor_conflicts: bool = False,
    mass_ppm: float = DEFAULT_MASS_PPM,
    minimum_intensity: float = DEFAULT_MINIMUM_INTENSITY,
    signal_to_noise: float = DEFAULT_SIGNAL_TO_NOISE,
    minimum_scan_points: int = DEFAULT_MINIMUM_SCAN_POINTS,
    ms2_precursor_ppm: float = DEFAULT_MS2_PRECURSOR_PPM,
    diagnostic_tolerance_da: float = DEFAULT_DIAGNOSTIC_TOLERANCE_DA,
    maximum_gap_ratio: float = DEFAULT_MAXIMUM_GAP_RATIO,
    maximum_fwhm_seconds: float = DEFAULT_MAXIMUM_FWHM_SECONDS,
    minimum_relative_valley_depth: float = (
        DEFAULT_MINIMUM_RELATIVE_VALLEY_DEPTH
    ),
    openms_features: Sequence[Mapping[str, Any]] | None = None,
    openms_rt_tolerance_minutes: float = (
        DEFAULT_OPENMS_RT_TOLERANCE_MINUTES
    ),
) -> dict[str, Any]:
    """Run the complete read-only expected-compounds workflow for one mzML."""

    path = Path(file_path).expanduser().resolve()
    if path.suffix.lower() != ".mzml":
        return unavailable_payload(
            status="unavailable_requires_mzml_with_ms1",
            input_path=path,
            warning=(
                "Expected Compounds requires continuous survey MS1 scans in "
                "mzML. MGF/HTML peak-list inputs do not contain the required "
                "chromatographic MS1 trace."
            ),
        )
    mode = _normalise_polarity(polarity)
    if mode not in {"POS", "NEG"}:
        return unavailable_payload(
            status="unavailable_unknown_or_mixed_polarity",
            input_path=path,
            warning=(
                "Expected-ion generation requires one resolved POS or NEG "
                "acquisition polarity."
            ),
        )
    selected_sources = list(dict.fromkeys(
        str(item).strip()
        for item in (dataset_sources or [])
        if str(item).strip()
    ))
    if not selected_sources:
        return unavailable_payload(
            status="unavailable_no_selected_target_dataset",
            input_path=path,
            warning=(
                "No sample-specific Expected Compounds catalogue was "
                "selected or inferred. The target layer was not run; an "
                "empty selection never searches all catalogues."
            ),
        )
    targets = load_expected_compound_catalog(
        db_path, dataset_sources=selected_sources
    )
    if not targets:
        return unavailable_payload(
            status="unavailable_no_presence_only_targets",
            input_path=path,
            warning=(
                "The configured TargetCompounds catalog has no selected "
                "literature_presence_only rows."
            ),
        )
    polarity_targets = [
        target
        for target in targets
        if str(target.get("target_polarity_scope") or "").strip().upper()
        in {"", mode}
    ]
    hypotheses = build_expected_ion_hypotheses(
        polarity_targets,
        mode,
        include_safe_common_adducts=include_safe_common_adducts,
        include_reported_precursor_conflicts=(
            include_reported_precursor_conflicts
        ),
        ambiguity_ppm=mass_ppm,
    )
    catalog_quality = target_catalog_quality_summary(polarity_targets)
    if not hypotheses:
        return unavailable_payload(
            status="unavailable_no_valid_expected_ions",
            input_path=path,
            warning=(
                "No selected target row had a usable reported precursor or "
                "formula/exact-mass ion hypothesis for this polarity."
            ),
            catalog={
                "source_tier": "literature_presence_only",
                "target_rows": len(polarity_targets),
                "selected_dataset_sources": selected_sources,
                "expected_ion_hypotheses": 0,
                "quality_audit": catalog_quality,
            },
            summary_audit={
                "target_catalog_rows": len(polarity_targets),
                "active_target_rows": catalog_quality[
                    "active_target_rows"
                ],
                "quarantined_target_rows": catalog_quality[
                    "quarantined_target_rows"
                ],
                "target_catalog_quality_audit": catalog_quality,
                "selected_dataset_sources": selected_sources,
            },
        )
    from pyteomics import mzml

    with mzml.read(str(path)) as reader:
        analysed = analyse_expected_compound_spectra(
            reader,
            hypotheses,
            ms2_rows=ms2_rows,
            mass_ppm=mass_ppm,
            minimum_intensity=minimum_intensity,
            signal_to_noise=signal_to_noise,
            minimum_scan_points=minimum_scan_points,
            ms2_precursor_ppm=ms2_precursor_ppm,
            diagnostic_tolerance_da=diagnostic_tolerance_da,
            maximum_gap_ratio=maximum_gap_ratio,
            maximum_fwhm_seconds=maximum_fwhm_seconds,
            minimum_relative_valley_depth=(
                minimum_relative_valley_depth
            ),
        )
    openms_validation = (
        link_expected_features_to_openms(
            analysed["features"],
            openms_features,
            mass_ppm=mass_ppm,
            rt_tolerance_minutes=openms_rt_tolerance_minutes,
        )
        if openms_features is not None
        else {
            "method": "targeted XIC to OpenMS feature cross-check (m/z + RT)",
            "status": "not_available",
            "targeted_feature_count": len(analysed["features"]),
            "openms_feature_count": 0,
            "linked_targeted_features": 0,
            "ambiguous_targeted_features": 0,
            "unlinked_targeted_features": len(analysed["features"]),
            "identity_confidence_increment": 0,
        }
    )
    analysed["summary"]["openms_feature_validation"] = openms_validation
    payload = {
        "schema_version": 1,
        "method": "Expected Compounds / continuous-MS1 targeted XIC",
        "method_reference": COMPOUND_DISCOVERER_METHOD_REFERENCE,
        "status": analysed["status"],
        "requested": True,
        "input": str(path),
        "target_catalog_db": str(Path(db_path).expanduser().resolve()),
        "parameters": {
            "mass_tolerance_ppm": float(mass_ppm),
            "minimum_intensity": float(minimum_intensity),
            "signal_to_noise": float(signal_to_noise),
            "minimum_scan_points": int(minimum_scan_points),
            "maximum_gap_ratio": float(maximum_gap_ratio),
            "maximum_fwhm_seconds": float(maximum_fwhm_seconds),
            "minimum_relative_valley_depth": float(
                minimum_relative_valley_depth
            ),
            "include_safe_common_adducts": bool(
                include_safe_common_adducts
            ),
            "include_reported_precursor_conflicts": bool(
                include_reported_precursor_conflicts
            ),
            "safe_common_adducts": list(SAFE_COMMON_ADDUCTS[mode]),
            "ms2_precursor_tolerance_ppm": float(ms2_precursor_ppm),
            "diagnostic_ion_tolerance_da": float(
                diagnostic_tolerance_da
            ),
            "high_accuracy_fragment_tolerance_da": (
                DEFAULT_HIGH_ACCURACY_FRAGMENT_TOLERANCE_DA
            ),
            "minimum_high_accuracy_fragment_pairs": (
                DEFAULT_MINIMUM_HIGH_ACCURACY_FRAGMENT_PAIRS
            ),
            "minimum_high_accuracy_library_intensity_coverage": (
                DEFAULT_MINIMUM_HIGH_ACCURACY_LIBRARY_COVERAGE
            ),
            "minimum_diagnostic_ions": (
                DEFAULT_MINIMUM_DIAGNOSTIC_IONS
            ),
            "minimum_diagnostic_coverage": (
                DEFAULT_MINIMUM_DIAGNOSTIC_COVERAGE
            ),
            "openms_rt_tolerance_minutes": float(
                openms_rt_tolerance_minutes
            ),
        },
        "catalog": {
            "source_tier": "literature_presence_only",
            "target_rows": len(polarity_targets),
            "dataset_sources": sorted({
                str(item.get("dataset_source") or "")
                for item in polarity_targets
                if str(item.get("dataset_source") or "")
            }),
            "selected_dataset_sources": sorted({
                str(item) for item in selected_sources if str(item)
            }),
            "expected_ion_hypotheses": len(hypotheses),
            "quality_audit": catalog_quality,
        },
        "summary": {
            **analysed["summary"],
            "target_catalog_rows": len(polarity_targets),
            "active_target_rows": catalog_quality["active_target_rows"],
            "quarantined_target_rows": catalog_quality[
                "quarantined_target_rows"
            ],
            "target_catalog_quality_audit": catalog_quality,
            "selected_dataset_sources": sorted({
                str(item) for item in selected_sources if str(item)
            }),
        },
        "features": analysed["features"],
        "scientific_boundaries": {
            "presence_only_rows_are_not_reference_spectra": True,
            "exact_mass_is_not_isomer_specific": True,
            "same_formula_adduct_mz_names_are_ambiguity_groups": True,
            "experimental_library_confirmed_requires": (
                "exact full InChIKey + High match + "
                "experimental/merged/curated/massbank_supplement measured "
                "library spectrum + at least three direct one-to-one fragment "
                "pairs within 0.0025 Da covering at least 50% of reference "
                "spectrum intensity"
            ),
            "experimental_library_candidate_support_is_standard_confirmation": (
                False
            ),
            "legacy_evidence_level_codes_retained_for_api_compatibility": True,
            "legacy_diagnostic_threshold_is_v15_fragment_calibrated": False,
            "ordinary_library_matched_count_unchanged": True,
            "ordinary_library_unique_compound_count_unchanged": True,
            "openms_feature_support_is_chromatographic_not_identity_evidence": (
                True
            ),
            "implementation_is_open_and_not_thermo_compound_discoverer": True,
        },
    }
    if analysed["status"] != "ok":
        payload["summary"]["warning"] = (
            "The mzML contains no usable survey MS1 scans, so Expected "
            "Compounds was not executed."
        )
    return payload


def write_expected_compounds_artifact(
    path: str | Path,
    payload: Mapping[str, Any],
) -> None:
    """Atomically persist one expected-compounds artifact."""

    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary_name = tempfile.mkstemp(
        prefix=f".{destination.name}.",
        suffix=".tmp",
        dir=str(destination.parent),
    )
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            json.dump(
                payload,
                handle,
                ensure_ascii=False,
                # These machine-consumed audit payloads can contain tens of
                # thousands of nested chromatographic features.  Pretty
                # printing adds hundreds of megabytes per QE task without
                # improving the API/export contract.
                separators=(",", ":"),
                default=str,
            )
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary_name, destination)
    except Exception:
        try:
            os.unlink(temporary_name)
        except OSError:
            pass
        raise
