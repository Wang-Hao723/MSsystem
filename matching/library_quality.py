"""Quality gates for user-contributed and literature-derived library records.

Literature target tables often contain only a short list of reported fragment
m/z values.  Assigning an arbitrary intensity to those ions does not turn the
list into an experimental spectrum.  This module validates chemical metadata,
separates measured spectra from presence-only target leads, and quarantines
records whose precursor identity is internally inconsistent.
"""

from __future__ import annotations

import math
import re
from collections import Counter, defaultdict
from copy import deepcopy
from functools import lru_cache
from typing import Any, Mapping

from utils.adduct import adduct_polarity, neutral_mass


TRUSTED_SOURCE_TIERS = {
    "experimental", "merged", "curated", "massbank_supplement",
}
_RECOGNIZED_DECLARED_SOURCE_TIERS = TRUSTED_SOURCE_TIERS | {
    "predicted", "in_silico", "unknown", "user_unverified",
    "metadata_conflict", "quarantined",
}
# The legacy Small database predates provenance storage.  It remains usable as
# an explicitly unresolved compatibility tier, with a confidence cap, until it
# is rebuilt from raw sources.  New imports may never assign this tier.
SEARCHABLE_SOURCE_TIERS = TRUSTED_SOURCE_TIERS | {"legacy_unresolved"}
PRESENCE_ONLY_SOURCE_TIER = "literature_presence_only"
QUARANTINED_SOURCE_TIER = "quarantined"

# Monoisotopic masses for elements commonly encountered in small-molecule
# libraries.  Unknown elements make the formula unverifiable rather than being
# silently treated as zero mass.
MONOISOTOPIC_MASS = {
    "H": 1.00782503223,
    "B": 11.00930536,
    "C": 12.0,
    "N": 14.00307400443,
    "O": 15.99491461957,
    "F": 18.99840316273,
    "Na": 22.9897692820,
    "Mg": 23.985041697,
    "P": 30.97376199842,
    "S": 31.9720711744,
    "Cl": 34.968852682,
    "K": 38.9637064864,
    "Ca": 39.962590863,
    "Fe": 55.93493633,
    "Br": 78.9183376,
    "I": 126.904468,
    "Si": 27.97692653465,
    "Se": 79.9165218,
}

_FORMULA_TOKEN = re.compile(r"([A-Z][a-z]?)(\d*)")
_FORMULA_CARET_CHARGE = re.compile(r"^(.*)\^(\d*)([+-])$")
_FORMULA_TRAILING_CHARGE = re.compile(r"^(.*)([+-])(\d*)$")
_FULL_INCHIKEY = re.compile(r"^[A-Z]{14}-[A-Z]{10}-[A-Z]$")
_MISSING = {"", "nan", "none", "null", "na", "n/a", "-"}
_SUBSCRIPT_TRANSLATION = str.maketrans("₀₁₂₃₄₅₆₇₈₉", "0123456789")
_CHARGE_TRANSLATION = str.maketrans({"−": "-", "–": "-", "＋": "+"})


def _clean(value: Any) -> str:
    text = "" if value is None else str(value).strip()
    return "" if text.lower() in _MISSING else text


def normalize_inchikey(value: Any) -> str:
    """Return one canonical full InChIKey, or ``""`` when it is invalid.

    Connectivity-only 14-character blocks and malformed public-library labels
    are useful provenance, but they are not stable compound identities.  The
    search/runtime layer therefore only groups spectra when all three standard
    InChIKey blocks are present.
    """
    text = _clean(value).upper()
    if text.startswith("INCHIKEY="):
        text = text.split("=", 1)[1]
    text = re.sub(r"\s+", "", text)
    return text if _FULL_INCHIKEY.fullmatch(text) else ""


def _finite_float(value: Any) -> float | None:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def normalize_formula_text(formula: Any) -> str:
    """Return a simple Hill-order formula, or ``""`` when it cannot be parsed."""
    text = _clean(formula)
    if not text:
        return ""
    # Keep this implementation local.  ``library_import`` imports this module
    # for its quality gate, so importing it back here creates a fragile circular
    # dependency during CLI/database migrations.
    text = text.translate(_SUBSCRIPT_TRANSLATION)
    text = re.split(r"[·.]", text, maxsplit=1)[0]
    text = re.sub(r"\s+", "", text)
    counts: dict[str, int] = {}
    pos = 0
    for match in _FORMULA_TOKEN.finditer(text):
        if match.start() != pos:
            return ""
        element, raw_count = match.groups()
        counts[element] = counts.get(element, 0) + (int(raw_count) if raw_count else 1)
        pos = match.end()
    if pos != len(text) or not counts:
        return ""
    ordered: list[tuple[str, int]] = []
    if "C" in counts:
        ordered.append(("C", counts.pop("C")))
        if "H" in counts:
            ordered.append(("H", counts.pop("H")))
    ordered.extend((element, counts[element]) for element in sorted(counts))
    return "".join(
        f"{element}{count if count != 1 else ''}" for element, count in ordered
    )


def _formula_counts(formula: Any) -> dict[str, int] | None:
    """Parse one uncharged elemental formula into positive integer counts."""
    text = normalize_formula_text(formula)
    if not text:
        return None
    counts: dict[str, int] = {}
    pos = 0
    for match in _FORMULA_TOKEN.finditer(text):
        if match.start() != pos:
            return None
        element, raw_count = match.groups()
        counts[element] = counts.get(element, 0) + (
            int(raw_count) if raw_count else 1
        )
        pos = match.end()
    return counts if pos == len(text) and counts else None


def _counts_formula(counts: Mapping[str, int]) -> str:
    """Render positive elemental counts in Hill order."""
    positive = {
        str(element): int(count)
        for element, count in counts.items()
        if int(count) > 0
    }
    if not positive:
        return ""
    ordered: list[tuple[str, int]] = []
    if "C" in positive:
        ordered.append(("C", positive.pop("C")))
        if "H" in positive:
            ordered.append(("H", positive.pop("H")))
    ordered.extend((element, positive[element]) for element in sorted(positive))
    return "".join(
        f"{element}{count if count != 1 else ''}" for element, count in ordered
    )


def _charged_formula_parts(formula: Any) -> tuple[str, str, int] | None:
    """Return ``(ion_formula, sign, magnitude)`` for an explicitly charged formula.

    Library exports commonly store a deprotonated structure as
    ``C21H19O12-`` while another source stores the neutral parent
    ``C21H20O12``.  The terminal sign is not part of an elemental formula and
    was previously discarded as an unparsable value, silently disabling
    same-neutral-formula isomer checks.
    """
    text = _clean(formula)
    if not text:
        return None
    text = text.translate(_SUBSCRIPT_TRANSLATION).translate(_CHARGE_TRANSLATION)
    text = re.split(r"[·.]", text, maxsplit=1)[0]
    text = re.sub(r"\s+", "", text)
    match = _FORMULA_CARET_CHARGE.fullmatch(text)
    if match:
        body, raw_magnitude, sign = match.groups()
    else:
        match = _FORMULA_TRAILING_CHARGE.fullmatch(text)
        if not match:
            return None
        body, sign, raw_magnitude = match.groups()
    magnitude = int(raw_magnitude) if raw_magnitude else 1
    if magnitude < 1 or not normalize_formula_text(body):
        return None
    return body, sign, magnitude


_ION_TO_NEUTRAL_DELTAS: dict[str, dict[str, int]] = {
    # delta is added to the explicitly charged ion composition to obtain M.
    "[M]+": {},
    "M+": {},
    "[M]-": {},
    "M-": {},
    "[M+H]+": {"H": -1},
    "M+H": {"H": -1},
    "[M-H]-": {"H": 1},
    "M-H": {"H": 1},
    "[M+NA]+": {"Na": -1},
    "M+NA": {"Na": -1},
    "[M+K]+": {"K": -1},
    "M+K": {"K": -1},
    "[M+NH4]+": {"N": -1, "H": -4},
    "M+NH4": {"N": -1, "H": -4},
    "[M+CL]-": {"Cl": -1},
    "M+CL": {"Cl": -1},
    "[M+FA-H]-": {"C": -1, "H": -1, "O": -2},
    "[M+HCOO]-": {"C": -1, "H": -1, "O": -2},
    "[M+HAC-H]-": {"C": -2, "H": -3, "O": -2},
    "[M+CH3COO]-": {"C": -2, "H": -3, "O": -2},
    "[M+AC-H]-": {"C": -2, "H": -3, "O": -2},
    "[M+H-H2O]+": {"H": 1, "O": 1},
    "[M-H2O+H]+": {"H": 1, "O": 1},
    "[M+H-2H2O]+": {"H": 3, "O": 2},
    "[M+H-NH3]+": {"N": 1, "H": 2},
    "[M-H-H2O]-": {"H": 3, "O": 1},
    "[M-H-CO2]-": {"C": 1, "H": 1, "O": 2},
}


def neutral_formula_text(formula: Any, adduct: Any = "") -> str:
    """Return the neutral-parent formula used for *comparison* and exact mass.

    Uncharged library formulae are already conventionally molecular formulae
    and are returned unchanged.  An explicitly charged ion formula is reversed
    only when a chemically specific, supported adduct is present; unknown or
    contradictory annotations fail closed instead of guessing a neutral
    composition.
    """
    neutral = normalize_formula_text(formula)
    if neutral:
        return neutral
    charged = _charged_formula_parts(formula)
    if charged is None:
        return ""
    ion_formula, sign, magnitude = charged
    if magnitude != 1:
        return ""
    normalized_adduct = re.sub(r"\s+", "", _clean(adduct)).upper()
    delta = _ION_TO_NEUTRAL_DELTAS.get(normalized_adduct)
    if delta is None:
        return ""
    expected_sign = (
        "+" if normalized_adduct.endswith("+")
        else "-" if normalized_adduct.endswith("-")
        else ""
    )
    if expected_sign and expected_sign != sign:
        return ""
    counts = _formula_counts(ion_formula)
    if counts is None:
        return ""
    for element, change in delta.items():
        counts[element] = counts.get(element, 0) + change
        if counts[element] < 0:
            return ""
    return _counts_formula(counts)


def compound_identity_key(item: dict[str, Any] | None) -> str:
    """Build a deterministic identity key without collapsing stereoisomers.

    Full InChIKey is authoritative.  Records without a valid key fall back to
    an exact structure string, then normalized name+formula, and finally a
    spectrum identifier.  A record that *has* a valid full InChIKey never falls
    back to its name, so two stereochemical keys cannot be merged merely
    because the source labels happen to match.
    """
    if not isinstance(item, dict):
        return ""
    key = normalize_inchikey(item.get("inchikey"))
    if key:
        return f"ik::{key}"
    smiles = _clean(item.get("smiles"))
    if smiles:
        return f"sm::{smiles}"
    name = re.sub(
        r"[^a-z0-9]+", "", _clean(item.get("compound_name")).lower()
    )
    formula = normalize_formula_text(item.get("formula") or item.get("molecular_formula"))
    if name or formula:
        return f"nf::{name}|{formula}"
    spectrum_id = _clean(
        item.get("lib_spectrum_id") or item.get("spectrum_id")
    )
    return f"sid::{spectrum_id}" if spectrum_id else ""


def formula_exact_mass(formula: Any) -> float | None:
    text = normalize_formula_text(formula)
    if not text:
        return None
    pos = 0
    total = 0.0
    for match in _FORMULA_TOKEN.finditer(text):
        if match.start() != pos:
            return None
        element, raw_count = match.groups()
        mass = MONOISOTOPIC_MASS.get(element)
        if mass is None:
            return None
        count = int(raw_count) if raw_count else 1
        total += mass * count
        pos = match.end()
    return total if pos == len(text) and total > 0 else None


def smiles_metadata(smiles: Any) -> dict[str, Any]:
    text = _clean(smiles)
    if not text:
        return {"status": "missing"}
    try:
        from rdkit import Chem, rdBase
        from rdkit.Chem import Descriptors, rdMolDescriptors

        # Invalid user-supplied structures are expected quality findings, not
        # server errors.  Capture them in the audit without flooding logs with
        # hundreds of RDKit parser diagnostics during a workbook preview.
        with rdBase.BlockLogs():
            mol = Chem.MolFromSmiles(text)
        if mol is None:
            return {"status": "invalid"}
        return {
            "status": "valid",
            "canonical_smiles": Chem.MolToSmiles(mol, isomericSmiles=True),
            "formula": rdMolDescriptors.CalcMolFormula(mol),
            "exact_mass": float(Descriptors.ExactMolWt(mol)),
            "inchikey": Chem.MolToInchiKey(mol),
            "inchi": Chem.MolToInchi(mol),
        }
    except Exception:
        return {"status": "unverifiable"}


@lru_cache(maxsize=200_000)
def resolve_structure_formula(
    smiles: str = "", formula: str = "", adduct: str = ""
) -> str:
    """Return a normalized stated formula or a structure-derived fallback.

    Several public-library rows have a valid SMILES but an empty molecular
    formula.  Leaving those rows empty in Top-K silently disables same-formula
    isomer checks.  Derivation is deliberately lazy and cached: it runs only
    for candidates that survive spectral scoring/presentation, not for all
    millions of library spectra.

    A non-empty stated formula is preserved as a neutral composition so the
    independent structure validator can still detect and quarantine a
    SMILES/formula conflict.
    """
    stated = neutral_formula_text(formula, adduct)
    if stated:
        return stated
    if _charged_formula_parts(formula) is not None:
        # An explicitly charged formula without a supported adduct cannot be
        # replaced by a structure-derived guess without hiding the conflict.
        return ""
    metadata = smiles_metadata(smiles)
    if metadata.get("status") != "valid":
        return ""
    return neutral_formula_text(metadata.get("formula"), adduct)


def inchi_metadata(inchi: Any) -> dict[str, Any]:
    """Parse an InChI as structure evidence when a SMILES is unavailable."""
    text = _clean(inchi)
    if not text:
        return {"status": "missing"}
    if not text.startswith("InChI="):
        return {"status": "invalid"}
    try:
        from rdkit import Chem, rdBase
        from rdkit.Chem import Descriptors, rdMolDescriptors

        with rdBase.BlockLogs():
            mol = Chem.MolFromInchi(text)
        if mol is None:
            return {"status": "invalid"}
        return {
            "status": "valid",
            "canonical_smiles": Chem.MolToSmiles(mol, isomericSmiles=True),
            "formula": rdMolDescriptors.CalcMolFormula(mol),
            "exact_mass": float(Descriptors.ExactMolWt(mol)),
            "inchikey": Chem.MolToInchiKey(mol),
            "inchi": Chem.MolToInchi(mol),
        }
    except Exception:
        return {"status": "unverifiable"}


@lru_cache(maxsize=200_000)
def validate_structure_metadata(
    smiles: str, formula: str = "", adduct: str = ""
) -> tuple[bool, str]:
    """Validate only candidates that survive spectral scoring (cached by metadata)."""
    metadata = smiles_metadata(smiles)
    if metadata.get("status") != "valid":
        return False, str(metadata.get("status") or "invalid")
    stated = neutral_formula_text(formula, adduct)
    derived = neutral_formula_text(metadata.get("formula"), adduct)
    if _charged_formula_parts(formula) is not None and not stated:
        return False, "charged_formula_neutralization_unresolved"
    if _charged_formula_parts(metadata.get("formula")) is not None and not derived:
        return False, "charged_structure_neutralization_unresolved"
    if stated and derived and stated != derived:
        return False, "smiles_formula_mismatch"
    return True, "validated"


def source_tier_is_searchable(source_tier: Any) -> bool:
    """Ordinary Top-K search accepts measured/curated spectra only.

    Predicted spectra and literature diagnostic-ion lists have separate uses;
    allowing them to compete with experimental references silently changes the
    meaning of a library hit even if the UI later caps its confidence.
    """
    return _clean(source_tier).lower() in SEARCHABLE_SOURCE_TIERS


def source_tier_confidence_cap(source_tier: Any) -> str | None:
    """Maximum compound-name confidence allowed for non-experimental evidence."""
    tier = _clean(source_tier).lower()
    if tier == PRESENCE_ONLY_SOURCE_TIER:
        return "Low"
    if tier in {
        "predicted", "in_silico", "unknown", "user_unverified",
        "legacy_unresolved",
    }:
        return "Medium"
    if tier == "metadata_conflict":
        return "Unknown"
    if tier == QUARANTINED_SOURCE_TIER:
        return "Unknown"
    return None


def _extract_peaks(record: dict[str, Any]) -> list[list[float]]:
    peaks = record.get("peaks") or []
    output: list[list[float]] = []
    for peak in peaks:
        try:
            mz = float(peak[0])
            intensity = float(peak[1])
        except (TypeError, ValueError, IndexError):
            continue
        if math.isfinite(mz) and math.isfinite(intensity) and mz > 0 and intensity > 0:
            output.append([mz, intensity])
    if output:
        return output

    fragments = record.get("fragments") or []
    intensities = record.get("fragment_intensities") or []
    origin = _clean(record.get("intensity_origin")).lower()
    for index, value in enumerate(fragments):
        mz = _finite_float(value)
        if mz is None or mz <= 0:
            continue
        intensity = _finite_float(intensities[index]) if index < len(intensities) else None
        if intensity is None or intensity <= 0:
            intensity = 1000.0
            if not origin:
                origin = "default"
        output.append([mz, intensity])
    if origin:
        record["intensity_origin"] = origin
    return output


def _repair_adduct_from_formula_precursor(
    precursor: float,
    exact_mass: float,
    declared_adduct: Any,
    max_precursor_ppm: float,
) -> tuple[str, float] | None:
    """Repair a clearly mistyped proton polarity using internal mass balance.

    The literature workbook contains rows whose m/z precisely fits [M+H]+ but
    is labelled [M-H]- (and vice versa).  Only proton polarity and a
    same-polarity permanent molecular ion are considered; broader adduct
    guessing could accidentally compensate for an incorrect formula.
    """
    polarity = adduct_polarity(declared_adduct)
    candidates = ["[M+H]+", "[M-H]-"]
    if polarity == "POS":
        candidates.append("[M]+")
    elif polarity == "NEG":
        candidates.append("[M]-")
    ranked: list[tuple[float, str, float]] = []
    for candidate in candidates:
        candidate_mass = neutral_mass(precursor, candidate)
        if candidate_mass is None:
            continue
        ppm = (candidate_mass - exact_mass) / exact_mass * 1e6
        ranked.append((abs(ppm), candidate, ppm))
    if not ranked:
        return None
    _, candidate, ppm = min(ranked)
    if abs(ppm) <= max_precursor_ppm and candidate != _clean(declared_adduct):
        return candidate, ppm
    return None


def assess_record(
    record: dict[str, Any],
    max_precursor_ppm: float = 50.0,
    *,
    require_experimental_provenance: bool = False,
) -> tuple[dict[str, Any], dict[str, Any]]:
    """Sanitize one record and return ``(record, audit)``.

    Invalid or formula-inconsistent SMILES are stripped so they cannot generate
    an InChIKey that merges a different named compound.  A precursor/formula
    conflict is a hard quarantine because the row's molecular identity is then
    ambiguous even if its four reported fragments look plausible.
    """
    clean = deepcopy(record)
    errors: list[str] = []
    warnings: list[str] = []

    # An InChIKey is a hash/identifier, not a structure that can be validated
    # on its own.  Keep the submitted value for provenance, but only expose an
    # ``inchikey`` as structural identity after deriving it from a parsed
    # structure below.  This distinction is important for imports: an
    # unverified user key must neither merge a compound nor lend a structure to
    # a name/formula row.
    submitted_inchikey = _clean(
        clean.get("reported_inchikey") or clean.get("inchikey")
    ).upper()
    clean["reported_inchikey"] = submitted_inchikey or None
    clean["structure_identity_verified"] = False
    clean["structure_identity_source"] = (
        "reported_unverified" if submitted_inchikey else "none"
    )

    precursor = _finite_float(clean.get("precursor_mz"))
    peaks = _extract_peaks(clean)
    clean["peaks"] = peaks
    if precursor is None or precursor <= 0:
        errors.append("missing_or_invalid_precursor_mz")
    if len(peaks) < 3:
        errors.append("fewer_than_3_fragment_peaks")

    adduct = clean.get("adduct")
    formula = neutral_formula_text(clean.get("formula"), adduct)
    smiles_structure = smiles_metadata(clean.get("smiles"))
    inchi_structure = inchi_metadata(clean.get("inchi"))
    structure_source = (
        "smiles"
        if smiles_structure.get("status") == "valid"
        else "inchi"
    )
    structure = (
        smiles_structure if structure_source == "smiles" else inchi_structure
    )
    if (
        submitted_inchikey
        and inchi_structure.get("status") == "valid"
        and _clean(inchi_structure.get("inchikey")).upper()
        == submitted_inchikey
        and (
            smiles_structure.get("status") != "valid"
            or _clean(smiles_structure.get("inchikey")).upper()
            != submitted_inchikey
        )
    ):
        # MassBank SMILES can omit stereochemistry while its InChI and
        # InChIKey retain it.  Prefer the independently matching structure
        # instead of falsely quarantining a genuine record.
        structure_source = "inchi"
        structure = inchi_structure
        if smiles_structure.get("status") == "valid":
            warnings.append(
                "smiles_stereochemistry_incomplete_inchi_identity_preferred"
            )
    if structure.get("status") == "valid":
        structure_formula = neutral_formula_text(structure.get("formula"), adduct)
        if formula and structure_formula and formula != structure_formula:
            warnings.append(f"{structure_source}_formula_mismatch_structure_removed")
            for key in ("smiles", "inchi", "inchikey"):
                clean[key] = None
        else:
            if not formula:
                formula = structure_formula
                warnings.append("formula_filled_from_smiles")
            clean["smiles"] = structure.get("canonical_smiles")
            clean["inchi"] = structure.get("inchi")
            supplied_key = submitted_inchikey
            derived_key = _clean(structure.get("inchikey")).upper()
            if supplied_key and supplied_key != derived_key:
                warnings.append(f"inchikey_repaired_from_{structure_source}")
            clean["inchikey"] = derived_key or None
            clean["structure_identity_verified"] = bool(
                derived_key and _FULL_INCHIKEY.fullmatch(derived_key)
            )
            clean["structure_identity_source"] = (
                f"{structure_source}_derived"
                if clean["structure_identity_verified"]
                else "none"
            )
    elif structure.get("status") == "invalid":
        if _clean(clean.get(structure_source)):
            warnings.append(f"invalid_{structure_source}_structure_removed")
        for key in ("smiles", "inchi", "inchikey"):
            clean[key] = None
    elif structure.get("status") == "unverifiable":
        warnings.append(f"{structure_source}_not_verified")

    if not clean.get("structure_identity_verified"):
        # Formula/name metadata may still be useful, and a literature row may
        # still be screened as a target.  The submitted key remains available
        # in ``reported_inchikey`` and the audit, but it is deliberately not
        # stored in the safe structural identity field.
        clean["inchikey"] = None
        if submitted_inchikey:
            warnings.append("unverified_inchikey_removed_from_structure_identity")

    clean["formula"] = formula or None
    exact_mass = formula_exact_mass(formula)
    if exact_mass is None and structure.get("status") == "valid" and clean.get("smiles"):
        exact_mass = _finite_float(structure.get("exact_mass"))
    clean["exact_mass"] = exact_mass

    precursor_ppm = None
    original_adduct = None
    if precursor and exact_mass:
        observed_neutral = neutral_mass(precursor, adduct)
        if observed_neutral is None:
            warnings.append("unrecognized_or_missing_adduct_precursor_not_verified")
        else:
            precursor_ppm = (observed_neutral - exact_mass) / exact_mass * 1e6
            if abs(precursor_ppm) > max_precursor_ppm:
                repaired = _repair_adduct_from_formula_precursor(
                    precursor, exact_mass, adduct, max_precursor_ppm
                )
                if repaired:
                    original_adduct = _clean(adduct)
                    clean["adduct"] = repaired[0]
                    precursor_ppm = repaired[1]
                    warnings.append("adduct_repaired_from_formula_precursor")
                else:
                    errors.append("formula_adduct_precursor_mismatch")
    elif not formula:
        warnings.append("missing_formula_precursor_not_verified")

    intensity_origin = _clean(clean.get("intensity_origin")).lower()
    if not intensity_origin:
        intensity_origin = "measured" if record.get("peaks") else "default"
    clean["intensity_origin"] = intensity_origin
    has_measured_intensity = intensity_origin in {
        "measured", "measured_consensus", "reported", "reported_relative",
        "tracefinder_measured",
    }
    if not has_measured_intensity:
        warnings.append("fragment_intensities_missing_presence_only")

    declared_tier = _clean(clean.get("source_tier")).lower()
    if not has_measured_intensity:
        source_tier = PRESENCE_ONLY_SOURCE_TIER
    elif declared_tier in _RECOGNIZED_DECLARED_SOURCE_TIERS:
        # Measured-looking intensities do not prove that a spectrum is
        # experimental.  Preserve explicit predicted/in-silico provenance.
        source_tier = declared_tier
        if require_experimental_provenance and source_tier in TRUSTED_SOURCE_TIERS:
            provenance_fields = (
                "instrument_type", "instrument", "instrument_model",
                "acquisition_method", "provenance", "dataset_accession",
                "spectrum_accession", "source_accession", "doi",
            )
            provenance_present = bool(clean.get("source_provenance_verified")) or any(
                _clean(clean.get(field)) for field in provenance_fields
            )
            if not provenance_present:
                source_tier = "user_unverified"
                warnings.append("experimental_provenance_missing_tier_downgraded")
    else:
        # Numeric fragment intensities prove only that numbers were supplied;
        # they do not establish instrument acquisition or library provenance.
        source_tier = "user_unverified"
        warnings.append("source_provenance_missing_not_experimental")
    accepted = not errors
    if not accepted:
        source_tier = QUARANTINED_SOURCE_TIER
    clean["source_tier"] = source_tier
    clean["quality_flags"] = sorted(set(errors + warnings))

    audit = {
        "name": clean.get("name"),
        "source_file": clean.get("source_file"),
        "accepted": accepted,
        "source_tier": source_tier,
        "errors": sorted(set(errors)),
        "warnings": sorted(set(warnings)),
        "precursor_ppm": round(precursor_ppm, 3) if precursor_ppm is not None else None,
        "n_peaks": len(peaks),
        "has_structure": bool(clean.get("smiles")),
        "reported_inchikey": clean.get("reported_inchikey"),
        "structure_identity_verified": bool(clean.get("structure_identity_verified")),
        "structure_identity_source": clean.get("structure_identity_source"),
        "formula": clean.get("formula"),
        "adduct": clean.get("adduct"),
        "original_adduct": original_adduct,
        "intensity_origin": intensity_origin,
    }
    return clean, audit


def assess_records(
    records: list[dict[str, Any]],
    *,
    require_experimental_provenance: bool = True,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    sanitized: list[dict[str, Any]] = []
    audits: list[dict[str, Any]] = []
    signatures: dict[tuple[Any, ...], list[int]] = defaultdict(list)

    for record in records:
        clean, audit = assess_record(
            record,
            require_experimental_provenance=require_experimental_provenance,
        )
        index = len(sanitized)
        sanitized.append(clean)
        audits.append(audit)
        signature = (
            round(float(clean.get("precursor_mz") or 0), 4),
            _clean(clean.get("adduct")),
            tuple(round(float(p[0]), 4) for p in clean.get("peaks") or []),
        )
        if signature[0] > 0 and signature[2]:
            signatures[signature].append(index)

    ambiguous_groups = 0
    for indices in signatures.values():
        names = {_clean(sanitized[i].get("name")).casefold() for i in indices}
        names.discard("")
        if len(names) <= 1:
            continue
        ambiguous_groups += 1
        for i in indices:
            flag = "shared_precursor_and_fragment_list_across_compound_names"
            sanitized[i].setdefault("quality_flags", []).append(flag)
            audits[i].setdefault("warnings", []).append(flag)

    tier_counts = Counter(audit["source_tier"] for audit in audits)
    error_counts = Counter(error for audit in audits for error in audit["errors"])
    warning_counts = Counter(warning for audit in audits for warning in audit["warnings"])
    accepted_records = [record for record, audit in zip(sanitized, audits) if audit["accepted"]]
    rejected_records = [record for record, audit in zip(sanitized, audits) if not audit["accepted"]]
    # Import/preview responses are returned through the web API.  Do not echo
    # every rejected source row (including long peak lists and structures) back
    # to the browser: large literature workbooks otherwise produce multi-MB
    # JSON responses.  The per-row audit samples below retain enough context for
    # diagnosis while the aggregate counters remain exact.
    rejected_samples = [
        {
            "name": record.get("name", ""),
            "source_file": record.get("source_file", ""),
            "source_tier": record.get("source_tier", ""),
            "quality_flags": list(record.get("quality_flags") or []),
        }
        for record in rejected_records[:50]
    ]
    return accepted_records, {
        "total": len(records),
        "accepted": len(accepted_records),
        "rejected": len(rejected_records),
        "source_tier_counts": dict(tier_counts),
        "error_counts": dict(error_counts),
        "warning_counts": dict(warning_counts),
        "ambiguous_signature_groups": ambiguous_groups,
        "provenance_held": sum(
            1
            for audit in audits
            if audit["accepted"] and audit["source_tier"] in {
                "user_unverified", "unknown", "metadata_conflict",
                "predicted", "in_silico",
            }
        ),
        "audits": audits[:200],
        "rejected_record_samples": rejected_samples,
    }


def _is_measured_intensity_origin(record: Mapping[str, Any]) -> bool:
    origin = _clean(record.get("intensity_origin")).lower()
    if not origin:
        return bool(record.get("peaks"))
    return origin in {
        "measured", "measured_consensus", "reported", "reported_relative",
        "tracefinder_measured",
    }


def is_expected_target_record(record: Mapping[str, Any]) -> bool:
    """Return whether a row belongs to the mass-list/target catalogue.

    Missing fragment intensities are the decisive boundary. A list of
    literature fragment m/z values with synthetic equal intensities is not an
    experimental spectrum and must never enter cosine Top-K search.
    """
    declared = _clean(record.get("source_tier")).lower()
    return (
        declared == PRESENCE_ONLY_SOURCE_TIER
        or not _is_measured_intensity_origin(record)
    )


def assess_target_record(
    record: dict[str, Any],
    max_precursor_ppm: float = 50.0,
) -> tuple[dict[str, Any], dict[str, Any]]:
    """Quality-assess one Expected Compound independently of spectral gates.

    A valid target may have no MS2 ions and may lack a reported precursor when
    its neutral formula/structure can generate theoretical ion m/z values.
    Formula/adduct/precursor conflicts remain visible quality findings but do
    not erase the compound from the catalogue; the runtime uses the formula
    mass and labels generated ion hypotheses separately.
    """
    clean, audit = assess_record(
        record,
        max_precursor_ppm=max_precursor_ppm,
        require_experimental_provenance=False,
    )
    clean["reported_precursor_mz"] = _finite_float(record.get("precursor_mz"))
    clean["reported_adduct"] = _clean(record.get("adduct")) or None
    clean["reported_retention_time"] = _finite_float(
        record.get("retention_time")
    )
    clean["reported_retention_time_unit"] = (
        _clean(record.get("retention_time_unit")) or None
    )
    collision_raw = (
        record.get("collision_energy_raw")
        if record.get("collision_energy_raw") is not None
        else record.get("collision_energy")
    )
    clean["reported_collision_energy_raw"] = _clean(collision_raw) or None
    clean["instrument"] = _clean(record.get("instrument")) or None

    spectral_only_errors = {
        "missing_or_invalid_precursor_mz",
        "fewer_than_3_fragment_peaks",
        "formula_adduct_precursor_mismatch",
    }
    original_errors = set(audit.get("errors") or [])
    remaining_errors = sorted(original_errors - spectral_only_errors)
    target_flags = set(clean.get("quality_flags") or [])
    target_flags.difference_update(spectral_only_errors)
    if "missing_or_invalid_precursor_mz" in original_errors:
        target_flags.add("reported_precursor_missing_theoretical_ion_required")
    if "fewer_than_3_fragment_peaks" in original_errors:
        target_flags.add("diagnostic_ions_below_spectral_library_minimum")
    if "formula_adduct_precursor_mismatch" in original_errors:
        target_flags.add("reported_ion_conflict_theoretical_ion_used")

    name = _clean(clean.get("name"))
    has_formula_mass = _finite_float(clean.get("exact_mass")) is not None
    has_reported_ion = (
        _finite_float(clean.get("reported_precursor_mz")) is not None
        and bool(_clean(clean.get("reported_adduct")))
    )
    has_identity = bool(
        name
        or clean.get("formula")
        or clean.get("smiles")
        or clean.get("inchi")
        or clean.get("reported_inchikey")
    )
    if not has_identity:
        remaining_errors.append("missing_target_identity")
        target_flags.add("missing_target_identity")

    if has_formula_mass:
        searchability = "theoretical_ion_ready"
    elif has_reported_ion:
        searchability = "reported_ion_only"
    else:
        searchability = "identity_only_unsearchable"
        target_flags.add("target_ion_mass_unavailable")
    if len(clean.get("peaks") or []) < 2:
        target_flags.add("ms1_only_expected_target")

    accepted = not remaining_errors
    clean["target_searchability"] = searchability
    clean["source_tier"] = (
        PRESENCE_ONLY_SOURCE_TIER if accepted else QUARANTINED_SOURCE_TIER
    )
    clean["quality_flags"] = sorted(target_flags | set(remaining_errors))
    audit.update({
        "accepted": accepted,
        "source_tier": clean["source_tier"],
        "errors": remaining_errors,
        "warnings": sorted(target_flags - set(remaining_errors)),
        "target_searchability": searchability,
        "n_diagnostic_ions": len(clean.get("peaks") or []),
    })
    return clean, audit


def assess_target_records(
    records: list[dict[str, Any]],
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Assess a batch of mass-list targets without MS2 spectral requirements."""
    cleaned: list[dict[str, Any]] = []
    audits: list[dict[str, Any]] = []
    for record in records:
        item, audit = assess_target_record(record)
        cleaned.append(item)
        audits.append(audit)

    accepted = [
        record for record, audit in zip(cleaned, audits) if audit["accepted"]
    ]
    rejected = [
        record for record, audit in zip(cleaned, audits) if not audit["accepted"]
    ]
    return accepted, {
        "total": len(records),
        "accepted": len(accepted),
        "rejected": len(rejected),
        "source_tier_counts": dict(Counter(
            audit["source_tier"] for audit in audits
        )),
        "error_counts": dict(Counter(
            error for audit in audits for error in audit.get("errors", [])
        )),
        "warning_counts": dict(Counter(
            warning for audit in audits for warning in audit.get("warnings", [])
        )),
        "target_searchability_counts": dict(Counter(
            audit.get("target_searchability") for audit in audits
        )),
        "ambiguous_signature_groups": 0,
        "provenance_held": 0,
        "audits": audits[:200],
        "rejected_record_samples": [
            {
                "name": record.get("name", ""),
                "source_file": record.get("source_file", ""),
                "source_tier": record.get("source_tier", ""),
                "quality_flags": list(record.get("quality_flags") or []),
            }
            for record in rejected[:50]
        ],
    }


def assess_import_records(
    records: list[dict[str, Any]],
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Route spectral records and Expected Compounds through separate gates."""
    target_input = [
        record for record in records if is_expected_target_record(record)
    ]
    spectral_input = [
        record for record in records if not is_expected_target_record(record)
    ]
    targets, target_quality = assess_target_records(target_input)
    spectra, spectral_quality = assess_records(spectral_input)

    audits = (
        list(target_quality.get("audits") or [])
        + list(spectral_quality.get("audits") or [])
    )
    rejected_samples = (
        list(target_quality.get("rejected_record_samples") or [])
        + list(spectral_quality.get("rejected_record_samples") or [])
    )
    return spectra + targets, {
        "total": len(records),
        "accepted": len(spectra) + len(targets),
        "rejected": (
            int(target_quality.get("rejected") or 0)
            + int(spectral_quality.get("rejected") or 0)
        ),
        "source_tier_counts": dict(
            Counter(target_quality.get("source_tier_counts") or {})
            + Counter(spectral_quality.get("source_tier_counts") or {})
        ),
        "error_counts": dict(
            Counter(target_quality.get("error_counts") or {})
            + Counter(spectral_quality.get("error_counts") or {})
        ),
        "warning_counts": dict(
            Counter(target_quality.get("warning_counts") or {})
            + Counter(spectral_quality.get("warning_counts") or {})
        ),
        "target_searchability_counts": target_quality.get(
            "target_searchability_counts", {}
        ),
        "ambiguous_signature_groups": int(
            target_quality.get("ambiguous_signature_groups") or 0
        ) + int(spectral_quality.get("ambiguous_signature_groups") or 0),
        "provenance_held": int(
            target_quality.get("provenance_held") or 0
        ) + int(spectral_quality.get("provenance_held") or 0),
        "audits": audits[:200],
        "rejected_record_samples": rejected_samples[:50],
        "record_types": {
            "expected_target_input": len(target_input),
            "experimental_spectrum_input": len(spectral_input),
            "expected_targets_accepted": len(targets),
            "experimental_spectra_accepted": len(spectra),
        },
    }
