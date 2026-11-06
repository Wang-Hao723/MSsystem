#!/usr/bin/env python3
"""Audit, repair, deduplicate, and incrementally import the CIeaD library.

The public CIeaD MSP export has two properties that make a generic filename-
driven import unsafe:

* records from positive/negative and CID/EAD files overlap heavily;
* MSP records omit structures, instrument provenance, ion mode, and a stable
  source accession.

This importer therefore treats the four MSP files as redundant peak evidence,
not as authoritative mode labels.  It joins them to the official merged MGF
(which carries ion mode and instrument) and to the peer-reviewed supporting
table (which carries structures, formulas, and NPClassifier classes).

Dry-run is the default.  ``--apply`` creates a consistent SQLite backup and
writes the accepted records in one transaction while holding the library lock.
"""
from __future__ import annotations

import argparse
from bisect import bisect_left, insort
from collections import Counter, defaultdict
from contextlib import nullcontext
from dataclasses import asdict, dataclass, field
from datetime import datetime
from difflib import SequenceMatcher
import gzip
from hashlib import sha256
import io
import json
import math
import os
from pathlib import Path
import re
import shutil
import sqlite3
import sys
import unicodedata
from typing import Any, Iterable, Iterator
import zipfile

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from config.settings import DB_PATH, LIBRARY_WRITE_ENABLED  # noqa: E402
from matching.library_quality import source_tier_is_searchable  # noqa: E402
from mssystem.classification import (  # noqa: E402
    RuleBasedClassifier,
    normalize_filename_class,
)
from rdkit import Chem, rdBase  # noqa: E402
from rdkit.Chem import Descriptors, rdMolDescriptors  # noqa: E402
from utils.adduct import (  # noqa: E402
    adduct_polarity,
    ion_mz,
)
from utils.spectrum import clean_peaks  # noqa: E402
from utils.storage_guard import (  # noqa: E402
    estimate_incremental_import_bytes,
    exclusive_file_lock,
    require_storage_capacity,
)


DOI = "10.1021/jasms.5c00329"
SUPPLEMENT_DOI = "10.1021/jasms.5c00329.s002"
LICENSE = "CC BY-NC 4.0"
DATASET_SOURCE = "CIeaD"
SOURCE_TIER = "experimental"
INSTRUMENT = "SCIEX ZenoTOF 7600"
INSTRUMENT_TYPE = "QTOF"
MGF_MEMBER_DEFAULT = "Library .mgf/CleaD.mgf"
PRECURSOR_TOLERANCE_DA = 0.015
PRECURSOR_MAX_PPM = 25.0
MIN_CLEAN_PEAKS = 3
TOP_N_PEAKS = 50
RELATIVE_INTENSITY_THRESHOLD = 0.01
DB_BATCH = 600

DEFAULT_STAGING = Path(
    "/Users/zizi/Desktop/Mssystem_shared_data/import_staging/"
    "ciead_official_20260729"
)
DEFAULT_MGF_ZIP = DEFAULT_STAGING / "Library.mgf.zip"
DEFAULT_METADATA = DEFAULT_STAGING / "ciead_metadata.json"
DEFAULT_REPAIRS = DEFAULT_STAGING / "pubchem_parent_repairs.json"
DEFAULT_MSP_FILES = (
    Path(
        "/Users/zizi/Desktop/MsQCModel_v2/_local_only/"
        "mgf_msp_database/MSP_data/Neg CID MSP.msp"
    ),
    Path(
        "/Users/zizi/Desktop/MsQCModel_v2/_local_only/"
        "mgf_msp_database/MSP_data/Neg EAD MSP.msp"
    ),
    Path(
        "/Users/zizi/Desktop/MsQCModel_v2/_local_only/"
        "mgf_msp_database/MSP_data/Pos CID MSP.msp"
    ),
    Path(
        "/Users/zizi/Desktop/MsQCModel_v2/_local_only/"
        "mgf_msp_database/MSP_data/Pos EAD MSP.msp"
    ),
)

_GREEK = str.maketrans(
    {
        "α": "alpha",
        "β": "beta",
        "γ": "gamma",
        "δ": "delta",
        "ε": "epsilon",
        "κ": "kappa",
        "λ": "lambda",
        "ω": "omega",
    }
)
_ZERO_WIDTH = re.compile(r"[\u200B-\u200D\u2060\uFEFF]")
_DASHES = str.maketrans("‐‑‒–—−", "------")
_FORMULA_TOKEN = re.compile(r"([A-Z][a-z]?)(\d*)")
_FULL_INCHIKEY = re.compile(r"^[A-Z]{14}-[A-Z]{10}-[A-Z]$")
_MISSING = {"", "n/a", "na", "none", "null", "-", "nan"}

# These public-library labels are mixtures, polymers, or unresolved extracts;
# assigning one convenient PubChem structure would create false identity.
UNRESOLVED_MIXTURES = {
    "Casanthranol",
    "Enocyanin",
    "Isomalt",
    "Silymarin",
    "Xylan",
}

POSITIVE_ADDUCTS = (
    "[M+H]+",
    "[M]+",
    "[M+Na]+",
    "[M+NH4]+",
    "[M+K]+",
    "[M-H2O+H]+",
    "[M+H-NH3]+",
)
NEGATIVE_ADDUCTS = (
    "[M-H]-",
    "[M]-",
    "[M+Cl]-",
    "[M+FA-H]-",
    "[M+Hac-H]-",
    "[M-H-H2O]-",
)


@dataclass
class StructureRecord:
    name: str
    aliases: tuple[str, ...]
    source_formula: str
    formula: str
    smiles: str
    inchi: str
    inchikey: str
    exact_mass: float
    compound_class_raw: str
    superclass_raw: str
    classification_v3: str
    classification_source: str
    classification_confidence: float
    structure_source: str
    pubchem_cid: str = ""


@dataclass
class MspEvidence:
    source_memberships: set[str] = field(default_factory=set)
    formulas: Counter[str] = field(default_factory=Counter)
    collision_energies: set[float] = field(default_factory=set)
    rows: int = 0


@dataclass
class Candidate:
    source_name: str
    structure: StructureRecord
    precursor_mz: float
    polarity: str
    adduct: str
    charge: int
    collision_energy: float | None
    collision_energy_raw: str
    peaks: list[list[float]]
    peak_fingerprint: str
    source_accession: str
    spectrum_id: str
    activation_method: str
    source_memberships: tuple[str, ...]
    quality_flags: list[str]
    metadata_match_method: str
    formula_evidence: str


@dataclass
class CompoundUpdate:
    compound_id: int
    old_values: dict[str, Any]
    new_values: dict[str, Any]


def clean_text(value: Any) -> str:
    text = str(value or "").strip()
    return "" if text.casefold() in _MISSING else text


def normalized_name(value: Any) -> str:
    text = unicodedata.normalize("NFKC", clean_text(value))
    text = _ZERO_WIDTH.sub("", text).translate(_DASHES).translate(_GREEK)
    text = re.sub(r"\s+", " ", text)
    return text.strip().casefold()


def simple_name(value: Any) -> str:
    text = normalized_name(value)
    text = re.sub(r"\b(alpha|beta|gamma|delta)\b", r"\1", text)
    return re.sub(r"[^a-z0-9]+", "", text)


def normalized_formula(value: Any) -> str:
    text = clean_text(value)
    if not text:
        return ""
    text = unicodedata.normalize("NFKC", text)
    text = re.sub(r"\s+", "", text)
    # Simple neutral and permanent-ion formulae only.  Salts/hydrates are
    # retained as source_formula but are not used as canonical identities.
    if "." in text or "/" in text or "x" in text:
        return ""
    terminal_charge = ""
    if text.endswith(("+", "-")):
        terminal_charge = text[-1]
        text = text[:-1]
    pos = 0
    tokens: list[tuple[str, int]] = []
    for match in _FORMULA_TOKEN.finditer(text):
        if match.start() != pos:
            return ""
        tokens.append(
            (match.group(1), int(match.group(2)) if match.group(2) else 1)
        )
        pos = match.end()
    if not tokens or pos != len(text):
        return ""
    counts: Counter[str] = Counter()
    for element, count in tokens:
        counts[element] += count
    ordered: list[str] = []
    if "C" in counts:
        count = counts.pop("C")
        ordered.append(f"C{count if count != 1 else ''}")
        if "H" in counts:
            count = counts.pop("H")
            ordered.append(f"H{count if count != 1 else ''}")
    for element in sorted(counts):
        count = counts[element]
        ordered.append(f"{element}{count if count != 1 else ''}")
    return "".join(ordered) + terminal_charge


def neutral_formula_for_comparison(value: Any) -> str:
    return normalized_formula(value).rstrip("+-")


def chunks(values: list[Any], size: int = DB_BATCH) -> Iterable[list[Any]]:
    for start in range(0, len(values), size):
        yield values[start : start + size]


def hash_file(path: Path, block_size: int = 8 * 1024 * 1024) -> str:
    digest = sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(block_size), b""):
            digest.update(block)
    return digest.hexdigest()


def raw_peak_fingerprint(peaks: Iterable[Iterable[Any]]) -> str:
    digest = sha256()
    found = False
    for peak in peaks:
        try:
            mz = float(peak[0])
            intensity = float(peak[1])
        except (TypeError, ValueError, IndexError):
            continue
        if not (
            math.isfinite(mz)
            and math.isfinite(intensity)
            and mz > 0
            and intensity > 0
        ):
            continue
        digest.update(f"{mz:.12g}\t{intensity:.12g}\n".encode())
        found = True
    return digest.hexdigest() if found else ""


def normalized_peak_fingerprint(peaks: Iterable[Iterable[Any]]) -> str:
    parsed: list[tuple[float, float]] = []
    for peak in peaks:
        try:
            mz = float(peak[0])
            intensity = float(peak[1])
        except (TypeError, ValueError, IndexError):
            continue
        if (
            math.isfinite(mz)
            and math.isfinite(intensity)
            and mz > 0
            and intensity > 0
        ):
            parsed.append((mz, intensity))
    if not parsed:
        return ""
    maximum = max(intensity for _, intensity in parsed)
    if maximum <= 0:
        return ""
    payload = "\n".join(
        f"{mz:.6f}:{intensity / maximum:.6f}"
        for mz, intensity in sorted(parsed)
    )
    return sha256(payload.encode()).hexdigest()


def source_key(name: Any, precursor_mz: Any, peaks: list[list[float]]) -> tuple:
    try:
        precursor = round(float(precursor_mz), 4)
    except (TypeError, ValueError):
        precursor = 0.0
    return normalized_name(name), precursor, raw_peak_fingerprint(peaks)


def iter_msp_records(path: Path) -> Iterator[dict[str, Any]]:
    fields: dict[str, str] = {}
    peaks: list[list[float]] = []

    def finish() -> dict[str, Any] | None:
        if not fields:
            return None
        record = dict(fields)
        record["peaks"] = list(peaks)
        return record

    with Path(path).open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            text = line.strip()
            if text.startswith("Name:") and fields:
                record = finish()
                if record is not None:
                    yield record
                fields = {}
                peaks = []
            if not text:
                continue
            if ":" in text and text[0] not in "+-.0123456789":
                key, value = text.split(":", 1)
                fields[key.strip()] = value.strip()
                continue
            parts = text.split()
            if len(parts) < 2:
                continue
            try:
                mz = float(parts[0])
                intensity = float(parts[1])
            except (TypeError, ValueError):
                continue
            if (
                math.isfinite(mz)
                and math.isfinite(intensity)
                and mz > 0
                and intensity > 0
            ):
                peaks.append([mz, intensity])
    record = finish()
    if record is not None:
        yield record


def iter_mgf_records(zip_path: Path, member: str) -> Iterator[dict[str, Any]]:
    with zipfile.ZipFile(zip_path) as archive:
        with archive.open(member) as binary:
            handle = io.TextIOWrapper(binary, encoding="utf-8", errors="replace")
            headers: dict[str, str] = {}
            peaks: list[list[float]] = []
            active = False
            for line in handle:
                text = line.strip()
                if text == "BEGIN IONS":
                    headers = {}
                    peaks = []
                    active = True
                    continue
                if text == "END IONS":
                    if headers:
                        yield {"headers": dict(headers), "peaks": peaks}
                    headers = {}
                    peaks = []
                    active = False
                    continue
                if not active or not text:
                    continue
                if "=" in text and text[0] not in "+-.0123456789":
                    key, value = text.split("=", 1)
                    headers[key.strip().upper()] = value.strip()
                    continue
                parts = text.split()
                if len(parts) < 2:
                    continue
                try:
                    mz = float(parts[0])
                    intensity = float(parts[1])
                except (TypeError, ValueError):
                    continue
                if (
                    math.isfinite(mz)
                    and math.isfinite(intensity)
                    and mz > 0
                    and intensity > 0
                ):
                    peaks.append([mz, intensity])


def load_pubchem_repairs(path: Path | None) -> dict[str, dict[str, Any]]:
    if path is None or not Path(path).exists():
        return {}
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    output: dict[str, dict[str, Any]] = {}
    for item in payload:
        properties = item.get("properties") or []
        if item.get("status") != "ok" or not properties:
            continue
        first = properties[0]
        if not first.get("SMILES"):
            continue
        output[str(item.get("source_name") or "")] = {
            **first,
            "url": item.get("url"),
            "query": item.get("query"),
        }
    return output


def _structure_from_smiles(
    smiles: str,
) -> tuple[Any, str, str, str, str, float] | None:
    with rdBase.BlockLogs():
        molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        return None
    canonical = Chem.MolToSmiles(molecule, isomericSmiles=True)
    inchi = Chem.MolToInchi(molecule)
    inchikey = Chem.MolToInchiKey(molecule)
    formula = rdMolDescriptors.CalcMolFormula(molecule)
    exact_mass = float(Descriptors.ExactMolWt(molecule))
    if not _FULL_INCHIKEY.fullmatch(inchikey):
        return None
    return molecule, canonical, inchi, inchikey, formula, exact_mass


def _split_aliases(value: Any) -> list[str]:
    text = clean_text(value)
    if not text:
        return []
    return [
        part.strip()
        for part in re.split(r"\s*[;|]\s*", text)
        if part.strip()
    ]


def load_structures(
    metadata_path: Path,
    repair_path: Path | None,
) -> tuple[list[StructureRecord], dict[str, Any], list[dict[str, Any]]]:
    payload = json.loads(Path(metadata_path).read_text(encoding="utf-8"))
    repairs = load_pubchem_repairs(repair_path)
    classifier = RuleBasedClassifier()
    records: list[StructureRecord] = []
    quarantine: list[dict[str, Any]] = []
    stats: Counter[str] = Counter()

    for row in payload.get("records") or []:
        name = clean_text(row.get("name"))
        aliases = tuple(_split_aliases(row.get("synonyms")))
        source_formula = clean_text(row.get("formula"))
        smiles = clean_text(row.get("smiles"))
        source = "CIeaD_Table_S2"
        pubchem_cid = ""

        parsed = _structure_from_smiles(smiles)
        formula_matches = False
        if parsed is not None:
            formula_matches = (
                neutral_formula_for_comparison(parsed[4])
                == neutral_formula_for_comparison(source_formula)
                and bool(neutral_formula_for_comparison(source_formula))
            )
        if name in UNRESOLVED_MIXTURES:
            stats["quarantined_unresolved_mixture_or_polymer"] += 1
            quarantine.append(
                {
                    "stage": "structure",
                    "reason": "unresolved_mixture_or_polymer",
                    "name": name,
                    "formula": source_formula,
                }
            )
            continue

        if parsed is None or not formula_matches:
            repair = repairs.get(name)
            if repair:
                repaired = _structure_from_smiles(clean_text(repair.get("SMILES")))
            else:
                repaired = None
            if repaired is None:
                reason = (
                    "invalid_smiles"
                    if parsed is None
                    else "smiles_formula_mismatch"
                )
                stats[f"quarantined_{reason}"] += 1
                quarantine.append(
                    {
                        "stage": "structure",
                        "reason": reason,
                        "name": name,
                        "formula": source_formula,
                    }
                )
                continue
            parsed = repaired
            source = "PubChem_parent_repair"
            pubchem_cid = str(repair.get("CID") or "")
            stats["pubchem_parent_repairs"] += 1

        (
            _molecule,
            canonical,
            inchi,
            inchikey,
            derived_formula,
            exact_mass,
        ) = parsed
        canonical_formula = neutral_formula_for_comparison(derived_formula)
        if not canonical_formula:
            stats["quarantined_unusable_derived_formula"] += 1
            quarantine.append(
                {
                    "stage": "structure",
                    "reason": "unusable_derived_formula",
                    "name": name,
                    "formula": derived_formula,
                }
            )
            continue

        np_class = (
            normalize_filename_class(row.get("compound_class"))
            or normalize_filename_class(row.get("superclass"))
        )
        if np_class:
            classification = np_class
            classification_source = "curated_ciead_npclassifier"
            classification_confidence = 0.97
        else:
            structural = classifier.predict(canonical)
            if structural != "Unknown":
                classification = structural
                classification_source = "smarts_ciead_structure"
                classification_confidence = 0.98
            else:
                classification = "Unknown"
                classification_source = "quality_reaudit_unresolved_v6"
                classification_confidence = 0.0

        records.append(
            StructureRecord(
                name=name,
                aliases=aliases,
                source_formula=source_formula,
                formula=canonical_formula,
                smiles=canonical,
                inchi=inchi,
                inchikey=inchikey,
                exact_mass=exact_mass,
                compound_class_raw=clean_text(row.get("compound_class")),
                superclass_raw=clean_text(row.get("superclass")),
                classification_v3=classification,
                classification_source=classification_source,
                classification_confidence=classification_confidence,
                structure_source=source,
                pubchem_cid=pubchem_cid,
            )
        )
        stats["accepted_structures"] += 1
        stats[f"classification:{classification}"] += 1
        stats[f"structure_source:{source}"] += 1

    stats["input_structures"] = len(payload.get("records") or [])
    stats["unique_inchikeys"] = len({record.inchikey for record in records})
    return records, dict(stats), quarantine


class StructureResolver:
    def __init__(self, records: list[StructureRecord]) -> None:
        self.records = records
        self.by_name: dict[str, list[int]] = defaultdict(list)
        self.by_simple: dict[str, list[int]] = defaultdict(list)
        self.by_formula: dict[str, list[int]] = defaultdict(list)
        for index, record in enumerate(records):
            names = (record.name, *record.aliases)
            for name in names:
                self.by_name[normalized_name(name)].append(index)
                self.by_simple[simple_name(name)].append(index)
            for formula in {
                neutral_formula_for_comparison(record.source_formula),
                neutral_formula_for_comparison(record.formula),
            }:
                if formula:
                    self.by_formula[formula].append(index)

    def resolve(
        self,
        source_name: str,
        formula: str,
    ) -> tuple[StructureRecord | None, str, dict[str, Any]]:
        name_key = normalized_name(source_name)
        simple_key = simple_name(source_name)
        formula_key = neutral_formula_for_comparison(formula)

        exact = list(dict.fromkeys(self.by_name.get(name_key, [])))
        if formula_key:
            exact = [
                index
                for index in exact
                if index in self.by_formula.get(formula_key, [])
            ]
        if len(exact) == 1:
            return self.records[exact[0]], "exact_name_formula", {}

        simple = list(dict.fromkeys(self.by_simple.get(simple_key, [])))
        if formula_key:
            simple = [
                index
                for index in simple
                if index in self.by_formula.get(formula_key, [])
            ]
        if len(simple) == 1:
            return self.records[simple[0]], "normalized_name_formula", {}

        pool = list(dict.fromkeys(self.by_formula.get(formula_key, [])))
        if not pool:
            return None, "no_formula_candidate", {}
        ranked: list[tuple[float, int, str]] = []
        for index in pool:
            record = self.records[index]
            candidates = (record.name, *record.aliases)
            best = max(
                (
                    SequenceMatcher(
                        None,
                        simple_key,
                        simple_name(candidate),
                    ).ratio()
                    for candidate in candidates
                    if simple_name(candidate)
                ),
                default=0.0,
            )
            ranked.append((best, index, record.name))
        ranked.sort(key=lambda item: (-item[0], item[2]))
        best = ranked[0]
        second_score = ranked[1][0] if len(ranked) > 1 else 0.0
        acceptable = (
            best[0] >= 0.90
            and (best[0] - second_score >= 0.035 or best[0] >= 0.975)
        ) or (len(pool) == 1 and best[0] >= 0.65)
        details = {
            "best_score": round(best[0], 4),
            "second_score": round(second_score, 4),
            "candidates": [
                {"name": name, "score": round(score, 4)}
                for score, _index, name in ranked[:5]
            ],
        }
        if not acceptable:
            return None, "ambiguous_fuzzy_name_formula", details
        return self.records[best[1]], "fuzzy_name_formula", details


def build_msp_index(
    paths: list[Path],
) -> tuple[dict[tuple, MspEvidence], dict[str, Any]]:
    index: dict[tuple, MspEvidence] = {}
    stats: Counter[str] = Counter()
    by_file: dict[str, Counter[str]] = defaultdict(Counter)
    for path in paths:
        label = path.stem
        for record in iter_msp_records(path):
            stats["raw_rows"] += 1
            by_file[label]["raw_rows"] += 1
            try:
                precursor = float(record.get("PrecursorMZ"))
            except (TypeError, ValueError):
                stats["invalid_precursor"] += 1
                continue
            peaks = record.get("peaks") or []
            key = source_key(record.get("Name"), precursor, peaks)
            if not key[0] or not key[2]:
                stats["invalid_source_key"] += 1
                continue
            evidence = index.setdefault(key, MspEvidence())
            evidence.source_memberships.add(label)
            formula = clean_text(record.get("Formula"))
            if formula:
                evidence.formulas[formula] += 1
            try:
                energy = float(record.get("Collision_energy"))
            except (TypeError, ValueError):
                energy = None
            if energy is not None and math.isfinite(energy):
                evidence.collision_energies.add(energy)
            evidence.rows += 1
            by_file[label]["indexed_rows"] += 1
        by_file[label]["unique_memberships"] = sum(
            label in evidence.source_memberships for evidence in index.values()
        )
    stats["unique_source_records"] = len(index)
    stats["records_repeated_across_files"] = sum(
        len(evidence.source_memberships) > 1 for evidence in index.values()
    )
    stats["records_in_all_four_files"] = sum(
        len(evidence.source_memberships) == 4 for evidence in index.values()
    )
    return index, {
        "counts": dict(stats),
        "by_file": {
            name: dict(counts) for name, counts in sorted(by_file.items())
        },
    }


def infer_adduct(
    exact_mass: float,
    precursor_mz: float,
    polarity: str,
) -> tuple[str | None, float | None]:
    candidates = POSITIVE_ADDUCTS if polarity == "POS" else NEGATIVE_ADDUCTS
    ranked: list[tuple[float, float, str]] = []
    for adduct in candidates:
        theoretical = ion_mz(exact_mass, adduct)
        if theoretical is None or theoretical <= 0:
            continue
        error_da = precursor_mz - theoretical
        ppm = error_da / theoretical * 1e6
        ranked.append((abs(ppm), abs(error_da), adduct))
    if not ranked:
        return None, None
    ppm_abs, error_abs, adduct = min(ranked)
    tolerance_da = max(
        0.003,
        precursor_mz * PRECURSOR_MAX_PPM / 1e6,
    )
    if ppm_abs > PRECURSOR_MAX_PPM or error_abs > tolerance_da:
        return None, ppm_abs
    return adduct, ppm_abs


def _activation_method(memberships: set[str]) -> str:
    methods = {
        method
        for name in memberships
        for method in ("CID", "EAD")
        if method.casefold() in name.casefold()
    }
    if len(methods) == 1:
        return next(iter(methods))
    return "unresolved"


def _choose_formula(evidence: MspEvidence) -> tuple[str, bool]:
    if not evidence.formulas:
        return "", False
    ordered = evidence.formulas.most_common()
    conflict = len(ordered) > 1
    return ordered[0][0], conflict


def _collision_energy(
    evidence: MspEvidence,
) -> tuple[float | None, str, bool]:
    if not evidence.collision_energies:
        return None, "", False
    raw = sorted(evidence.collision_energies)
    absolute = {abs(value) for value in raw}
    if len(absolute) != 1:
        return None, ",".join(f"{value:g}" for value in raw), True
    return next(iter(absolute)), ",".join(f"{value:g}" for value in raw), False


def prepare_candidates(
    mgf_zip: Path,
    mgf_member: str,
    msp_index: dict[tuple, MspEvidence],
    resolver: StructureResolver,
) -> tuple[list[Candidate], list[dict[str, Any]], dict[str, Any]]:
    candidates: list[Candidate] = []
    quarantine: list[dict[str, Any]] = []
    stats: Counter[str] = Counter()
    match_methods: Counter[str] = Counter()
    mode_counts: Counter[str] = Counter()

    for record_number, record in enumerate(
        iter_mgf_records(mgf_zip, mgf_member),
        start=1,
    ):
        stats["mgf_rows"] += 1
        headers = record["headers"]
        source_name = clean_text(headers.get("NAME"))
        instrument = clean_text(headers.get("SOURCE_INSTRUMENT"))
        if not instrument.startswith("ZenoTOF"):
            stats["excluded_non_zenotof"] += 1
            continue
        stats["zenotof_rows"] += 1
        try:
            precursor_mz = float(clean_text(headers.get("PEPMASS")).split()[0])
        except (TypeError, ValueError, IndexError):
            stats["quarantined_invalid_precursor"] += 1
            quarantine.append(
                {
                    "stage": "spectrum",
                    "reason": "invalid_precursor",
                    "record": record_number,
                    "name": source_name,
                }
            )
            continue
        ion_mode = clean_text(headers.get("IONMODE")).upper()
        polarity = (
            "POS"
            if ion_mode.startswith("POS")
            else "NEG"
            if ion_mode.startswith("NEG")
            else ""
        )
        if not polarity:
            stats["quarantined_invalid_ion_mode"] += 1
            quarantine.append(
                {
                    "stage": "spectrum",
                    "reason": "invalid_ion_mode",
                    "record": record_number,
                    "name": source_name,
                }
            )
            continue
        mode_counts[polarity] += 1
        key = source_key(source_name, precursor_mz, record["peaks"])
        evidence = msp_index.get(key)
        if evidence is None:
            stats["quarantined_not_cross_validated_in_msp"] += 1
            quarantine.append(
                {
                    "stage": "spectrum",
                    "reason": "not_cross_validated_in_four_msp_files",
                    "record": record_number,
                    "name": source_name,
                    "precursor_mz": precursor_mz,
                    "polarity": polarity,
                }
            )
            continue
        stats["msp_cross_validated"] += 1
        source_formula, formula_conflict = _choose_formula(evidence)
        structure, match_method, match_details = resolver.resolve(
            source_name,
            source_formula,
        )
        match_methods[match_method] += 1
        if structure is None:
            stats["quarantined_structure_unresolved"] += 1
            quarantine.append(
                {
                    "stage": "spectrum",
                    "reason": "structure_unresolved",
                    "record": record_number,
                    "name": source_name,
                    "formula": source_formula,
                    "precursor_mz": precursor_mz,
                    "metadata_match_method": match_method,
                    "metadata_match_details": match_details,
                }
            )
            continue

        adduct, precursor_ppm = infer_adduct(
            structure.exact_mass,
            precursor_mz,
            polarity,
        )
        if not adduct:
            stats["quarantined_formula_precursor_mismatch"] += 1
            quarantine.append(
                {
                    "stage": "spectrum",
                    "reason": "formula_adduct_precursor_mismatch",
                    "record": record_number,
                    "name": source_name,
                    "resolved_name": structure.name,
                    "formula": structure.formula,
                    "precursor_mz": precursor_mz,
                    "polarity": polarity,
                    "best_ppm_abs": (
                        round(precursor_ppm, 3)
                        if precursor_ppm is not None
                        else None
                    ),
                }
            )
            continue
        if adduct_polarity(adduct) != polarity:
            raise RuntimeError("internal adduct polarity mismatch")

        cleaned = clean_peaks(
            record["peaks"],
            precursor_mz,
            top_n=TOP_N_PEAKS,
            rel_threshold=RELATIVE_INTENSITY_THRESHOLD,
            precursor_charge=1,
        )
        if len(cleaned) < MIN_CLEAN_PEAKS:
            stats["quarantined_fewer_than_3_clean_peaks"] += 1
            quarantine.append(
                {
                    "stage": "spectrum",
                    "reason": "fewer_than_3_clean_peaks",
                    "record": record_number,
                    "name": source_name,
                    "resolved_name": structure.name,
                    "precursor_mz": precursor_mz,
                    "clean_peaks": len(cleaned),
                }
            )
            continue

        fingerprint = normalized_peak_fingerprint(cleaned)
        collision_energy, collision_raw, ce_conflict = _collision_energy(
            evidence
        )
        activation = _activation_method(evidence.source_memberships)
        flags = [
            "ciead_official_experimental_library",
            f"doi:{DOI}",
            f"license:{LICENSE}",
            "msp_four_file_content_deduplicated",
            "polarity_from_official_mgf",
            "adduct_from_structure_precursor_mass_balance",
            f"structure_source:{structure.structure_source}",
            f"metadata_match:{match_method}",
        ]
        if structure.superclass_raw:
            flags.append(
                "ciead_npclassifier_superclass:"
                f"{structure.superclass_raw}"
            )
        if structure.compound_class_raw:
            flags.append(
                "ciead_npclassifier_class:"
                f"{structure.compound_class_raw}"
            )
        if activation == "unresolved":
            flags.append("activation_method_unresolved_cid_ead_export_overlap")
        else:
            flags.append(f"activation_method:{activation}")
        if formula_conflict:
            flags.append("msp_formula_membership_conflict")
        if ce_conflict:
            flags.append("collision_energy_membership_conflict")
        if collision_energy is not None:
            flags.append("collision_energy_unit_unreported")
        if structure.pubchem_cid:
            flags.append(f"pubchem_cid:{structure.pubchem_cid}")
        if precursor_ppm is not None:
            flags.append(f"precursor_ppm:{precursor_ppm:.3f}")

        accession_payload = "|".join(
            (
                structure.inchikey,
                f"{precursor_mz:.6f}",
                polarity,
                adduct,
                fingerprint,
            )
        )
        accession_hash = sha256(accession_payload.encode()).hexdigest()
        accession = f"CIeaD:{accession_hash[:24]}"
        candidates.append(
            Candidate(
                source_name=source_name,
                structure=structure,
                precursor_mz=precursor_mz,
                polarity=polarity,
                adduct=adduct,
                charge=1 if polarity == "POS" else -1,
                collision_energy=collision_energy,
                collision_energy_raw=collision_raw,
                peaks=cleaned,
                peak_fingerprint=fingerprint,
                source_accession=accession,
                spectrum_id=accession,
                activation_method=activation,
                source_memberships=tuple(
                    sorted(evidence.source_memberships)
                ),
                quality_flags=sorted(set(flags)),
                metadata_match_method=match_method,
                formula_evidence=source_formula,
            )
        )
        stats["quality_accepted_before_dedup"] += 1

    # One identity/acquisition/content signature is one spectrum.  Merge public
    # export repetitions without collapsing different collision energies or
    # different fragment-intensity patterns.
    selected_by_signature: dict[tuple[Any, ...], Candidate] = {}
    for candidate in candidates:
        signature = (
            candidate.structure.inchikey,
            round(candidate.precursor_mz, 5),
            candidate.polarity,
            candidate.adduct,
            candidate.collision_energy,
            candidate.peak_fingerprint,
        )
        previous = selected_by_signature.get(signature)
        if previous is None:
            selected_by_signature[signature] = candidate
            continue
        stats["deduplicated_clean_identity_spectrum"] += 1
        previous.quality_flags = sorted(
            set(previous.quality_flags) | set(candidate.quality_flags)
        )
        previous.source_memberships = tuple(
            sorted(
                set(previous.source_memberships)
                | set(candidate.source_memberships)
            )
        )

    selected = list(selected_by_signature.values())
    broad_groups: dict[tuple[Any, ...], set[str]] = defaultdict(set)
    for candidate in selected:
        broad_groups[
            (
                round(candidate.precursor_mz, 5),
                candidate.polarity,
                candidate.peak_fingerprint,
            )
        ].add(candidate.structure.inchikey)
    ambiguous_groups = {
        key: identities
        for key, identities in broad_groups.items()
        if len(identities) > 1
    }
    for candidate in selected:
        key = (
            round(candidate.precursor_mz, 5),
            candidate.polarity,
            candidate.peak_fingerprint,
        )
        if key in ambiguous_groups:
            candidate.quality_flags = sorted(
                set(candidate.quality_flags)
                | {"shared_spectrum_across_structural_isomers"}
            )
    stats["shared_spectrum_isomer_groups"] = len(ambiguous_groups)
    stats["rows_in_shared_spectrum_isomer_groups"] = sum(
        key
        in ambiguous_groups
        for candidate in selected
        for key in [
            (
                round(candidate.precursor_mz, 5),
                candidate.polarity,
                candidate.peak_fingerprint,
            )
        ]
    )
    stats["selected_after_source_dedup"] = len(selected)
    stats["quarantined_rows"] = len(quarantine)
    selected_mode_counts = Counter(
        candidate.polarity for candidate in selected
    )
    selected_activation_counts = Counter(
        candidate.activation_method for candidate in selected
    )
    selected_adduct_counts = Counter(
        candidate.adduct for candidate in selected
    )
    selected_class_counts = Counter(
        candidate.structure.classification_v3 for candidate in selected
    )
    return selected, quarantine, {
        "counts": dict(stats),
        "metadata_match_methods": dict(match_methods),
        "mode_counts": dict(mode_counts),
        "selected_mode_counts": dict(selected_mode_counts),
        "selected_activation_counts": dict(selected_activation_counts),
        "selected_adduct_counts": dict(selected_adduct_counts),
        "selected_class_counts": dict(selected_class_counts),
        "selected_unique_inchikeys": len(
            {candidate.structure.inchikey for candidate in selected}
        ),
    }


def readonly_connection(path: Path) -> sqlite3.Connection:
    connection = sqlite3.connect(
        f"file:{Path(path).resolve()}?mode=ro",
        uri=True,
        timeout=60,
    )
    connection.row_factory = sqlite3.Row
    connection.execute("PRAGMA query_only=ON")
    connection.execute("PRAGMA busy_timeout=60000")
    return connection


def _valid_db_structure(smiles: Any, inchikey: str) -> bool:
    text = clean_text(smiles)
    if not text:
        return False
    parsed = _structure_from_smiles(text)
    return bool(parsed and parsed[3] == inchikey)


def load_database_state(
    connection: sqlite3.Connection,
    candidates: list[Candidate],
) -> tuple[
    dict[str, int],
    dict[int, sqlite3.Row],
    dict[tuple[str, str, str], list[float]],
    dict[str, int],
]:
    inchikeys = sorted({candidate.structure.inchikey for candidate in candidates})
    rows_by_key: dict[str, list[sqlite3.Row]] = defaultdict(list)
    rows_by_id: dict[int, sqlite3.Row] = {}
    for batch in chunks(inchikeys):
        placeholders = ",".join("?" for _ in batch)
        query = f"""
            SELECT id, compound_name, molecular_formula, smiles, inchi,
                   inchikey, exact_mass, classification, compound_class,
                   classification_v3, classification_source,
                   classification_confidence, reported_inchikey,
                   structure_identity_verified, structure_identity_source
            FROM Compounds
            WHERE upper(trim(inchikey)) IN ({placeholders})
        """
        for row in connection.execute(query, batch):
            key = clean_text(row["inchikey"]).upper()
            rows_by_key[key].append(row)
            rows_by_id[int(row["id"])] = row

    preferred: dict[str, int] = {}
    for key, rows in rows_by_key.items():
        ranked = sorted(
            rows,
            key=lambda row: (
                -int(bool(row["structure_identity_verified"])),
                -int(_valid_db_structure(row["smiles"], key)),
                -int(bool(clean_text(row["classification_v3"]))),
                int(row["id"]),
            ),
        )
        preferred[key] = int(ranked[0]["id"])

    content_index: dict[tuple[str, str, str], list[float]] = defaultdict(list)
    relevant_ids = sorted(rows_by_id)
    stats = {
        "candidate_inchikeys": len(inchikeys),
        "matched_inchikeys": len(preferred),
        "relevant_compound_rows": len(relevant_ids),
        "relevant_spectra": 0,
        "searchable_spectra": 0,
    }
    for batch in chunks(relevant_ids, 350):
        placeholders = ",".join("?" for _ in batch)
        query = f"""
            SELECT m.compound_id, s.id, s.precursor_mz, s.adduct, s.charge,
                   s.peaks_json, s.source_tier
            FROM Mapping m
            JOIN Spectra s ON s.id=m.spectrum_id
            WHERE m.compound_id IN ({placeholders})
        """
        seen: set[tuple[int, int]] = set()
        for row in connection.execute(query, batch):
            compound_id = int(row["compound_id"])
            spectrum_id = int(row["id"])
            if (compound_id, spectrum_id) in seen:
                continue
            seen.add((compound_id, spectrum_id))
            stats["relevant_spectra"] += 1
            tier = clean_text(row["source_tier"]).lower()
            if not source_tier_is_searchable(tier):
                continue
            try:
                precursor = float(row["precursor_mz"])
                peaks = json.loads(row["peaks_json"] or "[]")
            except (TypeError, ValueError, json.JSONDecodeError):
                continue
            if not math.isfinite(precursor) or precursor <= 0:
                continue
            polarity = adduct_polarity(row["adduct"])
            if polarity not in {"POS", "NEG"}:
                try:
                    charge = int(row["charge"])
                except (TypeError, ValueError):
                    charge = 0
                polarity = "NEG" if charge < 0 else "POS" if charge > 0 else None
            if polarity not in {"POS", "NEG"}:
                continue
            cleaned = clean_peaks(
                peaks,
                precursor,
                top_n=TOP_N_PEAKS,
                rel_threshold=RELATIVE_INTENSITY_THRESHOLD,
                precursor_charge=row["charge"],
            )
            fingerprint = normalized_peak_fingerprint(cleaned)
            if not fingerprint:
                continue
            inchikey = clean_text(rows_by_id[compound_id]["inchikey"]).upper()
            content_index[(inchikey, polarity, fingerprint)].append(precursor)
            stats["searchable_spectra"] += 1
    for values in content_index.values():
        values.sort()
    return preferred, rows_by_id, content_index, stats


def _has_precursor_match(
    values: list[float] | None,
    target: float,
) -> bool:
    if not values:
        return False
    position = bisect_left(values, target - PRECURSOR_TOLERANCE_DA)
    return (
        position < len(values)
        and abs(values[position] - target) <= PRECURSOR_TOLERANCE_DA
    )


def select_database_novel(
    candidates: list[Candidate],
    content_index: dict[tuple[str, str, str], list[float]],
) -> tuple[list[Candidate], dict[str, Any]]:
    selected: list[Candidate] = []
    stats: Counter[str] = Counter()
    seen_accessions: set[str] = set()
    for candidate in candidates:
        key = (
            candidate.structure.inchikey,
            candidate.polarity,
            candidate.peak_fingerprint,
        )
        if _has_precursor_match(
            content_index.get(key),
            candidate.precursor_mz,
        ):
            stats["duplicate_existing_database_content"] += 1
            continue
        if candidate.source_accession in seen_accessions:
            stats["duplicate_source_accession"] += 1
            continue
        selected.append(candidate)
        seen_accessions.add(candidate.source_accession)
        insort(content_index[key], candidate.precursor_mz)
        stats["selected_for_import"] += 1
    return selected, dict(stats)


def _classification_trusted(value: Any) -> bool:
    source = clean_text(value).lower()
    return source.startswith(
        (
            "curated",
            "manual",
            "smarts",
            "structure",
            "rule",
            "existing",
            "database_curated",
        )
    )


def compound_updates_for_existing(
    rows_by_id: dict[int, sqlite3.Row],
    preferred: dict[str, int],
    candidates: list[Candidate],
) -> list[CompoundUpdate]:
    structures = {
        candidate.structure.inchikey: candidate.structure
        for candidate in candidates
    }
    updates: list[CompoundUpdate] = []
    for inchikey, compound_id in preferred.items():
        if inchikey not in structures:
            continue
        row = rows_by_id[compound_id]
        structure = structures[inchikey]
        old: dict[str, Any] = {}
        new: dict[str, Any] = {}

        fill_fields = {
            "smiles": structure.smiles,
            "inchi": structure.inchi,
            "molecular_formula": structure.formula,
            "exact_mass": structure.exact_mass,
            "reported_inchikey": structure.inchikey,
            "structure_identity_verified": 1,
            "structure_identity_source": (
                f"{structure.structure_source}_smiles_derived"
            ),
        }
        for column, value in fill_fields.items():
            current = row[column]
            if column in {
                "structure_identity_verified",
                "structure_identity_source",
            }:
                if column == "structure_identity_verified" and bool(current):
                    continue
                if column == "structure_identity_source" and clean_text(current):
                    continue
            elif clean_text(current):
                continue
            old[column] = current
            new[column] = value

        current_class = clean_text(row["classification_v3"])
        current_source = clean_text(row["classification_source"])
        if (
            structure.classification_v3 != "Unknown"
            and (
                not current_class
                or current_class == "Unknown"
                or not _classification_trusted(current_source)
            )
        ):
            for column, value in (
                ("classification", structure.classification_v3),
                ("compound_class", structure.classification_v3),
                ("classification_v3", structure.classification_v3),
                ("classification_source", structure.classification_source),
                (
                    "classification_confidence",
                    structure.classification_confidence,
                ),
            ):
                old[column] = row[column]
                new[column] = value
        if new:
            updates.append(
                CompoundUpdate(
                    compound_id=compound_id,
                    old_values=old,
                    new_values=new,
                )
            )
    return updates


def write_json_atomic(path: Path, payload: Any) -> None:
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name(
        f".{destination.name}.{os.getpid()}.tmp"
    )
    try:
        with temporary.open("w", encoding="utf-8") as handle:
            json.dump(payload, handle, ensure_ascii=False, indent=2)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, destination)
    finally:
        temporary.unlink(missing_ok=True)


def write_jsonl_gzip(path: Path, rows: Iterable[dict[str, Any]]) -> None:
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name(
        f".{destination.name}.{os.getpid()}.tmp"
    )
    try:
        with gzip.open(temporary, "wt", encoding="utf-8") as handle:
            for row in rows:
                handle.write(json.dumps(row, ensure_ascii=False) + "\n")
        os.replace(temporary, destination)
    finally:
        temporary.unlink(missing_ok=True)


def consistent_backup(db_path: Path) -> Path:
    stamp = datetime.now().astimezone().strftime("%Y%m%dT%H%M%S%z")
    backup = db_path.with_name(f"{db_path.name}.pre_ciead_{stamp}.bak")
    source = sqlite3.connect(
        f"file:{db_path.resolve()}?mode=ro",
        uri=True,
        timeout=60,
    )
    destination = sqlite3.connect(str(backup))
    try:
        source.backup(destination, pages=8192)
        destination.commit()
    finally:
        destination.close()
        source.close()
    return backup


def _new_compound_values(structure: StructureRecord) -> tuple[Any, ...]:
    return (
        structure.name,
        structure.formula,
        structure.smiles,
        structure.inchi,
        structure.inchikey,
        structure.exact_mass,
        structure.classification_v3,
        structure.classification_v3,
        structure.classification_v3,
        structure.classification_source,
        structure.classification_confidence,
        None,
        structure.inchikey,
        1,
        f"{structure.structure_source}_smiles_derived",
    )


def apply_import(
    db_path: Path,
    candidates: list[Candidate],
    preferred: dict[str, int],
    rows_by_id: dict[int, sqlite3.Row],
    compound_updates: list[CompoundUpdate],
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    if not LIBRARY_WRITE_ENABLED:
        raise RuntimeError("MSSYSTEM_LIBRARY_WRITE_ENABLED forbids library writes")
    connection = sqlite3.connect(str(db_path), timeout=60)
    connection.row_factory = sqlite3.Row
    connection.execute("PRAGMA foreign_keys=ON")
    connection.execute("PRAGMA busy_timeout=60000")
    before_spectra = int(
        connection.execute("SELECT COUNT(*) FROM Spectra").fetchone()[0]
    )
    before_compounds = int(
        connection.execute("SELECT COUNT(*) FROM Compounds").fetchone()[0]
    )
    ledger_rows: list[dict[str, Any]] = []
    try:
        connection.execute("BEGIN IMMEDIATE")
        structure_ids = dict(preferred)
        unique_structures = {
            candidate.structure.inchikey: candidate.structure
            for candidate in candidates
        }
        new_compounds = 0
        for inchikey, structure in sorted(unique_structures.items()):
            if inchikey in structure_ids:
                continue
            cursor = connection.execute(
                """
                INSERT INTO Compounds (
                    compound_name, molecular_formula, smiles, inchi, inchikey,
                    exact_mass, classification, compound_class,
                    classification_v3, classification_source,
                    classification_confidence, fragmentation_pathway,
                    reported_inchikey, structure_identity_verified,
                    structure_identity_source
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                """,
                _new_compound_values(structure),
            )
            compound_id = int(cursor.lastrowid)
            structure_ids[inchikey] = compound_id
            new_compounds += 1
            ledger_rows.append(
                {
                    "action": "insert_compound",
                    "compound_id": compound_id,
                    "inchikey": inchikey,
                    "name": structure.name,
                }
            )

        applied_updates = 0
        for update in compound_updates:
            if not update.new_values:
                continue
            assignments = ", ".join(
                f"{column}=?" for column in update.new_values
            )
            values = list(update.new_values.values()) + [update.compound_id]
            connection.execute(
                f"UPDATE Compounds SET {assignments} WHERE id=?",
                values,
            )
            applied_updates += 1
            ledger_rows.append(
                {
                    "action": "update_compound",
                    **asdict(update),
                }
            )

        inserted_spectra = 0
        for candidate in candidates:
            flags = list(candidate.quality_flags)
            flags.extend(
                f"source_membership:{name}"
                for name in candidate.source_memberships
            )
            cursor = connection.execute(
                """
                INSERT INTO Spectra (
                    spectrum_id, source_file, precursor_mz, retention_time,
                    collision_energy, charge, adduct, peaks_json,
                    dataset_source, source_tier, intensity_origin,
                    quality_flags_json, collision_energy_raw,
                    collision_energy_unit, collision_energy_mode,
                    collision_energy_min, collision_energy_max,
                    instrument_type, instrument, ion_mode, source_accession
                ) VALUES (
                    ?, ?, ?, NULL, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?,
                    ?, ?, ?, ?
                )
                """,
                (
                    candidate.spectrum_id,
                    "|".join(candidate.source_memberships),
                    candidate.precursor_mz,
                    candidate.collision_energy,
                    candidate.charge,
                    candidate.adduct,
                    json.dumps(candidate.peaks, ensure_ascii=False),
                    DATASET_SOURCE,
                    SOURCE_TIER,
                    "measured",
                    json.dumps(sorted(set(flags)), ensure_ascii=False),
                    candidate.collision_energy_raw,
                    "",
                    "single" if candidate.collision_energy is not None else "",
                    candidate.collision_energy,
                    candidate.collision_energy,
                    INSTRUMENT_TYPE,
                    INSTRUMENT,
                    "Positive" if candidate.polarity == "POS" else "Negative",
                    candidate.source_accession,
                ),
            )
            spectrum_row_id = int(cursor.lastrowid)
            compound_id = structure_ids[candidate.structure.inchikey]
            connection.execute(
                "INSERT INTO Mapping (compound_id, spectrum_id, score) "
                "VALUES (?, ?, NULL)",
                (compound_id, spectrum_row_id),
            )
            inserted_spectra += 1
            ledger_rows.append(
                {
                    "action": "insert_spectrum",
                    "spectrum_row_id": spectrum_row_id,
                    "spectrum_id": candidate.spectrum_id,
                    "compound_id": compound_id,
                    "inchikey": candidate.structure.inchikey,
                }
            )

        after_spectra = int(
            connection.execute("SELECT COUNT(*) FROM Spectra").fetchone()[0]
        )
        after_compounds = int(
            connection.execute("SELECT COUNT(*) FROM Compounds").fetchone()[0]
        )
        if after_spectra - before_spectra != inserted_spectra:
            raise RuntimeError("spectra count invariant failed")
        if after_compounds - before_compounds != new_compounds:
            raise RuntimeError("compound count invariant failed")
        foreign_key_errors = list(connection.execute("PRAGMA foreign_key_check"))
        if foreign_key_errors:
            raise RuntimeError(
                f"foreign_key_check failed: {foreign_key_errors[:5]}"
            )
        connection.commit()
        return (
            {
                "before_spectra": before_spectra,
                "after_spectra": after_spectra,
                "inserted_spectra": inserted_spectra,
                "before_compounds": before_compounds,
                "after_compounds": after_compounds,
                "inserted_compounds": new_compounds,
                "updated_compounds": applied_updates,
            },
            ledger_rows,
        )
    except Exception:
        connection.rollback()
        raise
    finally:
        connection.close()


def validate_database(db_path: Path, spectrum_ids: list[str]) -> dict[str, Any]:
    connection = readonly_connection(db_path)
    try:
        quick_check = [
            str(row[0])
            for row in connection.execute("PRAGMA quick_check")
        ]
        foreign_keys = list(connection.execute("PRAGMA foreign_key_check"))
        found = 0
        tier_counts: Counter[str] = Counter()
        mode_counts: Counter[str] = Counter()
        for batch in chunks(spectrum_ids):
            placeholders = ",".join("?" for _ in batch)
            query = f"""
                SELECT source_tier, ion_mode
                FROM Spectra
                WHERE spectrum_id IN ({placeholders})
            """
            for row in connection.execute(query, batch):
                found += 1
                tier_counts[clean_text(row["source_tier"])] += 1
                mode_counts[clean_text(row["ion_mode"])] += 1
        return {
            "quick_check": quick_check,
            "foreign_key_errors": len(foreign_keys),
            "requested_spectrum_ids": len(spectrum_ids),
            "found_spectrum_ids": found,
            "source_tiers": dict(tier_counts),
            "ion_modes": dict(mode_counts),
        }
    finally:
        connection.close()


def default_report_path(db_path: Path) -> Path:
    stamp = datetime.now().astimezone().strftime("%Y%m%dT%H%M%S%z")
    return (
        db_path.parent
        / "migration_ledgers"
        / f"ciead_import_{stamp}.json"
    )


def run(args: argparse.Namespace) -> None:
    db_path = Path(args.db).expanduser().resolve()
    mgf_zip = Path(args.mgf_zip).expanduser().resolve()
    metadata_path = Path(args.metadata).expanduser().resolve()
    repair_path = (
        Path(args.repairs).expanduser().resolve()
        if args.repairs
        else None
    )
    msp_paths = [Path(path).expanduser().resolve() for path in args.msp]
    required = [db_path, mgf_zip, metadata_path, *msp_paths]
    if repair_path is not None:
        required.append(repair_path)
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise FileNotFoundError(f"missing required files: {missing}")

    report_path = (
        Path(args.report).expanduser().resolve()
        if args.report
        else default_report_path(db_path)
    )
    quarantine_path = report_path.with_name(
        f"{report_path.stem}.quarantine.jsonl.gz"
    )
    rollback_path = report_path.with_name(
        f"{report_path.stem}.rollback.jsonl.gz"
    )
    source_hashes = {
        str(path): {
            "sha256": hash_file(path),
            "bytes": path.stat().st_size,
        }
        for path in [mgf_zip, metadata_path, *msp_paths]
        + ([repair_path] if repair_path is not None else [])
    }

    structures, structure_stats, structure_quarantine = load_structures(
        metadata_path,
        repair_path,
    )
    resolver = StructureResolver(structures)
    msp_index, msp_stats = build_msp_index(msp_paths)
    candidates, spectrum_quarantine, source_stats = prepare_candidates(
        mgf_zip,
        args.mgf_member,
        msp_index,
        resolver,
    )
    quarantine = structure_quarantine + spectrum_quarantine

    with readonly_connection(db_path) as connection:
        preferred, rows_by_id, content_index, db_stats = load_database_state(
            connection,
            candidates,
        )
    novel, db_selection_stats = select_database_novel(
        candidates,
        content_index,
    )
    compound_updates = compound_updates_for_existing(
        rows_by_id,
        preferred,
        novel,
    )

    report: dict[str, Any] = {
        "type": "ciead_incremental_import",
        "schema_version": 1,
        "created_at": datetime.now().astimezone().isoformat(
            timespec="seconds"
        ),
        "apply_requested": bool(args.apply),
        "database": str(db_path),
        "doi": DOI,
        "supplement_doi": SUPPLEMENT_DOI,
        "license": LICENSE,
        "rdkit_version": rdBase.rdkitVersion,
        "source_hashes": source_hashes,
        "structure_audit": structure_stats,
        "msp_audit": msp_stats,
        "source_join_audit": source_stats,
        "database_audit": db_stats,
        "database_selection": db_selection_stats,
        "planned_new_spectra": len(novel),
        "planned_new_compounds": len(
            {
                candidate.structure.inchikey
                for candidate in novel
                if candidate.structure.inchikey not in preferred
            }
        ),
        "planned_compound_updates": len(compound_updates),
        "quarantine_rows": len(quarantine),
        "quarantine_path": str(quarantine_path),
        "rollback_path": str(rollback_path),
    }
    write_jsonl_gzip(quarantine_path, quarantine)

    print(json.dumps(report, ensure_ascii=False, indent=2))
    if not args.apply:
        write_json_atomic(report_path, report)
        print(f"\nDRY RUN: report={report_path}")
        print(f"quarantine={quarantine_path}")
        return

    lock_path = db_path.with_suffix(db_path.suffix + ".write.lock")
    planned_bytes = estimate_incremental_import_bytes(
        [db_path],
        source_bytes=sum(path.stat().st_size for path in msp_paths)
        + mgf_zip.stat().st_size,
    )
    with exclusive_file_lock(lock_path):
        report["storage"] = require_storage_capacity(
            db_path,
            planned_bytes=planned_bytes,
            operation="CIeaD incremental import",
        )
        backup = consistent_backup(db_path)
        report["backup"] = str(backup)

        # Recheck content while holding the writer lock, so a concurrent import
        # cannot invalidate the dry-run selection.
        with readonly_connection(db_path) as connection:
            preferred, rows_by_id, fresh_content, db_stats_locked = (
                load_database_state(connection, candidates)
            )
        novel_locked, selection_locked = select_database_novel(
            candidates,
            fresh_content,
        )
        updates_locked = compound_updates_for_existing(
            rows_by_id,
            preferred,
            novel_locked,
        )
        report["database_audit_locked"] = db_stats_locked
        report["database_selection_locked"] = selection_locked
        report["planned_new_spectra_locked"] = len(novel_locked)

        results, ledger_rows = apply_import(
            db_path,
            novel_locked,
            preferred,
            rows_by_id,
            updates_locked,
        )
        header = {
            "header": {
                "type": "ciead_import_rollback",
                "schema_version": 1,
                "database": str(db_path),
                "backup": str(backup),
                "created_at": datetime.now().astimezone().isoformat(
                    timespec="seconds"
                ),
                "rows": len(ledger_rows),
            }
        }
        write_jsonl_gzip(rollback_path, [header, *ledger_rows])
        report["result"] = results
        report["validation"] = validate_database(
            db_path,
            [candidate.spectrum_id for candidate in novel_locked],
        )
        report["completed_at"] = datetime.now().astimezone().isoformat(
            timespec="seconds"
        )
        write_json_atomic(report_path, report)
        print(json.dumps(report["result"], ensure_ascii=False, indent=2))
        print(json.dumps(report["validation"], ensure_ascii=False, indent=2))
        print(f"report={report_path}")
        print(f"rollback={rollback_path}")
        print(f"backup={backup}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db", default=str(DB_PATH))
    parser.add_argument("--mgf-zip", default=str(DEFAULT_MGF_ZIP))
    parser.add_argument("--mgf-member", default=MGF_MEMBER_DEFAULT)
    parser.add_argument("--metadata", default=str(DEFAULT_METADATA))
    parser.add_argument("--repairs", default=str(DEFAULT_REPAIRS))
    parser.add_argument(
        "--msp",
        action="append",
        default=None,
        help="Repeat for each of the four public MSP files",
    )
    parser.add_argument("--report")
    parser.add_argument("--apply", action="store_true")
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    if args.msp is None:
        args.msp = [str(path) for path in DEFAULT_MSP_FILES]
    run(args)


if __name__ == "__main__":
    main()
