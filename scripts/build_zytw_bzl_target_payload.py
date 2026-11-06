#!/usr/bin/env python3
"""Build the corrected zytw + BZL literature-target migration payload.

The source workbooks contain compound-presence metadata, not experimental
MS/MS spectra.  This builder therefore emits only ``literature_presence_only``
target records.  Every accepted structure is re-derived with RDKit and must
pass the source molecular-formula and mass checks.  Ambiguous names remain
formula-only.

The BZL workbook is read from JSON exported with the bundled artifact-tool;
this avoids a second spreadsheet parser with different type/coercion rules.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
import zipfile
from collections import Counter
from pathlib import Path
from typing import Any

from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, rdMolDescriptors


PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from matching.library_quality import formula_exact_mass, normalize_formula_text


RDLogger.DisableLog("rdApp.*")

DATASET_SOURCE = "user_contrib:literature_template_zytw_bzl"
SOURCE_TIER = "literature_presence_only"
PROTON_MASS = 1.007276466621
PLANTS = {
    "zytw": "Ligularia cymbulifera",
    "bzl": "Scutellaria barbata",
}

COCONUT_2024_RECORD = "https://doi.org/10.5281/zenodo.13692394"
COCONUT_2024_FILE = "coconut-09-2024.csv"
COCONUT_2024_MD5 = "439fca469248b2eb5194ffb75ddf62ed"
COCONUT_2021_RECORD = "https://doi.org/10.5281/zenodo.5336447"
COCONUT_2021_FILE = "COCONUT4MetFrag_april.csv"
COCONUT_2021_MD5 = "c3925f2f8da49047101967438816c9db"

CURRENT_COCONUT_IDS = {
    "CNP0520159.1",
    "CNP0520159.2",
    "CNP0041905.1",
    "CNP0041905.2",
    "CNP0070290.1",
    "CNP0522150.1",
    "CNP0527537.1",
    "CNP0580670.1",
    "CNP0497085.1",
}

# These three versioned identifiers disappeared from the COCONUT 2.0 snapshot.
# The legacy official export retains their base identifiers and connectivity.
# The provenance explicitly records that narrower identity scope.
LEGACY_COCONUT_CONNECTIVITY = {
    "CNP0136943.1": {
        "base_identifier": "CNP0136943",
        "smiles": (
            "O=C(OC(CC(O)C(=C)C1CC2OC2(C)C(OC(=O)C(=CC)C)C1O)"
            "C(O)(C)C)C(=CC)C"
        ),
        "reported_inchikey": "XNXMMOGHWCPTQZ-UHFFFAOYSA-N",
    },
    "CNP0125049.1": {
        "base_identifier": "CNP0125049",
        "smiles": "O=C(C=CC1=CC=C(O)C=C1)C2=CC=C(O)C=C2OC",
        "reported_inchikey": "PACBGANPVNHGNP-UHFFFAOYSA-N",
    },
    "CNP0211261.1": {
        "base_identifier": "CNP0211261",
        "smiles": "O=C(C=CSC)NCCO",
        "reported_inchikey": "NBRCIHAGNJRRRI-UHFFFAOYSA-N",
    },
}

# Explicit nomenclature rewrites supplied to official OPSIN.  These do not
# invent a scaffold: each query is a deterministic expansion of the source
# systematic name, including the source-stated stereochemistry.
BZL_OPSIN_QUERIES = {
    0: "(2S)-2,3-dihydro-2-(4-hydroxyphenyl)-6-methyl-4H-pyran-4-one",
    5: (
        "7-[(2S,3R,4S,5S,6R)-4,5-dihydroxy-6-(hydroxymethyl)-3-"
        "[(2S,3R,4S,5R,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]"
        "oxyoxan-2-yl]oxy-5-hydroxy-2-(4-methoxyphenyl)chromen-4-one"
    ),
    88: "(2S)-2',7-dihydroxy-5,8-dimethoxyflavanone",
    90: "(2S)-4',5,6,7-tetrahydroxyflavanone",
    140: "2',3',5,7-tetrahydroxyflavone",
    141: "3-Methyl-6-(1-methylethyl)-cyclohexene",
    142: "5-hydroxy-2-(4-hydroxy-3,5-dimethoxyphenyl)-6,7-dimethoxychromen-4-one",
    144: "(6R,7aS)-3,6-dimethyl-5,6,7,7a-tetrahydro-4H-benzofuran-2-one",
    147: "5,7-dihydroxy-8-methoxy-2-(2-methoxyphenyl)-4H-chromen-4-one",
    148: "5-hydroxy-7,8-dimethoxy-2-(4-methoxyphenyl)chromone",
    159: (
        "(8S,9S,10R,13R,14S,17R)-10,13-dimethyl-17-[(E,1R,4R)-"
        "1,4,5-trimethylhex-2-enyl]-2,3,8,9,11,12,14,15,16,17-"
        "decahydro-1H-cyclopenta[a]phenanthrene"
    ),
    160: (
        "(3S,8S,9S,10R,13R,14S,17R)-10,13-dimethyl-17-[(E,1R,4R)-"
        "1,4,5-trimethylhex-2-enyl]-2,3,8,9,11,12,14,15,16,17-"
        "decahydro-1H-cyclopenta[a]phenanthren-3-ol"
    ),
    170: "Stigmastan-3,5,22-triene",
    188: "4,4-dimethylcholesta-6,22,24-triene",
}

FORMULA_CORRECTIONS = {
    # Name + reported neutral mass independently support these corrections.
    147: ("C17H14O6", "source_formula_C17H14O7_to_C17H14O6"),
    200: ("C8H16O", "source_formula_C8H6O_to_C8H16O"),
    # The source name is non-stereospecific; retain formula-only identity.
    207: ("C16H18O9", "formula_filled_for_chlorogenic_name_mass_only"),
}

MASS_CORRECTIONS = {
    186: "reported_average_mass_121.12_to_formula_consistent_122.12",
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _text(value: Any) -> str:
    if value is None:
        return ""
    return str(value).replace("\u00a0", " ").strip()


def _formula(value: Any) -> str:
    return normalize_formula_text(_text(value))


def _artifact_rows(path: Path) -> list[dict[str, Any]]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    sheets = payload.get("sheets") or []
    if len(sheets) != 1:
        raise ValueError(f"expected one worksheet in {path}, found {len(sheets)}")
    values = sheets[0].get("values") or []
    if len(values) < 2:
        raise ValueError(f"no data rows in {path}")
    header = [str(value) for value in values[0]]
    return [dict(zip(header, row)) for row in values[1:]]


def _structure(smiles: str, reported_key: str = "") -> dict[str, Any]:
    mol = Chem.MolFromSmiles(_text(smiles))
    if mol is None:
        raise ValueError(f"RDKit could not parse SMILES: {smiles[:120]}")
    canonical = Chem.MolToSmiles(mol, isomericSmiles=True)
    key = Chem.MolToInchiKey(mol)
    if reported_key and key != _text(reported_key).upper():
        raise ValueError(f"InChIKey mismatch: derived={key}, reported={reported_key}")
    return {
        "smiles": canonical,
        "inchi": Chem.MolToInchi(mol),
        "inchikey": key,
        "formula": rdMolDescriptors.CalcMolFormula(mol),
        "exact_mass": round(float(Descriptors.ExactMolWt(mol)), 6),
        "average_mass": round(float(Descriptors.MolWt(mol)), 6),
        "assigned_chiral_centers": sum(
            label != "?"
            for _, label in Chem.FindMolChiralCenters(
                mol, includeUnassigned=True, includeCIP=True
            )
        ),
        "unassigned_chiral_centers": sum(
            label == "?"
            for _, label in Chem.FindMolChiralCenters(
                mol, includeUnassigned=True, includeCIP=True
            )
        ),
    }


def _formula_gate(
    structure: dict[str, Any], expected_formula: str, reported_mass: float | None,
    *, mass_kind: str,
) -> tuple[bool, str]:
    if expected_formula and structure["formula"] != expected_formula:
        return False, (
            f"formula mismatch {structure['formula']} != {expected_formula}"
        )
    if reported_mass is not None:
        calculated = (
            structure["exact_mass"] if mass_kind == "exact"
            else structure["average_mass"]
        )
        tolerance = 0.05 if mass_kind == "exact" else 0.5
        if abs(float(calculated) - float(reported_mass)) > tolerance:
            return False, (
                f"{mass_kind} mass mismatch {calculated} vs {reported_mass}"
            )
    return True, "ok"


def _opsin_structure(query: str, cache: dict[str, Any]) -> dict[str, Any]:
    cached = cache.get(query)
    if cached:
        return _structure(cached["smiles"], cached.get("inchikey", ""))
    url = (
        "https://www.ebi.ac.uk/opsin/ws/"
        + urllib.parse.quote(query, safe="")
        + ".json"
    )
    request = urllib.request.Request(
        url, headers={"User-Agent": "MSsystem-V10-zytw-bzl-migration/1.0"}
    )
    try:
        with urllib.request.urlopen(request, timeout=40) as response:
            payload = json.load(response)
    except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError) as exc:
        raise RuntimeError(f"OPSIN failed for {query!r}: {exc}") from exc
    result = _structure(payload.get("smiles") or payload.get("canonicalSmiles") or "")
    cache[query] = {"smiles": result["smiles"], "inchikey": result["inchikey"]}
    time.sleep(0.2)
    return result


def _current_coconut_structures(snapshot_zip: Path) -> dict[str, dict[str, Any]]:
    found: dict[str, dict[str, Any]] = {}
    with zipfile.ZipFile(snapshot_zip) as archive:
        members = [
            name for name in archive.namelist()
            if name.endswith(".csv") and not name.startswith("__MACOSX/")
        ]
        if len(members) != 1:
            raise ValueError(f"expected one CSV in {snapshot_zip}, found {members}")
        with archive.open(members[0]) as raw:
            reader = csv.DictReader(io.TextIOWrapper(raw, encoding="utf-8-sig"))
            for row in reader:
                identifier = _text(row.get("identifier"))
                if identifier not in CURRENT_COCONUT_IDS:
                    continue
                structure = _structure(
                    _text(row.get("canonical_smiles")),
                    _text(row.get("standard_inchi_key")),
                )
                structure.update({
                    "reported_inchi": _text(row.get("standard_inchi")),
                    "identifier": identifier,
                })
                found[identifier] = structure
    missing = sorted(CURRENT_COCONUT_IDS - set(found))
    if missing:
        raise ValueError(f"COCONUT snapshot missing identifiers: {missing}")
    return found


def _base_record(
    *, dataset: str, source_index: int, name: str, formula: str,
    structure: dict[str, Any] | None, identity_source: str, confidence: str,
    reason: str, provenance_detail: dict[str, Any], source_mass: float | None,
    cas: str = "", dataset_accession: str = "",
) -> dict[str, Any]:
    exact_mass = formula_exact_mass(formula) if formula else None
    if structure:
        if structure["formula"] != formula:
            raise ValueError(
                f"final formula/structure mismatch for {dataset}[{source_index}] {name}: "
                f"{formula} != {structure['formula']}"
            )
        exact_mass = float(structure["exact_mass"])
    precursor = round(float(exact_mass) + PROTON_MASS, 6) if exact_mass else None
    provenance = {
        "dataset": dataset,
        "source_index": source_index,
        "reported_name": name,
        "reported_formula": formula,
        "reported_mass": source_mass,
        "cas": cas or None,
        **provenance_detail,
    }
    return {
        "dataset": dataset,
        "source_index": source_index,
        "name": name,
        "formula": formula or None,
        "exact_mass": round(float(exact_mass), 6) if exact_mass else None,
        "precursor_mz": precursor,
        "charge": 1 if precursor else None,
        "adduct": "[M+H]+" if precursor else None,
        "plant_source": PLANTS[dataset],
        "smiles": structure["smiles"] if structure else None,
        "inchi": structure["inchi"] if structure else None,
        "inchikey": structure["inchikey"] if structure else None,
        "structure_identity_verified": bool(structure),
        "structure_identity_source": identity_source if structure else "none",
        "confidence": confidence if structure else None,
        "resolve_reason": reason,
        "source_tier": SOURCE_TIER,
        "dataset_source": DATASET_SOURCE,
        "instrument": "Literature-curated",
        "intensity_origin": "default",
        "peaks": [],
        "source_mass": source_mass,
        "cas": cas or None,
        "dataset_accession": dataset_accession or None,
        "provenance": json.dumps(provenance, sort_keys=True, ensure_ascii=False),
    }


def build_zytw(
    zytw_path: Path, coconut: dict[str, dict[str, Any]],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    rows = list(csv.DictReader(zytw_path.open(encoding="utf-8-sig", newline="")))
    if len(rows) != 118:
        raise ValueError(f"expected 118 zytw rows, found {len(rows)}")
    output: list[dict[str, Any]] = []
    audits: list[dict[str, Any]] = []
    for index, row in enumerate(rows):
        name = _text(row.get("Name"))
        formula = _formula(row.get("Formula"))
        reported_mass = float(row["AvgMass_file"]) if _text(row.get("AvgMass_file")) else None
        structure: dict[str, Any] | None = None
        identity_source = ""
        confidence = ""
        reason = _text(row.get("ResolveReason")) or "unresolved"
        detail: dict[str, Any] = {
            "source_file": zytw_path.name,
            "source_structure_status": _text(row.get("StructureResolved")),
        }
        if _text(row.get("StructureResolved")).casefold() == "yes":
            structure = _structure(
                _text(row.get("SMILES")), _text(row.get("InChIKey"))
            )
            ok, message = _formula_gate(
                structure, formula, reported_mass, mass_kind="exact"
            )
            if not ok:
                raise ValueError(f"zytw[{index}] failed v3 structure gate: {message}")
            identity_source = _text(row.get("Provenance")) or "zytw_v3_rdkit_verified"
            confidence = _text(row.get("Confidence")) or "high"
            reason = "resolved_in_zytw_v3"
        elif name in CURRENT_COCONUT_IDS:
            structure = coconut[name]
            ok, message = _formula_gate(
                structure, formula, reported_mass, mass_kind="exact"
            )
            if not ok:
                raise ValueError(f"zytw[{index}] COCONUT gate failed: {message}")
            identity_source = "coconut_2_0_exact_identifier_rdkit"
            confidence = "high"
            reason = "resolved_from_coconut_2_0_exact_identifier"
            detail.update({
                "resolver": "COCONUT 2.0 official snapshot",
                "resolver_record": COCONUT_2024_RECORD,
                "resolver_file": COCONUT_2024_FILE,
                "resolver_file_md5": COCONUT_2024_MD5,
                "resolver_identifier": name,
            })
        elif name in LEGACY_COCONUT_CONNECTIVITY:
            evidence = LEGACY_COCONUT_CONNECTIVITY[name]
            structure = _structure(
                evidence["smiles"], evidence["reported_inchikey"]
            )
            ok, message = _formula_gate(
                structure, formula, reported_mass, mass_kind="exact"
            )
            if not ok:
                raise ValueError(f"zytw[{index}] legacy COCONUT gate failed: {message}")
            identity_source = "coconut_2021_base_identifier_connectivity_rdkit"
            confidence = "medium_connectivity"
            reason = "resolved_from_legacy_coconut_base_identifier"
            detail.update({
                "resolver": "COCONUT legacy official snapshot",
                "resolver_record": COCONUT_2021_RECORD,
                "resolver_file": COCONUT_2021_FILE,
                "resolver_file_md5": COCONUT_2021_MD5,
                "resolver_identifier": evidence["base_identifier"],
                "identity_scope": "base_identifier_connectivity; version suffix absent",
            })
        record = _base_record(
            dataset="zytw", source_index=index, name=name, formula=formula,
            structure=structure, identity_source=identity_source,
            confidence=confidence, reason=reason, provenance_detail=detail,
            source_mass=reported_mass, dataset_accession=(name if name.startswith("CNP") else ""),
        )
        output.append(record)
        audits.append({
            "dataset": "zytw", "source_index": index, "name": name,
            "formula": formula, "resolved": bool(structure),
            "identity_source": record["structure_identity_source"],
            "reason": reason,
        })
    return output, audits


def build_bzl(
    original_path: Path, enriched_path: Path, opsin_cache: dict[str, Any],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    original = _artifact_rows(original_path)
    enriched = _artifact_rows(enriched_path)
    if len(original) != 228 or len(enriched) != 228:
        raise ValueError(
            f"expected 228 BZL rows, found original={len(original)}, enriched={len(enriched)}"
        )
    output: list[dict[str, Any]] = []
    audits: list[dict[str, Any]] = []
    for index, (source, old) in enumerate(zip(original, enriched)):
        name = _text(source.get("Name"))
        source_formula = _formula(source.get("Formula"))
        formula = source_formula
        source_mass = float(source["Mass"]) if _text(source.get("Mass")) else None
        cas = _text(source.get("CAS"))
        corrections: list[str] = []
        if index in FORMULA_CORRECTIONS:
            formula, correction = FORMULA_CORRECTIONS[index]
            corrections.append(correction)
        if index in MASS_CORRECTIONS:
            corrections.append(MASS_CORRECTIONS[index])

        structure: dict[str, Any] | None = None
        identity_source = ""
        confidence = ""
        reason = "unresolved_formula_only" if formula else "unresolved_identity_only"
        demotion_reason = ""

        old_smiles = _text(old.get("SMILES"))
        if old_smiles:
            candidate = _structure(old_smiles, _text(old.get("InChIKey")))
            expected = formula or candidate["formula"]
            ok, message = _formula_gate(
                candidate, expected, source_mass, mass_kind="average"
            )
            if index == 186 and candidate["formula"] == "C7H6O2":
                ok, message = True, "source mass typo; name and formula agree"
            if ok:
                structure = candidate
                formula = formula or candidate["formula"]
                identity_source = _text(old.get("IdentitySource")) or "pubchem_rdkit_verified"
                confidence = "high" if cas else "medium_name_exact"
                reason = "retained_after_full_rdkit_formula_mass_revalidation"
            else:
                demotion_reason = message
                corrections.append(f"old_structure_demoted:{message}")

        if index in BZL_OPSIN_QUERIES:
            query = BZL_OPSIN_QUERIES[index]
            candidate = _opsin_structure(query, opsin_cache)
            expected = formula or candidate["formula"]
            ok, message = _formula_gate(
                candidate, expected, source_mass, mass_kind="average"
            )
            if not ok:
                raise ValueError(f"BZL[{index}] OPSIN gate failed: {message}")
            structure = candidate
            formula = candidate["formula"]
            identity_source = "opsin_ebi_systematic_name_rdkit"
            confidence = (
                "high" if candidate["unassigned_chiral_centers"] == 0
                else "medium_connectivity"
            )
            reason = "resolved_from_source_systematic_name_via_opsin"

        detail: dict[str, Any] = {
            "source_file": original_path.name,
            "enrichment_file": enriched_path.name,
            "source_formula_original": source_formula or None,
            "corrections": corrections,
        }
        if index in BZL_OPSIN_QUERIES:
            detail.update({
                "resolver": "OPSIN official web service",
                "resolver_url": "https://www.ebi.ac.uk/opsin/ws",
                "resolver_query": BZL_OPSIN_QUERIES[index],
            })
        if demotion_reason:
            detail["demoted_structure_reason"] = demotion_reason
        record = _base_record(
            dataset="bzl", source_index=index, name=name, formula=formula,
            structure=structure, identity_source=identity_source,
            confidence=confidence, reason=reason, provenance_detail=detail,
            source_mass=source_mass, cas=cas,
        )
        output.append(record)
        audits.append({
            "dataset": "bzl", "source_index": index, "name": name,
            "source_formula": source_formula or None, "final_formula": formula or None,
            "resolved": bool(structure),
            "identity_source": record["structure_identity_source"],
            "reason": reason, "corrections": corrections,
        })
    return output, audits


def _write_csv(path: Path, records: list[dict[str, Any]]) -> None:
    fragment_columns = [f"fragment_mz_{index}" for index in range(1, 16)]
    columns = [
        "spectrum_name", "precursor_mz", "charge", "collision_energy",
        "retention_time", "instrument", "adduct", "component_class",
        "plant_source", "smiles", "inchi", "inchikey", "formula",
        "source_tier", "dataset_source", "structure_identity_verified",
        "structure_identity_source", "confidence", "resolve_reason",
        "source_row_index", "source_dataset", "reported_mass", "cas",
        "provenance",
    ] + fragment_columns
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        for index, record in enumerate(records):
            row = {
                "spectrum_name": record["name"],
                "precursor_mz": record["precursor_mz"] or "",
                "charge": record["charge"] or "",
                "collision_energy": "",
                "retention_time": "",
                "instrument": record["instrument"],
                "adduct": record["adduct"] or "",
                "component_class": "",
                "plant_source": record["plant_source"],
                "smiles": record["smiles"] or "",
                "inchi": record["inchi"] or "",
                "inchikey": record["inchikey"] or "",
                "formula": record["formula"] or "",
                "source_tier": SOURCE_TIER,
                "dataset_source": DATASET_SOURCE,
                "structure_identity_verified": int(record["structure_identity_verified"]),
                "structure_identity_source": record["structure_identity_source"],
                "confidence": record["confidence"] or "",
                "resolve_reason": record["resolve_reason"],
                "source_row_index": index,
                "source_dataset": record["dataset"],
                "reported_mass": record["source_mass"] or "",
                "cas": record["cas"] or "",
                "provenance": record["provenance"],
            }
            row.update({column: "" for column in fragment_columns})
            writer.writerow(row)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--zytw-v3", type=Path, required=True)
    parser.add_argument("--bzl-original-json", type=Path, required=True)
    parser.add_argument("--bzl-enriched-json", type=Path, required=True)
    parser.add_argument("--coconut-snapshot-zip", type=Path, required=True)
    parser.add_argument("--opsin-cache", type=Path, required=True)
    parser.add_argument("--output-payload", type=Path, required=True)
    parser.add_argument("--output-audit", type=Path, required=True)
    parser.add_argument("--output-template-csv", type=Path, required=True)
    parser.add_argument("--output-zytw-csv", type=Path, required=True)
    parser.add_argument("--output-bzl-csv", type=Path, required=True)
    args = parser.parse_args()

    input_paths = [
        args.zytw_v3, args.bzl_original_json, args.bzl_enriched_json,
        args.coconut_snapshot_zip,
    ]
    for path in input_paths:
        if not path.is_file():
            raise FileNotFoundError(path)
    snapshot_md5 = hashlib.md5(args.coconut_snapshot_zip.read_bytes()).hexdigest()
    if snapshot_md5 != COCONUT_2024_MD5:
        raise ValueError(
            f"COCONUT snapshot MD5 mismatch: {snapshot_md5} != {COCONUT_2024_MD5}"
        )

    opsin_cache = (
        json.loads(args.opsin_cache.read_text(encoding="utf-8"))
        if args.opsin_cache.exists() else {}
    )
    coconut = _current_coconut_structures(args.coconut_snapshot_zip)
    zytw, zytw_audit = build_zytw(args.zytw_v3, coconut)
    bzl, bzl_audit = build_bzl(
        args.bzl_original_json, args.bzl_enriched_json, opsin_cache
    )
    records = zytw + bzl
    if len(records) != 346:
        raise ValueError(f"expected 346 combined rows, found {len(records)}")

    for row_index, record in enumerate(records):
        record["source_row_index"] = row_index
        record["source_file"] = args.output_template_csv.name

    input_hashes = {str(path): _sha256(path) for path in input_paths}
    payload = {
        "schema_version": 1,
        "dataset_source": DATASET_SOURCE,
        "source_tier": SOURCE_TIER,
        "expected_record_count": 346,
        "input_sha256": input_hashes,
        "records": records,
    }
    audits = zytw_audit + bzl_audit
    audit = {
        "schema_version": 1,
        "dataset_source": DATASET_SOURCE,
        "inputs": input_hashes,
        "summary": {
            "source_rows": len(records),
            "resolved_total": sum(bool(row["smiles"]) for row in records),
            "unresolved_total": sum(not row["smiles"] for row in records),
            "zytw_rows": len(zytw),
            "zytw_resolved": sum(bool(row["smiles"]) for row in zytw),
            "bzl_rows": len(bzl),
            "bzl_resolved": sum(bool(row["smiles"]) for row in bzl),
            "identity_source_counts": dict(Counter(
                row["structure_identity_source"] for row in records
            )),
            "formula_ready": sum(bool(row["formula"]) for row in records),
            "identity_only": sum(not row["formula"] for row in records),
            "metadata_correction_rows": sum(
                bool(row.get("corrections")) for row in audits
            ),
        },
        "rows": audits,
    }

    for path, data in (
        (args.output_payload, payload),
        (args.output_audit, audit),
    ):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(
            json.dumps(data, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
    args.opsin_cache.parent.mkdir(parents=True, exist_ok=True)
    args.opsin_cache.write_text(
        json.dumps(opsin_cache, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    _write_csv(args.output_template_csv, records)
    _write_csv(args.output_zytw_csv, zytw)
    _write_csv(args.output_bzl_csv, bzl)

    print(json.dumps(audit["summary"], ensure_ascii=False, indent=2))
    print(f"payload={args.output_payload}")
    print(f"audit={args.output_audit}")
    print(f"template={args.output_template_csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
