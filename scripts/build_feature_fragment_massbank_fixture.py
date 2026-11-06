#!/usr/bin/env python3
"""Build the explicit V15 feature-fragment fixture from MassBank records.

The script downloads only the 30 allow-listed accessions below.  It never
discovers or reads project data directories.  Each record retains repository
provenance, a raw-response hash, an acquisition-aware hash, and its full peak
list under the record's CC BY license.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
import urllib.request
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from matching.feature_fragments.schema import stable_compound_name
from utils.collision_energy import parse_collision_energy


DEFAULT_OUTPUT = (
    ROOT
    / "tests"
    / "fixtures"
    / "feature_fragments"
    / "massbank_antwerp_phospholipids_v1.json"
)
API_ROOT = "https://massbank.eu/MassBank-api/records"
FETCH_DATE = "2026-08-26"
EXPECTED_INSTRUMENT = "Agilent 6560 QTOF"
EXPECTED_ADDUCT = "[M+H]+"
CLASS_ACCESSIONS = {
    "Glycerophosphocholines": [
        "MSBNK-Antwerp_Univ-METOX_N105906_EF88",
        "MSBNK-Antwerp_Univ-METOX_N105906_F638",
        "MSBNK-Antwerp_Univ-METOX_N105906_FB57",
        "MSBNK-Antwerp_Univ-METOX_N107406_EF88",
        "MSBNK-Antwerp_Univ-METOX_N107406_F638",
        "MSBNK-Antwerp_Univ-METOX_N107406_FB57",
        "MSBNK-Antwerp_Univ-METOX_N107506_EF88",
        "MSBNK-Antwerp_Univ-METOX_N107506_F638",
        "MSBNK-Antwerp_Univ-METOX_N107506_FB57",
        "MSBNK-Antwerp_Univ-METOX_N107606_EF88",
        "MSBNK-Antwerp_Univ-METOX_N107606_F638",
        "MSBNK-Antwerp_Univ-METOX_N107606_FB57",
        "MSBNK-Antwerp_Univ-METOX_N107706_EF88",
        "MSBNK-Antwerp_Univ-METOX_N107706_F638",
        "MSBNK-Antwerp_Univ-METOX_N107706_FB57",
        "MSBNK-Antwerp_Univ-METOX_N107806_EF88",
        "MSBNK-Antwerp_Univ-METOX_N107806_F638",
        "MSBNK-Antwerp_Univ-METOX_N107806_FB57",
    ],
    "Glycerophosphoethanolamines": [
        "MSBNK-Antwerp_Univ-METOX_N106006_EF88",
        "MSBNK-Antwerp_Univ-METOX_N106006_F638",
        "MSBNK-Antwerp_Univ-METOX_N106006_FB57",
        "MSBNK-Antwerp_Univ-METOX_N106106_EF88",
        "MSBNK-Antwerp_Univ-METOX_N106106_F638",
        "MSBNK-Antwerp_Univ-METOX_N106106_FB57",
        "MSBNK-Antwerp_Univ-METOX_N107906_EF88",
        "MSBNK-Antwerp_Univ-METOX_N107906_F638",
        "MSBNK-Antwerp_Univ-METOX_N107906_FB57",
    ],
    "Glycerophosphoserines": [
        "MSBNK-Antwerp_Univ-METOX_N106206_EF88",
        "MSBNK-Antwerp_Univ-METOX_N106206_F638",
        "MSBNK-Antwerp_Univ-METOX_N106206_FB57",
    ],
}


def _canonical_bytes(value: Any) -> bytes:
    return json.dumps(
        value, ensure_ascii=False, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")


def _sha256(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def _tags(items: list[dict[str, Any]] | None) -> dict[str, str]:
    return {
        str(item.get("subtag")): str(item.get("value"))
        for item in (items or [])
        if item.get("subtag") and item.get("value") is not None
    }


def _links(items: list[dict[str, Any]] | None) -> dict[str, str]:
    return {
        str(item.get("database")): str(item.get("identifier"))
        for item in (items or [])
        if item.get("database") and item.get("identifier") is not None
    }


def _float_pair(text: str) -> tuple[float, float]:
    numbers = [float(token) for token in re.findall(r"\d+(?:\.\d+)?", text)]
    if len(numbers) < 2:
        raise ValueError(f"cannot parse mass range: {text!r}")
    return min(numbers), max(numbers)


def _fetch(accession: str) -> tuple[dict[str, Any], bytes]:
    url = f"{API_ROOT}/{accession}"
    request = urllib.request.Request(
        url,
        headers={"User-Agent": "Mssystem-V15-feature-fragment-audit/1.0"},
    )
    with urllib.request.urlopen(request, timeout=60) as response:
        raw = response.read()
    return json.loads(raw), raw


def _peak_rows(record: dict[str, Any]) -> list[list[float]]:
    raw_peaks = record.get("peak", {}).get("peak", {}).get("values", [])
    peaks: list[list[float]] = []
    for item in raw_peaks:
        mz = float(item["mz"])
        intensity = float(item["intensity"])
        relative = float(item.get("rel", 0.0))
        peaks.append([mz, intensity, relative])
    if not peaks:
        raise ValueError(f"{record.get('accession')}: missing peak list")
    return sorted(peaks, key=lambda peak: peak[0])


def _peak_hash(peaks: list[list[float]]) -> str:
    base = max(peak[1] for peak in peaks)
    normalized = [[round(peak[0], 5), round(peak[1] / base, 6)] for peak in peaks]
    return _sha256(_canonical_bytes(normalized))


def _chemont_class(links: dict[str, str]) -> str:
    value = links.get("ChemOnt", "")
    return value.rsplit(";", 1)[-1].strip()


def _convert_record(
    record: dict[str, Any],
    raw: bytes,
    expected_class: str,
) -> dict[str, Any]:
    accession = str(record.get("accession") or "")
    compound = record.get("compound", {})
    acquisition = record.get("acquisition", {})
    acquisition_ms = acquisition.get("mass_spectrometry", {})
    acquisition_tags = _tags(acquisition_ms.get("subtags"))
    focused_tags = _tags(record.get("mass_spectrometry", {}).get("focused_ion"))
    links = _links(compound.get("link"))
    names = [
        str(name).strip() for name in compound.get("names", []) if str(name).strip()
    ]
    if not names or not stable_compound_name(names[0]):
        raise ValueError(
            f"{accession}: first compound name is missing or identifier-only"
        )
    name = names[0]
    inchikey = links.get("INCHIKEY", "")
    if not re.fullmatch(r"[A-Z]{14}-[A-Z]{10}-[A-Z]", inchikey):
        raise ValueError(f"{accession}: invalid full InChIKey")
    if _chemont_class(links) != expected_class:
        raise ValueError(
            f"{accession}: ChemOnt class {_chemont_class(links)!r} != {expected_class!r}"
        )
    if acquisition.get("instrument") != EXPECTED_INSTRUMENT:
        raise ValueError(f"{accession}: unexpected instrument")
    if acquisition_tags.get("FRAGMENTATION_MODE") != "CID":
        raise ValueError(f"{accession}: unexpected activation")
    if focused_tags.get("PRECURSOR_TYPE") != EXPECTED_ADDUCT:
        raise ValueError(f"{accession}: unexpected adduct")
    ion_mode = str(acquisition_ms.get("ion_mode") or "").upper()
    if ion_mode != "POSITIVE":
        raise ValueError(f"{accession}: unexpected polarity {ion_mode!r}")
    parsed_ce = parse_collision_energy(acquisition_tags.get("COLLISION_ENERGY"))
    if parsed_ce["unit"] != "eV" or parsed_ce["mode"] != "single":
        raise ValueError(f"{accession}: CE must be a single reported eV value")
    precursor_mz = float(focused_tags["PRECURSOR_M/Z"])
    mass_range = _float_pair(acquisition_tags["MASS_RANGE_M/Z"])
    peaks = _peak_rows(record)
    peak_hash = _peak_hash(peaks)
    acquisition_identity = {
        "connectivity": inchikey[:14],
        "precursor_mz": round(precursor_mz, 6),
        "adduct": EXPECTED_ADDUCT,
        "polarity": "positive",
        "instrument": EXPECTED_INSTRUMENT,
        "activation": "beam_cid",
        "collision_cell_region": "axial_hexapole",
        "collision_energy": {
            "value": parsed_ce["value"],
            "unit": "eV",
            "mode": "fixed",
        },
        "peak_hash": peak_hash,
    }
    group_material = {
        "source": "MassBank:Antwerp_Univ",
        "accession": accession,
        "connectivity": inchikey[:14],
        "ce": parsed_ce["value"],
    }
    complete_mass_range = mass_range[0] <= 184.0733 and mass_range[1] >= precursor_mz
    return {
        "spectrum_id": f"MassBank:{accession}",
        "evidence_stage": "raw_msms_spectrum",
        "experimental_status": "public_experimental_spectrum",
        "reference_identity_status": "repository_reported_metabolomics_standard",
        "query_confirmation_status": "not_applicable",
        "source": {
            "source_id": "massbank_api_antwerp_metox_2022",
            "repository": "MassBank Europe",
            "repository_record": accession,
            "url": f"{API_ROOT}/{accession}",
            "license": str(record.get("license") or ""),
            "record_created": record.get("date", {}).get("created"),
            "record_updated": record.get("date", {}).get("updated"),
            "record_sha256": _sha256(raw),
        },
        "compound": {
            "name": name,
            "synonyms": names[1:],
            "inchi": compound.get("inchi"),
            "inchikey": inchikey,
            "connectivity": inchikey[:14],
            "smiles": compound.get("smiles"),
            "formula": compound.get("formula"),
            "chemical_class": expected_class,
            "class_provenance": {
                "system": "ChemOnt",
                "raw": links.get("ChemOnt"),
                "source_record": accession,
            },
        },
        "acquisition": {
            "instrument_vendor": "Agilent Technologies",
            "instrument_family": "6500 Series Q-TOF",
            "instrument_model": EXPECTED_INSTRUMENT,
            "instrument_type": acquisition.get("instrument_type"),
            "analyzer": "QTOF",
            "activation_raw": acquisition_tags.get("FRAGMENTATION_MODE"),
            "activation": "beam_cid",
            "collision_cell_region": "axial_hexapole",
            "instrument_domain_id": "agilent_6560_qtof_beam_cid_fixed_ev_reported",
            "collision_energy": {
                "raw": parsed_ce["raw"],
                "value": parsed_ce["value"],
                "unit": "eV",
                "mode": "fixed",
                "minimum": parsed_ce["minimum"],
                "maximum": parsed_ce["maximum"],
                "values": parsed_ce["values"],
            },
            "adduct": EXPECTED_ADDUCT,
            "polarity": "positive",
            "charge": 1,
            "precursor_mz": precursor_mz,
            "resolution_reported": acquisition_tags.get("RESOLUTION"),
            "mass_range_mz": [mass_range[0], mass_range[1]],
        },
        "independence": {
            "independent_fragmentation_group_id": _sha256(
                _canonical_bytes(group_material)
            )[:24],
            "peak_only_mirror_hash": peak_hash,
            "acquisition_aware_hash": _sha256(_canonical_bytes(acquisition_identity)),
            "technical_replicate_id": accession,
            "source_family": "MassBank:Antwerp_Univ",
        },
        "quality": {
            "complete_mass_range": complete_mass_range,
            "centroided": True,
            "sufficient_signal": max(peak[1] for peak in peaks) > 0,
            "precursor_isolated": acquisition_ms.get("ms_type") == "MS2",
            "collision_energy_verified": parsed_ce["value"] in {10.0, 20.0, 40.0},
            "instrument_domain_verified": True,
        },
        "peaks": peaks,
    }


def build_fixture() -> dict[str, Any]:
    spectra: list[dict[str, Any]] = []
    for chemical_class, accessions in CLASS_ACCESSIONS.items():
        for accession in accessions:
            record, raw = _fetch(accession)
            if record.get("accession") != accession:
                raise ValueError(f"API returned the wrong record for {accession}")
            converted = _convert_record(record, raw, chemical_class)
            if converted["source"]["license"] != "CC BY":
                raise ValueError(f"{accession}: expected CC BY license")
            spectra.append(converted)
    spectra.sort(key=lambda item: item["spectrum_id"])
    identities = {item["compound"]["connectivity"] for item in spectra}
    formulas = {item["compound"]["formula"] for item in spectra}
    fixture = {
        "schema_version": "1.0.0",
        "fixture_id": "massbank_antwerp_phospholipids_v1",
        "generated_on": FETCH_DATE,
        "builder": "scripts/build_feature_fragment_massbank_fixture.py",
        "purpose": "Public experimental pilot for class-conditional feature-fragment engineering validation; not cross-platform or exact-identity validation.",
        "record_count": len(spectra),
        "independent_2d_identity_count": len(identities),
        "normalized_formula_family_count": len(formulas),
        "spectra": spectra,
    }
    fixture["content_sha256"] = _sha256(_canonical_bytes(fixture))
    return fixture


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    fixture = build_fixture()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(fixture, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(
        json.dumps(
            {
                "output": str(args.output),
                "record_count": fixture["record_count"],
                "identity_count": fixture["independent_2d_identity_count"],
                "formula_family_count": fixture["normalized_formula_family_count"],
                "content_sha256": fixture["content_sha256"],
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
