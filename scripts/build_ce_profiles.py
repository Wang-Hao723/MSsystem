#!/usr/bin/env python3
"""Build an atomic, provenance-rich multi-source CE summary database.

The output is a compact sidecar.  It does not duplicate the source MSP or the
production spectrum library: only multi-energy profile metadata, per-energy
consensus fragments, and fragment trajectories are retained.  MSP and MGF
inputs from MassBank, MoNA, GNPS, community libraries, and local curated
libraries use the same conservative metadata and identity gates.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import sqlite3
import statistics
import sys
import uuid
from typing import Any, Iterable, Iterator

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from rdkit import Chem, rdBase  # noqa: E402

from matching.rule_quality import normalize_adduct, normalize_polarity  # noqa: E402
from utils.collision_energy import parse_collision_energy  # noqa: E402


_FULL_INCHIKEY = re.compile(r"^[A-Z]{14}-[A-Z]{10}-[A-Z]$")
_PEAK_LINE = re.compile(
    r"^\s*([-+]?(?:\d+(?:\.\d*)?|\.\d+))\s+"
    r"([-+]?(?:\d+(?:\.\d*)?|\.\d+))"
)
_SUPPORTED_INPUT_EXTENSIONS = {".msp", ".mgf"}

SCHEMA_VERSION = "2"


_METADATA_KEY_ALIASES = {
    "name": "name",
    "title": "name",
    "compoundname": "name",
    "inchikey": "inchikey",
    "inchiaux": "inchikey",
    "inchi": "inchi",
    "smiles": "smiles",
    "canonicalsmiles": "smiles",
    "formula": "formula",
    "molecularformula": "formula",
    "precursortype": "precursor_type",
    "adduct": "precursor_type",
    "iontype": "precursor_type",
    "precursormz": "precursor_mz",
    "pepmass": "precursor_mz",
    "parentmass": "precursor_mz",
    "ionmode": "ion_mode",
    "ionizationmode": "ion_mode",
    "polarity": "ion_mode",
    "collisionenergy": "collision_energy",
    "collisionenergyvalue": "collision_energy",
    "activationenergy": "collision_energy",
    "hcdenergy": "collision_energy",
    "ce": "collision_energy",
    "collisionenergyunit": "collision_energy_unit",
    "ceunit": "collision_energy_unit",
    "instrumenttype": "instrument_type",
    "massanalyzertype": "instrument_type",
    "instrument": "instrument",
    "instrumentmodel": "instrument",
    "spectrumtype": "spectrum_type",
    "mslevel": "ms_level",
    "numpeaks": "num_peaks",
    "db": "source_accession",
    "accession": "source_accession",
    "spectrumid": "source_accession",
    "usi": "source_accession",
    "scans": "source_accession",
    "dataset": "source_dataset",
    "datasetsource": "source_dataset",
    "library": "source_dataset",
    "libraryname": "source_dataset",
    "source": "source_dataset",
}


def _canonical_metadata_key(value: Any) -> str:
    compact = re.sub(r"[^a-z0-9]+", "", str(value or "").casefold())
    return _METADATA_KEY_ALIASES.get(compact, compact)


def _set_metadata(
    metadata: dict[str, Any],
    key: Any,
    value: Any,
) -> None:
    canonical = _canonical_metadata_key(key)
    if not canonical:
        return
    text = str(value or "").strip()
    if not text:
        return
    if canonical in {"comments", "synon", "synonym"}:
        metadata.setdefault("comments", []).append(text)
        return
    metadata.setdefault(canonical, text)


def _comment_metadata(comment: str) -> dict[str, str]:
    """Extract only CE/identity/provenance fields from heterogeneous comments."""
    output: dict[str, str] = {}
    aliases = "|".join(
        sorted(
            (
                "collision[_ ]?energy",
                "collisionenergy",
                "activation[_ ]?energy",
                "hcd[_ ]?energy",
                "ce",
                "collision[_ ]?energy[_ ]?unit",
                "ce[_ ]?unit",
                "inchikey",
                "inchiaux",
                "inchi",
                "smiles",
                "formula",
                "adduct",
                "precursor[_ ]?type",
                "ion[_ ]?mode",
                "polarity",
                "instrument[_ ]?type",
                "instrument",
                "spectrum[_ ]?id",
                "accession",
                "dataset[_ ]?source",
                "library[_ ]?name",
            ),
            key=len,
            reverse=True,
        )
    )
    pattern = re.compile(
        rf"(?i)(?:^|[;\t ])(?P<key>{aliases})\s*[:=]\s*"
        r'(?P<value>"[^"]*"|[^;\t]+?)(?='
        rf"\s+(?:{aliases})\s*[:=]|[;\t]|$)"
    )
    for match in pattern.finditer(comment or ""):
        value = match.group("value").strip().strip('"')
        if value:
            output[_canonical_metadata_key(match.group("key"))] = value
    return output


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def iter_msp_records(path: Path) -> Iterator[dict[str, Any]]:
    """Stream MSP records without loading a large file at once."""
    metadata: dict[str, Any] = {}
    peaks: list[list[float]] = []

    def emit() -> dict[str, Any] | None:
        if not metadata and not peaks:
            return None
        return {"metadata": dict(metadata), "peaks": list(peaks)}

    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\r\n")
            if not line:
                record = emit()
                if record is not None:
                    yield record
                metadata.clear()
                peaks.clear()
                continue

            peak_match = _PEAK_LINE.match(line)
            if peak_match and ("num_peaks" in metadata or peaks):
                try:
                    peaks.append(
                        [float(peak_match.group(1)), float(peak_match.group(2))]
                    )
                except (TypeError, ValueError):
                    pass
                continue

            if ":" not in line:
                continue
            key, value = line.split(":", 1)
            canonical = _canonical_metadata_key(key)
            if canonical == "comments":
                metadata.setdefault("comments", []).append(value.strip())
                for parsed_key, parsed_value in _comment_metadata(value).items():
                    metadata.setdefault(parsed_key, parsed_value)
            else:
                _set_metadata(metadata, key, value)

    record = emit()
    if record is not None:
        yield record


def iter_mgf_records(path: Path) -> Iterator[dict[str, Any]]:
    """Stream MGF records and preserve arbitrary source CE metadata."""
    metadata: dict[str, Any] = {}
    peaks: list[list[float]] = []
    inside = False

    def emit() -> dict[str, Any] | None:
        if not metadata and not peaks:
            return None
        return {"metadata": dict(metadata), "peaks": list(peaks)}

    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            upper = line.upper()
            if upper == "BEGIN IONS":
                metadata.clear()
                peaks.clear()
                inside = True
                continue
            if upper == "END IONS":
                record = emit()
                if record is not None:
                    yield record
                metadata.clear()
                peaks.clear()
                inside = False
                continue
            if not inside or not line:
                continue
            peak_match = _PEAK_LINE.match(line)
            if peak_match:
                try:
                    peaks.append(
                        [float(peak_match.group(1)), float(peak_match.group(2))]
                    )
                except (TypeError, ValueError):
                    pass
                continue
            separator = "=" if "=" in line else ":" if ":" in line else ""
            if not separator:
                continue
            key, value = line.split(separator, 1)
            canonical = _canonical_metadata_key(key)
            if canonical == "comments":
                metadata.setdefault("comments", []).append(value.strip())
                for parsed_key, parsed_value in _comment_metadata(value).items():
                    metadata.setdefault(parsed_key, parsed_value)
            else:
                _set_metadata(metadata, key, value)

    record = emit()
    if inside and record is not None:
        yield record


def iter_source_records(path: Path) -> Iterator[dict[str, Any]]:
    suffix = path.suffix.casefold()
    if suffix == ".msp":
        yield from iter_msp_records(path)
    elif suffix == ".mgf":
        yield from iter_mgf_records(path)
    else:
        raise ValueError(f"Unsupported CE source format: {path}")


def resolve_input_files(inputs: Iterable[Path]) -> list[Path]:
    """Resolve files/directories deterministically and remove path duplicates."""
    resolved: dict[str, Path] = {}
    for candidate in inputs:
        path = Path(candidate).expanduser().resolve()
        if path.is_dir():
            children = (
                child
                for child in path.rglob("*")
                if child.is_file()
                and child.suffix.casefold() in _SUPPORTED_INPUT_EXTENSIONS
            )
        elif path.is_file():
            if path.suffix.casefold() not in _SUPPORTED_INPUT_EXTENSIONS:
                raise ValueError(f"Unsupported CE source format: {path}")
            children = (path,)
        else:
            raise FileNotFoundError(path)
        for child in children:
            resolved[str(child.resolve())] = child.resolve()
    return [resolved[key] for key in sorted(resolved)]


def _finite_float(value: Any) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def _validated_inchikey(
    metadata: dict[str, Any],
    cache: dict[tuple[str, str, str], tuple[str | None, str]],
    *,
    allow_reported_inchikey: bool,
) -> tuple[str | None, str]:
    reported = str(metadata.get("inchikey") or "").strip().upper()
    if not _FULL_INCHIKEY.fullmatch(reported):
        return None, "identity_missing_or_invalid"
    inchi = str(metadata.get("inchi") or "").strip()
    smiles = str(metadata.get("smiles") or "").strip()
    key = (reported, inchi, smiles)
    if key in cache:
        return cache[key]
    derived = ""
    try:
        with rdBase.BlockLogs():
            molecule = (
                Chem.MolFromInchi(inchi)
                if inchi
                else Chem.MolFromSmiles(smiles) if smiles else None
            )
            if molecule is not None:
                derived = Chem.MolToInchiKey(molecule)
    except Exception:
        derived = ""
    if derived == reported:
        cache[key] = (reported, "structure_verified")
    elif derived:
        cache[key] = (None, "identity_structure_conflict")
    elif allow_reported_inchikey:
        cache[key] = (reported, "reported_full_inchikey")
    else:
        cache[key] = (None, "identity_not_structure_verified")
    return cache[key]


def _clean_peaks(
    peaks: Iterable[Iterable[Any]],
    precursor_mz: float | None,
    *,
    relative_floor: float = 0.01,
    top_n: int = 50,
) -> list[list[float]]:
    parsed: list[list[float]] = []
    for peak in peaks:
        try:
            mz, intensity = float(peak[0]), float(peak[1])
        except (TypeError, ValueError, IndexError):
            continue
        if not (math.isfinite(mz) and math.isfinite(intensity)):
            continue
        if mz <= 0 or intensity <= 0:
            continue
        if precursor_mz is not None and mz >= precursor_mz - 1.5:
            continue
        parsed.append([mz, intensity])
    if not parsed:
        return []
    maximum = max(item[1] for item in parsed)
    if maximum <= 0:
        return []
    normalized = [
        [mz, intensity / maximum]
        for mz, intensity in parsed
        if intensity / maximum >= relative_floor
    ]
    strongest = sorted(normalized, key=lambda item: item[1], reverse=True)[:top_n]
    return sorted(strongest, key=lambda item: item[0])


def _polarity(metadata: dict[str, Any], adduct: str) -> str:
    explicit = normalize_polarity(
        metadata.get("ion_mode") or metadata.get("ionmode")
    )
    if explicit in {"pos", "neg"}:
        return explicit
    if adduct.endswith("+"):
        return "pos"
    if adduct.endswith("-"):
        return "neg"
    return ""


def _profile_record(
    record: dict[str, Any],
    structure_cache: dict[
        tuple[str, str, str],
        tuple[str | None, str],
    ],
    *,
    source_dataset: str,
    allow_reported_inchikey: bool,
) -> tuple[dict[str, Any] | None, str]:
    metadata = record["metadata"]
    spectrum_type = str(metadata.get("spectrum_type") or "").strip().upper()
    ms_level = str(metadata.get("ms_level") or "").strip().casefold()
    if (
        spectrum_type
        and "MS2" not in spectrum_type
        and "MS/MS" not in spectrum_type
    ) or (ms_level and ms_level not in {"2", "ms2", "ms/ms"}):
        return None, "not_ms2"
    instrument_type = str(metadata.get("instrument_type") or "").strip()
    if re.search(r"(^|[- ])(?:EI|GC)(?:[- ]|$)", instrument_type, re.I):
        return None, "not_lc_esi_msms"

    inchikey, identity_verification = _validated_inchikey(
        metadata,
        structure_cache,
        allow_reported_inchikey=allow_reported_inchikey,
    )
    if not inchikey:
        return None, identity_verification
    adduct = normalize_adduct(
        metadata.get("precursor_type") or metadata.get("adduct")
    )
    polarity = _polarity(metadata, adduct)
    if not adduct or not polarity:
        return None, "missing_adduct_or_polarity"

    ce = parse_collision_energy(
        metadata.get("collision_energy"),
        metadata.get("collision_energy_unit"),
    )
    if ce["mode"] != "single" or ce["value"] is None:
        return None, f"ce_mode_{ce['mode']}"
    if not ce["unit"]:
        return None, "ce_unit_unknown"

    precursor_mz = _finite_float(
        metadata.get("precursormz") or metadata.get("precursor_mz")
    )
    peaks = _clean_peaks(record.get("peaks") or [], precursor_mz)
    if len(peaks) < 3:
        return None, "fewer_than_3_clean_fragments"

    return {
        "inchikey": inchikey,
        "inchikey14": inchikey[:14],
        "identity_verification": identity_verification,
        "compound_name": str(metadata.get("name") or "").strip(),
        "formula": str(metadata.get("formula") or "").strip(),
        "adduct": adduct,
        "ion_mode": polarity,
        "instrument_type": instrument_type or "unknown",
        "instrument_model": str(metadata.get("instrument") or "").strip()
        or "unknown",
        "ce_value": float(ce["value"]),
        "ce_unit": str(ce["unit"]),
        "ce_raw": str(ce["raw"]),
        "precursor_mz": precursor_mz,
        "source_dataset": str(
            metadata.get("source_dataset") or source_dataset
        ).strip(),
        "accession": str(metadata.get("source_accession") or "").strip(),
        "peaks": peaks,
    }, "accepted"


def _default_source_dataset(path: Path) -> str:
    parent = path.parent.name.strip()
    stem = path.stem.strip()
    return f"{parent}/{stem}" if parent else stem


def _spectrum_fingerprint(record: dict[str, Any]) -> str:
    """Identify exact/near-exact redistributed copies without using accession."""
    payload = {
        "inchikey": record["inchikey"],
        "adduct": record["adduct"],
        "ion_mode": record["ion_mode"],
        "instrument_type": record["instrument_type"].casefold(),
        "instrument_model": record["instrument_model"].casefold(),
        "ce_value": round(float(record["ce_value"]), 6),
        "ce_unit": record["ce_unit"],
        "precursor_mz": (
            round(float(record["precursor_mz"]), 5)
            if record["precursor_mz"] is not None
            else None
        ),
        "peaks": [
            [round(float(mz), 4), round(float(intensity), 4)]
            for mz, intensity in record["peaks"]
        ],
    }
    encoded = json.dumps(
        payload,
        ensure_ascii=True,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _mz_tolerance(mz: float) -> float:
    return max(0.01, abs(mz) * 20e-6)


def _cluster_spectral_peaks(
    spectra: list[dict[str, Any]],
) -> list[dict[str, float]]:
    """Create a reproducible median consensus with per-spectrum support."""
    clusters: list[dict[str, Any]] = []
    for spectrum_index, spectrum in enumerate(spectra):
        for mz, intensity in spectrum["peaks"]:
            best: dict[str, Any] | None = None
            best_delta = float("inf")
            for cluster in clusters:
                delta = abs(mz - cluster["center"])
                if delta <= _mz_tolerance(cluster["center"]) and delta < best_delta:
                    best, best_delta = cluster, delta
            if best is None:
                best = {"center": mz, "observations": {}}
                clusters.append(best)
            previous = best["observations"].get(spectrum_index)
            if previous is None or intensity > previous[1]:
                best["observations"][spectrum_index] = (mz, intensity)
            best["center"] = statistics.median(
                observation[0] for observation in best["observations"].values()
            )

    required = max(1, math.ceil(len(spectra) * 0.5))
    consensus: list[dict[str, float]] = []
    for cluster in clusters:
        observations = list(cluster["observations"].values())
        if len(observations) < required:
            continue
        frequency = len(observations) / len(spectra)
        consensus.append(
            {
                "mz": round(statistics.median(item[0] for item in observations), 6),
                "intensity": statistics.median(item[1] for item in observations)
                * frequency,
                "frequency": frequency,
            }
        )
    if not consensus:
        return []
    maximum = max(item["intensity"] for item in consensus)
    for item in consensus:
        item["intensity"] = round(item["intensity"] / maximum, 6)
        item["frequency"] = round(item["frequency"], 6)
    return sorted(consensus, key=lambda item: item["mz"])


def _fragment_trajectories(
    bins: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    tracks: list[dict[str, Any]] = []
    for energy_bin in bins:
        energy = float(energy_bin["ce_value"])
        for peak in energy_bin["consensus"]:
            best: dict[str, Any] | None = None
            best_delta = float("inf")
            for track in tracks:
                delta = abs(float(peak["mz"]) - track["center"])
                if delta <= _mz_tolerance(track["center"]) and delta < best_delta:
                    best, best_delta = track, delta
            if best is None:
                best = {
                    "center": float(peak["mz"]),
                    "intensity": {},
                    "frequency": {},
                    "mzs": [],
                }
                tracks.append(best)
            best["intensity"][energy] = float(peak["intensity"])
            best["frequency"][energy] = float(peak["frequency"])
            best["mzs"].append(float(peak["mz"]))
            best["center"] = statistics.median(best["mzs"])

    output: list[dict[str, Any]] = []
    for track in tracks:
        if len(track["intensity"]) < 2:
            continue
        points = sorted(track["intensity"].items())
        energies = [item[0] for item in points]
        intensities = [item[1] for item in points]
        slope = (intensities[-1] - intensities[0]) / max(
            energies[-1] - energies[0], 1e-9
        )
        span = max(intensities) - min(intensities)
        apex_index = intensities.index(max(intensities))
        if span < 0.15:
            trend = "stable"
        elif 0 < apex_index < len(intensities) - 1:
            trend = "mid_energy_apex"
        elif slope > 0:
            trend = "increasing"
        elif slope < 0:
            trend = "decreasing"
        else:
            trend = "mixed"
        output.append(
            {
                "fragment_mz": round(statistics.median(track["mzs"]), 6),
                "energy_intensity": {
                    str(key): round(value, 6)
                    for key, value in sorted(track["intensity"].items())
                },
                "energy_frequency": {
                    str(key): round(value, 6)
                    for key, value in sorted(track["frequency"].items())
                },
                "n_levels": len(points),
                "trend": trend,
                "apex_energy": energies[apex_index],
            }
        )
    return sorted(output, key=lambda item: item["fragment_mz"])


def _schema(connection: sqlite3.Connection) -> None:
    connection.executescript(
        """
        PRAGMA journal_mode=DELETE;
        PRAGMA synchronous=FULL;
        CREATE TABLE Meta (
            key TEXT PRIMARY KEY,
            value TEXT NOT NULL
        );
        CREATE TABLE CEProfiles (
            id INTEGER PRIMARY KEY,
            inchikey TEXT NOT NULL,
            inchikey14 TEXT NOT NULL,
            identity_verification TEXT NOT NULL,
            structure_verified_spectra INTEGER NOT NULL,
            compound_name TEXT,
            formula TEXT,
            adduct TEXT NOT NULL,
            ion_mode TEXT NOT NULL,
            instrument_type TEXT NOT NULL,
            instrument_model TEXT NOT NULL,
            ce_unit TEXT NOT NULL,
            source_dataset TEXT NOT NULL,
            n_energy_levels INTEGER NOT NULL,
            n_spectra INTEGER NOT NULL,
            ce_min REAL NOT NULL,
            ce_max REAL NOT NULL,
            quality_grade TEXT NOT NULL,
            source_file TEXT NOT NULL,
            source_sha256 TEXT NOT NULL,
            built_at TEXT NOT NULL,
            UNIQUE (
                inchikey, adduct, ion_mode, instrument_type,
                instrument_model, ce_unit, source_dataset
            )
        );
        CREATE TABLE SourceFiles (
            source_file TEXT PRIMARY KEY,
            source_dataset TEXT NOT NULL,
            source_format TEXT NOT NULL,
            source_sha256 TEXT NOT NULL,
            source_size_bytes INTEGER NOT NULL,
            source_records INTEGER NOT NULL,
            accepted_records INTEGER NOT NULL,
            duplicate_records INTEGER NOT NULL
        );
        CREATE TABLE ProfileSources (
            profile_id INTEGER NOT NULL,
            source_dataset TEXT NOT NULL,
            source_file TEXT NOT NULL,
            source_sha256 TEXT NOT NULL,
            n_spectra INTEGER NOT NULL,
            FOREIGN KEY(profile_id) REFERENCES CEProfiles(id),
            UNIQUE(profile_id, source_file)
        );
        CREATE TABLE CEBins (
            id INTEGER PRIMARY KEY,
            profile_id INTEGER NOT NULL,
            ce_value REAL NOT NULL,
            n_spectra INTEGER NOT NULL,
            median_precursor_mz REAL,
            consensus_peaks_json TEXT NOT NULL,
            accessions_json TEXT NOT NULL,
            FOREIGN KEY(profile_id) REFERENCES CEProfiles(id),
            UNIQUE(profile_id, ce_value)
        );
        CREATE TABLE CEFragmentTrajectories (
            id INTEGER PRIMARY KEY,
            profile_id INTEGER NOT NULL,
            fragment_mz REAL NOT NULL,
            energy_intensity_json TEXT NOT NULL,
            energy_frequency_json TEXT NOT NULL,
            n_levels INTEGER NOT NULL,
            trend TEXT NOT NULL,
            apex_energy REAL,
            FOREIGN KEY(profile_id) REFERENCES CEProfiles(id)
        );
        CREATE INDEX ix_ce_profiles_lookup
            ON CEProfiles(inchikey, adduct, ion_mode, ce_unit);
        CREATE INDEX ix_ce_profiles_source
            ON CEProfiles(source_dataset);
        CREATE INDEX ix_ce_profiles_connectivity
            ON CEProfiles(inchikey14);
        CREATE INDEX ix_ce_bins_profile_energy
            ON CEBins(profile_id, ce_value);
        CREATE INDEX ix_ce_fragments_profile
            ON CEFragmentTrajectories(profile_id);
        CREATE INDEX ix_ce_profile_sources_profile
            ON ProfileSources(profile_id);
        """
    )


def build_database(
    inputs: list[Path],
    output: Path,
    *,
    overwrite: bool = False,
    min_energy_levels: int = 2,
    dry_run: bool = False,
    allow_reported_inchikey: bool = True,
) -> dict[str, Any]:
    inputs = resolve_input_files(inputs)
    if not inputs:
        raise ValueError("No .msp or .mgf source files were found")
    if output.exists() and not overwrite and not dry_run:
        raise FileExistsError(
            f"{output} exists; pass --overwrite to replace it atomically"
        )
    free = shutil.disk_usage(
        output.parent if output.parent.exists() else PROJECT_ROOT
    ).free
    required = max(2 * 1024**3, sum(path.stat().st_size for path in inputs) * 4)
    if not dry_run and free < required:
        raise OSError(
            f"Insufficient free space: need at least {required / 1024**3:.1f} GiB"
        )

    structure_cache: dict[
        tuple[str, str, str],
        tuple[str | None, str],
    ] = {}
    groups: dict[tuple[str, ...], list[dict[str, Any]]] = defaultdict(list)
    fingerprints: dict[str, dict[str, Any]] = {}
    counters: Counter[str] = Counter()
    source_hashes = {str(path): _sha256(path) for path in inputs}
    source_stats: dict[str, Counter[str]] = {
        str(path): Counter() for path in inputs
    }

    for path in inputs:
        source_file = str(path)
        source_dataset = _default_source_dataset(path)
        for record in iter_source_records(path):
            counters["source_records"] += 1
            source_stats[source_file]["source_records"] += 1
            parsed, status = _profile_record(
                record,
                structure_cache,
                source_dataset=source_dataset,
                allow_reported_inchikey=allow_reported_inchikey,
            )
            counters[status] += 1
            source_stats[source_file][status] += 1
            if parsed is None:
                continue
            parsed["source_file"] = source_file
            parsed["source_sha256"] = source_hashes[source_file]
            parsed["source_files"] = {source_file}
            parsed["source_datasets"] = {parsed["source_dataset"]}
            parsed["source_hashes"] = {parsed["source_sha256"]}
            parsed["source_provenance"] = {
                source_file: {
                    "source_dataset": parsed["source_dataset"],
                    "source_sha256": parsed["source_sha256"],
                }
            }
            parsed["accessions"] = (
                {parsed["accession"]} if parsed["accession"] else set()
            )
            fingerprint = _spectrum_fingerprint(parsed)
            previous = fingerprints.get(fingerprint)
            if previous is not None:
                counters["duplicate_spectra_removed"] += 1
                source_stats[source_file]["duplicate_spectra_removed"] += 1
                previous["source_files"].add(source_file)
                previous["source_datasets"].add(parsed["source_dataset"])
                previous["source_hashes"].add(parsed["source_sha256"])
                previous["source_provenance"].update(
                    parsed["source_provenance"]
                )
                previous["accessions"].update(parsed["accessions"])
                continue
            fingerprints[fingerprint] = parsed
            key = (
                parsed["inchikey"],
                parsed["adduct"],
                parsed["ion_mode"],
                parsed["instrument_type"],
                parsed["instrument_model"],
                parsed["ce_unit"],
                parsed["source_dataset"],
            )
            groups[key].append(parsed)

    eligible: list[tuple[tuple[str, ...], list[dict[str, Any]]]] = []
    for key, records in groups.items():
        levels = {round(float(row["ce_value"]), 6) for row in records}
        if len(levels) >= min_energy_levels:
            eligible.append((key, records))
    counters["comparable_groups"] = len(groups)
    counters["multi_ce_profiles"] = len(eligible)
    counters["multi_ce_identities"] = len({key[0] for key, _ in eligible})

    summary = {
        "schema_version": SCHEMA_VERSION,
        "inputs": [str(path) for path in inputs],
        "source_sha256": source_hashes,
        "counters": dict(sorted(counters.items())),
        "source_files": {
            path: dict(sorted(counts.items()))
            for path, counts in source_stats.items()
        },
        "min_energy_levels": min_energy_levels,
        "allow_reported_inchikey": allow_reported_inchikey,
        "built_at": datetime.now(timezone.utc).isoformat(),
    }
    if dry_run:
        return summary

    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(f".{output.name}.{uuid.uuid4().hex}.tmp")
    connection = sqlite3.connect(str(temporary))
    try:
        _schema(connection)
        for path in inputs:
            stats = source_stats[str(path)]
            connection.execute(
                """
                INSERT INTO SourceFiles (
                    source_file, source_dataset, source_format, source_sha256,
                    source_size_bytes, source_records, accepted_records,
                    duplicate_records
                ) VALUES (?,?,?,?,?,?,?,?)
                """,
                (
                    str(path),
                    _default_source_dataset(path),
                    path.suffix.casefold().lstrip("."),
                    source_hashes[str(path)],
                    path.stat().st_size,
                    int(stats["source_records"]),
                    int(stats["accepted"]),
                    int(stats["duplicate_spectra_removed"]),
                ),
            )
        for key, records in sorted(eligible, key=lambda item: item[0]):
            energy_groups: dict[float, list[dict[str, Any]]] = defaultdict(list)
            for record in records:
                energy_groups[round(float(record["ce_value"]), 6)].append(record)
            bins: list[dict[str, Any]] = []
            for energy, energy_records in sorted(energy_groups.items()):
                bins.append(
                    {
                        "ce_value": energy,
                        "records": energy_records,
                        "consensus": _cluster_spectral_peaks(energy_records),
                    }
                )
            records_per_bin = [len(item["records"]) for item in bins]
            if len(bins) >= 3 and min(records_per_bin) >= 2:
                quality = "A"
            elif len(bins) >= 3:
                quality = "B"
            else:
                quality = "C"
            first = records[0]
            sources = sorted(
                {
                    source
                    for row in records
                    for source in row["source_files"]
                }
            )
            hashes = sorted(
                {
                    source_hash
                    for row in records
                    for source_hash in row["source_hashes"]
                }
            )
            verified_count = sum(
                row["identity_verification"] == "structure_verified"
                for row in records
            )
            if verified_count == len(records):
                identity_verification = "structure_verified"
            elif verified_count:
                identity_verification = "mixed"
            else:
                identity_verification = "reported_full_inchikey"
            connection.execute(
                """
                INSERT INTO CEProfiles (
                    inchikey, inchikey14, identity_verification,
                    structure_verified_spectra, compound_name, formula,
                    adduct, ion_mode, instrument_type, instrument_model,
                    ce_unit, source_dataset, n_energy_levels, n_spectra,
                    ce_min, ce_max, quality_grade, source_file,
                    source_sha256, built_at
                ) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)
                """,
                (
                    first["inchikey"],
                    first["inchikey14"],
                    identity_verification,
                    verified_count,
                    first["compound_name"],
                    first["formula"],
                    first["adduct"],
                    first["ion_mode"],
                    first["instrument_type"],
                    first["instrument_model"],
                    first["ce_unit"],
                    first["source_dataset"],
                    len(bins),
                    len(records),
                    min(energy_groups),
                    max(energy_groups),
                    quality,
                    json.dumps(sources, ensure_ascii=False),
                    json.dumps(hashes),
                    summary["built_at"],
                ),
            )
            profile_id = int(
                connection.execute("SELECT last_insert_rowid()").fetchone()[0]
            )
            source_counts: Counter[str] = Counter()
            source_provenance: dict[str, dict[str, str]] = {}
            for row in records:
                source_counts.update(row["source_files"])
                source_provenance.update(row["source_provenance"])
            for source_file, source_count in sorted(source_counts.items()):
                provenance = source_provenance[source_file]
                connection.execute(
                    """
                    INSERT INTO ProfileSources (
                        profile_id, source_dataset, source_file,
                        source_sha256, n_spectra
                    ) VALUES (?,?,?,?,?)
                    """,
                    (
                        profile_id,
                        provenance["source_dataset"],
                        source_file,
                        provenance["source_sha256"],
                        int(source_count),
                    ),
                )
            for energy_bin in bins:
                precursor_values = [
                    row["precursor_mz"]
                    for row in energy_bin["records"]
                    if row["precursor_mz"] is not None
                ]
                connection.execute(
                    """
                    INSERT INTO CEBins (
                        profile_id, ce_value, n_spectra, median_precursor_mz,
                        consensus_peaks_json, accessions_json
                    ) VALUES (?,?,?,?,?,?)
                    """,
                    (
                        profile_id,
                        energy_bin["ce_value"],
                        len(energy_bin["records"]),
                        (
                            statistics.median(precursor_values)
                            if precursor_values
                            else None
                        ),
                        json.dumps(
                            energy_bin["consensus"],
                            ensure_ascii=False,
                            separators=(",", ":"),
                        ),
                        json.dumps(
                            sorted(
                                {
                                    accession
                                    for row in energy_bin["records"]
                                    for accession in row["accessions"]
                                    if accession
                                }
                            ),
                            ensure_ascii=False,
                            separators=(",", ":"),
                        ),
                    ),
                )
            for trajectory in _fragment_trajectories(bins):
                connection.execute(
                    """
                    INSERT INTO CEFragmentTrajectories (
                        profile_id, fragment_mz, energy_intensity_json,
                        energy_frequency_json, n_levels, trend, apex_energy
                    ) VALUES (?,?,?,?,?,?,?)
                    """,
                    (
                        profile_id,
                        trajectory["fragment_mz"],
                        json.dumps(
                            trajectory["energy_intensity"],
                            separators=(",", ":"),
                        ),
                        json.dumps(
                            trajectory["energy_frequency"],
                            separators=(",", ":"),
                        ),
                        trajectory["n_levels"],
                        trajectory["trend"],
                        trajectory["apex_energy"],
                    ),
                )

        meta = {
            "schema_version": SCHEMA_VERSION,
            "summary_json": json.dumps(summary, ensure_ascii=False),
            "profile_count": str(len(eligible)),
            "identity_count": str(counters["multi_ce_identities"]),
        }
        connection.executemany(
            "INSERT INTO Meta(key, value) VALUES (?,?)", meta.items()
        )
        connection.commit()
        check = connection.execute("PRAGMA quick_check").fetchone()[0]
        if check != "ok":
            raise RuntimeError(f"sidecar quick_check failed: {check}")
        connection.close()
        connection = None
        os.replace(temporary, output)
    except Exception:
        if connection is not None:
            connection.close()
        temporary.unlink(missing_ok=True)
        raise
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        action="append",
        required=True,
        type=Path,
        help=(
            "MSP/MGF file or directory; repeat for multiple source roots. "
            "Directories are searched recursively."
        ),
    )
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--min-energy-levels", type=int, default=2)
    parser.add_argument(
        "--require-structure-verified-inchikey",
        action="store_true",
        help=(
            "Reject full InChIKeys that cannot be independently derived from "
            "a supplied SMILES/InChI. By default they are retained with a "
            "reported_full_inchikey provenance grade."
        ),
    )
    args = parser.parse_args()
    if args.min_energy_levels < 2:
        parser.error("--min-energy-levels must be at least 2")
    summary = build_database(
        [path.resolve() for path in args.input],
        args.output.resolve(),
        overwrite=args.overwrite,
        min_energy_levels=args.min_energy_levels,
        dry_run=args.dry_run,
        allow_reported_inchikey=not args.require_structure_verified_inchikey,
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
