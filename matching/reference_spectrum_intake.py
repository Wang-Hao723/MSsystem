"""Auditable intake for raw reference spectra and digitized spectrum images.

Raw MGF, MSP, and centroid mzML peak lists can become experimental-library records
only when compound identity and acquisition provenance are supplied.  A plot
image can be retained with a manually digitized peak list, but its intensities
remain ``digitized_image`` evidence and therefore cannot enter ordinary
spectral search through the existing library quality gate.
"""

from __future__ import annotations

import hashlib
import json
import math
import re
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from utils.collision_energy import parse_collision_energy


MANIFEST_SCHEMA_VERSION = "mssystem-reference-intake/1.0"
RAW_SPECTRUM_EXTENSIONS = {".mgf", ".msp", ".mzml", ".mzxml", ".xml"}
IMAGE_EXTENSIONS = {".png", ".jpg", ".jpeg", ".tif", ".tiff"}
REFERENCE_INTAKE_EXTENSIONS = (
    RAW_SPECTRUM_EXTENSIONS | IMAGE_EXTENSIONS | {".json"}
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _original_name(path: Path) -> str:
    match = re.fullmatch(r"[0-9a-fA-F]{8}_(.+)", path.name)
    return match.group(1) if match else path.name


def _clean(value: Any) -> str:
    return str(value or "").strip()


def _number(value: Any) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def _metadata_value(metadata: Mapping[str, Any], *aliases: str) -> Any:
    wanted = {re.sub(r"[^a-z0-9]", "", alias.lower()) for alias in aliases}
    for key, value in metadata.items():
        normalized = re.sub(r"[^a-z0-9]", "", str(key or "").lower())
        if normalized in wanted and value is not None and _clean(value):
            return value
    return None


def _peaks(value: Any) -> list[list[float]]:
    output: list[list[float]] = []
    if isinstance(value, str):
        rows: Sequence[Any] = value.replace(";", "\n").splitlines()
    elif isinstance(value, Sequence) and not isinstance(value, (bytes, bytearray)):
        rows = value
    else:
        rows = []
    for row in rows:
        if isinstance(row, str):
            parts = re.split(r"[\s,\t]+", row.strip())
        elif isinstance(row, Sequence):
            parts = list(row)
        else:
            continue
        if len(parts) < 2:
            continue
        mz, intensity = _number(parts[0]), _number(parts[1])
        if mz is not None and intensity is not None and mz > 0 and intensity > 0:
            output.append([mz, intensity])
    output.sort(key=lambda peak: peak[0])
    return output


def _manifest_payload(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError("reference intake manifest must be a JSON object")
    if payload.get("schema_version") != MANIFEST_SCHEMA_VERSION:
        raise ValueError(
            f"reference intake manifest requires schema_version={MANIFEST_SCHEMA_VERSION}"
        )
    if not isinstance(payload.get("spectra"), list):
        raise ValueError("reference intake manifest requires a spectra array")
    return payload


def _flatten_manifest_context(
    payload: Mapping[str, Any], entry: Mapping[str, Any]
) -> dict[str, Any]:
    dataset = payload.get("dataset") if isinstance(payload.get("dataset"), Mapping) else {}
    defaults = payload.get("defaults") if isinstance(payload.get("defaults"), Mapping) else {}
    merged = {**dict(dataset), **dict(defaults), **dict(entry)}
    compound = merged.get("compound") if isinstance(merged.get("compound"), Mapping) else {}
    instrument = merged.get("instrument_metadata")
    instrument = instrument if isinstance(instrument, Mapping) else {}
    standard = merged.get("reference_standard")
    standard = standard if isinstance(standard, Mapping) else {}
    provenance = merged.get("provenance")
    provenance = provenance if isinstance(provenance, Mapping) else {}
    collision = merged.get("collision_energy")
    collision = collision if isinstance(collision, Mapping) else {}
    return {
        **merged,
        "name": _clean(
            merged.get("name")
            or merged.get("compound_name")
            or compound.get("name")
        ),
        "formula": _clean(merged.get("formula") or compound.get("formula")),
        "smiles": _clean(merged.get("smiles") or compound.get("smiles")),
        "inchi": _clean(merged.get("inchi") or compound.get("inchi")),
        "reported_inchikey": _clean(
            merged.get("reported_inchikey")
            or merged.get("inchikey")
            or compound.get("inchikey")
        ),
        "instrument_type": _clean(
            merged.get("instrument_type") or instrument.get("analyzer")
        ),
        "instrument": _clean(
            merged.get("instrument") or instrument.get("model")
        ),
        "activation_method": _clean(
            merged.get("activation_method") or instrument.get("activation")
        ),
        "collision_cell_region": _clean(
            merged.get("collision_cell_region")
            or instrument.get("collision_cell_region")
        ),
        "instrument_domain_id": _clean(
            merged.get("instrument_domain_id")
            or instrument.get("instrument_domain_id")
        ),
        "collision_energy_value": (
            merged.get("collision_energy_value")
            if merged.get("collision_energy_value") is not None
            else collision.get("value")
        ),
        "collision_energy_unit": _clean(
            merged.get("collision_energy_unit") or collision.get("unit")
        ),
        "collision_energy_mode": _clean(
            merged.get("collision_energy_mode") or collision.get("mode")
        ),
        "collision_energy_min": (
            merged.get("collision_energy_min")
            if merged.get("collision_energy_min") is not None
            else collision.get("minimum")
        ),
        "collision_energy_max": (
            merged.get("collision_energy_max")
            if merged.get("collision_energy_max") is not None
            else collision.get("maximum")
        ),
        "reference_material_type": _clean(
            merged.get("reference_material_type")
            or standard.get("material_type")
        ),
        "standard_lot": _clean(
            merged.get("standard_lot") or standard.get("lot")
        ),
        "standard_purity": (
            merged.get("standard_purity")
            if merged.get("standard_purity") is not None
            else standard.get("purity_percent")
        ),
        "source_provenance_verified": bool(
            merged.get("source_provenance_verified")
            or provenance.get("verified")
        ),
        "doi": _clean(merged.get("doi") or provenance.get("doi")),
        "source_accession": _clean(
            merged.get("source_accession")
            or provenance.get("source_accession")
        ),
        "dataset_accession": _clean(
            merged.get("dataset_accession")
            or provenance.get("dataset_accession")
        ),
        "acquisition_method": _clean(
            merged.get("acquisition_method")
            or instrument.get("acquisition_method")
        ),
    }


def _manifest_entries(payloads: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    entries: list[dict[str, Any]] = []
    for payload in payloads:
        for raw in payload.get("spectra") or []:
            if not isinstance(raw, Mapping):
                raise ValueError("each manifest spectra item must be an object")
            entry = _flatten_manifest_context(payload, raw)
            entry["_manifest_id"] = _clean(payload.get("manifest_id"))
            entries.append(entry)
    return entries


def _entry_for_raw(
    source_name: str, spectrum_id: Any, entries: Sequence[Mapping[str, Any]]
) -> dict[str, Any]:
    source = Path(source_name).name.casefold()
    matching = [
        dict(entry) for entry in entries
        if Path(_clean(entry.get("source_file"))).name.casefold() == source
        and not _clean(entry.get("image_file"))
    ]
    if not matching:
        return {}
    exact = [
        entry for entry in matching
        if _clean(entry.get("spectrum_id") or entry.get("scan_number"))
        == _clean(spectrum_id)
    ]
    if len(exact) == 1:
        return exact[0]
    unscoped = [
        entry for entry in matching
        if not _clean(entry.get("spectrum_id") or entry.get("scan_number"))
    ]
    return unscoped[0] if len(unscoped) == 1 else {}


def _qualified_source_tier(record: Mapping[str, Any], *, raw_peak_file: bool) -> str:
    if not raw_peak_file:
        return "user_unverified"
    purity = _number(record.get("standard_purity"))
    authentic = (
        _clean(record.get("reference_material_type")).lower()
        in {"authentic_standard", "reference_standard", "certified_reference_material"}
        and purity is not None
        and 0 < purity <= 100
        and bool(_clean(record.get("standard_lot")))
        and bool(_clean(record.get("instrument")))
        and bool(_clean(record.get("acquisition_method")))
        and record.get("source_provenance_verified") is True
    )
    public = bool(
        _clean(record.get("doi"))
        or _clean(record.get("source_accession"))
        or _clean(record.get("dataset_accession"))
    ) and bool(_clean(record.get("instrument")))
    return "experimental" if authentic or public else "user_unverified"


def _auditable_identity_present(record: Mapping[str, Any]) -> bool:
    """Require a readable name plus independently interpretable chemistry.

    A library accession or an unverified InChIKey alone is provenance, not a
    compound identity that this importer may merge into the search library.
    Formula, SMILES, or InChI supplies the minimum chemical descriptor; the
    downstream quality gate validates structure/formula/precursor consistency.
    """

    return bool(_clean(record.get("name"))) and any(
        _clean(record.get(field)) for field in ("formula", "smiles", "inchi")
    )


def _apply_context(
    base: Mapping[str, Any], context: Mapping[str, Any], source_path: Path
) -> dict[str, Any]:
    record = {**dict(base), **{key: value for key, value in context.items() if value not in (None, "")}}
    record["name"] = _clean(record.get("name") or record.get("compound_name"))
    record["compound_name"] = record["name"]
    original_name = _original_name(source_path)
    record["source_file"] = original_name
    record["source_content_sha256"] = _sha256(source_path)
    record["source_asset_sha256"] = record["source_content_sha256"]
    record["dataset_source"] = _clean(
        record.get("dataset_source")
        or record.get("dataset_id")
        or f"user_contrib:{Path(original_name).stem}"
    )[:240]
    nested_collision = record.get("collision_energy")
    if isinstance(nested_collision, Mapping):
        record["collision_energy"] = None
    ce_raw = record.get("collision_energy_raw")
    if ce_raw in (None, ""):
        ce_raw = record.get("collision_energy_value")
    if ce_raw in (None, "") and not isinstance(nested_collision, Mapping):
        ce_raw = nested_collision
    parsed = parse_collision_energy(
        ce_raw, record.get("collision_energy_unit")
    )
    if record.get("collision_energy") is None:
        record["collision_energy"] = parsed.get("value")
    record["collision_energy_raw"] = _clean(
        record.get("collision_energy_raw") or parsed.get("raw")
    )
    record["collision_energy_unit"] = _clean(
        record.get("collision_energy_unit") or parsed.get("unit")
    )
    mode = _clean(record.get("collision_energy_mode"))
    record["collision_energy_mode"] = (
        "fixed" if mode.lower() == "single" else mode
    ) or ("fixed" if parsed.get("mode") == "single" else _clean(parsed.get("mode")))
    record["collision_energy_min"] = (
        _number(record.get("collision_energy_min"))
        if record.get("collision_energy_min") not in (None, "")
        else parsed.get("minimum")
    )
    record["collision_energy_max"] = (
        _number(record.get("collision_energy_max"))
        if record.get("collision_energy_max") not in (None, "")
        else parsed.get("maximum")
    )
    record["charge"] = record.get("charge") or record.get("precursor_charge")
    record["ion_mode"] = _clean(record.get("ion_mode") or record.get("polarity"))
    record["intensity_origin"] = _clean(
        record.get("intensity_origin") or "measured"
    )
    record["source_tier"] = _qualified_source_tier(
        record, raw_peak_file=source_path.suffix.lower() in RAW_SPECTRUM_EXTENSIONS
    )
    if not _auditable_identity_present(record):
        # Keep the record visible in preview, but never allow it into the
        # searchable experimental library. Confirmed import runs in strict
        # mode and rejects the batch with the same explicit reason.
        record["source_tier"] = "user_unverified"
        record.setdefault("reference_intake_flags", []).append(
            "missing_auditable_compound_identity"
        )
    record["source_provenance"] = {
        "manifest_id": _clean(record.get("_manifest_id")),
        "reference_material_type": _clean(record.get("reference_material_type")),
        "standard_lot": _clean(record.get("standard_lot")),
        "standard_purity": _number(record.get("standard_purity")),
        "source_provenance_verified": bool(record.get("source_provenance_verified")),
        "doi": _clean(record.get("doi")),
        "dataset_accession": _clean(record.get("dataset_accession")),
        "source_accession": _clean(record.get("source_accession")),
        "source_asset_sha256": record["source_asset_sha256"],
    }
    return record


def _parse_mgf(path: Path, entries: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    from matchms.importing import load_from_mgf

    records: list[dict[str, Any]] = []
    for index, spectrum in enumerate(load_from_mgf(str(path))):
        precursor = _number(spectrum.get("precursor_mz"))
        peaks = list(zip(spectrum.peaks.mz.tolist(), spectrum.peaks.intensities.tolist()))
        if precursor is None or precursor <= 0 or not peaks:
            continue
        metadata = dict(spectrum.metadata)
        spectrum_id = _metadata_value(
            metadata, "spectrum_id", "spectrumid", "scan", "scan_number", "title"
        ) or index
        ce_raw = _metadata_value(
            metadata, "collision_energy", "collisionenergy", "ce", "nce"
        )
        context = _entry_for_raw(_original_name(path), spectrum_id, entries)
        base = {
            "spectrum_id": _clean(spectrum_id),
            "source_row_index": index,
            "precursor_mz": precursor,
            "peaks": [[float(mz), float(intensity)] for mz, intensity in peaks],
            "name": _clean(_metadata_value(
                metadata, "compound_name", "compoundname", "name"
            )),
            "formula": _clean(_metadata_value(
                metadata, "formula", "molecular_formula", "molecularformula"
            )),
            "smiles": _clean(_metadata_value(
                metadata, "smiles", "canonical_smiles", "isomeric_smiles"
            )),
            "inchi": _clean(_metadata_value(metadata, "inchi")),
            "reported_inchikey": _clean(_metadata_value(
                metadata, "inchikey", "inchi_key"
            )),
            "adduct": _clean(_metadata_value(
                metadata, "adduct", "precursor_type", "precursortype", "iontype"
            )),
            "charge": _metadata_value(metadata, "charge"),
            "ion_mode": _clean(_metadata_value(metadata, "ionmode", "ion_mode")),
            "collision_energy_raw": ce_raw,
            "collision_energy_unit": _clean(_metadata_value(
                metadata, "collision_energy_unit", "ce_unit"
            )),
            "collision_energy_mode": _clean(_metadata_value(
                metadata, "collision_energy_mode", "ce_mode"
            )),
            "instrument_type": _clean(_metadata_value(
                metadata, "instrument_type", "mass_analyzer", "massanalyzer"
            )),
            "instrument": _clean(_metadata_value(
                metadata, "instrument", "instrument_model", "instrumentmodel"
            )),
            "activation_method": _clean(_metadata_value(
                metadata, "activation", "activation_method", "fragmentation_method"
            )),
            "source_accession": _clean(_metadata_value(
                metadata, "source_accession", "accession"
            )),
            "doi": _clean(_metadata_value(metadata, "doi")),
            "acquisition_method": _clean(_metadata_value(
                metadata, "acquisition_method", "method_id", "method"
            )),
            "reference_material_type": _clean(_metadata_value(
                metadata, "reference_material_type", "sample_type"
            )),
            "standard_lot": _clean(_metadata_value(metadata, "standard_lot", "lot")),
            "standard_purity": _metadata_value(
                metadata, "standard_purity", "purity_percent", "purity"
            ),
            "source_provenance_verified": str(_metadata_value(
                metadata, "source_provenance_verified"
            ) or "").strip().lower() in {"1", "true", "yes"},
            "intensity_origin": "measured",
        }
        records.append(_apply_context(base, context, path))
    return records


def _parse_msp(path: Path, entries: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    """Parse MSP through matchms while preserving acquisition provenance."""
    from matchms.importing import load_from_msp

    records: list[dict[str, Any]] = []
    for index, spectrum in enumerate(
        load_from_msp(str(path), metadata_harmonization=False)
    ):
        metadata = dict(spectrum.metadata)
        spectrum_id = _metadata_value(
            metadata,
            "spectrum_id",
            "spectrumid",
            "accession",
            "db#",
            "id",
            "name",
        ) or index
        context = _entry_for_raw(_original_name(path), spectrum_id, entries)
        precursor = _number(_metadata_value(
            metadata,
            "precursor_mz",
            "precursormz",
            "precursor m/z",
            "pepmass",
        ))
        if precursor is None:
            precursor = _number(context.get("precursor_mz"))
        peaks = list(zip(
            spectrum.peaks.mz.tolist(), spectrum.peaks.intensities.tolist()
        ))
        if precursor is None or precursor <= 0 or not peaks:
            continue
        ce_raw = _metadata_value(
            metadata, "collision_energy", "collisionenergy", "ce", "nce"
        )
        base = {
            "spectrum_id": _clean(spectrum_id),
            "source_row_index": index,
            "precursor_mz": precursor,
            "peaks": [[float(mz), float(intensity)] for mz, intensity in peaks],
            "name": _clean(_metadata_value(
                metadata, "compound_name", "compoundname", "name"
            )),
            "formula": _clean(_metadata_value(
                metadata, "formula", "molecular_formula", "molecularformula"
            )),
            "smiles": _clean(_metadata_value(
                metadata, "smiles", "canonical_smiles", "isomeric_smiles"
            )),
            "inchi": _clean(_metadata_value(metadata, "inchi")),
            "reported_inchikey": _clean(_metadata_value(
                metadata, "inchikey", "inchi_key"
            )),
            "adduct": _clean(_metadata_value(
                metadata, "adduct", "precursor_type", "precursortype", "iontype"
            )),
            "charge": _metadata_value(metadata, "charge"),
            "ion_mode": _clean(_metadata_value(
                metadata, "ionmode", "ion_mode", "ion mode"
            )),
            "collision_energy_raw": ce_raw,
            "collision_energy_unit": _clean(_metadata_value(
                metadata, "collision_energy_unit", "ce_unit"
            )),
            "collision_energy_mode": _clean(_metadata_value(
                metadata, "collision_energy_mode", "ce_mode"
            )),
            "instrument_type": _clean(_metadata_value(
                metadata, "instrument_type", "mass_analyzer", "massanalyzer"
            )),
            "instrument": _clean(_metadata_value(
                metadata, "instrument", "instrument_model", "instrumentmodel"
            )),
            "activation_method": _clean(_metadata_value(
                metadata, "activation", "activation_method", "fragmentation_method"
            )),
            "source_accession": _clean(_metadata_value(
                metadata, "source_accession", "accession", "db#", "id"
            )),
            "doi": _clean(_metadata_value(metadata, "doi")),
            "acquisition_method": _clean(_metadata_value(
                metadata, "acquisition_method", "method_id", "method"
            )),
            "reference_material_type": _clean(_metadata_value(
                metadata, "reference_material_type", "sample_type"
            )),
            "standard_lot": _clean(_metadata_value(
                metadata, "standard_lot", "lot"
            )),
            "standard_purity": _metadata_value(
                metadata, "standard_purity", "purity_percent", "purity"
            ),
            "source_provenance_verified": str(_metadata_value(
                metadata, "source_provenance_verified"
            ) or "").strip().lower() in {"1", "true", "yes"},
            "intensity_origin": "measured",
        }
        records.append(_apply_context(base, context, path))
    return records


def _parse_mzml(path: Path, entries: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    if path.suffix.lower() in {".mzxml", ".xml"}:
        raise ValueError("mzXML must be normalized to centroid mzML before reference import")
    from utils.mzml_io import parse_mzml_file

    parsed = parse_mzml_file(path, allow_unknown_representation=False)
    records: list[dict[str, Any]] = []
    for index, spectrum in enumerate(parsed):
        spectrum_id = spectrum.get("spectrum_id", index)
        context = _entry_for_raw(_original_name(path), spectrum_id, entries)
        if not context:
            raise ValueError(
                f"{path.name} spectrum {spectrum_id} lacks compound metadata in the reference manifest"
            )
        base = {
            **spectrum,
            "source_row_index": index,
            "collision_energy_raw": spectrum.get("collision_energy"),
            "intensity_origin": "measured",
        }
        records.append(_apply_context(base, context, path))
    return records


def _image_records(
    entries: Sequence[Mapping[str, Any]], files_by_name: Mapping[str, Path]
) -> tuple[list[dict[str, Any]], set[str]]:
    records: list[dict[str, Any]] = []
    used: set[str] = set()
    for index, context in enumerate(entries):
        image_name = Path(_clean(context.get("image_file"))).name
        embedded_peaks = _peaks(context.get("peaks"))
        if not image_name and not embedded_peaks:
            continue
        if not image_name:
            raise ValueError("embedded peak-list manifest record requires image_file")
        image = files_by_name.get(image_name.casefold())
        if image is None or image.suffix.lower() not in IMAGE_EXTENSIONS:
            raise ValueError(f"manifest image_file is missing from upload batch: {image_name}")
        if len(embedded_peaks) < 3:
            raise ValueError(f"digitized image requires at least three m/z-intensity pairs: {image_name}")
        digitization = context.get("digitization")
        digitization = digitization if isinstance(digitization, Mapping) else {}
        if not bool(digitization.get("reviewer_confirmed")):
            raise ValueError("digitized image peak list requires reviewer_confirmed=true")
        if not _clean(digitization.get("method")):
            raise ValueError("digitized image peak list requires a digitization method")
        base = {
            "spectrum_id": _clean(context.get("spectrum_id")) or f"image-{index}",
            "source_row_index": index,
            "precursor_mz": context.get("precursor_mz"),
            "peaks": embedded_peaks,
            "intensity_origin": "digitized_image",
            "digitization_method": _clean(digitization.get("method")),
            "digitization_reviewer_confirmed": True,
            "source_image_file": _original_name(image),
            "source_image_sha256": _sha256(image),
        }
        record = _apply_context(base, context, image)
        if not _auditable_identity_present(record):
            raise ValueError(
                f"digitized image lacks compound name plus formula/SMILES/InChI: {image_name}"
            )
        if (_number(record.get("precursor_mz")) or 0) <= 0:
            raise ValueError(
                f"digitized image lacks a positive precursor_mz: {image_name}"
            )
        # Image-derived intensity cannot become an ordinary experimental
        # library score without a separate, reviewed promotion workflow.
        record["source_tier"] = "user_unverified"
        records.append(record)
        used.add(_original_name(image).casefold())
    return records, used


def parse_reference_intake(
    file_paths: Sequence[Path], *, strict: bool = False, build_consensus: bool = True
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Parse one mixed raw-spectrum/manifest/image contribution batch."""
    paths = [Path(path) for path in file_paths]
    files_by_name = {_original_name(path).casefold(): path for path in paths}
    manifest_paths = [path for path in paths if path.suffix.lower() == ".json"]
    payloads = [_manifest_payload(path) for path in manifest_paths]
    entries = _manifest_entries(payloads)
    records: list[dict[str, Any]] = []
    stats: list[dict[str, Any]] = []

    image_records, used_images = _image_records(entries, files_by_name)
    records.extend(image_records)
    for path in paths:
        extension = path.suffix.lower()
        if extension == ".json":
            stats.append({
                "file": _original_name(path),
                "total": len(entries),
                "written": len(image_records),
                "format": "reference_manifest",
            })
            continue
        if extension in IMAGE_EXTENSIONS:
            if _original_name(path).casefold() not in used_images:
                message = "spectrum image has no reviewed digitized manifest record"
                stats.append({"file": _original_name(path), "error": message, "total": 0, "written": 0})
                if strict:
                    raise ValueError(f"{_original_name(path)}: {message}")
            continue
        try:
            if extension == ".mgf":
                parsed = _parse_mgf(path, entries)
                format_name = "raw_mgf_reference"
            elif extension == ".msp":
                parsed = _parse_msp(path, entries)
                format_name = "raw_msp_reference"
            elif extension in {".mzml", ".mzxml", ".xml"}:
                parsed = _parse_mzml(path, entries)
                format_name = "raw_mzml_reference"
            else:
                continue
            if strict and not parsed:
                raise ValueError("no valid MS/MS spectra were parsed")
            if strict:
                missing_identity = [
                    str(record.get("spectrum_id") or "unknown")
                    for record in parsed
                    if not _auditable_identity_present(record)
                ]
                if missing_identity:
                    raise ValueError(
                        "raw reference spectrum lacks compound name plus "
                        "formula/SMILES/InChI: " + ", ".join(missing_identity[:5])
                    )
            records.extend(parsed)
            stats.append({
                "file": _original_name(path),
                "total": len(parsed),
                "written": len(parsed),
                "skipped_no_mz": 0,
                "skipped_no_frag": 0,
                "format": format_name,
            })
        except Exception as exc:
            stats.append({
                "file": _original_name(path),
                "error": str(exc),
                "total": 0,
                "written": 0,
            })
            if strict:
                raise
    if build_consensus:
        from matching.reference_consensus import (
            replace_technical_replicates_with_consensus,
        )

        records, consensus = replace_technical_replicates_with_consensus(records)
        stats.append({
            "file": "same-condition authenticated-standard spectra",
            "total": consensus["member_records_replaced"],
            "written": consensus["consensus_groups_built"],
            "format": "technical_replicate_consensus",
            "consensus": consensus,
        })
    return records, stats


__all__ = [
    "IMAGE_EXTENSIONS",
    "MANIFEST_SCHEMA_VERSION",
    "RAW_SPECTRUM_EXTENSIONS",
    "REFERENCE_INTAKE_EXTENSIONS",
    "parse_reference_intake",
]
