"""Build one auditable reference spectrum from technical replicates.

The builder is intentionally strict. It only combines authenticated-standard
records with the same structure-derived identity, precursor/adduct, ion mode,
instrument, acquisition method, standard lot, and fixed collision energy.
Raw assets remain recoverable through their file hashes and member spectrum
identifiers in ``source_provenance``.
"""

from __future__ import annotations

import hashlib
import json
import math
import re
from collections import defaultdict
from collections.abc import Mapping, Sequence
from copy import deepcopy
from statistics import median
from typing import Any

from matching.library_quality import (
    inchi_metadata,
    normalize_formula_text,
    smiles_metadata,
)
from utils.collision_energy import canonical_collision_energy_unit


CONSENSUS_SCHEMA_VERSION = "mssystem-v15-reference-consensus/1.0"
_STANDARD_MATERIAL_TYPES = {
    "authentic_standard",
    "reference_standard",
    "certified_reference_material",
}
_FULL_INCHIKEY = re.compile(r"^[A-Z]{14}-[A-Z]{10}-[A-Z]$")


def _clean(value: Any) -> str:
    return str(value or "").strip()


def _finite_number(value: Any) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def _normalized_text(value: Any) -> str:
    return re.sub(r"\s+", " ", _clean(value)).casefold()


def _normalized_adduct(value: Any) -> str:
    return re.sub(r"\s+", "", _clean(value)).upper()


def _normalized_ion_mode(value: Any) -> str:
    compact = re.sub(r"[^a-z]", "", _clean(value).casefold())
    if compact in {"pos", "positive", "p"}:
        return "POS"
    if compact in {"neg", "negative", "n"}:
        return "NEG"
    return ""


def _structure_identity(record: Mapping[str, Any]) -> str:
    """Return only a structure-derived full InChIKey."""
    for field, parser in (("smiles", smiles_metadata), ("inchi", inchi_metadata)):
        if not _clean(record.get(field)):
            continue
        metadata = parser(record.get(field))
        key = _clean(metadata.get("inchikey")).upper()
        if metadata.get("status") == "valid" and _FULL_INCHIKEY.fullmatch(key):
            return key
    return ""


def _consensus_group_key(record: Mapping[str, Any]) -> tuple[Any, ...] | None:
    """Return a fail-closed technical-replicate key, or ``None``."""
    if _clean(record.get("source_tier")).casefold() != "experimental":
        return None
    if _clean(record.get("intensity_origin")).casefold() != "measured":
        return None
    if not bool(record.get("source_provenance_verified")):
        return None
    material = _clean(record.get("reference_material_type")).casefold()
    if material not in _STANDARD_MATERIAL_TYPES:
        return None

    structure_key = _structure_identity(record)
    formula = normalize_formula_text(record.get("formula"))
    adduct = _normalized_adduct(record.get("adduct"))
    ion_mode = _normalized_ion_mode(record.get("ion_mode"))
    instrument_type = _normalized_text(record.get("instrument_type"))
    instrument = _normalized_text(record.get("instrument"))
    activation = _normalized_text(record.get("activation_method"))
    method = _normalized_text(record.get("acquisition_method"))
    lot = _normalized_text(record.get("standard_lot"))
    precursor_mz = _finite_number(record.get("precursor_mz"))
    ce_value = _finite_number(record.get("collision_energy"))
    ce_unit = canonical_collision_energy_unit(record.get("collision_energy_unit"))
    ce_mode = _clean(record.get("collision_energy_mode")).casefold()
    charge = _finite_number(record.get("charge"))

    required = (
        structure_key,
        formula,
        adduct,
        ion_mode,
        instrument_type,
        instrument,
        activation,
        method,
        lot,
        ce_unit,
    )
    if not all(required):
        return None
    if precursor_mz is None or precursor_mz <= 0:
        return None
    if ce_value is None or ce_value < 0 or ce_mode not in {"fixed", "single"}:
        return None
    if charge is None or charge == 0:
        return None

    ce_min = _finite_number(record.get("collision_energy_min"))
    ce_max = _finite_number(record.get("collision_energy_max"))
    if ce_min is not None and abs(ce_min - ce_value) > 1e-9:
        return None
    if ce_max is not None and abs(ce_max - ce_value) > 1e-9:
        return None

    return (
        structure_key,
        formula,
        adduct,
        ion_mode,
        int(charge),
        round(precursor_mz, 4),
        instrument_type,
        instrument,
        activation,
        _normalized_text(record.get("collision_cell_region")),
        _normalized_text(record.get("instrument_domain_id")),
        method,
        material,
        lot,
        round(_finite_number(record.get("standard_purity")) or 0.0, 4),
        ce_unit,
        round(ce_value, 6),
    )


def _member_identity(record: Mapping[str, Any], index: int) -> tuple[Any, ...]:
    """Identify a scan without counting an uploaded duplicate as a replicate."""
    asset_hash = _clean(
        record.get("source_asset_sha256") or record.get("source_content_sha256")
    )
    return (
        asset_hash,
        record.get("source_row_index"),
        _clean(record.get("spectrum_id")),
        index if not asset_hash else None,
    )


def _normalized_peaks(
    record: Mapping[str, Any], mz_tolerance_da: float
) -> list[tuple[float, float]]:
    raw: list[tuple[float, float]] = []
    for peak in record.get("peaks") or []:
        if not isinstance(peak, Sequence) or len(peak) < 2:
            continue
        mz = _finite_number(peak[0])
        intensity = _finite_number(peak[1])
        if mz is not None and intensity is not None and mz > 0 and intensity > 0:
            raw.append((mz, intensity))
    if not raw:
        return []
    maximum = max(intensity for _, intensity in raw)
    scaled = sorted((mz, intensity / maximum * 100.0) for mz, intensity in raw)

    # Retain one contribution per replicate and tolerance window.
    collapsed: list[tuple[float, float]] = []
    for mz, intensity in scaled:
        if collapsed and abs(mz - collapsed[-1][0]) <= mz_tolerance_da:
            if intensity > collapsed[-1][1]:
                collapsed[-1] = (mz, intensity)
        else:
            collapsed.append((mz, intensity))
    return collapsed


def _cluster_replicate_peaks(
    records: Sequence[Mapping[str, Any]], mz_tolerance_da: float
) -> list[dict[str, Any]]:
    points: list[tuple[float, float, int]] = []
    for replicate_index, record in enumerate(records):
        points.extend(
            (mz, intensity, replicate_index)
            for mz, intensity in _normalized_peaks(record, mz_tolerance_da)
        )
    points.sort(key=lambda item: (item[0], item[2]))

    clusters: list[dict[str, Any]] = []
    for mz, intensity, replicate_index in points:
        candidates = [
            (abs(mz - cluster["center"]), cluster_index)
            for cluster_index, cluster in enumerate(clusters)
            if abs(mz - cluster["center"]) <= mz_tolerance_da
        ]
        if candidates:
            _, cluster_index = min(candidates)
            cluster = clusters[cluster_index]
        else:
            cluster = {"center": mz, "members": {}}
            clusters.append(cluster)
        existing = cluster["members"].get(replicate_index)
        if existing is None or intensity > existing[1]:
            cluster["members"][replicate_index] = (mz, intensity)
        cluster["center"] = median(
            value[0] for value in cluster["members"].values()
        )
    return clusters


def _build_consensus_record(
    records: Sequence[Mapping[str, Any]],
    group_key: tuple[Any, ...],
    *,
    mz_tolerance_da: float,
    minimum_peak_prevalence: float,
) -> dict[str, Any] | None:
    replicate_count = len(records)
    required_members = max(2, math.ceil(replicate_count * minimum_peak_prevalence))
    retained: list[dict[str, Any]] = []
    for cluster in _cluster_replicate_peaks(records, mz_tolerance_da):
        members = cluster["members"]
        if len(members) < required_members:
            continue
        retained.append({
            "mz": round(median(value[0] for value in members.values()), 6),
            "intensity": median(value[1] for value in members.values()),
            "replicate_count": len(members),
            "prevalence": len(members) / replicate_count,
        })
    if len(retained) < 3:
        return None

    max_intensity = max(item["intensity"] for item in retained)
    peaks = [
        [item["mz"], round(item["intensity"] / max_intensity * 100.0, 6)]
        for item in sorted(retained, key=lambda item: item["mz"])
    ]
    member_provenance = [
        {
            "source_file": _clean(record.get("source_file")),
            "source_asset_sha256": _clean(
                record.get("source_asset_sha256")
                or record.get("source_content_sha256")
            ),
            "spectrum_id": _clean(record.get("spectrum_id")),
            "source_row_index": record.get("source_row_index"),
        }
        for record in records
    ]
    hash_payload = {
        "schema_version": CONSENSUS_SCHEMA_VERSION,
        "group_key": group_key,
        "members": member_provenance,
        "peaks": peaks,
        "mz_tolerance_da": mz_tolerance_da,
        "minimum_peak_prevalence": minimum_peak_prevalence,
    }
    consensus_hash = hashlib.sha256(
        json.dumps(
            hash_payload,
            ensure_ascii=False,
            sort_keys=True,
            separators=(",", ":"),
            default=str,
        ).encode("utf-8")
    ).hexdigest()

    consensus = deepcopy(dict(records[0]))
    consensus.update({
        "spectrum_id": f"v15_consensus_{consensus_hash[:24]}",
        "source_file": f"V15_consensus_{consensus_hash[:16]}.msp",
        "source_row_index": 0,
        "source_content_sha256": consensus_hash,
        "source_asset_sha256": consensus_hash,
        "dataset_source": f"user_contrib_consensus:{consensus_hash[:24]}",
        "source_accession": f"V15-CONSENSUS-{consensus_hash[:16].upper()}",
        "peaks": peaks,
        "source_tier": "experimental",
        "intensity_origin": "measured_consensus",
        "source_provenance_verified": True,
        "consensus_member_count": replicate_count,
        "technical_replicates_count_as_independent_evidence": False,
        "independent_evidence_count": 1,
    })
    flags = list(consensus.get("reference_intake_flags") or [])
    consensus["reference_intake_flags"] = sorted(
        set(flags + ["technical_replicate_consensus"])
    )
    provenance = dict(consensus.get("source_provenance") or {})
    provenance.update({
        "consensus_schema_version": CONSENSUS_SCHEMA_VERSION,
        "consensus_algorithm": "base_peak_normalize_median_peak_consensus",
        "consensus_member_count": replicate_count,
        "technical_replicates_count_as_independent_evidence": False,
        "independent_evidence_count": 1,
        "mz_tolerance_da": mz_tolerance_da,
        "minimum_peak_prevalence": minimum_peak_prevalence,
        "minimum_peak_members": required_members,
        "members": member_provenance,
        "retained_peaks": [
            {
                "mz": item["mz"],
                "replicate_count": item["replicate_count"],
                "prevalence": round(item["prevalence"], 6),
            }
            for item in retained
        ],
        "consensus_sha256": consensus_hash,
    })
    consensus["source_provenance"] = provenance
    return consensus


def replace_technical_replicates_with_consensus(
    records: Sequence[Mapping[str, Any]],
    *,
    mz_tolerance_da: float = 0.01,
    minimum_replicates: int = 2,
    minimum_peak_prevalence: float = 0.5,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Replace qualifying technical-replicate groups with one consensus."""
    if not math.isfinite(mz_tolerance_da) or not 0 < mz_tolerance_da <= 0.1:
        raise ValueError("mz_tolerance_da must be > 0 and <= 0.1")
    if minimum_replicates < 2:
        raise ValueError("minimum_replicates must be at least 2")
    if not 0 < minimum_peak_prevalence <= 1:
        raise ValueError("minimum_peak_prevalence must be in (0, 1]")

    prepared = [deepcopy(dict(record)) for record in records]
    grouped: dict[tuple[Any, ...], list[tuple[int, dict[str, Any]]]] = defaultdict(list)
    eligible_records = 0
    for index, record in enumerate(prepared):
        key = _consensus_group_key(record)
        if key is not None:
            eligible_records += 1
            grouped[key].append((index, record))

    replacements: dict[int, dict[str, Any]] = {}
    replaced_indices: set[int] = set()
    skipped_insufficient_replicates = 0
    skipped_insufficient_reproducible_peaks = 0
    member_records_replaced = 0
    for key, indexed_records in grouped.items():
        unique: list[tuple[int, dict[str, Any]]] = []
        seen_members: set[tuple[Any, ...]] = set()
        for index, record in indexed_records:
            member_id = _member_identity(record, index)
            if member_id in seen_members:
                continue
            seen_members.add(member_id)
            unique.append((index, record))
        if len(unique) < minimum_replicates:
            skipped_insufficient_replicates += 1
            continue
        consensus = _build_consensus_record(
            [record for _, record in unique],
            key,
            mz_tolerance_da=mz_tolerance_da,
            minimum_peak_prevalence=minimum_peak_prevalence,
        )
        if consensus is None:
            skipped_insufficient_reproducible_peaks += 1
            continue
        first_index = min(index for index, _ in unique)
        replacements[first_index] = consensus
        replaced_indices.update(index for index, _ in unique)
        member_records_replaced += len(unique)

    output: list[dict[str, Any]] = []
    for index, record in enumerate(prepared):
        if index in replacements:
            output.append(replacements[index])
        elif index not in replaced_indices:
            output.append(record)

    summary = {
        "schema_version": CONSENSUS_SCHEMA_VERSION,
        "input_records": len(prepared),
        "eligible_records": eligible_records,
        "eligible_condition_groups": len(grouped),
        "consensus_groups_built": len(replacements),
        "member_records_replaced": member_records_replaced,
        "output_records": len(output),
        "skipped_insufficient_replicates": skipped_insufficient_replicates,
        "skipped_insufficient_reproducible_peaks": (
            skipped_insufficient_reproducible_peaks
        ),
        "mz_tolerance_da": mz_tolerance_da,
        "minimum_replicates": minimum_replicates,
        "minimum_peak_prevalence": minimum_peak_prevalence,
        "technical_replicates_count_as_independent_evidence": False,
    }
    return output, summary


__all__ = [
    "CONSENSUS_SCHEMA_VERSION",
    "replace_technical_replicates_with_consensus",
]
