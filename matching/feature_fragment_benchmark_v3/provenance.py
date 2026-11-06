"""Per-spectrum provenance, license, domain, and independence governance."""

from __future__ import annotations

import hashlib
import json
import math
import re
from collections.abc import Mapping, Sequence
from typing import Any

from matching.feature_fragment_benchmark.mirror import exact_peak_hash

from .parsing import ce_program_has_complete_vector, ce_program_key


_FULL_INCHIKEY = re.compile(r"^[A-Z]{14}-[A-Z]{10}-[A-Z]$")
_DB_IDENTIFIER_NAME = re.compile(
    r"^(?:NCGC[-_A-Z0-9]*|CCMSLIB[A-Z0-9_-]*|HMDB\d+|CHEBI:\d+|"
    r"CID:\d+|DTXSID\d+|MSBNK[-_A-Z0-9]+|LM[A-Z0-9]+)$",
    re.IGNORECASE,
)
_UNREPORTED = re.compile(r"(?:unreported|unknown|unresolved|not[_ -]?reported)", re.I)


class MetadataValidationError(ValueError):
    """Raised when a record lacks minimum audit provenance or is inconsistent."""

    def __init__(self, code: str, detail: str = "") -> None:
        self.code = code
        self.detail = detail
        super().__init__(f"{code}: {detail}" if detail else code)


def _canonical_bytes(value: Any) -> bytes:
    return json.dumps(
        value, ensure_ascii=False, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")


def _sha(value: Any, length: int | None = None) -> str:
    digest = hashlib.sha256(_canonical_bytes(value)).hexdigest()
    return digest[:length] if length else digest


def _present(value: Any) -> bool:
    return value is not None and str(value).strip() != ""


def _normalize_peaks(peaks: Sequence[Sequence[Any]]) -> list[list[float]]:
    cleaned: list[tuple[float, float]] = []
    for peak in peaks:
        if len(peak) < 2:
            continue
        try:
            mz, intensity = float(peak[0]), float(peak[1])
        except (TypeError, ValueError):
            continue
        if math.isfinite(mz) and math.isfinite(intensity) and mz > 0 and intensity > 0:
            cleaned.append((mz, intensity))
    if not cleaned:
        raise MetadataValidationError("positive_peak_list_missing")
    base = max(intensity for _, intensity in cleaned)
    return [
        [round(mz, 6), round(intensity / base, 8)]
        for mz, intensity in sorted(cleaned)
    ]


def license_decision(
    license_label: str, policy: Mapping[str, Any]
) -> dict[str, Any]:
    """Resolve a literal license label to the preregistered analysis tier."""
    decisions = policy.get("label_decisions") or {}
    row = decisions.get(license_label)
    if not isinstance(row, Mapping):
        return {
            "label": license_label,
            "tier": "excluded",
            "eligible_for_primary_analysis": False,
            "eligible_for_sensitivity_analysis": False,
            "reason": "license_missing_or_not_preregistered",
            "license_url": None,
        }
    tier = str(row.get("tier") or "excluded")
    return {
        "label": license_label,
        "tier": tier,
        "eligible_for_primary_analysis": tier == "primary",
        "eligible_for_sensitivity_analysis": tier in {"primary", "sensitivity"},
        "reason": row.get("reason"),
        "license_url": row.get("license_url"),
        "share_alike_review_complete": bool(row.get("share_alike_review_complete")),
        "commercial_use_allowed": row.get("commercial_use_allowed"),
    }


def resolve_exact_instrument_domain(
    acquisition: Mapping[str, Any],
    domains: Sequence[Mapping[str, Any]],
) -> dict[str, Any]:
    """Resolve an exact acquisition domain without unit or vendor conversion."""
    program = acquisition.get("collision_energy") or {}
    candidates: list[Mapping[str, Any]] = []
    for domain in domains:
        if acquisition.get("instrument_vendor") != domain.get("vendor"):
            continue
        if acquisition.get("instrument_model") not in (domain.get("models") or []):
            continue
        if acquisition.get("analyzer") != domain.get("analyzer"):
            continue
        if acquisition.get("activation") not in (domain.get("activation") or []):
            continue
        if acquisition.get("collision_cell_region") not in (
            domain.get("collision_cell_regions") or []
        ):
            continue
        if program.get("unit") not in (domain.get("ce_units") or []):
            continue
        if program.get("mode") not in (domain.get("ce_modes") or []):
            continue
        candidates.append(domain)
    if len(candidates) != 1:
        return {
            "status": "unresolved" if not candidates else "ambiguous",
            "instrument_domain_id": None,
            "candidate_domain_ids": sorted(
                str(row.get("domain_id")) for row in candidates
            ),
            "cross_domain_pooling_allowed": False,
        }
    domain = candidates[0]
    exact_fields = {
        "vendor": acquisition.get("instrument_vendor"),
        "model": acquisition.get("instrument_model"),
        "analyzer": acquisition.get("analyzer"),
        "activation": acquisition.get("activation"),
        "collision_cell_region": acquisition.get("collision_cell_region"),
        "ce_unit": program.get("unit"),
        "ce_mode": program.get("mode"),
        "polarity": acquisition.get("polarity"),
        "adduct": acquisition.get("adduct"),
        "charge": acquisition.get("charge"),
    }
    complete = all(_present(value) and not _UNREPORTED.search(str(value)) for value in exact_fields.values())
    return {
        "status": "resolved",
        "instrument_domain_id": domain["domain_id"],
        "family_id": domain.get("family_id"),
        "exact_condition_metadata_complete": complete,
        "exact_fields": exact_fields,
        "cross_domain_pooling_allowed": False,
    }


def make_acquisition_replication_group(metadata: Mapping[str, Any]) -> str:
    """Create the true acquisition replication group without peak content.

    Peak hashes belong only to mirror detection.  When scan, injection, or raw
    file identifiers are missing, all records sharing the same source,
    identity, domain, CE program, polarity, and adduct collapse into one group.
    This loses apparent sample size rather than creating pseudoreplicates.
    """
    source = metadata.get("source") or {}
    compound = metadata.get("compound") or {}
    acquisition = metadata.get("acquisition") or {}
    raw = metadata.get("raw_acquisition") or {}
    material = {
        "repository": source.get("repository"),
        "origin_dataset_id": source.get("origin_dataset_id"),
        "origin_study_id": source.get("origin_study_id"),
        "source_id": source.get("source_id"),
        "laboratory_id": source.get("laboratory_id"),
        "identity": compound.get("connectivity"),
        "instrument_domain_id": acquisition.get("instrument_domain_id"),
        "collision_energy_program": ce_program_key(
            acquisition.get("collision_energy") or {}
        ),
        "polarity": acquisition.get("polarity"),
        "adduct": acquisition.get("adduct"),
        "charge": acquisition.get("charge"),
        "raw_file_id": raw.get("raw_file_id"),
        "injection_id": raw.get("injection_id"),
        "acquisition_event_id": raw.get("acquisition_event_id"),
        "scan_id": raw.get("scan_id"),
    }
    return "acqrep:" + _sha(material, 24)


def make_fragment_pathway_group(
    *, chemical_class: str, stratum_id: str, identity_presence: Sequence[str]
) -> str:
    """Create a rule-level fragment-pathway group from discovery data only."""
    return "fragpath:" + _sha(
        {
            "chemical_class": chemical_class,
            "stratum_id": stratum_id,
            "identity_presence": sorted(set(identity_presence)),
        },
        24,
    )


def _missing_fields(record: Mapping[str, Any]) -> dict[str, bool]:
    source = record["source"]
    compound = record["compound"]
    acquisition = record["acquisition"]
    ce = acquisition["collision_energy"]
    scan = acquisition["scan_range"]
    isolation = acquisition["precursor_isolation"]
    fields = {
        "repository": _present(source.get("repository")),
        "repository_release": _present(source.get("repository_release")),
        "accession": _present(source.get("accession")),
        "record_url": _present(source.get("record_url")),
        "usi": _present(source.get("usi")),
        "origin_dataset_id": _present(source.get("origin_dataset_id")),
        "origin_study_id": _present(source.get("origin_study_id")),
        "source_id": _present(source.get("source_id")),
        "laboratory_id": _present(source.get("laboratory_id")),
        "license_label": _present(source.get("license", {}).get("label")),
        "license_url": _present(source.get("license", {}).get("license_url")),
        "attribution": _present(source.get("license", {}).get("attribution")),
        "rights_holder": _present(
            source.get("license", {}).get("rights_holder")
        ),
        "record_sha256": _present(source.get("record_canonical_json_sha256")),
        "source_file_sha256": _present(source.get("source_file_sha256")),
        "stable_compound_name": _present(compound.get("display_name")),
        "full_inchikey": bool(_FULL_INCHIKEY.fullmatch(str(compound.get("inchikey") or ""))),
        "formula": _present(compound.get("formula")),
        "identity_basis": _present(compound.get("identity_basis")),
        "vendor": _present(acquisition.get("instrument_vendor")),
        "model": _present(acquisition.get("instrument_model")),
        "analyzer": _present(acquisition.get("analyzer")),
        "activation": _present(acquisition.get("activation")),
        "collision_cell_region": _present(acquisition.get("collision_cell_region"))
        and not _UNREPORTED.search(str(acquisition.get("collision_cell_region"))),
        "collision_energy_program": ce.get("status") == "resolved"
        and ce_program_has_complete_vector(ce),
        "scan_range": scan.get("status") == "resolved",
        "polarity": acquisition.get("polarity") in {"positive", "negative"},
        "adduct": _present(acquisition.get("adduct")),
        "charge": isinstance(acquisition.get("charge"), int),
        "precursor_mz": isinstance(acquisition.get("precursor_mz"), (int, float)),
        "precursor_isolation": isolation.get("status") == "verified",
        "raw_file_id": _present(record.get("raw_acquisition", {}).get("raw_file_id")),
        "scan_id": _present(record.get("raw_acquisition", {}).get("scan_id")),
        "injection_id": _present(record.get("raw_acquisition", {}).get("injection_id")),
    }
    return fields


def normalize_public_spectrum(
    record: Mapping[str, Any],
    *,
    license_policy: Mapping[str, Any],
    instrument_domains: Sequence[Mapping[str, Any]],
) -> dict[str, Any]:
    """Normalize one adapter-produced public spectrum and apply failure gates."""
    spectrum_id = str(record.get("spectrum_id") or "").strip()
    if not spectrum_id:
        raise MetadataValidationError("spectrum_id_missing")
    source = dict(record.get("source") or {})
    if not all(_present(source.get(field)) for field in ("repository", "accession", "record_url", "source_id")):
        raise MetadataValidationError("minimum_source_trace_missing", spectrum_id)
    compound = dict(record.get("compound") or {})
    display_name = str(compound.get("display_name") or "").strip()
    if not display_name or _DB_IDENTIFIER_NAME.fullmatch(display_name):
        raise MetadataValidationError("stable_compound_name_missing", display_name)
    inchikey = str(compound.get("inchikey") or "")
    if not _FULL_INCHIKEY.fullmatch(inchikey):
        raise MetadataValidationError("full_inchikey_missing", inchikey)
    compound["connectivity"] = inchikey[:14]
    compound.setdefault("identity_basis", "repository_reported_structure")
    compound.setdefault("reference_material_claim", "not_reported")
    compound.setdefault("standard_confirmation_status", "not_confirmed_by_mssystem")

    peaks = _normalize_peaks(record.get("peaks") or [])
    acquisition = dict(record.get("acquisition") or {})
    domain = resolve_exact_instrument_domain(acquisition, instrument_domains)
    acquisition["instrument_domain_resolution"] = domain
    acquisition["instrument_domain_id"] = domain.get("instrument_domain_id")
    ce = acquisition.get("collision_energy") or {}
    if ce.get("conversion_applied") is not False:
        raise MetadataValidationError("collision_energy_conversion_forbidden")

    decision = license_decision(str(source.get("license_label") or ""), license_policy)
    license_data = {
        **decision,
        "attribution": source.pop("attribution", None),
        "rights_holder": source.pop("rights_holder", None),
        "record_literal_preserved": True,
    }
    source["license"] = license_data
    source.pop("license_label", None)

    normalized: dict[str, Any] = {
        "schema_version": "3.0.0",
        "spectrum_id": spectrum_id,
        "evidence_stage": "raw_msms_spectrum",
        "identification_scope": {
            "raw_ion_or_feature": "not_applicable_reference_spectrum",
            "ms1_candidate": False,
            "msms_library_candidate": False,
            "repository_reported_reference_identity": True,
            "isomer_group_resolved": False,
            "standard_confirmed_by_mssystem": False,
        },
        "source": source,
        "compound": compound,
        "acquisition": acquisition,
        "raw_acquisition": dict(record.get("raw_acquisition") or {}),
        "peaks": peaks,
    }
    peak_hash = exact_peak_hash(peaks)
    acquisition_group = make_acquisition_replication_group(normalized)
    scan = acquisition.get("scan_range") or {}
    precursor = acquisition.get("precursor_mz")
    range_contains_precursor = (
        scan.get("status") == "resolved"
        and isinstance(precursor, (int, float))
        and scan["minimum_mz"] <= precursor <= scan["maximum_mz"]
    )
    fields = _missing_fields(normalized)
    missing = sorted(field for field, present in fields.items() if not present)
    isolation = acquisition.get("precursor_isolation") or {}
    verified_single_compound = bool(
        compound.get("reference_material_claim") == "verified_single_compound"
        and compound.get("material_type")
        and compound.get("purity")
        and compound.get("batch_id")
    )
    precursor_selection_verified = bool(
        isolation.get("status") == "verified" or verified_single_compound
    )
    dia_or_in_source = ce.get("mode") in {"MSE", "bbCID", "all_ions"} or acquisition.get(
        "activation_location"
    ) == "in_source"
    common_quantitative = bool(
        decision["eligible_for_primary_analysis"]
        and domain.get("status") == "resolved"
        and domain.get("exact_condition_metadata_complete")
        and ce.get("status") == "resolved"
        and ce_program_has_complete_vector(ce)
        and fields["polarity"]
        and fields["adduct"]
        and fields["charge"]
        and fields["precursor_mz"]
        and range_contains_precursor
    )
    positive_eligible = common_quantitative and precursor_selection_verified
    negative_eligible = positive_eligible and isolation.get("status") == "verified" and not dia_or_in_source
    normalized["independence"] = {
        "two_dimensional_identity_id": compound["connectivity"],
        "formula_family_id": compound.get("formula"),
        "source_origin_id": source.get("origin_dataset_id") or source.get("source_id"),
        "laboratory_id": source.get("laboratory_id"),
        "exact_peak_hash": peak_hash,
        "mirror_cluster_id": "exact:" + peak_hash,
        "acquisition_replication_group_id": acquisition_group,
        "peak_hash_used_in_acquisition_replication_group": False,
        "fragment_pathway_group_id": None,
    }
    normalized["metadata_completeness"] = {
        "fields": fields,
        "missing_fields": missing,
        "complete_field_count": sum(fields.values()),
        "field_count": len(fields),
    }
    normalized["quantitative_eligibility"] = {
        "primary_license_eligible": decision["eligible_for_primary_analysis"],
        "range_contains_precursor": range_contains_precursor,
        "precursor_selection_verified": precursor_selection_verified,
        "verified_single_compound_reference_acquisition": verified_single_compound,
        "positive_lr_eligible": positive_eligible,
        "negative_lr_eligible": negative_eligible,
        "positive_only_composite_or_in_source": bool(dia_or_in_source and positive_eligible),
        "ineligible_reason_codes": sorted(
            code
            for code, failed in {
                "license_not_primary": not decision["eligible_for_primary_analysis"],
                "exact_domain_unresolved": domain.get("status") != "resolved",
                "exact_condition_metadata_incomplete": not domain.get(
                    "exact_condition_metadata_complete", False
                ),
                "collision_energy_unresolved": ce.get("status") != "resolved",
                "collision_energy_vector_incomplete": not ce_program_has_complete_vector(ce),
                "scan_range_missing_or_excludes_precursor": not range_contains_precursor,
                "precursor_selection_unverified": not precursor_selection_verified,
                "precursor_isolation_unverified_for_negative_evidence": isolation.get("status")
                != "verified",
                "composite_or_in_source_negative_evidence_forbidden": dia_or_in_source,
            }.items()
            if failed
        ),
    }
    return normalized
