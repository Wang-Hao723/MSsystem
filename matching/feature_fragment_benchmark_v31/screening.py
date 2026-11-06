"""Dataset-specific public-spectrum intake screening."""

from __future__ import annotations

from collections import Counter, defaultdict
from collections.abc import Mapping, Sequence
from typing import Any


def _present(value: Any) -> bool:
    return value is not None and bool(str(value).strip())


def screen_source_descriptor(descriptor: Mapping[str, Any]) -> dict[str, Any]:
    """Screen a prospective dataset descriptor before downloading spectra."""

    repository = str(descriptor.get("repository") or "")
    if repository == "NIST Tandem Mass Spectral Library" or descriptor.get(
        "commercial_distributor_content"
    ) is True:
        return {
            "status": "excluded",
            "reason_codes": ["commercial_content_outside_public_intake_scope"],
            "quantitative_use_allowed": False,
            "negative_evidence_allowed": False,
        }
    checks = {
        "dataset_specific_origin": _present(descriptor.get("origin_dataset_id")),
        "stable_per_spectrum_accession": descriptor.get(
            "stable_per_spectrum_accession"
        )
        is True,
        "stable_record_url": descriptor.get("stable_record_url") is True,
        "record_or_dataset_license_permits_analysis": descriptor.get(
            "license_permits_analysis"
        )
        is True,
        "reference_identity_with_full_inchikey_and_formula": descriptor.get(
            "reference_identity_complete"
        )
        is True,
        "source_origin_id": _present(descriptor.get("source_origin_id")),
        "laboratory_id": _present(descriptor.get("laboratory_id")),
        "exact_instrument_model": descriptor.get("exact_instrument_model") is True,
        "activation_literal": descriptor.get("activation_literal") is True,
        "native_collision_energy_program": descriptor.get(
            "native_collision_energy_program"
        )
        is True,
        "polarity_adduct_and_charge": descriptor.get(
            "polarity_adduct_and_charge"
        )
        is True,
        "scan_range_contains_precursor": descriptor.get(
            "scan_range_contains_precursor"
        )
        is True,
        "precursor_selection_verified": descriptor.get(
            "precursor_selection_verified"
        )
        is True,
        "raw_acquisition_mapping_or_verified_single_compound_material": bool(
            descriptor.get("raw_acquisition_mapping") is True
            or descriptor.get("verified_single_compound_material") is True
        ),
    }
    missing = sorted(name for name, passed in checks.items() if not passed)
    quantitative = not missing
    negative_checks = {
        "verified_precursor_isolation_window": descriptor.get(
            "verified_precursor_isolation_window"
        )
        is True,
        "not_all_ions_mse_bbcid_or_in_source": descriptor.get(
            "not_composite_or_in_source"
        )
        is True,
        "reported_scan_range_contains_precursor": descriptor.get(
            "scan_range_contains_precursor"
        )
        is True,
    }
    negative_missing = sorted(
        name for name, passed in negative_checks.items() if not passed
    )
    return {
        "status": "quantitative_candidate" if quantitative else "audit_only",
        "checks": checks,
        "missing_requirement_codes": missing,
        "negative_evidence_checks": negative_checks,
        "negative_evidence_missing_requirement_codes": negative_missing,
        "quantitative_use_allowed": quantitative,
        "negative_evidence_allowed": quantitative and not negative_missing,
    }


def _record_descriptor(record: Mapping[str, Any]) -> dict[str, Any]:
    source = record.get("source") or {}
    compound = record.get("compound") or {}
    acquisition = record.get("acquisition") or {}
    raw = record.get("raw_acquisition") or {}
    eligibility = record.get("quantitative_eligibility") or {}
    license_data = source.get("license") or {}
    isolation = acquisition.get("precursor_isolation") or {}
    ce = acquisition.get("collision_energy") or {}
    composite = ce.get("mode") in {"MSE", "bbCID", "all_ions"} or acquisition.get(
        "activation_location"
    ) == "in_source"
    return {
        "repository": source.get("repository"),
        "origin_dataset_id": source.get("origin_dataset_id"),
        "stable_per_spectrum_accession": _present(source.get("accession")),
        "stable_record_url": _present(source.get("record_url")),
        "license_permits_analysis": license_data.get(
            "eligible_for_primary_analysis"
        )
        is True,
        "reference_identity_complete": bool(
            _present(compound.get("inchikey")) and _present(compound.get("formula"))
        ),
        "source_origin_id": (record.get("independence") or {}).get(
            "source_origin_id"
        ),
        "laboratory_id": source.get("laboratory_id"),
        "exact_instrument_model": _present(acquisition.get("instrument_model")),
        "activation_literal": _present(acquisition.get("activation")),
        "native_collision_energy_program": bool(
            ce.get("status") == "resolved" and ce.get("conversion_applied") is False
        ),
        "polarity_adduct_and_charge": bool(
            acquisition.get("polarity") in {"positive", "negative"}
            and _present(acquisition.get("adduct"))
            and isinstance(acquisition.get("charge"), int)
        ),
        "scan_range_contains_precursor": eligibility.get(
            "range_contains_precursor"
        )
        is True,
        "precursor_selection_verified": eligibility.get(
            "precursor_selection_verified"
        )
        is True,
        "raw_acquisition_mapping": bool(
            _present(raw.get("raw_file_id")) and _present(raw.get("scan_id"))
        ),
        "verified_single_compound_material": eligibility.get(
            "verified_single_compound_reference_acquisition"
        )
        is True,
        "verified_precursor_isolation_window": isolation.get("status")
        == "verified",
        "not_composite_or_in_source": not composite,
    }


def screen_manifest_sources(
    records: Sequence[Mapping[str, Any]],
) -> dict[str, Any]:
    """Aggregate per-record gates without turning a repository into one source."""

    grouped: dict[str, list[Mapping[str, Any]]] = defaultdict(list)
    for record in records:
        source = record.get("source") or {}
        source_id = str(source.get("source_id") or "unreported_source")
        grouped[source_id].append(record)
    rows: list[dict[str, Any]] = []
    overall_status = Counter()
    for source_id, source_records in sorted(grouped.items()):
        decisions = [screen_source_descriptor(_record_descriptor(row)) for row in source_records]
        quantitative_count = sum(
            decision["quantitative_use_allowed"] for decision in decisions
        )
        negative_count = sum(
            decision["negative_evidence_allowed"] for decision in decisions
        )
        reasons = Counter(
            reason
            for decision in decisions
            for reason in decision.get("missing_requirement_codes", [])
        )
        negative_reasons = Counter(
            reason
            for decision in decisions
            for reason in decision.get(
                "negative_evidence_missing_requirement_codes", []
            )
        )
        first = source_records[0].get("source") or {}
        status = "quantitative_candidate" if quantitative_count else "audit_only"
        overall_status[status] += 1
        rows.append(
            {
                "source_id": source_id,
                "repository": first.get("repository"),
                "origin_dataset_ids": sorted(
                    {
                        str((row.get("source") or {}).get("origin_dataset_id"))
                        for row in source_records
                    }
                ),
                "laboratory_ids": sorted(
                    {
                        str((row.get("source") or {}).get("laboratory_id"))
                        for row in source_records
                    }
                ),
                "record_count": len(source_records),
                "identity_count": len(
                    {
                        str((row.get("compound") or {}).get("connectivity"))
                        for row in source_records
                    }
                ),
                "current_quantitative_candidate_record_count": quantitative_count,
                "current_negative_evidence_record_count": negative_count,
                "status": status,
                "missing_requirement_record_counts": dict(sorted(reasons.items())),
                "negative_missing_requirement_record_counts": dict(
                    sorted(negative_reasons.items())
                ),
            }
        )
    return {
        "source_count": len(rows),
        "status_counts": dict(sorted(overall_status.items())),
        "sources": rows,
        "repository_aggregate_counts_do_not_equal_independent_sources": True,
    }
