"""MassBank adapter and immutable V3 corpus construction."""

from __future__ import annotations

import hashlib
import json
import math
import re
from collections import Counter, defaultdict
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any, TextIO

from matching.feature_fragment_benchmark.corpus import (
    RecordReject,
    _charge_from_adduct,
    _chemont,
    _first_link,
    _links,
    _normalize_peaks,
    _precursor_mz,
    _structure_checks,
    _tags,
    canonical_bytes,
    file_sha256,
    iter_json_array,
    normalize_formula,
    sha256_bytes,
    stable_compound_name,
)

from .parsing import parse_collision_energy_program, parse_scan_range
from .provenance import MetadataValidationError, normalize_public_spectrum


def _write_jsonl(handle: TextIO, value: Mapping[str, Any]) -> bytes:
    line = canonical_bytes(value) + b"\n"
    handle.write(line.decode("utf-8"))
    return line


def _download(url: str, target: Path) -> None:
    import shutil
    import urllib.request

    target.parent.mkdir(parents=True, exist_ok=True)
    partial = target.with_suffix(target.suffix + ".partial")
    request = urllib.request.Request(
        url,
        headers={"User-Agent": "Mssystem-V15-feature-fragment-benchmark/3.0"},
    )
    with urllib.request.urlopen(request, timeout=240) as response, partial.open(
        "wb"
    ) as handle:
        shutil.copyfileobj(response, handle)
    partial.replace(target)


def _domain_metadata(
    *,
    source: Mapping[str, Any],
    domains_by_id: Mapping[str, Mapping[str, Any]],
    instrument: str,
    activation: str,
    ce: Mapping[str, Any],
) -> dict[str, Any]:
    candidates = []
    for domain_id in source.get("domain_ids") or []:
        domain = domains_by_id.get(str(domain_id))
        if domain is None:
            raise ValueError(f"unknown domain id in source registry: {domain_id}")
        if instrument not in (domain.get("models") or []):
            continue
        if activation not in (domain.get("activation") or []):
            continue
        if ce.get("unit") not in (domain.get("ce_units") or []):
            continue
        if ce.get("mode") not in (domain.get("ce_modes") or []):
            continue
        candidates.append(domain)
    if len(candidates) == 1:
        return dict(candidates[0])
    model_candidates = [
        domains_by_id[str(domain_id)]
        for domain_id in source.get("domain_ids") or []
        if instrument in (domains_by_id[str(domain_id)].get("models") or [])
    ]
    if len(model_candidates) == 1:
        return dict(model_candidates[0])
    return {
        "domain_id": None,
        "vendor": None,
        "family_id": None,
        "models": [instrument] if instrument else [],
        "analyzer": None,
        "activation": [activation] if activation else [],
        "collision_cell_regions": ["unreported"],
        "activation_location": "unreported",
    }


def _isolation_metadata(
    focused: Mapping[str, str], acquisition_tags: Mapping[str, str], precursor_mz: float
) -> dict[str, Any]:
    merged = {**acquisition_tags, **focused}
    for key in (
        "ISOLATION_WINDOW_WIDTH",
        "ISOLATION_WIDTH",
        "PRECURSOR_ISOLATION_WINDOW",
        "ISOLATION_WINDOW",
    ):
        raw = str(merged.get(key) or "").strip()
        if not raw:
            continue
        match = re.search(r"(?:\d+(?:\.\d*)?|\.\d+)", raw)
        if not match:
            continue
        width = float(match.group(0))
        if math.isfinite(width) and width > 0:
            return {
                "status": "verified",
                "raw": raw,
                "width_mz": width,
                "lower_mz": precursor_mz - width / 2.0,
                "upper_mz": precursor_mz + width / 2.0,
                "basis": f"MassBank {key} literal",
            }
    return {
        "status": "unreported",
        "raw": None,
        "width_mz": None,
        "lower_mz": None,
        "upper_mz": None,
        "basis": "no isolation-window field in repository record",
    }


def _raw_identifier(tags: Mapping[str, str], *names: str) -> str | None:
    for name in names:
        value = str(tags.get(name) or "").strip()
        if value:
            return value
    return None


def adapt_massbank_record(
    raw: Mapping[str, Any],
    *,
    source: Mapping[str, Any],
    repository: Mapping[str, Any],
    source_file_sha256: str,
    license_policy: Mapping[str, Any],
    instrument_domains: Sequence[Mapping[str, Any]],
) -> dict[str, Any]:
    """Convert one parsed MassBank API record to the public V3 schema."""
    accession = str(raw.get("accession") or "").strip()
    if not accession:
        raise MetadataValidationError("accession_missing")
    if raw.get("deprecated") not in ({}, None):
        raise MetadataValidationError("deprecated_record")
    compound_raw = raw.get("compound") or {}
    acquisition_raw = raw.get("acquisition") or {}
    acquisition_ms = acquisition_raw.get("mass_spectrometry") or {}
    acquisition_tags = _tags(acquisition_ms.get("subtags"))
    focused = _tags((raw.get("mass_spectrometry") or {}).get("focused_ion"))
    links = _links(compound_raw.get("link"))
    try:
        name, name_index = stable_compound_name(compound_raw.get("names"))
        formula = normalize_formula(compound_raw.get("formula"))
        chemont = _chemont(_first_link(links, "ChemOnt"))
    except RecordReject as exc:
        raise MetadataValidationError(exc.code, exc.detail) from exc
    inchikey = _first_link(links, "INCHIKEY")
    inchi = str(compound_raw.get("inchi") or "")
    try:
        structure_validation = _structure_checks(inchi, inchikey, formula)
    except RecordReject as exc:
        raise MetadataValidationError(exc.code, exc.detail) from exc

    activation = str(acquisition_tags.get("FRAGMENTATION_MODE") or "").strip()
    ce = parse_collision_energy_program(acquisition_tags.get("COLLISION_ENERGY"))
    instrument = str(acquisition_raw.get("instrument") or "").strip()
    domain = _domain_metadata(
        source=source,
        domains_by_id={str(row["domain_id"]): row for row in instrument_domains},
        instrument=instrument,
        activation=activation,
        ce=ce,
    )
    polarity_literal = str(acquisition_ms.get("ion_mode") or "").strip().upper()
    polarity = (
        polarity_literal.lower()
        if polarity_literal in {"POSITIVE", "NEGATIVE"}
        else "unresolved"
    )
    adduct = str(
        focused.get("PRECURSOR_TYPE") or focused.get("ION_TYPE") or ""
    ).strip()
    try:
        charge = _charge_from_adduct(adduct)
    except RecordReject:
        charge = None
    if charge is not None and polarity in {"positive", "negative"}:
        if (charge > 0) != (polarity == "positive"):
            raise MetadataValidationError(
                "adduct_polarity_inconsistent", f"{adduct}/{polarity}"
            )
    precursor_raw = focused.get("PRECURSOR_M/Z") or focused.get(
        "FULL_SCAN_FRAGMENT_ION_PEAK"
    )
    try:
        precursor_mz, precursor_provenance = _precursor_mz(
            precursor_raw,
            neutral_exact_mass=compound_raw.get("mass"),
            adduct=adduct,
            source=source,
        )
    except RecordReject as exc:
        raise MetadataValidationError(exc.code, exc.detail) from exc
    peaks = _normalize_peaks(raw)
    scan_range = parse_scan_range(acquisition_tags.get("MASS_RANGE_M/Z"))
    isolation = _isolation_metadata(focused, acquisition_tags, precursor_mz)
    authors = [
        str(row.get("name") or "").strip()
        for row in raw.get("authors") or []
        if str(row.get("name") or "").strip()
    ]
    record_url = f"{repository['api_root']}/records/{accession}"
    attribution = "; ".join(
        filter(
            None,
            [
                name,
                accession,
                ", ".join(authors),
                repository.get("name"),
                record_url,
            ],
        )
    )
    classes_text = " ".join(str(value) for value in compound_raw.get("classes") or [])
    material_claim = (
        "repository_mentions_reference_or_standard"
        if re.search(r"\b(?:reference|standard)\b", classes_text, flags=re.I)
        else "not_reported"
    )
    adapter_record = {
        "spectrum_id": f"MassBank:{accession}",
        "source": {
            "repository": repository["name"],
            "repository_release": repository["release"],
            "repository_release_timestamp": repository.get("release_timestamp"),
            "source_id": source["source_id"],
            "origin_dataset_id": source.get("origin_dataset_id")
            or source["source_id"],
            "origin_study_id": source.get("origin_study_id")
            or source.get("research_id"),
            "origin_study_doi": source.get("origin_study_doi"),
            "laboratory_id": source.get("laboratory_id"),
            "contributor": source.get("contributor"),
            "accession": accession,
            "usi": source.get("usi_template", "").format(accession=accession)
            if source.get("usi_template")
            else None,
            "record_url": record_url,
            "query_url": source["query_url"],
            "license_label": str(raw.get("license") or ""),
            "attribution": attribution,
            "rights_holder": None,
            "rights_holder_status": "not_reported_in_record",
            "creators": authors,
            "source_file_sha256": source_file_sha256,
            "record_canonical_json_sha256": sha256_bytes(canonical_bytes(raw)),
            "record_hash_scope": "canonical parsed MassBank API record",
            "record_created": (raw.get("date") or {}).get("created"),
            "record_updated": (raw.get("date") or {}).get("updated"),
            "retrieved_on": repository.get("retrieved_on"),
        },
        "compound": {
            "display_name": name,
            "name_provenance": {
                "source": "MassBank compound.names",
                "list_index": name_index,
                "accession": accession,
            },
            "database_identifiers": links,
            "inchi": inchi or None,
            "inchikey": inchikey,
            "formula": formula,
            "smiles": compound_raw.get("smiles"),
            "structure_validation": structure_validation,
            "chemical_class": chemont,
            "identity_basis": "repository_reported_structure_with_inchi_consistency_check",
            "reference_material_claim": material_claim,
            "standard_confirmation_status": "not_confirmed_by_mssystem",
            "material_type": None,
            "purity": None,
            "batch_id": None,
        },
        "acquisition": {
            "instrument_vendor": domain.get("vendor"),
            "instrument_family": domain.get("family_id"),
            "instrument_model": instrument,
            "instrument_type": acquisition_raw.get("instrument_type"),
            "analyzer": domain.get("analyzer"),
            "activation_raw": activation,
            "activation": activation,
            "activation_location": domain.get("activation_location"),
            "activation_evidence": "MassBank contributor literal; family architecture does not upgrade the record CV",
            "collision_cell_region": (domain.get("collision_cell_regions") or [None])[0],
            "collision_energy": ce,
            "scan_range": scan_range,
            "polarity": polarity,
            "adduct": adduct or None,
            "charge": charge,
            "precursor_mz": precursor_mz,
            "precursor_mz_provenance": precursor_provenance,
            "precursor_isolation": isolation,
            "ms_level": 2,
            "peak_representation": "centroid",
        },
        "raw_acquisition": {
            "raw_file_id": _raw_identifier(
                acquisition_tags, "RAW_FILE", "RAW_FILE_NAME", "FILE_NAME"
            ),
            "scan_id": _raw_identifier(
                acquisition_tags, "SCAN_NUMBER", "SCAN_ID", "NATIVE_ID"
            ),
            "injection_id": _raw_identifier(
                acquisition_tags, "INJECTION_ID", "SAMPLE_ID"
            ),
            "acquisition_event_id": _raw_identifier(
                acquisition_tags, "ACQUISITION_EVENT_ID", "EVENT_ID"
            ),
        },
        "peaks": peaks,
    }
    return normalize_public_spectrum(
        adapter_record,
        license_policy=license_policy,
        instrument_domains=instrument_domains,
    )


def build_massbank_corpus_v3(
    source_manifest_path: str | Path,
    domain_path: str | Path,
    license_policy_path: str | Path,
    *,
    cache_dir: str | Path,
    output_jsonl: str | Path,
    rejects_jsonl: str | Path,
    download_missing: bool = False,
    max_records: int | None = None,
) -> dict[str, Any]:
    """Build an auditable V3 MassBank corpus without touching V2 artifacts."""
    source_manifest_path = Path(source_manifest_path)
    domain_path = Path(domain_path)
    license_policy_path = Path(license_policy_path)
    manifest = json.loads(source_manifest_path.read_text(encoding="utf-8"))
    domain_document = json.loads(domain_path.read_text(encoding="utf-8"))
    license_policy = json.loads(license_policy_path.read_text(encoding="utf-8"))
    repository = manifest.get("repository") or {}
    sources = manifest.get("sources") or []
    instrument_domains = domain_document.get("instrument_domains") or []
    if not sources or not instrument_domains:
        raise ValueError("V3 corpus requires sources and instrument domains")

    cache_dir = Path(cache_dir)
    output_jsonl = Path(output_jsonl)
    rejects_jsonl = Path(rejects_jsonl)
    output_jsonl.parent.mkdir(parents=True, exist_ok=True)
    rejects_jsonl.parent.mkdir(parents=True, exist_ok=True)
    corpus_digest = hashlib.sha256()
    reject_digest = hashlib.sha256()
    raw_total = accepted = rejected = 0
    source_counts: dict[str, dict[str, int]] = {}
    source_file_hashes: dict[str, str] = {}
    reason_counts: Counter[str] = Counter()
    domain_counts: Counter[str] = Counter()
    license_tiers: Counter[str] = Counter()
    scan_status: Counter[str] = Counter()
    scan_contains_precursor: Counter[str] = Counter()
    isolation_status: Counter[str] = Counter()
    ce_unit_mode_counts: Counter[str] = Counter()
    exact_condition_completeness: Counter[str] = Counter()
    reference_material_claims: Counter[str] = Counter()
    positive_eligibility: Counter[str] = Counter()
    negative_eligibility: Counter[str] = Counter()
    missing_fields: Counter[str] = Counter()
    identities: set[str] = set()
    formulas: set[str] = set()
    acquisition_replication_groups: set[str] = set()

    with output_jsonl.open("w", encoding="utf-8") as out, rejects_jsonl.open(
        "w", encoding="utf-8"
    ) as rejects:
        for source in sources:
            source_id = str(source["source_id"])
            cache_path = cache_dir / str(source["cache_filename"])
            if not cache_path.exists():
                if not download_missing:
                    raise ValueError(f"missing explicit public-data cache: {cache_path}")
                _download(str(source["query_url"]), cache_path)
            file_hash = file_sha256(cache_path)
            expected_hash = source.get("expected_source_file_sha256")
            if expected_hash and expected_hash != file_hash:
                raise ValueError(
                    f"{source_id}: source hash differs from preregistered hash"
                )
            source_file_hashes[source_id] = file_hash
            local = Counter()
            local_identities: set[str] = set()
            local_formulas: set[str] = set()
            for raw in iter_json_array(cache_path):
                if max_records is not None and raw_total >= max_records:
                    break
                raw_total += 1
                local["raw"] += 1
                accession = str(raw.get("accession") or "") if isinstance(raw, Mapping) else ""
                try:
                    if not isinstance(raw, Mapping):
                        raise MetadataValidationError("source_record_not_object")
                    normalized = adapt_massbank_record(
                        raw,
                        source=source,
                        repository=repository,
                        source_file_sha256=file_hash,
                        license_policy=license_policy,
                        instrument_domains=instrument_domains,
                    )
                except (MetadataValidationError, RecordReject) as exc:
                    code = getattr(exc, "code", type(exc).__name__)
                    detail = getattr(exc, "detail", str(exc))
                    row = {
                        "source_id": source_id,
                        "accession": accession,
                        "reason_code": code,
                        "detail": detail,
                    }
                    line = _write_jsonl(rejects, row)
                    reject_digest.update(line)
                    rejected += 1
                    local["rejected"] += 1
                    reason_counts[code] += 1
                    continue
                line = _write_jsonl(out, normalized)
                corpus_digest.update(line)
                accepted += 1
                local["accepted"] += 1
                domain_counts[
                    normalized["acquisition"].get("instrument_domain_id")
                    or "unresolved"
                ] += 1
                license_tiers[normalized["source"]["license"]["tier"]] += 1
                scan_status[
                    normalized["acquisition"]["scan_range"]["status"]
                ] += 1
                scan_contains_precursor[
                    str(
                        normalized["quantitative_eligibility"][
                            "range_contains_precursor"
                        ]
                    ).lower()
                ] += 1
                isolation_status[
                    normalized["acquisition"]["precursor_isolation"]["status"]
                ] += 1
                ce = normalized["acquisition"]["collision_energy"]
                ce_unit_mode_counts[
                    f"{ce.get('unit') or 'unresolved'}|{ce.get('mode') or 'unresolved'}"
                ] += 1
                exact_condition_completeness[
                    str(
                        bool(
                            normalized["acquisition"][
                                "instrument_domain_resolution"
                            ].get("exact_condition_metadata_complete")
                        )
                    ).lower()
                ] += 1
                reference_material_claims[
                    normalized["compound"]["reference_material_claim"]
                ] += 1
                eligibility = normalized["quantitative_eligibility"]
                positive_eligibility[str(bool(eligibility["positive_lr_eligible"])).lower()] += 1
                negative_eligibility[str(bool(eligibility["negative_lr_eligible"])).lower()] += 1
                missing_fields.update(normalized["metadata_completeness"]["missing_fields"])
                acquisition_replication_groups.add(
                    normalized["independence"][
                        "acquisition_replication_group_id"
                    ]
                )
                identity = normalized["compound"]["connectivity"]
                formula = normalized["compound"].get("formula")
                identities.add(identity)
                local_identities.add(identity)
                if formula:
                    formulas.add(formula)
                    local_formulas.add(formula)
            if max_records is None and int(source["expected_record_count"]) != local["raw"]:
                raise ValueError(
                    f"{source_id}: expected {source['expected_record_count']} records, observed {local['raw']}"
                )
            source_counts[source_id] = {
                "raw_records": local["raw"],
                "accepted_records": local["accepted"],
                "rejected_records": local["rejected"],
                "two_dimensional_identities": len(local_identities),
                "formula_families": len(local_formulas),
            }
            if max_records is not None and raw_total >= max_records:
                break

    return {
        "schema_version": "3.0.0",
        "status": "built_research_only",
        "runtime_effect": "none_fail_closed",
        "source_manifest_sha256": file_sha256(source_manifest_path),
        "instrument_domains_sha256": file_sha256(domain_path),
        "license_policy_sha256": file_sha256(license_policy_path),
        "repository_release": repository.get("release"),
        "raw_record_count": raw_total,
        "accepted_audit_record_count": accepted,
        "rejected_record_count": rejected,
        "two_dimensional_identity_count": len(identities),
        "formula_family_count": len(formulas),
        "source_counts": source_counts,
        "domain_counts": dict(sorted(domain_counts.items())),
        "license_tier_counts": dict(sorted(license_tiers.items())),
        "scan_range_status_counts": dict(sorted(scan_status.items())),
        "scan_range_contains_precursor_counts": dict(
            sorted(scan_contains_precursor.items())
        ),
        "precursor_isolation_status_counts": dict(sorted(isolation_status.items())),
        "collision_energy_unit_mode_counts": dict(
            sorted(ce_unit_mode_counts.items())
        ),
        "exact_condition_metadata_complete_counts": dict(
            sorted(exact_condition_completeness.items())
        ),
        "reference_material_claim_counts": dict(
            sorted(reference_material_claims.items())
        ),
        "acquisition_replication_group_count": len(
            acquisition_replication_groups
        ),
        "positive_lr_eligibility_counts": dict(sorted(positive_eligibility.items())),
        "negative_lr_eligibility_counts": dict(sorted(negative_eligibility.items())),
        "metadata_missing_field_counts": dict(sorted(missing_fields.items())),
        "rejection_reason_counts": dict(sorted(reason_counts.items())),
        "source_file_sha256": source_file_hashes,
        "corpus_jsonl": str(output_jsonl),
        "corpus_jsonl_sha256": corpus_digest.hexdigest(),
        "rejects_jsonl": str(rejects_jsonl),
        "rejects_jsonl_sha256": reject_digest.hexdigest(),
    }
