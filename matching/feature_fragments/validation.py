"""Offline validation for the public V15 feature-fragment pilot."""

from __future__ import annotations

import hashlib
import json
import math
import re
from collections import defaultdict
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from .evidence import conditional_log_likelihood_ratio, mass_window
from .schema import CatalogValidationError, stable_compound_name, validate_catalog


_FORMULA_TOKEN = re.compile(r"([A-Z][a-z]?)(\d*)")
_INCHIKEY = re.compile(r"[A-Z]{14}-[A-Z]{10}-[A-Z]")
REQUIRED_VENDORS = {
    "Thermo Fisher Scientific",
    "Agilent Technologies",
    "SCIEX",
    "Bruker Daltonics",
    "Waters Corporation",
}


class FixtureValidationError(ValueError):
    """Raised when a public spectrum fixture is incomplete or inconsistent."""


def _canonical_bytes(value: Any) -> bytes:
    return json.dumps(
        value, ensure_ascii=False, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")


def _sha256(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def file_sha256(path: Path) -> str:
    hasher = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            hasher.update(block)
    return hasher.hexdigest()


def normalize_formula(formula: str) -> str:
    """Normalize a simple molecular formula to Hill order."""
    text = re.sub(r"[\s()+\-]", "", str(formula or ""))
    counts: dict[str, int] = defaultdict(int)
    consumed = ""
    for element, raw_count in _FORMULA_TOKEN.findall(text):
        counts[element] += int(raw_count or "1")
        consumed += f"{element}{raw_count}"
    if not counts or consumed != text:
        raise FixtureValidationError(f"unparseable molecular formula: {formula!r}")
    order: list[str] = []
    if "C" in counts:
        order.append("C")
    if "H" in counts:
        order.append("H")
    order.extend(sorted(element for element in counts if element not in {"C", "H"}))
    return "".join(
        element + (str(counts[element]) if counts[element] != 1 else "")
        for element in order
    )


def _peak_only_hash(peaks: Sequence[Sequence[float]]) -> str:
    base = max(float(peak[1]) for peak in peaks)
    normalized = [
        [round(float(peak[0]), 5), round(float(peak[1]) / base, 6)]
        for peak in sorted(peaks, key=lambda item: float(item[0]))
    ]
    return _sha256(_canonical_bytes(normalized))


def _acquisition_hash(record: Mapping[str, Any]) -> str:
    compound = record["compound"]
    acquisition = record["acquisition"]
    identity = {
        "connectivity": compound["connectivity"],
        "precursor_mz": round(float(acquisition["precursor_mz"]), 6),
        "adduct": acquisition["adduct"],
        "polarity": acquisition["polarity"],
        "instrument": acquisition["instrument_model"],
        "activation": acquisition["activation"],
        "collision_cell_region": acquisition["collision_cell_region"],
        "collision_energy": {
            "value": acquisition["collision_energy"]["value"],
            "unit": acquisition["collision_energy"]["unit"],
            "mode": acquisition["collision_energy"]["mode"],
        },
        "peak_hash": record["independence"]["peak_only_mirror_hash"],
    }
    return _sha256(_canonical_bytes(identity))


def validate_fixture(fixture: Mapping[str, Any]) -> dict[str, Any]:
    """Validate provenance, identity, acquisition, independence, and hashes."""
    spectra = fixture.get("spectra")
    if not isinstance(spectra, list) or not spectra:
        raise FixtureValidationError("fixture.spectra must be a non-empty list")
    without_hash = dict(fixture)
    recorded_hash = without_hash.pop("content_sha256", None)
    if recorded_hash != _sha256(_canonical_bytes(without_hash)):
        raise FixtureValidationError("fixture content_sha256 mismatch")
    if fixture.get("record_count") != len(spectra):
        raise FixtureValidationError("fixture record_count mismatch")

    source_records: set[str] = set()
    spectrum_ids: set[str] = set()
    acquisition_hashes: set[str] = set()
    identities: set[str] = set()
    formulas: set[str] = set()
    mirror_clusters: dict[str, list[str]] = defaultdict(list)
    exact_domain_groups: set[tuple[str, str, float, str]] = set()
    required_quality = {
        "complete_mass_range",
        "centroided",
        "sufficient_signal",
        "precursor_isolated",
        "collision_energy_verified",
        "instrument_domain_verified",
    }

    for index, record in enumerate(spectra):
        path = f"spectra[{index}]"
        for field in (
            "spectrum_id",
            "evidence_stage",
            "experimental_status",
            "source",
            "compound",
            "acquisition",
            "independence",
            "quality",
            "peaks",
        ):
            if record.get(field) in (None, "", [], {}):
                raise FixtureValidationError(f"{path}: missing {field}")
        if record["evidence_stage"] != "raw_msms_spectrum":
            raise FixtureValidationError(f"{path}: invalid evidence stage")
        if record["experimental_status"] != "public_experimental_spectrum":
            raise FixtureValidationError(f"{path}: predicted spectra are forbidden")
        if record.get("query_confirmation_status") != "not_applicable":
            raise FixtureValidationError(f"{path}: query confirmation was upgraded")

        source = record["source"]
        for field in ("repository_record", "url", "license", "record_sha256"):
            if not source.get(field):
                raise FixtureValidationError(f"{path}.source: missing {field}")
        if source["license"] != "CC BY":
            raise FixtureValidationError(f"{path}: incompatible record license")
        if not re.fullmatch(r"[0-9a-f]{64}", source["record_sha256"]):
            raise FixtureValidationError(f"{path}: invalid source record hash")
        if source["repository_record"] in source_records:
            raise FixtureValidationError(f"{path}: duplicate source record")
        source_records.add(source["repository_record"])

        compound = record["compound"]
        if not stable_compound_name(compound.get("name")):
            raise FixtureValidationError(f"{path}: database identifier used as name")
        if not _INCHIKEY.fullmatch(str(compound.get("inchikey") or "")):
            raise FixtureValidationError(f"{path}: invalid full InChIKey")
        if compound["connectivity"] != compound["inchikey"][:14]:
            raise FixtureValidationError(f"{path}: connectivity/InChIKey mismatch")
        if not compound.get("class_provenance", {}).get("raw"):
            raise FixtureValidationError(f"{path}: missing class provenance")
        formula = normalize_formula(compound["formula"])
        if formula != compound["formula"]:
            raise FixtureValidationError(f"{path}: formula is not normalized")
        identities.add(compound["connectivity"])
        formulas.add(formula)

        acquisition = record["acquisition"]
        for field in (
            "instrument_vendor",
            "instrument_family",
            "instrument_model",
            "analyzer",
            "activation",
            "collision_cell_region",
            "instrument_domain_id",
            "collision_energy",
            "adduct",
            "polarity",
            "charge",
            "precursor_mz",
        ):
            if acquisition.get(field) in (None, "", [], {}):
                raise FixtureValidationError(f"{path}.acquisition: missing {field}")
        ce = acquisition["collision_energy"]
        if ce.get("unit") != "eV" or ce.get("mode") != "fixed":
            raise FixtureValidationError(f"{path}: fixture CE domain changed")
        if ce.get("value") not in {10.0, 20.0, 40.0}:
            raise FixtureValidationError(f"{path}: unexpected CE")
        if acquisition["adduct"] != "[M+H]+" or acquisition["charge"] != 1:
            raise FixtureValidationError(f"{path}: fixture adduct/charge changed")

        peaks = record["peaks"]
        if not peaks or any(len(peak) < 3 for peak in peaks):
            raise FixtureValidationError(f"{path}: invalid peak list")
        independence = record["independence"]
        peak_hash = _peak_only_hash(peaks)
        if peak_hash != independence.get("peak_only_mirror_hash"):
            raise FixtureValidationError(f"{path}: peak-only hash mismatch")
        if _acquisition_hash(record) != independence.get("acquisition_aware_hash"):
            raise FixtureValidationError(f"{path}: acquisition-aware hash mismatch")
        mirror_clusters[peak_hash].append(record["spectrum_id"])
        acquisition_hashes.add(independence["acquisition_aware_hash"])
        group = (
            compound["connectivity"],
            acquisition["instrument_domain_id"],
            float(ce["value"]),
            acquisition["adduct"],
        )
        if group in exact_domain_groups:
            raise FixtureValidationError(
                f"{path}: repeated identity within exact acquisition domain"
            )
        exact_domain_groups.add(group)

        if not required_quality.issubset(record["quality"]):
            raise FixtureValidationError(f"{path}: incomplete quality gates")
        spectrum_id = record["spectrum_id"]
        if spectrum_id in spectrum_ids:
            raise FixtureValidationError(f"{path}: duplicate spectrum_id")
        spectrum_ids.add(spectrum_id)

    if fixture.get("independent_2d_identity_count") != len(identities):
        raise FixtureValidationError("fixture identity count mismatch")
    if fixture.get("normalized_formula_family_count") != len(formulas):
        raise FixtureValidationError("fixture formula-family count mismatch")
    mirrored = {key: value for key, value in mirror_clusters.items() if len(value) > 1}
    if mirrored:
        raise FixtureValidationError("fixture contains duplicated peak-only mirrors")
    return {
        "record_count": len(spectra),
        "independent_2d_identity_count": len(identities),
        "normalized_formula_family_count": len(formulas),
        "source_record_count": len(source_records),
        "acquisition_aware_hash_count": len(acquisition_hashes),
        "peak_only_mirror_duplicate_count": 0,
    }


def _rule_record_applicable(rule: Mapping[str, Any], record: Mapping[str, Any]) -> bool:
    acquisition = record["acquisition"]
    ce = acquisition["collision_energy"]
    rule_ce = rule["collision_energy"]
    return (
        acquisition["polarity"] == rule["polarity"]
        and acquisition["adduct"] in rule["adducts"]
        and acquisition["charge"] in rule["charges"]
        and acquisition["instrument_domain_id"] in rule["instrument_domain_ids"]
        and acquisition["activation"] == rule["activation"]
        and ce["unit"] == rule_ce["unit"]
        and ce["mode"] == rule_ce["mode"]
        and float(rule_ce["minimum"]) <= float(ce["value"]) <= float(rule_ce["maximum"])
    )


def _rule_present(rule: Mapping[str, Any], record: Mapping[str, Any]) -> bool:
    peaks = record["peaks"]
    base = max(float(peak[1]) for peak in peaks)
    if rule["observation_type"] == "product_ion":
        expected = float(rule["product_mz"])
    else:
        expected = float(record["acquisition"]["precursor_mz"]) - float(
            rule["neutral_loss_da"]
        )
    window = mass_window(expected, rule["mass_tolerance"])
    floor = float(rule["minimum_relative_intensity"])
    return any(
        abs(float(peak[0]) - expected) <= window and float(peak[1]) / base >= floor
        for peak in peaks
    )


def _deduplicated_rule_rows(
    rule: Mapping[str, Any], spectra: Sequence[Mapping[str, Any]]
) -> list[dict[str, Any]]:
    background_classes = set(rule["background_statistics"].get("classes", []))
    target_classes = set(rule["chemical_classes"])
    rows: dict[str, dict[str, Any]] = {}
    for record in spectra:
        chemical_class = record["compound"]["chemical_class"]
        if chemical_class not in target_classes | background_classes:
            continue
        if not _rule_record_applicable(rule, record):
            continue
        connectivity = record["compound"]["connectivity"]
        if connectivity in rows:
            raise FixtureValidationError(
                f"{rule['rule_id']}: more than one record per 2D identity/domain"
            )
        rows[connectivity] = {
            "connectivity": connectivity,
            "formula_family": normalize_formula(record["compound"]["formula"]),
            "target": chemical_class in target_classes,
            "present": _rule_present(rule, record),
            "source_record": record["source"]["repository_record"],
        }
    return sorted(rows.values(), key=lambda row: row["connectivity"])


def _auc(labels: Sequence[bool], scores: Sequence[float]) -> float:
    positives = [score for label, score in zip(labels, scores) if label]
    negatives = [score for label, score in zip(labels, scores) if not label]
    if not positives or not negatives:
        return float("nan")
    wins = 0.0
    for positive in positives:
        for negative in negatives:
            wins += 1.0 if positive > negative else 0.5 if positive == negative else 0.0
    return wins / (len(positives) * len(negatives))


def _sigmoid(value: float) -> float:
    if value >= 0:
        return 1.0 / (1.0 + math.exp(-value))
    exp_value = math.exp(value)
    return exp_value / (1.0 + exp_value)


def _heldout_metrics(rows: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    """Leave out every identity sharing the same normalized formula."""
    predictions: list[dict[str, Any]] = []
    formula_families = sorted({row["formula_family"] for row in rows})
    leakage_violations: list[dict[str, Any]] = []
    for formula in formula_families:
        test = [row for row in rows if row["formula_family"] == formula]
        train = [row for row in rows if row["formula_family"] != formula]
        train_formulae = {row["formula_family"] for row in train}
        test_formulae = {row["formula_family"] for row in test}
        train_ids = {row["connectivity"] for row in train}
        test_ids = {row["connectivity"] for row in test}
        if train_formulae & test_formulae or train_ids & test_ids:
            leakage_violations.append(
                {"formula_family": formula, "reason": "identity_or_formula_overlap"}
            )
            continue
        target_train = [row for row in train if row["target"]]
        background_train = [row for row in train if not row["target"]]
        if not target_train or not background_train:
            continue
        target_present = sum(bool(row["present"]) for row in target_train)
        background_present = sum(bool(row["present"]) for row in background_train)
        prior_log_odds = math.log(
            (len(target_train) + 0.5) / (len(background_train) + 0.5)
        )
        for row in test:
            lr = conditional_log_likelihood_ratio(
                target_present=target_present,
                target_total=len(target_train),
                background_present=background_present,
                background_total=len(background_train),
                observed_present=bool(row["present"]),
            )
            log_odds = prior_log_odds + lr["log_likelihood_ratio"]
            predictions.append(
                {
                    "connectivity": row["connectivity"],
                    "formula_family": row["formula_family"],
                    "target": bool(row["target"]),
                    "present": bool(row["present"]),
                    "training_target_identities": len(target_train),
                    "training_background_identities": len(background_train),
                    "log_likelihood_ratio": lr["log_likelihood_ratio"],
                    "posterior_probability_from_training_prior": _sigmoid(log_odds),
                }
            )
    labels = [bool(row["target"]) for row in predictions]
    scores = [
        float(row["posterior_probability_from_training_prior"]) for row in predictions
    ]
    epsilon = 1e-12
    brier = (
        sum((score - float(label)) ** 2 for label, score in zip(labels, scores))
        / len(predictions)
        if predictions
        else float("nan")
    )
    log_loss = (
        -sum(
            float(label) * math.log(max(score, epsilon))
            + (1.0 - float(label)) * math.log(max(1.0 - score, epsilon))
            for label, score in zip(labels, scores)
        )
        / len(predictions)
        if predictions
        else float("nan")
    )
    accuracy = (
        sum((score >= 0.5) == label for label, score in zip(labels, scores))
        / len(predictions)
        if predictions
        else float("nan")
    )
    return {
        "method": "leave_one_normalized_formula_family_out",
        "formula_family_count": len(formula_families),
        "heldout_identity_predictions": len(predictions),
        "identity_overlap_violations": 0,
        "formula_family_overlap_violations": len(leakage_violations),
        "roc_auc": _auc(labels, scores),
        "brier_score": brier,
        "log_loss": log_loss,
        "accuracy_at_0_5": accuracy,
        "predictions": predictions,
    }


def _instrument_source_audit(catalog: Mapping[str, Any]) -> dict[str, Any]:
    sources = {source["source_id"]: source for source in catalog["sources"]}
    by_vendor: dict[str, dict[str, Any]] = {}
    for vendor in sorted(REQUIRED_VENDORS):
        domains = [
            domain
            for domain in catalog["instrument_domains"]
            if domain["vendor"] == vendor
        ]
        authoritative = sorted(
            {
                source_id
                for domain in domains
                for source_id in domain["source_ids"]
                if sources[source_id]["source_type"]
                in {"official_vendor_document", "peer_reviewed_primary"}
            }
        )
        by_vendor[vendor] = {
            "instrument_domain_ids": [domain["domain_id"] for domain in domains],
            "authoritative_source_ids": authoritative,
            "passed": bool(domains and authoritative),
        }
    missing = sorted(
        vendor for vendor, result in by_vendor.items() if not result["passed"]
    )
    return {
        "required_vendors": sorted(REQUIRED_VENDORS),
        "vendor_results": by_vendor,
        "missing_or_unsourced_vendors": missing,
        "cross_unit_conversion_allowed": False,
        "passed": not missing,
    }


def validate_public_pilot(
    catalog: Mapping[str, Any], fixture: Mapping[str, Any]
) -> dict[str, Any]:
    """Validate schema, configured counts, heldout performance, and source coverage."""
    validate_catalog(catalog)
    fixture_summary = validate_fixture(fixture)
    rule_results: list[dict[str, Any]] = []
    for rule in catalog["rules"]:
        rows = _deduplicated_rule_rows(rule, fixture["spectra"])
        target = [row for row in rows if row["target"]]
        background = [row for row in rows if not row["target"]]
        observed_counts = {
            "target_present": sum(bool(row["present"]) for row in target),
            "target_total": len(target),
            "background_present": sum(bool(row["present"]) for row in background),
            "background_total": len(background),
        }
        configured = {
            "target_present": int(rule["occurrence_statistics"]["present"]),
            "target_total": int(rule["occurrence_statistics"]["total"]),
            "background_present": int(rule["background_statistics"]["present"]),
            "background_total": int(rule["background_statistics"]["total"]),
        }
        if observed_counts != configured:
            raise FixtureValidationError(
                f"{rule['rule_id']}: catalog counts {configured} != fixture {observed_counts}"
            )
        metrics = _heldout_metrics(rows)
        if metrics["formula_family_overlap_violations"]:
            raise FixtureValidationError(f"{rule['rule_id']}: formula-family leakage")
        rule_results.append(
            {
                "rule_id": rule["rule_id"],
                "independent_fragmentation_group": rule[
                    "independent_fragmentation_group"
                ],
                "collision_energy": rule["collision_energy"],
                "observed_counts": observed_counts,
                "configured_counts_match": True,
                "heldout_validation": metrics,
            }
        )
    instrument_audit = _instrument_source_audit(catalog)
    if not instrument_audit["passed"]:
        raise CatalogValidationError("instrument family source audit failed")
    ce_response: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for result in rule_results:
        counts = result["observed_counts"]
        ce_response[result["independent_fragmentation_group"]].append(
            {
                "energy": result["collision_energy"]["minimum"],
                "unit": result["collision_energy"]["unit"],
                "target_frequency": counts["target_present"] / counts["target_total"],
                "background_frequency": counts["background_present"]
                / counts["background_total"],
            }
        )
    return {
        "receipt_schema_version": "1.0.0",
        "module": "v15_feature_fragments",
        "mode": "shadow_only",
        "runtime_influence": False,
        "identity_rank_influence": False,
        "fixture_summary": fixture_summary,
        "catalog_summary": {
            "source_count": len(catalog["sources"]),
            "instrument_domain_count": len(catalog["instrument_domains"]),
            "rule_count": len(catalog["rules"]),
        },
        "instrument_family_source_audit": instrument_audit,
        "rule_results": rule_results,
        "collision_energy_response": {
            group: sorted(rows, key=lambda row: row["energy"])
            for group, rows in sorted(ce_response.items())
        },
        "evidence_boundaries": {
            "raw_ions": "Fixture peak lists remain raw experimental MS/MS observations.",
            "ms1_candidate": "No MS1 candidate is created or upgraded.",
            "msms_library_candidate": "No library-match status is created or upgraded.",
            "isomer_group": "Shared fragments cannot distinguish formula-sharing or structural isomers.",
            "standard_confirmation": "No query is standard-confirmed by this module.",
        },
        "claims": {
            "engineering_pipeline_executable": True,
            "public_experimental_counts_reproduced": True,
            "identity_and_formula_family_isolation_executed": True,
            "cross_platform_rule_validated": False,
            "exact_identity_validated": False,
            "isomer_discrimination_validated": False,
        },
        "limitations": [
            "The empirical pilot contains 30 public spectra, 10 2D identities, 9 formula families, one laboratory, one Agilent 6560 QTOF, one adduct, and three reported fixed eV settings.",
            "Instrument-family records are source-validated metadata constraints; Thermo, SCIEX, Bruker, and Waters do not yet have public spectra in this fixture.",
            "Heldout metrics assess class-conditional fragment occurrence, not candidate Top-1 identity accuracy.",
            "The catalog remains shadow-only until larger independent laboratories and exact instrument domains pass the same isolation checks.",
        ],
    }
