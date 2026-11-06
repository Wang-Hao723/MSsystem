"""Schema checks for auditable feature-fragment catalogs.

The schema is enforced in Python so both command-line validation and unit
tests use exactly the same fail-closed rules.  It rejects ambiguous collision
energy metadata and incomplete provenance before a rule can be evaluated.
"""

from __future__ import annotations

import math
import re
from collections.abc import Mapping, Sequence
from typing import Any


SCHEMA_VERSION = "1.0.0"
SOURCE_TYPES = {
    "public_experimental_spectrum",
    "official_repository_document",
    "peer_reviewed_primary",
    "official_vendor_document",
    "book",
    "manual_rule",
    "computational_inference",
    "unsourced_constant",
}
OBSERVATION_TYPES = {"product_ion", "neutral_loss"}
POLARITIES = {"positive", "negative"}
CE_UNITS = {"eV", "NCE", "V", "percent_nominal"}
CE_MODES = {
    "fixed",
    "stepped",
    "ramp",
    "rolling",
    "dynamic",
    "CES",
    "bbCID",
    "MSE",
}
EVIDENCE_STAGES = {"msms_rule_support"}
_DB_IDENTIFIER_NAME = re.compile(
    r"^(?:MSBNK|CCMSLIB|NCGC|HMDB|CHEBI|CID|LM[A-Z]{2})[-_: ]?[A-Z0-9_.-]+$",
    re.IGNORECASE,
)


class CatalogValidationError(ValueError):
    """Raised when a catalog cannot support an auditable evaluation."""


def _require(mapping: Mapping[str, Any], keys: Sequence[str], path: str) -> None:
    missing = [key for key in keys if mapping.get(key) in (None, "", [], {})]
    if missing:
        raise CatalogValidationError(f"{path}: missing {', '.join(missing)}")


def _finite_nonnegative(value: Any, path: str) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise CatalogValidationError(f"{path}: expected a number") from exc
    if not math.isfinite(result) or result < 0:
        raise CatalogValidationError(f"{path}: expected a finite non-negative value")
    return result


def stable_compound_name(name: Any) -> bool:
    """Return whether *name* is a human-readable chemical name, not an ID."""
    text = str(name or "").strip()
    return bool(text and not _DB_IDENTIFIER_NAME.fullmatch(text))


def validate_source(source: Mapping[str, Any], path: str = "source") -> None:
    _require(
        source,
        ("source_id", "source_type", "title", "locator", "access_level"),
        path,
    )
    if source["source_type"] not in SOURCE_TYPES:
        raise CatalogValidationError(f"{path}.source_type: unsupported value")
    if source["source_type"] in {
        "public_experimental_spectrum",
        "official_repository_document",
        "peer_reviewed_primary",
        "official_vendor_document",
    }:
        _require(source, ("url", "supported_claims", "limitations"), path)
    if source["source_type"] == "public_experimental_spectrum":
        _require(source, ("repository", "license"), path)


def validate_instrument_domain(
    domain: Mapping[str, Any], path: str = "instrument_domain"
) -> None:
    _require(
        domain,
        (
            "domain_id",
            "vendor",
            "family",
            "analyzer",
            "activation",
            "collision_cell_region",
            "ce_unit",
            "ce_mode",
            "source_ids",
            "limitations",
        ),
        path,
    )
    if domain["ce_unit"] not in CE_UNITS:
        raise CatalogValidationError(f"{path}.ce_unit: unsupported value")
    if domain["ce_mode"] not in CE_MODES:
        raise CatalogValidationError(f"{path}.ce_mode: unsupported value")
    if domain.get("cross_unit_conversion_allowed") is not False:
        raise CatalogValidationError(
            f"{path}.cross_unit_conversion_allowed: must be false"
        )


def validate_rule(
    rule: Mapping[str, Any],
    *,
    source_ids: set[str],
    domain_ids: set[str],
    path: str = "rule",
) -> None:
    _require(
        rule,
        (
            "rule_id",
            "rule_version",
            "status",
            "observation_type",
            "claim_scope",
            "chemical_classes",
            "class_provenance",
            "adducts",
            "charges",
            "polarity",
            "instrument_domain_ids",
            "activation",
            "collision_energy",
            "mass_tolerance",
            "minimum_relative_intensity",
            "occurrence_statistics",
            "background_statistics",
            "independent_fragmentation_group",
            "positive_evidence_limitations",
            "negative_evidence",
            "source_ids",
            "evidence_level",
            "validation",
            "output_semantics",
        ),
        path,
    )
    if rule["observation_type"] not in OBSERVATION_TYPES:
        raise CatalogValidationError(f"{path}.observation_type: unsupported value")
    if rule["polarity"] not in POLARITIES:
        raise CatalogValidationError(f"{path}.polarity: unsupported value")
    if not all(isinstance(charge, int) and charge != 0 for charge in rule["charges"]):
        raise CatalogValidationError(f"{path}.charges: expected non-zero integers")
    if rule.get("runtime_influence") is not False:
        raise CatalogValidationError(f"{path}.runtime_influence: must be false")
    if rule.get("identity_rank_influence") is not False:
        raise CatalogValidationError(f"{path}.identity_rank_influence: must be false")

    if rule["observation_type"] == "product_ion":
        _finite_nonnegative(rule.get("product_mz"), f"{path}.product_mz")
    else:
        _finite_nonnegative(rule.get("neutral_loss_da"), f"{path}.neutral_loss_da")

    ce = rule["collision_energy"]
    if not isinstance(ce, Mapping):
        raise CatalogValidationError(f"{path}.collision_energy: expected object")
    _require(
        ce, ("raw", "unit", "mode", "minimum", "maximum"), f"{path}.collision_energy"
    )
    if ce["unit"] not in CE_UNITS or ce["mode"] not in CE_MODES:
        raise CatalogValidationError(f"{path}.collision_energy: unsupported domain")
    minimum = _finite_nonnegative(ce["minimum"], f"{path}.collision_energy.minimum")
    maximum = _finite_nonnegative(ce["maximum"], f"{path}.collision_energy.maximum")
    if maximum < minimum:
        raise CatalogValidationError(f"{path}.collision_energy: reversed interval")
    if ce["mode"] == "fixed" and minimum != maximum:
        raise CatalogValidationError(
            f"{path}.collision_energy: fixed mode requires one value"
        )

    tolerance = rule["mass_tolerance"]
    if not isinstance(tolerance, Mapping):
        raise CatalogValidationError(f"{path}.mass_tolerance: expected object")
    _require(tolerance, ("absolute_da", "ppm"), f"{path}.mass_tolerance")
    absolute_da = _finite_nonnegative(
        tolerance["absolute_da"], f"{path}.mass_tolerance.absolute_da"
    )
    ppm = _finite_nonnegative(tolerance["ppm"], f"{path}.mass_tolerance.ppm")
    if absolute_da == 0 and ppm == 0:
        raise CatalogValidationError(f"{path}.mass_tolerance: both values are zero")

    for stat_name in ("occurrence_statistics", "background_statistics"):
        stats = rule[stat_name]
        _require(
            stats,
            ("present", "total", "frequency", "counting_unit"),
            f"{path}.{stat_name}",
        )
        present = int(
            _finite_nonnegative(stats["present"], f"{path}.{stat_name}.present")
        )
        total = int(_finite_nonnegative(stats["total"], f"{path}.{stat_name}.total"))
        frequency = float(stats["frequency"])
        if present > total or not 0 <= frequency <= 1:
            raise CatalogValidationError(f"{path}.{stat_name}: invalid counts")
        expected = present / total if total else 0.0
        if abs(expected - frequency) > 1e-9:
            raise CatalogValidationError(
                f"{path}.{stat_name}: frequency/count mismatch"
            )
        if stats["counting_unit"] != "independent_2d_identity":
            raise CatalogValidationError(
                f"{path}.{stat_name}.counting_unit: unsafe unit"
            )

    unknown_sources = set(rule["source_ids"]) - source_ids
    unknown_domains = set(rule["instrument_domain_ids"]) - domain_ids
    if unknown_sources:
        raise CatalogValidationError(
            f"{path}.source_ids: unknown {sorted(unknown_sources)}"
        )
    if unknown_domains:
        raise CatalogValidationError(
            f"{path}.instrument_domain_ids: unknown {sorted(unknown_domains)}"
        )

    negative = rule["negative_evidence"]
    _require(
        negative,
        ("absence_interpretable", "required_quality_fields", "limitations"),
        f"{path}.negative_evidence",
    )
    semantics = rule["output_semantics"]
    _require(
        semantics,
        (
            "evidence_stage",
            "exact_identity_confirmed",
            "standard_confirmed",
            "isomer_discriminating",
        ),
        f"{path}.output_semantics",
    )
    if semantics["evidence_stage"] not in EVIDENCE_STAGES:
        raise CatalogValidationError(f"{path}.output_semantics: invalid stage")
    if any(
        semantics[key] is not False
        for key in (
            "exact_identity_confirmed",
            "standard_confirmed",
            "isomer_discriminating",
        )
    ):
        raise CatalogValidationError(
            f"{path}.output_semantics: forbidden evidence upgrade"
        )


def validate_catalog(catalog: Mapping[str, Any]) -> None:
    """Validate the complete source/domain/rule catalog or raise."""
    _require(
        catalog, ("schema_version", "sources", "instrument_domains", "rules"), "catalog"
    )
    if catalog["schema_version"] != SCHEMA_VERSION:
        raise CatalogValidationError("catalog.schema_version: unsupported version")
    sources = catalog["sources"]
    domains = catalog["instrument_domains"]
    rules = catalog["rules"]
    if not all(isinstance(item, Mapping) for item in (*sources, *domains, *rules)):
        raise CatalogValidationError("catalog: arrays must contain objects")
    for index, source in enumerate(sources):
        validate_source(source, f"sources[{index}]")
    source_ids = [str(source["source_id"]) for source in sources]
    if len(source_ids) != len(set(source_ids)):
        raise CatalogValidationError("sources: duplicate source_id")
    for index, domain in enumerate(domains):
        validate_instrument_domain(domain, f"instrument_domains[{index}]")
    domain_ids = [str(domain["domain_id"]) for domain in domains]
    if len(domain_ids) != len(set(domain_ids)):
        raise CatalogValidationError("instrument_domains: duplicate domain_id")
    for domain in domains:
        unknown = set(domain["source_ids"]) - set(source_ids)
        if unknown:
            raise CatalogValidationError(
                f"instrument domain {domain['domain_id']}: unknown sources {sorted(unknown)}"
            )
    rule_ids: list[str] = []
    for index, rule in enumerate(rules):
        validate_rule(
            rule,
            source_ids=set(source_ids),
            domain_ids=set(domain_ids),
            path=f"rules[{index}]",
        )
        rule_ids.append(str(rule["rule_id"]))
    if len(rule_ids) != len(set(rule_ids)):
        raise CatalogValidationError("rules: duplicate rule_id")
