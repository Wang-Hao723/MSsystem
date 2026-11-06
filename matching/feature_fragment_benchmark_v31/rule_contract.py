"""Strict validation for every future V3.1 feature-fragment rule record."""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any


def _present(value: Any) -> bool:
    return value is not None and bool(str(value).strip())


def _frequency_errors(name: str, value: Any) -> list[str]:
    if not isinstance(value, Mapping):
        return [f"{name}_missing"]
    errors = []
    events = value.get("event_count")
    denominator = value.get("identity_denominator")
    frequency = value.get("frequency")
    if not isinstance(events, int) or events < 0:
        errors.append(f"{name}_event_count_invalid")
    if not isinstance(denominator, int) or denominator <= 0:
        errors.append(f"{name}_identity_denominator_invalid")
    if not isinstance(frequency, (int, float)) or not 0 <= float(frequency) <= 1:
        errors.append(f"{name}_frequency_invalid")
    if (
        isinstance(events, int)
        and isinstance(denominator, int)
        and denominator > 0
        and isinstance(frequency, (int, float))
        and abs(events / denominator - float(frequency)) > 1e-12
    ):
        errors.append(f"{name}_frequency_not_equal_to_counts")
    return errors


def validate_rule_record(
    rule: Mapping[str, Any], schema: Mapping[str, Any]
) -> dict[str, Any]:
    errors: list[str] = []
    for field in schema["required_fields"]:
        if field not in rule:
            errors.append(f"required_field_missing:{field}")
    for field in schema["forbidden_fields"]:
        if field in rule:
            errors.append(f"forbidden_fixed_score_field:{field}")
    if rule.get("rule_kind") not in schema["rule_kinds"]:
        errors.append("rule_kind_invalid")
    if rule.get("evidence_grade") not in schema["evidence_grades"]:
        errors.append("evidence_grade_invalid")
    sources = rule.get("source_evidence")
    if not isinstance(sources, list) or not sources:
        errors.append("source_evidence_missing")
    else:
        for index, source in enumerate(sources):
            if not isinstance(source, Mapping):
                errors.append(f"source_evidence_invalid:{index}")
                continue
            for field in schema["source_evidence_minimum"]:
                if not _present(source.get(field)):
                    errors.append(f"source_evidence_{field}_missing:{index}")
    chemical = rule.get("chemical_scope") or {}
    if not isinstance(chemical, Mapping) or not chemical.get("leaf_classes"):
        errors.append("chemical_scope_leaf_classes_missing")
    if not rule.get("adducts"):
        errors.append("adduct_scope_missing")
    if not rule.get("polarities"):
        errors.append("polarity_scope_missing")
    instrument = rule.get("instrument_scope") or {}
    for field in schema["instrument_scope_fields"]:
        if not instrument.get(field):
            errors.append(f"instrument_scope_{field}_missing")
    ce = rule.get("collision_energy_scope") or {}
    for field in schema["collision_energy_scope_fields"]:
        if field not in ce:
            errors.append(f"collision_energy_scope_{field}_missing")
    if ce.get("cross_unit_conversion_applied") is not False:
        errors.append("collision_energy_cross_unit_conversion_forbidden")
    tolerance = rule.get("mass_tolerance") or {}
    for field in schema["mass_tolerance_fields"]:
        if field not in tolerance:
            errors.append(f"mass_tolerance_{field}_missing")
    if not (
        isinstance(tolerance.get("mz_tolerance_da"), (int, float))
        and float(tolerance["mz_tolerance_da"]) > 0
    ):
        errors.append("mass_tolerance_da_invalid")
    if not (
        isinstance(tolerance.get("mz_tolerance_ppm"), (int, float))
        and float(tolerance["mz_tolerance_ppm"]) > 0
    ):
        errors.append("mass_tolerance_ppm_invalid")
    errors.extend(_frequency_errors("target_frequency", rule.get("target_frequency")))
    errors.extend(
        _frequency_errors("background_frequency", rule.get("background_frequency"))
    )
    if not _present(rule.get("fragment_pathway_group_id")):
        errors.append("fragment_pathway_group_id_missing")
    if not isinstance(rule.get("positive_evidence_limit"), Mapping):
        errors.append("positive_evidence_limit_missing")
    if not isinstance(rule.get("negative_evidence_limit"), Mapping):
        errors.append("negative_evidence_limit_missing")
    validation = rule.get("validation") or {}
    for axis in (
        "identity_isolated",
        "formula_isolated",
        "source_isolated",
        "laboratory_isolated",
        "mirror_isolated",
        "raw_acquisition_isolated",
    ):
        if validation.get(axis) is not True:
            errors.append(f"validation_{axis}_not_true")
    for field, required in schema["claim_limits"].items():
        if rule.get(field) is not required:
            errors.append(f"claim_limit_violation:{field}")
    if rule.get("rule_kind") == "atomic_joint_fragment_event":
        if rule.get("joint_likelihood_ratio_fitted") is not True:
            errors.append("joint_likelihood_ratio_not_fitted")
        if rule.get("query_time_component_score_combination") not in {None, False}:
            errors.append("joint_rule_component_score_combination_forbidden")
    return {
        "status": "valid" if not errors else "invalid",
        "errors": sorted(set(errors)),
        "runtime_export_allowed": not errors
        and rule.get("evidence_grade") in {"G3", "G4"},
    }
