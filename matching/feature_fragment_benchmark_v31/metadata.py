"""Apply narrowly scoped, source-cited static instrument metadata derivations."""

from __future__ import annotations

import re
from collections.abc import Mapping
from copy import deepcopy
from typing import Any


_UNREPORTED = re.compile(r"(?:unreported|unknown|unresolved|not[_ -]?reported)", re.I)


def _present(value: Any) -> bool:
    return value is not None and bool(str(value).strip()) and not _UNREPORTED.search(
        str(value)
    )


def _path_get(value: Mapping[str, Any], path: str) -> Any:
    cursor: Any = value
    for part in path.split("."):
        if not isinstance(cursor, Mapping):
            return None
        cursor = cursor.get(part)
    return cursor


def _path_set(value: dict[str, Any], path: str, replacement: Any) -> None:
    parts = path.split(".")
    cursor: dict[str, Any] = value
    for part in parts[:-1]:
        child = cursor.get(part)
        if not isinstance(child, dict):
            child = {}
            cursor[part] = child
        cursor = child
    cursor[parts[-1]] = replacement


def _matches(acquisition: Mapping[str, Any], condition: Mapping[str, Any]) -> bool:
    if acquisition.get("instrument_vendor") != condition.get("instrument_vendor"):
        return False
    if acquisition.get("instrument_model") not in condition.get(
        "instrument_models", []
    ):
        return False
    return all(
        acquisition.get(field) == condition.get(field)
        for field in ("analyzer", "activation")
        if condition.get(field) is not None
    )


def _refresh_completeness(record: dict[str, Any]) -> None:
    acquisition = record.get("acquisition") or {}
    resolution = acquisition.get("instrument_domain_resolution") or {}
    if resolution.get("status") == "resolved":
        exact_fields = dict(resolution.get("exact_fields") or {})
        exact_fields["collision_cell_region"] = acquisition.get(
            "collision_cell_region"
        )
        resolution["exact_fields"] = exact_fields
        resolution["exact_condition_metadata_complete"] = all(
            _present(item) for item in exact_fields.values()
        )
        acquisition["instrument_domain_resolution"] = resolution
    completeness = record.get("metadata_completeness") or {}
    fields = dict(completeness.get("fields") or {})
    if "collision_cell_region" in fields:
        fields["collision_cell_region"] = _present(
            acquisition.get("collision_cell_region")
        )
        completeness["fields"] = fields
        completeness["missing_fields"] = sorted(
            field for field, present in fields.items() if not present
        )
        completeness["complete_field_count"] = sum(bool(v) for v in fields.values())
        completeness["field_count"] = len(fields)
        record["metadata_completeness"] = completeness


def apply_static_metadata_derivations(
    record: Mapping[str, Any], contract: Mapping[str, Any]
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    """Return an audit copy with permitted static fields filled.

    Quantitative eligibility flags are never edited here.  Callers recompute
    readiness from the cited metadata copy and preserve the source V3 record.
    """

    result = deepcopy(dict(record))
    protected = tuple(str(path) for path in contract["protected_per_scan_fields"])
    before = {path: deepcopy(_path_get(result, path)) for path in protected}
    allowed = set(str(path) for path in contract["permitted_vendor_manual_derived_fields"])
    audit: list[dict[str, Any]] = []
    acquisition = result.get("acquisition") or {}
    for rule in contract.get("derivation_rules", []):
        if not rule.get("enabled") or not _matches(acquisition, rule.get("match") or {}):
            continue
        field = str(rule.get("field") or "")
        if field not in allowed:
            raise ValueError(f"derivation field is not permitted: {field}")
        current = _path_get(result, field)
        if _present(current):
            audit.append(
                {
                    "derivation_id": rule["derivation_id"],
                    "field": field,
                    "status": "record_literal_retained",
                    "original_value": current,
                }
            )
            continue
        _path_set(result, field, rule["value"])
        audit.append(
            {
                "derivation_id": rule["derivation_id"],
                "field": field,
                "status": "derived",
                "original_value": current,
                "derived_value": rule["value"],
                "value_origin": rule["value_origin"],
                "source": deepcopy(rule["source"]),
            }
        )
    _refresh_completeness(result)
    after = {path: _path_get(result, path) for path in protected}
    changed_protected = [path for path in protected if before[path] != after[path]]
    if changed_protected:
        raise RuntimeError(
            "static derivation changed protected per-scan fields: "
            + ", ".join(changed_protected)
        )
    if audit:
        result["v31_metadata_derivation"] = {
            "contract_id": contract["contract_id"],
            "events": audit,
        }
    return result, audit


def common_quantitative_metadata_ready(record: Mapping[str, Any]) -> bool:
    """Check non-selection fields shared by LR+ and LR- eligibility."""

    source = record.get("source") or {}
    acquisition = record.get("acquisition") or {}
    resolution = acquisition.get("instrument_domain_resolution") or {}
    ce = acquisition.get("collision_energy") or {}
    eligibility = record.get("quantitative_eligibility") or {}
    return bool(
        (source.get("license") or {}).get("eligible_for_primary_analysis") is True
        and resolution.get("status") == "resolved"
        and resolution.get("exact_condition_metadata_complete") is True
        and ce.get("status") == "resolved"
        and ce.get("conversion_applied") is False
        and ce.get("unit")
        and ce.get("mode")
        and acquisition.get("polarity") in {"positive", "negative"}
        and acquisition.get("adduct")
        and isinstance(acquisition.get("charge"), int)
        and isinstance(acquisition.get("precursor_mz"), (int, float))
        and eligibility.get("range_contains_precursor") is True
    )


def recompute_v31_eligibility(record: Mapping[str, Any]) -> dict[str, bool]:
    common = common_quantitative_metadata_ready(record)
    eligibility = record.get("quantitative_eligibility") or {}
    acquisition = record.get("acquisition") or {}
    ce = acquisition.get("collision_energy") or {}
    isolation = acquisition.get("precursor_isolation") or {}
    composite = ce.get("mode") in {"MSE", "bbCID", "all_ions"} or acquisition.get(
        "activation_location"
    ) == "in_source"
    positive = bool(common and eligibility.get("precursor_selection_verified") is True)
    negative = bool(
        positive and isolation.get("status") == "verified" and not composite
    )
    return {
        "common_metadata_ready": common,
        "positive_lr_eligible": positive,
        "negative_lr_eligible": negative,
        "positive_ceiling_if_precursor_selection_is_verified": common,
        "negative_ceiling_if_isolation_is_verified": bool(common and not composite),
    }
