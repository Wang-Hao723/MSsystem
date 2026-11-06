"""Frozen-release exact-domain capacity audit before rule discovery."""

from __future__ import annotations

from collections import Counter, defaultdict
from collections.abc import Mapping, Sequence
from typing import Any


def _domain_row(
    domain_id: str,
    rows: Sequence[Mapping[str, Any]],
    thresholds: Mapping[str, Any],
) -> dict[str, Any]:
    primary_license_rows = [
        row
        for row in rows
        if (row.get("quantitative_eligibility") or {}).get("primary_license_eligible")
        is True
    ]
    exact_metadata_rows = [
        row
        for row in primary_license_rows
        if ((row.get("acquisition") or {}).get("instrument_domain_resolution") or {}).get(
            "exact_condition_metadata_complete"
        )
        is True
    ]
    eligible = [
        row
        for row in rows
        if (row.get("quantitative_eligibility") or {}).get("positive_lr_eligible")
        is True
    ]
    identities = {
        (row.get("independence") or {}).get("two_dimensional_identity_id")
        for row in eligible
        if (row.get("independence") or {}).get("two_dimensional_identity_id")
    }
    formulas = {
        (row.get("independence") or {}).get("formula_family_id")
        for row in eligible
        if (row.get("independence") or {}).get("formula_family_id")
    }
    sources = {
        (row.get("independence") or {}).get("source_origin_id")
        for row in eligible
        if (row.get("independence") or {}).get("source_origin_id")
    }
    laboratories = {
        (row.get("independence") or {}).get("laboratory_id")
        for row in eligible
        if (row.get("independence") or {}).get("laboratory_id")
    }
    identity_minimum = int(thresholds["minimum_target_identities"]) + int(
        thresholds["minimum_background_identities"]
    )
    formula_minimum = max(
        int(thresholds["minimum_target_formula_families"]),
        int(thresholds["minimum_background_formula_families"]),
    )
    source_minimum = max(
        int(thresholds["minimum_target_sources"]),
        int(thresholds["minimum_background_sources"]),
    )
    laboratory_minimum = max(
        int(thresholds["minimum_target_laboratories"]),
        int(thresholds["minimum_background_laboratories"]),
    )
    gates = {
        "disjoint_target_background_identity_capacity": len(identities)
        >= identity_minimum,
        "formula_family_capacity": len(formulas) >= formula_minimum,
        "source_capacity": len(sources) >= source_minimum,
        "laboratory_capacity": len(laboratories) >= laboratory_minimum,
    }
    status = (
        "capacity_sufficient_requires_rule_specific_check"
        if all(gates.values())
        else "infeasible_for_frozen_release"
    )
    reason_counts = Counter(
        reason
        for row in rows
        for reason in (
            (row.get("quantitative_eligibility") or {}).get("ineligible_reason_codes")
            or []
        )
    )
    return {
        "instrument_domain_id": domain_id,
        "status": status,
        "rule_specific_target_background_distribution_assessed": False,
        "audit_record_count": len(rows),
        "primary_license_record_count": len(primary_license_rows),
        "exact_condition_metadata_complete_record_count": len(exact_metadata_rows),
        "positive_lr_eligible_record_count": len(eligible),
        "maximum_counts_within_frozen_release": {
            "spectra": len(eligible),
            "two_dimensional_identities": len(identities),
            "formula_families": len(formulas),
            "origin_sources": len(sources),
            "origin_laboratories": len(laboratories),
        },
        "necessary_capacity_thresholds": {
            "disjoint_target_background_identities": identity_minimum,
            "formula_families": formula_minimum,
            "origin_sources": source_minimum,
            "origin_laboratories": laboratory_minimum,
        },
        "necessary_capacity_gates": gates,
        "failed_capacity_gates": sorted(name for name, passed in gates.items() if not passed),
        "ineligibility_reason_counts": dict(sorted(reason_counts.items())),
    }


def build_domain_feasibility_report(
    records: Sequence[Mapping[str, Any]],
    *,
    instrument_domains: Sequence[Mapping[str, Any]],
    thresholds: Mapping[str, Any],
    frozen_release_id: str,
    frozen_corpus_sha256: str,
) -> dict[str, Any]:
    """Report necessary capacity per exact domain; never claim permanence."""
    by_domain: dict[str, list[Mapping[str, Any]]] = defaultdict(list)
    unresolved: list[Mapping[str, Any]] = []
    for row in records:
        domain_id = (row.get("acquisition") or {}).get("instrument_domain_id")
        if domain_id:
            by_domain[str(domain_id)].append(row)
        else:
            unresolved.append(row)
    rows = [
        _domain_row(str(domain["domain_id"]), by_domain[str(domain["domain_id"])], thresholds)
        for domain in instrument_domains
    ]
    status_counts = Counter(row["status"] for row in rows)
    unresolved_reasons = Counter(
        reason
        for row in unresolved
        for reason in (
            (row.get("quantitative_eligibility") or {}).get("ineligible_reason_codes")
            or []
        )
    )
    return {
        "schema_version": "3.1.0",
        "report_id": "v15-feature-fragment-domain-feasibility-v3.1-20260827",
        "frozen_release_id": frozen_release_id,
        "frozen_corpus_sha256": frozen_corpus_sha256,
        "interpretation": "necessary aggregate capacity only; a domain that passes still requires rule-specific target and background checks",
        "permanent_infeasibility_claim_allowed": False,
        "domain_count": len(rows),
        "status_counts": dict(sorted(status_counts.items())),
        "domains": rows,
        "unresolved_domain_audit_bucket": {
            "status": "not_an_exact_domain",
            "record_count": len(unresolved),
            "positive_lr_eligible_record_count": sum(
                (row.get("quantitative_eligibility") or {}).get("positive_lr_eligible")
                is True
                for row in unresolved
            ),
            "ineligibility_reason_counts": dict(sorted(unresolved_reasons.items())),
        },
        "partial_pooling": {
            "may_rescue_failed_exact_domain": False,
            "evaluation_status": "not_reached_until_each_contributing_exact_domain_passes_all_primary_gates",
        },
    }
