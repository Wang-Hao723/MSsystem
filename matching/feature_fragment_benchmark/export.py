"""Fail-closed export gate for benchmark-tested feature-fragment rules."""

from __future__ import annotations

import hashlib
import json
from collections import Counter
from collections.abc import Mapping, Sequence
from typing import Any


def _rule_rejection_reasons(rule: Mapping[str, Any]) -> list[str]:
    reasons = []
    if rule.get("status") != "sealed_test_evaluated_shadow":
        reasons.append("sealed_test_not_completed")
    if rule.get("evidence_grade") != "empirically_validated_exact_domain":
        reasons.append("empirical_exact_domain_grade_missing")
    if rule.get("empirical_validation_gate_unmet"):
        reasons.append("empirical_validation_gate_unmet")
    if rule.get("evidence_stage_output") != "class_level_msms_support_only":
        reasons.append("invalid_output_evidence_stage")
    for field in (
        "identity_upgrade_allowed",
        "formula_upgrade_allowed",
        "isomer_resolution_allowed",
        "standard_confirmation_allowed",
    ):
        if rule.get(field) is not False:
            reasons.append(f"{field}_must_be_false")
    return reasons


def export_empirically_validated_rules(
    tested_rules: Sequence[Mapping[str, Any]],
) -> dict[str, Any]:
    """Return only exact-domain rules that pass every sealed-test gate."""
    accepted = []
    rejected = []
    reason_counts: Counter[str] = Counter()
    for rule in tested_rules:
        reasons = _rule_rejection_reasons(rule)
        if reasons:
            rejected.append({"rule_id": rule.get("rule_id"), "reasons": reasons})
            reason_counts.update(reasons)
        else:
            accepted.append(dict(rule))
    material = [
        {
            "rule_id": rule["rule_id"],
            "frozen_parameter_sha256": rule["frozen_parameter_sha256"],
            "evidence_grade": rule["evidence_grade"],
        }
        for rule in sorted(accepted, key=lambda row: row["rule_id"])
    ]
    return {
        "schema_version": "2.0.0",
        "runtime_influence": False,
        "exported_rule_count": len(accepted),
        "rejected_rule_count": len(rejected),
        "accepted_rules": accepted,
        "rejected_rules": rejected,
        "rejection_reason_counts": dict(sorted(reason_counts.items())),
        "accepted_rule_set_sha256": hashlib.sha256(
            json.dumps(material, sort_keys=True, separators=(",", ":")).encode("utf-8")
        ).hexdigest(),
        "integration_status": (
            "eligible_rules_require_separate_reviewed_runtime_adapter"
            if accepted
            else "no_rule_eligible_for_runtime_adapter"
        ),
    }
