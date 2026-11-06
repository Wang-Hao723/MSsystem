"""Review-only V3 export; runtime export remains hard-coded off."""

from __future__ import annotations

import hashlib
import json
from collections.abc import Mapping, Sequence
from typing import Any


def export_review_candidates(rules: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    """Export passing candidates for human review without runtime activation."""
    candidates = [
        dict(rule)
        for rule in rules
        if rule.get("sealed_test_gate", {}).get("passed") is True
        and rule.get("evidence_grade") in {"G3", "G4"}
    ]
    material = [
        [row["rule_id"], row["frozen_parameter_sha256"], row["evidence_grade"]]
        for row in sorted(candidates, key=lambda item: item["rule_id"])
    ]
    return {
        "schema_version": "3.0.0",
        "status": "candidates_waiting_review" if candidates else "no_candidate_passed",
        "review_candidate_count": len(candidates),
        "review_candidates": candidates,
        "review_candidate_set_sha256": hashlib.sha256(
            json.dumps(material, sort_keys=True, separators=(",", ":")).encode("utf-8")
        ).hexdigest(),
        "runtime_exported_rule_count": 0,
        "runtime_catalog_modified": False,
        "runtime_gate_modified": False,
        "runtime_ranking_enabled": False,
        "score_contribution": 0.0,
        "next_state": "waiting_user_review",
    }
