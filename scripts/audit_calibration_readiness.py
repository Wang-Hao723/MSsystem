#!/usr/bin/env python3
"""Verify registered FDR evidence and publish P3 production-readiness truth."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
from datetime import datetime, timezone
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from calibration.project import load_calibration_method, method_semantic_sha256


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=PROJECT_ROOT / "docs" / "audits" / "V9_P3_CALIBRATION_READINESS.json",
    )
    args = parser.parse_args()
    registry_path = PROJECT_ROOT / "config" / "calibration_registry.json"
    registry = json.loads(registry_path.read_text(encoding="utf-8"))
    method = load_calibration_method(PROJECT_ROOT / "config" / "calibration_method.json")
    entries = []
    hashes_ok = True
    for entry in registry.get("entries") or []:
        path = PROJECT_ROOT / entry["artifact"]
        actual = hashlib.sha256(path.read_bytes()).hexdigest() if path.is_file() else None
        hash_ok = actual == entry.get("sha256")
        hashes_ok &= hash_ok
        entries.append(
            {
                **entry,
                "artifact_exists": path.is_file(),
                "actual_sha256": actual,
                "hash_verified": hash_ok,
            }
        )
    golden = json.loads(
        (PROJECT_ROOT / "config" / "golden_task_sets.json").read_text(encoding="utf-8")
    )
    golden_task_count = sum(
        len(group.get("task_ids") or [])
        for group in golden.get("cohorts") or []
    )
    payload = {
        "schema_version": "1.0.0",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "calibration_method": {
            "method_id": method.method_id,
            "method_semantic_sha256": method_semantic_sha256(method),
            "score_field": method.score_field,
            "unit_levels": [item.value for item in method.allowed_unit_levels],
            "claim_scopes": [item.value for item in method.allowed_claim_scopes],
            "predeclared_stratum_fallback": method.stratum_fallback,
            "tie_policy": method.winner_tie_policy,
            "pseudocount": method.pseudocount,
        },
        "registered_evidence": entries,
        "registry_hashes_verified": hashes_ok,
        "production_calibration_available": bool(
            registry.get("production_calibration_id")
        ),
        "golden_tasks_audited": golden_task_count,
        "golden_tasks_with_assignable_project_qvalues": 0,
        "qvalues_withheld_reason_codes": [
            "no_current_project_target_decoy_competition_artifact",
            "registered_mix33_evidence_is_development_informed",
            "registered_artifacts_do_not_contain_assignable_entity_level_calls",
            "registered_score_is_cosine_not_current_rank_score",
        ],
        "probability_calibration": {
            "available": False,
            "brier_score": None,
            "expected_calibration_error": None,
            "reason": "no_independent_per-call_probability_and_binary-label_set",
            "qvalues_treated_as_probabilities": False,
        },
        "policy_assertions": registry.get("policy") or {},
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.output.with_name(f".{args.output.name}.tmp-{os.getpid()}")
    with temporary.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, args.output)
    print(json.dumps(payload, ensure_ascii=False, indent=2))
    return 0 if hashes_ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
