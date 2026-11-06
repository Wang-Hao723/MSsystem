from __future__ import annotations

import json
from pathlib import Path

from matching.feature_fragment_benchmark.corpus import file_sha256
from scripts.audit_feature_fragment_benchmark_v3_1 import audit


ROOT = Path(__file__).resolve().parents[2]


def _json(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def test_committed_v31_phase1_receipt_recomputes() -> None:
    artifact_dir = ROOT / "docs/audits/feature_fragment_benchmark_v3_1"
    committed = _json(
        ROOT
        / "docs/audits/V15_FEATURE_FRAGMENT_BENCHMARK_V3_1_RECEIPT_20260827.json"
    )
    recomputed = audit(artifact_dir)
    assert committed["status"] == recomputed["status"] == "pass"
    assert committed["failures"] == recomputed["failures"] == []
    assert all(recomputed["checks"].values())
    assert recomputed["benchmark_state"] == "infeasible_for_frozen_release"
    assert recomputed["next_state"] == "waiting_user_review"


def test_runtime_gate_and_v3_freeze_remain_unchanged() -> None:
    receipt = _json(
        ROOT
        / "docs/audits/V15_FEATURE_FRAGMENT_BENCHMARK_V3_1_RECEIPT_20260827.json"
    )
    runtime_path = ROOT / "config/feature_fragments/runtime_gate_v2.json"
    assert file_sha256(runtime_path) == receipt["runtime"]["runtime_gate_sha256"]
    assert receipt["runtime"]["runtime_exported_rule_count"] == 0
    assert receipt["runtime"]["runtime_ranking_enabled"] is False
    assert receipt["runtime"]["score_contribution"] == 0.0
    assert receipt["rules"]["identity_upgrade_allowed"] is False
    assert receipt["rules"]["formula_upgrade_allowed"] is False
    assert receipt["rules"]["isomer_resolution_allowed"] is False
    assert receipt["rules"]["standard_confirmation_allowed"] is False
    assert all(
        row["matches"] for row in receipt["base_v3_named_hash_checks"].values()
    )
    assert all(
        row["matches"] for row in receipt["base_v3_artifact_hash_checks"].values()
    )
