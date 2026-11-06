from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]


def _json(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def test_committed_v3_phase1_receipt_recomputes_and_runtime_stays_closed(
    tmp_path: Path,
) -> None:
    output = tmp_path / "receipt.json"
    subprocess.run(
        [
            sys.executable,
            str(ROOT / "scripts/audit_feature_fragment_benchmark_v3.py"),
            "--output",
            str(output),
        ],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    recomputed = _json(output)
    committed = _json(
        ROOT
        / "docs/audits/V15_FEATURE_FRAGMENT_BENCHMARK_V3_RECEIPT_20260827.json"
    )
    assert recomputed == committed
    assert recomputed["status"] == "pass"
    assert recomputed["failures"] == []
    assert all(recomputed["checks"].values())
    assert recomputed["runtime"]["runtime_exported_rule_count"] == 0
    assert recomputed["runtime"]["runtime_ranking_enabled"] is False
    assert recomputed["runtime"]["score_contribution"] == 0.0
    assert recomputed["next_state"] == "waiting_user_review"


def test_v2_artifact_hashes_remain_identical() -> None:
    receipt = _json(
        ROOT
        / "docs/audits/V15_FEATURE_FRAGMENT_BENCHMARK_V3_RECEIPT_20260827.json"
    )
    assert receipt["checks"]["v2_artifacts_unchanged"] is True
    assert receipt["v2_artifact_sha256"]["record_manifest.jsonl.gz"] == (
        "8ab3a8d39884cee0b637039442ec8dfdc2cac9f492a5c1811b5073823185e50f"
    )
