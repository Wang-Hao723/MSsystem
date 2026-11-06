from __future__ import annotations

import json
import hashlib
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]


def _json(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def test_committed_v31_audit_recomputes_exactly_and_runtime_stays_closed(
    tmp_path: Path,
) -> None:
    output = tmp_path / "receipt.json"
    subprocess.run(
        [
            sys.executable,
            str(ROOT / "scripts/audit_feature_fragment_benchmark_v31.py"),
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
        / "docs/audits/V15_FEATURE_FRAGMENT_BENCHMARK_V31_RECEIPT_20260827.json"
    )
    # The V16 path resolver changes executable provenance, never the frozen
    # scientific receipt. Verify present code bytes independently and compare
    # every remaining field to the original receipt without rewriting it.
    current_code = recomputed.pop("code_sha256")
    old_code = committed.pop("code_sha256")
    changed_code = {
        "scripts/audit_feature_fragment_benchmark_v31.py",
        "scripts/frozen_artifact_paths.py",
    }
    assert set(current_code) == set(old_code) | {"scripts/frozen_artifact_paths.py"}
    for name, digest in current_code.items():
        assert digest == hashlib.sha256((ROOT / name).read_bytes()).hexdigest()
        if name not in changed_code:
            assert digest == old_code[name]
    assert recomputed == committed
    assert recomputed["status"] == "pass"
    assert recomputed["failures"] == []
    assert all(recomputed["checks"].values())
    assert recomputed["runtime"]["runtime_exported_rule_count"] == 0
    assert recomputed["runtime"]["runtime_ranking_enabled"] is False
    assert recomputed["runtime"]["score_contribution"] == 0.0


def test_committed_full_manifest_counts_and_gates() -> None:
    receipt = _json(
        ROOT / "docs/audits/feature_fragment_benchmark_v31/benchmark_receipt.json"
    )
    assert receipt["corpus"]["record_count"] == 33928
    assert receipt["corpus"]["resolved_exact_domain_record_count"] == 19675
    assert receipt["corpus"]["two_dimensional_identity_count"] == 3402
    assert receipt["corpus"]["formula_family_count"] == 2747
    assert receipt["metadata_derivation"]["derived_record_count"] == 2083
    assert receipt["eligibility"]["lane_record_counts"]["v31_positive"] == 0
    assert receipt["eligibility"]["lane_record_counts"]["v31_negative"] == 0
    assert receipt["rules"]["runtime_exported"] == 0
