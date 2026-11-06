from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]


def _json(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def test_committed_benchmark_receipt_recomputes_without_failures(tmp_path: Path) -> None:
    output = tmp_path / "receipt.json"
    subprocess.run(
        [
            sys.executable,
            str(ROOT / "scripts/audit_feature_fragment_benchmark_v2.py"),
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
        / "docs/audits/V15_FEATURE_FRAGMENT_BENCHMARK_V2_RECEIPT_20260826.json"
    )
    assert recomputed == committed
    assert recomputed["status"] == "pass"
    assert recomputed["failures"] == []
    assert all(recomputed["checks"].values())
    assert recomputed["rules"]["runtime_exported"] == 0
    assert recomputed["split"]["forbidden_overlap_total_recomputed"] == 0


def test_instrument_evidence_sources_resolve_and_ce_conversion_is_forbidden() -> None:
    ledger = _json(
        ROOT / "config/feature_fragment_benchmark_v2/instrument_evidence.json"
    )
    sources = {row["source_id"] for row in ledger["sources"]}
    assert len(sources) == len(ledger["sources"])
    assert {row["family_id"] for row in ledger["families"]} == {
        "thermo_qe_orbitrap_hcd",
        "agilent_qtof_collision_cell_cid",
        "sciex_qtof_q2_cid",
        "bruker_maxis_impact_qtof_cid",
        "waters_qtof_collision_cell_cid",
    }
    for family in ledger["families"]:
        assert family["source_ids"]
        assert set(family["source_ids"]) <= sources
    policy = ledger["record_level_policy"]
    assert policy["massbank_cid_is_not_automatically_psi_beam_type_cid"] is True
    assert policy["ce_unit_conversion_allowed"] is False
    assert policy["cross_vendor_equivalence_allowed"] is False
    assert policy["cross_model_equivalence_allowed"] is False
