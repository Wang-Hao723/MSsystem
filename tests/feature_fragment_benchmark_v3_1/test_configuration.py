from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest

from matching.feature_fragment_benchmark_v3_1.rules import discover_rules_with_policy


ROOT = Path(__file__).resolve().parents[2]


def _json(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def test_amendment_freezes_real_v3_inputs_and_runtime_boundary() -> None:
    amendment = _json(
        ROOT / "config/feature_fragment_benchmark_v3_1/amendment.json"
    )
    for path_key, hash_key in (
        ("machine_preregistration_path", "machine_preregistration_sha256"),
        ("human_preregistration_path", "human_preregistration_sha256"),
        ("phase1_receipt_path", "phase1_receipt_sha256"),
        ("benchmark_receipt_path", "benchmark_receipt_sha256"),
    ):
        path = ROOT / amendment["base_release"][path_key]
        assert hashlib.sha256(path.read_bytes()).hexdigest() == amendment[
            "base_release"
        ][hash_key]
    assert amendment["scope"]["runtime_effect"] == "none_fail_closed"
    assert amendment["scope"]["phase_2_allowed"] is False
    assert amendment["scope"]["standards_purchase_or_acquisition_allowed"] is False
    assert amendment["collision_energy"]["nce_to_ev_conversion_allowed"] is False


def test_qc_policy_is_external_and_mismatch_fails_closed(qc_policy) -> None:
    assert discover_rules_with_policy([], qc_policy=qc_policy) == []
    changed = json.loads(json.dumps(qc_policy))
    changed["peak_presence"]["minimum_relative_intensity"] = 0.02
    with pytest.raises(ValueError, match="qc_policy_mismatch"):
        discover_rules_with_policy([], qc_policy=changed)


def test_mona_is_not_integrated_or_counted_as_independent_source() -> None:
    registry = _json(
        ROOT
        / "config/feature_fragment_benchmark_v3_1/public_source_amendments.json"
    )
    mona = registry["sources"][0]
    assert mona["integration_status"] == "screened_not_integrated"
    assert mona["default_site_license_is_per_spectrum_proof"] is False
    assert mona["aggregator_counts_as_new_source_or_laboratory"] is False
    assert registry["new_public_spectra_integrated_in_this_amendment"] == 0
