from __future__ import annotations

import json
from pathlib import Path

from governance.authentic_anchor_registry import (
    audit_active_authentic_anchor_registry,
)


ROOT = Path(__file__).resolve().parents[1]


def _load(relative: str) -> dict:
    return json.loads((ROOT / relative).read_text(encoding="utf-8"))


def test_v14_authentic_reference_runtime_artifact_is_exactly_bound():
    audit = audit_active_authentic_anchor_registry(
        ROOT / "config/authentic_anchor_registry.json",
        project_root=ROOT,
    )

    assert audit["passed"] is True
    assert audit["entry_count"] == 77
    assert audit["actual_file_sha256"] == audit["expected_file_sha256"]
    assert audit["actual_content_sha256"] == audit["expected_content_sha256"]


def test_v14_effective_lab_gain_is_paired_additive_and_specific():
    result = _load(
        "eval/v14_validation/results/v14_evitwin_effective_evaluation.json"
    )
    display = result["headline"][
        "v14_core_to_v14_effective_candidate_display"
    ]
    safe = result["headline"]["v14_core_to_v14_effective_safe_scope"]

    assert display["total"] == 39
    assert display["before_successes"] == 25
    assert display["after_successes"] == 29
    assert display["gained"] == 4
    assert display["lost"] == 0
    assert safe["total"] == 45
    assert safe["before_successes"] == 36
    assert safe["after_successes"] == 38
    assert result["specificity_and_applicability"][
        "priority_predictions_outside_declared_composition"
    ] == 0
    assert result["headline"]["formal_level1_release_increment"] == 0


def test_v14_external_gate_abstains_cross_instrument_without_regression():
    result = _load(
        "eval/v14_validation/results/v14_evitwin_external_gate.json"
    )

    assert result["qe_orbitrap"]["cross_instrument_false_promotions"] == 0
    assert result["qe_orbitrap"]["v14_effective_candidate_recalled"] == 38
    assert result["qe_orbitrap"]["v14_effective_top1_2d_correct"] == 35
    assert result["mix33_qtof"]["qualified_authentic_priority_additions"] == 0
    assert result["mix33_qtof"]["v14_effective_top1_correct"] == 48
    assert result["no_regression_gate"]["passed"] is True
