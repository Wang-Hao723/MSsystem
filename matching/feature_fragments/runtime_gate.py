"""Fail-closed runtime gate for V15 feature-fragment evidence.

The research catalog and its deterministic 30-spectrum fixture remain useful
for schema, parser, and regression tests.  They are not a production rule set.
Only rules exported by the sealed V2 benchmark may enable runtime evaluation;
the current accepted set is empty.
"""

from __future__ import annotations

import json
from collections.abc import Mapping
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_RUNTIME_GATE = ROOT / "config" / "feature_fragments" / "runtime_gate_v2.json"
SCHEMA_VERSION = "mssystem-v15-feature-fragment-runtime-gate/1.0"
EMPTY_RULE_SET_SHA256 = (
    "4f53cda18c2baa0c0354bb5f9a3ecbe5ed12ab4d8e11ba873c2f11161202b945"
)


class FeatureFragmentRuntimeGateError(ValueError):
    """Raised when the deployment gate is missing, malformed, or inconsistent."""


def validate_runtime_gate(gate: Mapping[str, Any]) -> dict[str, Any]:
    """Return a validated copy of the deployment gate.

    V15 currently supports only the sealed fail-closed state.  Enabling an
    adapter requires a separately reviewed implementation and is deliberately
    rejected here instead of silently loading the research catalog.
    """
    if not isinstance(gate, Mapping):
        raise FeatureFragmentRuntimeGateError("runtime_gate_not_an_object")
    required = {
        "schema_version",
        "deployment_state",
        "runtime_evaluation_enabled",
        "runtime_adapter_available",
        "pilot_v1_runtime_allowed",
        "pilot_v1_rule_count",
        "exported_rule_count",
        "accepted_rule_set_sha256",
        "benchmark_runtime_export_gate",
        "reason",
        "runtime_influence",
        "identity_rank_influence",
        "score_contribution",
    }
    missing = sorted(key for key in required if key not in gate)
    if missing:
        raise FeatureFragmentRuntimeGateError(
            f"runtime_gate_missing_fields:{','.join(missing)}"
        )
    if gate["schema_version"] != SCHEMA_VERSION:
        raise FeatureFragmentRuntimeGateError("runtime_gate_schema_mismatch")
    if gate["deployment_state"] != "fail_closed":
        raise FeatureFragmentRuntimeGateError("runtime_gate_not_fail_closed")
    if any(
        gate[key] is not False
        for key in (
            "runtime_evaluation_enabled",
            "runtime_adapter_available",
            "pilot_v1_runtime_allowed",
            "runtime_influence",
            "identity_rank_influence",
        )
    ):
        raise FeatureFragmentRuntimeGateError("runtime_gate_unsafe_boolean")
    try:
        pilot_count = int(gate["pilot_v1_rule_count"])
        exported_count = int(gate["exported_rule_count"])
        score_contribution = float(gate["score_contribution"])
    except (TypeError, ValueError) as exc:
        raise FeatureFragmentRuntimeGateError("runtime_gate_invalid_numeric_field") from exc
    if pilot_count < 0 or exported_count != 0 or score_contribution != 0.0:
        raise FeatureFragmentRuntimeGateError("runtime_gate_unsafe_numeric_value")
    if gate["accepted_rule_set_sha256"] != EMPTY_RULE_SET_SHA256:
        raise FeatureFragmentRuntimeGateError("runtime_gate_nonempty_rule_set_hash")
    benchmark = gate["benchmark_runtime_export_gate"]
    if not isinstance(benchmark, Mapping):
        raise FeatureFragmentRuntimeGateError("runtime_gate_benchmark_not_an_object")
    if not str(benchmark.get("path") or "").strip():
        raise FeatureFragmentRuntimeGateError("runtime_gate_benchmark_path_missing")
    digest = str(benchmark.get("sha256") or "").strip().lower()
    if len(digest) != 64 or any(char not in "0123456789abcdef" for char in digest):
        raise FeatureFragmentRuntimeGateError("runtime_gate_benchmark_sha256_invalid")
    if gate["reason"] != "no_feature_fragment_rule_passed_sealed_runtime_gate":
        raise FeatureFragmentRuntimeGateError("runtime_gate_reason_mismatch")
    return dict(gate)


def load_runtime_gate(path: Path | str | None = None) -> dict[str, Any]:
    """Load the source-controlled runtime gate without reading research data."""
    gate_path = Path(path) if path is not None else DEFAULT_RUNTIME_GATE
    try:
        with gate_path.open("r", encoding="utf-8") as handle:
            payload = json.load(handle)
    except (OSError, json.JSONDecodeError) as exc:
        raise FeatureFragmentRuntimeGateError("runtime_gate_unreadable") from exc
    return validate_runtime_gate(payload)


def runtime_gate_audit(gate: Mapping[str, Any]) -> dict[str, Any]:
    """Expose non-ranking deployment metadata in every unavailable result."""
    validated = validate_runtime_gate(gate)
    benchmark = validated["benchmark_runtime_export_gate"]
    return {
        "schema_version": validated["schema_version"],
        "deployment_state": validated["deployment_state"],
        "runtime_evaluation_enabled": False,
        "runtime_adapter_available": False,
        "pilot_v1_runtime_allowed": False,
        "pilot_v1_rule_count": int(validated["pilot_v1_rule_count"]),
        "exported_rule_count": 0,
        "accepted_rule_set_sha256": validated["accepted_rule_set_sha256"],
        "benchmark_runtime_export_gate_path": benchmark["path"],
        "benchmark_runtime_export_gate_sha256": benchmark["sha256"],
        "research_catalog_evaluated": False,
    }
