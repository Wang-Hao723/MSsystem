from __future__ import annotations

import copy
import hashlib
import importlib
from pathlib import Path

import pytest

from matching.feature_fragments import (
    FeatureFragmentRuntimeGateError,
    load_runtime_gate,
    runtime_gate_audit,
    validate_runtime_gate,
)


search_module = importlib.import_module("matching.search")


def _gate():
    return load_runtime_gate()


def test_source_controlled_gate_exports_no_runtime_rules():
    gate = _gate()
    audit = runtime_gate_audit(gate)

    assert gate["runtime_evaluation_enabled"] is False
    assert gate["runtime_adapter_available"] is False
    assert gate["pilot_v1_runtime_allowed"] is False
    assert gate["exported_rule_count"] == 0
    assert audit["research_catalog_evaluated"] is False
    assert audit["benchmark_runtime_export_gate_sha256"] == (
        "ea0cb31811259a8a8da0d40fd50ee552a195f7dc147e0bb5df6896a02a916624"
    )


def test_source_controlled_gate_pins_the_sealed_export_artifact():
    root = Path(__file__).resolve().parents[1]
    gate = _gate()
    benchmark = gate["benchmark_runtime_export_gate"]
    artifact = root / benchmark["path"]

    assert artifact.is_file()
    assert hashlib.sha256(artifact.read_bytes()).hexdigest() == benchmark["sha256"]


@pytest.mark.parametrize(
    ("field", "value", "message"),
    [
        ("runtime_evaluation_enabled", True, "runtime_gate_unsafe_boolean"),
        ("pilot_v1_runtime_allowed", True, "runtime_gate_unsafe_boolean"),
        ("exported_rule_count", 1, "runtime_gate_unsafe_numeric_value"),
        ("deployment_state", "production", "runtime_gate_not_fail_closed"),
    ],
)
def test_gate_rejects_unsafe_state(field, value, message):
    gate = copy.deepcopy(_gate())
    gate[field] = value

    with pytest.raises(FeatureFragmentRuntimeGateError, match=message):
        validate_runtime_gate(gate)


def test_gate_loader_rejects_missing_or_invalid_json(tmp_path):
    missing = tmp_path / "missing.json"
    invalid = tmp_path / "invalid.json"
    invalid.write_text("{not-json", encoding="utf-8")

    with pytest.raises(FeatureFragmentRuntimeGateError, match="runtime_gate_unreadable"):
        load_runtime_gate(missing)
    with pytest.raises(FeatureFragmentRuntimeGateError, match="runtime_gate_unreadable"):
        load_runtime_gate(invalid)


def test_search_fails_closed_without_loading_research_catalog(monkeypatch):
    monkeypatch.setattr(
        search_module,
        "_feature_fragment_runtime_gate",
        lambda: (_ for _ in ()).throw(
            FeatureFragmentRuntimeGateError("runtime_gate_unreadable")
        ),
    )
    monkeypatch.setattr(
        search_module,
        "_feature_fragment_catalog",
        lambda: pytest.fail("research catalog must not be loaded"),
    )

    evidence = search_module._evaluate_feature_fragment_shadow(
        query_peaks=[[100.0, 10.0], [184.0733, 100.0], [250.0, 12.0]],
        query_pmz=493.57,
        candidate_class="Glycerophosphocholines",
        classification_reliable=True,
        classification_source="curated_massbank",
        query_polarity="POS",
        query_adduct="[M+H]+",
        query_adduct_reliable=True,
        query_charge=1,
        query_collision_energy={"value": 20.0, "unit": "eV"},
        query_collision_energy_mode="fixed",
        query_instrument_type="qtof",
        query_instrument="Agilent 6560 QTOF",
        query_activation="CID",
        query_instrument_domain_id=None,
        query_collision_cell_region=None,
        query_spectrum_representation="centroid",
        query_raw_peaks_complete=True,
        query_isolation_evidence={"status": "clean"},
    )

    assert evidence["status"] == "unavailable"
    assert evidence["reason"] == "feature_fragment_runtime_gate_invalid"
    assert evidence["runtime_influence"] is False
    assert evidence["runtime_export_gate"]["research_catalog_evaluated"] is False
