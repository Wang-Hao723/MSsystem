from __future__ import annotations

from pathlib import Path

import pytest
from fastapi.testclient import TestClient
from pydantic import ValidationError

from config import settings
from governance.module_registry import (
    ModuleDefinition,
    audit_module_registry,
    load_module_registry,
)
from schemas.evidence import ModuleState
from web import server


def test_shadow_module_cannot_declare_ranking_influence():
    with pytest.raises(ValidationError):
        ModuleDefinition(
            module_id="bad-shadow",
            display_name="Bad",
            state=ModuleState.SHADOW,
            scientific_role="test",
            allowed_effects=("ranking",),
            forbidden_effects=(),
        )


def test_module_registry_matches_current_runtime_and_artifacts():
    root = Path(__file__).resolve().parents[1]
    registry = load_module_registry(root / "config" / "module_governance.json")
    report = audit_module_registry(
        registry,
        project_root=root,
        runtime_settings=settings,
    )
    assert report["passed"] is True
    assert report["violation_count"] == 0
    states = {module.module_id: module.state for module in registry.modules}
    assert states["v11_flagship_chem_entity_ranker"] == ModuleState.PRODUCTION
    assert states["v12_flagship_evirank_mv"] == ModuleState.SHADOW
    assert states["v12_evirank_candidate_recall"] == ModuleState.SHADOW
    assert states["recap_formula_conditioned_adduct"] == ModuleState.SHADOW
    assert states["isopair_explicit_adduct_ranking"] == ModuleState.SHADOW
    assert states["isopair_selective_high_veto"] == ModuleState.QUALIFIED
    assert states["candidate_pathway_alignment"] == ModuleState.SHADOW
    assert states["candidate_independent_lockmass_recalibration"] == ModuleState.SHADOW
    assert states["structure_classifier_ml"] == ModuleState.OFF
    assert states["project_target_decoy_calibration"] == ModuleState.QUALIFIED
    assert states["eviacquire_active_evidence"] == ModuleState.SHADOW
    assert states["evograph_fragmentation"] == ModuleState.SHADOW
    assert states["v16_dual_polarity_evidence_graph"] == ModuleState.SHADOW
    assert states["formula_prediction_fiddle"] == ModuleState.SHADOW
    assert states["formula_conditioned_library_enumeration"] == ModuleState.SHADOW
    assert states["evidiscrim_candidate_contrast"] == ModuleState.SHADOW
    assert states["forward_spectrum_fiora"] == ModuleState.SHADOW
    assert states["forward_spectrum_iceberg"] == ModuleState.OFF
    assert states["de_novo_structure_diffms"] == ModuleState.OFF


def test_cross_occurrence_ranker_is_explicitly_qualified():
    root = Path(__file__).resolve().parents[1]
    registry = load_module_registry(root / "config" / "module_governance.json")
    states = {module.module_id: module.state for module in registry.modules}
    assert states["v16_cross_occurrence_consensus_ranker"] == ModuleState.QUALIFIED


def test_module_governance_api_exposes_validated_registry():
    response = TestClient(server.app).get("/api/modules/governance")

    assert response.status_code == 200
    payload = response.json()
    assert payload["audit"]["passed"] is True
    assert payload["audit"]["violation_count"] == 0
    assert len(payload["registry"]["modules"]) == 32
    assert any(
        module.get("module_id") == "evitwin_ms"
        and module.get("state") == "qualified"
        for module in payload["registry"]["modules"]
    )
    flagship = next(
        module
        for module in payload["registry"]["modules"]
        if module["module_id"] == "v11_flagship_chem_entity_ranker"
    )
    assert flagship["state"] == "production"
    assert "entity_identity_selection" in flagship["allowed_effects"]
    challenger = next(
        module
        for module in payload["registry"]["modules"]
        if module["module_id"] == "v12_flagship_evirank_mv"
    )
    assert challenger["state"] == "shadow"
    assert "production_ranking" in challenger["forbidden_effects"]
    assert payload["claim_compiler"] == "EviGate 2.0"
