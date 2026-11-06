from scripts.audit_feature_fragment_rules import build_audit


def test_legacy_audit_is_row_level_and_keeps_source_types_separate():
    audit = build_audit()
    summary = audit["summary"]

    assert summary["row_level_record_count"] == 2030
    assert summary["counts_by_system"] == {
        "kb_v3.empirical_rules": 379,
        "kb_v3.empirical_rules_candidates_v2": 1156,
        "matching.fragment_feature_library.COMPOUND_CLASS_PROFILES": 159,
        "matching.fragmentation_rules.COMMON_NEUTRAL_LOSSES": 35,
        "matching.fragmentation_rules.DIAGNOSTIC_IONS": 44,
        "matching.fragmentation_rules.FRAGMENTATION_RULES": 22,
        "matching.fragmentation_rules_curation": 235,
    }
    assert summary["records_without_atomic_source"] == 329
    assert summary["existing_static_rules_with_official_vendor_conditions"] == 0
    assert summary["production_rules_with_target_and_background_frequency"] == 0
    assert summary["production_rules_with_independent_fragmentation_groups"] == 0
    assert (
        summary["provenance_class_occurrences"]["real_experimental_spectrum_statistics"]
        == 1535
    )
    assert summary["provenance_class_occurrences"]["computational_inference"] == 1156


def test_legacy_audit_marks_fixed_gates_and_evidence_stage_errors():
    audit = build_audit()
    constants = {item["path"]: item for item in audit["fixed_constants"]}
    assert constants["matching/evidence_governance.py:38"]["values"] == {
        "expert_release_gate": 0.40
    }
    assert all(
        item["provenance_class"] == "unsourced_constant" for item in constants.values()
    )
    defects = {(item["severity"], item["path"]) for item in audit["known_defects"]}
    assert ("P0", "matching/expected_compounds.py:1828") in defects
    assert ("P0", "matching/evidence_governance.py:139") in defects
    expected_system = next(
        item
        for item in audit["knowledge_systems"]
        if item["system"] == "expected_compounds.experimental_library"
    )
    assert "not same-method standard confirmation" in expected_system["audit_result"]
