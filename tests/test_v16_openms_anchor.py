from __future__ import annotations

from copy import deepcopy

from matching.v16_openms_anchor import (
    compact_openms_validation_audit,
    validate_openms_adduct_group_manifest,
)


def _feature(
    feature_id: str,
    adduct: str,
    neutral_mass: float,
    rt_seconds: float,
    *,
    group_id: str = "AG000001",
) -> dict:
    return {
        "feature_id": feature_id,
        "adduct_group_id": group_id,
        "adduct": adduct,
        "neutral_mass": neutral_mass,
        "rt_seconds": rt_seconds,
    }


def _manifest(tmp_path, polarity: str = "POS") -> dict:
    feature_xml = tmp_path / f"{polarity}.featureXML"
    feature_xml.write_text("<featureMap/>", encoding="utf-8")
    features = [
        _feature("LF1", "[M+H]+", 300.0000, 60.0),
        _feature("LF2", "[M+Na]+", 300.0008, 63.0),
    ]
    return {
        "status": "ok",
        "backend": "PyOpenMS 3.5.0",
        "pipeline": [
            "per-spectrum PeakPickerHiRes/PeakPickerIterative for profile MS1 only",
            "MassTraceDetection",
            "ElutionPeakDetection",
            "FeatureFindingMetabo",
            "MetaboliteFeatureDeconvolution",
        ],
        "output_features": str(feature_xml),
        "feature_count": len(features),
        "features": features,
        "adduct_deconvolution_status": "ok",
        "parameters": {"group_adducts": True, "ion_mode": polarity},
    }


def test_valid_group_emits_feature_bound_candidate_independent_evidence(tmp_path):
    manifest = _manifest(tmp_path)
    validation = validate_openms_adduct_group_manifest(
        manifest,
        task_polarity="POS",
        manifest_sha256="a" * 64,
        feature_xml_sha256="b" * 64,
        feature_xml_path=manifest["output_features"],
    )

    assert validation["status"] == "validated"
    assert validation["eligible_group_count"] == 1
    assert validation["eligible_feature_count"] == 2
    assert validation["candidate_identity_used"] is False
    assert validation["production_runtime_influence"] is False
    assert set(validation["evidence_by_feature_id"]) == {"LF1", "LF2"}
    evidence = validation["evidence_by_feature_id"]["LF1"]
    assert evidence["status"] == "validated"
    assert evidence["neutral_mass_consistent"] is True
    assert evidence["member_feature_ids"] == ["LF1", "LF2"]
    assert evidence["member_feature_count"] == 2
    assert evidence["candidate_identity_used"] is False

    audit = compact_openms_validation_audit(validation)
    assert "evidence_by_feature_id" not in audit
    assert audit["eligible_group_count"] == 1


def test_mass_rt_adduct_and_polarity_fail_closed_per_group(tmp_path):
    manifest = _manifest(tmp_path)
    manifest["features"].extend(
        [
            _feature("LF3", "[M+H]+", 400.0, 100.0, group_id="AG-BADMASS"),
            _feature("LF4", "[M+Na]+", 400.02, 101.0, group_id="AG-BADMASS"),
            _feature("LF5", "[M+H]+", 500.0, 200.0, group_id="AG-BADRT"),
            _feature("LF6", "[M+Na]+", 500.0005, 220.1, group_id="AG-BADRT"),
            _feature("LF7", "[M+H]+", 600.0, 300.0, group_id="AG-SAME"),
            _feature("LF8", "[M+H]+", 600.0005, 301.0, group_id="AG-SAME"),
            _feature("LF9", "[M-H]-", 700.0, 400.0, group_id="AG-SIGN"),
            _feature("LF10", "[M+Na]+", 700.0005, 401.0, group_id="AG-SIGN"),
        ]
    )
    manifest["feature_count"] = len(manifest["features"])

    validation = validate_openms_adduct_group_manifest(
        manifest,
        task_polarity="POS",
        feature_xml_path=manifest["output_features"],
    )

    assert validation["eligible_group_count"] == 1
    assert validation["rejected_group_count"] == 4
    counts = validation["rejected_group_reason_counts"]
    assert counts["neutral_mass_spread_exceeds_limit"] == 1
    assert counts["retention_time_spread_exceeds_limit"] == 1
    assert counts["insufficient_distinct_adducts"] == 1
    assert counts["adduct_sign_mismatch"] == 1
    assert set(validation["evidence_by_feature_id"]) == {"LF1", "LF2"}


def test_manifest_provenance_and_pipeline_mismatches_reject_all_groups(tmp_path):
    manifest = _manifest(tmp_path)
    manifest["output_features"] = str(tmp_path / "different.featureXML")
    manifest["parameters"]["ion_mode"] = "NEG"
    manifest["pipeline"].remove("MetaboliteFeatureDeconvolution")

    validation = validate_openms_adduct_group_manifest(
        manifest,
        task_polarity="POS",
        feature_xml_path=tmp_path / "POS.featureXML",
    )

    assert validation["status"] == "rejected_manifest"
    assert validation["evidence_by_feature_id"] == {}
    assert set(validation["manifest_rejection_reasons"]) == {
        "feature_xml_binding_mismatch",
        "manifest_ion_mode_mismatch",
        "required_openms_pipeline_stage_missing",
    }


def test_duplicate_feature_ids_and_declared_count_mismatch_reject_manifest(tmp_path):
    manifest = _manifest(tmp_path)
    duplicate = deepcopy(manifest["features"][0])
    manifest["features"].append(duplicate)

    validation = validate_openms_adduct_group_manifest(
        manifest,
        task_polarity="POS",
        feature_xml_path=manifest["output_features"],
    )

    assert validation["status"] == "rejected_manifest"
    assert "duplicate_feature_id" in validation["manifest_rejection_reasons"]
    assert "feature_count_mismatch" in validation["manifest_rejection_reasons"]
