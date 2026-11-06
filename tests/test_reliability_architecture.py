from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import re

import pytest

from pipeline import compound_conclusions as compound_projection
from pipeline.compound_conclusions import build_compound_conclusions
from pipeline.reliability_annotations import (
    background_evidence,
    ce_comparability,
    fragmentation_annotation,
)
from web.server import AnalysisParams


ROOT = Path(__file__).resolve().parents[1]


def _match(name: str = "Library candidate") -> dict:
    return {
        "compound_name": name,
        "formula": "C10H12O3",
        "smiles": "CC1=CC=CC=C1O",
        "inchikey": "AAAAAAAAAAAAAA-UHFFFAOYSA-N",
        "confidence": "Medium",
        "exact_structure_confidence": "Medium",
        "two_dimensional_connectivity_confidence": "Medium",
        "identification_scope": "exact_2d_connectivity",
        "rank_score": 0.61,
        "final_score": 0.61,
        "dataset_source": "formal-library",
    }


def _row(scan: str, feature: str, rt: float, match: dict | None) -> dict:
    return {
        "spectrum_id": scan,
        "precursor_mz": 330.1696,
        "retention_time": rt,
        "retention_time_unit": "minute",
        "ion_mode": "POS",
        "lcms_feature_id": feature,
        "precursor_feature_id": f"PF-{scan}",
        "collision_energy": 15.0,
        "collision_energy_unit": "eV",
        "collision_energy_mode": "fixed",
        "match": match,
    }


def test_p0_default_is_target_free_and_runtime_provider_is_absent():
    assert AnalysisParams().enable_expected_compounds is False
    assert not (ROOT / "web/qft_expected_dataset.py").exists()
    product_files = [
        ROOT / "web/server.py",
        ROOT / "matching/expected_compounds.py",
        ROOT / "web/static/js/upload.js",
    ]
    assert all(
        "qft_expected_dataset" not in path.read_text(encoding="utf-8")
        for path in product_files
    )
    html = (ROOT / "web/static/index.html").read_text(encoding="utf-8")
    checkbox = html.split('id="enableExpectedCompounds"', 1)[1].split(">", 1)[0]
    assert re.search(r"\schecked(?:\s|=|$)", checkbox) is None


def test_r1_groups_repeated_dda_by_feature_and_keeps_distinct_rt_features():
    match = _match()
    payload = build_compound_conclusions(
        "anonymous-run",
        [
            _row("scan-1", "LF-1", 4.10, match),
            _row("scan-2", "LF-1", 4.12, match),
            _row("scan-3", "LF-2", 7.80, match),
        ],
    )
    assert payload["metadata"]["unique_report_rows"] == 2
    by_feature = {item["feature_id"]: item for item in payload["compounds"]}
    assert by_feature["lcms:LF-1"]["spectra_count"] == 2
    assert by_feature["lcms:LF-2"]["spectra_count"] == 1
    assert all(item["display"]["default_visible"] for item in by_feature.values())


def test_r2_ce_annotation_has_zero_rank_and_governance_influence():
    same = ce_comparability(
        {
            "collision_energy": 15,
            "collision_energy_unit": "eV",
            "collision_energy_mode": "fixed",
            "collision_energy_quantities": [
                {"term": "collision energy", "value": 15, "unit": "eV"}
            ],
        },
        {
            "lib_collision_energy": 15,
            "lib_collision_energy_unit": "eV",
            "lib_collision_energy_mode": "fixed",
        },
    )
    assert same["status"] == "same_recorded_condition"
    assert same["score_multiplier"] == 1.0
    assert same["ranking_influence"] is False
    assert same["governance_influence"] is False
    unknown = ce_comparability(
        {"collision_energy": 15, "collision_energy_unit": "eV"}, {}
    )
    assert unknown["status"] == "reference_unknown"


def test_r3_no_blank_retains_raw_evidence_without_downweighting():
    annotation = background_evidence({"blank_used": False})
    assert annotation["B3_blank_confirmed"] == "unknown"
    assert annotation["action"] == "retain"
    assert annotation["new_background_excluded_peak_count"] == 0
    assert annotation["new_background_soft_downweighted_peak_count"] == 0
    assert annotation["score_multiplier"] == 1.0
    assert annotation["raw_peaks_preserved"] is True


def test_r4_unmatched_feature_remains_visible_and_counted_in_partition():
    payload = build_compound_conclusions(
        "anonymous-run",
        [_row("scan-0", "LF-0", 2.5, None)],
    )
    item = payload["compounds"][0]
    projection = payload["metadata"]["feature_result_projection"]
    assert item["conclusion"]["key"] == "unannotated"
    assert item["display"] == {
        "tier": "unannotated",
        "primary_candidate_id": None,
        "default_visible": True,
        "warnings": ["no_library_candidate"],
    }
    assert projection["n_detected_ion_features"] is None
    assert projection["n_unannotated_features"] == 1
    assert projection["independent_inventory_available"] is False
    assert projection["partition_satisfied"] is None
    assert fragmentation_annotation({"query_peaks": []})["status"] == (
        "insufficient_fragments"
    )


def test_p3_release_review_and_unannotated_counts_follow_real_claim_decision(
    tmp_path: Path,
):
    released = _match("Released library identity")
    blocked = _match("Visible review candidate")
    blocked["confidence"] = "High"
    blocked["evidence_governance"] = {
        "policy_id": "evigate-ms-locked-v1-2026-08-02",
        "claim_status": "insufficient_for_release",
        "release_eligible_before_replication": False,
        "release_blockers": ["insufficient_independent_injection_groups"],
        "ranking_changed": False,
        "confidence_field_overwritten": False,
    }
    claim_path = tmp_path / "claim-decisions.json"
    claim_path.write_text(
        json.dumps({
            "decisions": [{
                "spectrum_id": "scan-released",
                "decision": {
                    "released": True,
                    "policy_id": "evigate-ms-locked-v1-2026-08-02",
                },
            }]
        }),
        encoding="utf-8",
    )

    payload = build_compound_conclusions(
        "p3-partition",
        [
            _row("scan-released", "LF-RELEASED", 1.0, released),
            _row("scan-blocked", "LF-REVIEW", 2.0, blocked),
        ],
        claim_decisions_path=claim_path,
        feature_inventory=[
            {"feature_id": "LF-RELEASED", "mz": 330.1696, "rt_minutes": 1.0},
            {"feature_id": "LF-REVIEW", "mz": 330.1696, "rt_minutes": 2.0},
            {"feature_id": "LF-UNANNOTATED", "mz": 342.1700, "rt_minutes": 3.0},
        ],
    )

    by_feature = {item["feature_id"]: item for item in payload["compounds"]}
    projection = payload["metadata"]["feature_result_projection"]
    assert by_feature["lcms:LF-RELEASED"]["display"]["tier"] == "library_released"
    assert by_feature["lcms:LF-REVIEW"]["display"] == {
        "tier": "library_review",
        "primary_candidate_id": "AAAAAAAAAAAAAA-UHFFFAOYSA-N",
        "default_visible": True,
        "warnings": [],
    }
    assert by_feature["lcms:LF-REVIEW"]["representative"]["match"][
        "evidence_governance"
    ]["claim_status"] == "insufficient_for_release"
    assert projection["n_detected_ion_features"] == 3
    assert projection["n_library_released_features"] == 1
    assert projection["n_library_review_features"] == 1
    assert projection["n_unannotated_features"] == 1
    assert projection["partition_satisfied"] is True


def test_independent_ms1_inventory_adds_features_without_dda_and_proves_partition(
    tmp_path: Path,
):
    inventory_source = tmp_path / "feature_inventory.json"
    inventory_source.write_text('{"features":[]}', encoding="utf-8")
    payload = build_compound_conclusions(
        "anonymous-run",
        [_row("scan-1", "LF-1", 4.1, _match())],
        feature_inventory=[
            {
                "feature_id": "LF-1",
                "mz": 330.1696,
                "rt_minutes": 4.1,
                "intensity": 2000.0,
                "charge": 1,
                "ion_mode": "POS",
            },
            {
                "feature_id": "LF-2",
                "mz": 342.1700,
                "rt_minutes": 7.8,
                "intensity": 1500.0,
                "charge": 1,
                "ion_mode": "POS",
            },
        ],
        feature_inventory_source=inventory_source,
    )
    by_feature = {item["feature_id"]: item for item in payload["compounds"]}
    assert set(by_feature) == {"lcms:LF-1", "lcms:LF-2"}
    assert by_feature["lcms:LF-2"]["conclusion"]["key"] == "unannotated"
    assert by_feature["lcms:LF-2"]["unresolved_event_reason"] == (
        "no_associated_ms2"
    )
    assert by_feature["lcms:LF-2"]["spectra_count"] == 0
    assert by_feature["lcms:LF-2"]["representative"][
        "lcms_feature_inventory_only"
    ] is True

    projection = payload["metadata"]["feature_result_projection"]
    assert projection["independent_inventory_available"] is True
    assert projection["n_detected_ion_features"] == 2
    assert projection["n_visible_lcms_features"] == 2
    assert projection["result_lcms_feature_count_before_inventory_projection"] == 1
    assert projection["missing_inventory_feature_count"] == 0
    assert projection["unexpected_lcms_feature_count"] == 0
    assert projection["partition_satisfied"] is True
    assert payload["metadata"]["feature_inventory_source_sha256"] == (
        hashlib.sha256(inventory_source.read_bytes()).hexdigest()
    )


def test_independent_ms1_inventory_rejects_duplicate_feature_ids():
    with pytest.raises(ValueError, match="Duplicate feature inventory identifier"):
        build_compound_conclusions(
            "anonymous-run",
            [],
            feature_inventory=[
                {"feature_id": "LF-1", "mz": 100.0},
                {"feature_id": "LF-1", "mz": 101.0},
            ],
        )


def test_inventory_content_change_invalidates_artifact_even_with_same_stat(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    source = tmp_path / "scan-results.json"
    source.write_text("[]\n", encoding="utf-8")
    inventory = tmp_path / "feature-inventory.json"
    inventory.write_bytes(b"AAAAAAAA")
    inventory_stat = inventory.stat()
    artifact = tmp_path / "compounds.json"
    metadata = {
        "algorithm_id": compound_projection.COMPOUND_CONCLUSION_ALGORITHM_ID,
        "report_role_catalog": compound_projection.report_role_catalog_summary(),
        "source_artifact_size_bytes": source.stat().st_size,
        "source_artifact_mtime_ns": source.stat().st_mtime_ns,
        "feature_inventory_source": str(inventory.resolve()),
        "feature_inventory_source_size_bytes": inventory_stat.st_size,
        "feature_inventory_source_mtime_ns": inventory_stat.st_mtime_ns,
        "feature_inventory_source_sha256": hashlib.sha256(
            inventory.read_bytes()
        ).hexdigest(),
        "primary_identity_artifact": None,
    }
    artifact.write_text(
        json.dumps({"metadata": metadata}), encoding="utf-8"
    )
    monkeypatch.setattr(
        compound_projection, "load_v15_flagship_policy", lambda: {}
    )
    assert compound_projection.artifact_is_current(
        artifact, source, feature_inventory_artifact=inventory
    )

    inventory.write_bytes(b"BBBBBBBB")
    os.utime(
        inventory,
        ns=(inventory_stat.st_atime_ns, inventory_stat.st_mtime_ns),
    )
    assert inventory.stat().st_size == inventory_stat.st_size
    assert inventory.stat().st_mtime_ns == inventory_stat.st_mtime_ns
    assert not compound_projection.artifact_is_current(
        artifact, source, feature_inventory_artifact=inventory
    )


def test_r4_feature_identifier_and_unannotated_rows_reach_ui_and_exports():
    export_js = (ROOT / "web/static/js/export.js").read_text(encoding="utf-8")
    overview_js = (ROOT / "web/static/js/sample-overview.js").read_text(
        encoding="utf-8"
    )
    report_js = (ROOT / "web/static/js/report-generator.js").read_text(
        encoding="utf-8"
    )
    report_py = (ROOT / "web/report_common.py").read_text(encoding="utf-8")
    assert "'Feature ID', 'Result ID', '结果状态'" in export_js
    assert "matched_only: 'false'" in export_js
    assert "unannotated_ion_feature" in export_js
    assert "Feature ID ${esc(featureId)}" in overview_js
    assert "3.9 未注释离子特征" in report_js
    assert '"Record_Type": "unannotated_feature"' in report_py
    assert 'base["Feature_ID"] = compound.get("feature_id", "")' in report_py
