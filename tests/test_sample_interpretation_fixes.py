"""Regression coverage for the V10 sample-interpretation corrections."""

import sqlite3
import re
from pathlib import Path

from matching.compound_annotations import annotate_compound_dict
from matching.fragmentation_rules import normalize_classification
from matching.honest_calibration import annotate
from matching.screening_catalog import (
    audit_library_coverage,
    catalog_summary,
    load_screening_catalog,
    screening_annotation_for,
    validate_screening_catalog,
)
from mssystem.classification import RuleBasedClassifier


ROOT = Path(__file__).resolve().parents[1]
REPORT = (ROOT / "web/static/js/report-generator.js").read_text(encoding="utf-8")
OVERVIEW = (ROOT / "web/static/js/sample-overview.js").read_text(encoding="utf-8")
RESULTS = (ROOT / "web/static/js/results.js").read_text(encoding="utf-8")
INDEX = (ROOT / "web/static/index.html").read_text(encoding="utf-8")


def test_diarylheptanoid_structure_rule_is_narrow_and_normalized():
    classifier = RuleBasedClassifier()

    assert (
        classifier.predict("COc1cc(CCC(=O)CCCCc2ccccc2)ccc1O")
        == "diarylheptanoid"
    )
    assert normalize_classification("diarylheptanoid") == "diarylheptanoid"
    assert normalize_classification("二芳基庚烷") == "diarylheptanoid"
    assert classifier.predict("COc1cc(/C=C/C(=O)O)ccc1O") == "phenolic_acid"
    assert (
        classifier.predict("C[C@H]1OC(=O)c2c(O)ccc(O)c2[C@H]1O")
        == "coumarin"
    )
    assert classifier.predict("COC(=O)c1ccccc1O") == "Unknown"
    assert classifier.predict("COC(=O)/C=C/c1ccc(O)c(OC)c1") != "phenolic_acid"
    assert classifier.predict("CCCCOC(=O)c1ccccc1C(=O)OCCCC") == "Unknown"


def test_glycerophospholipid_rule_requires_core_acyl_and_carbon_context():
    classifier = RuleBasedClassifier()
    lpa_18_2 = "CCCCCC=CCC=CCCCCCCCC(=O)OCC(O)COP(=O)(O)O"

    assert classifier.predict(lpa_18_2) == "lipid"
    assert classifier.predict_detailed(lpa_18_2)["source"] == "smarts"
    assert classifier.predict("OCC(O)COP(=O)(O)O") == "Unknown"
    assert classifier.predict("NC1=NC=NC2=C1N=CN2C3OC(COP(=O)(O)O)C(O)C3O") != "lipid"
    assert classifier.predict("CCCCCCCCCCCCCCCC(=O)O") != "lipid"


def test_lower_ranked_class_cannot_be_projected_onto_unverified_top1():
    top1 = {
        "compound_name": "unverified candidate",
        "inchikey": "AAAAAAAAAAAAAA-BBBBBBBBBB-C",
        "formula": "C10H12O2",
        "smiles": "CCCCCCCCCC",
        "has_structure": True,
        "cosine": 0.82,
        "final_score": 0.52,
        "rank_score": 0.52,
        "reverse_match": 0.8,
        "mass_ppm": 2.0,
        "classification": "lipid",
        "classification_source": "ml",
        "source_tier": "experimental",
    }
    alternative = {
        **top1,
        "compound_name": "trusted glycoside alternative",
        "inchikey": "CCCCCCCCCCCCCC-DDDDDDDDDD-E",
        "formula": "C11H14O3",
        "classification": "glycoside",
        "classification_source": "curated",
        "cosine": 0.78,
        "final_score": 0.49,
        "rank_score": 0.49,
    }
    top1["top_candidates"] = [top1.copy(), alternative]

    result = annotate(top1)

    assert result["class_consensus"] == "other"
    assert result["class_agreement"] == 0.0
    assert (
        result["class_consensus_status"]
        == "top1_unverified_no_cross_candidate_transfer"
    )
    assert result["class_candidate_hypotheses"] == {"glycoside": 1}


def test_reliable_runtime_class_is_projected_without_rewriting_raw_library_value():
    match = {
        "classification": "Unknown",
        "expert_details": {
            "classification_provenance": {
                "reliable": True,
                "effective_classification": "phenolic_acid",
                "reliability_reason": "smiles_structure_rule",
                "input_confidence": 0.2,
            }
        },
    }

    annotate_compound_dict(match)

    assert match["classification"] == "Unknown"
    assert match["display_class"] == "phenolic_acid"
    assert match["classification_reliable"] is True
    assert match["classification_projection_preserved_raw"] is True
    assert match["display_classification_confidence"] == 1.0


def test_yakuchinone_identity_annotation_repairs_display_class_only():
    match = {
        "inchikey": "TXELARZTKDBEKS-UHFFFAOYSA-N",
        "compound_name": "Yakuchinone A",
        "classification": "Unknown",
    }

    annotate_compound_dict(match)

    assert match["classification"] == "Unknown"
    assert match["display_class"] == "diarylheptanoid"
    assert match["compound_role"] == "curated_natural"
    assert match["natural_product_candidate"] is True


def test_dihydroisocoumarin_annotation_overrides_stale_phenolic_acid_projection():
    match = {
        "inchikey": "UKOUEEGEHXWHHO-MOFOKWOHSA-N",
        "classification": "Unknown",
        "expert_details": {
            "classification_provenance": {
                "reliable": True,
                "effective_classification": "phenolic_acid",
                "reliability_reason": "smiles_structure_rule",
            }
        },
    }

    annotate_compound_dict(match)

    assert match["classification"] == "Unknown"
    assert match["display_class"] == "coumarin"
    assert match["display_classification_source"] == "curated_identity_annotation"


def test_screening_catalog_is_exact_identity_first_and_structure_validated():
    payload = load_screening_catalog()
    summary = catalog_summary()

    assert validate_screening_catalog(payload, validate_structures=True) == []
    assert summary["entries"] == 18
    assert summary["indexed_full_inchikeys"] == 20
    assert summary["regulatory_claims_encoded"] is False

    atrazine = screening_annotation_for("MXWJVTOOROXGIU-UHFFFAOYSA-N")
    assert atrazine["canonical_name"] == "Atrazine"
    assert atrazine["compound_role"] == "pesticide_screening_target"
    assert atrazine["regulatory_status"] == "not_encoded"
    assert screening_annotation_for("Atrazine") == {}
    assert screening_annotation_for("MXWJVTOOROXGIU") == {}

    thiamethoxam_alias = screening_annotation_for(
        "NWWZPOKUUAIXIW-UHFFFAOYSA-N"
    )
    assert thiamethoxam_alias["canonical_name"] == "Thiamethoxam"
    assert (
        thiamethoxam_alias["screening_identity_scope"]
        == "connectivity_alias_full_key"
    )


def test_screening_coverage_audit_uses_the_runtime_source_tier_gate(tmp_path):
    database = tmp_path / "catalog-audit.db"
    connection = sqlite3.connect(database)
    connection.executescript(
        """
        CREATE TABLE Compounds (
            id INTEGER PRIMARY KEY,
            compound_name TEXT,
            inchikey TEXT,
            structure_identity_verified INTEGER
        );
        CREATE TABLE Spectra (
            id INTEGER PRIMARY KEY,
            source_tier TEXT,
            ion_mode TEXT
        );
        CREATE TABLE Mapping (compound_id INTEGER, spectrum_id INTEGER);
        INSERT INTO Compounds VALUES (
            1, 'Atrazine', 'MXWJVTOOROXGIU-UHFFFAOYSA-N', 1
        );
        INSERT INTO Spectra VALUES
            (1, 'experimental', 'POS'),
            (2, 'curated', 'POS'),
            (3, 'massbank_supplement', 'NEG'),
            (4, 'merged', 'POS'),
            (5, 'legacy_unresolved', NULL),
            (6, 'reference_standard', 'POS'),
            (7, 'quarantined', 'POS');
        INSERT INTO Mapping VALUES
            (1, 1), (1, 2), (1, 3), (1, 4), (1, 5), (1, 6), (1, 7);
        """
    )
    connection.commit()
    connection.close()

    audit = audit_library_coverage(database)
    atrazine = next(
        item for item in audit["identities"]
        if item["inchikey"] == "MXWJVTOOROXGIU-UHFFFAOYSA-N"
    )

    assert atrazine["spectrum_count"] == 7
    assert atrazine["searchable_spectrum_count"] == 5


def test_integrated_report_uses_local_candidate_events_and_explains_question_marks():
    assert "matchedRowsFromGroups(groups)" in REPORT
    assert "const status = group?.status || SampleResultData.evidenceStatus(rows)" in REPORT
    assert "const clusterToleranceMinutes = 0.025" in REPORT
    assert "不会再用化合物的全样本唯一代表 RT 隐藏其他事件" in REPORT
    assert "附近未采 MS2" in REPORT
    assert "已采该峰离子但无候选过门槛" in REPORT
    assert "仅 Low/Unknown" in REPORT
    assert "rankScore: Number(match.flagship_entity_score" in REPORT


def test_expert_score_is_auditable_but_has_zero_flagship_identity_influence():
    assert "规则证据评估" in OVERVIEW
    assert "match.expert_score_source" in OVERVIEW
    assert "分类、生源、规则证据与碎裂解释" in REPORT
    assert "不直接参与身份排序" in REPORT
    assert "Legacy_rule_aggregate_uncalibrated_audit" in REPORT
    assert "Rule_aggregate_source" in REPORT
    assert "Legacy_rule_weighted_term_audit" in REPORT
    assert "confidence-guide.js?v=14.0.4" in INDEX
    assert re.search(r'<script\s+src="/static/js/sample-overview\.js\?v=[^"\s]+"', INDEX)
    assert re.search(r'<script\s+src="/static/js/report-generator\.js\?v=[^"\s]+"', INDEX)
    assert "实测 MS1 特征质心" in REPORT
    assert "ms2_feature_link_tolerance_ppm" in REPORT


def test_candidate_inferred_adduct_is_labeled_as_nonindependent_evidence():
    assert "候选离子假设" in RESULTS
    assert "非样品独立证据" in RESULTS
    assert "candidateIonFormHypothesis" in OVERVIEW
    assert "候选离子假设" in REPORT
    assert "非独立证据" in REPORT
