"""Regression coverage for the result-page fixes requested on 2026-08-12."""

from pathlib import Path

from matching.compound_annotations import annotate_compound_dict
from mssystem.classification import RuleBasedClassifier


ROOT = Path(__file__).resolve().parents[1]
OVERVIEW = (ROOT / "web/static/js/sample-overview.js").read_text(
    encoding="utf-8"
)
RESULTS = (ROOT / "web/static/js/results.js").read_text(encoding="utf-8")
STYLE = (ROOT / "web/static/css/style.css").read_text(encoding="utf-8")


def test_result_home_loads_and_renders_the_selected_targeted_database_layer():
    assert "/expected-compounds/groups?limit=1000" in OVERVIEW
    assert 'id="sampleTargetedResults"' in OVERVIEW
    assert "Results.expectedCompoundsCard(payload)" in OVERVIEW
    assert "expectedCompoundsCard: _expectedCompoundsCard" in RESULTS
    assert "靶向数据库结果（Expected Compounds" in RESULTS
    assert "selected_dataset_sources" in RESULTS
    assert "已应用数据集" in RESULTS


def test_blank_status_uses_a_compact_wrapping_label_and_full_hover_explanation():
    assert "未评估（无空白证据）" in OVERVIEW
    assert "未提供可用空白，背景尚未评估；" in OVERVIEW
    assert "已执行背景排除；本条保留" in OVERVIEW
    assert "任务共扣除" in OVERVIEW
    assert "当前这条结果未被排除" in OVERVIEW
    assert 'class="sample-key-score-note"' in OVERVIEW
    assert ".sample-key-scores > .sample-key-score-note" in STYLE
    assert "overflow-wrap: anywhere" in STYLE


def test_high_identity_without_a_class_is_described_as_two_independent_axes():
    expected = "尚无可靠类别归属（谱库候选证据为 High）"
    assert expected in OVERVIEW
    assert expected in RESULTS
    assert "分类可靠性与谱库身份匹配等级独立评估" in OVERVIEW
    assert "当前类别体系未覆盖" not in OVERVIEW
    assert "当前类别体系未覆盖" not in RESULTS


def test_organic_sulfonic_acid_structure_rule_repairs_cached_unknown_class():
    classifier = RuleBasedClassifier()
    decyl = "CCCCCCCCCCC1=CC=CC=C1S(=O)(=O)O"
    dodecyl = "CCCCCCCCCCCCC1=CC=C(C=C1)S(O)(=O)=O"

    assert classifier.predict(decyl) == "sulfonate"
    assert classifier.predict(dodecyl) == "sulfonate"
    assert classifier.predict("CS(=O)(=O)NC1=CC=CC=C1") != "sulfonate"

    match = {
        "classification": "Unknown",
        "display_class": "other",
        "classification_source": "other",
        "smiles": dodecyl,
        "structure_quality": "validated",
    }
    annotate_compound_dict(match)

    assert match["classification"] == "Unknown"
    assert match["display_class"] == "sulfonate"
    assert match["classification_reliable"] is True
    assert match["display_classification_source"].startswith(
        "runtime_structure_rule:"
    )


def test_lauramine_oxide_exact_identity_has_reliable_specific_class():
    match = {
        "inchikey": "SYELZBGXAIXKHU-UHFFFAOYSA-N",
        "compound_name": "N,N-Dimethyldodecylamine N-oxide",
        "classification": "other",
        "classification_source": "other",
        "smiles": "CCCCCCCCCCCC[N+](C)(C)[O-]",
        "structure_quality": "validated",
    }

    annotate_compound_dict(match)

    assert match["compound_name"] == "Dodecyldimethylamine oxide"
    assert match["display_class"] == "amine_oxide"
    assert match["classification_source"] == "curated_identity_annotation"
    assert match["classification_reliable"] is True
    assert match["classification_confidence"] == 1.0


def test_untrusted_ml_sugar_label_is_revalidated_from_open_chain_structure():
    match = {
        "compound_name": "L-Mannose",
        "classification": "sugar",
        "display_class": "sugar",
        "classification_source": "ml",
        "display_classification_source": "ml",
        "classification_reliable": False,
        "smiles": "OC[C@H](O)[C@@H](O)[C@H](O)[C@H](O)C=O",
        "structure_quality": "validated",
    }

    annotate_compound_dict(match)

    assert match["classification"] == "sugar"
    assert match["display_class"] == "sugar"
    assert match["classification_reliable"] is True
    assert match["display_classification_source"].startswith(
        "runtime_structure_rule:"
    )
