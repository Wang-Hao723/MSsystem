from pathlib import Path

from web.report_common import (
    classification_evidence_status,
    fragmentation_evidence_status,
    row_to_dict,
)


ROOT = Path(__file__).resolve().parents[1]


def test_model_display_class_is_not_promoted_by_presence_or_legacy_flag():
    status = classification_evidence_status(
        {
            "classification": "other",
            "classification_source": "other",
            "display_class": "flavonoid",
            "display_classification_source": "msclass_spectrum",
            "classification_reliable": True,
        }
    )

    assert status["label"] == "flavonoid"
    assert status["status"] == "tentative"
    assert status["reliable"] is False


def test_curated_and_runtime_structure_classes_are_reliable_but_existing_is_not():
    curated = classification_evidence_status(
        {
            "classification": "alkaloid",
            "classification_source": "curated_identity_annotation",
        }
    )
    runtime = classification_evidence_status(
        {
            "display_class": "terpenoid",
            "display_classification_source": "runtime_structure_rule:smiles",
            "classification_reliable": True,
        }
    )
    legacy = classification_evidence_status(
        {"classification": "flavonoid", "classification_source": "existing"}
    )

    assert curated["reliable"] is True
    assert runtime["reliable"] is True
    assert legacy["status"] == "tentative"


def test_fragmentation_export_names_evidence_and_never_confirms_mechanism():
    match = {
        "pathway_json": {"neutral_losses": [{"loss_da": 18.0106}]},
        "expert_details": {
            "neutral_loss_analysis": {"best_pathway": [{"neutral_loss": "H2O"}]},
            "pathway_align": {
                "eligible_for_fusion": False,
                "ranking_changed": False,
                "score_contribution": 0.0,
            },
        },
    }
    status = fragmentation_evidence_status(match)
    exported = row_to_dict({"spectrum_id": "scan=1", "match": match})

    assert status["evidence_type"] == (
        "library_derived_fragmentation_signature;empirical_neutral_loss_sequence"
    )
    assert status["mechanism_confirmed"] is False
    assert exported["Candidate_Pathway_Eligible_For_Fusion"] is False
    assert exported["Candidate_Pathway_Runtime_Influence"] is False


def test_frontend_source_contains_no_presence_equals_reliable_shortcut():
    results_source = (ROOT / "web" / "static" / "js" / "results.js").read_text(
        encoding="utf-8"
    )
    sample_source = (
        ROOT / "web" / "static" / "js" / "sample-overview.js"
    ).read_text(encoding="utf-8")

    assert "return { label: display, reliable: true" not in results_source
    assert "source === 'curated' || source === 'existing'" not in results_source
    assert "source === 'curated' || source === 'existing'" not in sample_source
    assert "Empirical Fragment Evidence" in results_source
