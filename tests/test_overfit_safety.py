"""Cross-module fail-closed defaults introduced by the overfit audit."""

from config import settings
from matching.classification_provenance import classification_source_trust
from scripts.mine_pathways import mine_compound


def _spectra():
    return [
        {
            "precursor_mz": 300.0,
            "peaks": [(153.0188, 100.0), (137.0239, 40.0), (121.029, 20.0)],
            "adduct": "[M+H]+",
            "charge": 1,
            "polarity": "pos",
            "spectrum_id": "source-a:1",
            "dataset_source": "source-a",
        },
        {
            "precursor_mz": 300.0,
            "peaks": [(153.0189, 100.0), (137.0240, 35.0), (121.0291, 18.0)],
            "adduct": "[M+H]+",
            "charge": 1,
            "polarity": "pos",
            "spectrum_id": "source-b:1",
            "dataset_source": "source-b",
        },
    ]


def test_unvalidated_production_influence_defaults_are_off():
    assert settings.USE_PATHWAY_ALIGN is False
    assert settings.USE_PATHWAY_CONFIDENCE is False
    assert settings.USE_DIAGNOSTIC_RESCUE is False
    assert settings.USE_MARGIN_UPGRADE is False
    assert settings.ISOPAIR_RANK_INFLUENCE_ENABLED is False


def test_existing_and_ml_class_sources_are_not_independent_evidence():
    assert classification_source_trust("existing")[0] is False
    assert classification_source_trust("ml")[0] is False
    assert classification_source_trust("ml_research_unvalidated")[0] is False
    assert classification_source_trust("smarts_structure_rule")[0] is True


def test_mined_pathway_records_candidate_conditioning_and_weak_label_source():
    pathway = mine_compound(
        1,
        _spectra(),
        classification="flavonoid",
        classification_source="ml",
    )

    provenance = pathway["provenance"]
    assert provenance["candidate_conditioned"] is True
    assert provenance["independent_evidence"] is False
    assert provenance["eligible_for_production_score"] is False
    assert provenance["classification_source_trusted"] is False
    assert provenance["dataset_sources"] == ["source-a", "source-b"]
