import json
from pathlib import Path

from paper_v10.software_validation import (
    audit_saved_task,
    build_software_article_validation,
    render_markdown,
)


ROOT = Path(__file__).resolve().parents[1]


def test_retrospective_validation_keeps_claim_boundary_and_primary_unit():
    result = build_software_article_validation(
        ROOT,
        bootstrap_iterations=1000,
    )

    assert result["retrospective_evigate"]["confirmatory"] is False
    full = result["retrospective_evigate"]["arms"]["full_governed"]
    assert full["reported_identities"] == 37
    assert full["scope_compatible_correct"] == 18
    assert full["incompatible_or_incorrect"] == 19

    expert = result["paired_module_ablation"]["comparisons"]["expert_increment"]
    assert expert["n_paired_standards"] == 33
    assert expert["baseline_correct"] == 19
    assert expert["variant_correct"] == 17
    assert expert["positive_gain_supported"] is False
    assert result["claim_decisions"]["expert_rules_improve_top1"]["status"] == "unsupported"
    assert result["claim_decisions"]["journal_confirmatory_performance"]["status"] == "pending_independent_v10_lockbox"

    markdown = render_markdown(result)
    assert "development-informed" in markdown
    assert "one standard, one call" in markdown
    assert "Any-scan results are exploratory only" in markdown


def test_explicit_saved_task_audit_streams_only_named_artifacts(tmp_path):
    task_id = "tiny"
    task_dir = tmp_path / "output" / "tasks"
    task_dir.mkdir(parents=True)
    rows = [
        {"spectrum_id": "unmatched", "match": None},
        {
            "spectrum_id": "a-low",
            "match": {
                "inchikey": "AAAAAAAAAAAAAA-BBBBBBBBBB-C",
                "rank_score": 0.2,
                "classification": "other",
                "display_class": "flavonoid",
                "display_classification_source": "msclass_spectrum",
                "expert_score": 0.2,
                "expert_details": {
                    "neutral_loss_analysis": {"best_pathway": [{"loss": 18.01}]},
                    "score_fusion": {
                        "selected_source": "base_8d",
                        "components": {"v3_kb": {"score": 0.0}},
                    },
                },
            },
        },
        {
            "spectrum_id": "a-high",
            "match": {
                "inchikey": "AAAAAAAAAAAAAA-DDDDDDDDDD-E",
                "rank_score": 0.8,
                "classification": "flavonoid",
                "classification_source": "curated_identity_annotation",
                "expert_score": 0.4,
                "pathway_json": {"neutral_losses": [{"loss_da": 18.01}]},
                "expert_details": {
                    "classification_provenance": {
                        "analysis_class": "flavonoid",
                        "reliable": True,
                    },
                    "score_fusion": {
                        "selected_source": "v3_kb",
                        "components": {"v3_kb": {"score": 0.6}},
                    },
                },
            },
        },
        {
            "spectrum_id": "b",
            "match": {
                "inchikey": "CCCCCCCCCCCCCC-DDDDDDDDDD-E",
                "rank_score": 0.5,
                "classification": "other",
                "classification_source": "existing",
                "expert_score": 0.0,
                "expert_details": {"score_fusion": {"selected_source": "base_8d"}},
            },
        },
    ]
    (task_dir / f"{task_id}_results.json").write_text(
        json.dumps(rows), encoding="utf-8"
    )
    for suffix in (
        "claim_decisions.json",
        "evograph_fragmentation.json",
        "identification_records.json",
    ):
        (task_dir / f"{task_id}_{suffix}").write_text(
            json.dumps({"metadata": {"record_count": 2}}), encoding="utf-8"
        )

    audit = audit_saved_task(tmp_path, task_id)

    assert audit["spectra"] == 4
    assert audit["matched_spectra"] == 3
    assert audit["representative_candidate_identities"] == 3
    assert audit["representative_2d_connectivities"] == 2
    assert audit["classification"]["counts"] == {
        "tentative": 1,
        "reliable": 1,
        "unresolved": 1,
    }
    assert audit["fragmentation_evidence"]["library_derived_signature"] == 1
    assert audit["fragmentation_evidence"]["candidate_specific_pathway_eligible_for_fusion"] == 0
    assert audit["expert_evidence"]["positive_v3_knowledge_base_score"] == 1
