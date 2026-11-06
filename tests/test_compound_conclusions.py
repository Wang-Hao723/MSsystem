from __future__ import annotations

import json

from fastapi.testclient import TestClient

from pipeline.compound_conclusions import build_compound_conclusions
from matching.v15_identity_resolution import apply_v15_identity_resolution
from web import server


def _match(
    name: str,
    *,
    inchikey: str,
    formula: str = "C2H6O",
    confidence: str = "High",
    isomer: bool = False,
) -> dict:
    return {
        "compound_name": name,
        "inchikey": inchikey,
        "formula": formula,
        "smiles": "CCO",
        "has_structure": True,
        "dataset_source": "synthetic-test",
        "confidence": confidence,
        "isomer_group_confidence": "Medium" if isomer else "Unknown",
        "isomer_ambiguous": isomer,
        "identification_scope": "isomer_group" if isomer else "2d_connectivity",
        "rank_score": 0.91,
        "final_score": 0.90,
        "cosine": 0.93,
        "mass_ppm": 1.2,
        "classification": "flavonoid",
        "classification_source": "curated",
        "classification_reliable": True,
        "expert_score": 0.3,
        "expert_score_source": "rule_engine",
        "expert_details": {
            "reasoning": ["生源：莽草酸途径候选", "酚羟基规则提供正向支持"],
            "neutral_loss_analysis": {
                "best_pathway": [
                    {
                        "type": "neutral_loss",
                        "from_mz": 317.065,
                        "to_mz": 302.042,
                        "neutral_loss": "CH3",
                    }
                ],
                "neutral_losses": [{"formula": "CH3"}],
            },
        },
    }


def _row(
    spectrum_id: str,
    match: dict,
    *,
    rt: float,
    rt_unit: str = "minute",
    blank_category: str = "sample_enriched",
) -> dict:
    return {
        "spectrum_id": spectrum_id,
        "precursor_mz": 317.065,
        "retention_time": rt,
        "retention_time_unit": rt_unit,
        "tic": 1000.0,
        "blank_used": True,
        "blank_evidence_category": blank_category,
        "match": match,
    }


def _event_row(
    spectrum_id: str,
    match: dict,
    *,
    event_id: str,
    collision_energy: float,
    rt: float = 1.0,
) -> dict:
    row = _row(spectrum_id, match, rt=rt)
    row.update({
        "precursor_feature_id": event_id,
        "collision_energy": collision_energy,
        "collision_energy_unit": "eV",
    })
    return row


def _primary_artifact(
    task_id: str,
    candidates: list[dict],
    member_indices: list[int],
) -> dict:
    payload = {
        "metadata": {
            "task_id": task_id,
            "production_identity_selector": True,
            "algorithm_id": "synthetic-v15-primary",
            "policy_sha256": "synthetic-policy",
            "entity_semantic_sha256": "synthetic-entity-digest",
        },
        "policy": {
            "identity_resolution": {
                "enabled": True,
                "algorithm_id": (
                    "mssystem-v15-identity-scope-resolution-1.2.0"
                ),
                "competition_set_rule": (
                    "all_non_vetoed_v15_ranked_shortlist_candidates"
                ),
                "maximum_ranked_shortlist_size": 5,
            }
        },
        "entities": [
            {
                "entity_id": f"entity:{task_id}:1",
                "member_result_indices": member_indices,
                "member_spectrum_ids": [f"scan-{index}" for index in member_indices],
                "distinct_collision_energy_bins": 2,
                "winner": candidates[0],
                "leading_candidate_when_abstained": None,
                "candidates": candidates,
                "decision": {
                    "status": "selected",
                    "scope": "2d_connectivity_candidate",
                    "blockers": [],
                },
            }
        ],
    }
    return apply_v15_identity_resolution(payload)


def test_repeated_scans_collapse_to_one_compound_and_one_peak_region():
    match = _match(
        "Ethanol",
        inchikey="LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
    )
    payload = build_compound_conclusions(
        "dedup",
        [
            {
                **_row("scan-1", match, rt=60.0, rt_unit="second"),
                "precursor_feature_id": "PF000001",
            },
            {
                **_row("scan-2", match, rt=1.2),
                "precursor_feature_id": "PF000001",
            },
        ],
    )

    assert payload["metadata"]["unique_report_rows"] == 1
    assert payload["metadata"]["high_confidence_exact_structure_count"] == 1
    assert payload["metadata"]["high_confidence_isomer_group_count"] == 0
    assert payload["metadata"]["high_confidence_scope_group_count"] == 1
    compound = payload["compounds"][0]
    assert compound["conclusion"]["key"] == "high"
    assert compound["spectra_count"] == 2
    assert compound["peak_cluster_count"] == 1
    assert compound["occurrences"][0]["retention_time_unit"] == "minute"


def test_evidence_release_blocker_never_hides_or_reorders_saved_candidate():
    match = _match(
        "Sinomenine",
        inchikey="INYYVPJSBIVGPH-QHRIQVFBSA-N",
        confidence="Medium",
    )
    match.update({
        "cosine": 0.35,
        "rank_score": 0.71,
        "final_score": 0.39,
        "evidence_governance": {
            "policy_id": "evigate-ms-locked-v1-2026-08-02",
            "claim_status": "insufficient_for_release",
            "release_eligible_before_replication": False,
            "release_blockers": ["cosine_below_locked_release_threshold"],
            "ranking_changed": False,
            "confidence_field_overwritten": False,
        },
    })

    payload = build_compound_conclusions(
        "release-blocker-visible",
        [_row("scan-visible", match, rt=2.0)],
    )

    compound = payload["compounds"][0]
    assert compound["display"]["default_visible"] is True
    assert compound["display_name"] == "Sinomenine"
    assert compound["representative"]["match"]["rank_score"] == 0.71
    assert payload["metadata"]["main_result_count"] == 1


def test_language_alias_never_replaces_or_decorates_database_compound_name():
    match = _match(
        "Dodecyldimethylamine oxide",
        inchikey="WMSZOKIFGYFJHL-UHFFFAOYSA-N",
        formula="C14H31NO",
    )
    match["compound_name_cn"] = "十二烷基二甲基氧化胺"
    match["source_compound_name"] = "Dodecyldimethylamine oxide"

    payload = build_compound_conclusions(
        "language-boundary",
        [_row("scan-name", match, rt=9.778)],
    )

    compound = payload["compounds"][0]
    projected = compound["representative"]["match"]
    assert compound["display_name"] == "Dodecyldimethylamine oxide"
    assert projected["compound_name"] == "Dodecyldimethylamine oxide"
    assert projected["compound_name_cn"] == "十二烷基二甲基氧化胺"
    assert "十二烷基" not in compound["display_name"]


def test_library_pathway_alone_is_not_reported_as_sample_reconstruction():
    match = _match(
        "Reference-only pathway",
        inchikey="LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
    )
    match["expert_details"] = {}
    match["pathway_json"] = {
        "n_spectra": 4,
        "diagnostic_ions": [
            {"mz": 31.018, "freq": 0.75, "recurring": True}
        ],
        "neutral_losses": [],
    }
    payload = build_compound_conclusions(
        "reference-only",
        [_row("scan-reference", match, rt=1.0)],
    )
    compound = payload["compounds"][0]
    fragmentation = compound["auxiliary_interpretation"]["fragmentation"]
    flags = compound["representative"]["match"]["display_evidence_flags"]

    assert fragmentation["status"] == "reference_library_pathway_only"
    assert fragmentation["step_count"] == 0
    assert fragmentation["product_structure_count"] == 0
    assert fragmentation["reason"] == (
        "reference_pathway_available_but_no_sample_pathway_reconstructed"
    )
    assert flags["has_fragmentation"] is False
    assert flags["has_library_pathway_reference"] is True

    compact_flags = server._summary_evidence_flags(match)
    assert compact_flags["has_fragmentation"] is False
    assert compact_flags["has_library_pathway_reference"] is True
    assert compact_flags["fragmentation_step_count"] == 0


def test_candidate_inferred_adduct_is_preserved_but_not_independent_sample_evidence():
    match = _match(
        "Candidate-conditioned ion form",
        inchikey="LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
    )
    match["adduct"] = ""
    match["candidate_adduct_hypothesis"] = "[M-H]-"
    match["candidate_adduct_hypothesis_source"] = (
        "structure_formula_library_precursor"
    )
    row = _row("scan-adduct-hypothesis", match, rt=1.0)
    row["ion_mode"] = "NEG"

    payload = build_compound_conclusions(
        "candidate-adduct-hypothesis", [row]
    )
    adduct = payload["compounds"][0]["auxiliary_interpretation"][
        "adduct_assignment"
    ]

    assert adduct["sample_ion_form"] == "[M-H]-"
    assert adduct["sample_ion_form_source"] == (
        "candidate_structure_formula_precursor_hypothesis"
    )
    assert adduct["candidate_independent"] is False
    assert adduct["sample_ion_form_reliable"] is False
    assert adduct["candidate_ion_form_hypothesis"] == "[M-H]-"
    assert adduct["library_adduct"] is None


def test_old_task_roles_are_reprojected_and_external_surfactant_is_not_first():
    ordinary = _match(
        "Ordinary sample candidate",
        inchikey="AAAAAAAAAAAAAA-UHFFFAOYSA-N",
    )
    ordinary["rank_score"] = 0.80
    surfactant = _match(
        "Decylbenzenesulfonic acid",
        inchikey="UAZLASMTBCLJKO-UHFFFAOYSA-N",
        formula="C16H26O3S",
    )
    surfactant["smiles"] = "CCCCCCCCCCC1=CC=CC=C1S(=O)(=O)O"
    surfactant["rank_score"] = 0.99
    # Simulate the stale fields stored before the V12 exact-role supplement.
    surfactant["compound_role"] = "unclassified"
    surfactant["potential_contaminant"] = False

    payload = build_compound_conclusions(
        "role-refresh",
        [
            _row("surfactant", surfactant, rt=2.0),
            _row("ordinary", ordinary, rt=1.0),
        ],
    )

    assert payload["compounds"][0]["display_name"] == "Ordinary sample"
    refreshed = payload["compounds"][1]["representative"]["match"]
    assert refreshed["compound_role"] == "external_contaminant"
    assert refreshed["potential_contaminant"] is True
    assert refreshed["contaminant_category"] == (
        "alkylbenzene_sulfonate_surfactant"
    )


def test_same_formula_names_at_distinct_unlinked_events_stay_separate():
    methylquercetin = _match(
        "3'-O-Methylquercetin",
        inchikey="AAAAAAAAAAAAAA-UHFFFAOYSA-N",
        formula="C16H12O7",
        confidence="Medium",
        isomer=True,
    )
    nepetin = _match(
        "Nepetin",
        inchikey="BBBBBBBBBBBBBB-UHFFFAOYSA-N",
        formula="C16H12O7",
        confidence="Medium",
        isomer=True,
    )
    payload = build_compound_conclusions(
        "isomers",
        [_row("scan-1", methylquercetin, rt=14.05), _row("scan-2", nepetin, rt=14.20)],
    )

    assert payload["metadata"]["unique_report_rows"] == 2
    assert len({item["feature_id"] for item in payload["compounds"]}) == 2


def test_high_isomer_group_is_counted_at_high_reportable_scope_not_as_exact_structure():
    isomer = _match(
        "High isomer group",
        inchikey="ZZZZZZZZZZZZZZ-UHFFFAOYSA-N",
        formula="C16H12O7",
        confidence="Medium",
        isomer=True,
    )
    isomer["isomer_group_confidence"] = "High"

    payload = build_compound_conclusions(
        "high-isomer-group", [_row("scan-high-isomer", isomer, rt=7.5)]
    )

    conclusion = payload["compounds"][0]["conclusion"]
    assert conclusion["key"] == "isomer"
    assert conclusion["isomer_group_confidence"] == "High"
    assert conclusion["high_confidence_group_detected"] is True
    assert payload["metadata"]["high_confidence_detected_count"] == 0
    assert payload["metadata"]["high_confidence_exact_structure_count"] == 0
    assert payload["metadata"]["high_confidence_isomer_group_count"] == 1
    assert payload["metadata"]["high_confidence_scope_group_count"] == 1


def test_same_connectivity_variants_report_one_high_2d_connectivity_not_isomer_group():
    match = _match(
        "Sorbitol",
        inchikey="FBPFZTCFMRRESA-JGWLITMVSA-N",
        formula="C6H14O6",
        confidence="Medium",
        isomer=True,
    )
    match.update({
        "exact_structure_confidence": "Medium",
        "two_dimensional_connectivity_confidence": "High",
        "isomer_group_confidence": "High",
        "isomer_ambiguity_type": "same_connectivity_variant_group",
        "identification_scope": "exact_2d_connectivity",
    })

    payload = build_compound_conclusions(
        "same-connectivity", [_row("scan-sorbitol", match, rt=3.2)]
    )

    compound = payload["compounds"][0]
    assert compound["identity_scope"] == "2d_connectivity"
    assert compound["conclusion"]["key"] == "high"
    assert compound["conclusion"]["high_confidence_2d_connectivity_detected"] is True
    assert compound["conclusion"]["high_confidence_exact_structure_detected"] is False
    assert payload["metadata"]["high_confidence_2d_connectivity_count"] == 1
    assert payload["metadata"]["high_confidence_exact_structure_count"] == 0
    assert payload["metadata"]["high_confidence_isomer_group_count"] == 0
    assert payload["metadata"]["high_confidence_scope_group_count"] == 1


def test_connectivity_dedup_never_merges_incompatible_neutral_formulas():
    first = _match(
        "Entry A",
        inchikey="AAAAAAAAAAAAAA-UHFFFAOYSA-N",
        formula="C10H20O2",
    )
    second = _match(
        "Entry B",
        inchikey="AAAAAAAAAAAAAA-UHFFFAOYSA-N",
        formula="C12H24O2",
    )
    # Simulate a bad library row that reused the same structure string.
    second["smiles"] = first["smiles"]

    payload = build_compound_conclusions(
        "formula-guard",
        [_row("scan-a", first, rt=1.0), _row("scan-b", second, rt=2.0)],
    )

    assert payload["metadata"]["unique_report_rows"] == 2
    assert {item["formula"] for item in payload["compounds"]} == {
        "C10H20O2", "C12H24O2"
    }


def test_one_precursor_event_with_consistent_cross_energy_identity_counts_once():
    match = _match(
        "Consistent database name",
        inchikey="AAAAAAAAAAAAAA-UHFFFAOYSA-N",
        formula="C10H12O3",
    )
    payload = build_compound_conclusions(
        "event-consistent",
        [
            _event_row("ce20", match, event_id="PF000001", collision_energy=20),
            _event_row("ce40", match, event_id="PF000001", collision_energy=40),
        ],
    )

    assert payload["metadata"]["unique_report_rows"] == 1
    assert payload["metadata"]["main_result_count"] == 1
    assert payload["metadata"]["precursor_event_aggregation"] == {
        "policy": "one_report_row_per_ion_feature",
        "matched_event_count": 1,
        "multi_spectrum_event_count": 1,
        "consistent_identity_event_count": 1,
        "same_formula_unresolved_event_count": 0,
        "cross_formula_unresolved_event_count": 0,
        "incomplete_formula_unresolved_event_count": 0,
    }
    compound = payload["compounds"][0]
    assert compound["precursor_event_ids"] == ["PF000001"]
    assert compound["precursor_event_resolution_statuses"] == [
        "consistent_identity"
    ]


def test_one_precursor_event_with_same_formula_structures_becomes_one_isomer_group():
    first = _match(
        "Anisole",
        inchikey="AAAAAAAAAAAAAA-UHFFFAOYSA-N",
        formula="C7H8O",
    )
    second = _match(
        "Benzyl alcohol",
        inchikey="BBBBBBBBBBBBBB-UHFFFAOYSA-N",
        formula="C7H8O",
    )
    second["smiles"] = "CCOC"

    payload = build_compound_conclusions(
        "event-same-formula",
        [
            _event_row("ce20", first, event_id="PF000002", collision_energy=20),
            _event_row("ce40", second, event_id="PF000002", collision_energy=40),
        ],
    )

    assert payload["metadata"]["unique_report_rows"] == 1
    assert payload["metadata"]["main_result_count"] == 1
    assert payload["metadata"]["precursor_event_aggregation"][
        "same_formula_unresolved_event_count"
    ] == 1
    compound = payload["compounds"][0]
    assert compound["identity_scope"] == "isomer_group"
    assert compound["formula"] == "C7H8O"
    assert compound["conclusion"]["key"] == "isomer"
    assert compound["precursor_event_ids"] == ["PF000002"]
    assert compound["precursor_event_resolution_statuses"] == [
        "same_formula_isomer_group"
    ]
    assert {member["compound_name"] for member in compound["candidate_members"]} >= {
        "Anisole",
        "Benzyl alcohol",
    }


def test_one_precursor_event_with_different_formulas_is_not_a_compound_count():
    first = _match(
        "Acetophenone",
        inchikey="AAAAAAAAAAAAAA-UHFFFAOYSA-N",
        formula="C8H8O",
    )
    second = _match(
        "Ethyl benzoate",
        inchikey="BBBBBBBBBBBBBB-UHFFFAOYSA-N",
        formula="C9H10O2",
    )
    second["smiles"] = "CCOC"

    payload = build_compound_conclusions(
        "event-cross-formula",
        [
            _event_row("ce20", first, event_id="PF000003", collision_energy=20),
            _event_row("ce40", second, event_id="PF000003", collision_energy=40),
        ],
    )

    assert payload["metadata"]["unique_report_rows"] == 1
    assert payload["metadata"]["main_result_count"] == 0
    assert payload["metadata"]["unresolved_precursor_event_count"] == 1
    assert payload["metadata"]["precursor_event_aggregation"][
        "cross_formula_unresolved_event_count"
    ] == 1
    event = payload["compounds"][0]
    assert event["identity_scope"] == "unresolved_precursor_event"
    assert event["display_name"] == "Unresolved precursor event"
    assert event["formula"] == ""
    assert event["candidate_formulas"] == ["C8H8O", "C9H10O2"]
    assert event["conclusion"]["key"] == "unresolved_event"
    assert event["conclusion"]["counted_as_sample_compound"] is False
    assert event["precursor_event_ids"] == ["PF000003"]
    assert event["precursor_event_resolution_statuses"] == [
        "cross_formula_conflict"
    ]
    assert {member["compound_name"] for member in event["candidate_members"]} >= {
        "Acetophenone",
        "Ethyl benzoate",
    }


def test_one_precursor_event_with_reused_identity_but_different_formulas_is_unresolved():
    first = _match(
        "Bad metadata entry A",
        inchikey="AAAAAAAAAAAAAA-UHFFFAOYSA-N",
        formula="C10H20O2",
    )
    second = _match(
        "Bad metadata entry B",
        inchikey="AAAAAAAAAAAAAA-UHFFFAOYSA-N",
        formula="C12H24O2",
    )
    second["smiles"] = first["smiles"]

    payload = build_compound_conclusions(
        "event-cross-formula-reused-identity",
        [
            _event_row("ce20", first, event_id="PF000004", collision_energy=20),
            _event_row("ce40", second, event_id="PF000004", collision_energy=40),
        ],
    )

    assert payload["metadata"]["unique_report_rows"] == 1
    assert payload["metadata"]["main_result_count"] == 0
    assert payload["metadata"]["unresolved_precursor_event_count"] == 1
    event = payload["compounds"][0]
    assert event["identity_scope"] == "unresolved_precursor_event"
    assert event["candidate_formulas"] == ["C10H20O2", "C12H24O2"]
    assert event["conclusion"]["key"] == "unresolved_event"
    assert event["conclusion"]["counted_as_sample_compound"] is False


def test_one_precursor_event_with_structural_disagreement_and_missing_formula_is_separately_unresolved():
    first = _match(
        "Candidate with formula",
        inchikey="AAAAAAAAAAAAAA-UHFFFAOYSA-N",
        formula="C8H8O",
    )
    second = _match(
        "Candidate missing formula",
        inchikey="BBBBBBBBBBBBBB-UHFFFAOYSA-N",
        formula="",
    )
    second["smiles"] = "CCOC"

    payload = build_compound_conclusions(
        "event-incomplete-formula",
        [
            _event_row("ce20", first, event_id="PF000005", collision_energy=20),
            _event_row("ce40", second, event_id="PF000005", collision_energy=40),
        ],
    )

    aggregation = payload["metadata"]["precursor_event_aggregation"]
    assert aggregation["cross_formula_unresolved_event_count"] == 0
    assert aggregation["incomplete_formula_unresolved_event_count"] == 1
    event = payload["compounds"][0]
    assert event["identity_scope"] == "unresolved_precursor_event"
    assert event["candidate_formulas"] == ["C8H8O"]
    assert event["unresolved_event_reason"] == (
        "incomplete_candidate_formula_assignment"
    )
    assert event["precursor_event_resolution_statuses"] == [
        "formula_assignment_incomplete"
    ]
    assert event["conclusion"]["counted_as_sample_compound"] is False
    assert "至少一个候选缺少可用分子式" in event["conclusion"]["note"]
    missing = next(
        member for member in event["candidate_members"]
        if member["relationship"] == "candidate_missing_formula_assignment"
    )
    assert "formula" not in missing


def test_background_high_is_not_counted_and_auxiliary_outputs_are_retained():
    sample = _match(
        "Sample compound",
        inchikey="CCCCCCCCCCCCCC-UHFFFAOYSA-N",
    )
    background = _match(
        "Background compound",
        inchikey="DDDDDDDDDDDDDD-UHFFFAOYSA-N",
        formula="C2H7N",
    )
    background["smiles"] = "CCN"
    payload = build_compound_conclusions(
        "background",
        [
            _row("sample", sample, rt=2.0),
            _row("blank", background, rt=3.0, blank_category="blank_dominant"),
        ],
    )
    by_name = {item["display_name"]: item for item in payload["compounds"]}

    assert by_name["Sample compound"]["conclusion"]["key"] == "high"
    assert by_name["Background compound"]["conclusion"]["key"] == "background"
    auxiliary = by_name["Sample compound"]["auxiliary_interpretation"]
    assert auxiliary["classification"]["status"] == "available"
    assert auxiliary["biosource"]["status"] == "class_context"
    assert auxiliary["biosource"]["evidence_scope"] == (
        "chemical_class_context"
    )
    assert auxiliary["biosource"]["candidate_specific"] is False
    # A legacy aggregate and free-text reasoning cannot substitute for an
    # auditable per-dimension evidence record.
    assert auxiliary["expert"]["positive_support"] is False
    assert auxiliary["expert"]["status"] == "unavailable"
    assert auxiliary["expert"]["availability_reason"] == (
        "no_structured_dimension_assessments"
    )
    assert auxiliary["fragmentation"]["steps"][0]["neutral_loss"] == "CH3"
    assert all(
        section.get("identity_rank_influence") is False
        for section in (
            auxiliary["biosource"],
            auxiliary["expert"],
            auxiliary["fragmentation"],
        )
    )


def test_untyped_polarity_only_blank_warns_but_does_not_remove_sample_candidate():
    match = _match(
        "Potential sample compound",
        inchikey="EEEEEEEEEEEEEE-UHFFFAOYSA-N",
    )
    row = _row(
        "unverified-blank", match, rt=3.0,
        blank_category="blank_associated",
    )
    row["blank_evidence_role"] = "unspecified"

    payload = build_compound_conclusions("unverified-blank", [row])

    compound = payload["compounds"][0]
    assert compound["conclusion"]["key"] == "high"
    assert compound["conclusion"]["unverified_blank_association"] is True
    assert "不据此剔除" in compound["conclusion"]["blank_warning"]
    assert payload["metadata"]["background_related_count"] == 0
    assert payload["metadata"]["unverified_blank_warning_count"] == 1
    assert payload["metadata"]["unverified_blank_warning_high_count"] == 1


def test_lauramine_oxide_projection_separates_identity_adduct_and_fragment_hypotheses():
    match = {
        "compound_name": "N,N-Dimethyldodecylamine N-oxide",
        "inchikey": "SYELZBGXAIXKHU-UHFFFAOYSA-N",
        "formula": "C14H31NO",
        "smiles": "CCCCCCCCCCCC[N+](C)(C)[O-]",
        "structure_quality": "validated",
        "adduct": "[M+H]+",
        "dataset_source": "MassBank",
        "confidence": "High",
        "final_score": 0.7005,
        "rank_score": 0.7005,
        "cosine": 0.9564,
        "mass_ppm": 4.04,
        "classification": "Unknown",
        "expert_score": 0.225,
        "expert_details": {
            "dimension_assessments": {
                "neutral_loss": {
                    "status": "scored", "score": 0.2,
                    "evidence_count": 1,
                },
                "rule_coverage": {
                    "status": "not_applicable", "score": 0,
                    "evidence_count": 0,
                },
            },
            "neutral_loss_analysis": {
                "best_pathway": [{
                    "from_mz": 230.24877,
                    "to_mz": 213.24250,
                    "mass_difference": 17.00627,
                    "neutral_loss": "氨",
                    "neutral_loss_en": "Ammonia",
                    "formula": "NH3",
                }],
                "neutral_losses": [{
                    "mass_difference": 17.00627,
                    "neutral_loss": "氨",
                    "formula": "NH3",
                }],
            },
            "structure_fragment_analysis": {
                "plausible_fragments": [{
                    "observed_mz": 213.24250,
                    "fragment_smiles": "CCCCCCCCCCCC[N+](C)C",
                    "expected_mz": 213.24507,
                    "delta": 0.00257,
                    "adduct": "[M]+",
                    "assignment_type": "candidate_structure_hypothesis",
                }]
            },
        },
        "pathway_json": {
            "n_spectra": 10,
            "diagnostic_ions": [{
                "mz": 58.064, "freq": 0.8, "recurring": True,
            }],
            "neutral_losses": [{
                "loss_da": 17.00627,
                "annotation": "Ammonia",
                "recurring": True,
            }],
        },
    }
    row = _row("586678", match, rt=9.778)
    row.update({
        "precursor_mz": 230.24877,
        "ion_mode": "POS",
        "lcms_adduct": "[M+H]+",
        "lcms_neutral_mass": 229.24145,
    })

    payload = build_compound_conclusions("lauramine-oxide", [row])
    compound = payload["compounds"][0]
    projected = compound["representative"]["match"]
    auxiliary = compound["auxiliary_interpretation"]

    assert compound["display_name"] == "Dodecyldimethylamine oxide"
    assert projected["compound_name"] == "Dodecyldimethylamine oxide"
    assert projected["display_class"] == "amine_oxide"
    assert projected["contaminant_category"] == "amine_oxide_surfactant"
    assert auxiliary["biosource"]["status"] == (
        "not_applicable_external_background"
    )
    assert auxiliary["biosource"]["reason"] == (
        "natural_product_biosynthesis_not_applicable_to_external_background"
    )
    assert auxiliary["biosource"]["pathways"] == []
    assert auxiliary["biosource"]["knowledge_coverage"] == "not_applicable"
    adduct = auxiliary["adduct_assignment"]
    assert adduct["sample_ion_form"] == "[M+H]+"
    assert adduct["sample_ion_form_source"] == "lcms_ion_identity"
    assert adduct["candidate_independent"] is True
    assert adduct["library_adduct"] == "[M+H]+"
    assert auxiliary["expert"]["coverage"] == {
        "scored": 1,
        "no_evidence": 0,
        "not_applicable": 1,
        "unavailable": 0,
        "total": 2,
    }
    fragmentation = auxiliary["fragmentation"]
    assert fragmentation["status"] == "product_structure_candidates_available"
    assert fragmentation["product_structure_count"] == 1
    assert fragmentation["steps"] == []
    assert fragmentation["rejected_legacy_assignments"][0][
        "reason"
    ] == "accurate_mass_mismatch"
    assert fragmentation["candidate_specific_library_pathway"][
        "recurring_losses"
    ][0]["annotation_status"] == (
        "legacy_annotation_rejected_by_current_policy"
    )
    assert fragmentation["mechanism_confirmed"] is False


def test_only_formally_released_claim_is_called_confirmed(tmp_path):
    match = _match(
        "Released compound",
        inchikey="EEEEEEEEEEEEEE-UHFFFAOYSA-N",
    )
    claim_path = tmp_path / "claims.json"
    claim_path.write_text(
        json.dumps({
            "decisions": [
                {"spectrum_id": "released", "decision": {"released": True}}
            ]
        }),
        encoding="utf-8",
    )

    payload = build_compound_conclusions(
        "released",
        [_row("released", match, rt=1.0)],
        claim_decisions_path=claim_path,
    )

    assert payload["metadata"]["confirmed_present_count"] == 1
    assert payload["compounds"][0]["conclusion"]["key"] == "confirmed"


def test_primary_v15_identity_is_atomic_and_candidate_order_is_inherited():
    legacy = _match(
        "Legacy scan title",
        inchikey="LLLLLLLLLLLLLL-UHFFFAOYSA-N",
        formula="C9H10O3",
    )
    rows = [
        _event_row(
            f"scan-{index}",
            legacy,
            event_id="PF-PRIMARY-1",
            collision_energy=ce,
        )
        for index, ce in enumerate((20.0, 40.0))
    ]
    candidates = [
        {
            "identity_key": "std2d:AAAAAAAAAAAAAA",
            "entity_rank": 1,
            "entity_score": 0.81,
            "collision_energy_bins_supported": 2,
            "collision_energy_bins_total": 2,
            "vetoed": False,
            "candidate": {
                "compound_name": "Database Lead Name",
                "formula": "C7H8O",
                "inchikey": "AAAAAAAAAAAAAA-UHFFFAOYSA-N",
                "smiles": "COC1=CC=CC=C1",
                "dataset_source": "database-a",
            },
        },
        {
            "identity_key": "std2d:BBBBBBBBBBBBBB",
            "entity_rank": 2,
            "entity_score": 0.66,
            "collision_energy_bins_supported": 2,
            "collision_energy_bins_total": 2,
            "vetoed": False,
            "candidate": {
                "compound_name": "Database Alternative Name",
                "formula": "C7H8O",
                "inchikey": "BBBBBBBBBBBBBB-UHFFFAOYSA-N",
                "smiles": "OCC1=CC=CC=C1",
                "dataset_source": "database-b",
            },
        },
    ]
    primary = _primary_artifact("primary-atomic", candidates, [0, 1])

    payload = build_compound_conclusions(
        "primary-atomic",
        rows,
        primary_identity_artifact=primary,
    )

    compound = payload["compounds"][0]
    representative = compound["representative"]["match"]
    assert compound["identity_scope"] == "isomer_group"
    assert compound["display_name"] == "Database Lead Name"
    assert [member["compound_name"] for member in compound["candidate_members"]] == [
        "Database Lead Name",
        "Database Alternative Name",
    ]
    assert [member["flagship_entity_rank"] for member in compound["candidate_members"]] == [
        1,
        2,
    ]
    assert {
        field: representative[field]
        for field in ("compound_name", "formula", "inchikey", "smiles")
    } == {
        "compound_name": "Database Lead Name",
        "formula": "C7H8O",
        "inchikey": "AAAAAAAAAAAAAA-UHFFFAOYSA-N",
        "smiles": "COC1=CC=CC=C1",
    }
    assert "Legacy scan title" not in compound["candidate_names"]
    assert payload["metadata"]["atomic_candidate_identity_invariant"][
        "satisfied"
    ] is True
    assert payload["metadata"]["primary_identity_projection"][
        "candidate_order_source"
    ] == "v15_entity_rank"


def test_primary_database_formula_is_preserved_verbatim_for_charged_record():
    legacy = _match(
        "Legacy formula",
        inchikey="LLLLLLLLLLLLLL-UHFFFAOYSA-N",
        formula="C23H42N",
    )
    rows = [
        _event_row(
            f"scan-{index}",
            legacy,
            event_id="PF-PRIMARY-CHARGED",
            collision_energy=ce,
        )
        for index, ce in enumerate((20.0, 40.0))
    ]
    candidates = [
        {
            "identity_key": "std2d:WNBGYVXHFTYOBY",
            "entity_rank": 1,
            "entity_score": 0.91,
            "collision_energy_bins_supported": 2,
            "collision_energy_bins_total": 2,
            "vetoed": False,
            "candidate": {
                "compound_name": "Benzyltetradecyldimethylammonium",
                "formula": "C23H42N+",
                "inchikey": "WNBGYVXHFTYOBY-UHFFFAOYSA-N",
                "smiles": "CCCC[N+](C)(C)CC1=CC=CC=C1",
                "dataset_source": "database-charged-record",
            },
        }
    ]
    primary = _primary_artifact("primary-charged-formula", candidates, [0, 1])

    payload = build_compound_conclusions(
        "primary-charged-formula",
        rows,
        primary_identity_artifact=primary,
    )

    compound = payload["compounds"][0]
    assert compound["lead_candidate"]["formula"] == "C23H42N+"
    assert compound["representative"]["match"]["formula"] == "C23H42N+"
    assert payload["metadata"]["atomic_candidate_identity_invariant"][
        "satisfied"
    ] is True


def test_merged_connectivity_keeps_representative_candidate_fields_atomic():
    task_id = "primary-representative-atomic"
    legacy = _match(
        "Legacy connectivity title",
        inchikey="AZAZKLKDEOMJBJ-LEGACYTEST-N",
        formula="C22H33NO3",
    )
    rows = [
        _event_row(
            "scan-low-release",
            legacy,
            event_id="PF-ATOMIC-1",
            collision_energy=20.0,
        ),
        _event_row(
            "scan-high-release",
            legacy,
            event_id="PF-ATOMIC-2",
            collision_energy=40.0,
        ),
    ]
    higher_score_candidate_only = {
        "identity_key": "std2d:AZAZKLKDEOMJBJ",
        "entity_rank": 1,
        "entity_score": 0.90,
        "collision_energy_bins_supported": 1,
        "collision_energy_bins_total": 1,
        "vetoed": False,
        "candidate": {
            "compound_name": "Database connectivity lead",
            "formula": "C22H33NO3",
            "inchikey": "AZAZKLKDEOMJBJ-UHFFFAOYSA-N",
            "smiles": "OC1C(=C)C2CC13C4CC5C6(C)CN(CC)C4C5(C(O)CC6)C3CC2O",
        },
    }
    lower_score_released = {
        **higher_score_candidate_only,
        "entity_score": 0.75,
        "collision_energy_bins_supported": 2,
        "collision_energy_bins_total": 2,
        "candidate": {
            **higher_score_candidate_only["candidate"],
            "smiles": "C=C1C2CC3(C1O)C1CC4C5(C)CCC(O)C4(C1N(CC)C5)C3CC2O",
        },
    }
    primary = _primary_artifact(
        task_id,
        [higher_score_candidate_only],
        [0],
    )
    first = primary["entities"][0]
    first["decision"]["status"] = "selected_candidate_only"
    first["decision"]["scope"] = "single_collision_energy_candidate"
    first["identity_resolution"]["upstream_decision_status"] = (
        "selected_candidate_only"
    )
    first["identity_resolution"]["upstream_decision_scope"] = (
        "single_collision_energy_candidate"
    )
    second_payload = _primary_artifact(
        task_id,
        [lower_score_released],
        [1],
    )
    second = second_payload["entities"][0]
    second["entity_id"] = f"entity:{task_id}:2"
    primary["entities"].append(second)

    payload = build_compound_conclusions(
        task_id,
        rows,
        primary_identity_artifact=primary,
    )

    compound = payload["compounds"][0]
    representative = compound["representative"]["match"]
    assert representative["confidence"] == "High"
    assert compound["lead_candidate"]["smiles"] == representative["smiles"]
    assert compound["lead_candidate"]["smiles"] == (
        lower_score_released["candidate"]["smiles"]
    )
    assert payload["metadata"]["atomic_candidate_identity_invariant"][
        "satisfied"
    ] is True


def test_primary_cross_formula_competition_is_excluded_from_compound_count():
    legacy = _match(
        "Legacy selected title",
        inchikey="LLLLLLLLLLLLLL-UHFFFAOYSA-N",
        formula="C12H18O4",
    )
    rows = [
        _event_row(
            f"scan-{index}",
            legacy,
            event_id="PF-PRIMARY-2",
            collision_energy=ce,
        )
        for index, ce in enumerate((20.0, 40.0))
    ]
    candidates = [
        {
            "identity_key": "std2d:AAAAAAAAAAAAAA",
            "entity_rank": 1,
            "entity_score": 0.80,
            "collision_energy_bins_supported": 2,
            "collision_energy_bins_total": 2,
            "vetoed": False,
            "candidate": {
                "compound_name": "Formula-A ranking lead",
                "formula": "C12H18O4",
                "inchikey": "AAAAAAAAAAAAAA-UHFFFAOYSA-N",
                "smiles": "CCO",
            },
        },
        {
            "identity_key": "std2d:BBBBBBBBBBBBBB",
            "entity_rank": 2,
            "entity_score": 0.55,
            "collision_energy_bins_supported": 2,
            "collision_energy_bins_total": 2,
            "vetoed": False,
            "candidate": {
                "compound_name": "Formula-B competitor",
                "formula": "C12H16O3",
                "inchikey": "BBBBBBBBBBBBBB-UHFFFAOYSA-N",
                "smiles": "CCC",
            },
        },
    ]
    primary = _primary_artifact("primary-formula-conflict", candidates, [0, 1])

    payload = build_compound_conclusions(
        "primary-formula-conflict",
        rows,
        primary_identity_artifact=primary,
    )

    event = payload["compounds"][0]
    assert payload["metadata"]["main_result_count"] == 0
    assert event["conclusion"]["key"] == "unresolved_event"
    assert event["candidate_formulas"] == ["C12H16O3", "C12H18O4"]
    assert [member["compound_name"] for member in event["candidate_members"]] == [
        "Formula-A ranking lead",
        "Formula-B competitor",
    ]


def test_compound_api_materializes_current_sidecar_and_pages(tmp_path, monkeypatch):
    task_id = "compound-api"
    match = _match(
        "API compound",
        inchikey="FFFFFFFFFFFFFF-UHFFFAOYSA-N",
    )
    source = tmp_path / f"{task_id}_results.json"
    source.write_text(json.dumps([_row("scan-1", match, rt=1.0)]), encoding="utf-8")
    primary = _primary_artifact(
        task_id,
        [
            {
                "identity_key": "std2d:FFFFFFFFFFFFFF",
                "entity_rank": 1,
                "entity_score": 0.91,
                "collision_energy_bins_supported": 2,
                "collision_energy_bins_total": 2,
                "vetoed": False,
                "candidate": {
                    "compound_name": "API compound",
                    "formula": "C2H6O",
                    "inchikey": "FFFFFFFFFFFFFF-UHFFFAOYSA-N",
                    "smiles": "CCO",
                    "dataset_source": "synthetic-test",
                },
            }
        ],
        [0],
    )
    from matching.v15_flagship_ranker import load_v15_flagship_policy
    from pipeline.active_artifact_binding import active_implementation_sha256
    from pipeline.derived_artifact_paths import binding_metadata
    active_policy = load_v15_flagship_policy()
    primary['metadata'].update(binding_metadata(active_policy,
        actual_implementation_sha256=active_implementation_sha256(active_policy)))
    server._flagship_artifact_path(tmp_path, task_id).write_text(
        json.dumps(primary), encoding="utf-8"
    )
    monkeypatch.setattr(server, "TASK_DIR", tmp_path)
    monkeypatch.setattr(
        server,
        "_tasks",
        {task_id: {"id": task_id, "status": "completed", "filename": "sample.mzML"}},
    )
    monkeypatch.setattr(server, "_save_tasks", lambda: None)
    monkeypatch.setattr(server, "_flagship_artifact_is_current", lambda _path: True)

    response = TestClient(server.app).get(
        f"/api/tasks/{task_id}/compounds?offset=0&limit=1"
    )

    assert response.status_code == 200
    payload = response.json()
    assert payload["result_unit"] == (
        "deduplicated_compound_unresolved_isomer_group_or_"
        "unresolved_precursor_event"
    )
    assert payload["confirmed_is_level1_only"] is True
    assert payload["auxiliary_rank_influence"] is False
    assert payload["compounds_total"] == 1
    assert payload["compounds"][0]["display_name"] == "API compound"
    assert server._compound_conclusion_artifact_path(tmp_path, task_id).is_file()
    assert server._tasks[task_id]["compound_conclusion_summary"]["status"] == "available"
