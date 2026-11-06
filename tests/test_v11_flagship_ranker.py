from __future__ import annotations

import copy
import json

import pytest

from matching.flagship_ranker import build_flagship_result
from matching.library_quality import formula_exact_mass
from pipeline.flagship import read_flagship_artifact, write_flagship_artifact
from utils.adduct import ion_mz


IK_A = "AAAAAAAAAAAAAA-BBBBBBBBBB-C"
IK_A_STEREO = "AAAAAAAAAAAAAA-CCCCCCCCCC-D"
IK_B = "DDDDDDDDDDDDDD-EEEEEEEEEE-F"


def _candidate(
    name: str,
    inchikey: str,
    strength: float,
    *,
    formula: str = "",
    adduct: str = "",
    classification: str = "class-a",
) -> dict:
    return {
        "compound_name": name,
        "inchikey": inchikey,
        "smiles": "",
        "formula": formula,
        "adduct": adduct,
        "cosine": strength,
        "entropy": strength,
        "reverse_match": strength,
        "fragment_score": strength,
        "intensity_score": strength,
        "mass_ppm": 0.0,
        "final_score": strength,
        "rank_score": strength,
        "expert_score": 0.99,
        "classification": classification,
        "classification_source": "synthetic-test",
        "dataset_source": "synthetic-test",
        "has_structure": bool(formula),
    }


def _row(
    ce: float,
    candidates: list[dict],
    *,
    feature: str = "PF-1",
    precursor_mz: float = 300.0,
    spectrum_id: str | None = None,
    raw_peaks: list[list[float]] | None = None,
    ion_mode: str = "POS",
    query_adduct: str = "",
    query_adduct_reliable: bool = False,
    lcms_adduct: str = "",
    lcms_neutral_mass: float | None = None,
    lcms_adduct_group_id: str = "",
    ion_identity: str = "",
    ion_identity_status: str = "",
    ion_identity_network_id: str = "",
    ion_identity_neutral_mass: float | None = None,
    ion_identity_evidence: dict | None = None,
    search_precursor_mz: float | None = None,
    search_precursor_mz_source: str = "",
) -> dict:
    match = copy.deepcopy(candidates[0])
    match["top_candidates"] = copy.deepcopy(candidates)
    return {
        "spectrum_id": spectrum_id or f"scan-{feature}-{ce}",
        "precursor_feature_id": feature,
        "lcms_feature_id": None,
        "precursor_mz": precursor_mz,
        "search_precursor_mz": search_precursor_mz,
        "search_precursor_mz_source": search_precursor_mz_source,
        "lcms_feature_used_for_search": search_precursor_mz is not None,
        "precursor_charge": 1,
        "retention_time": 1.0,
        "retention_time_unit": "minutes",
        "collision_energy": ce,
        "collision_energy_unit": "eV",
        "ion_mode": ion_mode,
        "query_adduct": query_adduct,
        "query_adduct_reliable": query_adduct_reliable,
        "lcms_adduct": lcms_adduct,
        "lcms_neutral_mass": lcms_neutral_mass,
        "lcms_adduct_group_id": lcms_adduct_group_id,
        "ion_identity": ion_identity,
        "ion_identity_status": ion_identity_status,
        "ion_identity_network_id": ion_identity_network_id,
        "ion_identity_neutral_mass": ion_identity_neutral_mass,
        "ion_identity_evidence": ion_identity_evidence or {},
        "raw_peaks": raw_peaks or [[100.0, 100.0], [150.0, 50.0]],
        "query_peaks": raw_peaks or [[100.0, 100.0], [150.0, 50.0]],
        "match": match,
    }


def _winner(payload: dict) -> dict:
    assert len(payload["entities"]) == 1
    winner = payload["entities"][0]["winner"]
    assert winner is not None
    return winner


def test_distinct_collision_energies_overturn_one_scan_winner():
    candidate_a = _candidate("stable-a", IK_A, 0.78)
    candidate_b_hot = _candidate("one-shot-b", IK_B, 0.98)
    candidate_b_weak = _candidate("one-shot-b", IK_B, 0.20)
    rows = [
        _row(10, [candidate_b_hot, candidate_a]),
        _row(20, [candidate_a, candidate_b_weak]),
        _row(30, [candidate_a, candidate_b_weak]),
    ]

    result = build_flagship_result("multi-ce", rows)
    entity = result["entities"][0]

    assert _winner(result)["identity_key"] == "ik14:AAAAAAAAAAAAAA"
    assert entity["distinct_collision_energy_bins"] == 3
    assert entity["ranking_changed_from_scan_top1"] is True
    assert entity["decision"]["scope"] == "2d_connectivity_candidate"


def test_missing_from_truncated_topk_is_not_imputed_as_zero_support():
    supported = _candidate("one-energy-evidence", IK_A, 0.82)
    other = _candidate("other", IK_B, 0.82)
    rows = [
        _row(10, [supported]),
        _row(20, [other]),
        _row(30, [other]),
    ]

    result = build_flagship_result("topk-missingness", rows)
    candidate = next(
        item for item in result["entities"][0]["candidates"]
        if item["identity_key"] == "ik14:AAAAAAAAAAAAAA"
    )

    assert candidate["score_components"]["mean_collision_energy_support"] > 0.5
    assert candidate["score_components"]["topk_absence_imputed_as_zero"] is False
    assert result["invariants"]["missing_from_truncated_topk_is_negative_evidence"] is False


def test_same_energy_technical_repeats_are_one_condition_record_not_independent_samples():
    candidate_a = _candidate("stable-a", IK_A, 0.82)
    candidate_b_hot = _candidate("repeat-heavy-b", IK_B, 0.99)
    candidate_b_weak = _candidate("repeat-heavy-b", IK_B, 0.18)
    rows = [
        *[
            _row(
                10,
                [candidate_b_hot, candidate_a],
                spectrum_id=f"duplicate-{index}",
            )
            for index in range(8)
        ],
        _row(20, [candidate_a, candidate_b_weak]),
        _row(30, [candidate_a, candidate_b_weak]),
    ]

    result = build_flagship_result("repeat-collapse", rows)
    entity = result["entities"][0]

    assert _winner(result)["identity_key"] == "ik14:AAAAAAAAAAAAAA"
    ten_ev = next(item for item in entity["collision_energy_bins"] if item["value"] == 10.0)
    assert ten_ev["technical_replicates"] == 8
    assert ten_ev["condition_evidence_count"] == 1
    assert ten_ev["independent_sample_vote_count"] == 0
    assert result["invariants"]["distinct_collision_energies_are_independent_samples"] is False


def test_expert_class_and_pathway_fields_cannot_change_ranking():
    rows = [
        _row(10, [_candidate("a", IK_A, 0.81), _candidate("b", IK_B, 0.75)]),
        _row(20, [_candidate("a", IK_A, 0.80), _candidate("b", IK_B, 0.74)]),
    ]
    perturbed = copy.deepcopy(rows)
    for row in perturbed:
        row["match"]["expert_score"] = 0.0
        row["match"]["classification"] = "opposite-class"
        row["match"]["pathway_json"] = {"invented": True}
        for index, candidate in enumerate(row["match"]["top_candidates"]):
            candidate["expert_score"] = 1.0 if index else 0.0
            candidate["classification"] = f"changed-{index}"
            candidate["pathway_json"] = {"arbitrary": index}

    baseline = build_flagship_result("forbidden-inputs", rows)["entities"][0]
    changed = build_flagship_result("forbidden-inputs", perturbed)["entities"][0]

    assert baseline["winner"]["identity_key"] == changed["winner"]["identity_key"]
    assert baseline["winner"]["entity_score"] == changed["winner"]["entity_score"]
    assert baseline["decision"] == changed["decision"]


def test_explicit_formula_adduct_contradiction_is_vetoed():
    formula = "C8H10N4O2"
    precursor = ion_mz(formula_exact_mass(formula), "[M+H]+")
    assert precursor is not None
    impossible = _candidate(
        "impossible-high",
        IK_B,
        0.98,
        formula="C20H40",
        adduct="[M+H]+",
    )
    viable = _candidate(
        "viable-lower",
        IK_A,
        0.78,
        formula=formula,
        adduct="[M+H]+",
    )

    result = build_flagship_result(
        "chemistry-veto",
        [_row(
            20,
            [impossible, viable],
            precursor_mz=precursor,
            query_adduct="[M+H]+",
            query_adduct_reliable=True,
        )],
    )
    entity = result["entities"][0]

    assert _winner(result)["identity_key"] == "ik14:AAAAAAAAAAAAAA"
    vetoed = next(item for item in entity["candidates"] if item["identity_key"] == "ik14:DDDDDDDDDDDDDD")
    assert vetoed["vetoed"] is True
    assert vetoed["chemistry_gate"]["status"] == "contradicted"
    assert vetoed["chemistry_gate"]["candidate_conditioned"] is False


def test_candidate_library_adduct_mismatch_cannot_veto_cross_adduct_identity():
    formula = "C42H72O14"
    precursor = ion_mz(formula_exact_mass(formula), "[M+FA-H]-")
    assert precursor is not None
    candidate = _candidate(
        "same-neutral-identity",
        IK_A,
        0.90,
        formula=formula,
        adduct="[M-H]-",
    )

    result = build_flagship_result(
        "cross-adduct",
        [_row(
            20,
            [candidate],
            precursor_mz=precursor,
            ion_mode="NEG",
        )],
    )
    chemistry = _winner(result)["chemistry_gate"]

    assert chemistry["veto"] is False
    assert chemistry["status"] == "consistent"
    assert chemistry["adduct"] == "[M+FA-H]-"
    assert chemistry["adduct_source"] == "candidate_formula_mass_unique_inference"
    assert chemistry["reported_library_adduct"] == "[M-H]-"
    assert chemistry["library_adduct_audit"] == "mass_mismatch_not_vetoed"
    assert chemistry["library_adduct_can_veto"] is False


def test_flagship_chemistry_uses_ms1_search_mass_and_audits_selected_ion():
    formula = "C6H12O6"
    measured = ion_mz(formula_exact_mass(formula), "[M+H]+")
    assert measured is not None
    candidate = _candidate(
        "glucose",
        IK_A,
        0.90,
        formula=formula,
        adduct="[M+H]+",
    )

    result = build_flagship_result(
        "ms1-mass-axis",
        [_row(
            20,
            [candidate],
            precursor_mz=181.0689,
            search_precursor_mz=measured,
            search_precursor_mz_source=(
                "measured_continuous_ms1_feature_centroid"
            ),
            query_adduct="[M+H]+",
            query_adduct_reliable=True,
        )],
    )
    entity = result["entities"][0]
    chemistry = _winner(result)["chemistry_gate"]

    assert chemistry["status"] == "consistent"
    assert abs(chemistry["mass_error_ppm"]) < 1.0
    assert entity["representative_precursor_mz"] == pytest.approx(measured)
    assert entity["representative_search_precursor_mz"] == pytest.approx(measured)
    assert entity["representative_selected_ion_mz"] == 181.0689
    assert entity["lcms_feature_used_for_search"] is True


def test_missing_chemistry_is_neutral_to_entity_score():
    complete = _candidate("a", IK_A, 0.80, formula="C8H10N4O2", adduct="")
    missing = _candidate("a", IK_A, 0.80, formula="", adduct="")
    result_complete = build_flagship_result("neutral-complete", [_row(10, [complete])])
    result_missing = build_flagship_result("neutral-missing", [_row(10, [missing])])

    assert _winner(result_complete)["entity_score"] == _winner(result_missing)["entity_score"]
    assert _winner(result_missing)["chemistry_gate"]["veto"] is False


def test_stereochemical_variants_collapse_to_connectivity():
    variant_1 = _candidate("stereo-one", IK_A, 0.85)
    variant_2 = _candidate("stereo-two", IK_A_STEREO, 0.83)
    result = build_flagship_result("stereo-collapse", [_row(10, [variant_1, variant_2])])
    entity = result["entities"][0]

    assert entity["candidate_count"] == 1
    assert entity["winner"]["identity_key"] == "ik14:AAAAAAAAAAAAAA"
    assert entity["winner"]["stereochemistry_claimed"] is False


def test_reference_source_families_are_counted_without_accession_pseudoreplication():
    source_a = _candidate("stereo-one", IK_A, 0.85)
    source_a["dataset_source"] = "MassBank_NIST-Pos"
    source_a["lib_spectrum_id"] = "source-a-1"
    source_a_repeat = copy.deepcopy(source_a)
    source_a_repeat["lib_spectrum_id"] = "source-a-2"
    source_b = _candidate("stereo-two", IK_A_STEREO, 0.83)
    source_b["dataset_source"] = "MoNA-export"

    result = build_flagship_result(
        "source-family-collapse",
        [_row(10, [source_a, source_a_repeat, source_b])],
    )
    winner = _winner(result)

    assert winner["reference_source_family_count"] == 2
    assert winner["reference_source_families"] == ["massbank", "mona"]
    assert winner["reference_source_accessions_count_as_independent"] is False
    assert winner["score_components"]["reference_source_diversity_affects_score"] is False


def test_tautomeric_structure_encodings_collapse_to_standardized_connectivity():
    keto_enol_variant = _candidate(
        "quercetin-systematic",
        "UMVWYQXKBPJMOF-UHFFFAOYSA-N",
        0.88,
        formula="C15H10O7",
    )
    keto_enol_variant["smiles"] = (
        "O=C1C(=O)C(OC=2C=C(O)C=C(O)C12)C3=CC=C(O)C(O)=C3"
    )
    canonical_variant = _candidate(
        "Quercetin",
        "REFJWTPEDVJJIY-UHFFFAOYSA-N",
        0.86,
        formula="C15H10O7",
    )
    canonical_variant["smiles"] = "O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12"

    result = build_flagship_result(
        "tautomer-collapse",
        [_row(20, [keto_enol_variant, canonical_variant])],
    )
    entity = result["entities"][0]

    assert entity["candidate_count"] == 1
    assert entity["winner"]["identity_key"] == "std2d:REFJWTPEDVJJIY"
    assert entity["winner"]["candidate"]["inchikey"].startswith("REFJWTPEDVJJIY-")
    assert entity["winner"]["candidate"]["standardized_connectivity_inchikey"] == "REFJWTPEDVJJIY"


def test_frozen_name_formula_registry_repairs_conflicting_library_structure():
    mislabeled = _candidate(
        "Glycyrrhizic acid",
        "QCQGNTXXPCJRFZ-UHFFFAOYSA-N",
        0.90,
        formula="C42H62O16",
    )
    mislabeled["smiles"] = (
        "CC1(C(=O)O)CCC2(C)CCC3(C)C(=C2C1)CC(=O)C1C2(C)CCC"
        "(OC4OC(C(=O)O)C(O)C(O)C4OC4OC(C(=O)O)C(O)C(O)C4O)"
        "C(C)(C)C2CCC13C"
    )

    result = build_flagship_result(
        "registry-repair",
        [_row(20, [mislabeled], precursor_mz=823.4111)],
    )
    winner = _winner(result)
    display = winner["candidate"]

    assert winner["identity_key"] == "std2d:LPLVUJXQOOQHMX"
    assert winner["identity_source"] == "pubchem_exact_name_formula_registry"
    assert display["inchikey"] == "LPLVUJXQOOQHMX-QWBHMCJMSA-N"
    assert display["reported_inchikey_before_identity_resolution"].startswith("QCQGNTXXPCJRFZ-")
    assert display["canonical_identity_resolution"]["ground_truth_inputs_used"] is False
    assert result["invariants"]["identity_registry_uses_ground_truth_inputs"] is False


def test_v14_registry_preserves_licoricesaponin_g2_connectivity():
    candidate = _candidate(
        "Licoricesaponin G2",
        "WBQVRPYEEYUEBQ-OJVDLISWSA-N",
        0.94,
        formula="C42H62O17",
    )
    candidate["smiles"] = (
        "C[C@]12CC[C@@](C)(C[C@H]1C1=CC(=O)[C@@H]3[C@@]4(C)CC[C@H]"
        "(O[C@@H]5O[C@@H]([C@@H](O)[C@H](O)[C@H]5O[C@@H]5O[C@@H]"
        "([C@@H](O)[C@H](O)[C@H]5O)C(O)=O)C(O)=O)[C@](C)(CO)[C@@H]4"
        "CC[C@@]3(C)[C@]1(C)CC2)C(O)=O"
    )

    result = build_flagship_result(
        "v14-licoricesaponin-g2",
        [_row(30, [candidate], precursor_mz=839.40598)],
    )
    winner = _winner(result)
    display = winner["candidate"]

    assert winner["identity_key"] == "std2d:WBQVRPYEEYUEBQ"
    assert winner["identity_source"] == "pubchem_exact_name_formula_registry"
    assert display["inchikey"] == "WBQVRPYEEYUEBQ-OJVDLISWSA-N"
    assert display["standardized_connectivity_inchikey"] == "WBQVRPYEEYUEBQ"
    assert display[
        "ranker_standardized_connectivity_before_identity_resolution"
    ] == "WIGUUKSQPAOWFS"
    assert display["canonical_identity_resolution"]["pubchem_cid"] == 14891565
    assert display["canonical_identity_resolution"]["ground_truth_inputs_used"] is False


def test_same_connectivity_display_keeps_stronger_resolved_library_candidate():
    glycyrrhizic = _candidate(
        "Glycyrrhizic acid",
        "LPLVUJXQOOQHMX-QWBHMCJMSA-N",
        0.94,
        formula="C42H62O16",
    )
    uralsaponin = _candidate(
        "Uralsaponin A",
        "LPLVUJXQOOQHMX-MOGLOQIBSA-N",
        0.90,
        formula="C42H62O16",
    )

    result = build_flagship_result(
        "v14-same-connectivity-display",
        [_row(30, [glycyrrhizic, uralsaponin], precursor_mz=823.4111)],
    )
    entity = result["entities"][0]
    winner = _winner(result)

    assert entity["candidate_count"] == 1
    assert winner["identity_key"] == "std2d:LPLVUJXQOOQHMX"
    assert winner["candidate"]["compound_name"] == "Glycyrrhizic acid"
    assert winner["candidate"]["inchikey"] == "LPLVUJXQOOQHMX-QWBHMCJMSA-N"


def test_entity_ranking_is_input_order_invariant():
    rows = [
        _row(10, [_candidate("a", IK_A, 0.82), _candidate("b", IK_B, 0.71)]),
        _row(20, [_candidate("a", IK_A, 0.79), _candidate("b", IK_B, 0.73)]),
        _row(30, [_candidate("a", IK_A, 0.84), _candidate("b", IK_B, 0.70)]),
    ]
    forward = build_flagship_result("order", rows)["entities"][0]
    reverse = build_flagship_result("order", list(reversed(rows)))["entities"][0]

    assert forward["winner"]["identity_key"] == reverse["winner"]["identity_key"]
    assert forward["winner"]["entity_score"] == reverse["winner"]["entity_score"]
    assert forward["decision"]["scope"] == reverse["decision"]["scope"]


def test_materialization_does_not_mutate_legacy_rows(tmp_path):
    rows = [_row(10, [_candidate("a", IK_A, 0.82)])]
    before = copy.deepcopy(rows)
    artifact = tmp_path / "task_v11_flagship.json"

    summary = write_flagship_artifact("artifact", rows, artifact)
    page = read_flagship_artifact(artifact, offset=0, limit=10)

    assert rows == before
    assert summary["legacy_result_mutated"] is False
    assert summary["production_identity_selector"] is True
    assert page["entities_total"] == 1
    assert page["entities"][0]["winner"]["identity_key"] == "ik14:AAAAAAAAAAAAAA"
    assert not list(tmp_path.glob(".*.tmp-*"))
    json.loads(artifact.read_text(encoding="utf-8"))


def test_member_spectrum_ids_remain_positionally_bound_to_result_indices():
    candidate = _candidate("a", IK_A, 0.82)
    rows = [
        _row(10, [candidate], spectrum_id="z-last-lexically"),
        _row(20, [candidate], spectrum_id="a-first-lexically"),
    ]

    entity = build_flagship_result("member-binding", rows)["entities"][0]

    assert entity["member_result_indices"] == [0, 1]
    assert entity["member_spectrum_ids"] == [
        "z-last-lexically",
        "a-first-lexically",
    ]


def test_fragmentation_output_never_contains_mass_only_mechanism():
    formula = "C6H12O6"
    precursor = ion_mz(formula_exact_mass(formula), "[M+H]+")
    glucose = _candidate(
        "formula-test",
        IK_A,
        0.90,
        formula=formula,
        adduct="[M+H]+",
    )
    peaks = [[163.0601, 100.0], [145.0495, 60.0], [127.0390, 30.0]]

    result = build_flagship_result(
        "fragmentation",
        [_row(20, [glucose], precursor_mz=precursor, raw_peaks=peaks)],
    )
    explanation = _winner(result)["fragmentation_explanation"]

    assert explanation["formula_conserving_only"] is True
    assert explanation["mass_only_mechanisms_emitted"] == 0
    assert all(edge["kind"] == "formula_conserving" for edge in explanation.get("edges", []))
    assert all(edge["formula_conservation"] is True for edge in explanation.get("edges", []))


def test_policy_rejects_probability_semantics():
    result = build_flagship_result("semantics", [_row(10, [_candidate("a", IK_A, 0.82)])])

    assert result["metadata"]["score_is_probability"] is False
    assert result["invariants"]["expert_score_affects_ranking"] is False
    assert result["invariants"]["mass_only_neutral_loss_called_mechanism"] is False
