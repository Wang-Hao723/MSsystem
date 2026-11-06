from __future__ import annotations

import copy

import pytest

from matching.conditioned_multilibrary import (
    aggregate_references,
    reference_condition,
    reference_content_fingerprint,
)
from matching.flagship_ranker import build_flagship_result
from matching.library_quality import formula_exact_mass
from matching.v15_flagship_ranker import (
    build_v15_flagship_result,
    load_v15_flagship_policy,
)
from utils.adduct import ion_mz


def _reference(
    score: float,
    source: str,
    peaks: list[list[float]],
    *,
    spectrum_id: str,
    instrument_type: str = "LC-ESI-QTOF",
    instrument: str = "Agilent 6546 Q-TOF",
    collision_energy: float = 20.0,
    collision_energy_unit: str = "eV",
    collision_energy_mode: str = "fixed",
    activation_method: str = "CID",
) -> dict:
    return {
        "spectral_score": score,
        "score_components": {"cosine": score},
        "dataset_source": source,
        "source_tier": "experimental",
        "lib_spectrum_id": spectrum_id,
        "lib_source_accession": f"ACC-{spectrum_id}",
        "lib_peaks": peaks,
        "lib_instrument_type": instrument_type,
        "lib_instrument": instrument,
        "lib_collision_energy": collision_energy,
        "lib_collision_energy_unit": collision_energy_unit,
        "lib_collision_energy_mode": collision_energy_mode,
        "lib_activation_method": activation_method,
    }


POLICY = {
    "best_reference_fraction": 0.5,
    "consensus_fraction": 0.5,
    "minimum_provenance_strata": 2,
    "minimum_unique_content_profiles": 2,
    "near_collision_energy_delta": {"eV": 5.0, "NCE": 5.0},
}


def test_peak_identical_cross_library_mirrors_create_one_profile_not_two_votes():
    query = {
        "instrument_type": "qtof",
        "instrument": "Agilent 6546 Q-TOF",
        "collision_energy": 20,
        "collision_energy_unit": "eV",
        "collision_energy_mode": "fixed",
        "activation_method": "CID",
    }
    mirror = [[50.0, 10.0], [75.0, 100.0], [90.0, 25.0]]
    distinct = [[50.0, 20.0], [75.0, 60.0], [91.0, 100.0]]
    result = aggregate_references(
        query,
        [
            _reference(0.92, "MassBank_NIST", mirror, spectrum_id="m1"),
            _reference(0.91, "ALL_GNPS", mirror, spectrum_id="g1"),
            _reference(0.70, "MoNA-export", distinct, spectrum_id="o1"),
        ],
        POLICY,
    )

    assert result["reference_record_count"] == 3
    assert result["unique_content_profile_count"] == 2
    assert result["mirror_duplicate_record_count"] == 1
    assert result["provenance_stratum_count"] == 2
    assert result["conditioned_provenance_stratum_count"] == 2
    assert result["conditioned_multilibrary_score"] == pytest.approx(0.865)
    assert result["mirrored_sources_count_as_independent_evidence"] is False


def test_unconditioned_multilibrary_median_is_display_only_and_cannot_rerank():
    query = {
        "instrument_type": "qtof",
        "instrument": "Agilent 6546 Q-TOF",
        "collision_energy": 20,
        "collision_energy_unit": "eV",
        "collision_energy_mode": "fixed",
        "activation_method": "CID",
    }
    references = [
        _reference(
            0.90,
            "MassBank_NIST",
            [[50.0, 10.0], [75.0, 100.0], [90.0, 25.0]],
            spectrum_id="far-1",
            collision_energy=40.0,
        ),
        _reference(
            0.50,
            "ALL_GNPS",
            [[51.0, 15.0], [76.0, 100.0], [92.0, 20.0]],
            spectrum_id="far-2",
            collision_energy=45.0,
        ),
    ]

    result = aggregate_references(query, references, POLICY)

    assert result["unconditioned_corroboration_available"] is True
    assert result["display_only_unconditioned_consensus_score"] == pytest.approx(0.70)
    assert result["unconditioned_consensus_affects_ranking"] is False
    assert result["rank_influence_allowed"] is False
    assert result["conditioned_multilibrary_score"] == pytest.approx(0.90)
    assert result["selected_consensus_basis"] == (
        "legacy_best_reference_no_condition_qualified_corroboration"
    )


def test_library_matrix_separates_hit_no_hit_and_unknown_search_state():
    query = {
        "instrument_type": "qtof",
        "instrument": "Agilent 6546 Q-TOF",
        "collision_energy": 20,
        "collision_energy_unit": "eV",
        "collision_energy_mode": "fixed",
        "activation_method": "CID",
        "searched_source_families": ["massbank", "gnps", "mona"],
    }
    policy = {
        **POLICY,
        "report_source_families": ["massbank", "gnps", "mona", "nist"],
    }
    result = aggregate_references(
        query,
        [
            _reference(
                0.82,
                "MassBank_NIST",
                [[50.0, 10.0], [75.0, 100.0], [90.0, 25.0]],
                spectrum_id="massbank-hit",
                collision_energy=40.0,
            ),
            _reference(
                0.79,
                "ALL_GNPS",
                [[51.0, 15.0], [76.0, 100.0], [92.0, 20.0]],
                spectrum_id="gnps-hit",
            ),
        ],
        policy,
    )
    matrix = {
        row["source_family"]: row
        for row in result["library_evidence_matrix"]
    }

    assert matrix["massbank"]["search_status"] == (
        "matched_condition_incompatible_or_incomplete"
    )
    assert matrix["gnps"]["search_status"] == "matched_conditioned"
    assert matrix["mona"]["search_status"] == "no_candidate_reference"
    assert matrix["nist"]["search_status"] == "search_state_unknown"
    assert matrix["mona"]["absence_is_negative_evidence"] is False
    assert matrix["gnps"]["representative_source_accession"] == (
        "ACC-gnps-hit"
    )


def test_ce_is_compared_only_within_unit_and_same_analyzer_family():
    query = {
        "instrument_type": "qtof",
        "instrument": "Agilent 6546 Q-TOF",
        "collision_energy": 20,
        "collision_energy_unit": "eV",
        "collision_energy_mode": "fixed",
        "activation_method": "CID",
    }
    near = reference_condition(
        query,
        _reference(
            0.8,
            "MassBank_NIST",
            [[50, 1], [80, 2]],
            spectrum_id="near",
            collision_energy=24,
            collision_energy_unit="eV",
        ),
        near_delta_by_unit={"eV": 5.0, "NCE": 5.0},
    )
    unlike_unit = reference_condition(
        query,
        _reference(
            0.8,
            "MassBank_NIST",
            [[50, 1], [80, 2]],
            spectrum_id="nce",
            collision_energy=20,
            collision_energy_unit="NCE",
        ),
        near_delta_by_unit={"eV": 5.0, "NCE": 5.0},
    )
    cross_platform = reference_condition(
        query,
        _reference(
            0.8,
            "MassBank_NIST",
            [[50, 1], [80, 2]],
            spectrum_id="orbitrap",
            instrument_type="Orbitrap",
            instrument="Q Exactive HF",
            collision_energy=20,
            collision_energy_unit="eV",
        ),
        near_delta_by_unit={"eV": 5.0, "NCE": 5.0},
    )

    assert near["tier"] == "same_model_near_ce"
    assert near["eligible_for_conditioned_consensus"] is True
    assert unlike_unit["collision_energy_delta"] is None
    assert unlike_unit["cross_unit_conversion_applied"] is False
    assert unlike_unit["eligible_for_conditioned_consensus"] is False
    assert cross_platform["tier"] == "cross_analyzer_same_ce_unit"
    assert cross_platform["eligible_for_conditioned_consensus"] is False

    cross_vendor_qtof = reference_condition(
        query,
        _reference(
            0.8,
            "MassBank_NIST",
            [[50, 1], [80, 2]],
            spectrum_id="waters",
            instrument_type="QTOF",
            instrument="Waters SYNAPT G2-S",
            collision_energy=20,
            collision_energy_unit="eV",
        ),
        near_delta_by_unit={"eV": 5.0, "NCE": 5.0},
    )
    assert cross_vendor_qtof["same_instrument_family"] is True
    assert cross_vendor_qtof["platform_anchor_eligible"] is False
    assert cross_vendor_qtof["eligible_for_conditioned_consensus"] is False


def test_stepped_and_fixed_ce_are_not_reduced_to_one_number():
    query = {
        "instrument_type": "qtof",
        "instrument": "Agilent 6546 Q-TOF",
        "activation_method": "CID",
        "collision_energy": 20,
        "collision_energy_unit": "eV",
        "collision_energy_mode": "stepped",
    }
    condition = reference_condition(
        query,
        _reference(
            0.8,
            "MassBank_NIST",
            [[50, 1], [80, 2]],
            spectrum_id="fixed",
            collision_energy=20,
            collision_energy_unit="eV",
            collision_energy_mode="fixed",
        ),
        near_delta_by_unit={"eV": 5.0},
    )

    assert condition["collision_energy_delta"] == 0
    assert condition["collision_energy_mode_compatible"] is False
    assert condition["eligible_for_conditioned_consensus"] is False
    assert condition["composite_energy_reduced_to_midpoint"] is False


def _candidate(
    identity: str,
    score: float,
    source: str,
    suffix: int,
) -> dict:
    return {
        "compound_name": f"display-{identity}-{suffix}",
        "inchikey": f"{identity}-UHFFFAOYSA-N",
        "raw_cosine": score,
        "raw_entropy": score,
        "raw_reverse_match": score,
        "raw_fragment_score": score,
        "raw_intensity_score": score,
        "mass_score": score,
        "dataset_source": source,
        "source_tier": "experimental",
        "lib_spectrum_id": f"{identity}-{suffix}",
        "lib_source_accession": f"ACC-{identity}-{suffix}",
        "lib_peaks": [
            [50.0 + suffix * 0.01, 10.0 + suffix],
            [70.0 + suffix * 0.01, 100.0],
            [90.0 + suffix * 0.01, 20.0 + suffix],
        ],
        "lib_instrument_type": "qtof",
        "lib_instrument": "Agilent 6546 Q-TOF",
        "lib_collision_energy": 20.0,
        "lib_collision_energy_unit": "eV",
        "lib_collision_energy_mode": "fixed",
        "lib_activation_method": "CID",
    }


def test_v15_can_rerank_from_robust_references_and_retains_v14_audit_score():
    identity_a = "AAAAAAAAAAAAAA"
    identity_b = "BBBBBBBBBBBBBB"
    match = _candidate(identity_a, 0.90, "MassBank_NIST", 1)
    match["top_candidates"] = [
        _candidate(identity_a, 0.40, "ALL_GNPS", 2),
        _candidate(identity_a, 0.40, "MoNA-export", 3),
        _candidate(identity_b, 0.84, "MassBank_NIST", 4),
        _candidate(identity_b, 0.84, "ALL_GNPS", 5),
        _candidate(identity_b, 0.84, "MoNA-export", 6),
    ]
    row = {
        "spectrum_id": "query-1",
        "lcms_feature_id": "feature-1",
        "precursor_mz": 100.0,
        "search_precursor_mz": 100.0,
        "ion_mode": "POS",
        "instrument_type": "qtof",
        "instrument": "Agilent 6546 Q-TOF",
        "collision_energy": 20.0,
        "collision_energy_unit": "eV",
        "collision_energy_mode": "fixed",
        "activation_method": "CID",
        "raw_peaks": [[50, 10], [70, 100], [90, 20]],
        "match": match,
    }

    v14 = build_flagship_result("audit", [copy.deepcopy(row)])
    v15 = build_v15_flagship_result("production", [copy.deepcopy(row)])
    before = v14["entities"][0]
    after = v15["entities"][0]

    assert before["winner"]["identity_key"].endswith(identity_a)
    assert after["winner"]["identity_key"].endswith(identity_b)
    assert after["ranking_changed_from_v14"] is True
    assert after["v14_audit_baseline"]["identity_key"].endswith(identity_a)
    winner = after["winner"]
    assert winner["legacy_v14_entity_score"] != winner["entity_score"]
    evidence = winner["collision_energy_evidence"][0][
        "conditioned_multilibrary_evidence"
    ]
    assert evidence["reference_record_count"] == 3
    assert evidence["selected_consensus_basis"] == (
        "compatible_instrument_domain_near_ce_stratum_median"
    )
    assert all(
        "compound_name" not in reference
        for reference in evidence["reference_records"]
    )
    assert v15["metadata"]["module_id"] == (
        "v15_conditioned_multilibrary_flagship"
    )
    assert v15["metadata"]["entities_reranked_from_v14"] == 1


def test_reference_fingerprint_is_scale_invariant():
    first = [[50, 1], [75, 10], [90, 2]]
    second = [[50, 10], [75, 100], [90, 20]]
    assert reference_content_fingerprint(first) == reference_content_fingerprint(second)


def test_feature_fragment_payload_is_display_only_and_cannot_change_v15_rank():
    identity_a = "AAAAAAAAAAAAAA"
    identity_b = "BBBBBBBBBBBBBB"
    match = _candidate(identity_a, 0.86, "MassBank_NIST", 1)
    match["top_candidates"] = [
        _candidate(identity_b, 0.82, "ALL_GNPS", 2),
    ]
    row = {
        "spectrum_id": "query-feature-shadow",
        "lcms_feature_id": "feature-shadow",
        "precursor_mz": 100.0,
        "search_precursor_mz": 100.0,
        "ion_mode": "POS",
        "instrument_type": "qtof",
        "instrument": "Agilent 6546 Q-TOF",
        "collision_energy": 20.0,
        "collision_energy_unit": "eV",
        "collision_energy_mode": "fixed",
        "activation_method": "CID",
        "raw_peaks": [[50, 10], [70, 100], [90, 20]],
        "match": match,
    }
    baseline = build_v15_flagship_result("ff-baseline", [copy.deepcopy(row)])

    row["match"]["feature_fragment_evidence"] = {
        "status": "available",
        "likelihood_ratio": 1e200,
        "combined_log_likelihood_ratio": 460.0,
        "runtime_influence": False,
        "identity_rank_influence": False,
    }
    row["match"]["top_candidates"][0]["feature_fragment_evidence"] = {
        "status": "available",
        "likelihood_ratio": 1e-200,
        "combined_log_likelihood_ratio": -460.0,
        "runtime_influence": False,
        "identity_rank_influence": False,
    }
    changed = build_v15_flagship_result("ff-mutated", [copy.deepcopy(row)])

    base_entity = baseline["entities"][0]
    changed_entity = changed["entities"][0]
    assert changed_entity["winner"]["identity_key"] == base_entity["winner"]["identity_key"]
    assert [item["identity_key"] for item in changed_entity["candidates"]] == [
        item["identity_key"] for item in base_entity["candidates"]
    ]
    assert [item["entity_score"] for item in changed_entity["candidates"]] == [
        item["entity_score"] for item in base_entity["candidates"]
    ]
    assert changed["invariants"]["feature_fragment_evidence_affects_ranking"] is False


def _chemistry_rows(
    formula: str,
    precursor_mz: float,
    *,
    ion_identity: str = "",
    ion_identity_neutral_mass: float | None = None,
    ion_identity_networked: bool = False,
    lcms_adduct_group_id: str = "",
    lcms_adduct_group_evidence: dict | None = None,
) -> list[dict]:
    candidate = _candidate("CHEMISTRYABCDE", 0.92, "MassBank_NIST", 30)
    candidate.update({
        "formula": formula,
        "adduct": "[M+FA-H]-",
        "has_structure": True,
    })
    rows = []
    for ce in (35.0, 40.0, 45.0):
        match = copy.deepcopy(candidate)
        match["top_candidates"] = [copy.deepcopy(candidate)]
        rows.append({
            "spectrum_id": f"chemistry-{ce:g}",
            "lcms_feature_id": "feature-chemistry",
            "precursor_mz": precursor_mz,
            "search_precursor_mz": precursor_mz,
            "ion_mode": "NEG",
            "collision_energy": ce,
            "collision_energy_unit": "eV",
            "lcms_adduct": "[M-H]-",
            "lcms_neutral_mass": precursor_mz + 1.00727646677,
            "lcms_adduct_group_id": lcms_adduct_group_id,
            "lcms_adduct_group_evidence": (
                copy.deepcopy(lcms_adduct_group_evidence) or {}
            ),
            "ion_identity": ion_identity,
            "ion_identity_status": (
                "networked" if ion_identity_networked else "unlinked"
            ),
            "ion_identity_network_id": (
                "IIN-TEST-1" if ion_identity_networked else ""
            ),
            "ion_identity_neutral_mass": ion_identity_neutral_mass,
            "ion_identity_evidence": ({
                "status": "networked",
                "identity_conflict": False,
                "runtime_influence": True,
            } if ion_identity_networked else {}),
            "raw_peaks": [[80.0, 20.0], [120.0, 100.0], [180.0, 40.0]],
            "match": match,
        })
    return rows


def test_v15_unlinked_default_adduct_cannot_hard_veto_formate_candidate():
    formula = "C42H72O14"
    precursor = ion_mz(formula_exact_mass(formula), "[M+FA-H]-")
    assert precursor is not None
    rows = _chemistry_rows(formula, precursor)

    legacy = build_flagship_result("legacy-default-adduct", copy.deepcopy(rows))
    v15 = build_v15_flagship_result("v15-default-adduct", copy.deepcopy(rows))

    assert legacy["entities"][0]["decision"]["status"] == "abstained"
    assert v15["entities"][0]["decision"]["status"] == "selected"
    chemistry = v15["entities"][0]["winner"]["chemistry_gate"]
    assert chemistry["veto"] is False
    assert chemistry["adduct"] == "[M+FA-H]-"
    assert chemistry["adduct_source"] == (
        "candidate_formula_mass_unique_inference"
    )
    assert chemistry["unlinked_lcms_feature_adduct_ignored"] is True
    assert chemistry["unlinked_lcms_feature_neutral_mass_ignored"] is True


def test_v15_networked_ion_identity_overrides_conflicting_feature_default():
    formula = "C27H42O20"
    neutral = formula_exact_mass(formula)
    precursor = ion_mz(neutral, "[M+FA-H]-")
    assert neutral is not None and precursor is not None
    rows = _chemistry_rows(
        formula,
        precursor,
        ion_identity="[M+FA-H]-",
        ion_identity_neutral_mass=neutral,
        ion_identity_networked=True,
    )

    v15 = build_v15_flagship_result("v15-ion-network", rows)
    chemistry = v15["entities"][0]["winner"]["chemistry_gate"]

    assert chemistry["veto"] is False
    assert chemistry["adduct"] == "[M+FA-H]-"
    assert chemistry["adduct_source"] == "ion_identity_network_adduct"
    assert chemistry["ion_identity_network_evidence_used"] is True
    assert chemistry["candidate_conditioned"] is False


def test_v15_unvalidated_adduct_group_cannot_hard_veto():
    formula = "C42H72O14"
    precursor = ion_mz(formula_exact_mass(formula), "[M+FA-H]-")
    assert precursor is not None
    rows = _chemistry_rows(
        formula,
        precursor,
        lcms_adduct_group_id="ADDUCT-GROUP-1",
    )

    v15 = build_v15_flagship_result("v15-explicit-group", rows)
    entity = v15["entities"][0]

    assert entity["decision"]["status"] == "selected"
    chemistry = entity["winner"]["chemistry_gate"]
    assert chemistry["veto"] is False
    assert chemistry["unvalidated_lcms_adduct_group_ignored"] is True
    assert chemistry["lcms_validated_adduct_group_evidence_available"] is False


def test_v15_validated_adduct_group_remains_eligible_for_hard_veto():
    formula = "C42H72O14"
    precursor = ion_mz(formula_exact_mass(formula), "[M+FA-H]-")
    assert precursor is not None
    rows = _chemistry_rows(
        formula,
        precursor,
        lcms_adduct_group_id="ADDUCT-GROUP-1",
        lcms_adduct_group_evidence={
            "status": "validated",
            "neutral_mass_consistent": True,
            "identity_conflict": False,
            "member_feature_count": 2,
            "runtime_influence": True,
        },
    )

    v15 = build_v15_flagship_result("v15-validated-group", rows)
    entity = v15["entities"][0]

    assert entity["decision"]["status"] == "abstained"
    candidate = entity["candidates"][0]
    assert candidate["chemistry_gate"]["veto"] is True
    assert candidate["chemistry_gate"][
        "lcms_validated_adduct_group_evidence_available"
    ] is True


def _neutral_loss_ambiguity_rows(
    *,
    candidate_adduct: str,
) -> list[dict]:
    formula = "C17H26O4"
    neutral = formula_exact_mass(formula)
    precursor = ion_mz(neutral, candidate_adduct)
    assert neutral is not None and precursor is not None
    network_neutral = precursor - 1.007276
    candidate = _candidate(
        "GINGEROLABCDEF", 0.95, "MassBank_NIST", 41
    )
    candidate.update({
        "formula": formula,
        "adduct": candidate_adduct,
        "has_structure": True,
    })
    rows = []
    for ce in (10.0, 20.0, 30.0):
        match = copy.deepcopy(candidate)
        match["top_candidates"] = [copy.deepcopy(candidate)]
        rows.append({
            "spectrum_id": f"neutral-loss-{candidate_adduct}-{ce:g}",
            "lcms_feature_id": "feature-neutral-loss",
            "precursor_mz": precursor,
            "search_precursor_mz": precursor,
            "ion_mode": "POS",
            "collision_energy": ce,
            "collision_energy_unit": "eV",
            "ion_identity": "[M+H]+",
            "ion_identity_status": "networked",
            "ion_identity_network_id": "IIN-NEUTRAL-LOSS-1",
            "ion_identity_neutral_mass": network_neutral,
            "ion_identity_evidence": {
                "status": "networked",
                "identity_conflict": False,
                "runtime_influence": True,
            },
            "raw_peaks": [[80.0, 20.0], [120.0, 100.0]],
            "match": match,
        })
    return rows


def test_v15_mass_consistent_library_neutral_loss_preserves_candidate_but_demotes_claim():
    result = build_v15_flagship_result(
        "v15-neutral-loss-ambiguity",
        _neutral_loss_ambiguity_rows(candidate_adduct="[M-H2O+H]+"),
    )
    entity = result["entities"][0]
    chemistry = entity["winner"]["chemistry_gate"]

    assert entity["winner"]["vetoed"] is False
    assert chemistry["status"] == "ambiguous"
    assert chemistry["precursor_transformation_ambiguity"] is True
    assert chemistry["candidate_independent_contradiction_preserved"] is True
    assert chemistry["precursor_transformation_identity_confirmed"] is False
    assert entity["decision"]["status"] == "selected_candidate_only"
    assert entity["decision"]["scope"] == (
        "precursor_transformation_ambiguity_candidate"
    )


def test_v15_ordinary_adduct_conflict_cannot_use_transformation_exception():
    result = build_v15_flagship_result(
        "v15-ordinary-adduct-conflict",
        _neutral_loss_ambiguity_rows(candidate_adduct="[M+Na]+"),
    )
    entity = result["entities"][0]
    candidate = entity["candidates"][0]

    assert entity["winner"] is None
    assert candidate["vetoed"] is True
    assert candidate["chemistry_gate"].get(
        "precursor_transformation_ambiguity", False
    ) is False


def _abstention_rows(
    runner_score: float,
    *,
    leader_score: float = 0.95,
) -> list[dict]:
    leader = _candidate(
        "AAAAAAAAAAAAAA", leader_score, "MassBank_NIST", 51
    )
    runner = _candidate("BBBBBBBBBBBBBB", runner_score, "MassBank_NIST", 52)
    rows = []
    for ce in (20.0, 30.0):
        match = copy.deepcopy(leader)
        match["top_candidates"] = [
            copy.deepcopy(leader),
            copy.deepcopy(runner),
        ]
        rows.append(
            {
                "spectrum_id": f"abstention-{ce:g}",
                "lcms_feature_id": "feature-abstention",
                "precursor_mz": 300.0,
                "search_precursor_mz": 300.0,
                "ion_mode": "POS",
                "collision_energy": ce,
                "collision_energy_unit": "eV",
                "raw_peaks": [[80.0, 20.0], [120.0, 100.0], [180.0, 40.0]],
                "match": match,
            }
        )
    return rows


def test_cross_platform_abstention_demotes_release_without_changing_rank_or_score():
    rows = _abstention_rows(0.20, leader_score=0.40)
    ungated_policy = load_v15_flagship_policy()
    ungated_policy["cross_platform_abstention_gate"]["enabled"] = False
    ungated = build_flagship_result("ungated-abstention", rows, policy=ungated_policy)
    gated = build_v15_flagship_result("gated-abstention", rows)
    before = ungated["entities"][0]
    after = gated["entities"][0]

    assert before["decision"]["status"] == "selected"
    assert after["decision"]["status"] == "selected_candidate_only"
    assert after["decision"]["scope"] == "cross_platform_abstention_candidate"
    assert after["decision"]["cross_platform_abstention_gate"]["status"] == "demoted"
    assert after["winner"]["identity_key"] == before["winner"]["identity_key"]
    assert after["winner"]["entity_score"] == before["winner"]["entity_score"]
    assert [item["identity_key"] for item in after["candidates"]] == [
        item["identity_key"] for item in before["candidates"]
    ]


def test_cross_platform_abstention_retains_well_separated_multi_ce_release():
    result = build_v15_flagship_result(
        "passed-abstention",
        _abstention_rows(0.40),
    )
    entity = result["entities"][0]

    assert entity["decision"]["status"] == "selected"
    assert entity["decision"]["scope"] == "2d_connectivity_candidate"
    assert entity["decision"]["cross_platform_abstention_gate"]["status"] == "passed"


def test_query_spectrum_gate_demotes_two_ce_release_with_too_few_fragments():
    rows = _abstention_rows(0.40)
    for row in rows:
        row["raw_peaks"] = [[80.0, 20.0], [120.0, 100.0]]
    policy = load_v15_flagship_policy()
    policy["query_spectrum_release_gate"]["affects_release_scope"] = True
    result = build_flagship_result(
        "low-fragment-evidence", rows, policy=policy
    )
    entity = result["entities"][0]
    audit = entity["decision"]["query_spectrum_release_gate"]

    assert entity["decision"]["status"] == "selected_candidate_only"
    assert entity["decision"]["scope"] == "query_spectrum_quality_candidate"
    assert audit["status"] == "demoted"
    assert audit["observed_quality_qualified_ce_bins"] == 0
    assert entity["winner"]["identity_key"] == "ik14:AAAAAAAAAAAAAA"
    assert result["metadata"][
        "entities_demoted_by_query_spectrum_release_gate"
    ] == 1


def test_query_spectrum_gate_requires_accepted_deconvolution_for_mixed_scan():
    rows = _abstention_rows(0.40)
    for row in rows:
        row["match"]["chimeric_deconvolution"] = {
            "isolation_status": "mixed",
            "status": "poor_reconstruction",
        }
    policy = load_v15_flagship_policy()
    policy["query_spectrum_release_gate"]["affects_release_scope"] = True
    demoted = build_flagship_result(
        "mixed-not-deconvolved", rows, policy=policy
    )
    demoted_entity = demoted["entities"][0]

    assert demoted_entity["decision"]["scope"] == (
        "query_spectrum_quality_candidate"
    )
    assert demoted_entity["decision"]["query_spectrum_release_gate"][
        "failure_reason_counts"
    ]["mixed_or_severe_isolation_without_accepted_deconvolution"] == 2

    for row in rows:
        row["match"]["chimeric_deconvolution"]["status"] = (
            "accepted_score_applied"
        )
    passed = build_flagship_result(
        "mixed-qualified-deconvolution", rows, policy=policy
    )
    passed_entity = passed["entities"][0]
    assert passed_entity["decision"]["status"] == "selected"
    assert passed_entity["decision"]["query_spectrum_release_gate"][
        "status"
    ] == "passed"


def test_query_spectrum_gate_reads_precursor_reliability_from_mass_audit():
    rows = _abstention_rows(0.40)
    rows[0]["precursor_mass_audit"] = {
        "search_precursor_mz_reliable": False
    }
    rows[1]["precursor_mass_audit"] = {
        "search_precursor_mz_reliable": True
    }
    policy = load_v15_flagship_policy()
    policy["query_spectrum_release_gate"]["affects_release_scope"] = True
    result = build_flagship_result("one-unreliable-ce", rows, policy=policy)
    entity = result["entities"][0]
    audit = entity["decision"]["query_spectrum_release_gate"]

    assert entity["decision"]["scope"] == "query_spectrum_quality_candidate"
    assert audit["observed_quality_qualified_ce_bins"] == 1
    first = audit["collision_energy_evidence"][0]["scans"][0]
    assert first["search_precursor_mz_reliability_source"] == (
        "precursor_mass_audit"
    )


def test_production_query_quality_gate_is_audit_only_after_transfer_failure():
    rows = _abstention_rows(0.40)
    for row in rows:
        row["match"]["chimeric_deconvolution"] = {
            "isolation_status": "severe",
            "status": "poor_reconstruction",
        }
    result = build_v15_flagship_result("audit-only-query-quality", rows)
    entity = result["entities"][0]
    audit = entity["decision"]["query_spectrum_release_gate"]

    assert entity["decision"]["status"] == "selected"
    assert entity["decision"]["scope"] == "2d_connectivity_candidate"
    assert audit["status"] == "failed_audit_only"
    assert audit["affects_release_scope"] is False
    assert result["metadata"][
        "entities_demoted_by_query_spectrum_release_gate"
    ] == 0
    assert result["metadata"][
        "entities_failed_query_spectrum_audit_only"
    ] == 1
