"""Regression for real Unknown-to-Medium inflation and discarded scan scores."""
from copy import deepcopy

import pytest

from pipeline.compound_conclusions import (
    _project_primary_identity, _selected_spectrum_evidence,
)


def fixture(grade="Unknown"):
    reference = {
        "compound_name": "Y-33075", "inchikey": "JTVBXQAYBIJXRP-SNVBAGLBSA-N",
        "lib_spectrum_id": "ALL_GNPS_140106", "formula": "C16H16N4O",
        "smiles": "C[C@@H](N)c1ccc(C(=O)Nc2c3cc[nH]c3ncc2)cc1",
        "confidence": grade, "cosine": .364, "mass_ppm": 5.05,
        "entropy": .48, "final_score": .42, "rank_score": .42,
    }
    row = {"spectrum_id": "649148", "precursor_mz": 281.141418457031,
           "match": {**reference, "top_candidates": [deepcopy(reference)]}}
    candidate = {"identity_key": "std2d:JTVBXQAYBIJXRP", "entity_rank": 1,
                 "entity_score": .640674, "candidate": {
                     k: reference[k] for k in ["compound_name", "inchikey", "formula",
                                                "smiles", "lib_spectrum_id"]}}
    entity = {"entity_id": "entity:case", "winner": candidate, "candidates": [candidate],
              "decision": {"status": "selected_candidate_only",
                           "scope": "single_collision_energy_candidate"}}
    return row, entity


@pytest.mark.parametrize("raw,expected", [("Unknown", "Unknown"), ("Low", "Low"),
                                          ("Medium", "Medium"), ("High", "Medium")])
def test_candidate_only_never_promotes_saved_spectral_grade(raw, expected):
    row, entity = fixture(raw)
    before = deepcopy((row, entity))
    result = _project_primary_identity(row, entity, result_index=400)
    match = result["match"]
    assert match["confidence"] == expected
    assert match["cosine"] == .364 and match["mass_ppm"] == 5.05
    assert match["rank_score"] == .640674  # unchanged entity score
    binding = match["spectral_evidence_binding"]
    assert binding["metrics"]["rank_score"] == .42  # separate native scan score
    assert binding["spectrum_id"] == "649148" and binding["result_index"] == 400
    assert binding["saved_confidence"] == raw
    assert (row, entity) == before


@pytest.mark.parametrize("field,value", [("inchikey", "OTHERABCDEFGHIJ-UHFFFAOYSA-N"),
                                         ("lib_spectrum_id", "another_reference")])
def test_different_candidate_or_reference_cannot_supply_scores(field, value):
    row, entity = fixture()
    entity["winner"]["candidate"][field] = value
    result = _project_primary_identity(row, entity, result_index=400)["match"]
    assert "cosine" not in result and "mass_ppm" not in result
    assert result["confidence"] == "Unknown"
    assert result["spectral_evidence_binding"]["status"] == "unavailable"


@pytest.mark.parametrize("field,value", [("spectrum_id", "wrong_scan"), ("result_index", 401)])
def test_explicit_cross_scan_reference_is_rejected(field, value):
    row, entity = fixture()
    row["match"][field] = value
    row["match"]["top_candidates"][0][field] = value
    source, binding = _selected_spectrum_evidence(row, entity["winner"]["candidate"], 400)
    assert source is None and binding["status"] == "unavailable"


def test_selected_candidate_in_same_scan_shortlist_recovers_own_evidence_only():
    row, entity = fixture("Low")
    actual = deepcopy(row["match"]["top_candidates"][0])
    row["match"].update(compound_name="Different Top1", inchikey="OTHERABCDEFGHIJ-UHFFFAOYSA-N",
                        lib_spectrum_id="other", cosine=.99, mass_ppm=.01,
                        confidence="High", evidence_tier="Strong Candidate")
    row["match"]["top_candidates"] = [actual]
    result = _project_primary_identity(row, entity, result_index=400)["match"]
    assert result["cosine"] == .364 and result["mass_ppm"] == 5.05
    assert result["confidence"] == "Low" and result["evidence_tier"] == ""
    assert result["spectral_evidence_binding"]["source_field"] == "match.top_candidates[0]"


def test_disagreeing_duplicate_reference_metrics_are_not_cherry_picked():
    row, entity = fixture()
    candidate = deepcopy(row["match"]["top_candidates"][0])
    other = {**candidate, "cosine": .8}
    row["match"].update(inchikey="OTHERABCDEFGHIJ-UHFFFAOYSA-N", lib_spectrum_id="other")
    row["match"]["top_candidates"] = [candidate, other]
    source, binding = _selected_spectrum_evidence(row, entity["winner"]["candidate"], 400)
    assert source is None and binding["reason"] == "conflicting_saved_reference_metrics"


def test_bound_scores_do_not_override_existing_high_evidence_or_primary_ranking():
    row, entity = fixture("High")
    proof = {"scope": "2d_connectivity", "supporting_spectra": [{"spectrum_id": "649148"}]}
    projected = _project_primary_identity(row, entity, proof, result_index=400)
    assert projected["match"]["confidence"] == "High"
    assert projected["match"]["flagship_entity_score"] == .640674
    assert projected["match"]["entity_high_evidence"] == proof
