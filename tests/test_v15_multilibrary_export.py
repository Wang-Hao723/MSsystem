from __future__ import annotations

import json

from web.v15_multilibrary_export import (
    MULTILIBRARY_COLUMNS,
    multilibrary_evidence_rows,
    read_multilibrary_evidence,
)


def _entity():
    evidence = {
        "legacy_best_reference_score": 0.91,
        "provenance_consensus_score": 0.82,
        "conditioned_multilibrary_score": 0.865,
        "display_only_unconditioned_consensus_score": 0.80,
        "unconditioned_consensus_affects_ranking": False,
        "rank_influence_allowed": True,
        "rank_influence_basis": (
            "compatible_instrument_domain_near_ce_stratum_median"
        ),
        "score_is_probability": False,
        "library_evidence_matrix": [
            {
                "source_family": "massbank",
                "search_status": "matched_conditioned",
                "source_participation_explicit": True,
                "reference_record_count": 2,
                "condition_compatible_record_count": 1,
                "unique_content_profile_count": 2,
                "independent_single_source_profile_count": 1,
                "best_score": 0.91,
                "best_condition_compatible_score": 0.88,
                "representative_spectrum_id": "MB-1",
                "representative_source_accession": "ACCESSION-1",
                "rank_influence_eligible": True,
                "representative_condition": {
                    "tier": "same_model_near_ce",
                    "reference_instrument_family": "qtof",
                    "reference_activation": "beam_cid",
                    "reference_collision_energy": 20.0,
                    "collision_energy_unit": "eV",
                    "collision_energy_delta": 0.0,
                },
            }
        ],
    }
    candidate = {
        "identity_key": "ik14:AAAAAAAAAAAAAA",
        "entity_rank": 1,
        "entity_score": 0.86,
        "legacy_v14_entity_score": 0.84,
        "candidate": {
            "compound_name": "Candidate A",
            "formula": "C10H10O2",
            "inchikey": "AAAAAAAAAAAAAA-BBBBBBBBBB-C",
        },
        "collision_energy_evidence": [
            {
                "collision_energy_bin": "eV:20",
                "conditioned_multilibrary_evidence": evidence,
            }
        ],
    }
    return {
        "entity_id": "entity-1",
        "decision": {"status": "selected", "scope": "2d_connectivity_candidate"},
        "winner": candidate,
        "candidates": [candidate],
    }


def test_multilibrary_export_is_one_row_per_entity_candidate_ce_library():
    rows = multilibrary_evidence_rows([_entity()])

    assert len(rows) == 1
    assert list(rows[0]) == MULTILIBRARY_COLUMNS
    assert rows[0]["Library_Source_Family"] == "massbank"
    assert rows[0]["Representative_Source_Accession"] == "ACCESSION-1"
    assert rows[0]["Rank_Influence_Eligible"] is True
    assert rows[0]["Unconditioned_Consensus_Affects_Ranking"] is False
    assert rows[0]["Score_Is_Probability"] is False


def test_multilibrary_json_reader_retains_explicit_non_probability_semantics(
    tmp_path,
):
    path = tmp_path / "task_v15_flagship.json"
    path.write_text(json.dumps({"entities": [_entity()]}), encoding="utf-8")

    payload = read_multilibrary_evidence(path)

    assert payload["schema_version"].endswith("/1.0")
    assert len(payload["rows"]) == 1
    assert payload["score_is_probability"] is False
    assert payload["unconditioned_consensus_affects_ranking"] is False
