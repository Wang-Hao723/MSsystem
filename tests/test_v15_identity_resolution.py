from __future__ import annotations

import copy

from matching.v15_identity_resolution import apply_v15_identity_resolution


def _candidate(
    identity: str,
    formula: str,
    score: float,
    rank: int,
    *,
    ce: int = 2,
) -> dict:
    return {
        "identity_key": identity,
        "entity_rank": rank,
        "entity_score": score,
        "collision_energy_bins_supported": ce,
        "collision_energy_bins_total": ce,
        "vetoed": False,
        "candidate": {
            "compound_name": f"database-name-{identity}",
            "formula": formula,
            "inchikey": identity,
            "smiles": f"structure-{identity}",
        },
    }


def _payload(candidates: list[dict], *, bins: int = 2) -> dict:
    winner = copy.deepcopy(candidates[0])
    return {
        "metadata": {
            "task_id": "identity-test",
            "production_identity_selector": True,
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
                "entity_id": "entity-1",
                "distinct_collision_energy_bins": bins,
                "winner": winner,
                "leading_candidate_when_abstained": None,
                "candidates": copy.deepcopy(candidates),
                "decision": {
                    "status": "selected",
                    "scope": "2d_connectivity_candidate",
                    "blockers": [],
                },
            }
        ],
    }


def _rank_snapshot(payload: dict) -> list[tuple[str, int, float]]:
    return [
        (
            candidate["identity_key"],
            candidate["entity_rank"],
            candidate["entity_score"],
        )
        for candidate in payload["entities"][0]["candidates"]
    ]


def test_same_formula_competitive_connectivities_release_only_as_isomer_group():
    payload = _payload(
        [
            _candidate("AAAAAAAAAAAAAA", "C10H12O2", 0.90, 1),
            _candidate("BBBBBBBBBBBBBB", "C10H12O2", 0.68, 2),
        ]
    )
    before = _rank_snapshot(payload)

    resolved = apply_v15_identity_resolution(payload)
    entity = resolved["entities"][0]

    assert _rank_snapshot(resolved) == before
    assert entity["winner"]["identity_key"] == "AAAAAAAAAAAAAA"
    assert entity["decision"]["status"] == "selected_isomer_group"
    assert entity["decision"]["scope"] == "molecular_formula_isomer_group"
    assert entity["decision"]["reportable_identity_key"] is None
    assert entity["identity_resolution"]["candidate_formulae"] == ["C10H12O2"]


def test_cross_formula_competition_is_unresolved_and_not_silently_ranked_as_identity():
    payload = _payload(
        [
            _candidate("AAAAAAAAAAAAAA", "C12H18O4", 0.80, 1),
            _candidate("BBBBBBBBBBBBBB", "C12H16O3", 0.55, 2),
        ]
    )
    before = _rank_snapshot(payload)

    resolved = apply_v15_identity_resolution(payload)
    entity = resolved["entities"][0]

    assert _rank_snapshot(resolved) == before
    assert entity["decision"]["status"] == "unresolved_formula_event"
    assert entity["decision"]["scope"] == "unresolved_precursor_event"
    assert entity["decision"]["reportable_identity_key"] is None
    assert entity["winner"]["reporting_semantics"] == (
        "ranking_lead_not_reportable_compound_identity"
    )


def test_missing_formula_in_competitive_set_is_a_distinct_unresolved_state():
    payload = _payload(
        [
            _candidate("AAAAAAAAAAAAAA", "C8H8O", 0.82, 1),
            _candidate("BBBBBBBBBBBBBB", "", 0.70, 2),
        ]
    )

    resolved = apply_v15_identity_resolution(payload)
    resolution = resolved["entities"][0]["identity_resolution"]

    assert resolution["status"] == "unresolved_incomplete_formula_event"
    assert resolution["missing_formula_candidate_identity_keys"] == [
        "BBBBBBBBBBBBBB"
    ]


def test_every_non_vetoed_shortlist_entry_enters_scope_resolution():
    payload = _payload(
        [
            _candidate("AAAAAAAAAAAAAA", "C8H8O", 0.90, 1),
            _candidate("BBBBBBBBBBBBBB", "C9H10O2", 0.40, 2),
        ]
    )

    resolved = apply_v15_identity_resolution(payload)
    entity = resolved["entities"][0]

    assert entity["decision"]["status"] == "unresolved_formula_event"
    assert entity["decision"]["scope"] == "unresolved_precursor_event"
    assert entity["decision"]["reportable_identity_key"] is None
    assert entity["identity_resolution"]["competitive_candidate_count"] == 2
    assert entity["identity_resolution"]["thresholds"][
        "fixed_score_threshold_used"
    ] is False


def test_ranking_lead_with_one_ce_remains_in_cross_formula_resolution():
    payload = _payload(
        [
            _candidate("AAAAAAAAAAAAAA", "C13H18O5S", 0.80, 1, ce=1),
            _candidate("BBBBBBBBBBBBBB", "C16H14O5", 0.65, 2, ce=2),
            _candidate("CCCCCCCCCCCCCC", "C16H14O5", 0.58, 3, ce=2),
        ]
    )

    resolved = apply_v15_identity_resolution(payload)
    entity = resolved["entities"][0]
    resolution = entity["identity_resolution"]

    assert entity["winner"]["identity_key"] == "AAAAAAAAAAAAAA"
    assert entity["decision"]["status"] == "unresolved_formula_event"
    assert resolution["status"] == "unresolved_cross_formula_event"
    assert resolution["competitive_candidate_identity_keys"] == [
        "AAAAAAAAAAAAAA",
        "BBBBBBBBBBBBBB",
        "CCCCCCCCCCCCCC",
    ]
    assert resolution["thresholds"][
        "ranking_lead_always_included_in_conflict_resolution"
    ] is True


def test_candidate_names_cannot_change_formula_scope_resolution():
    candidates = [
        _candidate("AAAAAAAAAAAAAA", "C10H12O2", 0.90, 1),
        _candidate("BBBBBBBBBBBBBB", "C10H12O2", 0.68, 2),
    ]
    first = apply_v15_identity_resolution(_payload(candidates))
    renamed = copy.deepcopy(candidates)
    renamed[0]["candidate"]["compound_name"] = "completely different title"
    renamed[1]["candidate"]["compound_name"] = "另一个名称"
    second = apply_v15_identity_resolution(_payload(renamed))

    first_resolution = first["entities"][0]["identity_resolution"]
    second_resolution = second["entities"][0]["identity_resolution"]
    assert first_resolution["status"] == second_resolution["status"]
    assert first_resolution["competitive_candidate_identity_keys"] == (
        second_resolution["competitive_candidate_identity_keys"]
    )
