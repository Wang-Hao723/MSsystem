from __future__ import annotations

import copy

import pytest

from matching.v16_dual_polarity import (
    ALGORITHM_ID,
    MODULE_ID,
    build_v16_dual_polarity_shadow,
    load_v16_dual_polarity_policy,
    rt_alignment_semantic_sha256,
    semantic_sha256,
    sidecar_semantic_sha256,
)


KEY_A = "AAAAAAAAAAAAAA"
KEY_B = "BBBBBBBBBBBBBB"
KEY_C = "CCCCCCCCCCCCCC"


def _task(
    task_id: str,
    polarity: str,
    *,
    status: str = "completed",
    filename: str | None = None,
) -> dict:
    return {
        "id": task_id,
        "status": status,
        "ion_mode": polarity,
        "filename": filename or f"{task_id}.mzML",
        "result_run_id": f"run-{task_id}",
        "untargeted_workflow": {
            "feature_xml": f"/shadow/{task_id}.featureXML"
        },
    }


def _candidate(
    connectivity: str | None,
    *,
    formula: str = "C10H12O2",
    name: str = "Database candidate",
    source: str = "test-library",
    potential_contaminant: bool = False,
    compound_role: str = "",
) -> dict:
    candidate = {
        "compound_name": name,
        "formula": formula,
        "dataset_source": source,
        "lib_source_accession": f"ACC-{connectivity or 'NONE'}",
        "final_score": 0.9,
    }
    if connectivity:
        candidate["inchikey"] = f"{connectivity}-UHFFFAOYSA-N"
    if potential_contaminant:
        candidate["potential_contaminant"] = True
    if compound_role:
        candidate["compound_role"] = compound_role
    return candidate


def _row(
    polarity: str,
    group_id: str,
    *,
    neutral_mass: float = 300.0,
    rt: float = 5.0,
    feature_id: str | None = None,
    connectivity: str | None = KEY_A,
    formula: str = "C10H12O2",
    name: str = "Database candidate",
    anchor_source: str = "iimn",
    conflict: bool = False,
    member_count: int | None = 2,
    member_ids: list[str] | None = None,
    adduct: str | None = None,
    top40: list[dict] | None = None,
) -> dict:
    explicit_adduct = adduct
    if explicit_adduct is None:
        explicit_adduct = "[M+H]+" if polarity == "POS" else "[M-H]-"
    row = {
        "spectrum_id": f"scan-{polarity}-{group_id}",
        "lcms_feature_id": feature_id or f"feature-{polarity}-{group_id}",
        "precursor_feature_id": f"event-{polarity}-{group_id}",
        "precursor_mz": neutral_mass
        + (1.007276 if polarity == "POS" else -1.007276),
        "retention_time": rt,
        "retention_time_unit": "minutes",
        "ion_mode": polarity,
        "match": _candidate(connectivity, formula=formula, name=name),
    }
    if anchor_source == "iimn":
        evidence = {
            "status": "networked",
            "network_id": group_id,
            "neutral_mass": neutral_mass,
            "identity_conflict": conflict,
            "runtime_influence": True,
        }
        if member_count is not None:
            evidence["member_feature_count"] = member_count
        if member_ids is not None:
            evidence["member_feature_ids"] = list(member_ids)
        row.update(
            {
                "ion_identity_status": "networked",
                "ion_identity_network_id": group_id,
                "ion_identity": explicit_adduct,
                "ion_identity_neutral_mass": neutral_mass,
                "ion_identity_evidence": evidence,
            }
        )
    elif anchor_source in {"openms", "openms_unvalidated"}:
        evidence = {
            "status": "validated" if anchor_source == "openms" else "hypothesis",
            "neutral_mass_consistent": True,
            "identity_conflict": conflict,
            "runtime_influence": True,
        }
        if member_count is not None:
            evidence["member_feature_count"] = member_count
        if member_ids is not None:
            evidence["member_feature_ids"] = list(member_ids)
        row.update(
            {
                "lcms_adduct_group_id": group_id,
                "lcms_adduct": explicit_adduct,
                "lcms_neutral_mass": neutral_mass,
                "lcms_adduct_group_evidence": evidence,
            }
        )
    elif anchor_source == "unlinked":
        row.update(
            {
                "ion_identity_status": "unlinked",
                "lcms_neutral_mass": neutral_mass,
                "lcms_adduct": explicit_adduct,
            }
        )
    else:
        raise AssertionError(anchor_source)
    if top40:
        row["candidate_recall_support"] = {
            "expanded_candidates": copy.deepcopy(top40)
        }
    return row


def _baseline(
    *,
    positive_rows: list[dict] | None = None,
    negative_rows: list[dict] | None = None,
) -> tuple[dict, list[dict], dict, list[dict]]:
    return (
        _task("task-pos", "POS"),
        positive_rows or [_row("POS", "IIN-POS")],
        _task("task-neg", "NEG"),
        negative_rows or [_row("NEG", "IIN-NEG")],
    )


def _build(
    *,
    positive_rows: list[dict] | None = None,
    negative_rows: list[dict] | None = None,
    **kwargs,
) -> dict:
    pos_task, pos_rows, neg_task, neg_rows = _baseline(
        positive_rows=positive_rows,
        negative_rows=negative_rows,
    )
    return build_v16_dual_polarity_shadow(
        pos_task,
        pos_rows,
        neg_task,
        neg_rows,
        **kwargs,
    )


def _flagship_entity(
    entity_id: str,
    row_indices: list[int],
    *,
    winner: dict | None = None,
    candidates: list[dict] | None = None,
    expanded: list[dict] | None = None,
) -> dict:
    entity = {
        "entity_id": entity_id,
        "member_result_indices": list(row_indices),
        "winner": copy.deepcopy(winner),
        "candidates": copy.deepcopy(candidates or []),
    }
    if expanded:
        entity["candidate_recall_support"] = {
            "expanded_candidates": copy.deepcopy(expanded)
        }
    return entity


def test_policy_is_locked_shadow_and_freezes_new_safety_rules():
    policy = load_v16_dual_polarity_policy()

    assert policy["module_id"] == MODULE_ID
    assert policy["algorithm_id"] == ALGORITHM_ID
    assert policy["execution"]["mode"] == "shadow_only"
    assert policy["cross_polarity_edge"]["maximum_rt_delta_seconds"] == 20.0
    assert policy["cross_polarity_edge"]["pairing_strategy"] == (
        "degree_one_unique_only"
    )
    assert policy["cross_polarity_edge"]["one_to_many_action"] == (
        "ambiguous_no_merge"
    )
    assert policy["candidate_independent_anchor"]["iimn_required"][
        "minimum_member_feature_count"
    ] == 2
    assert policy["candidate_support"][
        "candidate_attachment_occurs_after_graph_freeze"
    ] is True

    changed = copy.deepcopy(policy)
    changed["cross_polarity_edge"]["maximum_rt_delta_seconds"] = 21.0
    with pytest.raises(ValueError, match="20 seconds"):
        _build(policy=changed)


def test_pairing_requires_distinct_completed_explicit_pos_and_neg_tasks():
    pos = _task("explicit-pos", "POS", filename="looks-like-NEG.mzML")
    neg = _task("explicit-neg", "NEG", filename="looks-like-POS.mzML")
    payload = build_v16_dual_polarity_shadow(
        pos, [_row("POS", "P")], neg, [_row("NEG", "N")]
    )
    assert payload["metadata"]["task_pair"] == {
        "POS": "explicit-pos",
        "NEG": "explicit-neg",
        "explicit_selection_required": True,
        "filename_used_for_pairing": False,
    }

    with pytest.raises(ValueError, match="distinct explicit task IDs"):
        build_v16_dual_polarity_shadow(
            _task("same", "POS"),
            [_row("POS", "P")],
            _task("same", "NEG"),
            [_row("NEG", "N")],
        )
    with pytest.raises(ValueError, match="completed"):
        build_v16_dual_polarity_shadow(
            _task("pending", "POS", status="running"),
            [_row("POS", "P")],
            _task("neg", "NEG"),
            [_row("NEG", "N")],
        )
    with pytest.raises(ValueError, match="exactly one POS and one NEG"):
        build_v16_dual_polarity_shadow(
            _task("p1", "POS"),
            [_row("POS", "P1")],
            _task("p2", "POS"),
            [_row("POS", "P2")],
        )


def test_task_ion_mode_is_the_only_role_source_and_rows_only_reject_conflicts():
    legacy_only = _task("legacy-pos", "POS")
    legacy_only.pop("ion_mode")
    legacy_only["polarity"] = "POS"
    with pytest.raises(ValueError, match=r"task\['ion_mode'\]"):
        build_v16_dual_polarity_shadow(
            legacy_only,
            [_row("POS", "P")],
            _task("neg", "NEG"),
            [_row("NEG", "N")],
        )

    conflicting = _row("POS", "P")
    conflicting["ion_mode"] = "NEG"
    with pytest.raises(ValueError, match="result-row polarity conflicts"):
        build_v16_dual_polarity_shadow(
            _task("pos", "POS"),
            [conflicting],
            _task("neg", "NEG"),
            [_row("NEG", "N")],
        )


@pytest.mark.parametrize(
    ("positive_rows", "expected"),
    [
        ([_row("POS", "P", anchor_source="unlinked")], 0),
        ([_row("POS", "P", anchor_source="openms_unvalidated")], 0),
        ([_row("POS", "P", anchor_source="openms")], 1),
        ([_row("POS", "P", member_count=None, member_ids=None)], 0),
        (
            [
                _row(
                    "POS",
                    "P",
                    member_count=None,
                    member_ids=["f1", "f2"],
                )
            ],
            1,
        ),
        (
            [
                _row("POS", "P", feature_id="f1", member_count=None),
                _row("POS", "P", feature_id="f2", member_count=None),
            ],
            1,
        ),
    ],
)
def test_only_reliable_multi_feature_iimn_or_validated_openms_creates_anchor(
    positive_rows,
    expected,
):
    payload = _build(positive_rows=positive_rows)
    assert payload["summary"]["neutral_entity_count"] == expected
    assert payload["summary"]["compound_count_increment"] == 0


@pytest.mark.parametrize("anchor_source", ["iimn", "openms"])
@pytest.mark.parametrize("adduct", ["", "[M-H]-"])
def test_anchor_requires_nonempty_explicit_adduct_with_task_compatible_sign(
    anchor_source,
    adduct,
):
    payload = _build(
        positive_rows=[
            _row("POS", "P", anchor_source=anchor_source, adduct=adduct)
        ]
    )
    assert payload["summary"]["neutral_entity_count"] == 0
    assert payload["summary"]["reliable_anchor_counts"]["POS"] == 0


def test_candidate_identity_conflict_is_recorded_after_graph_and_never_blocks_iimn():
    clean = _build(positive_rows=[_row("POS", "P", conflict=False)])
    conflict = _build(positive_rows=[_row("POS", "P", conflict=True)])

    assert conflict["summary"]["neutral_entity_count"] == 1
    assert clean["candidate_independent_graph"]["semantic_sha256"] == conflict[
        "candidate_independent_graph"
    ]["semantic_sha256"]
    evidence = conflict["neutral_entities"][0]["precursor_event_evidence"]
    assert next(item for item in evidence if item["polarity"] == "POS")[
        "candidate_identity_conflict_observed"
    ] is True
    assert conflict["invariants"][
        "candidate_identity_conflict_used_for_graph_membership"
    ] is False


@pytest.mark.parametrize(
    ("negative_row", "expected"),
    [
        (_row("NEG", "N", neutral_mass=300.0031, rt=5.0), 0),
        (_row("NEG", "N", neutral_mass=300.0029, rt=5.0), 1),
        (_row("NEG", "N", rt=5.0 + 20.0 / 60.0), 1),
        (_row("NEG", "N", rt=5.0 + 20.001 / 60.0), 0),
    ],
)
def test_neutral_mass_and_fixed_twenty_second_rt_gate(negative_row, expected):
    payload = _build(negative_rows=[negative_row])
    assert payload["summary"]["neutral_entity_count"] == expected
    if expected:
        edge = payload["neutral_entities"][0]["candidate_independent_edge"]
        assert edge["rt_tolerance_seconds"] == 20.0
        assert edge["candidate_identity_used"] is False
        assert edge["spectrum_similarity_used"] is False


def test_one_to_many_edges_are_all_ambiguous_and_never_nearest_selected():
    payload = _build(
        positive_rows=[_row("POS", "P", neutral_mass=300.0, rt=5.0)],
        negative_rows=[
            _row("NEG", "N1", neutral_mass=300.0001, rt=5.0),
            _row("NEG", "N2", neutral_mass=300.0002, rt=5.01),
        ],
    )

    assert payload["summary"]["candidate_independent_compatible_edge_count"] == 2
    assert payload["summary"]["candidate_independent_accepted_edge_count"] == 0
    assert payload["summary"]["ambiguous_edge_count"] == 2
    assert payload["summary"]["ambiguous_anchor_count"] == 3
    assert payload["neutral_entities"] == []
    assert payload["ambiguity"]["automatic_merge_allowed"] is False
    assert {
        edge["pairing_status"] for edge in payload["ambiguity"]["edges"]
    } == {"ambiguous_no_merge"}
    assert {
        edge["positive_anchor_degree"] for edge in payload["ambiguity"]["edges"]
    } == {2}


def test_input_order_and_candidate_display_names_do_not_change_frozen_graph():
    pos_task, pos_rows, neg_task, neg_rows = _baseline()
    forward = build_v16_dual_polarity_shadow(pos_task, pos_rows, neg_task, neg_rows)
    reverse = build_v16_dual_polarity_shadow(neg_task, neg_rows, pos_task, pos_rows)
    assert forward == reverse

    renamed_pos = copy.deepcopy(pos_rows)
    renamed_neg = copy.deepcopy(neg_rows)
    renamed_pos[0]["match"]["compound_name"] = "完全不同的显示名称"
    renamed_neg[0]["match"]["compound_name"] = "another display label"
    renamed = build_v16_dual_polarity_shadow(
        pos_task, renamed_pos, neg_task, renamed_neg
    )
    assert renamed["candidate_independent_graph"] == forward[
        "candidate_independent_graph"
    ]
    assert renamed["summary"] == forward["summary"]
    assert renamed["invariants"][
        "candidate_display_name_used_for_structure_identity"
    ] is False


def test_optional_rt_alignment_is_hash_bound_and_can_supply_aligned_rt():
    pos_task, pos_rows, neg_task, neg_rows = _baseline(
        positive_rows=[_row("POS", "P", rt=5.0)],
        negative_rows=[_row("NEG", "N", rt=5.5)],
    )
    raw = build_v16_dual_polarity_shadow(pos_task, pos_rows, neg_task, neg_rows)
    assert raw["summary"]["neutral_entity_count"] == 0
    alignment = {
        "status": "ok",
        "task_ids": ["task-pos", "task-neg"],
        "source_results_sha256": raw["metadata"][
            "source_results_semantic_sha256"
        ],
        "source_task_metadata_sha256": raw["metadata"][
            "source_task_metadata_sha256"
        ],
        "alignment_status": {"task-pos": "reference", "task-neg": "aligned"},
        "aligned_rt_minutes": {
            "task-pos": {"P": 5.10},
            "task-neg": {"N": 5.12},
        },
    }
    alignment["semantic_sha256"] = rt_alignment_semantic_sha256(alignment)
    aligned = build_v16_dual_polarity_shadow(
        pos_task, pos_rows, neg_task, neg_rows, rt_alignment=alignment
    )
    assert aligned["summary"]["neutral_entity_count"] == 1
    assert {
        item["rt_basis"]
        for item in aligned["neutral_entities"][0]["precursor_event_evidence"]
    } == {"hash_bound_rt_alignment"}

    wrong_source = copy.deepcopy(alignment)
    wrong_source.pop("semantic_sha256")
    wrong_source["source_results_sha256"]["task-neg"] = "0" * 64
    wrong_source["semantic_sha256"] = rt_alignment_semantic_sha256(wrong_source)
    with pytest.raises(ValueError, match="different source results"):
        build_v16_dual_polarity_shadow(
            pos_task, pos_rows, neg_task, neg_rows, rt_alignment=wrong_source
        )


def test_only_stable_exact_2d_connectivity_can_receive_dual_support():
    exact = _build()
    support = exact["neutral_entities"][0]["candidate_support"]
    assert support[0]["connectivity"] == KEY_A
    assert support[0]["dual_polarity_support"] is True

    same_name_different_structure = _build(
        positive_rows=[_row("POS", "P", connectivity=KEY_A, name="Same name")],
        negative_rows=[_row("NEG", "N", connectivity=KEY_B, name="Same name")],
    )
    assert all(
        item["dual_polarity_support"] is False
        for item in same_name_different_structure["neutral_entities"][0][
            "candidate_support"
        ]
    )

    formula_conflict = _build(
        positive_rows=[_row("POS", "P", connectivity=KEY_A, formula="C10H12O2")],
        negative_rows=[_row("NEG", "N", connectivity=KEY_A, formula="C11H14O2")],
    )
    conflict = formula_conflict["neutral_entities"][0]["candidate_support"][0]
    assert conflict["formula_conflict"] is True
    assert conflict["dual_polarity_support"] is False
    assert conflict["status"] == "cross_polarity_formula_conflict"


def test_single_polarity_top40_candidate_stays_support_only():
    expanded = {
        "connectivity_identity": f"std2d:{KEY_B}",
        "formula": "C10H12O2",
        "compound_name": "Expanded candidate",
        "dataset_source": "rrf40",
    }
    payload = _build(positive_rows=[_row("POS", "P", top40=[expanded])])
    support = {
        item["connectivity"]: item
        for item in payload["neutral_entities"][0]["candidate_support"]
    }
    assert support[KEY_B]["status"] == "single_polarity_support"
    assert support[KEY_B]["polarities_observed"] == ["POS"]
    assert support[KEY_B]["can_replace_winner"] is False
    assert support[KEY_B]["can_receive_2d_connectivity_decision"] is False
    assert support[KEY_B]["compound_count_increment"] == 0


def test_v15_flagship_candidates_attach_by_member_result_indices_after_graph_freeze():
    artifacts = {
        "task-pos": {
            "entities": [
                _flagship_entity(
                    "fp",
                    [0],
                    winner=_candidate(KEY_A, name="POS A"),
                    expanded=[
                        {
                            "identity_key": f"std2d:{KEY_B}",
                            "candidate": _candidate(KEY_B, name="POS B"),
                        }
                    ],
                )
            ]
        },
        "task-neg": {
            "entities": [
                _flagship_entity(
                    "fn", [0], winner=_candidate(KEY_A, name="NEG A")
                )
            ]
        },
    }
    original = copy.deepcopy(artifacts)
    payload = _build(candidate_support_by_task=artifacts)
    support = {
        item["connectivity"]: item
        for item in payload["neutral_entities"][0]["candidate_support"]
    }

    assert support[KEY_A]["dual_polarity_support"] is True
    assert support[KEY_B]["status"] == "single_polarity_support"
    assert support[KEY_B]["polarities_observed"] == ["POS"]
    assert payload["neutral_entities"][0]["candidate_support_attachment"] == {
        "state": "attached_after_candidate_independent_graph_freeze",
        "source": "v15_flagship_member_result_indices",
        "binding_field": "member_result_indices",
        "candidate_fields_used_for_graph": False,
    }
    provenance = support[KEY_A]["polarity_provenance"]["POS"]["provenance"]
    assert provenance[0]["flagship_entity_id"] == "fp"
    assert provenance[0]["member_result_indices"] == [0]
    assert artifacts == original


def test_contaminant_and_exogenous_flagship_candidates_are_qc_only():
    clean = _candidate(KEY_A, name="Plant-compatible candidate")
    contaminant = _candidate(
        KEY_B,
        name="Plasticizer record",
        potential_contaminant=True,
    )
    exogenous = _candidate(
        KEY_C,
        name="Industrial record",
        compound_role="exogenous_chemical",
    )
    artifacts = {
        "task-pos": {
            "entities": [
                _flagship_entity(
                    "fp",
                    [0],
                    winner=clean,
                    candidates=[contaminant, exogenous],
                )
            ]
        },
        "task-neg": {
            "entities": [_flagship_entity("fn", [0], winner=clean)]
        },
    }
    payload = _build(candidate_support_by_task=artifacts)
    entity = payload["neutral_entities"][0]

    assert {item["connectivity"] for item in entity["candidate_support"]} == {KEY_A}
    qc = entity["quality_control_excluded_candidates"]
    assert {item["connectivity"] for item in qc} == {KEY_B, KEY_C}
    reasons = {reason for item in qc for reason in item["exclusion_reasons"]}
    assert "potential_contaminant=true" in reasons
    assert "compound_role=exogenous_chemical" in reasons
    assert all(item["retained_for_quality_control_only"] is True for item in qc)
    assert payload["summary"]["quality_control_excluded_candidate_count"] == 2


def test_perfect_ms2_similarity_is_auxiliary_and_cannot_create_edge_or_support():
    similarity = [
        {
            "positive_group_id": "P",
            "negative_group_id": "N",
            "modified_cosine": 1.0,
            "neutral_losses_cosine": 1.0,
            "ms2deepscore": 1.0,
        }
    ]
    mismatched = _build(
        positive_rows=[_row("POS", "P", neutral_mass=300.0, rt=5.0)],
        negative_rows=[_row("NEG", "N", neutral_mass=301.0, rt=7.0)],
        auxiliary_spectral_similarity=similarity,
    )
    assert mismatched["summary"]["neutral_entity_count"] == 0

    matched = _build(
        positive_rows=[_row("POS", "P")],
        negative_rows=[_row("NEG", "N")],
        auxiliary_spectral_similarity=similarity,
    )
    auxiliary = matched["neutral_entities"][0][
        "auxiliary_cross_polarity_ms2_similarity"
    ]
    assert auxiliary["state"] == "auxiliary_only"
    assert auxiliary["can_create_candidate_independent_edge"] is False
    assert auxiliary["can_create_dual_connectivity_support"] is False
    assert auxiliary["can_change_ranking"] is False
    assert auxiliary["score_is_probability"] is False


def test_sidecar_is_deterministic_hash_bound_read_only_and_never_double_counts():
    pos_task, pos_rows, neg_task, neg_rows = _baseline()
    inputs_before = copy.deepcopy((pos_task, pos_rows, neg_task, neg_rows))

    first = build_v16_dual_polarity_shadow(pos_task, pos_rows, neg_task, neg_rows)
    second = build_v16_dual_polarity_shadow(pos_task, pos_rows, neg_task, neg_rows)

    assert first == second
    assert (pos_task, pos_rows, neg_task, neg_rows) == inputs_before
    assert first["metadata"]["source_results_semantic_sha256"] == first[
        "metadata"
    ]["source_results_semantic_sha256_after"]
    assert first["metadata"]["sidecar_semantic_sha256"] == sidecar_semantic_sha256(
        first
    )
    graph = copy.deepcopy(first["candidate_independent_graph"])
    graph_hash = graph.pop("semantic_sha256")
    assert graph_hash == semantic_sha256(graph)
    assert graph_hash == first["metadata"][
        "candidate_independent_graph_sha256_before_candidate_attachment"
    ]
    assert graph_hash == first["metadata"][
        "candidate_independent_graph_sha256_after_candidate_attachment"
    ]
    assert first["invariants"][
        "candidate_independent_graph_unchanged_after_candidate_attachment"
    ] is True
    assert first["summary"]["status"] == "evaluated_with_edges"
    assert first["summary"]["accepted_neutral_entity_count"] == 1
    assert first["summary"]["ambiguous_edge_count"] == 0
    assert first["summary"]["ambiguous_anchor_count"] == 0
    assert first["summary"]["not_evaluable_anchor_count"] == 0
    assert first["summary"]["dual_polarity_supported_connectivity_count"] == 1
    assert first["summary"]["precursor_event_evidence_record_count"] == 2
    assert len(first["neutral_entities"][0]["precursor_event_evidence"]) == 2
    assert first["summary"]["compound_count_increment"] == 0
    assert all(
        first["invariants"][field] is False
        for field in (
            "production_top1_changed",
            "production_top5_changed",
            "production_score_changed",
            "production_decision_changed",
            "production_confidence_changed",
            "production_compound_count_changed",
        )
    )
