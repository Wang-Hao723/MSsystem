"""A truncated shortlist is unobserved evidence, not a measured zero match."""
from __future__ import annotations

import copy

import pytest

from matching.flagship_ranker import build_flagship_result
from matching.v15_flagship_ranker import load_v15_flagship_policy
from tests.test_v11_flagship_ranker import IK_A, IK_B, _candidate, _row


RULE = "observed_candidates_only_v16"


def policy(*, corrected=True, negative_absence=False):
    value = copy.deepcopy(load_v15_flagship_policy())
    value["candidate_pool"]["same_energy_replicate_support"] = RULE if corrected else "legacy_zero_fill"
    value["candidate_pool"]["missing_from_truncated_topk_is_negative_evidence"] = negative_absence
    return value


def candidate(payload, identity):
    return next(item for item in payload["entities"][0]["candidates"] if item["identity_key"] == identity)


@pytest.mark.parametrize("missing", [1, 2, 8])
def test_absent_same_energy_replicates_do_not_erase_measured_support(missing):
    a, b = _candidate("a", IK_A, .90), _candidate("b", IK_B, .60)
    one = [_row(10, [a, b], spectrum_id="observed")]
    rows = one + [_row(10, [b], spectrum_id=f"missing-{i}") for i in range(missing)]
    expected = candidate(build_flagship_result("fixed", one, policy=policy()), "ik14:AAAAAAAAAAAAAA")
    actual = candidate(build_flagship_result("fixed", rows, policy=policy()), "ik14:AAAAAAAAAAAAAA")
    assert actual["score_components"]["mean_collision_energy_support"] == expected["score_components"]["mean_collision_energy_support"]
    evidence = actual["collision_energy_evidence"][0]
    assert evidence["technical_replicates"] == missing + 1
    assert evidence["replicates_with_candidate"] == 1
    assert evidence["replicates_missing_from_truncated_topk"] == missing
    assert evidence["truncated_topk_absence_imputed_as_zero"] is False
    assert actual["collision_energy_bins_total"] == 1


def test_legacy_policy_remains_reproducible_and_v14_lane_is_not_rewritten():
    a, b = _candidate("a", IK_A, .90), _candidate("b", IK_B, .60)
    rows = [_row(10, [a, b]), _row(10, [b], spectrum_id="b2"), _row(10, [b], spectrum_id="b3")]
    old = candidate(build_flagship_result("fixed", rows, policy=policy(corrected=False)), "ik14:AAAAAAAAAAAAAA")
    new = candidate(build_flagship_result("fixed", rows, policy=policy()), "ik14:AAAAAAAAAAAAAA")
    assert old["entity_score"] == 0
    assert new["entity_score"] > .9
    assert new["legacy_v14_entity_score"] == old["legacy_v14_entity_score"] == 0


def test_explicit_negative_absence_policy_is_still_obeyed():
    a, b = _candidate("a", IK_A, .90), _candidate("b", IK_B, .60)
    rows = [_row(10, [a, b]), _row(10, [b]), _row(10, [b])]
    result = build_flagship_result("negative", rows, policy=policy(negative_absence=True))
    value = candidate(result, "ik14:AAAAAAAAAAAAAA")
    assert value["entity_score"] == 0
    assert value["collision_energy_evidence"][0]["truncated_topk_absence_imputed_as_zero"] is True


def test_complete_candidate_vectors_keep_existing_scores_order_and_decision():
    rows = [_row(ce, [_candidate("a", IK_A, .82), _candidate("b", IK_B, .60)], spectrum_id=f"s-{ce}-{i}")
            for ce in (10, 20, 30) for i in range(3)]
    before = build_flagship_result("complete", rows, policy=policy(corrected=False))["entities"][0]
    after = build_flagship_result("complete", rows, policy=policy())["entities"][0]
    project = lambda entity: [(c["identity_key"], c["entity_score"], c["legacy_v14_entity_score"]) for c in entity["candidates"]]
    assert project(after) == project(before)
    assert after["decision"] == before["decision"]


def test_true_low_measurements_still_lower_the_median():
    a = [_candidate("a", IK_A, strength) for strength in (.90, .10, .10)]
    rows = [_row(10, [item], spectrum_id=f"s-{index}") for index, item in enumerate(a)]
    result = build_flagship_result("observed-low", rows, policy=policy())
    value = candidate(result, "ik14:AAAAAAAAAAAAAA")
    assert value["collision_energy_evidence"][0]["support"] < .5
    assert value["collision_energy_evidence"][0]["replicates_missing_from_truncated_topk"] == 0


def test_candidate_names_and_input_order_cannot_control_the_fix():
    a, b = _candidate("a", IK_A, .90), _candidate("b", IK_B, .60)
    rows = [_row(10, [a, b], spectrum_id="a"), _row(10, [b], spectrum_id="b"), _row(10, [b], spectrum_id="c")]
    modified = copy.deepcopy(list(reversed(rows)))
    for row in modified:
        for item in [row["match"], *row["match"]["top_candidates"]]:
            item.update(compound_name="unrelated renamed candidate", expert_score=1.0, truth_connectivity="UNUSED")
    project = lambda payload: [(c["identity_key"], c["entity_score"]) for c in payload["entities"][0]["candidates"]]
    assert project(build_flagship_result("same", modified, policy=policy())) == project(build_flagship_result("same", rows, policy=policy()))


def test_one_ce_does_not_become_multiple_independent_conditions():
    a, b = _candidate("a", IK_A, .90), _candidate("b", IK_B, .60)
    rows = [_row(10, [a, b])] + [_row(10, [b], spectrum_id=f"s{i}") for i in range(10)]
    entity = build_flagship_result("one-ce", rows, policy=policy())["entities"][0]
    assert entity["distinct_collision_energy_bins"] == 1
    assert entity["decision"]["scope"] != "2d_connectivity_candidate"


def test_v16_result_path_keeps_old_v15_artifact(tmp_path):
    from pipeline.v15_flagship import v15_flagship_artifact_path
    old = tmp_path / "example_v15_flagship.json"
    old.write_text("existing result, do not replace")
    current = v15_flagship_artifact_path(tmp_path, "example")
    assert current.name.startswith("example_derived_")
    assert current.name.endswith("_v16_flagship.json")
    assert current != old
    assert old.read_text() == "existing result, do not replace"


def test_observed_audit_does_not_promote_a_regressing_primary_ranker():
    from pipeline.v16_candidate_missingness import attach_replicate_missingness_audit
    primary_policy = load_v15_flagship_policy()
    assert primary_policy["candidate_pool"]["same_energy_replicate_support"] == "legacy_zero_fill"
    a, b = _candidate("a", IK_A, .90), _candidate("b", IK_B, .60)
    rows = [_row(10, [a, b]), _row(10, [b]), _row(10, [b])]
    before = build_flagship_result("audit", rows, policy=primary_policy)
    frozen = copy.deepcopy(before)
    after = attach_replicate_missingness_audit(before, rows)
    assert before == frozen
    old_entity, new_entity = before["entities"][0], after["entities"][0]
    audit = new_entity.pop("replicate_missingness_audit")
    assert new_entity == old_entity
    assert audit["counterfactual_winner_differs"] is True
    assert audit["primary_ranking_changed"] is False
    assert audit["observed_only_candidates"][0]["observed_only_entity_score"] > .9


def test_missingness_audit_rejects_rank_promotion_or_changed_source():
    from pipeline.v16_candidate_missingness import attach_replicate_missingness_audit
    a, b = _candidate("a", IK_A, .90), _candidate("b", IK_B, .60)
    rows = [_row(10, [a, b]), _row(10, [b])]
    before = build_flagship_result("audit", rows, policy=load_v15_flagship_policy())
    with pytest.raises(ValueError, match="different source"):
        attach_replicate_missingness_audit(before, rows[:1])
    before["policy"]["replicate_missingness_audit"]["affects_primary_ranking"] = True
    with pytest.raises(ValueError, match="cannot modify"):
        attach_replicate_missingness_audit(before, rows)


def test_no_replicates_does_not_run_extra_ranking(monkeypatch):
    import pipeline.v16_candidate_missingness as audit
    a = _candidate("a", IK_A, .90)
    rows = [_row(10, [a]), _row(20, [a])]
    before = build_flagship_result("single", rows, policy=load_v15_flagship_policy())
    def forbidden(*args, **kwargs):
        raise AssertionError("No repeated-query correction is needed")
    monkeypatch.setattr(audit, "build_flagship_result", forbidden)
    after = audit.attach_replicate_missingness_audit(before, rows)
    assert after["entities"] == before["entities"]
    assert after["metadata"]["entities_with_incomplete_replicate_candidates"] == 0
