from itertools import permutations

from matching.evidence_governance import (
    MIN_SCOPED_EXPERT_SCORE,
    POLICY_ID,
    annotate_evidence_passports,
    build_evidence_passport,
    build_feature_consensus,
    build_replicate_consensus,
    summarize_evidence_governance,
)


def _candidate(
    identity="AAAAAAAAAAAAAA-BBBBBBBBBB-C",
    name="Candidate A",
    *,
    cosine=0.93,
    final=0.72,
    expert=0.55,
    scoped=True,
    formula="C10H12O2",
):
    return {
        "compound_name": name,
        "inchikey": identity,
        "smiles": "CC1=CC=CC=C1O",
        "formula": formula,
        "adduct": "[M+H]+",
        "cosine": cosine,
        "entropy": 0.88,
        "reverse_match": 0.91,
        "final_score": final,
        "rank_score": final,
        "expert_score": expert,
        "expert_score_source": "v3",
        "expert_details": {
            "dimension_assessments": {
                "neutral_loss": {
                    "status": "scored" if scoped else "not_applicable",
                    "reason": "matched_curated_neutral_loss" if scoped else "class_unconfirmed",
                },
                "diagnostic_ion": {
                    "status": "not_applicable",
                    "reason": "class_unconfirmed",
                },
            },
            "score_fusion": {"selected_source": "v3"},
        },
        "mass_ppm": 2.1,
        "mass_metadata_consistent": True,
        "source_tier": "experimental",
        "dataset_source": "unit-test-library",
        "intensity_origin": "experimental",
        "quality_flags": [],
        "has_structure": True,
        "confidence": "High",
        "confidence_class": "High",
        "isomer_ambiguous": False,
        "top_candidates": [],
    }


def _row(
    scan,
    candidate=None,
    *,
    feature="PF000001",
    rt=5.0,
    mz=165.091,
    energy=20.0,
):
    return {
        "spectrum_id": str(scan),
        "precursor_mz": mz,
        "retention_time": rt,
        "ion_mode": "POS",
        "collision_energy": energy,
        "precursor_feature_id": feature,
        "precursor_feature_spectra": 3,
        "precursor_feature_collision_energies": [20.0, 40.0, 60.0],
        "precursor_feature_identity_consensus": 1.0,
        "precursor_feature_identity_consistent": True,
        "blank_status": "not_requested",
        "blank_used": False,
        "blank_flagged": False,
        "match": candidate or _candidate(),
    }


def _node(passport, key):
    return next(item for item in passport["nodes"] if item["key"] == key)


def test_passport_separates_applicability_independence_and_claim_scope():
    passport = build_evidence_passport(_row("1"))

    assert passport["policy_id"] == POLICY_ID
    assert passport["ranking_changed"] is False
    assert passport["confidence_field_overwritten"] is False
    assert passport["release_eligible_before_replication"] is True
    assert passport["claim_status"] == "provisional_needs_independent_injection"
    assert passport["coverage_guarantee"] is False

    spectral = _node(passport, "spectral_similarity")
    expert = _node(passport, "scoped_expert_reasoning")
    within = _node(passport, "within_injection_stability")
    blank = _node(passport, "blank_contrast")
    assert spectral["independent_measurement"] is False
    assert expert["independent_measurement"] is False
    assert "ranking" in expert["used_for"]
    assert "expert_score_already_contributes_to_final_score_do_not_count_as_independent" in expert["warnings"]
    assert within["independent_measurement"] is False
    assert "pseudoreplicate_guard_applied" in within["warnings"]
    assert blank["status"] == "not_run"
    assert blank["direction"] == "neutral"


def test_not_applicable_expert_rules_do_not_create_pseudo_evidence():
    row = _row("1", _candidate(scoped=False, expert=0.99))
    passport = build_evidence_passport(row)

    expert = _node(passport, "scoped_expert_reasoning")
    assert expert["status"] == "not_applicable"
    assert passport["counterfactuals"]["expert_unscoped_reports"] is True
    assert passport["counterfactuals"]["expert_scoped_reports"] is False
    assert "expert_not_applicable" in passport["release_blockers"]


def test_blank_role_and_same_formula_competition_limit_the_claim():
    candidate = _candidate()
    candidate["isomer_ambiguous"] = True
    candidate["isomer_group_size"] = 2
    passport = build_evidence_passport(_row("1", candidate))
    assert passport["claim_status"] == "isomer_group_only"
    assert passport["claim_scope"] == "molecular_formula_isomer_group"

    blank_passport = build_evidence_passport(
        _row("1", _candidate()), sample_role="blank"
    )
    assert blank_passport["claim_status"] == "blank_background_context_only"
    assert "declared_blank_query" in blank_passport["release_blockers"]
    assert _node(blank_passport, "blank_contrast")["status"] == "role_veto"


def test_unspecified_soft_blank_flag_is_context_not_a_hard_veto():
    row = _row("1")
    row.update({
        "blank_used": True,
        "blank_flagged": True,
        "blank_clean_view_excluded": False,
        "blank_status": "blank_associated",
        "blank_evidence_role": "unspecified",
    })
    passport = build_evidence_passport(row)

    assert passport["release_eligible_before_replication"] is True
    assert "blank_background_conflict" not in passport["release_blockers"]
    assert _node(passport, "blank_contrast")["status"] == "context_warning"

    row["blank_clean_view_excluded"] = True
    vetoed = build_evidence_passport(row)
    assert vetoed["release_eligible_before_replication"] is False
    assert "blank_background_conflict" in vetoed["release_blockers"]


def test_feature_consensus_collapses_stepped_energy_before_voting():
    candidate_a = _candidate()
    candidate_b = _candidate(
        "CCCCCCCCCCCCCC-DDDDDDDDDD-E", "Candidate B", final=0.75
    )
    rows = []
    for scan, energy, top, alternate in (
        ("1", 20.0, candidate_a, candidate_b),
        ("2", 40.0, candidate_a, candidate_b),
        ("3", 60.0, candidate_b, candidate_a),
    ):
        top = dict(top)
        top["top_candidates"] = [dict(alternate)]
        rows.append(_row(scan, top, energy=energy))

    features = build_feature_consensus(rows, "task-one")

    assert len(features) == 1
    feature = features[0]
    assert feature["n_spectra"] == 3
    assert feature["pseudoreplicates_collapsed"] == 2
    assert feature["consensus_lead"]["compound_name"] == "Candidate A"
    assert feature["consensus_lead"]["top1_votes"] == 2
    assert feature["shadow_consensus_only"] is True
    assert feature["ranking_changed"] is False


def _three_energy_task(rt=5.0, mz=165.091, candidate=None):
    candidate = candidate or _candidate()
    return [
        _row(str(index), dict(candidate), rt=rt, mz=mz, energy=energy)
        for index, energy in enumerate((20.0, 40.0, 60.0), start=1)
    ]


def test_replicate_consensus_requires_independent_groups_and_reports_ablation():
    results = {
        "inj-1": _three_energy_task(rt=5.000),
        "inj-2": _three_energy_task(rt=5.025, mz=165.0917),
    }
    consensus = build_replicate_consensus(results)
    cluster = consensus["clusters"][0]

    assert cluster["action"] == "release_replicate_stable_candidate"
    assert cluster["independent_injection_groups"] == 2
    assert cluster["qualified_injection_groups"] == 2
    assert consensus["summary"]["released_unique_identities"] == 1
    assert consensus["summary"]["pseudoreplicate_scans_collapsed"] == 4
    assert consensus["summary"]["counterfactual_unique_identities"] == {
        "spectral_only": 1,
        "expert_unscoped": 1,
        "expert_scoped": 1,
        "full_governed": 1,
    }
    assert consensus["summary"]["ranking_changed"] is False

    same_batch = build_replicate_consensus(
        results,
        independence_groups={"inj-1": "batch-a", "inj-2": "batch-a"},
    )
    assert same_batch["clusters"][0]["action"] == (
        "needs_more_independent_injections"
    )


def test_blank_recurrence_vetoes_analyte_release_without_deleting_the_match():
    results = {
        "inj-1": _three_energy_task(rt=5.000),
        "inj-2": _three_energy_task(rt=5.020, mz=165.0915),
        "blank-1": _three_energy_task(rt=5.010, mz=165.0912),
    }
    consensus = build_replicate_consensus(
        results, task_roles={"blank-1": "blank"}
    )
    cluster = consensus["clusters"][0]

    assert cluster["action"] == "abstain_blank_recurrent"
    assert cluster["blank_conflict"] is True
    assert cluster["ranking_changed"] is False
    assert consensus["summary"]["released_unique_identities"] == 0


def test_cohort_alignment_is_invariant_to_task_submission_order():
    # Dense neighbouring features previously exposed an order-dependent greedy
    # assignment.  The result must be identical for all six task permutations.
    task_positions = {
        "inj-a": (5.00, 5.05),
        "inj-b": (5.00, 5.05),
        "inj-c": (5.05, 5.10),
    }
    base = {
        task_id: [
            _row(
                f"{task_id}-{index}",
                feature=f"PF-{task_id}-{index}",
                rt=rt,
                mz=165.091,
            )
            for index, rt in enumerate(retention_times)
        ]
        for task_id, retention_times in task_positions.items()
    }

    signatures = set()
    for order in permutations(base):
        consensus = build_replicate_consensus(
            {task_id: base[task_id] for task_id in order},
            rt_window_min=0.10,
        )
        assert consensus["task_ids"] == sorted(base)
        signatures.add(tuple(
            sorted(
                tuple(sorted(
                    (member["task_id"], member["feature_id"])
                    for member in cluster["member_features"]
                ))
                for cluster in consensus["clusters"]
            )
        ))

    assert len(signatures) == 1


def test_annotation_summary_is_idempotent_and_identity_aware():
    rows = [_row("1"), _row("2")]
    annotate_evidence_passports(rows)
    first = summarize_evidence_governance(rows)
    annotate_evidence_passports(rows)
    second = summarize_evidence_governance(rows)

    assert first == second
    assert first["matched_spectra_audited"] == 2
    assert first["unique_identities_audited"] == 1
    assert first["release_eligible_before_replication_spectra"] == 2
    assert MIN_SCOPED_EXPERT_SCORE == 0.40
