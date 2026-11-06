import json
from pathlib import Path

import pytest

from scripts.audit_feature_fragment_atomicization import (
    EXPECTED_ACCEPTED_SHADOW_RULE_IDS,
    audit_curation,
)
from scripts.curate_feature_fragment_legacy import (
    DEFAULT_LEDGER,
    DEFAULT_OUTPUT_DIR,
    EXPECTED_LEDGER_RECORDS,
    EXPECTED_LEDGER_SHA256,
    FORMULA_MASS_QUARANTINE_IDS,
    TERMINAL_RELATIONSHIPS,
    _alias_graph_is_acyclic,
    build_assets,
    formula_mz,
    load_ledger,
    sha256_bytes,
)


def _jsonl(name: str) -> list[dict]:
    path = DEFAULT_OUTPUT_DIR / name
    return [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines()]


def test_frozen_ledger_sha_and_2030_row_coverage():
    ledger_bytes = DEFAULT_LEDGER.read_bytes()
    ledger = load_ledger(DEFAULT_LEDGER)
    aliases = _jsonl("legacy_aliases_v1.jsonl")

    assert sha256_bytes(ledger_bytes) == EXPECTED_LEDGER_SHA256
    assert len(ledger) == EXPECTED_LEDGER_RECORDS == 2030
    assert len(aliases) == 2030
    assert {item["legacy_item_id"] for item in aliases} == {
        item["record"]["item_id"] for item in ledger
    }


def test_every_item_has_one_mutually_exclusive_terminal_relationship():
    aliases = _jsonl("legacy_aliases_v1.jsonl")
    item_ids = [item["legacy_item_id"] for item in aliases]

    assert len(item_ids) == len(set(item_ids)) == 2030
    assert all(isinstance(item["relationship"], str) for item in aliases)
    assert all(item["relationship"] in TERMINAL_RELATIONSHIPS for item in aliases)
    assert all("terminal_relationships" not in item for item in aliases)
    assert all(item["runtime_influence"] is False for item in aliases)

    receipt = json.loads(
        (DEFAULT_OUTPUT_DIR / "acceptance_receipt_v1.json").read_text(encoding="utf-8")
    )
    assert receipt["coverage"]["terminal_decision_count"] == 2030
    assert receipt["coverage"]["unique_terminal_item_count"] == 2030
    assert receipt["coverage"]["orphan_ledger_item_count"] == 0
    assert receipt["coverage"]["unknown_overlay_item_count"] == 0
    assert receipt["coverage"]["terminal_relationships_mutually_exclusive"] is True


def test_alias_and_supersession_graph_is_a_dag():
    aliases = _jsonl("legacy_aliases_v1.jsonl")
    acyclic, edge_count = _alias_graph_is_acyclic(aliases)

    assert acyclic is True
    assert edge_count == 15


def test_alias_dag_checker_rejects_a_cycle():
    aliases = [dict(item) for item in _jsonl("legacy_aliases_v1.jsonl")]
    aliases[0]["alias_of_item_id"] = aliases[1]["legacy_item_id"]
    aliases[1]["alias_of_item_id"] = aliases[0]["legacy_item_id"]

    with pytest.raises(ValueError, match="alias DAG cycle"):
        _alias_graph_is_acyclic(aliases)


def test_nine_formula_mass_inconsistent_profiles_are_quarantined():
    quarantine = _jsonl("quarantine_v1.jsonl")
    aliases = {
        item["legacy_item_id"]: item for item in _jsonl("legacy_aliases_v1.jsonl")
    }
    formula_rows = [
        item for item in quarantine if item["category"] == "formula_mass_inconsistent"
    ]

    assert len(formula_rows) == 9
    assert {item["legacy_item_id"] for item in formula_rows} == set(
        FORMULA_MASS_QUARANTINE_IDS
    )
    for item in formula_rows:
        assert aliases[item["legacy_item_id"]]["relationship"] == "quarantined_formula_mass"
        assert abs(item["mass_error_da"]) > 0.01
        assert formula_mz(item["formula"]) == pytest.approx(
            item["calculated_mz"], abs=1e-6
        )
        assert item["correction_policy"].startswith("do_not_guess")
        assert item["runtime_influence"] is False


def test_six_public_pilot_rules_are_ci_only_and_not_runtime_export_eligible():
    proposals = _jsonl("atomicization_queue_v1.jsonl")
    accepted = [item for item in proposals if item["status"] == "accepted_shadow"]

    assert len(proposals) == 40
    assert len(accepted) == 6
    assert {item["existing_runtime_rule_id"] for item in accepted} == (
        EXPECTED_ACCEPTED_SHADOW_RULE_IDS
    )
    assert all(item["export_eligible"] is False for item in proposals)
    for proposal in proposals:
        assert proposal["runtime_influence"] is False
        assert proposal["identity_rank_influence"] is False
        assert proposal["research_audit_only"] is True
        assert proposal["empirical_validation_gate_passed"] is False
        assert abs(proposal["observation"]["formula_mass_error_da"]) <= 0.01
        assert proposal["output_semantics"] == {
            "evidence_stage": "msms_rule_support",
            "exact_identity_confirmed": False,
            "standard_confirmed": False,
            "isomer_discriminating": False,
        }
    for proposal in accepted:
        assert proposal["benchmark_role"] == "deterministic_ci_fixture_only"
        assert proposal["acquisition_domain"]["activation"] == "CID"
        assert (
            proposal["acquisition_domain"]["psi_beam_type_cid_confirmed"]
            is False
        )
        assert "contributor_text" in proposal["acquisition_domain"][
            "instrument_domain_ids"
        ][0]


def test_checked_in_assets_are_deterministic_generator_output():
    expected_assets = build_assets(DEFAULT_LEDGER)

    assert set(expected_assets) == {
        "ledger_lock_v1.json",
        "legacy_aliases_v1.jsonl",
        "atomicization_queue_v1.jsonl",
        "quarantine_v1.jsonl",
        "source_claims_v1.jsonl",
        "statistics_plans_v1.json",
        "acceptance_receipt_v1.json",
    }
    for name, expected in expected_assets.items():
        assert (DEFAULT_OUTPUT_DIR / name).read_bytes() == expected


def test_full_atomicization_audit_and_runtime_isolation_pass():
    result = audit_curation()

    assert result["passed"] is True
    assert result["ledger"]["record_count"] == 2030
    assert result["ledger"]["source_file_count"] == 5
    assert result["coverage"]["terminal_states_mutually_exclusive"] is True
    assert result["alias_graph"]["acyclic"] is True
    assert result["formula_mass_quarantine"]["count"] == 9
    assert result["runtime_isolation"] == {
        "loader_path": "matching/feature_fragments/catalog.py",
        "explicit_runtime_files": [
            "instrument_domains.json",
            "rules_v1.json",
            "source_catalog.json",
        ],
        "curation_reference_present": False,
        "directory_discovery_present": False,
    }
