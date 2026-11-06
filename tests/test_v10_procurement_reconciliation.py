import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DATA = (
    ROOT
    / "outputs"
    / "019fc5bf-15dc-76a1-9d64-57597ecfcda1"
    / "inventory_reconcile_corrected"
    / "八个方_库存采购核对.json"
)


def load_payload():
    return json.loads(DATA.read_text(encoding="utf-8"))


def test_formula_inventory_reconciliation_is_complete_and_nonduplicated():
    payload = load_payload()
    summary = payload["summary"]
    rows = payload["procurement"]

    assert summary["frozen_standard_entries"] == 88
    assert summary["unique_formula_targets"] == 56
    assert summary["source_inventory"] == "/Users/zizi/Desktop/标品.xlsx"
    assert summary["inventory_source_rows"] == 91
    assert summary["inventory_count"] == 76
    assert summary["existing_formula_targets"] == 55
    assert summary["inventory_matched_formula_targets"] == 55
    assert summary["identity_conflict_targets"] == 1
    assert summary["unresolved_formula_targets"] == 1
    assert summary["missing_formula_targets"] == 0
    assert len(rows) == 56
    assert len({row["name_zh"] for row in rows}) == 56


def test_purchase_phases_separate_core_safety_and_c18_deferred_items():
    payload = load_payload()
    rows = payload["procurement"]
    by_name = {row["name_zh"]: row for row in rows}

    assert sum(row["priority"] in {"P0", "P1", "P2"} for row in rows) == 0
    assert by_name["石斛酚"]["priority"] == "待真值确认"
    assert sum(row["priority"] == "P1" for row in rows) == 0
    assert sum(row["priority"] == "P2" for row in rows) == 0
    assert by_name["马兜铃酸A"]["priority"] == "已有"
    assert by_name["马兜铃酸A"]["handling_class"] == "已有-安全管控"
    for name in {"L-亮氨酸", "精氨酸", "缬氨酸", "甘油三油酸酯"}:
        assert by_name[name]["priority"] == "已有"
        assert by_name[name]["handling_class"] == "已有-方法域受限"
    assert len(payload["hard_companions"]) == 6
    assert all(row["inventory_status"] == "缺少" for row in payload["hard_companions"])


def test_existing_but_failed_targets_are_not_reordered_as_new_purchases():
    payload = load_payload()
    cannot_solve = {row["name_zh"] for row in payload["purchase_cannot_solve"]}
    rows = {row["name_zh"]: row for row in payload["procurement"]}

    assert cannot_solve == {
        "甘草酸铵",
        "梓醇",
        "芥子碱硫氰酸盐",
        "蒿本内酯",
        "白桦脂酸",
        "马兜铃酸A",
        "苦杏仁苷",
    }
    assert rows["甘草酸铵"]["priority"] == "已有"
    assert rows["梓醇"]["priority"] == "已有"
    assert all(not str(row["identifier_review"]).startswith("unresolved") for row in rows.values())


def test_inventory_hits_use_exact_curated_names_not_fuzzy_substrings():
    payload = load_payload()
    existing = {row["name_zh"] for row in payload["procurement"] if row["priority"] == "已有"}

    assert len(existing) == 55
    assert existing == {row["name_zh"] for row in payload["procurement"]} - {"石斛酚"}
    catechin = next(row for row in payload["procurement"] if row["name_zh"] == "儿茶素")
    assert catechin["inventory_ids"] == "A1668"
    assert catechin["inventory_cas_recorded"] == "154-23-4"


def test_common_name_collision_does_not_create_a_false_inventory_hit():
    payload = load_payload()
    row = next(row for row in payload["procurement"] if row["name_zh"] == "石斛酚")

    assert row["inventory_status"] == "已有货号-真值结构冲突"
    assert row["inventory_ids"] == "A1558"
    assert row["name_en"] == "Moscatilin"
    assert row["cas_candidates"] == "108853-14-1"
    assert row["inchikey"] == "YTRAYUIKLRABOQ-UHFFFAOYSA-N"
    assert row["inventory_english_recorded"] == "Gigantol"
    assert row["inventory_cas_recorded"] == "83088-28-2"
    assert row["inventory_known_inchikeys"] == "SDXKZPQOVUDXIY-UHFFFAOYSA-N"
    assert row["order_decision"].startswith("确认前不采购")
    erratum = payload["truth_identity_errata"][0]
    assert erratum["standard_code"] == "A1558"
    assert erratum["legacy_structure"]["name"] == "Moscatilin"
    assert erratum["authoritative_inventory_record"]["name"] == "Gigantol"
    assert erratum["performance_handling"].startswith("exclude this entity")


def test_inventory_rows_are_deduplicated_and_excel_date_cas_values_are_recovered():
    payload = load_payload()
    inventory = {row["inventory_id"]: row for row in payload["inventory"]}

    assert len(inventory) == 76
    assert inventory["A0776"]["source_cells"] == "Sheet1!15；Sheet1!58；Sheet2!12"
    assert inventory["A0307"]["aliases_zh"] == "桂皮醛；肉桂醛"
    assert inventory["A0109"]["cas_recorded"] == "5373-11-5"
    assert inventory["A0564"]["cas_recorded"] == "2141-09-5"
