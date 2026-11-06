#!/usr/bin/env python3
"""Audit the V15 feature-fragment curation overlay without runtime loading.

The audit is deliberately independent of the production catalog loader.  It
reads only the frozen provenance ledger, the seven explicitly named curation
assets, five frozen source files, and the loader source used for the isolation
assertion.  It does not inspect data/, output/, or any database.
"""

from __future__ import annotations

import argparse
import ast
import json
import math
import sys
from collections import Counter
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.curate_feature_fragment_legacy import (
    DEFAULT_LEDGER,
    DEFAULT_OUTPUT_DIR,
    EXPECTED_LEDGER_RECORDS,
    EXPECTED_LEDGER_SHA256,
    FORMULA_MASS_QUARANTINE_IDS,
    PROPOSAL_STATUSES,
    SOURCE_FILE_HASHES,
    TERMINAL_RELATIONSHIPS,
    _alias_graph_is_acyclic,
    formula_mz,
    load_ledger,
    sha256_bytes,
    sha256_file,
)


REQUIRED_ASSETS = (
    "ledger_lock_v1.json",
    "legacy_aliases_v1.jsonl",
    "atomicization_queue_v1.jsonl",
    "quarantine_v1.jsonl",
    "source_claims_v1.jsonl",
    "statistics_plans_v1.json",
    "acceptance_receipt_v1.json",
)
EXPECTED_ACCEPTED_SHADOW_RULE_IDS = {
    "FF-PC-184-POS-H-AG6560-10EV",
    "FF-PC-184-POS-H-AG6560-20EV",
    "FF-PC-184-POS-H-AG6560-40EV",
    "FF-PE-NL141-POS-H-AG6560-10EV",
    "FF-PE-NL141-POS-H-AG6560-20EV",
    "FF-PE-NL141-POS-H-AG6560-40EV",
}


class AtomicizationAuditError(ValueError):
    """Raised when an overlay invariant fails."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise AtomicizationAuditError(message)


def _read_json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise AtomicizationAuditError(f"cannot read JSON asset {path}: {exc}") from exc
    _require(isinstance(value, dict), f"{path}: expected a JSON object")
    return value


def _read_jsonl(path: Path) -> list[dict[str, Any]]:
    try:
        raw_lines = path.read_bytes().splitlines()
    except OSError as exc:
        raise AtomicizationAuditError(f"cannot read JSONL asset {path}: {exc}") from exc
    records: list[dict[str, Any]] = []
    for line_number, raw_line in enumerate(raw_lines, start=1):
        try:
            value = json.loads(raw_line)
        except json.JSONDecodeError as exc:
            raise AtomicizationAuditError(
                f"{path}:{line_number}: invalid JSON: {exc}"
            ) from exc
        _require(isinstance(value, dict), f"{path}:{line_number}: expected object")
        records.append(value)
    return records


def _check_loader_isolation(loader_path: Path) -> dict[str, Any]:
    source = loader_path.read_text(encoding="utf-8")
    tree = ast.parse(source, filename=str(loader_path))
    string_values = {
        node.value
        for node in ast.walk(tree)
        if isinstance(node, ast.Constant) and isinstance(node.value, str)
    }
    expected_runtime_files = {
        "source_catalog.json",
        "instrument_domains.json",
        "rules_v1.json",
    }
    _require(
        expected_runtime_files.issubset(string_values),
        "runtime loader no longer names exactly the three expected catalog files",
    )
    _require("curation" not in source, "runtime loader references the curation overlay")
    forbidden_calls = {"glob", "rglob", "walk", "iterdir"}
    observed_calls = {
        node.func.attr
        for node in ast.walk(tree)
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute)
    }
    _require(
        not (forbidden_calls & observed_calls),
        "runtime loader performs directory discovery instead of explicit file loading",
    )
    return {
        "loader_path": str(loader_path.relative_to(ROOT)),
        "explicit_runtime_files": sorted(expected_runtime_files),
        "curation_reference_present": False,
        "directory_discovery_present": False,
    }


def audit_curation(
    *,
    ledger_path: Path = DEFAULT_LEDGER,
    curation_dir: Path = DEFAULT_OUTPUT_DIR,
    verify_source_files: bool = True,
) -> dict[str, Any]:
    for name in REQUIRED_ASSETS:
        _require((curation_dir / name).is_file(), f"missing curation asset: {name}")

    ledger_rows = load_ledger(ledger_path)
    ledger_ids = [entry["record"]["item_id"] for entry in ledger_rows]
    ledger_by_id = {entry["record"]["item_id"]: entry for entry in ledger_rows}
    _require(len(ledger_rows) == EXPECTED_LEDGER_RECORDS, "ledger row count drift")
    _require(len(set(ledger_ids)) == EXPECTED_LEDGER_RECORDS, "ledger item IDs not unique")

    lock = _read_json(curation_dir / "ledger_lock_v1.json")
    _require(lock.get("runtime_influence") is False, "ledger lock may not influence runtime")
    locked_ledger = lock.get("ledger") or {}
    _require(
        locked_ledger.get("sha256") == EXPECTED_LEDGER_SHA256,
        "ledger lock SHA-256 mismatch",
    )
    _require(
        locked_ledger.get("record_count") == EXPECTED_LEDGER_RECORDS,
        "ledger lock row count mismatch",
    )
    _require(lock.get("source_files") == SOURCE_FILE_HASHES, "source-file lock mismatch")
    source_hash_results: dict[str, str] = {}
    if verify_source_files:
        for relative, expected_hash in SOURCE_FILE_HASHES.items():
            path = ROOT / relative
            _require(path.is_file(), f"locked source file missing: {relative}")
            actual_hash = sha256_file(path)
            _require(
                actual_hash == expected_hash,
                f"locked source file drift: {relative} {actual_hash} != {expected_hash}",
            )
            source_hash_results[relative] = actual_hash

    aliases = _read_jsonl(curation_dir / "legacy_aliases_v1.jsonl")
    alias_ids = [item.get("legacy_item_id") for item in aliases]
    _require(len(aliases) == EXPECTED_LEDGER_RECORDS, "terminal overlay is not 2,030 rows")
    _require(len(set(alias_ids)) == len(alias_ids), "multiple terminal rows for one item")
    _require(set(alias_ids) == set(ledger_ids), "terminal overlay does not exactly cover ledger")
    for alias in aliases:
        item_id = alias["legacy_item_id"]
        relationship = alias.get("relationship")
        _require(
            isinstance(relationship, str) and relationship in TERMINAL_RELATIONSHIPS,
            f"{item_id}: invalid or non-scalar terminal relationship",
        )
        _require(
            "terminal_relationships" not in alias,
            f"{item_id}: multiple terminal relationship container is forbidden",
        )
        _require(alias.get("runtime_influence") is False, f"{item_id}: runtime influence")
        ledger_entry = ledger_by_id[item_id]
        _require(
            alias.get("ledger_line") == ledger_entry["line"],
            f"{item_id}: ledger line mismatch",
        )
        _require(
            alias.get("ledger_row_sha256") == ledger_entry["row_sha256"],
            f"{item_id}: ledger row hash mismatch",
        )
    acyclic, alias_edge_count = _alias_graph_is_acyclic(aliases)
    _require(acyclic, "alias graph is cyclic")

    proposals = _read_jsonl(curation_dir / "atomicization_queue_v1.jsonl")
    proposal_ids = [item.get("proposal_id") for item in proposals]
    _require(len(set(proposal_ids)) == len(proposal_ids), "duplicate proposal IDs")
    proposal_by_id = {item["proposal_id"]: item for item in proposals}
    alias_by_id = {item["legacy_item_id"]: item for item in aliases}
    source_claims = _read_jsonl(curation_dir / "source_claims_v1.jsonl")
    source_claim_ids = [item.get("source_claim_id") for item in source_claims]
    _require(len(set(source_claim_ids)) == len(source_claim_ids), "duplicate source claim IDs")
    source_claim_id_set = set(source_claim_ids)
    statistics_document = _read_json(curation_dir / "statistics_plans_v1.json")
    plans = statistics_document.get("plans") or []
    plan_ids = [item.get("statistics_plan_id") for item in plans]
    _require(len(set(plan_ids)) == len(plan_ids), "duplicate statistics plan IDs")
    plan_id_set = set(plan_ids)

    accepted_rule_ids: set[str] = set()
    for proposal in proposals:
        proposal_id = proposal["proposal_id"]
        status = proposal.get("status")
        _require(status in PROPOSAL_STATUSES, f"{proposal_id}: invalid status")
        _require(proposal.get("runtime_influence") is False, f"{proposal_id}: runtime influence")
        _require(
            proposal.get("identity_rank_influence") is False,
            f"{proposal_id}: identity rank influence",
        )
        _require(
            proposal.get("export_eligible") is False,
            f"{proposal_id}: curation proposal may not be runtime-export eligible",
        )
        _require(
            proposal.get("research_audit_only") is True,
            f"{proposal_id}: curation proposal is not restricted to research audit",
        )
        _require(
            proposal.get("empirical_validation_gate_passed") is False,
            f"{proposal_id}: curation proposal claims empirical validation",
        )
        _require(
            bool(proposal.get("empirical_validation_gate_unmet")),
            f"{proposal_id}: missing empirical validation failures",
        )
        observation = proposal.get("observation") or {}
        _require(
            abs(float(observation.get("formula_mass_error_da", math.inf))) <= 0.01,
            f"{proposal_id}: atomic proposal formula/mass mismatch",
        )
        _require(
            proposal.get("statistics_plan_id") in plan_id_set,
            f"{proposal_id}: unknown statistics plan",
        )
        _require(
            set(proposal.get("source_claim_ids") or []) <= source_claim_id_set,
            f"{proposal_id}: unknown source claim",
        )
        legacy_links = proposal.get("legacy_item_ids") or []
        _require(legacy_links, f"{proposal_id}: no legacy links")
        _require(set(legacy_links) <= set(ledger_ids), f"{proposal_id}: unknown legacy link")
        for legacy_item_id in legacy_links:
            _require(
                proposal_id in alias_by_id[legacy_item_id].get("canonical_proposal_ids", []),
                f"{proposal_id}: non-reciprocal link to {legacy_item_id}",
            )
        semantics = proposal.get("output_semantics") or {}
        for field in (
            "exact_identity_confirmed",
            "standard_confirmed",
            "isomer_discriminating",
        ):
            _require(semantics.get(field) is False, f"{proposal_id}: forbidden evidence upgrade")
        if status == "accepted_shadow":
            rule_id = proposal.get("existing_runtime_rule_id")
            _require(isinstance(rule_id, str) and rule_id, f"{proposal_id}: missing rule ID")
            accepted_rule_ids.add(rule_id)
            acquisition = proposal.get("acquisition_domain") or {}
            _require(
                acquisition.get("instrument_domain_ids")
                == ["agilent_6560_qtof_cid_contributor_text_fixed_ev"],
                f"{proposal_id}: accepted shadow has wrong domain",
            )
            _require(acquisition.get("activation") == "CID", f"{proposal_id}: activation")
            _require(
                acquisition.get("activation_evidence") == "MassBank contributor text",
                f"{proposal_id}: activation evidence",
            )
            _require(
                acquisition.get("psi_beam_type_cid_confirmed") is False,
                f"{proposal_id}: contributor CID was upgraded to PSI beam-type CID",
            )
            ce = acquisition.get("collision_energy") or {}
            _require(
                ce.get("unit") == "eV" and ce.get("mode") == "fixed",
                f"{proposal_id}: accepted CE is not exact fixed eV",
            )
            stats = proposal.get("existing_statistics") or {}
            for group in ("target", "background"):
                values = stats.get(group) or {}
                present = values.get("present")
                total = values.get("total")
                frequency = values.get("frequency")
                _require(
                    isinstance(present, int)
                    and isinstance(total, int)
                    and 0 <= present <= total
                    and total > 0,
                    f"{proposal_id}: invalid {group} counts",
                )
                _require(
                    abs(float(frequency) - present / total) <= 1e-12,
                    f"{proposal_id}: {group} frequency/count mismatch",
                )
                _require(
                    values.get("counting_unit") == "independent_2d_identity",
                    f"{proposal_id}: unsafe {group} counting unit",
                )
            _require(
                proposal.get("missing_schema_fields") == [],
                f"{proposal_id}: accepted rule still has missing fields",
            )
        else:
            _require(
                proposal.get("existing_runtime_rule_id") is None,
                f"{proposal_id}: nonaccepted proposal binds a runtime rule",
            )
            _require(
                bool(proposal.get("missing_schema_fields")),
                f"{proposal_id}: nonaccepted proposal hides missing validation",
            )
    _require(
        accepted_rule_ids == EXPECTED_ACCEPTED_SHADOW_RULE_IDS,
        "accepted-shadow set is not the six frozen public-pilot rules",
    )
    for alias in aliases:
        _require(
            set(alias.get("canonical_proposal_ids") or []) <= set(proposal_ids),
            f"{alias['legacy_item_id']}: unknown proposal link",
        )
        for proposal_id in alias.get("canonical_proposal_ids") or []:
            _require(
                alias["legacy_item_id"]
                in proposal_by_id[proposal_id].get("legacy_item_ids", []),
                f"{alias['legacy_item_id']}: non-reciprocal link to {proposal_id}",
            )

    quarantine = _read_jsonl(curation_dir / "quarantine_v1.jsonl")
    formula_quarantine = [
        item for item in quarantine if item.get("category") == "formula_mass_inconsistent"
    ]
    formula_ids = {item.get("legacy_item_id") for item in formula_quarantine}
    _require(len(formula_quarantine) == 9, "formula/mass quarantine count is not nine")
    _require(
        formula_ids == set(FORMULA_MASS_QUARANTINE_IDS),
        "formula/mass quarantine item set mismatch",
    )
    for item in formula_quarantine:
        item_id = item["legacy_item_id"]
        _require(
            alias_by_id[item_id]["relationship"] == "quarantined_formula_mass",
            f"{item_id}: terminal state is not formula quarantine",
        )
        recalculated = formula_mz(str(item["formula"]))
        _require(
            abs(recalculated - float(item["calculated_mz"])) <= 1e-6,
            f"{item_id}: calculated formula m/z drift",
        )
        _require(abs(float(item["mass_error_da"])) > 0.01, f"{item_id}: not out of tolerance")
        _require(item.get("runtime_influence") is False, f"{item_id}: quarantine influences runtime")

    receipt = _read_json(curation_dir / "acceptance_receipt_v1.json")
    coverage = receipt.get("coverage") or {}
    _require(coverage.get("ledger_item_count") == 2030, "receipt ledger count mismatch")
    _require(coverage.get("terminal_decision_count") == 2030, "receipt decision count mismatch")
    _require(coverage.get("unique_terminal_item_count") == 2030, "receipt unique count mismatch")
    _require(coverage.get("orphan_ledger_item_count") == 0, "receipt reports orphan items")
    _require(coverage.get("unknown_overlay_item_count") == 0, "receipt reports unknown items")
    actual_relationship_counts = dict(sorted(Counter(item["relationship"] for item in aliases).items()))
    _require(
        coverage.get("terminal_relationship_counts") == actual_relationship_counts,
        "receipt relationship counts are stale",
    )
    expected_hashes = receipt.get("asset_sha256") or {}
    _require(
        set(expected_hashes) == set(REQUIRED_ASSETS) - {"acceptance_receipt_v1.json"},
        "receipt asset hash manifest has the wrong members",
    )
    for name, expected_hash in expected_hashes.items():
        actual_hash = sha256_bytes((curation_dir / name).read_bytes())
        _require(actual_hash == expected_hash, f"receipt hash mismatch: {name}")
    runtime_boundary = receipt.get("runtime_boundary") or {}
    _require(
        runtime_boundary.get("curation_directory_loaded_by_runtime_catalog") is False,
        "receipt permits runtime curation loading",
    )
    _require(receipt.get("runtime_influence") is False, "receipt influences runtime")
    _require(
        (receipt.get("proposals") or {}).get("runtime_export_eligible_rule_ids") == [],
        "receipt exposes curation proposals for runtime export",
    )

    loader_isolation = _check_loader_isolation(
        ROOT / "matching" / "feature_fragments" / "catalog.py"
    )
    return {
        "passed": True,
        "ledger": {
            "record_count": len(ledger_rows),
            "sha256": EXPECTED_LEDGER_SHA256,
            "source_hashes_verified": verify_source_files,
            "source_file_count": len(source_hash_results),
        },
        "coverage": {
            "terminal_decisions": len(aliases),
            "unique_items": len(set(alias_ids)),
            "relationship_counts": actual_relationship_counts,
            "terminal_states_mutually_exclusive": True,
        },
        "alias_graph": {"acyclic": True, "edge_count": alias_edge_count},
        "formula_mass_quarantine": {
            "count": len(formula_quarantine),
            "all_outside_0_01_da": True,
        },
        "atomicization": {
            "proposal_count": len(proposals),
            "status_counts": dict(sorted(Counter(item["status"] for item in proposals).items())),
            "accepted_shadow_rule_ids": sorted(accepted_rule_ids),
            "runtime_export_eligible_rule_ids": [],
        },
        "runtime_isolation": loader_isolation,
        "receipt_hashes_verified": True,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ledger", type=Path, default=DEFAULT_LEDGER)
    parser.add_argument("--curation-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--skip-source-file-hashes", action="store_true")
    args = parser.parse_args()
    try:
        result = audit_curation(
            ledger_path=args.ledger,
            curation_dir=args.curation_dir,
            verify_source_files=not args.skip_source_file_hashes,
        )
    except (AtomicizationAuditError, ValueError, OSError) as exc:
        print(json.dumps({"passed": False, "error": str(exc)}, ensure_ascii=False, indent=2))
        return 1
    print(json.dumps(result, ensure_ascii=False, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
