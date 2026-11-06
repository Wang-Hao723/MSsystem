#!/usr/bin/env python3
"""Verify committed V15 feature-fragment benchmark artifacts and invariants."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Iterable

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from matching.feature_fragment_benchmark.corpus import canonical_bytes, file_sha256


DB_IDENTIFIER_NAME = re.compile(
    r"^(?:NCGC[-_A-Z0-9]*|CCMSLIB[A-Z0-9_-]*|HMDB\d+|CHEBI:\d+|"
    r"CID:\d+|DTXSID\d+|MSBNK[-_A-Z0-9]+|LM[A-Z0-9]+)$",
    re.IGNORECASE,
)
SPLITS = ("discovery", "calibration", "test")


def _read_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def _iter_gzip_jsonl(path: Path) -> Iterable[dict[str, Any]]:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            value = json.loads(line)
            if not isinstance(value, dict):
                raise ValueError(f"{path}:{line_number}: JSONL value is not an object")
            yield value


def _decompressed_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with gzip.open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _split_seals(rows: list[dict[str, Any]]) -> tuple[str, str]:
    seal_material = [
        [row["spectrum_id"], row["split"], row["component_id"]]
        for row in sorted(rows, key=lambda item: item["spectrum_id"])
    ]
    test_material = sorted(
        row["spectrum_id"] for row in rows if row["split"] == "test"
    )
    return (
        hashlib.sha256(canonical_bytes(seal_material)).hexdigest(),
        hashlib.sha256(canonical_bytes(test_material)).hexdigest(),
    )


def audit(
    artifact_dir: Path,
    curation_receipt_path: Path,
    source_manifest_path: Path,
    instrument_evidence_path: Path,
) -> dict[str, Any]:
    corpus = _read_json(artifact_dir / "corpus_receipt.json")
    benchmark = _read_json(artifact_dir / "benchmark_receipt.json")
    export_gate = _read_json(artifact_dir / "runtime_export_gate.json")
    tested_rules = _read_json(artifact_dir / "sealed_test_rules.json")
    curation = _read_json(curation_receipt_path)
    instrument_evidence = _read_json(instrument_evidence_path)

    artifact_hash_checks = {}
    for filename, expected in benchmark["artifact_sha256"].items():
        artifact_hash_checks[filename] = {
            "expected": expected,
            "observed": file_sha256(artifact_dir / filename),
        }
        artifact_hash_checks[filename]["matches"] = (
            artifact_hash_checks[filename]["observed"] == expected
        )

    record_rows = list(_iter_gzip_jsonl(artifact_dir / "record_manifest.jsonl.gz"))
    split_rows = list(_iter_gzip_jsonl(artifact_dir / "split_manifest.jsonl.gz"))
    mirror_rows = list(_iter_gzip_jsonl(artifact_dir / "mirror_clusters.jsonl.gz"))
    reject_rows = list(_iter_gzip_jsonl(artifact_dir / "rejects.jsonl.gz"))

    record_ids = [row["spectrum_id"] for row in record_rows]
    split_ids = [row["spectrum_id"] for row in split_rows]
    split_counts = Counter(row["split"] for row in split_rows)
    identities_by_split: dict[str, set[str]] = defaultdict(set)
    formulas_by_split: dict[str, set[str]] = defaultdict(set)
    mirrors_by_split: dict[str, set[str]] = defaultdict(set)
    domains: set[str] = set()
    sources: set[str] = set()
    labs: set[str] = set()
    names_that_are_database_ids: list[str] = []
    cross_ce_conversion_count = 0
    qtof_beam_upgrade_count = 0
    invalid_input_stage_count = 0
    invalid_query_confirmation_count = 0
    for row in record_rows:
        split = row["split"]["split"]
        identities_by_split[split].add(row["compound"]["connectivity"])
        formulas_by_split[split].add(row["compound"]["formula"])
        mirrors_by_split[split].add(row["independence"]["mirror_cluster_id"])
        domains.add(row["acquisition"]["instrument_domain_id"])
        sources.add(row["source"]["source_id"])
        labs.add(row["source"]["laboratory_id"])
        if DB_IDENTIFIER_NAME.fullmatch(row["compound"]["display_name"]):
            names_that_are_database_ids.append(row["spectrum_id"])
        cross_ce_conversion_count += bool(
            row["acquisition"]["collision_energy"].get("conversion_applied")
        )
        if row["acquisition"]["analyzer"] == "QTOF":
            qtof_beam_upgrade_count += bool(
                row["acquisition"].get("psi_beam_type_cid_confirmed")
            )
        invalid_input_stage_count += row.get("evidence_stage") != "raw_msms_spectrum"
        invalid_query_confirmation_count += (
            row.get("query_confirmation_status") != "not_applicable"
        )

    overlap_counts: dict[str, dict[str, int]] = {}
    for label, values in (
        ("identity", identities_by_split),
        ("formula", formulas_by_split),
        ("mirror", mirrors_by_split),
    ):
        overlap_counts[label] = {}
        for index, left in enumerate(SPLITS):
            for right in SPLITS[index + 1 :]:
                overlap_counts[label][f"{left}__{right}"] = len(
                    values[left] & values[right]
                )
    forbidden_overlap_total = sum(
        count for dimension in overlap_counts.values() for count in dimension.values()
    )
    split_manifest_sha256, sealed_test_ids_sha256 = _split_seals(split_rows)

    mirror_duplicate_edges = sum(row["member_count"] - 1 for row in mirror_rows)
    tested_grade_counts = Counter(rule["evidence_grade"] for rule in tested_rules)
    tested_domain_counts = Counter(
        rule["applicability"]["instrument_domain_id"] for rule in tested_rules
    )
    tested_class_counts = Counter(rule["chemical_class"]["name"] for rule in tested_rules)
    upgrade_violation_count = sum(
        bool(rule.get(field))
        for rule in tested_rules
        for field in (
            "identity_upgrade_allowed",
            "formula_upgrade_allowed",
            "isomer_resolution_allowed",
            "standard_confirmation_allowed",
        )
    )

    curation_hash_checks = {}
    for filename, expected in curation["asset_sha256"].items():
        path = curation_receipt_path.parent / filename
        curation_hash_checks[filename] = {
            "expected": expected,
            "observed": file_sha256(path),
        }
        curation_hash_checks[filename]["matches"] = (
            curation_hash_checks[filename]["observed"] == expected
        )

    checks = {
        "benchmark_artifact_hashes_match": all(
            row["matches"] for row in artifact_hash_checks.values()
        ),
        "curation_asset_hashes_match": all(
            row["matches"] for row in curation_hash_checks.values()
        ),
        "source_manifest_hash_matches": (
            file_sha256(source_manifest_path) == corpus["manifest_sha256"]
        ),
        "record_count_matches": len(record_rows) == corpus["accepted_record_count"],
        "reject_count_matches": len(reject_rows) == corpus["rejected_record_count"],
        "raw_count_balances": (
            len(record_rows) + len(reject_rows) == corpus["raw_record_count"]
        ),
        "record_ids_unique": len(record_ids) == len(set(record_ids)),
        "split_ids_unique": len(split_ids) == len(set(split_ids)),
        "record_and_split_ids_match": set(record_ids) == set(split_ids),
        "split_counts_match": dict(sorted(split_counts.items()))
        == benchmark["split"]["split_record_counts"],
        "split_manifest_seal_matches": split_manifest_sha256
        == benchmark["split"]["split_manifest_sha256"],
        "sealed_test_id_seal_matches": sealed_test_ids_sha256
        == benchmark["sealed_test_guard"]["sealed_test_spectrum_ids_sha256"],
        "identity_formula_mirror_overlap_is_zero": forbidden_overlap_total == 0,
        "source_count_matches": len(sources) == benchmark["input"]["source_count"],
        "laboratory_count_matches": len(labs)
        == benchmark["input"]["laboratory_count"],
        "instrument_domain_count_matches": len(domains)
        == benchmark["input"]["instrument_domain_count"],
        "exact_mirror_edge_count_matches": mirror_duplicate_edges
        == benchmark["mirror_audit"]["exact_duplicate_edge_count"],
        "database_identifiers_not_used_as_names": not names_that_are_database_ids,
        "input_stage_is_raw_msms_only": invalid_input_stage_count == 0,
        "query_spectra_not_claimed_as_confirmations": invalid_query_confirmation_count == 0,
        "cross_ce_conversion_count_is_zero": cross_ce_conversion_count == 0,
        "qtof_contributor_cid_not_upgraded_to_psi_beam_cid": qtof_beam_upgrade_count == 0,
        "tested_rule_count_matches": len(tested_rules)
        == benchmark["rule_counts"]["sealed_test_evaluated"],
        "tested_rule_grades_match": dict(sorted(tested_grade_counts.items()))
        == benchmark["rule_counts"]["sealed_test_evidence_grades"],
        "tested_rules_do_not_upgrade_identity_formula_isomer_or_standard": (
            upgrade_violation_count == 0
        ),
        "runtime_export_is_empty": export_gate["exported_rule_count"] == 0,
        "legacy_ledger_has_terminal_decision_for_every_item": (
            curation["coverage"]["ledger_item_count"]
            == curation["coverage"]["terminal_decision_count"]
            == curation["coverage"]["unique_terminal_item_count"]
            == 2030
        ),
        "legacy_curation_runtime_export_is_empty": (
            curation["proposals"].get("runtime_export_eligible_rule_ids") == []
        ),
        "instrument_evidence_covers_five_families": len(
            instrument_evidence["families"]
        )
        >= 5,
        "instrument_evidence_forbids_cross_domain_ce_equivalence": (
            instrument_evidence["record_level_policy"][
                "cross_vendor_equivalence_allowed"
            ]
            is False
            and instrument_evidence["record_level_policy"][
                "cross_model_equivalence_allowed"
            ]
            is False
            and instrument_evidence["record_level_policy"][
                "ce_unit_conversion_allowed"
            ]
            is False
        ),
    }
    failures = sorted(name for name, passed in checks.items() if not passed)
    return {
        "schema_version": "2.0.0",
        "receipt_id": "v15-feature-fragment-benchmark-v2-20260826",
        "generated_on": "2026-08-26",
        "status": "pass" if not failures else "fail",
        "scientific_scope": "public_repository_experimental_msms_spectra",
        "runtime_effect": "none_shadow_benchmark_only",
        "checks": checks,
        "failures": failures,
        "corpus": {
            "repository_release": corpus["repository_release"],
            "raw_record_count": corpus["raw_record_count"],
            "accepted_record_count": corpus["accepted_record_count"],
            "rejected_record_count": corpus["rejected_record_count"],
            "independent_2d_identity_count": corpus[
                "independent_2d_identity_count"
            ],
            "formula_family_count": corpus["formula_family_count"],
            "source_count": len(sources),
            "laboratory_count": len(labs),
            "instrument_domain_count": len(domains),
            "corpus_jsonl_sha256": corpus["corpus_jsonl_sha256"],
            "raw_peak_corpus_committed": False,
            "record_manifest_sha256": file_sha256(
                artifact_dir / "record_manifest.jsonl.gz"
            ),
            "rejects_uncompressed_sha256": _decompressed_sha256(
                artifact_dir / "rejects.jsonl.gz"
            ),
        },
        "split": {
            "record_counts": dict(sorted(split_counts.items())),
            "identity_observation_counts": benchmark[
                "identity_observation_counts"
            ],
            "forbidden_overlap_counts_recomputed": overlap_counts,
            "forbidden_overlap_total_recomputed": forbidden_overlap_total,
            "source_holdout_enforced": benchmark["split"]["source_isolation"][
                "enforced"
            ],
            "source_overlap_counts_reported": benchmark["split"][
                "isolation_audit"
            ]["source_overlap_counts_reported_not_forbidden"],
        },
        "rules": {
            "discovered": benchmark["rule_counts"]["discovered"],
            "calibrated": benchmark["rule_counts"]["calibrated_total"],
            "selected_for_sealed_test": benchmark["rule_counts"][
                "selected_for_sealed_test"
            ],
            "sealed_test_evaluated": len(tested_rules),
            "sealed_test_grade_counts": dict(sorted(tested_grade_counts.items())),
            "sealed_test_domain_counts": dict(sorted(tested_domain_counts.items())),
            "sealed_test_class_counts": dict(sorted(tested_class_counts.items())),
            "runtime_exported": export_gate["exported_rule_count"],
            "runtime_rejected": export_gate["rejected_rule_count"],
            "identity_formula_isomer_standard_upgrade_violation_count": (
                upgrade_violation_count
            ),
        },
        "collision_energy": {
            "series_count": benchmark["ce_dependence"]["series_count"],
            "status_counts": benchmark["ce_dependence"]["status_counts"],
            "minimum_shared_test_identities": benchmark["ce_dependence"][
                "minimum_shared_test_identities"
            ],
            "cross_unit_conversion_applied": benchmark["ce_dependence"][
                "cross_unit_conversion_applied"
            ],
        },
        "legacy_curation": {
            "ledger_sha256": curation["ledger"]["sha256"],
            "ledger_item_count": curation["coverage"]["ledger_item_count"],
            "terminal_relationship_counts": curation["coverage"][
                "terminal_relationship_counts"
            ],
            "proposal_status_counts": curation["proposals"]["status_counts"],
            "runtime_export_eligible_rule_ids": curation["proposals"][
                "runtime_export_eligible_rule_ids"
            ],
            "quarantined_formula_mass_count": curation["formula_mass_quarantine"][
                "actual_count"
            ],
        },
        "artifact_hash_checks": artifact_hash_checks,
        "curation_hash_checks": curation_hash_checks,
        "instrument_evidence": {
            "path": str(instrument_evidence_path.relative_to(ROOT)),
            "sha256": file_sha256(instrument_evidence_path),
            "source_count": len(instrument_evidence["sources"]),
            "family_count": len(instrument_evidence["families"]),
        },
        "limitations": benchmark["limitations"]
        + [
            "the committed record manifest omits raw peaks; the full normalized peak corpus is identified by its SHA-256 and retained outside git",
            "the current corpus uses the recorded MassBank 2025.10 API release and is not relabelled as a later repository release",
            "all sealed-test rules remain provisional and none may enter runtime scoring"
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--artifact-dir",
        type=Path,
        default=ROOT / "docs/audits/feature_fragment_benchmark_v2",
    )
    parser.add_argument(
        "--curation-receipt",
        type=Path,
        default=ROOT
        / "config/feature_fragments/curation/acceptance_receipt_v1.json",
    )
    parser.add_argument(
        "--source-manifest",
        type=Path,
        default=ROOT
        / "config/feature_fragment_benchmark_v2/public_massbank_sources.json",
    )
    parser.add_argument(
        "--instrument-evidence",
        type=Path,
        default=ROOT
        / "config/feature_fragment_benchmark_v2/instrument_evidence.json",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT
        / "docs/audits/V15_FEATURE_FRAGMENT_BENCHMARK_V2_RECEIPT_20260826.json",
    )
    args = parser.parse_args()
    receipt = audit(
        args.artifact_dir,
        args.curation_receipt,
        args.source_manifest,
        args.instrument_evidence,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(receipt, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(receipt, ensure_ascii=False, indent=2, sort_keys=True))
    return 0 if receipt["status"] == "pass" else 1


if __name__ == "__main__":
    raise SystemExit(main())
