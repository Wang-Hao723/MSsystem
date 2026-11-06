#!/usr/bin/env python3
"""Create a row-level provenance audit of existing V15 fragment knowledge.

Only explicit source files are read.  Project ``data/`` and ``output/``
directories are never discovered or scanned.  The resulting ledger separates
experimental statistics, papers/books, vendor documents, manual rules,
computational inference, and constants without an atomic source.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from collections import Counter
from datetime import date
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from matching import fragment_feature_library as feature_library
from matching import fragmentation_rules


DEFAULT_OUTPUT = (
    ROOT
    / "docs"
    / "audits"
    / "V15_FEATURE_FRAGMENT_LEGACY_PROVENANCE_AUDIT_20260826.json"
)
DEFAULT_LEDGER_OUTPUT = (
    ROOT
    / "docs"
    / "audits"
    / "V15_FEATURE_FRAGMENT_LEGACY_PROVENANCE_LEDGER_20260826.jsonl"
)


def _load(relative: str) -> dict[str, Any]:
    with (ROOT / relative).open("r", encoding="utf-8") as handle:
        return json.load(handle)


def _sha256(path: Path) -> str:
    hasher = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            hasher.update(block)
    return hasher.hexdigest()


def _reference_classes(references: list[str] | tuple[str, ...]) -> list[str]:
    cleaned = [str(value).strip() for value in references if str(value).strip()]
    return ["paper_or_book"] if cleaned else ["manual_rule", "unsourced_constant"]


def _base_record(
    *,
    item_id: str,
    system: str,
    path: str,
    kind: str,
    provenance_classes: list[str],
    references: list[str] | tuple[str, ...] = (),
    executable: bool | None = None,
    runtime_influence: str = "unknown",
    limitations: list[str] | None = None,
) -> dict[str, Any]:
    return {
        "item_id": item_id,
        "system": system,
        "path": path,
        "kind": kind,
        "provenance_classes": provenance_classes,
        "references": list(references),
        "executable": executable,
        "runtime_influence": runtime_influence,
        "instrument_conditioned": False,
        "activation_conditioned": False,
        "collision_energy_conditioned": False,
        "occurrence_frequency_available": False,
        "background_frequency_available": False,
        "independent_fragmentation_group_available": False,
        "negative_evidence_limit_available": False,
        "limitations": limitations or [],
    }


def _static_rule_records() -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    for name, loss in fragmentation_rules.COMMON_NEUTRAL_LOSSES.items():
        item = _base_record(
            item_id=f"static_neutral_loss:{name}",
            system="matching.fragmentation_rules.COMMON_NEUTRAL_LOSSES",
            path="matching/fragmentation_rules.py:128",
            kind="neutral_loss",
            provenance_classes=["manual_rule", "unsourced_constant"],
            executable=True,
            runtime_influence="expert_evidence_and_release_gate_indirectly",
            limitations=[
                "No atomic reference",
                "Priority is an unvalidated integer",
                "No target/background frequency or instrument/activation domain",
            ],
        )
        item.update(
            {
                "formula": loss.formula,
                "mass": loss.mass,
                "priority": loss.priority,
                "polarity": loss.polarity,
                "adducts": list(loss.allowed_adducts),
                "structure_constraints": list(loss.required_smarts),
                "collision_energy": {
                    "minimum": loss.collision_energy_min,
                    "maximum": loss.collision_energy_max,
                    "unit": loss.collision_energy_unit,
                },
            }
        )
        records.append(item)

    for name, rule in fragmentation_rules.FRAGMENTATION_RULES.items():
        references = list(rule.references)
        item = _base_record(
            item_id=f"static_fragmentation_rule:{name}",
            system="matching.fragmentation_rules.FRAGMENTATION_RULES",
            path="matching/fragmentation_rules.py:801",
            kind="fragmentation_rule",
            provenance_classes=_reference_classes(references),
            references=references,
            executable=True,
            runtime_influence="expert_evidence_and_release_gate_indirectly",
            limitations=[
                "References do not provide occurrence/background counts in the runtime domain",
                "Confidence is not calibrated on identity- or formula-isolated spectra",
            ],
        )
        item.update(
            {
                "rule_id": rule.rule_id,
                "mechanism_type": rule.mechanism_type,
                "confidence": rule.confidence_score,
                "polarity": rule.charge_preference,
                "adducts": list(rule.allowed_adducts),
                "structure_constraints": list(rule.required_smarts),
                "evidence_level": rule.evidence_level,
                "collision_energy": {
                    "minimum": rule.collision_energy_min,
                    "maximum": rule.collision_energy_max,
                    "unit": rule.collision_energy_unit,
                },
            }
        )
        item["collision_energy_conditioned"] = bool(rule.collision_energy_unit)
        records.append(item)

    for name, ion in fragmentation_rules.DIAGNOSTIC_IONS.items():
        references = list(ion.references)
        item = _base_record(
            item_id=f"static_diagnostic_ion:{name}",
            system="matching.fragmentation_rules.DIAGNOSTIC_IONS",
            path="matching/fragmentation_rules.py:1463",
            kind="diagnostic_ion",
            provenance_classes=_reference_classes(references),
            references=references,
            executable=True,
            runtime_influence="expert_evidence_and_release_gate_indirectly",
            limitations=[
                "No target/background frequency",
                "No independent fragmentation group",
                "Shared fragments cannot establish a unique class or isomer",
            ],
        )
        item.update(
            {
                "mz": ion.mz,
                "formula": ion.formula,
                "polarity": ion.polarity,
                "adducts": list(ion.allowed_adducts),
                "structure_constraints": list(ion.required_smarts),
                "evidence_level": ion.evidence_level,
                "collision_energy": {
                    "minimum": ion.collision_energy_min,
                    "maximum": ion.collision_energy_max,
                    "unit": ion.collision_energy_unit,
                },
            }
        )
        item["collision_energy_conditioned"] = bool(ion.collision_energy_unit)
        records.append(item)
    return records


def _feature_library_records() -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    for class_name, profile in feature_library.COMPOUND_CLASS_PROFILES.items():
        for ordinal, fragment in enumerate(profile.characteristic_fragments, start=1):
            references = list(fragment.references)
            item = _base_record(
                item_id=f"feature_profile:{class_name}:{ordinal:03d}",
                system="matching.fragment_feature_library.COMPOUND_CLASS_PROFILES",
                path="matching/fragment_feature_library.py:104",
                kind="class_feature_fragment",
                provenance_classes=_reference_classes(references),
                references=references,
                executable=True,
                runtime_influence="expert_evidence_and_fragment_only_classification",
                limitations=[
                    "Diagnostic value and profile confidence are not calibrated",
                    "The same peak can match repeated or near-duplicate entries",
                    "No occurrence/background frequency or independent fragmentation group",
                ],
            )
            item.update(
                {
                    "chemical_class": class_name,
                    "ordinal": ordinal,
                    "mz": fragment.mz,
                    "formula": fragment.formula,
                    "diagnostic_value": fragment.diagnostic_value,
                    "profile_confidence": profile.confidence_score,
                    "polarity": fragment.polarity,
                    "adduct_preference": fragment.adduct_preference,
                    "allowed_adducts": list(fragment.allowed_adducts),
                    "structure_constraints": list(fragment.required_smarts)
                    + list(fragment.required_smarts_all),
                    "evidence_level": fragment.evidence_level,
                    "collision_energy": {
                        "minimum": fragment.collision_energy_min,
                        "maximum": fragment.collision_energy_max,
                        "unit": fragment.collision_energy_unit,
                    },
                }
            )
            item["collision_energy_conditioned"] = bool(fragment.collision_energy_unit)
            records.append(item)
    return records


def _legacy_curation_records() -> list[dict[str, Any]]:
    payload = _load("matching/fragmentation_rules_curation.json")
    records: list[dict[str, Any]] = []
    for raw in payload["records"]:
        references = list(raw.get("reference_identifiers") or [])
        status = str(raw.get("status") or "")
        item = _base_record(
            item_id=f"legacy_rule:{raw['occurrence_id']}",
            system="matching.fragmentation_rules_curation",
            path="matching/fragmentation_rules_curation.json",
            kind="legacy_fragmentation_rule_occurrence",
            provenance_classes=_reference_classes(references),
            references=references,
            executable=status == "promoted_runtime",
            runtime_influence=(
                "promoted_runtime_only" if status == "promoted_runtime" else "none"
            ),
            limitations=[
                str(raw.get("status_rationale") or ""),
                "Legacy record lacks per-rule experimental target/background statistics",
            ],
        )
        item.update(
            {
                "source_rule_id": raw.get("source_rule_id"),
                "chemical_class": raw.get("class_en") or raw.get("class"),
                "status": status,
                "candidate_fragment_mz": raw.get("candidate_fragment_mz"),
                "candidate_polarity": raw.get("candidate_polarity"),
                "runtime_links": raw.get("runtime_links", []),
            }
        )
        records.append(item)
    return records


def _empirical_records() -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    old = _load("kb_v3/empirical_rules.json")
    for chemical_class, polarities in old.items():
        for polarity, group in polarities.items():
            for kind, value_key in (
                ("diagnostic_ions", "mz"),
                ("neutral_losses", "loss_da"),
            ):
                for ordinal, raw in enumerate(group.get(kind, []), start=1):
                    item = _base_record(
                        item_id=f"empirical_v1:{chemical_class}:{polarity}:{kind}:{ordinal:03d}",
                        system="kb_v3.empirical_rules",
                        path="kb_v3/empirical_rules.json",
                        kind=(
                            "diagnostic_ion"
                            if kind == "diagnostic_ions"
                            else "neutral_loss"
                        ),
                        provenance_classes=["real_experimental_spectrum_statistics"],
                        executable=True,
                        runtime_influence="expert_v3_fixed_weight_component",
                        limitations=[
                            "Counting unit is spectra, so technical repeats increase support",
                            "Instrument, activation, CE, and adduct conditions are pooled",
                            "No class-background frequency or identity/formula holdout",
                        ],
                    )
                    item.update(
                        {
                            "chemical_class": chemical_class,
                            "polarity": polarity,
                            value_key: raw.get(value_key),
                            "frequency": raw.get("freq"),
                            "supporting_spectra": raw.get("n_support"),
                            "total_spectra": group.get("n_spectra"),
                            "annotation": raw.get("annotation"),
                        }
                    )
                    item["occurrence_frequency_available"] = True
                    records.append(item)

    contrastive = _load("kb_v3/empirical_rules_candidates_v2.json")
    for chemical_class, polarities in contrastive["classes"].items():
        for polarity, group in polarities.items():
            for kind, value_key in (
                ("diagnostic_ions", "mz"),
                ("neutral_losses", "loss_da"),
            ):
                for ordinal, raw in enumerate(group.get(kind, []), start=1):
                    item = _base_record(
                        item_id=f"contrastive_v2:{chemical_class}:{polarity}:{kind}:{ordinal:03d}",
                        system="kb_v3.empirical_rules_candidates_v2",
                        path="kb_v3/empirical_rules_candidates_v2.json",
                        kind=(
                            "diagnostic_ion"
                            if kind == "diagnostic_ions"
                            else "neutral_loss"
                        ),
                        provenance_classes=[
                            "real_experimental_spectrum_statistics",
                            "computational_inference",
                        ],
                        executable=False,
                        runtime_influence="shadow_only",
                        limitations=[
                            "All instruments, activation methods, CE settings, and adduct spectra for an identity are unioned",
                            "Discovery counts and LR come from the same identities",
                            "Source-family labels do not guarantee mirror independence",
                        ],
                    )
                    item.update(
                        {
                            "chemical_class": chemical_class,
                            "polarity": polarity,
                            value_key: raw.get(value_key),
                            "occurrence_frequency": raw.get("prevalence"),
                            "background_frequency": raw.get("background_prevalence"),
                            "likelihood_ratio": raw.get("likelihood_ratio"),
                            "target_supporting_identities": raw.get(
                                "n_independent_compounds"
                            ),
                            "background_supporting_identities": raw.get(
                                "n_background_compounds"
                            ),
                            "source_family_count": raw.get("n_source_families"),
                            "source_support": raw.get("source_support", {}),
                            "evidence_level": raw.get("evidence_level"),
                        }
                    )
                    item["occurrence_frequency_available"] = True
                    item["background_frequency_available"] = True
                    records.append(item)
    return records


def _knowledge_systems() -> list[dict[str, Any]]:
    return [
        {
            "system": "expected_compounds.DiagnosticIons",
            "path": "matching/library_import.py:813; matching/expected_compounds.py:1963",
            "source_class": "manual_rule_or_literature_target",
            "row_level_inventory": "not_available_without_scanning_project_database",
            "runtime_effect": "Two ions and 25 percent coverage can upgrade a target to MS1_MS2_DIAGNOSTIC_PROBABLE",
            "audit_result": "P0: table rows have only target_id, fragment_mz, and ordinal; no per-ion source, polarity, adduct, instrument, activation, CE, frequency, background, or independence group",
        },
        {
            "system": "expected_compounds.experimental_library",
            "path": "matching/expected_compounds.py:1828",
            "source_class": "real_experimental_spectrum_match",
            "runtime_effect": "Fixed fragment-count and intensity-coverage gates label public library support as confirmed",
            "audit_result": "P0: a public MS/MS library candidate is not same-method standard confirmation",
        },
        {
            "system": "CEProfiles",
            "path": "scripts/build_ce_profiles.py; matching/ce_profiles.py",
            "source_class": "real_experimental_spectrum_statistics",
            "runtime_effect": "shadow",
            "audit_result": "Retains CE units and source records but omits activation and compares broad qtof/orbitrap families instead of exact vendor/model/cell domains",
        },
        {
            "system": "EvoGraph fragmentation",
            "path": "matching/evograph_fragmentation.py",
            "source_class": "computational_inference",
            "runtime_effect": "shadow",
            "audit_result": "Formula-conservation graph; not experimental fragment occurrence",
        },
        {
            "system": "CFM fragment prediction",
            "path": "utils/cfm_fragment.py",
            "source_class": "computational_inference",
            "runtime_effect": "report_only",
            "audit_result": "Predicted spectrum; cache key omits model hash, instrument, and CE",
        },
        {
            "system": "EviDiscrim",
            "path": "matching/evidiscrim.py; config/evidiscrim_policy.json",
            "source_class": "computational_inference",
            "runtime_effect": "shadow",
            "audit_result": "Candidate-pair counterfactual; default single-model HCD runner must abstain and cannot be counted as experimental support",
        },
    ]


def _fixed_constants() -> list[dict[str, Any]]:
    return [
        {
            "path": "matching/expert_scorer.py:174",
            "values": [0.22, 0.22, 0.18, 0.15, 0.10, 0.08, 0.03, 0.02],
            "role": "eight expert dimension weights",
            "provenance_class": "unsourced_constant",
        },
        {
            "path": "matching/expert_v3.py:19",
            "values": {
                "mass_tolerance_da": 0.05,
                "ion_weight": 0.6,
                "loss_weight": 0.4,
                "minimum_idf": 0.3,
            },
            "role": "v3 fragment matching and fusion",
            "provenance_class": "unsourced_constant",
        },
        {
            "path": "matching/fragment_feature_library.py:2101",
            "values": {
                "exclusivity_multiplier": 1.5,
                "minimum_score": 0.15,
                "ratio_threshold": 1.3,
            },
            "role": "fragment-only class assignment",
            "provenance_class": "unsourced_constant",
        },
        {
            "path": "matching/expected_compounds.py:47",
            "values": {
                "minimum_diagnostic_ions": 2,
                "minimum_coverage": 0.25,
                "library_fragment_pairs": 3,
                "library_intensity_coverage": 0.5,
            },
            "role": "evidence-stage upgrades",
            "provenance_class": "unsourced_constant",
        },
        {
            "path": "matching/evidence_governance.py:38",
            "values": {"expert_release_gate": 0.40},
            "role": "release gate",
            "provenance_class": "unsourced_constant",
        },
        {
            "path": "matching/expert_loss.py:560",
            "values": [0.25, 0.25, 0.15, 0.20, 0.15],
            "role": "training loss weights",
            "provenance_class": "unsourced_constant",
        },
    ]


def _known_defects() -> list[dict[str, Any]]:
    return [
        {
            "severity": "P0",
            "path": "matching/expert_scorer.py:1329",
            "defect": "Aconitine-like ester cascades can run for broad alkaloid/diterpene classes without an aconitine scaffold, instrument, activation, or CE guard.",
        },
        {
            "severity": "P0",
            "path": "matching/evidence_governance.py:139",
            "defect": "An uncalibrated expert score still gates release even though ranking contribution is zero.",
        },
        {
            "severity": "P0",
            "path": "matching/expected_compounds.py:1963",
            "defect": "Fixed diagnostic-ion thresholds upgrade an MS1 target without atomic provenance or acquisition conditions.",
        },
        {
            "severity": "P0",
            "path": "matching/expected_compounds.py:1828",
            "defect": "A public library match is labeled experimental-library confirmed although it is not same-method standard confirmation.",
        },
        {
            "severity": "P1",
            "path": "matching/fragmentation_rules.py:128",
            "defect": "Exact duplicate neutral losses HCOOH/H2O_CO and CH3COOH/C2H4O2_CGLY can double-represent one event.",
        },
        {
            "severity": "P1",
            "path": "matching/fragmentation_rules.py:1968",
            "defect": "Class signature characteristic_losses arrays contain values used elsewhere as product ions.",
        },
        {
            "severity": "P1",
            "path": "matching/fragment_feature_library.py:104",
            "defect": "Repeated and near-duplicate fragments allow one query peak to satisfy multiple entries.",
        },
        {
            "severity": "P1",
            "path": "matching/ce_profiles.py:361",
            "defect": "Instrument compatibility collapses all vendor QTOFs to qtof and all Orbitrap-family instruments to orbitrap.",
        },
        {
            "severity": "P2",
            "path": "matching/fragment_feature_library.py",
            "defect": "Nine formula/mass-inconsistent entries remain in source arrays even though a runtime gate rejects them.",
        },
        {
            "severity": "P2",
            "path": "kb_v3/build_kb_v3.py:301",
            "defect": "Common-loss annotations include formula/mass mismatches and an unresolved formula with a question mark.",
        },
    ]


def build_audit() -> dict[str, Any]:
    records = (
        _static_rule_records()
        + _feature_library_records()
        + _legacy_curation_records()
        + _empirical_records()
    )
    system_counts = Counter(record["system"] for record in records)
    kind_counts = Counter(record["kind"] for record in records)
    provenance_counts = Counter(
        provenance for record in records for provenance in record["provenance_classes"]
    )
    missing_atomic_sources = sum(
        "unsourced_constant" in record["provenance_classes"] for record in records
    )
    return {
        "schema_version": "1.0.0",
        "audit_id": "v15_feature_fragment_legacy_provenance_20260826",
        "generated_on": date.today().isoformat(),
        "scope": {
            "included": [
                "matching/fragmentation_rules.py",
                "matching/fragment_feature_library.py",
                "matching/fragmentation_rules_curation.json",
                "kb_v3/empirical_rules.json",
                "kb_v3/empirical_rules_candidates_v2.json",
                "explicit runtime and shadow call paths",
            ],
            "excluded": [
                "data/",
                "output/",
                "row-level project databases",
            ],
        },
        "summary": {
            "row_level_record_count": len(records),
            "counts_by_system": dict(sorted(system_counts.items())),
            "counts_by_kind": dict(sorted(kind_counts.items())),
            "provenance_class_occurrences": dict(sorted(provenance_counts.items())),
            "records_without_atomic_source": missing_atomic_sources,
            "existing_static_rules_with_official_vendor_conditions": 0,
            "production_rules_with_target_and_background_frequency": 0,
            "production_rules_with_independent_fragmentation_groups": 0,
        },
        "source_files": {
            relative: _sha256(ROOT / relative)
            for relative in (
                "matching/fragmentation_rules.py",
                "matching/fragment_feature_library.py",
                "matching/fragmentation_rules_curation.json",
                "kb_v3/empirical_rules.json",
                "kb_v3/empirical_rules_candidates_v2.json",
            )
        },
        "knowledge_systems": _knowledge_systems(),
        "fixed_constants": _fixed_constants(),
        "known_defects": _known_defects(),
        "records": records,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--ledger-output", type=Path, default=DEFAULT_LEDGER_OUTPUT)
    args = parser.parse_args()
    audit = build_audit()
    records = audit.pop("records")
    ledger_text = "".join(
        json.dumps(record, ensure_ascii=False, sort_keys=True) + "\n"
        for record in records
    )
    args.ledger_output.parent.mkdir(parents=True, exist_ok=True)
    args.ledger_output.write_text(ledger_text, encoding="utf-8")
    audit["row_level_ledger"] = {
        "path": str(args.ledger_output),
        "format": "JSON Lines; one existing rule occurrence per line",
        "record_count": len(records),
        "sha256": hashlib.sha256(ledger_text.encode("utf-8")).hexdigest(),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(audit, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(
        json.dumps(
            {
                "output": str(args.output),
                "ledger_output": str(args.ledger_output),
                **audit["summary"],
            },
            ensure_ascii=False,
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
