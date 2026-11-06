#!/usr/bin/env python3
"""Build a deterministic, one-record-per-rule curation inventory.

The legacy JSON is valuable source material but is not executable: most entries
mix precursor-relative notation, narrative mechanisms and instrument-dependent
intensity claims.  This script never promotes a rule by text parsing alone.  It
assigns a stable occurrence ID, identifies duplicate/misplaced/incorrect records,
and records the small set that has been independently reimplemented with
polarity, adduct and RDKit structure guards in the production rule modules.
"""

from __future__ import annotations

import hashlib
import json
import re
from collections import Counter
from pathlib import Path
from typing import Any, Iterator


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "matching" / "fragmentation_rules_library.json"
OUTPUT = ROOT / "matching" / "fragmentation_rules_curation.json"


PROMOTED = {
    "FLAV-020": ["GLYCOSIDIC_CLEAVAGE"],
    "FLAV-021": ["C_GLYCOSIDE_CROSS_RING"],
    "LIP-001": ["phosphocholine_184"],
    "PL-001": ["phosphocholine_184"],
    "PL-002": ["PHOSPHATIDYLETHANOLAMINE_141_LOSS", "phosphoethanolamine_196"],
    "SPH-002": ["phosphocholine_184", "sphingoid_d18_1_264"],
    "SUL-001": ["SULFATE_SO3_LOSS"],
}


def iter_rules(value: Any, context: dict[str, str] | None = None) -> Iterator[tuple[dict, dict]]:
    context = dict(context or {})
    if isinstance(value, dict):
        for key in ("class", "class_en", "subclass"):
            if value.get(key):
                context[key] = str(value[key])
        if value.get("rule_id"):
            yield context, value
        for child in value.values():
            yield from iter_rules(child, context)
    elif isinstance(value, list):
        for child in value:
            yield from iter_rules(child, context)


def compact_fingerprint(rule: dict) -> str:
    payload = json.dumps(rule, ensure_ascii=False, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()[:16]


def infer_polarity(text: str) -> str:
    lower = text.lower()
    positive = bool(re.search(r"正离子|positive|\[m\+h[^]]*\]\+", lower))
    negative = bool(re.search(r"负离子|negative|\[m-h[^]]*\]-", lower))
    if positive and negative:
        return "both_candidate_only"
    if positive:
        return "positive_candidate_only"
    if negative:
        return "negative_candidate_only"
    return "unknown"


def infer_ionization(text: str) -> list[str]:
    found = []
    for label, pattern in (
        ("ESI", r"\besi\b|electrospray"),
        ("APCI", r"\bapci\b"),
        ("EI", r"\bei\b|electron ionization"),
        ("CI", r"\bci\b|chemical ionization"),
    ):
        if re.search(pattern, text, flags=re.I):
            found.append(label)
    return found


def reference_ids(reference: str) -> list[str]:
    patterns = (
        r"PMCID\s*:?\s*PMC\d+",
        r"PMID\s*:?\s*\d+",
        r"DOI\s*:?\s*10\.\d{4,9}/[^;\s]+",
    )
    found = []
    for pattern in patterns:
        found.extend(re.findall(pattern, reference or "", flags=re.I))
    return sorted({re.sub(r"\s+", "", item) for item in found})


def direct_fragment_mz(value: Any) -> float | None:
    text = str(value or "").strip()
    match = re.fullmatch(r"(?:m/z\s*)?(\d+(?:\.\d+)?)", text, flags=re.I)
    return float(match.group(1)) if match else None


def classify_status(
    rule: dict,
    context: dict[str, str],
    fingerprint: str,
    seen_fingerprints: set[str],
) -> tuple[str, list[str], str]:
    rule_id = str(rule["rule_id"])
    text = json.dumps(rule, ensure_ascii=False)
    if rule_id == "ALK-009b" and context.get("class_en") == "Terpenoids":
        return "quarantined_misplaced", [], "乌头碱规则误置于萜类章节；保留生物碱章节来源"
    if fingerprint in seen_fingerprints:
        return "quarantined_exact_duplicate", [], "内容完全重复，不得重复计分"
    if rule_id == "LIP-002":
        return (
            "rejected_incorrect",
            [],
            "将PC与PE的141.019 Da头基丢失混写；141.019 Da生产规则仅对PE结构启用",
        )
    if rule_id == "NUC-001" and "脱氧核糖" not in text:
        return "superseded_source_record", [], "由同ID的核糖/脱氧核糖完整记录取代"
    if rule_id == "NUC-001" and "脱氧核糖" in text:
        return (
            "promoted_runtime",
            [
                "RIBONUCLEOSIDE_GLYCOSIDIC_CLEAVAGE",
                "DEOXYRIBONUCLEOSIDE_GLYCOSIDIC_CLEAVAGE",
                "nucleoside_base_diagnostic_ions",
            ],
            "已拆成核糖132.0423、脱氧核糖116.0473及碱基离子原子规则",
        )
    if rule_id == "NUC-004":
        return (
            "promoted_runtime",
            ["RIBONUCLEOSIDE_GLYCOSIDIC_CLEAVAGE", "nucleoside_base_diagnostic_ions"],
            "已用核苷结构守卫、[M+H]+及精确碱基离子重建",
        )
    if rule_id == "ALK-009b":
        return (
            "promoted_runtime",
            ["ChemicalExpertScorer._detect_alkaloid_ester_loss_cascade"],
            "已作为乌头碱乙酸/苯甲酸连续丢失路径实现",
        )
    if rule_id in {"FLAV-023", "FLAV-032"}:
        return (
            "held_collision_energy_specific",
            [],
            "强度比或竞争裂解依赖碰撞能/仪器；缺少可转移的仪器与CE单位元数据",
        )
    if rule_id == "FLAV-027" and "SO3" in text:
        return "promoted_runtime", ["SULFATE_SO3_LOSS"], "已用硫酸酯结构与负离子适用域重建"
    if rule_id == "GLY-001" and context.get("class_en") == "General Glycoside Rules":
        return "promoted_runtime", ["GLYCOSIDIC_CLEAVAGE"], "已用O-糖苷结构守卫重新实现"
    if rule_id in PROMOTED:
        return "promoted_runtime", PROMOTED[rule_id], "已用结构化适用域重新实现"
    return (
        "pending_atomicization",
        [],
        "尚需拆分数值产物/中性丢失并补齐结构、极性、加合物和证据适用域",
    )


def build_inventory() -> dict:
    data = json.loads(SOURCE.read_text(encoding="utf-8"))
    source_records = list(iter_rules(data.get("classes", [])))
    occurrences: Counter[str] = Counter()
    seen_fingerprints: set[str] = set()
    records = []

    for context, rule in source_records:
        rule_id = str(rule["rule_id"])
        occurrences[rule_id] += 1
        occurrence_id = f"{rule_id}#{occurrences[rule_id]:02d}"
        fingerprint = compact_fingerprint(rule)
        status, runtime_links, rationale = classify_status(
            rule, context, fingerprint, seen_fingerprints
        )
        seen_fingerprints.add(fingerprint)
        text = json.dumps(rule, ensure_ascii=False)
        records.append({
            "occurrence_id": occurrence_id,
            "source_rule_id": rule_id,
            "source_fingerprint": fingerprint,
            "class": context.get("class", ""),
            "class_en": context.get("class_en", ""),
            "subclass": context.get("subclass", ""),
            "status": status,
            "status_rationale": rationale,
            "runtime_links": runtime_links,
            "candidate_fragment_mz": direct_fragment_mz(rule.get("fragment_mz")),
            "candidate_polarity": infer_polarity(text),
            "candidate_ionization_methods": infer_ionization(text),
            "reference_identifiers": reference_ids(str(rule.get("reference") or "")),
            "source_confidence_label": rule.get("confidence"),
            "execution_gate": {
                "atomic_numeric_rule": status == "promoted_runtime",
                "polarity_curated": status == "promoted_runtime",
                "adduct_curated": status == "promoted_runtime",
                "structure_guard_curated": status == "promoted_runtime",
                "collision_energy_scope_curated": False,
                "production_enabled": status == "promoted_runtime",
            },
        })

    status_counts = Counter(record["status"] for record in records)
    return {
        "meta": {
            "schema_version": "1.0",
            "source": str(SOURCE.relative_to(ROOT)),
            "record_count": len(records),
            "unique_source_rule_ids": len({record["source_rule_id"] for record in records}),
            "status_counts": dict(sorted(status_counts.items())),
            "policy": (
                "Only promoted_runtime records are represented by independently curated "
                "production rules; inferred text fields never activate scoring."
            ),
        },
        "records": records,
    }


def main() -> None:
    inventory = build_inventory()
    OUTPUT.write_text(
        json.dumps(inventory, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(inventory["meta"], ensure_ascii=False, indent=2))
    print(f"wrote {OUTPUT}")


if __name__ == "__main__":
    main()
