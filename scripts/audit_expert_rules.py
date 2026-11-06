#!/usr/bin/env python3
"""Reproducible audit of V6 expert-rule data and optional cached task triggers."""

from __future__ import annotations

import argparse
import json
import re
import sqlite3
import sys
import urllib.parse
import urllib.request
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from matching.expert_scorer import ChemicalExpertScorer  # noqa: E402
from matching.fragment_feature_library import (  # noqa: E402
    COMPOUND_CLASS_PROFILES,
    get_fragments_by_class,
)
from matching.fragmentation_rules import (  # noqa: E402
    COMMON_NEUTRAL_LOSSES,
    DIAGNOSTIC_IONS,
    FRAGMENTATION_RULES,
    FUNCTIONAL_GROUP_SIGNATURES,
)
from matching.rule_quality import formula_mz_error, normalize_polarity  # noqa: E402


def _flatten_rule_dicts(value: Any) -> list[dict[str, Any]]:
    found: list[dict[str, Any]] = []
    if isinstance(value, dict):
        if value.get("rule_id"):
            found.append(value)
        for child in value.values():
            found.extend(_flatten_rule_dicts(child))
    elif isinstance(value, list):
        for child in value:
            found.extend(_flatten_rule_dicts(child))
    return found


def _json_library_audit() -> dict[str, Any]:
    path = ROOT / "matching" / "fragmentation_rules_library.json"
    data = json.loads(path.read_text(encoding="utf-8"))
    rules = _flatten_rule_dicts(data.get("classes", []))
    ids = [str(rule.get("rule_id")) for rule in rules]
    duplicate_ids = sorted(key for key, count in Counter(ids).items() if count > 1)

    direct_numeric = 0
    contains_numeric = 0
    explicit_polarity = 0
    mode_counts: Counter[str] = Counter()
    for rule in rules:
        fragment = rule.get("fragment_mz")
        text_fragment = "" if fragment is None else str(fragment).strip()
        if re.fullmatch(r"(?:m/z\s*)?\d+(?:\.\d+)?", text_fragment, flags=re.I):
            direct_numeric += 1
        if re.search(r"\d", text_fragment):
            contains_numeric += 1
        if any(key in rule for key in ("polarity", "ion_mode", "ionization_mode")):
            explicit_polarity += 1

        text = json.dumps(rule, ensure_ascii=False).lower()
        if "apci" in text:
            mode_counts["APCI"] += 1
        elif re.search(r"\b(ei|electron ionization)\b", text):
            mode_counts["EI"] += 1
        elif re.search(r"\b(ci|chemical ionization)\b", text):
            mode_counts["CI"] += 1
        elif "esi" in text or "electrospray" in text:
            mode_counts["ESI"] += 1
        else:
            mode_counts["unspecified"] += 1

    accessor_names = (
        "load_compound_class_rules",
        "get_compound_class_fragments",
        "get_compound_class_neutral_losses",
        "get_sugar_loss_reference",
        "get_aglycone_reference",
    )
    external_consumers = []
    for py_file in ROOT.rglob("*.py"):
        if (
            py_file.name == "fragmentation_rules.py"
            or py_file == Path(__file__).resolve()
            or py_file.name == "curate_legacy_rule_inventory.py"
        ):
            continue
        try:
            source = py_file.read_text(encoding="utf-8")
        except (OSError, UnicodeDecodeError):
            continue
        if "fragmentation_rules_library.json" in source or any(
            re.search(rf"\b{re.escape(name)}\b", source) for name in accessor_names
        ):
            external_consumers.append(str(py_file.relative_to(ROOT)))

    return {
        "path": str(path.relative_to(ROOT)),
        "classes": len(data.get("classes", [])),
        "rules": len(rules),
        "unique_rule_ids": len(set(ids)),
        "duplicate_rule_ids": duplicate_ids,
        "direct_numeric_fragment_mz": direct_numeric,
        "fragment_mz_with_any_number": contains_numeric,
        "non_numeric_or_missing_fragment_mz": len(rules) - contains_numeric,
        "explicit_structured_polarity": explicit_polarity,
        "ionization_mode_markers": dict(mode_counts),
        "external_runtime_consumers": sorted(external_consumers),
    }


def _runtime_audit() -> dict[str, Any]:
    diagnostic_bad = []
    for key, ion in DIAGNOSTIC_IONS.items():
        error = formula_mz_error(ion.formula, ion.mz)
        if error is None or abs(error) > 0.01:
            diagnostic_bad.append({
                "id": key, "mz": ion.mz, "formula": ion.formula, "error_da": error,
            })

    feature_bad = []
    class_rows = []
    for class_name, profile in COMPOUND_CLASS_PROFILES.items():
        for feature in profile.characteristic_fragments:
            error = formula_mz_error(feature.formula, feature.mz)
            if error is None or abs(error) > 0.01:
                feature_bad.append({
                    "class": class_name,
                    "mz": feature.mz,
                    "formula": feature.formula,
                    "error_da": error,
                    "name": feature.name,
                })
        raw = len(profile.characteristic_fragments)
        effective = len(get_fragments_by_class(class_name))
        class_rows.append({
            "class": class_name,
            "raw": raw,
            "effective": effective,
            "quarantined": raw - effective,
            "positive_effective": len(get_fragments_by_class(class_name, polarity="POS")),
            "negative_effective": len(get_fragments_by_class(class_name, polarity="NEG")),
        })

    return {
        "neutral_losses": len(COMMON_NEUTRAL_LOSSES),
        "mechanism_rules": len(FRAGMENTATION_RULES),
        "mechanism_rule_polarity": dict(Counter(
            rule.charge_preference for rule in FRAGMENTATION_RULES.values()
        )),
        "mechanism_rules_structure_guarded": sum(
            bool(rule.required_smarts or rule.required_smarts_all)
            for rule in FRAGMENTATION_RULES.values()
        ),
        "mechanism_rules_adduct_guarded": sum(
            bool(rule.allowed_adducts) for rule in FRAGMENTATION_RULES.values()
        ),
        "mechanism_rules_ce_guarded": sum(
            rule.collision_energy_min is not None or rule.collision_energy_max is not None
            for rule in FRAGMENTATION_RULES.values()
        ),
        "diagnostic_ions": len(DIAGNOSTIC_IONS),
        "diagnostic_ion_polarity": dict(Counter(
            ion.polarity for ion in DIAGNOSTIC_IONS.values()
        )),
        "diagnostic_formula_mz_failures": diagnostic_bad,
        "diagnostic_ions_structure_guarded": sum(
            bool(ion.required_smarts or ion.required_smarts_all)
            for ion in DIAGNOSTIC_IONS.values()
        ),
        "diagnostic_ions_adduct_guarded": sum(
            bool(ion.allowed_adducts) for ion in DIAGNOSTIC_IONS.values()
        ),
        "functional_group_signatures": len(FUNCTIONAL_GROUP_SIGNATURES),
        "feature_profiles": len(COMPOUND_CLASS_PROFILES),
        "feature_records_raw": sum(row["raw"] for row in class_rows),
        "feature_records_effective": sum(row["effective"] for row in class_rows),
        "feature_records_structure_guarded": sum(
            bool(feature.required_smarts or feature.required_smarts_all)
            for profile in COMPOUND_CLASS_PROFILES.values()
            for feature in profile.characteristic_fragments
        ),
        "feature_records_adduct_guarded": sum(
            bool(feature.allowed_adducts)
            for profile in COMPOUND_CLASS_PROFILES.values()
            for feature in profile.characteristic_fragments
        ),
        "feature_formula_mz_quarantined": feature_bad,
        "feature_classes": class_rows,
    }


def _curation_audit() -> dict[str, Any]:
    path = ROOT / "matching" / "fragmentation_rules_curation.json"
    if not path.exists():
        return {"available": False}
    data = json.loads(path.read_text(encoding="utf-8"))
    records = data.get("records", [])
    return {
        "available": True,
        "path": str(path.relative_to(ROOT)),
        "record_count": len(records),
        "unique_occurrence_ids": len({r.get("occurrence_id") for r in records}),
        "status_counts": dict(Counter(r.get("status") for r in records)),
        "production_enabled": sum(
            bool(r.get("execution_gate", {}).get("production_enabled")) for r in records
        ),
        "runtime_links": sorted({
            link for record in records for link in record.get("runtime_links", [])
        }),
    }


def _empirical_candidate_audit() -> dict[str, Any]:
    path = ROOT / "kb_v3" / "empirical_rules_candidates_v2.json"
    if not path.exists():
        return {"available": False}
    data = json.loads(path.read_text(encoding="utf-8"))
    class_rows = []
    diagnostic_total = neutral_loss_total = 0
    for class_name, polarities in data.get("classes", {}).items():
        for polarity, payload in polarities.items():
            diagnostic_count = len(payload.get("diagnostic_ions", []))
            neutral_loss_count = len(payload.get("neutral_losses", []))
            diagnostic_total += diagnostic_count
            neutral_loss_total += neutral_loss_count
            class_rows.append({
                "class": class_name,
                "polarity": polarity,
                "n_independent_compounds": payload.get("n_independent_compounds", 0),
                "diagnostic_candidates": diagnostic_count,
                "neutral_loss_candidates": neutral_loss_count,
            })
    return {
        "available": True,
        "path": str(path.relative_to(ROOT)),
        "runtime_enabled": bool(data.get("meta", {}).get("runtime_enabled")),
        "rows_seen": data.get("meta", {}).get("rows_seen", 0),
        "rows_used": data.get("meta", {}).get("rows_used", 0),
        "unique_identities_seen": data.get("meta", {}).get("unique_identities_seen", 0),
        "classes": len(data.get("classes", {})),
        "diagnostic_candidates": diagnostic_total,
        "neutral_loss_candidates": neutral_loss_total,
        "class_rows": class_rows,
        "collision_energy_policy": data.get("meta", {}).get("collision_energy_policy", ""),
    }


def _library_anchor_audit() -> list[dict[str, Any]]:
    db_path = ROOT / "data" / "local_library.db"
    if not db_path.exists():
        return []
    anchors = (
        ("adenosine", "lower(c.compound_name)='adenosine'", 268.104, "RIBONUCLEOSIDE_GLYCOSIDIC_CLEAVAGE"),
        ("2-deoxyadenosine", "lower(c.compound_name) like '%2''-deoxyadenosine%'", 252.109, "DEOXYRIBONUCLEOSIDE_GLYCOSIDIC_CLEAVAGE"),
        ("5-deoxyadenosine", "lower(c.compound_name) like '%5''-deoxyadenosine%'", 252.109, "DEOXYRIBONUCLEOSIDE_GLYCOSIDIC_CLEAVAGE"),
    )
    scorer = ChemicalExpertScorer()
    output = []
    connection = sqlite3.connect(f"file:{db_path}?mode=ro", uri=True)
    connection.row_factory = sqlite3.Row
    try:
        for label, name_where, expected_precursor, expected_rule in anchors:
            row = connection.execute(f"""
                SELECT c.compound_name, c.smiles, c.molecular_formula,
                       c.classification_v3, s.precursor_mz, s.peaks_json,
                       s.adduct, s.charge, s.dataset_source
                FROM Compounds c
                JOIN Mapping m ON m.compound_id=c.id
                JOIN Spectra s ON s.id=m.spectrum_id
                WHERE {name_where}
                  AND s.source_tier='experimental'
                  AND s.adduct='[M+H]+'
                  AND ABS(s.precursor_mz-?) < 0.02
                  AND s.peaks_json LIKE '%136.06%'
                ORDER BY length(s.peaks_json) DESC
                LIMIT 1
            """, (expected_precursor,)).fetchone()
            if row is None:
                output.append({"anchor": label, "status": "not_found"})
                continue
            details = scorer.evaluate(
                query_peaks=json.loads(row["peaks_json"]),
                query_pmz=row["precursor_mz"],
                candidate_smiles=row["smiles"] or "",
                candidate_formula=row["molecular_formula"] or "",
                candidate_classification=row["classification_v3"] or "Unknown",
                query_polarity="POS",
                query_adduct=row["adduct"] or "",
                query_charge=row["charge"] or 1,
            )
            covered = [
                item["rule_id"]
                for item in details.get("rule_coverage_analysis", {}).get("covered_rules", [])
            ]
            output.append({
                "anchor": label,
                "status": "passed" if expected_rule in covered else "failed",
                "compound_name": row["compound_name"],
                "dataset_source": row["dataset_source"],
                "source_classification": row["classification_v3"],
                "effective_classification": details.get("effective_classification"),
                "covered_rules": covered,
                "diagnostic_ions": [
                    item.get("ion_name")
                    for item in details.get("diagnostic_ion_analysis", {}).get("matched_ions", [])
                ],
            })
    finally:
        connection.close()
    return output


def _get_json(url: str) -> Any:
    with urllib.request.urlopen(url, timeout=60) as response:
        return json.load(response)


def _old_diagnostic_cross_polarity(details: dict[str, Any], polarity: str) -> int:
    wanted = normalize_polarity(polarity)
    total = 0
    for match in details.get("diagnostic_ion_analysis", {}).get("matched_ions", []):
        expected = float(match.get("expected_mz", -999))
        entries = [ion for ion in DIAGNOSTIC_IONS.values() if abs(ion.mz - expected) <= 0.05]
        if entries and all(normalize_polarity(ion.polarity) not in (None, wanted) for ion in entries):
            total += 1
    return total


def _old_feature_cross_polarity(
    details: dict[str, Any], class_name: str, polarity: str
) -> int:
    wanted = normalize_polarity(polarity)
    profile = COMPOUND_CLASS_PROFILES.get(str(class_name or "").lower())
    if not profile:
        return 0
    total = 0
    for match in details.get("feature_library_analysis", {}).get("matched_features", []):
        expected = float(match.get("expected_mz", -999))
        entries = [
            feature for feature in profile.characteristic_fragments
            if abs(feature.mz - expected) <= 0.05
        ]
        if entries and all(normalize_polarity(feature.polarity) not in (None, wanted) for feature in entries):
            total += 1
    return total


def _task_trigger_audit(api_base: str, task_id: str) -> dict[str, Any]:
    task = _get_json(f"{api_base}/api/tasks/{urllib.parse.quote(task_id)}")
    payload = _get_json(
        f"{api_base}/api/tasks/{urllib.parse.quote(task_id)}/results"
        "?matched_only=true&limit=10000"
    )
    rows = payload.get("results", [])
    scorer = ChemicalExpertScorer()
    counts: Counter[str] = Counter()
    old_di_names: Counter[str] = Counter()
    new_di_names: Counter[str] = Counter()

    for row in rows:
        match = row.get("match") or {}
        details = match.get("expert_details") or {}
        polarity = row.get("ion_mode") or ""
        analysis_class = details.get("analysis_class") or match.get("classification") or "Unknown"
        peaks = row.get("query_peaks") or []
        precursor = float(row.get("precursor_mz") or 0)
        fragment_mzs = [
            float(peak[0]) for peak in peaks
            if peak and len(peak) >= 2 and float(peak[0]) < precursor - 1.5
        ]

        old_di = details.get("diagnostic_ion_analysis", {}).get("matched_ions", [])
        old_ff = details.get("feature_library_analysis", {}).get("matched_features", [])
        counts["old_diagnostic_hits"] += len(old_di)
        counts["old_feature_hits"] += len(old_ff)
        counts["old_cross_polarity_diagnostic_hits"] += _old_diagnostic_cross_polarity(
            details, polarity
        )
        counts["old_cross_polarity_feature_hits"] += _old_feature_cross_polarity(
            details, analysis_class, polarity
        )
        for item in old_di:
            old_di_names[str(item.get("ion_name") or item.get("expected_mz"))] += 1

        new_di = scorer._analyze_diagnostic_ions(
            fragment_mzs,
            classification=analysis_class,
            polarity=polarity,
            candidate_smiles=match.get("smiles") or "",
            raw_classification=match.get("raw_classification") or match.get("classification") or "",
            query_adduct=row.get("adduct") or match.get("query_adduct") or "",
            query_collision_energy=row.get("collision_energy"),
        )
        new_ff = scorer._analyze_feature_library(
            fragment_mzs,
            classification=analysis_class,
            polarity=polarity,
            candidate_smiles=match.get("smiles") or "",
            raw_classification=match.get("raw_classification") or match.get("classification") or "",
        )
        counts["new_diagnostic_hits"] += new_di.get("matched_count", 0)
        counts["new_feature_hits"] += new_ff.get("matched_count", 0)
        counts["rows_old_diagnostic_nonzero"] += int(bool(old_di))
        counts["rows_new_diagnostic_nonzero"] += int(bool(new_di.get("matched_ions")))
        counts["rows_old_feature_nonzero"] += int(bool(old_ff))
        counts["rows_new_feature_nonzero"] += int(bool(new_ff.get("matched_features")))
        for item in new_di.get("matched_ions", []):
            new_di_names[str(item.get("ion_name") or item.get("expected_mz"))] += 1

    return {
        "task_id": task_id,
        "filename": task.get("filename") or task.get("name") or "",
        "polarity": task.get("polarity") or (rows[0].get("ion_mode") if rows else ""),
        "matched_rows": len(rows),
        **dict(counts),
        "old_top_diagnostic_ions": old_di_names.most_common(10),
        "new_top_diagnostic_ions": new_di_names.most_common(10),
    }


def _markdown(report: dict[str, Any]) -> str:
    runtime = report["runtime"]
    json_lib = report["json_library"]
    curation = report.get("curation", {})
    empirical = report.get("empirical_candidates", {})
    lines = [
        "# V6 专家规则库数据审计",
        "",
        f"生成时间：{report['generated_at']}",
        "",
        "## 结论",
        "",
        "- 专家知识可以辅助鉴别，但只能作为有适用域、可追溯的支持证据，不能替代标准品谱图或用规则数量直接制造 High。",
        "- 生产规则此前存在极性未执行、HRMS 容差过宽、公式与 m/z 冲突、规则缺失时仍给正分等问题；本审计对应的质量门已在运行代码中落实。",
        "- 综合 JSON 规则库目前适合作为待策展语料，不适合直接接入评分：它缺少结构化极性/加合物/碰撞能字段，并混合 ESI、EI、CI、APCI 规则。",
        "",
        "## 生产规则库",
        "",
        f"- 中性丢失：{runtime['neutral_losses']} 条",
        f"- 机理规则：{runtime['mechanism_rules']} 条；极性分布 {runtime['mechanism_rule_polarity']}",
        f"- 机理规则适用域：结构守卫 {runtime['mechanism_rules_structure_guarded']} 条；加合物守卫 {runtime['mechanism_rules_adduct_guarded']} 条；CE守卫 {runtime['mechanism_rules_ce_guarded']} 条",
        f"- 诊断离子：{runtime['diagnostic_ions']} 条；极性分布 {runtime['diagnostic_ion_polarity']}",
        f"- 诊断离子适用域：结构守卫 {runtime['diagnostic_ions_structure_guarded']} 条；加合物守卫 {runtime['diagnostic_ions_adduct_guarded']} 条",
        f"- 官能团签名：{runtime['functional_group_signatures']} 类",
        f"- 碎片类别档案：{runtime['feature_profiles']} 类，原始 {runtime['feature_records_raw']} 条，生产有效 {runtime['feature_records_effective']} 条，隔离 {len(runtime['feature_formula_mz_quarantined'])} 条",
        f"- 碎片特征适用域：结构守卫 {runtime['feature_records_structure_guarded']} 条；加合物守卫 {runtime['feature_records_adduct_guarded']} 条",
        f"- 诊断离子公式/m/z 失败：{len(runtime['diagnostic_formula_mz_failures'])} 条",
        "",
        "### 被隔离的碎片条目",
        "",
        "| 类别 | m/z | 分子式 | 误差(Da) | 名称 |",
        "|---|---:|---|---:|---|",
    ]
    for item in runtime["feature_formula_mz_quarantined"]:
        error = item["error_da"]
        error_text = "不可计算" if error is None else f"{error:.6f}"
        lines.append(
            f"| {item['class']} | {item['mz']:.4f} | {item['formula']} | {error_text} | {item['name']} |"
        )

    lines.extend([
        "",
        "## 综合 JSON 规则库",
        "",
        f"- 类别：{json_lib['classes']}；规则记录：{json_lib['rules']}；唯一 ID：{json_lib['unique_rule_ids']}",
        f"- 重复 ID：{', '.join(json_lib['duplicate_rule_ids']) or '无'}",
        f"- 可直接解析为单一数值 m/z：{json_lib['direct_numeric_fragment_mz']}；仅包含某些数字：{json_lib['fragment_mz_with_any_number']}；非数值/缺失：{json_lib['non_numeric_or_missing_fragment_mz']}",
        f"- 结构化极性字段：{json_lib['explicit_structured_polarity']} 条",
        f"- 电离方式文本标记：{json_lib['ionization_mode_markers']}",
        f"- 生产代码外部消费者：{json_lib['external_runtime_consumers'] or '无（当前未接入生产评分）'}",
    ])

    if curation.get("available"):
        lines.extend([
            "",
            "## 逐条策展台账",
            "",
            f"- 台账记录：{curation['record_count']}；唯一 occurrence ID：{curation['unique_occurrence_ids']}",
            f"- 状态分布：{curation['status_counts']}",
            f"- 已由结构化运行时规则替代：{curation['production_enabled']} 条来源记录",
            f"- 运行时链接：{', '.join(curation['runtime_links'])}",
            "- 台账中的文本推断字段不直接参与评分；只有 promoted_runtime 对应的独立实现可以进入生产。",
        ])

    if empirical.get("available"):
        lines.extend([
            "",
            "## 对照式经验规则候选",
            "",
            f"- 读取谱图：{empirical['rows_seen']:,}；可用：{empirical['rows_used']:,}；独立 InChIKey/结构：{empirical['unique_identities_seen']:,}",
            f"- 覆盖类别：{empirical['classes']}；诊断离子候选：{empirical['diagnostic_candidates']}；中性丢失候选：{empirical['neutral_loss_candidates']}",
            f"- 运行时启用：{empirical['runtime_enabled']}（候选未经结构/公式/基准验证，不自动接入）",
            f"- 碰撞能策略：{empirical['collision_energy_policy']}",
            "",
            "| 类别 | 极性 | 独立化合物 | 诊断离子候选 | 中性丢失候选 |",
            "|---|---|---:|---:|---:|",
        ])
        for row in empirical["class_rows"]:
            lines.append(
                f"| {row['class']} | {row['polarity']} | {row['n_independent_compounds']} | "
                f"{row['diagnostic_candidates']} | {row['neutral_loss_candidates']} |"
            )

    if report.get("library_anchor_audits"):
        lines.extend(["", "## 共享谱库锚点验证", ""])
        for anchor in report["library_anchor_audits"]:
            lines.append(
                f"- {anchor.get('anchor')}: {anchor.get('status')}；"
                f"原分类={anchor.get('source_classification', 'n/a')}；"
                f"纠正后={anchor.get('effective_classification', 'n/a')}；"
                f"规则={anchor.get('covered_rules', [])}；"
                f"诊断离子={anchor.get('diagnostic_ions', [])}；"
                f"来源={anchor.get('dataset_source', 'n/a')}"
            )

    if report.get("task_trigger_audits"):
        lines.extend(["", "## 缓存真实任务触发审计", ""])
        for task in report["task_trigger_audits"]:
            lines.extend([
                f"### {task['task_id']} {task.get('filename', '')}",
                "",
                f"- 匹配候选：{task['matched_rows']}",
                f"- 旧诊断离子命中：{task.get('old_diagnostic_hits', 0)}；其中跨极性：{task.get('old_cross_polarity_diagnostic_hits', 0)}",
                f"- 新诊断离子命中：{task.get('new_diagnostic_hits', 0)}（跨极性规则已在入口过滤）",
                f"- 旧碎片特征命中：{task.get('old_feature_hits', 0)}；其中跨极性：{task.get('old_cross_polarity_feature_hits', 0)}",
                f"- 新碎片特征命中：{task.get('new_feature_hits', 0)}（含公式质量门、极性门和窄类结构守卫）",
                f"- 旧 Top 诊断离子：{task.get('old_top_diagnostic_ions', [])}",
                f"- 新 Top 诊断离子：{task.get('new_top_diagnostic_ions', [])}",
                "",
            ])

    lines.extend([
        "## 后续策展要求",
        "",
        "1. 每条可执行规则必须有极性、加合物、碰撞能适用范围、精确式/m/z、类别或 SMARTS 守卫、证据等级和可核查来源。",
        "2. JSON 文本规则需先拆为原子化数值规则、消除重复 ID、区分 ESI 与 EI/CI/APCI，再通过标准谱和真实样本回归，方可接入评分。",
        "3. High 的提升应来自独立证据汇合（谱库相似度、前体质量、结构一致性、同位素、诊断碎片/中性丢失、RT/标准品），而不是单个专家规则反复计分。",
        "",
    ])
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--api-base", default="http://127.0.0.1:8506")
    parser.add_argument("--task", action="append", default=[])
    parser.add_argument("--output", type=Path)
    parser.add_argument("--json-output", type=Path)
    args = parser.parse_args()

    report = {
        "generated_at": datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds"),
        "runtime": _runtime_audit(),
        "json_library": _json_library_audit(),
        "curation": _curation_audit(),
        "empirical_candidates": _empirical_candidate_audit(),
        "library_anchor_audits": _library_anchor_audit(),
        "task_trigger_audits": [
            _task_trigger_audit(args.api_base.rstrip("/"), task_id)
            for task_id in args.task
        ],
    }
    markdown = _markdown(report)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(markdown, encoding="utf-8")
    else:
        print(markdown)
    if args.json_output:
        args.json_output.parent.mkdir(parents=True, exist_ok=True)
        args.json_output.write_text(
            json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
