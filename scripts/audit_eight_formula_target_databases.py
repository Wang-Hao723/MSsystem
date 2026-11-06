#!/usr/bin/env python3
"""Reproducible, read-only audit of the eight formula targeted catalogues.

The catalogue rows are literature-presence targets, not experimental spectra
and not confirmed sample identifications.  This audit checks whether those
rows can be executed safely by the MS1 targeted-screening layer while keeping
source defects visible.  It never rewrites ``local_library.db``.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import sqlite3
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Mapping

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from matching.expected_compounds import (
    _normalise_adduct,
    build_expected_ion_hypotheses,
    load_expected_compound_catalog,
    target_catalog_quality_summary,
)
from matching.library_quality import (
    formula_exact_mass,
    normalize_formula_text,
    normalize_inchikey,
)
from utils.adduct import ion_mz


DEFAULT_DATABASE = Path(
    "/Users/zizi/Desktop/Mssystem_shared_data/local_library.db"
)
DEFAULT_OUTPUT = Path(__file__).resolve().parents[1] / "docs" / "audits"

FORMULA_DATASETS = {
    "F001": "F001_独活寄生汤_targeted_20260817",
    "F002": "F002_四神煎_targeted_20260817",
    "F003": "F003_乌头汤_targeted_20260817",
    "F004": "F004_双合汤_targeted_20260817",
    "F005": "XBT_targeted_20260817",
    "F006": "STZYT_targeted_20260817",
    "F007": "F007_HQGZWW_targeted_20260817",
    "F008": "F008_QHSS_targeted_20260817",
}

FORMULA_NAMES = {
    "F001": "独活寄生汤",
    "F002": "四神煎",
    "F003": "乌头汤",
    "F004": "双合汤",
    "F005": "宣痹汤",
    "F006": "身痛逐瘀汤",
    "F007": "黄芪桂枝五物汤",
    "F008": "羌活胜湿汤",
}

FULL_INCHIKEY = re.compile(r"^[A-Z]{14}-[A-Z]{10}-[A-Z]$")
HEX64 = re.compile(r"^[0-9a-fA-F]{64}$")
EXAMPLE_LIMIT = 20


def _text(value: Any) -> str:
    return str(value or "").strip()


def _number(value: Any) -> float | None:
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        return None
    return parsed if math.isfinite(parsed) else None


def _ratio(value: int, total: int) -> float:
    return round(value / total, 6) if total else 0.0


def _parse_flags(value: Any) -> list[str]:
    try:
        parsed = json.loads(_text(value) or "[]")
    except json.JSONDecodeError:
        return ["invalid_quality_flags_json"]
    if isinstance(parsed, list):
        return [str(item) for item in parsed]
    if isinstance(parsed, dict):
        flags = parsed.get("flags")
        result = [str(item) for item in flags] if isinstance(flags, list) else []
        result.extend(
            f"metadata:{key}" for key in sorted(parsed)
            if key != "flags"
        )
        return result
    return ["invalid_quality_flags_json"]


def _add_example(examples: dict[str, list[dict[str, Any]]], key: str, row: Mapping[str, Any], **extra: Any) -> None:
    if len(examples[key]) >= EXAMPLE_LIMIT:
        return
    examples[key].append({
        "target_id": _text(row["target_id"]),
        "compound_name": _text(row["compound_name"]),
        "formula": _text(row["molecular_formula"]),
        **extra,
    })


def _structure_assessor():
    """Return a cached SMILES assessor, or None if RDKit is unavailable."""

    try:
        from rdkit import Chem, RDLogger
        from rdkit.Chem import rdMolDescriptors
    except ImportError:
        return None
    RDLogger.DisableLog("rdApp.error")

    cache: dict[str, dict[str, Any]] = {}

    def assess(smiles: str) -> dict[str, Any]:
        if smiles in cache:
            return cache[smiles]
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            result = {"valid": False, "formula": "", "inchikey": ""}
        else:
            try:
                inchikey = Chem.MolToInchiKey(mol)
            except Exception:
                inchikey = ""
            result = {
                "valid": True,
                "formula": normalize_formula_text(
                    rdMolDescriptors.CalcMolFormula(mol)
                ),
                "inchikey": normalize_inchikey(inchikey),
            }
        cache[smiles] = result
        return result

    return assess


def _open_read_only(path: Path) -> sqlite3.Connection:
    connection = sqlite3.connect(f"file:{path.resolve()}?mode=ro", uri=True)
    connection.row_factory = sqlite3.Row
    return connection


def _query_rows(connection: sqlite3.Connection, dataset: str) -> Iterable[sqlite3.Row]:
    return connection.execute(
        """
        SELECT t.*,
               COALESCE(d.diagnostic_ion_count, 0) AS diagnostic_ion_count
        FROM TargetCompounds t
        LEFT JOIN (
            SELECT target_id, COUNT(*) AS diagnostic_ion_count
            FROM DiagnosticIons
            GROUP BY target_id
        ) d ON d.target_id=t.id
        WHERE t.dataset_source=?
        ORDER BY t.id
        """,
        (dataset,),
    )


def _audit_dataset(
    connection: sqlite3.Connection,
    formula_id: str,
    dataset: str,
    assess_structure,
    runtime_targets: list[dict[str, Any]],
) -> dict[str, Any]:
    counts: Counter[str] = Counter()
    source_tiers: Counter[str] = Counter()
    quality_flags: Counter[str] = Counter()
    searchability_labels: Counter[str] = Counter()
    canonical_adducts: Counter[str] = Counter()
    quarantined_adduct_values: Counter[str] = Counter()
    examples: dict[str, list[dict[str, Any]]] = defaultdict(list)
    entity_keys: Counter[tuple[str, str, str]] = Counter()
    source_row_keys: Counter[tuple[str, str, str]] = Counter()
    fingerprint = hashlib.sha256()

    for row in _query_rows(connection, dataset):
        counts["rows"] += 1
        source_tiers[_text(row["source_tier"]) or "missing"] += 1
        searchability_labels[_text(row["target_searchability"]) or "missing"] += 1
        for flag in _parse_flags(row["quality_flags_json"]):
            quality_flags[flag] += 1

        stable_values = [
            row["target_id"], row["compound_name"], row["molecular_formula"],
            row["smiles"], row["inchikey"], row["exact_mass"],
            row["precursor_mz"], row["adduct"], row["reported_precursor_mz"],
            row["reported_adduct"], row["structure_identity_verified"],
            row["source_content_sha256"], row["source_row_index"],
            row["diagnostic_ion_count"],
        ]
        fingerprint.update(
            ("\x1f".join(_text(value) for value in stable_values) + "\n").encode(
                "utf-8"
            )
        )

        name = _text(row["compound_name"])
        formula_raw = _text(row["molecular_formula"])
        formula = normalize_formula_text(formula_raw)
        formula_mass = formula_exact_mass(formula)
        exact_mass = _number(row["exact_mass"])
        smiles = _text(row["smiles"])
        inchikey_raw = _text(row["inchikey"]).upper()
        inchikey = normalize_inchikey(inchikey_raw)
        verified = bool(row["structure_identity_verified"])

        for field, present in (
            ("name", bool(name)),
            ("formula", bool(formula)),
            ("smiles", bool(smiles)),
            ("inchikey", bool(inchikey)),
            ("exact_mass", exact_mass is not None),
            ("precursor_mz", _number(row["precursor_mz"]) is not None),
            ("reported_precursor_mz", _number(row["reported_precursor_mz"]) is not None),
            ("adduct", bool(_text(row["adduct"]))),
            ("source_hash", bool(HEX64.fullmatch(_text(row["source_content_sha256"])))),
            ("source_row_index", row["source_row_index"] is not None),
            ("provenance", bool(_text(row["provenance"]))),
            ("external_accession", bool(
                _text(row["doi"]) or _text(row["dataset_accession"])
                or _text(row["spectrum_accession"])
            )),
            ("diagnostic_ions", int(row["diagnostic_ion_count"] or 0) > 0),
            ("structure_identity_verified", verified),
        ):
            if present:
                counts[f"present_{field}"] += 1

        if formula_raw and not formula:
            counts["invalid_formula"] += 1
            _add_example(examples, "invalid_formula", row)
        elif formula:
            counts["valid_formula"] += 1
        if formula_mass is not None and exact_mass is not None and formula_mass > 0:
            ppm = (exact_mass - formula_mass) / formula_mass * 1e6
            counts["formula_exact_mass_comparable"] += 1
            if abs(ppm) <= 10:
                counts["formula_exact_mass_within_10ppm"] += 1
            elif abs(ppm) > 20:
                counts["formula_exact_mass_conflict_over_20ppm"] += 1
                _add_example(
                    examples, "formula_exact_mass_conflict", row,
                    exact_mass=exact_mass,
                    calculated_mass=round(formula_mass, 8),
                    error_ppm=round(ppm, 3),
                )

        if inchikey_raw and not FULL_INCHIKEY.fullmatch(inchikey_raw):
            counts["malformed_inchikey"] += 1
            _add_example(
                examples, "malformed_inchikey", row,
                reported_inchikey=inchikey_raw,
            )

        if smiles:
            if assess_structure is None:
                counts["smiles_not_checked_rdkit_unavailable"] += 1
            else:
                structure = assess_structure(smiles)
                if structure["valid"]:
                    counts["valid_smiles"] += 1
                    if formula and structure["formula"]:
                        counts["smiles_formula_comparable"] += 1
                        if structure["formula"] == formula:
                            counts["smiles_formula_consistent"] += 1
                        else:
                            counts["smiles_formula_mismatch"] += 1
                            if verified:
                                counts["verified_smiles_formula_mismatch"] += 1
                            _add_example(
                                examples, "smiles_formula_mismatch", row,
                                structure_formula=structure["formula"],
                                identity_verified=verified,
                            )
                    if inchikey and structure["inchikey"]:
                        counts["smiles_inchikey_comparable"] += 1
                        if structure["inchikey"] == inchikey:
                            counts["smiles_inchikey_consistent"] += 1
                        elif structure["inchikey"][:14] == inchikey[:14]:
                            # Same 2-D connectivity but the saved SMILES often
                            # omits stereochemistry that is present in the
                            # reported full InChIKey.  This is an
                            # underspecification, not a connectivity error.
                            counts[
                                "smiles_inchikey_connectivity_consistent_full_key_mismatch"
                            ] += 1
                            if verified:
                                counts[
                                    "verified_smiles_stereochemistry_or_protonation_underspecified"
                                ] += 1
                            _add_example(
                                examples,
                                "smiles_stereochemistry_or_protonation_underspecified",
                                row,
                                structure_inchikey=structure["inchikey"],
                                reported_inchikey=inchikey,
                                identity_verified=verified,
                            )
                        else:
                            counts["smiles_inchikey_connectivity_mismatch"] += 1
                            if verified:
                                counts[
                                    "verified_smiles_inchikey_connectivity_mismatch"
                                ] += 1
                            _add_example(
                                examples, "smiles_inchikey_connectivity_mismatch", row,
                                structure_inchikey=structure["inchikey"],
                                reported_inchikey=inchikey,
                                identity_verified=verified,
                            )
                else:
                    counts["invalid_smiles"] += 1
                    if verified:
                        counts["verified_invalid_smiles"] += 1
                    _add_example(
                        examples, "invalid_smiles", row,
                        identity_verified=verified,
                    )

        raw_active_adduct = _text(row["adduct"])
        raw_reported_adduct = _text(row["reported_adduct"]) or raw_active_adduct
        canonical_adduct = _normalise_adduct(raw_active_adduct)
        if canonical_adduct:
            counts["canonical_supported_adduct"] += 1
            canonical_adducts[canonical_adduct] += 1
        elif raw_reported_adduct:
            counts["quarantined_reported_adduct"] += 1
            quarantined_adduct_values[raw_reported_adduct] += 1
            _add_example(
                examples, "quarantined_reported_adduct", row,
                reported_adduct=raw_reported_adduct,
            )
        else:
            counts["adduct_not_reported"] += 1

        neutral_mass = formula_mass if formula_mass is not None else exact_mass
        precursor = (
            _number(row["reported_precursor_mz"])
            if _number(row["reported_precursor_mz"]) is not None
            else _number(row["precursor_mz"])
        )
        theoretical = (
            ion_mz(neutral_mass, canonical_adduct)
            if neutral_mass is not None and canonical_adduct else None
        )
        if precursor is not None and theoretical is not None and theoretical > 0:
            ppm = (precursor - theoretical) / theoretical * 1e6
            counts["precursor_theoretical_comparable"] += 1
            if abs(ppm) <= 10:
                counts["precursor_theoretical_within_10ppm"] += 1
            elif abs(ppm) > 20:
                counts["precursor_theoretical_conflict_over_20ppm"] += 1
                _add_example(
                    examples, "precursor_theoretical_conflict", row,
                    precursor_mz=precursor,
                    calculated_mz=round(theoretical, 8),
                    canonical_adduct=canonical_adduct,
                    error_ppm=round(ppm, 3),
                )

        # Runtime-safe means at least one formula/exact-mass common-adduct
        # hypothesis can be generated, or a canonical reported ion is usable.
        mass_available = neutral_mass is not None and neutral_mass > 0
        reported_ion_usable = (
            precursor is not None and precursor > 0 and bool(canonical_adduct)
        )
        if mass_available or reported_ion_usable:
            counts["runtime_searchable"] += 1
        else:
            counts["runtime_unsearchable"] += 1
            _add_example(examples, "runtime_unsearchable", row)
        if counts["quarantined_reported_adduct"] and mass_available:
            # Per-row condition, intentionally independent of the running count.
            if raw_reported_adduct and not canonical_adduct:
                counts["quarantined_adduct_recovered_by_formula_mass"] += 1

        identity_key = (
            name.casefold(), formula, inchikey or _text(row["reported_inchikey"])
        )
        entity_keys[identity_key] += 1
        source_hash = _text(row["source_content_sha256"])
        source_index = _text(row["source_row_index"])
        if source_hash and source_index:
            source_row_keys[(
                source_hash,
                source_index,
                _text(row["source_file"]),
            )] += 1

    total = counts["rows"]
    duplicated_entities = [count for count in entity_keys.values() if count > 1]
    duplicated_source_rows = [count for count in source_row_keys.values() if count > 1]
    counts["duplicate_entity_groups"] = len(duplicated_entities)
    counts["rows_in_duplicate_entity_groups"] = sum(duplicated_entities)
    counts["duplicate_source_row_target_groups"] = len(duplicated_source_rows)

    critical = (
        counts["runtime_unsearchable"]
        + counts["verified_invalid_smiles"]
        + counts["verified_smiles_formula_mismatch"]
        + counts["verified_smiles_inchikey_connectivity_mismatch"]
    )
    high = (
        counts["precursor_theoretical_conflict_over_20ppm"]
        + counts["formula_exact_mass_conflict_over_20ppm"]
        + counts["quarantined_reported_adduct"]
    )

    coverage_fields = (
        "name", "formula", "smiles", "inchikey", "exact_mass",
        "precursor_mz", "reported_precursor_mz", "adduct", "source_hash",
        "source_row_index", "provenance", "external_accession",
        "diagnostic_ions", "structure_identity_verified",
    )
    coverage = {
        field: {
            "present": int(counts[f"present_{field}"]),
            "missing": int(total - counts[f"present_{field}"]),
            "fraction": _ratio(counts[f"present_{field}"], total),
        }
        for field in coverage_fields
    }

    runtime_quality = target_catalog_quality_summary(runtime_targets)

    def runtime_key(target: Mapping[str, Any]) -> str:
        return str(
            target.get("catalog_row_id")
            or target.get("target_id")
            or ""
        )

    active_runtime_keys = {
        runtime_key(target)
        for target in runtime_targets
        if target.get("runtime_active") is not False and runtime_key(target)
    }
    invalid_adduct_keys = {
        runtime_key(target)
        for target in runtime_targets
        if target.get("runtime_active") is not False
        and target.get("adduct_validation_status") == "quarantined"
        and runtime_key(target)
    }
    execution: dict[str, Any] = {}
    represented_any: set[str] = set()
    conflict_any: set[str] = set()
    for polarity in ("POS", "NEG"):
        hypotheses = build_expected_ion_hypotheses(
            runtime_targets,
            polarity,
            include_safe_common_adducts=True,
            include_reported_precursor_conflicts=False,
        )
        represented = {
            str(member.get("catalog_row_id") or member.get("target_id") or "")
            for hypothesis in hypotheses
            for member in hypothesis.get("members") or []
            if str(member.get("target_id") or member.get("catalog_row_id") or "")
        }
        conflict_members = {
            str(member.get("catalog_row_id") or member.get("target_id") or "")
            for hypothesis in hypotheses
            if "reported_precursor_conflicts_with_theoretical_ion"
            in (hypothesis.get("ion_quality_flags") or [])
            for member in hypothesis.get("members") or []
            if str(member.get("target_id") or member.get("catalog_row_id") or "")
        }
        represented_any.update(represented)
        conflict_any.update(conflict_members)
        execution[polarity] = {
            "executable_ion_groups": len(hypotheses),
            "active_catalog_rows_represented": len(represented),
            "active_catalog_rows_not_represented": len(
                active_runtime_keys - represented
            ),
            "formula_theoretical_primary_groups": sum(
                bool(item.get("formula_theoretical_primary"))
                for item in hypotheses
            ),
            "reported_precursor_only_groups": sum(
                bool(item.get("reported_precursor_priority"))
                for item in hypotheses
            ),
            "source_precursor_conflict_groups_using_theoretical_mass": sum(
                "reported_precursor_conflicts_with_theoretical_ion"
                in (item.get("ion_quality_flags") or [])
                for item in hypotheses
            ),
            "source_precursor_conflict_hypotheses_executed": 0,
        }

    runtime_execution = {
        **runtime_quality,
        "safe_common_adduct_policy": {
            "POS": ["[M+H]+", "[M+Na]+", "[M+NH4]+"],
            "NEG": ["[M-H]-"],
            "formate_or_acetate_defaulted_when_unreported": False,
        },
        "polarity_execution": execution,
        "active_catalog_rows_represented_in_any_polarity": len(
            active_runtime_keys & represented_any
        ),
        "active_catalog_rows_without_any_executable_ion": len(
            active_runtime_keys - represented_any
        ),
        "invalid_reported_adduct_rows_recovered_by_safe_formula_ions": len(
            invalid_adduct_keys & represented_any
        ),
        "source_precursor_conflict_rows_using_theoretical_mass": len(
            conflict_any
        ),
        "reported_precursor_conflict_hypotheses_executed_by_default": 0,
    }
    return {
        "formula_id": formula_id,
        "formula_name": FORMULA_NAMES[formula_id],
        "dataset_source": dataset,
        "catalog_semantics": "literature_presence_targets_not_spectral_library_identifications",
        "row_fingerprint_sha256": fingerprint.hexdigest(),
        "row_count": int(total),
        "coverage": coverage,
        "runtime_execution": runtime_execution,
        "metrics": {key: int(value) for key, value in sorted(counts.items())},
        "source_tiers": dict(source_tiers.most_common()),
        "target_searchability_labels": dict(searchability_labels.most_common()),
        "canonical_adducts": dict(canonical_adducts.most_common()),
        "quarantined_adduct_values_top20": dict(
            quarantined_adduct_values.most_common(20)
        ),
        "quality_flags_top30": dict(quality_flags.most_common(30)),
        "issue_summary": {
            "critical_row_issue_events": int(critical),
            "high_priority_issue_events": int(high),
            "runtime_searchable_fraction": _ratio(
                counts["runtime_searchable"], total
            ),
            "interpretation": (
                "Counts are issue events and can overlap within a row; they are "
                "not counts of confirmed false identifications."
            ),
        },
        "examples": dict(examples),
    }


def _markdown(payload: Mapping[str, Any]) -> str:
    lines = [
        "# MSsystem V12 八个复方靶向数据库审计",
        "",
        f"生成时间（UTC）：{payload['generated_at_utc']}",
        "",
        "本审计直接读取 SQLite `TargetCompounds`/`DiagnosticIons`，全程只读。"
        "这些行是文献存在性靶标，不是实验谱库，也不是样本中已经确认的化合物。",
        "",
        "## 总览",
        "",
        "| 方剂 | 数据集 | 源库行 | 进入计算 | 运行时隔离 | POS 离子组 | NEG 离子组 | 任一极性可执行行 | 无可执行离子的有效行 | 异常加合物安全恢复 |",
        "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for item in payload["datasets"]:
        metrics = item["metrics"]
        rows = item["row_count"]
        runtime = item["runtime_execution"]
        polarity = runtime["polarity_execution"]

        def pct(value: int) -> str:
            return f"{value}/{rows} ({value / rows * 100:.1f}%)" if rows else "0"

        lines.append(
            "| {formula_id} {formula_name} | `{dataset_source}` | {rows} | {active} | "
            "{quarantined} | {pos} | {neg} | {represented} | {unexecutable} | "
            "{adduct_recovered} |".format(
                formula_id=item["formula_id"],
                formula_name=item["formula_name"],
                dataset_source=item["dataset_source"],
                rows=rows,
                active=pct(runtime["active_target_rows"]),
                quarantined=runtime["quarantined_target_rows"],
                pos=polarity["POS"]["executable_ion_groups"],
                neg=polarity["NEG"]["executable_ion_groups"],
                represented=pct(
                    runtime["active_catalog_rows_represented_in_any_polarity"]
                ),
                unexecutable=runtime[
                    "active_catalog_rows_without_any_executable_ion"
                ],
                adduct_recovered=runtime[
                    "invalid_reported_adduct_rows_recovered_by_safe_formula_ions"
                ],
            )
        )

    lines.extend([
        "",
        "## 已实施的安全修复",
        "",
        "- 加合物现在使用严格白名单规范化。m/z 前缀、残缺括号、多个离子形式拼接、"
        "中性丢失备注和二聚体文本不会再作为主动提取假设。原始文本仍保留用于审计。",
        "- 异常加合物被隔离后，只要有可解析分子式或精确质量，POS 仅生成 "
        "`[M+H]+/[M+Na]+/[M+NH4]+`，NEG 仅生成 `[M-H]-`；没有来源依据时不默认"
        "甲酸根或乙酸根。",
        "- 分子式理论离子为主；来源前体与理论值相差超过 20 ppm 时，默认不执行冲突的"
        "来源 m/z，仅保留理论离子和审计记录。",
        "- 未对 5 GB 共享数据库做破坏性批量改写。源数据缺失、结构冲突和异常字符串均在"
        "JSON 明细中保留，避免丢失原始谱系。",
        "",
        "## 结论边界",
        "",
        "运行时可检索不等同于化合物身份正确；MS1 精确质量只能形成候选。异构体级或具体"
        "连接结构结论仍需同方法标准品 RT/MS/MS 或经验证的实验参考谱。结构已验证标志与"
        "SMILES/分子式/InChIKey 冲突的行应优先回源校正，不能由软件静默猜测。",
        "",
        "完整逐库指标、异常示例、数据集指纹和质量标记见同名 JSON。",
        "",
    ])
    return "\n".join(lines)


def run(database: Path, output_dir: Path) -> dict[str, Any]:
    if not database.is_file():
        raise FileNotFoundError(database)
    output_dir.mkdir(parents=True, exist_ok=True)
    assess_structure = _structure_assessor()
    with _open_read_only(database) as connection:
        available = {
            str(row[0])
            for row in connection.execute(
                "SELECT DISTINCT dataset_source FROM TargetCompounds"
            )
        }
        missing = [
            dataset for dataset in FORMULA_DATASETS.values()
            if dataset not in available
        ]
        if missing:
            raise RuntimeError(f"Missing targeted datasets: {missing}")
        datasets = [
            _audit_dataset(
                connection,
                formula_id,
                dataset,
                assess_structure,
                load_expected_compound_catalog(
                    database, dataset_sources=[dataset]
                ),
            )
            for formula_id, dataset in FORMULA_DATASETS.items()
        ]
        schema_version = int(
            connection.execute("PRAGMA schema_version").fetchone()[0]
        )

    generated = datetime.now(timezone.utc).replace(microsecond=0).isoformat()
    payload = {
        "schema_version": "1.1.0",
        "audit_id": "mssystem-v12-eight-formula-runtime-audit-20260820",
        "generated_at_utc": generated,
        "read_only": True,
        "database": {
            "path": str(database.resolve()),
            "size_bytes": database.stat().st_size,
            "mtime_ns": database.stat().st_mtime_ns,
            "sqlite_schema_version": schema_version,
            "full_database_hash_omitted_reason": (
                "The shared database is approximately 5 GB; every audited row "
                "has instead been streamed into a per-dataset SHA-256 fingerprint."
            ),
        },
        "runtime_remediation": {
            "invalid_adduct_policy": "quarantine_raw_text_do_not_execute",
            "safe_formula_adduct_fallback": True,
            "reported_precursor_conflict_threshold_ppm": 20.0,
            "database_rows_rewritten": False,
        },
        "datasets": datasets,
    }
    stem = "MSSYSTEM_V12_EIGHT_FORMULA_RUNTIME_AUDIT_20260820"
    json_path = output_dir / f"{stem}.json"
    md_path = output_dir / f"{stem}.md"
    json_path.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    md_path.write_text(_markdown(payload), encoding="utf-8")
    return {"json": str(json_path), "markdown": str(md_path), "payload": payload}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--database", type=Path, default=DEFAULT_DATABASE)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    result = run(args.database, args.output_dir)
    print(json.dumps({key: value for key, value in result.items() if key != "payload"}, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
