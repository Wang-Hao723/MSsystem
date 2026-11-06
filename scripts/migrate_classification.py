#!/usr/bin/env python3
"""就地迁移质谱库: 清洗被污染的 classification + 新增 dataset_source 列。

背景: build_from_mgf_msp.py 旧版用 file_path.stem 当 classification, 导致 86% 的
Compounds.classification 是数据集文件名 (MoNA-export-All_Spectra 等)。本脚本:
  1. Spectra 加 dataset_source 列 = normalize_source_file(source_file)  (数据集名[-极性])
  2. Compounds.classification 重算 = resolve_compound_class(smiles, raw, source_file)
     (SMARTS 结构规则 → 文件名归一化 → 原 raw 若非数据集名 → 'other')

安全: 先 cp 整库到 .bak; 幂等(检查列是否存在 + WHERE 条件); 失败可从 .bak 恢复。

Usage:
    python scripts/migrate_classification.py                # 迁移 data/ 下两个库
    python scripts/migrate_classification.py --dry-run      # 只统计不写
    python scripts/migrate_classification.py --db data/local_library.db
"""
from __future__ import annotations

import argparse
import shutil
import sqlite3
import sys
import time
from collections import Counter
from pathlib import Path

# 让脚本能 import 项目模块
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from mssystem.classification import (  # noqa: E402
    resolve_compound_class,
    normalize_source_file,
    is_source_or_generic,
)


def _column_exists(conn: sqlite3.Connection, table: str, col: str) -> bool:
    rows = conn.execute(f"PRAGMA table_info({table})").fetchall()
    return any(r[1] == col for r in rows)


def _backup(db_path: Path) -> Path:
    bak = db_path.with_suffix(db_path.suffix + ".bak")
    if not bak.exists():
        shutil.copy2(db_path, bak)
        print(f"  ✓ 备份 → {bak.name}")
    else:
        print(f"  · 备份已存在, 跳过 ({bak.name})")
    return bak


def _compound_source_map(conn: sqlite3.Connection) -> dict[int, str]:
    """每个 compound 的众数 source_file (一条 GROUP BY 查询)。"""
    cur = conn.execute(
        """
        SELECT mp.compound_id, sp.source_file, COUNT(*) c
        FROM Mapping mp
        JOIN Spectra sp ON sp.id = mp.spectrum_id
        WHERE sp.source_file IS NOT NULL AND sp.source_file != ''
        GROUP BY mp.compound_id, sp.source_file
        """
    )
    best: dict[int, tuple[int, str]] = {}
    for cid, sf, cnt in cur:
        prev = best.get(cid)
        if prev is None or cnt > prev[0]:
            best[cid] = (cnt, sf)
    return {cid: sf for cid, (_, sf) in best.items()}


def migrate(db_path: Path, dry_run: bool = False) -> None:
    print(f"\n{'='*60}\n处理: {db_path}")
    if not db_path.exists():
        print(f"  ✗ 不存在, 跳过")
        return

    if not dry_run:
        _backup(db_path)

    conn = sqlite3.connect(str(db_path))
    conn.row_factory = sqlite3.Row

    # ── 迁移前统计 ──
    total = conn.execute("SELECT COUNT(*) FROM Compounds").fetchone()[0]
    polluted = conn.execute(
        "SELECT COUNT(*) FROM Compounds WHERE classification IS NOT NULL AND "
        "(classification LIKE 'MoNA%' OR classification LIKE 'MSMS%' OR "
        "classification LIKE 'Spectraverse%' OR classification LIKE 'MassBank%' OR "
        "classification LIKE 'GNPS%' OR classification LIKE 'MoNA_export%')"
    ).fetchone()[0]
    print(f"  迁移前: {total:,} 化合物, 其中 ~{polluted:,} ({100*polluted/max(total,1):.1f}%) classification 是数据集名")

    # ── Step 1: Spectra.dataset_source ──
    has_col = _column_exists(conn, "Spectra", "dataset_source")
    if not has_col:
        if dry_run:
            print("  [dry-run] 将 ALTER TABLE Spectra ADD COLUMN dataset_source TEXT")
        else:
            conn.execute("ALTER TABLE Spectra ADD COLUMN dataset_source TEXT")
            print("  ✓ ALTER TABLE Spectra ADD COLUMN dataset_source")
    else:
        print("  · Spectra.dataset_source 列已存在")

    if not dry_run:
        t0 = time.time()
        # SQLite 不能在 UPDATE 里直接调 Python 函数 → 取出再批量写
        rows = conn.execute(
            "SELECT id, source_file FROM Spectra WHERE source_file IS NOT NULL "
            "AND (dataset_source IS NULL OR dataset_source = '')"
        ).fetchall()
        conn.executemany(
            "UPDATE Spectra SET dataset_source = ? WHERE id = ?",
            [(normalize_source_file(r["source_file"]), r["id"]) for r in rows],
        )
        conn.commit()
        print(f"  ✓ Spectra.dataset_source backfill: {len(rows):,} 行 ({time.time()-t0:.1f}s)")

    # ── Step 2: Compounds.classification 重算 ──
    print("  ⏳ 重算 Compounds.classification (SMARTS+归一化) ...")
    src_map = _compound_source_map(conn) if not dry_run else {}
    t0 = time.time()
    rows = conn.execute(
        "SELECT id, smiles, classification FROM Compounds"
    ).fetchall()
    updates = []
    dist = Counter()
    via = Counter()
    from mssystem.classification import _get_rule_clf, normalize_filename_class
    rule_clf = _get_rule_clf()
    for r in rows:
        cid = r["id"]
        smi = r["smiles"]
        raw = r["classification"]
        src = src_map.get(cid)
        new_cls = resolve_compound_class(smiles=smi, raw=raw, source_file=src)
        dist[new_cls] += 1
        # 判类来源统计 (与 resolve 优先级一致)
        if smi and rule_clf.predict(smiles=smi) not in ("Unknown", None):
            via["SMARTS结构"] += 1
        elif normalize_filename_class(src):
            via["文件名"] += 1
        elif raw and normalize_filename_class(raw):
            via["raw归一化"] += 1
        else:
            via["other"] += 1
        if not dry_run and new_cls != (raw or ""):
            updates.append((new_cls, cid))
        if not dry_run and len(updates) >= 5000:
            conn.executemany("UPDATE Compounds SET classification = ? WHERE id = ?", updates)
            conn.commit()
            updates.clear()
            print(f"    ... {cid:,}/{total:,}", end="\r")
    if not dry_run and updates:
        conn.executemany("UPDATE Compounds SET classification = ? WHERE id = ?", updates)
        conn.commit()
    print(f"  ✓ classification 重算完成 ({time.time()-t0:.1f}s)            ")

    # ── 迁移后统计 ──
    oth = dist.get("other", 0)
    print(f"\n  迁移后 classification 分布 (clean-class 覆盖 {total-oth:,}/{total:,} = {100*(total-oth)/max(total,1):.1f}%):")
    for k, v in dist.most_common(16):
        print(f"    {k:18s} {v:>8,}  ({100*v/max(total,1):.1f}%)")
    print(f"  判类来源: {dict(via)}")
    if not dry_run and _column_exists(conn, "Spectra", "dataset_source"):
        ds_dist = conn.execute(
            "SELECT COUNT(DISTINCT dataset_source) n FROM Spectra WHERE dataset_source IS NOT NULL AND dataset_source!=''"
        ).fetchone()[0]
        print(f"  dataset_source distinct 数: {ds_dist}")

    conn.close()
    print(f"  {'[dry-run] 未写入' if dry_run else '✓ 已写入'} {db_path.name}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--db", action="append", default=None, help="指定单个库路径(可重复); 默认两个都跑")
    ap.add_argument("--dry-run", action="store_true", help="只统计, 不写库")
    args = ap.parse_args()

    project = Path(__file__).resolve().parent.parent
    dbs = [project / p for p in (args.db or ["data/local_library.db", "data/local_library_small.db"])]
    for db in dbs:
        migrate(db, dry_run=args.dry_run)
    print("\n完成。出错可从 .bak 恢复: cp <db>.bak <db>")


if __name__ == "__main__":
    main()
