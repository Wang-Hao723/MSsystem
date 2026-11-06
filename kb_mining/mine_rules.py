"""
mssystem v3 Phase 3.1 — 经验挖掘裂解规则 (数据驱动, 全数值化)
===============================================================
从 v3 DB 的真实谱图(按 classification_v3 分类)统计:
  - 诊断离子 (diagnostic ions): 高频碎片 m/z (按 ±0.025Da bin, 跨谱计数)
  - 中性丢失 (neutral losses): 高频 (precursor - fragment) 质量差

每条规则带**数值 m/z/loss + 出现频率 + 支持谱图数**, 天然可机器匹配。
根治 v2 KB "符号 m/z 无法匹配真实谱" 的致命缺陷。

输出: kb_v3/empirical_rules.json  {class: {pos/neg: {diagnostic_ions, neutral_losses, n_spectra}}}

cd /Users/zizi/Desktop/Mssystem/mssystem_v3
python3 kb_mining/mine_rules.py            # small DB
python3 kb_mining/mine_rules.py --full     # full DB
"""
from __future__ import annotations

import argparse
import json
import sqlite3
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

BIN = 0.025           # m/z 分箱 (Da) — 诊断离子聚类粒度
MIN_SPECTRA = 50      # 类-极性 最少谱图数才挖
FREQ_MIN = 0.15       # 频率阈值: 出现在 ≥15% 该类谱图中才算规律
TOPK = 25             # 每类每极性 最多保留 Top-K 规则
MZ_MIN, MZ_MAX = 30.0, 1000.0
LOSS_MIN, LOSS_MAX = 1.0, 400.0   # 中性丢失质量范围
OUT = ROOT / "kb_v3" / "empirical_rules.json"


def polarity_of(adduct: str, charge) -> str:
    """从 adduct/charge 推极性。"""
    a = (adduct or "").strip()
    if a:
        if "+" in a and "-" not in a:
            return "pos"
        if "-" in a and "+" not in a:
            return "neg"
    try:
        c = int(charge) if charge not in (None, "") else 0
        if c > 0:
            return "pos"
        if c < 0:
            return "neg"
    except Exception:
        pass
    return "unknown"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--full", action="store_true")
    args = ap.parse_args()
    db = str(ROOT / "data" / ("local_library.db" if args.full else "local_library_small.db"))
    OUT.parent.mkdir(parents=True, exist_ok=True)

    print(f"=== Phase 3.1 经验挖掘 ({'full' if args.full else 'small'} DB) ===")
    con = sqlite3.connect(db); con.row_factory = sqlite3.Row
    cur = con.cursor()
    # 检查 classification_v3 列
    has_v3 = any(r[1] == "classification_v3" for r in cur.execute("PRAGMA table_info(Compounds)").fetchall())
    cls_col = "c.classification_v3" if has_v3 else "c.classification"
    print(f"  用分类列: {cls_col}")

    rows = cur.execute(f"""
        SELECT s.precursor_mz, s.peaks_json, s.adduct, s.charge, {cls_col} AS cls
        FROM Spectra s JOIN Mapping m ON m.spectrum_id=s.id
        JOIN Compounds c ON m.compound_id=c.id
        WHERE {cls_col} IS NOT NULL AND {cls_col}!='' AND {cls_col}!='other'
    """).fetchall()
    print(f"  谱图(有分类): {len(rows):,}")

    # 累积: diag[(cls,pol)][bin] = set(spectrum_idx);  loss[(cls,pol)][bin] = set
    # 为控内存, 用 dict[int,int] 计数 (每个 bin 一个 counter, 每谱只 +1 一次)
    diag = defaultdict(lambda: defaultdict(int))   # (cls,pol) -> bin -> n_spectra_containing
    loss = defaultdict(lambda: defaultdict(int))
    n_spec = defaultdict(int)                       # (cls,pol) -> total spectra
    n_skipped = 0

    for idx, r in enumerate(rows):
        cls = r["cls"]
        if not cls:
            continue
        try:
            pmz = float(r["precursor_mz"])
            peaks = json.loads(r["peaks_json"] or "[]")
        except Exception:
            n_skipped += 1; continue
        if not peaks or pmz <= 0:
            continue
        pol = polarity_of(r["adduct"], r["charge"])
        if pol == "unknown":
            continue
        key = (cls, pol)
        n_spec[key] += 1
        seen_diag = set()
        seen_loss = set()
        for p in peaks:
            try:
                mz, inten = float(p[0]), float(p[1])
            except Exception:
                continue
            if inten <= 0:
                continue
            if MZ_MIN <= mz <= MZ_MAX:
                seen_diag.add(int(round(mz / BIN)))
            d = pmz - mz
            if LOSS_MIN <= d <= LOSS_MAX:
                seen_loss.add(int(round(d / BIN)))
        dmap = diag[key]; lmap = loss[key]
        for b in seen_diag:
            dmap[b] += 1
        for b in seen_loss:
            lmap[b] += 1
        if idx % 80000 == 0 and idx:
            print(f"    ... {idx:,} 谱处理", flush=True)

    print(f"  跳过(解析失败): {n_skipped}; 类×极性 组合: {len(n_spec)}")

    # 汇总规则
    result = {}
    for (cls, pol), ntot in sorted(n_spec.items()):
        if ntot < MIN_SPECTRA:
            continue
        thr = max(FREQ_MIN * ntot, 3)  # 至少 3 谱支持
        # 诊断离子: bin -> (center_mz, freq, n)
        di = sorted(
            ((b * BIN + BIN / 2, cnt / ntot, cnt) for b, cnt in diag[(cls, pol)].items() if cnt >= thr),
            key=lambda x: -x[1])[:TOPK]
        nl = sorted(
            ((b * BIN + BIN / 2, cnt / ntot, cnt) for b, cnt in loss[(cls, pol)].items() if cnt >= thr),
            key=lambda x: -x[1])[:TOPK]
        result.setdefault(cls, {})[pol] = {
            "n_spectra": ntot,
            "diagnostic_ions": [{"mz": round(m, 4), "freq": round(f, 3), "n_support": int(c)} for m, f, c in di],
            "neutral_losses": [{"loss_da": round(m, 4), "freq": round(f, 3), "n_support": int(c)} for m, f, c in nl],
        }

    # 常见中性丢失注释 (用精确质量匹配常见生物丢失)
    NL_ANNOT = {
        18.0106: "H2O", 28.0101: "CO", 35.9768: "HCl", 43.9898: "CO2",
        46.0055: "HCOOH", 57.0215: "CHNCO/HCNO", 60.0211: "CH2CO+OH?",
        71.0371: "C3H3O2?", 79.9568: "SO3", 87.0320: "C4H7NO?",
        90.0317: "(脱氧)糖? ", 120.0423: "C7H4O2?",
        146.0579: "脱氧己糖 deoxyhexose", 162.0528: "己糖 hexose",
        176.0321: "葡萄糖醛酸 glucuronic acid", 180.0634: "己糖+H2O?",
        246.0618: "乙酰己糖?", 308.0905: "芸香糖 rutinose",
    }
    def annot(loss_da):
        for k, v in NL_ANNOT.items():
            if abs(loss_da - k) < 0.03:
                return v
        return None

    # 补注释
    for cls, pols in result.items():
        for pol, d in pols.items():
            for nl in d["neutral_losses"]:
                nl["annotation"] = annot(nl["loss_da"])

    with open(OUT, "w", encoding="utf-8") as f:
        json.dump(result, f, ensure_ascii=False, indent=2)

    # 摘要
    print(f"\n  挖掘完成, {len(result)} 类有规律 → {OUT}")
    for cls in sorted(result):
        pols = result[cls]
        parts = []
        for pol in ("pos", "neg"):
            if pol in pols:
                d = pols[pol]
                top_nl = ", ".join(f"{nl['loss_da']}({nl['annotation'] or '?'})" for nl in d["neutral_losses"][:3])
                parts.append(f"{pol}:{d['n_spectra']}谱/{len(d['diagnostic_ions'])}诊断离子/{len(d['neutral_losses'])}丢失[top:{top_nl}]")
        print(f"    {cls:16s} " + " | ".join(parts))


if __name__ == "__main__":
    main()
