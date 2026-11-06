"""
mssystem v3 — P2·2b 逐化合物裂解路径自动挖掘
================================================
遍历库内**每个化合物自己的馆藏谱图**(Spectra.peaks_json + precursor_mz + adduct/charge),
构建【前体 → 中性丢失链 → 诊断离子】的碎裂路径, 序列化为 JSON 存入
Compounds 新增列 fragmentation_pathway (纯增量, 不动任何匹配/评分代码)。

与 kb_mining/mine_rules.py 的区别:
  - mine_rules: 按类跨化合物统计 (类级模板)
  - 本脚本:     按化合物跨它自己的多张谱聚合 (化合物级路径) ← 把 "0% 路径标注" 填满

化学一致性: 复用 matching.fragmentation_rules 的
  COMMON_NEUTRAL_LOSSES (数值质量 + 官能团 + 注释) 与 DIAGNOSTIC_IONS (精确 m/z + 来源)
  做匹配标注; "复现" = (a) 命中上述策展表 (化学上已知的丢失/诊断离子), 或
  (b) 该化合物 ≥2 张谱中相同 bin 重复出现 (数据驱动的新丢失/诊断离子)。

资源安全: 按 compound_id 流式 (cursor ORDER BY id, 不 fetchall),
  按 BATCH_COMPOUNDS 批写; 续跑点存 data/pathway_mining_progress.json (按 last_id 续算)。
  peaks_json 逐条 try/except, 统计跳过数。

用法:
  cd /Users/zizi/Desktop/Mssystem/mssystem_v3
  # 原型: 小库 P2 副本前 1000 个化合物
  python3 scripts/mine_pathways.py --db data/local_library_p2.db --limit 1000
  # 小库全量 (默认 P2 副本)
  python3 scripts/mine_pathways.py --db data/local_library_p2.db
  # 全库 (先 cp data/local_library.db data/local_library_p2_full.db)
  python3 scripts/mine_pathways.py --db data/local_library_p2_full.db
  # 续跑 (默认即续跑; --reset 清空进度从头来)
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sqlite3
import sys
import time
from collections import defaultdict
from pathlib import Path
from statistics import median

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from matching.fragmentation_rules import (  # noqa: E402
    COMMON_NEUTRAL_LOSSES, DIAGNOSTIC_IONS,
    normalize_classification, _is_diagnostic_ion_compatible,
)
from matching.classification_provenance import (  # noqa: E402
    classification_source_trust,
)

# ──────────────────────────────────────────────────────────────────────
# 参数 (与 kb_mining/mine_rules.py 一致的 bin 粒度, 便于跨工具可比)
# ──────────────────────────────────────────────────────────────────────
BIN = 0.025                 # m/z / loss 分箱 (Da)
MZ_MIN, MZ_MAX = 30.0, 1200.0          # 诊断离子 m/z 范围
LOSS_MIN, LOSS_MAX = 1.0, 600.0        # 中性丢失质量范围
NL_TOL = 0.05               # 命中策展中性丢失/诊断离子的容差 (Da) — 混合分辨率库, 稳健取 0.05
INT_FLOOR_PCT = 1.0         # 链构建时只考虑相对强度 ≥1% base peak 的产物离子
TOPK_LOSS = 40              # 每化合物最多保留中性丢失条目 (完整粒度上限)
TOPK_DIAG = 40              # 每化合物最多保留诊断离子条目
TOPK_CHAIN = 5             # 最多保留候选丢失链
MAX_CHAIN_DEPTH = 4         # 链最大深度 (前体 → 最多 4 次丢失)
MAX_CHAIN_NODES = 256       # DFS 节点预算 (防爆)
SCHEMA_VERSION = 2
BATCH_COMPOUNDS = 2000      # 批写化合物数

PROGRESS_DEFAULT = ROOT / "data" / "pathway_mining_progress.json"

# 策展表 → 按质量排序的查找结构 (一次性构建)
_NL_LIST = sorted(
    ((nl.mass, key, nl.formula, nl.name_en, nl.functional_groups)
     for key, nl in COMMON_NEUTRAL_LOSSES.items()),
    key=lambda x: x[0])
_NL_MASSES = [x[0] for x in _NL_LIST]
_DIAG_LIST = sorted(
    ((di.mz, key, di.formula, di.name, di.polarity, di.functional_group)
     for key, di in DIAGNOSTIC_IONS.items()),
    key=lambda x: x[0])
_DIAG_MZS = [x[0] for x in _DIAG_LIST]


def _bisect_match(values, target, tol):
    """在有序 values 中找 |v-target|<tol 的第一个, 返回 index 或 None。"""
    import bisect
    lo = bisect.bisect_left(values, target - tol)
    hi = bisect.bisect_right(values, target + tol)
    if lo >= hi:
        return None
    # 取最接近的
    best, best_d = None, tol
    for i in range(lo, hi):
        d = abs(values[i] - target)
        if d <= best_d:
            best, best_d = i, d
    return best


def annotate_neutral_loss(loss_da: float) -> dict | None:
    i = _bisect_match(_NL_MASSES, loss_da, NL_TOL)
    if i is None:
        return None
    _mass, key, formula, name_en, fgs = _NL_LIST[i]
    return {"key": key, "formula": formula, "name_en": name_en}


def annotate_diagnostic_ion(mz: float, classification: str = "") -> dict | None:
    """诊断离子标注 + **分类门控**: 只接受与化合物分类化学兼容的离子。
    (诊断离子是类特异结构报告离子; 不门控会把随机峰误标成 '黄酮A环' 等。)
    用 bisect 预筛质量窗口内的候选(通常 0-2 个), 再逐个过 _is_diagnostic_ion_compatible。"""
    import bisect
    eff = normalize_classification(classification or "")
    lo = bisect.bisect_left(_DIAG_MZS, mz - NL_TOL)
    hi = bisect.bisect_right(_DIAG_MZS, mz + NL_TOL)
    best, best_d = None, NL_TOL
    for i in range(lo, hi):
        _mz, key, formula, name, pol, fg = _DIAG_LIST[i]
        d = abs(_mz - mz)
        if d <= best_d and _is_diagnostic_ion_compatible(DIAGNOSTIC_IONS[key], eff):
            best, best_d = (key, formula, name, pol, fg), d
    if best is None:
        return None
    key, formula, name, pol, fg = best
    return {"key": key, "formula": formula, "name": name, "polarity": pol, "fg": fg}


def polarity_of(adduct, charge) -> str:
    """稳健极性推断 (库 adduct 字段脏: '[M+H]+' / 'M+H' / 'Positive' / 'negative' 混存)。"""
    a = (adduct or "").strip()
    if a:
        if "+" in a and "-" not in a:
            return "pos"
        if "-" in a and "+" not in a:
            return "neg"
        low = a.lower()
        if low.startswith("pos"):
            return "pos"
        if low.startswith("neg"):
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


def _norm_peaks(peaks_raw):
    """解析 peaks_json 列表 → [(mz, rel_pct)], 相对 base peak 归一, 丢零强度/非法。"""
    out = []
    if not peaks_raw:
        return out
    max_i = 0.0
    parsed = []
    for p in peaks_raw:
        try:
            mz, it = float(p[0]), abs(float(p[1]))
        except Exception:
            continue
        if mz <= 0 or it <= 0:
            continue
        parsed.append((mz, it))
        if it > max_i:
            max_i = it
    if max_i <= 0:
        return out
    for mz, it in parsed:
        out.append((mz, 100.0 * it / max_i))
    return out


def mine_compound(
    cid: int,
    spectra: list[dict],
    classification: str = "",
    classification_source: str | None = None,
) -> dict | None:
    """对单化合物的多张谱构建 fragmentation_pathway dict。

    spectra: [{precursor_mz, peaks:[(mz,rel_pct)], adduct, charge}, ...] (已解析)
    classification: classification_v3。仅当 classification_source 明确来自
      策展/结构规则时才可用于诊断离子门控；ML/未知来源 fail closed。
    返回 pathway dict 或 None (无可用谱)。

    设计:
      - 中性丢失 = 宽指纹 (质量差稳健): 按 (复现,强度) 取 Top-K, 策展标注(类无关)。
      - 诊断离子 = 紧标注 (类特异): 仅保留 "类门控策展命中" 或 "跨谱复现≥2"。
    """
    valid = [s for s in spectra if s["precursor_mz"] and s["precursor_mz"] > 0 and s["peaks"]]
    if not valid:
        return None
    n_spec = len(valid)
    pmz_rep = median(s["precursor_mz"] for s in valid)   # 代表前体 (跨谱中位, 抗离群)
    pols = sorted({s["polarity"] for s in valid if s["polarity"] != "unknown"})
    adducts = sorted({s["adduct"] for s in valid if s["adduct"]})
    n_peaks_total = sum(len(s["peaks"]) for s in valid)
    class_trusted, class_trust_reason = classification_source_trust(
        classification_source,
        allow_legacy_missing=False,
    )
    effective_classification = classification if class_trusted else ""

    # 跨谱聚合: 每个 bin 在多少张谱里出现 (复现计数用 bin),
    # 但**策展标注用该 bin 最强峰的真实 m/z/loss 值** (bin 中心会丢精度, 误丢边界 H2O 等)。
    diag_cnt = defaultdict(int)      # bin -> n_spectra_containing
    loss_cnt = defaultdict(int)
    diag_inten = defaultdict(float)  # bin -> max rel% (排序用)
    loss_inten = defaultdict(float)
    diag_real = {}                   # bin -> 该 bin 最强峰的真实 m/z
    loss_real = {}                   # lb -> 该 lb 最强峰的真实 loss (pmz - mz)
    # 链构建用的产物离子集合 (跨谱并集, 强度≥floor): bin -> (center_mz, max rel%)
    prod_bin = {}

    for s in valid:
        pmz = s["precursor_mz"]
        seen_d, seen_l = set(), set()
        for mz, rel in s["peaks"]:
            b = int(round(mz / BIN))
            if MZ_MIN <= mz <= MZ_MAX:
                if b not in seen_d:
                    diag_cnt[b] += 1
                    seen_d.add(b)
                if rel >= diag_inten[b]:
                    diag_inten[b] = rel
                    diag_real[b] = mz
            d = pmz - mz
            if LOSS_MIN <= d <= LOSS_MAX:
                lb = int(round(d / BIN))
                if lb not in seen_l:
                    loss_cnt[lb] += 1
                    seen_l.add(lb)
                if rel >= loss_inten[lb]:
                    loss_inten[lb] = rel
                    loss_real[lb] = d
            # 产物离子 (链校验用)
            if rel >= INT_FLOOR_PCT and MZ_MIN <= mz <= MZ_MAX:
                cur = prod_bin.get(b)
                if cur is None or rel > cur[1]:
                    prod_bin[b] = (b * BIN + BIN / 2, rel)

    # ── 中性丢失条目 (策展命中 + 跨谱复现新丢失; 类无关标注) ──
    # 不收录 "单谱任意质量差" — 那是峰列表噪声, 非裂解路径。
    losses = []
    for lb, cnt in loss_cnt.items():
        real = loss_real.get(lb, lb * BIN + BIN / 2)
        ann = annotate_neutral_loss(real)         # 真实 loss 标注, 抗 bin 精度损失
        recurring = cnt >= 2                       # ≥2 张谱复现
        if not ann and not recurring:
            continue                  # 单谱且非策展 → 视为噪声丢弃
        frag_mz = pmz_rep - real
        losses.append({
            "loss_da": round(real, 4),
            "freq": round(cnt / n_spec, 3),
            "n_spectra": int(cnt),
            "fragment_mz": round(frag_mz, 4) if frag_mz > 0 else None,
            "rel_intensity": round(loss_inten[lb], 1),
            "annotation": ann["name_en"] if ann else None,
            "nl_key": ann["key"] if ann else None,
            "recurring": bool(recurring),
        })
    losses.sort(key=lambda x: (not x["recurring"], -(x["rel_intensity"]), -x["freq"]))
    losses = losses[:TOPK_LOSS]

    # ── 诊断离子条目 (紧标注: 类特异, 须类门控命中或跨谱复现) ──
    # 单谱且无类兼容策展命中 → 丢弃 (避免随机峰误标 '黄酮A环' 等)。
    diags = []
    for b, cnt in diag_cnt.items():
        real = diag_real.get(b, b * BIN + BIN / 2)
        ann = annotate_diagnostic_ion(real, effective_classification)
        recurring = cnt >= 2
        if not ann and not recurring:
            continue
        diags.append({
            "mz": round(real, 4),
            "freq": round(cnt / n_spec, 3),
            "n_spectra": int(cnt),
            "rel_intensity": round(diag_inten[b], 1),
            "annotation": ann["name"] if ann else None,
            "diag_key": ann["key"] if ann else None,
            "recurring": bool(recurring),
        })
    diags.sort(key=lambda x: (not x["recurring"], -(x["rel_intensity"]), -x["freq"]))
    diags = diags[:TOPK_DIAG]

    # ── 候选丢失链: 前体逐次减去策展丢失, 每步产物须为观察峰 ──
    chains = _build_chains(pmz_rep, prod_bin)
    source_ids = sorted({
        str(s.get("spectrum_id") or "").strip()
        for s in valid
        if str(s.get("spectrum_id") or "").strip()
    })
    source_datasets = sorted({
        str(s.get("dataset_source") or "").strip()
        for s in valid
        if str(s.get("dataset_source") or "").strip()
    })
    spectrum_set_sha256 = hashlib.sha256(
        "\n".join(source_ids).encode("utf-8")
    ).hexdigest()

    return {
        "schema_version": SCHEMA_VERSION,
        "n_spectra": n_spec,
        "precursor_mz": round(pmz_rep, 4),
        "polarities": pols,
        "adducts": adducts,
        "n_peaks_avg": round(n_peaks_total / n_spec, 1),
        "neutral_losses": losses,
        "diagnostic_ions": diags,
        "loss_chains": chains,
        "provenance": {
            "candidate_conditioned": True,
            "independent_evidence": False,
            "eligible_for_production_score": False,
            "query_exclusion_policy": "not_applied_full_library_summary",
            "source_exclusion_policy": "not_applied_full_library_summary",
            "spectrum_set_sha256": spectrum_set_sha256,
            "source_spectrum_count": len(source_ids),
            "dataset_sources": source_datasets,
            "classification_source": classification_source or "unknown",
            "classification_source_trusted": bool(class_trusted),
            "classification_trust_reason": class_trust_reason,
        },
    }


def _build_chains(pmz: float, prod_bin: dict) -> list[list[float]]:
    """贪心 DFS: 从前体 pmz 出发, 每步减一个策展中性丢失, 要求产物落在观察峰 bin 内。
    返回 Top-K 条 [pmz, p1, p2, ...] (每步均为观察到的产物离子 m/z)。
    深度在递归**前**判 (MAX_CHAIN_DEPTH 次丢失), 避免稠密谱生成 20+ 步无意义长链。"""
    obs_set = set(round(c / BIN) for (c, _r) in prod_bin.values())

    def near_observed(m):
        # ±1 bin 容差: 链步用策展 NL 精确质量算产物, 与观察峰(可能差 ~0.02Da)会跨 bin 边界。
        b = int(round(m / BIN))
        return b in obs_set or (b - 1) in obs_set or (b + 1) in obs_set

    chains = []
    nodes = [0]

    def dfs(current_mass, path):
        if len(chains) >= TOPK_CHAIN * 4 or nodes[0] > MAX_CHAIN_NODES:
            return
        nodes[0] += 1
        depth = len(path) - 1            # 已发生的丢失次数
        if depth >= MAX_CHAIN_DEPTH:     # 深度上限 → 记录并回溯 (在前, 不再下探)
            chains.append(path)
            return
        extended = False
        for mass, _key, _formula, _name, _fgs in _NL_LIST:
            prod = current_mass - mass
            if prod < MZ_MIN:
                continue
            if near_observed(prod):
                extended = True
                dfs(prod, path + [round(prod, 4)])
        if not extended and len(path) > 1:
            chains.append(path)          # 无法继续 → 叶子链

    dfs(pmz, [round(pmz, 4)])

    # 去重 + 排序 (长链优先, 再按终点强度), 截 TopK
    uniq = []
    seen = set()
    for c in chains:
        t = tuple(c)
        if t in seen:
            continue
        seen.add(t)
        uniq.append(c)
    # 终点 m/z 在观察池里的强度 (prod_bin.values() = [(center, rel), ...])
    inten_of = {round(c, 4): r for (c, r) in prod_bin.values()}

    def term_inten(chain):
        return inten_of.get(chain[-1], 0.0)

    uniq.sort(key=lambda c: (-len(c), -term_inten(c)))
    return uniq[:TOPK_CHAIN]


# ──────────────────────────────────────────────────────────────────────
# 主流程: 流式 + 续跑
# ──────────────────────────────────────────────────────────────────────
def add_column(cur):
    cols = {r[1] for r in cur.execute("PRAGMA table_info(Compounds)")}
    if "fragmentation_pathway" not in cols:
        cur.execute("ALTER TABLE Compounds ADD COLUMN fragmentation_pathway TEXT")
        print("  ALTER TABLE Compounds ADD COLUMN fragmentation_pathway")


def load_progress(path: Path, db_name: str) -> int:
    if not path.exists():
        return 0
    try:
        d = json.loads(path.read_text())
        if d.get("db") and d["db"] != db_name:
            print(f"  ⚠ 进度文件 db={d['db']} ≠ 当前 {db_name} → 从头开始 (避免误续到别的库)")
            return 0
        return int(d.get("last_id", 0))
    except Exception:
        return 0


def save_progress(path: Path, last_id: int, db: str, n_done: int):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps({"last_id": last_id, "db": db, "n_done": n_done}))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--db", default=str(ROOT / "data" / "local_library_p2.db"),
                    help="目标库 (默认 P2 副本)")
    ap.add_argument("--limit", type=int, default=0, help="只处理前 N 个化合物 (0=全量)")
    ap.add_argument("--reset", action="store_true", help="清空续跑进度从头来")
    ap.add_argument("--progress", default=str(PROGRESS_DEFAULT))
    args = ap.parse_args()

    db = str(Path(args.db).resolve())
    progress_path = Path(args.progress)
    if args.reset and progress_path.exists():
        progress_path.unlink()

    print(f"=== P2·2b 逐化合物裂解路径挖掘 ===")
    print(f"  DB: {db}")
    con = sqlite3.connect(db)
    con.execute("PRAGMA journal_mode=WAL")   # 允许读写并发 (P1 可能同时读 live 库)
    cur = con.cursor()
    add_column(cur)
    con.commit()

    last_id = load_progress(progress_path, Path(db).name)
    n_total = cur.execute("SELECT COUNT(*) FROM Compounds").fetchone()[0]
    n_with_path = cur.execute(
        "SELECT COUNT(*) FROM Compounds WHERE fragmentation_pathway IS NOT NULL").fetchone()[0]
    print(f"  化合物: {n_total:,} | 已有 pathway: {n_with_path:,} | 续跑起点 id>{last_id}")

    # 按 compound_id **区间批处理**: 每批取一个 id 区间的全部 (化合物×谱) 行,
    # fetchall 完成后游标关闭 → 再写库+commit。避免 "游标迭代中同连接 commit 撕毁游标"。
    compound_columns = {
        row[1] for row in cur.execute("PRAGMA table_info(Compounds)")
    }
    spectra_columns = {
        row[1] for row in cur.execute("PRAGMA table_info(Spectra)")
    }
    class_source_expr = (
        "c.classification_source"
        if "classification_source" in compound_columns
        else "NULL"
    )
    spectrum_id_expr = (
        "s.spectrum_id"
        if "spectrum_id" in spectra_columns
        else "CAST(s.id AS TEXT)"
    )
    dataset_source_expr = (
        "s.dataset_source" if "dataset_source" in spectra_columns else "NULL"
    )
    instrument_expr = (
        "s.instrument_type" if "instrument_type" in spectra_columns else "NULL"
    )
    sql = f"""
        SELECT c.id, c.classification_v3, {class_source_expr},
               s.precursor_mz, s.peaks_json, s.adduct, s.charge,
               {spectrum_id_expr}, {dataset_source_expr}, {instrument_expr}
        FROM Compounds c
        JOIN Mapping m ON m.compound_id = c.id
        JOIN Spectra  s ON m.spectrum_id = s.id
        WHERE c.id > ? AND c.id <= ?
        ORDER BY c.id
    """
    max_id = cur.execute("SELECT MAX(id) FROM Compounds").fetchone()[0] or 0
    ID_STEP = 4000   # id 区间步长 (控内存: ~4000 化合物 × ~2.6 谱 ≈ 1万行/批)

    t0 = time.time()
    n_done = 0
    n_with_pathway = 0
    n_skipped_spec = 0
    n_skipped_parse = 0
    rows_seen = 0
    stopped = False
    pos = last_id
    last_written_id = last_id

    def flush(batch):
        if batch:
            cur.executemany("UPDATE Compounds SET fragmentation_pathway=? WHERE id=?", batch)
            con.commit()

    def parse_spectrum(
        pmz, peaks_json, adduct, charge,
        spectrum_id=None, dataset_source=None, instrument_type=None,
    ):
        """逐条容错解析一行谱 → dict 或 None(完全无效)。"""
        nonlocal n_skipped_spec, n_skipped_parse
        try:
            pmz_f = float(pmz) if pmz not in (None, "") else None
        except Exception:
            pmz_f = None
            n_skipped_spec += 1
        try:
            raw = json.loads(peaks_json) if peaks_json else []
        except Exception:
            raw = []
            n_skipped_parse += 1
        return {"precursor_mz": pmz_f,
                "peaks": _norm_peaks(raw),
                "adduct": adduct, "charge": charge,
                "polarity": polarity_of(adduct, charge),
                "spectrum_id": spectrum_id,
                "dataset_source": dataset_source,
                "instrument_type": instrument_type}

    batch = []
    while pos < max_id:
        hi = pos + ID_STEP
        rows = cur.execute(sql, (pos, hi)).fetchall()
        rows_seen += len(rows)
        # 按 cid 分组 (保序)
        order, by_cid = [], {}
        for (
            cid, cls_v3, cls_source, pmz, peaks_json, adduct, charge,
            spectrum_id, dataset_source, instrument_type,
        ) in rows:
            if cid not in by_cid:
                by_cid[cid] = {
                    "cls": cls_v3 or "",
                    "cls_source": cls_source,
                    "spectra": [],
                }
                order.append(cid)
            sp = parse_spectrum(
                pmz, peaks_json, adduct, charge,
                spectrum_id, dataset_source, instrument_type,
            )
            if sp is not None:
                by_cid[cid]["spectra"].append(sp)
        # 逐化合物建 pathway
        for cid in order:
            rec = by_cid[cid]
            pathway = mine_compound(
                cid,
                rec["spectra"],
                rec["cls"],
                rec["cls_source"],
            )
            if pathway is None:
                batch.append((None, cid))
            else:
                batch.append((json.dumps(pathway, ensure_ascii=False), cid))
                n_with_pathway += 1
            n_done += 1
            last_written_id = cid
            if n_done % BATCH_COMPOUNDS == 0:
                flush(batch); batch = []
                save_progress(progress_path, cid, Path(db).name, n_done)
                rate = n_done / (time.time() - t0 + 1e-6)
                print(f"    ... {n_done:,} 化合物 | {n_with_pathway:,} 有路径 | "
                      f"{rate:.0f} 个/s | last_id={cid}", flush=True)
            if args.limit and n_done >= args.limit:
                stopped = True
                break
        flush(batch); batch = []
        save_progress(progress_path, hi, Path(db).name, n_done)  # 区间已耗尽
        pos = hi
        if stopped:
            break

    elapsed = time.time() - t0
    print(f"\n  完成: 处理 {n_done:,} 化合物 | 有 pathway {n_with_pathway:,} "
          f"({100*n_with_pathway/max(n_done,1):.1f}%)")
    print(f"  跳过谱: 解析失败 {n_skipped_parse}, 前体非法 {n_skipped_spec} (共扫 {rows_seen:,} 行)")
    print(f"  耗时 {elapsed:.0f}s | 续跑点 last_id={last_written_id}")
    con.close()


if __name__ == "__main__":
    main()
