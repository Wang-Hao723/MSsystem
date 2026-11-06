"""
mssystem v3 Phase 3 — v3 裂解规则库加载器 (供 Phase 4 专家系统消费)
===============================================================
统一 API 读取 kb_v3/fragmentation_rules_v3.json:
  - get_class_info(cls)          → 骨架/生源/亚型/官能团/机制/科属
  - get_diagnostic_ions(cls, pol, tol) → 数值诊断离子 (带 m/z + freq)
  - get_neutral_losses(cls, pol) → 数值中性丢失
  - get_family_chemistry(family) → 科属化学特征
  - lookup_neutral_loss(loss_da, tol) → 通用中性丢失注释

所有规则数值化, 可与真实谱图 ±tol 匹配 (Phase 4 expert_scorer 用)。
"""
from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
KB_PATH = ROOT / "kb_v3" / "fragmentation_rules_v3.json"


@lru_cache(maxsize=1)
def load_kb() -> dict:
    with open(KB_PATH, encoding="utf-8") as f:
        return json.load(f)


def classes() -> list[str]:
    return sorted(load_kb()["classes"].keys())


def get_class_info(cls: str) -> dict | None:
    """返回该类全部策展信息 (骨架/生源/亚型/官能团/机制/科属)。"""
    return load_kb()["classes"].get(cls)


def get_diagnostic_ions(cls: str, polarity: str = "both") -> list[dict]:
    """该类该极性的经验诊断离子 [{mz, freq, n_support}]。polarity: pos/neg/both。"""
    info = get_class_info(cls)
    if not info:
        return []
    emp = info.get("empirical", {})
    out = []
    pols = [polarity] if polarity != "both" else ["pos", "neg"]
    for pol in pols:
        for ion in emp.get(pol, {}).get("diagnostic_ions", []):
            out.append({**ion, "polarity": pol})
    return sorted(out, key=lambda x: -x.get("freq", 0))


def get_neutral_losses(cls: str, polarity: str = "both") -> list[dict]:
    """该类该极性的经验中性丢失 [{loss_da, freq, annotation, polarity}]。"""
    info = get_class_info(cls)
    if not info:
        return []
    emp = info.get("empirical", {})
    out = []
    pols = [polarity] if polarity != "both" else ["pos", "neg"]
    for pol in pols:
        for nl in emp.get(pol, {}).get("neutral_losses", []):
            out.append({**nl, "polarity": pol})
    return sorted(out, key=lambda x: -x.get("freq", 0))


def get_mechanisms(cls: str) -> list[dict]:
    info = get_class_info(cls)
    return info.get("fragmentation_mechanisms", []) if info else []


def get_family_chemistry(family: str) -> dict | None:
    """科属 → {特征类, 标志化合物, 裂解备注}。支持模糊匹配 (含中文/英文)。"""
    fams = load_kb()["plant_family_chemistry"]
    if family in fams:
        return fams[family]
    low = family.lower()
    for k, v in fams.items():
        if family in k or low in k.lower():
            return v
    return None


def families_for_class(cls: str) -> list[str]:
    info = get_class_info(cls)
    return info.get("plant_families", []) if info else []


def lookup_neutral_loss(loss_da: float, tol: float = 0.02) -> str | None:
    """通用中性丢失表反查 (hexose/deoxyhexose/CO2...)。"""
    cnl = load_kb()["common_neutral_losses"]
    best = None
    for k, v in cnl.items():
        if abs(loss_da - float(k)) <= tol:
            return v
    return best


def match_diagnostic_ions(query_peaks, cls: str, polarity: str = "both",
                          mz_tol: float = 0.05) -> list[dict]:
    """把查询谱的碎片与该类诊断离子匹配 → 命中的诊断离子 (含 freq)。Phase 4 用。"""
    hits = []
    ions = get_diagnostic_ions(cls, polarity)
    qmz = [p[0] for p in query_peaks]
    for ion in ions:
        if any(abs(q - ion["mz"]) <= mz_tol for q in qmz):
            hits.append(ion)
    return hits


if __name__ == "__main__":
    # 自检
    kb = load_kb()
    print(f"v3 KB: {kb['meta']['n_classes']} 类, schema={kb['meta']['version']}")
    print(f"\n类: {classes()}")
    # flavonoid 抽查
    print(f"\n[flavonoid] 科属: {families_for_class('flavonoid')}")
    print(f"[flavonoid] 诊断离子(neg, top5): {[(i['mz'],i['freq']) for i in get_diagnostic_ions('flavonoid','neg')[:5]]}")
    print(f"[flavonoid] 中性丢失: {[(n['loss_da'],n.get('annotation')) for n in get_neutral_losses('flavonoid')]}")
    # 科属反查
    fc = get_family_chemistry("毛茛科")
    print(f"\n[毛茛科] 特征: {fc['classes']} 标志物: {fc['markers']}")
    # 通用NL反查
    print(f"\nlookup(162.05) = {lookup_neutral_loss(162.05)}")
    print(f"lookup(43.99)  = {lookup_neutral_loss(43.99)}")
