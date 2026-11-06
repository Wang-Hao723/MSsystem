"""
mssystem v3 Phase 4 — v3 裂解规则库专家加分 + 结构推理
=====================================================
消费 kb_v3 数值规则 (经验挖掘的诊断离子/中性丢失, 带频率), 给候选化合物打"化学一致性"分,
并生成可解释的结构推理链 (脱糖→苷元→母核→诊断离子)。

与原 expert_scorer 并存 (纯增量): 取 max(base_expert, v3_score) 融合, 让原本喂符号规则
而基本失效的专家权重 (0.05) 真正承载化学信号。config.use_expert_v3 开关 (默认开, 关=v2 行为)。

API:
  score_and_reason(query_peaks, query_pmz, candidate_class, polarity) ->
      {v3_score, matched_ions, matched_losses, reasoning}
"""
from __future__ import annotations

import math
from typing import Any

MZ_TOL = 0.05  # 诊断离子 / 中性丢失 匹配容差 (Da)
IDF_BIN = 0.05  # 跨类 DF 分箱粒度 (用于特异性 IDF)

# 强度地板: <1% base peak 的峰视为噪声, 不参与诊断离子/中性丢失匹配。
# 此前 expert_v3 取 p[0] 丢强度, 0.1% 噪声峰与基峰等权 → 专家分过噪 (是其被降到 0.05 的根因之一)。
try:
    from config.settings import EXPERT_INTENSITY_FLOOR as _INT_FLOOR
except Exception:  # 独立调用 / 单测时无 config
    _INT_FLOOR = 0.01

# 规则缺失/无效类时返回零分, 不影响打分
_NON_CLASSES = {"", "unknown", "other", "unclassified", "none", "nan"}

# ── 跨类特异性 IDF (TF-IDF 思路) ──
# 一个离子若被多个类同时判为"诊断"(如 m/z 91/107/121 在萜/甾体/黄酮都有),
# 则特异性低、IDF 低 → 降权; 只被少数类特有的离子 IDF 高 → 升权。
# 修复 quercetin 误得高分 vs steroid 的 bug (萜/甾体的泛芳香碎片被黄酮也产生)。
_ION_DF: dict[int, set] = {}   # bin -> set(classes)
_LOSS_DF: dict[int, set] = {}
_N_CLASSES = 1
_IDF_BUILT = False


def _ensure_idf() -> None:
    global _IDF_BUILT, _N_CLASSES
    if _IDF_BUILT:
        return
    try:
        from kb_v3 import loader as kbv3
    except Exception:
        _IDF_BUILT = True
        return
    ion_df: dict[int, set] = {}
    loss_df: dict[int, set] = {}
    classes = kbv3.classes()
    for cls in classes:
        for ion in kbv3.get_diagnostic_ions(cls, "both"):
            b = int(round(ion["mz"] / IDF_BIN))
            ion_df.setdefault(b, set()).add(cls)
        for nl in kbv3.get_neutral_losses(cls, "both"):
            b = int(round(nl["loss_da"] / IDF_BIN))
            loss_df.setdefault(b, set()).add(cls)
    _ION_DF.update(ion_df)
    _LOSS_DF.update(loss_df)
    _N_CLASSES = max(len(classes), 1)
    _IDF_BUILT = True


def _ion_idf(mz: float) -> float:
    """离子特异性: log(N / (1+DF))。跨类越常见 IDF 越低 (min 钳到 0.3 避免全压制)。"""
    b = int(round(mz / IDF_BIN))
    df = len(_ION_DF.get(b, ()))
    return max(math.log(_N_CLASSES / (1 + df)), 0.3)


def _loss_idf(loss_da: float) -> float:
    b = int(round(loss_da / IDF_BIN))
    df = len(_LOSS_DF.get(b, ()))
    return max(math.log(_N_CLASSES / (1 + df)), 0.3)


def _polarity(adduct: str | None, declared_polarity: str | None = None) -> str:
    # 统一走 adduct_polarity() (审计 #4): 旧"含-不含+"启发式对 [M+HCOO]- 等含 '+' 负加合物
    # 返回 "both" → Expert V3 对负离子化合物应用正负双套规则。现 POS→pos / NEG→neg / None→both。
    declared = str(declared_polarity or "").strip().lower()
    if declared in {"pos", "positive"}:
        return "pos"
    if declared in {"neg", "negative"}:
        return "neg"
    from utils.adduct import adduct_polarity as _ap
    return {"POS": "pos", "NEG": "neg"}.get(_ap(adduct), "both")


def score_and_reason(
    query_peaks: list[list[float]],
    query_pmz: float,
    candidate_class: str,
    adduct: str | None = None,
    query_polarity: str | None = None,
    query_charge: int | None = None,
    *,
    classification_reliable: bool = False,
    classification_reason: str | None = None,
) -> dict[str, Any]:
    """v3 KB 专家加分 + 结构推理链。

    v3_score = 0.6·诊断离子匹配率(频率加权) + 0.4·中性丢失匹配率; 无规则→0。
    """
    cls = (candidate_class or "").strip()
    if classification_reliable is not True:
        return _empty(
            cls,
            status="not_applicable",
            reason=classification_reason or "classification_unconfirmed",
            classification_reliable=False,
        )
    if not query_peaks or cls.lower() in _NON_CLASSES:
        return _empty(
            cls,
            reason="missing_query_or_classification",
            classification_reliable=True,
        )

    try:
        from kb_v3 import loader as kbv3
    except Exception:
        return _empty(
            cls,
            status="unavailable",
            reason="v3_knowledge_base_unavailable",
            classification_reliable=True,
        )

    pol = _polarity(adduct, query_polarity)
    # 强度感知: 按 base peak 归一化后过滤 < _INT_FLOOR 的噪声峰 (修复此前丢强度致专家分过噪)
    _iv = [float(p[1]) for p in query_peaks if p and len(p) >= 2]
    _imax = max(_iv) if _iv else 0.0
    _thr = _INT_FLOOR * _imax
    qmz = [float(p[0]) for p in query_peaks
           if p and len(p) >= 2 and (_imax <= 0 or float(p[1]) >= _thr)]
    if not qmz:
        return _empty(
            cls,
            reason="no_query_peaks_above_intensity_floor",
            classification_reliable=True,
        )

    _ensure_idf()

    # ── 1. 诊断离子匹配 (TF·IDF 加权: 频率 × 跨类特异性) ──
    ions = kbv3.get_diagnostic_ions(cls, pol)
    matched_ions = []
    for ion in ions:
        if any(abs(q - ion["mz"]) <= MZ_TOL for q in qmz):
            matched_ions.append(ion)
    if ions:
        num = sum(i.get("freq", 0) * _ion_idf(i["mz"]) for i in matched_ions)
        den = sum(i.get("freq", 0) * _ion_idf(i["mz"]) for i in ions)
        ion_score = num / den if den > 0 else 0.0
    else:
        ion_score = 0.0

    # ── 2. 中性丢失匹配 (precursor - fragment, 同样 TF·IDF) ──
    try:
        multiply_charged = abs(int(query_charge)) > 1 if query_charge is not None else False
    except (TypeError, ValueError):
        multiply_charged = False
    losses = [] if multiply_charged else kbv3.get_neutral_losses(cls, pol)
    matched_losses = []
    for nl in losses:
        lda = nl["loss_da"]
        if any(abs((query_pmz - q) - lda) <= MZ_TOL for q in qmz):
            matched_losses.append(nl)
    if losses:
        num = sum(l.get("freq", 0) * _loss_idf(l["loss_da"]) for l in matched_losses)
        den = sum(l.get("freq", 0) * _loss_idf(l["loss_da"]) for l in losses)
        loss_score = num / den if den > 0 else 0.0
    else:
        loss_score = 0.0

    # ── 3. 综合 v3 分 (有任一类规则才给分; 否则 0 不惩罚) ──
    if ions or losses:
        v3_score = 0.6 * ion_score + 0.4 * loss_score
    else:
        v3_score = 0.0

    return {
        "v3_score": round(min(v3_score, 1.0), 4),
        "matched_ions": [{"mz": i["mz"], "freq": i.get("freq"), "polarity": i.get("polarity")} for i in matched_ions[:8]],
        "matched_losses": [{"loss_da": l["loss_da"], "annotation": l.get("annotation"), "freq": l.get("freq")} for l in matched_losses[:8]],
        "reasoning": _reason(cls, query_pmz, matched_losses, matched_ions, kbv3),
        "status": "scored" if (ions or losses) else "not_applicable",
        "reason": "trusted_class_rules_evaluated" if (ions or losses) else "no_rules_for_class",
        "classification_reliable": True,
        "classification_reason": classification_reason or "independently_confirmed",
    }


def _reason(cls: str, pmz: float, matched_losses: list, matched_ions: list, kbv3) -> list[str]:
    """构建可解释结构推理链。"""
    info = kbv3.get_class_info(cls) or {}
    chain = [f"分类: {cls} ({info.get('class_en', '')})", f"母核: {info.get('core_skeleton', '未知')}"]
    # 脱糖→苷元 推理
    agly = pmz
    sugar_lost = False
    for nl in matched_losses:
        ann = nl.get("annotation") or kbv3.lookup_neutral_loss(nl["loss_da"]) or f"{nl['loss_da']:.1f}Da"
        chain.append(f"检测到中性丢失 {nl['loss_da']:.2f} ({ann})")
        if "hexose" in ann or "deoxyhexose" in ann or "rutinose" in ann or "pentose" in ann or "glucuronic" in ann:
            agly -= nl["loss_da"]
            sugar_lost = True
    if sugar_lost:
        chain.append(f"→ 推测脱糖得苷元 m/z ≈ {agly:.1f}")
    # 诊断离子
    if matched_ions:
        top = matched_ions[0]
        chain.append(f"诊断离子 m/z {top['mz']:.2f} ({cls} 特征碎片, 出现率 {top.get('freq', 0):.0%})")
    # 生源途径 (辅助)
    if info.get("biosynthesis_pathway"):
        chain.append(f"生源: {info['biosynthesis_pathway']}")
    return chain


def _empty(
    cls: str,
    *,
    status: str = "not_applicable",
    reason: str = "no_rules_for_class",
    classification_reliable: bool = False,
) -> dict[str, Any]:
    return {
        "v3_score": 0.0,
        "matched_ions": [],
        "matched_losses": [],
        "reasoning": (
            []
            if cls.lower() in _NON_CLASSES
            else [f"分类 {cls} 未进入生产 v3 规则评分: {reason}"]
        ),
        "status": status,
        "reason": reason,
        "classification_reliable": bool(classification_reliable),
    }


EXPERT_V3_AVAILABLE = True
