"""
mssystem v3 P1 — 分类一致性过滤 (Classification Consistency Filter)
=================================================================
把 DL 分类器从"展示徽章"升级为"交叉校验器": 解析 query 谱图的可靠化合物大类,
与匹配候选的库标签做跨骨架冲突检查, 冲突 → 软降权 (生产) / 硬否决 (opt-in)。

为何需要: 两个 DL 分类器此前从不参与决策 (search.py 在 best 选完后才调用, 仅用于
display_class)。跨类假阳性 (如谱图是糖苷却匹配到无糖甾体) 当前无任何机制否决 ——
这是最大未开发杠杆。本模块让"预测类 vs 命中类"的冲突成为评分信号。

化学细节 (避免假冲突, 保护 11/14 红线):
  - glycoside 是正交糖基修饰, 与任何骨架类兼容 (黄酮苷 / 甾体皂苷都合理) → 永不冲突。
  - {terpenoid, steroid, saponin} 同属萜类生源家族, 互不冲突 (甾体=三萜衍生; 皂苷=三萜/甾体糖苷)。
  - {lipid, fatty_acid, fatty_acyl} 同属脂质家族 (DB 已把 fatty_acid 合并进 lipid; DL 用 fatty_acyl)。
  - 其余骨架类 (flavonoid/alkaloid/phenolic_acid/coumarin/lignan/quinone/amino_acid/xanthone/
    stilbene/benzenoid/sulfonate) 互斥 — 跨之即冲突。

可靠性门控 (resolve_query_class): 仅当 query 类来自 (a) DL 且该类 per-class 准确率 ≥ 阈值
+ 置信度 ≥ 阈值, 或 (b) 无模型碎片推断且 confident, 才采信; 否则返回 None (不过滤)。
terpenoid(.60) / fatty_acyl(.63) 因 per-class 准确率低自动排除, 避免误伤。
"""
from __future__ import annotations
from typing import Optional, Tuple

# 萜类生源家族 (甾体=三萜衍生; 皂苷=三萜/甾体糖苷)
_TERP_FAMILY = {"terpenoid", "steroid", "saponin"}
# 脂质家族
_LIPID_FAMILY = {"lipid", "fatty_acid", "fatty_acyl"}
# 糖基正交修饰 — 与任何骨架类兼容
_ORTHOGONAL = {"glycoside"}
# 无信号占位
_NON = {"", "unknown", "other", "unclassified", "none", "nan"}


def _alias(c: str) -> str:
    """脂质家族折叠: fatty_acid / fatty_acyl → lipid (统一比较键)。"""
    if c in ("fatty_acid", "fatty_acyl"):
        return "lipid"
    return c


def _norm(raw) -> str:
    """规范化类别标签到比较键 (复用 fragmentation_rules.normalize_classification + 别名折叠)。

    DL class_names 用 'fatty_acyl', 而 normalize_classification 不识别它 (会回 Unknown),
    故先对原始小写串做别名命中, 再 normalize, 再别名折叠。
    """
    if not raw:
        return "Unknown"
    s = str(raw).strip().lower()
    if s in _NON:
        return "Unknown"
    if s in ("fatty_acid", "fatty_acyl"):  # 别名直接命中 (绕过 normalize 的盲区)
        return "lipid"
    try:
        from .fragmentation_rules import normalize_classification
        c = normalize_classification(raw)
    except Exception:
        c = s
    if not c or c == "Unknown":
        return "Unknown"
    return _alias(c)


def _same_family(a: str, b: str) -> bool:
    if a == b:
        return True
    if a in _TERP_FAMILY and b in _TERP_FAMILY:
        return True
    if a in _LIPID_FAMILY and b in _LIPID_FAMILY:
        return True
    return False


def is_conflict(query_class, candidate_class) -> bool:
    """两个类是否跨骨架冲突 (True = 应软降权 / 硬否决)。

    保守策略: 任一方 Unknown/无信号 → False (不过滤); 任一方 glycoside → False;
    同家族 → False。仅明确的跨骨架冲突才 True。
    """
    qc = _norm(query_class)
    cc = _norm(candidate_class)
    if qc == "Unknown" or cc == "Unknown":
        return False
    if qc in _ORTHOGONAL or cc in _ORTHOGONAL:
        return False
    if _same_family(qc, cc):
        return False
    return True


def resolve_query_class(
    query_peaks,
    query_pmz: float,
    clf_wrapper=None,
    min_dl_conf: float = 0.5,
    min_class_acc: float = 0.85,
    reliable_whitelist: str = "",
) -> Tuple[Optional[str], str]:
    """解析 query 谱图的可靠化合物大类 (per-query 算一次, 复用于所有候选)。

    优先 DL (仅当预测类 per-class 准确率 + 置信度均达标), 否则无模型碎片推断回退
    (仅 confident 且 conf≥0.5), 否则 (None, 'none') — None 表示无可靠信号, 不过滤。

    Args:
        query_peaks: [[mz, intensity], ...]
        query_pmz: precursor m/z
        clf_wrapper: MS2ClfWrapper 实例 (有 .predict / .class_acc); None 则跳过 DL
        min_dl_conf: DL 置信度下限
        min_class_acc: DL 该类 per-class 准确率下限 (动态读 ckpt, 自动排除弱类)
        reliable_whitelist: 可选显式白名单 (逗号分隔类名), 命中则视为可靠 (绕过 class_acc)

    Returns:
        (class_name | None, source)  source ∈ {'dl', 'fragment', 'none'}
    """
    # 1. DL primary — 仅可靠类才采信
    if clf_wrapper is not None:
        try:
            r = clf_wrapper.predict(query_peaks, query_pmz)
            c = r.get("class") if isinstance(r, dict) else None
            conf = float(r.get("confidence", 0.0)) if isinstance(r, dict) else 0.0
            if c and conf >= min_dl_conf:
                acc = 0.0
                ca = getattr(clf_wrapper, "class_acc", None)
                if isinstance(ca, dict):
                    acc = float(ca.get(c, 0.0))
                wl = {x.strip().lower() for x in (reliable_whitelist or "").split(",") if x.strip()}
                reliable = (acc >= min_class_acc) or (c.lower() in wl)
                if reliable:
                    return c, "dl"
        except Exception:
            pass

    # 2. 无模型碎片推断回退 (廉价; 仅 confident 才用)
    try:
        from .fragment_feature_library import classify_compound_from_fragments
        frag_mzs = [float(p[0]) for p in query_peaks
                    if p and len(p) >= 2 and float(p[0]) < query_pmz - 1.5]
        if frag_mzs:
            r = classify_compound_from_fragments(frag_mzs)
            if r and r.get("note") == "confident" and float(r.get("confidence", 0.0)) >= 0.5:
                cn = r.get("class_name")
                if cn and str(cn).strip().lower() not in _NON:
                    return cn, "fragment"
    except Exception:
        pass

    return None, "none"


CONSISTENCY_FILTER_AVAILABLE = True
