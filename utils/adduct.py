"""加合物检测模块 — 自动标注 [M+H]+, [M+Na]+, [M+NH4]+ 等常见加合物

基于质量差匹配:
1. 已知加合物质量偏移表 (正/负离子模式)
2. 对 precursor_mz 尝试反推中性分子质量 M
3. 多个候选加合物时取质量误差最小的

参考: MS-DIAL, MZmine 加合物检测策略
"""
import logging
from typing import Optional

logger = logging.getLogger(__name__)

# ── 常见加合物定义 (基于 IUPAC 标准原子量) ──────────────────────────
# 格式: (名称, 质量偏移, 电荷)
# 质量偏移 = 加合物质量 - 中性分子质量
# 例如 [M+H]+: 偏移 = 1.007276 (质子质量)
#      [M+Na]+: 偏移 = 22.989218 - 中性 = +22.989218 (不对)
# 更准确: [M+H]+ → M + 1.007276, 所以 precursor = M + 1.007276
#         M = precursor - 1.007276

POSITIVE_ADDUCTS = [
    # (名称,               加合物离子质量偏移,      电荷)
    ("[M+H]+",             1.007276,               1),
    ("[M+Na]+",           22.989218,               1),
    ("[M+K]+",            38.963158,               1),
    ("[M+NH4]+",          18.033822,               1),
    ("[M-H2O+H]+",       -17.003288,               1),  # 脱水+质子
    ("[M+H-H2O]+",        -17.003288217451,               1),  # = M + 1.007276 - 18.010565
    ("[2M+H]+",            1.007276,               1),  # 二聚体 (需特殊处理)
    ("[M+2H]2+",           2.014552,               2),  # 双电荷
    ("[M+H]+Na",          23.996494,               1),  # M+H+Na 同位素
]

NEGATIVE_ADDUCTS = [
    ("[M-H]-",            -1.007276,               1),
    ("[M+Cl]-",           34.969401261909,               1),
    ("[M+FA-H]-",         44.998201,               1),  # 甲酸
    ("[M+Hac-H]-",        59.013851,               1),  # 乙酸
    ("[M-H-H2O]-",       -19.017841150609,               1),
    ("[2M-H]-",           -1.007276,               1),
]

# 最常见加合物优先级 (检测时优先匹配)
COMMON_POSITIVE = ["[M+H]+", "[M+Na]+", "[M+NH4]+", "[M+K]+", "[M-H2O+H]+"]
COMMON_NEGATIVE = ["[M-H]-", "[M+Cl]-", "[M+FA-H]-"]


# ── 中性质量归一化 (跨加合物匹配的核心) ──────────────────────────────
# 加合物名 → shift, 使 neutral_mass = precursor_mz*|charge| - shift
# 修正表 (原 [M+H-H2O]+ 值有误; [M+H-H2O]+ = [M-H2O+H]+, shift 应为 -17.003289)
NEUTRAL_SHIFTS = {
    # 正离子
    "[M]+":         (0.0, 1),       # 永久正离子（如胆碱）
    "M+":           (0.0, 1),
    "[M+H]+":       (1.007276, 1),
    "M+H":          (1.007276, 1),
    "[M+Na]+":      (22.989218, 1),
    "M+Na":         (22.989218, 1),
    "[M+K]+":       (38.963158, 1),
    "M+K":          (38.963158, 1),
    "[M+NH4]+":     (18.033822, 1),
    "M+NH4":        (18.033822, 1),
    "[M+H-H2O]+":   (-17.003289, 1),   # 失水质子化
    "[M-H2O+H]+":   (-17.003289, 1),
    "[M+H-NH3]+":   (-16.018724, 1),   # 失氨
    # 负离子
    "[M]-":         (0.0, 1),       # 永久负离子
    "M-":           (0.0, 1),
    "[M-H]-":       (-1.007276, 1),
    "M-H":          (-1.007276, 1),
    "[M+Cl]-":      (34.969401261909, 1),
    "M+Cl":         (34.969401261909, 1),
    "[M+FA-H]-":    (44.998201, 1),    # 甲酸根
    "[M+HCOO]-":    (44.998201, 1),
    "[M+Hac-H]-":   (59.013851, 1),    # 乙酸根
    "[M+CH3COO]-":  (59.013851, 1),
    "[M-H-H2O]-":   (-19.017841150609, 1),
}


def parse_adduct_shift(adduct_str, polarity=None):
    """加合物字符串 → (shift, charge), 用于中性质量计算。未知返回 None。

    对库/样本里非标准标注 ('Positive'/'Negative'/空) 按极性默认 [M+H]+/[M-H]-。
    polarity: 'POS'/'NEG'/None。当 adduct_str 仅指示极性时用 polarity 兜底。
    """
    a = (adduct_str or "").strip()
    if a in NEUTRAL_SHIFTS:
        return NEUTRAL_SHIFTS[a]
    au = a.upper().replace(" ", "")
    if au in {k.upper().replace(" ", "") for k in NEUTRAL_SHIFTS}:
        for k, v in NEUTRAL_SHIFTS.items():
            if k.upper().replace(" ", "") == au:
                return v
    # 启发式
    if "+NA" in au:
        return (22.989218, 1)
    if "+NH4" in au or "NH4" in au:
        return (18.033822, 1)
    if "+K]" in au or au.endswith("+K") or au == "M+K":
        return (38.963158, 1)
    if "H2O" in au and "+" in au:
        return (-17.003289, 1)
    if "+CL" in au:
        return (34.969401261909, 1)
    if "FA" in au and "-" in au:
        return (44.998201, 1)
    if "HAC" in au and "-" in au:
        return (59.013851, 1)
    if ("M+H" in au or "[M+H]" in au) and "-" not in au:
        return (1.007276, 1)
    if ("M-H" in au or "[M-H]" in au) and "+" not in au:
        return (-1.007276, 1)
    # 仅极性信息 → 默认最常见加合物
    if au in ("POSITIVE", "POS", "+") or polarity == "POS":
        return (1.007276, 1)
    if au in ("NEGATIVE", "NEG", "-") or polarity == "NEG":
        return (-1.007276, 1)
    return None


def adduct_polarity(adduct_str):
    """库谱图 adduct 标注 → 'POS' / 'NEG' / None (无法判定)。

    离子模式过滤用: 把库的加合物标注映射到极性, 与 query 极性比对剔除跨模式候选。

    判定优先级 (越靠前越可靠):
      1. 末尾电荷符号 — `[M+H-H2O]+` 是正、`[M-H-H2O]-` 是负, 括号内同时含 +/- 时
         以末尾电荷为准 (脱水/脱氨这类中性丢失不改变极性)。
      2. 显式极性串 — 'Positive'/'POSITIVE'/'POS' → POS; 'Negative'/'NEGATIVE'/'NEG' → NEG
         (库里有 ~6.7万条 adduct 字段只写了极性串, 无具体加合物)。
      3. 启发式 — +H/+Na/+K/+NH4 → POS; -H/HCOO/+Cl/FA- → NEG。

    返回 None 的条目在过滤时保留 (adduct 空/异常, 不误杀)。
    """
    a = (adduct_str or "").strip()
    if not a:
        return None
    au = a.upper().replace(" ", "")
    # 1) 末尾电荷符号 (最可靠)
    if au.endswith("+"):
        return "POS"
    if au.endswith("-"):
        return "NEG"
    # 2) 显式极性串
    if au in ("POSITIVE", "POS", "+") or au.startswith("POS"):
        return "POS"
    if au in ("NEGATIVE", "NEG", "-") or au.startswith("NEG"):
        return "NEG"
    # 3) 启发式 (无括号电荷符号的简写, 如 'M+H'/'M-H'/'M+Na')
    if "+H" in au or "+NA" in au or "+K" in au or "+NH4" in au or "M+H" in au:
        return "POS"
    if "-H" in au or "HCOO" in au or "+CL" in au or ("FA" in au and "-" in au) or "M-H" in au:
        return "NEG"
    if "+" in au and "-" not in au:
        return "POS"
    if "-" in au and "+" not in au:
        return "NEG"
    return None


def neutral_mass(precursor_mz, adduct_str=None, polarity=None):
    """从 precursor m/z + 加合物标注反推中性分子质量 M。未知加合物返回 None。"""
    if precursor_mz is None or precursor_mz <= 0:
        return None
    sc = parse_adduct_shift(adduct_str, polarity)
    if sc is None:
        return None
    shift, charge = sc
    nm = precursor_mz * abs(charge) - shift
    return nm if 0 < nm < 5000 else None


def ion_mz(neutral_exact_mass, adduct_str=None):
    """Calculate theoretical ion m/z from neutral exact mass and an adduct.

    Unlike :func:`neutral_mass`, this helper does not guess from a polarity-only
    label.  Callers must provide a chemically specific adduct supported by the
    local shift table.
    """
    try:
        mass = float(neutral_exact_mass)
    except (TypeError, ValueError):
        return None
    if not (0 < mass < 5000):
        return None
    normalized = str(adduct_str or "").strip().upper().replace(" ", "")
    parsed = next(
        (
            value for key, value in NEUTRAL_SHIFTS.items()
            if key.upper().replace(" ", "") == normalized
        ),
        None,
    )
    if parsed is None:
        return None
    shift, charge = parsed
    charge = abs(int(charge))
    if charge < 1:
        return None
    value = (mass + shift) / charge
    return value if 0 < value < 5000 else None



def detect_adduct(
    precursor_mz: float,
    peaks: list = None,
    tol_da: float = 0.05,
    mode: str = "positive",
    neutral_mass_reference: float | None = None,
) -> Optional[dict]:
    """Evaluate adduct hypotheses for one precursor.

    A single precursor m/z is mathematically compatible with several
    ``M + shift`` hypotheses.  Without an independent neutral-mass reference,
    this function returns an explicit ambiguous hypothesis set and never labels
    the first common adduct as a high-confidence identification.
    """
    if precursor_mz <= 0:
        return None
    adducts = POSITIVE_ADDUCTS if mode == "positive" else NEGATIVE_ADDUCTS
    hypotheses = []
    for name, shift, charge in adducts:
        # Dimers need a separate 2M stoichiometric expression.
        if name.startswith("[2M"):
            continue
        neutral = precursor_mz * abs(charge) - shift
        if not 0 < neutral <= 2000:
            continue
        error = (
            abs(neutral - float(neutral_mass_reference))
            if neutral_mass_reference is not None
            else None
        )
        hypotheses.append(
            {
                "adduct": name,
                "neutral_mass": round(neutral, 6),
                "mass_error_da": round(error, 6) if error is not None else None,
                "charge": charge,
            }
        )

    if neutral_mass_reference is None:
        return {
            "status": "ambiguous",
            "adduct": None,
            "neutral_mass": None,
            "mass_error_da": None,
            "confidence": "unknown",
            "reason": (
                "A single precursor m/z cannot identify a specific adduct "
                "without an independent neutral-mass or ion-pair constraint"
            ),
            "hypotheses": hypotheses,
        }

    within = [
        hypothesis
        for hypothesis in hypotheses
        if hypothesis["mass_error_da"] is not None
        and hypothesis["mass_error_da"] <= tol_da
    ]
    if not within:
        return {
            "status": "no_match",
            "adduct": None,
            "neutral_mass": None,
            "mass_error_da": None,
            "confidence": "unknown",
            "hypotheses": sorted(
                hypotheses,
                key=lambda item: item["mass_error_da"]
                if item["mass_error_da"] is not None
                else float("inf"),
            )[:5],
        }
    best = min(within, key=lambda item: item["mass_error_da"])
    return {
        "status": "assigned",
        **best,
        "confidence": "high" if best["mass_error_da"] <= tol_da / 5 else "medium",
        "hypotheses": sorted(within, key=lambda item: item["mass_error_da"]),
        "evidence_source": "independent_neutral_mass",
    }


def detect_adduct_pair(
    mz1: float,
    mz2: float,
    tol_da: float = 0.05,
    mode: str = "positive",
) -> Optional[dict]:
    """
    通过两个可能的加合物峰推断加合物类型和分子量。
    例如: m/z 301([M+H]+) 和 323([M+Na]+) → 分子量 300。

    Args:
        mz1, mz2: 两个前体离子 m/z
        tol_da: 质量容差

    Returns:
        dict or None:
            {
                "adduct1": "[M+H]+",
                "adduct2": "[M+Na]+",
                "neutral_mass": 300.0,
                "mass_diff_da": 21.99,
                "mass_diff_theory": 21.9819,
                "error_da": 0.008
            }
    """
    if mz1 <= 0 or mz2 <= 0:
        return None

    diff = abs(mz2 - mz1)
    if diff < 0.1:  # 太接近，无意义
        return None

    common_names = COMMON_POSITIVE if mode == "positive" else COMMON_NEGATIVE
    definitions = {
        name: (shift, abs(charge))
        for name, shift, charge in (
            POSITIVE_ADDUCTS if mode == "positive" else NEGATIVE_ADDUCTS
        )
        if name in common_names and not name.startswith("[2M")
    }
    best = None
    best_error = tol_da
    names = list(definitions)
    for index, first in enumerate(names):
        first_shift, first_charge = definitions[first]
        if first_charge != 1:
            continue
        for second in names[index + 1:]:
            second_shift, second_charge = definitions[second]
            if second_charge != 1:
                continue
            theoretical = abs(second_shift - first_shift)
            error = abs(diff - theoretical)
            if error >= best_error:
                continue
            if first_shift <= second_shift:
                lower_name, lower_shift = first, first_shift
                upper_name, upper_shift = second, second_shift
            else:
                lower_name, lower_shift = second, second_shift
                upper_name, upper_shift = first, first_shift
            if mz1 <= mz2:
                name1, shift1 = lower_name, lower_shift
                name2, shift2 = upper_name, upper_shift
            else:
                name1, shift1 = upper_name, upper_shift
                name2, shift2 = lower_name, lower_shift
            neutral1 = mz1 - shift1
            neutral2 = mz2 - shift2
            best_error = error
            best = {
                "status": "assigned_pair",
                "adduct1": name1,
                "adduct2": name2,
                "neutral_mass": round((neutral1 + neutral2) / 2.0, 6),
                "neutral_mass_spread_da": round(abs(neutral1 - neutral2), 6),
                "mass_diff_da": round(diff, 6),
                "mass_diff_theory": round(theoretical, 6),
                "error_da": round(error, 6),
                "confidence": "high" if error <= tol_da / 5 else "medium",
                "evidence_source": "co_observed_adduct_pair",
            }

    return best


def batch_detect_adducts(precursor_mz_list: list[float],
                         tol_da: float = 0.05,
                         mode: str = "positive") -> list[dict]:
    """
    批量检测**有离子对证据**的加合物组。

    孤立 precursor 一律保持 ambiguous。只有两个共现 m/z 的差值与同一极性
    下已知加合物偏移差一致时，才进行具体标注。

    Args:
        precursor_mz_list: 前体离子 m/z 列表
        tol_da: 质量容差
        mode: "positive" 或 "negative"

    Returns:
        list[dict]: 每个输入 mz 对应的加合物标注
            {
                "precursor_mz": 301.1234,
                "adduct": "[M+H]+",
                "neutral_mass": 300.1161,
                "group_id": 0,        # 同一分子不同加合物共享 group_id
                "is_primary": True,   # 组内主加合物
            }
    """
    if not precursor_mz_list:
        return []

    parent = list(range(len(precursor_mz_list)))

    def find(index):
        while parent[index] != index:
            parent[index] = parent[parent[index]]
            index = parent[index]
        return index

    def union(left, right):
        left_root, right_root = find(left), find(right)
        if left_root != right_root:
            parent[right_root] = left_root

    evidence = [[] for _ in precursor_mz_list]
    for left in range(len(precursor_mz_list)):
        for right in range(left + 1, len(precursor_mz_list)):
            pair = detect_adduct_pair(
                precursor_mz_list[left],
                precursor_mz_list[right],
                tol_da=tol_da,
                mode=mode,
            )
            if not pair:
                continue
            union(left, right)
            evidence[left].append(
                (pair["adduct1"], pair["neutral_mass"], pair["confidence"])
            )
            evidence[right].append(
                (pair["adduct2"], pair["neutral_mass"], pair["confidence"])
            )

    group_ids = {}
    results = []
    for index, precursor_mz in enumerate(precursor_mz_list):
        items = evidence[index]
        labels = {item[0] for item in items}
        masses = [item[1] for item in items]
        if items and len(labels) == 1 and max(masses) - min(masses) <= tol_da:
            root = find(index)
            if root not in group_ids:
                group_ids[root] = len(group_ids)
            results.append(
                {
                    "precursor_mz": precursor_mz,
                    "status": "assigned_pair",
                    "adduct": next(iter(labels)),
                    "neutral_mass": round(sum(masses) / len(masses), 6),
                    "group_id": group_ids[root],
                    "confidence": (
                        "high"
                        if all(item[2] == "high" for item in items)
                        else "medium"
                    ),
                    "evidence_count": len(items),
                    "evidence_source": "co_observed_adduct_pair",
                }
            )
        else:
            ambiguous = detect_adduct(
                precursor_mz,
                tol_da=tol_da,
                mode=mode,
            )
            results.append(
                {
                    "precursor_mz": precursor_mz,
                    **(ambiguous or {}),
                    "group_id": None,
                }
            )
    return results
