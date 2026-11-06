"""
质谱碎片特征标准库 — Fragment Feature Standard Library
========================================================

基于三份权威资料构建的质谱碎片特征标准库，编码化合物类别 → 特征碎片 →
裂解路径 → 诊断离子的完整映射关系。

知识来源:
  1. 《质谱分析技术原理与应用》台湾质谱学会 (2019) — MS原理与各领域应用
  2. 《中药质谱分析》刘淑莹等 (2012) — 中药复杂体系ESI-MS裂解特征
  3. Spectraverse (Gupta et al., Anal. Chem. 2026) — 公共MS/MS谱库策展标准

核心功能:
  1. 化合物类别 → 特征碎片精确质量数据库
  2. 多级裂解路径 (逐级糖丢失 → 苷元裂解 → 小分子丢失)
  3. 结构解析推理规则 (从碎片逆向推导结构)
  4. 异构体区分离子对
  5. 加合物与电离模式特征

Author: Enhanced Expert System Module
Version: 3.0.0
"""

from __future__ import annotations

import json
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Any, Dict, FrozenSet, List, Optional, Set, Tuple

from .rule_quality import (
    collision_energy_value_unit,
    formula_mz_consistent,
    normalize_adduct,
    normalize_polarity,
    smiles_matches_all,
    smiles_matches_any,
)

# ─────────────────────────────────────────────────────────────────
# 基础数据结构
# ─────────────────────────────────────────────────────────────────


@dataclass(frozen=True)
class FragmentFeature:
    """单个碎片特征条目"""
    mz: float                      # 精确 m/z
    formula: str                   # 分子式
    name: str                      # 中文名称
    name_en: str                   # 英文名称
    intensity_rank: str            # "base" / "major" / "medium" / "minor"
    structural_origin: str         # 结构来源描述
    fragmentation_path: str        # 裂解路径 (如: "糖苷键断裂→脱己糖")
    diagnostic_value: float = 0.5  # 诊断价值 (0-1, 越高越具特征性)
    polarity: str = "positive"     # positive/negative/both
    adduct_preference: str = "[M+H]+"  # 常见加合物类型
    compatible_class_terms: Tuple[str, ...] = ()
    required_smarts: Tuple[str, ...] = ()
    required_smarts_all: Tuple[str, ...] = ()
    allowed_adducts: Tuple[str, ...] = ()
    collision_energy_min: Optional[float] = None
    collision_energy_max: Optional[float] = None
    collision_energy_unit: str = ""
    references: Tuple[str, ...] = ()
    evidence_level: str = "legacy_curated"


@dataclass(frozen=True)
class CompoundClassProfile:
    """化合物类别的完整质谱裂解特征档案"""
    class_name: str                    # 类别名称
    class_name_cn: str                 # 中文名称
    parent_class: str                  # 父类别
    characteristic_fragments: Tuple[FragmentFeature, ...]  # 特征碎片
    diagnostic_ion_pairs: Tuple[Tuple[float, float], ...]  # 诊断离子对 (用于异构体区分)
    fragmentation_pathway: str         # 完整裂解路径描述
    sequential_losses: Tuple[str, ...] # 顺序中性丢失序列
    isomer_differentiation_rules: Tuple[str, ...]  # 异构体区分规则
    reference_compounds: Tuple[str, ...]  # 代表性化合物
    literature_refs: Tuple[str, ...]   # 文献来源
    confidence_score: float = 0.85     # 该特征档案的可靠性


@dataclass(frozen=True)
class AdductProfile:
    """加合物离子特征"""
    adduct_name: str           # [M+H]+, [M+Na]+, etc.
    mass_shift: float          # 相对于中性分子的质量偏移
    polarity: str              # positive/negative
    prevalence: float          # 在代谢组学中的普遍性 (0-1)
    typical_compounds: str     # 典型形成该加合物的化合物类型
    fragmentation_behavior: str  # 对该加合物的碎裂行为描述


# ─────────────────────────────────────────────────────────────────
# 第一部分: 化合物类别 → 特征碎片数据库
# ─────────────────────────────────────────────────────────────────

# ═══════════════════════════════════════════════════════════════
# 1. 黄酮类 (Flavonoids) — 最详细的特征库
# ═══════════════════════════════════════════════════════════════

FLAVONOID_FRAGMENTS: Tuple[FragmentFeature, ...] = (
    # 注: 糖氧碳鎓离子 (m/z 163, 147, 133) 已移至 sugar 类别专属碎片
    # 黄酮苷元的糖丢失通过 COMPOUND_CLASS_PROFILES["flavonoid"].sequential_losses 记录
    # m/z 109.0290 (三重重叠 flavonoid/phenolic_acid/sugar) 已移除以消除歧义


    # ── 黄酮苷元 A环 rDA 裂解 (Retro-Diels-Alder of C-ring) ──
    FragmentFeature(
        mz=152.0473, formula="C8H8O3+", name="^1,3A+ (5,7-二OH-A环)",
        name_en="^1,3A+ (5,7-diOH A-ring)", intensity_rank="base",
        structural_origin="5,7-二羟基黄酮A环rDA裂解产物",
        fragmentation_path="C环rDA→A环碎片(1,3键断裂)", diagnostic_value=0.98,
    ),
    FragmentFeature(
        mz=153.0188, formula="C7H5O4+", name="^1,3A+ (5,7-二OH-4-酮-A环,黄酮醇)",
        name_en="^1,3A+ (flavonol A-ring with 4-keto)",
        intensity_rank="base", structural_origin="黄酮醇(槲皮素/山奈酚类)A环rDA裂解产物",
        fragmentation_path="C环rDA→A环(黄酮醇特征)", diagnostic_value=0.98,
    ),
    FragmentFeature(
        mz=165.0188, formula="C8H5O4+", name="^1,3A+ (5,6,7-三OH-A环)",
        name_en="^1,3A+ (5,6,7-triOH A-ring)",
        intensity_rank="major", structural_origin="多羟基黄酮A环rDA裂解",
        fragmentation_path="C环rDA→多OH-A环", diagnostic_value=0.92,
    ),
    FragmentFeature(
        mz=121.0290, formula="C7H5O2+", name="^1,3B+ (无OH-B环)",
        name_en="^1,3B+ (B-ring, 0 OH)", intensity_rank="major",
        structural_origin="无羟基取代B环的rDA裂解产物",
        fragmentation_path="C环rDA→B环碎片", diagnostic_value=0.90,
    ),
    FragmentFeature(
        mz=137.0239, formula="C7H5O3+", name="^1,3B+ (4'-OH-B环或3',4'-diOH)",
        name_en="^1,3B+ (mono/di-OH B-ring)", intensity_rank="major",
        structural_origin="单/二羟基B环的rDA裂解产物",
        fragmentation_path="C环rDA→OH-B环", diagnostic_value=0.95,
    ),
    FragmentFeature(
        mz=153.0188, formula="C7H5O4+", name="^1,3B+ (3',4',5'-triOH-B环)",
        name_en="^1,3B+ (triOH B-ring, 如myricetin)",
        intensity_rank="major", structural_origin="三羟基B环rDA裂解产物",
        fragmentation_path="C环rDA→triOH-B环", diagnostic_value=0.93,
    ),

    # ── 黄酮二糖苷特征碎片 (Diglycoside-specific) ──
    FragmentFeature(
        mz=309.0974, formula="C14H17O8+", name="黄酮单糖苷碎片(1→2连接)",
        name_en="Flavonoid monoglycoside fragment (1→2 linkage)",
        intensity_rank="medium", structural_origin="二糖苷中1→2糖连接处断裂",
        fragmentation_path="二糖苷→单糖苷", diagnostic_value=0.88,
    ),
    FragmentFeature(
        mz=471.1288, formula="C23H27O11+", name="黄酮-葡萄糖醛酸苷",
        name_en="Flavonoid monoglucuronide",
        intensity_rank="major", structural_origin="葡萄糖醛酸苷元离子",
        fragmentation_path="脱葡萄糖醛酸(176Da)", diagnostic_value=0.90,
    ),

    # ── 异黄酮特征 (Isoflavone-specific) ──
    FragmentFeature(
        mz=118.0419, formula="C8H6O+", name="异黄酮B环碎片",
        name_en="Isoflavone B-ring fragment",
        intensity_rank="major", structural_origin="异黄酮B环(3位连接)的特征裂解",
        fragmentation_path="异黄酮rDA→B环(不同于黄酮)", diagnostic_value=0.92,
    ),
    FragmentFeature(
        mz=132.0575, formula="C9H8O+", name="异黄酮B环(4'-OCH3)",
        name_en="Isoflavone B-ring (4'-methoxy)",
        intensity_rank="medium", structural_origin="甲基化异黄酮B环裂解",
        fragmentation_path="异黄酮rDA→OCH3-B环", diagnostic_value=0.88,
    ),

    # ── 小分子中性丢失 (后续裂解) ──
    FragmentFeature(
        mz=134.0368, formula="C8H6O2+", name="A环脱H2O(5-OH)",
        name_en="A-ring after H2O loss (5-OH)",
        intensity_rank="medium", structural_origin="152碎片脱H2O (5-OH邻位效应)",
        fragmentation_path="1,3A+→-H2O", diagnostic_value=0.80,
    ),
    FragmentFeature(
        mz=124.0158, formula="C6H4O3+", name="A环脱CO",
        name_en="A-ring after CO loss",
        intensity_rank="minor", structural_origin="152碎片脱CO产生",
        fragmentation_path="1,3A+→-CO", diagnostic_value=0.70,
    ),
)

# ═══════════════════════════════════════════════════════════════
# 2. 皂苷类 (Saponins) — 三萜皂苷 + 甾体皂苷
# ═══════════════════════════════════════════════════════════════

SAPONIN_FRAGMENTS: Tuple[FragmentFeature, ...] = (
    # ── 糖链顺序丢失 (Sequential Sugar Loss) ──
    FragmentFeature(
        mz=162.0528, formula="C6H10O5", name="己糖残基(葡萄糖/半乳糖)",
        name_en="Hexose residue (Glc/Gal)",
        intensity_rank="base", structural_origin="糖链末端或内部己糖丢失",
        fragmentation_path="糖苷键顺序断裂→脱己糖", diagnostic_value=0.97,
    ),
    FragmentFeature(
        mz=146.0579, formula="C6H10O4", name="脱氧己糖(鼠李糖/夫糖)",
        name_en="Deoxyhexose (Rha/Fuc)",
        intensity_rank="base", structural_origin="脱氧糖糖苷键断裂",
        fragmentation_path="脱氧糖苷键断裂", diagnostic_value=0.96,
    ),
    FragmentFeature(
        mz=132.0423, formula="C5H8O4", name="戊糖(阿拉伯糖/木糖)",
        name_en="Pentose (Ara/Xyl)",
        intensity_rank="major", structural_origin="戊糖苷键断裂",
        fragmentation_path="戊糖苷键断裂", diagnostic_value=0.94,
    ),
    FragmentFeature(
        mz=176.0321, formula="C6H8O6", name="葡萄糖醛酸",
        name_en="Glucuronic acid (GlcA)",
        intensity_rank="major", structural_origin="葡萄糖醛酸苷键断裂",
        fragmentation_path="脱GlcA", diagnostic_value=0.93,
    ),

    # ── 四环三萜皂苷特征 (Tetracyclic Triterpenoid) ──
    FragmentFeature(
        mz=109.1017, formula="C8H13+", name="三萜骨架碎片(C8H13+)",
        name_en="Triterpenoid skeleton C8H13+ fragment",
        intensity_rank="major", structural_origin="四环三萜C环裂解",
        fragmentation_path="三萜骨架→C环开环", diagnostic_value=0.88,
    ),
    FragmentFeature(
        mz=189.1643, formula="C14H21+", name="三萜D/E环碎片",
        name_en="Triterpenoid D/E ring fragment",
        intensity_rank="medium", structural_origin="达玛烷型三萜D环侧链裂解",
        fragmentation_path="侧链裂解→D/E环", diagnostic_value=0.85,
    ),
    FragmentFeature(
        mz=203.1800, formula="C15H23+", name="三萜四环骨架",
        name_en="Triterpenoid tetracyclic core",
        intensity_rank="medium", structural_origin="四环三萜骨架(失去全部糖和侧链)",
        fragmentation_path="完全脱糖→骨架离子", diagnostic_value=0.87,
    ),

    # 注: 甾体碎片 (m/z 255, 271, 253) 已移至独立的 steroid 类别档案
    # 注: 薯蓣皂苷元碎片 (m/z 415, 397, 282) 已移至 steroid 类别档案
    # 甾体碎片是苷元特征，不属于皂苷(糖苷)的专有碎片
    # 皂苷专有碎片仅包含糖链丢失 + 三萜骨架碎片
)

# ═══════════════════════════════════════════════════════════════
# 3. 生物碱类 (Alkaloids) — 乌头碱型 + 原小檗碱型
# ═══════════════════════════════════════════════════════════════

ALKALOID_FRAGMENTS: Tuple[FragmentFeature, ...] = (
    # ── 乌头碱型 C20-二萜生物碱 (Aconitine-type) ──
    FragmentFeature(
        mz=586.3016, formula="C32H44NO9+", name="乌头碱 [M+H]+",
        name_en="Aconitine [M+H]+",
        intensity_rank="base", structural_origin="乌头碱准分子离子",
        fragmentation_path="前体离子", diagnostic_value=1.0,
    ),
    FragmentFeature(
        mz=526.2808, formula="C30H40NO7+", name="乌头碱脱乙酸",
        name_en="Aconitine - CH3COOH (deacetate)",
        intensity_rank="major", structural_origin="C8位乙酰氧基丢失",
        fragmentation_path="脱乙酸(60Da)→去乙酰乌头碱", diagnostic_value=0.95,
    ),
    FragmentFeature(
        mz=494.2543, formula="C29H36NO6+", name="乌头碱脱乙酸-甲醇",
        name_en="Aconitine - HOAc - CH3OH",
        intensity_rank="major", structural_origin="8位脱乙酸→再脱甲醇",
        fragmentation_path="脱乙酸→脱甲醇", diagnostic_value=0.93,
    ),
    FragmentFeature(
        mz=462.2278, formula="C28H32NO5+", name="乌头碱多重丢失",
        name_en="Aconitine - HOAc - CH3OH - CH3OH (bis-demethanol)",
        intensity_rank="medium", structural_origin="连续丢失乙酸和两分子甲醇",
        fragmentation_path="脱乙酸→脱甲醇→脱甲醇", diagnostic_value=0.88,
    ),
    FragmentFeature(
        mz=105.0335, formula="C7H5O+", name="苯甲酰基离子",
        name_en="Benzoyl ion (C6H5CO+)",
        intensity_rank="major", structural_origin="C14位苯甲酰酯基裂解",
        fragmentation_path="C14酯基裂解→PhCO+", diagnostic_value=0.92,
    ),
    FragmentFeature(
        mz=122.0964, formula="C8H12N+", name="二萜骨架含氮碎片",
        name_en="Diterpenoid N-containing fragment",
        intensity_rank="medium", structural_origin="乌头碱骨架N-乙基裂解",
        fragmentation_path="骨架裂解→含N碎片", diagnostic_value=0.85,
    ),

    # ── 原小檗碱型生物碱 (Protoberberine-type) ──
    FragmentFeature(
        mz=336.1236, formula="C20H18NO4+", name="小檗碱 [M]+",
        name_en="Berberine [M]+ (quaternary ammonium)",
        intensity_rank="base", structural_origin="小檗碱季铵盐分子离子",
        fragmentation_path="前体离子(无加合)", diagnostic_value=1.0,
        adduct_preference="[M]+",
    ),
    FragmentFeature(
        mz=321.1001, formula="C19H15NO4+", name="小檗碱脱CH3",
        name_en="Berberine - CH3 (demethyl)",
        intensity_rank="major", structural_origin="甲氧基脱甲基",
        fragmentation_path="脱CH3→去甲小檗碱", diagnostic_value=0.94,
    ),
    FragmentFeature(
        mz=320.0923, formula="C19H14NO4+", name="小檗碱脱CH4",
        name_en="Berberine - CH4 (demethyl + H)",
        intensity_rank="major", structural_origin="甲氧基脱甲基重排",
        fragmentation_path="脱CH4→重排产物", diagnostic_value=0.92,
    ),
    FragmentFeature(
        mz=306.0766, formula="C18H12NO4+", name="小檗碱脱2×CH3",
        name_en="Berberine - 2CH3 (bisdemethyl)",
        intensity_rank="medium", structural_origin="两个甲氧基脱甲基",
        fragmentation_path="脱2×CH3", diagnostic_value=0.88,
    ),
    FragmentFeature(
        mz=292.0974, formula="C18H14NO3+", name="小檗碱脱C2H4O",
        name_en="Berberine - C2H4O (亚甲二氧基开环)",
        intensity_rank="medium", structural_origin="亚甲二氧基(-O-CH2-O-)开环",
        fragmentation_path="亚甲二氧基开环→脱CH2O", diagnostic_value=0.90,
    ),
    FragmentFeature(
        mz=278.0817, formula="C18H12NO2+", name="小檗碱多步碎片",
        name_en="Berberine multi-step fragment",
        intensity_rank="minor", structural_origin="脱2CH3+脱CO",
        fragmentation_path="脱2CH3→脱CO", diagnostic_value=0.78,
    ),

    # ── 黄嘌呤类生物碱 (Xanthine alkaloids) ──
    FragmentFeature(
        mz=138.0667, formula="C6H8N3O+", name="咖啡因特征碎片",
        name_en="Caffeine rDA fragment (imidazole opening)",
        intensity_rank="base", structural_origin="黄嘌呤咪唑环rDA开环→脱CH3NCO",
        fragmentation_path="咪唑rDA→脱CH3NCO(57Da)", diagnostic_value=0.97,
    ),
    FragmentFeature(
        mz=110.0713, formula="C5H8N3+", name="咖啡因二级碎片",
        name_en="Caffeine 138→ -CO secondary fragment",
        intensity_rank="major", structural_origin="138碎片脱CO产生",
        fragmentation_path="138→脱CO(28Da)", diagnostic_value=0.93,
    ),
    FragmentFeature(
        mz=83.0604, formula="C4H7N2+", name="咖啡因三级碎片(咪唑/嘧啶)",
        name_en="Caffeine tertiary fragment (imidazole/pyrimidine)",
        intensity_rank="medium", structural_origin="110碎片脱HCN产生",
        fragmentation_path="110→脱HCN(27Da)", diagnostic_value=0.88,
    ),
    FragmentFeature(
        mz=67.0296, formula="C3H3N2+", name="咖啡因四级碎片",
        name_en="Caffeine quaternary fragment (pyrazole)",
        intensity_rank="minor", structural_origin="83碎片继续裂解",
        fragmentation_path="83→脱CH4", diagnostic_value=0.75,
    ),

    # ── 喹嗪啶类生物碱 (Quinolizidine) ──
    FragmentFeature(
        mz=150.1283, formula="C10H16N+", name="喹嗪啶rDA开环特征离子",
        name_en="Quinolizidine rDA opening ion",
        intensity_rank="base", structural_origin="喹嗪啶双环rDA开环产物",
        fragmentation_path="喹嗪啶rDA开环", diagnostic_value=0.95,
    ),
    FragmentFeature(
        mz=136.1126, formula="C9H14N+", name="喹嗪啶脱CH2",
        name_en="Quinolizidine - CH2 fragment",
        intensity_rank="major", structural_origin="150碎片脱CH2",
        fragmentation_path="150→脱CH2", diagnostic_value=0.85,
    ),
)

# ═══════════════════════════════════════════════════════════════
# 4. 酚酸类 (Phenolic Acids) — 丹参素、迷迭香酸等
# ═══════════════════════════════════════════════════════════════

PHENOLIC_ACID_FRAGMENTS: Tuple[FragmentFeature, ...] = (
    FragmentFeature(
        mz=163.0395, formula="C9H7O3+", name="对-香豆酸 [M+H]+",
        name_en="p-Coumaric acid [M+H]+",
        intensity_rank="base", structural_origin="对-香豆酸准分子离子",
        fragmentation_path="前体离子", diagnostic_value=1.0,
    ),
    FragmentFeature(
        mz=119.0497, formula="C8H7O+", name="对-香豆酸脱CO2",
        name_en="p-Coumaric acid - CO2 (decarboxylation)",
        intensity_rank="major", structural_origin="羧基脱CO2产生乙烯基苯酚离子",
        fragmentation_path="脱CO2(44Da)", diagnostic_value=0.95,
    ),
    FragmentFeature(
        mz=147.0446, formula="C9H7O2+", name="肉桂酸脱H2O",
        name_en="Cinnamic acid - H2O",
        intensity_rank="medium", structural_origin="酚酸脱水",
        fragmentation_path="脱H2O(18Da)", diagnostic_value=0.82,
    ),
    FragmentFeature(
        mz=109.0290, formula="C6H5O2+", name="二羟基苯离子",
        name_en="Dihydroxybenzene ion (catechol type)",
        intensity_rank="major", structural_origin="咖啡酸类脱侧链→二羟基苯",
        fragmentation_path="脱丙烯酸侧链→儿茶酚离子", diagnostic_value=0.90,
    ),
    FragmentFeature(
        mz=197.0450, formula="C9H9O5+", name="丹酚酸类特征碎片",
        name_en="Salvianolic acid fragment (danshensu type)",
        intensity_rank="major", structural_origin="丹参素脱羧产物",
        fragmentation_path="脱CO2→丹酚酸碎片", diagnostic_value=0.88,
    ),
    FragmentFeature(
        mz=179.0344, formula="C9H7O4+", name="咖啡酸脱H2O",
        name_en="Caffeic acid - H2O (lactone formation)",
        intensity_rank="medium", structural_origin="咖啡酸邻位羟基脱水形成内酯",
        fragmentation_path="邻位效应→内酯化脱H2O", diagnostic_value=0.85,
    ),
)

# ═══════════════════════════════════════════════════════════════
# 5. 香豆素类 (Coumarins)
# ═══════════════════════════════════════════════════════════════

COUMARIN_FRAGMENTS: Tuple[FragmentFeature, ...] = (
    FragmentFeature(
        mz=147.0446, formula="C9H7O2+", name="香豆素 [M+H]+",
        name_en="Coumarin [M+H]+",
        intensity_rank="base", structural_origin="简单香豆素准分子离子",
        fragmentation_path="前体离子", diagnostic_value=1.0,
    ),
    FragmentFeature(
        mz=103.0548, formula="C8H7+", name="香豆素脱CO2",
        name_en="Coumarin - CO2 (styrene type ion)",
        intensity_rank="major", structural_origin="内酯环脱CO2→苯乙烯离子",
        fragmentation_path="内酯开环→脱CO2(44Da)", diagnostic_value=0.95,
    ),
    FragmentFeature(
        mz=91.0548, formula="C7H7+", name="莼离子(苯乙烯→䓬)",
        name_en="Tropylium (from styrene rearrangement)",
        intensity_rank="major", structural_origin="苯乙烯离子重排→䓬离子",
        fragmentation_path="103→重排→䓬(91Da)", diagnostic_value=0.90,
    ),
    FragmentFeature(
        mz=119.0497, formula="C8H7O+", name="香豆素脱CO",
        name_en="Coumarin - CO (benzofuran type)",
        intensity_rank="medium", structural_origin="内酯环脱CO→苯并呋喃型",
        fragmentation_path="内酯环→脱CO(28Da)", diagnostic_value=0.85,
    ),
    FragmentFeature(
        mz=131.0497, formula="C9H7O+", name="甲基香豆素脱CH4",
        name_en="Methylcoumarin - CH4",
        intensity_rank="medium", structural_origin="甲基香豆素脱甲烷",
        fragmentation_path="脱CH4→呋喃并苯", diagnostic_value=0.80,
    ),
)

# ═══════════════════════════════════════════════════════════════
# 6. 木脂素类 (Lignans) — 五味子等
# ═══════════════════════════════════════════════════════════════

LIGNAN_FRAGMENTS: Tuple[FragmentFeature, ...] = (
    FragmentFeature(
        mz=385.1651, formula="C22H25O6+", name="五味子甲素 [M+H]+",
        name_en="Schisandrin A [M+H]+",
        intensity_rank="base", structural_origin="联苯环辛烯型木脂素准分子离子",
        fragmentation_path="前体离子", diagnostic_value=1.0,
    ),
    FragmentFeature(
        mz=355.1545, formula="C21H23O5+", name="五味子甲素脱OCH2",
        name_en="Schisandrin A - OCH2 (demethylenedioxy)",
        intensity_rank="major", structural_origin="亚甲二氧基脱OCH2",
        fragmentation_path="脱OCH2→开环产物", diagnostic_value=0.90,
    ),
    FragmentFeature(
        mz=331.1545, formula="C20H21O4+", name="五味子甲素脱2xOCH3",
        name_en="Schisandrin A - 2xOCH3 (bisdemethoxy)",
        intensity_rank="medium", structural_origin="两甲氧基丢失",
        fragmentation_path="脱2xOCH3", diagnostic_value=0.85,
    ),
    FragmentFeature(
        mz=315.1596, formula="C19H23O4+", name="五味子醇甲脱H2O",
        name_en="Schisandrol A - H2O",
        intensity_rank="major", structural_origin="醇羟基脱水生成共轭二烯",
        fragmentation_path="醇脱水→共轭烯", diagnostic_value=0.87,
    ),
)

# ═══════════════════════════════════════════════════════════════
# 7. 醌类 (Quinones) — 丹参酮等
# ═══════════════════════════════════════════════════════════════

QUINONE_FRAGMENTS: Tuple[FragmentFeature, ...] = (
    FragmentFeature(
        mz=277.1223, formula="C19H17O2+", name="丹参酮IIA脱水碎片",
        name_en="Tanshinone IIA [M+H-H2O]+",
        intensity_rank="base", structural_origin="丹参酮IIA质子化分子脱水",
        fragmentation_path="[M+H]+→-H2O", diagnostic_value=1.0,
        references=("PMCID:PMC6268525", "PMCID:PMC6269062"),
        evidence_level="curated_primary",
    ),
    FragmentFeature(
        mz=249.1274, formula="C18H17O+", name="丹参酮IIA脱水脱CO碎片",
        name_en="Tanshinone IIA [M+H-H2O-CO]+",
        intensity_rank="major", structural_origin="脱水离子进一步失去醌羰基CO",
        fragmentation_path="295→-H2O→277→-CO→249", diagnostic_value=0.92,
        references=("PMCID:PMC6268525", "PMCID:PMC6269062"),
        evidence_level="curated_primary",
    ),
    FragmentFeature(
        mz=235.0754, formula="C16H11O2+", name="丹参酮IIA环重排碎片",
        name_en="Tanshinone IIA [277-C3H6]+",
        intensity_rank="medium", structural_origin="m/z 277非芳香A环重排并失去C3H6",
        fragmentation_path="277→-C3H6→235", diagnostic_value=0.85,
        references=("DOI:10.3389/fchem.2022.964744",),
        evidence_level="curated_primary",
    ),
    FragmentFeature(
        mz=221.1325, formula="C17H17+", name="丹参酮IIA多重脱羰",
        name_en="Tanshinone IIA - 2CO",
        intensity_rank="medium", structural_origin="两分子CO丢失",
        fragmentation_path="脱2×CO", diagnostic_value=0.88,
        references=("PMCID:PMC6269062", "PMCID:PMC7183602"),
        evidence_level="curated_primary",
    ),
)

# ═══════════════════════════════════════════════════════════════
# 8. 糖类 (Sugars/Monosaccharides) — 单糖、双糖及其裂解特征
# ═══════════════════════════════════════════════════════════════

SUGAR_FRAGMENTS: Tuple[FragmentFeature, ...] = (
    # ── 己糖 (Hexose) 特征碎片 ──
    FragmentFeature(
        mz=163.0606, formula="C6H11O5+", name="己糖氧碳鎓离子",
        name_en="Hexose oxocarbenium ion",
        intensity_rank="base", structural_origin="己糖(葡萄糖/半乳糖)糖苷键断裂或环裂解",
        fragmentation_path="己糖环裂解/C3-C4断裂", diagnostic_value=0.97, polarity="positive",
    ),
    FragmentFeature(
        mz=145.0500, formula="C6H9O4+", name="己糖脱H2O氧碳鎓离子",
        name_en="Hexose - H2O oxocarbenium ion",
        intensity_rank="major", structural_origin="己糖氧碳鎓离子脱水产物",
        fragmentation_path="己糖脱H2O", diagnostic_value=0.95, polarity="positive",
    ),
    FragmentFeature(
        mz=127.0395, formula="C6H7O3+", name="己糖脱2H2O碎片",
        name_en="Hexose - 2H2O fragment",
        intensity_rank="medium", structural_origin="己糖连续脱水",
        fragmentation_path="己糖→脱H2O→脱H2O", diagnostic_value=0.90, polarity="positive",
    ),
    FragmentFeature(
        mz=115.0395, formula="C5H7O3+", name="戊糖氧碳鎓离子",
        name_en="Pentose oxocarbenium ion (C5H7O3+)",
        intensity_rank="medium", structural_origin="己糖环C-C键断裂产生戊糖碎片",
        fragmentation_path="己糖环断裂→C5碎片", diagnostic_value=0.88, polarity="positive",
    ),
    FragmentFeature(
        mz=109.0290, formula="C6H5O2+", name="己糖脱3H2O碎片",
        name_en="Hexose - 3H2O fragment",
        intensity_rank="medium", structural_origin="己糖深度脱水(呋喃型碎片)",
        fragmentation_path="己糖→3xH2O", diagnostic_value=0.85, polarity="positive",
    ),
    FragmentFeature(
        mz=103.0395, formula="C4H7O3+", name="C4糖碎片",
        name_en="C4 sugar fragment",
        intensity_rank="minor", structural_origin="己糖环C2-C3断裂",
        fragmentation_path="己糖环断裂→C4碎片", diagnostic_value=0.80, polarity="positive",
    ),
    FragmentFeature(
        mz=85.0290, formula="C4H5O2+", name="丁糖碎片",
        name_en="Tetrose fragment",
        intensity_rank="minor", structural_origin="己糖环裂解后脱水",
        fragmentation_path="己糖环断裂→C4脱H2O", diagnostic_value=0.75, polarity="positive",
    ),
    FragmentFeature(
        mz=73.0290, formula="C3H5O2+", name="丙糖碎片",
        name_en="Triose fragment (C3H5O2+)",
        intensity_rank="minor", structural_origin="糖环裂解产生三碳碎片",
        fragmentation_path="己糖→丙糖碎片", diagnostic_value=0.70, polarity="positive",
    ),
    FragmentFeature(
        mz=71.0134, formula="C3H3O2+", name="丙糖醛酸碎片",
        name_en="Triose acid fragment",
        intensity_rank="minor", structural_origin="糖氧化裂解",
        fragmentation_path="糖→氧化裂解", diagnostic_value=0.65, polarity="positive",
    ),
    FragmentFeature(
        mz=60.0211, formula="C2H4O2+", name="乙酸碎片(糖重排)",
        name_en="Acetic acid fragment (sugar rearrangement)",
        intensity_rank="minor", structural_origin="糖的McLafferty重排产物",
        fragmentation_path="糖McLafferty重排", diagnostic_value=0.60, polarity="positive",
    ),

    # ── 脱氧己糖 (Deoxyhexose/Rhamnose) 特征碎片 ──
    FragmentFeature(
        mz=147.0657, formula="C6H11O4+", name="脱氧己糖氧碳鎓离子",
        name_en="Deoxyhexose oxocarbenium ion (Rhamnose/Fucose)",
        intensity_rank="base", structural_origin="鼠李糖/夫糖等脱氧糖的糖苷键断裂或环裂解",
        fragmentation_path="脱氧己糖环裂解", diagnostic_value=0.96, polarity="positive",
    ),
    FragmentFeature(
        mz=129.0552, formula="C6H9O3+", name="脱氧己糖脱H2O碎片",
        name_en="Deoxyhexose - H2O fragment",
        intensity_rank="major", structural_origin="脱氧己糖脱水产物",
        fragmentation_path="脱氧己糖→脱H2O", diagnostic_value=0.93, polarity="positive",
    ),
    FragmentFeature(
        mz=111.0445, formula="C6H7O2+", name="脱氧己糖脱2H2O碎片",
        name_en="Deoxyhexose - 2H2O fragment",
        intensity_rank="medium", structural_origin="脱氧己糖连续脱水",
        fragmentation_path="脱氧己糖→2xH2O", diagnostic_value=0.88, polarity="positive",
    ),
    FragmentFeature(
        mz=101.0603, formula="C5H9O2+", name="脱氧戊糖碎片",
        name_en="Deoxypentose fragment",
        intensity_rank="medium", structural_origin="脱氧己糖环断裂产生脱氧戊糖碎片",
        fragmentation_path="脱氧己糖→脱氧戊糖碎片", diagnostic_value=0.82, polarity="positive",
    ),
    FragmentFeature(
        mz=87.0446, formula="C4H7O2+", name="脱氧丁糖碎片",
        name_en="Deoxytetrose fragment",
        intensity_rank="minor", structural_origin="脱氧己糖环裂解",
        fragmentation_path="脱氧己糖→C4碎片", diagnostic_value=0.75, polarity="positive",
    ),

    # ── 戊糖 (Pentose) 特征碎片 ──
    FragmentFeature(
        mz=133.0501, formula="C5H9O4+", name="戊糖氧碳鎓离子",
        name_en="Pentose oxocarbenium ion",
        intensity_rank="base", structural_origin="阿拉伯糖/木糖等戊糖的环裂解",
        fragmentation_path="戊糖环裂解", diagnostic_value=0.94, polarity="positive",
    ),
    FragmentFeature(
        mz=115.0395, formula="C5H7O3+", name="戊糖脱H2O碎片",
        name_en="Pentose - H2O fragment",
        intensity_rank="major", structural_origin="戊糖脱水产物",
        fragmentation_path="戊糖→脱H2O", diagnostic_value=0.90, polarity="positive",
    ),
    FragmentFeature(
        mz=97.0290, formula="C5H5O2+", name="戊糖脱2H2O碎片",
        name_en="Pentose - 2H2O fragment",
        intensity_rank="medium", structural_origin="戊糖连续脱水",
        fragmentation_path="戊糖→2xH2O", diagnostic_value=0.85, polarity="positive",
    ),

    # ── 糖醇/还原糖特征 (负离子模式) ──
    FragmentFeature(
        mz=161.0453, formula="C6H9O5-", name="己糖去质子化碎片",
        name_en="Hexose deprotonated fragment [M-H]-",
        intensity_rank="major", structural_origin="己糖在负离子模式下的去质子化碎片",
        fragmentation_path="己糖去质子化", diagnostic_value=0.85, polarity="negative",
    ),
    FragmentFeature(
        mz=179.0561, formula="C6H11O6-", name="己糖酸去质子化",
        name_en="Hexonic acid [M-H]-",
        intensity_rank="major", structural_origin="己糖酸负离子",
        fragmentation_path="己糖酸去质子化", diagnostic_value=0.80, polarity="negative",
    ),
    FragmentFeature(
        mz=89.0244, formula="C3H5O3-", name="乳酸碎片(糖裂解)",
        name_en="Lactate fragment (sugar cleavage)",
        intensity_rank="minor", structural_origin="糖深度裂解产生的小分子酸",
        fragmentation_path="糖→乳酸碎片", diagnostic_value=0.65, polarity="negative",
    ),
    FragmentFeature(
        mz=71.0139, formula="C3H3O2-", name="丙酮酸碎片",
        name_en="Pyruvate fragment",
        intensity_rank="minor", structural_origin="糖氧化裂解",
        fragmentation_path="糖→丙酮酸碎片", diagnostic_value=0.60, polarity="negative",
    ),
)

# ═══════════════════════════════════════════════════════════════
# 9. 有机酸/脂肪酸 (Organic/Fatty Acids)
# ═══════════════════════════════════════════════════════════════

FATTY_ACID_FRAGMENTS: Tuple[FragmentFeature, ...] = (
    FragmentFeature(
        mz=60.0211, formula="CH3COOH2+", name="乙酸特征(rearr)",
        name_en="Acetic acid (McLafferty from butyric+)",
        intensity_rank="major", structural_origin="脂肪酸McLafferty重排→m/z 60",
        fragmentation_path="脂肪酸McLafferty重排→60Da", diagnostic_value=0.93,
    ),
    FragmentFeature(
        mz=73.0290, formula="C3H5O2+", name="丙酸特征",
        name_en="Propionic acid fragment",
        intensity_rank="medium", structural_origin="脂肪酸酯裂解",
        fragmentation_path="酯裂解→C3H5O2+", diagnostic_value=0.80,
    ),
    FragmentFeature(
        mz=88.0524, formula="C4H8O2+", name="丁酸McLafferty产物",
        name_en="Butyric acid McLafferty product",
        intensity_rank="major", structural_origin="丁酸酯McLafferty重排",
        fragmentation_path="McLafferty→C4H8O2+", diagnostic_value=0.88,
    ),
)

# ═══════════════════════════════════════════════════════════════
# 10. 萜类 (Terpenoids) — 单萜、倍半萜、二萜、三萜
# ═══════════════════════════════════════════════════════════════

TERPENOID_FRAGMENTS: Tuple[FragmentFeature, ...] = (
    # ── 单萜 (Monoterpenes, C10) ──
    FragmentFeature(
        mz=137.1328, formula="C10H17+", name="单萜特征离子(失去H2O/H)",
        name_en="Monoterpene characteristic ion (C10H17+)",
        intensity_rank="base", structural_origin="单萜骨架(柠檬烯/蒎烯)脱水或脱氢",
        fragmentation_path="单萜环裂解→C10H17+", diagnostic_value=0.95, polarity="positive",
    ),
    FragmentFeature(
        mz=121.1017, formula="C9H13+", name="单萜脱CH4碎片",
        name_en="Monoterpene - CH4 fragment",
        intensity_rank="major", structural_origin="单萜骨架脱甲基",
        fragmentation_path="单萜→脱CH4", diagnostic_value=0.85, polarity="positive",
    ),
    FragmentFeature(
        mz=93.0704, formula="C7H9+", name="环庚三烯碎片(单萜环裂)",
        name_en="Cycloheptatriene fragment (monoterpene ring cleavage)",
        intensity_rank="major", structural_origin="单萜七元环中间体",
        fragmentation_path="单萜→环裂解→C7碎片", diagnostic_value=0.88, polarity="positive",
    ),
    FragmentFeature(
        mz=81.0704, formula="C6H9+", name="环己烯碎片(萜类通用)",
        name_en="Cyclohexenyl fragment (universal terpene)",
        intensity_rank="base", structural_origin="萜类环结构断裂产物",
        fragmentation_path="萜环裂解→C6碎片", diagnostic_value=0.90, polarity="positive",
    ),

    # ── 倍半萜 (Sesquiterpenes, C15) ──
    FragmentFeature(
        mz=149.1328, formula="C11H17+", name="倍半萜特征离子",
        name_en="Sesquiterpene characteristic ion",
        intensity_rank="major", structural_origin="倍半萜骨架碎裂",
        fragmentation_path="倍半萜→骨架碎片", diagnostic_value=0.90, polarity="positive",
    ),
    FragmentFeature(
        mz=119.0861, formula="C9H11+", name="倍半萜脱C2H6碎片",
        name_en="Sesquiterpene - C2H6 fragment",
        intensity_rank="medium", structural_origin="倍半萜侧链丢失",
        fragmentation_path="倍半萜→侧链丢失", diagnostic_value=0.82, polarity="positive",
    ),

    # ── 二萜 (Diterpenes, C20) ──
    FragmentFeature(
        mz=191.1790, formula="C14H23+", name="二萜骨架特征碎片",
        name_en="Diterpene skeleton fragment",
        intensity_rank="major", structural_origin="二萜(丹参酮/穿心莲内酯)骨架裂解",
        fragmentation_path="二萜骨架→C14碎片", diagnostic_value=0.90, polarity="positive",
    ),
    FragmentFeature(
        mz=173.1684, formula="C14H21+", name="二萜脱水碎片",
        name_en="Diterpene - H2O fragment",
        intensity_rank="medium", structural_origin="二萜骨架脱水",
        fragmentation_path="二萜→脱H2O", diagnostic_value=0.82, polarity="positive",
    ),

    # ── 三萜 (Triterpenes, C30) ──
    FragmentFeature(
        mz=203.1800, formula="C15H23+", name="三萜四环骨架(达玛烷型)",
        name_en="Triterpenoid tetracyclic core (dammarane-type)",
        intensity_rank="major", structural_origin="达玛烷型三萜骨架(失去全部糖和侧链)",
        fragmentation_path="完全脱糖→骨架离子", diagnostic_value=0.87, polarity="positive",
    ),
    FragmentFeature(
        mz=189.1643, formula="C14H21+", name="三萜D/E环碎片",
        name_en="Triterpenoid D/E ring fragment",
        intensity_rank="medium", structural_origin="达玛烷型三萜D环侧链裂解",
        fragmentation_path="侧链裂解→D/E环", diagnostic_value=0.85, polarity="positive",
    ),
    FragmentFeature(
        mz=109.1017, formula="C8H13+", name="萜类链碎片(C8H13+)",
        name_en="Terpenoid chain fragment C8H13+",
        intensity_rank="major", structural_origin="萜类侧链或环裂解碎片",
        fragmentation_path="萜类→C8碎片", diagnostic_value=0.80, polarity="positive",
    ),
    FragmentFeature(
        mz=95.0861, formula="C7H11+", name="萜类环碎片(C7H11+)",
        name_en="Terpenoid ring fragment C7H11+",
        intensity_rank="medium", structural_origin="倍半萜环裂解",
        fragmentation_path="萜环裂解→C7碎片", diagnostic_value=0.78, polarity="positive",
    ),
)

# ═══════════════════════════════════════════════════════════════
# 11. 甾体类 (Steroids) — 甾体骨架、甾体皂苷元
# ═══════════════════════════════════════════════════════════════

STEROID_FRAGMENTS: Tuple[FragmentFeature, ...] = (
    # ── 甾体骨架碎片 (Steroid skeleton) ──
    FragmentFeature(
        mz=255.2113, formula="C19H27+", name="甾体D环碎片(脱侧链)",
        name_en="Steroid D-ring fragment (side chain loss)",
        intensity_rank="base", structural_origin="甾体D环侧链(如C8H17)丢失后的环系碎片",
        fragmentation_path="甾体→侧链丢失→D环", diagnostic_value=0.96, polarity="positive",
    ),
    FragmentFeature(
        mz=271.2063, formula="C19H27O+", name="甾体D环碎片(含3-OH)",
        name_en="Steroid D-ring fragment with 3-OH",
        intensity_rank="major", structural_origin="含3-羟基甾体侧链丢失",
        fragmentation_path="甾体→侧链丢失→3-OH-D环", diagnostic_value=0.94, polarity="positive",
    ),
    FragmentFeature(
        mz=273.2219, formula="C19H29O+", name="甾体饱和D环碎片",
        name_en="Steroid saturated D-ring fragment",
        intensity_rank="major", structural_origin="甾体饱和环系碎片",
        fragmentation_path="甾体→侧链丢失→饱和D环", diagnostic_value=0.92, polarity="positive",
    ),
    FragmentFeature(
        mz=213.1643, formula="C16H21+", name="甾体A/B/C环碎片",
        name_en="Steroid A/B/C ring fragment",
        intensity_rank="medium", structural_origin="甾体A/B/C三环碎片",
        fragmentation_path="甾体D环断裂→A/B/C环", diagnostic_value=0.85, polarity="positive",
    ),

    # ── 薯蓣皂苷元 (Diosgenin, 甾体皂苷元) ──
    FragmentFeature(
        mz=415.3212, formula="C27H43O3+", name="薯蓣皂苷元离子(25R)",
        name_en="Diosgenin ion (25R configuration)",
        intensity_rank="base", structural_origin="25R-薯蓣皂苷脱全部糖",
        fragmentation_path="完全脱糖→25R苷元", diagnostic_value=0.96, polarity="positive",
    ),
    FragmentFeature(
        mz=397.3107, formula="C27H41O2+", name="薯蓣皂苷元脱水离子",
        name_en="Diosgenin dehydration ion",
        intensity_rank="major", structural_origin="苷元3-OH位脱水",
        fragmentation_path="苷元→-H2O", diagnostic_value=0.85, polarity="positive",
    ),
    FragmentFeature(
        mz=282.2168, formula="C19H30O+", name="薯蓣皂苷元E/F环裂解",
        name_en="Diosgenin E/F ring opening",
        intensity_rank="major", structural_origin="螺甾烷E/F开环",
        fragmentation_path="E/F环开环→碎片", diagnostic_value=0.89, polarity="positive",
    ),

    # ── 通用甾体中性丢失后碎片 ──
    FragmentFeature(
        mz=147.1170, formula="C11H15+", name="甾体环碎片(C11H15+)",
        name_en="Steroid ring fragment C11H15+",
        intensity_rank="medium", structural_origin="甾体多环裂解产物",
        fragmentation_path="甾体→环裂解→C11碎片", diagnostic_value=0.75, polarity="positive",
    ),
)

# ═══════════════════════════════════════════════════════════════
# 12. 氨基酸/肽类 (Amino Acids/Peptides)
# ═══════════════════════════════════════════════════════════════

AMINO_ACID_FRAGMENTS: Tuple[FragmentFeature, ...] = (
    # ── 亚铵离子 (Immonium ions, AA鉴定特征) ──
    FragmentFeature(
        mz=86.0964, formula="C5H12N+", name="亮氨酸/异亮氨酸亚铵离子",
        name_en="Leu/Ile immonium ion",
        intensity_rank="base", structural_origin="亮氨酸/异亮氨酸侧链特征",
        fragmentation_path="AA→亚铵离子(α-β断裂)", diagnostic_value=0.96, polarity="positive",
        required_smarts=(
            "CC(C)CC(N)C(=O)[O;H1,-1]",
            "CCC(C)C(N)C(=O)[O;H1,-1]",
        ),
        allowed_adducts=("[M+H]+",), references=("PMID:18564857",),
        evidence_level="curated_primary",
    ),
    FragmentFeature(
        mz=120.0813, formula="C8H10N+", name="苯丙氨酸亚铵离子",
        name_en="Phe immonium ion",
        intensity_rank="base", structural_origin="苯丙氨酸苄基侧链特征",
        fragmentation_path="Phe→亚铵离子", diagnostic_value=0.97, polarity="positive",
        required_smarts=("N[C](Cc1ccccc1)C(=O)[O;H1,-1]",),
        allowed_adducts=("[M+H]+",), references=("PMID:18564857",),
        evidence_level="curated_primary",
    ),
    FragmentFeature(
        mz=136.0762, formula="C8H10NO+", name="酪氨酸亚铵离子",
        name_en="Tyr immonium ion",
        intensity_rank="base", structural_origin="酪氨酸羟基苄基特征",
        fragmentation_path="Tyr→亚铵离子(+O)", diagnostic_value=0.96, polarity="positive",
        required_smarts=("N[C](Cc1ccc([O;H1,-1])cc1)C(=O)[O;H1,-1]",),
        allowed_adducts=("[M+H]+",), references=("PMID:18564857",),
        evidence_level="curated_primary",
    ),
    FragmentFeature(
        mz=110.0718, formula="C5H8N3+", name="组氨酸亚铵离子(咪唑环)",
        name_en="His immonium ion (imidazole)",
        intensity_rank="major", structural_origin="组氨酸咪唑环特征",
        fragmentation_path="His→咪唑亚铵", diagnostic_value=0.94, polarity="positive",
        required_smarts=("N[C](Cc1cnc[nH]1)C(=O)[O;H1,-1]",),
        allowed_adducts=("[M+H]+",), references=("PMID:18564857",),
        evidence_level="curated_primary",
    ),
    FragmentFeature(
        mz=159.0922, formula="C10H11N2+", name="色氨酸亚铵离子",
        name_en="Trp immonium ion",
        intensity_rank="major", structural_origin="色氨酸吲哚环特征",
        fragmentation_path="Trp→吲哚亚铵离子", diagnostic_value=0.93, polarity="positive",
        required_smarts=("N[C](Cc1c[nH]c2ccccc12)C(=O)[O;H1,-1]",),
        allowed_adducts=("[M+H]+",), references=("PMID:18564857",),
        evidence_level="curated_primary",
    ),
    FragmentFeature(
        mz=129.1140, formula="C5H13N4+", name="精氨酸亚铵离子",
        name_en="Arg immonium ion (guanidino)",
        intensity_rank="major", structural_origin="精氨酸胍基特征",
        fragmentation_path="Arg→胍基亚铵", diagnostic_value=0.92, polarity="positive",
        required_smarts=("N[C](CCCNC(=N)N)C(=O)[O;H1,-1]",),
        allowed_adducts=("[M+H]+",), references=("PMID:18564857",),
        evidence_level="curated_primary",
    ),
    FragmentFeature(
        mz=101.1079, formula="C5H13N2+", name="赖氨酸亚铵离子",
        name_en="Lys immonium ion",
        intensity_rank="major", structural_origin="赖氨酸侧链亚铵离子",
        fragmentation_path="Lys→侧链亚铵", diagnostic_value=0.91, polarity="positive",
        required_smarts=("N[C](CCCCN)C(=O)[O;H1,-1]",),
        allowed_adducts=("[M+H]+",), references=("PMID:18564857",),
        evidence_level="curated_primary",
    ),
    FragmentFeature(
        mz=72.0813, formula="C4H10N+", name="缬氨酸亚铵离子",
        name_en="Val immonium ion",
        intensity_rank="major", structural_origin="缬氨酸异丙基侧链亚铵离子",
        fragmentation_path="Val→亚铵离子", diagnostic_value=0.93, polarity="positive",
        required_smarts=("CC(C)C(N)C(=O)[O;H1,-1]",),
        allowed_adducts=("[M+H]+",), references=("PMID:18564857",),
        evidence_level="curated_primary",
    ),
    FragmentFeature(
        mz=70.0656, formula="C4H8N+", name="脯氨酸亚铵离子",
        name_en="Pro immonium ion (pyrrolidine)",
        intensity_rank="major", structural_origin="脯氨酸吡咯烷环特征",
        fragmentation_path="Pro→吡咯亚铵", diagnostic_value=0.90, polarity="positive",
        required_smarts=("O=C([O;H1,-1])C1CCCN1",),
        allowed_adducts=("[M+H]+",), references=("PMID:18564857",),
        evidence_level="curated_primary",
    ),
    FragmentFeature(
        mz=74.0242, formula="C2H4NO2-", name="甘氨酸碎片(负离子)",
        name_en="Glycine fragment [M-H]-",
        intensity_rank="major", structural_origin="甘氨酸脱羧特征碎片",
        fragmentation_path="Gly→去质子化离子", diagnostic_value=0.85, polarity="negative",
        required_smarts=("[N;X3,X4+]-[CH2;D2]-[C](=O)[O;H1,-1]",),
        allowed_adducts=("[M-H]-",), evidence_level="curated_measured",
    ),

    # ── 肽键断裂特征 ──
    FragmentFeature(
        mz=175.1190, formula="C7H15N4O+", name="精氨酸侧链(b/y离子特征)",
        name_en="Arg side chain (b/y-ion marker)",
        intensity_rank="major", structural_origin="肽段中精氨酸残基特征",
        fragmentation_path="肽键断裂→Arg侧链", diagnostic_value=0.88, polarity="positive",
    ),
)

# ═══════════════════════════════════════════════════════════════
# 13. 核苷/脱氧核苷 (Nucleosides)
# ═══════════════════════════════════════════════════════════════

NUCLEOSIDE_FRAGMENTS: Tuple[FragmentFeature, ...] = (
    FragmentFeature(
        mz=112.0505, formula="C4H6N3O+", name="质子化胞嘧啶",
        name_en="Protonated cytosine", intensity_rank="base",
        structural_origin="胞苷N-糖苷键断裂后的碱基离子",
        fragmentation_path="胞苷→失核糖/脱氧核糖→[cytosine+H]+",
        diagnostic_value=0.98, polarity="positive",
        required_smarts=("Nc1ccn([C;R])c(=O)n1",),
        allowed_adducts=("[M+H]+",),
        references=("DOI:10.1021/acs.chemrestox.3c00041", "PMCID:PMC11587168"),
        evidence_level="curated_primary",
    ),
    FragmentFeature(
        mz=113.0346, formula="C4H5N2O2+", name="质子化尿嘧啶",
        name_en="Protonated uracil", intensity_rank="base",
        structural_origin="尿苷N-糖苷键断裂后的碱基离子",
        fragmentation_path="尿苷→失核糖→[uracil+H]+",
        diagnostic_value=0.98, polarity="positive",
        required_smarts=("O=c1ccn([C;R])c(=O)[nH]1",),
        allowed_adducts=("[M+H]+",),
        references=("DOI:10.1021/acs.chemrestox.3c00041", "PMCID:PMC11587168"),
        evidence_level="curated_primary",
    ),
    FragmentFeature(
        mz=127.0502, formula="C5H7N2O2+", name="质子化胸腺嘧啶",
        name_en="Protonated thymine", intensity_rank="base",
        structural_origin="胸苷N-糖苷键断裂后的碱基离子",
        fragmentation_path="胸苷→失脱氧核糖→[thymine+H]+",
        diagnostic_value=0.98, polarity="positive",
        required_smarts=("Cc1cn([C;R])c(=O)[nH]c1=O",),
        allowed_adducts=("[M+H]+",), references=("PMCID:PMC11587168",),
        evidence_level="curated_primary",
    ),
    FragmentFeature(
        mz=136.0618, formula="C5H6N5+", name="质子化腺嘌呤",
        name_en="Protonated adenine", intensity_rank="base",
        structural_origin="腺苷N-糖苷键断裂后的碱基离子",
        fragmentation_path="腺苷→失核糖/脱氧核糖→[adenine+H]+",
        diagnostic_value=0.99, polarity="positive",
        required_smarts=("Nc1ncnc2n([C;R])cnc12",),
        allowed_adducts=("[M+H]+",),
        references=("DOI:10.1021/acs.chemrestox.3c00041", "PMCID:PMC7249026"),
        evidence_level="curated_primary",
    ),
    FragmentFeature(
        mz=152.0567, formula="C5H6N5O+", name="质子化鸟嘌呤",
        name_en="Protonated guanine", intensity_rank="base",
        structural_origin="鸟苷N-糖苷键断裂后的碱基离子",
        fragmentation_path="鸟苷→失核糖/脱氧核糖→[guanine+H]+",
        diagnostic_value=0.99, polarity="positive",
        required_smarts=("Nc1nc2c(ncn2[C;R])c(=O)[nH]1",),
        allowed_adducts=("[M+H]+",),
        references=("DOI:10.1021/acs.chemrestox.3c00041", "PMCID:PMC7249026"),
        evidence_level="curated_primary",
    ),
)

# ═══════════════════════════════════════════════════════════════
# 14. 脂质类 (Lipids) — 甘油磷脂、鞘脂、甘油酯
# ═══════════════════════════════════════════════════════════════

LIPID_FRAGMENTS: Tuple[FragmentFeature, ...] = (
    # ── 磷脂头基碎片 (Phospholipid head groups) ──
    FragmentFeature(
        mz=184.0739, formula="C5H15NO4P+", name="磷脂酰胆碱头基(PC)",
        name_en="Phosphocholine head group [M+H]+",
        intensity_rank="base", structural_origin="PC/SM类磷脂的极性头基",
        fragmentation_path="PC→头基断裂→184Da", diagnostic_value=0.98, polarity="positive",
        compatible_class_terms=("phosphatidylcholine", "lysophosphatidylcholine", "sphingomyelin"),
        required_smarts=(
            "P(=O)(O)OCC[N+](C)(C)C",
            "P(=O)([O-])OCC[N+](C)(C)C",
            "P(=O)(O)OCCN(C)(C)C",
        ),
        allowed_adducts=("[M+H]+",),
        references=("PMID:32423109", "PMCID:PMC7287777"),
        evidence_level="curated_primary",
    ),
    FragmentFeature(
        mz=104.1075, formula="C5H14NO+", name="胆碱碎片",
        name_en="Choline fragment",
        intensity_rank="major", structural_origin="磷脂胆碱基团特征",
        fragmentation_path="PC→胆碱碎片", diagnostic_value=0.90, polarity="positive",
        compatible_class_terms=("phosphatidylcholine", "lysophosphatidylcholine", "sphingomyelin"),
        required_smarts=(
            "P(=O)(O)OCC[N+](C)(C)C",
            "P(=O)([O-])OCC[N+](C)(C)C",
            "P(=O)(O)OCCN(C)(C)C",
        ),
        allowed_adducts=("[M+H]+",),
        references=("PMCID:PMC12113132",), evidence_level="curated_primary",
    ),
    FragmentFeature(
        mz=86.0964, formula="C5H12N+", name="脱水胆碱碎片",
        name_en="Dehydrated choline fragment (注意:与Leu/Ile亚铵重叠)",
        intensity_rank="major", structural_origin="胆碱脱H2O",
        fragmentation_path="胆碱→脱H2O", diagnostic_value=0.70, polarity="positive",
        compatible_class_terms=("phosphatidylcholine", "lysophosphatidylcholine", "sphingomyelin"),
        required_smarts=(
            "P(=O)(O)OCC[N+](C)(C)C",
            "P(=O)([O-])OCC[N+](C)(C)C",
            "P(=O)(O)OCCN(C)(C)C",
        ),
        allowed_adducts=("[M+H]+",),
        references=("PMCID:PMC12113132",), evidence_level="curated_primary",
    ),
    FragmentFeature(
        mz=264.2686, formula="C18H34N+", name="d18:1鞘氨醇双脱水离子",
        name_en="d18:1 sphingoid-base double-dehydration ion",
        intensity_rank="major", structural_origin="d18:1鞘氨醇长链碱骨架",
        fragmentation_path="ceramide/sphingolipid→long-chain base→2H2O loss",
        diagnostic_value=0.95, polarity="positive",
        compatible_class_terms=("sphingolipid", "ceramide", "sphingomyelin", "hexosylceramide"),
        required_smarts=("NCC(O)C=C", "NC(CO)C(O)C=C"),
        allowed_adducts=("[M+H]+",),
        references=("PMID:11070861", "PMCID:PMC2724058"),
        evidence_level="curated_primary",
    ),

    # ── 负离子磷脂类别头基 ──
    FragmentFeature(
        mz=224.0693, formula="C7H15NO5P-", name="PC/LPC甘油磷酸胆碱负离子",
        name_en="Glycerophosphocholine anion (PC/LPC)",
        intensity_rank="major", structural_origin="PC/LPC去甲基化极性头基",
        fragmentation_path="PC/LPC→demethylated glycerophosphocholine",
        diagnostic_value=0.95, polarity="negative",
        compatible_class_terms=("phosphatidylcholine", "lysophosphatidylcholine"),
        required_smarts=(
            "P(=O)(O)OCC[N+](C)(C)C",
            "P(=O)([O-])OCC[N+](C)(C)C",
        ),
        allowed_adducts=("[M+FA-H]-", "[M+Hac-H]-"),
        references=("PMID:32423109", "PMCID:PMC7287777"),
        evidence_level="curated_primary",
    ),
    FragmentFeature(
        mz=196.0380, formula="C5H11NO5P-", name="PE/LPE磷酸乙醇胺负离子",
        name_en="Glycerophosphoethanolamine anion (PE/LPE)",
        intensity_rank="major", structural_origin="PE/LPE磷酸乙醇胺极性头基",
        fragmentation_path="PE/LPE→glycerophosphoethanolamine",
        diagnostic_value=0.95, polarity="negative",
        compatible_class_terms=("phosphatidylethanolamine", "lysophosphatidylethanolamine"),
        required_smarts=("P(=O)(O)OCC[NX3;H1,H2]", "P(=O)([O-])OCC[NX3;H1,H2]"),
        allowed_adducts=("[M-H]-",),
        references=("PMID:32423109", "DOI:10.1002/rcm.9527"),
        evidence_level="curated_primary",
    ),
    FragmentFeature(
        mz=241.0117, formula="C6H10O8P-", name="PI肌醇磷酸负离子",
        name_en="Inositol phosphate anion (PI)",
        intensity_rank="major", structural_origin="PI肌醇磷酸极性头基",
        fragmentation_path="PI→inositol phosphate",
        diagnostic_value=0.95, polarity="negative",
        compatible_class_terms=("phosphatidylinositol", "lysophosphatidylinositol"),
        required_smarts=(
            "P(=O)(O)OC1C(O)C(O)C(O)C(O)C1O",
            "P(=O)([O-])OC1C(O)C(O)C(O)C(O)C1O",
        ),
        allowed_adducts=("[M-H]-",),
        references=("PMID:32423109", "PMCID:PMC7287777"),
        evidence_level="curated_primary",
    ),

    # ── 脂肪酸链碎片 (负离子模式) ──
    FragmentFeature(
        mz=171.1386, formula="C9H15O2-", name="脂肪酸C9碎片(CRF)",
        name_en="Fatty acid C9 fragment (CRF)",
        intensity_rank="major", structural_origin="长链脂肪酸电荷远程裂解",
        fragmentation_path="CRF→C9碎片", diagnostic_value=0.85, polarity="negative",
    ),
    FragmentFeature(
        mz=255.2329, formula="C16H31O2-", name="棕榈酸根离子 [M-H]-",
        name_en="Palmitate [M-H]- (C16:0)",
        intensity_rank="base", structural_origin="棕榈酸去质子化",
        fragmentation_path="Palmitate去质子化", diagnostic_value=0.95, polarity="negative",
    ),
    FragmentFeature(
        mz=281.2486, formula="C18H33O2-", name="油酸根离子 [M-H]-",
        name_en="Oleate [M-H]- (C18:1)",
        intensity_rank="base", structural_origin="油酸去质子化",
        fragmentation_path="Oleate去质子化", diagnostic_value=0.95, polarity="negative",
    ),
    FragmentFeature(
        mz=279.2329, formula="C18H31O2-", name="亚油酸根离子 [M-H]-",
        name_en="Linoleate [M-H]- (C18:2)",
        intensity_rank="base", structural_origin="亚油酸去质子化",
        fragmentation_path="Linoleate去质子化", diagnostic_value=0.95, polarity="negative",
    ),

    # ── 甘油骨架碎片 ──
    FragmentFeature(
        mz=152.9963, formula="C3H6O5P-", name="甘油磷酸碎片",
        name_en="Glycerophosphate fragment",
        intensity_rank="major", structural_origin="磷脂中甘油磷酸骨架",
        fragmentation_path="甘油磷酸断裂", diagnostic_value=0.88, polarity="negative",
    ),
)

# ═══════════════════════════════════════════════════════════════
# 14. 苯环类 (Benzenoids) — 简单芳香化合物
# ═══════════════════════════════════════════════════════════════

BENZENOID_FRAGMENTS: Tuple[FragmentFeature, ...] = (
    # ── 苯环骨架碎片 ──
    FragmentFeature(
        mz=77.0391, formula="C6H5+", name="苯基阳离子(Phenyl)",
        name_en="Phenyl cation",
        intensity_rank="base", structural_origin="苯环直接离子化或䓬离子脱乙炔",
        fragmentation_path="苯环→C6H5+", diagnostic_value=0.90, polarity="positive",
    ),
    FragmentFeature(
        mz=91.0548, formula="C7H7+", name="䓬离子(Tropylium)",
        name_en="Tropylium ion",
        intensity_rank="base", structural_origin="苄基裂解后重排产生的七元环芳香阳离子",
        fragmentation_path="苄基裂解→重排→䓬", diagnostic_value=0.92, polarity="positive",
    ),
    FragmentFeature(
        mz=93.0340, formula="C6H5O+", name="苯氧基阳离子(Phenoxy)",
        name_en="Phenoxy cation",
        intensity_rank="major", structural_origin="酚类苄基位裂解",
        fragmentation_path="酚→苯氧碎片", diagnostic_value=0.90, polarity="positive",
    ),
    FragmentFeature(
        mz=105.0335, formula="C7H5O+", name="苯甲酰基阳离子(Benzoyl)",
        name_en="Benzoyl cation (C6H5CO+)",
        intensity_rank="major", structural_origin="苯甲酸酯或芳香酮的α-裂解",
        fragmentation_path="芳香酮→α裂解→苯甲酰", diagnostic_value=0.93, polarity="positive",
    ),
    FragmentFeature(
        mz=107.0491, formula="C7H7O+", name="羟基苄基阳离子(Hydroxybenzyl)",
        name_en="Hydroxybenzyl cation",
        intensity_rank="major", structural_origin="酚类/苄醇类的苄基位裂解",
        fragmentation_path="苄醇→苄基裂解", diagnostic_value=0.88, polarity="positive",
    ),

    # ── 取代苯特征 ──
    FragmentFeature(
        mz=121.0290, formula="C7H5O2+", name="苯甲酸特征碎片",
        name_en="Benzoic acid fragment (苯甲酸特征)",
        intensity_rank="major", structural_origin="苯甲酸/水杨酸衍生物",
        fragmentation_path="苯甲酸→CO2丢失前体", diagnostic_value=0.85, polarity="positive",
    ),
    FragmentFeature(
        mz=79.0548, formula="C6H7+", name="环己二烯阳离子",
        name_en="Cyclohexadienyl cation",
        intensity_rank="medium", structural_origin="苯环部分氢化碎片",
        fragmentation_path="苯环→部分还原", diagnostic_value=0.72, polarity="positive",
    ),
    FragmentFeature(
        mz=65.0391, formula="C5H5+", name="环戊二烯基阳离子",
        name_en="Cyclopentadienyl cation",
        intensity_rank="minor", structural_origin="苯环深度碎裂",
        fragmentation_path="苯环→C2H2丢失→C5H5+", diagnostic_value=0.65, polarity="positive",
    ),
)

# ═══════════════════════════════════════════════════════════════
#

ADDUCT_PROFILES: Dict[str, AdductProfile] = {
    "[M+H]+": AdductProfile(
        adduct_name="[M+H]+", mass_shift=1.0078, polarity="positive",
        prevalence=0.98, typical_compounds="含碱性基团(N)的化合物(生物碱、胺类)",
        fragmentation_behavior="质子化促进诱导裂解和重排(质子移动模型)",
    ),
    "[M+Na]+": AdductProfile(
        adduct_name="[M+Na]+", mass_shift=22.9898, polarity="positive",
        prevalence=0.85, typical_compounds="多OH化合物(糖类、黄酮苷、皂苷)",
        fragmentation_behavior="Na+定位促进特定位置裂解; 通常需要更高CID能量",
    ),
    "[M+K]+": AdductProfile(
        adduct_name="[M+K]+", mass_shift=38.9637, polarity="positive",
        prevalence=0.40, typical_compounds="含多醚/多OH的大环化合物",
        fragmentation_behavior="类似[M+Na]+但K+结合更强,碎片较少",
    ),
    "[M+NH4]+": AdductProfile(
        adduct_name="[M+NH4]+", mass_shift=18.0344, polarity="positive",
        prevalence=0.55, typical_compounds="醛类、含羰基化合物",
        fragmentation_behavior="NH4+易丢失NH3(17Da), 产生[M+H]+等价离子",
    ),
    "[M]+": AdductProfile(
        adduct_name="[M]+", mass_shift=0.0, polarity="positive",
        prevalence=0.15, typical_compounds="季铵盐类生物碱(小檗碱型)、自由基阳离子",
        fragmentation_behavior="无质子化位点, 裂解由电荷位置驱动",
    ),
    "[M-H]-": AdductProfile(
        adduct_name="[M-H]-", mass_shift=-1.0078, polarity="negative",
        prevalence=0.90, typical_compounds="含酸性基团(COOH, 酚OH)化合物",
        fragmentation_behavior="去质子化促进电荷远程裂解(CRF); 黄酮C环rDA在负离子模式下也有特征",
    ),
    "[M+Cl]-": AdductProfile(
        adduct_name="[M+Cl]-", mass_shift=34.9689, polarity="negative",
        prevalence=0.25, typical_compounds="含多卤素或强酸性化合物",
        fragmentation_behavior="Cl-加合后中性丢失HCl(36Da)常见",
    ),
    "[M+HCOOH-H]-": AdductProfile(
        adduct_name="[M+HCOOH-H]-", mass_shift=44.9977, polarity="negative",
        prevalence=0.35, typical_compounds="弱酸性化合物(动/植物代谢物)",
        fragmentation_behavior="甲酸加合物→易丢失HCOOH(46Da)",
    ),
    "[M+CH3COOH-H]-": AdductProfile(
        adduct_name="[M+CH3COOH-H]-", mass_shift=59.0133, polarity="negative",
        prevalence=0.20, typical_compounds="黄酮、酚类等弱酸性天然产物",
        fragmentation_behavior="乙酸加合物→易丢失CH3COOH(60Da)",
    ),
}

# ═══════════════════════════════════════════════════════════════
# 第三部分: 化合物类别完整档案 (Compound Class Profiles)
# ═══════════════════════════════════════════════════════════════

COMPOUND_CLASS_PROFILES: Dict[str, CompoundClassProfile] = {
    "flavonoid": CompoundClassProfile(
        class_name="flavonoid", class_name_cn="黄酮类",
        parent_class="polyphenol",
        characteristic_fragments=FLAVONOID_FRAGMENTS,
        diagnostic_ion_pairs=(
            (152.0473, 120.0575),   # A环(5,7-diOH) + B环(无OH)
            (153.0188, 137.0239),   # A环(5,7-diOH-4-keto) + B环(3',4'-diOH) 槲皮素
            (165.0188, 153.0188),   # A环(多OH) + B环(triOH) 杨梅素
            (152.0473, 137.0239),   # A环(5,7-diOH) + B环(diOH) 木犀草素
        ),
        fragmentation_pathway=(
            "糖苷键顺序断裂(O-苷) → 脱全部糖基 → 苷元离子 → "
            "C环rDA裂解 → A环1,3A+离子 + B环1,3B+离子 → "
            "后续小分子丢失(H2O/CO/CO2/CH3)"
        ),
        sequential_losses=(
            "C6H10O5(己糖162Da)", "C6H10O4(脱氧己糖146Da)",
            "C5H8O4(戊糖132Da)", "C6H8O6(葡萄糖醛酸176Da)",
            "CH3(15Da)", "H2O(18Da)", "CO(28Da)", "CO2(44Da)",
        ),
        isomer_differentiation_rules=(
            "黄酮 vs 异黄酮: B环rDA碎片m/z不同(黄酮B环连接C2,异黄酮连接C3)",
            "3-OH vs 无3-OH: A环1,3A+离子差16Da(152→153, 含4-酮基)",
            "4'-OH vs 3',4'-diOH: B环1,3B+离子差16Da(121→137)",
            "O-苷 vs C-苷: C-苷不发生糖苷键直接断裂,出现糖环交叉裂解(120/90Da)",
            "1→2 vs 1→6二糖连接: 中间碎片m/z不同(1→2可见单糖苷碎片309Da)",
        ),
        reference_compounds=(
            "槲皮素(Quercetin)", "木犀草素(Luteolin)", "山奈酚(Kaempferol)",
            "芹菜素(Apigenin)", "杨梅素(Myricetin)", "黄芩素(Baicalein)",
            "大豆苷元(Daidzein)", "葛根素(Puerarin)", "染料木素(Genistein)",
        ),
        literature_refs=(
            "Cuyckens & Claeys, JMS 2004, 39, 1-15",
            "Ma et al., JASMS 1997, 8, 1180-1187",
            "Fabre et al., JASMS 2001, 12, 707-715",
            "刘淑莹等, 中药质谱分析, 第4章, 2012",
        ),
        confidence_score=0.97,
    ),
    "glycoside": CompoundClassProfile(
        class_name="glycoside", class_name_cn="通用糖苷类",
        parent_class="glycoconjugate",
        # 通用糖苷只能安全复用糖残基/氧碳鎓离子证据；苷元碎片必须交给
        # 更具体的黄酮、皂苷、香豆素等档案或候选专属经验通路解释。
        characteristic_fragments=SUGAR_FRAGMENTS,
        diagnostic_ion_pairs=(),
        fragmentation_pathway=(
            "糖苷键断裂 → 糖残基/氧碳鎓离子与苷元离子 → "
            "苷元按其自身骨架继续裂解"
        ),
        sequential_losses=(
            "C6H10O5(己糖162Da)",
            "C6H10O4(脱氧己糖146Da)",
            "C5H8O4(戊糖132Da)",
            "C6H8O6(葡萄糖醛酸176Da)",
        ),
        isomer_differentiation_rules=(
            "O-苷与C-苷需结合直接脱糖和糖环交叉裂解证据判断",
            "糖连接位点和端基构型通常不能仅凭通用糖残基离子确认",
        ),
        reference_compounds=(
            "Aucubin",
            "Arbutin",
            "Salicin",
            "Flavonoid O-glycosides",
        ),
        literature_refs=(
            "Cuyckens & Claeys, J Mass Spectrom 2004, 39, 1-15",
            "刘淑莹等, 中药质谱分析, 糖苷裂解章节, 2012",
        ),
        confidence_score=0.88,
    ),
    "saponin": CompoundClassProfile(
        class_name="saponin", class_name_cn="皂苷类",
        parent_class="terpenoid_glycoside",
        characteristic_fragments=SAPONIN_FRAGMENTS,
        diagnostic_ion_pairs=(
            (162.0528, 146.0579),  # 己糖 + 脱氧己糖 (糖链组成)
            (109.1017, 189.1643),  # 三萜C8碎片 + D/E环碎片 (骨架确认)
            (203.1800, 189.1643),  # 四环核心 vs D/E环 (达玛烷 vs 齐墩果烷)
            (132.0423, 162.0528),  # 戊糖 + 己糖 (糖链组合)
        ),
        fragmentation_pathway=(
            "糖链顺序丢失(从外端开始) → 逐级脱糖(己糖162/脱氧己糖146/戊糖132/GlcA176) → "
            "苷元离子 → 苷元骨架裂解 → "
            "三萜: A/B/C/D环碎片 → 侧链丢失(C3H6/C4H8) → "
            "甾体皂苷的苷元碎片见 steroid 类别档案"
        ),
        sequential_losses=(
            "C6H10O5(己糖162Da)", "C6H10O4(脱氧己糖146Da)",
            "C5H8O4(戊糖132Da)", "C6H8O6(GlcA176Da)",
            "H2O(18Da)", "CH3COOH(60Da)", "CH3OH(32Da)",
        ),
        isomer_differentiation_rules=(
            "25R vs 25S(薯蓣皂苷): 高分辨MS区分—苷元脱水离子丰度不同",
            "四环三萜 vs 五环三萜: 骨架碎片m/z分布差异显著",
            "达玛烷 vs 齐墩果烷: E环碎片m/z 203 vs 189",
            "糖基数确定: 母离子-[M+H]+与苷元离子质量差除以糖单位质量",
            "糖连接顺序: 逐级MSn,Yn-1离子丰度推理连接位点",
        ),
        reference_compounds=(
            "人参皂苷Rb1/Rg1/Re(Ginsenosides)", "薯蓣皂苷(Dioscin)",
            "三七皂苷R1(Notoginsenoside R1)", "绞股蓝皂苷(Gypenosides)",
            "黄芪甲苷(Astragaloside IV)", "柴胡皂苷(Saikosaponins)",
        ),
        literature_refs=(
            "刘淑莹等, 中药质谱分析, 第3章, 2012",
            "Fang et al., JMS 2019, 54, 587-599",
            "Li et al., Phytochem Anal 2020, 31, 45-57",
        ),
        confidence_score=0.94,
    ),
    "alkaloid": CompoundClassProfile(
        class_name="alkaloid", class_name_cn="生物碱类",
        parent_class="nitrogenous_compound",
        characteristic_fragments=ALKALOID_FRAGMENTS,
        diagnostic_ion_pairs=(
            (586.3016, 526.2808),  # 乌头碱→脱乙酸
            (526.2808, 494.2543),  # 去乙酰乌头碱→脱甲醇
            (336.1236, 321.1001),  # 小檗碱→脱CH3
            (138.0667, 110.0600),  # 咖啡因特征碎片对
            (150.1283, 136.1126),  # 喹嗪啶特征碎片对
        ),
        fragmentation_pathway=(
            "乌头碱型: [M+H]+ → 脱乙酸(60Da) → 脱甲醇(32Da) → 脱苯甲酸(122Da) → 骨架裂解 → "
            "原小檗碱型: [M]+ → 脱CH3(15Da) → 脱CH4(16Da) → 脱2xCH3 → 脱CO(28Da) → "
            "黄嘌呤型: [M+H]+ → 咪唑rDA→脱CH3NCO(57Da) → 脱CO(28Da) → 脱HCN(27Da)"
        ),
        sequential_losses=(
            "CH3COOH(乙酸60Da)", "CH3OH(甲醇32Da)",
            "C6H5COOH(苯甲酸122Da)", "CH3(甲基15Da)",
            "CH3NCO(异氰酸甲酯57Da)", "HNCO(异氰酸43Da)",
            "CO(28Da)", "HCN(氢氰酸27Da)", "CH4(甲烷16Da)",
        ),
        isomer_differentiation_rules=(
            "乌头碱 vs 次乌头碱: C3位OH→H差异(差16Da)导致碎片模式偏移",
            "乌头碱 vs 中乌头碱: C3位OH vs OCH3→脱CH3 vs 脱CH3OH优先",
            "小檗碱 vs 巴马汀: 亚甲二氧基 vs 甲氧基→脱CH2O(30Da) vs 脱CH3(15Da)差异",
            "原小檗碱 vs 二氢原小檗碱: [M]+ vs [M+H]+ 加合模式不同",
        ),
        reference_compounds=(
            "乌头碱(Aconitine)", "次乌头碱(Hypaconitine)",
            "中乌头碱(Mesaconitine)", "苯甲酰乌头原碱(Benzoylaconine)",
            "小檗碱(Berberine)", "巴马汀(Palmatine)", "黄连碱(Coptisine)",
            "咖啡因(Caffeine)", "苦参碱(Matrine)",
        ),
        literature_refs=(
            "刘淑莹等, 中药质谱分析, 第5章, 2012",
            "Yue et al., JMS 2009, 44, 1379-1385",
            "Wang et al., J Pharm Biomed Anal 2010, 52, 438-445",
            "台湾质谱学会, 质谱分析技术原理与应用, 第2章, 2019",
        ),
        confidence_score=0.95,
    ),
    "phenolic_acid": CompoundClassProfile(
        class_name="phenolic_acid", class_name_cn="酚酸类",
        parent_class="polyphenol",
        characteristic_fragments=PHENOLIC_ACID_FRAGMENTS,
        diagnostic_ion_pairs=(
            (163.0395, 119.0497),  # 对香豆酸→脱CO2
            (179.0344, 135.0446),  # 咖啡酸→脱CO2
            (197.0450, 179.0344),  # 丹参素→脱H2O
            (147.0446, 103.0548),  # 肉桂酸→脱CO2
        ),
        fragmentation_pathway=(
            "酚酸 [M-H]- or [M+H]+ → 脱CO2(44Da, 最特征) → "
            "脱H2O(18Da, 邻位效应) → 脱CO(28Da) → "
            "脱HCOOH(46Da, 邻位双OH/COOH) → 苯环骨架离子"
        ),
        sequential_losses=(
            "CO2(44Da, 最特征)", "H2O(18Da)", "CO(28Da)",
            "HCOOH(甲酸46Da)", "CH3(15Da)",
        ),
        isomer_differentiation_rules=(
            "咖啡酸 vs 异咖啡酸: 脱水难易不同(邻位vs对位)",
            "丹参素 vs 迷迭香酸: 特征碎片m/z 197 vs 359",
            "绿原酸异构体: 奎宁酸部分连接位点不同→MS3裂解路径不同",
        ),
        reference_compounds=(
            "咖啡酸(Caffeic acid)", "阿魏酸(Ferulic acid)",
            "迷迭香酸(Rosmarinic acid)", "丹酚酸B(Salvianolic acid B)",
            "绿原酸(Chlorogenic acid)", "没食子酸(Gallic acid)",
        ),
        literature_refs=(
            "刘淑莹等, 中药质谱分析, 第6章, 2012",
            "Clifford et al., J Sci Food Agric 2003, 83, 1097-1104",
        ),
        confidence_score=0.90,
    ),
    "coumarin": CompoundClassProfile(
        class_name="coumarin", class_name_cn="香豆素类",
        parent_class="polyphenol",
        characteristic_fragments=COUMARIN_FRAGMENTS,
        diagnostic_ion_pairs=(
            (147.0446, 103.0548),  # 香豆素→脱CO2
            (103.0548, 91.0548),   # 苯乙烯→莼离子
        ),
        fragmentation_pathway=(
            "香豆素 [M+H]+ → 内酯环脱CO2(44Da, 最特征) → "
            "苯乙烯型离子(103Da) → 重排→䓬离子(91Da) → "
            "脱CO(28Da) → 脱C2H2(26Da)"
        ),
        sequential_losses=("CO2(44Da)", "CO(28Da)", "CH4(16Da)", "C2H2(26Da)"),
        isomer_differentiation_rules=(
            "简单香豆素 vs 呋喃香豆素: 呋喃香豆素先脱呋喃环片段(42/56Da)",
            "6,7- vs 7,8-呋喃香豆素: 呋喃环开环路径不同→中间碎片差异",
            "香豆素木脂素: 联苯片段特征碎片m/z 137/152",
        ),
        reference_compounds=(
            "香豆素(Coumarin)", "伞形花内酯(Umbelliferone)",
            "补骨脂素(Psoralen)", "欧前胡素(Imperatorin)",
            "白花前胡甲素(Praeruptorin A)",
        ),
        literature_refs=(
            "Tine et al., J Nat Prod 2017, 80, 1239-1247",
            "台湾质谱学会, 质谱分析技术原理与应用, 2019",
        ),
        confidence_score=0.88,
    ),
    "lignan": CompoundClassProfile(
        class_name="lignan", class_name_cn="木脂素类",
        parent_class="polyphenol",
        characteristic_fragments=LIGNAN_FRAGMENTS,
        diagnostic_ion_pairs=(
            (385.1651, 355.1545),  # 五味子甲素→脱OCH2
            (355.1545, 331.1545),  # 脱OCH2→脱OCH3
        ),
        fragmentation_pathway=(
            "木脂素 [M+H]+ → 脱OCH2(30Da,亚甲二氧基) → "
            "脱OCH3(31Da,甲氧基) → 脱H2O(18Da,醇OH) → "
            "联苯键断裂→单体碎片"
        ),
        sequential_losses=("OCH2(30Da)", "OCH3(31Da)", "H2O(18Da)", "CH3(15Da)"),
        isomer_differentiation_rules=(
            "联苯环辛烯 vs 二苄基丁烷型: 骨架结构不同→裂解模式迥异",
            "五味子甲素 vs 五味子乙素: 亚甲二氧基数不同→脱OCH2个数差异",
        ),
        reference_compounds=(
            "五味子甲素(Schisandrin A)", "五味子乙素(Schisandrin B)",
            "五味子醇甲(Schisandrol A)", "厚朴酚(Magnolol)",
            "和厚朴酚(Honokiol)",
        ),
        literature_refs=(
            "刘淑莹等, 中药质谱分析, 第6章, 2012",
            "Huang et al., J Pharm Biomed Anal 2008, 46, 274-281",
        ),
        confidence_score=0.85,
    ),
    "quinone": CompoundClassProfile(
        class_name="quinone", class_name_cn="醌类",
        parent_class="polyketide",
        characteristic_fragments=QUINONE_FRAGMENTS,
        diagnostic_ion_pairs=(
            (277.1434, 249.1485),  # 丹参酮IIA→脱CO
            (249.1485, 235.1481),  # 脱CO→脱CH2
        ),
        fragmentation_pathway=(
            "二萜醌 [M+H]+ → 脱CO(28Da, 邻醌特征) → "
            "脱CH2(14Da) → 脱H2O(18Da) → 脱2xCO → 萘环骨架"
        ),
        sequential_losses=("CO(28Da)", "CH2(14Da)", "H2O(18Da)", "CH3(15Da)"),
        isomer_differentiation_rules=(
            "丹参酮IIA vs 隐丹参酮: 后者有呋喃环→额外脱CO差异",
            "邻醌 vs 对醌: 脱CO难易不同(邻醌易脱一个CO,对醌更难)",
        ),
        reference_compounds=(
            "丹参酮IIA(Tanshinone IIA)", "隐丹参酮(Cryptotanshinone)",
            "丹参酮I(Tanshinone I)", "紫草素(Shikonin)",
            "大黄素(Emodin)", "大黄酸(Rhein)",
        ),
        literature_refs=(
            "刘淑莹等, 中药质谱分析, 第6章, 2012",
            "Yang et al., J Pharm Biomed Anal 2006, 41, 1042-1048",
        ),
        confidence_score=0.87,
    ),
    "sugar": CompoundClassProfile(
        class_name="sugar", class_name_cn="糖类/单糖",
        parent_class="carbohydrate",
        characteristic_fragments=SUGAR_FRAGMENTS,
        diagnostic_ion_pairs=(
            (163.0606, 145.0500),  # 己糖 → 己糖脱H2O
            (147.0657, 129.0552),  # 脱氧己糖 → 脱氧己糖脱H2O
            (133.0501, 115.0395),  # 戊糖 → 戊糖脱H2O
        ),
        fragmentation_pathway=(
            "糖环裂解(交叉环断裂) → 氧碳鎓离子(特征m/z 163/147/133) → "
            "连续脱水(H2O, 18Da) → 碎片离子(C3-C5碎片) → "
            "小分子碎片(乳酸m/z 89, 丙酮酸m/z 71) → "
            "脱氧糖: 特征脱氧己糖碎片m/z 147(比己糖163少16Da)"
        ),
        sequential_losses=(
            "H2O(水18Da, 最常见)", "H2O(连续脱水, 糖类最显著特征)",
            "CH2O(甲醛30Da)", "C2H4O2(乙酸60Da)",
            "C3H6O3(丙三醇片段90Da)", "CO(28Da)", "CO2(44Da)",
        ),
        isomer_differentiation_rules=(
            "葡萄糖 vs 半乳糖: 碎片离子丰度分布不同,葡萄糖m/z 163丰度更高",
            "鼠李糖 vs 岩藻糖: 立体异构体,MS/MS需结合保留时间区分",
            "己糖 vs 脱氧己糖: 特征离子m/z 163 vs 147(差16Da,脱氧特征)",
            "戊糖 vs 己糖: 特征离子m/z 133 vs 163(差30Da,少一个CH2O)",
            "吡喃糖 vs 呋喃糖: 环裂解碎片分布不同",
        ),
        reference_compounds=(
            "葡萄糖(Glucose)", "半乳糖(Galactose)", "甘露糖(Mannose)",
            "鼠李糖(Rhamnose)", "岩藻糖(Fucose)", "阿拉伯糖(Arabinose)",
            "木糖(Xylose)", "果糖(Fructose)", "核糖(Ribose)",
        ),
        literature_refs=(
            "台湾质谱学会, 质谱分析技术原理与应用, 第3章, 2019",
            "Domon & Costello, Glycobiology 1988, 8, 557-568",
            "Harvey, D.J., JASMS 2000, 11, 900-915",
            "刘淑莹等, 中药质谱分析, 第4章 糖苷裂解, 2012",
        ),
        confidence_score=0.93,
    ),
    "fatty_acid": CompoundClassProfile(
        class_name="fatty_acid", class_name_cn="脂肪酸/有机酸",
        parent_class="organic_acid",
        characteristic_fragments=FATTY_ACID_FRAGMENTS,
        diagnostic_ion_pairs=((60.0211, 73.0290),),
        fragmentation_pathway=(
            "脂肪酸 [M-H]- → 脱CO2(44Da) → "
            "电荷远程裂解(CRF) → CH2=CH2序列丢失(28Da间隔) → "
            "不饱和脂肪酸: 烯丙基裂解优先"
        ),
        sequential_losses=("CO2(44Da)", "H2O(18Da)", "C2H4(28Da)xN"),
        isomer_differentiation_rules=(
            "双键位置: CRF系列碎片分布模式确定不饱和位点",
            "支链 vs 直链: 支链位CRF碎片强度异常",
        ),
        reference_compounds=(
            "棕榈酸(Palmitic acid)", "油酸(Oleic acid)",
            "亚油酸(Linoleic acid)", "α-亚麻酸(α-Linolenic acid)",
        ),
        literature_refs=(
            "Gross, JASMS 1992, 3, 105-116",
            "Murphy & Axelsen, Mass Spectrom Rev 2011, 30, 579-601",
        ),
        confidence_score=0.82,
    ),
    # ═══════════════════════════════════════════════════════════════
    # 新增5个类别档案 (v3.3.0)
    # ═══════════════════════════════════════════════════════════════
    "terpenoid": CompoundClassProfile(
        class_name="terpenoid", class_name_cn="萜类",
        parent_class="isoprenoid",
        characteristic_fragments=TERPENOID_FRAGMENTS,
        diagnostic_ion_pairs=(
            (137.1328, 81.0704),   # 单萜: C10H17+ + 环己烯碎片
            (149.1328, 93.0704),   # 倍半萜: C11H17+ + 环庚三烯
            (191.1790, 173.1684),  # 二萜: C14H23+ + 脱水碎片
            (203.1800, 189.1643),  # 三萜: 四环核心 + D/E环
        ),
        fragmentation_pathway=(
            "萜类环裂解 → 烯丙基断裂(优先) → "
            "异戊二烯单元丢失(68Da) → "
            "单萜→m/z 137; 倍半萜→m/z 149; 二萜→m/z 191; 三萜→m/z 203"
        ),
        sequential_losses=(
            "C5H8(异戊二烯68Da)", "C3H6(丙烯42Da)",
            "C4H8(丁烯56Da)", "H2O(18Da)", "CO(28Da)",
        ),
        isomer_differentiation_rules=(
            "单萜 vs 倍半萜: 特征离子m/z 137 vs 149(差12Da, C1)",
            "倍半萜 vs 二萜: 特征离子m/z 149 vs 191(差42Da, C3H6)",
            "达玛烷 vs 齐墩果烷: E环碎片m/z 203 vs 189",
            "含氧萜 vs 非含氧萜: 脱H2O碎片的存在性",
        ),
        reference_compounds=(
            "柠檬烯(Limonene)", "α-蒎烯(α-Pinene)", "青蒿素(Artemisinin)",
            "穿心莲内酯(Andrographolide)", "人参皂苷苷元(Protopanaxadiol)",
        ),
        literature_refs=(
            "台湾质谱学会, 质谱分析技术原理与应用, 第4章, 2019",
            "Demarque et al., Phytochem Anal 2016, 27, 129-138",
        ),
        confidence_score=0.85,
    ),
    "steroid": CompoundClassProfile(
        class_name="steroid", class_name_cn="甾体类",
        parent_class="terpenoid",
        characteristic_fragments=STEROID_FRAGMENTS,
        diagnostic_ion_pairs=(
            (255.2113, 271.2063),  # D环(无OH) vs D环(含3-OH)
            (415.3212, 397.3107),  # 薯蓣皂苷元 vs 脱水(25R/25S区分)
            (273.2219, 213.1643),  # 饱和D环 vs A/B/C环
        ),
        fragmentation_pathway=(
            "甾体骨架裂解 → 侧链丢失(丁烯56/异丁烯) → "
            "D环碎片(m/z 255/271/273) → A/B/C环碎片(m/z 213) → "
            "螺甾烷: E/F环开环(m/z 282) → 苷元脱水(m/z 397)"
        ),
        sequential_losses=(
            "C8H16(辛烯112Da)", "C3H6(丙烯42Da)",
            "C4H8(丁烯56Da)", "H2O(18Da)", "CO(28Da)",
        ),
        isomer_differentiation_rules=(
            "25R vs 25S(薯蓣皂苷元): 苷元脱水离子丰度差异",
            "胆甾烷 vs 豆甾烷: 侧链质量差异(C2H4=28Da)",
            "3-OH vs 3-酮: D环碎片m/z 271 vs 255(差16Da)",
        ),
        reference_compounds=(
            "薯蓣皂苷元(Diosgenin)", "胆甾醇(Cholesterol)",
            "豆甾醇(Stigmasterol)", "β-谷甾醇(β-Sitosterol)",
        ),
        literature_refs=(
            "刘淑莹等, 中药质谱分析, 第3章, 2012",
            "Fang et al., JMS 2019, 54, 587-599",
        ),
        confidence_score=0.88,
    ),
    "amino_acid": CompoundClassProfile(
        class_name="amino_acid", class_name_cn="氨基酸/肽类",
        parent_class="nitrogenous_compound",
        characteristic_fragments=AMINO_ACID_FRAGMENTS,
        diagnostic_ion_pairs=(
            (86.0964, 110.0718),   # Leu/Ile vs His 亚铵离子
            (120.0813, 136.0762),  # Phe vs Tyr 亚铵离子
            (159.0922, 129.1140),  # Trp vs Arg 亚铵离子
            (101.1079, 72.0813),   # Lys vs Val 亚铵离子
        ),
        fragmentation_pathway=(
            "肽键断裂(b/y离子) → 氨基酸残基亚铵离子形成 → "
            "侧链特征碎片(每个AA有独特亚铵离子) → "
            "脱NH3(17Da) → 脱CO2(44Da) → 脱H2O(18Da)"
        ),
        sequential_losses=(
            "H2O(18Da)", "CO2(44Da)", "NH3(17Da)",
            "HCOOH(甲酸46Da)", "CO(28Da)",
        ),
        isomer_differentiation_rules=(
            "亮氨酸 vs 异亮氨酸: 亚铵离子m/z相同(86), 需保留时间区分",
            "苯丙氨酸 vs 酪氨酸: 亚铵离子差16Da(120 vs 136, 苄基 vs 羟苄基)",
            "精氨酸 vs 赖氨酸: Arg亚铵m/z 129.1140, Lys亚铵m/z 101.1079",
        ),
        reference_compounds=(
            "甘氨酸(Glycine)", "亮氨酸(Leucine)", "苯丙氨酸(Phenylalanine)",
            "色氨酸(Tryptophan)", "精氨酸(Arginine)", "组氨酸(Histidine)",
        ),
        literature_refs=(
            "Hohmann et al., Anal Chem 2008, PMID:18564857",
            "J Mass Spectrom 2024, PMID:39563088",
            "Steen & Mann, Nat Rev Mol Cell Biol 2004, 5, 699-711",
        ),
        confidence_score=0.94,
    ),
    "nucleoside": CompoundClassProfile(
        class_name="nucleoside", class_name_cn="核苷/脱氧核苷",
        parent_class="nucleic_acid_component",
        characteristic_fragments=NUCLEOSIDE_FRAGMENTS,
        diagnostic_ion_pairs=(),
        fragmentation_pathway=(
            "质子化核苷N-糖苷键断裂 → 核糖核苷失132.0423 Da / "
            "脱氧核苷失116.0473 Da → 形成对应质子化核碱基"
        ),
        sequential_losses=(
            "C5H8O4(核糖132.0423Da)",
            "C5H8O3(脱氧核糖116.0473Da)",
        ),
        isomer_differentiation_rules=(
            "核糖核苷 vs 脱氧核苷: N-糖苷键裂解分别失132.0423 vs 116.0473 Da",
            "碱基身份由质子化胞嘧啶/尿嘧啶/胸腺嘧啶/腺嘌呤/鸟嘌呤精确质量约束",
        ),
        reference_compounds=(
            "腺苷(Adenosine)", "鸟苷(Guanosine)", "胞苷(Cytidine)",
            "尿苷(Uridine)", "胸苷(Thymidine)",
        ),
        literature_refs=(
            "Chem Res Toxicol 2023, DOI:10.1021/acs.chemrestox.3c00041",
            "PMCID:PMC7249026", "PMCID:PMC11587168",
        ),
        confidence_score=0.97,
    ),
    "lipid": CompoundClassProfile(
        class_name="lipid", class_name_cn="脂质/磷脂类",
        parent_class="lipid",
        characteristic_fragments=LIPID_FRAGMENTS,
        diagnostic_ion_pairs=(
            (184.0739, 104.1075),  # PC头基 + 胆碱碎片
            (255.2329, 281.2486),  # 棕榈酸根 + 油酸根 (FA链组成)
            (279.2329, 281.2486),  # 亚油酸根 + 油酸根 (不饱和度区分)
        ),
        fragmentation_pathway=(
            "磷脂头基断裂(正离子) → PC: m/z 184(磷脂酰胆碱头基) → "
            "脂肪酸链断裂(负离子) → CRF系列碎片(14Da间隔) → "
            "sn-1 vs sn-2 脂肪酸链区分(中性丢失差异)"
        ),
        sequential_losses=(
            "C2H4(28Da)xN(链长)", "H2O(18Da)", "CO2(44Da)",
            "R-COOH(脂肪酸链)", "C3H9N(三甲胺59Da)",
        ),
        isomer_differentiation_rules=(
            "PC vs PE: 头基碎片m/z 184(PC) vs m/z 196(PE-H2O)",
            "sn-1 vs sn-2 脂肪酸: 中性丢失丰度比差异",
            "饱和 vs 不饱和: CRF碎片间隔模式不同",
        ),
        reference_compounds=(
            "PC(16:0/18:1)", "PE(18:0/20:4)", "SM(d18:1/16:0)",
            "TG(16:0/18:1/18:2)", "棕榈酸(Palmitic acid)",
        ),
        literature_refs=(
            "Gross, JASMS 1992, 3, 105-116",
            "Murphy & Axelsen, Mass Spectrom Rev 2011, 30, 579-601",
        ),
        confidence_score=0.85,
    ),
    "benzenoid": CompoundClassProfile(
        class_name="benzenoid", class_name_cn="苯环类",
        parent_class="aromatic_compound",
        characteristic_fragments=BENZENOID_FRAGMENTS,
        diagnostic_ion_pairs=(
            (77.0391, 91.0548),    # 苯基 + 䓬离子
            (93.0340, 105.0335),   # 苯氧 + 苯甲酰
            (107.0491, 121.0290),  # 羟基苄基 + 苯甲酸碎片
        ),
        fragmentation_pathway=(
            "苄基裂解(优先) → 䓬离子(m/z 91) → 脱乙炔(C2H2, 26Da) → "
            "苯基阳离子(m/z 77) → 环戊二烯基(m/z 65) → "
            "取代苯: 取代基特征碎片(苯甲酰105, 苯氧93等)"
        ),
        sequential_losses=(
            "C2H2(乙炔26Da)", "CO(28Da)", "H2O(18Da)",
            "CO2(44Da)", "CH3(15Da)",
        ),
        isomer_differentiation_rules=(
            "邻位 vs 对位取代: 邻位效应(脱H2O/CO2)仅在邻位发生",
            "甲氧基 vs 羟基: 脱CH3(15Da) vs 脱H2O(18Da)",
            "苯甲醇 vs 苯甲醛: m/z 107(苄基+O) vs 105(苯甲酰)",
        ),
        reference_compounds=(
            "苯(Benzene)", "甲苯(Toluene)", "苯酚(Phenol)",
            "苯甲酸(Benzoic acid)", "水杨酸(Salicylic acid)",
        ),
        literature_refs=(
            "McLafferty, Interpretation of Mass Spectra, 4th Ed., Ch.3",
            "台湾质谱学会, 质谱分析技术原理与应用, 第4章, 2019",
        ),
        confidence_score=0.80,
    ),
}

# ═══════════════════════════════════════════════════════════════
# 第四部分: 中药炮制/配伍特征 (TCM Processing & Compatibility)
# ═══════════════════════════════════════════════════════════════

TCM_PROCESSING_FEATURES: Dict[str, Any] = {
    "aconitum_processing": {
        "description": "乌头类(附子/川乌/草乌)炮制过程中的化学变化",
        "raw_features": {
            "diester_diterpenoid_alkaloids": {
                "compounds": ["乌头碱", "次乌头碱", "中乌头碱"],
                "characteristic": "[M+H]+ ± 1 Da, 观察乙酸丢失(60Da)和苯甲酸丢失(122Da)",
                "toxicity": "高毒性(双酯型)",
            },
        },
        "processed_features": {
            "monoester_diterpenoid_alkaloids": {
                "compounds": ["苯甲酰乌头原碱", "苯甲酰次乌头原碱"],
                "characteristic": "C8位脱乙酸→单酯型, 仍可见苯甲酰离子m/z 105",
                "toxicity": "中低毒性(单酯型)",
            },
            "amino_alcohol_diterpenoid_alkaloids": {
                "compounds": ["乌头原碱", "次乌头原碱"],
                "characteristic": "完全脱酰基→胺醇型, 不见m/z 105和乙酸丢失",
                "toxicity": "极低毒性(胺醇型)",
            },
        },
        "processing_markers": {
            "diester→monoester": "质量减少60Da(脱乙酸), 碎片模式中苯甲酰离子保留",
            "monoester→amino_alcohol": "质量减少122Da(脱苯甲酸), m/z 105消失",
            "Maillard_reaction": "红参加工中出现Melanoidin产物(棕色), 特征吸收及m/z改变",
        },
    },
    "ginseng_processing": {
        "description": "人参→红参加工过程中的Maillard反应",
        "maillard_stage1": {
            "reaction": "还原糖 + 氨基化合物 → Schiff碱 → Amadori重排",
            "ms_features": "Amadori产物特征: [M+H]+ 含C6H10O5加成, 水合茚三酮反应阳性减弱",
        },
        "maillard_stage2": {
            "reaction": "Amadori产物 → 脱水 → 裂解 → 羟甲基糠醛(HMF)等",
            "ms_features": "HMF: m/z 127.0395 [M+H]+, 紫外284nm; 还原酮类碎片",
        },
        "ginsenoside_changes": {
            "decrease": "丙二酰基-人参皂苷(malonyl-ginsenosides)含量下降→脱丙二酰(86Da)",
            "increase": "脱酰基人参皂苷增加; Rh系列皂苷(次级苷)增加",
            "new_compounds": "梅拉德反应产物-人参皂苷加合物可被高分辨MS检测",
        },
    },
    "compatibility_chemistry": {
        "aconitum_licorice": {
            "formula": "附子-甘草配伍(四逆汤关键药对)",
            "chemical_change": "甘草酸与乌头碱形成复合物→降低乌头碱溶出率",
            "ms_evidence": "配伍后乌头碱类[M+H]+峰强度降低50-80%",
        },
        "aconitum_ginger": {
            "formula": "附子-干姜配伍",
            "chemical_change": "干姜促进乌头碱脱乙酰→加速解毒",
            "ms_evidence": "单酯型和胺醇型峰比例增加,双酯型峰降低",
        },
    },
}

# ═══════════════════════════════════════════════════════════════
# 第五部分: 质谱数据质量控制标准 (Spectraverse-inspired)
# ═══════════════════════════════════════════════════════════════

DATA_QUALITY_STANDARDS: Dict[str, Any] = {
    "spectrum_level": {
        "min_fragment_ions": 3,       # 最少碎片离子数
        "max_fragment_ions": 4096,    # 最多保留碎片数(top-N)
        "min_relative_intensity": 0.001,  # 最小相对强度(0.1%)
        "min_mz": 10.0,               # 最小m/z
        "max_mz": 1000.0,             # 最大m/z
        "precursor_ppm_tolerance": 10.0,  # 前体质量精度阈值(ppm)
        "require_high_resolution": True,   # 要求高分辨数据
    },
    "metadata_level": {
        "required_fields": [
            "precursor_mz", "ionization_mode", "adduct",
            "collision_energy", "smiles", "inchikey",
        ],
        "recommended_fields": [
            "instrument_type", "instrument_name", "compound_name",
            "formula", "retention_time", "publication_doi",
        ],
        "adduct_validation": "parent_mass_from_smiles must match precursor_mz - adduct_mass within 0.1Da",
    },
    "chemical_structure_level": {
        "smiles_must_be_valid": True,
        "remove_stereochemistry": True,
        "standardize_tautomers": True,
        "neutralize_charges_except": ["zwitterions", "quaternary_ammonium"],
        "remove_disconnected_fragments": True,
        "remove_metals": True,
    },
    "duplicate_detection": {
        "inchikey_match_length": 14,  # 前14位匹配
        "cosine_threshold": 0.99,     # 余弦阈值
        "require_same_polarity": True,
        "keep_most_informative": True,  # 保留碎片数最多的
    },
    "quality_flags": {
        "UNIFORM_INTENSITIES": "去除所有碎片强度相同的伪谱",
        "ZERO_INTENSITY_FRAGMENTS": "去除零强度碎片",
        "FRAGMENTS_ABOVE_PRECURSOR": "去除所有碎片m/z > 前体m/z的谱图",
        "SINGLE_FRAGMENT_NEAR_PRECURSOR": "去除仅有一个碎片且在前体±1.6Da内的谱图",
        "LOW_RESOLUTION_QQQ": "去除三重四极杆低分辨数据",
        "MSN_ABOVE_MS2": "去除MS3及以上谱图",
        "INCOHERENT_ANNOTATION": "去除母离子质量与SMILES计算质量不匹配的谱图",
    },
}

# ═══════════════════════════════════════════════════════════════
# 第六部分: 结构解析推理引擎辅助数据
# ═══════════════════════════════════════════════════════════════

STRUCTURE_ELUCIDATION_RULES: Dict[str, Any] = {
    "substructure_diagnosis": {
        "benzene_ring": {
            "evidence": ["m/z 77 (C6H5+)", "m/z 91 (C7H7+ tropylium)",
                        "系列间隔26Da丢失(C2H2)"],
            "confidence": 0.95,
        },
        "phenolic_oh": {
            "evidence": ["脱H2O(18Da)碎片", "m/z 93 (C6H5O+)",
                        "负离子模式响应强", "邻位效应(脱H2O+CO2)"],
            "confidence": 0.90,
        },
        "methoxy": {
            "evidence": ["脱CH3(15Da)", "脱CH3OH(32Da)",
                        "脱CH2O(30Da)", "同位素峰M+2较明显(含O)"],
            "confidence": 0.88,
        },
        "methylenedioxy": {
            "evidence": ["脱OCH2(30Da)", "脱CH2O(30Da) 与甲氧基区分",
                        "特征m/z 135/147含OCH2O碎片"],
            "confidence": 0.85,
        },
        "carboxyl": {
            "evidence": ["脱CO2(44Da, 最特征)", "脱HCOOH(46Da)",
                        "负离子模式[M-H]-强响应", "脱H2O(18Da) 邻位效应"],
            "confidence": 0.92,
        },
        "carbonyl": {
            "evidence": ["α-裂解碎片", "McLafferty重排产物(diff=28/42/56Da)",
                        "脱CO(28Da)", "酰基阳离子RCO+"],
            "confidence": 0.90,
        },
        "amino": {
            "evidence": ["α-裂解→亚胺离子m/z 30", "脱NH3(17Da)",
                        "含氮碎片偶数质量(氮规则)", "正离子模式响应强"],
            "confidence": 0.93,
        },
        "glycoside": {
            "evidence": ["脱162Da(己糖)", "脱146Da(脱氧己糖)", "脱132Da(戊糖)",
                        "糖基氧碳鎓离子m/z 163/147/133", "逐级脱糖模式"],
            "confidence": 0.96,
        },
    },
    "mass_spectral_rules": {
        "nitrogen_rule": "含奇数N→分子量为奇数(EI); 含偶数N→分子量为偶数",
        "ring_double_bond_rule": "不饱和度 = (2C + 2 + N - H - X) / 2",
        "isotope_patterns": {
            "Cl": "M:M+2 ≈ 3:1",
            "Br": "M:M+2 ≈ 1:1",
            "S": "M+2/M ≈ 4.4%",
            "Si": "M+2/M ≈ 3.4%",
        },
    },
}

# ─────────────────────────────────────────────────────────────────
# 查询与检索接口
# ─────────────────────────────────────────────────────────────────


def get_profile(class_name: str) -> Optional[CompoundClassProfile]:
    """获取化合物类别的完整质谱特征档案"""
    return COMPOUND_CLASS_PROFILES.get(class_name.lower())


def get_fragments_by_class(
    class_name: str,
    polarity: str = "",
    validated_only: bool = True,
    candidate_smiles: str = "",
    raw_classification: str = "",
    enforce_scope: bool = False,
    query_adduct: str = "",
    query_collision_energy: Any = None,
) -> Tuple[FragmentFeature, ...]:
    """获取特定类别的可适用碎片。

    Production callers default to excluding records whose formula and m/z fail
    the 0.01 Da quality gate.  The raw profile remains available for audit with
    ``validated_only=False``.
    """
    profile = get_profile(class_name)
    if not profile:
        return ()
    wanted_polarity = normalize_polarity(polarity)
    fragments = []
    for fragment in profile.characteristic_fragments:
        fragment_polarity = normalize_polarity(fragment.polarity)
        if wanted_polarity and fragment_polarity and fragment_polarity != wanted_polarity:
            continue
        if validated_only and not formula_mz_consistent(fragment.formula, fragment.mz):
            continue
        if enforce_scope and (
            fragment.compatible_class_terms
            or fragment.required_smarts
            or fragment.required_smarts_all
        ):
            raw_lower = str(raw_classification or "").lower()
            class_match = bool(
                fragment.compatible_class_terms
                and any(term.lower() in raw_lower for term in fragment.compatible_class_terms)
            )
            structure_match = bool(
                fragment.required_smarts
                and smiles_matches_any(str(candidate_smiles or ""), fragment.required_smarts)
            )
            structure_all_match = bool(
                fragment.required_smarts_all
                and smiles_matches_all(
                    str(candidate_smiles or ""), fragment.required_smarts_all
                )
            )
            if fragment.required_smarts_all:
                structure_match = structure_match and structure_all_match if fragment.required_smarts else structure_all_match
            if not (class_match or structure_match):
                continue
        if enforce_scope and fragment.allowed_adducts:
            observed_adduct = normalize_adduct(query_adduct)
            allowed_adducts = {normalize_adduct(item) for item in fragment.allowed_adducts}
            if not observed_adduct or observed_adduct not in allowed_adducts:
                continue
        if enforce_scope and (
            fragment.collision_energy_min is not None
            or fragment.collision_energy_max is not None
        ):
            value, unit = collision_energy_value_unit(query_collision_energy)
            expected_unit = str(fragment.collision_energy_unit or "").strip()
            if value is None or not unit or not expected_unit or unit.lower() != expected_unit.lower():
                continue
            if fragment.collision_energy_min is not None and value < fragment.collision_energy_min:
                continue
            if fragment.collision_energy_max is not None and value > fragment.collision_energy_max:
                continue
        fragments.append(fragment)
    return tuple(fragments)


def get_adduct_profile(adduct_name: str) -> Optional[AdductProfile]:
    """获取加合物特征"""
    return ADDUCT_PROFILES.get(adduct_name)


def match_fragment_to_class(
    mz: float, tolerance: float = 0.05,
    polarity: str = "",
    validated_only: bool = True,
) -> List[Dict[str, Any]]:
    """
    将观察到的碎片m/z匹配到最可能的化合物类别

    Args:
        mz: 观察到的碎片 m/z
        tolerance: 匹配容差 (Da)

    Returns:
        [{"class_name": ..., "fragment": ..., "diagnostic_value": ...}, ...]
    """
    matches = []
    for class_name, profile in COMPOUND_CLASS_PROFILES.items():
        # The generic glycoside profile deliberately reuses sugar-residue ions
        # for candidate validation.  Those ions alone cannot distinguish a
        # free sugar from a glycoside, so this profile must not vote in the
        # fragment-only classifier (which would create a circular/tied class).
        if class_name == "glycoside":
            continue
        for frag in get_fragments_by_class(
            class_name, polarity=polarity, validated_only=validated_only
        ):
            if abs(frag.mz - mz) <= tolerance:
                matches.append({
                    "class_name": class_name,
                    "class_name_cn": profile.class_name_cn,
                    "fragment": frag,
                    "diagnostic_value": frag.diagnostic_value,
                    "fragmentation_path": frag.fragmentation_path,
                })
    matches.sort(key=lambda x: x["diagnostic_value"], reverse=True)
    return matches


def classify_compound_from_fragments(
    fragment_mzs: List[float],
    tolerance: float = 0.05,
    polarity: str = "",
) -> Dict[str, Any]:
    """
    基于碎片特征对化合物进行分类 (v3.3.0 Enhanced)

    改进算法:
    1. 排他性加权 (Exclusivity Weighting):
       - 碎片若仅匹配单一类别 → 权重 ×1.5 (排他性碎片)
       - 碎片若匹配多个类别 → 按 diagnostic_value 占比分配
    2. 类大小归一化 (Class Size Normalization):
       - 将原始得分除以该类碎片总数 → 防止碎片多的类(如糖类)占优势
    3. 置信度阈值 (Confidence Threshold):
       - 最高分 < 0.15 → 返回 Unknown
       - 最高分/第二名 < 1.3 → 标记 low_confidence
    4. 默认容差从 0.02 Da 放宽到 0.05 Da

    Returns:
        {"class_name": "flavonoid", "confidence": 0.85, "all_scores": {"flavonoid": 0.85, ...}}
        或 {} (无法分类时)
    """
    # ── Step 1: 收集每个碎片匹配到的所有类别 ──
    fragment_matches: List[List[Dict[str, Any]]] = []
    for mz in fragment_mzs:
        matches = match_fragment_to_class(
            mz, tolerance, polarity=polarity, validated_only=True
        )
        if matches:
            fragment_matches.append(matches)

    if not fragment_matches:
        return {}

    # ── Step 2: 排他性加权评分 ──
    raw_scores = defaultdict(float)
    for matches in fragment_matches:
        if len(matches) == 1:
            # 排他性碎片: 仅匹配一个类别 → 加权 ×1.5
            match = matches[0]
            raw_scores[match["class_name"]] += match["diagnostic_value"] * 1.5
        else:
            # 多类别碎片: 按 diagnostic_value 占比分配
            total_dv = sum(m["diagnostic_value"] for m in matches)
            for match in matches:
                weight = match["diagnostic_value"] / total_dv
                raw_scores[match["class_name"]] += match["diagnostic_value"] * weight

    if not raw_scores:
        return {}

    # ── Step 3: 类大小归一化 ──
    # 碎片多的类(如 sugar 有22个碎片) vs 碎片少的类(如 quinone 有4个碎片)
    # 归一化: score / sqrt(n_fragments_in_class) 防止碎片数量偏差
    normalized = {}
    for class_name, raw_score in raw_scores.items():
        profile = COMPOUND_CLASS_PROFILES.get(class_name)
        if profile:
            n_frags = max(len(get_fragments_by_class(
                class_name, polarity=polarity, validated_only=True
            )), 1)
            # 使用 sqrt 归一化: 减少但不完全消除碎片数量的影响
            normalized[class_name] = raw_score / (n_frags ** 0.5)
        else:
            normalized[class_name] = raw_score

    # 再次归一化到 [0, 1]
    total_norm = sum(normalized.values())
    if total_norm > 0:
        normalized = {k: round(v / total_norm, 4) for k, v in normalized.items()}
    else:
        return {}

    sorted_scores = dict(sorted(normalized.items(), key=lambda x: x[1], reverse=True))

    # ── Step 4: 置信度阈值检查 ──
    items = list(sorted_scores.items())
    best_class, best_score = items[0]
    second_score = items[1][1] if len(items) > 1 else 0.0

    # 最低置信度阈值
    MIN_CONFIDENCE = 0.15
    if best_score < MIN_CONFIDENCE:
        return {
            "class_name": "Unknown",
            "confidence": round(best_score, 4),
            "all_scores": sorted_scores,
            "note": "below_min_confidence_threshold",
        }

    # 置信度评级
    ratio = best_score / second_score if second_score > 0 else float("inf")
    if ratio < 1.3 and len(items) > 1:
        note = "low_confidence_close_scores"
    else:
        note = "confident"

    return {
        "class_name": best_class,
        "confidence": round(best_score, 4),
        "all_scores": sorted_scores,
        "confidence_ratio": round(ratio, 2),
        "note": note,
    }


def get_diagnostic_ion_pairs(
    class_name: str,
    polarity: str = "",
) -> Tuple[Tuple[float, float], ...]:
    """获取经极性与质量一致性校验的异构体诊断离子对。"""
    profile = get_profile(class_name)
    if not profile:
        return ()
    if not normalize_polarity(polarity):
        return profile.diagnostic_ion_pairs
    applicable_mzs = {
        round(fragment.mz, 4)
        for fragment in get_fragments_by_class(
            class_name, polarity=polarity, validated_only=True
        )
    }
    return tuple(
        pair for pair in profile.diagnostic_ion_pairs
        if round(pair[0], 4) in applicable_mzs and round(pair[1], 4) in applicable_mzs
    )


def get_fragmentation_pathway(class_name: str) -> str:
    """获取完整裂解路径描述"""
    profile = get_profile(class_name)
    return profile.fragmentation_pathway if profile else ""


def get_isomer_rules(class_name: str) -> Tuple[str, ...]:
    """获取异构体区分规则"""
    profile = get_profile(class_name)
    return profile.isomer_differentiation_rules if profile else ()


def list_all_classes() -> List[str]:
    """列出所有已编码的化合物类别"""
    return sorted(COMPOUND_CLASS_PROFILES.keys())


def get_data_quality_standards() -> Dict[str, Any]:
    """获取数据质量控制标准"""
    return DATA_QUALITY_STANDARDS


def get_tcm_processing_features() -> Dict[str, Any]:
    """获取中药炮制/配伍质谱特征"""
    return TCM_PROCESSING_FEATURES


def get_structure_elucidation_rules() -> Dict[str, Any]:
    """获取结构解析推理规则"""
    return STRUCTURE_ELUCIDATION_RULES


def search_fragment_by_formula(formula: str) -> List[Dict[str, Any]]:
    """按分子式搜索碎片"""
    results = []
    for class_name, profile in COMPOUND_CLASS_PROFILES.items():
        for frag in profile.characteristic_fragments:
            if formula.upper() in frag.formula.upper():
                results.append({
                    "class_name": class_name,
                    "class_name_cn": profile.class_name_cn,
                    "fragment": frag,
                })
    return results


def export_library(filepath: str) -> None:
    """导出完整特征库为JSON"""
    data = {
        "compound_classes": {},
        "adduct_profiles": {},
        "data_quality_standards": DATA_QUALITY_STANDARDS,
        "tcm_processing": TCM_PROCESSING_FEATURES,
        "structure_elucidation_rules": STRUCTURE_ELUCIDATION_RULES,
    }
    for name, profile in COMPOUND_CLASS_PROFILES.items():
        data["compound_classes"][name] = {
            "class_name_cn": profile.class_name_cn,
            "parent_class": profile.parent_class,
            "fragmentation_pathway": profile.fragmentation_pathway,
            "sequential_losses": list(profile.sequential_losses),
            "isomer_differentiation_rules": list(profile.isomer_differentiation_rules),
            "reference_compounds": list(profile.reference_compounds),
            "literature_refs": list(profile.literature_refs),
            "confidence_score": profile.confidence_score,
            "characteristic_fragments": [
                {
                    "mz": f.mz, "formula": f.formula, "name": f.name,
                    "name_en": f.name_en, "intensity_rank": f.intensity_rank,
                    "structural_origin": f.structural_origin,
                    "diagnostic_value": f.diagnostic_value,
                    "fragmentation_path": f.fragmentation_path,
                }
                for f in profile.characteristic_fragments
            ],
            "diagnostic_ion_pairs": [
                [pair[0], pair[1]] for pair in profile.diagnostic_ion_pairs
            ],
        }

    for name, adduct in ADDUCT_PROFILES.items():
        data["adduct_profiles"][name] = {
            "mass_shift": adduct.mass_shift,
            "polarity": adduct.polarity,
            "prevalence": adduct.prevalence,
            "typical_compounds": adduct.typical_compounds,
            "fragmentation_behavior": adduct.fragmentation_behavior,
        }

    with open(filepath, "w", encoding="utf-8") as f:
        json.dump(data, f, ensure_ascii=False, indent=2)


# ─────────────────────────────────────────────────────────────────
# 统计摘要
# ─────────────────────────────────────────────────────────────────


def get_library_summary() -> Dict[str, Any]:
    """获取特征库统计摘要"""
    total_fragments = sum(
        len(p.characteristic_fragments)
        for p in COMPOUND_CLASS_PROFILES.values()
    )
    total_ion_pairs = sum(
        len(p.diagnostic_ion_pairs)
        for p in COMPOUND_CLASS_PROFILES.values()
    )
    return {
        "total_compound_classes": len(COMPOUND_CLASS_PROFILES),
        "total_characteristic_fragments": total_fragments,
        "total_diagnostic_ion_pairs": total_ion_pairs,
        "total_adduct_profiles": len(ADDUCT_PROFILES),
        "total_quality_standards": len(DATA_QUALITY_STANDARDS),
        "total_tcm_features": len(TCM_PROCESSING_FEATURES),
        "total_elucidation_rules": len(STRUCTURE_ELUCIDATION_RULES),
        "classes_covered": list(COMPOUND_CLASS_PROFILES.keys()),
        "average_confidence": round(
            sum(p.confidence_score for p in COMPOUND_CLASS_PROFILES.values())
            / len(COMPOUND_CLASS_PROFILES), 3
        ),
        "knowledge_sources": [
            "质谱分析技术原理与应用 (台湾质谱学会, 2019)",
            "中药质谱分析 (刘淑莹等, 2012)",
            "Spectraverse (Gupta et al., Anal. Chem. 2026)",
        ],
    }
