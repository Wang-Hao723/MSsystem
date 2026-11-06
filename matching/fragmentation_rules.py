"""
质谱裂解规则专家知识库 — Expert Fragmentation Rule Knowledge Base
====================================================================

编码了来自安捷伦(Agilent)、沃特世(Waters)、赛默飞(Thermo Fisher)等主流质谱厂商
以及质谱学经典文献中的裂解规则。覆盖 ESI-MS/MS 中常见的碎裂机理。

参考文献:
  - McLafferty, F.W. "Interpretation of Mass Spectra", 4th Ed.
  - Steen, H. & Mann, M. "The ABC's (and XYZ's) of peptide sequencing" Nat Rev Mol Cell Biol, 2004
  - Demarque, D.P. et al. "Fragmentation reactions using electrospray ionization mass spectrometry" 2016
  - Agilent MassHunter Qualitative Analysis - Fragmentation Library
  - Waters MassLynx - Fragmentation Rules Database
  - NIST MS/MS Library fragmentation patterns

Author: Expert System Module
Version: 2.0.0
"""

from __future__ import annotations

import json
import logging
import re
from collections import defaultdict
from dataclasses import dataclass, field
from itertools import combinations_with_replacement
from pathlib import Path
from typing import Any, Dict, FrozenSet, List, Optional, Set, Tuple

from .rule_quality import (
    collision_energy_value_unit,
    formula_mz_consistent,
    normalize_adduct,
    normalize_polarity,
    smiles_matches_all,
    smiles_matches_any,
)

logger = logging.getLogger(__name__)

# ─────────────────────────────────────────────────────────────────
# 基础数据结构
# ─────────────────────────────────────────────────────────────────


@dataclass(frozen=True)
class NeutralLoss:
    """中性丢失描述"""
    formula: str          # 分子式，如 "H2O", "CO2"
    mass: float           # 精确质量(Da)
    name: str             # 中文名
    name_en: str          # 英文名
    functional_groups: FrozenSet[str]  # 关联官能团
    mechanism: str        # 机理类型: rearrangement/cleavage/elimination
    priority: int = 0     # 优先级(越高越常见)
    polarity: str = "both"
    allowed_adducts: Tuple[str, ...] = ()
    required_smarts: Tuple[str, ...] = ()
    required_smarts_all: Tuple[str, ...] = ()
    collision_energy_min: Optional[float] = None
    collision_energy_max: Optional[float] = None
    collision_energy_unit: str = ""


@dataclass(frozen=True)
class FragmentationRule:
    """裂解规则"""
    rule_id: str                    # 唯一标识
    name: str                       # 规则名称(中文)
    name_en: str                    # 英文名称
    mechanism_type: str             # α-cleavage / β-cleavage / McLafferty / rDA / etc.
    description: str                # 详细描述
    structural_requirement: str     # 结构要求(SMARTS 或文字描述)
    smarts_pattern: str             # SMARTS 反应模式 (反应物>>产物)
    functional_groups: FrozenSet[str]  # 涉及的官能团
    typical_neutral_losses: Tuple[NeutralLoss, ...]  # 典型中性丢失
    charge_preference: str          # "positive" / "negative" / "both"
    energy_level: str               # "low" / "medium" / "high" (碰撞能量)
    references: Tuple[str, ...]     # 文献来源
    confidence_score: float = 0.8   # 该规则的可靠性评分(0-1)
    ionization_methods: Tuple[str, ...] = ("ESI",)
    allowed_adducts: Tuple[str, ...] = ()
    required_smarts: Tuple[str, ...] = ()
    required_smarts_all: Tuple[str, ...] = ()
    collision_energy_min: Optional[float] = None
    collision_energy_max: Optional[float] = None
    collision_energy_unit: str = ""
    evidence_level: str = "legacy_curated"


@dataclass(frozen=True)
class DiagnosticIon:
    """诊断离子(特征碎片)"""
    mz: float                       # m/z 值
    formula: str                    # 分子式
    name: str                       # 名称
    structural_origin: str          # 结构来源描述
    functional_group: str           # 关联官能团
    rule_id: str                    # 关联裂解规则ID
    polarity: str = "positive"      # positive/negative
    intensity_characteristic: str = "medium"  # low/medium/high
    compatible_class_terms: Tuple[str, ...] = ()  # 原始分类名约束(任一命中)
    required_smarts: Tuple[str, ...] = ()     # 候选结构约束(任一命中)
    required_smarts_all: Tuple[str, ...] = () # 候选结构约束(全部命中)
    allowed_adducts: Tuple[str, ...] = ()     # 明确适用的前体离子
    collision_energy_min: Optional[float] = None
    collision_energy_max: Optional[float] = None
    collision_energy_unit: str = ""
    references: Tuple[str, ...] = ()          # 可追溯文献/数据库标识
    evidence_level: str = "curated"           # curated/tentative/legacy


@dataclass
class BondDissociationEnergy:
    """化学键解离能"""
    bond_type: str              # 键类型，如 "C-C", "C-O", "C-N"
    environment: str            # 化学环境，如 "aliphatic", "allylic", "benzylic"
    energy_kj_per_mol: float    # 解离能 (kJ/mol)
    energy_kcal_per_mol: float  # 解离能 (kcal/mol)
    fragmentation_likelihood: float  # 裂解倾向分数(0-1, 越高越容易断)


# ─────────────────────────────────────────────────────────────────
# 常用中性丢失数据库
# ─────────────────────────────────────────────────────────────────

COMMON_NEUTRAL_LOSSES: Dict[str, NeutralLoss] = {
    # ── 小分子丢失(Low-mass loss) ──
    "H2O": NeutralLoss(
        formula="H2O", mass=18.0106, name="水", name_en="Water",
        functional_groups=frozenset({"alcohol", "carboxylic_acid", "aldehyde"}),
        mechanism="elimination", priority=10,
    ),
    "NH3": NeutralLoss(
        formula="NH3", mass=17.0265, name="氨", name_en="Ammonia",
        functional_groups=frozenset({"primary_amine", "amide"}),
        mechanism="elimination", priority=9,
        # NH3 loss is not a generic 17 Da label.  Require an explicit
        # primary amine/primary amide so N-oxides, alcohols and arbitrary
        # unknown candidates cannot acquire this chemically incompatible
        # explanation merely because a peak falls inside a loose window.
        required_smarts=(
            "[NX3;H2;+0]",
            "[CX3](=O)[NX3;H2;+0]",
        ),
    ),
    "CO": NeutralLoss(
        formula="CO", mass=27.9949, name="一氧化碳", name_en="Carbon monoxide",
        functional_groups=frozenset({"aldehyde", "ketone", "phenol", "quinone"}),
        mechanism="cleavage", priority=9,
    ),
    "CO2": NeutralLoss(
        formula="CO2", mass=43.9898, name="二氧化碳", name_en="Carbon dioxide",
        functional_groups=frozenset({"carboxylic_acid", "carbamate", "carbonate"}),
        mechanism="elimination", priority=10,
    ),
    "HCOOH": NeutralLoss(
        formula="CH2O2", mass=46.0055, name="甲酸", name_en="Formic acid",
        functional_groups=frozenset({"formate_ester", "carboxylic_acid_ortho"}),
        mechanism="rearrangement", priority=7,
    ),
    "CH3OH": NeutralLoss(
        formula="CH4O", mass=32.0262, name="甲醇", name_en="Methanol",
        functional_groups=frozenset({"methyl_ester", "methyl_ether"}),
        mechanism="elimination", priority=8,
    ),
    "CH3COOH": NeutralLoss(
        formula="C2H4O2", mass=60.0211, name="乙酸", name_en="Acetic acid",
        functional_groups=frozenset({"acetate_ester"}),
        mechanism="elimination", priority=7,
    ),
    "C6H5COOH": NeutralLoss(
        # 苯甲酸 (Benzoic acid) — 苯甲酰酯去酯基特征丢失。
        # 乌头碱型 C19-二萜生物碱(附子/川乌/草乌, 治类风关寒湿痹证核心药)的
        # 苯甲酰酯取代基丢失, 区分毒性构型: 双酯型DDA→脱苯甲酰基→乌头原碱。
        # 也见于苯甲酸/苯甲酰酯类(benzenoid/phenolic_acid)。
        # functional_groups 含 alkaloid/terpenoid 使 strict_fg 阶段对生物碱/
        # 萜类匹配(与 HCN/CH3CN 的"类名入 FG"既定模式一致, 见上 HCN 条目)。
        formula="C7H6O2", mass=122.0368, name="苯甲酸", name_en="Benzoic acid",
        functional_groups=frozenset({"benzoic_acid", "benzoyl", "benzoate_ester",
                                     "alkaloid", "terpenoid"}),
        mechanism="elimination", priority=7,
    ),
    "H2S": NeutralLoss(
        formula="H2S", mass=33.9877, name="硫化氢", name_en="Hydrogen sulfide",
        functional_groups=frozenset({"thiol", "thioether"}),
        mechanism="elimination", priority=6,
    ),
    "SO3": NeutralLoss(
        formula="SO3", mass=79.9568, name="三氧化硫", name_en="Sulfur trioxide",
        functional_groups=frozenset({"sulfonic_acid", "sulfate"}),
        mechanism="elimination", priority=7,
    ),
    "H3PO4": NeutralLoss(
        formula="H3PO4", mass=97.9769, name="磷酸", name_en="Phosphoric acid",
        functional_groups=frozenset({"phosphate", "phosphonate"}),
        mechanism="elimination", priority=7,
    ),

    # ── 常见碎片离子 ──
    "CH3": NeutralLoss(
        formula="CH3", mass=15.0235, name="甲基自由基", name_en="Methyl radical",
        functional_groups=frozenset({"methyl", "gem_dimethyl"}),
        mechanism="cleavage", priority=6,
    ),
    "C2H4": NeutralLoss(
        formula="C2H4", mass=28.0313, name="乙烯", name_en="Ethylene",
        functional_groups=frozenset({"ethyl_ester", "ethyl_ether"}),
        mechanism="McLafferty", priority=8,
    ),
    "HCN": NeutralLoss(
        formula="HCN", mass=27.0109, name="氢氰酸", name_en="Hydrogen cyanide",
        functional_groups=frozenset({"nitrile", "aromatic_nitrogen", "alkaloid", "purine"}),
        mechanism="cleavage", priority=7,
    ),

    # ── 含氮化合物特征丢失(N-containing specific) ──
    "CH3NCO": NeutralLoss(
        formula="CH3NCO", mass=57.0215, name="异氰酸甲酯", name_en="Methyl isocyanate",
        functional_groups=frozenset({"methylxanthine", "alkaloid", "purine", "caffeine"}),
        mechanism="retro_Diels_Alder", priority=9,
    ),
    "HNCO": NeutralLoss(
        formula="HNCO", mass=43.0058, name="异氰酸", name_en="Isocyanic acid",
        functional_groups=frozenset({"purine", "pyrimidine", "alkaloid", "urea", "carbamate"}),
        mechanism="retro_Diels_Alder", priority=8,
    ),
    "CH3CN": NeutralLoss(
        formula="C2H3N", mass=41.0265, name="乙腈", name_en="Acetonitrile",
        functional_groups=frozenset({"N_heterocycle", "alkaloid", "nitrile"}),
        mechanism="cleavage", priority=6,
    ),
    "CH2O": NeutralLoss(
        formula="CH2O", mass=30.0106, name="甲醛", name_en="Formaldehyde",
        functional_groups=frozenset({"methoxy", "methylenedioxy", "N_methyl"}),
        mechanism="elimination", priority=7,
    ),
    "C2H2O": NeutralLoss(
        formula="C2H2O", mass=42.0106, name="乙烯酮", name_en="Ketene",
        functional_groups=frozenset({"acetyl", "acetate", "ketone"}),
        mechanism="elimination", priority=7,
    ),
    "C3H6": NeutralLoss(
        formula="C3H6", mass=42.0470, name="丙烯", name_en="Propylene",
        functional_groups=frozenset({"isoprenyl", "terpenoid", "prenyl"}),
        mechanism="McLafferty", priority=7,
    ),

    # ── 糖类特征丢失 ──
    "C6H10O5": NeutralLoss(
        formula="C6H10O5", mass=162.0528, name="己糖残基", name_en="Hexose residue",
        functional_groups=frozenset({"glycoside", "sugar"}),
        mechanism="cleavage", priority=10,
    ),
    "C6H10O4": NeutralLoss(
        formula="C6H10O4", mass=146.0579, name="脱氧己糖", name_en="Deoxyhexose",
        functional_groups=frozenset({"rhamnoside", "deoxysugar"}),
        mechanism="cleavage", priority=9,
    ),
    "C5H8O4": NeutralLoss(
        formula="C5H8O4", mass=132.0423, name="戊糖残基", name_en="Pentose residue",
        functional_groups=frozenset({"pentoside"}),
        mechanism="cleavage", priority=9,
    ),
    "C6H8O6": NeutralLoss(
        formula="C6H8O6", mass=176.0321, name="葡萄糖醛酸", name_en="Glucuronic acid",
        functional_groups=frozenset({"glucuronide"}),
        mechanism="cleavage", priority=9,
    ),
    "C5H8O3": NeutralLoss(
        formula="C5H8O3", mass=116.0473, name="脱氧核糖", name_en="Deoxyribose",
        functional_groups=frozenset({"nucleoside", "deoxyribose"}),
        mechanism="N_glycosidic_cleavage", priority=10,
        polarity="positive", allowed_adducts=("[M+H]+",),
        required_smarts=(
            "[n;R]-[C;R]1[O;R][C;R]([CH2][O;H1,-1])[C;R]([O;H1,-1])[CH2;R]1",
            "[n;R]-[C;R]1[O;R][C;R]([CH2][O;H1,-1])[CH2;R][C;R]1[O;H1,-1]",
            "[n;R]-[C;R]1[O;R][C;R]([CH3])[C;R]([O;H1,-1])[C;R]1[O;H1,-1]",
        ),
    ),
    "H2O_CO": NeutralLoss(
        formula="CH2O2", mass=46.0055, name="水和一氧化碳",
        name_en="Water plus carbon monoxide",
        functional_groups=frozenset({"amino_acid", "alpha_amino_acid"}),
        mechanism="sequential_cleavage", priority=9,
        polarity="positive", allowed_adducts=("[M+H]+",),
        required_smarts=("[N;X3,X4+]-[C;X4]-[C](=O)[O;H1,-1]",),
    ),
    "C2H8NO4P": NeutralLoss(
        formula="C2H8NO4P", mass=141.0191, name="磷酸乙醇胺头基",
        name_en="Phosphoethanolamine head group",
        functional_groups=frozenset({"lipid", "phosphatidylethanolamine", "phosphate"}),
        mechanism="headgroup_cleavage", priority=10,
        polarity="positive", allowed_adducts=("[M+H]+",),
        required_smarts=("[P](=O)(O)OCC[NX3;H1,H2]", "[P](=O)([O-])OCC[NX3;H1,H2]"),
    ),
    "C4H8O4": NeutralLoss(
        formula="C4H8O4", mass=120.0423, name="糖环交叉裂解120",
        name_en="Sugar cross-ring loss 120",
        functional_groups=frozenset({"glycoside", "C_glycoside", "sugar"}),
        mechanism="cross_ring_cleavage", priority=9,
        polarity="negative", allowed_adducts=("[M-H]-",),
        required_smarts=("[c]-[C;R]1[O;R][C;R][C;R][C;R][C;R]1",),
    ),
    "C3H6O3": NeutralLoss(
        formula="C3H6O3", mass=90.0317, name="糖环交叉裂解90",
        name_en="Sugar cross-ring loss 90",
        functional_groups=frozenset({"glycoside", "C_glycoside", "sugar"}),
        mechanism="cross_ring_cleavage", priority=9,
        polarity="negative", allowed_adducts=("[M-H]-",),
        required_smarts=("[c]-[C;R]1[O;R][C;R][C;R][C;R][C;R]1",),
    ),
    "C2H4O2_CGLY": NeutralLoss(
        formula="C2H4O2", mass=60.0211, name="糖环交叉裂解60",
        name_en="Sugar cross-ring loss 60",
        functional_groups=frozenset({"glycoside", "C_glycoside", "sugar"}),
        mechanism="cross_ring_cleavage", priority=8,
        polarity="negative", allowed_adducts=("[M-H]-",),
        required_smarts=("[c]-[C;R]1[O;R][C;R][C;R][C;R][C;R]1",),
    ),

    # ── 黄酮类特征丢失 ──
    "C6H5O": NeutralLoss(
        formula="C6H5O", mass=93.0340, name="苯氧基", name_en="Phenoxy radical",
        functional_groups=frozenset({"flavonoid", "phenol"}),
        mechanism="cleavage", priority=7,
    ),
    "C8H8O": NeutralLoss(
        formula="C8H8O", mass=120.0575, name="B环C2H2O片段", name_en="B-ring fragment",
        functional_groups=frozenset({"flavonoid"}),
        mechanism="rDA", priority=8,
    ),
    "C7H5O2": NeutralLoss(
        formula="C7H5O2", mass=121.0290, name="B环裂解", name_en="B-ring cleavage",
        functional_groups=frozenset({"flavonoid", "isoflavonoid"}),
        mechanism="rDA", priority=8,
    ),

    # ── 生物碱类 ──
    "CH3N": NeutralLoss(
        formula="CH5N", mass=31.0422, name="甲胺", name_en="Methylamine",
        functional_groups=frozenset({"tertiary_amine", "alkaloid"}),
        mechanism="cleavage", priority=6,
    ),

    # ── 萜类/甾体 ──
    "C4H8": NeutralLoss(
        formula="C4H8", mass=56.0626, name="丁烯/异丁烯", name_en="Butene/isobutylene",
        functional_groups=frozenset({"terpenoid", "steroid", "isoprenyl"}),
        mechanism="cleavage", priority=7,
    ),
}


# ─────────────────────────────────────────────────────────────────
# 分类标准化规则 (v3.4.0 重构)
# 数据库的 classification 字段有 ~560 种取值: 约 88% 是来源文件名 (MoNA-export-.../
# MassBank_NIST/MSMS-.../YDF 等垃圾值), 其余是真实化学类别名 (含子类骨架、复数、
# 中英文混用)。这里用 *有序* 关键词规则把任意标签归一到规范类别之一。
# 顺序即优先级 —— 子类→父类、高特异→低特异, 解决重叠冲突:
#   "Bisbenzylisoquinoline" 命中 isoquinoline→alkaloid
#   "Flavonoid-3-O-glycosides" 命中 flavonoid→flavonoid (先于 glycoside)
#   "Steroids/Glycoside" 命中 steroid→steroid (先于 glycoside)
#   "Phenylpropanoid" 命中 phenylpropanoid→lignan (不再被旧的 "phenyl"→benzenoid 误伤)
# 无法识别的标签 (垃圾值/用途类别) 返回 "Unknown" —— 保守优于错误归类。
# ─────────────────────────────────────────────────────────────────

# 中文关键词 (对原始小写字符串做子串匹配; 顺序即优先级)
_CN_CLASS_RULES: List[tuple] = [
    ("胺氧化物", "amine_oxide"), ("氧化胺", "amine_oxide"),
    ("三萜皂苷", "saponin"), ("皂苷", "saponin"),                 # 皂苷优先于萜/甾
    ("甾体", "steroid"),
    ("花色", "flavonoid"), ("异黄酮", "flavonoid"), ("二氢黄酮", "flavonoid"),
    ("紫檀烷", "flavonoid"), ("黄烷醇", "flavonoid"), ("黄酮", "flavonoid"),
    ("香豆素", "coumarin"),
    ("醌", "quinone"),
    ("山酮", "benzenoid"),
    ("二芳基庚烷", "diarylheptanoid"),
    ("木脂素", "lignan"),
    ("姜辣素", "phenolic_acid"), ("原花青素", "phenolic_acid"), ("酚酸", "phenolic_acid"),
    ("萜", "terpenoid"), ("挥发油", "terpenoid"), ("倍半萜", "terpenoid"),
    ("二萜", "terpenoid"), ("三萜", "terpenoid"), ("单萜", "terpenoid"),
    ("核苷酸", "nucleoside"), ("核苷", "nucleoside"),
    ("氨基酸", "amino_acid"),
    ("生物碱", "alkaloid"),
    ("单糖", "sugar"), ("二糖", "sugar"), ("三糖", "sugar"), ("四糖", "sugar"),
    ("多糖", "sugar"), ("碳水化合物", "sugar"), ("糖", "sugar"),
    ("糖苷", "glycoside"),
    ("脂肪酸", "fatty_acid"),
    ("磷脂", "lipid"), ("脂质", "lipid"),
    ("苯并", "benzenoid"), ("苯环", "benzenoid"),
    ("磺酸", "sulfonate"),
]

# 英文关键词 (对去分隔符后的小写字符串做子串匹配; 顺序即优先级)
_EN_CLASS_RULES: List[tuple] = [
    # Tertiary-amine N-oxides are a distinct zwitterionic class.  Keep this
    # ahead of the broad amine/alkaloid fallbacks.
    ("amine_oxide", "amine_oxide"), ("amineoxide", "amine_oxide"),
    ("aminoxide", "amine_oxide"), ("n_oxide", "amine_oxide"),
    # 皂苷 (三萜/甾体糖苷, 优先于 terpene/steroid/glycoside)
    ("saponin", "saponin"),
    # 甾体
    ("sterol", "steroid"), ("steroid", "steroid"), ("steroidal", "steroid"),
    ("withanolide", "steroid"), ("bufanolide", "steroid"), ("spirosolane", "steroid"),
    ("cerveratrum", "steroid"), ("phytosterol", "steroid"), ("cholestane", "steroid"),
    ("epiminocholestane", "steroid"), ("cholesterol", "steroid"),
    # 黄酮 (优先于 phenol/glycoside)
    ("flavonoid", "flavonoid"), ("flavone", "flavonoid"), ("flavonol", "flavonoid"),
    ("flavanone", "flavonoid"), ("flavanonol", "flavonoid"), ("flavanol", "flavonoid"),
    ("flavan", "flavonoid"), ("isoflav", "flavonoid"), ("anthocyan", "flavonoid"),
    ("chalcone", "flavonoid"), ("aurone", "flavonoid"), ("catechin", "flavonoid"),
    ("pterocarpan", "flavonoid"), ("biflavonoid", "flavonoid"), ("flavylium", "flavonoid"),
    # 香豆素 (优先于 phenol/lactone)
    ("coumarin", "coumarin"),
    # 醌
    ("quinone", "quinone"), ("anthraquinone", "quinone"), ("naphthoquinone", "quinone"),
    ("benzoquinone", "quinone"),
    # 山酮是二苯并吡喃酮，不是醌；在没有专属运行时档案时保守落到芳香
    # 苯并类规则，而不能让结构分类器的 xanthone 输出退回 Unknown。
    ("xanthone", "benzenoid"),
    # 二芳基庚烷是独立的 Ar-C7-Ar 结构类，先于木脂素/苯丙素。
    ("diarylheptanoid", "diarylheptanoid"),
    ("diphenylheptanoid", "diarylheptanoid"),
    # 木脂素 / 二苯乙烯 / 苯丙素 (优先于 phenolic_acid)
    ("lignan", "lignan"), ("lignols", "lignan"), ("neolignan", "lignan"),
    ("stilbene", "lignan"), ("stilbenoid", "lignan"), ("podophyllotoxin", "lignan"),
    ("dibenzylbutyr", "lignan"), ("dibenzylbutane", "lignan"),
    ("dibenzocycloheptene", "lignan"), ("arylnaphthalene", "lignan"),
    ("epoxylignan", "lignan"), ("phenylpropanoid", "lignan"), ("phenylpropan", "lignan"),
    ("phenylethanoid", "lignan"), ("phenylethanol", "lignan"),
    ("cinnamyl", "lignan"), ("cinnamaldehyde", "lignan"),
    # 酚酸 / 多酚
    ("phenolic", "phenolic_acid"), ("polyphenol", "phenolic_acid"), ("phenol", "phenolic_acid"),
    ("hydroxycinnamic", "phenolic_acid"), ("hydroxybenzoic", "phenolic_acid"),
    ("caffeic", "phenolic_acid"), ("ferulic", "phenolic_acid"), ("coumaric", "phenolic_acid"),
    ("sinapic", "phenolic_acid"), ("sinapoyl", "phenolic_acid"), ("gallic", "phenolic_acid"),
    ("cinnamic", "phenolic_acid"), ("quinic", "phenolic_acid"), ("chlorogenic", "phenolic_acid"),
    ("phenolamide", "phenolic_acid"), ("salicylic", "phenolic_acid"), ("tannin", "phenolic_acid"),
    ("ellagic", "phenolic_acid"), ("sinapinic", "phenolic_acid"), ("sinapic", "phenolic_acid"),
    # 萜类 (甾体之后)
    ("terpene", "terpenoid"), ("terpenoid", "terpenoid"), ("triterpene", "terpenoid"),
    ("diterpene", "terpenoid"), ("sesquiterpene", "terpenoid"), ("monoterpene", "terpenoid"),
    ("iridoid", "terpenoid"), ("carotenoid", "terpenoid"), ("carotene", "terpenoid"),
    ("crocin", "terpenoid"), ("ginkgolide", "terpenoid"), ("steviol", "terpenoid"),
    ("retinoid", "terpenoid"), ("prenyl", "terpenoid"), ("prenol", "terpenoid"),
    ("isoprenoid", "terpenoid"), ("megastigmane", "terpenoid"), ("atisane", "terpenoid"),
    # 氨基酸 / 肽 (具体名称置于生物碱前; 裸 "amino" 移到生物碱之后, 以免
    # "Aminoquinolines" 等含生物碱骨架名的标签被误判为氨基酸)
    ("amino_acid", "amino_acid"), ("aminoacid", "amino_acid"), ("amino_acids", "amino_acid"),
    ("amino acid", "amino_acid"),
    ("peptide", "amino_acid"), ("dipeptide", "amino_acid"), ("tripeptide", "amino_acid"),
    ("peptidomimetic", "amino_acid"), ("phenylalanine", "amino_acid"),
    ("leucine", "amino_acid"), ("isoleucine", "amino_acid"), ("valine", "amino_acid"),
    ("histidine", "amino_acid"), ("glutamic", "amino_acid"), ("tyrosine", "amino_acid"),
    ("tryptophan", "amino_acid"), ("lysine", "amino_acid"), ("kainoid", "amino_acid"),
    # 核苷/核苷酸必须先于 purine→alkaloid 与 ribose→sugar。
    ("nucleoside", "nucleoside"), ("nucleotide", "nucleoside"),
    ("adenosine", "nucleoside"), ("guanosine", "nucleoside"),
    ("uridine", "nucleoside"), ("cytidine", "nucleoside"),
    ("thymidine", "nucleoside"), ("inosine", "nucleoside"),
    # 生物碱 (含大量子类骨架名: ClassyFire/NPClassifier 层级 + 俗名)
    ("alkaloid", "alkaloid"), ("morphinan", "alkaloid"), ("morphine", "alkaloid"),
    ("isoquinoline", "alkaloid"), ("benzylisoquinoline", "alkaloid"),
    ("aporphine", "alkaloid"), ("protoberberine", "alkaloid"),
    ("protopine", "alkaloid"), ("phthalideisoquinoline", "alkaloid"),
    ("tetrahydroisoquinoline", "alkaloid"), ("dihydroisoquinoline", "alkaloid"),
    ("quinoline", "alkaloid"), ("quinolone", "alkaloid"), ("quinolizidine", "alkaloid"),
    ("pyridine", "alkaloid"), ("piperidine", "alkaloid"), ("pyrrolidine", "alkaloid"),
    ("pyrrolizidine", "alkaloid"), ("indolizidine", "alkaloid"), ("quinuclidine", "alkaloid"),
    ("indole", "alkaloid"), ("indolic", "alkaloid"), ("carboline", "alkaloid"),
    ("ergoline", "alkaloid"), ("yohimbine", "alkaloid"), ("sarpagine", "alkaloid"),
    ("corynanthean", "alkaloid"), ("ibogan", "alkaloid"), ("eburnan", "alkaloid"),
    ("aspidosperm", "alkaloid"), ("macroline", "alkaloid"), ("camptothecin", "alkaloid"),
    ("cinchona", "alkaloid"), ("harmala", "alkaloid"), ("gelsemium", "alkaloid"),
    ("aconit", "alkaloid"), ("tropane", "alkaloid"), ("ecgonine", "alkaloid"),
    ("cytisine", "alkaloid"), ("sparteine", "alkaloid"), ("matrine", "alkaloid"),
    ("emetine", "alkaloid"), ("nicotine", "alkaloid"), ("veratrum", "alkaloid"),
    ("quinazoline", "alkaloid"), ("acridone", "alkaloid"), ("colchicine", "alkaloid"),
    ("purine", "alkaloid"), ("xanthine", "alkaloid"), ("ephedra", "alkaloid"),
    # 裸 "amino" (兜底: 如标签仅为 "Amino") 置于生物碱骨架名之后, 不再抢占
    # "Aminoquinolines" 这类含生物碱关键词的标签
    ("amino", "amino_acid"),
    # 有机酸 (兜底: OrganicAcid/OrganicAcids)
    ("organicacid", "fatty_acid"),
    # 糖类
    ("sugar", "sugar"), ("carbohydrate", "sugar"), ("saccharide", "sugar"),
    ("monosaccharide", "sugar"), ("disaccharide", "sugar"), ("oligosaccharide", "sugar"),
    ("polysaccharide", "sugar"), ("hexose", "sugar"), ("pentose", "sugar"),
    ("glucose", "sugar"), ("fructose", "sugar"), ("galactose", "sugar"), ("mannose", "sugar"),
    ("arabinose", "sugar"), ("xylose", "sugar"), ("ribose", "sugar"), ("fucose", "sugar"),
    ("rhamnose", "sugar"), ("glucuronic", "sugar"), ("inositol", "sugar"),
    # 糖苷 (通用, 置于各类特异糖苷之后)
    ("glycoside", "glycoside"),
    # 脂肪酸 / 有机酸
    ("fatty_acid", "fatty_acid"), ("fattyacid", "fatty_acid"), ("fatty acid", "fatty_acid"),
    ("fatty", "fatty_acid"), ("fahfa", "fatty_acid"), ("prostaglandin", "fatty_acid"),
    ("leukotriene", "fatty_acid"), ("organic_acid", "fatty_acid"), ("organic acid", "fatty_acid"),
    # 脂质
    ("lipid", "lipid"), ("phospholipid", "lipid"), ("phosphatidyl", "lipid"),
    ("glycerophospholipid", "lipid"), ("sphingolipid", "lipid"), ("ceramide", "lipid"),
    ("cer_ns", "lipid"), ("cer_np", "lipid"), ("cer_ap", "lipid"), ("cer_as", "lipid"),
    ("triglyceride", "lipid"), ("saccharolipid", "lipid"),
    # 苯环类 (通用芳香, 置于最后; 不再含误伤性 "phenyl" 关键词)
    ("benzenoid", "benzenoid"), ("benzopyran", "benzenoid"), ("benzofuran", "benzenoid"),
    ("benzaldehyde", "benzenoid"), ("biphenyl", "benzenoid"), ("phthalide", "benzenoid"),
    ("phenanthrene", "benzenoid"), ("naphthalene", "benzenoid"), ("anthracene", "benzenoid"),
    ("tetralin", "benzenoid"), ("indane", "benzenoid"), ("benzene", "benzenoid"),
    # 磺酸
    ("sulfonate", "sulfonate"), ("sulfonic", "sulfonate"), ("sulfonyl", "sulfonate"),
    ("sulfate", "sulfonate"), ("thioether", "sulfonate"), ("sulfoxide", "sulfonate"),
]


def normalize_classification(raw: str) -> str:
    """将化合物分类名称标准化为统一键名 (v3.4.0 有序规则重构)

    处理以下变体:
    - 大小写、空格、连字符、下划线
    - 中文类别名 (酚酸类→phenolic_acid, 生物碱类→alkaloid ...)
    - 复数形式 (flavonoids→flavonoid)
    - 子类/骨架名 (Morphinans/Aporphines/Benzylisoquinoline→alkaloid ...)
    - 复合名称子串匹配 (Steroids/Glycoside→steroid, Flavonoid-3-O-glycosides→flavonoid)
    - 数据库来源文件名 (MoNA-export-.../MassBank_NIST/YDF/MSMS-... → Unknown)

    使用 *有序* 关键词规则 (_CN_CLASS_RULES / _EN_CLASS_RULES): 第一个命中的关键词
    决定类别, 顺序解决子类→父类与跨类冲突。无法识别的标签返回 "Unknown"
    (保守策略, 优于错误归类)。

    Returns:
        标准化分类键名 (如 "phenolic_acid"), 无匹配返回 "Unknown"
    """
    if not raw or not raw.strip():
        return "Unknown"

    s = raw.lower().strip()

    # 空值 / 占位符 / 用途类别 (非化学分类)
    if s in ("n/a", "unknown", "", "other", "none", "-", "na"):
        return "Unknown"

    # ── 中文关键词 (对原始字符串) ──
    for cn_kw, std in _CN_CLASS_RULES:
        if cn_kw in s:
            return std

    # 统一分隔符: 空格/连字符 → 下划线, 便于英文子串匹配
    s_norm = s.replace(" ", "_").replace("-", "_")

    # ── 英文关键词 (对标准化字符串) ──
    for en_kw, std in _EN_CLASS_RULES:
        kw = en_kw.replace(" ", "_")  # "fatty acid" → "fatty_acid"
        if kw in s_norm:
            return std

    return "Unknown"


def neutral_loss_applicable(
    neutral_loss: NeutralLoss,
    polarity: str = "",
    candidate_smiles: str = "",
    query_adduct: str = "",
    query_collision_energy: Any = None,
) -> tuple[bool, str]:
    """Apply structured scope guards before a neutral loss can explain a peak."""
    wanted_polarity = normalize_polarity(polarity)
    loss_polarity = normalize_polarity(neutral_loss.polarity)
    if wanted_polarity and loss_polarity and wanted_polarity != loss_polarity:
        return False, "polarity_mismatch"
    if neutral_loss.allowed_adducts:
        observed_adduct = normalize_adduct(query_adduct)
        allowed_adducts = {
            normalize_adduct(item) for item in neutral_loss.allowed_adducts
        }
        if not observed_adduct:
            return False, "adduct_unknown"
        if observed_adduct not in allowed_adducts:
            return False, "adduct_mismatch"
    if neutral_loss.required_smarts:
        if not candidate_smiles:
            return False, "structure_unknown"
        if not smiles_matches_any(candidate_smiles, neutral_loss.required_smarts):
            return False, "structure_mismatch"
    if neutral_loss.required_smarts_all:
        if not candidate_smiles:
            return False, "structure_unknown"
        if not smiles_matches_all(candidate_smiles, neutral_loss.required_smarts_all):
            return False, "structure_mismatch"
    if (
        neutral_loss.collision_energy_min is not None
        or neutral_loss.collision_energy_max is not None
    ):
        value, unit = collision_energy_value_unit(query_collision_energy)
        expected_unit = str(neutral_loss.collision_energy_unit or "").strip()
        if value is None or not unit or not expected_unit:
            return False, "collision_energy_not_comparable"
        if unit.lower() != expected_unit.lower():
            return False, "collision_energy_unit_mismatch"
        if neutral_loss.collision_energy_min is not None and value < neutral_loss.collision_energy_min:
            return False, "collision_energy_below_range"
        if neutral_loss.collision_energy_max is not None and value > neutral_loss.collision_energy_max:
            return False, "collision_energy_above_range"
    return True, "applicable"


def get_neutral_loss_by_mass(mass: float, tolerance: float = 0.5,
                             classification: str = "",
                             strict_fg: bool = False,
                             polarity: str = "",
                             candidate_smiles: str = "",
                             query_adduct: str = "",
                             query_collision_energy: Any = None) -> List[NeutralLoss]:
    """按质量查找可能的中性丢失

    Args:
        mass: 质量差异(Da)
        tolerance: 质量匹配容差(Da)
        classification: 化合物分类(用于过滤不适用的中性丢失)
        strict_fg: 是否启用官能团严格过滤

    When strict_fg=True, only return NL entries whose functional_groups
    are compatible with the compound's classification.
    """
    # 标准化分类并修复 truthiness Bug (空字符串 "" 不再 bypass 过滤)
    effective_cls = normalize_classification(classification)

    results = []
    for nl in COMMON_NEUTRAL_LOSSES.values():
        if abs(nl.mass - mass) <= tolerance:
            applicable, _ = neutral_loss_applicable(
                nl,
                polarity=polarity,
                candidate_smiles=candidate_smiles,
                query_adduct=query_adduct,
                query_collision_energy=query_collision_energy,
            )
            if not applicable:
                continue
            # ── 分类兼容性过滤 ──
            # 即使 effective_cls == "Unknown"，也使用保守过滤（黑名单机制）
            if strict_fg:
                if not _is_nl_compatible(nl, effective_cls):
                    continue
            results.append(nl)
    # Accurate-mass agreement is the primary assignment criterion; priority
    # only breaks ties between chemically applicable candidates.
    return sorted(results, key=lambda x: (abs(x.mass - mass), -x.priority, x.formula))


def get_composite_neutral_loss_by_mass(
    mass: float,
    tolerance: float = 0.05,
    classification: str = "",
    strict_fg: bool = True,
) -> List[Dict[str, Any]]:
    """Return conservative two-step cumulative neutral-loss assignments.

    A precursor-to-fragment difference may represent two sequential losses
    (for example hexose followed by water), rather than one missing rule.  Only
    pairs of high-priority, class-compatible curated losses are considered and
    a tight mass tolerance is used to avoid explaining arbitrary differences.
    """
    effective_cls = normalize_classification(classification)
    # Only losses with established sequential/cumulative use are allowed.
    # Class-specific rearrangements (HNCO, CH3NCO, rDA fragments, etc.) must
    # not be freely combined merely because their summed mass fits.
    cumulative_keys = {
        "H2O", "NH3", "CO", "CO2", "CH3OH",
        "C6H10O5", "C6H10O4", "C5H8O4", "C6H8O6",
    }
    pool = [
        nl for key, nl in COMMON_NEUTRAL_LOSSES.items()
        if key in cumulative_keys
        and (not strict_fg or _is_nl_compatible(nl, effective_cls))
    ]
    results = []
    for first, second in combinations_with_replacement(pool, 2):
        combined_mass = first.mass + second.mass
        delta = abs(combined_mass - mass)
        if delta > tolerance:
            continue
        results.append({
            "name": f"{first.name} + {second.name}",
            "name_en": f"{first.name_en} + {second.name_en}",
            "formula": f"{first.formula} + {second.formula}",
            "mass": round(combined_mass, 4),
            "delta": round(delta, 4),
            "mechanism": "sequential/cumulative",
            "priority": min(first.priority, second.priority) - 1,
            "components": [first.name, second.name],
        })
    results.sort(key=lambda item: (item["delta"], -item["priority"], item["name"]))
    return results


# ── 分类到官能团映射(用于过滤中性丢失匹配) ──
_CLASS_TO_FGS: Dict[str, Set[str]] = {
    "amine_oxide": {"amine_oxide", "N_oxide", "tertiary_amine", "amine",
                     "N_methyl", "long_chain", "alkyl"},
    "flavonoid": {"flavonoid", "phenol", "alcohol", "glycoside", "ketone",
                  "carboxylic_acid", "aldehyde", "ether", "methyl_ester",
                  "methoxy", "methylenedioxy", "formate_ester", "carbonate",
                  "flavone", "flavonol"},
    "alkaloid": {"alkaloid", "N_heterocycle", "nitrile", "amine", "primary_amine",
                 "amide", "tertiary_amine", "purine", "methylxanthine", "urea",
                 "carbamate", "pyrimidine", "ketone", "alcohol", "carboxylic_acid",
                 "aldehyde"},
    "terpenoid": {"terpenoid", "isoprenyl", "prenyl", "steroid", "alkene",
                  "alcohol", "ketone", "carboxylic_acid", "ether",
                  "cyclohexene"},
    "phenolic_acid": {"carboxylic_acid", "phenol", "alcohol", "aldehyde",
                      "formate_ester", "carbonate", "carbamate", "aromatic",
                      "benzoic_acid", "benzyl"},
    "glycoside": {"glycoside", "sugar", "rhamnoside", "deoxysugar", "pentoside",
                  "glucuronide", "alcohol", "ether"},
    "sugar": {"sugar", "glycoside", "deoxysugar", "rhamnoside", "alcohol",
              "carboxylic_acid", "aldehyde", "ketone", "ether"},
    "fatty_acid": {"fatty_acid", "lipid", "long_chain", "carboxylic_acid",
                   "alcohol", "alkene"},
    "coumarin": {"alcohol", "carboxylic_acid", "ketone", "aldehyde",
                 "lactone", "phenol", "methyl_ester", "methyl_ether"},
    "quinone": {"ketone", "phenol", "quinone"},
    "steroid": {"steroid", "terpenoid", "alcohol", "ketone", "alkene"},
    "amino_acid": {"amine", "primary_amine", "amide", "carboxylic_acid",
                   "alcohol", "thiol"},
    "benzenoid": {"aromatic", "benzyl", "alkylbenzene", "phenol", "ether",
                  "alcohol", "carboxylic_acid", "benzoic_acid", "aniline"},
    "sulfonate": {"sulfonic_acid", "sulfate", "thiol", "thioether"},
    "lipid": {"lipid", "fatty_acid", "long_chain", "carboxylic_acid",
              "alcohol", "alkene", "phosphate"},
    "lignan": {"phenol", "alcohol", "methoxy", "methylenedioxy",
               "ether", "aromatic", "benzyl"},
    "saponin": {"glycoside", "sugar", "deoxysugar", "rhamnoside", "pentoside",
                "glucuronide", "alcohol", "ether", "terpenoid", "steroid"},
    "nucleoside": {"nucleoside", "ribose", "deoxyribose", "purine", "pyrimidine",
                    "phosphate", "primary_amine", "amide"},
}


def _is_nl_compatible(nl: NeutralLoss, classification: str) -> bool:
    """检查中性丢失是否与化合物分类化学兼容

    核心规则: NL 的 functional_groups 必须与化合物分类对应的
    合理官能团有交集,否则该中性丢失在化学上不合理。

    例如:
    - 乙腈(CH3CN) 需要 N_heterocycle/alkaloid/nitrile → 酚酸类不应匹配
    - 磷酸(H3PO4) 需要 phosphate/phosphonate → 大多数天然产物不应匹配
    - 异氰酸甲酯(CH3NCO) 需要 methylxanthine/alkaloid → 黄酮类不应匹配
    """
    cls_lower = classification.lower().strip()

    # 映射分类名到标准键
    cls_key = None
    for key in _CLASS_TO_FGS:
        if key in cls_lower:
            cls_key = key
            break

    if cls_key is None:
        # 未知分类: 使用保守过滤,排除化学上高度特异性的NL
        # 这些NL只在特定的化合物类型中出现,不应泛化匹配
        _highly_specific_fgs = {
            "phosphate", "phosphonate", "sulfonic_acid", "sulfate",
            "nitrile", "methylxanthine", "caffeine", "purine", "pyrimidine",
            "thiol", "thioether", "nitro",
            "N_heterocycle",  # 仅含氮杂环才可能丢失乙腈
            # 新增 — 乙烯酮 (C2H2O) 的 acetyl/acetate 需要酯/乙酰基上下文
            "acetyl", "acetate",
            # 新增 — 高度类特异标签
            "alkaloid", "glycoside",
        }
        if nl.functional_groups & _highly_specific_fgs:
            return False
        # 即使是通用NL,也检查质量是否在合理范围
        # NL质量 > 100 Da 的通常是糖类或黄酮特定片段,未知分类不应匹配
        if nl.mass > 100 and nl.name not in {"水", "一氧化碳", "二氧化碳"}:
            return False
        return True

    allowed_fgs = _CLASS_TO_FGS[cls_key]

    # NL 的至少一个 functional_group 必须在允许列表中
    if nl.functional_groups & allowed_fgs:
        return True

    # 一些通用NL(水、CO、CO2等)对几乎所有类别都合理
    _universal_nl_names = {"水", "一氧化碳", "二氧化碳", "甲基自由基", "甲醛"}
    if nl.name in _universal_nl_names:
        return True

    return False


# ─────────────────────────────────────────────────────────────────
# 裂解规则数据库
# ─────────────────────────────────────────────────────────────────

FRAGMENTATION_RULES: Dict[str, FragmentationRule] = {
    # ═══════════════════════════════════════════════════════════
    # 1. α-裂解 (Alpha-Cleavage) — 最常见的自由基引发裂解
    # ═══════════════════════════════════════════════════════════
    "ALPHA_CLEAVAGE_AMINE": FragmentationRule(
        rule_id="ALPHA_CLEAVAGE_AMINE",
        name="胺类α-裂解",
        name_en="Alpha-cleavage of amines",
        mechanism_type="alpha_cleavage",
        description=(
            "含氮化合物中，N原子的孤对电子不稳定化相邻的C-C键，"
            "导致α位C-C键断裂，产生亚胺阳离子。这是胺类和生物碱类"
            "化合物最特征的裂解方式。单电子转移至N原子→α位断键。"
        ),
        structural_requirement="含有胺基(伯/仲/叔胺)或含N杂环",
        smarts_pattern="[N;!R:1]-[C:2]-[C:3]>>[N:1]=[C:2]",
        functional_groups=frozenset({"amine", "alkaloid", "N_heterocycle"}),
        typical_neutral_losses=(
            COMMON_NEUTRAL_LOSSES.get("CH3NCO", NeutralLoss("CH3NCO", 57.0, "", "", frozenset(), "", 0)),
            COMMON_NEUTRAL_LOSSES.get("HNCO", NeutralLoss("HNCO", 43.0, "", "", frozenset(), "", 0)),
            COMMON_NEUTRAL_LOSSES.get("CH3", NeutralLoss("CH3", 15.0, "", "", frozenset(), "", 0)),
        ),
        charge_preference="positive",
        energy_level="low",
        references=("McLafferty 4th Ed. Ch.4", "Agilent MSD Productivity ChemStation Manual"),
        confidence_score=0.95,
    ),

    "ALPHA_CLEAVAGE_ETHER": FragmentationRule(
        rule_id="ALPHA_CLEAVAGE_ETHER",
        name="醚类α-裂解",
        name_en="Alpha-cleavage of ethers",
        mechanism_type="alpha_cleavage",
        description=(
            "醚的O原子孤对电子引发α-C-C或α-C-O键断裂，"
            "产生氧鎓离子(oxonium ion)。C-O键为主要的断裂位点。"
        ),
        structural_requirement="含有醚键(C-O-C)或醇-OH",
        smarts_pattern="[O:1]-[C:2]-[C:3]>>[O+:1]=[C:2]",
        functional_groups=frozenset({"ether", "alcohol"}),
        typical_neutral_losses=(),
        charge_preference="positive",
        energy_level="low",
        references=("McLafferty 4th Ed. Ch.7", "Waters QTOF Fragmentation Guide"),
        confidence_score=0.90,
    ),

    "ALPHA_CLEAVAGE_CARBONYL": FragmentationRule(
        rule_id="ALPHA_CLEAVAGE_CARBONYL",
        name="羰基α-裂解",
        name_en="Alpha-cleavage of carbonyl compounds",
        mechanism_type="alpha_cleavage",
        description=(
            "羰基(C=O)化合物中，O的孤对电子引发α位C-C键断裂，"
            "生成酰基阳离子(RC≡O+)。酮和醛类化合物的特征裂解。"
            "若两侧均有α-H，优先断裂取代基较多一侧(McLafferty规则)。"
        ),
        structural_requirement="含有羰基(酮/醛/酯/酰胺)",
        smarts_pattern="[C:1](=[O:2])-[C:3]>>[C+:1]#[O:2]",
        functional_groups=frozenset({"ketone", "aldehyde", "ester", "amide", "carboxylic_acid"}),
        typical_neutral_losses=(),
        charge_preference="positive",
        energy_level="low",
        references=("McLafferty 4th Ed. Ch.5", "Agilent 5977 MSD Theory Manual"),
        confidence_score=0.95,
    ),

    # ═══════════════════════════════════════════════════════════
    # 2. McLafferty 重排 — 最经典的重排反应
    # ═══════════════════════════════════════════════════════════
    "MCLAFFERTY_REARRANGEMENT": FragmentationRule(
        rule_id="MCLAFFERTY_REARRANGEMENT",
        name="McLafferty重排",
        name_en="McLafferty rearrangement",
        mechanism_type="McLafferty",
        description=(
            "含有γ-H的不饱和官能团(C=O, C=N, C=C等)经过六元环过渡态，"
            "γ-H转移到不饱和原子上，伴随β-键断裂。产生中性烯烃碎片+"
            "自由基阳离子烯醇或类似物。"
            "特征：Δm/z = 28 (C2H4), 42 (C3H6), 56 (C4H8) 等。"
            "必要条件：γ位有H，且六元环过渡态空间上可行。"
        ),
        structural_requirement="含有不饱和基团(C=O/C=N/C=C)且γ位有H",
        smarts_pattern="[C:1](=[O:2])-[C:3]-[C:4]-[H:5]>>[C:1]([O:2][H:5])=[C:3]",
        functional_groups=frozenset(
            {"ketone", "aldehyde", "ester", "carboxylic_acid", "amide"}
        ),
        typical_neutral_losses=(
            COMMON_NEUTRAL_LOSSES.get("C2H4", NeutralLoss("C2H4", 28.0, "", "", frozenset(), "", 0)),
        ),
        charge_preference="positive",
        energy_level="medium",
        references=(
            "McLafferty, F.W. Anal. Chem. 1959, 31, 82",
            "Agilent MassHunter Qualitative Analysis Guide",
        ),
        confidence_score=0.98,
    ),

    # ═══════════════════════════════════════════════════════════
    # 3. Retro-Diels-Alder (rDA) — 环状化合物的特征裂解
    # ═══════════════════════════════════════════════════════════
    "RETRO_DIELS_ALDER": FragmentationRule(
        rule_id="RETRO_DIELS_ALDER",
        name="逆Diels-Alder反应",
        name_en="Retro-Diels-Alder (rDA) reaction",
        mechanism_type="retro_Diels_Alder",
        description=(
            "含环己烯结构的化合物(如黄酮、萜类)经历[4+2]环加成逆反应，"
            "C环打开产生二烯+亲二烯体两个碎片。"
            "黄酮类化合物特征：产生m/z 152(A环)+m/z 120(B环)的互补离子对。"
            "黄酮C环的1,3和0,4键断裂。"
        ),
        structural_requirement="含有环己烯或类环己烯结构(如黄酮C环、萜类骨架)",
        smarts_pattern="[C:1]1[C:2][C:3]=[C:4][C:5][C:6]1>>[C:1]=[C:2].[C:3]=[C:4].[C:5]=[C:6]",
        functional_groups=frozenset({"flavonoid", "terpenoid", "steroid", "cyclohexene"}),
        typical_neutral_losses=(
            COMMON_NEUTRAL_LOSSES.get("C8H8O", NeutralLoss("C8H8O", 120.0, "", "", frozenset(), "", 0)),
            COMMON_NEUTRAL_LOSSES.get("C7H5O2", NeutralLoss("C7H5O2", 121.0, "", "", frozenset(), "", 0)),
            COMMON_NEUTRAL_LOSSES.get("CH3NCO", NeutralLoss("CH3NCO", 57.0, "", "", frozenset(), "", 0)),
            COMMON_NEUTRAL_LOSSES.get("HNCO", NeutralLoss("HNCO", 43.0, "", "", frozenset(), "", 0)),
        ),
        charge_preference="both",
        energy_level="medium",
        references=(
            "Ma, Y.L. et al. J. Am. Soc. Mass Spectrom. 1997, 8, 1180",
            "Fabre, N. et al. J. Am. Soc. Mass Spectrom. 2001, 12, 707",
            "Agilent Natural Product MS Library Guide",
        ),
        confidence_score=0.95,
    ),

    # ═══════════════════════════════════════════════════════════
    # 4. 诱导裂解 (Inductive Cleavage / i-cleavage)
    # ═══════════════════════════════════════════════════════════
    "INDUCTIVE_CLEAVAGE": FragmentationRule(
        rule_id="INDUCTIVE_CLEAVAGE",
        name="诱导裂解",
        name_en="Inductive cleavage (i-cleavage)",
        mechanism_type="inductive_cleavage",
        description=(
            "电荷定位在杂原子上后，通过诱导效应拉电子，"
            "导致相邻σ键极化并断裂。常见于卤代物、醚类和胺类。"
            "产生稳定的碳正离子或杂原子正离子。"
            "与α-裂解不同，i-cleavage主要靠诱导效应而非自由基。"
        ),
        structural_requirement="含有电负性杂原子(O/N/S/卤素)",
        smarts_pattern="[O,N,S,F,Cl,Br,I:1]-[C:2]>>[*+:1].[C:2]",
        functional_groups=frozenset({"halide", "ether", "alcohol", "amine", "thioether"}),
        typical_neutral_losses=(),
        charge_preference="positive",
        energy_level="low",
        references=("McLafferty 4th Ed. Ch.6", "Waters Synapt G2-Si Theory Guide"),
        confidence_score=0.85,
    ),

    # ═══════════════════════════════════════════════════════════
    # 5. 苄基裂解 (Benzylic Cleavage)
    # ═══════════════════════════════════════════════════════════
    "BENZYLIC_CLEAVAGE": FragmentationRule(
        rule_id="BENZYLIC_CLEAVAGE",
        name="苄基位裂解",
        name_en="Benzylic cleavage",
        mechanism_type="benzylic_cleavage",
        description=(
            "苄基位的C-C键由于形成的苄基阳离子共振稳定性极高而优先断裂。"
            "产生特征性的m/z 91 (C7H7+, 䓬离子/tropylium)。"
            "取代苄基产生C7H6R+系列碎片。"
        ),
        structural_requirement="含有苄基结构(芳香环-CH2-X)",
        smarts_pattern="[c:1]1[c:2][c:3][c:4][c:5][c:6]1-[C:7]-[X:8]>>[c:1]1[c:2][c:3][c:4][c:5][c:6]1=[C+:7]",
        functional_groups=frozenset({"aromatic", "benzyl", "alkylbenzene"}),
        typical_neutral_losses=(),
        charge_preference="positive",
        energy_level="low",
        references=("McLafferty 4th Ed. Ch.3", "NIST MS/MS Interpretation Guide"),
        confidence_score=0.95,
    ),

    # ═══════════════════════════════════════════════════════════
    # 6. 烯丙基裂解 (Allylic Cleavage)
    # ═══════════════════════════════════════════════════════════
    "ALLYLIC_CLEAVAGE": FragmentationRule(
        rule_id="ALLYLIC_CLEAVAGE",
        name="烯丙基位裂解",
        name_en="Allylic cleavage",
        mechanism_type="allylic_cleavage",
        description=(
            "双键的烯丙位C-C键由于共振稳定阳离子(allyl cation)而易于断裂。"
            "常见于萜类和脂肪酸。产生特征的烯丙基阳离子。"
        ),
        structural_requirement="含有烯丙基结构(C=C-CH2-X)",
        smarts_pattern="[C:1]=[C:2]-[C:3]-[X:4]>>[C:1]=[C:2]-[C+:3]",
        functional_groups=frozenset({"alkene", "terpenoid", "fatty_acid"}),
        typical_neutral_losses=(),
        charge_preference="positive",
        energy_level="low",
        references=("McLafferty 4th Ed. Ch.3", "Agilent GC/MS ChemStation Guide"),
        confidence_score=0.90,
    ),

    # ═══════════════════════════════════════════════════════════
    # 7. 消除反应 (Elimination)
    # ═══════════════════════════════════════════════════════════
    "DEHYDRATION": FragmentationRule(
        rule_id="DEHYDRATION",
        name="脱水反应",
        name_en="Dehydration (loss of H2O)",
        mechanism_type="elimination",
        description=(
            "醇羟基或羧基通过1,3-或1,4-消除丢失H2O。"
            "醇类经过碳正离子中间体消除水，羧酸邻位有H时更易发生。"
            "特征中性丢失：18 Da (H2O)。"
            "注意区分：热脱水 vs 碰撞诱导脱水。"
        ),
        structural_requirement="含有羟基(-OH)且β位有H",
        smarts_pattern="[C:1]-[C:2](-[O:3][H:4])-[H:5]>>[C:1]=[C:2]",
        functional_groups=frozenset({"alcohol", "carboxylic_acid", "phenol"}),
        typical_neutral_losses=(COMMON_NEUTRAL_LOSSES["H2O"],),
        charge_preference="positive",
        energy_level="low",
        references=("McLafferty 4th Ed. Ch.7", "Waters MassLynx Fragmentation Guide"),
        confidence_score=0.95,
    ),

    "DECARBOXYLATION": FragmentationRule(
        rule_id="DECARBOXYLATION",
        name="脱羧反应",
        name_en="Decarboxylation (loss of CO2)",
        mechanism_type="elimination",
        description=(
            "羧酸或羧酸酯通过消除CO2产生碳负离子/自由基。"
            "芳香族邻位有-OH/NH2等供电子基团时更易发生。"
            "特征中性丢失：44 Da (CO2)。"
        ),
        structural_requirement="含有羧基(-COOH)或羧酸酯",
        smarts_pattern="[C:1](=[O:2])-[O:3][H:4]>>[C:1][H:4]",
        functional_groups=frozenset({"carboxylic_acid", "carboxylate", "carbamate"}),
        typical_neutral_losses=(COMMON_NEUTRAL_LOSSES["CO2"],),
        charge_preference="both",
        energy_level="low",
        references=("McLafferty 4th Ed. Ch.8", "Agilent LC/MS TOF Manual"),
        confidence_score=0.92,
    ),

    # ═══════════════════════════════════════════════════════════
    # 8. 糖苷键裂解 (Glycosidic Cleavage)
    # ═══════════════════════════════════════════════════════════
    "GLYCOSIDIC_CLEAVAGE": FragmentationRule(
        rule_id="GLYCOSIDIC_CLEAVAGE",
        name="糖苷键裂解",
        name_en="Glycosidic bond cleavage",
        mechanism_type="cleavage",
        description=(
            "糖苷键(缩醛键)在酸性/碰撞条件下优先断裂，"
            "产生糖基正离子(oxocarbenium ion)和苷元。"
            "黄酮苷中糖基连接O-苷键断裂，逐一丢失糖单元。"
            "特征：丢失162 Da(己糖)、146 Da(脱氧己糖)、132 Da(戊糖)等。"
        ),
        structural_requirement="含有糖苷键(O-glycoside或C-glycoside)",
        smarts_pattern="[O:1]-[C:2]1-[O:3]-[C:4]-[C:5]-[C:6]1>>[O+:1]=[C:2]",
        functional_groups=frozenset({"glycoside", "sugar", "flavonoid_glycoside"}),
        typical_neutral_losses=(
            COMMON_NEUTRAL_LOSSES["C6H10O5"],
            COMMON_NEUTRAL_LOSSES["C6H10O4"],
            COMMON_NEUTRAL_LOSSES["C5H8O4"],
            COMMON_NEUTRAL_LOSSES["C6H8O6"],
        ),
        charge_preference="both",
        energy_level="low",
        references=(
            "Cuyckens, F. & Claeys, M. J. Mass Spectrom. 2004, 39, 1",
            "Agilent Flavonoid MS/MS Library",
        ),
        confidence_score=0.98,
        allowed_adducts=("[M+H]+", "[M-H]-"),
        required_smarts=(
            "[!#1]-[O;X2]-[C;R]1[O;R][C;R][C;R][C;R][C;R]1",
            "[!#1]-[O;X2]-[C;R]1[O;R][C;R][C;R][C;R]1",
        ),
        evidence_level="curated_primary",
    ),

    "C_GLYCOSIDE_CROSS_RING": FragmentationRule(
        rule_id="C_GLYCOSIDE_CROSS_RING",
        name="C-糖苷糖环交叉裂解",
        name_en="C-glycoside sugar cross-ring cleavage",
        mechanism_type="cross_ring_cleavage",
        description=(
            "C-糖苷不易直接断开C-C糖苷键，负离子ESI-CID更常见糖环交叉裂解，"
            "形成60.0211、90.0317和120.0423 Da中性丢失；这些质量必须与"
            "芳香碳-C糖环结构同时出现，不能当作通用脱糖证据。"
        ),
        structural_requirement="芳香碳直接连接吡喃糖环(C-glycosidic C-C bond)",
        smarts_pattern="[c:1]-[C:2]1[O:3][C:4][C:5][C:6][C:7]1",
        functional_groups=frozenset({"glycoside", "C_glycoside", "flavonoid"}),
        typical_neutral_losses=(
            COMMON_NEUTRAL_LOSSES["C2H4O2_CGLY"],
            COMMON_NEUTRAL_LOSSES["C3H6O3"],
            COMMON_NEUTRAL_LOSSES["C4H8O4"],
        ),
        charge_preference="negative",
        energy_level="medium",
        references=("PMID:15389856", "PMID:35088488"),
        confidence_score=0.95,
        allowed_adducts=("[M-H]-",),
        required_smarts=("[c]-[C;R]1[O;R][C;R][C;R][C;R][C;R]1",),
        evidence_level="curated_primary",
    ),

    "ALPHA_AMINO_ACID_46_LOSS": FragmentationRule(
        rule_id="ALPHA_AMINO_ACID_46_LOSS",
        name="质子化α-氨基酸失水/一氧化碳",
        name_en="Protonated alpha-amino-acid loss of H2O plus CO",
        mechanism_type="sequential_cleavage",
        description=(
            "质子化α-氨基酸CID可连续失去H2O和CO(总质量46.0055 Da)，"
            "并形成与侧链对应的亚铵离子。该规则不外推到β/γ/ε-氨基酸。"
        ),
        structural_requirement="α-氨基羧酸骨架",
        smarts_pattern="[N:1]-[C:2]-[C:3](=[O:4])[O:5]",
        functional_groups=frozenset({"amino_acid", "alpha_amino_acid"}),
        typical_neutral_losses=(COMMON_NEUTRAL_LOSSES["H2O_CO"],),
        charge_preference="positive",
        energy_level="medium",
        references=("PMID:39563088",),
        confidence_score=0.96,
        allowed_adducts=("[M+H]+",),
        required_smarts=("[N;X3,X4+]-[C;X4]-[C](=O)[O;H1,-1]",),
        evidence_level="curated_primary",
    ),

    "RIBONUCLEOSIDE_GLYCOSIDIC_CLEAVAGE": FragmentationRule(
        rule_id="RIBONUCLEOSIDE_GLYCOSIDIC_CLEAVAGE",
        name="核糖核苷N-糖苷键裂解",
        name_en="Ribonucleoside N-glycosidic-bond cleavage",
        mechanism_type="N_glycosidic_cleavage",
        description="质子化核糖核苷CID失去132.0423 Da核糖并形成质子化碱基。",
        structural_requirement="含N-糖苷键和2'-羟基的核糖核苷",
        smarts_pattern="[n:1]-[C:2]1[O:3][C:4][C:5][C:6]1",
        functional_groups=frozenset({"nucleoside", "ribose", "purine", "pyrimidine"}),
        typical_neutral_losses=(COMMON_NEUTRAL_LOSSES["C5H8O4"],),
        charge_preference="positive",
        energy_level="low",
        references=("PMCID:PMC11587168", "PMCID:PMC7249026", "DOI:10.1021/acs.chemrestox.3c00041"),
        confidence_score=0.98,
        allowed_adducts=("[M+H]+",),
        required_smarts=(
            "[n;R]-[C;R]1[O;R][C;R]([CH2][O;H1,-1])[C;R]([O;H1,-1])[C;R]1[O;H1,-1]",
        ),
        evidence_level="curated_primary",
    ),

    "DEOXYRIBONUCLEOSIDE_GLYCOSIDIC_CLEAVAGE": FragmentationRule(
        rule_id="DEOXYRIBONUCLEOSIDE_GLYCOSIDIC_CLEAVAGE",
        name="单脱氧核苷N-糖苷键裂解",
        name_en="Monodeoxynucleoside N-glycosidic-bond cleavage",
        mechanism_type="N_glycosidic_cleavage",
        description=(
            "质子化2'、3'或5'单脱氧核苷CID失去116.0473 Da脱氧戊糖并形成质子化碱基。"
        ),
        structural_requirement="含N-糖苷键和单脱氧戊糖的核苷",
        smarts_pattern="[n:1]-[C:2]1[O:3][C:4][C:5][CH2:6]1",
        functional_groups=frozenset({"nucleoside", "deoxyribose", "purine", "pyrimidine"}),
        typical_neutral_losses=(COMMON_NEUTRAL_LOSSES["C5H8O3"],),
        charge_preference="positive",
        energy_level="low",
        references=("PMCID:PMC11587168", "PMCID:PMC7249026", "DOI:10.1021/acs.chemrestox.3c00041"),
        confidence_score=0.98,
        allowed_adducts=("[M+H]+",),
        required_smarts=(
            "[n;R]-[C;R]1[O;R][C;R]([CH2][O;H1,-1])[C;R]([O;H1,-1])[CH2;R]1",
            "[n;R]-[C;R]1[O;R][C;R]([CH2][O;H1,-1])[CH2;R][C;R]1[O;H1,-1]",
            "[n;R]-[C;R]1[O;R][C;R]([CH3])[C;R]([O;H1,-1])[C;R]1[O;H1,-1]",
        ),
        evidence_level="curated_primary",
    ),

    "PHOSPHATIDYLETHANOLAMINE_141_LOSS": FragmentationRule(
        rule_id="PHOSPHATIDYLETHANOLAMINE_141_LOSS",
        name="磷脂酰乙醇胺头基丢失",
        name_en="Phosphatidylethanolamine headgroup loss",
        mechanism_type="headgroup_cleavage",
        description=(
            "质子化磷脂酰乙醇胺在正离子ESI低能CID中可失去"
            "C2H8NO4P(141.0191 Da)头基；仅对PE/LPE结构与[M+H]+适用。"
        ),
        structural_requirement="磷酸乙醇胺头基",
        smarts_pattern="[P:1](=O)(O)OCC[N:2]",
        functional_groups=frozenset({"lipid", "phosphatidylethanolamine", "phosphate"}),
        typical_neutral_losses=(COMMON_NEUTRAL_LOSSES["C2H8NO4P"],),
        charge_preference="positive",
        energy_level="low",
        references=("PMID:19277613", "DOI:10.1007/s00216-009-2717-9"),
        confidence_score=0.97,
        allowed_adducts=("[M+H]+",),
        required_smarts=("P(=O)(O)OCC[NX3;H1,H2]", "P(=O)([O-])OCC[NX3;H1,H2]"),
        evidence_level="curated_primary",
    ),

    "SULFATE_SO3_LOSS": FragmentationRule(
        rule_id="SULFATE_SO3_LOSS",
        name="硫酸酯失去SO3",
        name_en="Sulfate ester loss of SO3",
        mechanism_type="sulfate_cleavage",
        description="去质子化硫酸酯在负离子ESI-CID中失去79.9568 Da SO3。",
        structural_requirement="有机硫酸酯或硫酸根结构",
        smarts_pattern="[O:1]-[S:2](=O)(=O)-[O:3]",
        functional_groups=frozenset({"sulfonate", "sulfate"}),
        typical_neutral_losses=(COMMON_NEUTRAL_LOSSES["SO3"],),
        charge_preference="negative",
        energy_level="low",
        references=("PMCID:PMC6639915",),
        confidence_score=0.95,
        allowed_adducts=("[M-H]-",),
        required_smarts=("[O]-S(=O)(=O)[O;H1,-1]",),
        evidence_level="curated_primary",
    ),

    # ═══════════════════════════════════════════════════════════
    # 9. 黄酮类裂解特征
    # ═══════════════════════════════════════════════════════════
    "FLAVONOID_RDA_C_RING": FragmentationRule(
        rule_id="FLAVONOID_RDA_C_RING",
        name="黄酮C环rDA裂解",
        name_en="Flavonoid C-ring Retro-Diels-Alder",
        mechanism_type="retro_Diels_Alder",
        description=(
            "黄酮/黄酮醇的C环通过rDA反应裂解，产生特征性的A环和B环碎片离子。"
            "裂解位置在C环的1/3键和0/4键。"
            "产生的^1,3A+和^1,3B+离子可用于确定A环和B环的取代模式。"
            "黄酮母核在m/z 152, 120处产生特征碎片。"
        ),
        structural_requirement="黄酮C环(2-苯基色原酮骨架)",
        smarts_pattern=(
            "[O:1]=[C:2]1[C:3]=[C:4]([c:5]2[c:6][c:7][c:8][c:9][c:10]2)"
            "[O:11][c:12]2[c:13][c:14][c:15][c:16][c:17]2-1"
        ),
        functional_groups=frozenset({"flavonoid", "flavone", "flavonol"}),
        typical_neutral_losses=(
            COMMON_NEUTRAL_LOSSES.get("C8H8O", NeutralLoss("C8H8O", 120.0, "", "", frozenset(), "", 0)),
            COMMON_NEUTRAL_LOSSES.get("C7H5O2", NeutralLoss("C7H5O2", 121.0, "", "", frozenset(), "", 0)),
        ),
        charge_preference="both",
        energy_level="medium",
        references=(
            "Ma et al. JASMS 1997",
            "Cuyckens & Claeys, JMS 2004",
            "Waters Flavonoid Identification Application Note",
        ),
        confidence_score=0.97,
    ),

    # ═══════════════════════════════════════════════════════════
    # 10. 电荷远程裂解 (Charge-Remote Fragmentation, CRF)
    # ═══════════════════════════════════════════════════════════
    "CHARGE_REMOTE_FRAGMENTATION": FragmentationRule(
        rule_id="CHARGE_REMOTE_FRAGMENTATION",
        name="电荷远程裂解",
        name_en="Charge-remote fragmentation (CRF)",
        mechanism_type="charge_remote",
        description=(
            "在高能CID条件下，电荷固定在特定位置（如羧酸根），"
            "远离电荷的C-C键发生1,4-消除，产生系列烷烃/烯烃碎片。"
            "脂肪酸的负离子模式特征：CH2=CH2间隔的系列碎片(相差14 Da)。"
            "用于脂肪酸双键位置和支链位置的确认。"
        ),
        structural_requirement="长链脂肪酸/脂类(负离子模式)",
        smarts_pattern="[C:1]-[C:2]-[C:3]-[C:4]>>[C:1]=[C:2].[C:3]=[C:4]",
        functional_groups=frozenset({"fatty_acid", "lipid", "long_chain"}),
        typical_neutral_losses=(),
        charge_preference="negative",
        energy_level="high",
        references=(
            "Gross, M.L. J. Am. Soc. Mass Spectrom. 1992, 3, 105",
            "Waters High Energy CID Application Guide",
        ),
        confidence_score=0.80,
    ),

    "PHOSPHOLIPID_HEADGROUP_CLEAVAGE": FragmentationRule(
        rule_id="PHOSPHOLIPID_HEADGROUP_CLEAVAGE",
        name="磷脂极性头基裂解",
        name_en="Phospholipid polar-headgroup cleavage",
        mechanism_type="charge_directed_cleavage",
        description=(
            "ESI-CID/HCD 中磷脂极性头基保留电荷，形成类别特异的磷酸胆碱、"
            "磷酸乙醇胺或磷酸肌醇产物离子；必须结合离子模式与候选结构使用。"
        ),
        structural_requirement="含相应磷脂极性头基，且离子模式与产物离子一致",
        smarts_pattern="[P:1](=O)(O)O-[C:2]>>[P:1](=O)(O)O.[C:2]",
        functional_groups=frozenset({"lipid", "phospholipid", "phosphate"}),
        typical_neutral_losses=(),
        charge_preference="both",
        energy_level="medium",
        references=(
            "PMID:32423109",
            "PMCID:PMC7287777",
            "PMCID:PMC5697860",
        ),
        confidence_score=0.97,
    ),

    "SPHINGOID_BASE_CLEAVAGE": FragmentationRule(
        rule_id="SPHINGOID_BASE_CLEAVAGE",
        name="鞘氨醇长链碱裂解",
        name_en="Sphingoid-base cleavage",
        mechanism_type="charge_directed_cleavage",
        description=(
            "正离子 ESI-MS/MS 中神经酰胺、糖鞘脂等可形成反映长链碱组成的"
            "双脱水碳正离子；d18:1 常见于 m/z 264.268。链长或不饱和度改变时"
            "该离子会相应位移，不能作为全部鞘脂的固定通用峰。"
        ),
        structural_requirement="鞘氨醇/神经酰胺长链碱骨架，正离子模式",
        smarts_pattern="[N:1]-[C:2]-[C:3]([O:4])-[C:5]=[C:6]>>[N:1]-[C:2]-[C:3]=[C:5]-[C:6]",
        functional_groups=frozenset({"lipid", "sphingolipid", "ceramide"}),
        typical_neutral_losses=(COMMON_NEUTRAL_LOSSES["H2O"],),
        charge_preference="positive",
        energy_level="medium",
        references=("PMID:11070861", "PMCID:PMC2724058"),
        confidence_score=0.95,
    ),

    # ═══════════════════════════════════════════════════════════
    # 11. ortho效应 (Ortho Effect)
    # ═══════════════════════════════════════════════════════════
    "ORTHO_EFFECT": FragmentationRule(
        rule_id="ORTHO_EFFECT",
        name="邻位效应",
        name_en="Ortho effect",
        mechanism_type="rearrangement",
        description=(
            "芳香环上邻位取代基之间的相互作用导致的特征裂解。"
            "邻位有-OH和-COOH时，通过分子内氢键消除H2O/CO2。"
            "邻位有-NH2和-COOR时，消除ROH。"
            "这是区分位置异构体的重要工具。"
        ),
        structural_requirement="芳香环邻位有两个能相互作用的官能团",
        smarts_pattern="[c:1](-[X:2])[c:3](-[Y:4])>>[c:1][c:3]",
        functional_groups=frozenset({"aromatic", "phenol", "benzoic_acid", "aniline"}),
        typical_neutral_losses=(
            COMMON_NEUTRAL_LOSSES["H2O"],
            COMMON_NEUTRAL_LOSSES["CO2"],
        ),
        charge_preference="both",
        energy_level="low",
        references=("McLafferty 4th Ed. Ch.9", "Agilent Ortho Effect Application Brief"),
        confidence_score=0.85,
    ),
}


def get_rules_by_mechanism(mechanism_type: str) -> List[FragmentationRule]:
    """按机理类型获取规则"""
    return [r for r in FRAGMENTATION_RULES.values()
            if r.mechanism_type == mechanism_type]


def get_rules_by_functional_group(fg: str) -> List[FragmentationRule]:
    """按官能团获取相关规则"""
    results = []
    for rule in FRAGMENTATION_RULES.values():
        if fg in rule.functional_groups:
            results.append(rule)
    return sorted(results, key=lambda x: x.confidence_score, reverse=True)


def get_all_mechanism_types() -> List[str]:
    """获取所有机理类型"""
    return sorted(set(r.mechanism_type for r in FRAGMENTATION_RULES.values()))


# ─────────────────────────────────────────────────────────────────
# 化学键解离能数据库
# ─────────────────────────────────────────────────────────────────

BOND_DISSOCIATION_ENERGIES: Dict[str, BondDissociationEnergy] = {
    # C-C 键
    "C-C_aliphatic": BondDissociationEnergy(
        "C-C", "aliphatic", 347.0, 83.0, 0.5
    ),
    "C-C_allylic": BondDissociationEnergy(
        "C-C", "allylic", 251.0, 60.0, 0.85
    ),
    "C-C_benzylic": BondDissociationEnergy(
        "C-C", "benzylic", 264.0, 63.0, 0.80
    ),
    "C-C_alpha_carbonyl": BondDissociationEnergy(
        "C-C", "alpha_to_carbonyl", 293.0, 70.0, 0.70
    ),
    "C-C_aromatic": BondDissociationEnergy(
        "C-C", "aromatic", 418.0, 100.0, 0.05
    ),

    # C-O 键
    "C-O_aliphatic_ether": BondDissociationEnergy(
        "C-O", "aliphatic_ether", 358.0, 85.6, 0.45
    ),
    "C-O_glycosidic": BondDissociationEnergy(
        "C-O", "glycosidic", 293.0, 70.0, 0.75
    ),
    "C-O_ester_alkyl": BondDissociationEnergy(
        "C-O", "ester_alkyl", 360.0, 86.0, 0.40
    ),
    "C-O_phenolic": BondDissociationEnergy(
        "C-O", "phenolic", 418.0, 100.0, 0.10
    ),

    # C-N 键
    "C-N_aliphatic_amine": BondDissociationEnergy(
        "C-N", "aliphatic_amine", 305.0, 73.0, 0.65
    ),
    "C-N_amide": BondDissociationEnergy(
        "C-N", "amide", 355.0, 85.0, 0.35
    ),
    "C-N_alpha_amine": BondDissociationEnergy(
        "C-N", "alpha_to_amine", 272.0, 65.0, 0.80
    ),

    # C-H 键
    "C-H_aliphatic": BondDissociationEnergy(
        "C-H", "aliphatic", 413.0, 98.8, 0.10
    ),
    "C-H_allylic": BondDissociationEnergy(
        "C-H", "allylic", 368.0, 88.0, 0.35
    ),
    "C-H_benzylic": BondDissociationEnergy(
        "C-H", "benzylic", 355.0, 85.0, 0.40
    ),
    "C-H_aldehyde": BondDissociationEnergy(
        "C-H", "aldehyde", 364.0, 87.0, 0.40
    ),

    # O-H 键
    "O-H_alcohol": BondDissociationEnergy(
        "O-H", "alcohol", 439.0, 105.0, 0.05
    ),
    "O-H_carboxylic": BondDissociationEnergy(
        "O-H", "carboxylic_acid", 460.0, 110.0, 0.05
    ),
    "O-H_phenolic": BondDissociationEnergy(
        "O-H", "phenol", 368.0, 88.0, 0.30
    ),
}


def get_bde_by_type(bond_type: str) -> List[BondDissociationEnergy]:
    """获取特定键类型的解离能"""
    return [bde for bde in BOND_DISSOCIATION_ENERGIES.values()
            if bde.bond_type == bond_type]


def get_most_labile_bonds() -> List[BondDissociationEnergy]:
    """获取最易断裂的键(按裂解倾向降序)"""
    return sorted(BOND_DISSOCIATION_ENERGIES.values(),
                  key=lambda x: x.fragmentation_likelihood, reverse=True)


# ─────────────────────────────────────────────────────────────────
# 诊断离子数据库
# ─────────────────────────────────────────────────────────────────

DIAGNOSTIC_IONS: Dict[str, DiagnosticIon] = {
    # ── 通用诊断离子 ──
    "tropylium": DiagnosticIon(
        mz=91.0548, formula="C7H7+", name="䓬离子(Tropylium)",
        structural_origin="苄基裂解后重排产生的七元环芳香阳离子",
        functional_group="benzyl", rule_id="BENZYLIC_CLEAVAGE",
        polarity="positive", intensity_characteristic="high",
    ),
    "phenyl": DiagnosticIon(
        mz=77.0391, formula="C6H5+", name="苯基阳离子(Phenyl)",
        structural_origin="苯环直接离子化或䓬离子脱乙炔",
        functional_group="aromatic", rule_id="BENZYLIC_CLEAVAGE",
        polarity="positive", intensity_characteristic="medium",
    ),
    "acylium": DiagnosticIon(
        mz=43.0184, formula="C2H3O+", name="乙酰基阳离子(Acetyl)",
        structural_origin="羰基α-裂解产生",
        functional_group="ketone", rule_id="ALPHA_CLEAVAGE_CARBONYL",
        polarity="positive", intensity_characteristic="high",
    ),
    "iminium_m1": DiagnosticIon(
        mz=30.0344, formula="CH4N+", name="亚胺阳离子",
        structural_origin="胺α-裂解产生的CH2=NH2+",
        functional_group="amine", rule_id="ALPHA_CLEAVAGE_AMINE",
        polarity="positive", intensity_characteristic="medium",
    ),

    # ── 黄酮类诊断离子 ──
    "flavonoid_A_ring_152": DiagnosticIon(
        mz=152.0473, formula="C8H8O3+", name="黄酮A环碎片(^1,3A+,5,7-diOH)",
        structural_origin="5,7-二羟基黄酮A环rDA裂解产物",
        functional_group="flavonoid", rule_id="FLAVONOID_RDA_C_RING",
        polarity="positive", intensity_characteristic="high",
    ),
    "flavonoid_A_ring_153": DiagnosticIon(
        mz=153.0188, formula="C7H5O4+", name="黄酮醇A环碎片(^1,3A+,5,7-diOH-4-keto)",
        structural_origin="黄酮醇(如槲皮素)A环rDA裂解产物",
        functional_group="flavonoid", rule_id="FLAVONOID_RDA_C_RING",
        polarity="positive", intensity_characteristic="high",
    ),
    "flavonoid_B_ring_120": DiagnosticIon(
        mz=120.0575, formula="C8H8O+", name="黄酮B环碎片(^1,3B+,无OH)",
        structural_origin="B环rDA裂解产物(无羟基取代时)",
        functional_group="flavonoid", rule_id="FLAVONOID_RDA_C_RING",
        polarity="positive", intensity_characteristic="medium",
    ),
    "flavonoid_B_ring_137": DiagnosticIon(
        mz=137.0239, formula="C7H5O3+", name="黄酮B环碎片(^1,3B+,3\',4\'-diOH)",
        structural_origin="B环rDA裂解产物(3',4'-二羟基取代,如槲皮素/木犀草素)",
        functional_group="flavonoid", rule_id="FLAVONOID_RDA_C_RING",
        polarity="positive", intensity_characteristic="high",
    ),

    # ── 生物碱类诊断离子 ──
    "quinolizidine_m150": DiagnosticIon(
        mz=150.1283, formula="C10H16N+", name="喹嗪啶特征离子",
        structural_origin="喹嗪啶类生物碱rDA开环",
        functional_group="alkaloid", rule_id="RETRO_DIELS_ALDER",
        polarity="positive", intensity_characteristic="high",
    ),
    "caffeine_m138": DiagnosticIon(
        mz=138.0667, formula="C6H8N3O+", name="咖啡因特征碎片",
        structural_origin="黄嘌呤生物碱失去CH3NCO后的碎片(咪唑环rDA开环产物)",
        functional_group="alkaloid", rule_id="RETRO_DIELS_ALDER",
        polarity="positive", intensity_characteristic="high",
    ),
    "caffeine_m110": DiagnosticIon(
        mz=110.0713, formula="C5H8N3+", name="咖啡因二级碎片",
        structural_origin="138碎片进一步丢失CO产生",
        functional_group="alkaloid", rule_id="ALPHA_CLEAVAGE_CARBONYL",
        polarity="positive", intensity_characteristic="medium",
        references=("BMRB FIO00571: caffeine LC-ESI-QTOF MS2",),
    ),
    "caffeine_m83": DiagnosticIon(
        mz=83.0604, formula="C4H7N2+", name="咖啡因三级碎片",
        structural_origin="110碎片失去HCN产生(咪唑/嘧啶特征)",
        functional_group="alkaloid", rule_id="ALPHA_CLEAVAGE_AMINE",
        polarity="positive", intensity_characteristic="low",
    ),

    # ── 萜类/甾体诊断离子 ──
    "steroid_m255": DiagnosticIon(
        mz=255.2113, formula="C19H27+", name="甾体D环裂解",
        structural_origin="甾体D环侧链丢失",
        functional_group="steroid", rule_id="ALLYLIC_CLEAVAGE",
        polarity="positive", intensity_characteristic="high",
    ),

    # ── 糖类 ──
    "oxocarbenium_163": DiagnosticIon(
        mz=163.0606, formula="C6H11O5+", name="己糖氧碳鎓离子",
        structural_origin="糖苷键断裂产生的糖基正离子",
        functional_group="glycoside", rule_id="GLYCOSIDIC_CLEAVAGE",
        polarity="positive", intensity_characteristic="high",
    ),

    # ── 核苷N-糖苷键裂解后质子化碱基 ──
    "nucleoside_cytosine_112": DiagnosticIon(
        mz=112.0505, formula="C4H6N3O+", name="质子化胞嘧啶",
        structural_origin="胞苷N-糖苷键断裂后形成的质子化胞嘧啶",
        functional_group="nucleoside", rule_id="RIBONUCLEOSIDE_GLYCOSIDIC_CLEAVAGE",
        polarity="positive", intensity_characteristic="high",
        required_smarts=("Nc1ccn([C;R])c(=O)n1",),
        allowed_adducts=("[M+H]+",),
        references=("DOI:10.1021/acs.chemrestox.3c00041", "PMCID:PMC11587168"),
        evidence_level="curated_primary",
    ),
    "nucleoside_uracil_113": DiagnosticIon(
        mz=113.0346, formula="C4H5N2O2+", name="质子化尿嘧啶",
        structural_origin="尿苷N-糖苷键断裂后形成的质子化尿嘧啶",
        functional_group="nucleoside", rule_id="RIBONUCLEOSIDE_GLYCOSIDIC_CLEAVAGE",
        polarity="positive", intensity_characteristic="high",
        required_smarts=("O=c1ccn([C;R])c(=O)[nH]1",),
        allowed_adducts=("[M+H]+",),
        references=("DOI:10.1021/acs.chemrestox.3c00041", "PMCID:PMC11587168"),
        evidence_level="curated_primary",
    ),
    "nucleoside_thymine_127": DiagnosticIon(
        mz=127.0502, formula="C5H7N2O2+", name="质子化胸腺嘧啶",
        structural_origin="胸苷N-糖苷键断裂后形成的质子化胸腺嘧啶",
        functional_group="nucleoside", rule_id="DEOXYRIBONUCLEOSIDE_GLYCOSIDIC_CLEAVAGE",
        polarity="positive", intensity_characteristic="high",
        required_smarts=("Cc1cn([C;R])c(=O)[nH]c1=O",),
        allowed_adducts=("[M+H]+",),
        references=("PMCID:PMC11587168",), evidence_level="curated_primary",
    ),
    "nucleoside_adenine_136": DiagnosticIon(
        mz=136.0618, formula="C5H6N5+", name="质子化腺嘌呤",
        structural_origin="腺苷N-糖苷键断裂后形成的质子化腺嘌呤",
        functional_group="nucleoside", rule_id="RIBONUCLEOSIDE_GLYCOSIDIC_CLEAVAGE",
        polarity="positive", intensity_characteristic="high",
        required_smarts=("Nc1ncnc2n([C;R])cnc12",),
        allowed_adducts=("[M+H]+",),
        references=("DOI:10.1021/acs.chemrestox.3c00041", "PMCID:PMC11587168"),
        evidence_level="curated_primary",
    ),
    "nucleoside_guanine_152": DiagnosticIon(
        mz=152.0567, formula="C5H6N5O+", name="质子化鸟嘌呤",
        structural_origin="鸟苷N-糖苷键断裂后形成的质子化鸟嘌呤",
        functional_group="nucleoside", rule_id="RIBONUCLEOSIDE_GLYCOSIDIC_CLEAVAGE",
        polarity="positive", intensity_characteristic="high",
        required_smarts=("Nc1nc2c(ncn2[C;R])c(=O)[nH]1",),
        allowed_adducts=("[M+H]+",),
        references=("DOI:10.1021/acs.chemrestox.3c00041", "PMCID:PMC11587168"),
        evidence_level="curated_primary",
    ),

    # ── 通用碎片离子(覆盖更广的化合物类型) ──
    "protonated_water_19": DiagnosticIon(
        mz=19.0184, formula="H3O+", name="水合氢离子",
        structural_origin="脱水或质子化水",
        functional_group="alcohol", rule_id="DEHYDRATION",
        polarity="positive", intensity_characteristic="low",
    ),
    "formate_45": DiagnosticIon(
        mz=44.9982, formula="CHO2-", name="甲酸根离子",
        structural_origin="甲酸酯或羧酸在负离子模式下产生",
        functional_group="carboxylic_acid", rule_id="DECARBOXYLATION",
        polarity="negative", intensity_characteristic="high",
    ),
    "acetate_59": DiagnosticIon(
        mz=59.0139, formula="C2H3O2-", name="乙酸根离子",
        structural_origin="乙酸酯裂解产生",
        functional_group="acetate_ester", rule_id="ALPHA_CLEAVAGE_CARBONYL",
        polarity="negative", intensity_characteristic="high",
    ),
    "benzoate_121": DiagnosticIon(
        mz=121.0297, formula="C7H5O2-", name="苯甲酸根离子",
        structural_origin="苯甲酸酯类裂解",
        functional_group="benzoate", rule_id="DECARBOXYLATION",
        polarity="negative", intensity_characteristic="high",
    ),
    "sulfonate_80": DiagnosticIon(
        mz=79.9574, formula="SO3-", name="硫酸根碎片",
        structural_origin="磺酸/硫酸酯裂解",
        functional_group="sulfonic_acid", rule_id="DEHYDRATION",
        polarity="negative", intensity_characteristic="high",
    ),
    "sulfonate_81": DiagnosticIon(
        mz=80.9652, formula="HSO3-", name="亚硫酸氢根",
        structural_origin="磺酸类化合物特征碎片",
        functional_group="sulfonic_acid", rule_id="DEHYDRATION",
        polarity="negative", intensity_characteristic="high",
    ),
    "deprotonated_phenol_93": DiagnosticIon(
        mz=93.0346, formula="C6H5O-", name="酚氧负离子",
        structural_origin="酚类去质子化或苯氧基裂解",
        functional_group="phenol", rule_id="BENZYLIC_CLEAVAGE",
        polarity="negative", intensity_characteristic="high",
    ),
    "deprotonated_cinnamic_147": DiagnosticIon(
        mz=147.0446, formula="C9H7O2-", name="肉桂酸根离子",
        structural_origin="酚酸/肉桂酸衍生物",
        functional_group="phenolic_acid", rule_id="DECARBOXYLATION",
        polarity="negative", intensity_characteristic="high",
    ),
    "deprotonated_caffeic_179": DiagnosticIon(
        mz=179.0350, formula="C9H7O4-", name="咖啡酸根离子",
        structural_origin="咖啡酸/绿原酸去质子化",
        functional_group="phenolic_acid", rule_id="DECARBOXYLATION",
        polarity="negative", intensity_characteristic="high",
    ),
    "oxalate_89": DiagnosticIon(
        mz=88.9920, formula="C2HO4-", name="草酸根离子",
        structural_origin="草酸酯或二元羧酸裂解",
        functional_group="carboxylic_acid", rule_id="DECARBOXYLATION",
        polarity="negative", intensity_characteristic="medium",
    ),
    "nitro_46": DiagnosticIon(
        mz=45.9929, formula="NO2-", name="硝酸根碎片",
        structural_origin="硝基化合物裂解",
        functional_group="nitro", rule_id="ALPHA_CLEAVAGE_CARBONYL",
        polarity="negative", intensity_characteristic="high",
    ),
    "phosphate_97": DiagnosticIon(
        mz=96.9696, formula="H2PO4-", name="磷酸二氢根",
        structural_origin="磷酸酯/磷脂裂解",
        functional_group="phosphate", rule_id="DEHYDRATION",
        polarity="negative", intensity_characteristic="high",
    ),
    "amino_acid_decarb_74": DiagnosticIon(
        mz=74.0242, formula="C2H4NO2-", name="甘氨酸碎片",
        structural_origin="氨基酸脱羧/脱水特征碎片",
        functional_group="amino_acid", rule_id="DECARBOXYLATION",
        polarity="negative", intensity_characteristic="medium",
    ),
    "fatty_acid_fragment_171": DiagnosticIon(
        mz=171.1386, formula="C10H19O2-", name="脂肪酸C10碎片",
        structural_origin="长链脂肪酸电荷远程裂解",
        functional_group="fatty_acid", rule_id="CHARGE_REMOTE_FRAGMENTATION",
        polarity="negative", intensity_characteristic="medium",
    ),

    # ── 经文献核对的磷脂/鞘脂类别离子（必须通过结构或原始分类守卫） ──
    "phosphocholine_184": DiagnosticIon(
        mz=184.0733, formula="C5H15NO4P+", name="磷酸胆碱头基离子",
        structural_origin="PC/LPC/鞘磷脂的磷酸胆碱极性头基",
        functional_group="lipid", rule_id="PHOSPHOLIPID_HEADGROUP_CLEAVAGE",
        polarity="positive", intensity_characteristic="high",
        compatible_class_terms=("phosphatidylcholine", "lysophosphatidylcholine", "sphingomyelin"),
        required_smarts=(
            "P(=O)(O)OCC[N+](C)(C)C",
            "P(=O)([O-])OCC[N+](C)(C)C",
            "P(=O)(O)OCCN(C)(C)C",
        ),
        allowed_adducts=("[M+H]+",),
        references=("PMID:32423109", "PMCID:PMC7287777", "PMCID:PMC5697860"),
        evidence_level="curated_primary",
    ),
    "sphingoid_d18_1_264": DiagnosticIon(
        mz=264.2686, formula="C18H34N+", name="d18:1鞘氨醇双脱水离子",
        structural_origin="d18:1鞘氨醇长链碱骨架双脱水形成的碳正离子",
        functional_group="lipid", rule_id="SPHINGOID_BASE_CLEAVAGE",
        polarity="positive", intensity_characteristic="high",
        compatible_class_terms=("sphingolipid", "ceramide", "sphingomyelin", "hexosylceramide"),
        required_smarts=("NCC(O)C=C", "NC(CO)C(O)C=C"),
        allowed_adducts=("[M+H]+",),
        references=("PMID:11070861", "PMCID:PMC2724058"),
        evidence_level="curated_primary",
    ),
    "glycerophosphocholine_224": DiagnosticIon(
        mz=224.0693, formula="C7H15NO5P-", name="甘油磷酸胆碱负离子",
        structural_origin="PC/LPC负离子模式的去甲基化甘油磷酸胆碱头基",
        functional_group="lipid", rule_id="PHOSPHOLIPID_HEADGROUP_CLEAVAGE",
        polarity="negative", intensity_characteristic="high",
        compatible_class_terms=("phosphatidylcholine", "lysophosphatidylcholine"),
        required_smarts=(
            "P(=O)(O)OCC[N+](C)(C)C",
            "P(=O)([O-])OCC[N+](C)(C)C",
        ),
        allowed_adducts=("[M+FA-H]-", "[M+Hac-H]-"),
        references=("PMID:32423109", "PMCID:PMC7287777"),
        evidence_level="curated_primary",
    ),
    "phosphoethanolamine_196": DiagnosticIon(
        mz=196.0380, formula="C5H11NO5P-", name="甘油磷酸乙醇胺负离子",
        structural_origin="PE/LPE负离子模式的磷酸乙醇胺头基",
        functional_group="lipid", rule_id="PHOSPHOLIPID_HEADGROUP_CLEAVAGE",
        polarity="negative", intensity_characteristic="high",
        compatible_class_terms=("phosphatidylethanolamine", "lysophosphatidylethanolamine"),
        required_smarts=("P(=O)(O)OCC[NX3;H1,H2]", "P(=O)([O-])OCC[NX3;H1,H2]"),
        allowed_adducts=("[M-H]-",),
        references=("PMID:32423109", "DOI:10.1002/rcm.9527"),
        evidence_level="curated_primary",
    ),
    "inositol_phosphate_241": DiagnosticIon(
        mz=241.0117, formula="C6H10O8P-", name="肌醇磷酸头基负离子",
        structural_origin="PI负离子模式的肌醇磷酸极性头基",
        functional_group="lipid", rule_id="PHOSPHOLIPID_HEADGROUP_CLEAVAGE",
        polarity="negative", intensity_characteristic="high",
        compatible_class_terms=("phosphatidylinositol", "lysophosphatidylinositol"),
        required_smarts=(
            "P(=O)(O)OC1C(O)C(O)C(O)C(O)C1O",
            "P(=O)([O-])OC1C(O)C(O)C(O)C(O)C1O",
        ),
        allowed_adducts=("[M-H]-",),
        references=("PMID:32423109", "PMCID:PMC7287777"),
        evidence_level="curated_primary",
    ),
    "steroid_m271": DiagnosticIon(
        mz=271.2113, formula="C19H27O-", name="甾体D环碎片(负离子)",
        structural_origin="甾体D环侧链丢失(负离子模式)",
        functional_group="steroid", rule_id="ALLYLIC_CLEAVAGE",
        polarity="negative", intensity_characteristic="high",
    ),
    "steroid_m273": DiagnosticIon(
        mz=273.2269, formula="C19H29O-", name="甾体特征碎片",
        structural_origin="甾体环裂解产物",
        functional_group="steroid", rule_id="ALLYLIC_CLEAVAGE",
        polarity="positive", intensity_characteristic="high",
    ),
    "terpene_m137": DiagnosticIon(
        mz=137.1328, formula="C10H17+", name="单萜特征碎片",
        structural_origin="单萜环裂解产物",
        functional_group="terpenoid", rule_id="ALLYLIC_CLEAVAGE",
        polarity="positive", intensity_characteristic="high",
    ),
    "terpene_m81": DiagnosticIon(
        mz=81.0704, formula="C6H9+", name="环己烯特征碎片",
        structural_origin="萜类/甾体环裂解碎片",
        functional_group="terpenoid", rule_id="ALLYLIC_CLEAVAGE",
        polarity="positive", intensity_characteristic="high",
    ),
    "terpene_m95": DiagnosticIon(
        mz=95.0861, formula="C7H11+", name="萜类环碎片",
        structural_origin="倍半萜环裂解",
        functional_group="terpenoid", rule_id="ALLYLIC_CLEAVAGE",
        polarity="positive", intensity_characteristic="medium",
    ),
    "terpene_m109": DiagnosticIon(
        mz=109.1017, formula="C8H13+", name="萜类链碎片",
        structural_origin="萜类侧链裂解",
        functional_group="terpenoid", rule_id="ALLYLIC_CLEAVAGE",
        polarity="positive", intensity_characteristic="medium",
    ),
}


def get_diagnostic_ions_by_fg(functional_group: str) -> List[DiagnosticIon]:
    """获取特定官能团的诊断离子"""
    return [ion for ion in DIAGNOSTIC_IONS.values()
            if ion.functional_group == functional_group]


def _diagnostic_ion_scope_matches(
    ion: DiagnosticIon,
    candidate_smiles: str = "",
    raw_classification: str = "",
    query_adduct: str = "",
    query_collision_energy: Any = None,
) -> bool:
    """Require class/structure support for narrow-scope curated ions."""
    raw_lower = str(raw_classification or "").lower()
    class_match = bool(
        ion.compatible_class_terms
        and any(term.lower() in raw_lower for term in ion.compatible_class_terms)
    )
    structure_match = bool(
        ion.required_smarts
        and smiles_matches_any(str(candidate_smiles or ""), ion.required_smarts)
    )
    if ion.required_smarts_all:
        structure_all_match = smiles_matches_all(
            str(candidate_smiles or ""), ion.required_smarts_all
        )
        structure_match = (
            structure_match and structure_all_match
            if ion.required_smarts else structure_all_match
        )
    if (ion.compatible_class_terms or ion.required_smarts or ion.required_smarts_all) and not (
        class_match or structure_match
    ):
        return False
    if ion.allowed_adducts:
        observed_adduct = normalize_adduct(query_adduct)
        allowed_adducts = {normalize_adduct(item) for item in ion.allowed_adducts}
        if not observed_adduct or observed_adduct not in allowed_adducts:
            return False
    if ion.collision_energy_min is not None or ion.collision_energy_max is not None:
        value, unit = collision_energy_value_unit(query_collision_energy)
        expected_unit = str(ion.collision_energy_unit or "").strip()
        if value is None or not unit or not expected_unit or unit.lower() != expected_unit.lower():
            return False
        if ion.collision_energy_min is not None and value < ion.collision_energy_min:
            return False
        if ion.collision_energy_max is not None and value > ion.collision_energy_max:
            return False
    return True


def match_diagnostic_ion(mz: float, tolerance: float = 0.02,
                        classification: str = "",
                        strict_fg: bool = False,
                        polarity: str = "",
                        candidate_smiles: str = "",
                        raw_classification: str = "",
                        query_adduct: str = "",
                        query_collision_energy: Any = None) -> List[DiagnosticIon]:
    """匹配诊断离子

    Args:
        mz: 观测到的碎片 m/z
        tolerance: 质量匹配容差(Da)
        classification: 化合物分类(用于过滤不相关的诊断离子)
        strict_fg: 是否启用官能团严格过滤
        polarity: 查询谱图离子模式；已知时严格排除相反极性规则
        candidate_smiles: 候选结构，用于类别特异规则的 SMARTS 守卫
        raw_classification: 未归一化分类，用于磷脂等子类规则守卫

    When strict_fg=True, only return diagnostic ions whose functional_group
    is compatible with the compound's classification.
    """
    # 标准化分类并修复 truthiness Bug
    effective_cls = normalize_classification(classification)
    wanted_polarity = normalize_polarity(polarity)

    results = []
    for ion in DIAGNOSTIC_IONS.values():
        ion_polarity = normalize_polarity(ion.polarity)
        if wanted_polarity and ion_polarity and ion_polarity != wanted_polarity:
            continue
        # 数据质量门：公式和 m/z 自相矛盾的条目不得参与生产评分。
        if not formula_mz_consistent(ion.formula, ion.mz):
            continue
        if not _diagnostic_ion_scope_matches(
            ion,
            candidate_smiles=candidate_smiles,
            raw_classification=raw_classification,
            query_adduct=query_adduct,
            query_collision_energy=query_collision_energy,
        ):
            continue
        if abs(ion.mz - mz) <= tolerance:
            # 即使 effective_cls == "Unknown"，也使用保守过滤
            if strict_fg:
                if not _is_diagnostic_ion_compatible(ion, effective_cls):
                    continue
            results.append(ion)
    return results


def _is_diagnostic_ion_compatible(ion: DiagnosticIon, classification: str) -> bool:
    """检查诊断离子是否与化合物分类兼容

    核心规则: 诊断离子标记的功能团必须与化合物的已知分类有合理的化学关系。
    例如: 黄酮B环碎片不应出现在酚酸类化合物的分析中。
    """
    cls_lower = classification.lower().strip()
    ion_fg = ion.functional_group.lower()

    # 分类兼容映射
    # 每个 classification 键允许的 diagnostic ion functional_group 值
    _compatibility = {
        "flavonoid": {"flavonoid", "phenol", "aromatic", "benzyl", "glycoside",
                       "ketone", "alcohol"},
        "alkaloid": {"alkaloid", "amine", "N_heterocycle", "purine", "aromatic"},
        "terpenoid": {"terpenoid", "steroid", "alkene", "alcohol", "ketone"},
        "steroid": {"steroid", "terpenoid", "alkene", "alcohol", "ketone"},
        "phenolic_acid": {"phenolic_acid", "carboxylic_acid", "phenol", "aromatic",
                          "benzyl", "benzoate", "aniline"},
        "glycoside": {"glycoside", "sugar", "alcohol", "phenol", "carboxylic_acid",
                       "flavonoid"},
        "sugar": {"glycoside", "sugar", "alcohol", "carboxylic_acid"},
        "fatty_acid": {"fatty_acid", "lipid", "carboxylic_acid"},
        "lipid": {"lipid", "fatty_acid", "carboxylic_acid", "phosphate"},
        "coumarin": {"coumarin", "carboxylic_acid", "phenol", "aromatic", "benzyl"},
        "quinone": {"quinone", "ketone", "phenol"},
        "lignan": {"phenol", "aromatic", "benzyl", "alcohol"},
        "amino_acid": {"amine", "amino_acid", "amide", "carboxylic_acid"},
        "benzenoid": {"aromatic", "benzyl", "phenol", "aniline", "benzoate"},
        "sulfonate": {"sulfonic_acid"},
        "saponin": {"glycoside", "sugar", "terpenoid", "steroid", "alcohol"},
        "nucleoside": {"nucleoside", "purine", "pyrimidine"},
    }

    # 找到匹配的分类
    allowed_fgs = None
    for cls_key, fgs in _compatibility.items():
        if cls_key in cls_lower:
            allowed_fgs = fgs
            break

    if allowed_fgs is None:
        # 未知分类: 保守策略 — 阻断高度类特异的诊断离子
        # 这些离子只在特定化合物类别中出现，不应泛化匹配
        _class_specific_fgs = {
            "flavonoid", "alkaloid", "steroid", "glycoside",
            "sugar", "fatty_acid", "lipid", "purine",
            "pyrimidine", "nucleoside", "caffeine", "terpenoid",
        }
        if ion_fg in _class_specific_fgs:
            return False
        return True

    # 检查诊断离子的 functional_group 是否在允许列表中
    if ion_fg in allowed_fgs:
        return True

    return False


# ─────────────────────────────────────────────────────────────────
# 官能团特征裂解模式映射
# ─────────────────────────────────────────────────────────────────

FUNCTIONAL_GROUP_SIGNATURES: Dict[str, Dict[str, Any]] = {
    "flavonoid": {
        "characteristic_losses": [162.0528, 146.0579, 132.0423, 120.0575, 152.0473, 153.0188, 137.0239],
        "diagnostic_ions": [152.0473, 153.0188, 137.0239, 120.0575, 121.0290],
        "rules": ["FLAVONOID_RDA_C_RING", "GLYCOSIDIC_CLEAVAGE", "RETRO_DIELS_ALDER",
                  "C_GLYCOSIDE_CROSS_RING"],
        "typical_neutral_losses": ["C6H10O5", "C8H8O", "CO", "H2O", "CH3", "CO2"],
        "fragmentation_pathway": "glycosidic_cleavage → rDA_C_ring → A_ring/B_ring_ions → CO/H2O_losses",
    },
    "alkaloid": {
        "characteristic_losses": [57.0215, 43.0058, 31.0422, 27.0109, 17.0265, 15.0235, 41.0265],
        "diagnostic_ions": [30.0344, 58.0657, 138.0667, 110.0713, 150.1283],
        "rules": ["ALPHA_CLEAVAGE_AMINE", "RETRO_DIELS_ALDER"],
        "typical_neutral_losses": ["CH3NCO", "HNCO", "CH3N", "HCN", "NH3", "CH3", "CH3CN"],
        "fragmentation_pathway": "imidazole/purine_opening → CH3NCO/HNCO_loss → HCN_loss → stable_cation",
    },
    "terpenoid": {
        "characteristic_losses": [56.0626, 42.0470, 28.0313],
        "diagnostic_ions": [81.0704, 95.0861, 109.1017, 255.2113],
        "rules": ["ALLYLIC_CLEAVAGE", "RETRO_DIELS_ALDER", "MCLAFFERTY_REARRANGEMENT"],
        "typical_neutral_losses": ["C4H8", "C3H6", "C2H4", "H2O"],
        "fragmentation_pathway": "ring_opening → sequential_alkene_loss → stable_cation",
    },
    "phenolic_acid": {
        "characteristic_losses": [44.0, 18.0106, 46.0055],
        "diagnostic_ions": [109.0290, 121.0290, 163.0395],
        "rules": ["DECARBOXYLATION", "DEHYDRATION", "ORTHO_EFFECT"],
        "typical_neutral_losses": ["CO2", "H2O", "HCOOH", "CO"],
        "fragmentation_pathway": "decarboxylation → dehydration → CO_loss",
    },
    "glycoside": {
        "characteristic_losses": [162.0528, 146.0579, 132.0423, 176.0321,
                                  60.0211, 90.0317, 120.0423],
        "diagnostic_ions": [163.0606, 147.0657, 133.0501],
        "rules": ["GLYCOSIDIC_CLEAVAGE", "C_GLYCOSIDE_CROSS_RING"],
        "typical_neutral_losses": ["C6H10O5", "C6H10O4", "C5H8O4", "C6H8O6"],
        "fragmentation_pathway": "sequential_sugar_loss → aglycone → aglycone_fragmentation",
    },
    "sugar": {
        "characteristic_losses": [18.0106, 36.0212, 30.0106],
        "diagnostic_ions": [163.0606, 147.0657, 133.0501, 145.0500, 129.0552, 115.0395],
        "rules": ["DEHYDRATION", "INDUCTIVE_CLEAVAGE", "MCLAFFERTY_REARRANGEMENT"],
        "typical_neutral_losses": ["H2O", "2H2O", "CH2O", "C2H4O2", "CO", "CO2"],
        "fragmentation_pathway": "ring_opening → oxocarbenium_ion → sequential_H2O_loss → C3-C5_fragments",
    },
    "fatty_acid": {
        "characteristic_losses": [18.0106, 44.0],
        "diagnostic_ions": [],
        "rules": ["CHARGE_REMOTE_FRAGMENTATION", "DEHYDRATION", "DECARBOXYLATION"],
        "typical_neutral_losses": ["H2O", "CO2"],
        "fragmentation_pathway": "CO2_loss_neg → CRF_series → alkene_ladder",
    },
    "coumarin": {
        "characteristic_losses": [44.0, 28.0, 18.0106],
        "diagnostic_ions": [147.0446, 119.0497, 91.0548],
        "rules": ["DECARBOXYLATION", "MCLAFFERTY_REARRANGEMENT"],
        "typical_neutral_losses": ["CO2", "CO", "H2O"],
        "fragmentation_pathway": "lactone_opening → CO2_loss → CO_loss",
    },
    "quinone": {
        "characteristic_losses": [28.0, 56.0],
        "diagnostic_ions": [],
        "rules": ["ALPHA_CLEAVAGE_CARBONYL", "MCLAFFERTY_REARRANGEMENT"],
        "typical_neutral_losses": ["CO", "2CO"],
        "fragmentation_pathway": "sequential_CO_loss",
    },
    "steroid": {
        "characteristic_losses": [18.0106, 42.0470, 56.0626, 28.0313],
        "diagnostic_ions": [255.2113, 271.2113, 273.2269, 81.0704, 95.0861, 109.1017],
        "rules": ["ALLYLIC_CLEAVAGE", "DEHYDRATION", "MCLAFFERTY_REARRANGEMENT"],
        "typical_neutral_losses": ["H2O", "C3H6", "C4H8", "C2H4"],
        "fragmentation_pathway": "side_chain_loss → ring_opening → sequential_alkene_loss",
    },
    "amino_acid": {
        "characteristic_losses": [18.0106, 44.0095, 17.0265, 46.0055],
        "diagnostic_ions": [70.0657, 72.0813, 86.0970, 101.1079, 110.0718,
                             120.0813, 129.1140, 136.0762, 159.0922],
        "rules": ["DECARBOXYLATION", "DEHYDRATION", "ALPHA_CLEAVAGE_AMINE",
                  "ALPHA_AMINO_ACID_46_LOSS"],
        "typical_neutral_losses": ["H2O", "CO2", "NH3", "HCOOH"],
        "fragmentation_pathway": "deamination → decarboxylation → dehydration",
    },
    "lipid": {
        "characteristic_losses": [18.0106, 43.9898, 141.0191],
        "diagnostic_ions": [171.1386],
        "rules": ["CHARGE_REMOTE_FRAGMENTATION", "DEHYDRATION", "DECARBOXYLATION",
                  "PHOSPHATIDYLETHANOLAMINE_141_LOSS"],
        "typical_neutral_losses": ["H2O", "CO2"],
        "fragmentation_pathway": "CRF_series → alkene_ladder → CO2_loss",
    },
    "benzenoid": {
        "characteristic_losses": [28.0, 18.0106, 44.0095],
        "diagnostic_ions": [91.0548, 77.0391, 93.0346],
        "rules": ["BENZYLIC_CLEAVAGE", "DEHYDRATION", "DECARBOXYLATION"],
        "typical_neutral_losses": ["CO", "H2O", "CO2"],
        "fragmentation_pathway": "benzylic_cleavage → CO_loss → ring_fragmentation",
    },
    "sulfonate": {
        "characteristic_losses": [79.9568, 80.0646, 18.0106],
        "diagnostic_ions": [79.9574, 80.9652],
        "rules": ["DEHYDRATION", "SULFATE_SO3_LOSS"],
        "typical_neutral_losses": ["SO3", "H2SO4", "H2O"],
        "fragmentation_pathway": "SO3_loss → H2O_loss",
    },
    "lignan": {
        "characteristic_losses": [30.0106, 31.0422, 18.0106, 15.0235],
        "diagnostic_ions": [137.0239, 151.0395],
        "rules": ["DEHYDRATION", "INDUCTIVE_CLEAVAGE", "ORTHO_EFFECT"],
        "typical_neutral_losses": ["OCH2", "OCH3", "H2O", "CH3", "CO"],
        "fragmentation_pathway": "OCH2/OCH3_loss → H2O_loss → biphenyl_cleavage",
    },
    "saponin": {
        "characteristic_losses": [162.0528, 146.0579, 132.0423, 176.0321, 18.0106],
        "diagnostic_ions": [109.1017, 189.1643, 203.1800],
        "rules": ["GLYCOSIDIC_CLEAVAGE", "DEHYDRATION"],
        "typical_neutral_losses": ["C6H10O5", "C6H10O4", "C5H8O4", "C6H8O6", "H2O"],
        "fragmentation_pathway": "sequential_sugar_loss → aglycone_ion → aglycone_fragmentation",
    },
    "nucleoside": {
        "characteristic_losses": [132.0423, 116.0473],
        "diagnostic_ions": [112.0505, 113.0346, 127.0502, 136.0618, 152.0567],
        "rules": ["RIBONUCLEOSIDE_GLYCOSIDIC_CLEAVAGE",
                  "DEOXYRIBONUCLEOSIDE_GLYCOSIDIC_CLEAVAGE"],
        "typical_neutral_losses": ["C5H8O4", "C5H8O3"],
        "fragmentation_pathway": "N_glycosidic_bond_cleavage → protonated_nucleobase",
    },
}


def get_fg_signature(fg: str) -> Optional[Dict[str, Any]]:
    """获取官能团的特征裂解签名"""
    return FUNCTIONAL_GROUP_SIGNATURES.get(fg)


# ─────────────────────────────────────────────────────────────────
# 工具函数
# ─────────────────────────────────────────────────────────────────

def compute_possible_neutral_losses(
    precursor_mz: float, fragment_mzs: List[float], tolerance: float = 0.5,
    classification: str = "", strict_fg: bool = False,
    polarity: str = "", candidate_smiles: str = "",
    query_adduct: str = "", query_collision_energy: Any = None,
) -> List[Tuple[float, Optional[NeutralLoss]]]:
    """
    从前体离子与碎片离子差值推算可能的中性丢失

    Args:
        precursor_mz: 前体离子 m/z
        fragment_mzs: 碎片离子 m/z 列表
        tolerance: 质量匹配容差(Da)
        classification: 化合物分类(用于过滤不兼容的NL)
        strict_fg: 是否启用官能团严格过滤

    Returns:
        [(质量差, 匹配的中性丢失), ...]
    """
    results = []
    for frag_mz in fragment_mzs:
        diff = precursor_mz - frag_mz
        if diff <= 0:
            continue
        matches = get_neutral_loss_by_mass(
            diff, tolerance,
            classification=classification, strict_fg=strict_fg,
            polarity=polarity, candidate_smiles=candidate_smiles,
            query_adduct=query_adduct,
            query_collision_energy=query_collision_energy,
        )
        if matches:
            results.append((diff, matches[0]))
        else:
            results.append((diff, None))
    return results


def fragmentation_rule_applicable(
    rule: FragmentationRule,
    polarity: str = "",
    candidate_smiles: str = "",
    query_adduct: str = "",
    query_collision_energy: Any = None,
    ionization_method: str = "ESI",
) -> tuple[bool, str]:
    """Apply the structured domain guard for an executable rule.

    Collision-energy ranges are activated only when both rule and query carry
    compatible units.  Bare database numbers are deliberately not treated as
    eV or NCE because public libraries mix both conventions.
    """
    wanted_polarity = normalize_polarity(polarity)
    rule_polarity = normalize_polarity(rule.charge_preference)
    if wanted_polarity and rule_polarity and wanted_polarity != rule_polarity:
        return False, "polarity_mismatch"

    method = str(ionization_method or "").strip().upper()
    allowed_methods = {str(item).strip().upper() for item in rule.ionization_methods if item}
    if method and allowed_methods and method not in allowed_methods:
        return False, "ionization_method_mismatch"

    if rule.allowed_adducts:
        observed_adduct = normalize_adduct(query_adduct)
        allowed_adducts = {normalize_adduct(item) for item in rule.allowed_adducts}
        if not observed_adduct:
            return False, "adduct_unknown"
        if observed_adduct not in allowed_adducts:
            return False, "adduct_mismatch"

    if rule.required_smarts:
        if not candidate_smiles:
            return False, "structure_unknown"
        if not smiles_matches_any(candidate_smiles, rule.required_smarts):
            return False, "structure_mismatch"
    if rule.required_smarts_all:
        if not candidate_smiles:
            return False, "structure_unknown"
        if not smiles_matches_all(candidate_smiles, rule.required_smarts_all):
            return False, "structure_mismatch"

    has_energy_range = (
        rule.collision_energy_min is not None or rule.collision_energy_max is not None
    )
    if has_energy_range:
        value, unit = collision_energy_value_unit(query_collision_energy)
        rule_unit = str(rule.collision_energy_unit or "").strip()
        if value is None or not unit or not rule_unit:
            return False, "collision_energy_not_comparable"
        if unit.lower() != rule_unit.lower():
            return False, "collision_energy_unit_mismatch"
        if rule.collision_energy_min is not None and value < rule.collision_energy_min:
            return False, "collision_energy_below_range"
        if rule.collision_energy_max is not None and value > rule.collision_energy_max:
            return False, "collision_energy_above_range"
    return True, "applicable"


def compute_rule_coverage(
    neutral_losses: List[Tuple[float, Optional[NeutralLoss]]],
    functional_group: str,
    polarity: str = "",
    mass_tolerance: float = 0.05,
    candidate_smiles: str = "",
    query_adduct: str = "",
    query_collision_energy: Any = None,
    ionization_method: str = "ESI",
) -> Dict[str, Any]:
    """
    计算裂解规则覆盖率

    给定观测到的中性丢失和化合物类型，评估有多少已知裂解规则被覆盖

    Returns:
        {
            "total_rules": 该官能团相关规则总数,
            "covered_rules": 被覆盖的规则ID列表,
            "coverage_rate": 覆盖率(0-1),
            "uncovered_rules": 未覆盖的规则ID列表,
            "matched_neutral_losses": 匹配到的中性丢失名称列表,
        }
    """
    sig = get_fg_signature(functional_group)
    if not sig:
        return {"total_rules": 0, "covered_rules": [], "coverage_rate": 0.0,
                "uncovered_rules": [], "matched_neutral_losses": []}

    rule_ids = []
    excluded_rules = {}
    for rule_id in sig["rules"]:
        rule = FRAGMENTATION_RULES.get(rule_id)
        if rule is None:
            continue
        applicable, reason = fragmentation_rule_applicable(
            rule,
            polarity=polarity,
            candidate_smiles=candidate_smiles,
            query_adduct=query_adduct,
            query_collision_energy=query_collision_energy,
            ionization_method=ionization_method,
        )
        if not applicable:
            excluded_rules[rule_id] = reason
            continue
        rule_ids.append(rule_id)
    expected_losses = set(sig["characteristic_losses"])

    matched_nls = set()
    covered_rules = set()

    for diff, nl in neutral_losses:
        if nl:
            matched_nls.add(nl.name)
            for rule_id in rule_ids:
                rule = FRAGMENTATION_RULES.get(rule_id)
                if rule and nl.name in {n.name for n in rule.typical_neutral_losses}:
                    covered_rules.add(rule_id)

    # 也通过质量匹配检测(针对没有精确NL条目的情况)
    for diff, _ in neutral_losses:
        for expected_mass in expected_losses:
            if abs(diff - expected_mass) <= mass_tolerance:
                # 大概匹配,查找相关规则
                for rule_id in rule_ids:
                    rule = FRAGMENTATION_RULES.get(rule_id)
                    if rule:
                        for nl in rule.typical_neutral_losses:
                            if abs(nl.mass - expected_mass) <= mass_tolerance:
                                covered_rules.add(rule_id)

    return {
        "total_rules": len(rule_ids),
        "covered_rules": list(covered_rules),
        "coverage_rate": len(covered_rules) / len(rule_ids) if rule_ids else 0.0,
        "uncovered_rules": [r for r in rule_ids if r not in covered_rules],
        "matched_neutral_losses": list(matched_nls),
        "excluded_rules": excluded_rules,
    }


def estimate_bond_lability(
    smiles: str
) -> Dict[str, Any]:
    """
    估算化合物中各键的裂解倾向性 — 基于BDE数据库做启发式估计

    这是一个基于SMARTS模式的启发式评估，不是量子化学计算。
    用于为碎片预测提供化学合理性先验。

    Returns:
        {
            "labile_bonds": [(bond_description, lability_score), ...],
            "most_labile_region": str,
            "estimated_fragmentation_sites": int,
        }
    """
    try:
        from rdkit import Chem
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {"labile_bonds": [], "most_labile_region": "unknown",
                    "estimated_fragmentation_sites": 0}
    except ImportError:
        return {"labile_bonds": [], "most_labile_region": "unknown",
                "estimated_fragmentation_sites": 0}

    labile_bonds = []
    bond_type_map = {
        Chem.BondType.SINGLE: "single",
        Chem.BondType.DOUBLE: "double",
        Chem.BondType.AROMATIC: "aromatic",
    }

    for bond in mol.GetBonds():
        bt = bond_type_map.get(bond.GetBondType(), "other")
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()

        # 简单启发式：杂原子旁的单键更易断
        lability = 0.5  # baseline for single C-C

        # 杂原子参与的键
        if a1.GetAtomicNum() not in (6, 1) or a2.GetAtomicNum() not in (6, 1):
            if bt == "single":
                lability = 0.75  # C-O, C-N, etc.

        # 苄基位
        if (a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 6 and bt == "single"):
            # Check if either is aromatic neighbor
            for atom in (a1, a2):
                for nbr in atom.GetNeighbors():
                    if nbr.GetIsAromatic():
                        lability = max(lability, 0.8)
                        break

        # α-羰基
        for atom in (a1, a2):
            for nbr in atom.GetNeighbors():
                if nbr.GetAtomicNum() == 8:
                    for nbr2 in nbr.GetNeighbors():
                        if nbr2.GetIdx() == atom.GetIdx():
                            b = mol.GetBondBetweenAtoms(nbr.GetIdx(), atom.GetIdx())
                            if b and b.GetBondType() == Chem.BondType.DOUBLE:
                                lability = max(lability, 0.7)

        # 烯丙基位
        for atom in (a1, a2):
            for nbr in atom.GetNeighbors():
                if nbr.GetAtomicNum() == 6:
                    for nbr2 in nbr.GetNeighbors():
                        if nbr2.GetIdx() != atom.GetIdx():
                            b2 = mol.GetBondBetweenAtoms(nbr.GetIdx(), nbr2.GetIdx())
                            if b2 and b2.GetBondType() == Chem.BondType.DOUBLE:
                                lability = max(lability, 0.85)

        if lability > 0.5:
            labile_bonds.append({
                "bond_idx": bond.GetIdx(),
                "atom1": f"{a1.GetSymbol()}{a1.GetIdx()}",
                "atom2": f"{a2.GetSymbol()}{a2.GetIdx()}",
                "bond_type": bt,
                "lability_score": lability,
            })

    labile_bonds.sort(key=lambda x: x["lability_score"], reverse=True)

    # Determine most labile region
    if labile_bonds:
        top = labile_bonds[0]
        most_labile = f"{top['atom1']}-{top['atom2']} (score={top['lability_score']:.2f})"
    else:
        most_labile = "no_significantly_labile_bonds"

    return {
        "labile_bonds": labile_bonds[:10],  # Top 10 most labile
        "most_labile_region": most_labile,
        "estimated_fragmentation_sites": len(labile_bonds),
    }


def get_rules_summary() -> Dict[str, Any]:
    """获取规则库摘要统计"""
    mechanism_counts = defaultdict(int)
    for rule in FRAGMENTATION_RULES.values():
        mechanism_counts[rule.mechanism_type] += 1

    return {
        "total_rules": len(FRAGMENTATION_RULES),
        "total_neutral_losses": len(COMMON_NEUTRAL_LOSSES),
        "total_diagnostic_ions": len(DIAGNOSTIC_IONS),
        "total_bde_entries": len(BOND_DISSOCIATION_ENERGIES),
        "functional_groups_covered": len(FUNCTIONAL_GROUP_SIGNATURES),
        "mechanism_types": dict(mechanism_counts),
        "rule_confidence_avg": sum(r.confidence_score for r in FRAGMENTATION_RULES.values())
                               / len(FRAGMENTATION_RULES),
    }


def load_compound_class_rules(filepath: str = None) -> Dict[str, Any]:
    """
    加载化合物类别质谱裂解规则库 (JSON格式)

    从 fragmentation_rules_library.json 加载基于母核结构的裂解规则，
    涵盖18+化合物类别、149条裂解规则、精确质量参考表。

    Args:
        filepath: JSON规则库文件路径。默认为当前目录下的 fragmentation_rules_library.json

    Returns:
        包含所有化合物类别裂解规则的字典
    """
    if filepath is None:
        filepath = Path(__file__).parent / "fragmentation_rules_library.json"

    filepath = Path(filepath)
    if not filepath.exists():
        logger.warning(f"Fragmentation rules library not found: {filepath}")
        return {}

    with open(filepath, "r", encoding="utf-8") as f:
        data = json.load(f)

    return data


def get_compound_class_fragments(class_name: str, rules_data: Dict = None) -> List[Dict]:
    """
    获取指定化合物类别的特征碎片离子规则

    Args:
        class_name: 化合物类别名(中文或英文)
        rules_data: 规则数据(如未提供则自动加载)

    Returns:
        特征碎片规则列表
    """
    if rules_data is None:
        rules_data = load_compound_class_rules()

    if not rules_data:
        return []

    fragments = []
    for cls in rules_data.get("classes", []):
        # Match by Chinese or English name
        if class_name.lower() in cls.get("class", "").lower() or \
           class_name.lower() in cls.get("class_en", "").lower():

            # Extract from subclasses
            for sub in cls.get("subclasses", []):
                for rule in sub.get("rules", []):
                    if rule.get("fragment_mz"):
                        fragments.append({
                            "class": cls["class"],
                            "subclass": sub.get("subclass", ""),
                            "fragment_mz": rule["fragment_mz"],
                            "description": rule.get("description", ""),
                            "mechanism": rule.get("mechanism", ""),
                            "neutral_loss": rule.get("neutral_loss"),
                            "confidence": rule.get("confidence", "中"),
                        })

            # Extract from direct rules
            for rule in cls.get("rules", []):
                if rule.get("fragment_mz"):
                    fragments.append({
                        "class": cls["class"],
                        "subclass": "general",
                        "fragment_mz": rule["fragment_mz"],
                        "description": rule.get("description", ""),
                        "mechanism": rule.get("mechanism", ""),
                        "neutral_loss": rule.get("neutral_loss"),
                        "confidence": rule.get("confidence", "中"),
                    })

    return fragments


def get_compound_class_neutral_losses(class_name: str, rules_data: Dict = None) -> List[Dict]:
    """
    获取指定化合物类别的特征中性丢失规则

    Args:
        class_name: 化合物类别名
        rules_data: 规则数据

    Returns:
        特征中性丢失规则列表
    """
    if rules_data is None:
        rules_data = load_compound_class_rules()

    if not rules_data:
        return []

    losses = []
    for cls in rules_data.get("classes", []):
        if class_name.lower() in cls.get("class", "").lower() or \
           class_name.lower() in cls.get("class_en", "").lower():

            for sub in cls.get("subclasses", []):
                for rule in sub.get("rules", []):
                    if rule.get("neutral_loss"):
                        losses.append({
                            "class": cls["class"],
                            "subclass": sub.get("subclass", ""),
                            "neutral_loss": rule["neutral_loss"],
                            "description": rule.get("description", ""),
                            "mechanism": rule.get("mechanism", ""),
                        })

            for rule in cls.get("rules", []):
                if rule.get("neutral_loss"):
                    losses.append({
                        "class": cls["class"],
                        "subclass": "general",
                        "neutral_loss": rule["neutral_loss"],
                        "description": rule.get("description", ""),
                        "mechanism": rule.get("mechanism", ""),
                    })

    return losses


def get_sugar_loss_reference(rules_data: Dict = None) -> Dict:
    """获取糖基中性丢失精确质量参考表"""
    if rules_data is None:
        rules_data = load_compound_class_rules()
    return rules_data.get("sugar_loss_reference", {})


def get_aglycone_reference(rules_data: Dict = None) -> Dict:
    """获取常见苷元质谱参考表"""
    if rules_data is None:
        rules_data = load_compound_class_rules()
    return rules_data.get("common_aglycone_reference", {})


def export_knowledge_base(filepath: str) -> None:
    """导出知识库为JSON"""
    data = {
        "neutral_losses": {
            k: {"formula": v.formula, "mass": v.mass, "name": v.name,
                "name_en": v.name_en, "mechanism": v.mechanism,
                "functional_groups": sorted(v.functional_groups)}
            for k, v in COMMON_NEUTRAL_LOSSES.items()
        },
        "fragmentation_rules": {
            k: {"name": v.name, "name_en": v.name_en,
                "mechanism_type": v.mechanism_type,
                "description": v.description,
                "functional_groups": sorted(v.functional_groups),
                "confidence_score": v.confidence_score}
            for k, v in FRAGMENTATION_RULES.items()
        },
        "diagnostic_ions": {
            k: {"mz": v.mz, "formula": v.formula, "name": v.name,
                "functional_group": v.functional_group,
                "structural_origin": v.structural_origin}
            for k, v in DIAGNOSTIC_IONS.items()
        },
        "bond_dissociation_energies": {
            k: {"bond_type": v.bond_type, "energy_kj_per_mol": v.energy_kj_per_mol,
                "fragmentation_likelihood": v.fragmentation_likelihood}
            for k, v in BOND_DISSOCIATION_ENERGIES.items()
        },
        "summary": get_rules_summary(),
    }
    with open(filepath, "w", encoding="utf-8") as f:
        json.dump(data, f, ensure_ascii=False, indent=2)
