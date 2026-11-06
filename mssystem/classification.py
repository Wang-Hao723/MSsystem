from __future__ import annotations

import re
import warnings
from collections import Counter
from typing import Any, Dict, Iterable, List, Optional, Sequence

import numpy as np
from matchms import Spectrum
from matchms.similarity import CosineGreedy
try:
    from matchms.similarity import ModifiedCosineGreedy as ModifiedCosine
except ImportError:
    from matchms.similarity import ModifiedCosine
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import AllChem
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report
from sklearn.model_selection import GroupShuffleSplit
from sklearn.neural_network import MLPClassifier


def _is_valid_smiles(smiles: str) -> bool:
    """Check if a string looks like a real SMILES (not 'N/A', empty, etc.)."""
    if not smiles or not isinstance(smiles, str):
        return False
    s = smiles.strip()
    if not s:
        return False
    # Exclude common non-SMILES placeholder strings
    if s.lower() in ("n/a", "na", "none", "null", "-", ""):
        return False
    return True


def smiles_to_morgan(smiles: str, radius: int = 2, n_bits: int = 2048) -> np.ndarray:
    if not _is_valid_smiles(smiles):
        return np.zeros(n_bits, dtype=np.float32)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return np.zeros(n_bits, dtype=np.float32)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(
            mol, radius=radius, nBits=n_bits
        )
    # Use direct ToBinary approach compatible with all RDKit versions
    arr = np.frombuffer(fp.ToBinary(), dtype=np.uint8).astype(np.float32)
    if arr.size < n_bits:
        arr = np.pad(arr, (0, n_bits - arr.size), constant_values=0)
    elif arr.size > n_bits:
        arr = arr[:n_bits]
    return arr


class RuleBasedClassifier:
    """基于 SMARTS 子结构匹配的规则分类器 (v3.4.0)

    数据库 classification 字段 ~88% 为来源文件名 (MoNA-export-.../MSMS-.../YDF 等)
    → normalize 后为 Unknown。当 SMILES 可用时, 本分类器作为级联第一层 (Layer 1),
    用 *经过验证的可靠子结构* 主动判类, 避免把判类交给易错的碎片相似度 (Layer 2)。

    规则按 "高特异/高置信 → 宽泛" 排序, 首个命中即返回:
      - alpha-氨基酸/肽 (N-Cα-COOH) 先于含氮环, 避免脯氨酸/组氨酸等被当生物碱
      - 糖环 (吡喃/呋喃) 先于生物碱, 避免氨基糖误判
      - 黄酮 (色原酮/黄烷酮/查耳酮) 覆盖 4 大亚型
      - 3,4-二氢异香豆素按其稠合内酯骨架归入香豆素大类
      - 生物碱: 先匹配明确骨架 (吲哚/嘌呤/四氢异喹啉), 最后用 "环内碱性氮 [NX3;R]"
        兜底吗啡烷/哌啶/莨菪烷等 (如青藤碱 Sinomenine)
    仅保留经验证 *正例命中且负例不误伤* 的 SMARTS；甾体/萜类仍无可靠通用
    SMARTS。二芳基庚烷只使用明确的 Ar-C7-Ar（含羰基）窄规则，避免把普通
    单芳环苯丙素或宽泛多酚误判为该类。
    """

    def __init__(self, smarts_rules: Optional[Dict[str, str]] = None) -> None:
        self.smarts_rules = smarts_rules or {
            # ── 0. 核苷/核苷酸: 含氮碱基N直接连接呋喃糖环 ──
            # 必须先于 purine→alkaloid；否则腺苷/鸟苷会被误标为生物碱。
            "[n;R]-[C;R]1[O;R][C;R][C;R][C;R]1": "nucleoside",
            # ── 1. 氨基酸 / 肽 (alpha-氨基羧基; N-Cα-COOH) ──
            "N[CX4][CX3](=O)O": "amino_acid",
            # ── 2. 香豆素 (苯并吡喃-2-酮) — 先于黄酮, 避免黄酮模式误抓香豆素内酯 ──
            "c1ccc2oc(=O)ccc2c1": "coumarin",
            # ── 2b. 3,4-二氢异香豆素 (benzofused six-membered lactone) ──
            # 区别于五元内酯 phthalide，也不能被羟基苯甲酸规则误抓。
            "O=C1O[C;R][C;R]c2ccccc12": "coumarin",
            # ── 3. 山酮 (xanthone, 二苯并-γ-吡喃酮) — 必须先于黄酮! ── v3 新增 (已验证)
            # 原因: 黄酮 chromone SMARTS(O=c1ccoc2ccccc12) 是 xanthone 的子结构, 若黄酮在前会把
            # 山酮(龙胆/芒果苷元)误判为黄酮。故 xanthone 优先。
            "O=c1c2ccccc2oc2ccccc12": "xanthone",
            # ── 3b. 二苯乙烯 (stilbene: 白藜芦醇/虎杖苷) ── v3 新增 (已验证, 容忍糖基化)
            "c1ccc(C=Cc2ccccc2)cc1": "stilbene",
            # ── 4. 黄酮 (色原酮/黄烷酮/查耳酮) — 不含宽松变体(会误判香豆素) ──
            "O=c1ccoc2ccccc12": "flavonoid",                # chromone: 黄酮/黄酮醇/异黄酮
            "O=C1CC(c2ccccc2)Oc2ccccc21": "flavonoid",       # flavanone (黄烷酮)
            "O=C(/C=C/c1ccccc1)/C=C/c1ccccc1": "flavonoid",  # chalcone (查耳酮)
            # ── 5. 醌 (9,10-蒽醌: 大黄素/芦荟大黄素等) ──
            "O=C1c2ccccc2C(=O)c2ccccc21": "quinone",
            # ── 5b. 二芳基庚烷 (diarylheptanoid: Ar-C7-Ar) ──
            # 七个非环碳连接两个芳环，第三个链碳为羰基；覆盖 Yakuchinone A/B
            # 及相同核心的 curcuminoid，且不会命中普通 C6-C3 苯丙素。
            "c1ccccc1~[C;!R]~[C;!R]~[CX3;!R](=O)~[C;!R]~[C;!R]~[C;!R]~[C;!R]~c2ccccc2": "diarylheptanoid",
            # ── 4b. 酚酸 (羟基苯甲酸: 没食子酸/原儿茶酸/水杨酸/香草酸) ──
            "O=C([$([OX2H1]),$([OX1-])])c1ccccc1O": "phenolic_acid",
            # ── 4c. 酚酸 (羟基肉桂酸: 咖啡酸/阿魏酸/芥子酸/p-香豆酸) ──
            "O=C([$([OX2H1]),$([OX1-])])/C=C/c1ccccc1": "phenolic_acid",
            # 注: 脂肪酸无可靠结构 SMARTS (长链+COOH 会误判三萜酸/胆汁酸), 由标签归一化处理
            # ── 5. 生物碱 — 明确骨架 (优先于宽泛 ring-N) ──
            "[nH]1cccc2ccccc12": "alkaloid",   # 吲哚
            "c1ncnc2c1ncn2": "alkaloid",       # 嘌呤
            "n1cccc2ccccc12": "alkaloid",      # 1,2,3,4-四氢异喹啉
            # ── 7. 生物碱 — 宽泛兜底: 环内碱性氮 (吗啡烷/哌啶/吡咯烷/莨菪烷/嘌呤等) ──
            "[NX3;R]": "alkaloid",
            # ── 8. 糖 (吡喃糖/呋喃糖环) — 置于苷元类之后, 避免黄酮苷/香豆素苷等被误判为糖 ──
            "C1OC(O)C(O)C(O)C1O": "sugar",
            "C1OC(O)C(O)C1O": "sugar",
            # ── 8a. 三级胺氧化物 ──
            # N 必须为四配位正电、与 O− 单键相连且另外三个取代基均为碳。
            # 该约束可覆盖烷基二甲基胺氧化物/TMAO，同时不会把硝基、羟胺、
            # 普通季铵盐或酰胺误归为胺氧化物。
            "[N+X4]([#6])([#6])([#6])[O-]": "amine_oxide",
            # ── 8b. 有机磺酸 — 要求 C-S(=O)2-OH/O−，不把硫酸酯或磺酰胺混入 ──
            # 置于黄酮/香豆素/糖等骨架之后，避免覆盖更具结构特异性的主类别。
            "[#6][SX4](=[OX1])(=[OX1])[$([OX2H1]),$([OX1-])]": "sulfonate",
            # ── 9. 低置信兜底 (苯丙素, 命中即用, 否则交碎片库) ──
            "Oc1ccc(cc1)C=C": "lignan",        # 苯丙素 (C6-C3 侧链)
            "c1cc2c(cc1)OCO2": "lignan",       # 亚甲二氧基苯丙素
        }
        self._compiled = {
            smarts: Chem.MolFromSmarts(smarts) for smarts in self.smarts_rules.keys()
        }
        # A narrow glycerophospholipid rule is kept separate from the generic
        # first-match table because it requires two independently necessary
        # motifs plus a minimum carbon count.  Phosphate alone would mislabel
        # nucleotides/sugar phosphates, while an acyl ester alone would
        # mislabel many natural products.
        self._glycerophosphate_core = Chem.MolFromSmarts(
            "[CH2X4][CHX4]([OX2,OX1-])[CH2X4][OX2]"
            "[PX4](=[OX1])([OX2,OX1-])[OX2,OX1-]"
        )
        self._acyl_glycerol_ester = Chem.MolFromSmarts(
            "[CX3](=[OX1])[OX2][CH2X4]"
        )

    def _looks_like_glycerophospholipid(self, mol) -> bool:
        """Recognize a conservative acyl-glycerophosphate lipid core.

        The conjunction deliberately covers lysophosphatidic acids and common
        glycerophospholipids while rejecting glycerol phosphate, nucleotides,
        phosphorylated sugars, ordinary fatty acids, and neutral glycerides.
        """
        try:
            carbon_count = sum(
                atom.GetAtomicNum() == 6 for atom in mol.GetAtoms()
            )
            return bool(
                carbon_count >= 10
                and self._glycerophosphate_core is not None
                and self._acyl_glycerol_ester is not None
                and mol.HasSubstructMatch(self._glycerophosphate_core)
                and mol.HasSubstructMatch(self._acyl_glycerol_ester)
            )
        except Exception:
            return False

    @staticmethod
    def _looks_like_glycoside(mol) -> bool:
        """吡喃/呋喃糖匹配后, 判断分子是否实为 *糖苷* (含显著苷元) 而非游离糖/寡糖。

        判据: 含芳环 (黄酮/香豆素/酚/蒽醌/香豆素苷) 或 大分子 (皂苷/萜苷, 重原子>28)。
        游离单糖/双糖 (重原子≤24, 无芳环) 不触发 → 仍归 sugar。
        """
        try:
            n_heavy = mol.GetNumHeavyAtoms()
            has_aromatic = any(a.GetIsAromatic() for a in mol.GetAtoms())
            return has_aromatic or n_heavy > 28
        except Exception:
            return False

    @staticmethod
    def _looks_like_acyclic_monosaccharide_or_alditol(mol) -> bool:
        """Conservatively recognize validated open-chain sugars/polyols.

        The ring SMARTS above cannot match the open-chain records commonly
        stored for glucose, mannose, sorbose and sorbitol.  Restrict this
        fallback to small acyclic C/O-only structures with many hydroxyls and
        at most one non-acid carbonyl so polycarboxylic acids and lipids do not
        become carbohydrates.
        """

        try:
            atoms = list(mol.GetAtoms())
            if any(atom.GetIsAromatic() or atom.IsInRing() for atom in atoms):
                return False
            if any(atom.GetAtomicNum() not in {6, 8} for atom in atoms):
                return False
            carbon_count = sum(atom.GetAtomicNum() == 6 for atom in atoms)
            oxygen_count = sum(atom.GetAtomicNum() == 8 for atom in atoms)
            if not (3 <= carbon_count <= 8 and oxygen_count >= carbon_count - 1):
                return False
            acid = Chem.MolFromSmarts("[CX3](=[OX1])[OX2H1,OX1-]")
            if acid is not None and mol.HasSubstructMatch(acid):
                return False
            carbonyl = Chem.MolFromSmarts("[CX3]=[OX1]")
            carbonyl_count = (
                len(mol.GetSubstructMatches(carbonyl))
                if carbonyl is not None else 0
            )
            if carbonyl_count > 1:
                return False
            hydroxyl = Chem.MolFromSmarts("[OX2H1,OX1-]")
            hydroxyl_count = (
                len(mol.GetSubstructMatches(hydroxyl))
                if hydroxyl is not None else 0
            )
            return hydroxyl_count >= max(3, carbon_count - 2)
        except Exception:
            return False

    def predict(self, smiles: Optional[str] = None, inchi: Optional[str] = None) -> str:
        if smiles and not _is_valid_smiles(smiles):
            smiles = None
        mol = None
        if smiles:
            mol = Chem.MolFromSmiles(smiles)
        elif inchi:
            mol = Chem.MolFromInchi(inchi)

        if mol is None:
            return "Unknown"

        if self._looks_like_glycerophospholipid(mol):
            return "lipid"

        if self._looks_like_acyclic_monosaccharide_or_alditol(mol):
            return "sugar"

        for smarts, patt in self._compiled.items():
            if patt is not None and mol.HasSubstructMatch(patt):
                cls = self.smarts_rules[smarts]
                # 糖 vs 糖苷: 吡喃/呋喃糖匹配但分子含苷元 → 归 glycoside (黄酮苷/皂苷/萜苷等)
                if cls == "sugar" and self._looks_like_glycoside(mol):
                    return "glycoside"
                return cls
        return "Unknown"

    def predict_detailed(self, smiles: Optional[str] = None, inchi: Optional[str] = None) -> dict:
        """v3 新增: 返回 {class, source, confidence, matched_smarts} 供自动分类 cascade 用。

        命中 SMARTS → source='smarts', confidence=1.0; 否则 class='Unknown' (交给 Tier2 ML)。
        """
        if smiles and not _is_valid_smiles(smiles):
            smiles = None
        mol = Chem.MolFromSmiles(smiles) if smiles else (Chem.MolFromInchi(inchi) if inchi else None)
        if mol is None:
            return {"class": "Unknown", "source": "none", "confidence": 0.0, "matched_smarts": None}
        if self._looks_like_glycerophospholipid(mol):
            return {
                "class": "lipid",
                "source": "smarts",
                "confidence": 1.0,
                "matched_smarts": "glycerophosphate_core+acyl_glycerol_ester+C>=10",
            }
        if self._looks_like_acyclic_monosaccharide_or_alditol(mol):
            return {
                "class": "sugar",
                "source": "smarts",
                "confidence": 1.0,
                "matched_smarts": (
                    "acyclic_C3-C8_C/O_only_polyhydroxy_"
                    "max_one_nonacid_carbonyl"
                ),
            }
        for smarts, patt in self._compiled.items():
            if patt is not None and mol.HasSubstructMatch(patt):
                cls = self.smarts_rules[smarts]
                if cls == "sugar" and self._looks_like_glycoside(mol):
                    return {"class": "glycoside", "source": "smarts", "confidence": 1.0,
                            "matched_smarts": smarts + "(→glycoside)"}
                return {"class": cls, "source": "smarts", "confidence": 1.0, "matched_smarts": smarts}
        return {"class": "Unknown", "source": "none", "confidence": 0.0, "matched_smarts": None}


# ────────────────────────────────────────────────────────────────────
# 数据集来源名识别 + 文件名归一化 + 级联判类
# 这些工具用于: (1) 库迁移时把被污染的 classification(86%是数据集文件名)清洗为真化合物类;
#              (2) 推理时对仍为 generic/other 的命中走 dl_classification 兜底。
# ────────────────────────────────────────────────────────────────────

# 已知"数据集/来源名"token —— 出现在 classification 里即判为非化合物类
_DATASET_NAME_TOKENS = (
    "MoNA", "MassBank", "GNPS", "Spectraverse", "MSMS", "RIKEN", "MS-DIAL", "MSDIAL",
    "HMDB", "FooDB", "MiMeDB", "VS19", "PFAS", "FiehnHILIC", "Vaniya", "bmdms",
    "experimentspectra", "LC_MS", "QTOF", "All_Spectra", "All_LC", "export",
    "mcescaf", "Metabolite_IDX", "NPL", "VF-NPL", "YDF", "FAHFA", "Dipeptide",
)
# 引用/版本号模式 —— 文献引用串或 _YYYYMMDD 版本后缀
_CITATION_RE = re.compile(r"(et al|Sci\.|Anal\.|Nat\.|Chem\.|Adv\.|Data|\b19\d{2}\b|\b20\d{2}\b)")
# 非类的占位词
_NONCLASS_LITERALS = {
    "other", "n/a", "na", "unknown", "none", "class1", "natural", "oxidized",
    "drug", "industrial_process", "pharmaceutical", "nutrimetabolomics", "",
}
# 把来源文件名/数据集名/原始 classification 归一化成 *规范化合物类* 的映射。
# 规范类与 MS2Classifier 14 类对齐。键含常见同义词/亚类/中文, 子串匹配。
_FILENAME_CLASS_MAP = {
    # 黄酮 (含亚类)
    "Flavonoid": "flavonoid", "Flavonoids": "flavonoid", "Polyphenol": "flavonoid",
    "Flavone": "flavonoid", "Flavonol": "flavonoid", "Flavonols": "flavonoid",
    "Anthocyanidin": "flavonoid", "Anthocyanins": "flavonoid", "Anthocyanin": "flavonoid",
    "Biflavonoid": "flavonoid", "Isoflavone": "flavonoid", "Chalcone": "flavonoid",
    "黄酮": "flavonoid",
    # 生物碱 (含亚类骨架)
    "Alkaloid": "alkaloid", "Alkaloids": "alkaloid",
    "Aporphine": "alkaloid", "Aporphines": "alkaloid",
    "Benzylisoquinoline": "alkaloid", "Protoberberine": "alkaloid",
    "Indole": "alkaloid", "Piperidine": "alkaloid", "Quinolizidine": "alkaloid",
    "生物碱": "alkaloid",
    # 萜类 (含亚类)
    "Terpenoid": "terpenoid", "Terpenoids": "terpenoid", "Terpene": "terpenoid",
    "Triterpene": "terpenoid", "Triterpenoid": "terpenoid", "Triterpenoids": "terpenoid",
    "Sesquiterpene": "terpenoid", "Diterpenoid": "terpenoid", "Monoterpenoid": "terpenoid",
    "Monoterpene": "terpenoid", "Iridoid": "terpenoid", "Iridoids": "terpenoid",
    "Carotenoid": "terpenoid", "Sesterterpene": "terpenoid",
    "萜": "terpenoid", "三萜": "terpenoid", "二萜": "terpenoid", "单萜": "terpenoid", "倍半萜": "terpenoid",
    # 甾体
    "Steroid": "steroid", "Steroids": "steroid", "Sterol": "steroid", "Sterols": "steroid",
    "Cardenolide": "steroid", "Saponin": "saponin", "Saponins": "saponin",
    "皂苷": "saponin", "三萜皂苷": "saponin", "三萜苷": "saponin",
    # 苯丙素 / 香豆素 / 木脂素
    "Coumarin": "coumarin", "Coumarins": "coumarin",
    "Lignan": "lignan", "Lignans": "lignan",
    "Diarylheptanoid": "diarylheptanoid", "Diarylheptanoids": "diarylheptanoid",
    "Phenylpropanoid": "phenolic_acid", "Phenylpropanoids": "phenolic_acid",
    "Phenylethanoid": "phenolic_acid",
    "香豆素": "coumarin", "木脂素": "lignan", "苯丙素": "phenolic_acid",
    # 醌
    "Quinone": "quinone", "Quinones": "quinone", "Anthraquinone": "quinone",
    "Naphthoquinone": "quinone", "Benzoquinone": "quinone",
    "醌": "quinone", "蒽醌": "quinone",
    # 结构上可独立于普通黄酮/生物碱报告的窄类
    "Xanthone": "xanthone", "Xanthones": "xanthone",
    "Stilbene": "stilbene", "Stilbenes": "stilbene", "Stilbenoid": "stilbene",
    "Xanthine": "xanthine", "Methylxanthine": "xanthine",
    # 酚酸
    "Phenolic": "phenolic_acid", "Phenol": "phenolic_acid", "OrganicAcid": "phenolic_acid",
    "Tannin": "phenolic_acid", "Tannins": "phenolic_acid",
    "酚酸": "phenolic_acid",
    # 氨基酸 / 肽 / 核苷
    "AminoAcid": "amino_acid", "Amino_Acid": "amino_acid", "Peptide": "amino_acid",
    "Peptides": "amino_acid", "Amino": "amino_acid",
    "Nucleoside": "nucleoside", "Nucleosides": "nucleoside",
    "Nucleotide": "nucleoside", "Nucleotides": "nucleoside",
    # 糖
    "Sugar": "sugar", "Carbohydrate": "sugar", "Carbohydrates": "sugar",
    "Glycoside": "glycoside", "Glycosides": "glycoside",
    "糖": "sugar", "糖苷": "glycoside",
    # 脂类 / 脂肪酸
    "FattyAcid": "fatty_acyl", "FattyAcyl": "fatty_acyl", "Fatty_Acid": "fatty_acyl",
    "Fatty": "fatty_acyl",
    "Lipid": "lipid", "FAHFA": "lipid", "Glycerophospholipid": "lipid",
    "Sphingolipid": "lipid", "Ceramide": "lipid", "Phospholipid": "lipid",
    "脂": "lipid", "脂肪酸": "fatty_acyl",
    # 三级胺氧化物 / 氧化胺型表面活性剂
    "AmineOxide": "amine_oxide", "Amine_Oxide": "amine_oxide",
    "Amine oxide": "amine_oxide", "胺氧化物": "amine_oxide",
    # 有机磺酸
    "Sulfonate": "sulfonate", "SulfonicAcid": "sulfonate", "Sulfonic_Acid": "sulfonate",
    "磺酸": "sulfonate",
}
# 规范类白名单 (与 MS2Classifier.CLASS_NAMES 对齐 + glycoside)
CANONICAL_CLASSES = {
    "flavonoid", "alkaloid", "terpenoid", "steroid", "saponin",
    "phenolic_acid", "coumarin", "lignan", "quinone",
    "sugar", "glycoside", "fatty_acyl", "lipid", "amino_acid", "other",
    "nucleoside", "xanthone", "stilbene", "xanthine",
    "diarylheptanoid", "sulfonate", "amine_oxide",
}


def is_source_or_generic(cls: Optional[str]) -> bool:
    """判断一个 classification 值是否是 *数据集/来源名或占位词* 而非真化合物类。
    库构建时 86% 的 classification 被写成了文件名 (MoNA-export-All_Spectra 等),
    迁移和推理时必须把这些识别出来并丢弃, 改用结构/ML 判类。
    """
    if cls is None:
        return True
    s = str(cls).strip()
    if not s:
        return True
    low = s.lower()
    if low in _NONCLASS_LITERALS:
        return True
    # 含任一已知数据集 token
    if any(tok.lower() in low for tok in _DATASET_NAME_TOKENS):
        return True
    # 文献引用 / 版本号
    if _CITATION_RE.search(s):
        return True
    # 文件名后缀
    if re.search(r"_20\d{6}|\.ms[pg]|\.mgf", low):
        return True
    return False


def normalize_filename_class(stem: Optional[str]) -> Optional[str]:
    """把来源文件名/目录名归一化成化合物类。返回 None 表示文件名不含类信息。
    例: 'FAHFA' → 'lipid', 'Alkaloids_QE_pos' → 'alkaloid', 'Terpenoid' → 'terpenoid'。
    """
    if not stem:
        return None
    s = str(stem).strip()
    # 取 basename (兼容 'Spectraverse/GNPS' 这类伪路径)
    if "/" in s or "\\" in s:
        s = re.split(r"[\\/]", s)[-1]
    # 去扩展名
    s = re.sub(r"\.(mgf|msp|mzml|mzxml)$", "", s, flags=re.I)
    if s in _FILENAME_CLASS_MAP:
        return _FILENAME_CLASS_MAP[s]
    low = s.lower()
    for key, val in _FILENAME_CLASS_MAP.items():
        if key.lower() in low:
            return val
    return None


def normalize_source_file(raw: Optional[str]) -> str:
    """把 Spectra.source_file 原始值归一化为 '数据集名[-极性]' token (dataset_source 列用)。

    规则:
      - 'Spectraverse/<X>' → 原样保留 (已是干净的 数据集/子库 形态)
      - 真实路径 → 取 basename 去扩展名, 抽极性 (-Pos-/-Neg-/_pos_/_neg_/-pos-/-neg-/POS_/NEG_),
        去版本日期 (_YYYYMMDD / YYYYMMDD_) 与冗余后缀 (_all_lib_MS2/_library/...), 拼成 '数据集[-极性]'
      - 极性前缀 POS_/NEG_ 也识别
    例:
      '.../MSP_data/MSMS-Pos-GNPS.msp'         → 'MSMS-GNPS-Pos'
      '.../MSP_data/MoNA-export-All_Spectra.msp' → 'MoNA-export-All_Spectra'
      'Spectraverse/GNPS'                       → 'Spectraverse/GNPS'
      '.../MGF_other_data/20231130_mcescaf_library_pos_all_lib_MS2.mgf' → 'mcescaf-Pos'
    """
    if not raw:
        return ""
    s = str(raw).strip()
    if not s:
        return ""
    # 1. Spectraverse 伪路径原样保留
    if s.startswith("Spectraverse/"):
        return s
    # 2. 取 basename + 去扩展名
    base = re.split(r"[\\/]", s)[-1]
    stem = re.sub(r"\.(mgf|msp|mzml|mzxml|mzxml)$", "", base, flags=re.I)
    # 3. 抽极性 (各种写法)
    polarity = None
    m = re.search(r"(?:^|[-_])(pos|positive|neg|negative)(?:[-_]|$)", stem, re.I)
    if m:
        polarity = "Pos" if m.group(1).lower().startswith("pos") else "Neg"
    # 4. 去极性段 + 去日期 + 去冗余后缀
    dataset = stem
    if polarity:
        # 删掉极性 token (含其前导分隔符)
        dataset = re.sub(
            r"(^|[-_])(pos|positive|neg|negative)([-_]|$)", r"\1", dataset,
            flags=re.I,
        )
    dataset = re.sub(r"_?\d{8}_?", "", dataset)           # 20200806 / 20231130_
    dataset = re.sub(
        r"(_all_lib_MS2|_all_lib|_library|_LC_MS_Spectra|_LC_MS|export_All_Spectra|All_Spectra)$",
        "", dataset, flags=re.I,
    )
    dataset = re.sub(r"[_\-]{2,}", "-", dataset)           # 合并连续分隔符
    dataset = dataset.strip("-_")
    if not dataset:
        dataset = stem
    # 5. 拼 token
    token = dataset
    if polarity:
        token = f"{dataset}-{polarity}"
    return token


# 模块级单例, 避免反复编译 SMARTS
_RULE_CLF: Optional[RuleBasedClassifier] = None


def _get_rule_clf() -> RuleBasedClassifier:
    global _RULE_CLF
    if _RULE_CLF is None:
        _RULE_CLF = RuleBasedClassifier()
    return _RULE_CLF


def resolve_compound_class(
    smiles: Optional[str] = None,
    raw: Optional[str] = None,
    source_file: Optional[str] = None,
) -> str:
    """级联判类, 返回统一的化合物种类标签 (黄酮/生物碱/萜类/…/other)。

    优先级:
      1. SMARTS 结构规则 (RuleBasedClassifier.predict, 确定性高精度, ~25-30% 覆盖)
      2. 来源文件名归一化 (normalize_filename_class, MGF_classify_data 目录文件按类命名)
      3. 原 classification 值 (若不是数据集名/占位词)
      4. 'other' (留待推理时 dl_classification 兜底)

    Args:
        smiles: 化合物 SMILES (结构判类)
        raw: 原 classification 字段值 (可能是数据集名, 需过滤)
        source_file: 来源文件路径/名 (用于文件名归一化)
    """
    # 1. 结构判类
    if smiles:
        cls = _get_rule_clf().predict(smiles=smiles)
        if cls and cls != "Unknown":
            return cls
    # 2. 文件名归一化
    nf = normalize_filename_class(source_file)
    if nf:
        return nf
    # 3. 原 raw 值归一化为规范类 (拒绝数据集名/药材名/占位词)
    if raw:
        mapped = normalize_filename_class(raw)
        if mapped:
            return mapped
    # 4. 兜底
    return "other"



class FingerprintClassifier:
    def __init__(self, model: Any | None = None) -> None:
        self.model = model

    def predict(self, smiles: str) -> str:
        if self.model is None:
            return "Unknown"
        if not _is_valid_smiles(smiles):
            return "Unknown"
        feat = smiles_to_morgan(smiles)
        pred = self.model.predict([feat])[0]
        return str(pred)


def train_classifier(
    library_data,
    model_type: str = "rf",
    test_size: float = 0.2,
    random_state: int = 42,
    min_samples_per_class: int = 2,
) -> Dict[str, Any]:
    if hasattr(library_data, "to_dict"):
        records = library_data.to_dict(orient="records")
    else:
        records = list(library_data)

    # Filter: must have valid SMILES, a classification label that is not "Unknown"
    valid = [
        r
        for r in records
        if _is_valid_smiles(r.get("smiles"))
        and r.get("classification")
        and str(r.get("classification")).strip()
        and r.get("classification") != "Unknown"
    ]

    if len(valid) < 20:
        raise ValueError(
            f"Not enough labeled samples to train classifier "
            f"(got {len(valid)}, need >=20)."
        )

    # Filter out classes with too few samples for stratified split
    class_counts = Counter(r["classification"] for r in valid)
    valid = [
        r for r in valid
        if class_counts[r["classification"]] >= min_samples_per_class
    ]

    if len(valid) < 20:
        raise ValueError(
            f"After filtering rare classes (<{min_samples_per_class} samples), "
            f"only {len(valid)} samples remain — cannot train."
        )

    # Warn if any classes were removed
    remaining_classes = Counter(r["classification"] for r in valid)
    if len(remaining_classes) < len(class_counts):
        dropped = set(class_counts.keys()) - set(remaining_classes.keys())
        import logging
        logging.getLogger(__name__).warning(
            "Dropped %d rare class(es) from training: %s",
            len(dropped), sorted(dropped)
        )

    # Canonicalize/deduplicate before splitting.  Conflicting duplicate labels
    # are discarded rather than allowing the same ECFP vector to teach and
    # test different targets.
    deduplicated: Dict[str, dict] = {}
    conflicts = set()
    for record in valid:
        molecule = Chem.MolFromSmiles(str(record["smiles"]))
        if molecule is None:
            continue
        canonical = Chem.MolToSmiles(
            molecule, canonical=True, isomericSmiles=True
        )
        existing = deduplicated.get(canonical)
        if existing and existing["classification"] != record["classification"]:
            conflicts.add(canonical)
            continue
        copied = dict(record)
        copied["smiles"] = canonical
        deduplicated[canonical] = copied
    valid = [
        record for canonical, record in deduplicated.items()
        if canonical not in conflicts
    ]
    if len(valid) < 20:
        raise ValueError("Too few unique, non-conflicting structures after deduplication")

    def scaffold_group(smiles: str) -> str:
        molecule = Chem.MolFromSmiles(smiles)
        scaffold = MurckoScaffold.GetScaffoldForMol(molecule)
        if scaffold is not None and scaffold.GetNumAtoms() > 0:
            generic = MurckoScaffold.MakeScaffoldGeneric(scaffold)
            return "murcko:" + Chem.MolToSmiles(generic, canonical=True)
        return "acyclic_formula:" + rdMolDescriptors.CalcMolFormula(molecule)

    X = np.array([smiles_to_morgan(r["smiles"]) for r in valid])
    y = np.array([r["classification"] for r in valid])
    groups = np.array([scaffold_group(r["smiles"]) for r in valid])
    all_classes = set(np.unique(y).tolist())
    split = None
    splitter = GroupShuffleSplit(
        n_splits=20,
        test_size=test_size,
        random_state=random_state,
    )
    for train_index, test_index in splitter.split(X, y, groups):
        if set(np.unique(y[train_index]).tolist()) == all_classes:
            split = (train_index, test_index)
            break
    if split is None:
        raise ValueError(
            "Could not create a scaffold-disjoint split with every class in training"
        )
    train_index, test_index = split
    X_train, X_test = X[train_index], X[test_index]
    y_train, y_test = y[train_index], y[test_index]

    if model_type.lower() == "mlp":
        model = MLPClassifier(
            hidden_layer_sizes=(512, 128),
            max_iter=300,
            random_state=random_state,
        )
    else:
        model = RandomForestClassifier(
            n_estimators=300,
            random_state=random_state,
            n_jobs=-1,
            class_weight="balanced",
        )

    model.fit(X_train, y_train)
    y_pred = model.predict(X_test)

    report = classification_report(y_test, y_pred, output_dict=True)
    return {
        "model": model,
        "report": report,
        "label_distribution": Counter(y.tolist()),
        "split_policy": "canonical_smiles_dedup_bemis_murcko_group_split",
        "split_audit": {
            "train_n": int(len(train_index)),
            "test_n": int(len(test_index)),
            "scaffold_overlap": len(
                set(groups[train_index]) & set(groups[test_index])
            ),
            "conflicting_duplicate_structures_discarded": len(conflicts),
        },
        "production_eligible": False,
    }
