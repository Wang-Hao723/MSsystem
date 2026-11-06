"""
测试分类标准化、中性丢失过滤、连通路径构建、诊断离子过滤
以及全类别分类准确率验证 (v3.3.0 expanded)
"""

import pytest
import sys
from pathlib import Path

# 确保可以导入 matching 模块
sys.path.insert(0, str(Path(__file__).parent.parent))

from matching.fragmentation_rules import (
    normalize_classification,
    get_neutral_loss_by_mass,
    get_composite_neutral_loss_by_mass,
    match_diagnostic_ion,
    _is_nl_compatible,
    _is_diagnostic_ion_compatible,
    COMMON_NEUTRAL_LOSSES,
    DIAGNOSTIC_IONS,
    FUNCTIONAL_GROUP_SIGNATURES,
    _CLASS_TO_FGS,
    NeutralLoss,
    DiagnosticIon,
)
from matching.expert_scorer import (
    ChemicalExpertScorer, enrich_neutral_loss_analysis, expert_evaluate,
)
from matching.fragment_feature_library import (
    classify_compound_from_fragments,
    COMPOUND_CLASS_PROFILES,
    get_profile,
    get_library_summary,
)
from mssystem.classification import RuleBasedClassifier


# ═══════════════════════════════════════════════════════════════
# 1. 分类标准化测试
# ═══════════════════════════════════════════════════════════════

class TestNormalizeClassification:
    """normalize_classification 函数测试"""

    def test_chinese_names(self):
        assert normalize_classification("酚酸类") == "phenolic_acid"
        assert normalize_classification("黄酮类") == "flavonoid"
        assert normalize_classification("生物碱类") == "alkaloid"
        assert normalize_classification("萜类") == "terpenoid"
        assert normalize_classification("香豆素类") == "coumarin"
        assert normalize_classification("醌类") == "quinone"
        assert normalize_classification("甾体类") == "steroid"
        assert normalize_classification("脂肪酸") == "fatty_acid"
        assert normalize_classification("氨基酸") == "amino_acid"
        assert normalize_classification("木脂素类") == "lignan"
        assert normalize_classification("皂苷类") == "saponin"
        assert normalize_classification("糖类") == "sugar"
        assert normalize_classification("单糖") == "sugar"
        assert normalize_classification("苯环类") == "benzenoid"

    def test_spaces_and_underscores(self):
        assert normalize_classification("phenolic acid") == "phenolic_acid"
        assert normalize_classification("Phenolic_Acid") == "phenolic_acid"
        assert normalize_classification("PHENOLIC ACID") == "phenolic_acid"
        assert normalize_classification("fatty acid") == "fatty_acid"
        assert normalize_classification("amino acid") == "amino_acid"

    def test_unknown_values(self):
        assert normalize_classification("") == "Unknown"
        assert normalize_classification("  ") == "Unknown"
        assert normalize_classification("n/a") == "Unknown"
        assert normalize_classification("Unknown") == "Unknown"
        assert normalize_classification("other") == "Unknown"
        assert normalize_classification("random_string_xyz") == "Unknown"

    def test_compound_names(self):
        assert normalize_classification("Steroids/Glycoside") == "steroid"

    def test_plurals(self):
        assert normalize_classification("flavonoids") == "flavonoid"
        assert normalize_classification("alkaloids") == "alkaloid"
        assert normalize_classification("coumarins") == "coumarin"
        assert normalize_classification("quinones") == "quinone"
        assert normalize_classification("steroids") == "steroid"
        assert normalize_classification("terpenoids") == "terpenoid"
        assert normalize_classification("saponins") == "saponin"
        assert normalize_classification("lignans") == "lignan"
        assert normalize_classification("sugars") == "sugar"
        assert normalize_classification("lipids") == "lipid"

    def test_subtypes(self):
        assert normalize_classification("isoflavonoid") == "flavonoid"
        assert normalize_classification("flavonol") == "flavonoid"
        assert normalize_classification("flavanone") == "flavonoid"
        assert normalize_classification("diterpenoid") == "terpenoid"
        assert normalize_classification("sesquiterpenoid") == "terpenoid"
        assert normalize_classification("triterpenoid") == "terpenoid"
        assert normalize_classification("monosaccharide") == "sugar"
        assert normalize_classification("carbohydrate") == "sugar"
        assert normalize_classification("peptide") == "amino_acid"
        assert normalize_classification("purine alkaloids") == "alkaloid"

    def test_case_insensitive(self):
        assert normalize_classification("FLAVONOID") == "flavonoid"
        assert normalize_classification("Alkaloid") == "alkaloid"
        assert normalize_classification("PHENOLIC_ACID") == "phenolic_acid"
        assert normalize_classification("STEROID") == "steroid"
        assert normalize_classification("Lipid") == "lipid"


class TestAmineOxideClassification:
    def test_lauramine_oxide_and_tmao_use_specific_structure_rule(self):
        classifier = RuleBasedClassifier()
        assert classifier.predict(
            "CCCCCCCCCCCC[N+](C)(C)[O-]"
        ) == "amine_oxide"
        assert classifier.predict("C[N+](C)(C)[O-]") == "amine_oxide"

    def test_rule_does_not_capture_nitro_or_quaternary_ammonium(self):
        classifier = RuleBasedClassifier()
        assert classifier.predict("C[N+](=O)[O-]") != "amine_oxide"
        assert classifier.predict("C[N+](C)(C)C") != "amine_oxide"


class TestOpenChainSugarClassification:
    def test_validated_open_chain_hexose_and_alditol_have_structure_rule(self):
        classifier = RuleBasedClassifier()
        assert classifier.predict(
            "OC[C@H](O)[C@@H](O)[C@H](O)[C@H](O)C=O"
        ) == "sugar"
        assert classifier.predict(
            "OC[C@H](O)[C@@H](O)[C@H](O)[C@H](O)CO"
        ) == "sugar"

    def test_polycarboxylic_acid_is_not_mislabeled_as_sugar(self):
        classifier = RuleBasedClassifier()
        assert classifier.predict("O=C(O)CC(O)(CC(=O)O)C(=O)O") != "sugar"


# ═══════════════════════════════════════════════════════════════
# 2. 中性丢失过滤测试
# ═══════════════════════════════════════════════════════════════

class TestNeutralLossFiltering:
    """中性丢失过滤逻辑测试"""

    def test_blocks_acetonitrile_for_phenolic_acid(self):
        """酚酸类不应匹配乙腈"""
        results = get_neutral_loss_by_mass(
            41.0265, 0.5,
            classification="phenolic_acid", strict_fg=True
        )
        names = [nl.name for nl in results]
        assert "乙腈" not in names, f"乙腈不应出现在酚酸类中性丢失中"

    def test_blocks_phosphoric_acid_for_phenolic_acid(self):
        """酚酸类不应匹配磷酸"""
        results = get_neutral_loss_by_mass(
            97.9769, 0.5,
            classification="phenolic_acid", strict_fg=True
        )
        names = [nl.name for nl in results]
        assert "磷酸" not in names

    def test_blocks_ketene_for_unknown(self):
        """未知分类不应匹配乙烯酮"""
        results = get_neutral_loss_by_mass(
            42.0106, 0.5,
            classification="", strict_fg=True
        )
        names = [nl.name for nl in results]
        assert "乙烯酮" not in names

    def test_allows_water_for_all(self):
        """水应始终通过过滤"""
        for cls in ["phenolic_acid", "flavonoid", "alkaloid", "", "Unknown", "sugar", "steroid", "amino_acid"]:
            results = get_neutral_loss_by_mass(
                18.0106, 0.5,
                classification=cls, strict_fg=True
            )
            names = [nl.name for nl in results]
            assert "水" in names, f"水应出现在分类 '{cls}' 的中性丢失中"

    def test_allows_co2_for_phenolic_acid(self):
        """CO2 是酚酸类的最特征丢失"""
        results = get_neutral_loss_by_mass(
            43.9898, 0.5,
            classification="phenolic_acid", strict_fg=True
        )
        names = [nl.name for nl in results]
        assert "二氧化碳" in names

    def test_ammonia_loss_requires_primary_amine_or_primary_amide(self):
        amine_oxide = "CCCCCCCCCCCC[N+](C)(C)[O-]"
        assert not get_neutral_loss_by_mass(
            17.0265,
            0.02,
            classification="amine_oxide",
            candidate_smiles=amine_oxide,
        )
        assert not get_neutral_loss_by_mass(
            230.24877 - 213.24250,
            0.02,
            classification="amine_oxide",
            candidate_smiles=amine_oxide,
        )
        assert any(
            item.formula == "NH3"
            for item in get_neutral_loss_by_mass(
                17.0265,
                0.02,
                classification="alkaloid",
                candidate_smiles="CCCCN",
            )
        )

    def test_allows_acetonitrile_for_alkaloid(self):
        """生物碱类应允许乙腈丢失"""
        results = get_neutral_loss_by_mass(
            41.0265, 0.5,
            classification="alkaloid", strict_fg=True
        )
        names = [nl.name for nl in results]
        assert "乙腈" in names

    def test_allows_sugar_loss_for_saponin(self):
        """皂苷类应允许己糖丢失"""
        results = get_neutral_loss_by_mass(
            162.0528, 0.5,
            classification="saponin", strict_fg=True
        )
        names = [nl.name for nl in results]
        assert "己糖残基" in names

    def test_allows_sugar_loss_for_flavonoid(self):
        """黄酮类应允许己糖丢失(糖苷)"""
        results = get_neutral_loss_by_mass(
            162.0528, 0.5,
            classification="flavonoid", strict_fg=True
        )
        names = [nl.name for nl in results]
        assert "己糖残基" in names

    def test_cumulative_hexose_plus_water_loss(self):
        """180.0634 Da 应识别为脱己糖后再脱水，而不是无规则候选。"""
        results = get_composite_neutral_loss_by_mass(
            180.0634, 0.05,
            classification="glycoside", strict_fg=True,
        )
        assert results
        assert set(results[0]["components"]) == {"己糖残基", "水"}

    def test_arbitrary_difference_is_not_forced_into_composite(self):
        """任意 86 Da 差值不能为消除提示而强行编造累积中性丢失。"""
        assert not get_composite_neutral_loss_by_mass(
            86.0, 0.05,
            classification="phenolic_acid", strict_fg=True,
        )

    def test_legacy_analysis_is_enriched_but_unmatched_stays_unassigned(self):
        analysis = {
            "neutral_losses": [
                {"mass_difference": 180.0634, "neutral_loss": None},
                {"mass_difference": 86.0, "neutral_loss": None},
            ],
            "sequential_explained": 0,
            "cascade_recovered": 0,
        }
        enrich_neutral_loss_analysis(analysis, "glycoside")
        assert analysis["neutral_losses"][0]["neutral_loss"].startswith("累积丢失：")
        assert analysis["neutral_losses"][0]["assignment_level"] == "tentative_composite"
        assert analysis["neutral_losses"][1]["neutral_loss"] is None
        assert analysis["neutral_losses"][1]["assignment_status"] == "unassigned_fragment_difference"


# ═══════════════════════════════════════════════════════════════
# 3. 诊断离子过滤测试
# ═══════════════════════════════════════════════════════════════

class TestDiagnosticIonFiltering:
    """诊断离子过滤逻辑测试"""

    def test_blocks_flavonoid_ion_for_phenolic_acid(self):
        """酚酸类不应匹配黄酮B环碎片"""
        results = match_diagnostic_ion(
            120.0575, 0.05,
            classification="phenolic_acid", strict_fg=True
        )
        flavonoid_matches = [ion for ion in results if ion.functional_group == "flavonoid"]
        assert len(flavonoid_matches) == 0

    def test_allows_flavonoid_ion_for_flavonoid(self):
        """黄酮类应允许黄酮B环碎片"""
        results = match_diagnostic_ion(
            120.0575, 0.05,
            classification="flavonoid", strict_fg=True
        )
        names = [ion.name for ion in results]
        assert any("B环" in n for n in names)

    def test_blocks_class_specific_for_unknown(self):
        """未知分类应阻断类特异诊断离子"""
        results = match_diagnostic_ion(
            120.0575, 0.05,
            classification="", strict_fg=True
        )
        flav_matches = [ion for ion in results if ion.functional_group == "flavonoid"]
        assert len(flav_matches) == 0

    def test_allows_steroid_ion_for_steroid(self):
        """甾体类应允许甾体D环碎片"""
        results = match_diagnostic_ion(
            255.2113, 0.05,
            classification="steroid", strict_fg=True
        )
        names = [ion.name for ion in results]
        assert any("甾体" in n for n in names)

    def test_blocks_steroid_ion_for_sugar(self):
        """糖类不应匹配甾体D环碎片"""
        results = match_diagnostic_ion(
            255.2113, 0.05,
            classification="sugar", strict_fg=True
        )
        steroid_matches = [ion for ion in results if ion.functional_group == "steroid"]
        assert len(steroid_matches) == 0

    def test_allows_terpene_ion_for_terpenoid(self):
        """萜类应允许萜类诊断离子"""
        results = match_diagnostic_ion(
            137.1328, 0.05,
            classification="terpenoid", strict_fg=True
        )
        terp_matches = [ion for ion in results if ion.functional_group == "terpenoid"]
        assert len(terp_matches) > 0


# ═══════════════════════════════════════════════════════════════
# 4. 连通路径构建测试
# ═══════════════════════════════════════════════════════════════

class TestConnectedPathways:
    """连通碎裂链构建测试"""

    def setup_method(self):
        self.scorer = ChemicalExpertScorer()

    def test_connected_chain_basic(self):
        """真实链应正确构建"""
        sorted_frags = [300.0, 282.0, 254.0]
        pathways = self.scorer._build_connected_pathways(
            sorted_frags, classification="phenolic_acid"
        )
        assert len(pathways) > 0
        best = pathways[0]
        assert len(best) >= 2

    def test_no_false_chains(self):
        """无关碎片不应产生伪路径"""
        sorted_frags = [500.0, 499.0, 498.5]
        pathways = self.scorer._build_connected_pathways(
            sorted_frags, classification="phenolic_acid"
        )
        for pathway in pathways:
            for step in pathway:
                assert step["neutral_loss"] not in {"乙腈", "乙烯酮", "磷酸"}

    def test_no_acetonitrile_in_phenolic_acid_pathway(self):
        """酚酸类碎裂路径不应包含乙腈"""
        sorted_frags = [300.0, 259.0, 217.0, 119.0]
        pathways = self.scorer._build_connected_pathways(
            sorted_frags, classification="phenolic_acid"
        )
        for pathway in pathways:
            for step in pathway:
                assert step["neutral_loss"] != "乙腈"
                assert step["neutral_loss"] != "磷酸"

    def test_max_depth_respected(self):
        """路径深度不应超过 max_depth"""
        sorted_frags = [400.0, 382.0, 354.0, 326.0, 298.0, 270.0]
        pathways = self.scorer._build_connected_pathways(
            sorted_frags, classification="flavonoid", max_depth=3
        )
        for pathway in pathways:
            assert len(pathway) <= 3


# ═══════════════════════════════════════════════════════════════
# 5. 完整报告生成测试
# ═══════════════════════════════════════════════════════════════

class TestReportGeneration:
    """报告生成集成测试"""

    def setup_method(self):
        self.scorer = ChemicalExpertScorer()

    def test_no_nonsensical_pathway_in_report(self):
        """酚酸类数据的报告不应包含乙腈/乙烯酮/磷酸"""
        peaks = [
            [163.0390, 500.0],
            [135.0440, 300.0],
            [117.0340, 200.0],
            [89.0390, 100.0],
        ]
        result = self.scorer.evaluate(
            query_peaks=peaks,
            query_pmz=181.05,
            candidate_smiles="",
            candidate_formula="",
            candidate_classification="phenolic_acid",
        )
        reasoning = result.get("chemical_reasoning", "")
        assert "乙腈" not in reasoning
        assert "乙烯酮" not in reasoning
        assert "磷酸" not in reasoning

    def test_effective_classification_in_result(self):
        """evaluate 结果应包含 effective_classification"""
        peaks = [[120.0, 500.0]]
        result = self.scorer.evaluate(
            query_peaks=peaks,
            query_pmz=200.0,
            candidate_smiles="",
            candidate_classification="Phenolic Acid",
        )
        assert result["effective_classification"] == "phenolic_acid"

    def test_chinese_classification_normalized(self):
        """中文分类名应被正确标准化"""
        peaks = [[120.0, 500.0]]
        result = self.scorer.evaluate(
            query_peaks=peaks,
            query_pmz=200.0,
            candidate_smiles="",
            candidate_classification="酚酸类",
        )
        assert result["effective_classification"] == "phenolic_acid"


# ═══════════════════════════════════════════════════════════════
# 6. 糖类分类与裂解测试
# ═══════════════════════════════════════════════════════════════

class TestSugarClassification:
    """测试单糖/糖类不会被误判为生物碱类"""

    def test_sugar_normalize_variants(self):
        """糖类各种名称变体都应标准化为 sugar"""
        sugar_variants = [
            "Sugar", "Monosaccharide", "monosaccharides",
            "Carbohydrate", "Hexose", "Rhamnose", "Glucose",
            "糖类", "单糖", "碳水化合物",
        ]
        for v in sugar_variants:
            assert normalize_classification(v) == "sugar", \
                f"'{v}' should normalize to 'sugar', got '{normalize_classification(v)}'"

    def test_sugar_fragment_classification(self):
        """Rhamnose碎片应被分类为 sugar"""
        rhamnose_fragments = [145.9, 129.0, 111.0, 87.0, 101.1, 73.0]
        result = classify_compound_from_fragments(rhamnose_fragments, tolerance=0.5)
        assert result, "Should return a classification result"
        assert result["class_name"] == "sugar", \
            f"Rhamnose should classify as sugar, got {result['class_name']}"

    def test_sugar_no_ammonia_loss(self):
        """糖类裂解路径中不应出现氨丢失"""
        peaks = [
            [145.9, 100.0],
            [129.0, 60.0],
            [111.0, 30.0],
        ]
        result = expert_evaluate(
            query_peaks=peaks, query_pmz=165.07,
            candidate_smiles="C[C@H](O1)[C@H](O)[C@@H](O)[C@@H](O)C(O)1",
            candidate_classification="sugar",
        )
        nl_analysis = result.get("neutral_loss_analysis", {})
        sequential = nl_analysis.get("sequential_losses", [])
        nl_names = [s.get("neutral_loss", "") for s in sequential]
        assert "氨" not in nl_names, \
            f"Sugar should NOT have 氨 loss, but found: {nl_names}"

    def test_sugar_feature_library_match(self):
        """糖类碎片应匹配糖类特征库"""
        peaks = [
            [145.9, 100.0],
            [129.0, 60.0],
            [111.0, 30.0],
        ]
        result = expert_evaluate(
            query_peaks=peaks, query_pmz=165.07,
            candidate_smiles="C[C@H](O1)[C@H](O)[C@@H](O)[C@@H](O)C(O)1",
            candidate_classification="sugar",
        )
        ff = result.get("feature_library_analysis", {})
        profile = ff.get("class_profile", "")
        assert "糖" in profile, \
            f"Feature library should match 糖类 profile, got '{profile}'"

    def test_sugar_reasoning_no_alkaloid(self):
        """糖类化学推理中不应出现生物碱"""
        peaks = [[145.9, 100.0], [129.0, 60.0], [111.0, 30.0]]
        result = expert_evaluate(
            query_peaks=peaks, query_pmz=165.07,
            candidate_smiles="C[C@H](O1)[C@H](O)[C@@H](O)[C@@H](O)C(O)1",
            candidate_classification="sugar",
        )
        reasoning = result.get("chemical_reasoning", "")
        assert "生物碱" not in reasoning


# ═══════════════════════════════════════════════════════════════
# 7. 全类别分类准确率测试 (v3.3.0 新增)
# ═══════════════════════════════════════════════════════════════

class TestAllCategoryClassification:
    """验证所有化合物类别的分类准确率"""

    # ── 每个类别的典型碎片数据 ──
    TEST_CASES = {
        "flavonoid": {
            "name": "槲皮素 (Quercetin)",
            "smiles": "O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc12",
            "fragments": [153.02, 137.02, 152.05, 121.03, 109.03],
            "pmz": 303.05,
        },
        "alkaloid": {
            "name": "咖啡因 (Caffeine)",
            "smiles": "Cn1c(=O)c2c(ncn2C)n(C)c1=O",
            "fragments": [138.07, 110.06, 83.06, 67.03],
            "pmz": 195.09,
        },
        "phenolic_acid": {
            "name": "咖啡酸 (Caffeic acid)",
            "smiles": "O=C(/C=C/c1ccc(O)c(O)c1)O",
            "fragments": [163.04, 119.05, 197.05, 179.03],
            "pmz": 181.05,
        },
        "coumarin": {
            "name": "香豆素 (Coumarin)",
            "smiles": "O=c1ccc2ccccc2o1",
            "fragments": [147.04, 103.05, 91.05, 119.05],
            "pmz": 147.04,
        },
        "quinone": {
            "name": "丹参酮IIA (Tanshinone IIA)",
            "smiles": "CC1=CC2C3CCC4=CC(=O)C=CC4(C)C3CC2C1=O",
            "fragments": [277.14, 249.15, 235.15, 221.13],
            "pmz": 295.17,
        },
        "lignan": {
            "name": "五味子甲素 (Schisandrin A)",
            "smiles": "COc1ccc2c(c1)CC1OC3c4ccc(OC)c(OC)c4CCC3c2c1",
            "fragments": [385.17, 355.15, 331.15, 315.16],
            "pmz": 385.17,
        },
        "sugar": {
            "name": "鼠李糖 (Rhamnose)",
            "smiles": "C[C@H](O1)[C@H](O)[C@@H](O)[C@@H](O)C(O)1",
            "fragments": [147.07, 129.06, 111.04, 87.04, 73.03],
            "pmz": 165.07,
        },
        "saponin": {
            "name": "人参皂苷Rb1 (Ginsenoside Rb1)",
            "smiles": "",  # 复杂结构，空SMILES
            "fragments": [162.05, 146.06, 132.04, 109.10, 203.18],
            "pmz": 1109.60,
        },
        "terpenoid": {
            "name": "柠檬烯 (Limonene)",
            "smiles": "CC1=CCC(CC1)=C",
            "fragments": [137.13, 121.10, 93.07, 81.07, 68.06],
            "pmz": 137.13,
        },
        "steroid": {
            "name": "胆甾醇 (Cholesterol)",
            "smiles": "CC(C)CCCC(C)C1CCC2C3CC=C4CC(O)CCC4(C)C3CCC12C",
            "fragments": [255.21, 271.21, 273.22, 213.16],
            "pmz": 386.35,
        },
        "amino_acid": {
            "name": "亮氨酸 (Leucine)",
            "smiles": "CC(C)C[C@H](N)C(=O)O",
            "fragments": [86.10, 84.04, 130.07, 101.07, 74.02],
            "pmz": 132.10,
        },
        "fatty_acid": {
            "name": "棕榈酸 (Palmitic acid)",
            "smiles": "CCCCCCCCCCCCCCCC(=O)O",
            "fragments": [255.23, 60.02, 73.03, 88.05],
            "pmz": 257.25,
        },
        "benzenoid": {
            "name": "苯甲酸 (Benzoic acid)",
            "smiles": "OC(=O)c1ccccc1",
            "fragments": [105.03, 77.04, 91.05, 51.02],
            "pmz": 123.04,
        },
        "lipid": {
            "name": "PC(16:0/18:1)",
            "smiles": "",  # 复杂磷脂
            "fragments": [184.07, 104.11, 255.23, 281.25],
            "pmz": 760.59,
        },
    }

    @pytest.mark.parametrize("class_key", list(TEST_CASES.keys()))
    def test_fragment_classification_correct(self, class_key):
        """每个类别的碎片应被正确分类到该类别"""
        case = self.TEST_CASES[class_key]
        result = classify_compound_from_fragments(case["fragments"], tolerance=0.5)
        assert result, f"{case['name']}: Should return a classification result"
        classified = result["class_name"]
        assert classified == class_key, \
            f"{case['name']}: Expected '{class_key}', got '{classified}' " \
            f"(scores: {result.get('all_scores', {})})"

    @pytest.mark.parametrize("class_key", list(TEST_CASES.keys()))
    def test_expert_classify_correct(self, class_key):
        """每个类别通过expert_evaluate应得到正确的effective_classification"""
        case = self.TEST_CASES[class_key]
        peaks = [[mz, 100.0] for mz in case["fragments"]]
        result = expert_evaluate(
            query_peaks=peaks,
            query_pmz=case["pmz"],
            candidate_smiles=case["smiles"],
            candidate_classification=class_key,
        )
        effective = result["effective_classification"]
        assert effective == class_key, \
            f"{case['name']}: Expected effective='{class_key}', got '{effective}'"


# ═══════════════════════════════════════════════════════════════
# 8. 碎片重叠与排他性测试 (v3.3.0 新增)
# ═══════════════════════════════════════════════════════════════

class TestFragmentExclusivity:
    """验证碎片重叠问题已修复"""

    def test_saponin_no_steroid_fragments(self):
        """皂苷类特征碎片不应包含甾体碎片"""
        profile = get_profile("saponin")
        assert profile is not None
        steroid_mzs = {415.3212, 397.3107, 282.2168, 255.2113, 271.2113}
        for frag in profile.characteristic_fragments:
            assert frag.mz not in steroid_mzs, \
                f"Saponin should not contain steroid fragment m/z {frag.mz} ({frag.name})"

    def test_steroid_has_diosgenin_fragments(self):
        """甾体类应包含薯蓣皂苷元碎片"""
        profile = get_profile("steroid")
        assert profile is not None
        steroid_mzs = {f.mz for f in profile.characteristic_fragments}
        assert 415.3212 in steroid_mzs, "Steroid should have diosgenin m/z 415.3212"
        assert 255.2113 in steroid_mzs, "Steroid should have D-ring m/z 255.2113"

    def test_no_flavonoid_sugar_fragments(self):
        """黄酮类不应包含糖类特征碎片(氧碳鎓离子)"""
        profile = get_profile("flavonoid")
        assert profile is not None
        sugar_oxocarbenium = {163.0606, 147.0657, 133.0501}
        for frag in profile.characteristic_fragments:
            assert frag.mz not in sugar_oxocarbenium, \
                f"Flavonoid should not contain sugar oxocarbenium m/z {frag.mz}"

    def test_class_profiles_count(self):
        """类别档案数量应至少14个"""
        summary = get_library_summary()
        n_classes = summary["total_compound_classes"]
        assert n_classes >= 14, \
            f"Should have at least 14 class profiles, got {n_classes}: {summary['classes_covered']}"

    def test_all_new_classes_have_profiles(self):
        """新增5个类别都应有完整档案"""
        for cls in ["terpenoid", "steroid", "amino_acid", "lipid", "benzenoid"]:
            profile = get_profile(cls)
            assert profile is not None, f"Missing profile for '{cls}'"
            assert len(profile.characteristic_fragments) >= 5, \
                f"'{cls}' should have at least 5 fragments, got {len(profile.characteristic_fragments)}"
            assert len(profile.diagnostic_ion_pairs) >= 1, \
                f"'{cls}' should have at least 1 diagnostic ion pair"

    def test_class_to_fgs_coverage(self):
        """_CLASS_TO_FGS 应覆盖所有类别"""
        expected_classes = {
            "flavonoid", "alkaloid", "terpenoid", "phenolic_acid",
            "glycoside", "sugar", "fatty_acid", "coumarin", "quinone",
            "steroid", "amino_acid", "benzenoid", "sulfonate",
            "lipid", "lignan", "saponin",
        }
        for cls in expected_classes:
            assert cls in _CLASS_TO_FGS, f"'{cls}' missing from _CLASS_TO_FGS"

    def test_fg_signatures_coverage(self):
        """FUNCTIONAL_GROUP_SIGNATURES 应覆盖所有类别"""
        expected_fgs = {
            "flavonoid", "alkaloid", "terpenoid", "phenolic_acid",
            "glycoside", "sugar", "fatty_acid", "coumarin", "quinone",
            "steroid", "amino_acid", "benzenoid", "sulfonate",
            "lipid", "lignan", "saponin",
        }
        for fg in expected_fgs:
            assert fg in FUNCTIONAL_GROUP_SIGNATURES, \
                f"'{fg}' missing from FUNCTIONAL_GROUP_SIGNATURES"


# ═══════════════════════════════════════════════════════════════
# 9. 三层级联分类测试 (v3.3.0 新增)
# ═══════════════════════════════════════════════════════════════

class TestThreeTierClassification:
    """测试三层级联分类: SMILES规则 → 碎片库 → Unknown"""

    def test_smiles_layer_classifies_flavonoid(self):
        """Layer 1: SMILES规则应正确分类黄酮"""
        # Flavonoid SMILES (valid for RDKit + RuleBasedClassifier)
        result = expert_evaluate(
            query_peaks=[[120.0, 100.0]],
            query_pmz=271.0,
            candidate_smiles="c1cc2c(cc1O)C(=O)C=C(O2)c3ccc(O)cc3",
            candidate_classification="Unknown",
        )
        effective = result["effective_classification"]
        assert effective != "Unknown", \
            f"SMILES-based classification should not be Unknown for flavonoid, got '{effective}'"

    def test_fragment_layer_classifies_sugar(self):
        """Layer 2: 碎片库应正确分类糖类"""
        result = expert_evaluate(
            query_peaks=[[145.9, 100.0], [129.0, 60.0], [111.0, 30.0]],
            query_pmz=165.07,
            candidate_smiles="",  # No SMILES → skip Layer 1
            candidate_classification="Unknown",
        )
        effective = result["effective_classification"]
        assert effective == "sugar", \
            f"Fragment-based should classify as sugar, got '{effective}'"

    def test_no_false_positive_for_random_fragments(self):
        """随机碎片不应产生高置信度分类"""
        # m/z values that don't match any known profile well
        result = classify_compound_from_fragments([500.0, 499.0, 498.0], tolerance=0.5)
        if result:
            assert result["confidence"] < 0.5, \
                f"Random fragments should have low confidence, got {result['confidence']}"


# ═══════════════════════════════════════════════════════════════
# 10. 结构式图片测试
# ═══════════════════════════════════════════════════════════════

class TestStructureImageGeneration:
    """测试结构式图片生成"""

    def test_sanitize_rhamnose_smiles(self):
        """鼠李糖 SMILES 应能成功规范化"""
        from web.structure import _sanitize_smiles
        smiles = "C[C@H](O1)[C@H](O)[C@@H](O)[C@@H](O)C(O)1"
        sanitized = _sanitize_smiles(smiles)
        assert sanitized is not None, "Rhamnose SMILES should sanitize successfully"
        assert len(sanitized) > 5

    def test_image_generation_rhamnose(self):
        """鼠李糖 SMILES 应能生成结构式图片 (PNG 优先, 无 Cairo 时 SVG 兜底, env 无关)"""
        from web.structure import smiles_to_data_uri
        smiles = "C[C@H](O1)[C@H](O)[C@@H](O)[C@@H](O)C(O)1"
        uri = smiles_to_data_uri(smiles)
        assert uri is not None, "Rhamnose structure image should be generated"
        assert uri.startswith("data:image/"), "Should be a valid image data URI"
        assert len(uri) > 100, "Image data URI should be substantial"

    def test_invalid_smiles_returns_none(self):
        """无效 SMILES 应返回 None"""
        from web.structure import smiles_to_image_base64
        assert smiles_to_image_base64("") is None
        assert smiles_to_image_base64(None) is None
        assert smiles_to_image_base64("N/A") is None
        assert smiles_to_image_base64("INVALID_XXXXXXX") is None


# ═══════════════════════════════════════════════════════════════
# 11. 糖苷糖基顺序丢失级联 (v3.7.0 新增)
# ═══════════════════════════════════════════════════════════════

class TestSugarLossCascade:
    """糖苷/皂苷糖基顺序丢失级联检测

    化学依据: Domon-Costello 糖链裂解命名 + 《中药质谱分析》糖苷类糖链顺序丢失规律
    核心改进: 从前体起链 (而非从最大碎片起链), 捕获完整 [M+H]→逐级脱糖→苷元 级联。
    """

    def setup_method(self):
        self.scorer = ChemicalExpertScorer()

    def test_saponin_hexose_ladder(self):
        """人参皂苷: 3 个连续己糖丢失 → 苷元"""
        fragments = [947.5, 785.4, 623.4, 423.4]
        cascade = self.scorer._detect_sugar_loss_cascade(
            1109.6, fragments, classification="saponin"
        )
        assert cascade is not None, "Saponin should yield a sugar-loss cascade"
        assert cascade["n_sugars"] >= 2, f"Expected >=2 hexose losses, got {cascade['n_sugars']}"
        sugar_types = [s["sugar_type"] for s in cascade["cascade"]]
        assert "hexose" in sugar_types
        # 苷元应在最小碎片附近
        assert abs(cascade["aglycone_mz"] - 623.4) < 1.0

    def test_rutin_deoxyhexose_then_hexose(self):
        """芦丁: 先脱脱氧己糖(鼠李糖 146) 再脱己糖(葡萄糖 162) → 槲皮素苷元 303"""
        fragments = [465.10, 303.05]
        cascade = self.scorer._detect_sugar_loss_cascade(
            611.16, fragments, classification="flavonoid"
        )
        assert cascade is not None
        assert cascade["n_sugars"] == 2
        sugar_types = [s["sugar_type"] for s in cascade["cascade"]]
        assert "deoxyhexose" in sugar_types and "hexose" in sugar_types
        assert abs(cascade["aglycone_mz"] - 303.05) < 1.0

    def test_cascade_starts_from_precursor(self):
        """级联首步必须从前体出发 (修复: 旧算法从最大碎片起链丢失首步)"""
        fragments = [947.5, 785.4]
        cascade = self.scorer._detect_sugar_loss_cascade(
            1109.6, fragments, classification="saponin"
        )
        assert cascade is not None
        assert cascade["cascade"][0]["from_mz"] == 1109.6, \
            "Cascade must start from the precursor [M+H]+"

    def test_sugar_composition(self):
        """糖组成统计应正确计数各糖型"""
        fragments = [465.10, 303.05]  # rutin: deoxyhexose + hexose
        cascade = self.scorer._detect_sugar_loss_cascade(
            611.16, fragments, classification="flavonoid"
        )
        assert cascade is not None
        comp = cascade["sugar_composition"]
        assert comp.get("deoxyhexose", 0) == 1
        assert comp.get("hexose", 0) == 1

    def test_aglycone_no_sugar_loss(self):
        """黄酮苷元 (槲皮素, 无糖) 不应产生糖基级联"""
        fragments = [285.04, 153.02, 137.02]
        cascade = self.scorer._detect_sugar_loss_cascade(
            303.05, fragments, classification="flavonoid"
        )
        assert cascade is None, "Quercetin aglycone should not yield a sugar cascade"

    def test_non_glycoside_requires_strong_evidence(self):
        """非糖苷类: 单次巧合 162 丢失不应误报 (需 ≥2 步连续糖丢失)"""
        fragments = [300.0, 137.9]  # 单次 162 巧合丢失
        cascade = self.scorer._detect_sugar_loss_cascade(
            300.0, fragments, classification="phenolic_acid"
        )
        # 仅 1 步且非糖苷类核心 → 不报告 (避免误报)
        assert cascade is None

    def test_cascade_in_expert_evaluate_result(self):
        """expert_evaluate 结果中 neutral_loss_analysis 应包含 sugar_loss_cascade"""
        peaks = [[p, 100.0] for p in [947.5, 785.4, 623.4]]
        result = expert_evaluate(
            query_peaks=peaks, query_pmz=1109.6,
            candidate_smiles="", candidate_classification="saponin",
        )
        nl = result["neutral_loss_analysis"]
        assert "sugar_loss_cascade" in nl
        assert nl["sugar_loss_cascade"] is not None
        assert nl["sugar_loss_cascade"]["n_sugars"] >= 2

    def test_reasoning_mentions_cascade(self):
        """化学推理应包含糖基级联行"""
        peaks = [[p, 100.0] for p in [947.5, 785.4, 623.4]]
        result = expert_evaluate(
            query_peaks=peaks, query_pmz=1109.6,
            candidate_smiles="", candidate_classification="saponin",
        )
        reasoning = result["chemical_reasoning"]
        assert "Sugar-loss cascade" in reasoning or "糖" in reasoning


# ═══════════════════════════════════════════════════════════════
# 12. 乌头碱型酯基顺序丢失级联 (v3.7.0 新增)
# ═══════════════════════════════════════════════════════════════

class TestAlkaloidEsterCascade:
    """乌头碱型 C19 二萜生物碱 酯基顺序丢失级联检测

    化学依据: 《中药质谱分析》第5章 — 乌头碱骨架酯取代基特征中性丢失
    (乙酰基 -AcOH / 苯甲酰基 -PhCOOH / 甲氧基 -CH3OH / 羟基 -H2O),
    据此区分双酯型 DDA (剧毒) → 单酯型 MDA (低毒) → 乌头原碱。
    """

    def setup_method(self):
        self.scorer = ChemicalExpertScorer()

    def test_aconitine_dda_cascade(self):
        """乌头碱(双酯型): -AcOH → -PhCOOH → -CH3OH 三步酯基级联"""
        fragments = [586.31, 464.28, 432.27]
        cascade = self.scorer._detect_alkaloid_ester_loss_cascade(
            646.32, fragments, classification="alkaloid"
        )
        assert cascade is not None, "Aconitine should yield an ester-loss cascade"
        assert cascade["n_losses"] >= 2
        loss_types = {s["loss_type"] for s in cascade["cascade"]}
        assert "acetate" in loss_types, "Should detect acetyl (AcOH) loss"
        assert "benzoate" in loss_types, "Should detect benzoyl (PhCOOH) loss"

    def test_dda_subtype_classification(self):
        """同时丢失乙酰基+苯甲酰基 → 双酯型 (毒性最高)"""
        fragments = [586.31, 464.28]  # -AcOH, -PhCOOH
        cascade = self.scorer._detect_alkaloid_ester_loss_cascade(
            646.32, fragments, classification="alkaloid"
        )
        assert cascade is not None
        assert "DDA" in cascade["aconitine_subtype"], \
            f"Acetyl+benzoyl losses → DDA, got {cascade['aconitine_subtype']}"

    def test_mda_subtype_classification(self):
        """仅苯甲酰基丢失 → 单酯型 (MDA)"""
        fragments = [524.28]  # only -PhCOOH from a 646 precursor (hypothetical MDA product)
        cascade = self.scorer._detect_alkaloid_ester_loss_cascade(
            646.32, fragments, classification="alkaloid"
        )
        assert cascade is not None
        assert "MDA" in cascade["aconitine_subtype"]

    def test_caffeine_no_ester_cascade(self):
        """咖啡因 (嘌呤生物碱, 无乌头碱型酯基) 不应触发酯基级联"""
        fragments = [138.07, 110.06, 83.06]
        cascade = self.scorer._detect_alkaloid_ester_loss_cascade(
            195.09, fragments, classification="alkaloid"
        )
        assert cascade is None, "Caffeine (no aconitine-type esters) should not trigger"

    def test_non_alkaloid_not_triggered(self):
        """非生物碱类不应触发酯基级联检测器"""
        fragments = [586.31, 464.28]  # would match ester losses numerically
        cascade = self.scorer._detect_alkaloid_ester_loss_cascade(
            646.32, fragments, classification="flavonoid"
        )
        assert cascade is None, "Non-alkaloid should not trigger alkaloid ester cascade"

    def test_ester_cascade_in_expert_evaluate(self):
        """expert_evaluate 结果应包含 alkaloid_ester_cascade"""
        peaks = [[p, 100.0] for p in [586.31, 464.28, 432.27]]
        result = expert_evaluate(
            query_peaks=peaks, query_pmz=646.32,
            candidate_smiles="", candidate_classification="alkaloid",
        )
        nl = result["neutral_loss_analysis"]
        assert "alkaloid_ester_cascade" in nl
        assert nl["alkaloid_ester_cascade"] is not None
        assert nl["alkaloid_ester_cascade"]["n_losses"] >= 2


# ═══════════════════════════════════════════════════════════════
# 13. 碎裂途径结构化呈现 (报告渲染烟测, v3.7.0 新增)
# ═══════════════════════════════════════════════════════════════

class TestFragmentationPathwayRendering:
    """验证结构化碎裂途径在 HTML 报告中正确渲染"""

    def _build_result(self, peaks, pmz, cls):
        ed = expert_evaluate(
            query_peaks=peaks, query_pmz=pmz,
            candidate_smiles="", candidate_classification=cls,
        )
        return {
            "precursor_mz": pmz,
            "matches": [{
                "compound_name": cls + " sample",
                "smiles": "", "classification": cls, "inchikey": "",
                "formula": "", "cosine": 0.8, "final_score": 0.7,
                "expert_score": ed["expert_total_score"], "mass_ppm": 5.0,
                "expert_details": ed, "chemical_reasoning": ed["chemical_reasoning"],
            }],
        }

    def test_report_renders_sugar_cascade(self, tmp_path):
        """报告 HTML 应包含糖基级联"""
        from web.expert_report import generate_expert_report
        peaks = [[p, 100.0] for p in [947.5, 785.4, 623.4]]
        results = [self._build_result(peaks, 1109.6, "saponin")]
        out = tmp_path / "r.html"
        generate_expert_report(results, str(out), "test.mgf")
        html = out.read_text(encoding="utf-8")
        assert "经验碎裂证据与中性丢失序列" in html
        assert "糖基顺序丢失级联" in html
        assert "苷元" in html

    def test_report_renders_ester_cascade(self, tmp_path):
        """报告 HTML 应包含乌头碱酯基级联 + DDA 亚型"""
        from web.expert_report import generate_expert_report
        peaks = [[p, 100.0] for p in [586.31, 464.28, 432.27]]
        results = [self._build_result(peaks, 646.32, "alkaloid")]
        out = tmp_path / "r.html"
        generate_expert_report(results, str(out), "test.mgf")
        html = out.read_text(encoding="utf-8")
        assert "乌头碱型酯基丢失级联" in html
        assert "DDA" in html

    def test_report_renders_pathway_and_diagnostic_tables(self, tmp_path):
        """报告 HTML 应包含碎裂步骤表与诊断离子表"""
        from web.expert_report import generate_expert_report
        peaks = [[303.05, 999], [153.02, 500], [137.02, 300]]
        results = [self._build_result(peaks, 303.05, "flavonoid")]
        out = tmp_path / "r.html"
        generate_expert_report(results, str(out), "test.mgf")
        html = out.read_text(encoding="utf-8")
        assert "碎裂步骤" in html or "诊断离子" in html or "中性丢失" in html


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
