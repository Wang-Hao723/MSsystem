"""
Matching 模块 — 质谱匹配引擎 + 化学专家系统 + 碎片特征标准库

核心组件:
  - search: 质谱库搜索与匹配引擎
  - scores: 评分函数(余弦/质量/碎片/强度/专家)
  - model_predictor: MS2→化学结构预测模型
  - blank_subtraction: 空白溶剂峰扣除

专家系统(v2.0):
  - fragmentation_rules: 裂解规则知识库
  - expert_scorer: 化学专家评分器
  - expert_metrics: 专家评估指标
  - expert_loss: 专家引导的损失函数
  - expert_logger: 专家评估日志系统

增强系统(v3.0):
  - fragment_feature_library: 碎片特征标准库(基于3份权威资料)
  - structure_elucidation: 结构解析推理引擎
  - data_quality: 质谱数据质量控制(Spectraverse标准)
"""

import logging as _logging
_logger = _logging.getLogger(__name__)

from .search import search, load_library, find_matching_fragments
from .scores import (
    cosine_score, mass_score, mass_ppm, fragment_score,
    intensity_score, composite_score, enhanced_composite_score, confidence,
)

# 空白扣除
from .blank_subtraction import (
    BlankProfile, BlankSpectrum, load_blank_mgf, load_blank_file, subtract_blank,
    create_blank_profile_from_list, BLANK_SUBTRACTION_AVAILABLE,
)

# ── v3.0 增强模块(优先加载，作为其他模块的依赖) ──
try:
    from .fragment_feature_library import (
        FragmentFeature, CompoundClassProfile, AdductProfile,
        COMPOUND_CLASS_PROFILES, ADDUCT_PROFILES,
        FLAVONOID_FRAGMENTS, SAPONIN_FRAGMENTS, ALKALOID_FRAGMENTS,
        PHENOLIC_ACID_FRAGMENTS, COUMARIN_FRAGMENTS,
        LIGNAN_FRAGMENTS, QUINONE_FRAGMENTS, FATTY_ACID_FRAGMENTS,
        DATA_QUALITY_STANDARDS, TCM_PROCESSING_FEATURES,
        STRUCTURE_ELUCIDATION_RULES,
        get_profile, get_fragments_by_class, get_adduct_profile,
        match_fragment_to_class, classify_compound_from_fragments,
        get_diagnostic_ion_pairs, get_fragmentation_pathway,
        get_isomer_rules, list_all_classes,
        get_data_quality_standards, get_tcm_processing_features,
        get_structure_elucidation_rules,
        search_fragment_by_formula, export_library, get_library_summary,
    )
    FEATURE_LIBRARY_AVAILABLE = True
except ImportError as e:
    FEATURE_LIBRARY_AVAILABLE = False
    logger = __import__("logging").getLogger(__name__)
    logger.warning(f"Fragment feature library not available: {e}")

try:
    from .structure_elucidation import (
        StructureElucidationEngine, ElucidationReport,
        ElucidationEvidence, FragmentationPath, StructureHypothesis,
        get_elucidation_engine, elucidate_structure, explain_match,
    )
    ELUCIDATION_AVAILABLE = True
except ImportError as e:
    ELUCIDATION_AVAILABLE = False

try:
    from .data_quality import (
        DataQualityController, QualityReport, BatchQualityReport,
        get_quality_controller, quick_quality_check,
    )
    DATA_QUALITY_AVAILABLE = True
except ImportError:
    DATA_QUALITY_AVAILABLE = False

# ── 专家系统模块(v2.0) ──
try:
    from .fragmentation_rules import (
        FRAGMENTATION_RULES, COMMON_NEUTRAL_LOSSES, DIAGNOSTIC_IONS,
        FUNCTIONAL_GROUP_SIGNATURES, BOND_DISSOCIATION_ENERGIES,
        FragmentationRule, NeutralLoss, DiagnosticIon, BondDissociationEnergy,
        get_neutral_loss_by_mass, get_rules_by_mechanism,
        get_rules_by_functional_group, match_diagnostic_ion,
        compute_possible_neutral_losses, compute_rule_coverage,
        estimate_bond_lability, get_fg_signature, get_rules_summary,
        normalize_classification,
    )
    EXPERT_KB_AVAILABLE = True
except ImportError:
    EXPERT_KB_AVAILABLE = False

try:
    from .expert_scorer import (
        ChemicalExpertScorer, get_expert_scorer, expert_evaluate,
    )
    EXPERT_SCORER_AVAILABLE = True
except ImportError:
    EXPERT_SCORER_AVAILABLE = False

try:
    from .expert_metrics import (
        compute_all_expert_metrics, ExpertMetricSet,
        fragmentation_rule_coverage, structural_reasoning_score,
        functional_group_diagnostic_score, neutral_loss_rationality,
        bond_energy_consistency, evaluate_model_performance,
        ExpertMetricCallback,
    )
    EXPERT_METRICS_AVAILABLE = True
except ImportError:
    EXPERT_METRICS_AVAILABLE = False

try:
    from .expert_loss import (
        CompositeExpertLoss, FragmentRuleLoss,
        NeutralLossConsistencyLoss, BondEnergyRegularizer,
        FunctionalGroupSignatureLoss, StructureFragmentAlignmentLoss,
        create_expert_training_bias,
    )
    EXPERT_LOSS_AVAILABLE = True
except ImportError:
    EXPERT_LOSS_AVAILABLE = False

try:
    from .expert_logger import (
        ExpertLogger, get_expert_logger, log_match,
    )
    EXPERT_LOGGER_AVAILABLE = True
except ImportError:
    EXPERT_LOGGER_AVAILABLE = False


# ─────────────────────────────────────────────────────────────────
# 系统状态摘要
# ─────────────────────────────────────────────────────────────────

def get_system_status() -> dict:
    """获取系统各模块可用状态"""
    return {
        "version": "3.0.0",
        "modules": {
            "feature_library": FEATURE_LIBRARY_AVAILABLE,
            "structure_elucidation": ELUCIDATION_AVAILABLE,
            "data_quality": DATA_QUALITY_AVAILABLE,
            "expert_kb": EXPERT_KB_AVAILABLE,
            "expert_scorer": EXPERT_SCORER_AVAILABLE,
            "expert_metrics": EXPERT_METRICS_AVAILABLE,
            "expert_loss": EXPERT_LOSS_AVAILABLE,
            "expert_logger": EXPERT_LOGGER_AVAILABLE,
            "blank_subtraction": BLANK_SUBTRACTION_AVAILABLE,
        },
        "knowledge_base": get_library_summary() if FEATURE_LIBRARY_AVAILABLE else {},
    }
