"""全局配置参数。"""
import os
from pathlib import Path

# ── 路径 ──
PROJECT_ROOT = Path(__file__).resolve().parent.parent


def _env_path(name: str, default: Path) -> Path:
    value = os.environ.get(name, "").strip()
    return Path(value).expanduser().resolve() if value else default.resolve()


def _env_bool(name: str, default: bool) -> bool:
    value = os.environ.get(name)
    if value is None:
        return default
    return value.strip().lower() in {"1", "true", "yes", "on"}


# Full/Small 谱库、映射和用户投稿属于跨版本共享的权威数据。默认放在代码仓库
# 之外，避免 V3/V5/V6 任一目录被删除或复制时连带复制/丢失数据库。环境变量可
# 覆盖此位置，便于将来迁移到外置磁盘。
_CANONICAL_SHARED_ROOT = PROJECT_ROOT.parent / "Mssystem_shared_data"
_LEGACY_DATA_ROOT = PROJECT_ROOT / "data"
SHARED_DATA_ROOT = _env_path(
    "MSSYSTEM_SHARED_DATA_ROOT",
    _CANONICAL_SHARED_ROOT if _CANONICAL_SHARED_ROOT.exists() else _LEGACY_DATA_ROOT,
)

# 分析任务、上传暂存、缓存和日志必须按版本隔离，不能放入共享谱库目录。
RUNTIME_ROOT = _env_path("MSSYSTEM_RUNTIME_ROOT", PROJECT_ROOT / "output")
VENDOR_CONVERSION_CACHE_ROOT = _env_path(
    "MSSYSTEM_VENDOR_CONVERSION_CACHE_ROOT",
    RUNTIME_ROOT / "vendor_conversion_cache",
)
DB_PATH = SHARED_DATA_ROOT / "local_library.db"
DB_PATH_SMALL = SHARED_DATA_ROOT / "local_library_small.db"
LIBRARY_MAPPING_PATH = SHARED_DATA_ROOT / "library_mapping.csv"
LIBRARY_MAPPING_SMALL_PATH = SHARED_DATA_ROOT / "library_mapping_small.csv"
USER_CONTRIB_DIR = SHARED_DATA_ROOT / "user_contrib"

# V16 与仍在运行的 V9/V10/V11/V12/V14/V15 共用权威谱库，因此默认只读，避免多个服务同时导库。
# 只有在旧写入实例已停用并由操作者显式确认后，才可设
# MSSYSTEM_LIBRARY_WRITE_ENABLED=1 将写权限移交给 V16。
LIBRARY_WRITE_ENABLED = _env_bool("MSSYSTEM_LIBRARY_WRITE_ENABLED", False)
MGF_DATABASE_DIR = PROJECT_ROOT / "mgf_msp_database"
DEFAULT_SAMPLES_DIR = PROJECT_ROOT.parent / "samples"
DEFAULT_OUTPUT_DIR = RUNTIME_ROOT
RT_STANDARDS_PATH = SHARED_DATA_ROOT / "rt_standards.csv"
# The existing filename is retained so the verified MassBank-only v1 sidecar
# remains readable.  A v2 rebuild can aggregate every MSP/MGF source into the
# same compact sidecar without duplicating the production spectral database.
CE_PROFILE_DB = _env_path(
    "MSSYSTEM_CE_PROFILE_DB",
    SHARED_DATA_ROOT / "massbank_ce_profiles.db",
)
MASSBANK_CE_PROFILE_DB = CE_PROFILE_DB  # backward-compatible name

# ── 匹配参数 ──
# v3.4.1 调优: 实验数据匹配标准品库
#   - TOP_N_PEAKS 50→15: 库谱平均6峰, 50峰ModifiedCosine被噪声稀释→cosine≈0.14;
#     15峰保留足够碎片信息且显著降低稀释(15vs6→cosine≈0.28)
#   - REL_THRESHOLD 0.01→0.05: 过滤低强度噪声峰
#   - COSINE_THRESHOLD 0.6→0.3: 实验谱≠标准品(仪器/条件差异), 需容忍
#   - FINAL_SCORE_THRESHOLD 0.4→0.25: 配合cosine放宽
#   实测: 4.4% → 31.6% 匹配率提升
# High-resolution precursor pre-filter.  The old 0.5 Da floor overrode the
# 10 ppm setting (e.g. a 500 m/z query searched roughly +/-1000 ppm), allowing
# nominal-mass neighbours to enter before MS/MS scoring.  Keep a small absolute
# floor for low-m/z/vendor rounding while retaining the ppm window at high m/z.
MZ_TOLERANCE_DA = 0.02
MZ_TOLERANCE_PPM = 10.0
COSINE_THRESHOLD = 0.3
FINAL_SCORE_THRESHOLD = 0.25
MZ_BUCKET_DALTON = 0.1
TOP_N_PEAKS = 15
REL_THRESHOLD = 0.05
MAX_CANDIDATES_PER_QUERY = 500
TOP_K_MATCHES = 1

# ── v6: 自适应谱图视图 ──────────────────────────────────────────────
# 主检索继续使用经过验证的 Top15/5% 视图；仅当主结果为空或证据较弱时，才用
# Top40/1% 做一次补充检索。这样既不会用更多噪声扰动已有高分命中，也能找回
# 被 Top15 截掉的低丰度诊断离子。补充结果必须显著优于主结果才会被采用。
USE_ADAPTIVE_PEAK_SEARCH = True
# Top40 对 MIX33 仅改变少量弱命中，却让运行时间增加约 2.7 倍；保留为研究开关，
# 生产默认关闭。低 m/z 干扰重试仍启用，它对木兰花碱有明确净收益。
USE_BROAD_PEAK_RETRY = False
ADAPTIVE_TOP_N_PEAKS = 40
ADAPTIVE_REL_THRESHOLD = 0.01
ADAPTIVE_RETRY_MAX_COSINE = 0.60
ADAPTIVE_MIN_COSINE_GAIN = 0.05
ADAPTIVE_MAX_MASS_PPM = 50.0

# 混标中单个低质量碎片可能因共碎裂成为基峰（MIX33 木兰花碱的 m/z 58.065
# 即为实测案例）。只在弱结果时尝试移除一个“低 m/z 且显著独占”的基峰，
# 并仍需通过上面的质量、reverse 和显著增益守卫；原始谱和主结果均保留。
USE_LOW_MZ_INTERFERENCE_RETRY = True
LOW_MZ_INTERFERENCE_MAX_MZ = 60.0
LOW_MZ_INTERFERENCE_DOMINANCE = 4.0

# ── v5: 可选的弱证据假阳性门槛 ─────────────────────────────────────
# FINAL_SCORE_THRESHOLD=0.25 用于常规探索性分析。严格报告模式可额外要求
# entropy / reverse / ppm 等证据，但不能在用户选择 0.30/0.25 后暗中再把 final
# 提高到 0.50，因此默认关闭并由每个分析任务显式启用。
USE_STRICT_EVIDENCE_GATE = False
STRICT_MIN_FINAL = 0.50
STRICT_MIN_ENTROPY = 0.45
STRICT_MIN_REVERSE_MATCH = 0.55
STRICT_MAX_MASS_PPM = 100.0
# raw_blank_analysis HTML 只导出 top-5，但其前体来自高分辨原始数据并保留 4 位小数。
# 对这种证据受限输入单独收紧质量门禁，避免 40–100 ppm 的同名义质量/异分子式条目
# 仅凭共享碎片进入 Top-K（例如 329.2346 查询混入 [M-H]- 329.2486 的 C22:5 DPA）。
LIMITED_FRAGMENT_STRICT_MAX_MASS_PPM = 20.0

# ── 空白扣除 ──
BLANK_MZ_TOL_DA = 0.01
BLANK_MZ_TOL_PPM = 10
# 特征级空白证据：同极性、前体质量与 RT 对齐后，样品/空白 3 倍作为
# “样品富集”标签阈值，不再是默认谱库检索门槛。旧的碎片强度减法参数
# 保留用于读取历史任务；新任务不跨针直接相减未归一化的 MS/MS 强度。
DEFAULT_BLANK_SAMPLE_RATIO = 3.0
DEFAULT_BLANK_RT_TOLERANCE_MINUTES = 0.30
DEFAULT_BLANK_MIN_DETECTION_FRACTION = 0.50
# V7 default: blank evidence annotates every searched feature instead of
# suppressing library identification before the compound identity is known.
# ``hard_legacy`` remains an explicit compatibility option for old workflows.
DEFAULT_BLANK_FILTER_POLICY = "soft_evidence"
DEFAULT_BLANK_CLEAN_MAX_RATIO = 1.0
DEFAULT_BLANK_CLEAN_MIN_REPLICATES = 2
DEFAULT_BLANK_MODE = "absolute"
DEFAULT_BLANK_INTENSITY_MODE = "max"
DEFAULT_BLANK_THRESHOLD = 0.1

# ── 评分权重 ──
# v3.5 重校准: 谱图相似度主导(cos+ent=0.70)。
#   实验谱 vs 标准品库存在仪器/校准/条件差异, 使 mass/frag/inten 子项系统性偏低;
#   谱图整体形貌相似度(cosine + entropy)是质谱鉴定的金标准(GNPS/MassBank 等均以此主导),
#   应主导综合分, 而非被噪声大的子项稀释。真匹配(cos>=0.7)final 0.43→0.59, 区分度+73%。
#   (Li et al. 2021 Nature Methods 验证 entropy 系统性优于 cosine, 故两者并重主导)
WEIGHT_COSINE = 0.40        # 0.20→0.40  谱图余弦相似度(主导)
WEIGHT_ENTROPY = 0.30       # 0.25→0.30  熵相似度(Li 2021 Nature Methods, 主导)
WEIGHT_MASS = 0.10          # 0.15→0.10  质量精度(仪器差异敏感, 降权)
WEIGHT_FRAGMENT = 0.10      # 0.15→0.10  碎片重叠
WEIGHT_INTENSITY = 0.05     # 0.10→0.05  强度分布(半数匹配峰<3得0, 降权)
WEIGHT_EXPERT = 0.05        # 仅保留旧结果可读性；V11 排序贡献由下方硬开关归零

# ── 专家系统配置 ──
EXPERT_SYSTEM_ENABLED = True
# V11 旗舰排序的硬边界：专家系统继续计算并展示，但不能改变候选入选、
# scan-level 分数或实体级身份。历史消融显示专家项使 MIX33 Top-1 从
# 19/33 降到 17/33，因此默认关闭其排序权限，而不是关闭解释本身。
EXPERT_RANK_INFLUENCE_REQUESTED = _env_bool(
    "MSSYSTEM_EXPERT_RANK_INFLUENCE_ENABLED", False
)
# Release invariant, not a configurable deployment preference.  The requested
# value is retained above for audit, but no environment variable can reactivate
# an empirically harmful rank path in V11.
EXPERT_RANK_INFLUENCE_ENABLED = False
EXPERT_NEUTRAL_LOSS_TOLERANCE = 0.5       # 中性丢失质量容差(Da)
EXPERT_DIAGNOSTIC_ION_TOLERANCE = 0.05    # 诊断离子质量容差(Da) — 放宽以覆盖更多真实质谱数据
EXPERT_MIN_FRAGMENTS_FOR_ANALYSIS = 3     # 最少碎片数才触发详细分析
EXPERT_CONFIDENCE_HIGH = 0.70             # 专家评分高置信度阈值
EXPERT_CONFIDENCE_MEDIUM = 0.45           # 专家评分中等置信度阈值
EXPERT_LOG_ENABLED = True                 # 是否启用专家评估日志
EXPERT_LOG_DIR = "logs"                   # 专家日志目录

# 无结构馆藏谱默认不进入生产化合物排序。它们仍保留在 SQLite 中供专项审计，
# 但与同名有结构记录竞争时会把已知化合物降成 Spectral Lead，并显著增加 Full
# 内存/密集窗开销。未匹配查询的类别线索由 MSClassify 单独提供。
LOAD_NAME_ONLY_LEADS = False
EXPERT_INTENSITY_FLOOR = 0.01             # v3 专家匹配强度地板 (<1% base peak 视为噪声, 不参与诊断离子/中性丢失匹配)
# ── v3 Phase 4: v3 KB 专家加分 (纯增量; 关掉=纯 v2 行为) ──
# 取 max(base_expert, v3_score) 融合进 expert_score, 让 0.05 专家权重承载真实化学信号。
USE_EXPERT_V3 = True                      # 启用 v3 KB 专家加分 + 结构推理
# ── v3 Phase 5: Top-K 候选 (纯增量; 让用户看到最优候选之外的近邻, 含异构体) ──
TOP_K_CANDIDATES = 20                     # V11 实体排序候选池；显示层仍可分页/截断

# ── V7: Top-K 多指标谱图审计（只读）────────────────────────────────
# 生产排序仍由既有 Modified Cosine + entropy + mass/fragment/intensity/expert
# 组合分决定。这里仅在候选已选定后计算 direct/modified Hungarian 和 neutral-loss
# cosine，帮助判断高相似的来源；未经 compound-split 真值标定前绝不参与排序/High。
USE_MULTIMETRIC_SHADOW = _env_bool(
    "MSSYSTEM_USE_MULTIMETRIC_SHADOW", True
)
MULTIMETRIC_FRAGMENT_TOL_DA = 0.10
MULTIMETRIC_STRICT_FRAGMENT_TOL_DA = 0.02

# ── v5/v6: 异构体专属碎片证据 ────────────────────────────────────
# 文献和规则库都显示: 位置/连接异构体有时能靠特征 product ions 或相对丰度区分,
# 但全谱 cosine/entropy 会被大量共享骨架碎片主导; 旧 expert 的
# isomer_differentiation_score 最终权重约 0.001, 基本无法影响 Top-1。
# 这里只在"同分子式 + 不同 InChIKey"候选同时入围时, 用库谱之间的候选专属峰检查
# query 是否支持该候选。MIX33 审计发现该信号既能纠正也会误推翻 cosine 胜出者，
# 且两类情形无法靠分差/类别阈值分离；因此 V6 默认只计算并展示证据，不改 Top-1。
# 立体异构/谱图高度重叠者仍需同方法 RT/CCS/标准品确认。
USE_ISOMER_FRAGMENT_RERANK = True
# 实验开关：仅供具有完整结构真值的独立验证使用。默认 False 是生产安全边界；
# 可通过 MSSYSTEM_USE_ISOMER_FRAGMENT_RANK_OVERRIDE=1 临时恢复旧排序行为。
USE_ISOMER_FRAGMENT_RANK_OVERRIDE = _env_bool(
    "MSSYSTEM_USE_ISOMER_FRAGMENT_RANK_OVERRIDE", False
)
# 旧版实际加分权重。默认 evidence-only 模式下只作为 potential bonus 报告，不写入 rank_score。
ISOMER_FRAGMENT_RERANK_WEIGHT = 0.06
ISOMER_FRAGMENT_MIN_REL_INTENSITY = 0.05
ISOMER_FRAGMENT_MIN_SHARED_CANDIDATES = 2
ISOMER_FRAGMENT_MAX_CE_DELTA = 10.0
# Diagnostic evidence gets a broader measured-spectrum view than production
# Top-15 scoring.  It remains shadow-only, so low-abundance candidate-specific
# ions can be audited without diluting or reranking the validated search score.
ISOMER_DIAGNOSTIC_QUERY_TOP_N = 100
ISOMER_DIAGNOSTIC_QUERY_MIN_REL_INTENSITY = 0.01
# Candidate-discriminating ions use the HRMS tolerance configured for the
# query rather than the broad 0.10-Da library-search tolerance.  A 0.10-Da
# window can turn nearby, chemically different ions into false support.
ISOMER_DIAGNOSTIC_MZ_TOL_DA = 0.02
ISOMER_DIAGNOSTIC_MIN_REFERENCE_REPLICATES = 2
ISOMER_DIAGNOSTIC_MIN_PROFILE_PREVALENCE = 0.50
ISOMER_TIE_FINAL_EPS = 0.002
# A same-formula entry is a chemically possible alternative, but a very weak
# entry anywhere in Top-K must not automatically suppress a clear winner.  An
# alternative remains competitive when its cosine is within 0.05; an
# independently higher production rank also remains competitive.  This is the
# most permissive boundary in the MIX33 threshold sweep that introduced no new
# false-High standard call, including the Atractylenolide-II/Costunolide pair.
ISOMER_AMBIGUITY_MAX_RANK_GAP = 0.00
ISOMER_AMBIGUITY_MAX_COSINE_GAP = 0.05

# ── V9: IsoPair-MS 同分子式异构体排序与 High 精度门 ─────────────
# 该模型只在当前生产 Top-1 所在的“同分子式 + 可审计加合物上下文 +
# 不同 2-D connectivity”候选组内工作，不参与跨分子式或立体异构层排序。
# 缺少显式查询加合物时，公式条件化上下文只作 shadow，不影响生产排名。
# 训练正样本为同 connectivity 的不同实验谱；硬负样本为同分子式/同加合物
# 的不同 connectivity。训练、验证、测试按完整分子式家族隔离。
ISOPAIR_MODEL_PATH = _env_path(
    "MSSYSTEM_ISOPAIR_MODEL_PATH",
    PROJECT_ROOT / "models" / "isopair_metric_v2_noprec" / "isopair_metric_v1.txt",
)
ISOPAIR_CONFIG_PATH = _env_path(
    "MSSYSTEM_ISOPAIR_CONFIG_PATH",
    PROJECT_ROOT / "models" / "isopair_metric_v2_noprec" / "runtime_config.json",
)
USE_ISOPAIR_MS = _env_bool(
    "MSSYSTEM_USE_ISOPAIR_MS",
    ISOPAIR_MODEL_PATH.exists() and ISOPAIR_CONFIG_PATH.exists(),
)
# The model is still executed for counterfactual monitoring, but the 500-family
# independent confirmation did not demonstrate a statistically reliable
# retrieval gain (CI crossed zero).  Production ranking is therefore fail-
# closed until a newly locked, identity-disjoint cohort qualifies it.
ISOPAIR_RANK_INFLUENCE_REQUESTED = _env_bool(
    "MSSYSTEM_ISOPAIR_RANK_INFLUENCE_ENABLED", False
)
ISOPAIR_RANK_INFLUENCE_ENABLED = False
# The separately evaluated selective-High branch met its predeclared precision
# floor and is retained only as a veto on an already-High MS/MS result.  It can
# never promote or rerank a candidate.
USE_ISOPAIR_HIGH_PRECISION_VETO = _env_bool(
    "MSSYSTEM_USE_ISOPAIR_HIGH_PRECISION_VETO", True
)
# 在验证集锁定：仅当原 Modified-Cosine connectivity 第一/第二名分差 ≤0.10
# 时允许重排。清晰胜出者保持原结果；模型 selective 分数只作已有
# MS/MS High 的精度门，绝不将近并列 Medium 自动升为 High。
ISOPAIR_MAX_BASELINE_MARGIN = 0.10
ISOPAIR_MIN_CONNECTIVITIES = 2
ISOPAIR_MAX_CONNECTIVITIES = 12
ISOPAIR_ALLOWED_SOURCE_TIERS = (
    "experimental", "curated", "massbank_supplement"
)

# ── Multi-source collision-energy profiles ──────────────────────────
# The sidecar stores source/structure/adduct/instrument/unit-isolated multi-CE
# consensus summaries from MSP/MGF libraries. Bare CE numbers are never treated
# as eV. Evidence remains shadow-only until a compound-split identifier
# benchmark approves any ranking contribution.
USE_CE_PROFILE_SHADOW = _env_bool(
    "MSSYSTEM_USE_CE_PROFILE_SHADOW",
    _env_bool("MSSYSTEM_USE_MASSBANK_CE_SHADOW", True),
)
CE_PROFILE_MAX_NEAR_DELTA = 5.0
USE_MASSBANK_CE_SHADOW = USE_CE_PROFILE_SHADOW  # compatibility
MASSBANK_CE_MAX_NEAR_DELTA = CE_PROFILE_MAX_NEAR_DELTA  # compatibility

# ── 准确度审计修复 #1: 加合物中性质量归一化匹配 (纯增量; 关掉=纯 precursor m/z 匹配) ──
# 实测 (eval/benchmark_stratified.py): [M+Na]+ 命中仅 29-34% (vs [M+H]+ 75%), 因按 precursor
# m/z±0.5Da 匹配时 Na(+22)/NH4(+17) 加合物查不到质子化参比。91% 的 Na 加合物化合物库内也有
# [M+H] 参比 → 按"中性质量"对齐即可跨加合物命中。开启后候选 = precursor 候选 ∪ 中性质量候选。
USE_NEUTRAL_MASS_MATCHING = True
NEUTRAL_MASS_TOL_DA = 0.02                # 中性质量匹配容差 (与 MZ_TOLERANCE_DA 对齐)
# 准确度修复 #4: 加合物未知时的多假设跨加合物匹配 (遍历该极性所有常见加合物)。
# 实测 (benchmark_stratified): 质子化主导样本上 blind-多假设/fallback 均 NET 负 (-0.4~-0.6pp,
# 弱质子化匹配被错质量候选翻盘); 仅"已知加合物单假设"(#1, query_neutral_mass) 净正 (+0.6pp,
# per-adduct [M+H-H2O]+ +40pp/[M+Na]+ +8pp, 质子化零退化)。故默认关, 仅供 Na-rich/加合物多样
# 样本 opt-in。已知加合物的查询不受此开关影响 (走 #1 单假设)。
# Blind multi-adduct guessing is retained only as a legacy configuration name.
# Production search no longer uses it: cross-adduct matching requires an
# explicit, parseable query adduct whose polarity agrees with the query.
ADDUCT_TOLERANT_MATCHING = False

# 离子模式过滤 (默认开): query 极性已知时, 剔除库 adduct 反推极性不一致的候选。
# precursor m/z 匹配物理上无法区分离子模式 → 同一 m/z 既可能是 [M-H]- 也可能是 [M+H]+,
# 负离子样本会查到正离子库谱图 (假阳性, cos 必然低 → 拉低 High 占比/污染弱匹配池)。
# 库 adduct 空/无法判定 → polarity=None → 保留(不误杀); 关闭或 query 极性读不到 → 不过滤(零行为变化)。
ENFORCE_POLARITY_FILTER = True

# ── P1: 分类一致性过滤 (query 可靠类 vs 候选类跨骨架冲突 → 软降权 / 硬否决) ──
# 纯增量; 生产默认关 (measurement-first: 先在 benchmark 量出收益再开)。
# 把 DL 分类器从"展示徽章"升级为"交叉校验器" — 当前两个分类器从不参与决策, 这是最大未开发杠杆。
# 例: query 谱图 DL 可靠判为 flavonoid, 却匹配到 steroid 候选 → 跨骨架冲突 → 软降权。
USE_CONSISTENCY_FILTER = False            # 生产开关 (benchmark 用 --consistency 独立控制)
CONSISTENCY_MODE = "soft"                 # "soft" = final *= CONSISTENCY_PENALTY | "hard" = continue 否决
CONSISTENCY_PENALTY = 0.5                 # 软降权系数 (冲突候选 final 乘以此值; sweep {0.3,0.5,0.7})
CONSISTENCY_MIN_DL_CONF = 0.5             # DL 预测置信度下限才采信 query 类
CONSISTENCY_MIN_CLASS_ACC = 0.85          # DL 该类 per-class 准确率下限才采信 (动态读 ckpt class_acc; terpenoid .60 / fatty_acyl .63 自动排除)
CONSISTENCY_RELIABLE_CLASSES = ""         # 可选显式白名单 (逗号分隔; 留空 = 纯靠 class_acc 动态判定)

# ── P-A: 化合物专属裂解通路对齐（生产禁用，仅 shadow）─────────
# 首次消费 Compounds.fragmentation_pathway (scripts/mine_pathways.py 挖出的逐化合物中性
# 丢失/诊断离子)。历史 504-query 基准用同一候选的预挖矿通路评估其馆藏 query，
# 没有证明 query/source 排除，因而 +0.59pp 不能视为泛化收益。该信号现在即使被环境
# 变量误开启也不得进入 expert_score、rerank 或 confidence；shadow 仅用于重做无泄漏验证。
USE_PATHWAY_ALIGN = False
USE_PATHWAY_ALIGN_SHADOW = _env_bool("MSSYSTEM_PATHWAY_ALIGN_SHADOW", True)
PATHWAY_ALIGN_MZ_TOL = 0.05               # 通路特征 m/z 匹配容差 (Da, 与 expert_v3 MZ_TOL 一致)
PATHWAY_ALIGN_MIN_SUPPORT = 0.5           # freq 下限: <0.5 的低频碎片不采信 (单/双谱噪声)
PATHWAY_ALIGN_REQUIRE_RECURRING = True    # 仅采信 recurring==True 或带 diag_key/nl_key 的离子

# ── P-A 可选: 通路对齐 rerank 乘子 (在 final 上额外小幅放大, 突破 expert 0.05 权重上限) ──
# 与 mass_error_penalty 同型 (乘性)。仅当 USE_PATHWAY_ALIGN 已开且 max-merge 杠杆不足时启用。
# 绝不无 max-merge 单开 (乘性放大会同时抬错配候选)。默认关。
USE_PATHWAY_RERANK = False
PATHWAY_RERANK_STRENGTH = 0.03            # rerank = min(1.0, 1.0 + strength * align_score)

# ── P-B: 诊断碎片抢救 (纯增量; 关闭=完全当前行为) ──
# clean_peaks 在评分前硬过滤 rel>=0.05 / top15; 本特性在候选循环前抢救"匹配全局诊断离子集"
# (kb_v3 12 类 ~134 m/z) 的低强度碎片, 以强度地板注入 query_clean, 让下游评分可见。
# 极性过滤 (POS/NEG); 0.05Da 容差 + 策展集缓解噪声。504-query 基准
# Top-1 71.63%→72.02%，与 pathway 合用 72.22%、Top-5 不退化，因此开启。
# 实现: matching/diagnostic_rescue.py + search.py:383 后。
# 该抢救集由类别规则选择峰，会对规则评分形成选择偏差；在独立验证完成前关闭。
USE_DIAGNOSTIC_RESCUE = False
DIAGNOSTIC_RESCUE_MZ_TOL = 0.05           # 抢救匹配容差 (Da)
DIAGNOSTIC_RESCUE_INTENSITY_FLOOR = 0.05  # 抢救峰相对强度地板 (base peak 的 5%, 不主导)

# ── P-E: pathway-aware 置信度提升 (纯增量; 关闭=当前行为) ──
# 用 P-A 通路对齐信号开一条 *不依赖 cos≥0.9* 的 High 通路, 直击"高置信度偏少": cos<0.9 的真匹配
# (跨仪器谱图常态)被卡 Medium, 若化合物专属通路对齐够强 + 质量符 + 非异构歧义 → 诚实提升 Medium→High。
# 关键守卫(与 P4 一致, 避免过度声称): not isomer_ambiguous(异构 MS 物理不可分, 提升=假High) +
# cls_agreement≥0.5(类共识一致, 对齐 P4 降级条件)。集成于 honest_calibration.annotate。
# 依赖 USE_PATHWAY_ALIGN(pathway_align 须先由 P-A 算出并存入 expert_details)。
# 14标品验证中异构守卫正确挡住黄酮/生物碱异构，仅非异构真匹配可升
# (如 5-O-Methylvisammioside pathw0.529/cos0.88)。依赖上面已验证净正的通路对齐；
# 只改变诚实呈现层级，不改变检索排序。
USE_PATHWAY_CONFIDENCE = False
PATHWAY_CONF_HIGH_ALIGN = 0.5     # High: pathway_align ≥ 此 (化合物专属诊断证据)
PATHWAY_CONF_HIGH_COS = 0.8       #      AND cos ≥ 此 (放宽自 0.9; 仍需谱图相似)
PATHWAY_CONF_HIGH_FINAL = 0.55    #      AND final ≥ 此
PATHWAY_CONF_HIGH_PPM = 50        #      AND |ppm| ≤ 此 (质量须吻合)

# ── margin-aware 置信度校准 (cosine top1−top2 分差; 纯增量, 默认关) ──
# 用"清晰胜出度"(top1 比 top2 高多少)校准 compound_confidence。依据
# eval/struct_inference/validate_margin_calibration.py (443q, 零外部模型):
#   High 档整体精度仅 0.792 (21% 过度声称); margin<0.10 的 High 仅 0.636, ≥0.10 者 0.843;
#   Med 档 margin≥0.10 的精度 0.84 (≥ 现 High 0.79)。
# - USE_MARGIN_DOWNGRADE: High + margin<MARGIN_DOWNGRADE_THRESH → Medium (防过度声称; High 精度↑, High 数↓)
# - USE_MARGIN_UPGRADE  : Med + margin≥MARGIN_UPGRADE_THRESH + 非异构 + 类共识 → High (缓解"高置信偏少")
# presentation-only (只改 compound_confidence 档, 不动 final_score/排序) → 红线 11/14 零风险。
# 集成于 honest_calibration.annotate; 升档守卫同 P4/P-E (异构 MS 物理不可分不升)。
USE_MARGIN_DOWNGRADE = False
MARGIN_DOWNGRADE_THRESH = 0.10   # High + margin<此 → Medium
USE_MARGIN_UPGRADE = False
MARGIN_UPGRADE_THRESH = 0.10     # Med + margin≥此 + 守卫 → High

# ── v5: RT/标准品证据 + 分层报告 ─────────────────────────────────────
# 只采信用同一 LC 方法/本地标准品策展出来的 RT 窗口, 不默认用公共库 RT 改排序。
# 原因: 公共库 retention_time 跨色谱柱/梯度/仪器不可比, 乱用会把准确率拉低。
# 当 data/rt_standards.csv 中存在匹配标准品且 query RT 落入窗口时:
#   - 可把异构歧义候选升级为化合物级 High (RT 作为正交证据)
#   - 标注 annotation_level="Level 1-like (MS/MS+RT standard)"
# RT 不在窗口内则只做提示/降级, 不改变排名。
USE_RT_STANDARD_CONFIDENCE = True
QUERY_RT_UNIT = "auto"          # "auto" | "minute" | "second"; auto: >80 视为秒, 否则分钟
RT_DEFAULT_TOL_SEC = 30.0       # 无 tolerance 列时默认 ±30 s
RT_STRONG_TOL_SEC = 12.0        # ≤12 s 视为 strong RT confirm
RT_CONF_HIGH_COS = 0.80         # RT confirmed 升 High 的最低谱图相似要求
RT_CONF_HIGH_FINAL = 0.50
RT_CONF_HIGH_PPM = 50.0

# Strong Candidate: 不直接改成 High, 但在报告里单独统计, 对应 MSI/Schymanski 的
# high-quality putative structure / Level 2 candidate。
USE_STRONG_CANDIDATE_TIER = True
STRONG_CANDIDATE_COS = 0.80
STRONG_CANDIDATE_FINAL = 0.50
STRONG_CANDIDATE_PPM = 50.0
STRONG_CANDIDATE_CLASS_AGREE = 0.50

# P-A/pathway 只作为 Strong Candidate 的辅助证据, 不默认升 High。
PATHWAY_STRONG_ALIGN = 0.50

# 样本内多谱图一致性证据: 同一化合物在相近 RT 窗口内被多张 MS2 谱图重复命中,
# 且谱图/质量/Top-K 类共识过线时, 报告为 Strong Candidate。
# 这是实际样本中的正交重复证据; 只提升 evidence_tier, 不把 Medium 冒充为 High。
USE_SAMPLE_CONSENSUS_TIER = True
SAMPLE_CONSENSUS_MIN_SPECTRA = 2
SAMPLE_CONSENSUS_RT_WINDOW_MIN = 0.35
SAMPLE_CONSENSUS_MIN_COS = 0.75
SAMPLE_CONSENSUS_MIN_FINAL = 0.55
SAMPLE_CONSENSUS_MAX_PPM = 50.0
SAMPLE_CONSENSUS_MIN_CLASS_AGREE = 0.50

# MS2 acquisition-event accounting.  This is intentionally narrower than a
# full LC-MS feature finder: it collapses repeated DDA / stepped-energy scans
# using only precursor m/z, RT and polarity, and is labelled as approximate in
# every API/UI summary.
PRECURSOR_FEATURE_MZ_PPM = 10.0
PRECURSOR_FEATURE_MZ_FLOOR_DA = 0.002
PRECURSOR_FEATURE_MAX_RT_GAP_MIN = 0.03
PRECURSOR_FEATURE_MAX_RT_SPAN_MIN = 0.35

# ── v6: DecoID-like chimeric MS/MS deconvolution ───────────────────
# A positive L1 penalty is mandatory.  Runtime influence is permitted only for
# a genuinely mixed precursor-isolation window, a multi-component stable fit,
# and a measurable improvement of the observed target allocation.  Raw scores
# are always retained and the gain is capped until an independent truth-set
# benchmark supports wider use.
USE_CHIMERIC_DECONVOLUTION = True
CHIMERIC_DECONVOLUTION_L1_PENALTY = 0.02
CHIMERIC_DECONVOLUTION_MAX_COMPONENTS = 10
CHIMERIC_DECONVOLUTION_MIN_COSINE_GAIN = 0.03
CHIMERIC_DECONVOLUTION_MIN_FINAL_GAIN = 0.015
CHIMERIC_DECONVOLUTION_MAX_FINAL_GAIN = 0.12
# DecoID's DDA unknown-library idea: only near-pure same-run spectra are allowed
# to act as components for co-isolated ions absent from the reference library.
CHIMERIC_RUN_COMPONENT_MIN_PURITY = 0.90
CHIMERIC_RUN_COMPONENT_MIN_TIC = 10000.0

# ── v6: Ion Identity Molecular Networking ──────────────────────────
# Uses continuous MS1 XICs, not the approximate MS2 event clusters.  Thresholds
# mirror the conservative IIMN evidence regime (RT, overlap and peak-shape
# correlation); inferred identities are cross-checked against query/library
# adduct metadata before High confidence is allowed.
USE_ION_IDENTITY_NETWORKING = True
ION_IDENTITY_RT_TOLERANCE_MIN = 0.10
ION_IDENTITY_MIN_SHAPE_CORRELATION = 0.85
ION_IDENTITY_MIN_OVERLAP = 0.50
