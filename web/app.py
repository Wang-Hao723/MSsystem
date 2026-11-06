"""
MSsystem Web — LC-MS/MS Spectral Matching
启动: streamlit run web/app.py
"""
import io, csv, tempfile, time, sys, base64, json, threading, uuid, logging, zipfile, math
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import streamlit as st
import pandas as pd
from config.settings import (
    COSINE_THRESHOLD, DB_PATH, DB_PATH_SMALL, FINAL_SCORE_THRESHOLD,
    MZ_TOLERANCE_DA, MZ_TOLERANCE_PPM, TOP_N_PEAKS,
    USE_STRICT_EVIDENCE_GATE, LIMITED_FRAGMENT_STRICT_MAX_MASS_PPM,
    USE_CHIMERIC_DECONVOLUTION, USE_ION_IDENTITY_NETWORKING,
    ION_IDENTITY_RT_TOLERANCE_MIN, ION_IDENTITY_MIN_SHAPE_CORRELATION,
    ION_IDENTITY_MIN_OVERLAP,
)
from matching.search import (
    build_deconvolution_run_pool, find_matching_fragments,
)
from pipeline.query_search import search
from runtime.library_metadata_repairs import load_library_with_verified_repairs as load_library
from pipeline.mzml_metadata import query_acquisition_projection
from runtime.legacy_library_session import (
    COORDINATOR as _library_coordinator, capture_session_state,
    task_registry as _task_registry, atomic_json_write as _atomic_json_write,
)
from matching.scores import confidence
from matching.honest_calibration import annotate as _honest_annotate
from pipeline.query_input_guards import annotate_precursor_features
from pipeline.ion_identity_evidence import annotate_ion_identity_networks
from utils.collision_energy import collision_energy_evidence
from web.report_common import (
    _resolved_adduct,
    _standard_compound_name,
    classification_evidence_status,
    exact_structure_confidence,
    final_claim_badge_html,
    identification_scope,
    isomer_candidate_projection,
    reportable_identity_confidence,
    reportable_identity_label,
    reportable_identity_scope_label,
    render_diagnostic_ion_audit,
    render_isomer_candidate_competition_html,
)
from web.raw_spectrum_locator import (
    raw_spectrum_locator,
    raw_spectrum_locator_text,
)
from web.governance_blocked_projection import (
    governance_blocked_export_rows,
    project_governance_blocked_candidates,
    render_governance_blocked_html,
    zh_aliases_for,
)

# ── Logging ──
logging.basicConfig(level=logging.INFO, format="%(asctime)s  %(levelname)-8s  %(message)s", datefmt="%H:%M:%S")
_logger = logging.getLogger("mssystem.streamlit")

# ── Paths ──
PROJECT_ROOT = Path(__file__).resolve().parent.parent
OUTPUT_DIR = PROJECT_ROOT / "output"
UPLOAD_DIR = OUTPUT_DIR / "tmp_uploads"
TASK_DIR = OUTPUT_DIR / "tasks"
TASKS_FILE = OUTPUT_DIR / "tasks_index.json"
for _d in [OUTPUT_DIR, UPLOAD_DIR, TASK_DIR]:
    _d.mkdir(parents=True, exist_ok=True)

# ── Thread lock for background analysis ──
_tasks_lock = _library_coordinator.lock

# 专家系统(可选)
try:
    from matching.fragmentation_rules import get_rules_summary
    EXPERT_KB_AVAILABLE = True
except ImportError:
    EXPERT_KB_AVAILABLE = False

# ═══════════════════════════════════════════════════════════════════
# Custom CSS Injection
# ═══════════════════════════════════════════════════════════════════

def _inject_custom_css():
    """Inject comprehensive custom CSS for GNPS-inspired modern design."""
    st.markdown("""
<style>
:root {
    --primary: #2C7A7B;
    --primary-light: #E6FFFA;
    --primary-dark: #1A5253;
    --secondary: #4A5568;
    --bg: #F7FAFC;
    --card-bg: #FFFFFF;
    --text: #1A202C;
    --text-secondary: #718096;
    --border: #E2E8F0;
    --success: #48BB78;
    --warning: #ECC94B;
    --danger: #ED8936;
    --unknown: #A0AEC0;
    --blue: #3182CE;
    --red: #E53E3E;
    --radius: 10px;
    --shadow: 0 1px 3px rgba(0,0,0,0.06);
    --shadow-md: 0 4px 6px rgba(0,0,0,0.07);
    --shadow-lg: 0 10px 15px -3px rgba(0,0,0,0.08);
}

/* ── Page Base ── */
.stApp {
    background: var(--bg);
}
.main > div {
    padding: 0 1.5rem;
    max-width: 1400px;
    margin: 0 auto;
}

/* ── Typography ── */
h1 { letter-spacing: -0.02em; color: var(--text); }
h2 { color: var(--text); font-weight: 600; }
h3 { color: var(--text); font-size: 1.1rem; }

/* ── Sidebar ── */
div[data-testid="stSidebarContent"] {
    padding: 1.5rem 1rem;
}
div[data-testid="stSidebarContent"] .stMarkdown {
    margin-bottom: 0.5rem;
}
.sidebar-section {
    margin-bottom: 1.5rem;
    padding-bottom: 0.75rem;
    border-bottom: 1px solid var(--border);
}
.sidebar-section:last-child {
    border-bottom: none;
}

/* ── Hero Section ── */
.hero-section {
    background: linear-gradient(135deg, var(--primary-dark) 0%, var(--primary) 60%, #389B9D 100%);
    color: white;
    padding: 2rem 2.25rem;
    border-radius: var(--radius);
    margin-bottom: 1.5rem;
    box-shadow: var(--shadow-md);
}
.hero-title {
    font-size: 1.75rem;
    font-weight: 700;
    letter-spacing: -0.02em;
    margin-bottom: 0.35rem;
}
.hero-subtitle {
    font-size: 0.9rem;
    opacity: 0.9;
    line-height: 1.5;
}

/* ── Metric Cards ── */
.metric-card {
    background: var(--card-bg);
    border: 1px solid var(--border);
    border-radius: var(--radius);
    padding: 1rem 1.25rem;
    text-align: center;
    box-shadow: var(--shadow);
    transition: box-shadow 0.2s ease, transform 0.15s ease;
}
.metric-card:hover {
    box-shadow: var(--shadow-md);
    transform: translateY(-1px);
}
.metric-value {
    font-size: 1.6rem;
    font-weight: 700;
    color: var(--primary);
    line-height: 1.2;
}
.metric-label {
    font-size: 0.7rem;
    color: var(--text-secondary);
    text-transform: uppercase;
    letter-spacing: 0.06em;
    margin-top: 0.2rem;
}
.metric-sub {
    font-size: 0.75rem;
    color: var(--text-secondary);
    margin-top: 0.15rem;
}

/* ── Result Cards ── */
.result-card {
    background: var(--card-bg);
    border: 1px solid var(--border);
    border-radius: var(--radius);
    box-shadow: var(--shadow);
    padding: 1.5rem;
    margin-bottom: 1rem;
    transition: box-shadow 0.2s ease, transform 0.15s ease;
}
.result-card:hover {
    box-shadow: var(--shadow-md);
}

/* ── Badges ── */
.badge {
    display: inline-flex;
    align-items: center;
    padding: 3px 12px;
    border-radius: 9999px;
    font-size: 0.72rem;
    font-weight: 600;
    letter-spacing: 0.02em;
    white-space: nowrap;
}
.badge-high { background: #C6F6D5; color: #276749; }
.badge-medium { background: #FEFCBF; color: #744210; }
.badge-low { background: #FEEBC8; color: #7B341E; }
.badge-unknown { background: #EDF2F7; color: #4A5568; }

/* ── Score Chips ── */
.score-chip {
    display: inline-flex;
    flex-direction: column;
    align-items: center;
    background: #F7FAFC;
    border: 1px solid var(--border);
    border-radius: 8px;
    padding: 0.4rem 0.75rem;
    text-align: center;
    min-width: 75px;
}
.score-val {
    font-size: 0.9rem;
    font-weight: 700;
    color: var(--text);
}
.score-lbl {
    font-size: 0.6rem;
    color: var(--text-secondary);
    text-transform: uppercase;
    letter-spacing: 0.04em;
}

/* ── Section Header ── */
.section-header {
    font-size: 1rem;
    font-weight: 600;
    color: var(--text);
    border-left: 3px solid var(--primary);
    padding-left: 0.75rem;
    margin: 1.25rem 0 0.75rem;
}

/* ── Export Cards ── */
.export-card {
    background: var(--card-bg);
    border: 1px solid var(--border);
    border-radius: var(--radius);
    padding: 1.25rem;
    text-align: center;
    box-shadow: var(--shadow);
    transition: box-shadow 0.2s ease;
}
.export-card:hover {
    box-shadow: var(--shadow-md);
}
.export-card-title {
    font-size: 0.9rem;
    font-weight: 600;
    color: var(--text);
}
.export-card-desc {
    font-size: 0.75rem;
    color: var(--text-secondary);
    margin: 0.4rem 0 0.75rem;
}

/* ── Tab Styling ── */
button[data-baseweb="tab"] {
    font-size: 0.85rem !important;
    font-weight: 500 !important;
    padding: 0.5rem 1.25rem !important;
    border-radius: 8px 8px 0 0 !important;
    transition: all 0.15s ease !important;
    color: var(--text-secondary) !important;
}
button[data-baseweb="tab"]:hover {
    background: var(--primary-light) !important;
    color: var(--primary) !important;
}
button[data-baseweb="tab"][aria-selected="true"] {
    color: var(--primary) !important;
    border-bottom: 2px solid var(--primary) !important;
    font-weight: 600 !important;
}

/* ── DataFrame Styling ── */
div[data-testid="stDataFrame"] {
    border: 1px solid var(--border);
    border-radius: var(--radius);
    overflow: hidden;
}
div[data-testid="stDataFrame"] table {
    font-size: 0.82rem;
}
div[data-testid="stDataFrame"] th {
    background: #F7FAFC !important;
    color: var(--secondary) !important;
    font-weight: 600 !important;
    font-size: 0.75rem !important;
    text-transform: uppercase;
    letter-spacing: 0.04em;
}

/* ── Expander Styling ── */
div[data-testid="stExpander"] {
    border: 1px solid var(--border) !important;
    border-radius: var(--radius) !important;
    box-shadow: var(--shadow);
}

/* ── Sidebar Selectbox / Button ── */
div[data-testid="stSidebarContent"] button[kind="primary"] {
    background: var(--primary) !important;
}
div[data-testid="stSidebarContent"] .stSelectbox label {
    font-size: 0.8rem;
    color: var(--text-secondary);
    font-weight: 500;
}

/* ── Empty State ── */
.empty-state {
    text-align: center;
    padding: 4rem 2rem;
    color: var(--text-secondary);
}

/* ── Upload Zone ── */
div[data-testid="stFileUploader"] {
    border: 2px dashed var(--border) !important;
    border-radius: var(--radius) !important;
    background: var(--card-bg);
    transition: border-color 0.2s ease;
}
div[data-testid="stFileUploader"]:hover {
    border-color: var(--primary) !important;
}

/* ── Animations ── */
@keyframes fadeInUp {
    from { opacity: 0; transform: translateY(12px); }
    to { opacity: 1; transform: translateY(0); }
}
@keyframes fadeIn {
    from { opacity: 0; }
    to { opacity: 1; }
}
.fade-in { animation: fadeInUp 0.4s ease-out; }
.fade-in-fast { animation: fadeIn 0.25s ease-out; }
.stagger-1 { animation-delay: 0.05s; animation-fill-mode: both; }
.stagger-2 { animation-delay: 0.1s; animation-fill-mode: both; }
.stagger-3 { animation-delay: 0.15s; animation-fill-mode: both; }
.stagger-4 { animation-delay: 0.2s; animation-fill-mode: both; }

/* ── Spinner ── */
div[data-testid="stSpinner"] > div {
    border-top-color: var(--primary) !important;
}

/* ── Info / Warning / Error box overrides ── */
div[data-testid="stNotification"] {
    border-radius: var(--radius) !important;
}

/* ── Divider ── */
hr {
    border-color: var(--border) !important;
    margin: 1rem 0 !important;
}

/* ── Sidebar Sliders ── */
div[data-testid="stSidebarContent"] div[data-testid="stSlider"] label {
    font-size: 0.8rem;
    color: var(--text-secondary);
}
</style>
""", unsafe_allow_html=True)


# ═══════════════════════════════════════════════════════════════════
# SVG Spectrum Charts (module-level)
# ═══════════════════════════════════════════════════════════════════

def _single_peak_svg(peaks, w=700, h=200):
    """Generate a single-spectrum SVG chart (query peaks only)."""
    if not peaks:
        return ""
    mx = max(p[1] for p in peaks)
    if mx <= 0:
        return ""
    mz_max = max(p[0] for p in peaks) or 1
    mz_min = min(p[0] for p in peaks) or 0
    rng = mz_max - mz_min or 1
    lm, rm, tm, bm = 55, 20, 15, 30

    lines = [
        f'<rect x="0" y="0" width="{w}" height="{h}" fill="#FAFBFC" rx="4"/>',
        f'<line x1="{lm}" y1="{h-bm}" x2="{w-rm}" y2="{h-bm}" stroke="#CBD5E0" stroke-width="0.8"/>',
        f'<line x1="{lm}" y1="{tm}" x2="{lm}" y2="{h-bm}" stroke="#CBD5E0" stroke-width="0.8"/>',
    ]
    # Grid lines
    for i in range(1, 4):
        y = tm + i * (h - bm - tm) / 4
        lines.append(f'<line x1="{lm}" y1="{y:.0f}" x2="{w-rm}" y2="{y:.0f}" stroke="#EDF2F7" stroke-width="0.5" stroke-dasharray="3,3"/>')

    # Peaks
    for mz_val, inten in sorted(peaks, key=lambda x: x[1])[:80]:
        x = lm + (mz_val - mz_min) / rng * (w - lm - rm)
        y_top = h - bm - (inten / mx) * (h - tm - bm)
        opacity = 0.75
        lines.append(f'<line x1="{x:.1f}" y1="{h-bm:.0f}" x2="{x:.1f}" y2="{y_top:.0f}" stroke="#3182CE" stroke-width="1.8" stroke-linecap="round" opacity="{opacity}"/>')
        if inten / mx > 0.6:
            lines.append(f'<circle cx="{x:.1f}" cy="{y_top:.0f}" r="2" fill="#3182CE" opacity="0.9"/>')

    # m/z ticks
    for i in range(7):
        mz_tick = mz_min + rng * i / 6
        x_tick = lm + (mz_tick - mz_min) / rng * (w - lm - rm)
        lines.append(f'<line x1="{x_tick:.0f}" y1="{h-bm}" x2="{x_tick:.0f}" y2="{h-bm+3}" stroke="#A0AEC0" stroke-width="0.5"/>')
        lines.append(f'<text x="{x_tick:.0f}" y="{h-10}" font-size="8" fill="#718096" text-anchor="middle">{mz_tick:.1f}</text>')

    # Labels
    lines.append(f'<text x="{lm + (w-lm-rm)/2}" y="{h-2}" font-size="9" fill="#718096" text-anchor="middle">m/z</text>')

    # viewBox + 内联 width:100%: 卡片现与总览等宽(colspan=13), 主列 ~1100px。
    # 必须内联 width:100% (最高优先级), 否则 SVG 自带 width 属性会抢占布局致右侧留白。
    return f'<svg width="{w}" height="{h}" viewBox="0 0 {w} {h}" preserveAspectRatio="xMidYMid meet" xmlns="http://www.w3.org/2000/svg" style="width:100%;height:auto;font-family:system-ui,sans-serif;">{"".join(lines)}</svg>'


def _dual_peak_svg(query_peaks, lib_peaks, w=700, h=240):
    """Enhanced dual-spectrum SVG: Query (blue, up) vs Library (red, down)."""
    if not query_peaks and not lib_peaks:
        return ""
    all_peaks = query_peaks + lib_peaks
    mx = max(p[1] for p in all_peaks)
    if mx <= 0:
        return ""
    mz_max = max(p[0] for p in all_peaks) or 1
    mz_min = min(p[0] for p in all_peaks) or 0
    rng = mz_max - mz_min or 1
    lm, rm, tm, bm = 55, 20, 20, 30
    mid_y = h // 2

    lines = [
        f'<rect x="0" y="0" width="{w}" height="{h}" fill="#FAFBFC" rx="4"/>',
        # Center line
        f'<line x1="{lm}" y1="{mid_y}" x2="{w-rm}" y2="{mid_y}" stroke="#E2E8F0" stroke-width="0.5"/>',
        # Vertical axis
        f'<line x1="{lm}" y1="{tm}" x2="{lm}" y2="{h-bm}" stroke="#CBD5E0" stroke-width="0.8"/>',
    ]

    # Grid lines (dashed)
    for i in range(1, 4):
        y_up = tm + i * (mid_y - tm) / 4
        y_down = mid_y + i * (h - bm - mid_y) / 4
        lines.append(f'<line x1="{lm}" y1="{y_up:.0f}" x2="{w-rm}" y2="{y_up:.0f}" stroke="#EDF2F7" stroke-width="0.5" stroke-dasharray="3,3"/>')
        lines.append(f'<line x1="{lm}" y1="{y_down:.0f}" x2="{w-rm}" y2="{y_down:.0f}" stroke="#EDF2F7" stroke-width="0.5" stroke-dasharray="3,3"/>')

    # Query peaks (blue, upward)
    for mz_val, inten in sorted(query_peaks, key=lambda x: x[1])[:80]:
        x = lm + (mz_val - mz_min) / rng * (w - lm - rm)
        y_top = mid_y - (inten / mx) * (mid_y - tm)
        lines.append(f'<line x1="{x:.1f}" y1="{mid_y:.0f}" x2="{x:.1f}" y2="{y_top:.0f}" stroke="#3182CE" stroke-width="1.8" stroke-linecap="round" opacity="0.8"/>')
        if inten / mx > 0.6:
            lines.append(f'<circle cx="{x:.1f}" cy="{y_top:.0f}" r="2" fill="#3182CE" opacity="0.9"/>')

    # Library peaks (red, downward)
    for mz_val, inten in sorted(lib_peaks, key=lambda x: x[1])[:80]:
        x = lm + (mz_val - mz_min) / rng * (w - lm - rm)
        y_bot = mid_y + (inten / mx) * (h - bm - mid_y)
        lines.append(f'<line x1="{x:.1f}" y1="{mid_y:.0f}" x2="{x:.1f}" y2="{y_bot:.0f}" stroke="#E53E3E" stroke-width="1.8" stroke-linecap="round" opacity="0.8"/>')
        if inten / mx > 0.6:
            lines.append(f'<circle cx="{x:.1f}" cy="{y_bot:.0f}" r="2" fill="#E53E3E" opacity="0.9"/>')

    # Labels
    lines.append(f'<text x="{lm-8}" y="{tm+12}" font-size="9" fill="#3182CE" text-anchor="end" font-weight="600">Query</text>')
    lines.append(f'<text x="{lm-8}" y="{h-bm-5}" font-size="9" fill="#E53E3E" text-anchor="end" font-weight="600">Library</text>')

    # m/z ticks
    for i in range(7):
        mz_tick = mz_min + rng * i / 6
        x_tick = lm + (mz_tick - mz_min) / rng * (w - lm - rm)
        lines.append(f'<line x1="{x_tick:.0f}" y1="{h-bm}" x2="{x_tick:.0f}" y2="{h-bm+3}" stroke="#A0AEC0" stroke-width="0.5"/>')
        lines.append(f'<text x="{x_tick:.0f}" y="{h-10}" font-size="8" fill="#718096" text-anchor="middle">{mz_tick:.1f}</text>')

    lines.append(f'<text x="{lm + (w-lm-rm)/2}" y="{h-2}" font-size="9" fill="#718096" text-anchor="middle">m/z</text>')

    # viewBox + 内联 width:100%: 卡片现与总览等宽(colspan=13), 主列 ~1100px。
    # 必须内联 width:100% (最高优先级), 否则 SVG 自带 width 属性会抢占布局致右侧留白。
    return f'<svg width="{w}" height="{h}" viewBox="0 0 {w} {h}" preserveAspectRatio="xMidYMid meet" xmlns="http://www.w3.org/2000/svg" style="width:100%;height:auto;font-family:system-ui,sans-serif;">{"".join(lines)}</svg>'


# ═══════════════════════════════════════════════════════════════════
# 化合物中文名别名（展示与搜索辅助层）：
# 仅用于报告显示/检索索引，不参与打分、候选生成与任何治理判定（G1-G3 通用性不变）。
# key 为英文标准名小写子串；查不准译名的化合物不入表。
# ═══════════════════════════════════════════════════════════════════
ZH_COMPOUND_ALIASES = {
    "sinomenine": ["青藤碱", "清风藤碱"],
    "magnoflorine": ["木兰碱", "广玉兰碱"],
    "berberine": ["小檗碱"],
    "palmatine": ["巴马汀"],
    "jatrorrhizine": ["药根碱"],
    "coptisine": ["黄连碱"],
    "epiberberine": ["表小檗碱"],
    "tetrandrine": ["粉防己碱", "汉防己甲素"],
    "fangchinoline": ["防己诺林碱"],
    "dauricine": ["蝙蝠葛碱"],
    "puerarin": ["葛根素"],
    "daidzein": ["大豆素"],
    "daidzin": ["大豆苷"],
    "baicalin": ["黄芩苷"],
    "baicalein": ["黄芩素"],
    "wogonin": ["汉黄芩素"],
    "wogonoside": ["汉黄芩苷"],
    "astragaloside": ["黄芪甲苷"],
    "glycyrrhizic acid": ["甘草酸"],
    "glycyrrhizin": ["甘草酸"],
    "liquiritin": ["甘草苷"],
    "paeoniflorin": ["芍药苷"],
    "albiflorin": ["白芍苷"],
    "aconitine": ["乌头碱"],
    "mesaconitine": ["新乌头碱", "中乌头碱"],
    "hypaconitine": ["次乌头碱"],
    "benzoylaconine": ["苯甲酰乌头原碱"],
    "benzoylmesaconine": ["苯甲酰新乌头原碱"],
    "matrine": ["苦参碱"],
    "oxymatrine": ["氧化苦参碱"],
    "sophocarpine": ["槐果碱"],
    "ginsenoside": ["人参皂苷"],
    "notoginsenoside": ["三七皂苷"],
    "salidroside": ["红景天苷"],
    "tetramethylpyrazine": ["川芎嗪"],
    "ligustrazine": ["川芎嗪"],
    "ferulic acid": ["阿魏酸"],
    "tetrahydropalmatine": ["延胡索乙素"],
    "corydaline": ["延胡索甲素", "紫堇碱"],
    "cryptotanshinone": ["隐丹参酮"],
    "tanshinone": ["丹参酮"],
    "danshensu": ["丹参素"],
    "salvianolic acid": ["丹酚酸"],
    "protocatechuic acid": ["原儿茶酸"],
    "hesperidin": ["橙皮苷"],
    "hesperetin": ["橙皮素"],
    "naringin": ["柚皮苷"],
    "neohesperidin": ["新橙皮苷"],
    "synephrine": ["辛弗林"],
    "ephedrine": ["麻黄碱"],
    "pseudoephedrine": ["伪麻黄碱"],
    "quercetin": ["槲皮素"],
    "kaempferol": ["山奈酚"],
    "rutin": ["芦丁"],
    "hyperoside": ["金丝桃苷"],
    "isoquercitrin": ["异槲皮苷"],
    "vitexin": ["牡荆素"],
    "emodin": ["大黄素"],
    "aloe-emodin": ["芦荟大黄素"],
    "rhein": ["大黄酸"],
    "chrysophanol": ["大黄酚"],
    "physcion": ["大黄素甲醚"],
    "ligustilide": ["藁本内酯"],
    "catalpol": ["梓醇"],
    "acteoside": ["毛蕊花糖苷"],
    "verbascoside": ["毛蕊花糖苷"],
    "echinacoside": ["松果菊苷"],
    "forsythiaside": ["连翘酯苷"],
    "phillyrin": ["连翘苷"],
    "arctiin": ["牛蒡苷"],
    "arctigenin": ["牛蒡苷元"],
    "chlorogenic acid": ["绿原酸"],
    "caffeic acid": ["咖啡酸"],
    "gallic acid": ["没食子酸"],
    "ellagic acid": ["鞣花酸"],
    "bergenin": ["岩白菜素"],
    "schisandrin": ["五味子素"],
    "deoxyschisandrin": ["五味子甲素"],
    "curcumin": ["姜黄素"],
    "gingerol": ["姜辣素"],
    "aucubin": ["桃叶珊瑚苷"],
    "geniposide": ["栀子苷"],
    "genipin": ["栀子苷元"],
    "ursolic acid": ["熊果酸"],
    "oleanolic acid": ["齐墩果酸"],
    "betulinic acid": ["白桦脂酸"],
    "sitosterol": ["谷甾醇"],
    "stigmasterol": ["豆甾醇"],
}


def _zh_alias_for(name):
    """返回化合物名匹配到的中文别名列表（全词匹配，展示/搜索辅助专用）。

    必须按词边界整词匹配：子串匹配会把 protosinomenine 误标为
    青藤碱、cryptotanshinone 误标为丹参酮（它们是不同化合物）。
    表内更长的 key 优先，避免同一名字叠加多个家族别名。
    """
    return zh_aliases_for(name)


# ═══════════════════════════════════════════════════════════════════
# Rich HTML Report Generator (for detailed format export)
# ═══════════════════════════════════════════════════════════════════

def _gen_rich(
    results,
    source_name=None,
    *,
    total_spectra=None,
    matched_spectra=None,
    high_supported_compounds=None,
    high_confidence_scope_groups=None,
    high_confidence_exact_structures=None,
    high_confidence_isomer_groups=None,
    governance_blocked_candidates=None,
):
    """生成 Detailed HTML 证据报告。

    ``results`` 可以是去重后的候选行；调用方可额外传入原始
    MS2 总数与命中谱图数，避免把“去重候选数”误报为谱图总数。
    """
    from html import escape as _html_escape
    from copy import deepcopy
    from matching.compound_annotations import annotate_match_compound_roles
    from matching.compound_names import normalize_match_compound_names
    from pipeline.compound_conclusions import _auxiliary_interpretation
    from web.structure import smiles_to_image_base64, smiles_to_svg
    from matching.search import find_matching_fragments as _find_matched_frags

    projected_results = []
    for source_row in results:
        row = deepcopy(source_row)
        for candidate in row.get("matches") or []:
            if isinstance(candidate, dict):
                normalize_match_compound_names(candidate)
                annotate_match_compound_roles(candidate)
                candidate["auxiliary_interpretation"] = _auxiliary_interpretation(
                    candidate, row
                )
        projected_results.append(row)
    matched = [r for r in projected_results if r.get("matches")]
    matched.sort(key=lambda r: -r["matches"][0].get("final_score", 0))
    total = int(total_spectra) if total_spectra is not None else len(results)
    matched_scan_count = (
        int(matched_spectra)
        if matched_spectra is not None else len(matched)
    )
    matched_scan_count = max(0, min(matched_scan_count, total))
    unique_candidate_count = len(matched)
    unmatched_scan_count = max(0, total - matched_scan_count)
    rate = matched_scan_count / total * 100 if total else 0
    source_file = raw_spectrum_locator({}, source_name)["source_file"]
    source_meta = _html_escape(source_file or "未记录")
    # High is counted and rendered at the result's declared identity scope.
    # Full identity remains separate, so a 2-D-connectivity High does not turn
    # into a misleading Medium merely because stereochemistry is unresolved.
    derived_high_exact = sum(
        reportable_identity_confidence(r["matches"][0]) == "High"
        and identification_scope(r["matches"][0]) != "isomer_group"
        for r in matched
    )
    derived_high_isomer = sum(
        reportable_identity_confidence(r["matches"][0]) == "High"
        and identification_scope(r["matches"][0]) == "isomer_group"
        for r in matched
    )
    n_high_exact = max(0, int(
        high_confidence_exact_structures
        if high_confidence_exact_structures is not None
        else (
            high_supported_compounds
            if high_supported_compounds is not None
            else derived_high_exact
        )
    ))
    n_high_isomer = max(0, int(
        high_confidence_isomer_groups
        if high_confidence_isomer_groups is not None else derived_high_isomer
    ))
    n_high_scope = max(0, int(
        high_confidence_scope_groups
        if high_confidence_scope_groups is not None
        else n_high_exact + n_high_isomer
    ))

    _CONF_COLORS = {  # (bg, fg) per tier
        "High": ("#dcfce7", "#166534"), "Medium": ("#fef3c7", "#92400e"),
        "Low": ("#f1f5f9", "#475569"), "Unknown": ("#f5f5f5", "#9ca3af"),
    }

    def _conf_badge(label, conf):
        bg, fg = _CONF_COLORS.get(conf, _CONF_COLORS["Unknown"])
        return (f'<span style="font-size:9px;padding:1px 6px;border-radius:3px;font-weight:600;'
                f'background:{bg};color:{fg};white-space:nowrap;border:1px solid {fg}22">{label} {conf}</span>')

    # ── helpers ─────────────────────────────────────────────────────
    def _score_color(s):
        if s >= 0.6: return "#1a5f2a"
        if s >= 0.4: return "#8b6914"
        return "#666666"

    def _ppm_color(ppm):
        if ppm is None or ppm > 20: return "#c62828"
        if ppm > 5: return "#e65100"
        return "#2e7d32"

    def _class_color(cls):
        if not cls: return "#555"
        c = cls.lower()
        if "alkaloid" in c: return "#7b2d8b"
        if "flavonoid" in c: return "#c4881a"
        if any(x in c for x in ("phenyl", "terpen", "lipid")): return "#2a6b8b"
        if any(x in c for x in ("phenol", "acid", "amino")): return "#1a6b2a"
        return "#333"

    def _fragmentation_svg(m, size=(260, 200)):
        """碎裂高亮: 从 neutral_loss_analysis 取 Δm, 枚举非环键模拟断键, 两片段精确质量
        匹配 Δm → 高亮这些键画 SVG (真正的碎裂位置推断, 不依赖预存 bond_idx)。"""
        from rdkit import Chem
        from rdkit.Chem import AllChem, Descriptors, Draw
        try:
            mol = Chem.MolFromSmiles(m.get("smiles", ""))
            if not mol:
                return "", 0
            ed = m.get("expert_details") or {}
            nl = ed.get("neutral_loss_analysis") or {}
            losses = []
            for step in (nl.get("best_pathway") or []):
                d = step.get("mass_difference")
                if d and float(d) > 0:
                    losses.append(round(float(d), 2))
            for item in (nl.get("neutral_losses") or [])[:12]:
                d = item.get("mass_difference")
                if d and float(d) > 0:
                    losses.append(round(float(d), 2))
            losses = list(dict.fromkeys(losses))[:15]
            if not losses:
                return "", 0
            total = Descriptors.ExactMolWt(mol)
            hl = []
            for bond in mol.GetBonds():
                if bond.IsInRing():
                    continue
                em = Chem.RWMol(mol)
                try:
                    em.RemoveBond(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
                    frags = Chem.GetMolFrags(em, asMols=True)
                except Exception:
                    continue
                hit = False
                for f in frags:
                    fm = Descriptors.ExactMolWt(f)
                    for L in losses:
                        if abs(fm - L) <= 0.5 or abs(total - fm - L) <= 0.5:
                            hit = True
                            break
                    if hit:
                        break
                if hit:
                    hl.append(bond.GetIdx())
            hl = list(dict.fromkeys(hl))
            AllChem.Compute2DCoords(mol)
            drawer = Draw.MolDraw2DSVG(size[0], size[1])
            drawer.DrawMolecule(mol, highlightAtoms=[], highlightBonds=hl)
            drawer.FinishDrawing()
            svg = drawer.GetDrawingText().replace('<?xml version=\'1.0\' encoding=\'iso-8859-1\'?>\n', '')
            svg = svg.replace('<svg ', '<svg style="width:100%;height:auto" ', 1)
            return svg, len(hl)
        except Exception:
            return "", 0

    def _pathway_reaction_svg(m, query_peaks=None, max_frags=4):
        """碎裂产物结构级联: CFM-ID 预测产物(读预计算缓存) → RDKit 断键推断兜底 (统一渲染器)。

        解耦版(问题3): CFM 跑出预测产物即展示, 不再以观测 m/z 0.5Da 门卡死。报告路径
        allow_docker=False(读缓存瞬时, 不卡服务器); 缓存由分析完成时后台 precompute_cfm_fragments
        预灌 → Docker 开启 + 分析完成后, 报告即显示 CFM 产物级联。"""
        try:
            from web.report_common import render_fragmentation_viz
            smiles = m.get("smiles", "")
            if not smiles:
                return ""
            adduct = str(m.get("adduct", "") or "")
            # 极性派发统一走 adduct_polarity() (修旧启发式: [M+HCOO]- 等含 '+' 负加合物被误判正离子→CFM 错模式)。
            # 此处只有 match m (无 result 级 query ion_mode), 故 adduct_polarity 兜底。
            from utils.adduct import adduct_polarity as _adduct_polarity
            ion = "-" if _adduct_polarity(adduct) == "NEG" else "+"
            ed = m.get("expert_details") or {}
            nl = ed.get("neutral_loss_analysis") or {}
            losses = []
            for s in (nl.get("best_pathway") or []) + (nl.get("sequential_losses") or []):
                try:
                    d = float(s.get("mass_difference", 0))
                    if d > 0:
                        losses.append(d)
                except Exception:
                    pass
            for it in (nl.get("neutral_losses") or []):
                try:
                    d = float(it.get("mass_difference", 0))
                    if d > 0:
                        losses.append(d)
                except Exception:
                    pass
            structure_candidates = (
                (ed.get("structure_fragment_analysis") or {}).get(
                    "plausible_fragments"
                ) or []
            )
            return render_fragmentation_viz(
                smiles,
                ion,
                losses,
                query_peaks,
                allow_docker=False,
                product_candidates=structure_candidates,
            )
        except Exception:
            return ""

    def _dimension_score(m, details, assessment, legacy_key):
        """Read the structured dimension ledger before legacy flat fields."""
        for value in (
            assessment.get("score"),
            m.get(legacy_key),
            details.get(legacy_key),
        ):
            if value is None or value == "":
                continue
            try:
                number = float(value)
            except (TypeError, ValueError):
                continue
            if number == number and abs(number) != float("inf"):
                return max(0.0, min(1.0, number))
        return 0.0

    def _expert_dims_html(m):
        """Generate expert 8-dimension mini bar chart for detailed report.
        Handles both data formats: sub-scores at top-level (batch) or inside expert_details (web).
        """
        expert = m.get("expert_score", 0)
        if expert <= 0:
            return ""
        # Extract from expert_details if sub-scores not at top level
        ed = m.get("expert_details") or {}
        dims = [
            ("neutral_loss", "Neutral Loss", "neutral_loss_score"),
            ("rule_coverage", "Rule Coverage", "rule_coverage_score"),
            ("fg_signature", "FG Signature", "fg_signature_score"),
            ("feature_library", "Class Feature Lib", "feature_library_score"),
            ("structure_fragment", "Structure-Frag", "structure_fragment_score"),
            ("diagnostic_ion", "Diagnostic Ion", "diagnostic_ion_score"),
            ("bond_energy", "Bond Energy", "bond_energy_score"),
        ]
        assessments = ed.get("dimension_assessments") or {}
        bars = ""
        for dim_key, dim_name, legacy_key in dims:
            assessment = assessments.get(dim_key) or {}
            dim_val = _dimension_score(m, ed, assessment, legacy_key)
            status = assessment.get("status") or (
                "scored" if dim_val > 0 else "no_evidence"
            )
            pct = int(dim_val * 100)
            bar_color = (
                "#CBD5E0" if status in ("not_applicable", "unavailable")
                else "#A0AEC0" if status == "no_evidence"
                else "#38A169" if dim_val >= 0.6
                else "#D69E2E" if dim_val >= 0.3
                else "#3182CE"
            )
            display_value = (
                "不适用" if status == "not_applicable"
                else "未建立" if status == "unavailable"
                else "未观察" if status == "no_evidence"
                else f"{dim_val:.2f}"
            )
            bars += (
                f'<div style="display:flex;align-items:center;gap:4px;margin:2px 0">'
                f'<span style="font-size:9px;color:#718096;width:80px;text-align:right">{dim_name}</span>'
                f'<div style="flex:1;background:#EDF2F7;border-radius:3px;height:8px">'
                f'<div style="width:{pct}%;background:{bar_color};height:100%;border-radius:3px"></div>'
                f'</div>'
                f'<span style="font-size:9px;font-weight:600;color:#2D3748;width:42px" '
                f'title="{assessment.get("reason", "")}">{display_value}</span>'
                f'</div>'
            )
        reasoning = m.get("chemical_reasoning", "") or ed.get("chemical_reasoning", "")
        if "异构" in str(reasoning) or "isomer" in str(reasoning).lower():
            reasoning = ""
        reasoning_html = ""
        if reasoning:
            reasoning_html = (
                f'<div style="margin-top:6px;padding:5px 8px;background:#FFFBEB;'
                f'border:1px solid #F6E05E;border-radius:4px;font-size:10px;color:#744210;'
                f'line-height:1.4;max-height:90px;overflow:auto">{reasoning}</div>'
            )
        # 结构化碎裂途径: 预存馆藏途径(主, 直接读 DB fragmentation_pathway) + 本谱观测对齐(辅, on-the-fly)
        pathway_html = ""
        try:
            from web.report_common import render_db_pathway
            _db_pw = render_db_pathway(m.get("pathway_json")) or ""
        except Exception:
            _db_pw = ""
        try:
            from web.expert_report import _render_fragmentation_pathway
            _obs_pw = _render_fragmentation_pathway(ed) or ""
        except Exception:
            _obs_pw = ""
        pathway_html = _db_pw + _obs_pw
        # v3 结构推理链 + Top-K候选 (独立卡片, 分区清晰不拥挤)
        v3 = ed.get("v3") or {}
        v3_lines = [
            str(x)[:80]
            for x in (v3.get("reasoning") or [])[:4]
            if "异构" not in str(x) and "isomer" not in str(x).lower()
        ]
        v3_score = v3.get("v3_score", 0)
        _seen = set()
        cand_rows = []
        for c in (m.get("top_candidates") or [])[:6]:
            _rn = (c.get("compound_name") or "").lower().split("(")[0].split("[")[0].strip()
            _nn = "".join(ch for ch in _rn if ch.isalnum())
            _kk = (
                str(c.get("inchikey") or "").strip().upper()
                if unresolved and c.get("inchikey")
                else f"{_nn}|{c.get('formula')}|{c.get('rank') if unresolved else ''}"
            )
            if _kk in _seen:
                continue
            _seen.add(_kk)
            candidate_name = _standard_compound_name(
                c.get("compound_name", ""), c
            )[:26]
            cand_rows.append(f'{candidate_name} <span style="color:#A0AEC0">({c.get("cosine",0):.3f})</span>')

        cards = []
        # 卡片1: 规则证据覆盖（紫灰）；旧聚合值仅作审计
        fusion = ed.get("score_fusion") or {}
        fusion_components = fusion.get("components") or {}
        component_bits = []
        for source_key in ("base_8d", "v3_kb", "compound_pathway"):
            component = fusion_components.get(source_key) or {}
            if not component.get("available"):
                continue
            selected = source_key == fusion.get("selected_source")
            component_bits.append(
                f'<span style="display:inline-block;padding:1px 5px;margin:0 3px 4px 0;'
                f'border:1px solid {"#805AD5" if selected else "#CBD5E0"};border-radius:999px;'
                f'font-size:8px;color:{"#553C9A" if selected else "#718096"}">'
                f'{component.get("label") or source_key}'
                f'{" · 记录来源" if selected else ""}</span>'
            )
        fusion_html = (
            "".join(component_bits)
            + '<div style="font-size:8px;color:#A0AEC0;margin-bottom:3px">'
              '历史组件值仅供审计，未校准且不是概率；下方逐维状态为当前主要呈现。</div>'
            if component_bits else ""
        )
        dimension_assessments = ed.get("dimension_assessments") or {}
        assessed_dimensions = sum(
            str(value.get("status") or "unavailable") != "unavailable"
            for value in dimension_assessments.values()
            if isinstance(value, dict)
        )
        evidence_dimensions = sum(
            str(value.get("status") or "") == "scored"
            and (
                int(value.get("evidence_count") or 0) > 0
                or float(value.get("score") or 0) > 0
            )
            for value in dimension_assessments.values()
            if isinstance(value, dict)
        )
        cards.append(
            f'<div style="padding:7px 10px;background:#FAFBFC;border:1px solid #E2E8F0;border-radius:6px">'
            f'<div style="font-size:11px;font-weight:700;color:#553C9A;margin-bottom:6px">'
            f'🧪 规则证据评估 <span style="font-weight:400;color:#A0AEC0;font-size:10px">'
            f'支持 {evidence_dimensions}/{assessed_dimensions or len(dimension_assessments) or 8} 维</span></div>'
            f'{fusion_html}{bars}</div>'
        )
        # 卡片2: 结构推理 (绿)
        if v3_lines:
            cards.append(
                f'<div style="padding:7px 10px;background:#F0FFF4;border:1px solid #9AE6B4;'
                f'border-left:3px solid #38A169;border-radius:6px">'
                f'<div style="font-size:11px;font-weight:700;color:#22543D;margin-bottom:5px">🧬 结构推理'
                f' <span style="font-weight:400;color:#48BB78;font-size:10px">一致性 {v3_score:.2f}</span></div>'
                f'<div style="font-size:10px;color:#234E52;line-height:1.75">'
                + "<br>".join("• " + x for x in v3_lines) + "</div></div>"
            )
        if reasoning_html:
            cards.append(reasoning_html)
        # 卡片3: 碎裂途径 (橙) — 级联反应图(前体高亮断键 + m/z→[-丢失Δm]→碎片) + 详细表折叠
        if pathway_html:
            _react = _pathway_reaction_svg(m)
            cards.append(
                f'<div style="padding:7px 10px;background:#FFFAF0;border:1px solid #FBD38D;'
                f'border-left:3px solid #ED8936;border-radius:6px">'
                f'<div style="font-size:11px;font-weight:700;color:#7B341E;margin-bottom:6px">⚗️ 经验碎裂证据</div>'
                f'{_react}'
                f'<details style="margin-top:6px"><summary style="font-size:9px;color:#7B341E;cursor:pointer">📋 详细碎裂步骤表</summary>{pathway_html}</details>'
                f'</div>'
            )
        # 卡片4: Top-K候选 (蓝)
        if cand_rows:
            cards.append(
                f'<div style="padding:7px 10px;background:#EBF8FF;border:1px solid #90CDF4;'
                f'border-left:3px solid #4299E1;border-radius:6px">'
                f'<div style="font-size:11px;font-weight:700;color:#2A4365;margin-bottom:5px">🏆 Top-K 候选'
                f' <span style="font-weight:400;color:#A0AEC0;font-size:10px">{len(cand_rows)} 个</span></div>'
                f'<div style="font-size:10px;color:#2C5282;line-height:1.85">' + " · ".join(cand_rows) + "</div></div>"
            )
        # 2x2 grid 占满宽度(同行等高对齐), 4卡片: 规则证据|结构推理 / 碎裂途径|Top-K
        return '<div style="margin-top:10px;display:grid;grid-template-columns:1fr 1fr;gap:8px;align-items:stretch">' + "".join(cards) + "</div>"

    # ── v5.1 detail-card section builders (teal, 2-col + full-width) ────────
    def _rc_bars(m):
        """规则证据 8 维横向条（放左栏，可读字号）。"""
        ed = m.get("expert_details") or {}
        dims = [
            ("neutral_loss", "neutral_loss_score", "中性丢失"),
            ("rule_coverage", "rule_coverage_score", "规则覆盖"),
            ("fg_signature", "fg_signature_score", "官能团特征"),
            ("feature_library", "feature_library_score", "大类特征库"),
            ("structure_fragment", "structure_fragment_score", "结构-碎片"),
            ("diagnostic_ion", "diagnostic_ion_score", "诊断离子"),
            ("bond_energy", "bond_energy_score", "键能"),
        ]
        assessments = ed.get("dimension_assessments") or {}
        out = ""
        for dim_key, k, lbl in dims:
            assessment = assessments.get(dim_key) or {}
            v = _dimension_score(m, ed, assessment, k)
            status = assessment.get("status") or (
                "scored" if v > 0 else "no_evidence"
            )
            col = (
                "#cbd5e0" if status in ("not_applicable", "unavailable")
                else "#a0aec0" if status == "no_evidence"
                else "#1a8a5a" if v >= 0.6
                else "#c08a1a" if v >= 0.3
                else "#3182ce"
            )
            evidence_count = int(assessment.get("evidence_count") or 0)
            display_value = (
                "不适用" if status == "not_applicable"
                else "未建立" if status == "unavailable"
                else "未观察" if status == "no_evidence"
                else f"{v:.2f}" + (f" · {evidence_count}项" if evidence_count else "")
            )
            out += (f'<div class="rc-bar-row"><span class="rc-bar-nm">{lbl}</span>'
                    f'<div class="rc-bar-track"><div class="rc-bar-fill" style="width:{int(v*100)}%;background:{col}"></div></div>'
                    f'<span class="rc-bar-val" title="{assessment.get("reason", "")}">{display_value}</span></div>')
        return out

    def _rc_reasoning(m):
        """结构推理链 + 化学推理 (通栏)"""
        ed = m.get("expert_details") or {}
        v3 = ed.get("v3") or {}
        lines = [
            line
            for line in (v3.get("reasoning") or [])
            if "异构" not in str(line) and "isomer" not in str(line).lower()
        ]
        score = float(v3.get("v3_score", 0) or 0)
        items = "".join(
            f"<li>{_html_escape(str(x))}</li>" for x in lines
        )
        out = ""
        if items:
            out = (f'<div class="rc-rsn-h">结构推理链 · 一致性 {score:.2f}</div>'
                   f'<ul class="rc-rsn-list">{items}</ul>')
        chem = ed.get("chemical_reasoning") or ""
        saved_class_reliable = (
            (ed.get("classification_provenance") or {}).get("reliable")
        )
        if bool(m.get("classification_reliable")) and saved_class_reliable is False:
            chem = ""
            out += (
                '<div class="rc-chem" style="background:#eff6ff;border-color:#93c5fd;color:#1e40af">'
                '当前结构分类已更新；旧任务中基于“类别证据不足”的压缩推理文字已停止展示。'
                '重新分析可刷新全部规则证据维度。</div>'
            )
        if "异构" in str(chem) or "isomer" in str(chem).lower():
            chem = ""
        if chem:
            out += f'<div class="rc-chem">{_html_escape(str(chem))}</div>'
        return out or '<div class="rc-muted">无结构推理数据</div>'

    def _evograph_formula_graph(graph):
        """Render candidate-conditioned observed-fragment composition hypotheses."""

        if not isinstance(graph, dict):
            return ""
        context = graph.get("context") or {}
        audit = graph.get("graph_audit") or {}
        nodes = graph.get("nodes") or []
        node_by_id = {
            str(node.get("node_id") or ""): node
            for node in nodes if isinstance(node, dict)
        }

        fragment_rows = []
        for node in nodes:
            if not isinstance(node, dict) or node.get("node_type") != "fragment":
                continue
            candidates = node.get("formula_candidates") or []
            best = candidates[0] if candidates else {}
            status = str(node.get("incoming_edge_status") or "unassigned")
            status_text = {
                "heuristically_supported": "启发支持",
                "ambiguous": "多解",
                "unassigned": "未归属",
            }.get(status, status)
            try:
                mz_text = f'{float(node.get("mz")):.5f}'
            except (TypeError, ValueError):
                mz_text = "—"
            try:
                ppm_text = f'{float(best.get("mass_error_ppm")):.2f}'
            except (TypeError, ValueError):
                ppm_text = "—"
            ambiguity = (
                f"{len(candidates)} 个候选组成" if len(candidates) > 1
                else "唯一保留组成" if best else "质量/母式约束下无解"
            )
            fragment_rows.append(
                f'<tr><td class="mono">{mz_text}</td>'
                f'<td class="mono">{_html_escape(str(best.get("formula") or "—"))}</td>'
                f'<td class="mono">{ppm_text}</td>'
                f'<td>{_html_escape(status_text)}</td>'
                f'<td>{_html_escape(ambiguity)}</td></tr>'
            )

        edge_rows = []
        for edge in graph.get("edges") or []:
            if not isinstance(edge, dict):
                continue
            if not edge.get("selected") or edge.get("kind") != "formula_conserving":
                continue
            parent = node_by_id.get(str(edge.get("parent_id") or ""), {})
            child = node_by_id.get(str(edge.get("child_id") or ""), {})
            prior = edge.get("common_loss_prior") or {}
            try:
                route_text = (
                    f'{float(parent.get("mz")):.5f} → '
                    f'{float(child.get("mz")):.5f}'
                )
            except (TypeError, ValueError):
                route_text = "—"
            try:
                delta_text = f'{float(edge.get("observed_loss_da")):.5f}'
            except (TypeError, ValueError):
                delta_text = "—"
            try:
                ppm_text = f'{float(edge.get("mass_error_ppm")):.2f}'
            except (TypeError, ValueError):
                ppm_text = "—"
            edge_rows.append(
                f'<tr><td class="mono">{route_text}</td>'
                f'<td class="mono">{_html_escape(str(edge.get("parent_formula") or "—"))} → '
                f'{_html_escape(str(edge.get("child_formula") or "—"))}</td>'
                f'<td class="mono">−{delta_text}</td>'
                f'<td class="mono">{_html_escape(str(edge.get("loss_formula") or "—"))}</td>'
                f'<td class="mono">{ppm_text}</td>'
                f'<td>{_html_escape(str(prior.get("name") or edge.get("selection_status") or "—"))}</td></tr>'
            )

        formula_available = graph.get("formula_graph_available") is True
        reasons = context.get("reasons") or []
        status_text = (
            f'已检查 {int(audit.get("fragment_node_count") or 0)} 个观测碎片；'
            f'保留 {int(audit.get("selected_edge_count") or 0)} 条碎裂边。'
            if formula_available else
            "未生成分子式图：" + str(
                "、".join(str(reason) for reason in reasons)
                or graph.get("reason") or "候选分子式/明确加合物/电荷条件不完整"
            )
        )
        try:
            selected_mz = f'{float(graph.get("selected_ion_mz")):.6f}'
        except (TypeError, ValueError):
            selected_mz = "—"
        try:
            search_mz = f'{float(graph.get("search_precursor_mz", graph.get("precursor_mz"))):.6f}'
        except (TypeError, ValueError):
            search_mz = "—"
        try:
            parent_ppm = f'{float(context.get("mass_consistency_error_ppm")):.2f}'
        except (TypeError, ValueError):
            parent_ppm = "—"

        parts = [
            '<div class="rc-frag-h">观测碎片组成图（EvoGraph-F）</div>',
            f'<div class="diag-note"><b>{"已完成候选条件下的分子式守恒检查" if formula_available else "当前组成解释不可用"}</b> · '
            f'{_html_escape(status_text)}<br>隔离目标 m/z {selected_mz} · '
            f'检索/计算质量 m/z {search_mz} · 母离子式质量误差 {parent_ppm} ppm</div>',
        ]
        if fragment_rows:
            parts.append(
                '<div class="rc-frag-h">观测产物离子与候选中性组成</div>'
                '<div class="diag-table-wrap"><table class="diag-table"><thead><tr>'
                '<th>观测 m/z</th><th>候选组成</th><th>误差 ppm</th><th>入边状态</th><th>组成歧义</th>'
                f'</tr></thead><tbody>{"".join(fragment_rows)}</tbody></table></div>'
            )
        if edge_rows:
            parts.append(
                '<div class="rc-frag-h">分子式守恒的碎裂边</div>'
                '<div class="diag-table-wrap"><table class="diag-table"><thead><tr>'
                '<th>观测路径</th><th>母式 → 子式</th><th>观测 Δm</th><th>丢失组成</th><th>误差 ppm</th><th>注释</th>'
                f'</tr></thead><tbody>{"".join(edge_rows)}</tbody></table></div>'
            )
        parts.append(
            '<div class="diag-footnote">边界：上述分子式是在当前谱库候选、明确加合物和电荷保持假设下的组成假设，'
            '不是已确认的产物结构或键断裂机理；启发权重不是概率。EvoGraph-F 不参与检索、排序、置信度或 High 升档。'
            'CFM/RDKit 二维产物结构属于独立附加输出。</div>'
        )
        return "".join(parts)

    def _rc_fragpath(m, query_peaks, graph=None):
        """完整展示查询谱碎裂与候选馆藏碎裂特征。"""
        cascade = _pathway_reaction_svg(m, query_peaks=query_peaks) or ""
        ed = m.get("expert_details") or {}
        nl = ed.get("neutral_loss_analysis") or {}
        current_fragmentation = (
            (m.get("auxiliary_interpretation") or {}).get("fragmentation") or {}
        )
        steps = current_fragmentation.get("steps") or []
        step_rows = []
        for s in steps:
            try:
                fm = float(s.get("from_mz", 0)); tm = float(s.get("to_mz", 0))
            except Exception:
                continue
            nlname = s.get("neutral_loss") or ""
            nlen = s.get("neutral_loss_en") or ""
            fmula = s.get("formula") or ""
            mech = s.get("mechanism") or ""
            try:
                d = float(s.get("mass_difference", 0))
            except Exception:
                d = 0.0
            step_rows.append(
                f"<tr><td class='mono'>{fm:.4f}</td><td class='mono'>{tm:.4f}</td>"
                f"<td class='mono'>{d:.4f}</td>"
                f"<td>{_html_escape(str(nlname or '—'))}</td>"
                f"<td>{_html_escape(str(nlen or '—'))}</td>"
                f"<td class='mono'>{_html_escape(str(fmula or '—'))}</td>"
                f"<td>{_html_escape(str(mech or '—'))}</td>"
                f"<td>{_html_escape(str(s.get('priority', '—')))}</td></tr>"
            )
        nl_list = current_fragmentation.get("neutral_loss_assignments") or []
        neutral_rows = []
        for item in nl_list:
            try:
                delta = float(item.get("mass_difference", 0) or 0)
                delta_text = f"{delta:.4f}"
            except (TypeError, ValueError):
                delta_text = "—"
            candidates = item.get("possible_candidates") or []
            candidate_text = "; ".join(
                str(candidate.get("name") or candidate.get("formula") or "")
                for candidate in candidates if isinstance(candidate, dict)
            )
            assignment = (
                item.get("neutral_loss") or candidate_text
                or item.get("assignment_status") or "未归属"
            )
            neutral_rows.append(
                f"<tr><td class='mono'>{delta_text}</td>"
                f"<td>{_html_escape(str(assignment))}</td>"
                f"<td class='mono'>{_html_escape(str(item.get('formula') or '—'))}</td>"
                f"<td>{_html_escape(str(item.get('mechanism') or '—'))}</td>"
                f"<td>{_html_escape(str(item.get('priority', '—')))}</td>"
                f"<td>{_html_escape(str(item.get('assignment_status') or ('explained' if item.get('neutral_loss') else 'unassigned')))}</td></tr>"
            )

        special_cascades = ""
        try:
            from web.expert_report import (
                _render_ester_cascade,
                _render_sugar_cascade,
            )
            special_cascades = (
                _render_sugar_cascade(nl.get("sugar_loss_cascade"))
                + _render_ester_cascade(nl.get("alkaloid_ester_cascade"))
            )
        except Exception:
            special_cascades = ""

        db_pathway = ""
        try:
            from web.report_common import render_db_pathway
            db_pathway = render_db_pathway(
                m.get("pathway_json"),
                query_peaks=query_peaks,
                max_items=None,
                candidate_smiles=str(m.get("smiles") or ""),
                classification=str(m.get("display_class") or m.get("classification") or ""),
                adduct=str(m.get("adduct") or ""),
            ) or ""
        except Exception:
            db_pathway = ""
        # 可视化: 优先 CFM-ID 产物结构级联(读预计算缓存); 无缓存时用 RDKit 质量匹配法断键推断(无需Docker)
        viz = cascade
        if not viz:
            rsvg, nhl = _fragmentation_svg(m)
            if rsvg and nhl > 0:
                viz = (f'<div class="rc-rdkit"><div class="rc-frag-h">断键位置推断 · RDKit 质量匹配法 '
                       f'(基于中性丢失 Δm, 标注前体上 {nhl} 个可能断裂的键)</div>'
                       f'<div class="rc-rdkit-box">{rsvg}</div></div>')
        sections = []
        formula_graph = _evograph_formula_graph(graph)
        if formula_graph:
            sections.append(formula_graph)
        if viz:
            sections.append(f'<div class="rc-frag-viz">{viz}</div>')
        if special_cascades:
            sections.append(
                '<div class="rc-frag-h">查询谱特征级联</div>'
                + special_cascades
            )
        if step_rows:
            sections.append(
                '<div class="rc-frag-h">查询谱碎裂步骤（全部已保存步骤）</div>'
                '<div class="diag-table-wrap"><table class="diag-table"><thead><tr>'
                '<th>From m/z</th><th>To m/z</th><th>Δm</th><th>中性丢失</th>'
                '<th>English</th><th>分子式</th><th>机理</th><th>优先级</th>'
                f'</tr></thead><tbody>{"".join(step_rows)}</tbody></table></div>'
            )
        if neutral_rows:
            sections.append(
                '<div class="rc-frag-h">前体→产物中性丢失（全部已保存条目）</div>'
                '<div class="diag-table-wrap"><table class="diag-table"><thead><tr>'
                '<th>Δm</th><th>归属/候选</th><th>分子式</th><th>机理</th><th>优先级</th><th>状态</th>'
                f'</tr></thead><tbody>{"".join(neutral_rows)}</tbody></table></div>'
            )
        rejected = current_fragmentation.get("rejected_legacy_assignments") or []
        if rejected:
            rejected_rows = "".join(
                f'<tr><td class="mono">{_html_escape(str(item.get("mass_difference") or "—"))}</td>'
                f'<td>{_html_escape(str(item.get("neutral_loss") or item.get("neutral_loss_en") or "—"))}</td>'
                f'<td class="mono">{_html_escape(str(item.get("formula") or "—"))}</td>'
                f'<td>{_html_escape(str(item.get("reason") or "current_policy_rejected"))}</td></tr>'
                for item in rejected
            )
            sections.append(
                '<div class="rc-frag-h">旧版归属的当前策略复核</div>'
                '<div class="diag-table-wrap"><table class="diag-table"><thead><tr>'
                '<th>Δm</th><th>旧归属</th><th>分子式</th><th>驳回原因</th>'
                f'</tr></thead><tbody>{rejected_rows}</tbody></table></div>'
            )
        if db_pathway:
            sections.append(
                '<div class="rc-frag-h">候选化合物馆藏预存碎裂指纹（全部条目）</div>'
                + db_pathway
            )
        if not sections:
            return '<div class="rc-muted">无碎裂数据</div>'
        return "".join(sections)

    def _rc_topk(m):
        """Render the complete unresolved set, otherwise the ordinary Top-K."""
        projection = isomer_candidate_projection(m)
        unresolved = projection.get("status") != "not_applicable"
        cands = (
            projection.get("candidate_members") or []
            if unresolved else m.get("top_candidates") or []
        )
        if not cands:
            return '<div class="rc-muted">无候选</div>'
        seen, cells, n = set(), "", 0
        for c in cands:
            _rn = (c.get("compound_name") or "").lower().split("(")[0].split("[")[0].strip()
            _nn = "".join(ch for ch in _rn if ch.isalnum())
            _kk = f"{_nn}|{c.get('formula')}"
            if _kk in seen:
                continue
            seen.add(_kk); n += 1
            candidate_rank = int(c.get("rank") or n)
            star = "★ 当前最佳候选 · " if candidate_rank == 1 else f"{candidate_rank}. 竞争候选 · "
            smi = c.get("smiles", "") or ""
            svg = ""
            if smi:
                try:
                    svg = smiles_to_svg(smi, size=(150, 110)) or ""
                    if svg:
                        svg = svg.replace('<?xml version=\'1.0\' encoding=\'iso-8859-1\'?>\n', '')
                        svg = svg.replace('<svg ', '<svg preserveAspectRatio="xMidYMid meet" ', 1)
                except Exception:
                    svg = ""
            win = " rc-cand-win" if candidate_rank == 1 else ""
            name = _html_escape(_standard_compound_name(
                c.get("compound_name", ""), c
            ))
            _zh = _zh_alias_for(str(c.get("compound_name") or ""))
            zh_html = (
                f' <span style="color:#0e7490;font-size:11px">{_html_escape("·".join(_zh))}</span>'
                if _zh else ""
            )
            cos = float(c.get("cosine", 0) or 0)
            final = float(c.get("final_score", 0) or 0)
            rank_score = float(c.get("rank_score", 0) or 0)
            formula = c.get("formula", "") or ""
            inchikey = _html_escape(str(c.get("inchikey") or "—"))
            source = _html_escape(str(c.get("dataset_source") or "—"))
            spectrum = _html_escape(str(c.get("lib_spectrum_id") or "—"))
            relationship = _html_escape(str(c.get("relationship") or "—"))
            cells += (f'<div class="rc-cand{win}"><div class="rc-cand-s">{svg}</div>'
                      f'<div class="rc-cand-n">{star}{name}{zh_html}</div>'
                      f'<div class="rc-cand-meta">cos <b>{cos:.4f}</b> · final <b>{final:.4f}</b> · '
                      f'rank <b>{rank_score:.4f}</b><br>{_html_escape(str(formula))} · '
                      f'InChIKey <span class="mono">{inchikey}</span><br>'
                      f'{source} / {spectrum}<br>{relationship}</div></div>')
        if not cells:
            return '<div class="rc-muted">无候选</div>'
        boundary = ""
        if unresolved:
            boundary = (
                '<div class="diag-note"><b>未解析异构体组：</b>'
                f'当前共 {int(projection.get("candidate_count") or len(cands))} 个候选；'
                '第 1 名仅表示本次保存分析的排序最佳候选，不表示已经唯一解析。'
                '以下逐项列出全部保留的竞争候选与候选级证据。</div>'
            )
        return boundary + f'<div class="rc-cand-row">{cells}</div>'

    def _build_overlay_svg(
        q_peaks,
        lib_peaks,
        matched_frags=None,
        *,
        precursor_mz=None,
        precursor_charge=None,
        query_saved_count=None,
        w=1120,
        h=300,
    ):
        """Render a true shared-axis query/reference mirror spectrum.

        The previous plot independently stretched the query and library m/z
        ranges, so vertical alignment did not represent equal m/z.  Both sides
        now use one x-axis; the complete saved peak lists remain in the tables.
        """
        q_saved = [
            (float(p[0]), max(0.0, float(p[1]))) for p in (q_peaks or [])
            if isinstance(p, (list, tuple)) and len(p) >= 2
        ]
        l_saved = [
            (float(p[0]), max(0.0, float(p[1]))) for p in (lib_peaks or [])
            if isinstance(p, (list, tuple)) and len(p) >= 2
        ]
        if not q_saved and not l_saved:
            return ""

        # 单电荷 MS/MS 的高于前体窗口记录会把真正的碎片区压成一条窄线。
        # 镜像图只聚焦化学上可解释的产品离子窗口；所有保存峰（含超窗记录）
        # 仍逐项保留在后面的完整峰表，不在报告中删除。
        focus_limit = None
        try:
            parsed_precursor = float(precursor_mz)
            parsed_charge = abs(float(precursor_charge))
            if parsed_precursor > 0 and parsed_charge == 1:
                focus_limit = parsed_precursor + 1.5
        except (TypeError, ValueError):
            focus_limit = None
        q_valid = (
            [peak for peak in q_saved if peak[0] <= focus_limit]
            if focus_limit is not None else list(q_saved)
        )
        l_valid = (
            [peak for peak in l_saved if peak[0] <= focus_limit]
            if focus_limit is not None else list(l_saved)
        )
        if not q_valid and not l_valid:
            q_valid, l_valid = list(q_saved), list(l_saved)
            focus_limit = None
        q_outside = len(q_saved) - len(q_valid)
        l_outside = len(l_saved) - len(l_valid)

        matched_query = {
            round(float(item.get("query_mz", 0)), 5)
            for item in (matched_frags or [])
        }
        matched_library = {
            round(float(item.get("lib_mz", 0)), 5)
            for item in (matched_frags or [])
        }

        def _is_hit(mz, values):
            return any(abs(float(mz) - value) <= 0.02 for value in values)

        all_mz = [mz for mz, _ in q_valid + l_valid]
        mz_min = max(0.0, min(all_mz) - 2.0)
        mz_max = max(all_mz) + 2.0
        mz_range = max(mz_max - mz_min, 1.0)
        left, right = 54.0, 18.0
        top, bottom = 38.0, 40.0
        baseline = h / 2.0
        upper_height = baseline - top - 18.0
        lower_height = h - bottom - baseline - 18.0
        plot_width = w - left - right
        q_max = max((intensity for _, intensity in q_valid), default=1.0) or 1.0
        l_max = max((intensity for _, intensity in l_valid), default=1.0) or 1.0
        raw_query_count = int(query_saved_count or len(q_saved))
        raw_query_unplotted = max(raw_query_count - len(q_valid), 0)
        query_count_label = (
            f'质控后用于匹配 {len(q_valid)} 个峰；原始工件保存 {raw_query_count} 个峰'
        )
        library_count_label = (
            f'图中显示 {len(l_valid)} 个峰；完整峰表保存 {len(l_saved)} 个峰'
            if l_outside else f'图中与完整峰表均为 {len(l_saved)} 个峰'
        )
        query_color = "#007A78"
        library_color = "#6F4E9B"
        shared_color = "#C43D55"
        parts = [
            f'<rect x="0" y="0" width="{w}" height="{h}" rx="8" fill="#fbfcfe"/>',
            f'<text x="{left}" y="20" font-size="13" font-weight="700" fill="{query_color}">'
            f'样品 MS/MS（Query，向上）· {query_count_label}</text>',
            f'<text x="{left}" y="{h - 11}" font-size="13" font-weight="700" fill="{library_color}">'
            f'参考库 MS/MS（Library，向下）· {library_count_label}</text>',
            f'<line x1="{left}" y1="{baseline:.1f}" x2="{w-right}" y2="{baseline:.1f}" '
            'stroke="#64748b" stroke-width="1.2"/>',
        ]

        for tick_index in range(9):
            fraction = tick_index / 8.0
            mz_tick = mz_min + fraction * mz_range
            x_tick = left + fraction * plot_width
            parts.extend([
                f'<line x1="{x_tick:.1f}" y1="{top}" x2="{x_tick:.1f}" y2="{h-bottom}" '
                'stroke="#e7edf3" stroke-width="1"/>',
                f'<text x="{x_tick:.1f}" y="{baseline + 15:.1f}" font-size="10" '
                f'fill="#64748b" text-anchor="middle">{mz_tick:.1f}</text>',
            ])
        parts.extend([
            f'<text x="{left - 10}" y="{top + 2}" font-size="10" fill="#64748b" text-anchor="end">100%</text>',
            f'<text x="{left - 10}" y="{baseline - 5:.1f}" font-size="10" fill="#64748b" text-anchor="end">0</text>',
            f'<text x="{left - 10}" y="{h-bottom}" font-size="10" fill="#64748b" text-anchor="end">100%</text>',
            f'<text x="{w-right}" y="{baseline + 29:.1f}" font-size="11" fill="#475569" text-anchor="end">m/z（共享横轴）</text>',
        ])

        def _render_peaks(peaks, maximum, upward, base_color, hit_values):
            peak_parts = []
            top_indices = set(
                sorted(range(len(peaks)), key=lambda idx: -peaks[idx][1])[:8]
            )
            label_candidates = {
                index for index, (mz, _) in enumerate(peaks)
                if index in top_indices or _is_hit(mz, hit_values)
            }
            # 命中峰优先，再按强度选择；同一半谱中横向过近的标签只显示
            # 一个，精确的每一对 m/z 仍在紧随其后的碎片对表中。
            label_indices = set()
            occupied_x = []
            for index in sorted(
                label_candidates,
                key=lambda idx: (
                    0 if _is_hit(peaks[idx][0], hit_values) else 1,
                    -peaks[idx][1],
                ),
            ):
                x_candidate = left + (
                    peaks[index][0] - mz_min
                ) / mz_range * plot_width
                if all(abs(x_candidate - previous) >= 34 for previous in occupied_x):
                    label_indices.add(index)
                    occupied_x.append(x_candidate)
                if len(label_indices) >= 12:
                    break
            for index, (mz, intensity) in enumerate(peaks):
                x = left + (mz - mz_min) / mz_range * plot_width
                scaled = intensity / maximum
                y_tip = (
                    baseline - scaled * upper_height
                    if upward else baseline + scaled * lower_height
                )
                hit = _is_hit(mz, hit_values)
                color = base_color
                peak_parts.append(
                    f'<line x1="{x:.2f}" y1="{baseline:.1f}" x2="{x:.2f}" y2="{y_tip:.2f}" '
                    f'stroke="{color}" stroke-width="{1.9 if hit else 1.35}"/>'
                )
                if hit:
                    peak_parts.append(
                        f'<circle cx="{x:.2f}" cy="{y_tip:.2f}" r="2.7" '
                        f'fill="{shared_color}" stroke="#fff" stroke-width="0.8"/>'
                    )
                if index in label_indices:
                    label_y = y_tip - 5 if upward else y_tip + 11
                    label_y = max(38.0, min(h - 34.0, label_y))
                    peak_parts.append(
                        f'<text x="{x:.2f}" y="{label_y:.2f}" font-size="8.5" '
                        f'fill="{color}" font-weight="{700 if hit else 400}" '
                        f'text-anchor="middle">{mz:.3f}</text>'
                    )
            return "".join(peak_parts)

        parts.append(_render_peaks(q_valid, q_max, True, query_color, matched_query))
        parts.append(_render_peaks(l_valid, l_max, False, library_color, matched_library))
        parts.append(
            f'<text x="{w-right}" y="22" font-size="10.5" fill="#475569" text-anchor="end">'
            f'红色端点：样品谱与参考谱的 m/z 差不超过 0.020 Da（{len(matched_frags or [])} 对）；峰棒保留各自谱图颜色</text>'
        )
        focus_note = ""
        if focus_limit is not None and (q_outside or l_outside):
            focus_note = (
                f'为避免前体范围外记录压缩产品离子区，本图仅显示 m/z ≤ {focus_limit:.4f}；'
                f'未显示的样品谱 {q_outside} 个峰、参考谱 {l_outside} 个峰仍逐项列在完整峰表中。'
            )
        query_scope_note = (
            f'未显示的样品原始峰 {raw_query_unplotted} 个仍保存在完整审计峰表中。'
            if raw_query_unplotted else ''
        )
        return (
            '<div class="mirror-spectrum">'
            f'<svg data-spectrum-layout="mirror-shared-axis" role="img" '
            f'data-query-saved-peaks="{raw_query_count}" data-query-plotted-peaks="{len(q_valid)}" '
            f'aria-label="Query 与 Library MS/MS 共享 m/z 轴镜像对照" '
            f'viewBox="0 0 {w} {h}" preserveAspectRatio="xMidYMid meet">'
            f'{"".join(parts)}</svg>'
            '<div class="mirror-caption">上下两张谱使用同一 m/z 横轴；每张谱的峰强度分别以其最强峰归一化为 100%。'
            f'{focus_note}{query_scope_note}质控后匹配峰与原始保存峰采用不同数据范围并分别标注；所有原始 m/z 和强度见下方折叠审计峰表。</div></div>'
        )

    # ── Build HTML ──────────────────────────────────────────────────
    now_str = __import__("datetime").datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    h = f'''<!DOCTYPE html><html lang="zh-CN"><head><meta charset="UTF-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>LC-MS/MS Batch Analysis Report — MSsystem</title>
<style>
*{{box-sizing:border-box;margin:0;padding:0}}
body{{font-family:-apple-system,"Segoe UI",Arial,sans-serif;background:#e9eef4;color:#222;line-height:1.6;font-size:14px}}
.container{{width:calc(100% - 56px);max-width:1520px;margin:28px auto 48px;padding:20px 16px 28px;background:#f8fafc;border:1px solid #d8e1e9;border-radius:14px;box-shadow:0 8px 30px rgba(15,76,129,.08)}}
h1{{background:linear-gradient(135deg,#0f4c81,#0f4c81);color:white;padding:20px 28px;border-radius:10px;margin-bottom:20px;font-size:20px;display:flex;justify-content:space-between;align-items:center;flex-wrap:wrap;gap:10px}}
h1 .sub{{font-size:13px;color:rgba(255,255,255,0.8);font-weight:400}}
.stats-row{{display:grid;grid-template-columns:repeat(6,1fr);gap:12px;margin-bottom:16px}}
.stat-card{{background:white;padding:14px 16px;border-radius:10px;text-align:center;box-shadow:0 1px 4px rgba(0,0,0,0.08);border:1px solid #e0e0e0}}
.stat-card .val{{font-size:26px;font-weight:bold;color:#0f4c81;line-height:1.1}}
.stat-card .lbl{{font-size:12px;color:#888;margin-top:3px}}
.stat-card.hi{{background:linear-gradient(135deg,#0f4c81,#0f4c81);color:white}}
.stat-card.hi .val{{color:white}}
.stat-card.hi .lbl{{color:rgba(255,255,255,0.8)}}
.report-guide{{background:#fff;border:1px solid #cbd8e3;border-left:4px solid #0f4c81;border-radius:8px;padding:11px 14px;margin-bottom:14px;color:#334155;font-size:12.5px;line-height:1.65}}
.report-guide b{{color:#0f4c81}}
.high-gate-note{{margin-top:7px;padding:6px 8px;border:1px solid #f0b44c;border-left:3px solid #c47b00;border-radius:5px;background:#fff8e8;color:#7a4b00;font-size:10.5px;line-height:1.5}}
.high-gate-note b{{color:#8a5300}}

/* Toolbar */
.toolbar{{background:white;padding:10px 16px;border-radius:10px;margin-bottom:14px;box-shadow:0 1px 4px rgba(0,0,0,0.08);border:1px solid #e0e0e0;display:flex;gap:10px;flex-wrap:wrap;align-items:center;font-size:13px}}
.toolbar label{{font-weight:600;color:#555;white-space:nowrap}}
.toolbar select,.toolbar input{{padding:5px 8px;border:1px solid #ddd;border-radius:4px;font-size:13px;font-family:inherit;background:white}}
.toolbar select:focus,.toolbar input:focus{{outline:none;border-color:#0f4c81}}
.toolbar .result-count{{margin-left:auto;color:#888;font-size:13px}}

/* Table */
.table-container{{background:white;border-radius:10px;padding:0;box-shadow:0 2px 8px rgba(0,0,0,0.08);border:1px solid #e0e0e0;overflow:hidden}}
.table-scroll{{overflow-x:auto;max-height:600px;overflow-y:auto}}
.main-table{{width:100%;min-width:1180px;table-layout:fixed;border-collapse:collapse;font-size:13px}}
.main-table thead tr{{background:#0f4c81;color:white;position:sticky;top:0;z-index:5}}
.main-table th{{padding:10px 7px;text-align:center;font-size:12px;font-weight:600}}
.main-table td{{padding:10px 9px;border-bottom:1px solid #e5eaf0;vertical-align:top;overflow-wrap:anywhere}}
.candidate-cell{{font-size:13px;line-height:1.45}}
.candidate-name{{font-size:14px;font-weight:700;color:#153e5c;word-break:normal;overflow-wrap:anywhere}}
.candidate-meta{{display:flex;gap:5px;flex-wrap:wrap;margin-top:6px}}
.meta-chip{{display:inline-block;padding:1px 6px;border-radius:999px;background:#edf4f8;color:#35596b;font-size:10px;line-height:1.5}}
.locator-cell{{font-family:ui-monospace,SFMono-Regular,Menlo,monospace;font-size:11px;line-height:1.55;color:#334155}}
.locator-cell .loc-file{{font-family:-apple-system,"Segoe UI",Arial,sans-serif;font-weight:700;color:#0f4c81;margin-bottom:3px;overflow-wrap:anywhere}}
.locator-cell .loc-key{{color:#64748b;font-family:-apple-system,"Segoe UI",Arial,sans-serif}}
.evidence-cell{{font-size:11px;line-height:1.5}}
.evidence-grid{{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:5px}}
.evidence-item{{background:#f4f8fb;border:1px solid #e2e8f0;border-radius:5px;padding:4px 5px;text-align:center}}
.evidence-item b{{display:block;color:#0f4c81;font-family:ui-monospace,SFMono-Regular,Menlo,monospace;font-size:12px}}
.fragment-preview{{display:flex;gap:4px;flex-wrap:wrap;align-content:flex-start}}
.fragment-chip{{display:inline-block;background:#fff7ed;border:1px solid #fed7aa;border-radius:4px;color:#9a3412;padding:2px 5px;font:10.5px ui-monospace,SFMono-Regular,Menlo,monospace}}
.fragment-more{{font-size:10px;color:#64748b;padding:2px 0}}
.source-structure-cell{{font-size:11px;line-height:1.5;color:#475569}}
.source-structure-cell .mono{{display:block;margin-top:4px;font-size:9.5px;word-break:break-all;color:#64748b;max-height:48px;overflow:auto}}
tr.main-row{{cursor:pointer;transition:background 0.15s}}
tr.main-row:hover{{background:#e3f2fd!important}}
tr.main-row td .expand-icon{{display:inline-block;width:22px;height:22px;line-height:22px;text-align:center;background:#e8e8e8;border-radius:4px;font-size:14px;font-weight:bold;color:#555;margin-right:4px;transition:all 0.2s}}
tr.main-row:hover td .expand-icon{{background:#0f4c81;color:white}}
tr.main-row.expanded td .expand-icon{{background:#0f4c81;color:white}}
.frag-row{{display:none}}
.frag-row.visible{{display:table-row!important}}
.frag-row td{{padding:0!important;background:#f8f9fa}}
.frag-detail{{padding:14px 18px}}
.frag-grid-3{{display:grid;grid-template-columns:230px 1fr 1fr;gap:16px;margin-bottom:12px;align-items:start}}
.frag-peaks-grid{{display:grid;grid-template-columns:1fr 1fr;gap:16px;margin:14px 0}}
@media(max-width:800px){{.frag-grid-3{{grid-template-columns:1fr}} .frag-peaks-grid{{grid-template-columns:1fr}}}}
.frag-structure-box{{text-align:center;display:flex;flex-direction:column;align-items:center}}
.frag-structure-box img{{width:100%;max-width:220px;height:auto;border:1px solid #ddd;border-radius:4px;background:white;padding:4px;object-fit:contain}}
.frag-structure-box svg{{width:100%;max-width:220px;height:auto}}
.frag-info-box h4{{font-size:12px;font-weight:bold;color:#666;margin-bottom:4px;letter-spacing:0.5px;text-transform:uppercase}}
.frag-info-box .mono{{font-family:monospace;font-size:12px;word-break:break-all;background:white;padding:6px 8px;border-radius:4px;border:1px solid #ddd;line-height:1.4}}
.frag-match-title{{font-size:13px;font-weight:bold;color:#333;margin:10px 0 6px}}
.frag-table{{width:100%;border-collapse:collapse;font-size:12px;margin-top:4px}}
.frag-table th{{background:#0f4c81;color:white;padding:6px 5px;text-align:center;font-size:11px}}
.frag-table td{{padding:4px 5px;border-bottom:1px solid #e0e0e0;text-align:center;font-family:monospace;font-size:12px}}
.frag-table tr:hover{{background:#e3f2fd}}

/* Score colors */
.sc-green{{color:#1a5f2a;font-weight:bold}}
.sc-amber{{color:#8b6914;font-weight:bold}}
.sc-gray{{color:#666}}
.sc-ppm-good{{color:#2e7d32}}
.sc-ppm-warn{{color:#e65100}}
.sc-ppm-bad{{color:#c62828}}

/* Print */
.btn-print{{background:rgba(255,255,255,0.2);color:white;border:1px solid rgba(255,255,255,0.4);padding:6px 16px;border-radius:6px;cursor:pointer;font-size:13px;font-weight:600;font-family:inherit}}
.btn-print:hover{{background:rgba(255,255,255,0.3)}}
.no-match{{background:white;padding:20px;border-radius:10px;color:#999;text-align:center;border:1px dashed #ddd;margin-top:16px}}
@media print{{
  body{{background:white;font-size:11px}}
  .toolbar{{display:none}}
  .btn-print{{display:none}}
  .container{{width:100%;max-width:none;margin:0;padding:0.3in;border:0;border-radius:0;box-shadow:none;background:white}}
  .frag-row{{display:table-row!important}}
  tr.main-row{{cursor:default}}
  .table-scroll{{max-height:none;overflow:visible}}
}}
.pathway-card{{margin-top:6px;padding:6px 8px;background:#FFFAF0;border:1px solid #EDD6B0;border-radius:5px;font-size:10px}}
.pathway-title{{font-size:11px;font-weight:700;color:#7B341E;margin-bottom:3px}}
.pathway-section{{margin-top:4px}}
.pathway-section>b{{font-size:10px;color:#9C4221;display:inline-block;margin-bottom:2px}}
.cascade-sub{{font-size:9px;color:#A0AEC0;font-weight:400;margin-left:3px}}
.cascade-box{{border-radius:4px;padding:4px 6px;margin:3px 0;font-size:10px}}
.cascade-box>b{{font-size:10px;color:#276749}}
.cascade-sugar{{background:#F0FFF4;border:1px solid #9AE6B4}}
.cascade-ester{{background:#FFF5F5;border:1px solid #FEB2B2}}
.cascade-ester>b{{color:#9B2C2C}}
.cascade-chain{{margin:3px 0;line-height:1.8;word-break:break-word}}
.cascade-chain .node{{display:inline-block;background:white;border:1px solid #CBD5E0;border-radius:3px;padding:1px 4px;font-weight:600;font-size:9px;color:#2D3748}}
.cascade-chain .node.aglycone{{background:#FEFCBF;border-color:#ECC94B;color:#744210}}
.cascade-chain .node-tag{{display:block;font-size:7px;font-weight:400;color:#975A16;line-height:1}}
.cascade-chain .arrow{{color:#4A5568;font-size:9px;margin:0 1px}}
.cascade-meta{{font-size:9px;color:#4A5568;margin-top:1px}}
.subtype-badge{{display:inline-block;background:#9B2C2C;color:white;padding:1px 6px;border-radius:9999px;font-size:8px;font-weight:600;margin-left:3px}}
.nl-explained{{color:#276749;font-weight:600}}
.nl-unexplained{{color:#A0AEC0;font-style:italic}}
.spec{{display:inline-block;padding:1px 4px;border-radius:3px;font-size:8px;font-weight:600}}
.spec-high{{background:#C6F6D5;color:#22543D}}
.spec-med{{background:#FEEBC8;color:#7B341E}}
.fragment-table{{width:100%;font-size:9px;border-collapse:collapse;margin:3px 0}}
.fragment-table th{{background:#7B341E;color:white;padding:3px 5px;text-align:left;font-weight:600;font-size:8px}}
.fragment-table td{{padding:2px 5px;border-bottom:1px solid #F5E6D3}}
.fragment-table tr:nth-child(even){{background:#FFF8F0}}

/* ===== v5.1 redesigned expanded detail card (rc-*) — teal, 2-col + full-width ===== */
.rc-card{{padding:0!important;max-width:none;width:100%;margin:0 auto}}
.rc-hd{{display:flex;align-items:center;gap:12px;flex-wrap:wrap;padding:15px 22px;background:#0f4c81;color:#fff}}
.rc-num{{width:27px;height:27px;border-radius:7px;background:#ffffff26;color:#fff;font-weight:700;font-size:13px;display:flex;align-items:center;justify-content:center;flex-shrink:0}}
.rc-ttl{{font-size:18px;font-weight:700;color:#fff;overflow-wrap:anywhere}}
.rc-cbadge{{display:inline-block;color:#fff;background:#ffffff26;padding:2px 10px;border-radius:999px;font-size:10.5px;font-weight:600}}
.rc-hd-meta{{font-size:11.5px;color:#cfe2f3}}
.rc-badges{{margin-left:auto;display:flex;gap:6px;flex-wrap:wrap}}
.rc-bd{{font-size:11px;font-weight:700;padding:3px 9px;border-radius:5px}}
.rc-bd-cos{{background:#ffffff26;color:#fff}}.rc-bd-fin{{background:#1a8a5a;color:#fff}}.rc-bd-ppm{{background:#ffffff26;color:#fff}}
.rc-grid2{{display:grid;grid-template-columns:286px minmax(0,1fr);gap:0;border-bottom:1px solid #eef2f6;align-items:stretch}}
.rc-side{{padding:14px 16px;border-right:1px solid #eef2f6;background:#fcfdfe;text-align:center}}
.rc-struct{{width:240px;max-width:100%;height:176px;margin:4px auto 10px;background:#fff;border:1px solid #e2e8f0;border-radius:7px;display:flex;align-items:center;justify-content:center;padding:5px}}
.rc-struct img,.rc-struct svg{{max-width:100%;max-height:100%;object-fit:contain}}
.rc-struct>div{{max-width:100%}}
.rc-kv{{font-size:12.5px;color:#5d6d7e;margin:4px 0;display:flex;justify-content:space-between;gap:6px;text-align:left}}
.rc-kv>span{{flex:0 0 72px}}
.rc-kv b{{color:#0f4c81;font-weight:600;text-align:right;overflow-wrap:anywhere;min-width:0}}
.rc-main{{padding:14px 18px}}
.rc-scores{{display:grid;grid-template-columns:repeat(6,1fr);gap:8px;margin-bottom:12px}}
.rc-sc{{background:#f4f8fb;border:1px solid #e2e8f0;border-radius:6px;padding:9px 2px;text-align:center}}
.rc-sc-v{{font-size:17px;font-weight:700;color:#0f4c81}}.rc-sc-l{{font-size:10.5px;color:#7b8a99;margin-top:2px}}
.rc-sc-hi{{background:#e8f4f8}}.rc-sc-ex{{background:#e8f5ee;border-color:#b6dec4}}
.rc-main-content-grid{{display:grid;grid-template-columns:minmax(0,1.3fr) minmax(300px,.7fr);gap:14px;align-items:start}}
.rc-evidence-stack{{min-width:0}}
.rc-expert-panel{{min-width:0;background:#fbfcfe;border:1px solid #dce6ed;border-radius:7px;padding:10px 12px}}
.rc-expert-panel .rc-lbl{{margin-bottom:8px}}
.rc-source-box{{background:#f7fafc;border:1px solid #dce6ed;border-radius:7px;padding:9px 10px;text-align:left}}
.rc-source-grid{{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:7px}}
.rc-source-item{{min-width:0;background:#fff;border:1px solid #e4ebf0;border-radius:5px;padding:6px 7px}}
.rc-source-item span{{display:block;color:#718096;font-size:9px;text-transform:uppercase;letter-spacing:.25px;margin-bottom:1px}}
.rc-source-item b{{display:block;color:#234e6f;font-size:10.5px;overflow-wrap:anywhere}}
.rc-source-note{{margin-top:8px;padding-top:7px;border-top:1px dashed #cbd5e1;color:#526775;font-size:10.5px;line-height:1.5}}
.rc-smi{{font-size:10.5px;background:#f4f8fb;border:1px solid #e2e8f0;border-radius:5px;padding:7px 9px;line-height:1.45;max-height:90px;overflow:auto;text-align:left;margin:0 0 10px;word-break:break-all}}
.rc-spec{{width:100%;max-width:1080px;margin:0 auto;border:1px solid #e2e8f0;border-radius:7px;background:#fafbfc;text-align:center;padding:6px}}
.rc-spec svg{{display:block;width:100%;height:auto;margin:0 auto}}
.mirror-spectrum{{width:100%;border:1px solid #d9e3ec;border-radius:8px;background:#fbfcfe;overflow:hidden}}
.mirror-spectrum svg{{display:block;width:100%;min-height:0;height:auto}}
.mirror-caption{{padding:7px 10px;border-top:1px solid #e2e8f0;background:#f8fafc;color:#586a78;font-size:11px;text-align:left}}
.rc-lbl{{font-size:12.5px;font-weight:700;color:#0f4c81;text-transform:uppercase;letter-spacing:.3px;margin:0 0 10px;padding-bottom:5px;border-bottom:2px solid #0f4c8133;text-align:left}}
.rc-bars{{display:flex;flex-direction:column;gap:5px;text-align:left}}
.rc-bar-row{{display:flex;align-items:center;gap:6px}}
.rc-bar-nm{{font-size:10px;color:#5d6d7e;width:74px;flex-shrink:0;line-height:1.1}}
.rc-bar-track{{flex:1;height:8px;background:#eef2f6;border-radius:4px;overflow:hidden}}
.rc-bar-fill{{height:100%;border-radius:4px}}
.rc-bar-val{{font-size:10px;font-weight:600;color:#0f4c81;width:72px;text-align:right;flex-shrink:0}}
.rc-fw{{padding:14px 22px}}
.rc-fw+.rc-fw{{border-top:1px dashed #e8eef5}}
.rc-rsn-h{{font-size:11.5px;color:#0f4c81;font-weight:600;margin-bottom:3px}}
.rc-rsn-list{{list-style:none;padding-left:0;text-align:left}}
.rc-rsn-list li{{font-size:11.5px;color:#34495e;padding:2px 0 2px 13px;position:relative}}
.rc-rsn-list li:before{{content:'▸';color:#0f4c81;position:absolute;left:0}}
.rc-chem{{font-size:11px;color:#7c5a10;background:#fef9e7;border:1px solid #f5d98b;border-radius:5px;padding:5px 8px;line-height:1.5;margin-top:6px;text-align:left}}
.rc-frag-h{{font-size:11.5px;font-weight:600;color:#9c4221;margin:6px 0 3px}}
.rc-frag-list{{list-style:none;padding-left:0;text-align:left}}
.rc-frag-list li{{font-size:11.5px;color:#4a3500;padding:3px 0 3px 14px;position:relative;line-height:1.5}}
.rc-frag-list li:before{{content:'⚡';position:absolute;left:0;font-size:10px}}
.rc-nl-chips{{display:flex;flex-wrap:wrap;gap:5px;margin-top:5px}}
.rc-nl-chip{{font-size:10.5px;background:#fef5e7;border:1px solid #f5cba7;color:#9c4221;padding:2px 8px;border-radius:4px;font-family:monospace}}
.rc-muted{{font-size:11.5px;color:#94a3b8;font-style:italic}}
.rc-frag-viz{{width:100%;text-align:center;margin:6px auto 12px}}
.rc-frag-viz>div{{margin-left:auto!important;margin-right:auto!important;max-width:1040px}}
.rc-frag-viz>div>div[style*="display:flex"]{{justify-content:center}}
.rc-frag-viz>div>div[style*="max-width"]{{margin-left:auto!important;margin-right:auto!important}}
.rc-rdkit-box{{width:100%;max-width:520px;margin:8px auto;border:1px solid #e2e8f0;border-radius:7px;background:#fff;padding:6px}}
.rc-rdkit-box svg{{display:block;width:100%;height:auto}}
.rc-rdkit{{text-align:center}}
.rc-cand-row{{display:flex;gap:10px;flex-wrap:wrap;justify-content:center}}
.rc-cand{{width:190px;border:1px solid #e2e8f0;border-radius:7px;padding:7px;background:#fff;text-align:center}}
.rc-cand-win{{border-color:#86efac;background:#f0fdf4}}
.rc-cand-s{{width:145px;height:96px;margin:0 auto 4px;display:flex;align-items:center;justify-content:center}}
.rc-cand-s svg{{max-width:100%;max-height:100%}}
.rc-cand-n{{font-size:10px;color:#34495e;font-weight:600;line-height:1.35;min-height:26px;overflow-wrap:anywhere}}
.rc-cand-meta{{font-size:9.5px;color:#7b8a99;margin-top:3px}}
.rc-cand-meta b{{color:#0f4c81}}
.rc-pk-grid{{display:grid;grid-template-columns:1fr 1fr;gap:14px;margin-top:8px}}
.rc-pk-h{{font-size:11px;font-weight:600;color:#0f4c81;margin-bottom:4px}}
.rc-locator{{margin:0;padding:11px 13px;border:1px solid #9cc9c6;border-left:4px solid #176b68;border-radius:6px;background:#f0f8f7;color:#214c50;font-size:11px;line-height:1.55}}
.locator-grid{{display:grid;grid-template-columns:repeat(5,minmax(0,1fr));gap:8px}}
.locator-item{{background:#fff;border:1px solid #cce2df;border-radius:5px;padding:7px 8px;min-width:0}}
.locator-item span{{display:block;color:#668187;font-size:9.5px;margin-bottom:2px;text-transform:uppercase;letter-spacing:.2px}}
.locator-item b{{display:block;color:#184b50;font:11px ui-monospace,SFMono-Regular,Menlo,monospace;overflow-wrap:anywhere}}
.locator-help{{margin-top:7px;color:#58767a;font-size:10.5px}}
.diagnostic-audit{{border:1px solid #d7e1e9;border-radius:8px;background:#fbfdff;padding:12px}}
.diag-summary{{display:flex;gap:10px;justify-content:space-between;align-items:center;flex-wrap:wrap;color:#0f4c81;font-size:13px}}
.diag-sub{{font-size:10px;color:#64748b;font-weight:400;margin-left:7px}}
.diag-status{{display:inline-block;padding:3px 8px;border-radius:999px;font-size:10.5px;font-weight:700}}
.diag-ok{{background:#dcfce7;color:#166534}}.diag-miss{{background:#fff7ed;color:#9a3412}}.diag-na{{background:#eef2f6;color:#475569}}
.diag-note,.diag-empty{{margin-top:8px;padding:7px 9px;border-radius:5px;background:#f8fafc;color:#536572;font-size:11px;line-height:1.55}}
.diag-block{{margin-top:11px}}.diag-block-title{{font-size:11.5px;font-weight:700;color:#334155;margin-bottom:5px}}
.diag-table-wrap{{width:100%;overflow-x:auto;border:1px solid #e2e8f0;border-radius:6px}}
.diag-table{{width:100%;border-collapse:collapse;font-size:10.5px;background:#fff}}
.diag-table th{{background:#eaf1f6;color:#24465b;padding:6px 7px;text-align:left;white-space:nowrap}}
.diag-table td{{padding:5px 7px;border-top:1px solid #edf1f4;vertical-align:top}}
.diag-table tr:nth-child(even){{background:#fbfcfd}}
.diag-hit{{color:#15803d;font-weight:700}}.diag-not-seen{{color:#9a3412;font-weight:600}}
.diag-footnote{{margin-top:10px;border-top:1px dashed #cbd5e1;padding-top:8px;color:#64748b;font-size:10px;line-height:1.55}}
details.rc-big>summary{{font-size:12px;font-weight:600;color:#0f4c81;cursor:pointer;padding:5px 0;list-style:none}}
details.rc-big>summary::-webkit-details-marker{{display:none}}
details.rc-big>summary:before{{content:'▶ ';color:#0f4c81;font-size:10px}}
details.rc-big[open]>summary:before{{content:'▼ '}}
@media(max-width:1100px){{.stats-row{{grid-template-columns:repeat(3,1fr)}}.locator-grid{{grid-template-columns:repeat(2,minmax(0,1fr))}}}}
@media(max-width:900px){{.rc-grid2{{grid-template-columns:1fr}}.rc-side{{border-right:none;border-bottom:1px solid #eef2f6}}.rc-scores{{grid-template-columns:repeat(3,1fr)}}.rc-main-content-grid{{grid-template-columns:1fr}}.rc-source-grid{{grid-template-columns:1fr}}.rc-pk-grid{{grid-template-columns:1fr}}}}
@media(max-width:800px){{.container{{width:calc(100% - 24px);margin:12px auto 24px;padding:12px 10px 20px;border-radius:10px}}}}
</style></head><body>
<div class="container">
<h1>LC-MS/MS 分析证据报告
  <span class="sub">样品：{source_meta} &nbsp;|&nbsp; 生成时间：{now_str} &nbsp;|&nbsp; 获得谱库候选的 MS2：{matched_scan_count}/{total} ({rate:.1f}%) &nbsp;|&nbsp; 去重候选：{unique_candidate_count} &nbsp;|&nbsp; 数据库：MSsystem</span>
  <button class="btn-print" onclick="window.print()">🖨️ 打印 / 另存为 PDF</button>
</h1>

<div class="stats-row">
  <div class="stat-card hi"><div class="val">{unique_candidate_count}</div><div class="lbl">去重化合物候选</div></div>
  <div class="stat-card"><div class="val">{total}</div><div class="lbl">原始 MS2 谱图</div></div>
  <div class="stat-card"><div class="val">{matched_scan_count}</div><div class="lbl">有谱库候选的 MS2</div></div>
  <div class="stat-card"><div class="val">{rate:.1f}%</div><div class="lbl">获得谱库候选的 MS2 比例</div></div>
  <div class="stat-card" title="按实际可报告身份范围去重；单一二维连接结构与未能区分的同分子式异构体组分别计数，Medium 不进入。">
    <div class="val">{n_high_scope}</div>
    <div class="lbl">High 级鉴定组</div>
    <div style="font-size:10px;color:#718096;margin-top:3px">单一二维连接：{n_high_exact} · 异构体组：{n_high_isomer}</div>
  </div>
  <div class="stat-card"><div class="val">{unmatched_scan_count}</div><div class="lbl">未找到谱库候选的 MS2</div></div>
</div>

<div class="report-guide"><b>阅读方式：</b>主表每行是一个去重化合物候选，候选名称来自最佳谱库记录，不代表唯一结构已确证。点击任意行可展开全宽证据卡：先按源文件 / RT / 前体 m/z / MS2 native ID 回查原谱，再核对共享 m/z 轴镜像谱、全部碎片对、诊断离子、经验中性丢失序列和馆藏碎裂签名。</div>

<script>
/* 顶部早加载：大文档(30MB+)解析到底部主脚本前，输入框已可安全触发——
   applyFilters 未就绪时轮询补执行；150ms 防抖。 */
function debounceFilter() {{
  if (window.__filterTimer) clearTimeout(window.__filterTimer);
  window.__filterTimer = setTimeout(function() {{
    if (typeof applyFilters === "function") {{ applyFilters(); }}
    else if (!window.__filterQueued) {{
      window.__filterQueued = true;
      var t0 = Date.now();
      var iv = setInterval(function() {{
        if (typeof applyFilters === "function") {{
          clearInterval(iv); applyFilters();
          var rc = document.getElementById("resultCount");
          if (rc) rc.textContent += "（页面仍在加载，已补执行搜索）";
        }} else if (Date.now() - t0 > 60000) {{ clearInterval(iv); }}
      }}, 800);
    }}
  }}, 150);
}}
</script>
<div class="toolbar" style="position:sticky;top:0;z-index:50;background:#fff;box-shadow:0 2px 8px rgba(0,0,0,0.12);padding:10px 12px;border-radius:8px;flex-wrap:wrap;gap:8px">
  <label>排序：</label>
  <select id="sortSelect" onchange="applyFilters()">
    <option value="final_desc">综合评分 ↓</option>
    <option value="cosine_desc">谱图相似度 ↓</option>
    <option value="ppm_asc">质量偏差 ppm ↑</option>
    <option value="frags_desc">相符碎片对 ↓</option>
    <option value="name_asc">化合物名称 A–Z</option>
  </select>
  <label style="margin-left:8px">化合物：</label>
  <input type="text" id="nameFilter" placeholder="名称 / 中文名 / m/z 搜索（含竞争候选）…" oninput="debounceFilter()" style="width:260px;font-size:14px;padding:6px 10px;border:2px solid #2C7A7B;border-radius:6px">
  <label style="margin-left:8px">置信度：</label>
  <select id="confFilter" onchange="applyFilters()">
    <option value="all">全部</option><option value="High">High</option><option value="Medium">Medium</option><option value="Low">Low</option>
  </select>
  <span class="result-count" id="resultCount"></span>
</div>

<div class="table-container">
<div class="table-scroll">
<table class="main-table" id="mainTable">
<colgroup><col style="width:44px"><col style="width:27%"><col style="width:23%"><col style="width:18%"><col style="width:18%"><col style="width:14%"></colgroup>
<thead><tr>
  <th>#</th>
  <th>化合物候选名称</th>
  <th>原始 MS2 定位</th>
  <th>样品谱与参考谱比较结果</th>
  <th>相符碎片 m/z（预览）</th>
  <th>参考库与结构</th>
</tr></thead>
<tbody id="tableBody">
'''

    # ── Build rows ──────────────────────────────────────────────────
    row_idx = 0
    for i, r in enumerate(matched):
        m = r["matches"][0]
        smi = str(m.get("smiles", "") or "")
        img_b64 = smiles_to_image_base64(smi, size=(220, 170)) if smi else None
        # Fall back to SVG when Cairo/PNG not available
        svg_str = None
        if smi and not img_b64:
            svg_str = smiles_to_svg(smi, size=(220, 170))
            if svg_str:
                # Strip XML declaration for inline HTML embedding
                svg_str = svg_str.replace('<?xml version=\'1.0\' encoding=\'iso-8859-1\'?>\n', '')
                # Set responsive width for the report
                svg_str = svg_str.replace('width="400"', 'width="220"', 1)
                svg_str = svg_str.replace('height="300"', 'height="165"', 1)
        # 同时保留两个口径：query_peaks 是生产检索视图；raw_peaks
        # 是结果工件保存的原始查询峰视图（历史上限 100）。报告只跳过
        # 不能解析为数值对的损坏记录，不再按前体质量或 Top-N 二次删峰。
        def _report_peak_pairs(values):
            pairs = []
            for peak in values or []:
                try:
                    mz_value = float(peak[0])
                    intensity_value = float(peak[1])
                except (TypeError, ValueError, IndexError, KeyError):
                    continue
                if not (mz_value == mz_value and intensity_value == intensity_value):
                    continue
                pairs.append((mz_value, intensity_value))
            return pairs

        qpeaks_full = _report_peak_pairs(
            r.get("raw_peaks") or r.get("query_peaks", []) or r.get("peaks", [])
        )
        qpeaks = _report_peak_pairs(
            r.get("query_peaks", []) or r.get("peaks", [])
        )
        lpeaks = _report_peak_pairs(m.get("lib_peaks", []))
        try:
            source_peak_count = int(r.get("peaks_count") or 0)
        except (TypeError, ValueError):
            source_peak_count = 0
        raw_complete_value = r.get("raw_peaks_complete")
        if raw_complete_value is None:
            raw_artifact_complete = bool(qpeaks_full) and (
                source_peak_count <= 0 or len(qpeaks_full) >= source_peak_count
            )
        else:
            raw_artifact_complete = bool(raw_complete_value)
        if raw_artifact_complete:
            raw_storage_status = "全部有效数值峰已保存"
            raw_peak_storage_note = (
                "此任务工件保存的原始查询峰已全部展示，报告层未做 Top-N 截断。"
            )
        elif source_peak_count > len(qpeaks_full):
            raw_storage_status = "历史工件存在保存上限"
            raw_peak_storage_note = (
                f"源记录为 {source_peak_count} 峰，当前历史任务工件仅保存 "
                f"{len(qpeaks_full)} 峰；报告已全部展示这些保存峰，但不能伪造未写入工件的"
                "其余峰。V11 新任务现按全部有效数值峰保存。"
            )
        else:
            raw_storage_status = "保存完整性未声明"
            raw_peak_storage_note = (
                "报告已展示任务工件中的全部保存峰，未做 Top-N 截断；"
                "该历史工件没有保存完整性声明。"
            )
        has_lib = bool(lpeaks)
        scoring_matched_frags = (
            _find_matched_frags(qpeaks, lpeaks, tol=0.02)
            if has_lib else []
        )
        review_matched_frags = (
            _find_matched_frags(qpeaks_full, lpeaks, tol=0.02)
            if has_lib else []
        )
        match_count = len(scoring_matched_frags)
        review_match_count = len(review_matched_frags)

        fs = m.get("final_score", 0)
        cs = m.get("cosine", 0)
        full_identity_conf = exact_structure_confidence(m)
        report_conf = reportable_identity_confidence(m)
        report_conf_label = reportable_identity_label(m, language="zh")
        report_scope_label = reportable_identity_scope_label(m, language="zh")
        report_scope = identification_scope(m)
        full_identity_chip = (
            f'<span class="meta-chip">完整身份 {full_identity_conf}</span>'
            if report_scope in {"exact_2d_connectivity", "isomer_group"}
            or report_conf != full_identity_conf
            else ""
        )
        final_claim_badge = final_claim_badge_html(m, language="zh")
        gate_required = bool(m.get("precursor_high_evidence_gate_required"))
        gate_reason_raw = str(
            m.get("precursor_high_evidence_gate_reason_zh")
            or m.get("precursor_high_evidence_gate_reason")
            or ""
        ).strip()
        if not gate_reason_raw and "precursor_high_evidence_gate_required" not in m:
            try:
                gate_annotation = _honest_annotate(m, m.get("top_candidates"))
                gate_required = bool(
                    gate_annotation.get("precursor_high_evidence_gate_required")
                )
                gate_reason_raw = str(
                    gate_annotation.get("precursor_high_evidence_gate_reason_zh")
                    or gate_annotation.get("precursor_high_evidence_gate_reason")
                    or ""
                ).strip()
            except Exception:
                gate_required = False
        high_gate_note = (
            '<div class="high-gate-note"><b>High 证据门控：</b>'
            '该候选已降为 Medium。'
            f'{_html_escape(gate_reason_raw)}</div>'
            if gate_required else ""
        )
        ppm_val = m.get("mass_ppm")
        try:
            ppm_number = float(ppm_val) if ppm_val is not None else None
        except (TypeError, ValueError):
            ppm_number = None
        ppm_str = f"{ppm_number:.1f}" if ppm_number is not None else "—"
        ppm_abs = abs(ppm_number) if ppm_number is not None else None
        ppm_cls = "sc-ppm-good" if (ppm_abs is not None and ppm_abs <= 5) else ("sc-ppm-warn" if (ppm_abs is not None and ppm_abs <= 20) else "sc-ppm-bad")
        compound_name_raw = (
            _standard_compound_name(m.get("compound_name", ""), m)
            or "Unnamed compound"
        )
        compound_name = _html_escape(compound_name_raw)
        class_evidence = classification_evidence_status(m)
        classification_raw = str(class_evidence["label"])
        if class_evidence["status"] == "tentative":
            classification_raw = f"待核查类别线索：{classification_raw}"
        classification = _html_escape(classification_raw)
        dataset_source_raw = str(m.get("dataset_source", "") or "")
        dataset_source = _html_escape(dataset_source_raw or "—")
        formula_raw = str(m.get("formula", "—") or "—")
        formula_str = _html_escape(formula_raw)
        resolved_adduct, resolved_adduct_source, _ = _resolved_adduct(m, r)
        adduct_raw = str(resolved_adduct or "Unresolved")
        adduct_html = _html_escape(adduct_raw)
        adduct_source_html = _html_escape(str(resolved_adduct_source or "unresolved"))
        library_adduct_raw = str(m.get("adduct", "") or "")
        library_adduct_html = _html_escape(library_adduct_raw or "—")
        candidate_adduct_raw = str(
            m.get("candidate_adduct_hypothesis") or ""
        )
        candidate_adduct_html = _html_escape(candidate_adduct_raw or "—")
        inchikey_raw = str(m.get("inchikey", "—") or "—")
        inchikey_html = _html_escape(inchikey_raw)
        source_compound_name_raw = str(
            m.get("source_compound_name") or m.get("compound_name") or ""
        )
        source_compound_name_html = _html_escape(source_compound_name_raw)
        lib_spectrum_id = _html_escape(str(m.get("lib_spectrum_id") or "—"))
        lib_accession = _html_escape(str(m.get("lib_source_accession") or "—"))
        source_tier = _html_escape(str(m.get("source_tier") or "unknown"))
        expert_score = float(m.get("expert_score", 0) or 0)
        expert_dimensions = (
            (m.get("expert_details") or {}).get("dimension_assessments") or {}
        )
        expert_dimension_total = len(expert_dimensions)
        expert_dimension_scored = sum(
            str(value.get("status") or "") == "scored"
            for value in expert_dimensions.values()
            if isinstance(value, dict)
        )
        expert_dimension_no_evidence = sum(
            str(value.get("status") or "") == "no_evidence"
            for value in expert_dimensions.values()
            if isinstance(value, dict)
        )
        expert_dimension_na = sum(
            str(value.get("status") or "") in {"not_applicable", "unavailable"}
            for value in expert_dimensions.values()
            if isinstance(value, dict)
        )
        locator = raw_spectrum_locator(r, source_name)
        rt_minutes = locator["retention_time_min"]
        rt_display = f"{rt_minutes:.3f} min" if rt_minutes is not None else "未记录"
        rt_clock = locator["retention_time_clock"] or "—"
        selected_ion_mz = locator.get("selected_ion_mz")
        search_precursor_mz = locator.get("search_precursor_mz")
        lcms_feature_used_for_search = bool(
            locator.get("lcms_feature_used_for_search")
        )
        selected_ion_mz_text = (
            f"{selected_ion_mz:.5f}"
            if selected_ion_mz is not None
            else "未记录"
        )
        search_precursor_mz_text = (
            f"{search_precursor_mz:.5f}"
            if search_precursor_mz is not None
            else "未记录"
        )
        if lcms_feature_used_for_search:
            precursor_header_text = (
                f"MS2 隔离目标 {selected_ion_mz_text} · "
                f"检索用实测 MS1 {search_precursor_mz_text}"
            )
            precursor_locator_rows = (
                '<div><span class="loc-key">MS2 隔离目标</span> '
                f'm/z {selected_ion_mz_text}</div>'
                '<div><span class="loc-key">检索用实测 MS1</span> '
                f'm/z {search_precursor_mz_text}</div>'
            )
            search_mass_grid_item = (
                '<div class="locator-item"><span>检索/计算质量（实测 MS1 特征质心）</span>'
                f'<b>m/z {search_precursor_mz_text}</b></div>'
            )
        else:
            precursor_header_text = f"前体 m/z {selected_ion_mz_text}"
            precursor_locator_rows = (
                '<div><span class="loc-key">前体 m/z</span> '
                f'{selected_ion_mz_text}</div>'
            )
            search_mass_grid_item = (
                '<div class="locator-item"><span>检索/计算质量</span>'
                f'<b>m/z {search_precursor_mz_text}（使用 MS2 隔离目标）</b></div>'
            )
        spectrum_id = _html_escape(locator["spectrum_id"] or "未记录")
        source_file_html = _html_escape(locator["source_file"] or "未记录")
        precursor_scan_html = _html_escape(locator["precursor_scan_ref"] or "—")
        ce_value = locator.get("collision_energy")
        ce_text = (
            f"{ce_value:.1f} {_html_escape(locator.get('collision_energy_unit') or '(单位未记录)')}"
            if ce_value is not None else "未记录"
        )
        ion_mode_html = _html_escape(str(r.get("ion_mode") or "未记录"))
        charge_html = _html_escape(str(r.get("precursor_charge") or "未记录"))
        tic_value = r.get("tic")
        base_peak_value = r.get("base_peak")
        try:
            tic_text = f"{float(tic_value):.6g}"
        except (TypeError, ValueError):
            tic_text = "未记录"
        try:
            base_peak_text = f"{float(base_peak_value):.6g}"
        except (TypeError, ValueError):
            base_peak_text = "未记录"
        locator_text = _html_escape(
            raw_spectrum_locator_text(r, source_name, language="zh")
        )

        # Structure image — prefer PNG, fall back to inline SVG, else placeholder
        if img_b64:
            struct_html = (
                f'<img src="data:image/png;base64,{img_b64}" style="max-width:100%;height:auto;object-fit:contain" alt="structure"/>'
            )
        elif svg_str:
            struct_html = (
                f'<div style="width:100%;max-width:220px">{svg_str}</div>'
            )
        else:
            struct_html = '<div style="width:100%;aspect-ratio:4/3;display:flex;align-items:center;justify-content:center;background:#f5f5f5;border:1px solid #ddd;border-radius:4px;color:#999;font-size:13px">No Structure</div>'

        # Fragment ion comparison table
        frag_comp_html = ""
        if scoring_matched_frags:
            frag_comp_html = '<div class="diag-table-wrap"><table class="diag-table"><thead><tr><th>#</th><th>Query m/z</th><th>Query Int</th><th>Library m/z</th><th>Library Int</th><th>Δm/z</th></tr></thead><tbody>'
            for jj, mf in enumerate(scoring_matched_frags):
                frag_comp_html += f'<tr><td>{jj+1}</td><td>{mf["query_mz"]:.4f}</td><td>{mf["query_intensity"]:.1f}</td><td>{mf["lib_mz"]:.4f}</td><td>{mf["lib_intensity"]:.1f}</td><td>{mf["diff"]:.4f}</td></tr>'
            frag_comp_html += '</tbody></table></div>'
        elif has_lib:
            frag_comp_html = '<div class="diag-empty">用于评分的质控后查询峰与库谱在 ±0.020 Da 内没有共同碎片对。</div>'
        else:
            frag_comp_html = '<div class="diag-empty">该候选没有可用的参考库谱峰表。</div>'

        # The primary mirror shows the exact processed query spectrum used for
        # scoring.  The complete raw peak table remains available below as an
        # audit artifact, but thousands of noise/profile points no longer
        # overwhelm the analytical comparison.
        svg_html = _build_overlay_svg(
            qpeaks,
            lpeaks,
            scoring_matched_frags,
            precursor_mz=search_precursor_mz,
            precursor_charge=r.get("precursor_charge"),
            query_saved_count=len(qpeaks_full),
        ) if (qpeaks or lpeaks) else ""

        matched_mz_values = sorted({
            round(float(item.get("query_mz", 0)), 4)
            for item in scoring_matched_frags
        })
        msms_display = "; ".join(f"{value:.4f}" for value in matched_mz_values)
        preview_values = matched_mz_values[:8]
        msms_preview_html = "".join(
            f'<span class="fragment-chip">{value:.4f}</span>'
            for value in preview_values
        )
        if len(matched_mz_values) > len(preview_values):
            msms_preview_html += (
                f'<span class="fragment-more">+{len(matched_mz_values) - len(preview_values)} 个，展开查看全部</span>'
            )
        if not msms_preview_html:
            msms_preview_html = '<span class="fragment-more">无共同碎片</span>'

        # Query peak table: every peak retained in the saved raw-query view.
        query_peaks_html = ""
        if qpeaks_full:
            qp_sorted = sorted(qpeaks_full, key=lambda p: float(p[0]))
            query_peaks_html = '<div class="diag-table-wrap"><table class="diag-table"><thead><tr><th>#</th><th>m/z</th><th>Intensity</th><th>Rel. Int. %</th><th>库谱共同峰</th></tr></thead><tbody>'
            max_int = max((float(p[1]) for p in qp_sorted), default=1) or 1
            for jj, (mz_val, inten) in enumerate(qp_sorted):
                rel = inten / max_int * 100 if max_int > 0 else 0
                is_shared = any(
                    abs(float(mz_val) - float(item.get("query_mz", 0))) <= 0.02
                    for item in review_matched_frags
                )
                shared_mark = '<span class="diag-hit">✓</span>' if is_shared else "—"
                query_peaks_html += (
                    f'<tr><td>{jj+1}</td><td class="mono">{mz_val:.5f}</td>'
                    f'<td>{inten:.3f}</td><td>{rel:.2f}%</td>'
                    f'<td>{shared_mark}</td></tr>'
                )
            query_peaks_html += '</tbody></table></div>'
        else:
            query_peaks_html = '<div class="diag-empty">无可用的查询谱峰数据。</div>'

        # Library peak table: every library peak saved with the match.
        lib_peaks_html = ""
        if lpeaks:
            lp_sorted = sorted(lpeaks, key=lambda p: float(p[0]))
            lib_peaks_html = '<div class="diag-table-wrap"><table class="diag-table"><thead><tr><th>#</th><th>m/z</th><th>Intensity</th><th>Rel. Int. %</th><th>查询谱共同峰</th></tr></thead><tbody>'
            max_int_l = max((float(p[1]) for p in lp_sorted), default=1) or 1
            for jj, (mz_val, inten) in enumerate(lp_sorted):
                rel = inten / max_int_l * 100 if max_int_l > 0 else 0
                is_shared = any(
                    abs(float(mz_val) - float(item.get("lib_mz", 0))) <= 0.02
                    for item in review_matched_frags
                )
                shared_mark = '<span class="diag-hit">✓</span>' if is_shared else "—"
                lib_peaks_html += (
                    f'<tr><td>{jj+1}</td><td class="mono">{mz_val:.5f}</td>'
                    f'<td>{inten:.3f}</td><td>{rel:.2f}%</td>'
                    f'<td>{shared_mark}</td></tr>'
                )
            lib_peaks_html += '</tbody></table></div>'
        else:
            lib_peaks_html = '<div class="diag-empty">无可用的参考库谱峰数据。</div>'

        diagnostic_html = render_diagnostic_ion_audit(
            m, query_peaks=qpeaks, tolerance_da=0.02
        )

        original_name_html = ""
        if source_compound_name_raw and source_compound_name_raw != compound_name_raw:
            original_name_html = (
                '<div class="rc-lbl" style="margin-top:10px">谱库原始名称 Library Record Name</div>'
                f'<div class="rc-smi">{source_compound_name_html}</div>'
            )

        row_id = f"r{i}"
        # 通用搜索索引（展示层辅助）：top-1 名称 + 竞争候选名称 + 中文名别名 + 前体 m/z。
        # 不改变打分/计数/治理判定；data-name 保持原语义供排序兼容。
        _search_terms = [compound_name_raw.lower()]
        _zh_chips = []
        for _zh in _zh_alias_for(compound_name_raw):
            _search_terms.append(_zh)
            _zh_chips.append(_zh)
        _tc = m.get("top_candidates") or []
        for _c in _tc[:20]:
            _cn = str(_c.get("compound_name") or "").strip().lower()
            if _cn:
                _search_terms.append(_cn)
                for _zh in _zh_alias_for(_cn):
                    _search_terms.append(_zh)
        # 主名被规范为系统名时，中文名 chip 回退用 top-1 候选俗名的别名
        if not _zh_chips and _tc:
            _zh_chips = _zh_alias_for(str(_tc[0].get("compound_name") or ""))
        _search_text = " ".join(dict.fromkeys(t for t in _search_terms if t))
        try:
            _mz_val = float(search_precursor_mz) if search_precursor_mz is not None else 0.0
        except (TypeError, ValueError):
            _mz_val = 0.0
        _zh_chip_html = (
            f'<span class="meta-chip zh-alias">{_html_escape("·".join(_zh_chips))}</span>'
            if _zh_chips else ""
        )
        _zh_title_suffix = (
            f' <span style="color:#0e7490;font-size:12px">{_html_escape("·".join(_zh_chips))}</span>'
            if _zh_chips else ""
        )
        h += f'''
<tr class="main-row" data-idx="{i}" data-cosine="{cs:.4f}" data-final="{fs:.4f}"
    data-ppm="{ppm_abs if ppm_abs is not None else 999:.2f}" data-frags="{match_count}"
    data-name="{_html_escape(compound_name_raw.lower(), quote=True)}" data-class="{_html_escape(classification_raw, quote=True)}" data-conf="{_html_escape(report_conf, quote=True)}" data-full-conf="{_html_escape(full_identity_conf, quote=True)}"
    data-search="{_html_escape(_search_text, quote=True)}" data-mz="{_mz_val:.5f}"
    onclick="toggleFragRow('{row_id}', this)">
  <td style="text-align:center"><span class="expand-icon" id="icon_{row_id}">➕</span></td>
  <td class="candidate-cell"><div class="candidate-name">{compound_name}</div>
    <div class="candidate-meta"><span class="meta-chip">{classification}</span><span class="meta-chip mono">{formula_str}</span><span class="meta-chip">{_html_escape(report_conf_label)}</span>{_zh_chip_html}{full_identity_chip}</div>
    <div style="display:flex;gap:4px;flex-wrap:wrap;margin-top:6px">{final_claim_badge}</div></td>
  <td class="locator-cell"><div class="loc-file">{source_file_html}</div>
    <div><span class="loc-key">RT</span> {rt_display} ({rt_clock})</div>
    {precursor_locator_rows}
    <div><span class="loc-key">MS2 native ID</span> {spectrum_id}</div>
    <div><span class="loc-key">前级扫描</span> {precursor_scan_html}</div>
    <div><span class="loc-key">峰数</span> 源记录 {r.get("peaks_count",0)} · 报告保存视图 {len(qpeaks_full)}</div></td>
  <td class="evidence-cell"><div class="evidence-grid">
    <div class="evidence-item"><b>{cs:.3f}</b>Cosine</div>
    <div class="evidence-item"><b>{fs:.3f}</b>Final</div>
    <div class="evidence-item"><b class="{ppm_cls}">{ppm_str}</b>ppm</div>
    <div class="evidence-item"><b>{match_count}</b>分析评分使用的相符碎片对</div>
  </div><div style="margin-top:6px;color:#64748b">按完整保存峰表复核的相符碎片对：<b>{review_match_count}</b></div>{high_gate_note}</td>
  <td><div class="fragment-preview" title="{_html_escape(msms_display, quote=True)}">{msms_preview_html}</div></td>
  <td class="source-structure-cell"><b>{dataset_source}</b><br>证据层：{source_tier}<br>库谱：<span class="mono">{lib_spectrum_id}</span><br>样品离子形式：{adduct_html}<small>（{adduct_source_html}）</small><br>谱库加合物：{library_adduct_html}{f'<br>候选离子假设：{candidate_adduct_html}<small>（非样品独立证据）</small>' if candidate_adduct_raw else ''}<br><span class="mono">SMILES: {_html_escape(smi)}</span></td>
</tr>
<tr class="frag-row" id="frag_{row_id}">
  <td colspan="6" style="padding:0">
    <div class="frag-detail rc-card">
      <div class="rc-hd">
        <span class="rc-num">{i+1}</span>
        <div style="display:flex;flex-direction:column;gap:3px">
          <span class="rc-ttl">{compound_name}{_zh_title_suffix}</span>
          <span class="rc-hd-meta">{formula_str} · {adduct_html} · {classification} · RT {rt_display} ({rt_clock}) · MS2 {spectrum_id} · {precursor_header_text} · PPM {ppm_str}</span>
          <div style="display:flex;gap:4px;margin-top:3px;flex-wrap:wrap">{_conf_badge(report_scope_label, report_conf)}{_conf_badge("完整身份", full_identity_conf) if report_scope in {"exact_2d_connectivity", "isomer_group"} or report_conf != full_identity_conf else ""}{final_claim_badge}</div>
          {high_gate_note}
        </div>
      </div>
      <div class="rc-fw">
        <div class="rc-lbl">原始 MS2 定位信息 Raw Spectrum Locator</div>
        <div class="rc-locator">
          <div class="locator-grid">
            <div class="locator-item"><span>源文件</span><b>{source_file_html}</b></div>
            <div class="locator-item"><span>MS2 谱图编号 / native ID</span><b>{spectrum_id}</b></div>
            <div class="locator-item"><span>保留时间</span><b>{rt_display} ({rt_clock})</b></div>
            <div class="locator-item"><span>MS2 隔离目标</span><b>m/z {selected_ion_mz_text}</b></div>
            {search_mass_grid_item}
            <div class="locator-item"><span>关联前级扫描</span><b>{precursor_scan_html}</b></div>
            <div class="locator-item"><span>碰撞能</span><b>{ce_text}</b></div>
            <div class="locator-item"><span>离子模式</span><b>{ion_mode_html}</b></div>
            <div class="locator-item"><span>前体电荷</span><b>{charge_html}</b></div>
            <div class="locator-item"><span>TIC</span><b>{tic_text}</b></div>
            <div class="locator-item"><span>基峰强度</span><b>{base_peak_text}</b></div>
            <div class="locator-item"><span>源记录峰数</span><b>{r.get("peaks_count",0)}</b></div>
            <div class="locator-item"><span>报告保存的全部样品峰</span><b>{len(qpeaks_full)}</b></div>
            <div class="locator-item"><span>原始峰保存策略</span><b>{raw_storage_status}</b></div>
            <div class="locator-item"><span>分析评分使用的样品峰</span><b>{len(qpeaks)}</b></div>
            <div class="locator-item"><span>参考库谱峰</span><b>{len(lpeaks)}</b></div>
          </div>
          <div class="locator-help">{locator_text}<br>优先按谱图编号/native ID 回查；查看软件不支持编号检索时，再联合 RT 与前体 m/z，并以下方查询碎片峰核对。</div>
        </div>
      </div>
      <div class="rc-grid2">
        <div class="rc-side">
          <div class="rc-lbl">结构式 Structure</div>
          <div class="rc-struct">{struct_html}</div>
          <div class="rc-kv"><span>Formula</span><b>{formula_str}</b></div>
          <div class="rc-kv"><span>InChIKey</span><b class="mono">{inchikey_html}</b></div>
          <div class="rc-kv"><span>Sample ion form</span><b>{adduct_html}</b></div>
          <div class="rc-kv"><span>Ion-form source</span><b>{adduct_source_html}</b></div>
          <div class="rc-kv"><span>Library adduct</span><b>{library_adduct_html}</b></div>
          {f'<div class="rc-kv"><span>Candidate ion-form hypothesis</span><b>{candidate_adduct_html} · candidate-dependent</b></div>' if candidate_adduct_raw else ''}
        </div>
        <div class="rc-main">
          <div class="rc-lbl">打分明细 Scores</div>
          <div class="rc-scores">
            <div class="rc-sc rc-sc-hi"><div class="rc-sc-v">{cs:.3f}</div><div class="rc-sc-l">Cosine</div></div>
            <div class="rc-sc rc-sc-hi"><div class="rc-sc-v">{fs:.3f}</div><div class="rc-sc-l">Final</div></div>
            <div class="rc-sc"><div class="rc-sc-v">{m.get("mass_score",0):.3f}</div><div class="rc-sc-l">Mass</div></div>
            <div class="rc-sc"><div class="rc-sc-v">{m.get("fragment_score",0):.3f}</div><div class="rc-sc-l">Frag</div></div>
            <div class="rc-sc"><div class="rc-sc-v">{m.get("intensity_score",0):.3f}</div><div class="rc-sc-l">Intensity</div></div>
            <div class="rc-sc rc-sc-ex"><div class="rc-sc-v">{expert_dimension_scored}/{expert_dimension_total or 8}</div><div class="rc-sc-l">Rule-evidence dimensions with support</div></div>
          </div>
          <div class="rc-main-content-grid">
            <div class="rc-evidence-stack">
              <div class="rc-lbl">SMILES</div>
              <div class="rc-smi mono">{_html_escape(smi) if smi else "—"}</div>
              {original_name_html}
              <div class="rc-lbl" style="margin-top:8px">谱库来源与证据边界</div>
              <div class="rc-source-box">
                <div class="rc-source-grid">
                  <div class="rc-source-item"><span>Dataset</span><b>{dataset_source}</b></div>
                  <div class="rc-source-item"><span>Source tier</span><b>{source_tier}</b></div>
                  <div class="rc-source-item"><span>Library spectrum</span><b class="mono">{lib_spectrum_id}</b></div>
                  <div class="rc-source-item"><span>Accession</span><b class="mono">{lib_accession}</b></div>
                  <div class="rc-source-item"><span>Sample ion form</span><b>{adduct_html} · {adduct_source_html}</b></div>
                  <div class="rc-source-item"><span>Library adduct</span><b>{library_adduct_html}</b></div>
                  <div class="rc-source-item"><span>Candidate class</span><b>{classification}</b></div>
                </div>
                <div class="rc-source-note">候选名称是当前保存分析中的最佳谱库候选，不代表唯一结构已确证；请联合镜像谱、共同碎片、诊断离子和鉴定结论标签解读。</div>
              </div>
            </div>
            <div class="rc-expert-panel">
              <div class="rc-lbl">规则证据评估（有支持证据 {expert_dimension_scored}/{expert_dimension_total or 8}；已评价但未观察 {expert_dimension_no_evidence}；不适用/未建立 {expert_dimension_na}）</div>
              <div class="rc-bars">{_rc_bars(m)}</div>
              <div class="diag-footnote">各维度表示当前保存证据的支持状态；“不适用”表示没有可执行规则。旧式聚合值仅保存在机器可读审计字段中，不是概率，也不参与身份排序。</div>
            </div>
          </div>
        </div>
      </div>
      <div class="rc-fw">
        <div class="rc-lbl">查询谱 vs 参考库谱 Query / Library MS/MS</div>
        <div class="rc-spec">{svg_html}</div>
      </div>
      <div class="rc-fw">
        <div class="rc-lbl">共同碎片离子对 Matched Fragment Pairs</div>
        <div style="font-size:11px;color:#64748b;margin-bottom:8px">下表只列分析评分实际使用的质控后样品峰：{match_count} 对与参考谱的 m/z 差不超过 0.020 Da。原始保存峰表的只读复核得到 {review_match_count} 对，仅作审计，不回写评分。</div>
        {frag_comp_html}
      </div>
      <div class="rc-fw">
        <div class="rc-lbl">诊断离子 Diagnostic Ions</div>
        {diagnostic_html}
      </div>
      <div class="rc-fw">
        <div class="rc-lbl">经验碎裂证据与中性丢失序列 Empirical Fragment Evidence / Neutral-Loss Sequences</div>
        {_rc_fragpath(m, qpeaks, r.get("evograph_fragmentation"))}
      </div>
      <div class="rc-fw">
        <div class="rc-lbl">结构推理 Structure Reasoning</div>
        {_rc_reasoning(m)}
      </div>
      <div class="rc-fw">
        <div class="rc-lbl">异构体竞争候选／其他谱库候选 Isomer Competition / Alternative Library Matches</div>
        {_rc_topk(m)}
      </div>
      <div class="rc-fw">
        <details class="rc-big"><summary>📋 完整峰表 Peak Tables (Query / Library) — 点击展开</summary>
          <div class="rc-pk-grid">
            <div><div class="rc-pk-h">样品 MS/MS 完整保存峰表（{len(qpeaks_full)} 峰；源记录 {source_peak_count} 峰）</div>{query_peaks_html}</div>
            <div><div class="rc-pk-h">Library 参考库谱峰（{len(lpeaks)} 峰）</div>{lib_peaks_html}</div>
          </div>
          <div class="diag-footnote">{raw_peak_storage_note}</div>
        </details>
      </div>
    </div>
  </td>
</tr>'''
        row_idx += 1

    h += '''</tbody></table></div></div>'''

    # Keep the complete candidate set in the offline file, with closed native
    # disclosures on screen and print-expanded content for paper/PDF output.
    h += render_isomer_candidate_competition_html(
        matched, language="zh"
    )

    unmatched_list = [r for r in results if not r.get("matches")]
    if unmatched_list:
        h += f'<div class="no-match"><strong>未找到满足当前筛选条件的谱库候选：</strong>{len(unmatched_list)} 张 MS2 谱图</div>'

    if governance_blocked_candidates is None:
        scan_rows = []
        for row in projected_results:
            matches = row.get("matches") or []
            scan_rows.append({
                **row,
                "match": matches[0] if matches else None,
            })
        governance_blocked_candidates = project_governance_blocked_candidates(
            scan_rows
        )
    h += render_governance_blocked_html(
        governance_blocked_candidates, language="zh"
    )

    h += '''</div>

<script>
function toggleFragRow(rowId, mainRow) {
  var frag = document.getElementById("frag_" + rowId);
  var icon = document.getElementById("icon_" + rowId);
  var isOpen = frag.classList.contains("visible");
  if (isOpen) {
    frag.classList.remove("visible");
    mainRow.classList.remove("expanded");
    icon.textContent = "➕";
  } else {
    frag.classList.add("visible");
    mainRow.classList.add("expanded");
    icon.textContent = "➖";
  }
}

function applyFilters() {
  var sortBy = document.getElementById("sortSelect").value;
  var nameFilter = document.getElementById("nameFilter").value.toLowerCase();
  var confFilter = document.getElementById("confFilter").value;
  var tbody = document.getElementById("tableBody");
  var rows = Array.from(tbody.querySelectorAll("tr.main-row"));

  // Build items with data
  var items = rows.map(function(r) {
    return {
      row: r,
      cosine: parseFloat(r.dataset.cosine),
      final: parseFloat(r.dataset.final),
      ppm: parseFloat(r.dataset.ppm),
      frags: parseInt(r.dataset.frags),
      name: r.dataset.name,
      // 通用搜索索引: top-1 名称+竞争候选名称+中文别名; data-mz 为前体 m/z
      search: (r.dataset.search || r.dataset.name || "").toLowerCase(),
      mz: parseFloat(r.dataset.mz || "0"),
      cls: r.dataset.class,
      idx: parseInt(r.dataset.idx)
    };
  });

  // 置信度取 honest_calibration 化合物级 (data-conf, 与统计卡同源); 不再用 final 阈值(旧版与之矛盾)
  items.forEach(function(item) {
    item.conf = item.row.dataset.conf || "Unknown";
  });

  // Filter
  var filtered = items;
  if (nameFilter) {
    if (/^[0-9]+([.][0-9]+)?$/.test(nameFilter)) {
      // 纯数字输入 = 按前体 m/z 搜索 (±0.5 Da)
      var qmz = parseFloat(nameFilter);
      filtered = filtered.filter(function(item) { return item.mz > 0 && Math.abs(item.mz - qmz) <= 0.5; });
    } else {
      // 名称/中文别名搜索 (含该行全部竞争候选名称)
      filtered = filtered.filter(function(item) { return item.search.indexOf(nameFilter) >= 0; });
    }
  }
  if (confFilter !== "all") filtered = filtered.filter(function(item) { return item.conf === confFilter; });

  // Sort（对全量排序；匹配行按序显示，不匹配隐藏——绝不删除 DOM 行：
  // 旧版空结果会移除全部行且无法恢复，导致之后任何搜索都得到空表）
  var cmp;
  if (sortBy === "final_desc") cmp = function(a,b) { return b.final - a.final; };
  else if (sortBy === "cosine_desc") cmp = function(a,b) { return b.cosine - a.cosine; };
  else if (sortBy === "ppm_asc") cmp = function(a,b) { return a.ppm - b.ppm; };
  else if (sortBy === "frags_desc") cmp = function(a,b) { return b.frags - a.frags; };
  else if (sortBy === "name_asc") cmp = function(a,b) { return a.name.localeCompare(b.name); };
  else cmp = function(a,b) { return 0; };
  items.sort(cmp);

  // Keep set for visibility
  var toKeep = {};
  filtered.forEach(function(item) { toKeep[item.idx] = true; });

  // Re-append all rows in sorted order; hide non-matching instead of removing
  items.forEach(function(item) {
    if (!item.frag) item.frag = item.row.nextElementSibling;
    var show = toKeep[item.idx] ? "" : "none";
    item.row.style.display = show;
    if (item.frag) item.frag.style.display = show;
    tbody.appendChild(item.row);
    if (item.frag) tbody.appendChild(item.frag);
  });

  // Update count
  document.getElementById("resultCount").textContent = "当前显示 " + filtered.length + " / " + items.length + " 个候选";
}

// Initial
applyFilters();
</script>
</body></html>'''
    return h


# ═══════════════════════════════════════════════════════════════════
# Task Persistence & Serialization
# ═══════════════════════════════════════════════════════════════════

def _now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def _load_tasks_from_disk(session_state):
    """Attach every session to the process's single current task mapping."""
    _task_registry(TASKS_FILE).attach(session_state)


def _save_tasks_to_disk(session_state):
    """Reject stale mappings and propagate disk failures to the caller."""
    _task_registry(TASKS_FILE).save(session_state.tasks)


def _save_task_results(task_id: str, results: list):
    """Publish complete results only while the shared task generation is valid."""
    with _tasks_lock:
        _task_registry(TASKS_FILE).assert_current()
        _atomic_json_write(TASK_DIR / f"{task_id}_results.json", results)


def _delete_task_from_session(task_id: str):
    """Reject active work before changing shared tasks or removing its files."""
    try:
        session_state = capture_session_state()
        with _tasks_lock:
            _library_coordinator.assert_task_idle(task_id)
            registry = _task_registry(TASKS_FILE)
            registry.assert_current(session_state.tasks)
            task = session_state.tasks.get(task_id)
            if task is None:
                return False
            if task.get("status") in {"queued", "running", "finalizing", "cancelling"}:
                raise RuntimeError("This task is marked active. Use the application that owns it before deleting it.")
            removed = session_state.tasks.pop(task_id)
            try:
                registry.save(session_state.tasks)
            except BaseException:
                session_state.tasks[task_id] = removed
                raise
            # Files are removed only after a successful atomic index update.
            # A deletion failure is reported; no background lease was invalidated.
            path = Path(task.get("file_path", ""))
            if path.is_file():
                path.unlink()
            result_path = TASK_DIR / f"{task_id}_results.json"
            result_path.unlink(missing_ok=True)
        return True
    except Exception as exc:
        st.error(str(exc))
        return False


def _load_task_results(task_id: str) -> list:
    """Load results from disk on demand."""
    path = TASK_DIR / f"{task_id}_results.json"
    if path.exists():
        try:
            with open(path, "r", encoding="utf-8") as f:
                return json.load(f)
        except Exception:
            pass
    return []


def _msms_fragments_str(query_peaks, lib_peaks, tol=0.02, top_n=20):
    """Build semicolon-separated MS/MS fragment m/z string from matched query peaks."""
    if not query_peaks or not lib_peaks:
        return ""
    lib_mz_list = [lp[0] for lp in lib_peaks]
    matched = []
    for qp in query_peaks:
        for lmz in lib_mz_list:
            if abs(qp[0] - lmz) <= tol:
                matched.append(f"{qp[0]:.4f}")
                break
    return "; ".join(matched[:top_n]) if matched else ""


def _load_library_into_session(db_path, db_name):
    """Switch the shared library only when no analysis holds a reservation."""
    def prepare():
        st.session_state.records = None
        st.session_state.idx = None
        st.session_state.lib_loaded = False
        st.session_state.lib_name = None
        st.session_state.lib_spectra_count = 0
        st.session_state.lib_load_time = 0.0
    started = time.time()
    with _tasks_lock:
        snapshot = _library_coordinator.replace(lambda: load_library(db_path), prepare)
        elapsed = time.time() - started
        st.session_state.records = snapshot.records
        st.session_state.idx = snapshot.index
        st.session_state.lib_name = db_name
        st.session_state.lib_spectra_count = len(snapshot.records)
        st.session_state.lib_load_time = round(elapsed, 1)
        st.session_state.lib_loaded = True
        return len(snapshot.records), elapsed


def _parse_spectrum_file(
    file_path: str, isolation_window_width_da: float | None = None,
) -> list:
    """Parse MGF / mzML / raw_blank_analysis HTML 报告 into spectrum list."""
    from utils.mgf_io import parse_mgf_file as _parse_mgf
    from pipeline.mzml_input import parse_mzml_file as _parse_mzml
    from utils.report_html_io import (
        parse_raw_report_html as _parse_html,
        is_raw_report_html as _is_report_html,
    )
    ext = Path(file_path).suffix.lower()
    if ext in (".mgf",):
        return _parse_mgf(file_path)
    elif ext in (".mzml", ".mzxml", ".xml"):
        return _parse_mzml(
            file_path, isolation_window_width_da=isolation_window_width_da,
        )
    elif ext in (".html", ".htm"):
        # 外部前处理工具 (raw_blank_analysis_v2) 的 HTML 峰值报告, 内嵌 MS2 谱图;
        # 先探测签名, 防止把 mssystem 自身导出的结果 HTML 当样本解析。
        if not _is_report_html(file_path):
            raise ValueError(
                "该 .html 不是可解析的分析报告 (缺 MS2 碎片/前体签名)。"
                "仅支持 raw_blank_analysis 风格的峰值报告。"
            )
        return _parse_html(file_path)
    else:
        raise ValueError(f"Unsupported file format: {ext}")


def _serialize_results(results: list) -> list:
    """Convert raw search results to JSON-serializable format (same as FastAPI)."""
    def _serialized_peak_pairs(peaks):
        saved = []
        for peak in peaks or []:
            try:
                mz_value = float(peak[0])
                intensity_value = float(peak[1])
            except (TypeError, ValueError, IndexError, KeyError):
                continue
            if not (
                math.isfinite(mz_value) and math.isfinite(intensity_value)
            ):
                continue
            saved.append([mz_value, intensity_value])
        return saved

    out = []
    for r in results:
        m = r.get("match")
        raw_input_peaks = r.get("raw_peaks") or []
        raw_saved_peaks = _serialized_peak_pairs(raw_input_peaks)
        serialized = {
            "spectrum_id": str(r.get("spectrum_id", "")),
            "precursor_mz": r["precursor_mz"],
            "peaks_count": r.get("peaks_count", 0),
            "retention_time": r.get("retention_time"),
            "retention_time_unit": r.get("retention_time_unit"),
            "ion_mode": r.get("ion_mode"),
            "precursor_charge": r.get("precursor_charge"),
            **query_acquisition_projection(r),
            "precursor_intensity": r.get("precursor_intensity"),
            "precursor_ms1_intensity": r.get("precursor_ms1_intensity"),
            "collision_energy": r.get("collision_energy"),
            "collision_energy_unit": r.get("collision_energy_unit", ""),
            "instrument_type": r.get("instrument_type", ""),
            "instrument": r.get("instrument", ""),
            "rt_method_id": r.get("rt_method_id", ""),
            "precursor_isolation_evidence": r.get("precursor_isolation_evidence"),
            "input_format": r.get("input_format"),
            "limited_fragment_list": bool(r.get("limited_fragment_list", False)),
            "input_evidence_note": str(r.get("input_evidence_note") or ""),
            "raw_peaks": raw_saved_peaks,
            "raw_peaks_input_count": len(raw_input_peaks),
            "raw_peaks_saved_count": len(raw_saved_peaks),
            "raw_peaks_complete": len(raw_saved_peaks) == len(raw_input_peaks),
            "raw_peaks_storage_policy": "all_valid_numeric_pairs",
        }
        if m:
            qp = r.get("query_peaks", [])
            lp = m.get("lib_peaks", [])
            msms_str = _msms_fragments_str(qp, lp)
            matched_count = len(msms_str.split("; ")) if msms_str else 0
            top_qp = sorted(qp, key=lambda p: -p[1])[:5] if qp else []
            top_qp_str = "; ".join(f"{p[0]:.4f}" for p in top_qp) if top_qp else ""

            ann = _honest_annotate(m, m.get("top_candidates"))
            class_evidence = classification_evidence_status(m)
            serialized["match"] = {
                **{key: m[key] for key in ("query_adduct_charge_evidence",) if key in m},
                "cosine": round(m.get("cosine", 0), 4),
                "mass_score": round(m.get("mass_score", 0), 4),
                # 0 ppm is valid evidence, not a missing value.
                "mass_ppm": round(m["mass_ppm"], 2) if m.get("mass_ppm") is not None else None,
                "fragment_score": round(m.get("fragment_score", 0), 4),
                "intensity_score": round(m.get("intensity_score", 0), 4),
                "isomer_fragment_score": round(m.get("isomer_fragment_score", 0) or 0, 4),
                "isomer_fragment_bonus": round(m.get("isomer_fragment_bonus", 0) or 0, 4),
                "isomer_fragment_detail": m.get("isomer_fragment_detail"),
                "multimetric_shadow": m.get("multimetric_shadow"),
                "final_score": round(m.get("final_score", 0), 4),
                "rank_score": round(m.get("rank_score", m.get("final_score", 0)) or 0, 4),
                "raw_cosine": m.get("raw_cosine"),
                "raw_final_score": m.get("raw_final_score"),
                "chimeric_deconvolution": m.get("chimeric_deconvolution"),
                "compound_name": str(m.get("compound_name", ""))[:200],
                "smiles": str(m.get("smiles", "")),
                "inchikey": str(m.get("inchikey", "")),
                "formula": str(m.get("formula", "")),
                "adduct": str(m.get("adduct", "")),
                "classification": str(class_evidence["label"]),
                "display_class": str(m.get("display_class") or ""),
                "raw_classification": str(m.get("classification", "")),
                "classification_source": str(class_evidence["source"]),
                "display_classification_source": str(
                    m.get("display_classification_source") or ""
                ),
                "classification_reliable": bool(class_evidence["reliable"]),
                "classification_status": str(class_evidence["status"]),
                "classification_reliability_reason": str(
                    class_evidence["reason"]
                ),
                "dataset_source": str(m.get("dataset_source", "")),
                "confidence": ann.get("compound_confidence", "Unknown"),
                "confidence_class": ann.get("class_confidence", "Unknown"),
                "exact_structure_confidence": ann.get(
                    "exact_structure_confidence", "Unknown"
                ),
                "isomer_group_confidence": ann.get(
                    "isomer_group_confidence", "Unknown"
                ),
                "identification_scope": ann.get("identification_scope"),
                "isomer_group_formula": ann.get("isomer_group_formula", ""),
                "isomer_group_size": ann.get("isomer_group_size", 0),
                "isomer_group_evidence": ann.get(
                    "isomer_group_evidence", {}
                ),
                "isomer_ambiguity_type": ann.get(
                    "isomer_ambiguity_type", "none"
                ),
                "isomer_relationship_counts": ann.get(
                    "isomer_relationship_counts", {}
                ),
                "isomer_variant_subtype_counts": ann.get(
                    "isomer_variant_subtype_counts", {}
                ),
                "isomer_resolution_status": ann.get(
                    "isomer_resolution_status", "not_applicable"
                ),
                "isomer_literature_evidence": ann.get(
                    "isomer_literature_evidence", {}
                ),
                "honesty_tier": ann.get("honesty_tier"),
                "honest_headline": ann.get("honest_headline"),
                "evidence_tier": ann.get("evidence_tier"),
                "annotation_level": ann.get("annotation_level"),
                "candidate_comparison_source": ann.get(
                    "candidate_comparison_source"
                ),
                "candidate_comparison_cosine": ann.get(
                    "candidate_comparison_cosine"
                ),
                "candidate_margin": ann.get("candidate_margin"),
                "source_tier": m.get("source_tier", "unknown"),
                "has_structure": m.get("has_structure"),
                "reverse_match": m.get("reverse_match"),
                "top_candidates": m.get("top_candidates") or [],
                "ms_ms_fragments": msms_str,
                "fragment_count": matched_count,
                "top_query_peaks": top_qp_str,
                "expert_score": round(m.get("expert_score", 0), 4),
                "expert_details": m.get("expert_details"),
                "dl_classification": str(m.get("dl_classification", "")),
                "dl_confidence": round(m.get("dl_confidence", 0), 4),
                "dl_top3": m.get("dl_top3", []),
                "structure_score": round(m.get("structure_score", 0), 4),
                "predicted_smiles": str(m.get("predicted_smiles", "")),
            }
            serialized["query_peaks"] = _serialized_peak_pairs(qp)
            serialized["lib_peaks"] = _serialized_peak_pairs(lp)
            serialized["deconvolved_query_peaks"] = [
                [float(p[0]), float(p[1])]
                for p in (m.get("deconvolved_query_peaks") or [])
            ]
        else:
            serialized["match"] = None
            serialized["query_peaks"] = _serialized_peak_pairs(
                r.get("query_peaks", [])
            )
            serialized["lib_peaks"] = []
            # 未命中谱图用 MSClassify 从谱图判类 (补无化合物命中的判类空白)
            try:
                from classifier.msclass_wrapper import get_msclass_wrapper
                _msc_w = get_msclass_wrapper()
                if _msc_w is not None:
                    _msc_r = _msc_w.predict(r.get("query_peaks", []))
                    if _msc_r and _msc_r["confidence"] >= 0.15:
                        serialized["display_class"] = _msc_r["class"]
                        serialized["msclass_source"] = "unmatched"
                        serialized["dl_confidence"] = round(_msc_r["confidence"], 3)
                        serialized["dl_top3"] = _msc_r["top3"]
            except Exception:
                pass
        out.append(serialized)
    return out


# ═══════════════════════════════════════════════════════════════════
# Background Analysis Thread
# ═══════════════════════════════════════════════════════════════════

def _start_analysis_in_thread(task_id: str, params=None):
    """Reserve a generation before constructing a thread, for every UI entry."""
    try:
        session_state = capture_session_state()
        with _tasks_lock:
            if session_state.analysis_running:
                raise RuntimeError("Another analysis is running in this session. Please wait.")
            if not session_state.lib_loaded:
                raise RuntimeError("Load the library before starting analysis.")
            _task_registry(TASKS_FILE).assert_current(session_state.tasks)
            task = session_state.tasks.get(task_id)
            if task is None:
                raise RuntimeError("The requested task no longer exists.")
            lease = _library_coordinator.reserve(task_id, session_state.records, session_state.idx)
            try:
                if params is not None:
                    task["params"] = dict(params)
                task["status"] = "queued"
                task["progress"] = 0
                task["error"] = None
                task["message"] = "Queued for analysis..."
                session_state.analysis_running = True
                session_state.active_task_id = task_id
                _save_tasks_to_disk(session_state)
                thread = threading.Thread(target=_run_analysis_in_thread,
                                          args=(task_id, session_state, lease), daemon=True)
                thread.start()
            except BaseException as exc:
                _library_coordinator.release(lease)
                session_state.analysis_running = False
                session_state.active_task_id = None
                task["status"] = "failed"
                task["error"] = str(exc)
                task["message"] = "Analysis thread could not be started. Retry the task."
                _save_tasks_to_disk(session_state)
                raise
        return True
    except Exception as exc:
        st.error(str(exc))
        return False


def _run_analysis_in_thread(task_id: str, session_state, lease):
    """Release a reservation even when setup, parsing or finalization fails."""
    try:
        _library_coordinator.validate(lease)
        _run_analysis_with_library(task_id, session_state, lease)
    except Exception as exc:
        _logger.exception("[task %s] FAILED before finalization: %s", task_id, exc)
        with _tasks_lock:
            task = session_state.tasks.get(task_id)
            if task is not None:
                task["status"] = "failed"
                task["error"] = str(exc)
                task["message"] = f"Error: {exc}"
    finally:
        with _tasks_lock:
            _library_coordinator.release(lease)
            if session_state.active_task_id == task_id:
                session_state.analysis_running = False
                session_state.active_task_id = None
            _save_tasks_to_disk(session_state)


def _run_analysis_with_library(task_id: str, session_state, lease):
    """Execute analysis in background thread (mirrors FastAPI's _run_analysis)."""
    with _tasks_lock:
        task = session_state.tasks.get(task_id)
        if not task:
            return
        task["status"] = "running"
        task["start_time"] = _now_iso()
        task["progress"] = 0
        task["message"] = "Starting analysis..."
    _save_tasks_to_disk(session_state)

    t_start = time.time()
    try:
        file_path = task["file_path"]
        params = task.get("params", {})

        # ── Phase 1: Parse spectra ──
        _logger.info("[task %s] Phase 1/3: Parsing file %s ...", task_id, Path(file_path).name)
        with _tasks_lock:
            session_state.tasks[task_id]["progress"] = 2
            session_state.tasks[task_id]["message"] = "Parsing spectrum file..."
        _save_tasks_to_disk(session_state)

        spectra = _parse_spectrum_file(
            file_path,
            isolation_window_width_da=params.get("isolation_window_width_da"),
        )
        total = len(spectra)

        if total == 0:
            raise ValueError("No valid spectra found in file.")

        with _tasks_lock:
            session_state.tasks[task_id]["total_spectra"] = total
            session_state.tasks[task_id]["progress"] = 5
            session_state.tasks[task_id]["message"] = f"Parsed {total} spectra. Starting matching..."
        _save_tasks_to_disk(session_state)

        # ── Phase 2: Spectral matching ──
        _logger.info("[task %s] Phase 2/3: Matching %d spectra...", task_id, total)
        deconvolution_run_pool = (
            build_deconvolution_run_pool(spectra)
            if params.get(
                "enable_chimeric_deconvolution", USE_CHIMERIC_DECONVOLUTION
            ) else ([], {})
        )
        all_results = []
        n_matched = 0
        last_progress_val = 5

        for i, spec in enumerate(spectra):
            _library_coordinator.validate(lease)
            match = search(
                spec["peaks"],
                spec["precursor_mz"],
                records=lease.snapshot.records,
                mz_index=lease.snapshot.index,
                cos_thresh=params.get("cosine_threshold", COSINE_THRESHOLD),
                final_thresh=params.get("final_score_threshold", FINAL_SCORE_THRESHOLD),
                mz_da=params.get("mz_tolerance_da", MZ_TOLERANCE_DA),
                mz_ppm=params.get("mz_tolerance_ppm", MZ_TOLERANCE_PPM),
                top_n=params.get("top_n_peaks", TOP_N_PEAKS),
                strict_evidence_gate=params.get(
                    "strict_evidence_gate", USE_STRICT_EVIDENCE_GATE
                ) or bool(spec.get("limited_fragment_list")),
                strict_max_mass_ppm=(
                    LIMITED_FRAGMENT_STRICT_MAX_MASS_PPM
                    if spec.get("limited_fragment_list") else None
                ),
                query_adduct=spec.get("adduct"),
                query_adduct_reliable=spec.get("adduct_reliable"),
                query_polarity=spec.get("polarity"),
                query_rt=spec.get("retention_time"),
                query_rt_unit=(
                    "seconds" if Path(file_path).suffix.lower() == ".mgf" else "minutes"
                ),
                query_charge=spec.get("precursor_charge"),
                query_collision_energy=collision_energy_evidence(spec),
                query_instrument_type=spec.get("instrument_type"),
                query_instrument=spec.get("instrument"),
                query_rt_method_id=params.get("rt_method_id"),
                query_isotope_evidence=spec.get("precursor_isotope_evidence"),
                isotope_evidence_status=spec.get("isotope_evidence_status"),
                query_isolation_evidence=spec.get("precursor_isolation_evidence"),
                enable_chimeric_deconvolution=params.get(
                    "enable_chimeric_deconvolution", USE_CHIMERIC_DECONVOLUTION,
                ),
                deconvolution_run_pool=deconvolution_run_pool,
            )
            all_results.append({
                "spectrum_id": spec.get("spectrum_id", i),
                "precursor_mz": spec["precursor_mz"],
                "peaks_count": len(spec["peaks"]),
                "retention_time": spec.get("retention_time"),
                "retention_time_unit": spec.get("retention_time_unit", "minutes"),
                "ion_mode": spec.get("polarity"),
                "precursor_charge": spec.get("precursor_charge"),
                **query_acquisition_projection(spec),
                "precursor_intensity": spec.get("precursor_intensity"),
                "precursor_ms1_intensity": spec.get("precursor_ms1_intensity"),
                "collision_energy": spec.get("collision_energy"),
                "collision_energy_unit": spec.get("collision_energy_unit", ""),
                "instrument_type": spec.get("instrument_type", ""),
                "instrument": spec.get("instrument", ""),
                "rt_method_id": params.get("rt_method_id"),
                "isotope_evidence_status": spec.get(
                    "isotope_evidence_status", "not_available"
                ),
                "precursor_isotope_evidence": spec.get("precursor_isotope_evidence"),
                "precursor_isolation_evidence": spec.get(
                    "precursor_isolation_evidence"
                ),
                "input_format": spec.get("input_format"),
                "limited_fragment_list": bool(spec.get("limited_fragment_list", False)),
                "input_evidence_note": spec.get("input_evidence_note", ""),
                "match": match,
                "query_peaks": spec["peaks"],
            })
            if match is not None:
                n_matched += 1

            # Update progress periodically
            new_progress = 5 + int((i + 1) / total * 90)
            now = time.time()
            if new_progress != last_progress_val:
                elapsed = now - t_start
                rate = (i + 1) / elapsed if elapsed > 0 else 0
                eta = (total - i - 1) / rate if rate > 0 else 0
                with _tasks_lock:
                    session_state.tasks[task_id]["progress"] = new_progress
                    session_state.tasks[task_id]["message"] = (
                        f"Matching... {i+1}/{total} spectra "
                        f"({n_matched} matched, {rate:.1f} spec/s, ETA {eta:.0f}s)"
                    )
                last_progress_val = new_progress
                _save_tasks_to_disk(session_state)

        # ── Phase 3: Serialize & Finalize ──
        _logger.info("[task %s] Phase 3/3: Serializing results...", task_id)
        with _tasks_lock:
            session_state.tasks[task_id]["progress"] = 96
            session_state.tasks[task_id]["message"] = f"Serializing {n_matched} results..."
        _save_tasks_to_disk(session_state)

        _library_coordinator.validate(lease)
        serialized = _serialize_results(all_results)
        annotate_precursor_features(serialized)
        if params.get(
            "enable_ion_identity_networking", USE_ION_IDENTITY_NETWORKING
        ):
            ion_identity_summary = annotate_ion_identity_networks(
                serialized,
                file_path=file_path,
                max_rt_delta_min=ION_IDENTITY_RT_TOLERANCE_MIN,
                min_shape_correlation=ION_IDENTITY_MIN_SHAPE_CORRELATION,
                min_overlap_fraction=ION_IDENTITY_MIN_OVERLAP,
            )
            for row in serialized:
                match = row.get("match")
                if not isinstance(match, dict):
                    continue
                ann = _honest_annotate(match, match.get("top_candidates"))
                match["confidence"] = ann["compound_confidence"]
                match["confidence_class"] = ann["class_confidence"]
                for key in (
                    "exact_structure_confidence", "isomer_group_confidence",
                    "identification_scope", "isomer_group_formula",
                    "isomer_group_size", "isomer_group_evidence",
                    "isomer_variant_subtype_counts",
                    "isomer_literature_evidence",
                    "honesty_tier", "honest_headline", "evidence_tier",
                    "annotation_level", "isomer_ambiguous", "isomer_n",
                    "isomer_set", "class_consensus", "class_agreement",
                    "candidate_comparison_source",
                    "candidate_comparison_cosine", "candidate_margin",
                ):
                    match[key] = ann.get(key)
        else:
            ion_identity_summary = {"status": "disabled", "network_count": 0}
        _save_task_results(task_id, serialized)

        t_total = time.time() - t_start
        match_rate = n_matched / total * 100 if total else 0

        with _tasks_lock:
            session_state.tasks[task_id]["status"] = "completed"
            session_state.tasks[task_id]["progress"] = 100
            session_state.tasks[task_id]["message"] = (
                f"Completed: {n_matched}/{total} matched ({match_rate:.1f}%) in {t_total:.0f}s"
            )
            session_state.tasks[task_id]["end_time"] = _now_iso()
            session_state.tasks[task_id]["matched_spectra"] = n_matched
            session_state.tasks[task_id]["match_rate"] = round(match_rate, 1)
            session_state.tasks[task_id]["total_spectra"] = total
            session_state.tasks[task_id]["ion_identity_summary"] = ion_identity_summary
            session_state.tasks[task_id]["ion_identity_networks"] = int(
                ion_identity_summary.get("network_count") or 0
            )
        _save_tasks_to_disk(session_state)
        _logger.info("[task %s] COMPLETED: %d/%d matched (%.1f%%) in %.1fs",
                     task_id, n_matched, total, match_rate, t_total)

    except Exception as e:
        _logger.exception("[task %s] FAILED: %s", task_id, e)
        with _tasks_lock:
            session_state.tasks[task_id]["status"] = "failed"
            session_state.tasks[task_id]["error"] = str(e)
            session_state.tasks[task_id]["message"] = f"Error: {e}"
        _save_tasks_to_disk(session_state)


# ═══════════════════════════════════════════════════════════════════
# Export Helpers (Excel, Batch ZIP)
# ═══════════════════════════════════════════════════════════════════

def _export_excel(task: dict, results: list) -> bytes:
    """Generate multi-sheet Excel workbook (Matched, Unmatched, Summary)."""
    matched = [r for r in results if r.get("match")]
    unmatched = [r for r in results if not r.get("match")]

    match_data = []
    for r in matched:
        m = r["match"]
        match_data.append({
            "Spectrum": r["spectrum_id"], "Precursor_mz": round(r["precursor_mz"], 4),
            "Peaks_Count": r.get("peaks_count", 0), "Retention_Time": r.get("retention_time") or "",
            "Compound_Name": m["compound_name"], "SMILES": m.get("smiles", ""),
            "InChIKey": m.get("inchikey", ""), "Formula": m.get("formula", ""),
            "Adduct": m.get("adduct", ""), "Classification": m.get("classification", ""),
            "Confidence": m.get("confidence", ""), "Cosine_Score": m["cosine"],
            "Final_Score": m["final_score"], "Mass_Score": m.get("mass_score", 0),
            "Mass_PPM": m.get("mass_ppm"), "Fragment_Score": m.get("fragment_score", 0),
            "Intensity_Score": m.get("intensity_score", 0),
            "MS_MS_Fragments": m.get("ms_ms_fragments", ""),
            "Fragment_Count": m.get("fragment_count", 0),
        })

    df_match = pd.DataFrame(match_data) if match_data else pd.DataFrame()
    df_unmatch = pd.DataFrame([{
        "Spectrum": r["spectrum_id"], "Precursor_mz": round(r["precursor_mz"], 4),
        "Peaks_Count": r.get("peaks_count", 0),
    }] for r in unmatched) if unmatched else pd.DataFrame()
    governance_blocked_rows = governance_blocked_export_rows(
        project_governance_blocked_candidates(results)
    )

    buf = io.BytesIO()
    with pd.ExcelWriter(buf, engine="openpyxl") as writer:
        if not df_match.empty:
            df_match.to_excel(writer, sheet_name="Matched", index=False)
        if not df_unmatch.empty:
            df_unmatch.to_excel(writer, sheet_name="Unmatched", index=False)
        if governance_blocked_rows:
            pd.DataFrame(governance_blocked_rows).to_excel(
                writer, sheet_name="治理未放行候选", index=False
            )
        summary = pd.DataFrame([{
            "Filename": task["filename"],
            "Total_Spectra": task.get("total_spectra", 0),
            "Matched": task.get("matched_spectra", 0),
            "Match_Rate_%": round(task.get("match_rate", 0), 1),
            "Upload_Time": task.get("upload_time", ""),
            "Analysis_Time": task.get("end_time", ""),
        }])
        summary.to_excel(writer, sheet_name="Summary", index=False)
    buf.seek(0)
    return buf.getvalue()


def _export_batch_zip(task_ids: list, session_state) -> bytes:
    """Package multiple task results into a ZIP file."""
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
        for tid in task_ids:
            task = session_state.tasks.get(tid)
            if not task:
                continue
            results = _load_task_results(tid)
            matched = [r for r in results if r.get("match")]
            if not matched:
                continue
            csv_buf = io.StringIO()
            w = csv.writer(csv_buf)
            w.writerow(["Spectrum", "Precursor_mz", "Compound_Name", "SMILES",
                         "InChIKey", "Formula", "Final_Score", "Cosine_Score",
                         "Confidence", "Classification", "MS_MS_Fragments", "Fragment_Count"])
            for r in matched:
                m = r["match"]
                w.writerow([
                    r["spectrum_id"], f"{r['precursor_mz']:.4f}",
                    m["compound_name"], m.get("smiles", ""),
                    m.get("inchikey", ""), m.get("formula", ""),
                    f"{m['final_score']:.4f}", f"{m['cosine']:.4f}",
                    m.get("confidence", ""), m.get("classification", ""),
                    m.get("ms_ms_fragments", ""), m.get("fragment_count", 0),
                ])
            safe_name = Path(task["filename"]).stem
            zf.writestr(f"{safe_name}_results.csv", csv_buf.getvalue())
    buf.seek(0)
    return buf.getvalue()


# ═══════════════════════════════════════════════════════════════════
# Page Setup (only runs via `streamlit run web/app.py`)
# ═══════════════════════════════════════════════════════════════════

if __name__ == "__main__":

    st.set_page_config(page_title="MSsystem · LC-MS/MS Spectral Matching", page_icon="🧪", layout="wide")
    _inject_custom_css()

    # ═══════════════════════════════════════════════════════════════════
    # Session State
    # ═══════════════════════════════════════════════════════════════════

    _SESSION_DEFAULTS = {
        "lib_loaded": False, "records": None, "idx": None,
        "lib_name": "", "lib_spectra_count": 0, "lib_load_time": 0.0,
        "tasks": {}, "tasks_loaded_from_disk": False,
        "analysis_running": False, "active_task_id": None,
        "selected_task_id": None, "selected_results": None,
    }
    for k, v in _SESSION_DEFAULTS.items():
        if k not in st.session_state:
            st.session_state[k] = v

    # Load tasks from disk on first run
    _load_tasks_from_disk(st.session_state)

    # ═══════════════════════════════════════════════════════════════════
    # Sidebar
    # ═══════════════════════════════════════════════════════════════════

    with st.sidebar:
        # ── Brand ──
        st.markdown("""
        <div style="display:flex;align-items:center;gap:10px;margin-bottom:1rem;">
            <div style="background:linear-gradient(135deg,#1A5253,#2C7A7B);color:white;width:38px;height:38px;
                        border-radius:8px;display:flex;align-items:center;justify-content:center;
                        font-weight:700;font-size:17px;">L</div>
            <div>
                <div style="font-weight:700;font-size:1.1rem;color:#1A202C;line-height:1.2;">MSsystem</div>
                <div style="font-size:0.7rem;color:#718096;">Spectral Matching</div>
            </div>
        </div>
        """, unsafe_allow_html=True)

        # ── Library Status ──
        if st.session_state.lib_loaded:
            n_spec = len(st.session_state.records)
            st.markdown(f"""<div style="display:flex;align-items:center;gap:8px;
                        background:#F0FFF4;border:1px solid #C6F6D5;border-radius:6px;
                        padding:8px 12px;margin-bottom:1rem;font-size:0.82rem;">
            <span style="color:#48BB78;font-size:1.1rem;">&#9679;</span>
            <span style="color:#276749;">{n_spec:,} spectra loaded</span>
            </div>""", unsafe_allow_html=True)
        else:
            st.markdown("""<div style="display:flex;align-items:center;gap:8px;
                        background:#FFF5F5;border:1px solid #FED7D7;border-radius:6px;
                        padding:8px 12px;margin-bottom:1rem;font-size:0.82rem;">
            <span style="color:#FC8181;font-size:1.1rem;">&#9679;</span>
            <span style="color:#9B2C2C;">Library not loaded</span>
            </div>""", unsafe_allow_html=True)

        # ── Database Section ──
        st.markdown("""<div class="sidebar-section">""", unsafe_allow_html=True)
        st.markdown('<div style="font-size:0.7rem;font-weight:600;text-transform:uppercase;letter-spacing:0.08em;color:#718096;margin-bottom:0.4rem;">Database</div>', unsafe_allow_html=True)

        db_mode = st.selectbox(
            "Library", ["Small (<800 Da)", "Full"],
            label_visibility="collapsed"
        )
        db_path = DB_PATH_SMALL if "Small" in db_mode else DB_PATH

        col_ld, col_re = st.columns([3, 1])
        with col_ld:
            if st.button("Load Library", use_container_width=True, type="primary"):
                with st.spinner("Loading library..."):
                    try:
                        n_spec, elapsed = _load_library_into_session(db_path, db_mode)
                    except Exception as exc:
                        st.error(str(exc))
                    else:
                        st.success(f"Loaded {n_spec:,} spectra in {elapsed:.1f}s")
        with col_re:
            if st.session_state.lib_loaded:
                if st.button("↻", use_container_width=True, help="Reload library"):
                    with st.spinner("Reloading..."):
                        try:
                            n_spec, elapsed = _load_library_into_session(db_path, db_mode)
                        except Exception as exc:
                            st.error(str(exc))
                        else:
                            st.success(f"Loaded {n_spec:,} spectra in {elapsed:.1f}s")

        st.markdown("</div>", unsafe_allow_html=True)

        # ── Settings Section ──
        st.markdown("""<div class="sidebar-section">""", unsafe_allow_html=True)
        st.markdown('<div style="font-size:0.7rem;font-weight:600;text-transform:uppercase;letter-spacing:0.08em;color:#718096;margin-bottom:0.4rem;">Thresholds</div>', unsafe_allow_html=True)

        cos_th = st.slider("Cosine ≥", 0.0, 1.0, COSINE_THRESHOLD, 0.05, help="Minimum cosine similarity for a match")
        fin_th = st.slider("Final Score ≥", 0.0, 1.0, FINAL_SCORE_THRESHOLD, 0.05, help="Minimum composite confidence score")
        mz_da = st.number_input("m/z tolerance (Da)", 0.01, 5.0, MZ_TOLERANCE_DA, 0.01)
        mz_ppm = st.number_input("m/z tolerance (ppm)", 1.0, 100.0, MZ_TOLERANCE_PPM, 1.0)

        with st.expander("Advanced"):
            top_n = st.slider("Top N Peaks", 5, 200, TOP_N_PEAKS, 5, help="Number of most intense peaks to use for matching")
            strict_gate = st.checkbox(
                "Strict evidence gate (final reports)",
                value=USE_STRICT_EVIDENCE_GATE,
                help="Adds Final>=0.50, Entropy>=0.45, Reverse>=0.55 and |ppm|<=100 requirements.",
            )
            enable_deconvolution = st.checkbox(
                "DecoID-like chimeric MS/MS deconvolution",
                value=USE_CHIMERIC_DECONVOLUTION,
                help="Only stable non-negative L1 fits from mixed/severe isolation windows can affect scores.",
            )
            enable_iimn = st.checkbox(
                "Ion Identity Molecular Networking",
                value=USE_ION_IDENTITY_NETWORKING,
                help="Merge adducts and in-source fragments using continuous-MS1 co-elution and neutral mass.",
            )
            isolation_width = st.number_input(
                "Isolation window full width (Da; 0 = read mzML only)",
                min_value=0.0, max_value=20.0, value=0.0, step=0.1,
            )

        st.markdown("</div>", unsafe_allow_html=True)

        # ── About ──
        st.markdown("---")
        n_tasks = len(st.session_state.tasks)
        st.caption(f"MSsystem v3.0 · {n_tasks} tasks · LC-MS/MS Spectral Matching")

    # ═══════════════════════════════════════════════════════════════════
    # Main Content — 4-Tab Layout
    # ═══════════════════════════════════════════════════════════════════

    tab_upload, tab_tasks, tab_results, tab_export = st.tabs([
        "📤 Upload & Analyze", "📋 Task Manager", "🔬 Results Viewer", "📥 Export Center"
    ])

    # ───────────────────────────────────────────────────────────────────
    # TAB 1: Upload & Analyze
    # ───────────────────────────────────────────────────────────────────
    with tab_upload:
        # ── Hero Section ──
        st.markdown("""
        <div class="hero-section fade-in">
            <div class="hero-title">LC-MS/MS Spectral Matching</div>
            <div class="hero-subtitle">
                Upload experimental mass spectra (.mgf, .mzML, .mzXML) or a raw_blank_analysis
                HTML report (.html) and identify compounds against a reference library. Powered by
                cosine similarity, fragment ion matching, and composite confidence scoring.
            </div>
        </div>
        """, unsafe_allow_html=True)

        # ── Upload Zone ──
        files = st.file_uploader(
            "Drop sample files here",
            type=["mgf", "mzML", "mzxml", "html"],
            accept_multiple_files=True,
            label_visibility="collapsed"
        )

        # ── Status Messages ──
        if not st.session_state.lib_loaded:
            if files:
                st.warning("⚠️ Load the reference library from the sidebar before running analysis.")
            else:
                st.info("💡 Load a reference library from the sidebar, then upload sample files to begin analysis.")
        elif not files:
            st.info("📂 Upload .mgf, .mzML, .mzXML, or .html (raw_blank_analysis 报告) files to begin spectral matching.")

        col_btn1, col_btn2, col_btn3 = st.columns([1, 2, 1])
        with col_btn2:
            run_disabled = not (files and st.session_state.lib_loaded)
            run_clicked = st.button(
                "🚀 Start Analysis",
                type="primary",
                use_container_width=True,
                disabled=run_disabled,
            )

        # ── Handle Analysis Start ──
        if run_clicked and files and st.session_state.lib_loaded:
            # Check if analysis is already running
            if st.session_state.analysis_running:
                st.warning("⚠️ Analysis is already running. Please wait for it to complete.")
            else:
                # Create task entries and start analysis for each file
                started_tasks = []
                for uf in files:
                    data = uf.read()
                    if len(data) == 0:
                        st.warning(f"Empty file: {uf.name}")
                        continue

                    safe_name = f"{uuid.uuid4().hex[:8]}_{uf.name}"
                    dest = UPLOAD_DIR / safe_name
                    with open(dest, "wb") as f:
                        f.write(data)

                    file_size_mb = len(data) / (1024 * 1024)
                    task_id = uuid.uuid4().hex[:12]
                    with _tasks_lock:
                        st.session_state.tasks[task_id] = {
                            "id": task_id,
                            "filename": uf.name,
                            "file_path": str(dest),
                            "file_size_mb": round(file_size_mb, 2),
                            "status": "pending",
                            "progress": 0,
                            "message": "File uploaded, queued for analysis",
                            "upload_time": _now_iso(),
                            "start_time": None, "end_time": None,
                            "total_spectra": 0, "matched_spectra": 0, "match_rate": 0.0,
                            "results": [], "error": None,
                            "params": {
                                "cosine_threshold": cos_th,
                                "final_score_threshold": fin_th,
                                "mz_tolerance_da": mz_da,
                                "mz_tolerance_ppm": mz_ppm,
                                "top_n_peaks": top_n,
                                "strict_evidence_gate": strict_gate,
                                "enable_chimeric_deconvolution": enable_deconvolution,
                                "enable_ion_identity_networking": enable_iimn,
                                "isolation_window_width_da": (
                                    isolation_width if isolation_width > 0 else None
                                ),
                            },
                        }
                    _save_tasks_to_disk(st.session_state)
                    started_tasks.append((task_id, uf.name))

                if started_tasks:
                    # Start first task in background thread
                    first_id = started_tasks[0][0]
                    if not _start_analysis_in_thread(first_id):
                        st.stop()

                    # Queue remaining tasks info
                    if len(started_tasks) > 1:
                        st.info(f"📋 {len(started_tasks)} files uploaded. First file analyzing now, "
                                f"remaining files are queued. Use Task Manager to start them.")

                    names = ", ".join(n for _, n in started_tasks)
                    st.success(f"✅ Analysis started for: {names}")

        # ── Progress Polling ──
        if st.session_state.analysis_running and st.session_state.active_task_id:
            tid = st.session_state.active_task_id
            task = st.session_state.tasks.get(tid, {})
            status = task.get("status", "pending")

            if status in {"queued", "running"}:
                progress = task.get("progress", 0)
                message = task.get("message", "Processing...")
                st.progress(progress / 100.0, text=message)
                time.sleep(1.5)
                st.rerun()
            elif status == "completed":
                st.success(f"✅ {task.get('filename', 'Analysis')} completed: "
                           f"{task.get('matched_spectra', 0)}/{task.get('total_spectra', 0)} matched "
                           f"({task.get('match_rate', 0):.1f}%)")
                st.session_state.analysis_running = False
                st.session_state.active_task_id = None
                st.balloons()
            elif status == "failed":
                st.error(f"❌ Analysis failed: {task.get('error', 'Unknown error')}")
                st.session_state.analysis_running = False
                st.session_state.active_task_id = None

        # ── Quick Stats ──
        tasks = st.session_state.tasks
        completed = [t for t in tasks.values() if t["status"] == "completed"]
        running = [t for t in tasks.values() if t["status"] == "running"]
        if completed or running:
            st.markdown("---")
            c1, c2, c3 = st.columns(3)
            with c1:
                st.markdown(f"""<div class="metric-card"><div class="metric-value">{len(completed)}</div>
                    <div class="metric-label">Completed</div></div>""", unsafe_allow_html=True)
            with c2:
                st.markdown(f"""<div class="metric-card"><div class="metric-value">{len(running)}</div>
                    <div class="metric-label">Running</div></div>""", unsafe_allow_html=True)
            with c3:
                total_matched = sum(t.get("matched_spectra", 0) for t in completed)
                st.markdown(f"""<div class="metric-card"><div class="metric-value">{total_matched}</div>
                    <div class="metric-label">Total Matches</div></div>""", unsafe_allow_html=True)

    # ───────────────────────────────────────────────────────────────────
    # TAB 2: Task Manager
    # ───────────────────────────────────────────────────────────────────
    with tab_tasks:
        st.markdown('<div class="section-header">Task Manager</div>', unsafe_allow_html=True)

        tasks = st.session_state.tasks
        task_list = list(tasks.values())

        # Filter
        col_f1, col_f2 = st.columns([3, 1])
        with col_f1:
            status_filter = st.selectbox("Filter by status",
                ["All", "pending", "queued", "running", "completed", "failed"],
                format_func=lambda x: {"All": "All", "pending": "⏳ Pending", "queued": "⏳ Queued",
                                        "running": "🔄 Running", "completed": "✅ Completed",
                                        "failed": "❌ Failed"}.get(x, x))
        with col_f2:
            if st.button("🔄 Refresh", use_container_width=True):
                _load_tasks_from_disk(st.session_state)
                st.rerun()

        filtered = task_list if status_filter == "All" else [t for t in task_list if t["status"] == status_filter]
        filtered.sort(key=lambda t: t.get("upload_time", ""), reverse=True)

        if not filtered:
            st.info("No tasks found. Upload files in the first tab to get started.")
        else:
            for t in filtered:
                status_emoji = {"pending": "⏳", "queued": "⏳", "running": "🔄", "completed": "✅", "failed": "❌"}.get(t["status"], "❓")
                status_color = {"pending": "#A0AEC0", "queued": "#3182CE", "running": "#3182CE", "completed": "#48BB78", "failed": "#E53E3E"}.get(t["status"], "#A0AEC0")

                with st.expander(f"{status_emoji} **{t['filename']}** — {t['status']} "
                                 f"({'matched_spectra' in t and t.get('matched_spectra', 0) or 0} matched)" if t["status"] == "completed"
                                 else f"{status_emoji} **{t['filename']}** — {t['status']}"):
                    # Task info
                    col_i1, col_i2 = st.columns([3, 2])
                    with col_i1:
                        st.markdown(f"""
                        <table style="font-size:0.82rem;width:100%;">
                            <tr><td style="color:#718096;">Task ID</td><td style="font-family:monospace;">{t['id']}</td></tr>
                            <tr><td style="color:#718096;">File</td><td>{t['filename']}</td></tr>
                            <tr><td style="color:#718096;">Size</td><td>{t.get('file_size_mb', 0):.2f} MB</td></tr>
                            <tr><td style="color:#718096;">Uploaded</td><td>{t.get('upload_time', '—')[:19]}</td></tr>
                        </table>
                        """, unsafe_allow_html=True)
                    with col_i2:
                        if t["status"] in ("completed", "running"):
                            mc1, mc2 = st.columns(2)
                            with mc1:
                                st.metric("Total Spectra", t.get("total_spectra", 0))
                            with mc2:
                                st.metric("Matched", f"{t.get('matched_spectra', 0)} ({t.get('match_rate', 0):.1f}%)")

                    # Progress bar for running tasks
                    if t["status"] in {"queued", "running"}:
                        st.progress(t.get("progress", 0) / 100.0, text=t.get("message", "Processing..."))

                    # Error message
                    if t["status"] == "failed":
                        st.error(f"Error: {t.get('error', 'Unknown')}")

                    # Action buttons
                    col_a1, col_a2, col_a3, col_a4 = st.columns(4)
                    with col_a1:
                        if t["status"] == "completed":
                            if st.button(f"🔍 View", key=f"view_{t['id']}"):
                                st.session_state.selected_task_id = t["id"]
                                st.info("👆 Switch to the **Results Viewer** tab to see results.")
                    with col_a2:
                        if t["status"] in ("completed", "failed", "pending"):
                            if st.button(f"🔄 Re-analyze", key=f"retry_{t['id']}"):
                                # Re-analyze with current sidebar thresholds
                                file_path = t.get("file_path", "")
                                if file_path and Path(file_path).exists():
                                    if _start_analysis_in_thread(t["id"], params={
                                            "cosine_threshold": cos_th,
                                            "final_score_threshold": fin_th,
                                            "mz_tolerance_da": mz_da,
                                            "mz_tolerance_ppm": mz_ppm,
                                            "top_n_peaks": top_n,
                                            "strict_evidence_gate": strict_gate,
                                            "enable_chimeric_deconvolution": enable_deconvolution,
                                            "enable_ion_identity_networking": enable_iimn,
                                            "isolation_window_width_da": (
                                                isolation_width if isolation_width > 0 else None
                                            ),
                                        }):
                                        st.success("Re-analysis started!")
                                        time.sleep(1)
                                        st.rerun()
                                else:
                                    st.error("Original file no longer exists.")
                    with col_a3:
                        if t["status"] == "pending":
                            if st.button(f"▶ Start", key=f"start_{t['id']}"):
                                if not st.session_state.lib_loaded:
                                    st.error("Load the library first!")
                                elif st.session_state.analysis_running:
                                    st.warning("Another analysis is running. Please wait.")
                                else:
                                    if _start_analysis_in_thread(t["id"], params={
                                        "cosine_threshold": cos_th,
                                        "final_score_threshold": fin_th,
                                        "mz_tolerance_da": mz_da,
                                        "mz_tolerance_ppm": mz_ppm,
                                        "top_n_peaks": top_n,
                                        "strict_evidence_gate": strict_gate,
                                        "enable_chimeric_deconvolution": enable_deconvolution,
                                        "enable_ion_identity_networking": enable_iimn,
                                        "isolation_window_width_da": (
                                            isolation_width if isolation_width > 0 else None
                                        ),
                                    }):
                                        st.success("Analysis started!")
                                        time.sleep(1)
                                        st.rerun()

                    with col_a4:
                        if st.button(f"🗑 Delete", key=f"del_{t['id']}"):
                            if _delete_task_from_session(t["id"]):
                                st.rerun()

    # ───────────────────────────────────────────────────────────────────
    # TAB 3: Results Viewer
    # ───────────────────────────────────────────────────────────────────
    with tab_results:
        st.markdown('<div class="section-header">Results Viewer</div>', unsafe_allow_html=True)

        # Task selector
        completed_tasks = [t for t in st.session_state.tasks.values() if t["status"] == "completed"]
        completed_tasks.sort(key=lambda t: t.get("upload_time", ""), reverse=True)

        if not completed_tasks:
            st.info("No completed tasks yet. Upload and analyze files first.")
        else:
            task_options = {t["id"]: f"{t['filename']} ({t.get('matched_spectra', 0)}/{t.get('total_spectra', 0)} matched, {t.get('match_rate', 0):.1f}%)" for t in completed_tasks}
            selected_id = st.selectbox("Select task to view",
                options=list(task_options.keys()),
                format_func=lambda x: task_options[x],
                index=0 if st.session_state.selected_task_id not in task_options else list(task_options.keys()).index(st.session_state.selected_task_id) if st.session_state.selected_task_id in task_options else 0,
            )

            # Load results on selection change
            if selected_id != st.session_state.selected_task_id:
                st.session_state.selected_task_id = selected_id
                st.session_state.selected_results = _load_task_results(selected_id)

            if st.session_state.selected_task_id and not st.session_state.selected_results:
                st.session_state.selected_results = _load_task_results(st.session_state.selected_task_id)

            results = st.session_state.selected_results or []
            if not results:
                st.warning("No results available for this task.")
            else:
                matched = [r for r in results if r.get("match")]
                matched.sort(key=lambda r: -r["match"]["final_score"])
                total = len(results)
                n_matched = len(matched)
                rate = n_matched / total * 100 if total else 0
                n_comp_high = sum(1 for r in matched if r["match"].get("confidence") == "High")

                # ── Summary Metric Cards ──
                c1, c2, c3, c4 = st.columns(4)
                metrics = [
                    ("Total", f"{total}", ""),
                    ("Matched", f"{n_matched}", f"{rate:.1f}%"),
                    (
                        "High Confidence",
                        f"{n_comp_high}",
                        "化合物鉴定置信度",
                    ),
                    ("Unmatched", f"{total - n_matched}", ""),
                ]
                for col, (label, value, sub) in zip([c1, c2, c3, c4], metrics):
                    with col:
                        st.markdown(f"""
                        <div class="metric-card fade-in">
                            <div class="metric-value">{value}</div>
                            <div class="metric-label">{label}</div>
                            {f'<div class="metric-sub">{sub}</div>' if sub else ''}
                        </div>
                        """, unsafe_allow_html=True)

                st.markdown("<br>", unsafe_allow_html=True)

                # ── Sub-tabs: Table / Detail ──
                sub_table, sub_detail = st.tabs(["📋 Table", "🔬 Detail"])

                # ── Sub-tab: Table ──
                with sub_table:
                    if matched:
                        # Filter options
                        col_f1, col_f2 = st.columns(2)
                        with col_f1:
                            show_matched_only = st.checkbox("Matched only", value=True)
                        with col_f2:
                            min_cosine = st.slider("Min cosine", 0.0, 1.0, 0.0, 0.05)

                        display_list = matched
                        if min_cosine > 0:
                            display_list = [r for r in display_list if r["match"]["cosine"] >= min_cosine]

                        rows = []
                        for i, r in enumerate(display_list):
                            m = r["match"]
                            rows.append({
                                "#": i + 1,
                                "m/z": round(r["precursor_mz"], 4),
                                "Compound": str(m["compound_name"])[:60],
                                "Formula": str(m.get("formula", "")),
                                "Cosine": round(m["cosine"], 3),
                                "Final Score": round(m["final_score"], 3),
                                "Class": str(m.get("classification", ""))[:30],
                            })
                        df = pd.DataFrame(rows)
                        st.dataframe(
                            df,
                            height=500,
                            use_container_width=True,
                            column_config={
                                "Cosine": st.column_config.NumberColumn(format="%.3f"),
                                "Final Score": st.column_config.NumberColumn(format="%.3f"),
                                "m/z": st.column_config.NumberColumn(format="%.4f"),
                            },
                            hide_index=True,
                        )
                        st.caption(f"Showing {len(display_list)} of {n_matched} matched spectra")
                    else:
                        st.info("No matches found. Try lowering the thresholds and re-analyzing.")

                # ── Sub-tab: Detail (PRESERVED FROM ORIGINAL) ──
                with sub_detail:
                    if not matched:
                        st.info("No matches to display.")
                    else:
                        idx_sel = st.selectbox(
                            "Select match to inspect",
                            range(len(matched)),
                            format_func=lambda x: f"#{x+1} · m/z {matched[x]['precursor_mz']:.2f} · {str(matched[x]['match']['compound_name'])[:50]} (final={matched[x]['match']['final_score']:.3f})"
                        )
                        m = matched[idx_sel]["match"]
                        r = matched[idx_sel]
                        # 置信度优先取 honest_calibration 化合物级(与详细报告/网页同源), 缺字段才回退旧函数
                        conf_label = m.get("confidence") or confidence(m["cosine"], m["final_score"])

                        # ── Main Card ──
                        badge_cls_map = {"High": "badge-high", "Medium": "badge-medium", "Low": "badge-low", "Unknown": "badge-unknown"}
                        badge_cls = badge_cls_map.get(conf_label, "badge-unknown")

                        st.markdown(f"""
                        <div class="result-card fade-in">
                        <div style="display:flex;justify-content:space-between;align-items:flex-start;margin-bottom:0.75rem;flex-wrap:wrap;gap:0.5rem;">
                            <div>
                                <div style="font-size:1.15rem;font-weight:600;color:#1A202C;">{str(m['compound_name'])[:70]}</div>
                                <div style="font-size:0.82rem;color:#718096;margin-top:0.25rem;">
                                    m/z {r['precursor_mz']:.4f} &nbsp;·&nbsp;
                                    Formula {m.get('formula', '—')} &nbsp;·&nbsp;
                                    Class {m.get('classification', '—')}
                                </div>
                            </div>
                            <span class="badge {badge_cls}" style="font-size:0.8rem;padding:4px 14px;">{conf_label} Confidence</span>
                        </div>
                        """, unsafe_allow_html=True)

                        # Scores row
                        score_expert_details = m.get("expert_details") or {}
                        score_assessments = (
                            score_expert_details.get("dimension_assessments") or {}
                        )
                        rule_assessed = sum(
                            str(value.get("status") or "unavailable") != "unavailable"
                            for value in score_assessments.values()
                            if isinstance(value, dict)
                        )
                        rule_supported = sum(
                            str(value.get("status") or "") == "scored"
                            and (
                                int(value.get("evidence_count") or 0) > 0
                                or float(value.get("score") or 0) > 0
                            )
                            for value in score_assessments.values()
                            if isinstance(value, dict)
                        )
                        score_items = [
                            ("Cosine", m.get("cosine", 0), 0.9, 0.8, 0.7),
                            ("Final", m.get("final_score", 0), 0.7, 0.5, 0.3),
                            ("Mass", m.get("mass_score", 0), 0.8, 0.6, 0.4),
                            ("Fragment", m.get("fragment_score", 0), None, None, None),
                            ("Intensity", m.get("intensity_score", 0), None, None, None),
                        ]
                        st.markdown('<div style="display:flex;gap:0.5rem;flex-wrap:wrap;margin-bottom:0.75rem;">', unsafe_allow_html=True)
                        for label, val, hi, med, lo in score_items:
                            if hi is not None:
                                cls = "badge-high" if val >= hi else ("badge-medium" if val >= med else ("badge-low" if val >= lo else "badge-unknown"))
                            else:
                                cls = "badge-unknown"
                            st.markdown(
                                f'<span class="badge {cls}" style="border-radius:8px;padding:0.4rem 0.75rem;text-align:center;min-width:75px;">'
                                f'<div style="font-size:0.6rem;opacity:0.75;">{label}</div>'
                                f'<div style="font-size:0.88rem;font-weight:600;">{val:.4f}</div></span>',
                                unsafe_allow_html=True
                            )
                        st.markdown(
                            '<span class="badge badge-unknown" style="border-radius:8px;'
                            'padding:0.4rem 0.75rem;text-align:center;min-width:100px;">'
                            '<div style="font-size:0.6rem;opacity:0.75;">Rule evidence</div>'
                            f'<div style="font-size:0.88rem;font-weight:600;">{rule_supported}/{rule_assessed or len(score_assessments) or 8}</div></span>',
                            unsafe_allow_html=True,
                        )
                        ppm_val = m.get("mass_ppm")
                        ppm_str = f"{ppm_val:.2f}" if ppm_val is not None else "—"
                        st.markdown(
                            f'<span class="badge badge-unknown" style="border-radius:8px;padding:0.4rem 0.75rem;text-align:center;min-width:75px;">'
                            f'<div style="font-size:0.6rem;opacity:0.75;">PPM</div>'
                            f'<div style="font-size:0.88rem;font-weight:600;">{ppm_str}</div></span>',
                            unsafe_allow_html=True
                        )
                        st.markdown('</div>', unsafe_allow_html=True)

                        # Two-column: info + structure
                        col_info, col_struct = st.columns([2, 1])
                        with col_info:
                            st.markdown(f"""
                            <table style="width:100%;font-size:0.82rem;border-collapse:collapse;">
                                <tr><td style="padding:6px 8px;color:#718096;font-weight:500;">SMILES</td>
                                    <td style="padding:6px 8px;font-family:monospace;word-break:break-all;font-size:0.78rem;">
                                        {str(m.get('smiles',''))[:120]}</td></tr>
                                <tr><td style="padding:6px 8px;color:#718096;font-weight:500;">InChIKey</td>
                                    <td style="padding:6px 8px;font-family:monospace;font-size:0.78rem;">
                                        {str(m.get('inchikey',''))[:27]}</td></tr>
                                <tr><td style="padding:6px 8px;color:#718096;font-weight:500;">Adduct</td>
                                    <td style="padding:6px 8px;">{m.get('adduct','—') or '—'}</td></tr>
                                <tr><td style="padding:6px 8px;color:#718096;font-weight:500;">Spectral ID</td>
                                    <td style="padding:6px 8px;font-family:monospace;font-size:0.78rem;">
                                        {r.get('spectrum_id','—')}</td></tr>
                            </table>
                            """, unsafe_allow_html=True)

                        with col_struct:
                            smi = m.get("smiles", "")
                            if smi:
                                from web.structure import smiles_to_image_base64
                                img_b64 = smiles_to_image_base64(smi, size=(350, 280))
                                if img_b64:
                                    st.markdown(
                                        f'<div style="text-align:center;padding:0.5rem;">'
                                        f'<img src="data:image/png;base64,{img_b64}" '
                                        f'style="max-width:200px;border:1px solid #E2E8F0;border-radius:6px;padding:6px;">'
                                        f'<div style="font-size:0.7rem;color:#718096;margin-top:4px;">2D Structure</div>'
                                        f'</div>',
                                        unsafe_allow_html=True
                                    )
                                else:
                                    st.caption("Structure not available")
                            else:
                                st.caption("No SMILES available")

                        st.markdown("</div>", unsafe_allow_html=True)

                        # ── Structured rule-evidence assessment ──
                        expert_details = m.get("expert_details")
                        if expert_details:
                            with st.expander("🧪 Structured Rule-Evidence Assessment", expanded=False):
                                ed = expert_details
                                reasoning = ed.get("chemical_reasoning", "No reasoning available")
                                assessments = ed.get("dimension_assessments") or {}
                                ec1, ec2, ec3, ec4, ec5, ec6 = st.columns(6)
                                expert_sub_scores = [
                                    ("neutral_loss", "Neutral Loss", ed.get("neutral_loss_score", 0)),
                                    ("diagnostic_ion", "Diagnostic Ion", ed.get("diagnostic_ion_score", 0)),
                                    ("rule_coverage", "Rule Coverage", ed.get("rule_coverage_score", 0)),
                                    ("fg_signature", "FG Signature", ed.get("fg_signature_score", 0)),
                                    ("structure_fragment", "Structure-Frag", ed.get("structure_fragment_score", 0)),
                                    ("bond_energy", "Bond Energy", ed.get("bond_energy_score", 0)),
                                ]
                                for col, (dimension_key, label, val) in zip([ec1, ec2, ec3, ec4, ec5, ec6], expert_sub_scores):
                                    with col:
                                        assessment = assessments.get(dimension_key) or {}
                                        status = assessment.get("status") or (
                                            "scored" if val > 0 else "no_evidence"
                                        )
                                        display_value = (
                                            "N/A" if status == "not_applicable"
                                            else "Unavailable" if status == "unavailable"
                                            else "No support" if status == "no_evidence"
                                            else f"{val:.3f}"
                                        )
                                        color = "#A0AEC0" if status != "scored" else "#48BB78" if val >= 0.7 else "#ECC94B" if val >= 0.45 else "#3182CE"
                                        st.markdown(f"""
                                        <div style="text-align:center;padding:0.5rem;background:#F7FAFC;
                                                    border-radius:8px;border:1px solid #E2E8F0;">
                                            <div style="font-size:0.82rem;font-weight:700;color:{color};">{display_value}</div>
                                            <div style="font-size:0.6rem;color:#718096;">{label}</div>
                                        </div>
                                        """, unsafe_allow_html=True)
                                st.markdown(f"""
                                <div style="margin-top:0.75rem;padding:0.75rem;background:#F0FFF4;
                                            border-radius:8px;border:1px solid #C6F6D5;">
                                    <strong>🧠 Chemical Reasoning:</strong><br>
                                    <span style="font-size:0.85rem;color:#276749;">{reasoning}</span>
                                </div>
                                """, unsafe_allow_html=True)
                                rc = ed.get("rule_coverage_analysis", {})
                                if rc.get("covered_rules"):
                                    st.markdown("**📏 Fragmentation Rules Activated:**")
                                    for rule in rc["covered_rules"][:5]:
                                        st.markdown(f"- **{rule['name']}** ({rule['mechanism_type']}, confidence={rule['confidence_score']:.0%})")
                                di = ed.get("diagnostic_ion_analysis", {})
                                matched_ions = di.get("matched_ions", [])
                                if matched_ions:
                                    st.markdown("**🔬 Diagnostic Ions Detected:**")
                                    for ion in matched_ions[:5]:
                                        st.markdown(f"- m/z {ion['observed_mz']:.4f} → **{ion['ion_name']}** ({ion.get('structural_origin', '')[:60]}...)")

                        # ── Spectrum Analysis ──
                        qpeaks = r.get("query_peaks", []) or r.get("peaks", [])
                        lpeaks = m.get("lib_peaks", [])
                        has_lib = bool(lpeaks)

                        if qpeaks or lpeaks:
                            st.markdown("---")
                            st.markdown('<div class="section-header">Spectrum Analysis</div>', unsafe_allow_html=True)
                            matched_frags = find_matching_fragments(qpeaks, lpeaks, tol=0.02) if has_lib else []

                            cq, cl, cm = st.columns(3)
                            with cq:
                                st.markdown(f"""<div class="metric-card fade-in">
                                    <div class="metric-value" style="color:#3182CE;">{len(qpeaks)}</div>
                                    <div class="metric-label">Query Peaks</div></div>""", unsafe_allow_html=True)
                            with cl:
                                st.markdown(f"""<div class="metric-card fade-in">
                                    <div class="metric-value" style="color:#E53E3E;">{len(lpeaks)}</div>
                                    <div class="metric-label">Library Peaks</div></div>""", unsafe_allow_html=True)
                            with cm:
                                st.markdown(f"""<div class="metric-card fade-in">
                                    <div class="metric-value" style="color:#48BB78;">{len(matched_frags)}</div>
                                    <div class="metric-label">Matched Fragments</div></div>""", unsafe_allow_html=True)

                            # ── Plotly Spectrum ──
                            try:
                                import plotly.graph_objects as go
                                fig = go.Figure()
                                if qpeaks:
                                    q_sorted = sorted(qpeaks, key=lambda x: x[1], reverse=True)[:60]
                                    fig.add_trace(go.Bar(x=[p[0] for p in q_sorted], y=[p[1] for p in q_sorted],
                                        name="Query", marker=dict(color="#3182CE", line=dict(width=0)),
                                        opacity=0.75, width=1.5,
                                        hovertemplate="m/z: %{x:.2f}<br>Intensity: %{y:.1f}<extra>Query</extra>"))
                                if lpeaks:
                                    l_sorted = sorted(lpeaks, key=lambda x: x[1], reverse=True)[:60]
                                    fig.add_trace(go.Bar(x=[p[0] for p in l_sorted], y=[-p[1] for p in l_sorted],
                                        name="Library", marker=dict(color="#E53E3E", line=dict(width=0)),
                                        opacity=0.75, width=1.5,
                                        hovertemplate="m/z: %{x:.2f}<br>Intensity: %{customdata:.1f}<extra>Library</extra>",
                                        customdata=[abs(p[1]) for p in l_sorted]))
                                fig.update_layout(barmode="overlay",
                                    title=dict(text="Spectral Comparison", font=dict(size=14), x=0.5),
                                    xaxis=dict(title="m/z", gridcolor="#EDF2F7"), yaxis=dict(title="Intensity", gridcolor="#EDF2F7", zeroline=True, zerolinecolor="#CBD5E0"),
                                    plot_bgcolor="white", paper_bgcolor="white", hovermode="x unified",
                                    height=380, margin=dict(l=40, r=20, t=50, b=40),
                                    legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1))
                                st.plotly_chart(fig, use_container_width=True, config={"displayModeBar": True, "modeBarButtonsToRemove": ["lasso2d", "select2d"], "displaylogo": False})
                            except ImportError:
                                pass

                            # ── Matched Fragments ──
                            if matched_frags:
                                with st.expander(f"Matched Fragment Ions ({len(matched_frags)} pairs)", expanded=True):
                                    frag_data = [{"Query m/z": f"{mf['query_mz']:.4f}", "Query Int": f"{mf['query_intensity']:.1f}",
                                                  "Lib m/z": f"{mf['lib_mz']:.4f}", "Lib Int": f"{mf['lib_intensity']:.1f}",
                                                  "Δ m/z": f"{mf['diff']:.4f}"} for mf in matched_frags[:30]]
                                    st.dataframe(pd.DataFrame(frag_data), height=250, use_container_width=True, hide_index=True)
                            elif has_lib:
                                st.info("No matching fragment ion pairs found.")

                            # ── Peak Tables ──
                            col_q, col_l = st.columns(2)
                            with col_q:
                                st.markdown('<div style="font-size:0.85rem;font-weight:600;color:#3182CE;margin-bottom:0.4rem;">Query Peaks (Top 20)</div>', unsafe_allow_html=True)
                                q_data = [{"m/z": f"{p[0]:.4f}", "Intensity": f"{p[1]:.1f}"} for p in sorted(qpeaks, key=lambda x: -x[1])[:20]]
                                st.dataframe(pd.DataFrame(q_data), height=250, use_container_width=True, hide_index=True)
                            with col_l:
                                if has_lib:
                                    st.markdown('<div style="font-size:0.85rem;font-weight:600;color:#E53E3E;margin-bottom:0.4rem;">Library Peaks (Top 20)</div>', unsafe_allow_html=True)
                                    l_data = [{"m/z": f"{p[0]:.4f}", "Intensity": f"{p[1]:.1f}"} for p in sorted(lpeaks, key=lambda x: -x[1])[:20]]
                                    st.dataframe(pd.DataFrame(l_data), height=250, use_container_width=True, hide_index=True)

    # ───────────────────────────────────────────────────────────────────
    # TAB 4: Export Center
    # ───────────────────────────────────────────────────────────────────
    with tab_export:
        st.markdown('<div class="section-header">Export Center</div>', unsafe_allow_html=True)

        # Task selector for export
        export_tasks = [t for t in st.session_state.tasks.values() if t["status"] == "completed"]
        export_tasks.sort(key=lambda t: t.get("upload_time", ""), reverse=True)

        if not export_tasks:
            st.info("No completed tasks to export. Run analyses first.")
        else:
            export_opts = {t["id"]: f"{t['filename']} ({t.get('matched_spectra', 0)}/{t.get('total_spectra', 0)} matched)" for t in export_tasks}
            export_id = st.selectbox("Select task to export", options=list(export_opts.keys()), format_func=lambda x: export_opts[x])

            if export_id:
                export_task = st.session_state.tasks.get(export_id, {})
                export_results = _load_task_results(export_id)
                export_matched = [r for r in export_results if r.get("match")]

                if not export_matched:
                    st.warning("No matched results to export for this task.")
                else:
                    # ── Export Format Grid ──
                    st.markdown('<div style="font-size:0.95rem;font-weight:600;margin-bottom:0.75rem;">Single-Task Exports</div>', unsafe_allow_html=True)

                    ec1, ec2, ec3 = st.columns(3)
                    ec4, ec5, ec6 = st.columns(3)

                    # Results CSV
                    with ec1:
                        st.markdown("""<div class="export-card fade-in">""", unsafe_allow_html=True)
                        buf1 = io.StringIO()
                        w1 = csv.writer(buf1)
                        w1.writerow(["Spectrum", "Precursor_mz", "Match_Name", "Match_SMILES", "Match_InChIKey", "Match_Formula",
                                      "Final_Score", "Cosine_Score", "Confidence", "Class", "MS_MS_Fragments", "Fragment_Count"])
                        for r in export_matched:
                            m = r["match"]
                            msms_s = _msms_fragments_str(r.get("query_peaks", []), m.get("lib_peaks", []))
                            frag_n = len(msms_s.split("; ")) if msms_s else 0
                            w1.writerow([r["spectrum_id"], f"{r['precursor_mz']:.4f}", m["compound_name"],
                                         m.get("smiles", ""), m.get("inchikey", ""), m.get("formula", ""),
                                         f"{m['final_score']:.4f}", f"{m['cosine']:.4f}", m.get("confidence", ""),
                                         m.get("classification", ""), msms_s, frag_n])
                        fname_stem = Path(export_task.get("filename", "analysis")).stem
                        st.markdown('<div class="export-card-title">📊 Results CSV</div>', unsafe_allow_html=True)
                        st.markdown('<div class="export-card-desc">Match data with scores</div>', unsafe_allow_html=True)
                        st.download_button("Download", buf1.getvalue(), f"{fname_stem}_results.csv", "text/csv", use_container_width=True, key="dl_csv_results")
                        st.markdown("</div>", unsafe_allow_html=True)

                    # Report CSV
                    with ec2:
                        st.markdown("""<div class="export-card fade-in">""", unsafe_allow_html=True)
                        buf2 = io.StringIO()
                        w2 = csv.writer(buf2)
                        w2.writerow(["Spectrum", "Precursor_mz", "Peaks_Count", "Compound_Name", "Classification",
                                      "Cosine_Score", "Final_Score", "Mass_Score", "Mass_PPM", "Fragment_Score",
                                      "Intensity_Score", "InChIKey", "SMILES", "MS_MS_Fragments", "Fragment_Count"])
                        for r in export_matched:
                            m = r["match"]
                            ppm = m.get("mass_ppm")
                            msms_s = _msms_fragments_str(r.get("query_peaks", []), m.get("lib_peaks", []))
                            frag_n = len(msms_s.split("; ")) if msms_s else 0
                            w2.writerow([r["spectrum_id"], f"{r['precursor_mz']:.4f}", r.get("peaks_count", 0),
                                         m["compound_name"], m.get("classification", ""),
                                         f"{m['cosine']:.4f}", f"{m['final_score']:.4f}", f"{m.get('mass_score', 0):.4f}",
                                         f"{ppm:.2f}" if ppm else "", f"{m.get('fragment_score', 0):.4f}",
                                         f"{m.get('intensity_score', 0):.4f}", m.get("inchikey", ""), m.get("smiles", ""),
                                         msms_s, frag_n])
                        st.markdown('<div class="export-card-title">📈 Report CSV</div>', unsafe_allow_html=True)
                        st.markdown('<div class="export-card-desc">Full scores & metadata</div>', unsafe_allow_html=True)
                        st.download_button("Download", buf2.getvalue(), f"{fname_stem}_report.csv", "text/csv", use_container_width=True, key="dl_csv_report")
                        st.markdown("</div>", unsafe_allow_html=True)

                    # Excel
                    with ec3:
                        st.markdown("""<div class="export-card fade-in">""", unsafe_allow_html=True)
                        try:
                            excel_bytes = _export_excel(export_task, export_results)
                            st.markdown('<div class="export-card-title">📗 Excel (.xlsx)</div>', unsafe_allow_html=True)
                            st.markdown('<div class="export-card-desc">3 sheets: Matched, Unmatched, Summary</div>', unsafe_allow_html=True)
                            st.download_button("Download", excel_bytes, f"{fname_stem}_analysis.xlsx",
                                "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                                use_container_width=True, key="dl_excel")
                        except Exception as e:
                            st.markdown('<div class="export-card-title">📗 Excel</div>', unsafe_allow_html=True)
                            st.error(f"openpyxl needed: {e}")
                        st.markdown("</div>", unsafe_allow_html=True)

                    # Standard HTML
                    with ec4:
                        st.markdown("""<div class="export-card fade-in">""", unsafe_allow_html=True)
                        st.markdown('<div class="export-card-title">🌐 Standard HTML</div>', unsafe_allow_html=True)
                        st.markdown('<div class="export-card-desc">Interactive table with structures</div>', unsafe_allow_html=True)
                        if st.button("Generate", use_container_width=True, key="btn_html_std"):
                            html_results = []
                            for r in export_results:
                                m = r.get("match")
                                html_results.append({
                                    "precursor_mz": r["precursor_mz"], "spectrum_id": r.get("spectrum_id", ""),
                                    "peaks_count": r.get("peaks_count", 0), "retention_time": r.get("retention_time"),
                                    "retention_time_unit": r.get("retention_time_unit"),
                                    "precursor_scan_ref": r.get("precursor_scan_ref"),
                                    "query_peaks": r.get("query_peaks", []),
                                    "matches": [m] if m else [], "best_cosine": m["cosine"] if m else 0.0,
                                })
                            with tempfile.NamedTemporaryFile(suffix=".html", delete=False, mode="w") as tmp:
                                from web.report import generate_html
                                generate_html(html_results, tmp.name, fname_stem)
                            with open(tmp.name) as f:
                                st.download_button("Download HTML", f.read(), f"{fname_stem}_report.html", "text/html",
                                    use_container_width=True, key="dl_html_std")
                        st.markdown("</div>", unsafe_allow_html=True)

                    # Detailed HTML
                    with ec5:
                        st.markdown("""<div class="export-card fade-in">""", unsafe_allow_html=True)
                        st.markdown('<div class="export-card-title">🔬 Detailed HTML</div>', unsafe_allow_html=True)
                        st.markdown('<div class="export-card-desc">Dual spectra + fragment matching</div>', unsafe_allow_html=True)
                        if st.button("Generate", use_container_width=True, key="btn_html_det"):
                            html_results = []
                            for r in export_results:
                                m = r.get("match")
                                if m:
                                    m["lib_peaks"] = r.get("lib_peaks", [])
                                    ed = m.get("expert_details") or {}
                                    for ek in ["neutral_loss_score", "diagnostic_ion_score", "rule_coverage_score",
                                                "fg_signature_score", "structure_fragment_score", "bond_energy_score",
                                                "feature_library_score", "isomer_differentiation_score", "chemical_reasoning"]:
                                        if ek in ed:
                                            m[ek] = ed[ek]
                                html_results.append({
                                    "precursor_mz": r["precursor_mz"], "spectrum_id": r.get("spectrum_id", ""),
                                    "peaks_count": r.get("peaks_count", 0), "retention_time": r.get("retention_time"),
                                    "retention_time_unit": r.get("retention_time_unit"),
                                    "precursor_scan_ref": r.get("precursor_scan_ref"),
                                    "query_peaks": r.get("query_peaks", []),
                                    "matches": [m] if m else [], "best_cosine": m["cosine"] if m else 0.0,
                                })
                            rich_html = _gen_rich(
                                html_results,
                                source_name=fname_stem,
                                governance_blocked_candidates=(
                                    project_governance_blocked_candidates(export_results)
                                ),
                            )
                            st.download_button("Download HTML", rich_html, f"{fname_stem}_detailed.html", "text/html",
                                use_container_width=True, key="dl_html_det")
                        st.markdown("</div>", unsafe_allow_html=True)

                    # Structure PNG
                    with ec6:
                        st.markdown("""<div class="export-card fade-in">""", unsafe_allow_html=True)
                        st.markdown('<div class="export-card-title">🧬 Structure Images</div>', unsafe_allow_html=True)
                        st.markdown('<div class="export-card-desc">Download compound structure images</div>', unsafe_allow_html=True)
                        struct_idx = st.number_input("Match #", 0, len(export_matched) - 1, 0, key="struct_idx")
                        if struct_idx < len(export_matched):
                            smi = export_matched[struct_idx]["match"].get("smiles", "")
                            if smi:
                                from web.structure import smiles_to_image_base64
                                img_bytes = smiles_to_image_base64(smi, size=(600, 450))
                                if img_bytes:
                                    st.image(base64.b64decode(img_bytes), use_column_width=False)
                                    st.download_button("Download PNG", base64.b64decode(img_bytes),
                                        f"structure_{struct_idx+1}.png", "image/png", use_container_width=True, key="dl_struct")
                        st.markdown("</div>", unsafe_allow_html=True)

                    # ── Batch ZIP Export ──
                    st.markdown("---")
                    st.markdown('<div style="font-size:0.95rem;font-weight:600;margin:1rem 0 0.75rem;">📦 Batch Export</div>', unsafe_allow_html=True)
                    st.markdown('<div style="font-size:0.85rem;color:#718096;margin-bottom:0.75rem;">Select multiple completed tasks and download as a ZIP file.</div>', unsafe_allow_html=True)

                    batch_ids = []
                    for t in export_tasks:
                        checked = st.checkbox(
                            f"{t['filename']} ({t.get('matched_spectra', 0)}/{t.get('total_spectra', 0)} matched)",
                            value=False, key=f"batch_{t['id']}"
                        )
                        if checked:
                            batch_ids.append(t["id"])

                    if batch_ids:
                        if st.button(f"📦 Export {len(batch_ids)} Tasks as ZIP", use_container_width=True, type="primary"):
                            zip_bytes = _export_batch_zip(batch_ids, st.session_state)
                            st.download_button("Download ZIP", zip_bytes, "batch_export.zip", "application/zip",
                                use_container_width=True, key="dl_batch_zip")
