"""质谱库搜索与匹配引擎"""
import json
import logging
import math
import os
import re
import sqlite3
from statistics import median
from pathlib import Path

import numpy as np

from utils.adduct import (
    NEUTRAL_SHIFTS,
    adduct_polarity as _adduct_polarity,
    ion_mz as _calc_ion_mz,
    neutral_mass as _calc_neutral_mass,
)

from config.settings import (
    DB_PATH, MZ_TOLERANCE_DA, MZ_TOLERANCE_PPM, COSINE_THRESHOLD,
    FINAL_SCORE_THRESHOLD, MZ_BUCKET_DALTON, MAX_CANDIDATES_PER_QUERY,
    WEIGHT_COSINE, WEIGHT_ENTROPY, WEIGHT_MASS, WEIGHT_FRAGMENT, WEIGHT_INTENSITY,
    WEIGHT_EXPERT, TOP_N_PEAKS, REL_THRESHOLD,
    EXPERT_SYSTEM_ENABLED, EXPERT_LOG_ENABLED,
    EXPERT_RANK_INFLUENCE_ENABLED,
    USE_EXPERT_V3, TOP_K_CANDIDATES,
    LOAD_NAME_ONLY_LEADS,
    USE_NEUTRAL_MASS_MATCHING, NEUTRAL_MASS_TOL_DA,
    ENFORCE_POLARITY_FILTER,
    USE_CONSISTENCY_FILTER, CONSISTENCY_MODE, CONSISTENCY_PENALTY,
    CONSISTENCY_MIN_DL_CONF, CONSISTENCY_MIN_CLASS_ACC, CONSISTENCY_RELIABLE_CLASSES,
    USE_PATHWAY_ALIGN, USE_PATHWAY_ALIGN_SHADOW,
    PATHWAY_ALIGN_MZ_TOL, PATHWAY_ALIGN_MIN_SUPPORT,
    PATHWAY_ALIGN_REQUIRE_RECURRING, USE_PATHWAY_RERANK, PATHWAY_RERANK_STRENGTH,
    USE_DIAGNOSTIC_RESCUE, DIAGNOSTIC_RESCUE_MZ_TOL, DIAGNOSTIC_RESCUE_INTENSITY_FLOOR,
    USE_RT_STANDARD_CONFIDENCE,
    USE_MULTIMETRIC_SHADOW,
    MULTIMETRIC_FRAGMENT_TOL_DA, MULTIMETRIC_STRICT_FRAGMENT_TOL_DA,
    USE_ISOMER_FRAGMENT_RERANK, USE_ISOMER_FRAGMENT_RANK_OVERRIDE,
    USE_ISOPAIR_MS,
    ISOMER_FRAGMENT_RERANK_WEIGHT,
    ISOMER_FRAGMENT_MIN_REL_INTENSITY,
    ISOMER_FRAGMENT_MIN_SHARED_CANDIDATES, ISOMER_FRAGMENT_MAX_CE_DELTA,
    ISOMER_DIAGNOSTIC_QUERY_TOP_N,
    ISOMER_DIAGNOSTIC_QUERY_MIN_REL_INTENSITY,
    ISOMER_DIAGNOSTIC_MZ_TOL_DA,
    ISOMER_DIAGNOSTIC_MIN_REFERENCE_REPLICATES,
    ISOMER_DIAGNOSTIC_MIN_PROFILE_PREVALENCE,
    ISOMER_TIE_FINAL_EPS,
    USE_ADAPTIVE_PEAK_SEARCH, USE_BROAD_PEAK_RETRY,
    ADAPTIVE_TOP_N_PEAKS, ADAPTIVE_REL_THRESHOLD,
    ADAPTIVE_RETRY_MAX_COSINE, ADAPTIVE_MIN_COSINE_GAIN, ADAPTIVE_MAX_MASS_PPM,
    USE_LOW_MZ_INTERFERENCE_RETRY, LOW_MZ_INTERFERENCE_MAX_MZ,
    LOW_MZ_INTERFERENCE_DOMINANCE,
    USE_STRICT_EVIDENCE_GATE, STRICT_MIN_FINAL, STRICT_MIN_ENTROPY,
    STRICT_MIN_REVERSE_MATCH, STRICT_MAX_MASS_PPM,
    USE_CE_PROFILE_SHADOW,
    USE_CHIMERIC_DECONVOLUTION,
    CHIMERIC_DECONVOLUTION_L1_PENALTY,
    CHIMERIC_DECONVOLUTION_MAX_COMPONENTS,
    CHIMERIC_DECONVOLUTION_MIN_COSINE_GAIN,
    CHIMERIC_DECONVOLUTION_MIN_FINAL_GAIN,
    CHIMERIC_DECONVOLUTION_MAX_FINAL_GAIN,
    CHIMERIC_RUN_COMPONENT_MIN_PURITY,
    CHIMERIC_RUN_COMPONENT_MIN_TIC,
)
from matching.reference_mass import (
    DERIVED_NEUTRAL_MASS_SOURCE, verified_neutral_mass_correction,
)
from matching.scores import (
    batch_modified_cosine_scores, cosine_score, entropy_score, mass_score, mass_ppm, fragment_score,
    intensity_score, composite_score, enhanced_composite_score, reverse_match_ratio,
    mass_error_penalty,
    confidence as score_confidence,
)
from matching.compound_names import (
    is_readable_compound_name,
    resolve_compound_name,
)
from matching.compound_annotations import annotate_compound_dict
from matching.honest_calibration import (
    candidate_comparison_view,
    is_spectrally_competitive_isomer,
)
from matching.multimetric_shadow import (
    PEAK_VIEW_PRESETS,
    build_multimetric_shadow,
    build_peak_view_shadow,
)
from matching.conditioned_multilibrary import (
    normalize_source_family,
    reference_content_fingerprint,
)
from matching.feature_fragments import (
    FeatureFragmentRuntimeGateError,
    evaluate_feature_fragment_evidence,
    feature_fragment_score_standard,
    load_catalog as load_feature_fragment_catalog,
    load_runtime_gate as load_feature_fragment_runtime_gate,
    resolve_instrument_domain,
    runtime_gate_audit,
)
from matching.isomer_resolution import summarize_isomer_relationships
from matching.rule_quality import collision_energy_value_unit
from matching.library_quality import (
    PRESENCE_ONLY_SOURCE_TIER,
    compound_identity_key,
    formula_exact_mass,
    normalize_inchikey,
    resolve_structure_formula,
    source_tier_is_searchable,
    validate_structure_metadata,
)
from matching.library_polarity import classify_spectrum_polarity
from matching.classification_provenance import trusted_analysis_class
from matching.retrieval_policy import build_retrieval_audit
from utils.spectrum import (
    PeakProcessingConfig,
    clean_peaks,
    process_peaks,
    resolve_peak_processing_config,
)
from utils.instrument_metadata import normalize_instrument_family
from matching.source_ms_stage import source_ms_stage, source_ms2_pair_status
from .blank_subtraction import subtract_blank, BLANK_SUBTRACTION_AVAILABLE


_FEATURE_FRAGMENT_CATALOG = None
_FEATURE_FRAGMENT_RUNTIME_GATE = None


def _feature_fragment_catalog():
    global _FEATURE_FRAGMENT_CATALOG
    if _FEATURE_FRAGMENT_CATALOG is None:
        _FEATURE_FRAGMENT_CATALOG = load_feature_fragment_catalog()
    return _FEATURE_FRAGMENT_CATALOG


def _feature_fragment_runtime_gate():
    global _FEATURE_FRAGMENT_RUNTIME_GATE
    if _FEATURE_FRAGMENT_RUNTIME_GATE is None:
        _FEATURE_FRAGMENT_RUNTIME_GATE = load_feature_fragment_runtime_gate()
    return _FEATURE_FRAGMENT_RUNTIME_GATE


def _feature_fragment_unavailable(
        reason, *, acquisition_domain=None, runtime_export_gate=None):
    output = {
        'module': 'v15_feature_fragments',
        'status': 'unavailable',
        'reason': str(reason),
        'evidence_stage': 'msms_rule_support',
        'runtime_influence': False,
        'identity_rank_influence': False,
        'score_contribution': 0.0,
        'exact_identity_confirmed': False,
        'standard_confirmed': False,
        'isomer_discriminating': False,
        'confirmation_status': 'not_confirmed_by_feature_fragments',
        'feature_fragment_score_standard': feature_fragment_score_standard(),
        'rule_evidence': [],
        'informative_independent_group_count': 0,
        'present_independent_fragment_count': 0,
        'eligible_absence_independent_group_count': 0,
        'cumulative_evidence_by_independent_count': [],
        'combined_log_likelihood_ratio': 0.0,
        'combined_score_nats': 0.0,
        'combined_score_bits': 0.0,
        'combined_score_decibans': 0.0,
        'likelihood_ratio': 1.0,
        'acquisition_domain': acquisition_domain or {},
    }
    if runtime_export_gate is not None:
        output['runtime_export_gate'] = runtime_export_gate
    return output


def _evaluate_feature_fragment_shadow(
        *, query_peaks, query_pmz, candidate_class,
        classification_reliable, classification_source,
        query_polarity, query_adduct, query_adduct_reliable,
        query_charge, query_collision_energy,
        query_collision_energy_mode, query_instrument_type,
        query_instrument, query_activation,
        query_instrument_domain_id, query_collision_cell_region,
        query_spectrum_representation, query_raw_peaks_complete,
        query_isolation_evidence):
    """Attach conditional class evidence without touching any search score."""
    try:
        runtime_gate = _feature_fragment_runtime_gate()
        gate_audit = runtime_gate_audit(runtime_gate)
        if runtime_gate.get('runtime_evaluation_enabled') is not True:
            return _feature_fragment_unavailable(
                runtime_gate.get('reason')
                or 'no_feature_fragment_rule_passed_sealed_runtime_gate',
                runtime_export_gate=gate_audit,
            )
        catalog = _feature_fragment_catalog()
        ce_value, ce_unit = collision_energy_value_unit(query_collision_energy)
        acquisition = resolve_instrument_domain(
            catalog=catalog,
            instrument=query_instrument,
            instrument_type=query_instrument_type,
            activation=query_activation,
            collision_energy_unit=ce_unit,
            collision_energy_mode=query_collision_energy_mode,
            collision_cell_region=query_collision_cell_region,
            explicit_domain_id=query_instrument_domain_id,
        )
        if acquisition.get('status') != 'resolved':
            return _feature_fragment_unavailable(
                acquisition.get('reason') or 'instrument_domain_unresolved',
                acquisition_domain=acquisition,
            )
        if ce_value is None:
            return _feature_fragment_unavailable(
                'collision_energy_value_missing', acquisition_domain=acquisition
            )
        polarity = str(query_polarity or '').strip().upper()
        polarity = 'positive' if polarity == 'POS' else (
            'negative' if polarity == 'NEG' else ''
        )
        isolation = query_isolation_evidence or {}
        normalized_query = [
            peak for peak in query_peaks or []
            if isinstance(peak, (list, tuple)) and len(peak) >= 2
        ]
        quality = {
            'complete_mass_range': query_raw_peaks_complete is True,
            'centroided': str(query_spectrum_representation or '').lower() == 'centroid',
            'sufficient_signal': len(normalized_query) >= 3,
            'precursor_isolated': str(isolation.get('status') or '').lower() == 'clean',
            'collision_energy_verified': bool(ce_unit and query_collision_energy_mode),
            'instrument_domain_verified': True,
        }
        evidence = evaluate_feature_fragment_evidence(
            peaks=query_peaks,
            precursor_mz=query_pmz,
            context={
                'chemical_class': str(candidate_class or ''),
                'classification_reliable': classification_reliable is True,
                'classification_source': str(classification_source or ''),
                'polarity': polarity,
                # A library-record adduct is not substituted for missing query
                # metadata: that would feed the candidate premise back into the
                # claimed query-spectrum support.
                'adduct': str(query_adduct or '')
                if query_adduct_reliable is True else '',
                'charge': query_charge,
                'instrument_domain_id': acquisition['instrument_domain_id'],
                'activation': acquisition['activation'],
                'collision_energy': {
                    'value': ce_value,
                    'unit': ce_unit,
                    'mode': acquisition['collision_energy_mode'],
                },
                'quality': quality,
            },
            catalog=catalog,
        )
        evidence['acquisition_domain'] = acquisition
        evidence['classification_premise'] = {
            'chemical_class': str(candidate_class or ''),
            'reliable': classification_reliable is True,
            'source': str(classification_source or ''),
        }
        evidence['quality_gate'] = quality
        return evidence
    except FeatureFragmentRuntimeGateError as exc:
        logging.getLogger(__name__).error(
            'V15 feature-fragment runtime gate failed closed: %s', exc
        )
        return _feature_fragment_unavailable(
            'feature_fragment_runtime_gate_invalid',
            runtime_export_gate={
                'deployment_state': 'fail_closed',
                'runtime_evaluation_enabled': False,
                'runtime_adapter_available': False,
                'exported_rule_count': 0,
                'research_catalog_evaluated': False,
                'validation_error': str(exc),
            },
        )
    except Exception as exc:
        logging.getLogger(__name__).warning(
            'V15 feature-fragment shadow evaluation failed: %s', exc
        )
        return _feature_fragment_unavailable('feature_fragment_evaluation_error')

try:
    from .chimeric_deconvolution import deconvolve_chimeric_spectrum
    CHIMERIC_DECONVOLUTION_AVAILABLE = True
except Exception:
    CHIMERIC_DECONVOLUTION_AVAILABLE = False

try:
    from .model_predictor import predict_from_spectrum, structure_similarity
    MODEL_PREDICTOR_AVAILABLE = True
except ImportError:
    MODEL_PREDICTOR_AVAILABLE = False

# 专家系统模块(可选导入)
try:
    from .expert_scorer import (
        ChemicalExpertScorer,
        expert_evaluate,
        finalize_expert_evidence,
    )
    EXPERT_SYSTEM_AVAILABLE = True
except ImportError:
    EXPERT_SYSTEM_AVAILABLE = False

# v3 Phase 4: v3 KB 专家加分模块 (可选导入; 失败则降级为纯 v2 专家)
try:
    from .expert_v3 import score_and_reason as _v3_score_and_reason
    EXPERT_V3_AVAILABLE = True
except Exception:
    EXPERT_V3_AVAILABLE = False

# Compound/background contrastive rules are exposed as shadow evidence only.
try:
    from .expert_contrastive import (
        evaluate_contrastive_shadow as _contrastive_shadow,
    )
    EXPERT_CONTRASTIVE_AVAILABLE = True
except Exception:
    EXPERT_CONTRASTIVE_AVAILABLE = False

# P1: 分类一致性过滤模块 (可选导入; 失败则降级为不过滤)
try:
    from .consistency_filter import resolve_query_class as _resolve_query_class, is_conflict as _is_class_conflict
    CONSISTENCY_FILTER_AVAILABLE = True
except Exception:
    CONSISTENCY_FILTER_AVAILABLE = False

# P-A: 化合物专属裂解通路对齐 (可选导入; 失败则降级为不用通路信号)
try:
    from .pathway_align import parse_pathway as _parse_pathway, pathway_align_score as _pathway_align_score
    PATHWAY_ALIGN_AVAILABLE = True
except Exception:
    PATHWAY_ALIGN_AVAILABLE = False

# P-B: 诊断碎片抢救 (可选导入; 失败则降级为不抢救)
try:
    from .diagnostic_rescue import rescue_diagnostic_peaks as _rescue_diagnostic_peaks
    DIAGNOSTIC_RESCUE_AVAILABLE = True
except Exception:
    DIAGNOSTIC_RESCUE_AVAILABLE = False

# v5: local RT standards (same-method retention time evidence).
try:
    from .rt_reference import match_rt_evidence as _match_rt_evidence
    RT_REFERENCE_AVAILABLE = True
except Exception:
    RT_REFERENCE_AVAILABLE = False

# Multi-source CE profiles are a read-only, shadow evidence layer.
try:
    from .ce_profiles import lookup_ce_profile_evidence as _lookup_ce_profile
    CE_PROFILE_AVAILABLE = True
except Exception:
    CE_PROFILE_AVAILABLE = False

try:
    from .isopair_ranker import apply_isopair_ranker as _apply_isopair_ranker
    ISOPAIR_RUNTIME_AVAILABLE = True
except Exception:
    ISOPAIR_RUNTIME_AVAILABLE = False

logger = logging.getLogger(__name__)

_POLARITY_ONLY_ADDUCTS = {"POSITIVE", "POS", "+", "NEGATIVE", "NEG", "-"}
_SUPPORTED_SPECIFIC_ADDUCTS = {
    re.sub(r'\s+', '', key).upper() for key in NEUTRAL_SHIFTS
}
_ADDUCT_IDENTITY = {
    re.sub(r'\s+', '', key).upper(): (
        round(float(shift), 6),
        abs(int(charge)),
        _adduct_polarity(key),
    )
    for key, (shift, charge) in NEUTRAL_SHIFTS.items()
}
_COMMON_STRUCTURE_ADDUCT_HYPOTHESES = {
    'POS': (
        '[M]+', '[M+H]+', '[M+Na]+', '[M+K]+', '[M+NH4]+',
        '[M+H-H2O]+', '[M+H-NH3]+',
    ),
    'NEG': (
        '[M]-', '[M-H]-', '[M+Cl]-', '[M+FA-H]-',
        '[M+Hac-H]-', '[M-H-H2O]-',
    ),
}

# MS2Classifier 深度学习分类器 — 懒加载（首次调用时初始化，避免拖慢启动）
_clf_wrapper = None
_clf_loaded = False

def _ensure_clf():
    """按需加载 MS2Classifier，仅首次调用时初始化"""
    global _clf_wrapper, _clf_loaded
    if _clf_loaded:
        return _clf_wrapper is not None
    _clf_loaded = True
    try:
        from classifier.msclass_wrapper import get_msclass_wrapper
        _clf_wrapper = get_msclass_wrapper()
        if _clf_wrapper is None:
            raise RuntimeError('MSClassify adapter is unavailable')
        logger.info('MS2Classifier lazy-loaded: %d classes', len(_clf_wrapper.class_names))
        return True
    except Exception as _e:
        _clf_wrapper = None
        logger.warning('MS2Classifier not available: %s', _e)
        return False

# 全局缓存 + 线程锁
import threading
_lib_records = None
_mz_index = None
_nm_index = None       # 准确度修复 #1: 中性质量桶索引 (跨加合物匹配)
_lib_db_path = None    # 当前已加载库的路径 (检测 small↔full 切换并强制重载)
_lib_source_families = ()
_lib_lock = threading.Lock()


def _ik14(value):
    """Return InChIKey connectivity block when usable."""
    key = normalize_inchikey(value)
    if key:
        return key.split('-', 1)[0]
    # Compatibility for synthetic tests/legacy in-memory candidates.  Database
    # loading already clears malformed full keys, so this fallback cannot make
    # an invalid persisted key participate in production identity grouping.
    raw = str(value or '').strip().upper()
    head = raw.split('-', 1)[0]
    return head if re.fullmatch(r'[A-Z]{14}', head) else ''


def _rerank_identity(result):
    """Exact runtime identity with compatibility for legacy in-memory tests.

    Persisted malformed keys are cleared by the database loader, so the
    connectivity-block fallback cannot merge malformed database records.
    """
    if normalize_inchikey(result.get('inchikey')):
        return compound_identity_key(result)
    legacy = _ik14(result.get('inchikey'))
    return f'legacy_ik14::{legacy}' if legacy else compound_identity_key(result)


def _pathway_polarity_status(pathway, query_polarity):
    """Return available/unknown_polarity/mismatch for a mined DB pathway."""
    if not isinstance(pathway, dict):
        return 'missing', []
    polarities = [str(p).strip().lower() for p in pathway.get('polarities', []) if p]
    normalized = {
        ('neg' if p.startswith('neg') else 'pos' if p.startswith('pos') else p)
        for p in polarities
    }
    wanted = str(query_polarity or '').strip().lower()
    wanted = 'neg' if wanted.startswith('neg') else 'pos' if wanted.startswith('pos') else ''
    if not normalized:
        return 'unknown_polarity', []
    if wanted and wanted not in normalized:
        return 'mismatch', sorted(normalized)
    return 'available', sorted(normalized)


def _pathway_adduct_status(pathway, query_adduct, query_adduct_reliable):
    """Return the applicability of a legacy, compound-level pathway by adduct.

    Schema-v1 pathways aggregate all library spectra for a compound before
    constructing absolute product-ion m/z values.  If several chemically
    distinct adducts were aggregated, those absolute coordinates cannot be
    assigned to one query adduct and must not influence scoring.  A known
    single-adduct mismatch is likewise outside the pathway's applicability
    domain.  Missing query metadata remains ``unknown_query_adduct`` so older
    files preserve their previous behavior while the uncertainty is exposed.
    """
    if not isinstance(pathway, dict):
        return 'missing', []
    labels = [
        str(value).strip()
        for value in (pathway.get('adducts') or [])
        if str(value or '').strip()
    ]
    identities = {
        _ADDUCT_IDENTITY.get(re.sub(r'\s+', '', label).upper())
        for label in labels
        if re.sub(r'\s+', '', label).upper() not in _POLARITY_ONLY_ADDUCTS
    }
    identities.discard(None)
    if not identities:
        return 'unknown_pathway_adduct', labels
    if len(identities) > 1:
        return 'mixed_pathway_adducts', labels
    if not query_adduct_reliable:
        return 'unknown_query_adduct', labels
    query_identity = _ADDUCT_IDENTITY.get(
        re.sub(r'\s+', '', str(query_adduct or '')).upper()
    )
    if query_identity is None:
        return 'unknown_query_adduct', labels
    if query_identity not in identities:
        return 'mismatch', labels
    return 'available', labels


def _peak_present(mz, peaks, tol):
    return any(abs(float(p[0]) - float(mz)) <= tol for p in (peaks or []))


def _nearest_peak(target_mz, peaks, tol):
    """Return the nearest observed m/z and signed error within tolerance."""
    nearest = None
    for peak in peaks or []:
        try:
            observed = float(peak[0])
        except (TypeError, ValueError, IndexError):
            continue
        error = observed - float(target_mz)
        if abs(error) <= tol and (
            nearest is None or abs(error) < abs(nearest[1])
        ):
            try:
                intensity = abs(float(peak[1]))
            except (TypeError, ValueError, IndexError):
                intensity = None
            nearest = (observed, error, intensity)
    return nearest


def _relative_peaks(peaks, min_rel=0.05):
    if not peaks:
        return []
    try:
        max_i = max(abs(float(p[1])) for p in peaks) or 1.0
    except Exception:
        return []
    out = []
    for mz, inten in peaks:
        try:
            rel = abs(float(inten)) / max_i
            if rel >= min_rel:
                out.append((float(mz), rel))
        except Exception:
            continue
    return out


def _specific_adduct_identity(value):
    compact = re.sub(r'\s+', '', str(value or '')).upper()
    if not compact or compact in _POLARITY_ONLY_ADDUCTS:
        return None
    return _ADDUCT_IDENTITY.get(compact)


def _collision_energy_value(value, unit=''):
    if unit:
        parsed = collision_energy_value_unit({"value": value, "unit": unit})
        if parsed[0] is not None:
            return parsed
    return collision_energy_value_unit(value)


def _instrument_comparability(left_type, left_model, right_type, right_model):
    left_family = normalize_instrument_family(left_type, left_model)
    right_family = normalize_instrument_family(right_type, right_model)
    if left_family and right_family:
        return (
            'same_instrument_family'
            if left_family == right_family
            else 'different_instrument_families'
        )
    left_label = re.sub(
        r'[^a-z0-9]+', '', str(left_type or left_model or '').lower()
    )
    right_label = re.sub(
        r'[^a-z0-9]+', '', str(right_type or right_model or '').lower()
    )
    if left_label and right_label and left_label == right_label:
        return 'same_instrument_label'
    return 'instrument_metadata_incomplete'


def _candidate_pair_comparability(left, right):
    """Return whether two library spectra can support an isomer contrast."""
    left_adduct = _specific_adduct_identity(left.get('adduct'))
    right_adduct = _specific_adduct_identity(right.get('adduct'))
    if left_adduct and right_adduct:
        adduct_status = (
            'same_adduct' if left_adduct == right_adduct else 'adduct_mismatch'
        )
    else:
        adduct_status = 'adduct_metadata_incomplete'

    left_ce, left_unit = _collision_energy_value(
        left.get('lib_collision_energy'),
        left.get('lib_collision_energy_unit'),
    )
    right_ce, right_unit = _collision_energy_value(
        right.get('lib_collision_energy'),
        right.get('lib_collision_energy_unit'),
    )
    ce_delta = None
    if left_ce is not None and right_ce is not None and left_unit and right_unit:
        if left_unit != right_unit:
            ce_status = 'collision_energy_unit_mismatch'
        else:
            ce_delta = abs(left_ce - right_ce)
            ce_status = (
                'collision_energy_comparable'
                if ce_delta <= ISOMER_FRAGMENT_MAX_CE_DELTA
                else 'collision_energy_mismatch'
            )
    else:
        ce_status = 'collision_energy_metadata_incomplete'

    instrument_status = _instrument_comparability(
        left.get('lib_instrument_type'),
        left.get('lib_instrument'),
        right.get('lib_instrument_type'),
        right.get('lib_instrument'),
    )

    stage_check = source_ms2_pair_status(
        left.get('lib_source_ms_level'), right.get('lib_source_ms_level')
    )
    known_mismatch = (
        not stage_check['comparable_for_ordinary_ms2']
        or adduct_status == 'adduct_mismatch'
        or ce_status in {
            'collision_energy_unit_mismatch',
            'collision_energy_mismatch',
        }
        or instrument_status == 'different_instrument_families'
    )
    return {
        'comparable': not known_mismatch,
        'adduct_status': adduct_status,
        'collision_energy_status': ce_status,
        'collision_energy_delta': (
            round(ce_delta, 4) if ce_delta is not None else None
        ),
        'instrument_status': instrument_status,
        'source_ms_stage_comparability': stage_check,
    }


def _query_candidate_comparability(
    result,
    query_collision_energy,
    query_instrument_type=None,
    query_instrument=None,
):
    query_adduct = (
        _specific_adduct_identity(result.get('query_adduct'))
        if result.get('query_adduct_reliable') else None
    )
    library_adduct = _specific_adduct_identity(result.get('adduct'))
    if query_adduct and library_adduct:
        adduct_status = (
            'same_adduct' if query_adduct == library_adduct else 'adduct_mismatch'
        )
    else:
        adduct_status = 'adduct_metadata_incomplete'

    query_ce, query_unit = collision_energy_value_unit(query_collision_energy)
    library_ce, library_unit = _collision_energy_value(
        result.get('lib_collision_energy'),
        result.get('lib_collision_energy_unit'),
    )
    ce_delta = None
    if (
        query_ce is not None and library_ce is not None
        and query_unit and library_unit
    ):
        if query_unit != library_unit:
            ce_status = 'collision_energy_unit_mismatch'
        else:
            ce_delta = abs(query_ce - library_ce)
            ce_status = (
                'collision_energy_comparable'
                if ce_delta <= ISOMER_FRAGMENT_MAX_CE_DELTA
                else 'collision_energy_mismatch'
            )
    else:
        ce_status = 'collision_energy_metadata_incomplete'
    return {
        'adduct_status': adduct_status,
        'collision_energy_status': ce_status,
        'collision_energy_delta': (
            round(ce_delta, 4) if ce_delta is not None else None
        ),
        'instrument_status': _instrument_comparability(
            query_instrument_type,
            query_instrument,
            result.get('lib_instrument_type'),
            result.get('lib_instrument'),
        ),
        'source_ms_stage': source_ms_stage(result.get('lib_source_ms_level')),
    }


def _condition_metadata_complete(check):
    return bool(
        check.get('comparable')
        and (check.get('source_ms_stage_comparability') or {}).get('complete')
        and check.get('adduct_status') == 'same_adduct'
        and check.get('collision_energy_status')
        == 'collision_energy_comparable'
        and check.get('instrument_status')
        in {'same_instrument_family', 'same_instrument_label'}
    )


def _query_metadata_complete(check):
    return bool(
        (check.get('source_ms_stage') or {}).get('explicit_ms2_only')
        and check.get('adduct_status') == 'same_adduct'
        and check.get('collision_energy_status')
        == 'collision_energy_comparable'
        and check.get('instrument_status')
        in {'same_instrument_family', 'same_instrument_label'}
    )


def _reference_peak_profile(rows, *, neutral_losses=False):
    """Build a replicate-aware peak/loss profile for one exact identity."""
    rows = list(rows or [])
    reference_count = len(rows)
    if not reference_count:
        return []
    observations = []
    for spectrum_index, row in enumerate(rows):
        try:
            precursor_mz = float(row.get('lib_precursor_mz'))
        except (TypeError, ValueError):
            precursor_mz = None
        for product_mz, rel in _relative_peaks(
            row.get('lib_peaks') or [],
            min_rel=ISOMER_FRAGMENT_MIN_REL_INTENSITY,
        ):
            coordinate = product_mz
            if neutral_losses:
                if precursor_mz is None:
                    continue
                coordinate = precursor_mz - product_mz
                if coordinate <= 0 or coordinate >= precursor_mz:
                    continue
            observations.append({
                'coordinate': float(coordinate),
                'relative_intensity': float(rel),
                'product_mz': float(product_mz),
                'spectrum_index': spectrum_index,
            })
    clusters = []
    for observation in sorted(
        observations, key=lambda item: item['coordinate']
    ):
        nearest = None
        for cluster in clusters:
            error = abs(
                observation['coordinate'] - cluster['center']
            )
            if error <= ISOMER_DIAGNOSTIC_MZ_TOL_DA and (
                nearest is None or error < nearest[0]
            ):
                nearest = (error, cluster)
        if nearest is None:
            cluster = {
                'center': observation['coordinate'],
                'by_spectrum': {
                    observation['spectrum_index']: observation,
                },
            }
            clusters.append(cluster)
            continue
        cluster = nearest[1]
        existing = cluster['by_spectrum'].get(
            observation['spectrum_index']
        )
        if (
            existing is None
            or observation['relative_intensity']
            > existing['relative_intensity']
        ):
            cluster['by_spectrum'][
                observation['spectrum_index']
            ] = observation
        cluster['center'] = median(
            item['coordinate']
            for item in cluster['by_spectrum'].values()
        )

    required_support = (
        1 if reference_count == 1 else max(
            ISOMER_DIAGNOSTIC_MIN_REFERENCE_REPLICATES,
            int(math.ceil(
                reference_count
                * ISOMER_DIAGNOSTIC_MIN_PROFILE_PREVALENCE
            )),
        )
    )
    profile = []
    for cluster in clusters:
        values = list(cluster['by_spectrum'].values())
        support_count = len(values)
        if support_count < required_support:
            continue
        profile.append({
            'coordinate': float(median(
                item['coordinate'] for item in values
            )),
            'product_mz': float(median(
                item['product_mz'] for item in values
            )),
            'relative_intensity': float(median(
                item['relative_intensity'] for item in values
            )),
            'reference_support_count': support_count,
            'reference_spectrum_count': reference_count,
            'reference_prevalence': support_count / reference_count,
        })
    return sorted(
        profile,
        key=lambda item: (
            -item['reference_prevalence'],
            -item['relative_intensity'],
            item['coordinate'],
        ),
    )


def _candidate_specific_profile(target_profile, competitor_profiles):
    output = []
    for feature in target_profile:
        overlaps = []
        for identity, profile in competitor_profiles.items():
            nearest = min(
                (
                    item for item in profile
                    if abs(
                        item['coordinate'] - feature['coordinate']
                    ) <= ISOMER_DIAGNOSTIC_MZ_TOL_DA
                ),
                key=lambda item: abs(
                    item['coordinate'] - feature['coordinate']
                ),
                default=None,
            )
            if nearest is not None:
                overlaps.append({
                    'identity': identity,
                    'prevalence': nearest['reference_prevalence'],
                })
        if overlaps:
            continue
        output.append(dict(feature))
    return output


def _query_losses(query_clean, query_pmz):
    try:
        precursor = float(query_pmz)
    except (TypeError, ValueError):
        return []
    if precursor <= 0:
        return []
    losses = []
    for peak in query_clean or []:
        try:
            product_mz = float(peak[0])
            intensity = float(peak[1])
        except (TypeError, ValueError, IndexError):
            continue
        loss = precursor - product_mz
        if 0 < loss < precursor:
            losses.append([loss, intensity, product_mz])
    return losses


def _observed_profile(
    profile,
    query_clean,
    query_pmz,
    query_max_intensity,
    *,
    neutral_losses=False,
    candidate_identity=None,
    candidate_name=None,
):
    observed_rows = []
    query_view = (
        _query_losses(query_clean, query_pmz)
        if neutral_losses else query_clean
    )
    for feature in profile:
        nearest = _nearest_peak(
            feature['coordinate'],
            query_view,
            ISOMER_DIAGNOSTIC_MZ_TOL_DA,
        )
        item = {
            'library_relative_intensity': round(
                feature['relative_intensity'], 5
            ),
            'reference_support_count': int(
                feature['reference_support_count']
            ),
            'reference_spectrum_count': int(
                feature['reference_spectrum_count']
            ),
            'reference_prevalence': round(
                feature['reference_prevalence'], 5
            ),
            'observed': nearest is not None,
            'mass_error_da': (
                round(nearest[1], 6) if nearest is not None else None
            ),
            'query_relative_intensity': (
                round(nearest[2] / query_max_intensity, 5)
                if (
                    nearest is not None
                    and nearest[2] is not None
                    and query_max_intensity > 0
                ) else None
            ),
        }
        if candidate_identity:
            item['candidate_identity'] = candidate_identity
        if candidate_name:
            item['candidate_name'] = candidate_name
        if neutral_losses:
            item.update({
                'library_loss_da': round(feature['coordinate'], 5),
                'library_product_mz': round(feature['product_mz'], 5),
                'query_loss_da': (
                    round(nearest[0], 5)
                    if nearest is not None else None
                ),
                'query_product_mz': (
                    round(float(query_pmz) - nearest[0], 5)
                    if nearest is not None else None
                ),
            })
        else:
            item.update({
                'library_mz': round(feature['coordinate'], 5),
                'query_mz': (
                    round(nearest[0], 5)
                    if nearest is not None else None
                ),
            })
        observed_rows.append(item)
    return observed_rows


def _independent_observed_feature_count(
    product_ions,
    neutral_losses,
):
    coordinates = []
    for item in product_ions or []:
        if item.get('observed') and item.get('query_mz') is not None:
            coordinates.append(float(item['query_mz']))
    for item in neutral_losses or []:
        if (
            item.get('observed')
            and item.get('query_product_mz') is not None
        ):
            coordinates.append(float(item['query_product_mz']))
    independent = []
    for coordinate in sorted(coordinates):
        if not any(
            abs(coordinate - existing)
            <= ISOMER_DIAGNOSTIC_MZ_TOL_DA
            for existing in independent
        ):
            independent.append(coordinate)
    return len(independent)


def _apply_isomer_fragment_rerank(
    results,
    query_clean,
    query_pmz=None,
    query_collision_energy=None,
    query_instrument_type=None,
    query_instrument=None,
    apply_rank_override=None,
):
    """Measure candidate-specific product ions in true same-formula groups.

    Full-spectrum similarity is excellent for library search, but it can be dominated by
    common scaffold/adduct fragments.  When several same-formula structures survive
    the normal gates, the useful signal is often the product
    ions that are present in one candidate's library spectrum and absent from competitors.

    This evidence is not independently calibrated to identify a concrete isomer.  In
    production it is therefore annotation-only: scores and a potential bonus are exposed,
    while ``rank_score`` is left unchanged.  The legacy rank override is retained behind
    an explicit experimental flag for reproducible validation on authoritative structure
    ground truth.
    """
    if not results or len(results) < ISOMER_FRAGMENT_MIN_SHARED_CANDIDATES:
        return
    if apply_rank_override is None:
        apply_rank_override = USE_ISOMER_FRAGMENT_RANK_OVERRIDE
    apply_rank_override = bool(apply_rank_override)

    groups = {}
    for r in results:
        formula = (r.get('formula') or '').strip()
        identity = _rerank_identity(r)
        if formula and identity:
            groups.setdefault(formula, []).append(r)

    for group in groups.values():
        identity_count = len({_rerank_identity(r) for r in group})
        if identity_count < ISOMER_FRAGMENT_MIN_SHARED_CANDIDATES:
            continue
        try:
            query_max_intensity = max(
                abs(float(peak[1])) for peak in (query_clean or [])
            )
        except (TypeError, ValueError, IndexError):
            query_max_intensity = 0.0
        for r in group:
            rik = _rerank_identity(r)
            same_identity_rows = [
                candidate for candidate in group
                if (
                    _rerank_identity(candidate) == rik
                    and (
                        candidate is r
                        or _candidate_pair_comparability(
                            r, candidate
                        )['comparable']
                    )
                )
            ]
            other_identities = [
                other for other in group if _rerank_identity(other) != rik
            ]
            pair_checks = [
                (other, _candidate_pair_comparability(r, other))
                for other in other_identities
            ]
            comparable_others = [
                other for other, check in pair_checks if check['comparable']
            ]
            comparison_r = candidate_comparison_view(r)
            try:
                selected_cosine = float(
                    comparison_r.get('cosine', 0.0) or 0.0
                )
            except (TypeError, ValueError):
                selected_cosine = 0.0
            try:
                selected_rank = float(
                    comparison_r.get(
                        'rank_score',
                        comparison_r.get('final_score', 0.0),
                    ) or 0.0
                )
            except (TypeError, ValueError):
                selected_rank = 0.0
            competitive_others = [
                other for other in other_identities
                if is_spectrally_competitive_isomer(
                    selected_cosine,
                    selected_rank,
                    other,
                )
            ]
            competitive_identity_set = {
                _rerank_identity(other)
                for other in competitive_others
            }
            competitive_comparable_others = [
                other for other in comparable_others
                if _rerank_identity(other) in competitive_identity_set
            ]
            competitive_comparable_identity_set = {
                _rerank_identity(other)
                for other in competitive_comparable_others
            }
            all_competitive_candidates_comparable = bool(
                competitive_identity_set
                and competitive_comparable_identity_set
                == competitive_identity_set
            )
            query_check = _query_candidate_comparability(
                r,
                query_collision_energy,
                query_instrument_type=query_instrument_type,
                query_instrument=query_instrument,
            )
            target_product_profile = _reference_peak_profile(
                same_identity_rows,
            )
            target_loss_profile = _reference_peak_profile(
                same_identity_rows,
                neutral_losses=True,
            )
            competitor_rows_by_identity = {}
            for other in comparable_others:
                competitor_rows_by_identity.setdefault(
                    _rerank_identity(other), []
                ).append(other)
            competitor_product_profiles = {
                identity: _reference_peak_profile(rows)
                for identity, rows in competitor_rows_by_identity.items()
            }
            competitor_loss_profiles = {
                identity: _reference_peak_profile(
                    rows, neutral_losses=True
                )
                for identity, rows in competitor_rows_by_identity.items()
            }
            unique_product_profile = _candidate_specific_profile(
                target_product_profile,
                competitor_product_profiles,
            ) if comparable_others else []
            unique_loss_profile = _candidate_specific_profile(
                target_loss_profile,
                competitor_loss_profiles,
            ) if comparable_others else []
            product_ions = _observed_profile(
                unique_product_profile,
                query_clean,
                query_pmz,
                query_max_intensity,
            )
            neutral_losses = _observed_profile(
                unique_loss_profile,
                query_clean,
                query_pmz,
                query_max_intensity,
                neutral_losses=True,
            )

            opposing_product_ions = []
            opposing_neutral_losses = []
            competitor_names = {
                _rerank_identity(other): str(
                    other.get('compound_name') or ''
                )
                for other in comparable_others
            }
            for identity, profile in competitor_product_profiles.items():
                competitor_specific = _candidate_specific_profile(
                    profile,
                    {rik: target_product_profile},
                )
                opposing_product_ions.extend(_observed_profile(
                    competitor_specific,
                    query_clean,
                    query_pmz,
                    query_max_intensity,
                    candidate_identity=identity,
                    candidate_name=competitor_names.get(identity),
                ))
            for identity, profile in competitor_loss_profiles.items():
                competitor_specific = _candidate_specific_profile(
                    profile,
                    {rik: target_loss_profile},
                )
                opposing_neutral_losses.extend(_observed_profile(
                    competitor_specific,
                    query_clean,
                    query_pmz,
                    query_max_intensity,
                    neutral_losses=True,
                    candidate_identity=identity,
                    candidate_name=competitor_names.get(identity),
                ))

            competitive_opposing_product_ions = [
                item for item in opposing_product_ions
                if item.get('candidate_identity')
                in competitive_identity_set
            ]
            competitive_opposing_neutral_losses = [
                item for item in opposing_neutral_losses
                if item.get('candidate_identity')
                in competitive_identity_set
            ]
            supported_product = sum(
                bool(item['observed']) for item in product_ions
            )
            supported_losses = sum(
                bool(item['observed']) for item in neutral_losses
            )
            independent_supported = _independent_observed_feature_count(
                product_ions, neutral_losses
            )
            all_same_formula_opposing_supported = (
                _independent_observed_feature_count(
                    opposing_product_ions,
                    opposing_neutral_losses,
                )
            )
            opposing_supported = _independent_observed_feature_count(
                competitive_opposing_product_ions,
                competitive_opposing_neutral_losses,
            )
            if not competitive_others:
                diagnostic_status = (
                    'no_spectrally_competitive_isomer'
                )
            elif not all_competitive_candidates_comparable:
                diagnostic_status = 'candidate_metadata_not_comparable'
            elif not product_ions and not neutral_losses:
                diagnostic_status = 'no_candidate_specific_features'
            elif independent_supported >= 2:
                diagnostic_status = 'multiple_candidate_specific_features_observed'
            elif independent_supported == 1:
                diagnostic_status = 'single_candidate_specific_feature_observed'
            else:
                diagnostic_status = 'candidate_specific_features_not_observed'

            support = sum(
                item['relative_intensity']
                for item in unique_product_profile
                if _peak_present(
                    item['coordinate'],
                    query_clean,
                    ISOMER_DIAGNOSTIC_MZ_TOL_DA,
                )
            )
            total = sum(
                item['relative_intensity']
                for item in unique_product_profile
            ) or 1.0
            score = max(0.0, min(1.0, support / total))
            if not competitive_others:
                score = 0.0
            all_pair_statuses = [check for _, check in pair_checks]
            competitive_rows_by_identity = {
                identity: rows
                for identity, rows
                in competitor_rows_by_identity.items()
                if identity in competitive_identity_set
            }
            competitor_profile_metadata_complete = bool(
                competitive_rows_by_identity
                and all(
                    _condition_metadata_complete(check)
                    for other, check in pair_checks
                    if (
                        _rerank_identity(other)
                        in competitive_identity_set
                        and check.get('comparable')
                    )
                )
            )
            same_identity_metadata_complete = all(
                _condition_metadata_complete(
                    _candidate_pair_comparability(r, candidate)
                )
                for candidate in same_identity_rows
                if candidate is not r
            )
            metadata_complete = bool(
                all_competitive_candidates_comparable
                and competitor_profile_metadata_complete
                and same_identity_metadata_complete
                and _query_metadata_complete(query_check)
            )
            replication_complete = bool(
                len(same_identity_rows)
                >= ISOMER_DIAGNOSTIC_MIN_REFERENCE_REPLICATES
                and competitive_rows_by_identity
                and all(
                    len(rows)
                    >= ISOMER_DIAGNOSTIC_MIN_REFERENCE_REPLICATES
                    for rows in competitive_rows_by_identity.values()
                )
            )
            relationship_summary = summarize_isomer_relationships(
                r,
                list({
                    _rerank_identity(candidate): candidate
                    for candidate in competitive_others
                }.values()),
            )
            variant_subtype_counts = relationship_summary[
                'same_connectivity_variant_subtype_counts'
            ]
            if (
                relationship_summary['ambiguity_type']
                == 'mixed_isomer_group'
            ):
                resolution_status = (
                    'mixed_isomer_group_requires_multi_axis_resolution'
                )
            elif variant_subtype_counts.get(
                'structure_definition_specificity_variant'
            ):
                resolution_status = (
                    'same_connectivity_structure_definition_specificity_ambiguity'
                )
            elif variant_subtype_counts.get(
                'protonation_or_charge_layer_variant'
            ):
                resolution_status = (
                    'same_connectivity_protonation_variant_requires_standardization'
                )
            elif (
                relationship_summary['ambiguity_type']
                in {
                    'same_connectivity_variant_group',
                    'mixed_isomer_group',
                }
            ):
                resolution_status = (
                    'same_connectivity_variant_not_resolved_by_msms'
                )
            elif not competitive_others:
                resolution_status = (
                    'no_spectrally_competitive_isomer'
                )
            elif not all_competitive_candidates_comparable:
                resolution_status = 'candidate_metadata_not_comparable'
            elif independent_supported and opposing_supported:
                if not replication_complete:
                    resolution_status = (
                        'conflicting_features_observed_but_unreplicated'
                    )
                elif not metadata_complete:
                    resolution_status = (
                        'conflicting_features_observed_but_metadata_incomplete'
                    )
                else:
                    resolution_status = 'conflicting_diagnostic_evidence'
            elif opposing_supported and not independent_supported:
                if not replication_complete:
                    resolution_status = (
                        'competitor_features_observed_but_unreplicated'
                    )
                elif not metadata_complete:
                    resolution_status = (
                        'competitor_features_observed_but_metadata_incomplete'
                    )
                else:
                    resolution_status = (
                        'competitor_supported_pending_validation'
                    )
            elif independent_supported >= 2:
                if not replication_complete:
                    resolution_status = (
                        'selected_candidate_supported_but_unreplicated'
                    )
                elif not metadata_complete:
                    resolution_status = (
                        'selected_candidate_supported_but_metadata_incomplete'
                    )
                else:
                    resolution_status = (
                        'selected_candidate_supported_pending_validation'
                    )
            elif independent_supported == 1:
                resolution_status = (
                    'single_selected_feature_insufficient'
                )
            else:
                resolution_status = 'diagnostically_uninformative'
            detail = {
                'schema_version': 3,
                'diagnostic_status': diagnostic_status,
                'unique_supported': supported_product,
                'unique_total': len(unique_product_profile),
                'neutral_loss_supported': supported_losses,
                'neutral_loss_total': len(unique_loss_profile),
                'independent_supported_features': independent_supported,
                'opposing_supported_features': opposing_supported,
                'all_same_formula_opposing_supported_features': (
                    all_same_formula_opposing_supported
                ),
                'evidence_balance': (
                    independent_supported - opposing_supported
                ),
                'candidate_specific_product_ions': product_ions[:12],
                'candidate_specific_neutral_losses': neutral_losses[:12],
                'opposing_candidate_specific_product_ions': (
                    competitive_opposing_product_ions[:12]
                ),
                'opposing_candidate_specific_neutral_losses': (
                    competitive_opposing_neutral_losses[:12]
                ),
                'same_formula_candidate_count': identity_count,
                'comparable_candidate_count': len({
                    _rerank_identity(other) for other in comparable_others
                }),
                'spectrally_competitive_candidate_count': len(
                    competitive_identity_set
                ),
                'spectrally_competitive_comparable_candidate_count': len({
                    _rerank_identity(other)
                    for other in competitive_comparable_others
                }),
                'candidate_pair_comparability': all_pair_statuses[:12],
                'query_candidate_comparability': query_check,
                'metadata_comparability': (
                    'complete' if metadata_complete else 'incomplete_or_mismatched'
                ),
                'reference_replication': {
                    'status': (
                        'complete'
                        if replication_complete else 'insufficient'
                    ),
                    'selected_reference_spectra': len(
                        same_identity_rows
                    ),
                    'competitor_reference_spectra': {
                        identity: len(rows)
                        for identity, rows
                        in competitive_rows_by_identity.items()
                    },
                    'minimum_required': (
                        ISOMER_DIAGNOSTIC_MIN_REFERENCE_REPLICATES
                    ),
                },
                'ambiguity_type': relationship_summary[
                    'ambiguity_type'
                ],
                'relationship_counts': relationship_summary[
                    'relationship_counts'
                ],
                'same_connectivity_variant_subtype_counts': (
                    variant_subtype_counts
                ),
                'msms_exact_resolution_capable': relationship_summary[
                    'msms_exact_resolution_capable'
                ],
                'resolution_status': resolution_status,
                'diagnostic_mz_tolerance_da': (
                    ISOMER_DIAGNOSTIC_MZ_TOL_DA
                ),
                'weight': ISOMER_FRAGMENT_RERANK_WEIGHT,
                'rank_override_enabled': apply_rank_override,
                'runtime_influence': apply_rank_override,
                'eligible_for_exact_identity_promotion': False,
                'identity_promotion_blocker': (
                    'diagnostic_evidence_not_compound_split_calibrated'
                ),
                'query_evidence_view': 'broad_measured_ms2_shadow',
                'query_evidence_peak_count': len(query_clean or []),
                'query_evidence_min_relative_intensity': (
                    ISOMER_DIAGNOSTIC_QUERY_MIN_REL_INTENSITY
                ),
            }
            bonus = ISOMER_FRAGMENT_RERANK_WEIGHT * score
            r['isomer_fragment_score'] = round(score, 4)
            r['isomer_fragment_potential_bonus'] = round(bonus, 4)
            r['isomer_fragment_bonus'] = (
                round(bonus, 4) if apply_rank_override else 0.0
            )
            r['isomer_fragment_detail'] = detail
            r['isomer_fragment_rank_applied'] = apply_rank_override
            if apply_rank_override:
                r['rank_score'] = min(
                    1.0, r.get('final_score', 0.0) + bonus
                )


def _apply_formula_isotope_rerank(
    results,
    query_isotope_evidence,
    max_adjustment=0.01,
):
    """Use an adjacent-MS1 isotope envelope as a bounded ranking tie-breaker.

    The adjustment applies to ``rank_score`` only.  It deliberately does not
    alter ``final_score`` or confidence calibration because isotope ratios can
    be distorted by co-elution and isolation interference.
    """
    if not results or not query_isotope_evidence:
        return
    if query_isotope_evidence.get('source') != 'adjacent_ms1_scan':
        return
    if int(query_isotope_evidence.get('n_isotopes_observed') or 0) < 1:
        return
    try:
        from utils.isotope import evaluate_formula_isotope_consistency
    except Exception:
        return

    bound = max(0.0, min(float(max_adjustment), 0.02))
    for result in results:
        try:
            evidence = evaluate_formula_isotope_consistency(
                result.get('formula', ''), query_isotope_evidence
            )
        except Exception:
            evidence = None
        result['formula_isotope_evidence'] = evidence
        if not evidence:
            continue
        score = max(0.0, min(1.0, float(evidence.get('score', 0.0))))
        adjustment = bound * (2.0 * score - 1.0)
        base_rank = float(result.get('rank_score', result.get('final_score', 0.0)) or 0.0)
        result['isotope_rank_adjustment'] = round(adjustment, 5)
        result['rank_score'] = max(0.0, min(1.0, base_rank + adjustment))


def _ranking_key(result):
    ppm = result.get('mass_ppm')
    try:
        ppm_abs = abs(float(ppm)) if ppm is not None else 999999.0
    except (TypeError, ValueError):
        ppm_abs = 999999.0
    comparison = candidate_comparison_view(result)
    return (
        float(comparison.get(
            'rank_score', comparison.get('final_score', 0.0)
        ) or 0.0),
        (
            float(result.get('isomer_fragment_score', 0.0) or 0.0)
            if USE_ISOMER_FRAGMENT_RANK_OVERRIDE else 0.0
        ),
        float(comparison.get('cosine', 0.0) or 0.0),
        float(comparison.get('reverse_match', 0.0) or 0.0),
        -ppm_abs,
    )


_SOURCE_TIER_PRIORITY = {
    'curated': 4,
    'experimental': 3,
    'massbank_supplement': 3,
    'merged': 2,
    'legacy_unresolved': 1,
}


def _stable_identity_suffix(result):
    """Final lexical tie-breakers for reproducible ranking across DB order."""
    return (
        compound_identity_key(result),
        str(result.get('dataset_source') or ''),
        str(result.get('lib_spectrum_id') or result.get('spectrum_id') or ''),
        str(result.get('compound_name') or ''),
    )


def _stable_ranking_sort_key(result):
    """Ascending total order whose first item is the production Top-1.

    Score semantics remain unchanged.  Structure/source preferences are used
    only after all scientific score and mass-error fields tie; exact identity
    plus accession then makes the result independent of SQLite row order.
    """
    ppm = result.get('mass_ppm')
    try:
        ppm_abs = abs(float(ppm)) if ppm is not None else 999999.0
    except (TypeError, ValueError):
        ppm_abs = 999999.0
    source_priority = _SOURCE_TIER_PRIORITY.get(
        str(result.get('source_tier') or '').strip().lower(), 0
    )
    has_structure = bool(result.get('has_structure') and result.get('smiles'))
    comparison = candidate_comparison_view(result)
    return (
        -float(comparison.get(
            'rank_score', comparison.get('final_score', 0.0)
        ) or 0.0),
        -(
            float(result.get('isomer_fragment_score', 0.0) or 0.0)
            if USE_ISOMER_FRAGMENT_RANK_OVERRIDE else 0.0
        ),
        -float(comparison.get('cosine', 0.0) or 0.0),
        -float(comparison.get('reverse_match', 0.0) or 0.0),
        ppm_abs,
        -int(has_structure),
        -source_priority,
        *_stable_identity_suffix(result),
    )


def _top_candidates_with_best(best, candidates, limit=TOP_K_CANDIDATES):
    """Return raw-spectrum Top-K while guaranteeing selected identity is first."""
    ordered = sorted(candidates, key=_stable_ranking_sort_key)
    if not best:
        return [
            candidate_comparison_view(candidate)
            for candidate in ordered[:limit]
        ]
    best_key = compound_identity_key(best)
    selected = next(
        (candidate for candidate in ordered
         if compound_identity_key(candidate) == best_key),
        None,
    )
    if selected is None:
        # A legacy/cache caller can provide a selected result that is absent
        # from the candidate map.  The public contract is still stronger than
        # the cache contents: include the selected raw-spectrum view first.
        selected_view = candidate_comparison_view(best)
        return [
            selected_view,
            *[
                candidate_comparison_view(candidate)
                for candidate in ordered
                if compound_identity_key(candidate) != best_key
            ],
        ][:limit]
    return [
        candidate_comparison_view(selected),
        *[
            candidate_comparison_view(candidate)
            for candidate in ordered
            if compound_identity_key(candidate) != best_key
        ],
    ][:limit]


def _isomer_tie_key(result):
    ppm = result.get('mass_ppm')
    try:
        ppm_abs = abs(float(ppm)) if ppm is not None else 999999.0
    except (TypeError, ValueError):
        ppm_abs = 999999.0
    comparison = candidate_comparison_view(result)
    return (
        (
            float(result.get('isomer_fragment_score', 0.0) or 0.0)
            if USE_ISOMER_FRAGMENT_RANK_OVERRIDE else 0.0
        ),
        float(comparison.get('cosine', 0.0) or 0.0),
        float(comparison.get('reverse_match', 0.0) or 0.0),
        -ppm_abs,
        float(comparison.get(
            'rank_score', comparison.get('final_score', 0.0)
        ) or 0.0),
    )


def _stable_isomer_sort_key(result):
    ppm = result.get('mass_ppm')
    try:
        ppm_abs = abs(float(ppm)) if ppm is not None else 999999.0
    except (TypeError, ValueError):
        ppm_abs = 999999.0
    comparison = candidate_comparison_view(result)
    return (
        -(
            float(result.get('isomer_fragment_score', 0.0) or 0.0)
            if USE_ISOMER_FRAGMENT_RANK_OVERRIDE else 0.0
        ),
        -float(comparison.get('cosine', 0.0) or 0.0),
        -float(comparison.get('reverse_match', 0.0) or 0.0),
        ppm_abs,
        -float(comparison.get(
            'rank_score', comparison.get('final_score', 0.0)
        ) or 0.0),
        *_stable_identity_suffix(result),
    )


def _select_best_result(results):
    """Choose best candidate using calibrated scores and deterministic tie-breaks.

    Candidate-specific isomer fragments may lead the near-tie sort only when the
    explicit experimental rank-override flag is enabled.  Production defaults to
    cosine/reverse/mass tie-breaks and reports the isomer evidence separately.
    """
    if not results:
        return None
    best_by_final = min(results, key=_stable_ranking_sort_key)
    best_comparison = candidate_comparison_view(best_by_final)
    best_final = float(best_comparison.get(
        'rank_score', best_comparison.get('final_score', 0.0)
    ) or 0.0)
    selected = best_by_final
    formula = (best_by_final.get('formula') or '').strip()
    if formula:
        tied = [
            r for r in results
            if (r.get('formula') or '').strip() == formula
            and abs(float(
                candidate_comparison_view(r).get(
                    'rank_score',
                    candidate_comparison_view(r).get('final_score', 0.0),
                ) or 0.0
            ) - best_final) <= ISOMER_TIE_FINAL_EPS
        ]
        tied_identities = {
            compound_identity_key(r) for r in tied if compound_identity_key(r)
        }
        if len(tied_identities) >= 2:
            selected = min(tied, key=_stable_isomer_sort_key)

    # A high-quality name-only library spectrum can sit a fraction above an otherwise
    # equivalent curated structure (e.g. GNPS ``naringenin [M-H]-`` 0.9234 versus
    # structured Naringenin 0.9227).  Treating the former as the headline downgrades a
    # known compound to "Spectral Lead / no structure" even though the same Top-K
    # contains its exact-name InChIKey/SMILES record.  Promote only exact normalized
    # names and only inside a tight score/ppm neighborhood; never propagate structure
    # across fuzzy names or generic annotations.
    if selected.get('name_only_lead'):
        raw_name = str(selected.get('compound_name') or '').lower().split('(')[0].split('[')[0]
        name_key = re.sub(r'[^a-z0-9]', '', raw_name)
        selected_rank = float(selected.get('rank_score', selected.get('final_score', 0.0)) or 0.0)
        selected_cosine = float(selected.get('cosine', 0.0) or 0.0)
        selected_ppm = selected.get('mass_ppm')
        structured = []
        if len(name_key) >= 5:
            for candidate in results:
                if not candidate.get('has_structure') or not candidate.get('smiles'):
                    continue
                candidate_name = str(candidate.get('compound_name') or '').lower().split('(')[0].split('[')[0]
                if re.sub(r'[^a-z0-9]', '', candidate_name) != name_key:
                    continue
                rank = float(candidate.get('rank_score', candidate.get('final_score', 0.0)) or 0.0)
                cosine = float(candidate.get('cosine', 0.0) or 0.0)
                candidate_ppm = candidate.get('mass_ppm')
                ppm_close = (
                    selected_ppm is None or candidate_ppm is None
                    or abs(float(candidate_ppm) - float(selected_ppm)) <= 10.0
                )
                if rank >= selected_rank - 0.01 and cosine >= selected_cosine - 0.02 and ppm_close:
                    structured.append(candidate)
        if structured:
            promoted = min(structured, key=_stable_ranking_sort_key)
            promoted['structure_promoted_over_name_only'] = True
            promoted['supporting_name_only_lead'] = {
                'compound_name': selected.get('compound_name'),
                'source_compound_name': selected.get('source_compound_name', ''),
                'dataset_source': selected.get('dataset_source', ''),
                'lib_spectrum_id': selected.get('lib_spectrum_id'),
                'final_score': selected.get('final_score'),
                'cosine': selected.get('cosine'),
                'mass_ppm': selected.get('mass_ppm'),
            }
            selected = promoted
    return selected


def _passes_strict_evidence_gate(final, entropy, reverse, ppm, enabled=None,
                                 max_mass_ppm=None):
    if enabled is None:
        enabled = USE_STRICT_EVIDENCE_GATE
    if not enabled:
        return True
    if final < STRICT_MIN_FINAL:
        return False
    if entropy < STRICT_MIN_ENTROPY:
        return False
    if reverse < STRICT_MIN_REVERSE_MATCH:
        return False
    if max_mass_ppm is None:
        max_mass_ppm = STRICT_MAX_MASS_PPM
    if ppm is not None:
        try:
            if abs(float(ppm)) > float(max_mass_ppm):
                return False
        except (TypeError, ValueError):
            pass
    return True


def _drop_dominant_low_mz_interference(peaks, precursor_mz):
    """在严格守卫下移除一个疑似共碎裂低 m/z 基峰，返回新峰表或 None。

    这不是通用去噪：只有低于阈值、且强度至少为第二强峰指定倍数时才触发。
    原检索结果始终保留，补充结果还需经过 `_adaptive_result_better`。
    """
    usable = []
    for p in peaks or []:
        try:
            mz, inten = float(p[0]), abs(float(p[1]))
        except (TypeError, ValueError, IndexError):
            continue
        if inten > 0 and mz <= float(precursor_mz) + 1.5 and abs(mz - float(precursor_mz)) > 0.5:
            usable.append((mz, inten))
    if len(usable) < 2:
        return None
    strongest = sorted(usable, key=lambda x: x[1], reverse=True)
    (base_mz, base_i), (_, second_i) = strongest[:2]
    if base_mz > LOW_MZ_INTERFERENCE_MAX_MZ:
        return None
    if second_i <= 0 or base_i / second_i < LOW_MZ_INTERFERENCE_DOMINANCE:
        return None
    removed = False
    out = []
    for p in peaks or []:
        try:
            same = not removed and abs(float(p[0]) - base_mz) <= 1e-6 and abs(abs(float(p[1])) - base_i) <= 1e-6
        except (TypeError, ValueError, IndexError):
            same = False
        if same:
            removed = True
            continue
        out.append(p)
    return out if removed else None


def _clean_field(v, default=''):
    """清理库字段值: None / float('nan') / 'nan'/'None' / 空串 → default。
    关键陷阱: float('nan') 是 truthy, `nan or x` 返回 nan 而非 x, 必须显式检查
    (否则 compound_name=NaN 会显示成 'nan', smiles=NaN 不会被 `or ''` 兜底)。
    库(MoNA/Spectraverse)部分条目存在 NaN/空 字段。"""
    if v is None:
        return default
    try:
        if v != v:  # float('nan') 唯一性质: nan != nan → True
            return default
    except Exception:
        pass
    s = str(v).strip()
    if s.lower() in ('', 'nan', 'none', 'null', 'na', 'n/a', 'n.a.'):
        return default
    return s


def _charge_polarity(value):
    """Return POS/NEG for numeric and common textual charge encodings."""
    if value is None:
        return None
    text = str(value).strip()
    if not text:
        return None
    if text.startswith('-') or text.endswith('-'):
        return 'NEG'
    if text.startswith('+') or text.endswith('+'):
        return 'POS'
    try:
        number = float(text)
    except (TypeError, ValueError):
        return None
    if number < 0:
        return 'NEG'
    if number > 0:
        return 'POS'
    return None


def _library_record_neutral_mass(record):
    """Use a verified monoisotopic correction before the original fallback."""
    corrected = verified_neutral_mass_correction(record)
    if corrected is not None:
        return corrected, DERIVED_NEUTRAL_MASS_SOURCE
    mass = record.get('exact_mass')
    if mass is not None and mass > 0:
        return mass, 'compound_exact_mass'
    mass = (
        _calc_neutral_mass(record['precursor_mz'], record['adduct'])
        if _specific_adduct_polarity(record.get('adduct')) else None
    )
    return mass, 'library_adduct' if mass is not None else None


def iter_normalized_library_records(db_path):
    """Stream the exact filtered runtime records, closing SQLite on early exit.

    Used both by the legacy in-memory loader and the bounded-memory offline
    record-store builder. The SQL row order, peak cleaning and quality gates
    are shared; no candidate eligibility or scoring threshold is changed.
    """
    db_path = Path(db_path).resolve()
    conn = sqlite3.connect(f'file:{db_path}?mode=ro', uri=True)
    try:
        conn.row_factory = sqlite3.Row

        # dataset/source quality fields由迁移脚本新增; 旧库降级为 unknown。
        _has_ds = any(
            row[1] == 'dataset_source'
            for row in conn.execute('PRAGMA table_info(Spectra)').fetchall()
        )
        _ds_col = 's.dataset_source,' if _has_ds else "'' AS dataset_source,"
        _has_source_tier = any(
            row[1] == 'source_tier'
            for row in conn.execute('PRAGMA table_info(Spectra)').fetchall()
        )
        _source_tier_col = 's.source_tier,' if _has_source_tier else "'unknown' AS source_tier,"
        _has_intensity_origin = any(
            row[1] == 'intensity_origin'
            for row in conn.execute('PRAGMA table_info(Spectra)').fetchall()
        )
        _intensity_origin_col = (
            's.intensity_origin,' if _has_intensity_origin
            else "'unknown' AS intensity_origin,"
        )
        _has_quality_flags = any(
            row[1] == 'quality_flags_json'
            for row in conn.execute('PRAGMA table_info(Spectra)').fetchall()
        )
        _quality_flags_col = (
            's.quality_flags_json,' if _has_quality_flags
            else "'[]' AS quality_flags_json,"
        )
        _has_collision_energy = any(
            row[1] == 'collision_energy'
            for row in conn.execute('PRAGMA table_info(Spectra)').fetchall()
        )
        _collision_energy_col = (
            's.collision_energy,' if _has_collision_energy
            else "NULL AS collision_energy,"
        )
        _spectrum_columns = {
            row[1] for row in conn.execute('PRAGMA table_info(Spectra)').fetchall()
        }
        _ce_metadata_cols = " ".join(
            (
                f"s.{name} AS {name},"
                if name in _spectrum_columns
                else f"NULL AS {name},"
            )
            for name in (
                "collision_energy_raw",
                "collision_energy_unit",
                "collision_energy_mode",
                "collision_energy_min",
                "collision_energy_max",
                "instrument_type",
                "instrument",
                "ion_mode",
                "source_accession",
                "activation_method",
                "instrument_domain_id",
                "collision_cell_region",
                "acquisition_method",
            )
        )
        # v3 Phase 2: 优先用自动分类的 classification_v3 (治 "other" 67%→9%); 无此列/NULL → 降级原 classification
        _has_cv3 = any(
            row[1] == 'classification_v3'
            for row in conn.execute('PRAGMA table_info(Compounds)').fetchall()
        )
        _cv3_col = 'c.classification_v3,' if _has_cv3 else "NULL AS classification_v3,"
        _compound_columns = {
            row[1] for row in conn.execute('PRAGMA table_info(Compounds)').fetchall()
        }
        _classification_source_col = (
            'c.classification_source,'
            if 'classification_source' in _compound_columns
            else "NULL AS classification_source,"
        )
        _classification_confidence_col = (
            'c.classification_confidence,'
            if 'classification_confidence' in _compound_columns
            else "NULL AS classification_confidence,"
        )
        # P-A: fragmentation_pathway 列 (P2 scripts/mine_pathways.py 挖出; 旧库无此列 → NULL 降级)
        _has_fp = any(
            row[1] == 'fragmentation_pathway'
            for row in conn.execute('PRAGMA table_info(Compounds)').fetchall()
        )
        _fp_col = 'c.fragmentation_pathway' if _has_fp else "NULL AS fragmentation_pathway"
        # 流式遍历游标，避免 Full 库 460 万个 sqlite3.Row 与 records 同时常驻。
        # 旧 fetchall 会在 name-only 扩容后把加载峰值推到约 29 GB 并触发 macOS 压缩。
        rows = conn.execute('''
            SELECT s.precursor_mz, s.peaks_json, s.spectrum_id,
                   s.adduct, s.retention_time, s.charge, %s %s %s %s %s %s
                   c.compound_name, c.smiles, c.inchikey,
                   c.molecular_formula, c.exact_mass, %s %s %s
                   c.classification, %s
            FROM Spectra s
            JOIN Mapping m ON m.spectrum_id = s.id
            JOIN Compounds c ON m.compound_id = c.id
        ''' % (
            _collision_energy_col, _ce_metadata_cols, _ds_col, _source_tier_col,
            _intensity_origin_col, _quality_flags_col, _cv3_col, _classification_source_col,
            _classification_confidence_col, _fp_col,
        ))

        loaded_count = 0
        _drop_no_precursor = 0
        _drop_low_information = 0
        _drop_reverse_only = 0
        _drop_no_structure = 0
        _name_only_leads = 0
        for r in rows:
            _source_tier = _clean_field(r['source_tier'], 'unknown').lower()
            if not source_tier_is_searchable(_source_tier):
                continue
            _quality_flags_raw = _clean_field(r['quality_flags_json'], '[]')
            try:
                _quality_flags = json.loads(_quality_flags_raw)
                if not isinstance(_quality_flags, list):
                    _quality_flags = []
            except (TypeError, ValueError):
                _quality_flags = []
            pmz = _safe_float(r['precursor_mz'])
            if pmz is None or pmz <= 0:
                _drop_no_precursor += 1
                continue
            try:
                peaks = json.loads(r['peaks_json'])
            except Exception:
                peaks = []
            if peaks:
                peaks = sorted(peaks, key=lambda p: p[0])
            cleaned = clean_peaks(peaks, pmz, precursor_charge=r['charge']) if peaks else []
            if len(cleaned) < 3:
                # The import quality gate rejects fewer than three fragments.
                # Enforce the same invariant while loading legacy/direct rows
                # so one- or two-peak records cannot manufacture high cosine.
                _drop_low_information += 1
                continue
            if 'flat_fragment_intensities_reverse_only' in _quality_flags:
                _drop_reverse_only += 1
                continue
            _search_precursor_mz, _precursor_mz_source = (
                _audited_search_precursor_from_flags(_quality_flags)
            )
            # 无结构谱图只保留为严格门控的 name-only spectral lead；它绝不能进入
            # 结构专家、裂解路径或结构置信度输出。纯登记号/噪声名仍在加载阶段丢弃。
            _smi = _clean_field(r['smiles'])
            _has_structure = bool(_smi)
            _source_name = _clean_field(r['compound_name'])
            _raw_inchikey = _clean_field(r['inchikey'])
            _normalized_inchikey = normalize_inchikey(_raw_inchikey)
            _raw_compound_name = (
                _source_name
                or ((_normalized_inchikey or _raw_inchikey).split('-')[0] or '')
                or f"Spectrum_{r['spectrum_id']}"
            )
            _display_compound_name, _name_status = resolve_compound_name(
                _raw_compound_name, _normalized_inchikey
            )
            if not _has_structure:
                if not LOAD_NAME_ONLY_LEADS:
                    _drop_no_structure += 1
                    continue
                if not _source_name or not is_readable_compound_name(_source_name):
                    _drop_no_structure += 1
                    continue
                _name_only_leads += 1
            _record = {
                'precursor_mz': pmz,
                'reported_precursor_mz': pmz,
                'search_precursor_mz': _search_precursor_mz,
                'precursor_mz_source': _precursor_mz_source,
                'peaks': peaks,
                'clean_peaks': cleaned,
                'spectrum_id': r['spectrum_id'],
                'adduct': _clean_field(r['adduct']),
                'retention_time': _safe_float(r['retention_time']),
                # Explicit adduct is preferred; DB charge safely fills legacy
                # rows whose adduct is missing or only descriptive text.
                'charge': r['charge'],
                'collision_energy': _clean_field(r['collision_energy']),
                'collision_energy_raw': _clean_field(r['collision_energy_raw']),
                'collision_energy_unit': _clean_field(r['collision_energy_unit']),
                'collision_energy_mode': _clean_field(r['collision_energy_mode']),
                'collision_energy_min': _safe_float(r['collision_energy_min']),
                'collision_energy_max': _safe_float(r['collision_energy_max']),
                'instrument_type': _clean_field(r['instrument_type']),
                'instrument': _clean_field(r['instrument']),
                'source_ion_mode': _clean_field(r['ion_mode']),
                'source_accession': _clean_field(r['source_accession']),
                'activation_method': _clean_field(r['activation_method']),
                'instrument_domain_id': _clean_field(r['instrument_domain_id']),
                'collision_cell_region': _clean_field(r['collision_cell_region']),
                'acquisition_method': _clean_field(r['acquisition_method']),
                'polarity': classify_spectrum_polarity(
                    _clean_field(r['adduct']),
                    r['charge'],
                    _clean_field(r['ion_mode']),
                ),
                # 无名: compound_name→InChIKey前14→谱图ID 三级兜底。NCGC 原始标识
                # 只用于追溯；展示名由完整 InChIKey 精确别名或去前缀规则生成。
                'compound_name': _display_compound_name,
                'source_compound_name': (
                    _raw_compound_name
                    if _display_compound_name != _raw_compound_name else ''
                ),
                'name_resolution_status': _name_status,
                'smiles': _smi,
                'has_structure': _has_structure,
                'name_only_lead': not _has_structure,
                # Only a canonical full InChIKey participates in identity
                # aggregation.  Malformed source labels remain available for
                # audit but can no longer merge unrelated candidates.
                'inchikey': _normalized_inchikey,
                'reported_inchikey': (
                    _raw_inchikey if _raw_inchikey != _normalized_inchikey else ''
                ),
                'formula': _clean_field(r['molecular_formula']),
                'exact_mass': _safe_float(r['exact_mass']),
                'classification': (
                    _clean_field(r['classification_v3'], None)
                    or _clean_field(r['classification'], 'Unknown')
                ) if _has_structure else 'Unknown',
                'classification_source': (
                    _clean_field(r['classification_source'], 'unknown')
                    if _has_structure else 'none'
                ),
                'classification_confidence': (
                    _safe_float(r['classification_confidence'])
                    if _has_structure else None
                ),
                'dataset_source': _clean_field(r['dataset_source']),
                'source_tier': _source_tier,
                'intensity_origin': _clean_field(r['intensity_origin'], 'unknown'),
                'quality_flags': _quality_flags,
                'pathway_json': _clean_field(r['fragmentation_pathway']) if _has_structure else '',
            }
            annotate_compound_dict(_record)
            nm, nm_source = _library_record_neutral_mass(_record)
            if nm is not None:
                _record['neutral_mass'] = nm
                _record['neutral_mass_source'] = nm_source
            loaded_count += 1
            yield _record


        logger.info(
            '  Loaded %d spectra (%d name-only leads; dropped: %d no-precursor, %d low-information, %d reverse-only, %d unresolved/noisy no-structure)',
            loaded_count, _name_only_leads, _drop_no_precursor,
            _drop_low_information, _drop_reverse_only, _drop_no_structure,
        )
    finally:
        conn.close()


def load_library(db_path=None):
    """加载质谱库到内存，构建 m/z 索引（线程安全）。

    缓存键为 db_path：在 small↔full 之间切换时会自动重载，
    而不是返回旧库 (修复"切到 Full 仍返回 Small"的缓存 bug)。
    """
    global _lib_records, _mz_index, _nm_index, _lib_db_path
    global _lib_source_families
    if db_path is None:
        db_path = DB_PATH
    db_path = Path(db_path).resolve()

    # Fast path: 同一个库已加载就直接返回 (避免每次分析重复加载)
    if _lib_records is not None and _lib_db_path == db_path:
        return _lib_records, _mz_index

    # Slow path: acquire lock, double-check
    with _lib_lock:
        if _lib_records is not None and _lib_db_path == db_path:
            return _lib_records, _mz_index

        if _lib_records is not None and _lib_db_path != db_path:
            logger.info('Switching library: %s -> %s, reloading', _lib_db_path, db_path)
        else:
            logger.info('Loading library: %s', db_path)
        from matching.library_record_store import try_load_runtime_record_store
        cached = try_load_runtime_record_store(db_path)
        if cached is not None:
            records, mz_index, nm_index = cached
            _lib_records, _mz_index, _nm_index = records, mz_index, nm_index
            _lib_db_path = db_path
            _lib_source_families = tuple(records.source_families)
            logger.info('  Loaded %d spectra via verified read-only record mmap', len(records))
            return records, mz_index

        records = list(iter_normalized_library_records(db_path))

        bucket_size = MZ_BUCKET_DALTON
        mz_index = {}
        for idx, rec in enumerate(records):
            indexed_mz = _record_index_precursor_mz(rec)
            if indexed_mz is None:
                continue
            bucket = int(indexed_mz / bucket_size)
            mz_index.setdefault(bucket, []).append(idx)
        for bucket in mz_index:
            mz_index[bucket].sort(
                key=lambda i: _record_index_precursor_mz(records[i])
            )

        # 中性质量桶索引 (跨加合物匹配的核心)。
        # 准确度修复(检索端): 优先用分子结构(SMILES)算化合物真实中性质量——它加合物无关。
        # 旧实现用 precursor_mz+adduct标签 推, 但 ~6.7万条库谱图的 adduct 字段是极性串
        # ("Positive"/"Negative")或简写("M+H"), 会被当 [M+H]+/[M-H]- 假设 → 真实加合物若是
        # [M+Na]+/[M+K]+ 则中性质量建错桶(差22-38Da), 跨加合物命中静默失效。
        # 实例: Geniposidic acid 谱图 [M+Na]+@397.11 标 "Positive" → 旧法得 396.10, 真实 374.12。
        # exact_mass 是库导入阶段写入的化合物中性精确质量; Web 加载阶段不能再对
        # 133万条谱图逐条跑 RDKit MolFromSmiles, 否则 full 库会卡死并刷大量 SMILES
        # 解析错误。回退链: exact_mass → precursor+adduct。
        nm_index = {}
        _n_exact = 0
        _n_adduct = 0
        for idx, rec in enumerate(records):
            nm, nm_source = _library_record_neutral_mass(rec)
            if nm_source in ('compound_exact_mass', DERIVED_NEUTRAL_MASS_SOURCE):
                _n_exact += 1
            elif nm_source == 'library_adduct':
                _n_adduct += 1
            if nm_source is not None:
                rec['neutral_mass_source'] = nm_source
            if nm is None:
                continue
            rec['neutral_mass'] = nm
            nm_index.setdefault(int(nm / bucket_size), []).append(idx)
        for bucket in nm_index:
            nm_index[bucket].sort(key=lambda i: records[i].get('neutral_mass', 0))
        logger.info('  Neutral-mass index: %d spectra (exact_mass: %d, adduct fallback: %d)',
                    sum(len(v) for v in nm_index.values()), _n_exact, _n_adduct)

        _lib_records = records
        _mz_index = mz_index
        _nm_index = nm_index
        _lib_db_path = db_path
        _lib_source_families = tuple(sorted({
            normalize_source_family(record.get('dataset_source'))
            for record in records
            if normalize_source_family(record.get('dataset_source'))
        }))
        return records, mz_index


def _safe_float(v):
    if v is None:
        return None
    try:
        f = float(v)
        return f if not (f != f) else None
    except (TypeError, ValueError):
        return None


def _specific_adduct_polarity(adduct):
    """Return POS/NEG only for a supported, chemically specific adduct.

    ``Positive`` and ``Negative`` are common legacy values in the library's
    adduct column, but they do not establish whether an ion was protonated,
    sodiated, formate-adducted, etc.  They are adequate for polarity filtering,
    never for neutral-mass conversion.
    """
    text = str(adduct or '').strip()
    normalized = re.sub(r'\s+', '', text).upper()
    if normalized in _POLARITY_ONLY_ADDUCTS or normalized not in _SUPPORTED_SPECIFIC_ADDUCTS:
        return None
    polarity = _adduct_polarity(text)
    if polarity not in ('POS', 'NEG'):
        return None
    # Use an arbitrary positive precursor only to verify that the local adduct
    # parser supports this exact notation.  The real precursor is converted by
    # ``_trusted_query_neutral_mass`` below.
    if _calc_neutral_mass(500.0, text) is None:
        return None
    return polarity


def _trusted_query_neutral_mass(query_pmz, query_adduct, query_polarity=None):
    """Return ``(neutral_mass, polarity)`` for trustworthy query metadata.

    Caller-supplied ``query_neutral_mass`` is deliberately insufficient: an
    explicit adduct is required so the conversion is reproducible and its ion
    polarity can be checked.  This prevents blind multi-hypothesis searches
    from turning every precursor into several unrelated neutral masses.
    """
    polarity = _specific_adduct_polarity(query_adduct)
    if polarity is None:
        return None, None
    declared = str(query_polarity or '').strip().upper()
    if declared in ('POS', 'NEG') and declared != polarity:
        return None, None
    neutral = _calc_neutral_mass(query_pmz, query_adduct)
    if neutral is None:
        return None, None
    return neutral, polarity


def _record_source_is_searchable(record):
    """Filter explicit non-reference tiers while keeping legacy test records."""
    tier = str(record.get('source_tier') or '').strip()
    return not tier or source_tier_is_searchable(tier)


def _audited_search_precursor_from_flags(quality_flags):
    prefix = 'theoretical_precursor_mz='
    for flag in quality_flags or []:
        text = str(flag or '').strip()
        if not text.startswith(prefix):
            continue
        value = _safe_float(text[len(prefix):])
        if value is not None and 0 < value < 5000:
            return value, 'structure_adduct_theoretical_index'
    return None, 'reported_precursor_mz'


def _record_index_precursor_mz(record):
    """Return the audited retrieval mass without hiding the reported value."""
    indexed = _safe_float(record.get('search_precursor_mz'))
    if indexed is not None and indexed > 0:
        return indexed
    return _safe_float(record.get('precursor_mz'))


def _structure_adduct_precursor_audit(
    record,
    *,
    query_polarity=None,
    maximum_ppm=50.0,
    absolute_floor_da=0.02,
):
    """Audit whether a measured library precursor can belong to its structure.

    A spectrum whose precursor is 134.9 m/z cannot identify a C30H48O
    ``[M-H]-`` structure merely because its fragments are similar.  Old public
    library rows with this kind of stale structure mapping survived import and
    were previously only confidence-capped after ranking.  This audit is a
    search eligibility gate for records carrying a validated full structural
    identity and an executable adduct.
    """

    observed = _safe_float(
        record.get('reported_precursor_mz')
        if record.get('reported_precursor_mz') is not None
        else record.get('precursor_mz')
    )
    adduct = str(record.get('adduct') or '').strip()
    inchikey = normalize_inchikey(record.get('inchikey'))
    if (
        observed is None
        or observed <= 0
        or not str(record.get('smiles') or '').strip()
        or not inchikey
    ):
        return {
            'status': 'not_applicable_unverified_structure',
            'hard_reject': False,
            'reported_precursor_mz': observed,
        }

    formula = resolve_structure_formula(
        str(record.get('smiles') or ''),
        str(record.get('formula') or ''),
        adduct,
    )
    neutral_exact_mass = formula_exact_mass(formula)
    if neutral_exact_mass is None:
        return {
            'status': 'not_applicable_formula_mass_unavailable',
            'hard_reject': False,
            'reported_precursor_mz': observed,
            'formula': formula,
        }

    specific_polarity = _specific_adduct_polarity(adduct)
    if not specific_polarity:
        mode = str(
            query_polarity or record.get('polarity') or ''
        ).strip().upper()
        mode = (
            'POS' if mode.startswith('POS')
            else 'NEG' if mode.startswith('NEG')
            else ''
        )
        if not mode:
            return {
                'status': 'not_applicable_missing_adduct_and_polarity',
                'hard_reject': False,
                'reported_precursor_mz': observed,
                'formula': formula,
            }
        hypotheses = []
        for candidate_adduct in _COMMON_STRUCTURE_ADDUCT_HYPOTHESES[mode]:
            candidate_mz = _calc_ion_mz(
                neutral_exact_mass, candidate_adduct
            )
            if candidate_mz is None or candidate_mz <= 0:
                continue
            tolerance_da = max(
                float(absolute_floor_da),
                candidate_mz * float(maximum_ppm) / 1e6,
            )
            error_da = observed - candidate_mz
            if abs(error_da) <= tolerance_da:
                hypotheses.append({
                    'adduct': candidate_adduct,
                    'theoretical_precursor_mz': candidate_mz,
                    'error_da': error_da,
                    'error_ppm': error_da / candidate_mz * 1e6,
                    'tolerance_da': tolerance_da,
                })
        if hypotheses:
            hypotheses.sort(key=lambda item: abs(item['error_da']))
            selected = hypotheses[0]
            return {
                'status': 'consistent_via_structure_formula_adduct_hypothesis',
                'hard_reject': False,
                'reported_precursor_mz': observed,
                'formula': formula,
                'reported_adduct': adduct or None,
                'inferred_adduct': selected['adduct'],
                'inferred_adduct_source': (
                    'candidate_structure_formula_precursor_unique_hypothesis'
                    if len(hypotheses) == 1
                    else 'candidate_structure_formula_precursor_best_hypothesis'
                ),
                'adduct_hypotheses': hypotheses,
                'theoretical_precursor_mz': selected[
                    'theoretical_precursor_mz'
                ],
                'error_da': selected['error_da'],
                'error_ppm': selected['error_ppm'],
                'tolerance_da': selected['tolerance_da'],
            }
        return {
            'status': 'gross_structure_formula_precursor_conflict_without_adduct',
            'hard_reject': True,
            'reported_precursor_mz': observed,
            'formula': formula,
            'reported_adduct': adduct or None,
            'tested_adducts': list(
                _COMMON_STRUCTURE_ADDUCT_HYPOTHESES[mode]
            ),
        }

    theoretical = _calc_ion_mz(neutral_exact_mass, adduct)
    if theoretical is None or theoretical <= 0:
        return {
            'status': 'not_applicable_formula_mass_unavailable',
            'hard_reject': False,
            'reported_precursor_mz': observed,
            'formula': formula,
        }

    delta_da = observed - theoretical
    delta_ppm = delta_da / theoretical * 1e6
    tolerance_da = max(
        float(absolute_floor_da),
        theoretical * float(maximum_ppm) / 1e6,
    )
    if abs(delta_da) <= tolerance_da:
        status = 'consistent'
        hard_reject = False
    else:
        # Some explicitly audited legacy rows retain a coarsely rounded source
        # precursor but carry a structure-derived index m/z.  Preserve that
        # traceable rescue (still confidence-capped downstream); do not grant
        # the same exception to an unflagged stale mapping.
        indexed = _safe_float(record.get('search_precursor_mz'))
        indexed_source = str(record.get('precursor_mz_source') or '')
        indexed_consistent = (
            indexed is not None
            and abs(indexed - theoretical) <= tolerance_da
            and indexed_source == 'structure_adduct_theoretical_index'
        )
        status = (
            'audited_theoretical_index_rescue'
            if indexed_consistent
            else 'gross_structure_adduct_precursor_conflict'
        )
        hard_reject = not indexed_consistent
    return {
        'status': status,
        'hard_reject': hard_reject,
        'reported_precursor_mz': observed,
        'theoretical_precursor_mz': theoretical,
        'formula': formula,
        'adduct': adduct,
        'error_da': delta_da,
        'error_ppm': delta_ppm,
        'tolerance_da': tolerance_da,
    }


def _candidate_mass_metrics(query_pmz, library_record, query_neutral_mass=None,
                            cross_adduct=False, ppm_tolerance=None,
                            da_tolerance=None):
    """Calculate mass score/ppm on the physically correct mass axis."""
    default_da_tolerance = (
        NEUTRAL_MASS_TOL_DA if cross_adduct else MZ_TOLERANCE_DA
    )
    if da_tolerance is None:
        da_tolerance = default_da_tolerance
    if cross_adduct:
        query_mass = _safe_float(query_neutral_mass)
        library_mass = _safe_float(library_record.get('neutral_mass'))
        basis = 'neutral_mass'
    else:
        query_mass = _safe_float(query_pmz)
        observed_library_mass = _safe_float(library_record.get('precursor_mz'))
        indexed_library_mass = _safe_float(
            library_record.get('search_precursor_mz')
        )
        library_mass = indexed_library_mass
        basis = 'precursor_mz'
        if indexed_library_mass is not None:
            basis = str(
                library_record.get('precursor_mz_source')
                or 'audited_search_precursor_mz'
            )
        # Public libraries often round the reported precursor m/z even when a
        # structure-derived neutral exact mass and a specific adduct are
        # available.  Scoring against that rounded value can manufacture
        # 50-100+ ppm "High" hits.  Prefer the chemically reconstructable ion;
        # incompatible exact-mass/adduct metadata then becomes a large error
        # instead of silently receiving a strong mass score.
        if library_mass is None and (
            library_record.get('neutral_mass_source') in (
                'compound_exact_mass', DERIVED_NEUTRAL_MASS_SOURCE
            )
            and _specific_adduct_polarity(library_record.get('adduct'))
        ):
            library_mass = _calc_ion_mz(
                library_record.get('neutral_mass'), library_record.get('adduct')
            )
            if library_mass is not None:
                # Use theoretical m/z only when it is compatible with the
                # reported library precursor within the same retrieval window.
                # A gross conflict usually means a wrong adduct label, salt or
                # stale compound mapping; keep the measured spectral lead but
                # surface the conflict so it cannot become exact-structure High.
                internal_tol = max(
                    float(MZ_TOLERANCE_DA),
                    library_mass * float(
                        MZ_TOLERANCE_PPM
                    ) / 1e6,
                )
                if (
                    observed_library_mass is not None
                    and abs(observed_library_mass - library_mass) > internal_tol
                ):
                    library_mass = observed_library_mass
                    basis = 'precursor_mz_metadata_conflict'
                else:
                    basis = 'theoretical_ion_mz'
        if library_mass is None:
            library_mass = observed_library_mass
    if query_mass is None or library_mass is None:
        return 0.0, None, basis, query_mass, library_mass
    return (
        mass_score(query_mass, library_mass, ppm_tol=ppm_tolerance,
                   da_tol=da_tolerance),
        mass_ppm(query_mass, library_mass),
        basis,
        query_mass,
        library_mass,
    )


def _query_candidates(pmz, records, mz_index, tol_da=None, tol_ppm=None):
    """按 m/z 容差查找候选（精确范围，避免无效遍历）。"""
    if tol_da is None:
        tol_da = MZ_TOLERANCE_DA
    if tol_ppm is None:
        tol_ppm = MZ_TOLERANCE_PPM
    # Combine the configured small absolute floor with the relative HRMS
    # window.  The production floor is 0.02 Da (not the historical 0.5 Da).
    tol_from_ppm = pmz * tol_ppm / 1e6 if pmz > 0 else tol_da
    effective_tol = max(tol_da, tol_from_ppm)
    if hasattr(mz_index, "range_query"):
        # searchsorted provides the coarse closed interval, but the legacy
        # bucket path applies one final Python-float ``abs(diff) <= tol``.
        # Reapply that exact predicate so values on a binary-representation
        # boundary cannot create a candidate-set difference at cutover.
        return [
            idx
            for idx in mz_index.range_query(pmz, effective_tol)
            if (
                (ref_mz := _record_index_precursor_mz(records[idx])) is not None
                and abs(pmz - ref_mz) <= effective_tol
            )
        ]
    # 只搜索容差范围内的 bucket（+1 余量）
    bucket_range = max(1, int(effective_tol / MZ_BUCKET_DALTON) + 1)
    bucket = int(pmz / MZ_BUCKET_DALTON)
    candidates = []
    for b in range(bucket - bucket_range, bucket + bucket_range + 1):
        if b not in mz_index:
            continue
        for idx in mz_index[b]:
            ref_mz = _record_index_precursor_mz(records[idx])
            if ref_mz is None:
                continue
            diff = abs(pmz - ref_mz)
            if diff <= effective_tol:
                candidates.append(idx)
    # Never cap here.  Source tier and polarity have not been checked yet, so
    # an early nearest-mass cap can be filled by predicted/wrong-mode spectra.
    # search() filters those dimensions first and only then applies the dense
    # MS/MS pre-screen budget.
    return candidates


def _library_peaks_for_config(
    record: dict,
    peak_config: PeakProcessingConfig,
) -> list:
    """Return a library view processed identically to the query view.

    The loader's cached ``clean_peaks`` is valid only for the calibrated
    Top15/5% default.  Reusing it for a user-selected Top40/1% query made the
    comparison asymmetric.  An explicitly empty cached view also stays empty;
    it must never fall back to unfiltered raw peaks.
    """

    if peak_config.is_legacy_default and "clean_peaks" in record:
        return record.get("clean_peaks") or []
    cleaned, _ = process_peaks(
        record.get("peaks") or [],
        record.get("precursor_mz"),
        config=peak_config,
        precursor_charge=record.get("charge") or record.get("precursor_charge"),
    )
    return cleaned


def _prescreen_dense_candidates(candidates, records, query_clean, query_pmz,
                                limit, neutral_mass_queries=None,
                                peak_config=None, query_polarity=None):
    """Select a diverse dense-window candidate set by spectral relevance.

    The old dense-window path chose the nearest precursor spectrum for each compound
    and truncated to ``limit`` *before* calculating MS2 similarity.  In a window with
    hundreds of collision-energy replicates/isobars this can discard the only useful
    library spectrum.  Score the full polarity-filtered window with ModifiedCosine,
    retain the best spectrum from each compound before filling the remaining
    budget in rounds with additional spectra from those compounds.  ``cosine_score`` calls the same matchms numba kernels
    directly, preserving exact scores without per-candidate Spectrum object creation.

    Returns ``(selected_indices, cosine_cache)``.  Sparse windows are intentionally
    handled by the caller without this pre-screen so their established behaviour is
    unchanged.
    """
    candidates = [
        idx for idx in candidates
        if not _structure_adduct_precursor_audit(
            records[idx], query_polarity=query_polarity
        ).get('hard_reject')
    ]
    neutral_mass_queries = neutral_mass_queries or []

    def compound_key(i):
        rec = records[i]
        return compound_identity_key(rec) or f"sid::{rec.get('spectrum_id')}"

    def mass_distance(i):
        rec = records[i]
        indexed_mz = _record_index_precursor_mz(rec)
        distance = (
            abs(indexed_mz - query_pmz)
            if indexed_mz is not None else float('inf')
        )
        candidate_nm = rec.get('neutral_mass')
        if candidate_nm and neutral_mass_queries:
            distance = min(
                distance,
                min(abs(candidate_nm - value) for value in neutral_mass_queries),
            )
        return distance

    resolved_peak_config = peak_config or resolve_peak_processing_config()
    library_peak_lists = [
        _library_peaks_for_config(records[idx], resolved_peak_config)
        for idx in candidates
    ]
    batch_cosines = batch_modified_cosine_scores(
        query_clean,
        library_peak_lists,
        query_pmz,
        [_record_index_precursor_mz(records[idx]) for idx in candidates],
    )
    cosine_cache = {}
    best_by_compound = {}
    pools = {}
    for idx, lib_clean, cosine in zip(candidates, library_peak_lists, batch_cosines):
        rec = records[idx]
        cosine = float(cosine) if lib_clean else 0.0
        cosine_cache[idx] = cosine
        source_priority = _SOURCE_TIER_PRIORITY.get(
            str(rec.get('source_tier') or '').strip().lower(), 0
        )
        # Ascending total order: best similarity/mass/source first, then a
        # stable accession.  This is independent of SQLite row iteration.
        rank = (
            -cosine,
            mass_distance(idx),
            -source_priority,
            str(rec.get('spectrum_id') or ''),
        )
        key = compound_key(idx)
        pools.setdefault(key, []).append((rank, idx))
        previous = best_by_compound.get(key)
        if previous is None or rank < previous[0]:
            best_by_compound[key] = (rank, idx)

    representatives = [item[1] for item in best_by_compound.values()]
    representatives.sort(
        key=lambda i: (
            -cosine_cache[i],
            mass_distance(i),
            -_SOURCE_TIER_PRIORITY.get(
                str(records[i].get('source_tier') or '').strip().lower(), 0
            ),
            compound_key(i),
            str(records[i].get('spectrum_id') or ''),
        ),
    )
    queues = [sorted(pools[compound_key(i)]) for i in representatives]
    selected = []
    for depth in range(max(map(len, queues), default=0)):
        for queue in queues:
            if depth < len(queue): selected.append(queue[depth][1])
            if len(selected) == limit: break
        if len(selected) == limit: break
    return selected, {idx: cosine_cache[idx] for idx in selected}


def _neutral_candidates(neutral_mass_q, records, nm_index, tol_da=None,
                        tol_ppm=None):
    """准确度修复 #1: 按"中性质量"找候选 (跨加合物)。
    让 [M+Na]+/[M+NH4]+ 查询命中库内 [M+H]+/[M-H]- 参比 (同中性质量, 不同加合物)。"""
    if neutral_mass_q is None or neutral_mass_q <= 0 or not nm_index:
        return []
    if tol_da is None:
        tol_da = NEUTRAL_MASS_TOL_DA
    if tol_ppm is None:
        tol_ppm = MZ_TOLERANCE_PPM
    effective_tol = max(tol_da, neutral_mass_q * tol_ppm / 1e6)
    if hasattr(nm_index, "range_query"):
        return [
            idx
            for idx in nm_index.range_query(neutral_mass_q, effective_tol)
            if abs(records[idx].get('neutral_mass', -999) - neutral_mass_q)
            <= effective_tol
        ]
    bucket_size = MZ_BUCKET_DALTON
    bucket_range = max(1, int(effective_tol / bucket_size) + 1)
    bucket = int(neutral_mass_q / bucket_size)
    out = []
    for b in range(bucket - bucket_range, bucket + bucket_range + 1):
        for idx in nm_index.get(b, []):
            if abs(records[idx].get('neutral_mass', -999) - neutral_mass_q) <= effective_tol:
                out.append(idx)
    return out


def build_deconvolution_run_pool(
    spectra,
    *,
    min_purity=CHIMERIC_RUN_COMPONENT_MIN_PURITY,
    min_tic=CHIMERIC_RUN_COMPONENT_MIN_TIC,
):
    """Build a DecoID-style on-the-fly component library from clean DDA scans.

    Membership uses only acquisition evidence (isolation purity, TIC, polarity,
    precursor m/z), never the eventual library answer.  A contamination ceiling
    of 10% mirrors the conservative DDA unknown-library principle.
    """
    records = []
    mz_index = {}
    for index, spectrum in enumerate(spectra or []):
        isolation = spectrum.get("precursor_isolation_evidence") or {}
        try:
            purity = float(isolation.get("purity"))
        except (TypeError, ValueError):
            continue
        if isolation.get("status") != "clean" or purity < float(min_purity):
            continue
        raw_peaks = spectrum.get("peaks") or []
        tic = sum(float(peak[1]) for peak in raw_peaks if len(peak) >= 2)
        if tic < float(min_tic):
            continue
        precursor_mz = _safe_float(spectrum.get("precursor_mz"))
        if precursor_mz is None or precursor_mz <= 0:
            continue
        peaks = clean_peaks(
            raw_peaks, precursor_mz, top_n=60, rel_threshold=0.005,
            precursor_charge=spectrum.get("precursor_charge"),
        )
        if len(peaks) < 3:
            continue
        record = {
            "id": f"run:{spectrum.get('spectrum_id', index)}",
            "spectrum_id": str(spectrum.get("spectrum_id", index)),
            "precursor_mz": precursor_mz,
            "peaks": peaks,
            "polarity": spectrum.get("polarity"),
            "purity": purity,
            "tic": tic,
            "role": "same_run_clean_spectrum",
        }
        pool_index = len(records)
        records.append(record)
        bucket = int(precursor_mz / MZ_BUCKET_DALTON)
        mz_index.setdefault(bucket, []).append(pool_index)
    return records, mz_index


def _deconvolution_component_candidates(
    query_peaks,
    query_pmz,
    best,
    isolation_evidence,
    records,
    mz_index,
    *,
    query_polarity,
    mz_da,
    mz_ppm,
    max_components,
    run_component_pool=None,
):
    """Build a small, diverse component set from genuine co-isolated MS1 ions."""
    target_id = f"target:{best.get('lib_spectrum_id') or 'selected'}"
    components = [{
        "id": target_id,
        "peaks": best.get("lib_peaks") or [],
        "compound_name": best.get("compound_name") or "selected target",
        "precursor_mz": best.get("lib_precursor_mz"),
        "role": "selected_target",
    }]
    target_identity = compound_identity_key(best)
    proposals = []
    for co_rank, peak in enumerate(
        isolation_evidence.get("strongest_coisolated_peaks") or []
    ):
        if not isinstance(peak, (list, tuple)) or len(peak) < 2:
            continue
        try:
            co_mz = float(peak[0])
            co_intensity = float(peak[1])
        except (TypeError, ValueError):
            continue
        scored = []
        for idx in _query_candidates(co_mz, records, mz_index, mz_da, mz_ppm):
            record = records[idx]
            if not _record_source_is_searchable(record):
                continue
            if (query_polarity in {"POS", "NEG"}
                    and record.get("polarity") not in {None, query_polarity}):
                continue
            identity = compound_identity_key(record)
            if identity and identity == target_identity:
                continue
            lib_peaks = record.get("clean_peaks") or record.get("peaks") or []
            if len(lib_peaks) < 2:
                continue
            # Zero precursor shift is intentional: contaminant fragments remain
            # at their observed masses inside the chimeric query spectrum.
            similarity = cosine_score(
                query_peaks, lib_peaks,
                query_pmz=query_pmz, lib_pmz=query_pmz,
            )
            scored.append((
                similarity, str(record.get("spectrum_id") or ""), idx,
                identity, lib_peaks,
            ))
        for similarity, spectrum_id, idx, identity, lib_peaks in sorted(
            scored, key=lambda item: (-item[0], item[1], item[2])
        )[:2]:
            record = records[idx]
            proposals.append({
                "id": f"co:{idx}:{spectrum_id or co_rank}",
                "peaks": lib_peaks,
                "compound_name": record.get("compound_name") or spectrum_id,
                "precursor_mz": record.get("precursor_mz"),
                "coisolated_mz": co_mz,
                "coisolated_ms1_intensity": co_intensity,
                "prescreen_cosine": similarity,
                "identity_key": identity or f"spectrum:{idx}",
                "role": "coisolated_candidate",
            })

        run_records, run_index = run_component_pool or ([], {})
        run_scored = []
        for run_idx in (
            _query_candidates(co_mz, run_records, run_index, mz_da, mz_ppm)
            if run_records and run_index else []
        ):
            run_record = run_records[run_idx]
            if (query_polarity in {"POS", "NEG"}
                    and run_record.get("polarity") not in {None, query_polarity}):
                continue
            similarity = cosine_score(
                query_peaks, run_record["peaks"],
                query_pmz=query_pmz, lib_pmz=query_pmz,
            )
            run_scored.append((similarity, run_record["id"], run_record))
        for similarity, run_id, run_record in sorted(
            run_scored, key=lambda item: (-item[0], item[1])
        )[:1]:
            proposals.append({
                "id": run_id,
                "peaks": run_record["peaks"],
                "compound_name": "same-run unknown component",
                "precursor_mz": run_record["precursor_mz"],
                "coisolated_mz": co_mz,
                "coisolated_ms1_intensity": co_intensity,
                "prescreen_cosine": similarity,
                "identity_key": run_id,
                "role": "same_run_clean_spectrum",
                "source_purity": run_record.get("purity"),
                "source_tic": run_record.get("tic"),
            })

    proposals.sort(key=lambda item: (
        -float(item["coisolated_ms1_intensity"]),
        -float(item["prescreen_cosine"]), item["id"],
    ))
    seen = {target_identity} if target_identity else set()
    for proposal in proposals:
        identity = proposal.pop("identity_key")
        if identity in seen:
            continue
        seen.add(identity)
        components.append(proposal)
        if len(components) >= max(2, int(max_components)):
            break
    return target_id, components


def _apply_chimeric_deconvolution(
    best,
    query_peaks,
    query_pmz,
    isolation_evidence,
    records,
    mz_index,
    *,
    query_polarity,
    query_charge,
    mz_da,
    mz_ppm,
    score_weights,
    enabled,
    run_component_pool=None,
):
    """Apply guarded deconvolution to the selected target without reranking."""
    base = {
        "method": "DecoID-like non-negative L1 sparse deconvolution",
        "accepted": False,
        "runtime_influence": False,
        "score_contribution": 0.0,
        "candidate_comparison_policy": "raw_measured_spectrum_only",
        "candidate_margin_influence": False,
        "selected_identity_before_deconvolution": compound_identity_key(best),
    }
    if not enabled:
        best["chimeric_deconvolution"] = {**base, "status": "disabled"}
        return
    if not CHIMERIC_DECONVOLUTION_AVAILABLE:
        best["chimeric_deconvolution"] = {**base, "status": "module_unavailable"}
        return
    if not isinstance(isolation_evidence, dict):
        best["chimeric_deconvolution"] = {
            **base, "status": "isolation_evidence_unavailable",
        }
        return
    isolation_status = str(isolation_evidence.get("status") or "")
    if isolation_status not in {"mixed", "severe"}:
        best["chimeric_deconvolution"] = {
            **base,
            "status": (
                "clean_isolation_no_deconvolution"
                if isolation_status == "clean" else "isolation_window_unavailable"
            ),
            "isolation_status": isolation_status or "unavailable",
            "isolation_purity": isolation_evidence.get("purity"),
        }
        return

    extended_query = clean_peaks(
        query_peaks, query_pmz, top_n=60, rel_threshold=0.005,
        precursor_charge=query_charge,
    )
    target_id, components = _deconvolution_component_candidates(
        extended_query, query_pmz, best, isolation_evidence, records, mz_index,
        query_polarity=query_polarity, mz_da=mz_da, mz_ppm=mz_ppm,
        max_components=CHIMERIC_DECONVOLUTION_MAX_COMPONENTS,
        run_component_pool=run_component_pool,
    )
    component_audit = [
        {key: value for key, value in component.items() if key != "peaks"}
        for component in components
    ]
    if len(components) < 2:
        best["chimeric_deconvolution"] = {
            **base,
            "status": "no_coisolated_library_component",
            "isolation_status": isolation_status,
            "isolation_purity": isolation_evidence.get("purity"),
            "component_candidates": component_audit,
        }
        return

    evidence = {
        **base,
        **deconvolve_chimeric_spectrum(
            extended_query, components, target_id,
            penalty=CHIMERIC_DECONVOLUTION_L1_PENALTY,
            tolerance_da=min(max(float(mz_da), 0.005), 0.02),
            tolerance_ppm=max(float(mz_ppm), 20.0),
        ),
    }
    evidence.update({
        "isolation_status": isolation_status,
        "isolation_purity": isolation_evidence.get("purity"),
        "isolation_window_source": isolation_evidence.get("isolation_window_source"),
        "component_candidates": component_audit,
    })
    if not evidence.get("accepted"):
        evidence["score_increase_gates"] = {
            **(evidence.get("quality_gates") or {}),
            "spectral_gain": None,
        }
        evidence["all_score_increase_gates_passed"] = False
        best["chimeric_deconvolution"] = evidence
        return

    purified = evidence.get("purified_target_peaks") or []
    lib_peaks = best.get("lib_peaks") or []
    new_cosine = cosine_score(
        purified, lib_peaks,
        query_pmz=query_pmz, lib_pmz=best.get("lib_precursor_mz"),
    )
    new_entropy = entropy_score(purified, lib_peaks)
    new_fragment = fragment_score(purified, lib_peaks)
    new_intensity = intensity_score(purified, lib_peaks)
    new_reverse = reverse_match_ratio(
        purified, lib_peaks,
        query_pmz=query_pmz, lib_pmz=best.get("lib_precursor_mz"),
    )
    raw = {
        "cosine": float(best.get("cosine") or 0.0),
        "entropy": float(best.get("entropy") or 0.0),
        "fragment_score": float(best.get("fragment_score") or 0.0),
        "intensity_score": float(best.get("intensity_score") or 0.0),
        "reverse_match": float(best.get("reverse_match") or 0.0),
        "final_score": float(best.get("final_score") or 0.0),
        "rank_score": float(best.get("rank_score", best.get("final_score")) or 0.0),
    }
    updated = {
        "cosine": new_cosine,
        "entropy": new_entropy,
        "fragment_score": new_fragment,
        "intensity_score": new_intensity,
        "reverse_match": new_reverse,
    }
    cosine_weight, _mass_weight, fragment_weight, intensity_weight = score_weights
    spectral_gain = (
        cosine_weight * (new_cosine - raw["cosine"])
        + WEIGHT_ENTROPY * (new_entropy - raw["entropy"])
        + fragment_weight * (new_fragment - raw["fragment_score"])
        + intensity_weight * (new_intensity - raw["intensity_score"])
    ) * mass_error_penalty(best.get("mass_ppm"))
    applied_gain = min(
        max(0.0, spectral_gain), CHIMERIC_DECONVOLUTION_MAX_FINAL_GAIN
    )
    score_guard = (
        new_cosine >= raw["cosine"] + CHIMERIC_DECONVOLUTION_MIN_COSINE_GAIN
        and new_reverse >= max(0.50, raw["reverse_match"] - 0.05)
        and applied_gain >= CHIMERIC_DECONVOLUTION_MIN_FINAL_GAIN
    )
    evidence["score_increase_gates"] = {
        **(evidence.get("quality_gates") or {}),
        "spectral_gain": bool(score_guard),
    }
    evidence["all_score_increase_gates_passed"] = bool(
        evidence.get("all_pre_score_gates_passed")
        and all(
            gate is True
            for gate in evidence["score_increase_gates"].values()
        )
    )
    evidence["raw_scores"] = {key: round(value, 6) for key, value in raw.items()}
    evidence["deconvolved_scores"] = {
        **{key: round(value, 6) for key, value in updated.items()},
        "projected_final_score": round(
            min(1.0, raw["final_score"] + applied_gain), 6
        ),
    }
    if not score_guard:
        evidence.update({
            "accepted": False,
            "status": "deconvolved_without_guarded_score_gain",
            "runtime_influence": False,
            "score_contribution": 0.0,
        })
        best["chimeric_deconvolution"] = evidence
        return

    if not evidence["all_score_increase_gates_passed"]:
        # Defensive fail-closed guard: score changes are forbidden unless every
        # serialized pre-score gate and the spectral-gain gate are true.
        evidence.update({
            "accepted": False,
            "status": "score_gate_invariant_failed",
            "runtime_influence": False,
            "score_contribution": 0.0,
        })
        best["chimeric_deconvolution"] = evidence
        return

    best.update(updated)
    for key, value in raw.items():
        best[f"raw_{key}"] = value
    # Raw measured-spectrum ordering stays authoritative.  Deconvolution can
    # improve evidence for selected Top-1 but cannot manufacture a wider margin.
    best["rank_score"] = raw["rank_score"]
    best["final_score"] = min(1.0, raw["final_score"] + applied_gain)
    best["deconvolved_query_peaks"] = purified
    best["search_pass"] = f"{best.get('search_pass', 'primary_top15')}+chimeric_l1"
    evidence.update({
        "status": "accepted_score_applied",
        "accepted": True,
        "runtime_influence": True,
        "score_contribution": round(applied_gain, 6),
        "reranked_top1": False,
    })
    best["chimeric_deconvolution"] = evidence


def search(query_peaks, query_pmz, records=None, mz_index=None,
           blank_profile=None, blank_mode='absolute',
           blank_threshold=0.1, blank_int_mode='max',
           cos_thresh=None, final_thresh=None,
           mz_da=None, mz_ppm=None,
           w_cos=None, w_mass=None, w_frag=None, w_inten=None,
           top_n=None, rel_thresh=None,
           query_adduct=None, query_neutral_mass=None,
           query_adduct_reliable=None,
           query_polarity=None, query_rt=None, enforce_polarity=None,
           strict_evidence_gate=None, strict_max_mass_ppm=None,
           adaptive_peak_search=None, query_rt_unit=None,
           query_charge=None, query_collision_energy=None,
           query_instrument_type=None, query_instrument=None,
           query_activation=None, query_collision_energy_mode=None,
           query_instrument_domain_id=None,
           query_collision_cell_region=None,
           query_spectrum_representation=None,
           query_raw_peaks_complete=None,
           query_rt_method_id=None,
           query_isotope_evidence=None, isotope_evidence_status=None,
           query_isolation_evidence=None,
           enable_chimeric_deconvolution=None,
           deconvolution_run_pool=None,
           peak_processing_config=None, capture_peak_evidence=True):
    """
    搜索最佳匹配

    Args:
        query_peaks: 查询峰列表 [[mz, intensity], ...]
        query_pmz: 查询 precursor m/z
        records, mz_index: 库（None 则用已加载的全局库）
        blank_profile: BlankProfile (None 则跳过扣除)
        query_adduct: Explicit MGF/vendor adduct annotation.  Polarity alone
            is not enough to activate neutral-mass matching.
        query_adduct_reliable: Optional parser quality flag.  An explicit
            ``False`` (e.g. CHARGE/ADDUCT conflict) disables neutral matching.

    Returns:
        dict or None: 最佳匹配结果
    """
    if records is None or mz_index is None:
        records, mz_index = _lib_records, _mz_index
    if records is None:
        raise RuntimeError('Library not loaded. Call load_library() first.')
    if records is _lib_records and _lib_source_families:
        _searched_source_families = list(_lib_source_families)
    else:
        _searched_source_families = sorted({
            normalize_source_family(record.get('dataset_source'))
            for record in records
            if normalize_source_family(record.get('dataset_source'))
        })

    # 空白扣除
    if blank_profile is not None and BLANK_SUBTRACTION_AVAILABLE:
        query_peaks = subtract_blank(
            query_peaks, blank_profile,
            mode=blank_mode,
            intensity_threshold_factor=blank_threshold,
            blank_intensity_mode=blank_int_mode,
            precursor_mz=query_pmz,
            retention_time=query_rt,
            retention_time_unit=query_rt_unit,
        )

    _mz_da = mz_da if mz_da is not None else MZ_TOLERANCE_DA
    _mz_ppm = mz_ppm if mz_ppm is not None else MZ_TOLERANCE_PPM
    _cos_thresh = cos_thresh if cos_thresh is not None else COSINE_THRESHOLD
    _final_thresh = final_thresh if final_thresh is not None else FINAL_SCORE_THRESHOLD
    _query_instrument_family = normalize_instrument_family(
        query_instrument_type, query_instrument
    )
    _query_is_high_resolution = _query_instrument_family in {
        'orbitrap', 'qtof', 'fticr'
    }
    _w = (
        w_cos if w_cos is not None else WEIGHT_COSINE,
        w_mass if w_mass is not None else WEIGHT_MASS,
        w_frag if w_frag is not None else WEIGHT_FRAGMENT,
        w_inten if w_inten is not None else WEIGHT_INTENSITY,
    )
    _weight_override = any(v is not None for v in (w_cos, w_mass, w_frag, w_inten))
    if isinstance(peak_processing_config, PeakProcessingConfig):
        _peak_config = peak_processing_config.validated()
    elif isinstance(peak_processing_config, dict):
        _peak_config = resolve_peak_processing_config(
            peak_processing_config.get("preset"),
            relative_intensity_threshold=peak_processing_config.get(
                "relative_intensity_threshold"
            ),
            precursor_exclusion_da=peak_processing_config.get(
                "precursor_exclusion_da"
            ),
            top_n=peak_processing_config.get("top_n"),
            upper_mass_bound_enabled=peak_processing_config.get(
                "upper_mass_bound_enabled"
            ),
        )
    else:
        _peak_config = resolve_peak_processing_config(
            "library_match",
            relative_intensity_threshold=rel_thresh,
            top_n=top_n,
        )
    _top_n = _peak_config.top_n
    _rel_thresh = _peak_config.relative_intensity_threshold

    _query_nm, _adduct_mode = _trusted_query_neutral_mass(
        query_pmz, query_adduct, query_polarity
    )
    if query_adduct_reliable is False:
        _query_nm, _adduct_mode = None, None
    if query_polarity not in ('POS', 'NEG') and _adduct_mode in ('POS', 'NEG'):
        query_polarity = _adduct_mode

    direct_candidates = _query_candidates(
        query_pmz, records, mz_index, _mz_da, _mz_ppm
    )
    direct_candidate_set = set(direct_candidates)
    candidates = list(direct_candidates)

    # Cross-adduct matching is a single reproducible neutral-mass hypothesis.
    # ``query_neutral_mass`` alone is intentionally ignored: without a
    # trustworthy query adduct it cannot be audited and could have been guessed.
    q_nm_list = [_query_nm] if _query_nm is not None else []
    nm_cands = []
    if (USE_NEUTRAL_MASS_MATCHING and _query_nm is not None
            and records is _lib_records and _nm_index):
        nm_cands = _neutral_candidates(
            _query_nm, records, _nm_index,
            _mz_da, _mz_ppm,
        )
        seen = set(candidates)
        for i in nm_cands:
            if i not in seen:
                candidates.append(i)
                seen.add(i)

    _retrieval_direct_count = len(direct_candidates)
    _retrieval_neutral_count = len(nm_cands)

    # Source quality and ion mode are semantic eligibility checks, so both
    # happen before any candidate budget/cap.  Missing source_tier is allowed
    # for backward-compatible in-memory/evaluation records; an explicit
    # predicted/presence/quarantine tier is not an experimental library hit.
    candidates = [i for i in candidates if _record_source_is_searchable(records[i])]
    _retrieval_after_source_gate = len(candidates)
    # 离子模式过滤 (默认开): query 极性已知时, 剔除库 adduct 反推极性不一致的候选。
    # precursor m/z 匹配物理上无法区分离子模式 → 同一 m/z 既可能是 [M-H]- 也可能是 [M+H]+,
    # 负离子样本会查到正离子库谱图 (假阳性)。库 adduct 空/无法判定 → polarity=None → 保留(不误杀)。
    # enforce_polarity=False 或 query 极性读不到(None) → 不过滤 (零行为变化, 安全降级)。
    if enforce_polarity is None:
        enforce_polarity = ENFORCE_POLARITY_FILTER
    if enforce_polarity and query_polarity in ('POS', 'NEG'):
        candidates = [i for i in candidates
                      if records[i].get('polarity') in (None, query_polarity)]
        if not candidates:
            return None

    _retrieval_after_polarity_gate = len(candidates)

    if not candidates:
        return None

    candidate_set = set(candidates)
    cross_adduct_candidates = (
        set(nm_cands).intersection(candidate_set).difference(direct_candidate_set)
    )

    query_clean, _ = process_peaks(
        query_peaks,
        query_pmz,
        config=_peak_config,
        precursor_charge=query_charge,
    )

    # 候选多样性 (审计 #7): 密集质量窗里同一化合物几十张重复谱会占满 500 预算挤掉异构体。
    # 超 budget 时先对全部同极性候选计算轻量 Modified Cosine，再按化合物保留谱图最相似的
    # collision-energy/来源代表并截到昂贵评分预算。旧实现先按 m/z 最近截断、后算 MS2，可能
    # 静默漏掉稍远但真正相似的谱。稀疏窗(≤budget)不过滤，维持既有行为。
    _prescreen_cosine = {}
    if len(candidates) > MAX_CANDIDATES_PER_QUERY:
        candidates, _prescreen_cosine = _prescreen_dense_candidates(
            candidates, records, query_clean, query_pmz,
            MAX_CANDIDATES_PER_QUERY,
            neutral_mass_queries=q_nm_list,
            peak_config=_peak_config,
            query_polarity=query_polarity,
        )

    # v6 轻量多视图：只预处理峰表一次，不重复运行整套候选/专家/分类流程。
    # 候选循环中仅对质量误差≤50 ppm 且主 cosine 较弱者计算补充 cosine。
    _adaptive_enabled = USE_ADAPTIVE_PEAK_SEARCH if adaptive_peak_search is None else bool(adaptive_peak_search)
    _adaptive_views = []
    if _adaptive_enabled and _top_n < ADAPTIVE_TOP_N_PEAKS:
        if USE_BROAD_PEAK_RETRY:
            _broad_clean = clean_peaks(
                query_peaks, query_pmz,
                top_n=ADAPTIVE_TOP_N_PEAKS,
                rel_threshold=ADAPTIVE_REL_THRESHOLD,
                precursor_charge=query_charge,
            )
            if len(_broad_clean) > len(query_clean):
                _adaptive_views.append((f'adaptive_top{ADAPTIVE_TOP_N_PEAKS}', _broad_clean))
        if USE_LOW_MZ_INTERFERENCE_RETRY:
            _decontaminated = _drop_dominant_low_mz_interference(query_peaks, query_pmz)
            if _decontaminated:
                _decontam_clean = clean_peaks(
                    _decontaminated, query_pmz,
                    top_n=_top_n,
                    rel_threshold=_rel_thresh,
                    precursor_charge=query_charge,
                )
                if _decontam_clean:
                    _adaptive_views.append(('adaptive_low_mz_interference', _decontam_clean))

    # P-B: 诊断碎片抢救 — clean_peaks 砍掉的低丰度诊断峰, 若匹配全局诊断离子集(kb_v3 ~134 m/z)
    # 则救回。纯增量; 关闭(USE_DIAGNOSTIC_RESCUE=False)→query_extended==query_clean, 零行为变化。
    #
    # ⚠ 选择性喂入 (分解实验证实): 抢救峰只进 query_extended, 喂给 *纯存在性检查* 的评分
    #   (expert / expert_v3 / pathway_align — 它们只查 query 里有没有某 m/z, 不惩罚多余峰);
    #   不进 cosine/entropy/fragment/intensity/reverse_match。
    #   原因: 这五者做谱图匹配, query 有 lib 无的峰会受惩罚 — fragment Dice 按峰数稀释
    #   (实测 frag 1.0→0.488); cosine 平方强度点积被拖累; entropy(Li2021) 反上调低强度峰
    #   权重使 query-only 峰被放大惩罚 (实测 ent 1.0→0.812)。故抢救峰只进存在性检查类。
    query_extended = query_clean
    if USE_DIAGNOSTIC_RESCUE and DIAGNOSTIC_RESCUE_AVAILABLE and query_peaks:
        _pol = {'POS': 'pos', 'NEG': 'neg'}.get(query_polarity, 'both')
        try:
            _rescued = _rescue_diagnostic_peaks(
                query_peaks, query_pmz, _pol,
                mz_tol=DIAGNOSTIC_RESCUE_MZ_TOL,
                intensity_floor=DIAGNOSTIC_RESCUE_INTENSITY_FLOOR,
            )
            if _rescued:
                # 仅注入 query_clean 尚不存在的 m/z (容差内去重), 拼接后按 m/z 重排
                _have = [p[0] for p in query_clean]
                _inj = [p for p in _rescued
                        if not any(abs(p[0] - h) <= DIAGNOSTIC_RESCUE_MZ_TOL for h in _have)]
                if _inj:
                    query_extended = sorted(query_clean + _inj, key=lambda p: p[0])
        except Exception:
            pass

    # P1: 分类一致性过滤 — 解析 query 可靠大类 (DL 优先, 仅可靠类; 否则无模型碎片回退)。
    # per-query 算一次, 复用于所有候选; 同时缓存 DL 预测避免 best[] 处重复 forward。
    # 关闭时 _query_class=None → 循环内零行为变化, best[] 处回落到原懒加载逻辑。
    _query_class = None
    _query_clf_cached = None
    if USE_CONSISTENCY_FILTER and CONSISTENCY_FILTER_AVAILABLE:
        _clf_ok = _ensure_clf()
        _cw = _clf_wrapper if _clf_ok else None
        try:
            _query_class, _qsrc = _resolve_query_class(
                query_peaks, query_pmz, _cw,
                min_dl_conf=CONSISTENCY_MIN_DL_CONF,
                min_class_acc=CONSISTENCY_MIN_CLASS_ACC,
                reliable_whitelist=CONSISTENCY_RELIABLE_CLASSES,
            )
            if _cw is not None:
                try:
                    _query_clf_cached = _cw.predict(query_peaks, query_pmz)
                except Exception:
                    _query_clf_cached = None
        except Exception:
            _query_class = None

    best = None
    _topk_by_comp = {}   # v3 Phase 5 / Path1: Top-K 候选, 按化合物去重 (ik14 连接层, 保留每化合物最高 final), 供前端展示近邻/异构体
    _reference_evidence_by_comp = {}
    _topk_source_by_comp = {}  # Same representatives with private cleaned peaks for V7 shadow audit.
    _accepted_results = []
    for idx in candidates:
        lib = records[idx]
        _precursor_metadata_audit = _structure_adduct_precursor_audit(
            lib,
            query_polarity=query_polarity,
        )
        if _precursor_metadata_audit.get('hard_reject'):
            continue
        lib_clean = _library_peaks_for_config(lib, _peak_config)
        if not lib_clean:
            continue
        lib_search_pmz = _record_index_precursor_mz(lib)

        _cross_adduct_match = idx in cross_adduct_candidates
        ms, ppm_err, _mass_match_basis, _mass_query_value, _mass_library_value = (
            _candidate_mass_metrics(
                query_pmz,
                lib,
                query_neutral_mass=_query_nm,
                cross_adduct=_cross_adduct_match,
                ppm_tolerance=_mz_ppm,
                da_tolerance=_mz_da,
            )
        )
        if ppm_err is None:
            continue
        # On an explicitly high-resolution acquisition, the configured ppm
        # window is a hard identity-search boundary.  The 0.02 Da setting is a
        # low-resolution/rounding fallback, not permission for a 33--150 ppm
        # QE/QTOF candidate to enter ranking.
        if (
            _query_is_high_resolution
            and abs(float(ppm_err)) > float(_mz_ppm)
        ):
            continue

        _score_query = query_clean
        _search_pass = 'primary_top15'
        cos = _prescreen_cosine.get(idx)
        if cos is None:
            cos = cosine_score(query_clean, lib_clean,
                               query_pmz=query_pmz,
                               lib_pmz=lib_search_pmz)
        _ppm_pre = ppm_err
        if (_adaptive_views and cos < ADAPTIVE_RETRY_MAX_COSINE
                and _ppm_pre is not None and abs(_ppm_pre) <= ADAPTIVE_MAX_MASS_PPM):
            for _view_name, _view_peaks in _adaptive_views:
                _view_cos = cosine_score(
                    _view_peaks, lib_clean,
                    query_pmz=query_pmz, lib_pmz=lib_search_pmz,
                )
                if _view_cos >= cos + ADAPTIVE_MIN_COSINE_GAIN:
                    cos = _view_cos
                    _score_query = _view_peaks
                    _search_pass = _view_name
        # 碎片数门槛(根治单碎片虚高cos): ModifiedCosine对稀疏谱图会给虚高cos——
        # 如5-Methoxyflavanone: query/lib各仅1个碎片(151.00 vs 151.04), 单碎片完美匹配
        # cos=1.0, final被拉到0.845(假高分)。仅对 min_pk<3(单/双碎片)降cos(floor 0.3),
        # 保留 3+碎片命中(避免误伤Aucubin等少碎片标品, 验证: 红线11/14不破)。
        _min_pk = min(len(_score_query), len(lib_clean))
        # 仅惩罚"少碎片(min_pk<3) + cos虚高(>0.8)"的组合 = 单/双碎片完美匹配的虚高特征
        # (如5-Methoxyflavanone min_pk=1 cos=1.0)。不误伤"少碎片+cos边缘(~0.4)"的真实命中
        # (如Aucubin, 红线11/14不破)。cos<=0.8的少碎片命中保留(信息虽少但是真实相似)。
        if _min_pk < 3 and cos > 0.8:
            # 三重条件: 少碎片 + cos虚高 + PPM>100(质量明显不符) = 假阳性(如5-Methoxyflavanone
            # min_pk=1/cos=1.0/ppm=494)。质量符(PPM小)的少碎片高cos保留(如Aucubin标品)。
            _ppm_ck = ppm_err
            if _ppm_ck is not None and _ppm_ck > 100:
                cos *= max(0.3, _min_pk / 3.0)
        if cos < _cos_thresh:
            continue

        ent = entropy_score(_score_query, lib_clean, mz_tol=0.1)  # P-B: entropy 不喂抢救峰(ms-entropy上调低强度峰权重→query有lib无的峰反被放大惩罚, 实测0.812<1.0); 抢救峰只喂存在性检查类(expert/pathway)

        frag = fragment_score(_score_query, lib_clean)
        inten = intensity_score(_score_query, lib_clean)
        # 反向匹配率(NIST reverse dot product思路): lib峰在query找到的比例, 低=假阳性信号
        rev = reverse_match_ratio(
            _score_query,
            lib_clean,
            query_pmz=query_pmz,
            lib_pmz=lib_search_pmz,
        )

        # 预存馆藏碎裂路径是极性相关证据。三库回填后部分化合物只有正离子
        # 路径；负离子样本不得拿它评分或展示为本次路径。
        _candidate_pathway_json = lib.get('pathway_json')
        _candidate_pathway = None
        _pathway_status = 'missing'
        _pathway_polarities = []
        _pathway_adduct_status_value = 'missing'
        _pathway_adducts = []
        if _candidate_pathway_json:
            try:
                if PATHWAY_ALIGN_AVAILABLE:
                    _candidate_pathway = _parse_pathway(_candidate_pathway_json)
                elif isinstance(_candidate_pathway_json, str):
                    _candidate_pathway = json.loads(_candidate_pathway_json)
                elif isinstance(_candidate_pathway_json, dict):
                    _candidate_pathway = _candidate_pathway_json
                _pathway_status, _pathway_polarities = _pathway_polarity_status(
                    _candidate_pathway, query_polarity
                )
                _pathway_adduct_status_value, _pathway_adducts = (
                    _pathway_adduct_status(
                        _candidate_pathway,
                        query_adduct,
                        _query_nm is not None,
                    )
                )
                if (
                    _pathway_status == 'mismatch'
                    or _pathway_adduct_status_value
                    in {'mismatch', 'mixed_pathway_adducts'}
                ):
                    _candidate_pathway_json = None
            except Exception:
                _candidate_pathway = None
                _candidate_pathway_json = None
                _pathway_status = 'invalid'
                _pathway_adduct_status_value = 'invalid'

        _candidate_has_structure = bool(lib.get('has_structure'))
        _candidate_formula = resolve_structure_formula(
            str(lib.get('smiles') or ''),
            str(lib.get('formula') or ''),
            str(lib.get('adduct') or ''),
        )
        _structure_quality = 'missing'
        if _candidate_has_structure:
            _candidate_has_structure, _structure_quality = validate_structure_metadata(
                str(lib.get('smiles') or ''),
                str(lib.get('formula') or ''),
                str(lib.get('adduct') or ''),
            )

        expert_score = 0.0
        expert_details = None
        (
            _trusted_v3_class,
            _trusted_v3_reliable,
            _trusted_v3_reason,
        ) = trusted_analysis_class(lib, allow_legacy_missing=False)
        if _candidate_has_structure and EXPERT_SYSTEM_ENABLED and EXPERT_SYSTEM_AVAILABLE:
            try:
                expert_result = expert_evaluate(
                    query_peaks=_score_query,   # base expert (RDKit SMARTS, 随峰数变贵) 用最终评分视图
                    query_pmz=query_pmz,
                    candidate_smiles=lib['smiles'],
                    candidate_formula=_candidate_formula,
                    candidate_classification=lib.get('classification', 'Unknown'),
                    candidate_inchikey=lib.get('inchikey', ''),
                    query_polarity=query_polarity or '',
                    query_adduct=query_adduct or '',
                    query_charge=query_charge,
                    query_collision_energy=query_collision_energy,
                    candidate_adduct=lib.get('adduct') or '',
                    candidate_classification_source=lib.get(
                        'classification_source'
                    ),
                    candidate_classification_confidence=lib.get(
                        'classification_confidence'
                    ),
                )
                expert_score = expert_result.get('expert_total_score', 0.0)
                expert_details = expert_result
                (
                    _trusted_v3_class,
                    _trusted_v3_reliable,
                    _trusted_v3_reason,
                ) = trusted_analysis_class(
                    {**lib, 'expert_details': expert_result},
                    allow_legacy_missing=False,
                )
            except Exception:
                expert_score = 0.0
                expert_details = None

        # v3 class rules may influence production only when the class was
        # independently confirmed by structure/curation.  Raw ML or spectral
        # labels are recorded as not-applicable rather than silently scoring.
        if _candidate_has_structure and USE_EXPERT_V3 and EXPERT_V3_AVAILABLE:
            try:
                _v3 = _v3_score_and_reason(
                    _score_query,
                    query_pmz,
                    _trusted_v3_class,
                    lib.get('adduct'),
                    query_polarity=query_polarity,
                    query_charge=query_charge,
                    classification_reliable=_trusted_v3_reliable,
                    classification_reason=_trusted_v3_reason,
                )
                if _v3.get('v3_score', 0.0) > expert_score:
                    expert_score = _v3['v3_score']
                if expert_details is None:
                    expert_details = {}
                expert_details['v3'] = _v3
            except Exception:
                pass

        # Candidate pathways are mined from the candidate's own library
        # spectra.  Query/source exclusion was not guaranteed in the historic
        # artifact, so this is counterfactual shadow evidence only.  The hard
        # no-merge rule below remains in force even if a legacy env/config flag
        # requests pathway influence.
        if ((USE_PATHWAY_ALIGN_SHADOW or USE_PATHWAY_ALIGN)
                and PATHWAY_ALIGN_AVAILABLE
                and _candidate_pathway is not None
                and _pathway_status != 'mismatch'
                and _pathway_adduct_status_value
                not in {'mismatch', 'mixed_pathway_adducts'}):
            try:
                if _candidate_pathway:
                    _align, _align_detail = _pathway_align_score(
                        query_extended, query_pmz, _candidate_pathway,
                        mz_tol=PATHWAY_ALIGN_MZ_TOL,
                        min_support=PATHWAY_ALIGN_MIN_SUPPORT,
                        require_recurring=PATHWAY_ALIGN_REQUIRE_RECURRING,
                    )
                    _align_detail.update({
                        'candidate_conditioned': True,
                        'independent_evidence': False,
                        'eligible_for_fusion': False,
                        'runtime_influence': False,
                        'score_contribution': 0.0,
                        'status': 'shadow_only_query_exclusion_not_proven',
                        'legacy_influence_requested': bool(USE_PATHWAY_ALIGN),
                    })
                    if expert_details is None:
                        expert_details = {}
                    expert_details['pathway_align'] = _align_detail
            except Exception:
                pass

        if expert_details is not None:
            expert_details = finalize_expert_evidence(
                expert_details,
                fused_score=expert_score,
            )
            expert_score = float(
                expert_details.get('fused_expert_score', 0.0) or 0.0
            )
        _expert_fusion = (
            expert_details.get('score_fusion')
            if isinstance(expert_details, dict) else {}
        ) or {}

        # V11 keeps expert evidence for explanation and audit, but the flagship
        # identity selector cannot inherit the empirically negative expert
        # ablation.  A dedicated hard flag makes the zero contribution visible
        # and machine-auditable instead of relying on an accidental weight.
        expert_for_ranking = (
            expert_score if EXPERT_RANK_INFLUENCE_ENABLED else 0.0
        )
        final = enhanced_composite_score(
            cos, ms, frag, inten, expert_for_ranking, ent
        )
        if _weight_override:
            # Preserve the calibrated entropy/expert terms while honoring the
            # four documented per-call overrides.  Previously these arguments
            # were parsed into an unused tuple, so API callers saw no effect.
            final += (
                (_w[0] - WEIGHT_COSINE) * cos
                + (_w[1] - WEIGHT_MASS) * ms
                + (_w[2] - WEIGHT_FRAGMENT) * frag
                + (_w[3] - WEIGHT_INTENSITY) * inten
            )
            final = max(0.0, min(1.0, final))

        # 质量误差惩罚: 让 final_score 诚实反映 ppm 偏差 (解决"ppm>200 分数仍虚高")。
        # |ppm|≤50 无惩罚(自检索真匹配 ppm≈0 → 红线 11/14 不受影响); >200 强惩罚。
        # 同时作用于排序(质量不符候选下移)与展示分值。
        final *= mass_error_penalty(ppm_err)

        # 无结构馆藏只能作为同前体质量的高相似谱图线索。跨加合物、少碎片、
        # 低反向覆盖或普通相似度均拒绝，避免把 name-only 谱图冒充结构鉴定。
        if not _candidate_has_structure and not (
            not _cross_adduct_match
            and _min_pk >= 3
            and cos >= 0.85
            and rev >= 0.60
            and frag >= 0.25
            and ppm_err is not None
            and abs(ppm_err) <= 50.0
        ):
            continue

        # A literature table with reported fragment m/z but no measured
        # intensities is a target lead, not an experimental reference spectrum.
        # Keep only strong presence/coverage matches and mark them for a hard
        # confidence cap in honest_calibration.
        if lib.get('source_tier') == PRESENCE_ONLY_SOURCE_TIER and not (
            _min_pk >= 3
            and cos >= 0.85
            and rev >= 0.75
            and frag >= 0.35
            and ppm_err is not None
            and abs(ppm_err) <= 50.0
        ):
            continue

        # Legacy pathway rerank is deliberately suppressed.  Its requested
        # state is retained in provenance so an accidental enable is visible.
        if expert_details and 'pathway_align' in expert_details:
            expert_details['pathway_align']['rerank_requested'] = bool(
                USE_PATHWAY_RERANK and USE_PATHWAY_ALIGN
            )
            expert_details['pathway_align']['rerank_applied'] = False

        # P1: 分类一致性过滤 — query 可靠类 vs 候选类跨骨架冲突 → 软降权 / 硬否决
        _class_conflict = False
        if _query_class is not None:
            try:
                if _is_class_conflict(_query_class, lib.get('classification', 'Unknown')):
                    _class_conflict = True
                    if CONSISTENCY_MODE == 'hard':
                        continue
                    final *= CONSISTENCY_PENALTY
            except Exception:
                pass

        if final < _final_thresh:
            continue
        if not _passes_strict_evidence_gate(
                final, ent, rev, ppm_err, enabled=strict_evidence_gate,
                max_mass_ppm=strict_max_mass_ppm):
            continue
        rt_evidence = None
        if USE_RT_STANDARD_CONFIDENCE and RT_REFERENCE_AVAILABLE and query_rt is not None:
            try:
                rt_evidence = _match_rt_evidence(
                    lib,
                    query_rt,
                    query_method_id=query_rt_method_id,
                )
            except Exception:
                rt_evidence = None
        ce_profile_evidence = None
        if USE_CE_PROFILE_SHADOW and CE_PROFILE_AVAILABLE:
            try:
                ce_profile_evidence = _lookup_ce_profile(
                    inchikey=lib.get('inchikey'),
                    adduct=lib.get('adduct'),
                    ion_mode=query_polarity,
                    query_collision_energy=query_collision_energy,
                    query_peaks=_score_query,
                    query_precursor_mz=query_pmz,
                    query_instrument_type=query_instrument_type,
                    query_instrument=query_instrument,
                )
            except Exception:
                ce_profile_evidence = {
                    'status': 'lookup_error',
                    'runtime_influence': False,
                    'score_contribution': 0.0,
                }
        result = {
            'cosine': cos,
            'entropy': ent,
            'mass_score': ms,
            'mass_ppm': ppm_err,
            'mass_match_basis': _mass_match_basis,
            'mass_metadata_consistent': (
                True if _precursor_metadata_audit.get('status') == 'consistent'
                else False if (
                    _precursor_metadata_audit.get('status')
                    in {
                        'audited_theoretical_index_rescue',
                        'gross_structure_adduct_precursor_conflict',
                        'consistent_via_structure_formula_adduct_hypothesis',
                    }
                    or _mass_match_basis in {
                    'precursor_mz_metadata_conflict',
                    'structure_adduct_theoretical_index',
                    }
                )
                else True if _mass_match_basis in {'theoretical_ion_mz', 'neutral_mass'}
                else None
            ),
            'precursor_metadata_audit': _precursor_metadata_audit,
            'candidate_adduct_hypothesis': str(
                _precursor_metadata_audit.get('inferred_adduct') or ''
            ),
            'candidate_adduct_hypothesis_source': str(
                _precursor_metadata_audit.get('inferred_adduct_source') or ''
            ),
            'mass_query_value': _mass_query_value,
            'mass_library_value': _mass_library_value,
            'cross_adduct_match': _cross_adduct_match,
            'query_adduct': str(query_adduct or ''),
            'query_adduct_reliable': _query_nm is not None,
            'query_instrument_type': str(
                query_instrument_type or ''
            ),
            'query_instrument': str(query_instrument or ''),
            'query_instrument_family': _query_instrument_family,
            'query_activation': str(query_activation or ''),
            'query_collision_energy_mode': str(
                query_collision_energy_mode or ''
            ),
            'query_instrument_domain_id': str(
                query_instrument_domain_id or ''
            ),
            'query_collision_cell_region': str(
                query_collision_cell_region or ''
            ),
            'query_rt_method_id': str(query_rt_method_id or ''),
            'query_neutral_mass': _query_nm,
            'lib_neutral_mass': _safe_float(lib.get('neutral_mass')),
            'lib_neutral_mass_source': str(lib.get('neutral_mass_source') or ''),
            'fragment_score': frag,
            'intensity_score': inten,
            'reverse_match': round(rev, 4),
            'final_score': final,
            'rank_score': final,
            'expert_score': expert_score,
            'expert_rank_influence': bool(EXPERT_RANK_INFLUENCE_ENABLED),
            'expert_rank_contribution': round(
                WEIGHT_EXPERT * expert_for_ranking, 8
            ),
            'expert_details': expert_details,
            'expert_score_source': _expert_fusion.get('selected_source', 'none'),
            'expert_score_components': _expert_fusion.get('components', {}),
            'compound_name': lib['compound_name'] or lib['spectrum_id'],
            'source_compound_name': lib.get('source_compound_name', ''),
            'name_resolution_status': lib.get('name_resolution_status', ''),
            'smiles': lib['smiles'],
            'has_structure': _candidate_has_structure,
            'structure_quality': _structure_quality,
            'name_only_lead': bool(lib.get('name_only_lead') or not _candidate_has_structure),
            'inchikey': lib['inchikey'],
            'reported_inchikey': lib.get('reported_inchikey', ''),
            'formula': _candidate_formula,
            'adduct': lib['adduct'] or '',
            'lib_charge': lib.get('charge'),
            'lib_collision_energy': lib.get('collision_energy'),
            'lib_collision_energy_raw': lib.get('collision_energy_raw'),
            'lib_collision_energy_unit': lib.get('collision_energy_unit'),
            'lib_collision_energy_mode': lib.get('collision_energy_mode'),
            'lib_instrument_type': lib.get('instrument_type'),
            'lib_instrument': lib.get('instrument'),
            'lib_source_ms_level': lib.get('source_ms_level'),
            'lib_activation_method': lib.get('activation_method'),
            'lib_instrument_domain_id': lib.get('instrument_domain_id'),
            'lib_collision_cell_region': lib.get('collision_cell_region'),
            'lib_acquisition_method': lib.get('acquisition_method'),
            'lib_source_accession': lib.get('source_accession'),
            'ce_profile_evidence': ce_profile_evidence,
            'query_rt': query_rt,
            'lib_retention_time': lib.get('retention_time'),
            'rt_evidence': rt_evidence,
            'classification': lib['classification'],
            'classification_source': lib.get('classification_source', 'unknown'),
            'classification_confidence': lib.get('classification_confidence'),
            'classification_reliable': _trusted_v3_reliable,
            'classification_reliability_reason': _trusted_v3_reason,
            'dataset_source': lib.get('dataset_source', ''),
            'source_tier': (
                lib.get('source_tier', 'unknown')
                if _candidate_has_structure or not lib.get('has_structure')
                else 'metadata_conflict'
            ),
            'intensity_origin': lib.get('intensity_origin', 'unknown'),
            'quality_flags': lib.get('quality_flags') or [],
            'compound_role': lib.get('compound_role', ''),
            'potential_contaminant': bool(lib.get('potential_contaminant')),
            'contaminant_category': lib.get('contaminant_category', ''),
            'contaminant_note': lib.get('contaminant_note', ''),
            'exclude_from_natural_product_count': bool(
                lib.get('exclude_from_natural_product_count')
            ),
            'natural_product_candidate': bool(
                lib.get('natural_product_candidate', False)
            ),
            # P2-show: 透传预存裂解途径 JSON (Compounds.fragmentation_pathway) 到展示层。
            # 此前只在 :557 算 pathway_align 分后丢弃; 现带进 match → 序列化 → 报告直接读预存途径。
            'pathway_json': _candidate_pathway_json,
            'pathway_polarity_status': _pathway_status,
            'pathway_polarities': _pathway_polarities,
            'pathway_adduct_status': _pathway_adduct_status_value,
            'pathway_adducts': _pathway_adducts,
            'lib_precursor_mz': lib['precursor_mz'],
            'lib_search_precursor_mz': lib_search_pmz,
            'lib_precursor_mz_basis': _mass_match_basis,
            'lib_spectrum_id': lib['spectrum_id'],
            'lib_peaks': lib['peaks'],
            '_score_query_peaks': _score_query,
            '_score_library_peaks': lib_clean,
            'structure_score': 0.0,
            'class_conflict': _class_conflict,
            'search_pass': _search_pass,
        }
        result['feature_fragment_evidence'] = _evaluate_feature_fragment_shadow(
            # Feature-rule absences must be assessed against the complete
            # centroid peak list, not the Top-N/relative-intensity view used
            # for library scoring.  Passing the scoring view while declaring
            # raw_peaks_complete=True would manufacture false absences.
            query_peaks=query_peaks,
            query_pmz=query_pmz,
            candidate_class=_trusted_v3_class,
            classification_reliable=_trusted_v3_reliable,
            classification_source=_trusted_v3_reason,
            query_polarity=query_polarity,
            query_adduct=query_adduct,
            query_adduct_reliable=_query_nm is not None,
            query_charge=query_charge,
            query_collision_energy=query_collision_energy,
            query_collision_energy_mode=query_collision_energy_mode,
            query_instrument_type=query_instrument_type,
            query_instrument=query_instrument,
            query_activation=query_activation,
            query_instrument_domain_id=query_instrument_domain_id,
            query_collision_cell_region=query_collision_cell_region,
            query_spectrum_representation=query_spectrum_representation,
            query_raw_peaks_complete=query_raw_peaks_complete,
            query_isolation_evidence=query_isolation_evidence,
        )
        _accepted_results.append(result)

    # v6: 计算异构/同前体候选专属碎片证据。生产默认只展示，不改 rank_score；
    # 只有显式实验开关才恢复旧 rerank 行为。
    if _accepted_results:
        if USE_ISOMER_FRAGMENT_RERANK:
            try:
                _isomer_query_view = clean_peaks(
                    query_peaks,
                    query_pmz,
                    top_n=ISOMER_DIAGNOSTIC_QUERY_TOP_N,
                    rel_threshold=(
                        ISOMER_DIAGNOSTIC_QUERY_MIN_REL_INTENSITY
                    ),
                    precursor_charge=query_charge,
                )
                _apply_isomer_fragment_rerank(
                    _accepted_results,
                    _isomer_query_view or query_clean,
                    query_pmz=query_pmz,
                    query_collision_energy=query_collision_energy,
                    query_instrument_type=query_instrument_type,
                    query_instrument=query_instrument,
                )
            except Exception:
                pass
        _apply_formula_isotope_rerank(
            _accepted_results, query_isotope_evidence,
        )
        best = _select_best_result(_accepted_results)

        # IsoPair-MS is a formula-conditional decision layer, not a global
        # search score.  It may only replace the provisional winner with a
        # different 2-D connectivity inside that winner's same-formula,
        # auditable-adduct group, and only inside its locked near-tie domain.
        # Missing query adduct metadata may be resolved only conditionally on
        # the retrieved formula; that mode is shadow-only and cannot rerank or
        # support High.
        if USE_ISOPAIR_MS and ISOPAIR_RUNTIME_AVAILABLE and best is not None:
            try:
                _query_ce_value, _query_ce_unit = collision_energy_value_unit(
                    query_collision_energy
                )
                best, _isopair_evidence = _apply_isopair_ranker(
                    best,
                    _accepted_results,
                    query_peaks,
                    query_precursor_mz=query_pmz,
                    query_adduct=query_adduct,
                    query_adduct_reliable=query_adduct_reliable,
                    query_polarity=query_polarity,
                    query_charge=query_charge,
                    query_collision_energy=_query_ce_value,
                    query_collision_energy_unit=_query_ce_unit,
                    query_instrument_family=query_instrument_type,
                )
                if best is not None:
                    best['isopair_evidence'] = _isopair_evidence
            except Exception as exc:
                if best is not None:
                    best['isopair_evidence'] = {
                        'schema_version': 1,
                        'status': 'evaluation_failed',
                        'error': str(exc)[:160],
                        'evaluated': False,
                        'rank_influence_eligible': False,
                        'runtime_influence': False,
                        'ranking_changed': False,
                        'high_precision_gate_applicable': False,
                        'high_precision_support': False,
                        'high_eligible': False,
                        'high_promotion_allowed': False,
                        'identity_scope': '2d_connectivity_only',
                        'stereochemistry_claimed': False,
                    }

        # DecoID-like deconvolution is intentionally post-selection: raw MS/MS
        # chooses Top-1; a stable sparse fit may strengthen that target's
        # evidence but cannot use its own purified spectrum to rerank candidates.
        try:
            _apply_chimeric_deconvolution(
                best,
                query_peaks,
                query_pmz,
                query_isolation_evidence,
                records,
                mz_index,
                query_polarity=query_polarity,
                query_charge=query_charge,
                mz_da=_mz_da,
                mz_ppm=_mz_ppm,
                score_weights=_w,
                enabled=(
                    USE_CHIMERIC_DECONVOLUTION
                    if enable_chimeric_deconvolution is None
                    else bool(enable_chimeric_deconvolution)
                ),
                run_component_pool=deconvolution_run_pool,
            )
        except Exception as exc:
            best["chimeric_deconvolution"] = {
                "method": "DecoID-like non-negative L1 sparse deconvolution",
                "status": "evaluation_failed",
                "error": str(exc),
                "accepted": False,
                "runtime_influence": False,
                "score_contribution": 0.0,
            }

        # v3 Phase 5 / Path1: 维护 Top-K 候选, 按化合物去重 (同物多谱只保留最高 final,
        # 让 Top-K 显示不同化合物 & 异构体可检出)。放在 rerank 之后, 保证展示顺序与 best 一致。
        for result in _accepted_results:
            # Top-K is a frozen raw-spectrum comparison table.  The main match
            # may contain accepted deconvolution evidence, but those values may
            # not enter candidate ordering, margins or isomer competition.
            comparison_result = candidate_comparison_view(result)
            _ck = compound_identity_key(result)
            _reference_evidence_by_comp.setdefault(_ck, []).append({
                # Compact V15 evidence record.  Peak arrays stay in memory; the
                # stable content fingerprint is sufficient for mirror collapse
                # and avoids multiplying result-file size.
                'spectrum_id': result.get('lib_spectrum_id'),
                'source_accession': result.get('lib_source_accession'),
                'dataset_source': result.get('dataset_source', ''),
                'source_tier': result.get('source_tier', 'unknown'),
                'library_peak_fingerprint': reference_content_fingerprint(
                    result.get('_score_library_peaks')
                    or result.get('lib_peaks')
                    or []
                ),
                'raw_cosine': comparison_result.get('cosine'),
                'raw_entropy': comparison_result.get('entropy'),
                'raw_reverse_match': comparison_result.get('reverse_match'),
                'mass_score': comparison_result.get('mass_score'),
                'mass_ppm': comparison_result.get('mass_ppm'),
                'raw_fragment_score': comparison_result.get('fragment_score'),
                'raw_intensity_score': comparison_result.get('intensity_score'),
                'lib_collision_energy': result.get('lib_collision_energy'),
                'lib_collision_energy_raw': result.get(
                    'lib_collision_energy_raw'
                ),
                'lib_collision_energy_unit': result.get(
                    'lib_collision_energy_unit'
                ),
                'lib_collision_energy_mode': result.get(
                    'lib_collision_energy_mode'
                ),
                'lib_instrument_type': result.get('lib_instrument_type'),
                'lib_instrument': result.get('lib_instrument'),
                'lib_source_ms_level': result.get('lib_source_ms_level'),
                'lib_activation_method': result.get('lib_activation_method'),
                'lib_instrument_domain_id': result.get(
                    'lib_instrument_domain_id'
                ),
                'lib_collision_cell_region': result.get(
                    'lib_collision_cell_region'
                ),
            })
            _tk_entry = {
                'compound_name': result['compound_name'], 'smiles': result['smiles'],
                'compound_name_cn': result.get('compound_name_cn', ''),
                'source_compound_name': result.get('source_compound_name', ''),
                'name_resolution_status': result.get('name_resolution_status', ''),
                'has_structure': bool(result.get('has_structure')),
                'structure_quality': result.get('structure_quality', ''),
                'name_only_lead': bool(result.get('name_only_lead')),
                'inchikey': result['inchikey'],
                'reported_inchikey': result.get('reported_inchikey', ''),
                'formula': result['formula'],
                'final_score': round(
                    comparison_result.get('final_score', 0.0), 4
                ),
                'rank_score': round(
                    comparison_result.get(
                        'rank_score',
                        comparison_result.get('final_score', 0.0),
                    ),
                    4,
                ),
                'cosine': round(
                    comparison_result.get('cosine', 0.0), 4
                ),
                'candidate_comparison_source': comparison_result.get(
                    'candidate_comparison_source', 'measured_spectrum'
                ),
                'mass_ppm': (
                    round(result['mass_ppm'], 2)
                    if result.get('mass_ppm') is not None else None
                ),
                'mass_match_basis': result.get('mass_match_basis', 'precursor_mz'),
                'mass_metadata_consistent': result.get('mass_metadata_consistent'),
                'cross_adduct_match': bool(result.get('cross_adduct_match')),
                'query_adduct': result.get('query_adduct', ''),
                'query_adduct_reliable': bool(result.get('query_adduct_reliable')),
                'query_instrument_type': result.get(
                    'query_instrument_type', ''
                ),
                'query_instrument': result.get('query_instrument', ''),
                'query_rt_method_id': result.get(
                    'query_rt_method_id', ''
                ),
                'query_neutral_mass': result.get('query_neutral_mass'),
                'lib_neutral_mass': result.get('lib_neutral_mass'),
                'lib_neutral_mass_source': result.get('lib_neutral_mass_source', ''),
                'mass_query_value': result.get('mass_query_value'),
                'mass_library_value': result.get('mass_library_value'),
                'classification': result['classification'],
                'classification_source': result.get(
                    'classification_source', 'unknown'
                ),
                'classification_confidence': result.get(
                    'classification_confidence'
                ),
                'classification_reliable': bool(
                    result.get('classification_reliable')
                ),
                'classification_reliability_reason': result.get(
                    'classification_reliability_reason', ''
                ),
                'dataset_source': result.get('dataset_source', ''),
                'source_tier': result.get('source_tier', 'unknown'),
                'intensity_origin': result.get('intensity_origin', 'unknown'),
                'quality_flags': result.get('quality_flags') or [],
                'compound_role': result.get('compound_role', ''),
                'compound_role_note': result.get('compound_role_note', ''),
                'potential_contaminant': bool(result.get('potential_contaminant')),
                'contaminant_category': result.get('contaminant_category', ''),
                'exclude_from_natural_product_count': bool(
                    result.get('exclude_from_natural_product_count')
                ),
                'natural_product_candidate': bool(
                    result.get('natural_product_candidate', False)
                ),
                'adduct': result.get('adduct', ''),
                'lib_precursor_mz': result.get('lib_precursor_mz'),
                'lib_search_precursor_mz': result.get(
                    'lib_search_precursor_mz'
                ),
                'lib_precursor_mz_basis': result.get(
                    'lib_precursor_mz_basis'
                ),
                'lib_spectrum_id': result.get('lib_spectrum_id'),
                'reverse_match': round(
                    comparison_result.get('reverse_match', 0.0) or 0.0, 4
                ),
                # Per-candidate similarity evidence, exposed so entropy reranking
                # (weak-match regime) and calibrated-probability models can be
                # measured without re-scoring.  Read-only views of already-computed
                # terms; they change no score or ranking.
                'entropy': round(
                    comparison_result.get('entropy', 0.0) or 0.0, 4
                ),
                'fragment_score': round(
                    comparison_result.get('fragment_score', 0.0) or 0.0, 4
                ),
                'intensity_score': round(
                    comparison_result.get('intensity_score', 0.0) or 0.0, 4
                ),
                'lib_collision_energy': result.get('lib_collision_energy'),
                'lib_collision_energy_raw': result.get('lib_collision_energy_raw'),
                'lib_collision_energy_unit': result.get('lib_collision_energy_unit'),
                'lib_collision_energy_mode': result.get('lib_collision_energy_mode'),
                'lib_instrument_type': result.get('lib_instrument_type'),
                'lib_instrument': result.get('lib_instrument'),
                'lib_source_ms_level': result.get('lib_source_ms_level'),
                'lib_activation_method': result.get('lib_activation_method'),
                'lib_instrument_domain_id': result.get(
                    'lib_instrument_domain_id'
                ),
                'lib_collision_cell_region': result.get(
                    'lib_collision_cell_region'
                ),
                'lib_source_accession': result.get('lib_source_accession'),
                'ce_profile_evidence': result.get('ce_profile_evidence'),
                'feature_fragment_evidence': result.get(
                    'feature_fragment_evidence'
                ),
                'rt_evidence': result.get('rt_evidence'),
                'isomer_fragment_score': result.get('isomer_fragment_score'),
                'isomer_fragment_bonus': result.get('isomer_fragment_bonus'),
                'isomer_fragment_potential_bonus': result.get(
                    'isomer_fragment_potential_bonus'
                ),
                'isomer_fragment_rank_applied': bool(
                    result.get('isomer_fragment_rank_applied')
                ),
                'isomer_fragment_detail': result.get('isomer_fragment_detail'),
                'isopair_probability': result.get('isopair_probability'),
                'isopair_fused_score': result.get('isopair_fused_score'),
                'isopair_rank_applied': bool(
                    result.get('isopair_rank_applied')
                ),
                'isopair_selected_connectivity': bool(
                    result.get('isopair_selected_connectivity')
                ),
                'isopair_high_eligible': bool(
                    result.get('isopair_high_eligible')
                ),
                'isopair_high_precision_support': bool(
                    result.get('isopair_high_precision_support')
                ),
                'isopair_evidence': result.get('isopair_evidence'),
                'formula_isotope_evidence': result.get('formula_isotope_evidence'),
                'isotope_rank_adjustment': result.get('isotope_rank_adjustment'),
            }
            # Exact full InChIKey is the primary identity.  Name-first grouping
            # previously collapsed stereoisomers and made Top-K/Top-1 depend on
            # source labels.  Invalid/missing keys use the conservative shared
            # fallback implemented by ``compound_identity_key``.
            _prev = _topk_by_comp.get(_ck)
            if (_prev is None
                    or _stable_ranking_sort_key(_tk_entry)
                    < _stable_ranking_sort_key(_prev)):
                _topk_by_comp[_ck] = _tk_entry
                _topk_source_by_comp[_ck] = result

        for _candidate_key, _candidate in _topk_by_comp.items():
            _candidate['reference_evidence_records'] = sorted(
                _reference_evidence_by_comp.get(_candidate_key) or [],
                key=lambda record: (
                    str(record.get('dataset_source') or ''),
                    str(record.get('source_accession') or ''),
                    str(record.get('spectrum_id') or ''),
                ),
            )

        # V7 multi-metric evidence is deliberately post-selection and Top-K
        # bounded.  It cannot affect retrieval, acceptance, final_score,
        # rank_score, Top-1, or confidence calibration.
        if USE_MULTIMETRIC_SHADOW and best is not None:
            _shadow_sources = _top_candidates_with_best(
                best,
                _topk_source_by_comp.values(),
                TOP_K_CANDIDATES,
            )
            _best_key = compound_identity_key(best)
            _peak_view_shadows = {}
            for _candidate in _shadow_sources:
                _comparison = candidate_comparison_view(_candidate)
                _shadow = build_multimetric_shadow(
                    _candidate.get('_score_query_peaks') or query_clean,
                    _candidate.get('_score_library_peaks')
                    or _candidate.get('lib_peaks')
                    or [],
                    query_precursor_mz=query_pmz,
                    library_precursor_mz=_candidate.get('lib_precursor_mz'),
                    production_modified_cosine=_comparison.get('cosine'),
                    spectral_entropy=_comparison.get('entropy'),
                    reverse_coverage=_comparison.get('reverse_match'),
                    fragment_tolerance_da=MULTIMETRIC_FRAGMENT_TOL_DA,
                    strict_fragment_tolerance_da=(
                        MULTIMETRIC_STRICT_FRAGMENT_TOL_DA
                    ),
                )
                _peak_view = build_peak_view_shadow(
                    query_peaks,
                    _candidate.get('lib_peaks') or [],
                    query_precursor_mz=query_pmz,
                    library_precursor_mz=_candidate.get('lib_precursor_mz'),
                    query_precursor_charge=query_charge,
                    library_precursor_charge=_candidate.get('lib_charge'),
                    fragment_tolerance_da=MULTIMETRIC_FRAGMENT_TOL_DA,
                )
                _shadow['peak_view_shadow'] = _peak_view
                _candidate_key = compound_identity_key(_candidate)
                _peak_view_shadows[_candidate_key] = {
                    'compound_name': _candidate.get('compound_name', ''),
                    'evidence': _peak_view,
                }
                if _candidate_key in _topk_by_comp:
                    _topk_by_comp[_candidate_key][
                        'multimetric_shadow'
                    ] = _shadow
                if _candidate_key == _best_key:
                    best['multimetric_shadow'] = _shadow

            # This compares only the candidate similarity order inside the
            # already selected public Top-K. It is not a retrieval benchmark
            # and never mutates rank_score, Top-1, or confidence.
            _view_order_summary = {}
            for _view_name in PEAK_VIEW_PRESETS:
                _view_rows = []
                for _candidate_key, _payload in _peak_view_shadows.items():
                    try:
                        _view_score = float(
                            _payload['evidence']['views'][_view_name]
                            ['modified_cosine_hungarian']['score']
                        )
                    except (KeyError, TypeError, ValueError):
                        continue
                    _view_rows.append((
                        -_view_score,
                        _candidate_key,
                        _payload.get('compound_name', ''),
                    ))
                _view_rows.sort()
                for _rank, (_negative, _candidate_key, _name) in enumerate(
                    _view_rows, start=1
                ):
                    _payload = _peak_view_shadows[_candidate_key]['evidence']
                    _payload.setdefault('similarity_rank_by_view', {})[
                        _view_name
                    ] = _rank
                if _view_rows:
                    _winner = _view_rows[0]
                    _view_order_summary[_view_name] = {
                        'winner_identity_key': _winner[1],
                        'winner_compound_name': _winner[2],
                        'winner_similarity': round(-_winner[0], 6),
                        'production_top1_preserved': _winner[1] == _best_key,
                    }
            _best_shadow = _peak_view_shadows.get(_best_key, {}).get('evidence')
            if isinstance(_best_shadow, dict):
                _best_shadow['topk_similarity_order_by_view'] = _view_order_summary
                _best_shadow['top1_similarity_stable_across_views'] = bool(
                    _view_order_summary
                    and all(
                        item['production_top1_preserved']
                        for item in _view_order_summary.values()
                    )
                )
                _best_shadow['counterfactual_similarity_order_changed'] = any(
                    not item['production_top1_preserved']
                    for item in _view_order_summary.values()
                )

    # 模型结构预测得分（MS2 → 化学结构，与最佳候选对比）
    if best is not None and MODEL_PREDICTOR_AVAILABLE:
        try:
            pred = predict_from_spectrum(query_peaks)
            if pred and pred.get('valid') and pred.get('smiles'):
                ss = structure_similarity(pred['smiles'], best['smiles'])
                best['structure_score'] = ss
                best['predicted_smiles'] = pred['smiles']
        except Exception:
            pass

    # MS2Classifier 深度学习分类预测 (懒加载; P1 一致性过滤开启时复用循环前已缓存的预测)
    if best is not None:
        clf_result = _query_clf_cached  # P1: 一致性过滤已在循环前算过则直接复用
        if clf_result is None and _ensure_clf() and _clf_wrapper is not None:
            try:
                clf_result = _clf_wrapper.predict(query_peaks, query_pmz)
            except Exception as e:
                clf_result = None
                logger.debug('MS2Classifier prediction failed: %s', e)
        if clf_result is not None:
            best['dl_classification'] = clf_result['class']
            best['dl_confidence'] = clf_result['confidence']
            best['dl_top3'] = clf_result['top3']
            best['dl_probs'] = clf_result.get('probs', {})

    # ── 展示用分类 (display_class) resolution ──
    # 库 classification 若是 *具体* 化合物类(非 other/generic) → 直接用;
    # 否则(other/generic, 占库 ~70%)且 ML dl_classification 非 other 且置信度够高 → 用 ML 兜底。
    # 统一所有报告/表格读 display_class, 避免 86% 数据集名污染 + 提升 other 的覆盖率。
    if best is not None:
        try:
            from mssystem.classification import is_source_or_generic, CANONICAL_CLASSES
            lib_cls = best.get('classification') or ''
            if lib_cls and lib_cls != 'other' and (
                lib_cls in CANONICAL_CLASSES or not is_source_or_generic(lib_cls)
            ):
                best['display_class'] = lib_cls
                best['display_classification_source'] = best.get(
                    'classification_source', 'unknown'
                )
                best['display_classification_confidence'] = best.get(
                    'classification_confidence'
                )
            else:
                dl = best.get('dl_classification')
                dl_conf = best.get('dl_confidence', 0.0) or 0.0
                if dl and dl != 'other' and dl_conf >= 0.45:
                    best['display_class'] = dl
                    best['display_classification_source'] = 'msclass_spectrum'
                    best['display_classification_confidence'] = dl_conf
                else:
                    best['display_class'] = lib_cls or 'other'
                    best['display_classification_source'] = best.get(
                        'classification_source', 'unknown'
                    )
                    best['display_classification_confidence'] = best.get(
                        'classification_confidence'
                    )
        except Exception:
            best['display_class'] = best.get('classification') or 'other'
            best['display_classification_source'] = best.get(
                'classification_source', 'unknown'
            )
            best['display_classification_confidence'] = best.get(
                'classification_confidence'
            )

    # ── 记录专家评估日志 ──
    if best is not None and EXPERT_SYSTEM_ENABLED and EXPERT_SYSTEM_AVAILABLE and EXPERT_LOG_ENABLED:
        try:
            from .expert_logger import get_expert_logger
            elog = get_expert_logger()
            elog.log_evaluation(
                query_info={
                    'spectrum_id': query_peaks[0][0] if query_peaks else 'unknown',
                    'precursor_mz': query_pmz,
                    'num_peaks': len(query_peaks),
                },
                match_result=best,
                expert_scores=best.get('expert_details') or {'expert_total_score': 0.0},
            )
        except Exception:
            pass

    # ── Query precursor/adduct/MS1 evidence ──
    if best is not None:
        if _query_nm is not None:
            # This is metadata-backed conversion, not a spectrum-only guess.
            best['adduct_detected'] = str(query_adduct or '')
            best['neutral_mass'] = _query_nm
            best['adduct_detection_source'] = 'query_metadata'
        else:
            # ``detect_adduct`` cannot distinguish [M+H]+ from [M+Na]+ using
            # one precursor alone; its former first-item fallback falsely
            # labelled unknown queries as high-confidence protonated ions.
            best['adduct_detection_source'] = 'not_available'
        best['precursor_charge'] = query_charge
        best['collision_energy'] = query_collision_energy
        best['precursor_isotope_evidence'] = query_isotope_evidence
        best['precursor_isolation_evidence'] = query_isolation_evidence
        best['isotope_evidence_status'] = str(isotope_evidence_status or 'not_available')
        try:
            from utils.isotope import evaluate_formula_isotope_consistency
            best['formula_isotope_evidence'] = evaluate_formula_isotope_consistency(
                best.get('formula', ''), query_isotope_evidence
            )
        except Exception:
            best['formula_isotope_evidence'] = None

    # v3 Phase 5 / Path1: 附 Top-K 候选 (按化合物去重后, 降序, best 在首位) 供前端展示近邻/异构体
    if best is not None:
        # ``_select_best_result`` may use several deterministic near-tie rules.
        # Keep the public Top-K contract true:
        # the production best is always position zero, followed by the remaining
        # distinct identities.  Confidence margin code and external evaluators
        # rely on this ordering.
        best['top_candidates'] = _top_candidates_with_best(
            best, _topk_by_comp.values(), TOP_K_CANDIDATES
        )
        deconvolution = best.get("chimeric_deconvolution")
        if isinstance(deconvolution, dict):
            candidates = best["top_candidates"]
            selected_before = deconvolution.get(
                "selected_identity_before_deconvolution"
            )
            current_identity = compound_identity_key(best)
            top1_preserved = bool(
                candidates
                and selected_before
                and compound_identity_key(candidates[0])
                == current_identity
                == selected_before
            )
            deconvolution.update({
                "top1_identity_preserved": top1_preserved,
                "candidate_order_source": "raw_measured_spectrum",
                "deconvolved_candidate_margin_used": False,
            })
            if len(candidates) >= 2:
                deconvolution["raw_candidate_cosine_margin"] = round(
                    float(candidates[0].get("cosine", 0.0) or 0.0)
                    - float(candidates[1].get("cosine", 0.0) or 0.0),
                    6,
                )
                deconvolution["raw_candidate_rank_margin"] = round(
                    float(candidates[0].get(
                        "rank_score", candidates[0].get("final_score", 0.0)
                    ) or 0.0)
                    - float(candidates[1].get(
                        "rank_score", candidates[1].get("final_score", 0.0)
                    ) or 0.0),
                    6,
                )
            if deconvolution.get("accepted") and not top1_preserved:
                # Fail closed if a future refactor ever violates the immutable
                # raw Top-1 contract.
                for field in (
                    "cosine", "entropy", "fragment_score", "intensity_score",
                    "reverse_match", "final_score", "rank_score",
                ):
                    raw_value = best.get(f"raw_{field}")
                    if raw_value is not None:
                        best[field] = raw_value
                best.pop("deconvolved_query_peaks", None)
                deconvolution.update({
                    "status": "candidate_order_invariant_failed",
                    "accepted": False,
                    "runtime_influence": False,
                    "score_contribution": 0.0,
                })
        # Compound role is orthogonal to the learned broad class and must be
        # applied last so a contaminant cannot be relabelled as a natural class
        # by the optional MS2 classifier.
        annotate_compound_dict(best)
        for _candidate in best['top_candidates']:
            annotate_compound_dict(_candidate)
            if EXPERT_CONTRASTIVE_AVAILABLE:
                try:
                    _shadow = _contrastive_shadow(
                        query_extended,
                        query_pmz,
                        _candidate.get('classification', 'Unknown'),
                        query_polarity=query_polarity,
                        candidate_adduct=query_adduct or _candidate.get('adduct'),
                        query_charge=query_charge,
                    )
                    # Top-K can contain many candidates across thousands of
                    # spectra.  Keep comparison fields without duplicating the
                    # full rule provenance payload for every candidate.
                    _candidate['contrastive_expert_shadow'] = {
                        'status': _shadow.get('status'),
                        'candidate_class': _shadow.get('candidate_class'),
                        'polarity': _shadow.get('polarity'),
                        'adduct_specific': bool(_shadow.get('adduct_specific')),
                        'runtime_influence': False,
                        'score_contribution': 0.0,
                        'descriptive_evidence_weight': _shadow.get(
                            'descriptive_evidence_weight'
                        ),
                        'matched_ions_n': len(_shadow.get('matched_ions') or []),
                        'matched_losses_n': len(_shadow.get('matched_losses') or []),
                        'rule_set': _shadow.get('rule_set'),
                    }
                except Exception:
                    _candidate['contrastive_expert_shadow'] = {
                        'status': 'evaluation_failed',
                        'runtime_influence': False,
                        'score_contribution': 0.0,
                    }
        if EXPERT_CONTRASTIVE_AVAILABLE:
            try:
                best['contrastive_expert_shadow'] = _contrastive_shadow(
                    query_extended,
                    query_pmz,
                    best.get('classification', 'Unknown'),
                    query_polarity=query_polarity,
                    candidate_adduct=query_adduct or best.get('adduct'),
                    query_charge=query_charge,
                )
            except Exception:
                best['contrastive_expert_shadow'] = {
                    'status': 'evaluation_failed',
                    'runtime_influence': False,
                    'score_contribution': 0.0,
                }

    if best is not None:
        best.setdefault('search_pass', 'primary_top15')
        if capture_peak_evidence:
            try:
                from matching.spectral_evidence import attach_spectral_evidence
                attach_spectral_evidence(
                    best, _topk_source_by_comp.values(), query_precursor_mz=query_pmz,
                )
            except Exception as exc:
                unavailable = {
                    'schema': 'mssystem.spectral_evidence.v1',
                    'status': 'unavailable', 'reason': 'scoring_view_capture_failed',
                    'error_type': type(exc).__name__, 'ranking_influence': False,
                }
                for candidate in [best, *(best.get('top_candidates') or [])]:
                    candidate['spectral_evidence'] = dict(unavailable)
        # This explicit inventory lets V15 distinguish "searched with no
        # candidate reference" from "source participation unknown".  Absence
        # is display state only and is never converted to negative evidence.
        best['searched_source_families'] = _searched_source_families
        best['retrieval_audit'] = build_retrieval_audit(
            direct_candidates=_retrieval_direct_count,
            neutral_candidates=_retrieval_neutral_count,
            candidates_after_source_gate=_retrieval_after_source_gate,
            candidates_after_polarity_gate=_retrieval_after_polarity_gate,
            candidates_scored=len(candidates),
            reliable_specific_adduct=bool(_query_nm is not None),
            query_polarity=query_polarity,
            candidate_budget=MAX_CANDIDATES_PER_QUERY,
            mmap_index_used=hasattr(mz_index, 'range_query'),
        )
    return best


def find_matching_fragments(query_peaks, lib_peaks, tol=0.02,
                            query_pmz=None, lib_pmz=None):
    """找出直接或 precursor-shifted 匹配的碎片峰。"""
    if not query_peaks or not lib_peaks:
        return []
    q_mz = np.array([p[0] for p in query_peaks])
    matches = []
    mass_shift = None
    try:
        if query_pmz is not None and lib_pmz is not None:
            candidate_shift = float(query_pmz) - float(lib_pmz)
            if abs(candidate_shift) > tol:
                mass_shift = candidate_shift
    except (TypeError, ValueError):
        mass_shift = None
    for lmz, lint in lib_peaks:
        direct = np.abs(q_mz - lmz)
        candidates = [(float(direct.min()), int(np.argmin(direct)), 'direct')]
        if mass_shift is not None:
            shifted = np.abs(q_mz - (lmz + mass_shift))
            candidates.append((float(shifted.min()), int(np.argmin(shifted)), 'precursor_shift'))
        min_diff, min_idx, match_mode = min(candidates, key=lambda item: item[0])
        if min_diff <= tol:
            matches.append({
                'lib_mz': lmz, 'lib_intensity': lint,
                'query_mz': query_peaks[min_idx][0],
                'query_intensity': query_peaks[min_idx][1],
                'diff': min_diff,
                'match_mode': match_mode,
                'precursor_shift': mass_shift if match_mode == 'precursor_shift' else 0.0,
            })
    return matches
