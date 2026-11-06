"""Parse explicit source MS-stage text for ordinary-MS2 fragment contrasts.

This describes the original source annotation. It does not construct a real
MSn acquisition tree or infer stage from instrument, peaks, names, or energy.
"""
from functools import lru_cache
import re


@lru_cache(maxsize=256)
def _parse(text):
    normalized = text.upper().replace('\u2013', '-').replace('\u2014', '-')
    single = re.fullmatch(r'MS\s*(\d+)', normalized)
    composite = re.fullmatch(r'MS\s*(\d+)\s*-\s*MS\s*(\d+)\s+COMPOSITE', normalized)
    if single:
        lower = upper = int(single[1])
        kind = 'single_stage'
    elif composite:
        lower, upper = int(composite[1]), int(composite[2])
        kind = 'composite'
        if lower > upper:
            return ('unrecognized', None, None, 'invalid_stage_range', False, False)
    else:
        return ('missing' if not text else 'unrecognized', None, None,
                None if not text else 'unrecognized_source_stage_text', False, False)
    if lower == upper == 2:
        return ('explicit_ms2_only', lower, upper, None, True, False)
    if upper > 2:
        anomaly = ('mixed_ms2_msn_composite' if kind == 'composite' and lower == 2
                   else 'higher_order_composite' if kind == 'composite'
                   else 'higher_order_spectrum')
        return ('explicit_non_ms2', lower, upper, anomaly, False, True)
    return ('explicit_non_ms2', lower, upper, 'source_stage_below_ms2', False, True)


def source_ms_stage(value):
    text = '' if value is None else str(value).strip()
    status, lower, upper, anomaly, explicit_ms2, disqualified = _parse(text)
    return {'source_ms_level': text, 'status': status, 'minimum_stage': lower,
            'maximum_stage': upper, 'source_anomaly_type': anomaly,
            'explicit_ms2_only': explicit_ms2,
            'ineligible_for_ordinary_ms2_contrast': disqualified}


def source_ms2_pair_status(left, right):
    a, b = source_ms_stage(left), source_ms_stage(right)
    if a['ineligible_for_ordinary_ms2_contrast'] or b['ineligible_for_ordinary_ms2_contrast']:
        status = 'explicit_non_ms2_reference'
    elif a['explicit_ms2_only'] and b['explicit_ms2_only']:
        status = 'both_references_explicit_ms2'
    else:
        status = 'source_ms_level_incomplete'
    return {'status': status, 'left': a, 'right': b,
            'comparable_for_ordinary_ms2': status != 'explicit_non_ms2_reference',
            'complete': status == 'both_references_explicit_ms2'}
