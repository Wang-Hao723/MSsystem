"""Query-only adduct/charge checks before neutral-mass candidate lookup.

Positive bare values are charge magnitudes in the historical query schema.
Only explicit +/- strings or negative numeric values supply signed evidence.
This module neither introduces adduct aliases nor changes library identities.
"""
from __future__ import annotations
import inspect
import math
import re

from utils.adduct import NEUTRAL_SHIFTS
from matching.library_polarity import polarity_evidence

_SPECIFIC = {re.sub(r'\s+', '', key).upper(): value
             for key, value in NEUTRAL_SHIFTS.items()}


def _observed_charge(value):
    if value is None or (isinstance(value, str) and not value.strip()):
        return None, 'missing'
    if isinstance(value, bool):
        return None, 'invalid'
    text = str(value).strip()
    suffix = re.fullmatch(r'(\d+)([+-])', text)
    if suffix:
        text = suffix[2] + suffix[1]
    try:
        number = float(text)
    except (TypeError, ValueError, OverflowError):
        return None, 'invalid'
    if not math.isfinite(number) or not number.is_integer():
        return None, 'invalid'
    if number == 0:
        return None, 'unknown_zero'
    # Same absolute-charge comparison as utils.mgf_io, with malformed values
    # explicit instead of the MGF helper's arbitrary-digit extraction.
    return abs(int(number)), 'explicit'


def query_adduct_charge_evidence(adduct=None, charge=None, polarity=None):
    normalized = re.sub(r'\s+', '', str(adduct or '')).upper()
    specific = _SPECIFIC.get(normalized)
    expected = abs(int(specific[1])) if specific is not None else None
    magnitude, charge_status = _observed_charge(charge)
    polarity_info = polarity_evidence(
        adduct=adduct, charge=None if charge_status == "unknown_zero" else charge,
        ion_mode=polarity)
    reasons = []
    if charge_status == 'invalid':
        reasons.append('invalid_query_charge')
    if polarity_info['conflict']:
        reasons.append('query_ion_polarity_conflict')
    if expected is not None and magnitude is not None and magnitude != expected:
        reasons.append('explicit_adduct_charge_magnitude_conflict')
    return {
        'rule_version': 'query-adduct-charge-before-retrieval-v3',
        'reported_adduct': adduct, 'reported_charge': charge,
        'reported_polarity': polarity,
        'specific_adduct_supported': specific is not None,
        'expected_charge_magnitude': expected,
        'observed_charge_magnitude': magnitude,
        'charge_status': charge_status,
        'charge_polarity': polarity_info['charge_polarity'],
        'polarity_conflict': polarity_info['conflict'],
        'conflict': bool(reasons), 'reasons': reasons,
    }


def guard_query_record(record):
    """Preserve raw fields and prior parser rejection; never invent an adduct."""
    evidence = query_adduct_charge_evidence(
        record.get('adduct'), record.get('reported_precursor_charge', record.get('precursor_charge')),
        record.get('polarity'))
    updated = dict(record, adduct_charge_evidence=evidence)
    if evidence['conflict']:
        updated['adduct_reliable'] = False
    return updated


def guard_search_arguments(function, args, kwargs):
    """Validate positional and keyword metadata before the underlying search."""
    bound = inspect.signature(function).bind_partial(*args, **kwargs)
    evidence = query_adduct_charge_evidence(
        bound.arguments.get('query_adduct'), bound.arguments.get('query_charge'),
        bound.arguments.get('query_polarity'))
    if evidence['conflict']:
        bound.arguments['query_adduct_reliable'] = False
    return bound.args, bound.kwargs, evidence
