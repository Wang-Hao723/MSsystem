"""Acquisition checks at query boundaries, independent of library candidates."""
from __future__ import annotations
from collections import defaultdict
import importlib
import json
import math
import re
import types

from matching.library_polarity import polarity_evidence

_recap = importlib.import_module('matching.recap_adduct')
_precursor = importlib.import_module('matching.precursor_features')
_isopair = importlib.import_module('matching.isopair_ranker')


def charge_key(value):
    if value is None or value == '':
        return 'unknown'
    try:
        number = float(value)
    except (TypeError, ValueError):
        return 'invalid:' + str(value)
    if isinstance(value, bool) or not math.isfinite(number) or not number.is_integer():
        return 'invalid:' + str(value)
    return int(number) if number else 'unknown'


def canonical_adduct_key(value):
    return _recap._ADDUCT_KEY_BY_ALIAS.get(_recap._compact(value), '')


def resolve_isopair_adduct_context(**kwargs):
    reported_charge = kwargs.get('query_charge')
    numeric_charge = reported_charge
    # A bare positive value is a historical charge magnitude for either mode.
    # Explicit MGF signs (1+, 1-, +1, -1) retain their polarity evidence.
    if isinstance(reported_charge, str):
        match = re.fullmatch(r'(\d+)([+-])', reported_charge.strip())
        if match:
            numeric_charge = match[2] + match[1]
    charge = charge_key(numeric_charge)
    raw = str(kwargs.get('query_adduct') or '').strip()
    polarity = _recap._query_polarity(kwargs.get('query_polarity'))
    charge_polarity = (polarity_evidence(charge=reported_charge)['charge_polarity']
                       if isinstance(charge, int) else None)
    reason = None
    if isinstance(charge, str) and charge.startswith('invalid:'):
        reason = 'invalid_query_charge'
    elif isinstance(charge, int) and polarity and charge_polarity and charge_polarity != polarity:
        reason = 'query_polarity_charge_conflict'
    elif raw and _recap._compact(raw) not in {'POSITIVE','POS','+','NEGATIVE','NEG','-'}:
        if not canonical_adduct_key(raw):
            reason = 'unsupported_explicit_query_adduct'
        elif isinstance(charge, int) and abs(charge) != 1:
            # This trained RECAP stratum enumerates singly charged ions only.
            reason = 'explicit_adduct_charge_conflict'
    if reason:
        return {**_recap._base_context('unresolved', reason),
                'reported_query_adduct': raw, 'reported_query_charge': kwargs.get('query_charge'),
                'query_polarity': polarity}
    # Do not let the legacy resolver infer POS from an unsigned magnitude.
    resolved_polarity = polarity or charge_polarity or ''
    context = _recap.resolve_isopair_adduct_context(**{
        **kwargs, 'query_polarity': resolved_polarity,
        'query_charge': charge if isinstance(charge, int) and resolved_polarity else None,
    })
    if context.get('ranking_eligible') and isinstance(charge, int) and abs(charge) != 1:
        return {**_recap._base_context('unresolved', 'formula_adduct_charge_conflict'),
                'reported_query_adduct': raw, 'reported_query_charge': charge,
                'rejected_formula_conditioned_context': context}
    return context


def _clone(function, replacements):
    copied = types.FunctionType(function.__code__, {**function.__globals__, **replacements},
                                function.__name__, function.__defaults__, function.__closure__)
    copied.__kwdefaults__ = function.__kwdefaults__
    return copied


apply_isopair_ranker = _clone(_isopair.apply_isopair_ranker, {
    'resolve_isopair_adduct_context': resolve_isopair_adduct_context,
    'canonical_adduct_key': canonical_adduct_key,
})


def energy_condition(row):
    value = row.get('collision_energy')
    unit = row.get('collision_energy_unit') or ''
    if isinstance(value, dict):
        unit = value.get('unit') or unit
        value = value.get('value')
    def finite(v):
        try:
            n = float(v)
        except (ValueError, TypeError):
            return None
        return n if math.isfinite(n) and n >= 0 else None
    number = finite(value)
    values = row.get('collision_energy_values') or []
    values = sorted(set(n for n in map(finite, values) if n is not None))
    if number is None and not values:
        return None
    return {'value': number, 'unit': str(unit).strip(),
            'activation': str(row.get('activation_method') or '').strip(),
            'mode': str(row.get('collision_energy_mode') or '').strip(), 'values': values}


def _energy_metadata(rows):
    conditions = {}
    comparable = defaultdict(set)
    for row in rows:
        condition = energy_condition(row)
        if condition is None:
            continue
        conditions[json.dumps(condition, sort_keys=True)] = condition
        if condition['unit']:
            key = (condition['unit'], condition['activation'], condition['mode'])
            comparable[key].add((condition['value'], tuple(condition['values'])))
    return {'precursor_feature_collision_energy_conditions': list(conditions.values()),
            'precursor_feature_condition_count': len(conditions),
            'precursor_feature_multi_energy': any(len(values) >= 2 for values in comparable.values())}


def annotate_precursor_features(results, **kwargs):
    partitions = defaultdict(list)
    for index, row in enumerate(results):
        if isinstance(row, dict):
            partitions[charge_key(row.get('precursor_charge', row.get('charge')))].append((index,row))
    groups = []
    for charge, indexed in partitions.items():
        _precursor.annotate_precursor_features([row for _,row in indexed], **kwargs)
        local = defaultdict(list)
        for index,row in indexed:
            local[row['precursor_feature_id']].append((index,row))
        groups.extend((min(i for i,_ in values),charge,values) for values in local.values())
    for ordinal,(_,charge,indexed) in enumerate(sorted(groups,key=lambda x:x[0]),1):
        metadata = _energy_metadata([row for _,row in indexed])
        for _,row in indexed:
            row.update(metadata,precursor_feature_id=f'PF{ordinal:06d}',
                       precursor_feature_charge=charge if isinstance(charge,int) else None,
                       precursor_feature_charge_status='explicit' if isinstance(charge,int) else charge)
    return precursor_feature_summary(results)


def precursor_feature_summary(results):
    summary = _precursor.precursor_feature_summary(results)
    groups = defaultdict(list)
    for index,row in enumerate(results or []):
        if not isinstance(row,dict):
            continue
        key = ('lcms',row['lcms_feature_id']) if row.get('lcms_feature_id') else ('fallback',row.get('precursor_feature_id',index))
        groups[key].append(row)
    summary['multi_energy_precursor_features'] = sum(_energy_metadata(rows)['precursor_feature_multi_energy'] for rows in groups.values())
    summary['multiple_acquisition_condition_features'] = sum(_energy_metadata(rows)['precursor_feature_condition_count'] >= 2 for rows in groups.values())
    if not summary.get('lcms_feature_finding_available'):
        summary['method'] = 'MS2 precursor-event clustering (m/z + RT + polarity + explicit charge)'
    return summary
