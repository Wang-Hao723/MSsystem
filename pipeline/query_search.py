"""Current query validation over the unchanged accepted library processing code."""
import functools
import importlib
import copy
from collections import OrderedDict
from pipeline.structure_prediction import predict_from_spectrum
from pipeline.query_input_guards import _clone, apply_isopair_ranker
from pipeline.query_charge_guard import guard_search_arguments
from pipeline.query_record_memo import QueryRecordMemo
from pipeline.query_expert_memo import query_expert_memo

_base = importlib.import_module('matching.search')


def _reuse_isomer_reference_profiles(*args, **kwargs):
    """Reuse identical work within one annotation call, retaining every reference.

    The underlying routine adds annotation/rank fields but never changes the
    identity fields, library peaks, or precursor masses read by these helpers.
    Ordered row identities preserve multiplicity, order and condition membership.
    A new cache for every call also isolates spectra, settings and worker threads.
    """
    identities = {}
    profiles = OrderedDict()

    def identity(row):
        key = id(row)
        if key not in identities:
            identities[key] = _base._rerank_identity(row)
        return identities[key]

    def profile(rows, *, neutral_losses=False):
        rows = list(rows or [])
        key = (neutral_losses, tuple(id(row) for row in rows))
        if key not in profiles:
            profiles[key] = _base._reference_peak_profile(
                rows, neutral_losses=neutral_losses,
            )
            if len(profiles) > 256:
                profiles.popitem(last=False)
        profiles.move_to_end(key)
        # Preserve the original fresh-result contract if a consumer adds fields.
        return copy.deepcopy(profiles[key])

    function = _clone(_base._apply_isomer_fragment_rerank, {
        '_rerank_identity': identity,
        '_reference_peak_profile': profile,
    })
    return function(*args, **kwargs)


@functools.wraps(_base.search)
def search(*args, **kwargs):
    # Capture the current loaded library and flags per call; do not copy stale
    # module globals at import time or modify another query's shared namespace.
    args, kwargs, charge_evidence = guard_search_arguments(_base.search, args, kwargs)
    positional = list(args)
    records = positional[2] if len(positional) > 2 else kwargs.get('records')
    if records is None:
        records = getattr(_base, '_lib_records', None)
        replace_explicit_records = False
    else:
        replace_explicit_records = True
    query_records = None
    if (records is not None and hasattr(records, 'read_uncached')
            and hasattr(records, 'cache_info')):
        query_records = QueryRecordMemo(records)
        if replace_explicit_records:
            if len(positional) > 2:
                positional[2] = query_records
            else:
                kwargs['records'] = query_records
        args = tuple(positional)
    overrides = {
        '_apply_isopair_ranker': apply_isopair_ranker,
        'predict_from_spectrum': predict_from_spectrum,
        '_apply_isomer_fragment_rerank': _reuse_isomer_reference_profiles,
    }
    expert = getattr(_base, "expert_evaluate", None)
    if callable(expert):
        overrides["expert_evaluate"] = query_expert_memo(expert)
    # Preserve every identity-sensitive branch in matching.search. Production
    # passes the loaded mmap store explicitly, while that function compares it
    # with its module-global store before enabling accepted neutral-mass paths.
    if query_records is not None and records is getattr(_base, '_lib_records', None):
        overrides['_lib_records'] = query_records
    function = _clone(_base.search, overrides)
    result = function(*args, **kwargs)
    if isinstance(result, dict):
        result["query_adduct_charge_evidence"] = charge_evidence
    return result
