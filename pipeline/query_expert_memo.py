"""Query-scoped reuse of identical expert-evidence calls."""
from collections import OrderedDict
from copy import deepcopy
from functools import wraps
import math


class _Uncacheable(TypeError):
    pass


def _freeze(value):
    kind = type(value)
    if value is None:
        return ("none",)
    if kind in (str, int, bool, bytes):
        return (kind.__name__, value)
    if kind is float:
        if not math.isfinite(value):
            raise _Uncacheable()
        return ("float", value.hex())
    if kind in (list, tuple):
        return (kind.__name__, tuple(_freeze(item) for item in value))
    if kind is dict and all(type(key) is str for key in value):
        return ("dict", tuple((key, _freeze(value[key])) for key in sorted(value)))
    raise _Uncacheable()


def query_expert_memo(function, *, max_entries=512):
    """Return a fresh bounded cache for one query's expert-evidence calls."""
    if type(max_entries) is not int or max_entries < 0:
        raise ValueError("max_entries must be a nonnegative integer")
    entries = OrderedDict()
    counts = {"hits": 0, "misses": 0, "bypassed": 0, "evictions": 0}

    @wraps(function)
    def call(*args, **kwargs):
        try:
            key = (_freeze(args), _freeze(kwargs))
        except (TypeError, ValueError, RecursionError):
            counts["bypassed"] += 1
            return function(*args, **kwargs)
        if key in entries:
            counts["hits"] += 1
            entries.move_to_end(key)
            return deepcopy(entries[key])
        counts["misses"] += 1
        result = function(*args, **kwargs)
        if isinstance(result, dict) and max_entries:
            entries[key] = deepcopy(result)
            if len(entries) > max_entries:
                entries.popitem(last=False)
                counts["evictions"] += 1
        return result

    call.cache_info = lambda: {**counts, "entries": len(entries), "max_entries": max_entries}
    return call
