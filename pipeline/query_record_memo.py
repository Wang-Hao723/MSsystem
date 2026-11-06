"""Query-local decoded-record reuse for immutable mmap library records."""

from __future__ import annotations

import operator
from collections import OrderedDict
from collections.abc import Sequence


DEFAULT_QUERY_MEMO_ENTRIES = 4096


class QueryRecordMemo(Sequence):
    """Keep bounded strong references to records revisited by one search.

    ``MMapRecordStore`` already has a process-wide byte-bounded LRU. Dense
    precursor windows can cycle through more decoded records than that cache
    retains before the same records are revisited by polarity, prescreen and
    scoring passes. This view changes no record and delegates every first read
    to the accepted store; it only keeps a bounded per-query reference cache.
    """

    def __init__(self, records, *, max_entries: int = DEFAULT_QUERY_MEMO_ENTRIES):
        if isinstance(max_entries, bool) or int(max_entries) < 1:
            raise ValueError("max_entries must be a positive integer")
        self.records = records
        self.max_entries = int(max_entries)
        self._cache = OrderedDict()
        self.hits = 0
        self.misses = 0
        self.evictions = 0

    def __len__(self):
        return len(self.records)

    def __getitem__(self, index):
        if isinstance(index, slice):
            return [self[i] for i in range(*index.indices(len(self)))]
        index = operator.index(index)
        if index < 0:
            index += len(self)
        if not 0 <= index < len(self):
            raise IndexError(index)
        if index in self._cache:
            self.hits += 1
            self._cache.move_to_end(index)
            return self._cache[index]
        # Delegate each first query-local read to the accepted process store.
        # This preserves its cross-query hot cache; repeated accesses within
        # this query then reuse the same immutable object without another LRU
        # lookup, decode, or deep-size traversal after process-cache eviction.
        record = self.records[index]
        self.misses += 1
        self._cache[index] = record
        if len(self._cache) > self.max_entries:
            self._cache.popitem(last=False)
            self.evictions += 1
        return record

    def __getattr__(self, name):
        return getattr(self.records, name)

    def memo_info(self) -> dict:
        return {
            "entries": len(self._cache),
            "max_entries": self.max_entries,
            "hits": self.hits,
            "misses": self.misses,
            "evictions": self.evictions,
        }

