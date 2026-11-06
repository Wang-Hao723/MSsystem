from pipeline.query_record_memo import QueryRecordMemo


class CountingStore:
    def __init__(self, rows):
        self.rows = rows
        self.cached_reads = 0
        self.uncached_reads = 0

    def __len__(self):
        return len(self.rows)

    def __getitem__(self, index):
        self.cached_reads += 1
        return self.rows[index]

    def read_uncached(self, index):
        self.uncached_reads += 1
        return self.rows[index]

    def cache_info(self):
        return {"cached_reads": self.cached_reads,
                "uncached_reads": self.uncached_reads}


def test_query_memo_reuses_same_record_object_without_copying():
    rows = [{"value": 1}, {"value": 2}]
    store = CountingStore(rows)
    memo = QueryRecordMemo(store, max_entries=2)

    assert memo[0] is rows[0]
    assert memo[0] is rows[0]
    assert store.cached_reads == 1
    assert store.uncached_reads == 0
    assert memo.memo_info() == {
        "entries": 1, "max_entries": 2, "hits": 1, "misses": 1,
        "evictions": 0,
    }


def test_query_memo_is_bounded_and_delegates_evicted_reads():
    store = CountingStore([{"value": i} for i in range(3)])
    memo = QueryRecordMemo(store, max_entries=2)

    assert [memo[i]["value"] for i in (0, 1, 2, 0)] == [0, 1, 2, 0]
    assert store.cached_reads == 4
    assert store.uncached_reads == 0
    assert memo.memo_info()["entries"] == 2
    assert memo.memo_info()["evictions"] == 2


def test_query_memo_falls_back_for_generic_sequences():
    rows = [{"value": 1}]
    memo = QueryRecordMemo(rows)

    assert memo[0] is rows[0]
    assert memo[0] is rows[0]
    assert memo.memo_info()["hits"] == 1


def test_each_query_delegates_first_read_to_shared_process_store():
    store = CountingStore([{"value": 1}])

    assert QueryRecordMemo(store)[0] is store.rows[0]
    assert QueryRecordMemo(store)[0] is store.rows[0]
    assert store.cached_reads == 2
    assert store.uncached_reads == 0

