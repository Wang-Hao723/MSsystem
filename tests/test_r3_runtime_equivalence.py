import datetime
import io
import json
import math

import pytest

from pipeline.query_expert_memo import query_expert_memo
from runtime.json_encoding import write_json_encoded


@pytest.mark.parametrize(
    ("value", "kwargs"),
    [
        ({"name": "白术内酯II", "escaped": "line\n\"\t", "x": -0.0}, {"ensure_ascii": False, "separators": (",", ":")}),
        ({"when": datetime.date(2026, 10, 24), "nested": [{1: True}, None]}, {"ensure_ascii": False, "separators": (",", ":"), "default": str}),
        ({7: "number", "z": ["α", 1e-12, 1e20]}, {"ensure_ascii": True, "indent": 2, "sort_keys": False}),
        ({"x": float("nan"), "y": float("inf")}, {}),
        ({"x": float("nan")}, {"allow_nan": False}),
        ({(1, 2): 3}, {}),
        ({(1, 2): 3, "ok": 4}, {"skipkeys": True}),
        (None, {}),
    ],
)
def test_single_write_encoder_matches_json_dump_value_or_exception(value, kwargs):
    if value is None:
        value = []
        value.append(value)
    outcomes = []
    for writer in (json.dump, write_json_encoded):
        handle = io.StringIO()
        try:
            writer(value, handle, **kwargs)
            outcomes.append(("ok", handle.getvalue()))
        except Exception as exc:
            outcomes.append((type(exc), str(exc)))
    assert outcomes[0] == outcomes[1]


def _expert_fixture():
    calls = []

    def evidence(*args, **kwargs):
        calls.append((args, kwargs))
        return {"score": 0.3, "annotation": {"tokens": ["original"]}}

    return calls, evidence


def test_expert_memo_isolates_caller_mutation():
    calls, evidence = _expert_fixture()
    cached = query_expert_memo(evidence)
    kwargs = {"query_peaks": [[55.0177, 1.0]], "candidate_smiles": "C=C"}
    first = cached(**kwargs)
    first["annotation"]["tokens"].append("mutation")
    assert cached(**kwargs)["annotation"]["tokens"] == ["original"]
    assert len(calls) == 1


def test_expert_memo_distinguishes_peak_and_collision_energy_inputs():
    calls, evidence = _expert_fixture()
    cached = query_expert_memo(evidence)
    cached(query_peaks=[[55.0177, 1.0]], query_collision_energy={"value": 40.0, "unit": "eV"})
    cached(query_peaks=[[55.0177, 1.0]], query_collision_energy={"value": 45.0, "unit": "eV"})
    cached(query_peaks=[[55.0177, 0.5]], query_collision_energy={"value": 40.0, "unit": "eV"})
    assert len(calls) == 3


def test_expert_memo_has_bounded_lru_eviction():
    calls, evidence = _expert_fixture()
    cached = query_expert_memo(evidence, max_entries=2)
    for value in (1.0, 2.0, 3.0, 1.0):
        cached(query_peaks=[[55.0, value]])
    assert len(calls) == 4
    assert cached.cache_info()["evictions"] == 2


def test_expert_memo_is_fresh_per_query():
    calls, evidence = _expert_fixture()
    query_expert_memo(evidence)(candidate_smiles="C=C")
    query_expert_memo(evidence)(candidate_smiles="C=C")
    assert len(calls) == 2


def test_expert_memo_bypasses_nonfinite_values():
    calls, evidence = _expert_fixture()
    cached = query_expert_memo(evidence)
    cached(score=float("nan"))
    cached(score=float("nan"))
    assert len(calls) == 2
    assert cached.cache_info()["bypassed"] == 2


def test_expert_memo_bypasses_unsupported_types():
    calls, evidence = _expert_fixture()
    cached = query_expert_memo(evidence)
    unsupported = object()
    cached(candidate_classification=unsupported)
    cached(candidate_classification=unsupported)
    assert len(calls) == 2


def test_expert_memo_never_caches_exceptions():
    calls = []

    def failing(**kwargs):
        calls.append(kwargs)
        raise RuntimeError("unchanged failure")

    cached = query_expert_memo(failing)
    for _ in range(2):
        with pytest.raises(RuntimeError, match="unchanged failure"):
            cached(candidate_smiles="C=C")
    assert len(calls) == 2


def test_query_search_installs_one_expert_cache_per_query(monkeypatch):
    from pipeline import query_search

    calls = []

    def evidence(**kwargs):
        calls.append(kwargs)
        return {"score": 0.4}

    namespace = {"expert_evaluate": evidence}
    exec(
        "def synthetic_search(*args, **kwargs):\n"
        "    return {'first': expert_evaluate(candidate_smiles='C=C'), "
        "'second': expert_evaluate(candidate_smiles='C=C')}\n",
        namespace,
    )
    monkeypatch.setattr(query_search._base, "search", namespace["synthetic_search"])
    monkeypatch.setattr(query_search._base, "expert_evaluate", evidence)
    monkeypatch.setattr(
        query_search,
        "guard_search_arguments",
        lambda function, args, kwargs: (args, kwargs, {"status": "synthetic"}),
    )

    first = query_search.search()
    second = query_search.search()
    assert first["first"] == first["second"] == {"score": 0.4}
    assert second["first"] == second["second"] == {"score": 0.4}
    assert len(calls) == 2
