"""Evidence round trips and score/order invariants; synthetic inputs only."""
import copy
import importlib
import json
from pathlib import Path

import numpy as np
import pytest

from matching.scores import cosine_score
from matching.spectral_evidence import capture_spectral_evidence, attach_spectral_evidence


def candidate(q, l, *, qpmz=500.0, lpmz=500.0, reference="reference-1"):
    return {"_score_query_peaks": q, "_score_library_peaks": l,
            "cosine": cosine_score(q, l, query_pmz=qpmz, lib_pmz=lpmz),
            "final_score": 0.7, "rank_score": 0.71,
            "lib_spectrum_id": reference, "lib_search_precursor_mz": lpmz,
            "dataset_source": "synthetic", "lib_source_accession": "test-accession",
            "inchikey": "AAAAAAAAAAAAAA-BBBBBBBBBB-C", "smiles": "CCO"}


@pytest.mark.parametrize("q,l,qpmz,lpmz", [
    ([[100.001, 10], [100.004, 20]], [[100.001, 30], [100.004, 40]], 500, 500),
    ([[120, 20], [100, 10]], [[110, 10], [130, 20]], 500, 510),
    ([[100, 1], [100, 2]], [[100, 2], [100, 1]], 500, 500),
    ([[100, 0]], [[100, 1]], 500, 500),
])
def test_pairs_replay_same_score_and_point_to_unique_saved_indices(q, l, qpmz, lpmz):
    source = candidate(q, l, qpmz=qpmz, lpmz=lpmz)
    before = copy.deepcopy(source)
    evidence = capture_spectral_evidence(source, query_precursor_mz=qpmz)
    assert source == before
    assert evidence["replayed_modified_cosine"] == pytest.approx(source["cosine"], abs=1e-12)
    assert evidence["ranking_influence"] is False
    assert evidence["fragment_atom_assignment"] is False
    qids, lids = [], []
    for pair in evidence["matched_pairs"]:
        qi, li = pair["query_peak_index"], pair["library_peak_index"]
        assert evidence["query_peaks"][qi][0] == pair["query_mz"]
        assert evidence["library_peaks"][li][0] == pair["lib_mz"]
        assert abs(pair["residual_da"]) <= evidence["tolerance_da"]
        qids.append(qi); lids.append(li)
    assert len(set(qids)) == len(qids) and len(set(lids)) == len(lids)
    assert len(evidence["view_sha256"]) == 64
    assert json.loads(json.dumps(evidence, allow_nan=False)) == evidence


def test_randomized_replay_does_not_modify_scores_or_input_order():
    rng = np.random.default_rng(744)
    for _ in range(32):
        mz = rng.uniform(60, 400, 10)
        q = np.column_stack([mz, rng.uniform(1, 100, 10)]).tolist()
        shift = float(rng.choice([0, 8, -12]))
        l = np.column_stack([mz - shift + rng.uniform(-0.05, 0.05, 10), rng.uniform(1, 100, 10)]).tolist()
        source = candidate(q, l, qpmz=500, lpmz=500-shift)
        evidence = capture_spectral_evidence(source, query_precursor_mz=500)
        assert evidence["replayed_modified_cosine"] == pytest.approx(source["cosine"], abs=1e-12)
        assert evidence["replay_matches_comparison"]
        assert source["_score_query_peaks"] == q


def test_adjusted_score_is_kept_separate_from_pair_replay():
    source = candidate([[100, 1]], [[100, 1]])
    source["cosine"] = 0.3
    e = capture_spectral_evidence(source, query_precursor_mz=500)
    assert e["comparison_cosine"] == 0.3
    assert e["replayed_modified_cosine"] == 1.0
    assert e["replay_matches_comparison"] is False
    assert source["cosine"] == 0.3


def test_same_compound_different_reference_does_not_get_best_reference_peaks():
    best = candidate([[100, 1]], [[100, 1]], reference="ref-best")
    other = candidate([[100, 1]], [[101, 1]], reference="ref-other")
    public = {k: copy.deepcopy(v) for k,v in other.items() if not k.startswith("_")}
    best["top_candidates"] = [public]
    attach_spectral_evidence(best, [best, other], query_precursor_mz=500)
    assert best["spectral_evidence"]["reference"]["spectrum_id"] == "ref-best"
    assert public["spectral_evidence"]["reference"]["spectrum_id"] == "ref-other"
    assert public["spectral_evidence"]["library_peaks"] == [[101, 1]]


def test_missing_view_remains_unavailable():
    e = capture_spectral_evidence({}, query_precursor_mz=500)
    assert e["status"] == "unavailable"
    assert "matched_pairs" not in e
    with pytest.raises(ValueError):
        capture_spectral_evidence(candidate([[100, float("nan")]], [[100,1]]), query_precursor_mz=500)


def test_search_and_serialization_preserve_full_candidate_scores(monkeypatch):
    spec = importlib.util.spec_from_file_location(
        "reference_fixture", Path(__file__).with_name("test_v15_search_reference_evidence.py"),
    )
    fixtures = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(fixtures)
    _disable_optional_layers, _record = fixtures._disable_optional_layers, fixtures._record
    QUERY, PRECURSOR = fixtures.QUERY, fixtures.PRECURSOR
    engine = importlib.import_module("matching.search")
    _disable_optional_layers(monkeypatch)
    records = [_record(0, "GNPS", QUERY, 20.0), _record(1, "MassBank", QUERY, 20.0)]
    bucket = int(PRECURSOR / engine.MZ_BUCKET_DALTON)
    kwargs = dict(records=records, mz_index={bucket:[0,1]}, query_polarity="POS",
                  cos_thresh=0, final_thresh=0, strict_evidence_gate=False,
                  enable_chimeric_deconvolution=False)
    baseline = engine.search(QUERY, PRECURSOR, capture_peak_evidence=False, **kwargs)
    captured = engine.search(QUERY, PRECURSOR, capture_peak_evidence=True, **kwargs)
    def strip_evidence(value):
        if isinstance(value, dict):
            return {k:strip_evidence(v) for k,v in value.items() if k != "spectral_evidence"}
        if isinstance(value, list): return [strip_evidence(v) for v in value]
        return value
    assert strip_evidence(captured) == baseline
    e = captured["spectral_evidence"]
    assert e["status"] == "captured"
    assert e["reference"]["spectrum_id"] == captured["lib_spectrum_id"]
    from web.server import _serialize_results
    row = {"spectrum_id":"test-query", "precursor_mz":PRECURSOR,
           "query_peaks":QUERY, "raw_peaks":QUERY, "match":captured}
    saved = _serialize_results([row])[0]
    assert saved["match"]["spectral_evidence"] == e
    assert saved["match"]["top_candidates"][0]["spectral_evidence"]["status"] == "captured"
    assert saved["query_peaks"] == QUERY


@pytest.mark.parametrize("failure", ["import", "capture"])
def test_optional_capture_failure_preserves_original_search(monkeypatch, failure):
    import builtins
    spec = importlib.util.spec_from_file_location(
        "failure_reference_fixture", Path(__file__).with_name("test_v15_search_reference_evidence.py"),
    )
    fixtures = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(fixtures)
    fixtures._disable_optional_layers(monkeypatch)
    engine = importlib.import_module("matching.search")
    records = [fixtures._record(0, "GNPS", fixtures.QUERY, 20.0),
               fixtures._record(1, "MassBank", fixtures.QUERY, 20.0)]
    bucket = int(fixtures.PRECURSOR / engine.MZ_BUCKET_DALTON)
    kwargs = dict(records=records, mz_index={bucket: [0, 1]}, query_polarity="POS",
                  cos_thresh=0, final_thresh=0, strict_evidence_gate=False,
                  enable_chimeric_deconvolution=False)
    baseline = engine.search(fixtures.QUERY, fixtures.PRECURSOR, capture_peak_evidence=False, **kwargs)
    if failure == "import":
        original_import = builtins.__import__
        def fail_import(name, *args, **kwargs):
            if name == "matching.spectral_evidence":
                raise ImportError("injected optional observer import failure")
            return original_import(name, *args, **kwargs)
        monkeypatch.setattr(builtins, "__import__", fail_import)
    else:
        observer = importlib.import_module("matching.spectral_evidence")
        def fail_capture(best, *args, **kwargs):
            best["spectral_evidence"] = {"status": "partial"}
            raise RuntimeError("injected optional observer failure")
        monkeypatch.setattr(observer, "attach_spectral_evidence", fail_capture)
    observed = engine.search(fixtures.QUERY, fixtures.PRECURSOR, capture_peak_evidence=True, **kwargs)
    def strip(value):
        if isinstance(value, dict):
            return {k: strip(v) for k, v in value.items() if k != "spectral_evidence"}
        if isinstance(value, list):
            return [strip(v) for v in value]
        return value
    assert strip(observed) == baseline
    for entry in [observed, *observed["top_candidates"]]:
        assert entry["spectral_evidence"]["status"] == "unavailable"
        assert entry["spectral_evidence"]["reason"] == "scoring_view_capture_failed"
        assert entry["spectral_evidence"]["error_type"] == (
            "ImportError" if failure == "import" else "RuntimeError")
