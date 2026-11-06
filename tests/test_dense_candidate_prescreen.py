from matching.search import _prescreen_dense_candidates


def _record(idx, precursor_mz, peaks, inchikey=""):
    return {
        "precursor_mz": precursor_mz,
        "peaks": peaks,
        "clean_peaks": peaks,
        "spectrum_id": f"spec-{idx}",
        "inchikey": inchikey,
        "smiles": "",
    }


def test_dense_prescreen_keeps_spectral_match_that_is_not_mass_nearest():
    query_pmz = 500.0
    query = [[100.0, 1.0], [200.0, 0.8], [300.0, 0.6]]
    records = [
        _record(i, query_pmz + i * 0.0001, [[50.0, 1.0], [60.0, 0.8], [70.0, 0.6]])
        for i in range(501)
    ]
    # This spectrum would be lost by the old nearest-mass-first 500-candidate cap.
    records.append(_record(501, 500.4, query))

    selected, cosine_cache = _prescreen_dense_candidates(
        list(range(len(records))), records, query, query_pmz, limit=500,
    )

    assert 501 in selected
    assert cosine_cache[501] > 0.99


def test_dense_prescreen_uses_best_spectrum_for_same_compound():
    query_pmz = 500.0
    query = [[100.0, 1.0], [200.0, 0.8], [300.0, 0.6]]
    ik = "ABCDEFGHIJKLMN-AAAAAA-A"
    records = [
        _record(0, 500.0, [[50.0, 1.0], [60.0, 0.8], [70.0, 0.6]], ik),
        _record(1, 500.2, query, ik),
    ]

    selected, _ = _prescreen_dense_candidates(
        [0, 1], records, query, query_pmz, limit=1,
    )

    assert selected == [1]


def test_dense_prescreen_keeps_precursor_shifted_fragment_match():
    query_pmz = 510.0
    query = [[110.0, 1.0], [210.0, 0.8], [310.0, 0.6]]
    records = [
        _record(i, query_pmz + i * 0.0001, [[50.0, 1.0], [60.0, 0.8], [70.0, 0.6]])
        for i in range(501)
    ]
    # Raw fragments differ by +10 Da, exactly matching the precursor delta.
    records.append(_record(501, 500.0, [[100.0, 1.0], [200.0, 0.8], [300.0, 0.6]]))

    selected, cosine_cache = _prescreen_dense_candidates(
        list(range(len(records))), records, query, query_pmz, limit=500,
    )

    assert 501 in selected
    assert cosine_cache[501] > 0.99
