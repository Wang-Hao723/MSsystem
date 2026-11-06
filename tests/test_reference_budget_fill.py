"""Dense reference windows retain additional spectra within the original budget."""
import importlib
import random
import pytest

search = importlib.import_module('matching.search')


@pytest.fixture
def controlled_records(monkeypatch):
    monkeypatch.setattr(search, 'compound_identity_key', lambda r: r['cid'])
    monkeypatch.setattr(search, '_record_index_precursor_mz', lambda r: r['precursor_mz'])
    monkeypatch.setattr(search, '_library_peaks_for_config', lambda r, c: r['peaks'])
    monkeypatch.setattr(search, 'batch_modified_cosine_scores',
                        lambda q, spectra, pmz, pmzs: [p[0][1] for p in spectra])
    monkeypatch.setattr(search, '_structure_adduct_precursor_audit',
                        lambda r, query_polarity=None: {'hard_reject': r.get('invalid', False)})
    return [{'cid': name, 'precursor_mz': 100.0, 'neutral_mass': 99.0,
             'source_tier': 'experimental', 'spectrum_id': f'{name}{i}',
             'peaks': [[10.0, value]]}
            for name, values in [('A', [.99, .98, .97, .96]), ('B', [.85, .84]), ('C', [.75])]
            for i, value in enumerate(values)]


def run(records, indices=None, limit=5):
    return search._prescreen_dense_candidates(
        list(range(len(records))) if indices is None else indices,
        records, [[10.0, 1.0]], 100.0, limit, query_polarity='POS')


def test_unused_budget_is_filled_after_all_compounds_receive_one_reference(controlled_records):
    records = controlled_records
    selected, scores = run(records)
    assert [records[i]['spectrum_id'] for i in selected] == ['A0', 'B0', 'C0', 'A1', 'B1']
    assert len(selected) == len(set(selected)) == 5
    assert set(scores) == set(selected)


def test_record_order_does_not_change_selected_references(controlled_records):
    expected = run(controlled_records)
    for seed in range(20):
        indices = list(range(len(controlled_records)))
        random.Random(seed).shuffle(indices)
        assert run(controlled_records, indices) == expected


def test_invalid_best_reference_cannot_displace_valid_same_compound_reference(controlled_records, monkeypatch):
    controlled_records[0]['invalid'] = True
    calls = []
    def audit(record, query_polarity=None):
        calls.append(query_polarity)
        return {'hard_reject': record.get('invalid', False)}
    monkeypatch.setattr(search, '_structure_adduct_precursor_audit', audit)
    selected, scores = run(controlled_records)
    assert 0 not in selected and len(selected) == 5
    assert controlled_records[selected[0]]['spectrum_id'] == 'A1'
    assert calls == ['POS'] * len(controlled_records)
    assert set(scores) == set(selected)


def test_all_eligible_references_retained_below_budget(controlled_records):
    selected, _ = run(controlled_records, limit=500)
    assert set(selected) == set(range(len(controlled_records)))


def test_all_invalid_references_return_empty_result(controlled_records):
    for record in controlled_records:
        record['invalid'] = True
    assert run(controlled_records) == ([], {})
