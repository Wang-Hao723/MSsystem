"""Portable output-projection test; actual database windows have separate QA."""
import copy
from matching.evirank_retrieval import SQLiteCandidateRecallRetriever, load_evirecall_policy


def test_sql_reference_outputs_keep_stage_and_original_primary_key_with_duplicate_accessions():
    retriever = object.__new__(SQLiteCandidateRecallRetriever)
    retriever.policy = load_evirecall_policy()
    retriever.library_snapshot_id = 'synthetic-fixed-database-snapshot'
    tier = retriever.policy['library_gate']['searchable_source_tiers'][0]
    peaks = [[15.02, 10], [17.02, 20], [18.03, 40], [19.01, 30],
             [20.04, 60], [21.02, 10], [22.02, 40], [24.02, 20],
             [25.02, 35], [27.02, 70], [29.04, 90], [31.02, 100]]
    rows = []
    for pk, identity, smiles, stage, source in [
        (101, 'LFQSCWFLJHTTHZ', 'CCO', 'MS2', 'gnps'),
        (102, 'LCGLNKUTAGEVQW', 'COC', 'MS2-MS3 Composite', 'massbank'),
    ]:
        rows.append({'database_spectrum_id': pk, 'spectrum_id': 'duplicate-accession',
                     'connectivity_identity': identity, 'identity_source': 'synthetic_verified_structure',
                     'compound_name': smiles, 'molecular_formula': 'C2H6O', 'smiles': smiles,
                     'inchikey': '', 'precursor_mz': 47.049, 'adduct': '[M+H]+',
                     'charge': 1, 'ion_mode': 'POSITIVE', 'library_polarity': 'POS',
                     'source_tier': tier, 'dataset_source': source, 'quality_flags': [],
                     'peaks': copy.deepcopy(peaks), 'processed_views': {},
                     'instrument': 'Thermo Q Exactive HF', 'instrument_type': 'Orbitrap',
                     'source_ms_level': stage, 'collision_energy': None})
    # Only the database window is replaced. The real spectral-view scoring and
    # the two reference output projections execute unchanged.
    retriever._fetch_window = lambda precursor, polarity: (rows, {'status': 'available'})
    output = retriever.retrieve_scan({'query_id': 'portable-stage-output',
                                     'precursor_mz': 47.049, 'polarity': 'POS', 'peaks': peaks})
    assert set(output['candidates']) == {row['connectivity_identity'] for row in rows}
    observed_primary_keys = set(); observations = 0
    for identity, candidate in output['candidates'].items():
        expected = next(row for row in rows if row['connectivity_identity'] == identity)
        references = list(candidate['reference_records'].values())
        references += [record for strata in candidate['source_reference_records'].values()
                       for record in strata.values()]
        assert references
        for reference in references:
            assert reference['spectrum_id'] == 'duplicate-accession'
            assert reference['database_spectrum_id'] == expected['database_spectrum_id']
            assert reference['source_ms_level'] == expected['source_ms_level']
            observed_primary_keys.add(reference['database_spectrum_id']); observations += 1
    assert observed_primary_keys == {101, 102} and observations >= 4
