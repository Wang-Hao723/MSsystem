"""Physical input distinctions and candidate independence for the new experiments."""
from pathlib import Path
import sqlite3
import sys
import numpy as np
import pytest
import torch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'eval/msystem_method_breakthrough_20260908'))
from bottleneck_training_inventory import observed_energy, peak_digest
from bottleneck_catalog_recall import enumerate_candidates, SHIFTS
from fit_bottleneck_local_pairs import LocalPairModel, graph_features, peak_features, pair_choice
from bottleneck_catalog_ranking import canonical_max, new_representations, top_key


@pytest.mark.parametrize('raw', ['30', 'Ramp 5-60 V', ['10 eV', '20 eV'], None, '10-40 eV'])
def test_ambiguous_or_multiple_energies_do_not_become_single_observations(raw):
    assert observed_energy(raw)[0] is None


def test_explicit_units_remain_distinct():
    assert observed_energy('20 eV')[0] == (20., 'eV')
    assert observed_energy('20 V')[0] == (20., 'V')
    assert observed_energy('20 NCE')[0] == (20., 'NCE')
    assert observed_energy(['20% (nominal)'])[0] == (20., '%')


@pytest.mark.parametrize('precursor', [None, 0, -1, float('nan')])
def test_missing_precursor_is_not_reconstructed_from_a_candidate(precursor):
    with pytest.raises(ValueError, match='measured positive precursor'):
        peak_features([[50., 1.]], precursor)


def test_float_noise_cannot_rank_an_identical_composition_control():
    assert pair_choice([.2, .2 + 2.98e-8]) == 0
    assert pair_choice([.2 + 2.98e-8, .2]) == 0
    assert pair_choice([.2, .21]) == 1


def test_spectral_deduplication_ignores_order_and_overall_intensity_scale():
    first = peak_digest([[50., 2.], [60., 1.]])
    assert first == peak_digest([[60., 10.], [50., 20.]])
    assert first != peak_digest([[50., 2.], [60., 2.]])


def test_candidate_admission_ignores_answer_and_preserves_old_candidates():
    db = sqlite3.connect(':memory:')
    db.execute('CREATE TABLE structures(raw14 TEXT, atomic_mass REAL, charge INTEGER)')
    mass = 100.
    db.executemany('INSERT INTO structures VALUES(?,?,?)', [('admitted', mass, 0), ('too_far', mass + .1, 0)])
    base = {'query_id': 'observed', 'domain': 'test', 'precursor_mz': mass + SHIFTS['[M+H]+'],
            'adduct': '[M+H]+', 'original_keys': ['old_retained']}
    output = list(enumerate_candidates(db, [base]))[0]
    poisoned = list(enumerate_candidates(db, [dict(base, target_raw14='too_far')]))[0]
    assert output == poisoned
    assert output['c1_keys'] == ['admitted', 'old_retained']
    assert output['added_keys'] == ['admitted']


def test_pool_expansion_keeps_original_representations_and_checks_each_new_mass():
    db = sqlite3.connect(':memory:')
    db.execute('CREATE TABLE structures(raw14 TEXT, canonical TEXT, atomic_mass REAL, charge INTEGER)')
    db.executemany('INSERT INTO structures VALUES(?,?,?,?)', [
        ('old', 'alternative old representation', 100., 0),
        ('new', 'eligible new representation', 100., 0),
        ('new', 'wrong mass representation', 100.1, 0)])
    observed = {'precursor_mz': 100. + SHIFTS['[M+H]+'], 'adduct': '[M+H]+'}
    assert new_representations(db, observed, ['old']) == [('new', 'eligible new representation')]


def test_paired_score_uses_same_canonical_before_raw14_max_and_retains_missing():
    values = canonical_max(np.asarray([[.9, .8], [.7, .1], [.4, .3]]), [0, 0, 1], 3)
    assert values[0] == pytest.approx([.9, .8, .6])
    assert values[1] == pytest.approx([.4, .3, .1])
    assert np.isnan(values[2]).all()
    assert top_key(['z', 'a', 'missing'], values[:, 0]) == 'z'
    assert top_key(['z', 'a'], np.asarray([.5, .5])) == 'a'
    assert top_key(['missing'], np.asarray([float('nan')])) is None


def test_graph_model_atom_order_invariant_and_composition_control_collapses_isomers():
    torch.set_num_threads(1)
    torch.manual_seed(12)
    model = LocalPairModel().eval()
    peak = torch.randn(1, 4, 12)
    weights = torch.full((1, 4), .25)
    condition = torch.zeros(1, 6)

    def score(smiles, composition=False, order=None):
        atom, adjacency = graph_features(smiles, composition)
        if order is not None:
            atom, adjacency = atom[order], adjacency[np.ix_(order, order)]
        with torch.no_grad():
            return model(torch.tensor(atom)[None], torch.tensor(adjacency)[None],
                         torch.ones(1, len(atom)), peak, weights, condition)

    assert torch.allclose(score('CCCO'), score('CCCO', order=[3, 1, 0, 2]), atol=1e-6)
    assert torch.allclose(score('CCCO', True), score('CC(O)C', True), atol=1e-6)
    assert not torch.allclose(score('CCCO'), score('CC(O)C'), atol=1e-8)
