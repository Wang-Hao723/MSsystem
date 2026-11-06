import builtins
import sys
import types

import numpy as np
import pytest

from classifier.msclass_wrapper import MSClassWrapper
from pipeline import structure_prediction as structure


def test_classification_is_invariant_to_peak_order_and_filters_before_budget():
    peaks = [(i + 10.2, i + 1.) for i in range(250)]
    invalid = [(np.nan, 1), (10, np.inf), (np.inf, 1), (3, -1)] * 70
    expected = MSClassWrapper._bin_spectrum(peaks)
    np.testing.assert_array_equal(expected, MSClassWrapper._bin_spectrum(np.array((invalid + peaks)[::-1])))
    assert expected[10] == 0 and expected[259] == 1


def test_classification_finite_after_extreme_intensities():
    result = MSClassWrapper._bin_spectrum([(100.1, 1e308), (100.2, 1e308), (200, 1e308)])
    assert np.isfinite(result).all()
    assert result[100] == 1 and result[200] == .5


def test_classification_requires_three_usable_peaks():
    wrapper = object.__new__(MSClassWrapper)
    assert wrapper.predict(np.array([[10, 1], [20, 2], [30, np.nan]])) is None


def test_missing_classification_assets_checked_before_torch(tmp_path, monkeypatch):
    original = builtins.__import__
    def guarded(name, *args, **kwargs):
        if name.startswith('torch'):
            raise AssertionError('Must check assets before importing torch')
        return original(name, *args, **kwargs)
    monkeypatch.setattr(builtins, '__import__', guarded)
    with pytest.raises(FileNotFoundError):
        MSClassWrapper(tmp_path, tmp_path / 'missing.pt')


def test_missing_structure_checkpoint_needs_no_model_imports(tmp_path, monkeypatch):
    monkeypatch.setattr(structure, '_model', None)
    monkeypatch.setattr(structure, '_ensure_imports', lambda: pytest.fail('Import before asset check'))
    assert structure.load_model(tmp_path / 'absent.pt') == (None, None, None, None)


@pytest.mark.parametrize('value', [0, -1, np.nan, np.inf, True])
def test_invalid_generation_temperature_rejected(value):
    with pytest.raises(ValueError):
        structure.predict_from_spectrum([], temperature=value)


def test_structure_failed_load_does_not_publish_partial_model(tmp_path, monkeypatch):
    for name in ('model.py', 'selfies_utils.py', 'utils.py', 'checkpoint.pt'):
        (tmp_path / name).write_text('')
    for name in ('_model', '_vocab', '_config', '_device', '_loaded_path'):
        monkeypatch.setattr(structure, name, None)
    monkeypatch.setattr(structure, '_MODEL_ROOT', tmp_path)
    monkeypatch.setattr(structure, '_ensure_imports', lambda: None)
    monkeypatch.setattr(structure, 'get_device', lambda: 'cpu')
    config = {'model': {'decoder': {'d_model': 64}}}
    checkpoint = {'config': config, 'vocab': {'token_to_idx': {'s': 0}, 'idx_to_token': {'0': 's'}}, 'model_state_dict': {}}
    fake_torch = types.SimpleNamespace(load=lambda *a, **k: checkpoint)
    monkeypatch.setitem(sys.modules, 'torch', fake_torch)
    class Vocab:
        def __len__(self): return len(self.token_to_idx)
    instances = []
    class Model:
        def __init__(self, **kwargs):
            assert kwargs['d_model'] == 64
            instances.append(self)
        def to(self, device): return self
        def state_dict(self): return {'encoder.weight': 1}
        def load_state_dict(self, state, strict):
            assert strict and state == {'encoder.weight': 1}
        def eval(self): pass
    monkeypatch.setitem(sys.modules, 'selfies_utils', types.SimpleNamespace(SelfiesVocabulary=Vocab))
    monkeypatch.setitem(sys.modules, 'model', types.SimpleNamespace(MS2SELFIESModel=Model))
    with pytest.raises(ValueError, match='missing required weights'):
        structure.load_model(tmp_path / 'checkpoint.pt')
    assert structure._model is None and structure._vocab is None and structure._config is None
    checkpoint['model_state_dict'] = {'encoder.weight': 1}
    loaded = structure.load_model(tmp_path / 'checkpoint.pt')
    assert loaded[0] is instances[-1] and len(instances) == 2
    assert structure.load_model(tmp_path / 'checkpoint.pt')[0] is loaded[0]
