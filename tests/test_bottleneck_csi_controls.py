"""Supplemental CSI controls cannot use incomplete query sets or stale gates."""
import importlib.util
import json
from pathlib import Path
import sys
import numpy as np
import pytest

CODE = Path(__file__).resolve().parents[1] / 'eval/msystem_method_breakthrough_20260908'
sys.path.insert(0, str(CODE))
import bottleneck_csi_controls as controls


def test_paired_macro_counts_each_structure_equally():
    truth = {'a1': 'A', 'a2': 'A', 'b': 'B'}
    full = {'a1': 'A', 'a2': 'A', 'b': None}
    baseline = {'a1': None, 'a2': None, 'b': 'B'}
    assert controls.metric(full, truth)['structure_macro'] == .5
    assert controls.paired_effect(full, baseline, truth)['gain'] == 0.
    with pytest.raises(RuntimeError, match='denominator'):
        controls.paired_effect(full, {'a1': None}, truth)


def test_duplicate_or_out_of_pool_predictions_do_not_pass(tmp_path):
    path = tmp_path / 'predictions.jsonl'
    held = {'query': {'keys': ['allowed']}}
    row = json.dumps({'query_id': 'query', 'key': 'allowed'}) + '\n'
    path.write_text(row + row)
    with pytest.raises(RuntimeError, match='Duplicate'):
        controls.predictions(path, held)
    path.write_text(json.dumps({'query_id': 'query', 'key': 'inserted answer'}) + '\n')
    with pytest.raises(RuntimeError, match='outside'):
        controls.predictions(path, held)


def test_stale_successful_extra_gate_is_rejected(tmp_path, monkeypatch):
    audit, fit = tmp_path / 'audit', tmp_path / 'fit'
    (audit / 'CSI_ADDITIONAL_CONTROLS/inner').mkdir(parents=True)
    (fit / 'inner').mkdir(parents=True)
    protocol = audit / 'protocol.json'
    protocol.write_text('{}')
    pred, decision = fit / 'inner/full_predictions.jsonl', fit / 'inner/DECISION.json'
    pred.write_text('original predictions')
    decision.write_text('{"passed":true}')
    report = {'passed': True, 'original_csi_gate_passed': True, 'additional_protocol_sha256': controls.sha(protocol),
              'csi_prediction_sha256': controls.sha(pred), 'original_csi_decision_sha256': controls.sha(decision)}
    (audit / 'CSI_ADDITIONAL_CONTROLS/inner/REPORT.json').write_text(json.dumps(report))
    monkeypatch.setattr(controls, 'OUT', audit)
    monkeypatch.setattr(controls, 'FIT', fit)
    monkeypatch.setattr(controls, 'PROTOCOL', protocol)
    controls.require_inner_gate()
    pred.write_text('replaced predictions')
    with pytest.raises(RuntimeError, match='predictions changed'):
        controls.require_inner_gate()


def test_instrument_prediction_must_match_verified_native_scores(tmp_path, monkeypatch):
    monkeypatch.setattr(controls, 'GRAFF', tmp_path)
    directory = tmp_path / 'train_inner'
    directory.mkdir()
    scores = np.asarray([[.8, .2, .6], [.6, .1, .5]])
    np.save(directory / 'SCORES.npy', scores)
    np.save(directory / 'INSTRUMENT_SCORES.npy', np.repeat(scores[:, None, :], 4, axis=1))
    receipt = {'positions': 2, 'score_sha256': controls.sha(directory / 'SCORES.npy'),
               'instrument_scores_sha256': controls.sha(directory / 'INSTRUMENT_SCORES.npy')}
    (directory / 'COMPLETE.json').write_text(json.dumps(receipt))
    held = {'q': {'keys': ['a', 'b'], 'offset': 0, 'n': 2}}
    values = {'native_graff_' + label + side: {'q': 'a'}
              for label in ['', 'instrument_0_', 'instrument_1_', 'instrument_2_', 'observed_or_author_default_']
              for side in ['full', 'root', 'paired']}
    controls.verify_controls('inner', held, values, {'trained_variants': []})
    values['native_graff_instrument_2_root']['q'] = 'b'
    with pytest.raises(RuntimeError, match='disagreement'):
        controls.verify_controls('inner', held, values, {'trained_variants': []})
