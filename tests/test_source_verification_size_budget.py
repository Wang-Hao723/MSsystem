"""Cold-drive deadlines scale with bound bytes; no evidence check is skipped."""
import math
import pytest
from runtime import source_compound_mass_core as core
from runtime import source_file_verification as verifier


@pytest.mark.parametrize('size,expected', [
    (0,120.), (1024,120.), (960*1024**2,120.),
    (4533916079,60+4533916079/(16*1024**2)),
    (5121690843,60+5121690843/(16*1024**2)),
    (1740*16*1024**2,1800.), (10**1000,1800.),
])
def test_registered_size_budget_is_finite_and_bounded(size, expected):
    actual=core.source_verification_timeout({'size_bytes':size})
    assert actual==expected and math.isfinite(actual) and 120<=actual<=1800


@pytest.mark.parametrize('size', [None,True,False,-1,1.0,'500',float('inf')])
def test_invalid_sizes_are_rejected_before_source_access(size):
    with pytest.raises(ValueError,match='source size'):
        core.source_verification_timeout({'size_bytes':size})


def test_every_whole_source_keeps_exact_content_identity_and_one_attempt(monkeypatch):
    files={'one':{'size_bytes':4533916079,'sha256':'a'*64},
           'two':{'size_bytes':5121690843,'sha256':'b'*64}}
    observed=[]
    def verify(name, expected, **kwargs):
        assert expected is files[name]
        observed.append((name,expected.copy(),kwargs))
    monkeypatch.setattr(verifier,'verify_source_file',verify)
    core.verify_source_evidence({'source_files':files},{})
    assert [x[0] for x in observed]==list(files)
    for name,expected,kwargs in observed:
        assert expected==files[name]
        assert kwargs=={'timeout_seconds':core.source_verification_timeout(expected),'attempts':1}


def test_source_failure_aborts_before_any_later_source(monkeypatch):
    visited=[]
    def verify(name, expected, **kwargs):
        visited.append(name)
        raise TimeoutError('still unreadable')
    monkeypatch.setattr(verifier,'verify_source_file',verify)
    with pytest.raises(TimeoutError,match='still unreadable'):
        core.verify_source_evidence({'source_files':{'one':{'size_bytes':0},'two':{'size_bytes':0}}},{})
    assert visited==['one']
