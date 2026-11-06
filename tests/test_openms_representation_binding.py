"""Exercise the installed PyOpenMS binding in the production isolated environment."""
import subprocess
import sys
from pathlib import Path

import pytest

from utils.openms_backend import _runtime_environment


def test_real_openms_estimator_and_explicit_representation_contract():
    root = Path(__file__).resolve().parents[1]
    code = r'''
import importlib.util
import numpy as np
import pyopenms as oms
from pathlib import Path
spec = importlib.util.spec_from_file_location('representation_runner', Path('scripts/openms_untargeted.py'))
runner = importlib.util.module_from_spec(spec)
spec.loader.exec_module(runner)

centroid = oms.MSSpectrum()
centroid.set_peaks((np.arange(100.,200.,10.), np.ones(10)))
assert centroid.getType() == runner.UNKNOWN
assert runner._spectrum_representation(centroid) == ('centroid','openms_signal_shape')

profile = oms.MSSpectrum()
mz = np.linspace(100.,110.,101)
profile.set_peaks((mz,np.exp(-((mz-105.)**2))))
assert runner._spectrum_representation(profile) == ('profile','openms_signal_shape')

sparse = oms.MSSpectrum()
sparse.set_peaks(([100.,150.,200.],[2.,4.,3.]))
assert runner._spectrum_representation(sparse) == ('unknown','openms_signal_shape')

# The decision must survive serialization; preserving UNKNOWN would cause the
# post-conversion mzML guard to reject the otherwise successful normalization.
centroid.setRT(123.4)
centroid.setNativeID('scan=17')
centroid.setMSLevel(2)
before = centroid.get_peaks()
experiment = oms.MSExperiment()
experiment.addSpectrum(centroid)
counts = runner._pick_profile_spectra(experiment,strategy='hires',signal_to_noise=1.)
converted = experiment[0]
assert converted.getType() == runner.CENTROID
assert converted.getRT() == 123.4 and converted.getNativeID() == 'scan=17'
assert all(np.array_equal(a,b) for a,b in zip(before,converted.get_peaks()))
assert counts['2']['inferred_centroid'] == 1 and counts['2']['profile_picked'] == 0

# Explicit metadata always wins; estimation must not relabel these spectra.
profile.setType(runner.CENTROID)
assert runner._spectrum_representation(profile) == ('centroid','mzml_cv')
centroid.setType(runner.PROFILE)
assert runner._spectrum_representation(centroid) == ('profile','mzml_cv')
assert sparse.getType() == runner.UNKNOWN
print('real OpenMS representation contracts passed')
'''
    result = subprocess.run([sys.executable, '-c', code], cwd=root,
                            env=_runtime_environment(), capture_output=True,
                            text=True, timeout=60)
    assert result.returncode == 0, result.stdout + result.stderr
    assert 'real OpenMS representation contracts passed' in result.stdout


def test_openms_empty_spectrum_keeps_following_precursor_and_scan_metadata(tmp_path):
    path = tmp_path / 'with_empty.mzML'
    code = r'''
import sys
import numpy as np
import pyopenms as oms
exp=oms.MSExperiment()
for i in (1,2,3):
    s=oms.MSSpectrum()
    s.setNativeID('scan='+str(i)); s.setMSLevel(2); s.setRT(60.*i)
    s.setType(oms.SpectrumSettings.SpectrumType.CENTROID)
    settings=s.getInstrumentSettings(); settings.setPolarity(oms.IonSource.Polarity.POSITIVE)
    s.setInstrumentSettings(settings)
    precursor=oms.Precursor(); precursor.setMZ(100.*i+1.); precursor.setCharge(1)
    s.setPrecursors([precursor])
    s.set_peaks((np.array([]),np.array([])) if i==2 else (np.array([50.,75.]),np.array([10.,20.])))
    exp.addSpectrum(s)
oms.MzMLFile().store(sys.argv[1],exp)
'''
    result = subprocess.run([sys.executable, '-c', code, str(path)],
                            env=_runtime_environment(), capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stdout + result.stderr
    from pipeline.mzml_input import parse_mzml_file
    from utils import mzml_io
    original_reader = mzml_io.load_from_mzml
    records = parse_mzml_file(path)
    assert [r['spectrum_id'] for r in records] == [1, 3]
    assert [r['precursor_mz'] for r in records] == [101., 301.]
    assert [r['retention_time'] for r in records] == [1., 3.]
    assert all(r['polarity'] == 'POS' for r in records)
    assert mzml_io.load_from_mzml is original_reader
    corrupted = tmp_path / 'missing_nonempty_arrays.mzML'
    corrupted.write_text(path.read_text().replace('defaultArrayLength="0"', 'defaultArrayLength="3"'))
    with pytest.raises(ValueError, match='missing required peak arrays'):
        parse_mzml_file(corrupted)
