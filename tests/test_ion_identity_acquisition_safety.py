"""IIMN must use observed peak shapes from same-polarity MS1 scans."""
import base64
import struct
from pathlib import Path

import numpy as np
import pytest
from pyteomics import mzml, mzxml

import matching.ion_identity as legacy
from pipeline.ion_identity_evidence import (
    annotate_ion_identity_networks,
    chromatographic_shape_evidence,
    extract_targeted_ms1_traces,
)


def _trace(values=(1, 5, 100, 5, 1)):
    return {'rt': [5 + i * .01 for i in range(len(values))], 'intensity': list(values)}


def _rows():
    return [dict(precursor_feature_id=key, precursor_mz=mz, retention_time=5.02,
                 ion_mode='POS', match=None)
            for key, mz in [('P1', 300 + legacy.PROTON), ('P2', 300 + 22.989218)]]


def _array(values, accession, name):
    data = base64.b64encode(struct.pack('<' + 'd' * len(values), *values)).decode()
    return (f'<binaryDataArray encodedLength="{len(data)}">'
            '<cvParam cvRef="MS" accession="MS:1000523" name="64-bit float"/>'
            '<cvParam cvRef="MS" accession="MS:1000576" name="no compression"/>'
            f'<cvParam cvRef="MS" accession="{accession}" name="{name}"/>'
            f'<binary>{data}</binary></binaryDataArray>')


def _mixed_mzml(path: Path, signal_polarity='NEG'):
    spectra = []
    shape = [1, 5, 20, 5, 1]
    for index in range(10):
        positive = index % 2 == 0
        has_signal = ('POS' if positive else 'NEG') == signal_polarity
        values = [shape[index // 2], 3 * shape[index // 2]] if has_signal else [0, 0]
        spectra.append(
            f'<spectrum index="{index}" id="scan={index + 1}" defaultArrayLength="2">'
            '<cvParam cvRef="MS" accession="MS:1000511" name="ms level" value="1"/>'
            f'<cvParam cvRef="MS" accession="MS:{1000130 if positive else 1000129}" '
            f'name="{"positive" if positive else "negative"} scan"/>'
            '<scanList count="1"><scan><cvParam cvRef="MS" accession="MS:1000016" '
            f'name="scan start time" value="{5 + index * .005}" unitCvRef="UO" '
            'unitAccession="UO:0000031" unitName="minute"/></scan></scanList>'
            '<binaryDataArrayList count="2">'
            + _array([300 + legacy.PROTON, 300 + 22.989218], 'MS:1000514', 'm/z array')
            + _array(values, 'MS:1000515', 'intensity array')
            + '</binaryDataArrayList></spectrum>')
    path.write_text('<?xml version="1.0"?><mzML xmlns="http://psi.hupo.org/ms/mzml" '
                    'version="1.1.0"><run id="synthetic"><spectrumList count="10">'
                    + ''.join(spectra) + '</spectrumList></run></mzML>')


def _targets():
    return [dict(id=row['precursor_feature_id'], mz=row['precursor_mz'],
                 rt=row['retention_time'], polarity=row['ion_mode']) for row in _rows()]


@pytest.mark.parametrize('values', [(0, 0, 100, 0, 0), (0, 1, 100, 1, 0),
                                   (1, 0, 100, 1, 1)])
def test_shared_spikes_and_missing_peak_flanks_cannot_define_network(values):
    trace = _trace(values)
    assert legacy.chromatographic_shape_evidence(trace, trace)['accepted']
    evidence = chromatographic_shape_evidence(trace, trace)
    assert evidence['accepted'] is False
    assert evidence['status'] == 'insufficient_nonzero_peak_shape_points'
    summary = annotate_ion_identity_networks(_rows(), traces={'P1': trace, 'P2': trace})
    assert summary['network_count'] == 0
    assert summary['runtime_influence'] is False


def test_observed_smooth_peak_is_accepted_and_legacy_module_is_not_mutated():
    trace = _trace()
    assert chromatographic_shape_evidence(trace, trace)['accepted']
    summary = annotate_ion_identity_networks(_rows(), traces={'P1': trace, 'P2': trace})
    assert summary['network_count'] == 1
    assert legacy.annotate_ion_identity_networks.__globals__['chromatographic_shape_evidence'] is legacy.chromatographic_shape_evidence
    assert legacy.annotate_ion_identity_networks.__globals__['extract_targeted_ms1_traces'] is legacy.extract_targeted_ms1_traces


@pytest.mark.parametrize('signal_polarity,expected_networks', [('NEG', 0), ('POS', 1)])
def test_polarity_switching_uses_only_same_polarity_ms1(tmp_path, signal_polarity, expected_networks):
    path = tmp_path / 'switching.mzML'
    _mixed_mzml(path, signal_polarity)
    with mzml.read(str(path)) as reader:
        scans = list(reader)
    assert sum('positive scan' in scan for scan in scans) == 5
    assert sum('negative scan' in scan for scan in scans) == 5
    traces, extraction = extract_targeted_ms1_traces(path, _targets())
    assert extraction['scan_polarity_counts'] == {'POS': 5, 'NEG': 5}
    assert traces['P1']['points'] == 5
    assert traces['P1']['rt'] == pytest.approx([5, 5.01, 5.02, 5.03, 5.04])
    expected = [1, 5, 20, 5, 1] if signal_polarity == 'POS' else [0] * 5
    assert traces['P1']['intensity'] == expected
    assert annotate_ion_identity_networks(_rows(), file_path=path)['network_count'] == expected_networks


def test_missing_scan_polarity_cannot_be_used_as_positive_evidence(tmp_path):
    path = tmp_path / 'unknown.mzML'
    _mixed_mzml(path, 'POS')
    import re
    path.write_text(re.sub(r'<cvParam cvRef="MS" accession="MS:10001(?:29|30)" name="(?:positive|negative) scan"/>', '', path.read_text()))
    traces, extraction = extract_targeted_ms1_traces(path, _targets())
    assert extraction['scan_polarity_counts'] == {'UNKNOWN': 10}
    assert all(trace['points'] == 0 for trace in traces.values())
    assert annotate_ion_identity_networks(_rows(), file_path=path)['network_count'] == 0


@pytest.mark.parametrize('extension', ['.mzXML', '.xml'])
def test_mzxml_reader_retains_real_ms1_arrays_and_retention_units(tmp_path, extension):
    scans = []
    for index, intensity in enumerate([1, 5, 20, 5, 1], start=1):
        encoded = base64.b64encode(struct.pack('>dd', 300 + legacy.PROTON, intensity)).decode()
        scans.append(f'<scan num="{index}" msLevel="1" peaksCount="1" polarity="+" '
                     f'retentionTime="PT{300 + index}S"><peaks precision="64" '
                     'byteOrder="network" pairOrder="m/z-int" compressionType="none">'
                     f'{encoded}</peaks></scan>')
    path = tmp_path / ('synthetic' + extension)
    path.write_text('<?xml version="1.0"?><mzXML xmlns="http://sashimi.sourceforge.net/schema_revision/mzXML_3.2" '
                    'version="3.2"><msRun scanCount="5">' + ''.join(scans) + '</msRun></mzXML>')
    with mzxml.read(str(path)) as reader:
        independent_scans = list(reader)
    assert len(independent_scans) == 5
    traces, extraction = extract_targeted_ms1_traces(path, _targets()[:1])
    assert extraction['status'] == 'ok'
    assert extraction['source_format'] == 'mzXML'
    assert extraction['ms1_scan_count'] == 5
    assert traces['P1']['intensity'] == [1, 5, 20, 5, 1]
    assert traces['P1']['rt'] == pytest.approx([(300 + i) / 60 for i in range(1, 6)])
