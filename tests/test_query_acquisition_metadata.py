"""Synthetic mzML contracts: explicit quantities, no reference library or labels."""
import base64
import struct

import pytest
from pyteomics.auxiliary import unitfloat

from pipeline.mzml_metadata import acquisition_metadata, read_instrument_header
from pipeline.mzml_input import parse_mzml_file


def cv(accession, name, value="", unit=""):
    return f'<cvParam cvRef="MS" accession="{accession}" name="{name}" value="{value}" {unit}/>'


def write_query(path, activation, *, model="Q Exactive", polarity="positive", second_model=None):
    mz = base64.b64encode(struct.pack('<2d', 100., 150.)).decode()
    intensity = base64.b64encode(struct.pack('<2d', 10., 20.)).decode()
    params = cv('MS:1001911', model)
    other = f'<instrumentConfiguration id="IC2">{cv("MS:1003028", second_model)}</instrumentConfiguration>' if second_model else ''
    scanref = 'instrumentConfigurationRef="IC2"' if second_model else ''
    path.write_text(f'''<?xml version="1.0"?><mzML xmlns="http://psi.hupo.org/ms/mzml" version="1.1.0">
    <referenceableParamGroupList count="1"><referenceableParamGroup id="common">{params}</referenceableParamGroup></referenceableParamGroupList>
    <instrumentConfigurationList count="{2 if second_model else 1}"><instrumentConfiguration id="IC1"><referenceableParamGroupRef ref="common"/>
    <componentList count="2"><source order="1">{cv('MS:1000073','electrospray ionization')}</source><analyzer order="2">{cv('MS:1000484','orbitrap')}</analyzer></componentList>
    </instrumentConfiguration>{other}</instrumentConfigurationList>
    <run id="test" defaultInstrumentConfigurationRef="IC1"><spectrumList count="1">
    <spectrum id="controllerType=0 controllerNumber=1 scan=7" index="0" defaultArrayLength="2">
    {cv('MS:1000511','ms level',2)}{cv('MS:1000127','centroid spectrum')}{cv('MS:1000130' if polarity == 'positive' else 'MS:1000129',polarity+' scan')}
    <scanList count="1"><scan {scanref}>{cv('MS:1000016','scan start time',120,'unitName="second" unitAccession="UO:0000010" unitCvRef="UO"')}
    {cv('MS:1000927','ion injection time',20,'unitName="millisecond" unitAccession="UO:0000028" unitCvRef="UO"')}
    {cv('MS:1000800','mass resolving power',17500)}</scan></scanList>
    <precursorList count="1"><precursor><selectedIonList count="1"><selectedIon>{cv('MS:1000744','selected ion m/z',300.1234)}{cv('MS:1000041','charge state',1)}</selectedIon></selectedIonList>
    <activation>{cv('MS:1000422','beam-type collision-induced dissociation')}{activation}</activation></precursor></precursorList>
    <binaryDataArrayList count="2"><binaryDataArray encodedLength="24">{cv('MS:1000523','64-bit float')}{cv('MS:1000576','no compression')}{cv('MS:1000514','m/z array')}<binary>{mz}</binary></binaryDataArray>
    <binaryDataArray encodedLength="24">{cv('MS:1000523','64-bit float')}{cv('MS:1000576','no compression')}{cv('MS:1000515','intensity array')}<binary>{intensity}</binary></binaryDataArray></binaryDataArrayList>
    </spectrum></spectrumList></run></mzML>''')


NCE = cv('MS:1000138','normalized collision energy',30,'unitName="percent" unitAccession="UO:0000187" unitCvRef="UO"')
EV = cv('MS:1000045','collision energy',18,'unitName="electronvolt" unitAccession="UO:0000266" unitCvRef="UO"')


@pytest.mark.parametrize('activation, expected', [(NCE,(30.,'NCE')),(EV,(18.,'eV')),(EV+NCE,(18.,'eV'))])
def test_typed_energy_survives_entire_query_parser(tmp_path, activation, expected):
    path = tmp_path/'query.mzML'
    write_query(path, activation)
    row, = parse_mzml_file(path)
    assert (row['collision_energy'],row['collision_energy_unit']) == expected
    assert row['collision_energy_mode'] == 'fixed'
    assert row['instrument'] == 'Q Exactive'
    assert row['instrument_type'] == 'orbitrap'
    assert row['instrument_configuration_id'] == 'IC1'
    assert row['ion_injection_time']['value'] == 20.
    assert row['ion_injection_time']['unit'] == 'millisecond'
    assert row['mass_resolving_power']['value'] == 17500.
    assert row['retention_time'] == 2.
    assert row['precursor_mz'] == 300.1234
    assert row['peaks'] == [(100.,10.),(150.,20.)]
    if activation == EV+NCE:
        assert [(q['value'],q['unit']) for q in row['collision_energy_quantities']] == [(18.,'eV'),(30.,'NCE')]


def test_generic_percentage_is_not_normalized_energy(tmp_path):
    path=tmp_path/'percent.mzML'
    write_query(path, cv('MS:1000045','collision energy',30,'unitName="percent" unitAccession="UO:0000187" unitCvRef="UO"'))
    row,=parse_mzml_file(path)
    assert row['collision_energy_unit'] == 'percent'


def test_explicit_nce_precedes_untyped_percent(tmp_path):
    path=tmp_path/'percent_nce.mzML'
    write_query(path,cv('MS:1000045','collision energy',30,'unitName="percent" unitAccession="UO:0000187" unitCvRef="UO"')+NCE)
    row,=parse_mzml_file(path)
    assert row['collision_energy_unit'] == 'NCE'


def test_multi_precursor_spectrum_cannot_be_labeled_single_fixed_condition():
    raw={'precursorList':{'precursor':[{'activation':{'collision energy':unitfloat(20,'electronvolt')}},
                                      {'activation':{'collision energy':unitfloat(40,'electronvolt')}}]}}
    with pytest.raises(ValueError,match='multi-precursor MS2'):
        acquisition_metadata(raw,{'default':'','configurations':{}})


def test_per_scan_instrument_and_negative_polarity(tmp_path):
    path=tmp_path/'multi.mzML'
    write_query(path, NCE, second_model='Orbitrap Exploris 240', polarity='negative')
    row,=parse_mzml_file(path)
    assert row['instrument'] == 'Orbitrap Exploris 240'
    assert row['instrument_configuration_id'] == 'IC2'
    assert row['polarity'] == 'NEG'


def test_unknown_scan_configuration_does_not_inherit_qe(tmp_path):
    path=tmp_path/'unknown.mzML'
    write_query(path,NCE,second_model='unreported analyzer')
    row,=parse_mzml_file(path)
    assert row['instrument_configuration_id'] == 'IC2'
    assert row['instrument'] == row['instrument_type'] == ''


def test_stepped_nce_is_preserved_without_a_synthetic_average(tmp_path):
    path=tmp_path/'stepped.mzML'
    write_query(path,NCE+cv('MS:1000138','normalized collision energy',50))
    row,=parse_mzml_file(path)
    assert row['collision_energy'] is None
    assert row['collision_energy_mode'] == 'stepped'
    assert [q['value'] for q in row['collision_energy_quantities']] == [30.,50.]


@pytest.mark.parametrize('bad',[float('nan'),float('inf'),-20])
def test_invalid_energies_do_not_become_conditions(bad):
    raw={'precursorList':{'precursor':[{'activation':{'normalized collision energy':unitfloat(bad,'percent')}}]}}
    result=acquisition_metadata(raw,{'default':'','configurations':{}})
    assert result['collision_energy'] is None
    assert result['collision_energy_mode'] == 'missing'


def test_header_stops_before_malformed_spectral_body(tmp_path):
    path=tmp_path/'header.mzML'
    write_query(path,NCE)
    text=path.read_text()
    path.write_text(text[:text.index('<spectrum id=')]+'incomplete binary payload')
    header=read_instrument_header(path)
    assert header['configurations']['IC1']['instrument'] == 'Q Exactive'


def test_query_adapter_does_not_modify_library_parser_globals(tmp_path):
    from utils import mzml_io
    path=tmp_path/'query.mzML';write_query(path,NCE)
    original={key:mzml_io.parse_mzml_file.__globals__[key] for key in ('load_from_mzml','_extract_mzml_precursor_evidence','extract_mzml_instrument_metadata')}
    parse_mzml_file(path)
    assert all(mzml_io.parse_mzml_file.__globals__[key] is value for key,value in original.items())


def test_analysis_projection_and_serialization_keep_conditions(tmp_path):
    from pipeline.mzml_metadata import query_acquisition_projection
    from web.server import _serialize_results
    path=tmp_path/'query.mzML';write_query(path,EV+NCE)
    row,=parse_mzml_file(path)
    projected=query_acquisition_projection(row)
    result,=_serialize_results([{**row, **projected, 'raw_peaks':row['peaks'], 'match':None}])
    assert all(result[key] == value for key,value in projected.items())
    assert result['collision_energy_unit'] == 'eV'
    assert result['activation_method'] == 'beam-type collision-induced dissociation'
    assert result['collision_energy_mode'] == 'fixed'
