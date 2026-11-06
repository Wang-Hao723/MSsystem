"""Declared charge evidence survives the actual web parsing and serialization paths."""
import ast
import copy
from pathlib import Path

import pytest
from pipeline.mzml_metadata import query_acquisition_projection
from pipeline.query_charge_guard import query_adduct_charge_evidence
from tests.test_query_acquisition_metadata import write_query, NCE
from tests.test_regression_smoke import _mock_row
from web import server
from web.report_common import classification_evidence_status

APP_PATH = Path(__file__).resolve().parents[1] / 'web/app.py'


def _app_function(name):
    # Execute this function's actual source without starting the Streamlit UI.
    nodes = [n for n in ast.parse(APP_PATH.read_text()).body
             if isinstance(n, ast.FunctionDef) and n.name in {name, '_msms_fragments_str'}]
    namespace = dict(server.__dict__, query_acquisition_projection=query_acquisition_projection,
                     classification_evidence_status=classification_evidence_status)
    exec(compile(ast.Module(body=nodes, type_ignores=[]), str(APP_PATH), 'exec'), namespace)
    return namespace[name]


@pytest.mark.parametrize('surface', ['fastapi', 'streamlit'])
@pytest.mark.parametrize('token', ['0', '-0', '-1'])
def test_declared_xml_charge_survives_json_projection(tmp_path, surface, token):
    path = tmp_path/'query.mzML'
    write_query(path, NCE, polarity='negative' if token.startswith('-') else 'positive')
    path.write_text(path.read_text().replace('name="charge state" value="1"',
                                             f'name="charge state" value="{token}"'))
    parser = server._parse_spectrum_file if surface == 'fastapi' else _app_function('_parse_spectrum_file')
    serialize = server._serialize_results if surface == 'fastapi' else _app_function('_serialize_results')
    record, = parser(str(path))
    result, = serialize([{**record, **query_acquisition_projection(record),
                          'raw_peaks':record['peaks'], 'query_peaks':record['peaks'], 'match':None}])
    assert result['reported_precursor_charge'] == token
    assert result['precursor_charge_status'] == ('unknown_zero' if int(token) == 0 else 'explicit')
    assert result['adduct_charge_evidence']['reported_charge'] == token
    assert result['adduct_charge_evidence']['conflict'] is False


@pytest.mark.parametrize('surface', ['fastapi', 'streamlit'])
def test_real_xml_noninteger_charge_cannot_become_integer(tmp_path, surface):
    path=tmp_path/'query.mzML';write_query(path,NCE)
    path.write_text(path.read_text().replace('name="charge state" value="1"',
                                             'name="charge state" value="1.5"'))
    parser=server._parse_spectrum_file if surface=='fastapi' else _app_function('_parse_spectrum_file')
    with pytest.raises(ValueError, match='integer literal'):
        parser(str(path))


@pytest.mark.parametrize('surface', ['fastapi', 'streamlit'])
def test_search_neutral_rejection_evidence_is_not_dropped(surface):
    row=_mock_row();evidence=query_adduct_charge_evidence('[M+H]+',2,'POS')
    assert evidence['conflict']
    row['match']['query_adduct_charge_evidence']=evidence
    serialize=server._serialize_results if surface=='fastapi' else _app_function('_serialize_results')
    result,=serialize([copy.deepcopy(row)])
    assert result['match']['query_adduct_charge_evidence']==evidence
