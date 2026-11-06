"""Generation checks must finish before batches or individual rows are exposed."""
import sqlite3
import importlib.util
from pathlib import Path
import pytest
from runtime import verified_sql_retrieval as M
from runtime import source_compound_mass_repairs as mass
from matching.evirank_retrieval import SQLiteCandidateRecallRetriever
def _sibling_fixture_module(name):
    spec = importlib.util.spec_from_file_location(name, Path(__file__).with_name(name + '.py'))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


bound = _sibling_fixture_module('test_verified_sql_retrieval').bound
_mass_fixtures = _sibling_fixture_module('test_source_compound_mass_repairs')
fixture, _state = _mass_fixtures.fixture, _mass_fixtures._state


def verified(bound):
    root, db, backing, registry, bindings, config, raw = bound
    return M.VerifiedSQLBindings(config, bindings, registry, backing, raw,
        project_root=root, dependencies=[db, root/'source.mgf'])


def raw_rows(bound):
    conn=bound[-1]._connect()
    try:
        return [dict(x) for x in conn.execute('''SELECT s.id AS database_spectrum_id,
            s.spectrum_id,s.source_file,s.precursor_mz,s.peaks_json,s.charge,s.adduct,
            s.ion_mode,s.source_accession,c.smiles,c.inchikey,c.exact_mass,
            c.molecular_formula FROM Spectra s JOIN Mapping m ON m.id=s.id
            JOIN Compounds c ON c.id=m.compound_id ORDER BY s.id''').fetchall()]
    finally: conn.close()


class Rows:
    def __init__(self, rows, callback=None):
        self.rows=list(rows);self.callback=callback
    def take(self,n):
        if self.callback:
            callback,self.callback=self.callback,None;callback()
        result=self.rows[:n];del self.rows[:n];return result
    def fetchall(self):return self.take(len(self.rows))
    def fetchmany(self,size=1):return self.take(size)
    def fetchone(self):
        result=self.take(1);return result[0] if result else None
    def __iter__(self):return iter(self.rows)


def mutation(bound,kind):
    if kind=='source':
        source=bound[0]/'source.mgf';source.write_text(source.read_text()+'; changed\n')
    else:
        with sqlite3.connect(bound[1]) as c:
            c.execute('UPDATE Compounds SET compound_name=? WHERE id=2',('changed',))


@pytest.mark.parametrize('operation',['fetchall','fetchmany','fetchone'])
@pytest.mark.parametrize('kind',['source','database'])
def test_batch_change_during_row_processing_rejected_before_return(bound,monkeypatch,operation,kind):
    state=verified(bound);rows=raw_rows(bound);original=state._apply_checked_row
    def changed(row):
        result=original(row);mutation(bound,kind);return result
    monkeypatch.setattr(state,'_apply_checked_row',changed)
    cursor=M._BoundCursor(Rows(rows),state);returned=[]
    with pytest.raises(ValueError,match='dependencies changed'):
        returned.append(getattr(cursor,operation)())
    assert returned==[]


@pytest.mark.parametrize('operation',['fetchall','fetchmany','fetchone'])
@pytest.mark.parametrize('empty',[False,True])
def test_generation_change_during_cursor_read_even_empty_rejected(bound,operation,empty):
    state=verified(bound);rows=[] if empty else raw_rows(bound)
    cursor=M._BoundCursor(Rows(rows,lambda:mutation(bound,'source')),state)
    with pytest.raises(ValueError,match='dependencies changed'):
        getattr(cursor,operation)()


@pytest.mark.parametrize('operation',['fetchall','fetchmany','fetchone'])
def test_batch_entry_rejects_changed_generation_without_processing(bound,monkeypatch,operation):
    state=verified(bound);cursor=M._BoundCursor(Rows(raw_rows(bound)),state);called=[]
    monkeypatch.setattr(state,'_apply_checked_row',lambda row:called.append(row))
    mutation(bound,'source')
    with pytest.raises(ValueError,match='dependencies changed'):getattr(cursor,operation)()
    assert called==[]


def test_iterator_rechecks_when_caller_resumes(bound):
    state=verified(bound);rows=raw_rows(bound);it=iter(M._BoundCursor(Rows(rows),state))
    first=next(it);assert first['charge']==-1
    mutation(bound,'source')
    with pytest.raises(ValueError,match='dependencies changed'):next(it)


@pytest.mark.parametrize('stream',[False,True])
def test_individual_row_change_during_processing_rejected_before_exposure(bound,monkeypatch,stream):
    state=verified(bound);rows=raw_rows(bound);original=state._apply_checked_row
    def changed(row):
        result=original(row);mutation(bound,'source');return result
    monkeypatch.setattr(state,'_apply_checked_row',changed)
    with pytest.raises(ValueError,match='dependencies changed'):
        next(iter(M._BoundCursor(Rows(rows),state))) if stream else state.apply(rows[0])


@pytest.mark.parametrize('prime_cache',[False,True])
@pytest.mark.parametrize('kind',['source','database'])
def test_window_final_check_covers_normalization_and_cache_return(bound,monkeypatch,prime_cache,kind):
    state=verified(bound);fixed=M.VerifiedSQLiteCandidateRecallRetriever(bound[1],library_snapshot_id='fixture-snapshot',bindings=state)
    if prime_cache:fixed._fetch_window(59.01,'NEG')
    original=SQLiteCandidateRecallRetriever._fetch_window;seen=[]
    def altered(self,*args):
        records,audit=original(self,*args);seen.append(audit.get('cache_hit',False));mutation(bound,kind);return records,audit
    monkeypatch.setattr(SQLiteCandidateRecallRetriever,'_fetch_window',altered)
    with pytest.raises(ValueError,match='dependencies changed'):fixed._fetch_window(59.01,'NEG')
    assert seen==[prime_cache]


def test_stale_held_connection_rejected_before_candidate_execute(bound):
    state=verified(bound)
    class Connection:
        called=False
        def execute(self,*args):self.called=True;raise AssertionError('stale SQL executed')
    c=Connection();wrapped=M._BoundConnection(c,state);mutation(bound,'source')
    with pytest.raises(ValueError,match='dependencies changed'):wrapped.execute(M._SELECT_PREFIX)
    assert c.called is False


def test_mass_standalone_rechecks_before_return(fixture,monkeypatch):
    f=fixture;state=_state(f,dependencies=[f.source,f.db]);original=state._apply_sql_checked_row
    def changed(row):
        result=original(row);f.source.write_text('modified');return result
    monkeypatch.setattr(state,'_apply_sql_checked_row',changed)
    with pytest.raises(ValueError,match='dependencies changed'):state.apply_sql(f.raw)


def test_batch_keeps_complete_original_row_digest(bound,monkeypatch):
    state=verified(bound);rows=raw_rows(bound);rows[0]['smiles']='CCO'
    with pytest.raises(ValueError,match='original candidate changed'):
        M._BoundCursor(Rows(rows),state).fetchall()


def test_batch_repair_preserves_original_order_values_and_only_checks_boundaries(bound,monkeypatch):
    state=verified(bound);rows=raw_rows(bound);expected=[state.apply(row) for row in rows];calls=[];original=state.assert_fresh
    def counted():calls.append(1);return original()
    monkeypatch.setattr(state,'assert_fresh',counted)
    cursor=M._BoundCursor(Rows(rows),state)
    assert cursor.fetchall()==expected;assert len(calls)==2
    assert cursor.fetchone() is None;assert len(calls)==4
