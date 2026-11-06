"""Real shared-coordinator races and the actual Streamlit worker/start sources."""
import ast
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import threading
import time
from types import SimpleNamespace
import pytest
from runtime import legacy_library_session as M

APP=Path(__file__).resolve().parents[1]/'web/app.py'


def coordinator():return M.LibraryCoordinator(runtime_binding=lambda: ())

def publish(c,records=None,index=None):
    records=['reference'] if records is None else records;index={1:[0]} if index is None else index
    return c.replace(lambda:(records,index),lambda:None)


def session(snapshot,task_id='task'):
    return SimpleNamespace(records=snapshot.records,idx=snapshot.index,lib_loaded=True,
        lib_name='Full',lib_spectra_count=len(snapshot.records),lib_load_time=1.,
        analysis_running=False,active_task_id=None,tasks={task_id:{'id':task_id,'status':'pending','file_path':'fixture.mzML','params':{'enable_ion_identity_networking':False,'enable_chimeric_deconvolution':False}}})


def app_scope(c,state):
    from web import server
    names={'_load_library_into_session','_start_analysis_in_thread','_run_analysis_in_thread','_run_analysis_with_library'}
    nodes=[n for n in ast.parse(APP.read_text()).body if isinstance(n,ast.FunctionDef) and n.name in names]
    assert len(nodes)==4
    errors=[];saved=[]
    scope=dict(server.__dict__,_logger=server.logger,_library_coordinator=c,_tasks_lock=c.lock,
        capture_session_state=lambda:state,st=SimpleNamespace(session_state=state,error=errors.append),
        _task_registry=lambda path:SimpleNamespace(assert_current=lambda tasks=None:None),
        TASKS_FILE=Path('isolated-library-coordinator-test.json'),
        _now_iso=lambda:'fixed',_save_tasks_to_disk=lambda state:None,
        _save_task_results=lambda task,rows:saved.append(rows),
        _parse_spectrum_file=lambda *a,**k:[{'peaks':[[10.,1.]],'precursor_mz':100.},{'peaks':[[20.,1.]],'precursor_mz':101.}],
        _serialize_results=lambda rows:rows,annotate_precursor_features=lambda rows:None,
        search=lambda *a,**k:None,threading=threading,time=time)
    exec(compile(ast.Module(body=nodes,type_ignores=[]),str(APP),'exec'),scope)
    return scope,errors,saved


def test_reservation_blocks_replacement_before_prepare_and_after_thread_reservation():
    c=coordinator();s=publish(c);lease=c.reserve('task',s.records,s.index);prepared=[]
    with pytest.raises(M.LibraryBusyError,match='queued or running'):
        c.replace(lambda:([],{}),lambda:prepared.append(True))
    assert not prepared;c.validate(lease);c.release(lease);assert c.active_count()==0
    assert publish(c).generation!=s.generation


def test_old_session_rejected_after_another_session_reloads():
    c=coordinator();old=publish(c);new=publish(c)
    with pytest.raises(M.LibraryBusyError,match='current library generation'):c.reserve('a',old.records,old.index)
    lease=c.reserve('a',new.records,new.index);c.release(lease)


def test_same_task_id_cannot_start_twice_across_sessions():
    c=coordinator();s=publish(c);a=c.reserve('same',s.records,s.index)
    with pytest.raises(M.LibraryBusyError,match='already queued'):c.reserve('same',s.records,s.index)
    b=c.reserve('other',s.records,s.index);assert c.active_count()==2;c.release(a);c.validate(b);c.release(b)


def test_coordinator_and_lock_survive_script_rerun_imports():
    namespaces=[{},{}]
    for ns in namespaces:exec('from runtime.legacy_library_session import COORDINATOR\nlock=COORDINATOR.lock',ns)
    assert namespaces[0]['COORDINATOR'] is namespaces[1]['COORDINATOR'] is M.COORDINATOR
    assert namespaces[0]['lock'] is namespaces[1]['lock']


def test_capture_returns_concrete_session_not_proxy(monkeypatch):
    import streamlit.runtime.scriptrunner as runner
    concrete=object();monkeypatch.setattr(runner,'get_script_run_ctx',lambda:SimpleNamespace(session_state=concrete));assert M.capture_session_state() is concrete
    monkeypatch.setattr(runner,'get_script_run_ctx',lambda:None)
    with pytest.raises(RuntimeError,match='Streamlit session'):M.capture_session_state()


def test_replace_and_reserve_use_one_lock_under_race():
    c=coordinator();old=publish(c);entered=threading.Event();finish=threading.Event();attempt=threading.Event()
    def loader():entered.set();assert finish.wait(2);return ['new'],{2:[0]}
    def reserve():attempt.set();return c.reserve('task',old.records,old.index)
    with ThreadPoolExecutor(max_workers=2) as pool:
        replacing=pool.submit(c.replace,loader,lambda:None);assert entered.wait(2)
        reserving=pool.submit(reserve);assert attempt.wait(2);assert not reserving.done()
        finish.set();replacing.result(timeout=2)
        with pytest.raises(M.LibraryBusyError,match='current library generation'):reserving.result(timeout=2)
    assert c.active_count()==0


def test_reserved_before_scheduled_worker_blocks_other_session():
    c=coordinator();snap=publish(c);state=session(snap);scope,errors,_=app_scope(c,state);created=[]
    class WaitingThread:
        def __init__(self,**kw):created.append(kw);assert c.active_count()==1
        def start(self):
            with pytest.raises(M.LibraryBusyError):c.replace(lambda:([],{}),lambda:None)
    scope['threading']=SimpleNamespace(Thread=WaitingThread)
    assert scope['_start_analysis_in_thread']('task');assert state.tasks['task']['status']=='queued' and c.active_count()==1 and not errors
    scope['_run_analysis_in_thread'](*created[0]['args']);assert c.active_count()==0 and not state.analysis_running


@pytest.mark.parametrize('where',['construct','start'])
def test_thread_creation_or_start_failure_releases_and_reports(where):
    c=coordinator();snap=publish(c);state=session(snap);scope,errors,_=app_scope(c,state)
    class FailedThread:
        def __init__(self,**kw):
            if where=='construct':raise RuntimeError('cannot construct')
        def start(self):raise RuntimeError('cannot start')
    scope['threading']=SimpleNamespace(Thread=FailedThread)
    assert scope['_start_analysis_in_thread']('task') is False
    assert c.active_count()==0 and not state.analysis_running and state.active_task_id is None
    assert state.tasks['task']['status']=='failed' and errors;publish(c)


@pytest.mark.parametrize('failure',['parse','setup','missing_task'])
def test_worker_every_exit_releases_reservation(failure):
    binding=[object()];c=M.LibraryCoordinator(runtime_binding=lambda:tuple(binding));snap=publish(c);state=session(snap);scope,_,_=app_scope(c,state);lease=c.reserve('task',snap.records,snap.index);state.analysis_running=True;state.active_task_id='task'
    if failure=='parse':
        def fail(*a,**k):raise ValueError('bad spectrum')
        scope['_parse_spectrum_file']=fail
    if failure=='setup':binding[0]=object()
    if failure=='missing_task':state.tasks.clear()
    scope['_run_analysis_in_thread']('task',state,lease)
    assert c.active_count()==0 and not state.analysis_running and state.active_task_id is None
    if failure!='missing_task':assert state.tasks['task']['status']=='failed'


def test_each_spectrum_uses_fixed_library_even_if_session_proxy_changes():
    c=coordinator();snap=publish(c);state=session(snap);scope,_,saved=app_scope(c,state);calls=[]
    def search(*args,**kw):
        calls.append((kw['records'],kw['mz_index']));state.records=['wrong'];state.idx={999:[0]};return None
    scope['search']=search;lease=c.reserve('task',snap.records,snap.index);state.analysis_running=True;state.active_task_id='task'
    scope['_run_analysis_in_thread']('task',state,lease)
    assert len(calls)==2 and all(a is snap.records and b is snap.index for a,b in calls)
    assert len(saved)==1 and state.tasks['task']['status']=='completed' and c.active_count()==0


def test_external_global_index_change_fails_before_next_spectrum():
    binding=[object()];c=M.LibraryCoordinator(runtime_binding=lambda:tuple(binding));snap=publish(c);state=session(snap);scope,_,saved=app_scope(c,state);calls=[]
    def search(*a,**k):calls.append(True);binding[0]=object();return None
    scope['search']=search;lease=c.reserve('task',snap.records,snap.index)
    scope['_run_analysis_in_thread']('task',state,lease)
    assert calls==[True] and state.tasks['task']['status']=='failed' and not saved and c.active_count()==0


def test_failed_reload_clears_generation_and_prevents_other_old_session_start():
    c=coordinator();snap=publish(c);state=session(snap);scope,_,_=app_scope(c,state)
    def fail(path):raise ValueError('source changed')
    scope['load_library']=fail
    with pytest.raises(ValueError):scope['_load_library_into_session']('db','Full')
    assert state.records is None and state.idx is None and not state.lib_loaded
    with pytest.raises(M.LibraryBusyError):c.reserve('other',snap.records,snap.index)


def test_running_lease_rejects_reload_without_clearing_other_session():
    c=coordinator();snap=publish(c);first=session(snap);other=session(snap,'other');scope,_,_=app_scope(c,other);lease=c.reserve('task',first.records,first.idx)
    def forbidden(path):raise AssertionError('must reject before loader')
    scope['load_library']=forbidden
    with pytest.raises(M.LibraryBusyError):scope['_load_library_into_session']('db','Full')
    assert other.records is snap.records and other.lib_loaded;c.release(lease)


def test_all_three_entrypoints_use_reservation_and_only_one_thread_constructor():
    tree=ast.parse(APP.read_text());starts=[n for n in ast.walk(tree) if isinstance(n,ast.Call) and isinstance(n.func,ast.Name) and n.func.id=='_start_analysis_in_thread'];threads=[n for n in ast.walk(tree) if isinstance(n,ast.Call) and isinstance(n.func,ast.Attribute) and isinstance(n.func.value,ast.Name) and n.func.value.id=='threading' and n.func.attr=='Thread']
    assert len(starts)==3 and len(threads)==1
    worker=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='_run_analysis_with_library')
    assert 'st.session_state' not in ast.unparse(worker)


def test_actual_background_thread_holds_generation_through_completion():
    c=coordinator();snap=publish(c);state=session(snap);scope,errors,saved=app_scope(c,state)
    entered=threading.Event();finish=threading.Event();threads=[];calls=[]
    def search(*a,**kw):
        calls.append((kw['records'],kw['mz_index']))
        if len(calls)==1:entered.set();assert finish.wait(3)
        return None
    def make_thread(**kw):
        t=threading.Thread(**kw);threads.append(t);return t
    scope['search']=search;scope['threading']=SimpleNamespace(Thread=make_thread)
    assert scope['_start_analysis_in_thread']('task');assert entered.wait(3)
    try:
        assert c.active_count()==1
        with pytest.raises(M.LibraryBusyError):c.replace(lambda:(['replacement'],{}),lambda:None)
        # A rerun replaces script globals but imports this exact stable lock.
        again={};exec('from runtime.legacy_library_session import COORDINATOR',again)
        assert again['COORDINATOR'] is M.COORDINATOR
    finally:
        finish.set();threads[0].join(timeout=3)
    assert not threads[0].is_alive() and c.active_count()==0 and not errors
    assert len(calls)==2 and all(a is snap.records and b is snap.index for a,b in calls)
    assert state.tasks['task']['status']=='completed' and len(saved)==1


def test_duplicate_other_session_cannot_change_task_params_or_status():
    c=coordinator();snap=publish(c);first=c.reserve('task',snap.records,snap.index);other=session(snap);scope,errors,_=app_scope(c,other)
    old=dict(other.tasks['task']);oldparams=dict(old['params'])
    assert scope['_start_analysis_in_thread']('task',params={'cosine_threshold':.999}) is False
    assert other.tasks['task']==old and other.tasks['task']['params']==oldparams
    assert other.analysis_running is False and errors;c.release(first)


def test_start_after_other_session_switch_requires_reload_without_task_mutation():
    c=coordinator();snap=publish(c);state=session(snap);scope,errors,_=app_scope(c,state);publish(c)
    assert scope['_start_analysis_in_thread']('task',params={'cosine_threshold':.999}) is False
    assert state.tasks['task']['status']=='pending' and 'cosine_threshold' not in state.tasks['task']['params'] and errors


@pytest.mark.parametrize('button',[0,1])
@pytest.mark.parametrize('outcome',['success','failure','busy'])
def test_both_actual_load_button_handlers_report_result_without_stale_publication(button,outcome):
    c=coordinator();snap=publish(c);state=session(snap);scope,errors,_=app_scope(c,state);messages=[];calls=[]
    scope['st'].success=messages.append;scope['db_path']='fixture.db';scope['db_mode']='Full'
    newrecords=['new','other'];newindex={22:[0,1]}
    def loader(path):
        calls.append(path)
        if outcome=='failure':raise ValueError('source verification rejected')
        return newrecords,newindex
    scope['load_library']=loader
    lease=c.reserve('running-elsewhere',snap.records,snap.index) if outcome=='busy' else None
    tree=ast.parse(APP.read_text());handlers=[]
    for node in ast.walk(tree):
        if isinstance(node,ast.With) and any(isinstance(item,ast.Try) and any(isinstance(x,ast.Call) and isinstance(x.func,ast.Name) and x.func.id=='_load_library_into_session' for x in ast.walk(item)) for item in node.body):handlers.append(node.body)
    assert len(handlers)==2
    exec(compile(ast.Module(body=handlers[button],type_ignores=[]),str(APP),'exec'),scope)
    if outcome=='success':
        assert state.records is newrecords and state.idx is newindex and state.lib_loaded and state.lib_spectra_count==2
        assert len(messages)==1 and not errors and len(calls)==1
    elif outcome=='failure':
        assert state.records is None and state.idx is None and not state.lib_loaded and state.lib_spectra_count==0
        assert not messages and len(errors)==1 and len(calls)==1
    else:
        assert state.records is snap.records and state.idx is snap.index and state.lib_loaded
        assert not messages and len(errors)==1 and not calls;c.release(lease)
