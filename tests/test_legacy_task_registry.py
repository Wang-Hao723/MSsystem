"""Shared sessions write actual temporary JSON without dropping newer tasks."""
import ast
import copy
import json
import os
import importlib.util
from pathlib import Path
from types import SimpleNamespace
import pytest
from runtime import legacy_library_session as M
_spec = importlib.util.spec_from_file_location(
    'legacy_concurrency_test_helpers', Path(__file__).with_name('test_legacy_library_concurrency.py'))
_helpers = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_helpers)
coordinator, publish, session = _helpers.coordinator, _helpers.publish, _helpers.session
app_scope, APP = _helpers.app_scope, _helpers.APP


def registry(tmp_path,c):
    path=tmp_path/'tasks.json';path.write_text(json.dumps({'old':{'id':'old','status':'completed','progress':100,'extra':'preserved'}}))
    return M.TaskRegistry(path,c.lock)


def test_two_sessions_share_new_records_progress_and_completion_without_stale_overwrite(tmp_path):
    c=coordinator();r=registry(tmp_path,c)
    a=SimpleNamespace(tasks={'staleA':{}},tasks_loaded_from_disk=False);b=SimpleNamespace(tasks={'staleB':{}},tasks_loaded_from_disk=False)
    r.attach(a);r.attach(b);assert a.tasks is b.tasks is r.tasks
    with c.lock:a.tasks['a']={'id':'a','status':'queued','progress':0};r.save(a.tasks)
    with c.lock:b.tasks['b']={'id':'b','status':'queued','progress':0};r.save(b.tasks)
    with c.lock:a.tasks['a'].update(status='running',progress=45);r.save(a.tasks)
    with c.lock:b.tasks['b'].update(status='completed',progress=100);r.save(b.tasks)
    with c.lock:a.tasks['a'].update(status='completed',progress=100);r.save(a.tasks)
    observed=json.loads(r.path.read_text());assert set(observed)=={'old','a','b'}
    assert observed['old']['extra']=='preserved' and all(observed[k]['status']=='completed' and observed[k]['progress']==100 for k in observed)
    r.attach(a);assert a.tasks is r.tasks


def test_stale_session_copy_is_rejected_before_replacing_disk(tmp_path):
    c=coordinator();r=registry(tmp_path,c);before=r.path.read_bytes()
    with pytest.raises(M.TaskRegistryChangedError,match='outdated task mapping'):r.save(copy.deepcopy(r.tasks))
    assert r.path.read_bytes()==before


def test_external_task_writer_is_detected_without_overwriting_its_bytes(tmp_path):
    c=coordinator();r=registry(tmp_path,c);outside={'outside':{'status':'running'}};r.path.write_text(json.dumps(outside));before=r.path.read_bytes()
    with pytest.raises(M.TaskRegistryChangedError,match='another writer'):r.save(r.tasks)
    with pytest.raises(M.TaskRegistryChangedError):r.attach(SimpleNamespace())
    assert r.path.read_bytes()==before


def test_atomic_index_replace_failure_keeps_previous_json_and_raises(tmp_path,monkeypatch):
    c=coordinator();r=registry(tmp_path,c);before=r.path.read_bytes();r.tasks['new']={'status':'queued'}
    def fail(*args):raise PermissionError('injected replace denial')
    with monkeypatch.context() as context:
        context.setattr(os,'replace',fail)
        with pytest.raises(PermissionError,match='replace denial'):r.save(r.tasks)
    assert r.path.read_bytes()==before and not list(tmp_path.glob('.tasks.json.*'))
    r.save(r.tasks);assert 'new' in json.loads(r.path.read_text())


def persistent_scope(tmp_path):
    c=coordinator();snap=publish(c);state=session(snap);r=registry(tmp_path,c);r.attach(state)
    state.tasks['task']={'id':'task','status':'pending','file_path':str(tmp_path/'upload.mzML'),'params':{'enable_ion_identity_networking':False,'enable_chimeric_deconvolution':False}}
    Path(state.tasks['task']['file_path']).write_text('original upload');r.save(state.tasks)
    scope,errors,saved=app_scope(c,state);scope.update(_task_registry=lambda path:r,_atomic_json_write=M.atomic_json_write,TASKS_FILE=r.path,TASK_DIR=tmp_path)
    names={'_load_tasks_from_disk','_save_tasks_to_disk','_save_task_results','_delete_task_from_session'}
    nodes=[n for n in ast.parse(APP.read_text()).body if isinstance(n,ast.FunctionDef) and n.name in names]
    assert len(nodes)==4;exec(compile(ast.Module(body=nodes,type_ignores=[]),str(APP),'exec'),scope)
    return c,snap,state,r,scope,errors


def test_real_worker_result_disk_failure_is_failed_not_completed(tmp_path,monkeypatch):
    c,snap,state,r,scope,errors=persistent_scope(tmp_path);result=tmp_path/'task_results.json';result.write_text('["previous result"]');before=result.read_bytes();replace=os.replace
    def fail_result(source,dest):
        if Path(dest)==result:raise PermissionError('result replace denied')
        return replace(source,dest)
    lease=c.reserve('task',snap.records,snap.index);state.analysis_running=True;state.active_task_id='task'
    with monkeypatch.context() as context:
        context.setattr(os,'replace',fail_result);scope['_run_analysis_in_thread']('task',state,lease)
    assert result.read_bytes()==before and state.tasks['task']['status']=='failed' and c.active_count()==0
    assert json.loads(r.path.read_text())['task']['status']=='failed'


def test_real_save_wrapper_propagates_external_generation_conflict(tmp_path):
    c,snap,state,r,scope,errors=persistent_scope(tmp_path);r.path.write_text('{"new_writer":{"status":"running"}}');before=r.path.read_bytes()
    with pytest.raises(M.TaskRegistryChangedError):scope['_save_tasks_to_disk'](state)
    assert r.path.read_bytes()==before


@pytest.mark.parametrize('stored_status',['pending','running','completed'])
def test_delete_rejects_other_session_active_lease_even_with_stale_status(tmp_path,stored_status):
    c,snap,state,r,scope,errors=persistent_scope(tmp_path);state.tasks['task']['status']=stored_status;r.save(state.tasks)
    result=tmp_path/'task_results.json';result.write_text('[]');before=r.path.read_bytes();upload=Path(state.tasks['task']['file_path']);lease=c.reserve('task',snap.records,snap.index)
    assert scope['_delete_task_from_session']('task') is False
    assert 'task' in state.tasks and upload.read_text()=='original upload' and result.read_text()=='[]' and r.path.read_bytes()==before and errors;c.release(lease)


def test_idle_delete_updates_shared_mapping_and_removes_only_selected_files(tmp_path):
    c,snap,state,r,scope,errors=persistent_scope(tmp_path);other=SimpleNamespace();r.attach(other);upload=Path(state.tasks['task']['file_path']);result=tmp_path/'task_results.json';result.write_text('[]')
    assert scope['_delete_task_from_session']('task') is True
    assert 'task' not in state.tasks and 'task' not in other.tasks and 'old' in other.tasks
    assert not upload.exists() and not result.exists() and set(json.loads(r.path.read_text()))=={'old'} and not errors


def test_delete_index_write_failure_restores_mapping_and_preserves_files(tmp_path,monkeypatch):
    c,snap,state,r,scope,errors=persistent_scope(tmp_path);upload=Path(state.tasks['task']['file_path']);result=tmp_path/'task_results.json';result.write_text('[]');before=r.path.read_bytes();old=copy.deepcopy(state.tasks)
    def fail(*a):raise PermissionError('index replace denied')
    with monkeypatch.context() as context:
        context.setattr(os,'replace',fail);assert scope['_delete_task_from_session']('task') is False
    assert state.tasks==old and upload.exists() and result.exists() and r.path.read_bytes()==before and errors


def test_external_index_change_blocks_start_before_params_status_or_lease_mutation(tmp_path):
    c,snap,state,r,scope,errors=persistent_scope(tmp_path);old=copy.deepcopy(state.tasks['task']);r.path.write_text('{"outside":{"status":"running"}}');before=r.path.read_bytes()
    assert scope['_start_analysis_in_thread']('task',params={'cosine_threshold':.99}) is False
    assert state.tasks['task']==old and c.active_count()==0 and not state.analysis_running and errors
    assert r.path.read_bytes()==before
