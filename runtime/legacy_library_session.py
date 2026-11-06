"""Process-wide library generations and analysis leases for Streamlit sessions.

This imported module survives script reruns. Library replacement and reservation
use the same lock, and every running or not-yet-started worker holds a lease.
"""
from dataclasses import dataclass
import importlib
import threading
import uuid


class LibraryBusyError(RuntimeError):
    pass


def _runtime_binding():
    module=importlib.import_module('matching.search')
    return tuple(getattr(module,name,None) for name in
                 ('_lib_records','_mz_index','_nm_index','_lib_db_path','_lib_source_families'))


@dataclass(frozen=True)
class LibrarySnapshot:
    generation: str
    records: object
    index: object
    runtime_binding: tuple


@dataclass(frozen=True)
class AnalysisLease:
    token: str
    task_id: str
    snapshot: LibrarySnapshot


class LibraryCoordinator:
    def __init__(self, runtime_binding=_runtime_binding):
        self.lock=threading.RLock()
        self._runtime_binding=runtime_binding
        self._current=None
        self._leases={}
        self._tasks={}

    def replace(self, loader, prepare):
        """Reject active workers before clearing any caller or shared state."""
        with self.lock:
            if self._leases:
                raise LibraryBusyError('An analysis is queued or running in this process. Wait before reloading the library.')
            self._current=None
            prepare()
            records,index=loader()
            if records is None or index is None:
                raise ValueError("Library loader returned no records or index.")
            self._current=LibrarySnapshot(uuid.uuid4().hex,records,index,self._runtime_binding())
            return self._current

    def reserve(self, task_id, records, index):
        with self.lock:
            task_id=str(task_id)
            snapshot=self._current
            if snapshot is None or records is not snapshot.records or index is not snapshot.index:
                raise LibraryBusyError('This session has no current library generation. Load the library before starting analysis.')
            if task_id in self._tasks:
                raise LibraryBusyError('This task is already queued or running in another thread or session.')
            self._validate_binding(snapshot)
            lease=AnalysisLease(uuid.uuid4().hex,str(task_id),snapshot)
            self._leases[lease.token]=lease
            self._tasks[lease.task_id]=lease.token
            return lease

    def _validate_binding(self,snapshot):
        current=self._runtime_binding()
        if len(current)!=len(snapshot.runtime_binding) or any(a is not b for a,b in zip(current,snapshot.runtime_binding)):
            raise LibraryBusyError('The library changed outside the shared coordinator. Stop and reload it before analysis.')

    def validate(self,lease):
        with self.lock:
            if self._leases.get(lease.token) is not lease or self._current is not lease.snapshot:
                raise LibraryBusyError('The analysis library lease is no longer current.')
            self._validate_binding(lease.snapshot)

    def release(self,lease):
        with self.lock:
            if self._leases.get(lease.token) is lease:
                del self._leases[lease.token]
                self._tasks.pop(lease.task_id,None)

    def assert_task_idle(self,task_id):
        with self.lock:
            if str(task_id) in self._tasks:
                raise LibraryBusyError('This task is queued or running in a thread or session. Wait before deleting it.')

    def active_count(self):
        with self.lock:return len(self._leases)


COORDINATOR=LibraryCoordinator()


def capture_session_state():
    # A background worker receives the actual session object, not the global
    # proxy which resolves against the calling thread's Streamlit context.
    from streamlit.runtime.scriptrunner import get_script_run_ctx
    context=get_script_run_ctx()
    if context is None:
        raise RuntimeError('A Streamlit session is required to start analysis.')
    return context.session_state


class TaskRegistryChangedError(RuntimeError):
    pass


def _file_revision(path):
    import hashlib
    try:return hashlib.sha256(path.read_bytes()).hexdigest()
    except FileNotFoundError:return None


def atomic_json_write(path,payload):
    """Write a complete JSON file, preserving its previous bytes on write failure."""
    import json,os,tempfile
    from pathlib import Path
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True);temporary=None
    encoded=json.dumps(payload,ensure_ascii=False,indent=2,default=str).encode('utf-8')
    try:
        with tempfile.NamedTemporaryFile(mode='wb',dir=path.parent,prefix='.'+path.name+'.',delete=False) as handle:
            temporary=Path(handle.name);handle.write(encoded);handle.flush();os.fsync(handle.fileno())
        os.replace(temporary,path)
    finally:
        if temporary is not None:temporary.unlink(missing_ok=True)
    import hashlib
    return hashlib.sha256(encoded).hexdigest()


class TaskRegistry:
    """One task mapping for every Streamlit session in this process.

    A changed external file is rejected before saving. The revision check and
    replace are not a transaction with non-cooperating external processes.
    """
    def __init__(self,path,lock):
        from pathlib import Path
        import json
        self.path=Path(path).resolve();self.lock=lock
        with lock:
            before=_file_revision(self.path)
            payload=json.loads(self.path.read_text()) if before is not None else {}
            if not isinstance(payload,dict) or any(not isinstance(row,dict) for row in payload.values()):
                raise ValueError('The task index must contain task records indexed by ID.')
            if before!=_file_revision(self.path):raise TaskRegistryChangedError('Task index changed during loading. Reload through the active application.')
            self.tasks=payload;self.revision=before

    def assert_current(self,tasks=None):
        with self.lock:
            if tasks is not None and tasks is not self.tasks:
                raise TaskRegistryChangedError('This session has an outdated task mapping. Refresh it before changing tasks.')
            if _file_revision(self.path)!=self.revision:
                raise TaskRegistryChangedError('Task index was changed by another writer. No stale task snapshot was saved; use the active application or restart this instance.')

    def attach(self,session_state):
        with self.lock:
            self.assert_current();session_state.tasks=self.tasks;session_state.tasks_loaded_from_disk=True

    def save(self,tasks):
        with self.lock:
            self.assert_current(tasks)
            payload={tid:{key:value for key,value in task.items() if key!='results'} for tid,task in self.tasks.items()}
            self.revision=atomic_json_write(self.path,payload)


_TASK_REGISTRIES={}


def task_registry(path):
    from pathlib import Path
    key=str(Path(path).resolve())
    with COORDINATOR.lock:
        if key not in _TASK_REGISTRIES:_TASK_REGISTRIES[key]=TaskRegistry(path,COORDINATOR.lock)
        return _TASK_REGISTRIES[key]
