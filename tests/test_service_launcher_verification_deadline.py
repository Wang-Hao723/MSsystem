"""Exercise the real launcher argument parser and shell-to-main dispatch."""
import importlib.util
import json
import os
from pathlib import Path
import shutil
import socket
import subprocess
import sys
import pytest

ROOT=Path(__file__).resolve().parents[1]


def launcher():
    spec=importlib.util.spec_from_file_location('tested_service_launcher',ROOT/'scripts/service_launcher.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module


@pytest.mark.parametrize('wait',[600,3600,1,30,300])
@pytest.mark.parametrize('mode',['fastapi','streamlit'])
def test_real_main_dispatches_finite_deadline(tmp_path,monkeypatch,wait,mode):
    m=launcher();calls=[]
    def start(*args):calls.append(args);return 0
    monkeypatch.setattr(m,'start_service',start)
    monkeypatch.setattr(sys,'argv',['service_launcher.py','--project-root',str(tmp_path),'--wait-seconds',str(wait),'--mode',mode])
    assert m.main()==0;assert len(calls)==1
    assert calls[0][0]==tmp_path and calls[0][-1]==float(wait) and calls[0][-2]==mode


@pytest.mark.parametrize('wait',['3601','NaN','inf','-inf','0','-1','3600.1'])
def test_real_main_rejects_invalid_wait_before_spawn(tmp_path,monkeypatch,wait):
    m=launcher();calls=[];monkeypatch.setattr(m,'start_service',lambda *a:calls.append(a))
    monkeypatch.setattr(sys,'argv',['service_launcher.py','--project-root',str(tmp_path),'--wait-seconds='+wait])
    with pytest.raises(SystemExit) as error:m.main()
    assert error.value.code==2 and calls==[]


@pytest.mark.parametrize('override,expected',[(None,1800),('0008',8),('3600',3600)])
def test_actual_bash_entrypoint_reaches_real_main_parser(tmp_path,override,expected):
    checkout=tmp_path/'checkout';checkout.mkdir();shutil.copyfile(ROOT/'run.sh',checkout/'run.sh')
    recorded=tmp_path/'dispatched.json';fake=tmp_path/'python'
    fake.write_text('#!'+sys.executable+'\n'+'''import importlib.util,json,os,pathlib,sys
if len(sys.argv)>1 and sys.argv[1].endswith('/scripts/service_launcher.py'):
    spec=importlib.util.spec_from_file_location('real_service_launcher',os.environ['TEST_REAL_LAUNCHER'])
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
    def dispatch(*args):
        pathlib.Path(os.environ['TEST_DISPATCH']).write_text(json.dumps(args,default=str))
        return 0
    m.start_service=dispatch
    sys.argv=sys.argv[1:]
    raise SystemExit(m.main())
''');fake.chmod(0o755)
    env=os.environ.copy();env.pop('MSSYSTEM_START_WAIT_SECONDS',None)
    env.update(MSSYSTEM_PYTHON=str(fake),TEST_REAL_LAUNCHER=str(ROOT/'scripts/service_launcher.py'),TEST_DISPATCH=str(recorded),MSSYSTEM_RUNTIME_ROOT=str(tmp_path/'runtime'),MSSYSTEM_VENDOR_CONVERSION_CACHE_ROOT=str(tmp_path/'vendor'),MPLCONFIGDIR=str(tmp_path/'mpl'),NUMBA_CACHE_DIR=str(tmp_path/'numba'))
    if override is not None:env['MSSYSTEM_START_WAIT_SECONDS']=override
    with socket.socket() as sock:sock.bind(('127.0.0.1',0));port=sock.getsockname()[1]
    result=subprocess.run(['bash',str(checkout/'run.sh'),'--no-open','--port',str(port)],env=env,text=True,capture_output=True,timeout=20)
    assert result.returncode==0,result.stdout+result.stderr
    args=json.loads(recorded.read_text());assert args[-1]==expected and args[-2]=='fastapi'
