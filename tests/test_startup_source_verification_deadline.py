"""The shell entrypoint must forward a sufficient finite cold-start deadline."""
import json
import os
from pathlib import Path
import shutil
import socket
import subprocess
import sys

import pytest


def _run_entrypoint(tmp_path, override):
    project = tmp_path / "checkout"
    project.mkdir()
    shutil.copyfile(Path(__file__).resolve().parents[1] / "run.sh", project / "run.sh")
    recorded = tmp_path / "launcher-arguments.json"
    fake_python = tmp_path / "python-with-recorded-launcher"
    fake_python.write_text(
        "#!" + sys.executable + "\n"
        "import json, os, pathlib, sys\n"
        "if len(sys.argv)>1 and sys.argv[1].endswith('/scripts/service_launcher.py'):\n"
        "    pathlib.Path(os.environ['TEST_LAUNCHER_ARGUMENTS']).write_text(json.dumps(sys.argv[2:]))\n"
    )
    fake_python.chmod(0o755)
    env = os.environ.copy()
    env.pop("MSSYSTEM_START_WAIT_SECONDS", None)
    env.update(MSSYSTEM_PYTHON=str(fake_python), TEST_LAUNCHER_ARGUMENTS=str(recorded),
               MSSYSTEM_RUNTIME_ROOT=str(tmp_path / "runtime"),
               MSSYSTEM_VENDOR_CONVERSION_CACHE_ROOT=str(tmp_path / "vendor"),
               MPLCONFIGDIR=str(tmp_path / "mpl"), NUMBA_CACHE_DIR=str(tmp_path / "numba"))
    if override is not None:
        env["MSSYSTEM_START_WAIT_SECONDS"] = override
    with socket.socket() as listener:
        listener.bind(("127.0.0.1", 0))
        port = listener.getsockname()[1]
    result = subprocess.run(["bash", str(project / "run.sh"), "--no-open", "--port", str(port)],
                            env=env, text=True, capture_output=True, timeout=20)
    return result, json.loads(recorded.read_text()) if recorded.exists() else None


@pytest.mark.parametrize("override,expected", [(None, "1800"), ("75", "75"), ("0008", "8"), ("3600", "3600")])
def test_finite_cold_start_deadline_reaches_actual_launcher(tmp_path, override, expected):
    result, arguments = _run_entrypoint(tmp_path, override)
    assert result.returncode == 0, result.stdout + result.stderr
    assert arguments[arguments.index("--wait-seconds") + 1] == expected


@pytest.mark.parametrize("override", ["0", "-1", "nan", "1.5", "3601", "999999999999999999"])
def test_invalid_deadline_stops_before_launch(tmp_path, override):
    result, arguments = _run_entrypoint(tmp_path, override)
    assert result.returncode == 2
    assert "must be an integer from 1 to 3600" in result.stderr
    assert arguments is None
