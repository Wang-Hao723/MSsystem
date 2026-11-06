import copy

import pytest

from eval.v16_evidence_validation import record_store_runtime_smoke as smoke


def healthy():
    return {"status": "ok", "project_root": str(smoke.ROOT), "service_version": "16.0.0",
        "active_port": 8516, "storage": {"library_write_enabled": False},
        "runtime_manifest": {"fresh": True, "release_ready": True, "restart_required": False},
        "library_loaded": True, "library_spectra": 3153872}


def test_live_smoke_accepts_only_the_expected_readonly_full_library():
    smoke.check_health(healthy(), 3153872, loaded=True)


@pytest.mark.parametrize("change", ["root", "version", "write_enabled", "stale", "count"])
def test_live_smoke_rejects_wrong_instance_staleness_or_library(change):
    payload = copy.deepcopy(healthy())
    if change == "root":
        payload["project_root"] = "/unrelated/project"
    elif change == "version":
        payload["service_version"] = "15.0.0"
    elif change == "write_enabled":
        payload["storage"]["library_write_enabled"] = True
    elif change == "stale":
        payload["runtime_manifest"]["fresh"] = False
    else:
        payload["library_spectra"] = 1
    with pytest.raises(ValueError):
        smoke.check_health(payload, 3153872, loaded=True)
