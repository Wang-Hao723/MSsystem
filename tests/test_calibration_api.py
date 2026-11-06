from fastapi.testclient import TestClient

from web.server import app


def test_calibration_policy_is_explicitly_fail_closed():
    response = TestClient(app).get("/api/calibration/policy")
    assert response.status_code == 200
    payload = response.json()
    assert payload["production_calibration_available"] is False
    assert payload["production_calibration_id"] is None
    assert payload["current_legacy_high_is_qvalue_tier"] is False
    assert payload["policy"]["raw_similarity_is_probability"] is False
    assert payload["policy"]["qvalue_is_probability"] is False
