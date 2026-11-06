"""Completion checks for required formal outputs, independent of match scores."""
from collections.abc import Mapping


def formal_output_check(summaries: Mapping, *, primary_artifact_current: bool) -> dict:
    required = ("primary_identity", "identification_records", "claim_contract", "compound_conclusions")
    failures = []
    for name in required:
        value = summaries.get(name)
        if not isinstance(value, Mapping):
            failures.append({"output": name, "error": "required output summary is absent"})
        elif value.get("status") != "available" or value.get("error"):
            failures.append({"output": name, "error": str(value.get("error") or value.get("status") or "unavailable")})
    value = summaries.get("primary_identity")
    primary = value if isinstance(value, Mapping) else {}
    if primary.get("production_identity_selector") is not True:
        failures.append({"output": "primary_identity", "error": "formal production identity selector is absent"})
    if not primary_artifact_current:
        failures.append({"output": "primary_identity", "error": "active formal identity artifact is absent or stale"})
    return {"status": "failed" if failures else "passed", "required_outputs": list(required),
            "primary_artifact_current": bool(primary_artifact_current), "failures": failures,
            "audit_artifact_substitution_allowed": False}


def assert_formal_outputs(check: Mapping) -> None:
    if check.get("status") != "passed":
        errors = "; ".join(f"{item['output']}: {item['error']}" for item in check.get("failures") or [])
        raise RuntimeError("Required formal result finalization failed: " + errors)
