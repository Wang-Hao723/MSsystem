import copy
import hashlib
import json
from pathlib import Path

import pytest

from matching.feature_fragments.catalog import load_catalog
from matching.feature_fragments.validation import (
    FixtureValidationError,
    validate_fixture,
    validate_public_pilot,
)
from scripts.build_feature_fragment_massbank_fixture import CLASS_ACCESSIONS


ROOT = Path(__file__).resolve().parents[1]
FIXTURE = (
    ROOT
    / "tests"
    / "fixtures"
    / "feature_fragments"
    / "massbank_antwerp_phospholipids_v1.json"
)


def _fixture():
    return json.loads(FIXTURE.read_text(encoding="utf-8"))


def _rehash(fixture):
    payload = dict(fixture)
    payload.pop("content_sha256", None)
    raw = json.dumps(
        payload, ensure_ascii=False, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")
    fixture["content_sha256"] = hashlib.sha256(raw).hexdigest()


def test_explicit_massbank_allowlist_and_fixture_are_complete():
    assert sum(len(accessions) for accessions in CLASS_ACCESSIONS.values()) == 30
    assert len({a for values in CLASS_ACCESSIONS.values() for a in values}) == 30
    summary = validate_fixture(_fixture())
    assert summary == {
        "record_count": 30,
        "independent_2d_identity_count": 10,
        "normalized_formula_family_count": 9,
        "source_record_count": 30,
        "acquisition_aware_hash_count": 30,
        "peak_only_mirror_duplicate_count": 0,
    }


def test_public_pilot_reproduces_counts_and_formula_isolated_metrics():
    receipt = validate_public_pilot(load_catalog(), _fixture())

    assert receipt["mode"] == "shadow_only"
    assert receipt["runtime_influence"] is False
    assert receipt["identity_rank_influence"] is False
    assert receipt["instrument_family_source_audit"]["passed"] is True
    assert len(receipt["rule_results"]) == 6
    assert all(result["configured_counts_match"] for result in receipt["rule_results"])
    assert all(
        result["heldout_validation"]["identity_overlap_violations"] == 0
        and result["heldout_validation"]["formula_family_overlap_violations"] == 0
        for result in receipt["rule_results"]
    )
    pc = [
        result
        for result in receipt["rule_results"]
        if result["independent_fragmentation_group"] == "phosphocholine_headgroup_184"
    ]
    assert {result["heldout_validation"]["roc_auc"] for result in pc} == {1.0}
    pe40 = next(
        result
        for result in receipt["rule_results"]
        if result["rule_id"] == "FF-PE-NL141-POS-H-AG6560-40EV"
    )
    assert pe40["observed_counts"]["target_present"] == 2
    assert pe40["observed_counts"]["target_total"] == 3
    assert pe40["heldout_validation"]["roc_auc"] == pytest.approx(2 / 3)
    assert receipt["claims"]["cross_platform_rule_validated"] is False
    assert receipt["claims"]["exact_identity_validated"] is False
    assert receipt["claims"]["isomer_discrimination_validated"] is False


def test_formula_sharing_regioisomers_are_held_out_together():
    receipt = validate_public_pilot(load_catalog(), _fixture())
    pc20 = next(
        result
        for result in receipt["rule_results"]
        if result["rule_id"] == "FF-PC-184-POS-H-AG6560-20EV"
    )
    same_formula = [
        prediction
        for prediction in pc20["heldout_validation"]["predictions"]
        if prediction["formula_family"] == "C42H82NO8P"
    ]
    assert len(same_formula) == 2
    assert len({row["connectivity"] for row in same_formula}) == 2
    assert {row["training_target_identities"] for row in same_formula} == {4}
    assert all(row["target"] for row in same_formula)


def test_fixture_rejects_identifier_name_and_hash_tampering():
    identifier_name = _fixture()
    identifier_name["spectra"][0]["compound"]["name"] = "MSBNK-TEST-123"
    _rehash(identifier_name)
    with pytest.raises(FixtureValidationError, match="identifier used as name"):
        validate_fixture(identifier_name)

    tampered = _fixture()
    tampered["spectra"][0]["peaks"][0][1] *= 2
    with pytest.raises(FixtureValidationError, match="content_sha256 mismatch"):
        validate_fixture(tampered)


def test_validation_rejects_catalog_counts_not_reproduced_by_public_records():
    catalog = load_catalog()
    catalog = copy.deepcopy(catalog)
    catalog["rules"][0]["occurrence_statistics"] = {
        "present": 5,
        "total": 6,
        "frequency": 5 / 6,
        "counting_unit": "independent_2d_identity",
    }
    with pytest.raises(FixtureValidationError, match="catalog counts"):
        validate_public_pilot(catalog, _fixture())
