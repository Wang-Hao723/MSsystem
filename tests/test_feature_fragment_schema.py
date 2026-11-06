import copy

import pytest

from matching.feature_fragments.catalog import load_catalog
from matching.feature_fragments.schema import (
    CatalogValidationError,
    stable_compound_name,
    validate_catalog,
)


def test_public_pilot_catalog_has_atomic_audit_fields_and_is_shadow_only():
    catalog = load_catalog()

    assert len(catalog["rules"]) == 6
    assert {domain["vendor"] for domain in catalog["instrument_domains"]} >= {
        "Thermo Fisher Scientific",
        "Agilent Technologies",
        "SCIEX",
        "Bruker Daltonics",
        "Waters Corporation",
    }
    for rule in catalog["rules"]:
        assert rule["runtime_influence"] is False
        assert rule["identity_rank_influence"] is False
        assert rule["source_ids"]
        assert rule["chemical_classes"]
        assert rule["adducts"]
        assert rule["charges"] == [1]
        assert rule["activation"] == "beam_cid"
        assert rule["collision_energy"]["unit"] == "eV"
        assert rule["mass_tolerance"]["absolute_da"] > 0
        assert rule["mass_tolerance"]["ppm"] > 0
        assert rule["independent_fragmentation_group"]
        assert rule["positive_evidence_limitations"]
        assert rule["negative_evidence"]["limitations"]
        assert rule["output_semantics"] == {
            "evidence_stage": "msms_rule_support",
            "exact_identity_confirmed": False,
            "standard_confirmed": False,
            "isomer_discriminating": False,
        }


def test_catalog_rejects_unknown_source_and_cross_unit_conversion():
    catalog = load_catalog()
    bad_source = copy.deepcopy(catalog)
    bad_source["rules"][0]["source_ids"].append("missing-source")
    with pytest.raises(CatalogValidationError, match="unknown"):
        validate_catalog(bad_source)

    bad_domain = copy.deepcopy(catalog)
    bad_domain["instrument_domains"][0]["cross_unit_conversion_allowed"] = True
    with pytest.raises(CatalogValidationError, match="must be false"):
        validate_catalog(bad_domain)


def test_catalog_rejects_evidence_stage_upgrade_and_fixed_energy_range():
    catalog = load_catalog()
    upgraded = copy.deepcopy(catalog)
    upgraded["rules"][0]["output_semantics"]["exact_identity_confirmed"] = True
    with pytest.raises(CatalogValidationError, match="forbidden evidence upgrade"):
        validate_catalog(upgraded)

    ranged = copy.deepcopy(catalog)
    ranged["rules"][0]["collision_energy"]["maximum"] = 20.0
    with pytest.raises(CatalogValidationError, match="fixed mode"):
        validate_catalog(ranged)


def test_database_accession_cannot_be_used_as_compound_name():
    assert stable_compound_name("1,2-Dioleoyl-sn-glycero-3-phosphocholine")
    assert not stable_compound_name("MSBNK-Antwerp_Univ-METOX_N107506_EF88")
    assert not stable_compound_name("HMDB0000123")
    assert not stable_compound_name("CCMSLIB000001")
