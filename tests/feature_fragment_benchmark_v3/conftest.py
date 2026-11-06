from __future__ import annotations

from copy import deepcopy

import pytest

from matching.feature_fragment_benchmark_v3.parsing import (
    parse_collision_energy_program,
    parse_scan_range,
)
from matching.feature_fragment_benchmark_v3.provenance import normalize_public_spectrum


@pytest.fixture
def v3_domain():
    return {
        "domain_id": "test_vendor_model_hcd_v3",
        "family_id": "test_family",
        "vendor": "Test Vendor",
        "models": ["Test Model"],
        "analyzer": "Orbitrap",
        "activation": ["HCD"],
        "activation_location": "HCD_cell",
        "collision_cell_regions": ["HCD_cell"],
        "ce_units": ["NCE", "eV"],
        "ce_modes": ["fixed", "stepped", "MSE"],
    }


@pytest.fixture
def v3_license_policy():
    return {
        "label_decisions": {
            "CC0": {
                "tier": "primary",
                "license_url": "https://creativecommons.org/publicdomain/zero/1.0/",
                "commercial_use_allowed": True,
                "share_alike_review_complete": True,
                "reason": "test",
            },
            "CC BY-SA": {
                "tier": "quarantine",
                "license_url": "https://creativecommons.org/licenses/by-sa/4.0/",
                "commercial_use_allowed": True,
                "share_alike_review_complete": False,
                "reason": "test",
            },
        }
    }


@pytest.fixture
def v3_record_factory(v3_domain, v3_license_policy):
    def factory(
        spectrum_id: str,
        *,
        connectivity: str = "ABCDEFGHIJKLMN",
        formula: str = "C10H20O2",
        source_id: str = "source-a",
        laboratory_id: str = "lab-a",
        peaks=None,
        isolation_status: str = "verified",
        collision_energy: str = "30 NCE",
        license_label: str = "CC0",
    ):
        record = {
            "spectrum_id": spectrum_id,
            "source": {
                "repository": "Public Test Repository",
                "repository_release": "1.0",
                "source_id": source_id,
                "origin_dataset_id": source_id,
                "origin_study_id": source_id + "-study",
                "laboratory_id": laboratory_id,
                "accession": spectrum_id,
                "usi": "mzspec:test:" + spectrum_id,
                "record_url": "https://example.test/" + spectrum_id,
                "query_url": "https://example.test/query",
                "license_label": license_label,
                "attribution": "Test Author; " + spectrum_id,
                "rights_holder": ["Test Author"],
                "source_file_sha256": "a" * 64,
                "record_canonical_json_sha256": "b" * 64,
            },
            "compound": {
                "display_name": "Test compound " + spectrum_id,
                "inchikey": connectivity + "-ABCDEFGHIJ-A",
                "formula": formula,
                "identity_basis": "repository_reported_structure",
                "reference_material_claim": "not_reported",
                "standard_confirmation_status": "not_confirmed_by_mssystem",
                "chemical_class": {
                    "ontology": "ChemOnt",
                    "ontology_id": "CHEMONTID:1",
                    "leaf_class": "Test class",
                    "path": ["Organic compounds", "Test class"],
                },
            },
            "acquisition": {
                "instrument_vendor": "Test Vendor",
                "instrument_model": "Test Model",
                "analyzer": "Orbitrap",
                "activation": "HCD",
                "activation_location": "HCD_cell",
                "collision_cell_region": "HCD_cell",
                "collision_energy": parse_collision_energy_program(collision_energy),
                "scan_range": parse_scan_range("50-1000"),
                "polarity": "positive",
                "adduct": "[M+H]+",
                "charge": 1,
                "precursor_mz": 300.0,
                "precursor_isolation": {
                    "status": isolation_status,
                    "width_mz": 1.0 if isolation_status == "verified" else None,
                },
            },
            "raw_acquisition": {
                "raw_file_id": source_id + ".mzML",
                "scan_id": spectrum_id,
                "injection_id": source_id + "-injection",
                "acquisition_event_id": spectrum_id,
            },
            "peaks": deepcopy(
                peaks
                or [[50.0, 100.0], [100.0, 50.0], [150.0, 25.0]]
            ),
        }
        return normalize_public_spectrum(
            record,
            license_policy=v3_license_policy,
            instrument_domains=[v3_domain],
        )

    return factory
