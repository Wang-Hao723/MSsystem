from __future__ import annotations

import json
import sys
from copy import deepcopy
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


@pytest.fixture
def qc_policy():
    return json.loads(
        (
            ROOT
            / "config/feature_fragment_benchmark_v3_1/spectral_qc_policy.json"
        ).read_text(encoding="utf-8")
    )


@pytest.fixture
def v31_record_factory():
    def factory(spectrum_id: str = "test-spectrum"):
        return {
            "schema_version": "3.0.0",
            "spectrum_id": spectrum_id,
            "evidence_stage": "raw_msms_spectrum",
            "identification_scope": {
                "standard_confirmed_by_mssystem": False,
            },
            "source": {
                "repository": "Repository A",
                "source_id": "source-a",
                "origin_dataset_id": "dataset-a",
                "laboratory_id": "lab-a",
                "usi": None,
            },
            "compound": {
                "display_name": "Test molecule",
                "connectivity": "AAAAAAAAAAAAAA",
                "formula": "C10H10O2",
            },
            "acquisition": {
                "instrument_domain_id": "domain-a",
                "instrument_domain_resolution": {
                    "status": "resolved",
                    "exact_condition_metadata_complete": True,
                },
                "activation_location": "collision_cell",
                "collision_energy": {
                    "status": "resolved",
                    "mode": "fixed",
                    "unit": "eV",
                    "values": [20.0],
                    "weights": [1.0],
                    "conversion_applied": False,
                },
            },
            "raw_acquisition": {
                "raw_file_id": None,
                "scan_id": None,
                "injection_id": None,
            },
            "independence": {
                "two_dimensional_identity_id": "AAAAAAAAAAAAAA",
                "formula_family_id": "C10H10O2",
                "source_origin_id": "dataset-a",
                "laboratory_id": "lab-a",
                "exact_peak_hash": "a" * 64,
                "mirror_cluster_id": "exact:" + "a" * 64,
                "acquisition_replication_group_id": "acq:test",
                "peak_hash_used_in_acquisition_replication_group": False,
            },
            "peaks": [[50.0, 1.0], [75.0, 0.2], [100.0, 0.009]],
            "quantitative_eligibility": {
                "primary_license_eligible": True,
                "positive_lr_eligible": True,
                "negative_lr_eligible": True,
                "positive_only_composite_or_in_source": False,
                "ineligible_reason_codes": [],
            },
        }

    return factory
