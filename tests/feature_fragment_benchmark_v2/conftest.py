from __future__ import annotations

from copy import deepcopy
from typing import Any

import pytest

from matching.feature_fragment_benchmark.mirror import exact_peak_hash


_DEFAULT_PEAKS = [[50.0, 1.0], [100.0, 0.5], [200.0, 0.25]]


@pytest.fixture
def benchmark_record_factory():
    def make_record(
        spectrum_id: str,
        *,
        identity: str | None = None,
        formula: str = "C10H20O2",
        chemical_class: str = "Target class",
        peaks: list[list[float]] | None = None,
        mirror_cluster_id: str | None = None,
        source_id: str = "source-a",
        laboratory_id: str = "lab-a",
        research_id: str = "study-a",
        instrument_domain_id: str = "agilent_qtof_cid_contributor_text",
        activation: str = "CID",
        collision_cell_region: str = "collision_cell",
        ce_unit: str = "eV",
        ce_mode: str = "fixed",
        ce_value: float | None = 20.0,
        ce_minimum: float = 20.0,
        ce_maximum: float = 20.0,
        adduct: str = "[M+H]+",
        polarity: str = "positive",
        charge: int = 1,
        precursor_mz: float = 500.0,
        reported_mass_range: dict[str, float] | None = None,
        quality_overrides: dict[str, Any] | None = None,
    ) -> dict[str, Any]:
        identity = identity or f"identity-{spectrum_id}"
        normalized_peaks = deepcopy(peaks if peaks is not None else _DEFAULT_PEAKS)
        peak_hash = exact_peak_hash(normalized_peaks)
        if reported_mass_range is None:
            reported_mass_range = {"minimum": 40.0, "maximum": precursor_mz}
        quality = {
            "centroid_peak_list": True,
            "positive_finite_peaks": True,
            "at_least_three_peaks": len(normalized_peaks) >= 3,
            "collision_energy_verified": True,
            "cross_ce_unit_conversion_applied": False,
            "instrument_domain_verified": True,
            "reported_mass_range_present": reported_mass_range is not None,
            "reported_mass_range_contains_precursor": reported_mass_range is not None,
            "precursor_isolation_window_reported": False,
        }
        quality.update(quality_overrides or {})
        return {
            "schema_version": "2.0.0",
            "spectrum_id": spectrum_id,
            "evidence_stage": "raw_msms_spectrum",
            "experimental_status": "public_repository_experimental_spectrum",
            "reference_identity_status": "repository_reported_reference_identity",
            "query_confirmation_status": "not_applicable",
            "source": {
                "source_id": source_id,
                "laboratory_id": laboratory_id,
                "research_id": research_id,
            },
            "compound": {
                "display_name": f"Compound {spectrum_id}",
                "name_provenance": {"source": "test fixture"},
                "inchikey": "AAAAAAAAAAAAAA-BBBBBBBBBB-C",
                "connectivity": identity,
                "formula": formula,
                "chemical_class": {
                    "ontology": "ChemOnt",
                    "ontology_id": "CHEMONTID:TEST",
                    "leaf_class": chemical_class,
                    "raw": chemical_class,
                },
            },
            "acquisition": {
                "instrument_domain_id": instrument_domain_id,
                "activation": activation,
                "activation_raw": activation,
                "activation_evidence": "synthetic test fixture contributor text",
                "psi_beam_type_cid_confirmed": False,
                "psi_hcd_cv_confirmed": False,
                "activation_architecture_interpretation": (
                    "synthetic test fixture; no controlled-vocabulary activation claim"
                ),
                "collision_cell_region": collision_cell_region,
                "collision_energy": {
                    "raw": f"{ce_value} {ce_unit}",
                    "value": ce_value,
                    "minimum": ce_minimum,
                    "maximum": ce_maximum,
                    "unit": ce_unit,
                    "mode": ce_mode,
                    "conversion_applied": False,
                    "unit_provenance": "synthetic test fixture literal",
                    "cross_vendor_or_cross_model_equivalence_assumed": False,
                },
                "adduct": adduct,
                "polarity": polarity,
                "charge": charge,
                "precursor_mz": precursor_mz,
                "reported_mass_range": deepcopy(reported_mass_range),
            },
            "independence": {
                "exact_peak_hash": peak_hash,
                "mirror_cluster_id": mirror_cluster_id or f"mirror:{spectrum_id}",
                "two_dimensional_identity_id": identity,
                "formula_family_id": formula,
                "source_origin_id": source_id,
                "laboratory_id": laboratory_id,
                "research_id": research_id,
                "independent_fragmentation_group_id": f"experiment:{spectrum_id}",
            },
            "quality": quality,
            "peaks": normalized_peaks,
        }

    return make_record
