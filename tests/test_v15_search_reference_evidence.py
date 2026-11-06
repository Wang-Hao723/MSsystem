from __future__ import annotations

import importlib

from matching.v15_flagship_ranker import build_v15_flagship_result


search_module = importlib.import_module("matching.search")

PRECURSOR = 493.57067832056
QUERY = [[80.0, 15.0], [100.0, 100.0], [120.0, 50.0], [140.0, 25.0]]
DISTINCT = [[80.0, 10.0], [100.0, 45.0], [125.0, 100.0], [145.0, 40.0]]
INCHIKEY = "AAAAAAAAAAAAAA-BBBBBBBBBB-C"


def _record(index: int, source: str, peaks: list[list[float]], ce: float) -> dict:
    return {
        "precursor_mz": PRECURSOR,
        "neutral_mass": None,
        "neutral_mass_source": "",
        "peaks": peaks,
        "clean_peaks": peaks,
        "spectrum_id": f"ref-{index}",
        "source_accession": f"ACC-{index}",
        "compound_name": "same compound",
        "source_compound_name": "",
        "name_resolution_status": "",
        "smiles": "C" * 35,
        "has_structure": True,
        "name_only_lead": False,
        "inchikey": INCHIKEY,
        "reported_inchikey": INCHIKEY,
        "formula": "C35H72",
        "adduct": "[M+H]+",
        "charge": 1,
        "collision_energy": ce,
        "collision_energy_raw": f"{ce:g} eV",
        "collision_energy_unit": "eV",
        "collision_energy_mode": "single",
        "activation_method": "CID",
        "instrument_domain_id": "synthetic_qtof_fixed_ev",
        "collision_cell_region": "collision_cell",
        "instrument_type": "qtof",
        "instrument": "synthetic QTOF",
        "retention_time": None,
        "classification": "other",
        "classification_source": "test",
        "classification_confidence": None,
        "dataset_source": source,
        "source_tier": "experimental",
        "intensity_origin": "measured",
        "quality_flags": [],
        "pathway_json": "",
        "polarity": "POS",
    }


def _disable_optional_layers(monkeypatch) -> None:
    monkeypatch.setattr(search_module, "EXPERT_SYSTEM_ENABLED", False)
    monkeypatch.setattr(search_module, "USE_EXPERT_V3", False)
    monkeypatch.setattr(search_module, "USE_PATHWAY_ALIGN", False)
    monkeypatch.setattr(search_module, "USE_CONSISTENCY_FILTER", False)
    monkeypatch.setattr(search_module, "USE_ADAPTIVE_PEAK_SEARCH", False)
    monkeypatch.setattr(search_module, "USE_ISOMER_FRAGMENT_RERANK", False)
    monkeypatch.setattr(search_module, "USE_ISOPAIR_MS", False)
    monkeypatch.setattr(search_module, "MODEL_PREDICTOR_AVAILABLE", False)
    monkeypatch.setattr(search_module, "_ensure_clf", lambda: False)


def test_search_retains_all_reference_records_and_v15_collapses_mirrors(monkeypatch):
    _disable_optional_layers(monkeypatch)
    records = [
        _record(0, "GNPS", QUERY, 20.0),
        _record(1, "MassBank", [[mz, intensity * 10] for mz, intensity in QUERY], 20.0),
        _record(2, "MoNA", DISTINCT, 22.0),
    ]
    bucket = int(PRECURSOR / search_module.MZ_BUCKET_DALTON)
    result = search_module.search(
        QUERY,
        PRECURSOR,
        records=records,
        mz_index={bucket: [0, 1, 2]},
        query_polarity="POS",
        query_charge=1,
        query_collision_energy={"value": 20.0, "unit": "eV"},
        query_collision_energy_mode="fixed",
        query_instrument_type="qtof",
        query_instrument="synthetic QTOF",
        query_activation="CID",
        query_instrument_domain_id="synthetic_qtof_fixed_ev",
        query_collision_cell_region="collision_cell",
        cos_thresh=0.0,
        final_thresh=0.0,
        strict_evidence_gate=False,
        enable_chimeric_deconvolution=False,
    )

    assert result is not None
    assert len(result["top_candidates"]) == 1
    references = result["top_candidates"][0]["reference_evidence_records"]
    assert len(references) == 3
    assert {record["dataset_source"] for record in references} == {
        "GNPS",
        "MassBank",
        "MoNA",
    }
    assert len({record["library_peak_fingerprint"] for record in references}) == 2

    row = {
        "spectrum_id": "query-1",
        "lcms_feature_id": "feature-1",
        "precursor_mz": PRECURSOR,
        "search_precursor_mz": PRECURSOR,
        "ion_mode": "POS",
        "collision_energy": 20.0,
        "collision_energy_unit": "eV",
        "collision_energy_mode": "fixed",
        "instrument_type": "qtof",
        "instrument": "synthetic QTOF",
        "activation_method": "CID",
        "instrument_domain_id": "synthetic_qtof_fixed_ev",
        "collision_cell_region": "collision_cell",
        "raw_peaks": QUERY,
        "query_peaks": QUERY,
        "lib_peaks": result["lib_peaks"],
        "match": result,
        "tic": 1000.0,
    }
    payload = build_v15_flagship_result("search-reference-test", [row])
    candidate = payload["entities"][0]["candidates"][0]
    evidence = candidate["collision_energy_evidence"][0][
        "conditioned_multilibrary_evidence"
    ]

    assert evidence["reference_record_count"] == 3
    assert evidence["unique_content_profile_count"] == 2
    assert evidence["mirror_duplicate_record_count"] == 1
    assert evidence["conditioned_provenance_stratum_count"] == 2
    assert evidence["selected_consensus_basis"] == (
        "compatible_instrument_domain_near_ce_stratum_median"
    )
    assert evidence["runtime_influence"] is True
    assert evidence["accessions_count_as_independent_evidence"] is False

    from web.server import _serialize_results

    serialized = _serialize_results([{
        **row,
        "query_peaks": QUERY,
        "raw_peaks": QUERY,
        "peaks_count": len(QUERY),
        "spectrum_representation": "centroid",
        "raw_peaks_complete": True,
    }])[0]
    assert serialized["collision_energy_mode"] == "fixed"
    assert serialized["activation_method"] == "CID"
    assert serialized["instrument_domain_id"] == "synthetic_qtof_fixed_ev"
    assert serialized["collision_cell_region"] == "collision_cell"
    assert serialized["match"]["lib_activation_method"] == "CID"
    assert serialized["match"]["lib_instrument_domain_id"] == (
        "synthetic_qtof_fixed_ev"
    )
    assert serialized["match"]["lib_collision_cell_region"] == "collision_cell"
