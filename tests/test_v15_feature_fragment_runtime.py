from __future__ import annotations

import importlib

import pytest

from matching.evidence_governance import build_evidence_passport
from matching.feature_fragments import load_catalog, resolve_instrument_domain
from utils.mgf_io import parse_mgf_file
from utils.mzml_io import _activation_method, _collision_energy_mode


search_module = importlib.import_module("matching.search")
PRECURSOR = 493.57067832056
PEAKS = [[100.0, 20.0], [184.0733, 100.0], [250.0, 15.0]]


def _domain(**kwargs):
    return resolve_instrument_domain(catalog=load_catalog(), **kwargs)


def test_qe_nce_and_ev_are_distinct_domains_without_conversion():
    nce = _domain(
        instrument="Thermo Q Exactive HF-X",
        instrument_type="Orbitrap",
        activation="HCD",
        collision_energy_unit="NCE",
        collision_energy_mode="fixed",
    )
    ev = _domain(
        instrument="Thermo Q Exactive HF-X",
        instrument_type="Orbitrap",
        activation="HCD",
        collision_energy_unit="eV",
        collision_energy_mode="fixed",
    )

    assert nce["status"] == "resolved"
    assert nce["instrument_domain_id"] == "thermo_qe_orbitrap_hcd_fixed_nce"
    assert ev["status"] == "resolved"
    assert ev["instrument_domain_id"] == "thermo_qe_orbitrap_hcd_fixed_ev"
    assert nce["cross_unit_conversion_applied"] is False
    assert ev["cross_unit_conversion_applied"] is False


def test_qtof_voltage_and_energy_domains_are_not_interchanged():
    ev = _domain(
        instrument="Agilent 6560 QTOF",
        instrument_type="QTOF",
        activation="CID",
        collision_energy_unit="eV",
        collision_energy_mode="fixed",
    )
    voltage = _domain(
        instrument="Agilent 6560 QTOF",
        instrument_type="QTOF",
        activation="CID",
        collision_energy_unit="V",
        collision_energy_mode="fixed",
    )

    assert ev["instrument_domain_id"] == (
        "agilent_6560_qtof_beam_cid_fixed_ev_reported"
    )
    assert voltage["instrument_domain_id"] == (
        "agilent_6500_qtof_beam_cid_fixed_v"
    )
    assert ev["cross_unit_conversion_applied"] is False
    assert voltage["cross_unit_conversion_applied"] is False


def test_waters_trap_and_transfer_require_collision_region():
    ambiguous = _domain(
        instrument="Waters SYNAPT G2-Si",
        instrument_type="QTOF",
        activation="beam-type collision-induced dissociation",
        collision_energy_unit="eV",
        collision_energy_mode="fixed",
    )
    trap = _domain(
        instrument="Waters SYNAPT G2-Si",
        instrument_type="QTOF",
        activation="beam-type collision-induced dissociation",
        collision_energy_unit="eV",
        collision_energy_mode="fixed",
        collision_cell_region="trap",
    )

    assert ambiguous["status"] == "unavailable"
    assert ambiguous["reason"] == "instrument_domain_ambiguous"
    assert trap["instrument_domain_id"] == "waters_synapt_trap_cid_fixed_ev"


def test_mgf_preserves_activation_ce_mode_and_peak_list_quality(tmp_path):
    path = tmp_path / "agilent_pc.mgf"
    path.write_text(
        "BEGIN IONS\nTITLE=pc\nPEPMASS=493.570678\nCHARGE=1+\n"
        "ADDUCT=[M+H]+\nIONMODE=Positive\nINSTRUMENT=Agilent 6560 QTOF\n"
        "ACTIVATION=CID\nCOLLISIONENERGY=20 eV\n"
        "100 20\n184.0733 100\n250 15\nEND IONS\n",
        encoding="utf-8",
    )

    record = parse_mgf_file(path)[0]

    assert record["activation_method"] == "CID"
    assert record["collision_energy_mode"] == "fixed"
    assert record["collision_energy_unit"] == "eV"
    assert record["spectrum_representation"] == "centroid"
    assert record["raw_peaks_complete"] is True


def test_mzml_empty_cv_flag_still_preserves_activation_method():
    activation = {
        "collision-induced dissociation": "",
        "collision energy": 20.0,
    }

    assert _activation_method(activation) == "collision-induced dissociation"
    assert _collision_energy_mode(activation, 20.0) == "fixed"


def _record():
    return {
        "precursor_mz": PRECURSOR,
        "neutral_mass": None,
        "neutral_mass_source": "",
        "peaks": PEAKS,
        "clean_peaks": PEAKS,
        "spectrum_id": "pc-reference",
        "source_accession": "TEST-PC-1",
        "compound_name": "test phosphocholine candidate",
        "source_compound_name": "",
        "name_resolution_status": "",
        "smiles": "C" * 35,
        "has_structure": True,
        "name_only_lead": False,
        "inchikey": "AAAAAAAAAAAAAA-BBBBBBBBBB-C",
        "reported_inchikey": "",
        "formula": "C35H72",
        "adduct": "[M+H]+",
        "charge": 1,
        "collision_energy": 20.0,
        "collision_energy_raw": "20 eV",
        "collision_energy_unit": "eV",
        "collision_energy_mode": "fixed",
        "instrument_type": "qtof",
        "instrument": "Agilent 6560 QTOF",
        "retention_time": None,
        "classification": "Glycerophosphocholines",
        "classification_source": "curated_massbank",
        "classification_confidence": 1.0,
        "dataset_source": "unit-test",
        "source_tier": "experimental",
        "intensity_origin": "measured",
        "quality_flags": [],
        "pathway_json": "",
        "polarity": "POS",
    }


def _disable_optional_layers(monkeypatch):
    monkeypatch.setattr(search_module, "EXPERT_SYSTEM_ENABLED", False)
    monkeypatch.setattr(search_module, "USE_EXPERT_V3", False)
    monkeypatch.setattr(search_module, "USE_PATHWAY_ALIGN", False)
    monkeypatch.setattr(search_module, "USE_CONSISTENCY_FILTER", False)
    monkeypatch.setattr(search_module, "USE_ADAPTIVE_PEAK_SEARCH", False)
    monkeypatch.setattr(search_module, "USE_ISOMER_FRAGMENT_RERANK", False)
    monkeypatch.setattr(search_module, "USE_ISOPAIR_MS", False)
    monkeypatch.setattr(search_module, "MODEL_PREDICTOR_AVAILABLE", False)
    monkeypatch.setattr(search_module, "_ensure_clf", lambda: False)


def _run_search(*, activation):
    bucket = int(PRECURSOR / search_module.MZ_BUCKET_DALTON)
    return search_module.search(
        PEAKS,
        PRECURSOR,
        records=[_record()],
        mz_index={bucket: [0]},
        query_adduct="[M+H]+",
        query_adduct_reliable=True,
        query_polarity="POS",
        query_charge=1,
        query_collision_energy={"value": 20.0, "unit": "eV"},
        query_collision_energy_mode="fixed",
        query_instrument_type="qtof",
        query_instrument="Agilent 6560 QTOF",
        query_activation=activation,
        query_spectrum_representation="centroid",
        query_raw_peaks_complete=True,
        query_isolation_evidence={"status": "clean", "purity": 0.95},
        cos_thresh=0.0,
        final_thresh=0.0,
        strict_evidence_gate=False,
        enable_chimeric_deconvolution=False,
    )


def test_feature_fragment_runtime_gate_is_auditable_and_rank_invariant(monkeypatch):
    _disable_optional_layers(monkeypatch)
    applicable = _run_search(activation="CID")
    unavailable = _run_search(activation="")

    assert applicable is not None and unavailable is not None
    evidence = applicable["feature_fragment_evidence"]
    assert evidence["status"] == "unavailable"
    assert evidence["reason"] == (
        "no_feature_fragment_rule_passed_sealed_runtime_gate"
    )
    assert evidence["rule_evidence"] == []
    assert evidence["combined_log_likelihood_ratio"] == 0.0
    assert evidence["runtime_influence"] is False
    assert evidence["identity_rank_influence"] is False
    assert evidence["exact_identity_confirmed"] is False
    assert evidence["standard_confirmed"] is False
    gate = evidence["runtime_export_gate"]
    assert gate["deployment_state"] == "fail_closed"
    assert gate["exported_rule_count"] == 0
    assert gate["pilot_v1_rule_count"] == 6
    assert gate["pilot_v1_runtime_allowed"] is False
    assert gate["research_catalog_evaluated"] is False

    assert unavailable["feature_fragment_evidence"]["status"] == "unavailable"
    assert unavailable["feature_fragment_evidence"]["reason"] == evidence["reason"]
    assert applicable["inchikey"] == unavailable["inchikey"]
    assert applicable["final_score"] == pytest.approx(unavailable["final_score"])
    assert applicable["rank_score"] == pytest.approx(unavailable["rank_score"])


def test_search_passes_complete_centroid_list_to_feature_rules(monkeypatch):
    _disable_optional_layers(monkeypatch)
    captured = {}

    def capture(**kwargs):
        captured["peaks"] = kwargs["query_peaks"]
        captured["complete"] = kwargs["query_raw_peaks_complete"]
        return search_module._feature_fragment_unavailable("captured")

    monkeypatch.setattr(search_module, "_evaluate_feature_fragment_shadow", capture)
    peaks = [
        [100.0, 100.0],
        [184.0733, 80.0],
        [250.0, 60.0],
        [333.3333, 0.1],
    ]
    bucket = int(PRECURSOR / search_module.MZ_BUCKET_DALTON)
    result = search_module.search(
        peaks,
        PRECURSOR,
        records=[_record()],
        mz_index={bucket: [0]},
        query_adduct="[M+H]+",
        query_adduct_reliable=True,
        query_polarity="POS",
        query_charge=1,
        query_collision_energy={"value": 20.0, "unit": "eV"},
        query_collision_energy_mode="fixed",
        query_instrument_type="qtof",
        query_instrument="Agilent 6560 QTOF",
        query_activation="CID",
        query_spectrum_representation="centroid",
        query_raw_peaks_complete=True,
        query_isolation_evidence={"status": "clean"},
        top_n=3,
        rel_thresh=0.05,
        cos_thresh=0.0,
        final_thresh=0.0,
        strict_evidence_gate=False,
        enable_chimeric_deconvolution=False,
    )

    assert result is not None
    assert captured["complete"] is True
    assert captured["peaks"] == peaks
    assert [333.3333, 0.1] in captured["peaks"]


def test_evidence_governance_records_feature_fragment_gate_as_unavailable():
    evidence = search_module._evaluate_feature_fragment_shadow(
        query_peaks=PEAKS,
        query_pmz=PRECURSOR,
        candidate_class="Glycerophosphocholines",
        classification_reliable=True,
        classification_source="curated_massbank",
        query_polarity="POS",
        query_adduct="[M+H]+",
        query_adduct_reliable=True,
        query_charge=1,
        query_collision_energy={"value": 20.0, "unit": "eV"},
        query_collision_energy_mode="fixed",
        query_instrument_type="qtof",
        query_instrument="Agilent 6560 QTOF",
        query_activation="CID",
        query_instrument_domain_id=None,
        query_collision_cell_region=None,
        query_spectrum_representation="centroid",
        query_raw_peaks_complete=True,
        query_isolation_evidence={"status": "clean"},
    )
    candidate = {
        **_record(),
        "cosine": 0.95,
        "entropy": 0.90,
        "reverse_match": 0.95,
        "final_score": 0.80,
        "rank_score": 0.80,
        "expert_score": 0.0,
        "expert_details": {},
        "mass_ppm": 1.0,
        "mass_metadata_consistent": True,
        "confidence": "Medium",
        "isomer_ambiguous": False,
        "top_candidates": [],
        "feature_fragment_evidence": evidence,
    }
    passport = build_evidence_passport({
        "spectrum_id": "query-1",
        "collision_energy": 20.0,
        "instrument": "Agilent 6560 QTOF",
        "blank_status": "not_requested",
        "match": candidate,
    })
    node = next(
        item for item in passport["nodes"]
        if item["key"] == "feature_fragment_rule_support"
    )

    assert node["status"] == "unavailable"
    assert node["reason"] == "no_feature_fragment_rule_passed_sealed_runtime_gate"
    assert node["used_for"] == ["explanation"]
    assert node["independent_measurement"] is False
    assert "release_gate" not in node["used_for"]
    assert "ranking" not in node["used_for"]
    assert node["value"]["informative_rule_count"] == 0
    assert node["value"]["runtime_export_gate"]["exported_rule_count"] == 0
    assert passport["ranking_changed"] is False
