import sqlite3

from matching.compound_names import alias_quality, resolve_compound_name
from matching.evidence_network import prepare_feature_nodes
from matching.library_polarity import (
    POLARITY_RULE_VERSION,
    polarity_evidence,
    recompute_polarity_snapshot,
)
from matching.untargeted_workflow import (
    attach_ms2_to_lcms_features,
    build_lcms_consensus_spectra,
    resolve_search_precursor_mz,
)
from scripts.openms_untargeted import (
    MAX_SAFE_ADDUCT_GROUPING_FEATURES,
    _adduct_grouping_decision,
    _effective_adduct_grouping_limit,
)
from utils.spectrum import (
    aggregate_peak_processing_audits,
    process_peaks,
    resolve_peak_processing_config,
)


def test_peak_processing_is_staged_configurable_and_auditable():
    config = resolve_peak_processing_config(
        "sensitive",
        relative_intensity_threshold=0.10,
        precursor_exclusion_da=0.5,
        top_n=3,
    )
    peaks, audit = process_peaks(
        [
            [50.0, 100.0],
            [60.0, 9.0],
            [70.0, 80.0],
            [99.8, 70.0],
            [250.0, 60.0],
            ["bad", 4.0],
        ],
        100.0,
        config=config,
        precursor_charge=1,
    )
    assert peaks == [[50.0, 100.0], [70.0, 80.0]]
    assert [stage["stage"] for stage in audit["stages"]] == [
        "valid_numeric_positive",
        "relative_intensity",
        "precursor_exclusion",
        "upper_mass_bound",
        "top_n",
    ]
    aggregate = aggregate_peak_processing_audits([audit, audit])
    assert aggregate["spectra_audited"] == 2
    assert aggregate["input_peaks"] == 12
    assert aggregate["output_peaks"] == 4


def test_peak_processing_rescues_products_hidden_by_dominant_precursor():
    config = resolve_peak_processing_config(
        "library_match",
        relative_intensity_threshold=0.05,
        precursor_exclusion_da=0.5,
        top_n=15,
    )
    peaks, audit = process_peaks(
        [
            [50.0, 100.0],
            [75.0, 80.0],
            [125.0, 60.0],
            [199.99, 10_000.0],
        ],
        200.0,
        config=config,
    )

    assert peaks == [[50.0, 100.0], [75.0, 80.0], [125.0, 60.0]]
    assert audit["schema_version"] == 2
    assert audit["precursor_dominance_rescue_triggered"] is True
    assert audit["precursor_dominance_rescue_reason"] == (
        "excluded_precursor_dominated_relative_threshold"
    )
    rescue = next(
        stage
        for stage in audit["stages"]
        if stage["stage"] == "precursor_dominance_rescue"
    )
    assert rescue["input_peaks"] == 0
    assert rescue["retained_peaks"] == 3


def test_peak_processing_does_not_invent_minimum_fragment_evidence():
    config = resolve_peak_processing_config(
        "library_match",
        relative_intensity_threshold=0.05,
        precursor_exclusion_da=0.5,
        top_n=15,
    )
    peaks, audit = process_peaks(
        [[50.0, 100.0], [75.0, 80.0], [199.99, 10_000.0]],
        200.0,
        config=config,
    )

    assert peaks == []
    assert audit["precursor_dominance_rescue_triggered"] is False
    assert all(
        stage["stage"] != "precursor_dominance_rescue"
        for stage in audit["stages"]
    )


def test_large_feature_map_skips_only_optional_adduct_ilp():
    assert _adduct_grouping_decision(5_959, True, 10_000) == (
        True,
        "running",
    )
    # Real POS sample regression: 7,656 features completed feature finding in
    # 2.24 s but exceeded the old 300 s total bound inside the optional ILP.
    assert _adduct_grouping_decision(7_656, True, 10_000) == (
        False,
        "skipped_feature_limit_preserved_feature_map",
    )
    assert _adduct_grouping_decision(29_507, True, 10_000) == (
        False,
        "skipped_feature_limit_preserved_feature_map",
    )
    assert _adduct_grouping_decision(29_507, False, 10_000) == (
        False,
        "disabled_preserved_feature_map",
    )
    assert MAX_SAFE_ADDUCT_GROUPING_FEATURES == 6_000
    assert _effective_adduct_grouping_limit(10_000) == 6_000
    assert _effective_adduct_grouping_limit(4_000) == 4_000
    assert _effective_adduct_grouping_limit(0) == 6_000


def test_library_polarity_rule_is_conflict_aware_and_snapshot_is_consistent(
    tmp_path,
):
    db = tmp_path / "library.db"
    connection = sqlite3.connect(db)
    connection.executescript(
        """
        CREATE TABLE Spectra (
            id INTEGER PRIMARY KEY, adduct TEXT, charge TEXT, ion_mode TEXT
        );
        CREATE TABLE Mapping (compound_id INTEGER, spectrum_id INTEGER);
        INSERT INTO Spectra VALUES
            (1, '[M+H]+', '-1', 'NEG'),
            (2, '[M-H]-', '+1', NULL),
            (3, NULL, '+1', NULL),
            (4, NULL, NULL, NULL);
        INSERT INTO Mapping VALUES (1,1), (1,3), (2,2), (3,4);
        """
    )
    connection.commit()
    connection.close()

    explicit = polarity_evidence("[M+H]+", "-1", "NEG")
    assert explicit["polarity"] == "NEG"
    assert explicit["source"] == "explicit_ion_mode"
    assert explicit["conflict"] is True
    assert polarity_evidence(None, 1, None)["polarity"] is None
    assert polarity_evidence(None, "+1", None)["polarity"] == "POS"

    snapshot = recompute_polarity_snapshot(db)
    assert snapshot["polarity_rule_version"] == POLARITY_RULE_VERSION
    assert snapshot["positive_spectra"] == 1
    assert snapshot["negative_spectra"] == 2
    assert snapshot["unknown_polarity_spectra"] == 1
    assert snapshot["polarity_metadata_conflicts"] == 2
    assert snapshot["both_polarities"] == 1
    assert snapshot["negative_only_compounds"] == 1
    assert snapshot["unknown_only_compounds"] == 1


def test_dangerous_aliases_are_rejected_without_rejecting_real_names():
    for unsafe in (
        "methyl",
        "[M-H]-",
        "MEG_1023",
        "unknown?",
        "candidate_major",
        "putative collision 20 eV",
    ):
        accepted, _ = alias_quality(unsafe)
        assert accepted is False
    assert alias_quality("Naringenin") == (True, "safe")
    assert alias_quality("3-O-Caffeoylquinic acid") == (True, "safe")
    assert alias_quality("systematic-" + "hydroxy-" * 40 + "acid") == (
        True,
        "safe",
    )
    assert alias_quality("Dereplicator Identification - candidate compound")[0] is False
    for unsafe in ("methyl", "[M-H]-", "candidate_major"):
        display, status = resolve_compound_name(unsafe)
        assert display.startswith("Unresolved")
        assert status.startswith("unsafe_source_name_")


def test_ms2_links_to_chromatographic_feature_and_builds_consensus():
    spectra = [
        {
            "spectrum_id": "scan=10",
            "precursor_mz": 300.0008,
            "retention_time": 5.02,
            "precursor_charge": 1,
            "peaks": [[50.0, 100.0], [75.0, 50.0], [120.0, 30.0]],
        },
        {
            "spectrum_id": "scan=11",
            "precursor_mz": 300.0010,
            "retention_time": 5.04,
            "precursor_charge": 1,
            "peaks": [[50.01, 80.0], [75.01, 60.0], [140.0, 25.0]],
        },
    ]
    features = [
        {
            "feature_id": "LF0000001",
            "mz": 300.0009,
            "rt_minutes": 5.03,
            "rt_start_minutes": 4.95,
            "rt_end_minutes": 5.10,
            "intensity": 100000.0,
            "charge": 1,
            "isotope_trace_count": 2,
            "adduct": "[M+H]+",
            "neutral_mass": 298.9936,
            "adduct_group_id": "AG000001",
        }
    ]
    attachment = attach_ms2_to_lcms_features(spectra, features)
    assert attachment["linked_ms2_spectra"] == 2
    assert {row["lcms_feature_id"] for row in spectra} == {"LF0000001"}
    mass_audit = resolve_search_precursor_mz(spectra[0])
    assert mass_audit["selected_ion_mz"] == 300.0008
    assert mass_audit["search_precursor_mz"] == 300.0009
    assert mass_audit["search_precursor_mz_source"] == (
        "measured_continuous_ms1_feature_centroid"
    )
    assert mass_audit["lcms_feature_used_for_search"] is True

    rows = [
        {
            **spectrum,
            "raw_peaks": spectrum["peaks"],
            "query_peaks": spectrum["peaks"],
            "ion_mode": "POS",
            "collision_energy": 15.0,
            "collision_energy_unit": "eV",
            "collision_energy_mode": "fixed",
            "activation_method": "beam-type collision-induced dissociation",
            "instrument_configuration_id": "IC1",
            "match": None,
        }
        for spectrum in spectra
    ]
    artifacts, summary = build_lcms_consensus_spectra(rows)
    assert summary["consensus_spectra"] == 0
    assert summary["multi_spectrum_consensus"] == 0
    assert summary["multi_spectrum_representative_groups"] == 1
    assert summary["representative_spectra"] == 1
    assert artifacts[0]["member_count"] == 2
    assert sum(bool(row.get("lcms_consensus_is_representative")) for row in rows) == 1

    nodes, preparation = prepare_feature_nodes(rows, top_n=70)
    assert len(nodes) == 1
    assert nodes[0]["feature_source"] == "openms_lcms_feature"
    assert nodes[0]["peak_source"] == "lcms_feature_representative"
    assert nodes[0]["peak_evidence"]["spectrum_id"] == artifacts[0]["representative_spectrum_id"]
    assert nodes[0]["usable_member_count"] == 1
    assert preparation["lcms_feature_nodes"] == 1
    assert preparation["fallback_precursor_event_nodes"] == 0


def test_search_precursor_mass_fails_closed_without_unique_compatible_feature():
    ambiguous = {
        "precursor_mz": 417.115112,
        "precursor_charge": -1,
        "lcms_feature_link_status": "ambiguous_best_feature",
        "lcms_feature_mz": 417.118881,
        "lcms_feature_charge": -1,
    }
    conflict = {
        **ambiguous,
        "lcms_feature_link_status": "linked_openms_feature",
        "lcms_feature_charge": 2,
    }

    ambiguous_audit = resolve_search_precursor_mz(ambiguous)
    conflict_audit = resolve_search_precursor_mz(conflict)

    assert ambiguous_audit["search_precursor_mz"] == 417.115112
    assert ambiguous_audit["lcms_feature_used_for_search"] is False
    assert ambiguous_audit["fallback_reason"] == "ambiguous_best_feature"
    assert conflict_audit["search_precursor_mz"] == 417.115112
    assert conflict_audit["fallback_reason"] == "linked_feature_charge_conflict"


def test_negative_mode_vendor_charge_magnitude_matches_signed_openms_charge():
    spectrum = {
        "precursor_mz": 417.115112,
        "precursor_charge": 1,
        "polarity": "NEG",
        "retention_time": 6.0413,
        "peaks": [[151.0, 100.0]],
    }
    feature = {
        "feature_id": "LF_NEG_1",
        "mz": 417.118881,
        "rt_minutes": 6.0587,
        "rt_start_minutes": 6.00,
        "rt_end_minutes": 6.10,
        "intensity": 100000.0,
        "charge": -1,
    }

    summary = attach_ms2_to_lcms_features(
        [spectrum], [feature], mz_tolerance_ppm=20.0
    )
    audit = resolve_search_precursor_mz(spectrum)

    assert summary["linked_ms2_spectra"] == 1
    assert summary["charge_conflict_candidates_rejected"] == 0
    assert audit["lcms_feature_used_for_search"] is True
    assert audit["search_precursor_mz"] == 417.118881


def test_network_keeps_low_information_spectra_as_explicit_singletons():
    rows = [
        {
            "spectrum_id": "scan=poor",
            "precursor_mz": 100.0,
            "raw_peaks": [[99.9, 100.0]],
            "query_peaks": [],
            "match": None,
        }
    ]
    nodes, preparation = prepare_feature_nodes(rows)
    assert len(nodes) == 1
    assert nodes[0]["peaks_count"] == 0
    assert (
        nodes[0]["network_eligibility"]
        == "singleton_insufficient_fragments"
    )
    assert preparation["spectra_in_nodes"] == 1
    assert preparation["input_node_coverage_pct"] == 100.0
    assert preparation["insufficient_fragment_nodes"] == 1
