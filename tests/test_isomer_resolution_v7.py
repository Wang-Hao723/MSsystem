"""Regression guards for calibrated isomer-resolution evidence."""

from __future__ import annotations

import copy
import csv
import json

import pytest

from matching.honest_calibration import (
    annotate,
    normalize_isomer_fragment_resolution_status,
)
from matching.isomer_literature import (
    build_isomer_literature_context,
    search_pubmed_isomer_evidence,
)
from matching import isomer_literature
from matching.isomer_resolution import (
    same_connectivity_variant_subtype,
)
from matching.rt_reference import load_rt_standards, match_rt_evidence
from matching.search import _apply_isomer_fragment_rerank
from utils.instrument_metadata import (
    extract_mzml_instrument_metadata,
    normalize_instrument_family,
)


IK_A = "AAAAAAAAAAAAAA-BBBBBBBBBB-C"
IK_A_VARIANT = "AAAAAAAAAAAAAA-CCCCCCCCCC-D"
IK_B = "DDDDDDDDDDDDDD-EEEEEEEEEE-F"


def _candidate(
    name: str,
    inchikey: str,
    peaks: list[list[float]],
    *,
    spectrum_id: str,
) -> dict:
    return {
        "compound_name": name,
        "inchikey": inchikey,
        "formula": "C10H12O2",
        "smiles": "CCOc1ccccc1",
        "has_structure": True,
        "classification": "flavonoid",
        "source_tier": "experimental",
        "cosine": 0.95,
        "final_score": 0.70,
        "rank_score": 0.70,
        "reverse_match": 0.90,
        "mass_ppm": 2.0,
        "mass_metadata_consistent": True,
        "adduct": "[M+H]+",
        "query_adduct": "[M+H]+",
        "query_adduct_reliable": True,
        "lib_precursor_mz": 200.0,
        "lib_collision_energy": 20.0,
        "lib_collision_energy_unit": "eV",
        "lib_instrument_type": "QTOF",
        "lib_instrument": "Agilent 6545 Q-TOF",
        "lib_spectrum_id": spectrum_id,
        "lib_source_ms_level": "MS2",
        "lib_peaks": peaks,
    }


def _replicated_constitutional_candidates() -> list[dict]:
    a_peaks = [[100.0, 1.0], [120.0, 0.8], [150.0, 0.5]]
    b_peaks = [[110.0, 1.0], [130.0, 0.8], [150.0, 0.5]]
    return [
        _candidate("isomer A", IK_A, a_peaks, spectrum_id="a-1"),
        _candidate("isomer A", IK_A, a_peaks, spectrum_id="a-2"),
        _candidate("isomer B", IK_B, b_peaks, spectrum_id="b-1"),
        _candidate("isomer B", IK_B, b_peaks, spectrum_id="b-2"),
    ]


def _run_isomer_evidence(
    candidates: list[dict],
    query_peaks: list[list[float]],
) -> dict:
    original_ranks = [candidate["rank_score"] for candidate in candidates]
    _apply_isomer_fragment_rerank(
        candidates,
        query_peaks,
        query_pmz=200.0,
        query_collision_energy={"value": 20.0, "unit": "eV"},
        query_instrument_type="Q-TOF",
    )
    assert [candidate["rank_score"] for candidate in candidates] == original_ranks
    assert all(
        candidate["isomer_fragment_rank_applied"] is False
        for candidate in candidates
    )
    return candidates[0]["isomer_fragment_detail"]


def test_replicated_bidirectional_evidence_supports_selected_without_promoting():
    detail = _run_isomer_evidence(
        _replicated_constitutional_candidates(),
        [[100.0, 1.0], [120.0, 0.8], [150.0, 0.5]],
    )

    assert detail["schema_version"] == 3
    assert detail["ambiguity_type"] == (
        "constitutional_or_positional_isomer_group"
    )
    assert detail["reference_replication"]["status"] == "complete"
    assert detail["metadata_comparability"] == "complete"
    assert detail["independent_supported_features"] == 2
    assert detail["opposing_supported_features"] == 0
    assert detail["resolution_status"] == (
        "selected_candidate_supported_pending_validation"
    )
    assert detail["eligible_for_exact_identity_promotion"] is False


def test_bidirectional_diagnostic_evidence_detects_conflict():
    detail = _run_isomer_evidence(
        _replicated_constitutional_candidates(),
        [
            [100.0, 1.0],
            [110.0, 0.9],
            [120.0, 0.8],
            [130.0, 0.7],
            [150.0, 0.5],
        ],
    )

    assert detail["independent_supported_features"] == 2
    assert detail["opposing_supported_features"] == 2
    assert detail["evidence_balance"] == 0
    assert detail["resolution_status"] == "conflicting_diagnostic_evidence"


def test_competitor_direction_is_not_claimed_without_replicates():
    candidates = _replicated_constitutional_candidates()[::2]
    detail = _run_isomer_evidence(
        candidates,
        [[110.0, 1.0], [130.0, 0.8], [150.0, 0.5]],
    )

    assert detail["independent_supported_features"] == 0
    assert detail["opposing_supported_features"] == 2
    assert detail["reference_replication"]["status"] == "insufficient"
    assert detail["resolution_status"] == (
        "competitor_features_observed_but_unreplicated"
    )


def test_competitor_direction_is_not_claimed_with_incomplete_metadata():
    candidates = _replicated_constitutional_candidates()
    for candidate in candidates:
        candidate["lib_instrument_type"] = ""
        candidate["lib_instrument"] = ""
    detail = _run_isomer_evidence(
        candidates,
        [[110.0, 1.0], [130.0, 0.8], [150.0, 0.5]],
    )

    assert detail["reference_replication"]["status"] == "complete"
    assert detail["metadata_comparability"] == "incomplete_or_mismatched"
    assert detail["resolution_status"] == (
        "competitor_features_observed_but_metadata_incomplete"
    )


def test_cached_directional_status_is_recalibrated_symmetrically():
    detail = {
        "resolution_status": "diagnostic_evidence_favors_competitor",
        "reference_replication": {"status": "insufficient"},
        "metadata_comparability": "complete",
    }
    assert normalize_isomer_fragment_resolution_status(detail) == (
        "competitor_features_observed_but_unreplicated"
    )
    assert detail["resolution_status"] == (
        "competitor_features_observed_but_unreplicated"
    )

    detail = {
        "resolution_status": "conflicting_diagnostic_evidence",
        "reference_replication": {"status": "complete"},
        "metadata_comparability": "incomplete_or_mismatched",
    }
    assert normalize_isomer_fragment_resolution_status(detail) == (
        "conflicting_features_observed_but_metadata_incomplete"
    )


def test_weak_same_formula_neighbor_cannot_create_competitor_warning():
    candidates = _replicated_constitutional_candidates()
    for candidate in candidates[2:]:
        candidate.update({
            "cosine": 0.50,
            "final_score": 0.40,
            "rank_score": 0.40,
        })

    detail = _run_isomer_evidence(
        candidates,
        [[110.0, 1.0], [130.0, 0.8], [150.0, 0.5]],
    )

    assert detail["spectrally_competitive_candidate_count"] == 0
    assert detail["ambiguity_type"] == "none"
    assert detail["opposing_supported_features"] == 0
    assert detail["all_same_formula_opposing_supported_features"] == 2
    assert detail["diagnostic_status"] == "no_spectrally_competitive_isomer"
    assert detail["resolution_status"] == "no_spectrally_competitive_isomer"


def test_diagnostic_mass_tolerance_does_not_accept_old_point_one_da_window():
    detail = _run_isomer_evidence(
        _replicated_constitutional_candidates(),
        [[100.03, 1.0], [150.0, 0.5]],
    )

    product = next(
        item
        for item in detail["candidate_specific_product_ions"]
        if item["library_mz"] == pytest.approx(100.0)
    )
    assert detail["diagnostic_mz_tolerance_da"] == pytest.approx(0.02)
    assert product["observed"] is False


def test_same_connectivity_variant_cannot_be_resolved_by_ordinary_msms():
    candidates = _replicated_constitutional_candidates()
    for candidate in candidates[2:]:
        candidate["inchikey"] = IK_A_VARIANT

    detail = _run_isomer_evidence(
        candidates,
        [[100.0, 1.0], [120.0, 0.8], [150.0, 0.5]],
    )

    assert detail["ambiguity_type"] == "same_connectivity_variant_group"
    assert detail["msms_exact_resolution_capable"] is False
    assert detail["resolution_status"] == (
        "same_connectivity_variant_not_resolved_by_msms"
    )
    assert detail["eligible_for_exact_identity_promotion"] is False


def test_structure_definition_specificity_is_not_mislabeled_as_chiral_pair():
    unspecified_key = "RBBVPNQTBKHOEQ-UHFFFAOYSA-O"
    specific_key = "RBBVPNQTBKHOEQ-BJQOMGFOSA-O"
    selected = _candidate(
        "Phellodendrine",
        unspecified_key,
        [[100.0, 1.0]],
        spectrum_id="phellodendrine",
    )
    competitor = _candidate(
        "Phellodendrine chloride",
        specific_key,
        [[100.0, 1.0]],
        spectrum_id="phellodendrine-chloride",
    )
    competitor.update({"cosine": 0.94, "final_score": 0.69, "rank_score": 0.69})
    selected["top_candidates"] = [copy.deepcopy(selected), competitor]

    result = annotate(selected)

    assert same_connectivity_variant_subtype(
        selected, competitor
    ) == "structure_definition_specificity_variant"
    assert result["isomer_variant_subtype_counts"] == {
        "structure_definition_specificity_variant": 1
    }
    assert result["isomer_resolution_status"] == (
        "same_connectivity_structure_definition_unresolved"
    )
    assert "chiral" not in (
        result["isomer_group_evidence"]["resolution_requirement"].lower()
    )


def test_complete_candidate_vs_underspecified_entry_does_not_imply_chiral_lc():
    specific = _candidate(
        "Phellodendrine defined stereochemistry",
        "RBBVPNQTBKHOEQ-BJQOMGFOSA-O",
        [[100.0, 1.0]],
        spectrum_id="specific",
    )
    unspecified = _candidate(
        "Phellodendrine unspecified",
        "RBBVPNQTBKHOEQ-UHFFFAOYSA-O",
        [[100.0, 1.0]],
        spectrum_id="unspecified",
    )
    unspecified.update({
        "cosine": 0.94,
        "final_score": 0.69,
        "rank_score": 0.69,
    })
    specific["top_candidates"] = [copy.deepcopy(specific), unspecified]

    result = annotate(specific)

    assert result["isomer_resolution_status"] == (
        "same_connectivity_structure_definition_specificity_ambiguity"
    )
    requirement = result["isomer_group_evidence"]["resolution_requirement"]
    assert "chiral method is required only if" in requirement


def test_mixed_connectivity_and_definition_group_reports_both_axes():
    selected = _candidate(
        "Albiflorin",
        "QQUHMASGPODSIW-ICECTASOSA-N",
        [[100.0, 1.0]],
        spectrum_id="albiflorin-specific",
    )
    underspecified = _candidate(
        "Albiflorin unspecified",
        "QQUHMASGPODSIW-UHFFFAOYSA-N",
        [[100.0, 1.0]],
        spectrum_id="albiflorin-unspecified",
    )
    constitutional = _candidate(
        "Paeoniflorin",
        "YKRGDOXKVOZESV-RMNZJSDUSA-N",
        [[100.0, 1.0]],
        spectrum_id="paeoniflorin",
    )
    for candidate, cosine in [
        (underspecified, 0.94),
        (constitutional, 0.93),
    ]:
        candidate.update({
            "cosine": cosine,
            "final_score": 0.69,
            "rank_score": 0.69,
        })
    selected["top_candidates"] = [
        copy.deepcopy(selected),
        underspecified,
        constitutional,
    ]

    result = annotate(selected)

    assert result["isomer_ambiguity_type"] == "mixed_isomer_group"
    assert result["isomer_resolution_status"] == (
        "mixed_isomer_group_requires_multi_axis_resolution"
    )
    requirement = result["isomer_group_evidence"]["resolution_requirement"]
    assert "different-connectivity candidates" in requirement
    assert "complete underspecified structures" in requirement


def test_protonation_layer_variant_requires_standardization_not_chiral_lc():
    left = _candidate(
        "parent ion A",
        "AAAAAAAAAAAAAA-BBBBBBBBBB-C",
        [[100.0, 1.0]],
        spectrum_id="charge-c",
    )
    right = _candidate(
        "parent ion A salt",
        "AAAAAAAAAAAAAA-BBBBBBBBBB-D",
        [[100.0, 1.0]],
        spectrum_id="charge-d",
    )
    right.update({"cosine": 0.94, "final_score": 0.69, "rank_score": 0.69})
    left["top_candidates"] = [copy.deepcopy(left), right]

    result = annotate(left)

    assert same_connectivity_variant_subtype(
        left, right
    ) == "protonation_or_charge_layer_variant"
    assert result["isomer_resolution_status"] == (
        "same_connectivity_protonation_variant_requires_standardization"
    )
    requirement = result["isomer_group_evidence"]["resolution_requirement"]
    assert "chiral method is not implied" in requirement


def test_curated_pair_specific_literature_keeps_zero_scoring_influence():
    selected = _candidate(
        "Ephedrine",
        IK_A,
        [[100.0, 1.0]],
        spectrum_id="ephedrine",
    )
    competitor = _candidate(
        "Pseudoephedrine",
        IK_A_VARIANT,
        [[110.0, 1.0]],
        spectrum_id="pseudoephedrine",
    )
    context = build_isomer_literature_context(
        selected,
        [competitor],
        ambiguity_type="same_connectivity_variant_group",
        variant_subtype_counts={
            "stereochemical_isotopic_or_tautomeric_variant": 1
        },
        isomer_detail={},
    )

    assert context["status"] == (
        "pair_specific_method_evidence_available"
    )
    assert context["pubmed_lookup_required"] is False
    assert context["used_for_scoring"] is False
    assert context["eligible_for_exact_identity_promotion"] is False
    assert "31250489" in {
        source.get("pmid") for source in context["sources"]
    }


def test_pubmed_lookup_returns_bibliographic_leads_only_and_caches(
    tmp_path,
    monkeypatch,
):
    search_payload = json.dumps({
        "esearchresult": {"idlist": ["31250489"]},
    }).encode()
    article_payload = b"""<?xml version="1.0"?>
    <PubmedArticleSet><PubmedArticle><MedlineCitation>
      <PMID>31250489</PMID><Article>
        <Journal><JournalIssue><PubDate><Year>2019</Year></PubDate></JournalIssue>
          <Title>Chirality</Title></Journal>
        <ArticleTitle>Ephedrine and pseudoephedrine chiral LC-MS separation</ArticleTitle>
        <Abstract><AbstractText>Reference standards were separated by a chiral column.</AbstractText></Abstract>
      </Article>
    </MedlineCitation><PubmedData><ArticleIdList>
      <ArticleId IdType="pubmed">31250489</ArticleId>
      <ArticleId IdType="doi">10.1002/chir.23104</ArticleId>
    </ArticleIdList></PubmedData></PubmedArticle></PubmedArticleSet>"""
    calls = []

    def fake_request(url, timeout):
        calls.append(url)
        return article_payload if "efetch.fcgi" in url else search_payload

    monkeypatch.setattr(isomer_literature, "_request_bytes", fake_request)
    first = search_pubmed_isomer_evidence(
        ["Ephedrine", "Pseudoephedrine"],
        formula="C10H15NO",
        cache_dir=tmp_path,
    )
    second = search_pubmed_isomer_evidence(
        ["Ephedrine", "Pseudoephedrine"],
        formula="C10H15NO",
        cache_dir=tmp_path,
    )

    assert first["status"] == "bibliographic_leads_found"
    assert first["used_for_scoring"] is False
    assert first["records"][0]["diagnostic_fragment_data_extracted"] is False
    assert first["records"][0]["evidence_scope"] == (
        "pair_specific_bibliographic_method_lead"
    )
    assert first["records"][0]["used_for_ranking"] is False
    assert "_abstract_for_scope_only" not in first["records"][0]
    assert second["cache_hit"] is True
    assert len(calls) == 2


def test_pubmed_single_readable_name_is_never_pair_specific(
    tmp_path,
    monkeypatch,
):
    search_payload = json.dumps({
        "esearchresult": {"idlist": ["38776711"]},
    }).encode()
    article_payload = b"""<?xml version="1.0"?>
    <PubmedArticleSet><PubmedArticle><MedlineCitation>
      <PMID>38776711</PMID><Article>
        <Journal><JournalIssue><PubDate><Year>2024</Year></PubDate></JournalIssue>
          <Title>Journal of Chromatography</Title></Journal>
        <ArticleTitle>Pseudoephedrine chiral LC-MS analysis</ArticleTitle>
        <Abstract><AbstractText>Pseudoephedrine was separated by a chiral column.</AbstractText></Abstract>
      </Article>
    </MedlineCitation><PubmedData><ArticleIdList>
      <ArticleId IdType="pubmed">38776711</ArticleId>
    </ArticleIdList></PubmedData></PubmedArticle></PubmedArticleSet>"""

    def fake_request(url, timeout):
        return article_payload if "efetch.fcgi" in url else search_payload

    monkeypatch.setattr(isomer_literature, "_request_bytes", fake_request)
    result = search_pubmed_isomer_evidence(
        ["Pseudoephedrine"],
        formula="C10H15NO",
        cache_dir=tmp_path,
    )

    assert result["status"] == "bibliographic_leads_found"
    assert result["records"][0]["candidate_name_hits"] == [
        "Pseudoephedrine"
    ]
    assert result["records"][0]["evidence_scope"] == (
        "single_candidate_or_general_bibliographic_lead"
    )


def _write_rt_standards(path, *, stereochemistry_capable: bool) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "compound_name",
                "inchikey",
                "formula",
                "rt_sec",
                "tolerance_sec",
                "method_id",
                "separation_mode",
                "stereochemistry_capable",
            ],
        )
        writer.writeheader()
        writer.writerow({
            "compound_name": "isomer A",
            "inchikey": IK_A,
            "formula": "C10H12O2",
            "rt_sec": "120",
            "tolerance_sec": "15",
            "method_id": "LC_METHOD_A",
            "separation_mode": (
                "chiral LC" if stereochemistry_capable else "reversed phase"
            ),
            "stereochemistry_capable": (
                "true" if stereochemistry_capable else "false"
            ),
        })


def test_rt_requires_exact_full_identity_and_same_explicit_method(tmp_path):
    path = tmp_path / "rt.csv"
    _write_rt_standards(path, stereochemistry_capable=False)
    load_rt_standards.cache_clear()
    selected = _candidate("isomer A", IK_A, [[100.0, 1.0]], spectrum_id="a")

    exact = match_rt_evidence(
        selected,
        2.0,
        path,
        query_method_id="LC METHOD A",
    )
    assert exact["matched"] is True
    assert exact["matched_by"] == "full_inchikey"
    assert exact["method_status"] == "same_method"
    assert exact["eligible_for_exact_identity"] is True

    connectivity_only = match_rt_evidence(
        {**selected, "inchikey": IK_A_VARIANT},
        2.0,
        path,
        query_method_id="LC METHOD A",
    )
    assert connectivity_only["matched"] is True
    assert connectivity_only["identity_scope"] == "connectivity"
    assert connectivity_only["eligible_for_exact_identity"] is False

    missing_method = match_rt_evidence(selected, 2.0, path)
    assert missing_method["matched"] is False
    assert missing_method["method_status"] == "method_metadata_incomplete"
    assert missing_method["eligible_for_exact_identity"] is False


def test_connectivity_only_rt_mismatch_is_audit_only(tmp_path):
    path = tmp_path / "rt.csv"
    _write_rt_standards(path, stereochemistry_capable=False)
    load_rt_standards.cache_clear()
    selected = _candidate(
        "isomer A variant",
        IK_A_VARIANT,
        [[100.0, 1.0]],
        spectrum_id="a-variant",
    )
    selected["rt_evidence"] = match_rt_evidence(
        selected,
        3.0,
        path,
        query_method_id="LC_METHOD_A",
    )
    selected["top_candidates"] = [copy.deepcopy(selected)]

    assert selected["rt_evidence"]["identity_scope"] == "connectivity"
    assert selected["rt_evidence"]["chromatographic_status"] == "mismatch"
    result = annotate(selected)
    assert result["rt_warning"] is None
    assert result["exact_structure_confidence"] == "High"


def test_achiral_rt_does_not_resolve_same_connectivity_variant(tmp_path):
    selected = _candidate("isomer A", IK_A, [[100.0, 1.0]], spectrum_id="a")
    competitor = _candidate(
        "isomer A stereo",
        IK_A_VARIANT,
        [[110.0, 1.0]],
        spectrum_id="a-stereo",
    )
    competitor.update({"cosine": 0.94, "final_score": 0.69, "rank_score": 0.69})

    achiral_path = tmp_path / "achiral.csv"
    _write_rt_standards(achiral_path, stereochemistry_capable=False)
    load_rt_standards.cache_clear()
    achiral = copy.deepcopy(selected)
    achiral["rt_evidence"] = match_rt_evidence(
        achiral,
        2.0,
        achiral_path,
        query_method_id="LC_METHOD_A",
    )
    achiral["top_candidates"] = [copy.deepcopy(achiral), competitor]
    achiral_result = annotate(achiral)
    assert achiral_result["isomer_ambiguous"] is True
    assert achiral_result["rt_confirmed"] is False
    assert achiral_result["exact_structure_confidence"] == "Medium"
    assert achiral_result["isomer_resolution_status"] == (
        "same_connectivity_requires_chiral_or_orthogonal_standard"
    )

    chiral_path = tmp_path / "chiral.csv"
    _write_rt_standards(chiral_path, stereochemistry_capable=True)
    load_rt_standards.cache_clear()
    chiral = copy.deepcopy(selected)
    chiral["rt_evidence"] = match_rt_evidence(
        chiral,
        2.0,
        chiral_path,
        query_method_id="LC_METHOD_A",
    )
    chiral["top_candidates"] = [copy.deepcopy(chiral), competitor]
    chiral_result = annotate(chiral)
    assert chiral_result["rt_confirmed"] is True
    assert chiral_result["exact_structure_confidence"] == "High"
    assert chiral_result["identification_scope"] == "exact_structure"


def test_instrument_family_normalization_and_mzml_extraction(tmp_path):
    assert normalize_instrument_family("Q-TOF") == "qtof"
    assert normalize_instrument_family("quadrupole time-of-flight") == "qtof"
    assert normalize_instrument_family("Orbitrap Exploris 240") == "orbitrap"

    mzml = tmp_path / "instrument.mzML"
    mzml.write_text(
        '<?xml version="1.0" encoding="UTF-8"?>'
        '<mzML xmlns="http://psi.hupo.org/ms/mzml">'
        '<instrumentConfigurationList count="1">'
        '<instrumentConfiguration id="IC1">'
        '<cvParam accession="MS:1000126" name="quadrupole time-of-flight"/>'
        '<cvParam accession="MS:1000000" name="Agilent 6545 Q-TOF"/>'
        '</instrumentConfiguration>'
        '</instrumentConfigurationList>'
        '</mzML>',
        encoding="utf-8",
    )
    metadata = extract_mzml_instrument_metadata(mzml)
    assert metadata["instrument_type"] == "qtof"
    assert metadata["instrument"] == "Agilent 6545 Q-TOF"


def test_replicated_ms2_fragments_with_unreported_source_stage_remain_incomplete():
    candidates = _replicated_constitutional_candidates()
    for candidate in candidates:
        candidate.pop("lib_source_ms_level")
    before_count = len(candidates)
    detail = _run_isomer_evidence(candidates, [[100.0, 1.0], [120.0, 0.8], [150.0, 0.5]])
    assert len(candidates) == before_count
    assert detail["reference_replication"]["status"] == "complete"
    assert detail["metadata_comparability"] == "incomplete_or_mismatched"
    assert detail["eligible_for_exact_identity_promotion"] is False
