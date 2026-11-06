import re
import shutil
import subprocess

import pytest

from web.report_en import generate_english_html
from web.server import _build_pdf_report


def _sample_results():
    pathway = {
        "n_spectra": 3,
        "polarities": ["neg"],
        "adducts": ["[M-H]-"],
        "diagnostic_ions": [
            {
                "mz": 151.0037,
                "diag_key": "flavonoid_a_ring",
                "annotation": "黄酮A环离子",
                "freq": 0.8,
                "rel_intensity": 72.0,
                "recurring": True,
            }
        ],
        "neutral_losses": [
            {
                "loss_da": 18.0106,
                "fragment_mz": 283.0450,
                "nl_key": "H2O",
                "annotation": "Water",
                "freq": 0.67,
                "rel_intensity": 35.0,
                "recurring": True,
            }
        ],
        "loss_chains": [[301.0556, 283.0450]],
    }
    match = {
        "compound_name": "Quercetin",
        "formula": "C15H10O7",
        "smiles": "O=c1cc(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12",
        "inchikey": "REFJWTPEDVJJIY-UHFFFAOYSA-N",
        "classification": "flavonoid",
        "dataset_source": "MassBank",
        "source_tier": "experimental",
        "cosine": 0.9132,
        "entropy": 0.8811,
        "reverse_match": 0.82,
        "mass_ppm": 2.1,
        "final_score": 0.7712,
        "rank_score": 0.773,
        "expert_score": 0.61,
        "confidence": "High",
        "confidence_class": "High",
        "evidence_tier": "High",
        "annotation_level": "Level 2 candidate",
        "adduct": "[M-H]-",
        "query_adduct": "[M-H]-",
        "query_adduct_reliable": True,
        "mass_match_basis": "neutral_mass",
        "pathway_json": pathway,
        "pathway_polarity_status": "available",
        "pathway_adduct_status": "available",
        "expert_details": {"pathway_align": {"score": 0.75}},
        "contrastive_expert_shadow": {
            "status": "descriptive_only",
            "matched_ions": [151.0037],
            "matched_losses": [18.0106],
        },
        "formula_isotope_evidence": {
            "status": "supporting",
            "score": 0.89,
            "observed": {"m1_over_m": 0.16, "m2_over_m": 0.03},
            "theoretical": {"m1_over_m": 0.165, "m2_over_m": 0.028},
        },
        "isotope_rank_adjustment": 0.0078,
        "ce_profile_evidence": {
            "status": "available_near_energy",
            "runtime_influence": False,
            "score_contribution": 0.0,
            "query_collision_energy": 35.0,
            "query_collision_energy_unit": "eV",
            "best_profile": {
                "closest_energy": 35.0,
                "collision_energy_unit": "eV",
                "energy_delta": 0.0,
                "instrument_type": "LC-ESI-QTOF",
                "instrument_model": "Test QTOF",
                "n_energy_levels": 5,
                "n_profile_spectra": 12,
                "quality_grade": "A",
                "consensus_cosine": 0.91,
                "consensus_entropy": 0.88,
                "matched_fragment_trajectories": [
                    {
                        "fragment_mz": 151.0037,
                        "trend": "increasing",
                        "n_levels": 4,
                        "apex_energy": 45.0,
                    }
                ],
            },
        },
        "isomer_ambiguous": True,
        "isomer_set": ["Quercetin", "Morin"],
        "isomer_fragment_score": 0.4,
        "isomer_fragment_detail": {
            "unique_supported": 1,
            "unique_total": 3,
            "rank_override_enabled": False,
        },
        "top_candidates": [
            {
                "compound_name": "Quercetin",
                "formula": "C15H10O7",
                "cosine": 0.9132,
                "entropy": 0.8811,
                "final_score": 0.7712,
                "mass_ppm": 2.1,
            },
            {
                "compound_name": "Morin",
                "formula": "C15H10O7",
                "cosine": 0.902,
                "entropy": 0.86,
                "final_score": 0.76,
                "mass_ppm": 2.5,
            },
        ],
    }
    return [
        {
            "spectrum_id": "scan=1",
            "precursor_mz": 301.0556,
            "query_peaks": [[151.0037, 100.0], [179.0, 44.0]],
            "lib_peaks": [[151.0040, 88.0], [180.0, 30.0]],
            "isotope_evidence_status": "observed",
            "precursor_isolation_evidence": {"purity": 0.91},
            "match": match,
        },
        {
            "spectrum_id": "scan=2",
            "precursor_mz": 250.0,
            "query_peaks": [[100.0, 10.0]],
            "match": None,
        },
    ]


def test_detailed_english_report_is_structured_and_system_prose_is_english():
    rendered = generate_english_html(
        _sample_results(), "sample.mzML", detailed=True
    )

    assert '<html lang="en">' in rendered
    assert "LC–MS/MS Identification Report" in rendered
    assert "Library-derived fragmentation signature" in rendered
    assert "not atom-mapped bond-cleavage mechanisms" in rendered
    assert "MS1 isotope evidence" in rendered
    assert "Adduct and neutral-mass evidence" in rendered
    assert "Multi-source CE evidence" in rendered
    assert "shadow-only" in rendered
    assert "LC-ESI-QTOF" in rendered
    assert "Isomer assessment" in rendered
    assert "Isomer-group confidence" in rendered
    assert "Unresolved isomer competition — complete offline candidate list" in rendered
    assert "Ranking lead: Quercetin" in rendered
    assert "Morin" in rendered
    assert "Flavonoid A Ring" in rendered
    assert "黄酮A环离子" not in rendered
    assert not re.search(r"[\u4e00-\u9fff]", rendered)


def test_standard_english_report_omits_heavy_evidence_panels():
    rendered = generate_english_html(
        _sample_results(), "sample.mzML", detailed=False
    )

    assert "LC–MS/MS Identification Report" in rendered
    assert "Compound evidence details" not in rendered
    assert "Library-derived fragmentation signature" not in rendered
    assert 'href="#compound-1"' not in rendered
    assert "Unresolved isomer competition — complete offline candidate list" in rendered
    assert "Ranking lead: Quercetin" in rendered
    assert "Morin" in rendered


def test_english_report_keeps_scan_counts_separate_from_unique_compounds():
    # The API can pass a deduplicated matched view plus every unmatched scan.
    # Explicit scan totals must therefore override len(rows).
    rendered = generate_english_html(
        _sample_results(),
        "sample.mzML",
        total_spectra=10,
        matched_spectra=4,
        high_supported_compounds=1,
        high_confidence_scope_groups=10,
        high_confidence_exact_structures=1,
        high_confidence_isomer_groups=9,
    )

    assert '<strong>10</strong><span>Total MS/MS spectra</span>' in rendered
    assert (
        '<strong>4</strong><span>MS/MS spectra with a library candidate '
        '(40.0%)</span>'
    ) in rendered
    assert '<strong>1</strong><span>Unique matched compounds</span>' in rendered
    assert '<strong>6</strong><span>Unmatched MS/MS spectra</span>' in rendered
    assert (
        '<strong>10</strong><span>High evidence groups<br>'
        '<small>single structures 1 · isomer groups 9</small></span>'
    ) in rendered


def test_english_pdf_is_generated_from_the_same_results(tmp_path):
    pdf = _build_pdf_report(
        {
            "filename": "sample.mzML",
            "compound_conclusion_summary": {
                "high_confidence_scope_group_count": 10,
                "high_confidence_exact_structure_count": 1,
                "high_confidence_isomer_group_count": 9,
            },
        },
        _sample_results(),
        language="en",
    )

    assert pdf.startswith(b"%PDF")
    assert len(pdf) > 1_000
    assert b"/FontFile2" in pdf, "English PDF must embed a visible TrueType font"

    pdftoppm = shutil.which("pdftoppm")
    if pdftoppm is None:
        pytest.skip("Poppler is unavailable for PDF appearance regression")
    pdf_path = tmp_path / "english-report.pdf"
    image_prefix = tmp_path / "english-report-page"
    pdf_path.write_bytes(pdf)
    subprocess.run(
        [
            pdftoppm,
            "-f", "1",
            "-singlefile",
            "-png",
            "-r", "96",
            str(pdf_path),
            str(image_prefix),
        ],
        check=True,
        capture_output=True,
    )
    from PIL import Image

    rendered = Image.open(image_prefix.with_suffix(".png")).convert("L")
    title_region = rendered.crop((0, 0, rendered.width, 150))
    visible_dark_pixels = sum(value < 225 for value in title_region.getdata())
    assert visible_dark_pixels > 100, "English PDF title rendered invisibly"


def test_export_ui_exposes_single_and_batch_english_reports():
    index = open("web/static/index.html", encoding="utf-8").read()
    script = open("web/static/js/export.js", encoding="utf-8").read()

    assert "English HTML Detailed Report" in index
    assert "English PDF Report" in index
    assert "downloadHTML('detailed','en')" in index
    assert 'value="html_detailed|en"' in index
    assert "batchExportSelected" in script
    assert "language=${encodeURIComponent(language)}" in script
    assert "const requestBody = { format, taskIds, language };" in script
