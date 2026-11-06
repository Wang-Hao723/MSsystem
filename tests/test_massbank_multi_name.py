"""Tests for the MassBank '|'-concatenated multi-field name handler.

MassBank imports join accession + several CH$NAME fields with '|'; the handler
in ``compound_names._extract_massbank_multi_name`` recovers the single most
specific readable name instead of showing the whole concatenation.
"""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from matching.compound_names import resolve_compound_name  # noqa: E402


def test_ets00001_picks_specific_name_and_flags_tentative():
    raw = ("MassbankEU:ETS00001 C10-LAS (TENTATIVE)|C10-linear alkylbenzyl "
           "sulfonate|4-(decan-4-yl)benzenesulfonic acid [M-H]-")
    display, status = resolve_compound_name(raw, "NANHIUZYPFDGJS-UHFFFAOYSA-N")
    assert display == "4-(decan-4-yl)benzenesulfonic acid"
    assert status == "putative_extracted"  # (TENTATIVE) must downgrade confidence


def test_trailing_adduct_is_stripped_from_display():
    raw = ("MassbankEU:ETS00001 C10-LAS (TENTATIVE)|C10-linear alkylbenzyl "
           "sulfonate|4-(decan-4-yl)benzenesulfonic acid [M-H]-")
    display, _ = resolve_compound_name(raw, "NANHIUZYPFDGJS-UHFFFAOYSA-N")
    assert "[M-H]-" not in display
    assert "|" not in display
    assert "Massbank" not in display


def test_sulfamethazine_pipe_keeps_systematic_name():
    raw = ("Massbank:AU100802 Sulfamethazine|4-amino-N-(4,6-dimethylpyrimidin-2-yl)"
           "benzenesulfonamide [M+H]+")
    display, status = resolve_compound_name(raw, "")
    assert display == "4-amino-N-(4,6-dimethylpyrimidin-2-yl)benzenesulfonamide"
    assert status == "massbank_multi_name"


def test_standard_mix_not_treated_as_putative():
    raw = ("MassbankEU:ETS00014 C11-LAS (STANDARD MIX)|C11-linear alkylbenzyl "
           "sulfonate|4-(undecan-5-yl)benzenesulfonic acid [M-H]-")
    display, status = resolve_compound_name(raw, "FERBTPHUEYEGDN-UHFFFAOYSA-N")
    assert display == "4-(undecan-5-yl)benzenesulfonic acid"
    assert status == "massbank_multi_name"  # STANDARD MIX is not tentative


def test_plain_name_unaffected():
    display, status = resolve_compound_name("Acetaminophen", "")
    assert display == "Acetaminophen"
    assert status == "source_name"


def test_pipe_name_without_massbank_prefix_not_caught():
    # GNPS duplicate-style or other '|'-joins without a MassBank accession prefix
    # must not be silently rewritten by the MassBank handler.
    display, status = resolve_compound_name("Cyclohomonervilasterol|Cyclohomonervilasterol", "")
    assert status != "massbank_multi_name"


def test_no_pipe_massbank_name_not_caught():
    # Scope is intentionally the pipe-concatenated format; a clean single
    # MassBank accession+name without '|' is left to the existing path.
    display, status = resolve_compound_name("Massbank:AU100802 Sulfamethazine", "")
    assert status != "massbank_multi_name"


def test_all_unsafe_segments_falls_through():
    # Every segment is a formula or identifier -> no readable name -> handler
    # returns None and the existing path takes over (whatever it produces).
    # The contract here is only that the MassBank handler does not fabricate a
    # massbank_multi_name result from unsafe segments.
    raw = "Massbank:XX0001 C9H18O2|12345|AKOS0000001"
    display, status = resolve_compound_name(raw, "")
    assert status != "massbank_multi_name"
