"""Lock the recovery of chemical names wrapped in library provenance.

These cases used to collapse to ``Unresolved compound (<InChIKey>)`` because the
readable name was embedded inside a dereplicator / putative / SIRIUS provenance
string.  They must now display the real name (or an honest computational
-candidate placeholder) without leaking registry hashes or unconfirmed guesses.
"""

from matching.compound_names import resolve_compound_name


def test_dereplicator_wrapped_name_is_recovered():
    # The user-reported case: name was present in the data all along, wrapped in
    # "Dereplicator Identification - X_X M-H" provenance.
    raw = "Dereplicator Identification - Cyclohomonervilasterol_Cyclohomonervilasterol M-H"
    display, status = resolve_compound_name(raw, "QECYIHDHCKWWPU-MJUSBNMBSA-N")
    assert display == "Cyclohomonervilasterol"
    assert "Unresolved" not in display
    assert status != "missing"


def test_putative_parenthetical_keeps_core_name():
    display, status = resolve_compound_name(
        "Sagitol B (Putative Analogue)  M+H", ""
    )
    assert display == "Sagitol B"
    assert status == "putative_extracted"


def test_putative_trailing_word_and_unk_are_stripped():
    display, _ = resolve_compound_name("Surfactin (Putative) Unk", "")
    assert display == "Surfactin"


def test_putative_observed_from_parenthetical_is_stripped():
    display, _ = resolve_compound_name(
        "Aflatoxin B1 (putative_Observed from A. flavus)  M+Na", ""
    )
    assert display == "Aflatoxin B1"


def test_sirius_candidate_is_never_presented_as_a_confirmed_name():
    # An unconfirmed SIRIUS formula annotation must NOT leak "Carnitine" as if
    # it were an identification.  It gets an honest placeholder instead.
    raw = (
        '"Candidate Carnitine (delta mass -1.9793; '
        'sirius atomic difference: 1C, 2H, -1O)" M+H'
    )
    display, status = resolve_compound_name(raw, "QECYIHDHCKWWPU-MJUSBNMBSA-N")
    assert status == "computational_candidate"
    assert "Carnitine" not in display
    assert "candidate" in display.lower()


def test_sirius_candidate_without_inchikey_stays_honest():
    raw = '[BA_core: 3,12a-(OH)2] Candidate Dihydroxylated bile acid M+H'
    display, status = resolve_compound_name(raw, "")
    assert status == "computational_candidate"
    assert display == "Computational candidate"


def test_dereplicator_duplicated_payload_collapses_once():
    display, _ = resolve_compound_name(
        "Dereplicator Identification - Citrusin_VI_Citrusin_VI M+H", ""
    )
    # Underscores become spaces and the duplicated payload is collapsed.
    assert display.lower() == "citrusin vi"
    assert display.lower().count("citrusin") == 1


def test_trailing_adduct_still_stripped_from_plain_name():
    # Regression guard: the existing adduct stripping must keep working.
    display, _ = resolve_compound_name("Naringenin M+H", "FTVWIRXFELQLPI-UHFFFAOYSA-N")
    assert display == "Naringenin"


def test_plain_name_ending_in_m_word_is_not_mistaken_for_an_adduct():
    display, status = resolve_compound_name("glycerol monolinoleate", "")
    assert display == "glycerol monolinoleate"
    assert status == "source_name"


def test_bracketed_and_multimer_adduct_suffixes_still_strip():
    assert resolve_compound_name("Naringenin [M+H]+", "")[0] == "Naringenin"
    assert resolve_compound_name("Naringenin 2M+Na", "")[0] == "Naringenin"
    assert resolve_compound_name("Magnoflorine M", "")[0] == "Magnoflorine"
    assert resolve_compound_name("Magnoflorine [M]+", "")[0] == "Magnoflorine"


def test_exact_inchikey_alias_cannot_leak_adduct_into_compound_name():
    # This exact cached alias previously rendered ``Schaftoside [M+H]+`` in a
    # NEG targeted result.  Ion form belongs to the adduct field, never the
    # chemical-name field, even when the alias came from an exact InChIKey.
    display, status = resolve_compound_name(
        "nan", "MMDUKUSNQNWVET-UHFFFAOYSA-N"
    )
    assert display == "Schaftoside"
    assert status == "exact_inchikey_alias"


def test_metlin_spectral_match_title_shows_only_chemical_name():
    display, status = resolve_compound_name(
        "Spectral Match to Glu Phe from METLIN M+H", ""
    )
    assert display == "Glu Phe"
    assert status == "metlin_spectral_match_title"
