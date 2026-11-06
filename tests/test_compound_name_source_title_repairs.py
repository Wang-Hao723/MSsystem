from matching.compound_names import resolve_compound_name


def test_exact_pubchem_verified_alias_repairs_cyclohomonervilasterol():
    raw = (
        "Dereplicator Identification - "
        "Cyclohomonervilasterol_Cyclohomonervilasterol M-H"
    )
    display, status = resolve_compound_name(
        raw, "QECYIHDHCKWWPU-CMDGGOBGSA-N"
    )
    assert display == "Cyclohomonervilasterol"
    assert status == "exact_inchikey_alias"


def test_dereplicator_title_recovers_embedded_name_but_keeps_tentative_status():
    display, status = resolve_compound_name(
        "Dereplicator Identification - Achyroclinopyrone_Achyroclinopyrone M+H",
        "AAAAAAAAAAAAAA-BBBBBBBBBB-C",
    )
    assert display == "Achyroclinopyrone"
    assert status == "dereplicator_annotation_title"


def test_dereplicator_derivative_sentence_is_not_promoted_to_parent_identity():
    display, status = resolve_compound_name(
        "Suspect related to Dereplicator Identification - rhapsamine "
        "(predicted molecular formula SIRIUS: C31H62N2O6) [M+H]+",
        "AAAAAAAAAAAAAA-BBBBBBBBBB-C",
    )
    assert display == "Unresolved compound (AAAAAAAAAAAAAA)"
    assert status.startswith("unsafe_source_name_")


def test_quoted_collision_energy_export_recovers_the_chemical_title():
    display, status = resolve_compound_name(
        '"(22S,23S)-28-Homo Brassinolide CollisionEnergy:102040" 2M+Na',
        "AAAAAAAAAAAAAA-BBBBBBBBBB-C",
    )
    assert display == "(22S,23S)-28-Homo Brassinolide"
    assert status == "collision_energy_title"


def test_collision_energy_parser_does_not_hide_corrupt_stereochemistry():
    display, status = resolve_compound_name(
        '"(3,4-dimethoxyphenyl)-N-(???methyl)carboxamide '
        'CollisionEnergy:102040" M+H',
        "AAAAAAAAAAAAAA-BBBBBBBBBB-C",
    )
    assert display == "Unresolved compound (AAAAAAAAAAAAAA)"
    assert status == (
        "unsafe_collision_title_unresolved_stereochemistry_or_corrupt_text"
    )
