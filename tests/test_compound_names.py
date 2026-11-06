import json

from matching.compound_annotations import annotate_compound_dict
from matching.compound_names import (
    display_compound_name,
    is_compact_code_name,
    is_identifier_name,
    is_molecular_formula_name,
    is_readable_compound_name,
    normalize_match_compound_names,
    resolve_compound_name,
)
from matching.search import (
    _pathway_adduct_status,
    _pathway_polarity_status,
    _select_best_result,
)


def test_exact_inchikey_alias_restores_pinoresinol_name():
    raw = (
        "NCGC00347700-02!2-[4-[(3S,3aR,6S,6aR)-6-[3-methoxy-4-"
        "[3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxyphenyl]"
    )
    assert display_compound_name(raw, "ZJSJQWDXAYNLNS-XJUVHSKASA-N") == (
        "Pinoresinol diglucoside"
    )


def test_exact_connectivity_alias_restores_generic_pinoresinol_name():
    display, status = resolve_compound_name(
        "nan", "HGXBRUKMWQGOIE-UHFFFAOYSA-N"
    )
    assert display == "Pinoresinol"
    assert status == "exact_inchikey_alias"


def test_targeted_and_ms2_names_align_by_exact_inchikey():
    assert display_compound_name(
        "ReSpect:PM019124 Isoliquiritin|library alias",
        "YNWXJFQOCHMPCK-LXGDFETPSA-N",
    ) == "Isoliquiritin"
    assert display_compound_name(
        "Formononetin (Biochanin B)",
        "HKQYGTCOTHHOMP-UHFFFAOYSA-N",
    ) == "Formononetin"


def test_diisobutyl_phthalate_alias_and_cached_headline_are_idempotent():
    inchikey = "MGWAVDBGNNKXQV-UHFFFAOYSA-N"
    display, status = resolve_compound_name("nan", inchikey)
    assert display == "Diisobutyl phthalate"
    assert status == "exact_inchikey_alias"

    nested = {
        "compound_name": "Unresolved compound (MGWAVDBGNNKXQV)",
        "source_compound_name": "MGWAVDBGNNKXQV",
        "inchikey": inchikey,
    }
    derived = {
        "compound_name": "Unresolved compound (MGWAVDBGNNKXQV)",
        "inchikey": inchikey,
    }
    match = {
        **nested,
        "honest_headline": (
            "化合物可能为 Unresolved compound "
            "(Unresolved compound (MGWAVDBGNNKXQV)) / Dibutyl phthalate"
        ),
        "top_candidates": [
            {
                **nested,
                "honest_headline": (
                    "候选 Unresolved compound "
                    "(Unresolved compound (MGWAVDBGNNKXQV))"
                ),
                "same_formula_alternatives": [{**derived}],
            }
        ],
        "isomer_set": [{**derived}],
        "same_formula_alternatives": [{**derived}],
    }

    normalize_match_compound_names(match)
    normalized_once = match["honest_headline"]
    normalize_match_compound_names(match)

    assert match["compound_name"] == "Diisobutyl phthalate"
    assert match["honest_headline"] == normalized_once
    assert "Unresolved compound" not in match["honest_headline"]
    assert match["top_candidates"][0]["compound_name"] == display
    assert (
        match["top_candidates"][0]["same_formula_alternatives"][0][
            "compound_name"
        ]
        == display
    )
    assert "Unresolved compound" not in (
        match["top_candidates"][0]["honest_headline"]
    )
    assert match["isomer_set"][0]["compound_name"] == display
    assert match["same_formula_alternatives"][0]["compound_name"] == display
    assert "Unresolved compound" not in json.dumps(match, ensure_ascii=False)


def test_isomucronulatol_glucoside_exact_alias_repairs_cas_only_history():
    from rdkit import Chem

    inchikey = "SXHOGLPTLQBGDO-ZPHFBNLKSA-N"
    display, status = resolve_compound_name("94367-43-8", inchikey)
    assert display == "Isomucronulatol 7-O-glucoside"
    assert status == "exact_inchikey_alias"

    match = {
        "compound_name": "Unresolved compound (SXHOGLPTLQBGDO)",
        "source_compound_name": "94367-43-8",
        "inchikey": inchikey,
        "smiles": (
            "COc1ccc(C2COc3cc(O[C@@H]4O[C@H](CO)[C@@H](O)"
            "[C@H](O)[C@H]4O)ccc3C2)c(O)c1OC"
        ),
        "formula": "C23H28O10",
        "confidence": "Unknown",
        "top_candidates": [
            {
                "compound_name": "Unresolved compound (SXHOGLPTLQBGDO)",
                "source_compound_name": "94367-43-8",
                "inchikey": inchikey,
            }
        ],
    }
    normalize_match_compound_names(match)
    annotate_compound_dict(match)

    molecule = Chem.MolFromSmiles(match["smiles"])
    assert molecule is not None
    assert Chem.InchiToInchiKey(Chem.MolToInchi(molecule)) == inchikey
    assert match["compound_name"] == display
    assert match["compound_name_cn"] == "黄芪异黄烷苷"
    assert match["source_compound_name"] == "94367-43-8"
    assert match["compound_role"] == "curated_natural"
    assert match["natural_product_candidate"] is True
    assert match["confidence"] == "Unknown"
    assert match["top_candidates"][0]["compound_name"] == display


def test_unknown_identifier_headline_normalization_is_idempotent():
    match = {
        "compound_name": "AAAAAAAAAAAAAA",
        "inchikey": "AAAAAAAAAAAAAA-BBBBBBBBBB-C",
        "honest_headline": "候选 AAAAAAAAAAAAAA",
    }
    normalize_match_compound_names(match)
    normalized_once = match["honest_headline"]
    normalize_match_compound_names(match)

    assert match["honest_headline"] == normalized_once
    assert normalized_once.count("Unresolved compound") == 1


def test_formula_only_source_titles_never_pose_as_compound_names():
    raw = "NCGC00380481-01_C23H28O13S_ M+NH4"
    display, status = resolve_compound_name(
        raw, "MZHPYFUVDWYNQP-MIDFUKFWSA-N"
    )
    assert display == "Unresolved compound (MZHPYFUVDWYNQP)"
    assert status == "identifier_unresolved"

    match = {
        "compound_name": raw,
        "inchikey": "MZHPYFUVDWYNQP-MIDFUKFWSA-N",
    }
    normalize_match_compound_names(match)
    assert match["compound_name"] == "Unresolved compound (MZHPYFUVDWYNQP)"
    assert match["source_compound_name"] == raw

    for formula in ("C23H28O13S", "C6H12O6", "NaCl", "C21H19O12-"):
        assert is_molecular_formula_name(formula)
        assert is_identifier_name(formula)
        assert not is_readable_compound_name(formula)

    for readable in ("Vitamin B12", "Pinoresinol", "3-Methylhistidine"):
        assert not is_molecular_formula_name(readable)


def test_ncg_identifier_is_stripped_when_no_exact_alias_exists():
    assert display_compound_name(
        "NCGC99999999-01!NAPROXEN", "AAAAAAAAAAAAAA-BBBBBBBBBB-C"
    ) == "NAPROXEN"


def test_mls_identifier_prefix_and_trailing_cas_never_leak():
    display, status = resolve_compound_name(
        '"MLS000863595-01!(E)-3-(4-hydroxyphenyl)prop-2-en-1-one"',
        "YQHMWTPYORBCMF-ZZXKWVIFSA-N",
    )
    assert display == "(E)-3-(4-hydroxyphenyl)prop-2-en-1-one"
    assert status == "registry_payload"
    assert "MLS" not in display

    assert display_compound_name(
        "MLS002207021-01!Sinomenine115-53-7",
        "AAAAAAAAAAAAAA-BBBBBBBBBB-C",
    ) == "Sinomenine"


def test_mls_identifier_or_cas_only_payload_uses_safe_fallback():
    for raw in ("MLS000574953-01", "MLS001077311-01!2062-78-4"):
        display, status = resolve_compound_name(raw, "AAAAAAAAAAAAAA-BBBBBBBBBB-C")
        assert display == "Unresolved compound (AAAAAAAAAAAAAA)"
        assert status == "identifier_unresolved"
        assert "MLS" not in display


def test_raw_library_name_is_preserved_for_traceability():
    match = {
        "compound_name": "NCGC99999999-01!Readable chemical name",
        "inchikey": "AAAAAAAAAAAAAA-BBBBBBBBBB-C",
        "top_candidates": [],
    }
    normalize_match_compound_names(match)
    assert match["compound_name"] == "Readable chemical name"
    assert match["source_compound_name"] == "NCGC99999999-01!Readable chemical name"


def test_nist_spectral_match_title_shows_only_chemical_name_and_keeps_source():
    raw = "Spectral Match to Dibutyl phthalate from NIST14 M+H"
    display, status = resolve_compound_name(
        raw, "DOIRQSBPFJWKBE-UHFFFAOYSA-N"
    )
    assert display == "Dibutyl phthalate"
    assert status == "nist_spectral_match_title"

    match = {
        "compound_name": raw,
        "inchikey": "DOIRQSBPFJWKBE-UHFFFAOYSA-N",
        "top_candidates": [],
    }
    normalize_match_compound_names(match)
    assert match["compound_name"] == "Dibutyl phthalate"
    assert match["source_compound_name"] == raw


def test_legacy_honest_headline_is_normalized_with_main_name():
    raw = '"MLS000863595-01!(E)-3-(4-hydroxyphenyl)prop-2-en-1-one"'
    match = {
        "compound_name": raw,
        "inchikey": "YQHMWTPYORBCMF-ZZXKWVIFSA-N",
        "honest_headline": f"flavonoid (类级High) · 化合物可能为 {raw} / Chalcone",
        "top_candidates": [],
        "isomer_set": [],
    }
    normalize_match_compound_names(match)
    assert "MLS" not in match["compound_name"]
    assert "MLS" not in match["honest_headline"]
    assert "(E)-3-(4-hydroxyphenyl)" in match["honest_headline"]


def test_registry_identifiers_never_leak_when_alias_is_missing():
    for raw in (
        "2034542-40-8", "DTXSID7020182", "HMDB0000122", "LMFA01010001", "CID 123",
        "22341", "Z1575913012", "AKOS034798514", "KWIUHFFTVRNATP",
    ):
        display, status = resolve_compound_name(raw, "AAAAAAAAAAAAAA-BBBBBBBBBB-C")
        assert display == "Unresolved compound (AAAAAAAAAAAAAA)"
        assert status == "identifier_unresolved"
        assert is_identifier_name(raw)
    assert display_compound_name(
        "9999999-99-9", "ZZZZZZZZZZZZZZ-UHFFFAOYSA-N"
    ) == "Unresolved compound (ZZZZZZZZZZZZZZ)"


def test_identifier_like_name_without_inchikey_uses_safe_generic_fallback():
    for raw in (
        "1", "42", "3.14", "981271", "Z1183983796", "AKOS027472506",
        "AAWZDTNXLSGCEK",
    ):
        display, status = resolve_compound_name(raw, "")
        assert display == "Unresolved library compound"
        assert status == "identifier_unresolved"


def test_exact_curated_alias_replaces_cas_and_chembl_identifiers():
    assert display_compound_name("2034542-40-8", "FJMZMSNRMNFIHS-UHFFFAOYSA-N").startswith("N-(2-")
    assert "Dihydroxyphenyl" in display_compound_name(
        "CHEMBL4752551 [M+FA]-", "LMUQMEDGLAEQEK-UHFFFAOYSA-N"
    )


def test_exact_ao4703_alias_replaces_source_code():
    display, status = resolve_compound_name(
        "AO4703", "VMZVBRIIHDRYGK-UHFFFAOYSA-N"
    )
    assert display == (
        "2,6-Di-tert-butyl-4-(dimethylaminomethyl)phenol (AO4703)"
    )
    assert status == "exact_inchikey_alias"

    match = {
        "compound_name": "AO4703",
        "inchikey": "VMZVBRIIHDRYGK-UHFFFAOYSA-N",
    }
    normalize_match_compound_names(match)
    assert match["compound_name"].startswith("2,6-Di-tert-butyl")
    assert match["source_compound_name"] == "AO4703"


def test_unverified_compact_codes_do_not_pose_as_chemical_names():
    for raw in ("AO4703", "AZD5153", "MMV123456", "PF-07321332"):
        display, status = resolve_compound_name(
            raw, "AAAAAAAAAAAAAA-BBBBBBBBBB-C"
        )
        assert display == "Unresolved compound (AAAAAAAAAAAAAA)"
        assert status == "code_identifier_unresolved"
        assert is_compact_code_name(raw)
        assert is_identifier_name(raw)

    # Avoid overmatching ordinary systematic names and familiar digit-leading
    # chemical abbreviations.
    for readable in (
        "ANILINE-2-SULFONATE",
        "BUTYL2-FURYLPROPIONATE",
        "CIS-3-HEXENOL",
        "ANTI-13-TERT-BUTYLPHENOL",
        "BETA-3-PHENYLTROPIDINE",
        "SYN-15-TERT-BUTYLPHENOL",
        "6PPD",
    ):
        assert not is_compact_code_name(readable)
        assert display_compound_name(
            readable, "AAAAAAAAAAAAAA-BBBBBBBBBB-C"
        ) == readable


def test_name_only_lead_quality_rejects_noise_and_accessions():
    assert is_readable_compound_name("Apigenin")
    assert not is_readable_compound_name("$:00 ms2")
    assert not is_readable_compound_name("CHEMBL1873725 CollisionEnergy:102040 M+Na")


def test_name_only_library_records_are_not_loaded_in_default_production_mode():
    from config.settings import LOAD_NAME_ONLY_LEADS

    assert LOAD_NAME_ONLY_LEADS is False


def test_prestore_pathway_must_match_query_polarity():
    positive_pathway = {"polarities": ["pos"], "diagnostic_ions": [{"mz": 137.06}]}
    assert _pathway_polarity_status(positive_pathway, "NEG") == ("mismatch", ["pos"])
    assert _pathway_polarity_status(positive_pathway, "POS") == ("available", ["pos"])
    assert _pathway_polarity_status({"diagnostic_ions": []}, "NEG") == (
        "unknown_polarity",
        [],
    )


def test_prestore_pathway_must_not_mix_or_cross_specific_adducts():
    protonated = {"adducts": ["[M+H]+"]}
    assert _pathway_adduct_status(protonated, "[M+H]+", True) == (
        "available",
        ["[M+H]+"],
    )
    assert _pathway_adduct_status(protonated, "[M+Na]+", True)[0] == "mismatch"
    assert _pathway_adduct_status(
        {"adducts": ["[M+H]+", "[M+Na]+"]}, "[M+H]+", True
    )[0] == "mixed_pathway_adducts"
    assert _pathway_adduct_status(protonated, "", False)[0] == "unknown_query_adduct"


def test_exact_name_structured_candidate_beats_near_tied_name_only_lead():
    lead = {
        "compound_name": "naringenin", "source_compound_name": "naringenin [M-H]-",
        "name_only_lead": True, "has_structure": False, "smiles": "",
        "rank_score": 0.9834, "final_score": 0.9234, "cosine": 0.9915,
        "mass_ppm": 3.14, "dataset_source": "ALL_GNPS", "lib_spectrum_id": "lead",
    }
    structured = {
        "compound_name": "Naringenin", "name_only_lead": False, "has_structure": True,
        "smiles": "O=C1CC(c2ccc(O)cc2)Oc2cc(O)cc(O)c12",
        "inchikey": "FTVWIRXFELQLPI-UHFFFAOYSA-N", "formula": "C15H12O5",
        "rank_score": 0.9827, "final_score": 0.9227, "cosine": 0.9930,
        "mass_ppm": 3.87, "dataset_source": "MassBank_NIST", "lib_spectrum_id": "structured",
    }

    selected = _select_best_result([lead, structured])

    assert selected is structured
    assert selected["structure_promoted_over_name_only"] is True
    assert selected["supporting_name_only_lead"]["lib_spectrum_id"] == "lead"


def test_name_only_lead_does_not_borrow_structure_from_different_name():
    lead = {
        "compound_name": "Naringenin", "name_only_lead": True, "has_structure": False,
        "rank_score": 0.90, "final_score": 0.90, "cosine": 0.99, "mass_ppm": 3.0,
    }
    chalcone = {
        "compound_name": "Naringenin chalcone", "name_only_lead": False,
        "has_structure": True, "smiles": "structured", "rank_score": 0.899,
        "final_score": 0.899, "cosine": 0.99, "mass_ppm": 3.0,
    }
    assert _select_best_result([lead, chalcone]) is lead
