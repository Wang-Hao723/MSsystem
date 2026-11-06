"""Quality-gate regressions for the expert fragmentation knowledge bases."""

from dataclasses import replace
import json
from pathlib import Path

from rdkit import Chem

from matching.expert_scorer import ChemicalExpertScorer
from matching.expert_v3 import _polarity as v3_polarity
from matching.fragment_feature_library import (
    COMPOUND_CLASS_PROFILES,
    get_fragments_by_class,
    match_fragment_to_class,
)
from matching.fragmentation_rules import (
    DIAGNOSTIC_IONS,
    FRAGMENTATION_RULES,
    compute_rule_coverage,
    fragmentation_rule_applicable,
    get_neutral_loss_by_mass,
    match_diagnostic_ion,
)
from matching.rule_quality import formula_mz_consistent, normalize_adduct
from mssystem.classification import RuleBasedClassifier
from kb_v3.loader import get_family_chemistry


PC_SMILES = "CCCCC(=O)OCC(COP(=O)([O-])OCC[N+](C)(C)C)OC(=O)CCCC"
ADENOSINE_SMILES = "Nc1ncnc2c1ncn2[C@@H]1O[C@H](CO)[C@@H](O)[C@H]1O"
DEOXYADENOSINE_SMILES = "Nc1ncnc2c1ncn2[C@@H]1O[C@H](CO)[C@@H](O)C1"
FIVE_PRIME_DEOXYADENOSINE_SMILES = "C[C@H]1O[C@H]([C@H](O)[C@@H]1O)N1C=NC2=C(N)N=CN=C12"
CAFFEINE_SMILES = "Cn1c(=O)c2c(ncn2C)n(C)c1=O"
PE_SMILES = "CCCC(=O)OCC(COP(=O)(O)OCCN)OC(=O)CCC"
C_GLYCOSIDE_SMILES = "c1ccccc1C2OC(CO)C(O)C(O)C2O"
O_GLYCOSIDE_SMILES = "c1ccccc1OC2OC(CO)C(O)C(O)C2O"


def test_every_runtime_diagnostic_ion_passes_formula_mass_gate():
    failures = {
        key: (ion.mz, ion.formula)
        for key, ion in DIAGNOSTIC_IONS.items()
        if not formula_mz_consistent(ion.formula, ion.mz)
    }
    assert failures == {}


def test_invalid_legacy_fragment_is_retained_for_audit_but_quarantined():
    raw = COMPOUND_CLASS_PROFILES["flavonoid"].characteristic_fragments
    assert any(abs(fragment.mz - 309.0974) < 1e-6 for fragment in raw)
    production = get_fragments_by_class("flavonoid")
    assert not any(abs(fragment.mz - 309.0974) < 1e-6 for fragment in production)


def test_diagnostic_ions_never_cross_query_polarity():
    positive = match_diagnostic_ion(
        93.0346, 0.05, classification="phenolic_acid",
        strict_fg=True, polarity="POS",
    )
    negative = match_diagnostic_ion(
        91.0548, 0.05, classification="benzenoid",
        strict_fg=True, polarity="NEG",
    )
    assert positive == []  # phenoxide is a negative-ion rule
    assert negative == []  # tropylium is a positive-ion rule


def test_phosphocholine_rule_requires_candidate_scope_support():
    unsupported = match_diagnostic_ion(
        184.0733, 0.01, classification="lipid",
        strict_fg=True, polarity="POS",
    )
    supported = match_diagnostic_ion(
        184.0733, 0.01, classification="lipid",
        strict_fg=True, polarity="POS", candidate_smiles=PC_SMILES,
        query_adduct="[M+H]+",
    )
    wrong_adduct = match_diagnostic_ion(
        184.0733, 0.01, classification="lipid",
        strict_fg=True, polarity="POS", candidate_smiles=PC_SMILES,
        query_adduct="[M+Na]+",
    )
    assert unsupported == []
    assert [ion.formula for ion in supported] == ["C5H15NO4P+"]
    assert wrong_adduct == []


def test_feature_library_filters_polarity_and_candidate_scope():
    assert match_fragment_to_class(255.2329, 0.01, polarity="POS") == []
    assert any(
        hit["class_name"] == "lipid"
        for hit in match_fragment_to_class(255.2329, 0.01, polarity="NEG")
    )

    generic_lipid = get_fragments_by_class(
        "lipid", polarity="POS", enforce_scope=True,
        raw_classification="lipid",
    )
    pc_lipid = get_fragments_by_class(
        "lipid", polarity="POS", enforce_scope=True,
        candidate_smiles=PC_SMILES, query_adduct="[M+H]+",
    )
    assert not any(abs(fragment.mz - 184.0739) < 0.001 for fragment in generic_lipid)
    assert any(abs(fragment.mz - 184.0739) < 0.001 for fragment in pc_lipid)


def test_rule_coverage_excludes_opposite_charge_rules_from_denominator():
    positive = compute_rule_coverage([], "fatty_acid", polarity="POS")
    negative = compute_rule_coverage([], "fatty_acid", polarity="NEG")
    assert "CHARGE_REMOTE_FRAGMENTATION" not in positive["uncovered_rules"]
    assert "CHARGE_REMOTE_FRAGMENTATION" in negative["uncovered_rules"]


def test_absence_of_isomer_rules_is_not_positive_evidence():
    result = ChemicalExpertScorer()._analyze_isomer_differentiation(
        [91.0548], "not_a_curated_class", polarity="POS"
    )
    assert result["score"] == 0.0


def test_hrms_neutral_loss_tolerance_and_multicharge_guard():
    scorer = ChemicalExpertScorer()
    # V12 publication policy uses a 0.02 Da HRMS gate; 0.05 Da was the
    # inherited legacy value and can turn unrelated accurate-mass differences
    # into named neutral losses.
    assert scorer.nl_tol == 0.02
    result = scorer._analyze_neutral_losses(
        500.0, [482.0, 450.0], classification="lipid", query_charge=2
    )
    assert result["explanation_rate"] == 0.0
    assert "multiply_charged" in result["note"]


def test_expert_v3_declared_query_polarity_overrides_missing_library_adduct():
    assert v3_polarity(None, "POS") == "pos"
    assert v3_polarity(None, "NEG") == "neg"


def test_every_runtime_smarts_guard_compiles():
    invalid = []
    for key, rule in FRAGMENTATION_RULES.items():
        for smarts in (*rule.required_smarts, *rule.required_smarts_all):
            if Chem.MolFromSmarts(smarts) is None:
                invalid.append(("rule", key, smarts))
    for key, ion in DIAGNOSTIC_IONS.items():
        for smarts in (*ion.required_smarts, *ion.required_smarts_all):
            if Chem.MolFromSmarts(smarts) is None:
                invalid.append(("ion", key, smarts))
    for class_name, profile in COMPOUND_CLASS_PROFILES.items():
        for feature in profile.characteristic_fragments:
            for smarts in (*feature.required_smarts, *feature.required_smarts_all):
                if Chem.MolFromSmarts(smarts) is None:
                    invalid.append(("feature", class_name, smarts))
    assert invalid == []


def test_nucleosides_are_not_misclassified_as_alkaloids():
    classifier = RuleBasedClassifier()
    assert classifier.predict(ADENOSINE_SMILES) == "nucleoside"
    assert classifier.predict(DEOXYADENOSINE_SMILES) == "nucleoside"
    assert classifier.predict(CAFFEINE_SMILES) == "alkaloid"


def test_ribose_and_deoxyribose_rules_are_structure_and_adduct_specific():
    ribose = FRAGMENTATION_RULES["RIBONUCLEOSIDE_GLYCOSIDIC_CLEAVAGE"]
    deoxyribose = FRAGMENTATION_RULES["DEOXYRIBONUCLEOSIDE_GLYCOSIDIC_CLEAVAGE"]
    common = {"polarity": "POS", "query_adduct": "[M+H]+"}
    assert fragmentation_rule_applicable(
        ribose, candidate_smiles=ADENOSINE_SMILES, **common
    ) == (True, "applicable")
    assert fragmentation_rule_applicable(
        ribose, candidate_smiles=DEOXYADENOSINE_SMILES, **common
    )[0] is False
    assert fragmentation_rule_applicable(
        deoxyribose, candidate_smiles=DEOXYADENOSINE_SMILES, **common
    ) == (True, "applicable")
    assert fragmentation_rule_applicable(
        deoxyribose, candidate_smiles=ADENOSINE_SMILES, **common
    )[0] is False
    assert fragmentation_rule_applicable(
        deoxyribose, candidate_smiles=FIVE_PRIME_DEOXYADENOSINE_SMILES, **common
    ) == (True, "applicable")
    assert fragmentation_rule_applicable(
        ribose, candidate_smiles=ADENOSINE_SMILES,
        polarity="POS", query_adduct="[M+Na]+",
    )[0] is False


def test_nucleobase_ion_requires_matching_structure_and_protonated_adduct():
    supported = match_diagnostic_ion(
        136.0618, 0.01, classification="nucleoside", strict_fg=True,
        polarity="POS", candidate_smiles=ADENOSINE_SMILES,
        query_adduct="[M+H]+",
    )
    no_structure = match_diagnostic_ion(
        136.0618, 0.01, classification="nucleoside", strict_fg=True,
        polarity="POS", query_adduct="[M+H]+",
    )
    wrong_adduct = match_diagnostic_ion(
        136.0618, 0.01, classification="nucleoside", strict_fg=True,
        polarity="POS", candidate_smiles=ADENOSINE_SMILES,
        query_adduct="[M+Na]+",
    )
    assert [ion.name for ion in supported] == ["质子化腺嘌呤"]
    assert no_structure == []
    assert wrong_adduct == []


def test_duplicate_peaks_do_not_duplicate_one_diagnostic_rule():
    result = ChemicalExpertScorer()._analyze_diagnostic_ions(
        [136.0618, 136.0617, 119.0350],
        classification="nucleoside",
        polarity="POS",
        candidate_smiles=ADENOSINE_SMILES,
        query_adduct="[M+H]+",
    )
    adenine_hits = [
        hit for hit in result["matched_ions"]
        if hit["ion_name"] == "质子化腺嘌呤"
    ]
    assert len(adenine_hits) == 1
    assert result["matched_count"] == 1


def test_collision_energy_conditioning_requires_explicit_comparable_unit():
    base = FRAGMENTATION_RULES["RIBONUCLEOSIDE_GLYCOSIDIC_CLEAVAGE"]
    conditioned = replace(
        base, collision_energy_min=20.0,
        collision_energy_max=40.0, collision_energy_unit="eV",
    )
    common = {
        "polarity": "POS", "candidate_smiles": ADENOSINE_SMILES,
        "query_adduct": "[M+H]+",
    }
    assert fragmentation_rule_applicable(
        conditioned, query_collision_energy=30, **common
    ) == (False, "collision_energy_not_comparable")
    assert fragmentation_rule_applicable(
        conditioned, query_collision_energy={"value": 30, "unit": "eV"}, **common
    ) == (True, "applicable")
    assert fragmentation_rule_applicable(
        conditioned, query_collision_energy={"value": 30, "unit": "NCE"}, **common
    ) == (False, "collision_energy_unit_mismatch")


def test_adduct_aliases_normalize_without_polarity_inference():
    assert normalize_adduct("[M+HCOO]-") == "[M+FA-H]-"
    assert normalize_adduct("[M+CH3COO]-") == "[M+Hac-H]-"
    assert normalize_adduct("POS") == "POS"


def test_narrow_neutral_losses_cannot_bypass_structure_and_adduct_guards():
    pe_supported = get_neutral_loss_by_mass(
        141.0191, 0.01, classification="lipid", strict_fg=True,
        polarity="POS", candidate_smiles=PE_SMILES, query_adduct="[M+H]+",
    )
    pc_rejected = get_neutral_loss_by_mass(
        141.0191, 0.01, classification="lipid", strict_fg=True,
        polarity="POS", candidate_smiles=PC_SMILES, query_adduct="[M+H]+",
    )
    unknown_structure_rejected = get_neutral_loss_by_mass(
        141.0191, 0.01, classification="lipid", strict_fg=True,
        polarity="POS", query_adduct="[M+H]+",
    )
    assert [loss.formula for loss in pe_supported] == ["C2H8NO4P"]
    assert pc_rejected == []
    assert unknown_structure_rejected == []

    c_glycoside = get_neutral_loss_by_mass(
        120.0423, 0.01, classification="glycoside", strict_fg=True,
        polarity="NEG", candidate_smiles=C_GLYCOSIDE_SMILES,
        query_adduct="[M-H]-",
    )
    o_glycoside = get_neutral_loss_by_mass(
        120.0423, 0.01, classification="glycoside", strict_fg=True,
        polarity="NEG", candidate_smiles=O_GLYCOSIDE_SMILES,
        query_adduct="[M-H]-",
    )
    assert [loss.name_en for loss in c_glycoside] == ["Sugar cross-ring loss 120"]
    assert o_glycoside == []


def test_legacy_curation_inventory_accounts_for_every_source_rule():
    root = Path(__file__).resolve().parents[1]
    inventory = json.loads(
        (root / "matching" / "fragmentation_rules_curation.json").read_text()
    )
    records = inventory["records"]
    assert len(records) == 235
    assert len({record["occurrence_id"] for record in records}) == 235
    assert inventory["meta"]["unique_source_rule_ids"] == 221
    assert any(record["status"] == "rejected_incorrect" for record in records)
    assert all(
        record["status"] == "promoted_runtime"
        for record in records
        if record["execution_gate"]["production_enabled"]
    )


def test_family_chemistry_fuzzy_lookup_does_not_raise_name_error():
    assert get_family_chemistry("Ranunculaceae") is not None
