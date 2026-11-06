#!/usr/bin/env python3
"""Build the fail-closed V15 feature-fragment legacy curation overlay.

The script reads one explicitly named 2,030-row ledger.  It never discovers
or scans data/, output/, or project databases, and it never edits runtime
feature-fragment catalogs.  Every ledger item receives exactly one terminal
curation relationship.  Only proposals already backed by the frozen public
Agilent 6560 fixture are marked ``accepted_shadow``; all other proposals stay
non-executable until acquisition-conditioned target/background counts exist.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import sys
from collections import Counter
from pathlib import Path
from typing import Any, Iterable


ROOT = Path(__file__).resolve().parents[1]
LEDGER_RELATIVE = Path(
    "docs/audits/V15_FEATURE_FRAGMENT_LEGACY_PROVENANCE_LEDGER_20260826.jsonl"
)
DEFAULT_LEDGER = ROOT / LEDGER_RELATIVE
DEFAULT_OUTPUT_DIR = ROOT / "config" / "feature_fragments" / "curation"
EXPECTED_LEDGER_SHA256 = (
    "926a85ca00abebf0775205007b727dac2b641c860ba6083fb0a706aaf2e81fa3"
)
EXPECTED_LEDGER_RECORDS = 2030
CURATION_SCHEMA_VERSION = "1.0.0"
CURATION_DATE = "2026-08-26"

TERMINAL_RELATIONSHIPS = {
    "exact_alias",
    "numeric_correction",
    "split_parent",
    "superseded",
    "quarantined_duplicate",
    "quarantined_formula_mass",
    "quarantined_misplaced",
    "rejected_incorrect",
    "defer",
}
PROPOSAL_STATUSES = {
    "accepted_shadow",
    "needs_recount",
    "literature_only",
    "rejected",
}

SOURCE_FILE_HASHES = {
    "kb_v3/empirical_rules.json": (
        "988e237dda3a21d5dab080ac5bde834f43fac13cd460605cc488d03e61a2c16a"
    ),
    "kb_v3/empirical_rules_candidates_v2.json": (
        "e59ae49089a003b4ed3e5eef2158101e02f12cbca356b406b561a1c98facadb7"
    ),
    "matching/fragment_feature_library.py": (
        "f8838ccee8f6e986622cc5974cd6ad2ba5afea534e17ed6b4ca1e4d2614b3b72"
    ),
    "matching/fragmentation_rules.py": (
        "a8c0682db115bbf0b8ef151bc9800372f587c43aad4b885c910a9ba523a21e01"
    ),
    "matching/fragmentation_rules_curation.json": (
        "5e36b20cce02223316f4ad9a2756f029c99ec374820625e16b019c196a576129"
    ),
}

FORMULA_MASS_QUARANTINE_IDS = (
    "feature_profile:flavonoid:007",
    "feature_profile:flavonoid:008",
    "feature_profile:alkaloid:012",
    "feature_profile:lignan:003",
    "feature_profile:fatty_acid:001",
    "feature_profile:terpenoid:008",
    "feature_profile:steroid:007",
    "feature_profile:amino_acid:011",
    "feature_profile:lipid:008",
)

ATOMIC_MASSES = {
    "H": 1.00782503223,
    "C": 12.0,
    "N": 14.00307400443,
    "O": 15.99491461957,
    "P": 30.97376199842,
    "S": 31.9720711744,
}
ELECTRON_MASS_DA = 0.000548579909065
FORMULA_TOKEN = re.compile(r"([A-Z][a-z]?)(\d*)")


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def sha256_file(path: Path) -> str:
    hasher = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            hasher.update(block)
    return hasher.hexdigest()


def json_bytes(value: Any) -> bytes:
    return (
        json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n"
    ).encode("utf-8")


def jsonl_bytes(values: Iterable[dict[str, Any]]) -> bytes:
    lines = [
        json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"))
        for value in values
    ]
    return (("\n".join(lines) + "\n") if lines else "").encode("utf-8")


def formula_mz(formula: str) -> float:
    """Calculate a simple monoisotopic mass/m/z for the curated formulas."""
    text = formula.strip()
    charge = 0
    match = re.search(r"([+-])(\d*)$", text)
    if match:
        magnitude = int(match.group(2) or "1")
        charge = magnitude if match.group(1) == "+" else -magnitude
        text = text[: match.start()]
    total = 0.0
    position = 0
    for token in FORMULA_TOKEN.finditer(text):
        if token.start() != position or token.group(1) not in ATOMIC_MASSES:
            raise ValueError(f"unsupported formula: {formula}")
        count = int(token.group(2) or "1")
        total += ATOMIC_MASSES[token.group(1)] * count
        position = token.end()
    if position != len(text) or not text:
        raise ValueError(f"unsupported formula: {formula}")
    if charge:
        return (total - charge * ELECTRON_MASS_DA) / abs(charge)
    return total


def load_ledger(path: Path) -> list[dict[str, Any]]:
    payload = path.read_bytes()
    actual_hash = sha256_bytes(payload)
    if actual_hash != EXPECTED_LEDGER_SHA256:
        raise ValueError(
            f"ledger SHA-256 mismatch: {actual_hash} != {EXPECTED_LEDGER_SHA256}"
        )
    raw_lines = payload.splitlines()
    if len(raw_lines) != EXPECTED_LEDGER_RECORDS:
        raise ValueError(
            f"ledger row count mismatch: {len(raw_lines)} != {EXPECTED_LEDGER_RECORDS}"
        )
    rows: list[dict[str, Any]] = []
    seen: set[str] = set()
    for line_number, raw_line in enumerate(raw_lines, start=1):
        row = json.loads(raw_line)
        item_id = str(row.get("item_id") or "")
        if not item_id or item_id in seen:
            raise ValueError(f"invalid or duplicate item_id at ledger line {line_number}")
        seen.add(item_id)
        rows.append(
            {
                "line": line_number,
                "row_sha256": sha256_bytes(raw_line),
                "record": row,
            }
        )
    return rows


def _observation(kind: str, value: float, formula: str) -> dict[str, Any]:
    calculated = formula_mz(formula)
    key = "product_mz" if kind == "product_ion" else "neutral_loss_da"
    formula_key = "product_formula" if kind == "product_ion" else "neutral_loss_formula"
    return {
        "type": kind,
        key: value,
        formula_key: formula,
        "formula_calculated_da": round(calculated, 7),
        "formula_mass_error_da": round(value - calculated, 7),
    }


def _proposal(
    proposal_id: str,
    *,
    observation_type: str,
    value: float,
    formula: str,
    classes: list[str],
    structural_scope: str,
    polarity: str,
    adducts: list[str],
    charges: list[int],
    group: str,
    legacy_item_ids: list[str],
    statistics_plan_id: str,
    source_claim_ids: list[str],
    status: str = "needs_recount",
    existing_rule_id: str | None = None,
    domain_id: str | None = None,
    activation: str | None = None,
    collision_energy: dict[str, Any] | None = None,
    statistics: dict[str, Any] | None = None,
    limitations: list[str] | None = None,
    missing_schema_fields: list[str] | None = None,
) -> dict[str, Any]:
    if status not in PROPOSAL_STATUSES:
        raise ValueError(f"unsupported proposal status: {status}")
    return {
        "schema_version": CURATION_SCHEMA_VERSION,
        "proposal_id": proposal_id,
        "status": status,
        "runtime_influence": False,
        "identity_rank_influence": False,
        "export_eligible": False,
        "research_audit_only": True,
        "benchmark_role": (
            "deterministic_ci_fixture_only"
            if status == "accepted_shadow"
            else "curation_candidate"
        ),
        "existing_runtime_rule_id": existing_rule_id,
        "observation": _observation(observation_type, value, formula),
        "chemical_scope": {
            "classes": classes,
            "structural_scope": structural_scope,
            "claim_scope": "chemical_class_support_only",
        },
        "ionization": {
            "polarity": polarity,
            "adducts": adducts,
            "charges": charges,
        },
        "acquisition_domain": {
            "instrument_domain_ids": [domain_id] if domain_id else [],
            "activation": activation,
            "activation_evidence": (
                "MassBank contributor text" if domain_id else None
            ),
            "psi_beam_type_cid_confirmed": False if domain_id else None,
            "collision_energy": collision_energy,
            "cross_domain_transfer_allowed": False,
        },
        "mass_tolerance": {
            "absolute_da": 0.01 if observation_type == "product_ion" else 0.015,
            "ppm": 20.0,
            "policy": "both recorded; implementation may use the larger physical window",
        },
        "statistics_plan_id": statistics_plan_id,
        "existing_statistics": statistics,
        "independent_fragmentation_group": group,
        "positive_evidence_limitations": limitations
        or [
            "Class-support evidence only",
            "Does not confirm exact identity, formula, standard identity, or isomer",
        ],
        "negative_evidence": {
            "absence_interpretable": status == "accepted_shadow",
            "required_quality_fields": [
                "complete_mass_range",
                "centroided",
                "sufficient_signal",
                "precursor_isolated",
                "collision_energy_verified",
                "instrument_domain_verified",
            ],
            "limitations": [
                "Absence is not usable outside a validated exact acquisition stratum"
            ],
        },
        "source_claim_ids": source_claim_ids,
        "evidence_level": (
            "public_experimental_ci_fixture_single_laboratory_not_empirically_validated"
            if status == "accepted_shadow"
            else "curation_candidate_not_validated"
        ),
        "empirical_validation_gate_passed": False,
        "empirical_validation_gate_unmet": (
            [
                "independent_sources_lt_3",
                "independent_laboratories_lt_2",
                "sealed_identity_formula_holdout_not_run_for_v2_export",
            ]
            if status == "accepted_shadow"
            else ["strict_recount_not_completed"]
        ),
        "legacy_item_ids": legacy_item_ids,
        "missing_schema_fields": missing_schema_fields
        if missing_schema_fields is not None
        else (
            []
            if status == "accepted_shadow"
            else [
                "exact_instrument_activation_ce_stratum",
                "target_occurrence_statistics",
                "background_occurrence_statistics",
                "identity_and_formula_isolated_validation",
                "negative_evidence_calibration",
            ]
        ),
        "output_semantics": {
            "evidence_stage": "msms_rule_support",
            "exact_identity_confirmed": False,
            "standard_confirmed": False,
            "isomer_discriminating": False,
        },
    }


def build_proposals() -> list[dict[str, Any]]:
    proposals: list[dict[str, Any]] = []

    pc_legacy = [
        "static_fragmentation_rule:PHOSPHOLIPID_HEADGROUP_CLEAVAGE",
        "static_diagnostic_ion:phosphocholine_184",
        "feature_profile:lipid:001",
        "legacy_rule:LIP-001#01",
        "legacy_rule:PL-001#01",
        "legacy_rule:SPH-002#01",
    ]
    pe_legacy = [
        "static_neutral_loss:C2H8NO4P",
        "static_fragmentation_rule:PHOSPHATIDYLETHANOLAMINE_141_LOSS",
        "legacy_rule:PL-002#01",
    ]
    for energy, pc_present, pe_present in ((10, 6, 3), (20, 6, 3), (40, 6, 2)):
        ce = {
            "raw": f"{energy} eV",
            "unit": "eV",
            "mode": "fixed",
            "minimum": float(energy),
            "maximum": float(energy),
        }
        proposals.append(
            _proposal(
                f"CUR-PC184-AG6560-{energy}EV",
                observation_type="product_ion",
                value=184.0733,
                formula="C5H15NO4P+",
                classes=["phosphocholine-containing lipids"],
                structural_scope="PC, LPC, SM, and other phosphocholine headgroups",
                polarity="positive",
                adducts=["[M+H]+"],
                charges=[1],
                group="phosphocholine_headgroup_184",
                legacy_item_ids=pc_legacy,
                statistics_plan_id="STAT-PC184-AG6560-PUBLIC-PILOT",
                source_claim_ids=["SRC-MASSBANK-AG6560-PC-PE-PILOT"],
                status="accepted_shadow",
                existing_rule_id=f"FF-PC-184-POS-H-AG6560-{energy}EV",
                domain_id="agilent_6560_qtof_cid_contributor_text_fixed_ev",
                activation="CID",
                collision_energy=ce,
                statistics={
                    "target": {
                        "present": pc_present,
                        "total": 6,
                        "frequency": pc_present / 6,
                        "counting_unit": "independent_2d_identity",
                    },
                    "background": {
                        "present": 0,
                        "total": 4,
                        "frequency": 0.0,
                        "counting_unit": "independent_2d_identity",
                    },
                },
                limitations=[
                    "Supports a phosphocholine headgroup only in the exact Agilent 6560 domain",
                    "Cannot distinguish PC/LPC from sphingomyelin",
                    "Cannot assign acyl positions, double bonds, stereochemistry, or exact identity",
                ],
            )
        )
        proposals.append(
            _proposal(
                f"CUR-PE-NL141-AG6560-{energy}EV",
                observation_type="neutral_loss",
                value=141.0191,
                formula="C2H8NO4P",
                classes=["glycerophosphoethanolamines"],
                structural_scope="protonated PE/LPE within the measured domain",
                polarity="positive",
                adducts=["[M+H]+"],
                charges=[1],
                group="phosphoethanolamine_headgroup_loss_141",
                legacy_item_ids=pe_legacy,
                statistics_plan_id="STAT-PE141-AG6560-PUBLIC-PILOT",
                source_claim_ids=["SRC-MASSBANK-AG6560-PC-PE-PILOT"],
                status="accepted_shadow",
                existing_rule_id=f"FF-PE-NL141-POS-H-AG6560-{energy}EV",
                domain_id="agilent_6560_qtof_cid_contributor_text_fixed_ev",
                activation="CID",
                collision_energy=ce,
                statistics={
                    "target": {
                        "present": pe_present,
                        "total": 3,
                        "frequency": pe_present / 3,
                        "counting_unit": "independent_2d_identity",
                    },
                    "background": {
                        "present": 0,
                        "total": 7,
                        "frequency": 0.0,
                        "counting_unit": "independent_2d_identity",
                    },
                },
                limitations=[
                    "Supports a PE/LPE-class hypothesis only in the exact Agilent 6560 domain",
                    "Does not assign acyl positions, ether/plasmalogen state, double bonds, or exact identity",
                    "The 40 eV public pilot contains one target absence",
                ],
            )
        )

    nuc_common = ["legacy_rule:NUC-001#02"]
    proposals.extend(
        [
            _proposal(
                "CUR-NUC-RIBOSE-NL132-POS-H",
                observation_type="neutral_loss",
                value=132.0423,
                formula="C5H8O4",
                classes=["ribonucleosides"],
                structural_scope="N-glycosidic cleavage of protonated ribonucleosides",
                polarity="positive",
                adducts=["[M+H]+"],
                charges=[1],
                group="nucleoside_sugar_cleavage",
                legacy_item_ids=[
                    "static_neutral_loss:C5H8O4",
                    "static_fragmentation_rule:RIBONUCLEOSIDE_GLYCOSIDIC_CLEAVAGE",
                    "legacy_rule:NUC-001#02",
                    "legacy_rule:NUC-004#01",
                ],
                statistics_plan_id="STAT-NUCLEOSIDE-POS-H",
                source_claim_ids=["SRC-NUCLEOSIDE-ORBITRAP-DIA-2023"],
                limitations=[
                    "A 132 Da loss is shared with pentose O-glycosides",
                    "Does not distinguish nucleosides from nucleotides or modified nucleosides",
                ],
            ),
            _proposal(
                "CUR-NUC-DEOXYRIBOSE-NL116-POS-H",
                observation_type="neutral_loss",
                value=116.0473,
                formula="C5H8O3",
                classes=["deoxyribonucleosides"],
                structural_scope="N-glycosidic cleavage of protonated deoxyribonucleosides",
                polarity="positive",
                adducts=["[M+H]+"],
                charges=[1],
                group="nucleoside_sugar_cleavage",
                legacy_item_ids=[
                    "static_neutral_loss:C5H8O3",
                    "static_fragmentation_rule:DEOXYRIBONUCLEOSIDE_GLYCOSIDIC_CLEAVAGE",
                    "legacy_rule:NUC-001#02",
                ],
                statistics_plan_id="STAT-NUCLEOSIDE-POS-H",
                source_claim_ids=["SRC-NUCLEOSIDE-ORBITRAP-DIA-2023"],
            ),
        ]
    )
    nuc_ions = [
        ("CYTOSINE", 112.0505, "C4H6N3O+", "cytosine", "001", "nucleoside_cytosine_112", "contrastive_v2:nucleoside:pos:diagnostic_ions:024"),
        ("URACIL", 113.0346, "C4H5N2O2+", "uracil", "002", "nucleoside_uracil_113", None),
        ("THYMINE", 127.0502, "C5H7N2O2+", "thymine", "003", "nucleoside_thymine_127", None),
        ("ADENINE", 136.0618, "C5H6N5+", "adenine", "004", "nucleoside_adenine_136", "contrastive_v2:nucleoside:pos:diagnostic_ions:004"),
        ("GUANINE", 152.0567, "C5H6N5O+", "guanine", "005", "nucleoside_guanine_152", None),
    ]
    for label, mz, formula, base, ordinal, diagnostic_id, pooled_seed_id in nuc_ions:
        legacy_ids = [
            f"static_diagnostic_ion:{diagnostic_id}",
            f"feature_profile:nucleoside:{ordinal}",
            *nuc_common,
            "legacy_rule:NUC-004#01",
        ]
        if pooled_seed_id:
            legacy_ids.append(pooled_seed_id)
        proposals.append(
            _proposal(
                f"CUR-NUC-{label}-{mz:.4f}-POS-H",
                observation_type="product_ion",
                value=mz,
                formula=formula,
                classes=[f"{base}-containing nucleosides/nucleotides"],
                structural_scope=f"protonated {base} product ion",
                polarity="positive",
                adducts=["[M+H]+"],
                charges=[1],
                group=f"nucleobase_product_{base}",
                legacy_item_ids=legacy_ids,
                statistics_plan_id="STAT-NUCLEOSIDE-POS-H",
                source_claim_ids=[
                    "SRC-NUCLEOSIDE-ORBITRAP-DIA-2023",
                    "SRC-NUCLEOSIDE-CID-HCD-DOMAIN-2024",
                ],
                limitations=[
                    f"Supports a {base}-containing precursor hypothesis",
                    "Does not distinguish nucleoside, nucleotide, modified base, or positional isomer",
                ],
            )
        )

    ogly_losses = [
        ("HEXOSE", 162.0528, "C6H10O5", "static_neutral_loss:C6H10O5"),
        ("DEOXYHEXOSE", 146.0579, "C6H10O4", "static_neutral_loss:C6H10O4"),
        ("PENTOSE", 132.0423, "C5H8O4", "static_neutral_loss:C5H8O4"),
        ("URONIC", 176.0321, "C6H8O6", "static_neutral_loss:C6H8O6"),
    ]
    for label, loss, formula, static_id in ogly_losses:
        for polarity, suffix, adduct in (
            ("positive", "POS-H", "[M+H]+"),
            ("negative", "NEG-H", "[M-H]-"),
        ):
            legacy_ids = [
                static_id,
                "static_fragmentation_rule:GLYCOSIDIC_CLEAVAGE",
                "legacy_rule:GLY-001#01",
            ]
            if label != "URONIC":
                legacy_ids.append("legacy_rule:FLAV-020#01")
            proposals.append(
                _proposal(
                    f"CUR-OGLY-{label}-NL{loss:.4f}-{suffix}",
                    observation_type="neutral_loss",
                    value=loss,
                    formula=formula,
                    classes=["O-glycosides"],
                    structural_scope=f"O-linked {label.lower()} residue cleavage",
                    polarity=polarity,
                    adducts=[adduct],
                    charges=[1 if polarity == "positive" else -1],
                    group=f"oglycoside_{label.lower()}_cleavage",
                    legacy_item_ids=legacy_ids,
                    statistics_plan_id="STAT-OGLYCOSIDE-BY-DOMAIN",
                    source_claim_ids=["SRC-OGLYCOSIDE-LEGACY-REFERENCES-PENDING"],
                    limitations=[
                        "The same nominal loss occurs outside O-glycosides",
                        "O- versus C-glycoside discrimination requires explicit backgrounds",
                    ],
                )
            )

    for loss, formula, static_id in (
        (60.0211, "C2H4O2", "static_neutral_loss:C2H4O2_CGLY"),
        (90.0317, "C3H6O3", "static_neutral_loss:C3H6O3"),
        (120.0423, "C4H8O4", "static_neutral_loss:C4H8O4"),
    ):
        proposals.append(
            _proposal(
                f"CUR-CGLY-XRING-NL{loss:.4f}-NEG-H",
                observation_type="neutral_loss",
                value=loss,
                formula=formula,
                classes=["phenolic and flavonoid C-glycosides"],
                structural_scope="C-linked sugar cross-ring fragmentation",
                polarity="negative",
                adducts=["[M-H]-"],
                charges=[-1],
                group="cglycoside_cross_ring_60_90_120",
                legacy_item_ids=[
                    static_id,
                    "static_fragmentation_rule:C_GLYCOSIDE_CROSS_RING",
                    "legacy_rule:FLAV-021#01",
                ],
                statistics_plan_id="STAT-CGLYCOSIDE-NEG-H",
                source_claim_ids=["SRC-CGLYCOSIDE-LEGACY-REFERENCES-PENDING"],
                limitations=[
                    "The 60/90/120 series is correlated and must contribute at most one independent event",
                    "Does not identify the aglycone or sugar stereochemistry",
                ],
            )
        )

    proposals.extend(
        [
            _proposal(
                "CUR-AA-ALPHA-NL46.0055-POS-H",
                observation_type="neutral_loss",
                value=46.0055,
                formula="CH2O2",
                classes=["alpha-amino acids"],
                structural_scope="combined H2O and CO loss from protonated alpha-amino acids",
                polarity="positive",
                adducts=["[M+H]+"],
                charges=[1],
                group="neutral_loss_46_0055_observed_mass",
                legacy_item_ids=[
                    "static_neutral_loss:H2O_CO",
                    "static_fragmentation_rule:ALPHA_AMINO_ACID_46_LOSS",
                    "empirical_v1:amino_acid:pos:neutral_losses:001",
                ],
                statistics_plan_id="STAT-AA-ALPHA-NL46-POS-H",
                source_claim_ids=["SRC-ALPHA-AA-NL46-2024"],
                limitations=[
                    "The observed CH2O2 mass difference cannot distinguish H2O+CO from formic-acid loss",
                    "Beta, gamma, and epsilon amino acids are required backgrounds",
                ],
            ),
            _proposal(
                "CUR-SULFATE-NL79.9568-NEG-H",
                observation_type="neutral_loss",
                value=79.9568,
                formula="SO3",
                classes=["sulfate esters"],
                structural_scope="deprotonated sulfate ester SO3 loss",
                polarity="negative",
                adducts=["[M-H]-"],
                charges=[-1],
                group="sulfate_headgroup_so3_80_97",
                legacy_item_ids=[
                    "static_neutral_loss:SO3",
                    "static_fragmentation_rule:SULFATE_SO3_LOSS",
                    "legacy_rule:FLAV-027#02",
                    "legacy_rule:SUL-001#01",
                ],
                statistics_plan_id="STAT-SULFATE-SO3-NEG-H",
                source_claim_ids=["SRC-SULFATE-SO3-SOURCE-PENDING"],
                limitations=[
                    "First curation proposal is negative [M-H]- only",
                    "Sulfates, sulfonates, phosphates, sulfites, and nonsulfated analogs require separate backgrounds",
                ],
            ),
        ]
    )

    lipid_specs = [
        (
            "CUR-SPH-D18_1-264.2686-POS-H",
            "product_ion",
            264.2686,
            "C18H34N+",
            ["d18:1-containing sphingolipids"],
            "positive",
            ["[M+H]+"],
            [1],
            "sphingoid_d18_1_product_264",
            [
                "static_fragmentation_rule:SPHINGOID_BASE_CLEAVAGE",
                "static_diagnostic_ion:sphingoid_d18_1_264",
                "feature_profile:lipid:004",
                "legacy_rule:SPH-002#01",
            ],
            ["SRC-SPHINGOLIPID-264-LEGACY"],
        ),
        (
            "CUR-PC-224.0693-NEG-FA",
            "product_ion",
            224.0693,
            "C7H15NO5P-",
            ["glycerophosphocholines"],
            "negative",
            ["[M+FA-H]-"],
            [-1],
            "glycerophosphocholine_headgroup_224",
            [
                "static_diagnostic_ion:glycerophosphocholine_224",
                "feature_profile:lipid:005",
            ],
            ["SRC-LIPID-HEADGROUP-AIF-2020"],
        ),
        (
            "CUR-PC-224.0693-NEG-HAC",
            "product_ion",
            224.0693,
            "C7H15NO5P-",
            ["glycerophosphocholines"],
            "negative",
            ["[M+Hac-H]-"],
            [-1],
            "glycerophosphocholine_headgroup_224",
            [
                "static_diagnostic_ion:glycerophosphocholine_224",
                "feature_profile:lipid:005",
            ],
            ["SRC-LIPID-HEADGROUP-AIF-2020"],
        ),
        (
            "CUR-PE-196.0380-NEG-H",
            "product_ion",
            196.0380,
            "C5H11NO5P-",
            ["glycerophosphoethanolamines"],
            "negative",
            ["[M-H]-"],
            [-1],
            "phosphoethanolamine_product_196",
            [
                "static_diagnostic_ion:phosphoethanolamine_196",
                "feature_profile:lipid:006",
                "legacy_rule:PL-002#01",
            ],
            ["SRC-PE-196-RCM-2023", "SRC-LIPID-HEADGROUP-AIF-2020"],
        ),
        (
            "CUR-PI-241.0117-NEG-H",
            "product_ion",
            241.0117,
            "C6H10O8P-",
            ["glycerophosphoinositols"],
            "negative",
            ["[M-H]-"],
            [-1],
            "inositol_phosphate_product_241",
            [
                "static_diagnostic_ion:inositol_phosphate_241",
                "feature_profile:lipid:007",
            ],
            ["SRC-LIPID-HEADGROUP-AIF-2020"],
        ),
    ]
    for spec in lipid_specs:
        proposal_id, kind, value, formula, classes, polarity, adducts, charges, group, legacy_ids, claims = spec
        proposals.append(
            _proposal(
                proposal_id,
                observation_type=kind,
                value=value,
                formula=formula,
                classes=classes,
                structural_scope="headgroup or sphingoid-chain support; no sn-position, double-bond, or isomer assignment",
                polarity=polarity,
                adducts=adducts,
                charges=charges,
                group=group,
                legacy_item_ids=legacy_ids,
                statistics_plan_id="STAT-LIPID-HEADGROUP-BY-DOMAIN",
                source_claim_ids=claims,
            )
        )

    immonium_specs = [
        ("LEU_ILE", 86.0964, "C5H12N+", "001"),
        ("PHE", 120.0813, "C8H10N+", "002"),
        ("TYR", 136.0762, "C8H10NO+", "003"),
        ("HIS", 110.0718, "C5H8N3+", "004"),
        ("TRP", 159.0922, "C10H11N2+", "005"),
        ("ARG", 129.1140, "C5H13N4+", "006"),
        ("LYS", 101.1079, "C5H13N2+", "007"),
        ("VAL", 72.0813, "C4H10N+", "008"),
        ("PRO", 70.0656, "C4H8N+", "009"),
    ]
    for label, mz, formula, ordinal in immonium_specs:
        legacy_ids = [f"feature_profile:amino_acid:{ordinal}"]
        if label == "PRO":
            legacy_ids.append("contrastive_v2:amino_acid:pos:diagnostic_ions:002")
        proposals.append(
            _proposal(
                f"CUR-PEPTIDE-IMMONIUM-{label}-{mz:.4f}-POS",
                observation_type="product_ion",
                value=mz,
                formula=formula,
                classes=["peptides containing the named residue"],
                structural_scope="peptide-residue immonium ion; not yet validated for free amino acids",
                polarity="positive",
                adducts=[],
                charges=[1],
                group=(
                    "shared_product_86_0964_leu_ile_lipid"
                    if label == "LEU_ILE"
                    else f"peptide_immonium_{label.lower()}"
                ),
                legacy_item_ids=legacy_ids,
                statistics_plan_id="STAT-PEPTIDE-IMMONIUM-QTOF-CID",
                source_claim_ids=["SRC-PEPTIDE-IMMONIUM-QSTAR-2008"],
                status="literature_only",
                limitations=[
                    "The cited study concerns peptide residue composition, not free amino-acid standards",
                    "Leu/Ile and other isomeric contexts cannot be resolved by this ion alone",
                ],
                missing_schema_fields=[
                    "exact_precursor_adduct_and_charge_scope",
                    "free_amino_acid_applicability_if_claimed",
                    "exact_ce_strata",
                    "target_and_background_identity_counts",
                    "identity_and_formula_isolated_validation",
                ],
            )
        )

    proposal_ids = [item["proposal_id"] for item in proposals]
    if len(proposal_ids) != len(set(proposal_ids)):
        raise ValueError("duplicate proposal_id")
    for proposal in proposals:
        if abs(float(proposal["observation"]["formula_mass_error_da"])) > 0.01:
            raise ValueError(
                f"proposal formula mismatch: {proposal['proposal_id']} "
                f"{proposal['observation']['formula_mass_error_da']}"
            )
    return sorted(proposals, key=lambda item: item["proposal_id"])


def _proposal_ids(proposals: list[dict[str, Any]], prefix: str) -> list[str]:
    return [item["proposal_id"] for item in proposals if item["proposal_id"].startswith(prefix)]


def build_aliases(
    ledger_rows: list[dict[str, Any]], proposals: list[dict[str, Any]]
) -> list[dict[str, Any]]:
    proposal_ids = {item["proposal_id"] for item in proposals}
    overrides: dict[str, dict[str, Any]] = {}

    def mark(
        item_ids: str | Iterable[str],
        relationship: str,
        *,
        canonical: Iterable[str] = (),
        reason: str,
        **extra: Any,
    ) -> None:
        ids = [item_ids] if isinstance(item_ids, str) else list(item_ids)
        canonical_ids = list(canonical)
        unknown = set(canonical_ids) - proposal_ids
        if unknown:
            raise ValueError(f"unknown canonical proposals: {sorted(unknown)}")
        for item_id in ids:
            if item_id in overrides:
                raise ValueError(f"duplicate curation override: {item_id}")
            overrides[item_id] = {
                "relationship": relationship,
                "canonical_proposal_ids": canonical_ids,
                "reason": reason,
                **extra,
            }

    pc = _proposal_ids(proposals, "CUR-PC184-")
    pe141 = _proposal_ids(proposals, "CUR-PE-NL141-")
    nuc = _proposal_ids(proposals, "CUR-NUC-")
    ogly = _proposal_ids(proposals, "CUR-OGLY-")
    cgly = _proposal_ids(proposals, "CUR-CGLY-")
    immonium = _proposal_ids(proposals, "CUR-PEPTIDE-IMMONIUM-")

    mark(
        "static_fragmentation_rule:PHOSPHOLIPID_HEADGROUP_CLEAVAGE",
        "split_parent",
        canonical=pc,
        reason="Broad phospholipid parent split into exact CE-domain PC184 observations.",
    )
    mark(
        "static_diagnostic_ion:phosphocholine_184",
        "split_parent",
        canonical=pc,
        reason="One domain-free observation maps to three validated CE strata.",
    )
    mark(
        "feature_profile:lipid:001",
        "numeric_correction",
        canonical=pc,
        reason="184.0739 is corrected to the charged-formula value 184.0733 before CE stratification.",
        alias_of_item_id="static_diagnostic_ion:phosphocholine_184",
        field_overrides={
            "mz": {
                "from": 184.0739,
                "to": 184.0733,
                "basis": "C5H15NO4P+ monoisotopic m/z",
            }
        },
    )
    mark(
        ["legacy_rule:LIP-001#01", "legacy_rule:PL-001#01"],
        "split_parent",
        canonical=pc,
        reason="Legacy domain-free PC184 entry maps to three Agilent 6560 CE strata.",
    )
    mark(
        "legacy_rule:SPH-002#01",
        "split_parent",
        canonical=[*pc, "CUR-SPH-D18_1-264.2686-POS-H"],
        reason="Bundled sphingomyelin record split into phosphocholine 184 and d18:1 264 observations.",
    )

    mark(
        [
            "static_neutral_loss:C2H8NO4P",
            "static_fragmentation_rule:PHOSPHATIDYLETHANOLAMINE_141_LOSS",
        ],
        "split_parent",
        canonical=pe141,
        reason="Domain-free PE loss maps to three validated Agilent 6560 CE strata.",
    )
    mark(
        "legacy_rule:PL-002#01",
        "split_parent",
        canonical=[*pe141, "CUR-PE-196.0380-NEG-H"],
        reason="Positive PE loss 141 and negative product 196 are distinct atomic observations.",
    )
    mark(
        "legacy_rule:LIP-002#01",
        "rejected_incorrect",
        reason="Record incorrectly assigns the PE 141.019 Da loss to PC and mixes 183/141 events.",
    )

    mark(
        "static_neutral_loss:C5H8O4",
        "split_parent",
        canonical=["CUR-NUC-RIBOSE-NL132-POS-H", *_proposal_ids(proposals, "CUR-OGLY-PENTOSE-")],
        reason="The same 132.0423 loss occurs in ribonucleosides and pentose O-glycosides.",
        duplicate_observation_group="neutral_loss_132_0423_shared",
    )
    mark(
        "static_fragmentation_rule:RIBONUCLEOSIDE_GLYCOSIDIC_CLEAVAGE",
        "exact_alias",
        canonical=["CUR-NUC-RIBOSE-NL132-POS-H"],
        reason="Chemistry and ionization match; experimental domain statistics remain missing.",
    )
    mark(
        "static_neutral_loss:C5H8O3",
        "exact_alias",
        canonical=["CUR-NUC-DEOXYRIBOSE-NL116-POS-H"],
        reason="Formula, loss, class guard, polarity, and adduct agree.",
    )
    mark(
        "static_fragmentation_rule:DEOXYRIBONUCLEOSIDE_GLYCOSIDIC_CLEAVAGE",
        "exact_alias",
        canonical=["CUR-NUC-DEOXYRIBOSE-NL116-POS-H"],
        reason="Chemistry and ionization match; experimental domain statistics remain missing.",
    )
    nuc_pairs = [
        ("nucleoside_cytosine_112", "001", "CUR-NUC-CYTOSINE-112.0505-POS-H", "contrastive_v2:nucleoside:pos:diagnostic_ions:024"),
        ("nucleoside_uracil_113", "002", "CUR-NUC-URACIL-113.0346-POS-H", None),
        ("nucleoside_thymine_127", "003", "CUR-NUC-THYMINE-127.0502-POS-H", None),
        ("nucleoside_adenine_136", "004", "CUR-NUC-ADENINE-136.0618-POS-H", "contrastive_v2:nucleoside:pos:diagnostic_ions:004"),
        ("nucleoside_guanine_152", "005", "CUR-NUC-GUANINE-152.0567-POS-H", None),
    ]
    for diagnostic_id, ordinal, proposal_id, pooled_seed_id in nuc_pairs:
        static_id = f"static_diagnostic_ion:{diagnostic_id}"
        mark(
            static_id,
            "exact_alias",
            canonical=[proposal_id],
            reason="Formula, product m/z, positive polarity, and [M+H]+ scope agree.",
        )
        mark(
            f"feature_profile:nucleoside:{ordinal}",
            "exact_alias",
            canonical=[proposal_id],
            reason="Duplicate class-profile representation of the static nucleobase ion.",
            alias_of_item_id=static_id,
        )
        if pooled_seed_id:
            mark(
                pooled_seed_id,
                "split_parent",
                canonical=[proposal_id],
                reason="Pooled multi-instrument discovery statistics identify a recount seed but cannot migrate as strict validation.",
            )
    mark(
        "legacy_rule:NUC-001#02",
        "split_parent",
        canonical=nuc,
        reason="Bundled sugar losses and base ions are represented as separate atomic proposals.",
    )
    mark(
        "legacy_rule:NUC-004#01",
        "split_parent",
        canonical=[item for item in nuc if "DEOXYRIBOSE" not in item],
        reason="Bundled ribose loss and base ions are represented separately; this occurrence does not support deoxyribose loss 116.",
    )
    mark(
        "legacy_rule:NUC-001#01",
        "superseded",
        reason="Replaced by the fuller NUC-001#02 occurrence.",
        superseded_by_item_id="legacy_rule:NUC-001#02",
    )

    for static_id, prefix in (
        ("static_neutral_loss:C6H10O5", "CUR-OGLY-HEXOSE-"),
        ("static_neutral_loss:C6H10O4", "CUR-OGLY-DEOXYHEXOSE-"),
        ("static_neutral_loss:C6H8O6", "CUR-OGLY-URONIC-"),
    ):
        mark(
            static_id,
            "split_parent",
            canonical=_proposal_ids(proposals, prefix),
            reason="Both-polarity legacy loss split into positive and negative adduct-specific proposals.",
        )
    mark(
        "static_fragmentation_rule:GLYCOSIDIC_CLEAVAGE",
        "split_parent",
        canonical=ogly,
        reason="Broad glycosidic rule split by residue mass, polarity, and adduct.",
    )
    mark(
        "legacy_rule:FLAV-020#01",
        "split_parent",
        canonical=[item for item in ogly if "URONIC" not in item],
        reason="Legacy flavonoid rule bundles hexose, deoxyhexose, and pentose losses.",
    )
    mark(
        "legacy_rule:GLY-001#01",
        "split_parent",
        canonical=ogly,
        reason="General O-glycoside rule split into eight residue/polarity proposals.",
    )
    mark(
        "legacy_rule:GLY-001#02",
        "defer",
        reason="Glycoglycerolipid sodium/lithium sequential losses are not aliases of GLY-001#01 and lack an atomic source.",
    )

    for static_id, proposal_id in (
        ("static_neutral_loss:C2H4O2_CGLY", "CUR-CGLY-XRING-NL60.0211-NEG-H"),
        ("static_neutral_loss:C3H6O3", "CUR-CGLY-XRING-NL90.0317-NEG-H"),
        ("static_neutral_loss:C4H8O4", "CUR-CGLY-XRING-NL120.0423-NEG-H"),
    ):
        mark(
            static_id,
            "exact_alias",
            canonical=[proposal_id],
            reason="Formula, mass, negative polarity, adduct, and C-glycoside structure guard agree.",
            duplicate_observation_group=(
                "neutral_loss_60_0211_observed_mass" if "C2H4O2" in static_id else None
            ),
        )
    mark(
        ["static_fragmentation_rule:C_GLYCOSIDE_CROSS_RING", "legacy_rule:FLAV-021#01"],
        "split_parent",
        canonical=cgly,
        reason="Correlated 60/90/120 cross-ring series split into atomic rules sharing one evidence group.",
    )

    mark(
        "static_neutral_loss:H2O_CO",
        "exact_alias",
        canonical=["CUR-AA-ALPHA-NL46.0055-POS-H"],
        reason="Formula, mass, class guard, positive polarity, and [M+H]+ agree.",
        duplicate_observation_group="neutral_loss_46_0055_observed_mass",
    )
    mark(
        "static_fragmentation_rule:ALPHA_AMINO_ACID_46_LOSS",
        "exact_alias",
        canonical=["CUR-AA-ALPHA-NL46.0055-POS-H"],
        reason="Primary literature supports the alpha-amino-acid 46 Da channel, pending domain counts.",
    )
    mark(
        "empirical_v1:amino_acid:pos:neutral_losses:001",
        "numeric_correction",
        canonical=["CUR-AA-ALPHA-NL46.0055-POS-H"],
        reason="46.0125 is a pooled bin center, not the CH2O2 exact mass; pooled spectra cannot supply strict counts.",
        field_overrides={
            "neutral_loss_da": {
                "from": 46.0125,
                "to": 46.0055,
                "basis": "CH2O2 monoisotopic neutral mass",
            }
        },
    )
    mark(
        "static_neutral_loss:HCOOH",
        "defer",
        reason="Observed mass is identical to H2O+CO; mechanism cannot be inferred from mass alone.",
        duplicate_observation_group="neutral_loss_46_0055_observed_mass",
    )
    mark(
        "static_neutral_loss:CH3COOH",
        "defer",
        reason="Observed mass is identical to the C-glycoside C2H4O2 proposal and lacks its structural context.",
        duplicate_observation_group="neutral_loss_60_0211_observed_mass",
    )

    mark(
        "static_neutral_loss:SO3",
        "split_parent",
        canonical=["CUR-SULFATE-NL79.9568-NEG-H"],
        reason="Both-polarity common loss is narrowed to the supported negative sulfate-ester branch; no positive rule is accepted.",
    )
    mark(
        [
            "static_fragmentation_rule:SULFATE_SO3_LOSS",
            "legacy_rule:FLAV-027#02",
            "legacy_rule:SUL-001#01",
        ],
        "exact_alias",
        canonical=["CUR-SULFATE-NL79.9568-NEG-H"],
        reason="Negative sulfate-ester SO3 observation retained; both-polarity legacy generalization is narrowed.",
    )

    mark(
        "static_fragmentation_rule:SPHINGOID_BASE_CLEAVAGE",
        "exact_alias",
        canonical=["CUR-SPH-D18_1-264.2686-POS-H"],
        reason="d18:1 product proposal retains the narrow structural scope.",
    )
    lipid_pairs = [
        ("sphingoid_d18_1_264", "004", ["CUR-SPH-D18_1-264.2686-POS-H"]),
        ("glycerophosphocholine_224", "005", ["CUR-PC-224.0693-NEG-FA", "CUR-PC-224.0693-NEG-HAC"]),
        ("phosphoethanolamine_196", "006", ["CUR-PE-196.0380-NEG-H"]),
        ("inositol_phosphate_241", "007", ["CUR-PI-241.0117-NEG-H"]),
    ]
    for diagnostic_id, ordinal, canonical in lipid_pairs:
        static_id = f"static_diagnostic_ion:{diagnostic_id}"
        relationship = "split_parent" if len(canonical) > 1 else "exact_alias"
        mark(
            static_id,
            relationship,
            canonical=canonical,
            reason="Static diagnostic ion mapped to explicit adduct proposal(s).",
        )
        mark(
            f"feature_profile:lipid:{ordinal}",
            relationship,
            canonical=canonical,
            reason="Duplicate profile representation; negative entries ignore the inherited positive adduct preference.",
            alias_of_item_id=static_id,
            field_overrides=(
                {
                    "adduct_preference": {
                        "from": "[M+H]+",
                        "to": "use allowed_adducts only",
                        "basis": "profile polarity/adduct inconsistency",
                    }
                }
                if ordinal in {"005", "006", "007"}
                else {}
            ),
        )

    for proposal in [item for item in proposals if item["proposal_id"] in immonium]:
        for item_id in proposal["legacy_item_ids"]:
            if item_id in overrides:
                continue
            mark(
                item_id,
                "split_parent",
                canonical=[proposal["proposal_id"]],
                reason="Existing free-amino-acid profile is respecified as a peptide-residue literature candidate; no scoring allowed.",
            )

    mark(
        "feature_profile:flavonoid:006",
        "quarantined_duplicate",
        reason="Exact duplicate of feature_profile:flavonoid:002 at m/z 153.0188.",
        alias_of_item_id="feature_profile:flavonoid:002",
    )
    mark(
        "feature_profile:glycoside:017",
        "quarantined_duplicate",
        reason="Exact duplicate of feature_profile:glycoside:004 at m/z 115.0395.",
        alias_of_item_id="feature_profile:glycoside:004",
    )
    mark(
        "legacy_rule:FLAV-026#03",
        "quarantined_duplicate",
        reason="Legacy occurrence is an exact duplicate and cannot contribute a second rule.",
        alias_of_item_id="legacy_rule:FLAV-026#02",
    )
    mark(
        "legacy_rule:FLAV-027#03",
        "quarantined_duplicate",
        reason="Legacy occurrence is an exact duplicate and cannot contribute a second rule.",
        alias_of_item_id="legacy_rule:FLAV-027#02",
    )
    mark(
        "legacy_rule:ALK-009b#02",
        "quarantined_misplaced",
        reason="Aconitine-like alkaloid rule is misplaced in the terpenoid section and remains overly broad.",
    )
    for item_id in FORMULA_MASS_QUARANTINE_IDS:
        mark(
            item_id,
            "quarantined_formula_mass",
            reason="Stored formula and m/z are inconsistent; no correction is guessed without the original spectrum/source.",
        )

    aliases: list[dict[str, Any]] = []
    ledger_ids = {entry["record"]["item_id"] for entry in ledger_rows}
    unknown_overrides = set(overrides) - ledger_ids
    if unknown_overrides:
        raise ValueError(f"curation override IDs absent from ledger: {sorted(unknown_overrides)}")
    for entry in ledger_rows:
        row = entry["record"]
        item_id = row["item_id"]
        decision = overrides.get(
            item_id,
            {
                "relationship": "defer",
                "canonical_proposal_ids": [],
                "reason": "No atomic, acquisition-conditioned migration decision has been accepted.",
            },
        )
        if decision["relationship"] not in TERMINAL_RELATIONSHIPS:
            raise ValueError(f"invalid terminal relationship for {item_id}")
        aliases.append(
            {
                "schema_version": CURATION_SCHEMA_VERSION,
                "legacy_item_id": item_id,
                "ledger_line": entry["line"],
                "ledger_row_sha256": entry["row_sha256"],
                "relationship": decision["relationship"],
                "canonical_proposal_ids": decision.get("canonical_proposal_ids", []),
                "alias_of_item_id": decision.get("alias_of_item_id"),
                "superseded_by_item_id": decision.get("superseded_by_item_id"),
                "duplicate_observation_group": decision.get("duplicate_observation_group"),
                "field_overrides": decision.get("field_overrides", {}),
                "reason": decision["reason"],
                "source_locator": {
                    "system": row.get("system"),
                    "path": row.get("path"),
                    "kind": row.get("kind"),
                },
                "runtime_influence": False,
            }
        )
    return aliases


def build_quarantine(
    ledger_rows: list[dict[str, Any]], aliases: list[dict[str, Any]]
) -> list[dict[str, Any]]:
    by_id = {entry["record"]["item_id"]: entry for entry in ledger_rows}
    alias_by_id = {item["legacy_item_id"]: item for item in aliases}
    records: list[dict[str, Any]] = []
    for item_id in FORMULA_MASS_QUARANTINE_IDS:
        entry = by_id[item_id]
        row = entry["record"]
        observed = float(row["mz"])
        calculated = formula_mz(str(row["formula"]))
        records.append(
            {
                "schema_version": CURATION_SCHEMA_VERSION,
                "quarantine_id": f"Q-FORMULA-{entry['line']:04d}",
                "legacy_item_id": item_id,
                "category": "formula_mass_inconsistent",
                "ledger_line": entry["line"],
                "formula": row["formula"],
                "observed_mz": observed,
                "calculated_mz": round(calculated, 7),
                "mass_error_da": round(observed - calculated, 7),
                "correction_policy": "do_not_guess; recover original spectrum and atomic source",
                "reason": alias_by_id[item_id]["reason"],
                "runtime_influence": False,
            }
        )
    for item_id, category in (
        ("feature_profile:flavonoid:006", "exact_duplicate"),
        ("feature_profile:glycoside:017", "exact_duplicate"),
        ("legacy_rule:FLAV-026#03", "exact_duplicate"),
        ("legacy_rule:FLAV-027#03", "exact_duplicate"),
        ("legacy_rule:ALK-009b#02", "misplaced_overbroad_rule"),
        ("legacy_rule:LIP-002#01", "incorrect_class_assignment"),
    ):
        entry = by_id[item_id]
        records.append(
            {
                "schema_version": CURATION_SCHEMA_VERSION,
                "quarantine_id": f"Q-{category.upper().replace('_', '-')}-{entry['line']:04d}",
                "legacy_item_id": item_id,
                "category": category,
                "ledger_line": entry["line"],
                "reason": alias_by_id[item_id]["reason"],
                "runtime_influence": False,
            }
        )
    return sorted(records, key=lambda item: item["quarantine_id"])


def build_source_claims() -> list[dict[str, Any]]:
    common = {
        "schema_version": CURATION_SCHEMA_VERSION,
        "runtime_influence": False,
        "citation_presence_is_not_validation": True,
    }
    claims = [
        {
            **common,
            "source_claim_id": "SRC-MASSBANK-AG6560-PC-PE-PILOT",
            "source_type": "public_experimental_spectrum",
            "citation": "MassBank Europe Antwerp METOX phospholipid records",
            "locator": "tests/fixtures/feature_fragments/massbank_antwerp_phospholipids_v1.json",
            "url": "https://massbank.eu/MassBank-api/",
            "verification_status": "verified_full_record_fixture",
            "supported_claim": "PC product 184.0733 and PE loss 141.0191 occurrence at reported 10/20/40 eV on Agilent 6560 QTOF records whose contributor activation text is CID, positive [M+H]+.",
            "supported_domain": "Agilent 6560 QTOF; contributor text CID; fixed reported eV; positive [M+H]+",
            "limitations": [
                "One laboratory",
                "Ten two-dimensional identities",
                "Deterministic CI fixture only",
                "No PSI beam-type-CID term in the records",
                "No cross-platform CE transfer",
                "Class support only"
            ],
        },
        {
            **common,
            "source_claim_id": "SRC-NUCLEOSIDE-ORBITRAP-DIA-2023",
            "source_type": "peer_reviewed_primary",
            "citation": "Chem Res Toxicol. 2023; DOI 10.1021/acs.chemrestox.3c00041",
            "locator": "main text and supporting LC-Orbitrap DIA fragment tables",
            "url": "https://pubs.acs.org/doi/10.1021/acs.chemrestox.3c00041",
            "verification_status": "verified_primary_domain_limited",
            "supported_claim": "RNA-related 132 Da loss, DNA-related 116 Da loss, and protonated nucleobase product ions were observed in the reported Orbitrap DIA workflow.",
            "supported_domain": "LC-Orbitrap DIA as reported by the paper",
            "limitations": ["Not generic CID evidence", "Does not provide strict class-vs-background counts"],
        },
        {
            **common,
            "source_claim_id": "SRC-NUCLEOSIDE-CID-HCD-DOMAIN-2024",
            "source_type": "peer_reviewed_primary",
            "citation": "Automated identification of modified nucleosides; PMCID PMC11587168",
            "locator": "CID 42% and HCD 20 a.u. comparison",
            "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC11587168/",
            "verification_status": "verified_primary_domain_limited",
            "supported_claim": "Nucleoside fragment observations depend on activation and instrument conditions.",
            "supported_domain": "paper-specific CID and HCD conditions",
            "limitations": ["Percentage CID and HCD settings are not interchangeable"],
        },
        {
            **common,
            "source_claim_id": "SRC-ALPHA-AA-NL46-2024",
            "source_type": "peer_reviewed_primary",
            "citation": "PMID 39563088",
            "locator": "abstract and reported alpha/beta/gamma/epsilon amino-acid fragmentation comparison",
            "url": "https://pubmed.ncbi.nlm.nih.gov/39563088/",
            "verification_status": "verified_primary_domain_limited",
            "supported_claim": "Protonated alpha-amino acids show a 46 Da H2O+CO channel, while positional amino-acid classes differ.",
            "supported_domain": "reported tandem-MS/MSn conditions",
            "limitations": ["Observed CH2O2 mass does not prove mechanism", "No strict public target/background recount yet"],
        },
        {
            **common,
            "source_claim_id": "SRC-PE-196-RCM-2023",
            "source_type": "peer_reviewed_primary",
            "citation": "Ventura et al. Rapid Commun Mass Spectrom. 2023; DOI 10.1002/rcm.9527",
            "locator": "m/z 196.038 [C5H11NO5P]-; negative HILIC-AIF and selected-precursor LIT-CID experiments",
            "url": "https://analyticalsciencejournals.onlinelibrary.wiley.com/doi/full/10.1002/rcm.9527",
            "verification_status": "verified_primary_domain_limited",
            "supported_claim": "The m/z 196.038 phosphoethanolamine-related product is observed under the paper's negative-ion conditions.",
            "supported_domain": "paper-specific negative HILIC-AIF and LIT-CID",
            "limitations": ["PE/L-NAPE and related isomers require RT/other ions", "AIF is not precursor-selected MS/MS"],
        },
        {
            **common,
            "source_claim_id": "SRC-LIPID-HEADGROUP-AIF-2020",
            "source_type": "peer_reviewed_primary",
            "citation": "PMID 32423109; PMCID PMC7287777",
            "locator": "negative-ion HILIC-ESI-FTMS all-ion fragmentation headgroup ions",
            "url": "https://pubmed.ncbi.nlm.nih.gov/32423109/",
            "verification_status": "verified_primary_domain_limited",
            "supported_claim": "Negative lipid headgroup ions including m/z 224, 196, and 241 were used in the reported AIF workflow.",
            "supported_domain": "paper-specific HILIC-ESI-FTMS AIF",
            "limitations": ["AIF cannot establish precursor-linked absence", "No cross-platform portability claim"],
        },
        {
            **common,
            "source_claim_id": "SRC-PEPTIDE-IMMONIUM-QSTAR-2008",
            "source_type": "peer_reviewed_primary",
            "citation": "Hohmann et al.; PMCID PMC2638499",
            "locator": "ABI QSTAR Pulsar i peptide QTOF-CID immonium-ion study",
            "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC2638499/",
            "verification_status": "verified_scope_mismatch_with_legacy",
            "supported_claim": "Selected immonium ions correlate with peptide residue composition under the reported QTOF-CID conditions.",
            "supported_domain": "peptide QTOF-CID",
            "limitations": ["Not a free amino-acid validation", "Leu/Ile cannot be distinguished", "Arg/Lys not supported by the cited quantitative set"],
        },
        {
            **common,
            "source_claim_id": "SRC-OGLYCOSIDE-LEGACY-REFERENCES-PENDING",
            "source_type": "peer_reviewed_primary",
            "citation": "Legacy PMIDs 15389856 and 16418439",
            "locator": "exact supporting passages not yet captured",
            "url": "https://pubmed.ncbi.nlm.nih.gov/15389856/",
            "verification_status": "pending_exact_claim_verification",
            "supported_claim": "No claim accepted until exact loss, polarity, adduct, activation, and CE passages are recorded.",
            "supported_domain": None,
            "limitations": ["Citation presence alone is insufficient"],
        },
        {
            **common,
            "source_claim_id": "SRC-CGLYCOSIDE-LEGACY-REFERENCES-PENDING",
            "source_type": "peer_reviewed_primary",
            "citation": "Legacy PMIDs 25981540 and 15389856",
            "locator": "exact 60/90/120 passages and domains not yet captured",
            "url": "https://pubmed.ncbi.nlm.nih.gov/25981540/",
            "verification_status": "pending_exact_claim_verification",
            "supported_claim": "No claim accepted until the cross-ring series and exact acquisition domain are verified.",
            "supported_domain": None,
            "limitations": ["The three losses are correlated, not three independent votes"],
        },
        {
            **common,
            "source_claim_id": "SRC-SULFATE-SO3-SOURCE-PENDING",
            "source_type": "manual_rule",
            "citation": "Legacy sulfate-esters rule",
            "locator": "no atomic primary source accepted",
            "url": None,
            "verification_status": "unverified_manual_rule",
            "supported_claim": "None; acquisition-conditioned public spectra and a primary source are required.",
            "supported_domain": None,
            "limitations": ["Current source is insufficient for scoring"],
        },
        {
            **common,
            "source_claim_id": "SRC-SPHINGOLIPID-264-LEGACY",
            "source_type": "peer_reviewed_primary",
            "citation": "Legacy PMID 11070861 and PMCID PMC2724058",
            "locator": "d18:1-related m/z 264 passages require domain-specific recount",
            "url": "https://pubmed.ncbi.nlm.nih.gov/11070861/",
            "verification_status": "bibliography_verified_statistics_missing",
            "supported_claim": "m/z 264 is a d18:1-related sphingoid product under reported conditions.",
            "supported_domain": "paper-specific positive-ion collision activation",
            "limitations": ["Not universal to all sphingolipids", "d18:0/d18:2 and other chain lengths are required backgrounds"],
        },
    ]
    return sorted(claims, key=lambda item: item["source_claim_id"])


def _statistics_plan(
    plan_id: str,
    *,
    status: str,
    target: str,
    background: list[str],
    special: list[str],
    fixture: str | None = None,
) -> dict[str, Any]:
    return {
        "statistics_plan_id": plan_id,
        "status": status,
        "runtime_influence": False,
        "fixture": fixture,
        "target_definition": target,
        "background_definitions": background,
        "stratify_exactly_by": [
            "vendor",
            "instrument_model",
            "analyzer",
            "activation",
            "collision_cell_region",
            "collision_energy_raw",
            "collision_energy_unit",
            "collision_energy_mode",
            "adduct",
            "polarity",
            "precursor_charge",
        ],
        "counting_unit": "independent_2d_identity",
        "identity_key": "standardized_2d_connectivity",
        "formula_family_key": "charge-neutralized Hill formula",
        "deduplication": {
            "mirror_key": "base-peak-normalized peak-content SHA-256",
            "technical_replicates": "collapse within identity x method x CE before identity counting",
            "source_family": "record repository lineage; mirrored exports do not create a second identity",
        },
        "holdout": {
            "discovery_validation_identity_overlap_allowed": False,
            "discovery_validation_formula_family_overlap_allowed": False,
            "discovery_validation_source_family_overlap_allowed": False,
        },
        "required_outputs": [
            "target_present",
            "target_total",
            "target_frequency",
            "background_present",
            "background_total",
            "background_frequency",
            "conditional_likelihood_ratio_with_interval",
            "identity_formula_source_split_manifest_sha256",
        ],
        "negative_evidence_quality_gates": [
            "complete_mass_range",
            "centroided",
            "sufficient_signal",
            "precursor_isolated",
            "collision_energy_verified",
            "instrument_domain_verified",
        ],
        "special_requirements": special,
        "fixed_fragment_point_values_allowed": False,
    }


def build_statistics_plans() -> dict[str, Any]:
    plans = [
        _statistics_plan(
            "STAT-PC184-AG6560-PUBLIC-PILOT",
            status="complete_public_pilot_single_laboratory",
            target="phosphocholine-containing target identities represented as PC in the fixture",
            background=["glycerophosphoethanolamines", "glycerophosphoserines"],
            special=["Do not claim specificity against sphingomyelin until SM is an explicit background"],
            fixture="tests/fixtures/feature_fragments/massbank_antwerp_phospholipids_v1.json",
        ),
        _statistics_plan(
            "STAT-PE141-AG6560-PUBLIC-PILOT",
            status="complete_public_pilot_single_laboratory",
            target="protonated glycerophosphoethanolamine identities",
            background=["glycerophosphocholines", "glycerophosphoserines"],
            special=["Report the 2/3 target frequency at 40 eV; absence is weak evidence"],
            fixture="tests/fixtures/feature_fragments/massbank_antwerp_phospholipids_v1.json",
        ),
        _statistics_plan(
            "STAT-NUCLEOSIDE-POS-H",
            status="planned_no_strict_counts",
            target="each base x ribose/deoxyribose class counted separately",
            background=["other nucleobases", "modified nucleosides", "nucleotides", "O-glycosides sharing loss 132"],
            special=["Separate loss 132, loss 116, and each protonated-base product", "Do not treat pooled contrastive_v2 discovery counts as validation"],
        ),
        _statistics_plan(
            "STAT-OGLYCOSIDE-BY-DOMAIN",
            status="planned_no_strict_counts",
            target="O-glycosides separated by sugar residue and polarity/adduct",
            background=["C-glycosides", "other sugar types", "non-glycosides matched by aglycone/formula where possible"],
            special=["Loss 132 shares an evidence group with ribonucleosides", "Do not infer that every residue loss occurs in both polarities"],
        ),
        _statistics_plan(
            "STAT-CGLYCOSIDE-NEG-H",
            status="planned_no_strict_counts",
            target="phenolic/flavonoid C-glycosides under negative [M-H]-",
            background=["O-glycosides", "non-glycosides", "same-aglycone analogs"],
            special=["60/90/120 losses share one independent fragmentation group"],
        ),
        _statistics_plan(
            "STAT-AA-ALPHA-NL46-POS-H",
            status="planned_no_strict_counts",
            target="protonated alpha-amino-acid identities",
            background=["beta-amino acids", "gamma-amino acids", "epsilon-amino acids", "other CH2O2-loss compounds"],
            special=["HCOOH and H2O+CO labels share one observed-mass event", "Record MSn stage"],
        ),
        _statistics_plan(
            "STAT-SULFATE-SO3-NEG-H",
            status="planned_no_strict_counts",
            target="deprotonated sulfate esters",
            background=["sulfonates", "phosphates", "sulfites", "nonsulfated analogs"],
            special=["SO3 loss and m/z 80/97 products share one headgroup evidence group"],
        ),
        _statistics_plan(
            "STAT-LIPID-HEADGROUP-BY-DOMAIN",
            status="planned_no_strict_counts",
            target="each d18:1, PC224-adduct, PE196, or PI241 structural class separately",
            background=["other glycerophospholipid classes", "NAPE/L-NAPE", "other sphingoid chain lengths", "phosphorylated sugars"],
            special=["Separate PC224 formate and acetate adducts", "All-ion-fragmentation records cannot establish precursor-linked absence"],
        ),
        _statistics_plan(
            "STAT-PEPTIDE-IMMONIUM-QTOF-CID",
            status="literature_scope_repair_required",
            target="peptides containing each named residue under exact QTOF-CID conditions",
            background=["peptides lacking the residue", "isomeric residue contexts", "non-peptide compounds sharing the product m/z"],
            special=["Do not reuse these statistics for free amino acids", "Leu/Ile remain an unresolved isomer group", "Arg/Lys require an additional primary source"],
        ),
    ]
    return {
        "schema_version": CURATION_SCHEMA_VERSION,
        "runtime_influence": False,
        "plans": sorted(plans, key=lambda item: item["statistics_plan_id"]),
    }


def build_ledger_lock() -> dict[str, Any]:
    return {
        "schema_version": CURATION_SCHEMA_VERSION,
        "lock_id": "v15-feature-fragment-legacy-ledger-lock-20260826",
        "runtime_influence": False,
        "ledger": {
            "path": LEDGER_RELATIVE.as_posix(),
            "sha256": EXPECTED_LEDGER_SHA256,
            "record_count": EXPECTED_LEDGER_RECORDS,
            "format": "JSON Lines; one immutable legacy occurrence per line",
        },
        "source_files": SOURCE_FILE_HASHES,
        "coverage_contract": {
            "each_ledger_item_has_exactly_one_terminal_relationship": True,
            "unknown_item_ids_allowed": False,
            "default_relationship": "defer",
            "runtime_catalog_auto_load": False,
        },
    }


def _alias_graph_is_acyclic(aliases: list[dict[str, Any]]) -> tuple[bool, int]:
    graph: dict[str, list[str]] = {}
    item_ids = {item["legacy_item_id"] for item in aliases}
    edge_count = 0
    for item in aliases:
        targets = [
            item.get("alias_of_item_id"),
            item.get("superseded_by_item_id"),
        ]
        clean = [target for target in targets if target]
        unknown = set(clean) - item_ids
        if unknown:
            raise ValueError(
                f"alias target absent from ledger for {item['legacy_item_id']}: {sorted(unknown)}"
            )
        graph[item["legacy_item_id"]] = clean
        edge_count += len(clean)

    visiting: set[str] = set()
    visited: set[str] = set()

    def visit(node: str) -> None:
        if node in visiting:
            raise ValueError(f"alias DAG cycle at {node}")
        if node in visited:
            return
        visiting.add(node)
        for target in graph.get(node, []):
            visit(target)
        visiting.remove(node)
        visited.add(node)

    for node in graph:
        visit(node)
    return True, edge_count


def build_assets(ledger: Path) -> dict[str, bytes]:
    ledger_rows = load_ledger(ledger)
    proposals = build_proposals()
    aliases = build_aliases(ledger_rows, proposals)
    quarantine = build_quarantine(ledger_rows, aliases)
    source_claims = build_source_claims()
    statistics_plans = build_statistics_plans()
    lock = build_ledger_lock()

    ledger_ids = {entry["record"]["item_id"] for entry in ledger_rows}
    alias_ids = [item["legacy_item_id"] for item in aliases]
    if len(alias_ids) != EXPECTED_LEDGER_RECORDS or set(alias_ids) != ledger_ids:
        raise ValueError("legacy alias overlay does not cover the frozen ledger exactly")
    if len(alias_ids) != len(set(alias_ids)):
        raise ValueError("legacy alias overlay contains duplicate terminal decisions")
    proposal_ids = {item["proposal_id"] for item in proposals}
    for alias in aliases:
        unknown = set(alias["canonical_proposal_ids"]) - proposal_ids
        if unknown:
            raise ValueError(f"unknown proposal link from {alias['legacy_item_id']}")
    for proposal in proposals:
        unknown = set(proposal["legacy_item_ids"]) - ledger_ids
        if unknown:
            raise ValueError(f"proposal references unknown ledger IDs: {sorted(unknown)}")
        for item_id in proposal["legacy_item_ids"]:
            alias = aliases[next(i for i, item in enumerate(aliases) if item["legacy_item_id"] == item_id)]
            if proposal["proposal_id"] not in alias["canonical_proposal_ids"]:
                raise ValueError(
                    f"proposal/alias link is not reciprocal: {proposal['proposal_id']} -> {item_id}"
                )
    proposal_by_id = {item["proposal_id"]: item for item in proposals}
    for alias in aliases:
        for proposal_id in alias["canonical_proposal_ids"]:
            if alias["legacy_item_id"] not in proposal_by_id[proposal_id]["legacy_item_ids"]:
                raise ValueError(
                    f"alias/proposal link is not reciprocal: {alias['legacy_item_id']} -> {proposal_id}"
                )
    acyclic, edge_count = _alias_graph_is_acyclic(aliases)
    formula_quarantine = [
        item for item in quarantine if item["category"] == "formula_mass_inconsistent"
    ]
    if {item["legacy_item_id"] for item in formula_quarantine} != set(
        FORMULA_MASS_QUARANTINE_IDS
    ):
        raise ValueError("formula/mass quarantine set mismatch")
    if any(abs(float(item["mass_error_da"])) <= 0.01 for item in formula_quarantine):
        raise ValueError("formula/mass quarantine contains an in-tolerance record")

    assets: dict[str, bytes] = {
        "ledger_lock_v1.json": json_bytes(lock),
        "legacy_aliases_v1.jsonl": jsonl_bytes(aliases),
        "atomicization_queue_v1.jsonl": jsonl_bytes(proposals),
        "quarantine_v1.jsonl": jsonl_bytes(quarantine),
        "source_claims_v1.jsonl": jsonl_bytes(source_claims),
        "statistics_plans_v1.json": json_bytes(statistics_plans),
    }
    relationship_counts = Counter(item["relationship"] for item in aliases)
    status_counts = Counter(item["status"] for item in proposals)
    receipt = {
        "schema_version": CURATION_SCHEMA_VERSION,
        "receipt_id": "v15-feature-fragment-curation-overlay-20260826",
        "generated_on": CURATION_DATE,
        "generated_by": "scripts/curate_feature_fragment_legacy.py",
        "runtime_influence": False,
        "ledger": lock["ledger"],
        "coverage": {
            "ledger_item_count": len(ledger_ids),
            "terminal_decision_count": len(alias_ids),
            "unique_terminal_item_count": len(set(alias_ids)),
            "orphan_ledger_item_count": len(ledger_ids - set(alias_ids)),
            "unknown_overlay_item_count": len(set(alias_ids) - ledger_ids),
            "terminal_relationship_counts": dict(sorted(relationship_counts.items())),
            "terminal_relationships_mutually_exclusive": True,
        },
        "alias_graph": {
            "edge_count": edge_count,
            "acyclic": acyclic,
        },
        "formula_mass_quarantine": {
            "required_count": 9,
            "actual_count": len(formula_quarantine),
            "legacy_item_ids": sorted(FORMULA_MASS_QUARANTINE_IDS),
            "all_outside_0_01_da": True,
        },
        "proposals": {
            "count": len(proposals),
            "status_counts": dict(sorted(status_counts.items())),
            "accepted_shadow_rule_ids": sorted(
                item["existing_runtime_rule_id"]
                for item in proposals
                if item["status"] == "accepted_shadow"
            ),
            "runtime_export_eligible_rule_ids": [],
            "export_policy": "no curation proposal is runtime-export eligible; six accepted_shadow entries are deterministic CI and research-audit records only",
        },
        "runtime_boundary": {
            "curation_directory_loaded_by_runtime_catalog": False,
            "runtime_catalog_files_modified": False,
            "exact_identity_upgrade_allowed": False,
            "standard_confirmation_upgrade_allowed": False,
            "isomer_discrimination_upgrade_allowed": False,
        },
        "asset_sha256": {
            name: sha256_bytes(payload) for name, payload in sorted(assets.items())
        },
    }
    assets["acceptance_receipt_v1.json"] = json_bytes(receipt)
    return assets


def write_assets(output_dir: Path, assets: dict[str, bytes]) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    for name, payload in sorted(assets.items()):
        (output_dir / name).write_bytes(payload)


def check_assets(output_dir: Path, assets: dict[str, bytes]) -> list[str]:
    errors: list[str] = []
    for name, expected in sorted(assets.items()):
        path = output_dir / name
        if not path.is_file():
            errors.append(f"missing asset: {path}")
            continue
        actual = path.read_bytes()
        if actual != expected:
            errors.append(
                f"stale asset: {path} sha256={sha256_bytes(actual)} expected={sha256_bytes(expected)}"
            )
    return errors


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ledger", type=Path, default=DEFAULT_LEDGER)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--write", action="store_true", help="write deterministic curation assets")
    mode.add_argument("--check", action="store_true", help="compare checked-in assets byte-for-byte")
    args = parser.parse_args()

    assets = build_assets(args.ledger)
    if args.write:
        write_assets(args.output_dir, assets)
    else:
        errors = check_assets(args.output_dir, assets)
        if errors:
            for error in errors:
                print(error, file=sys.stderr)
            return 1
    summary = {
        "mode": "write" if args.write else "check",
        "ledger_sha256": EXPECTED_LEDGER_SHA256,
        "ledger_record_count": EXPECTED_LEDGER_RECORDS,
        "asset_count": len(assets),
        "output_dir": str(args.output_dir),
    }
    print(json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
