"""RDKit-defined chemical identity boundaries for blind validation.

The benchmark unit is two-dimensional molecular connectivity.  Full
stereochemical identity remains recorded, but ordinary achiral LC-MS/MS is not
allowed to claim separation of structures that share the same connectivity
layer.
"""

from __future__ import annotations

import hashlib
from typing import Any

import rdkit
from rdkit import Chem
from rdkit.Chem import inchi, rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold

from matching.library_quality import normalize_formula_text, normalize_inchikey


def _fallback_connectivity_key(canonical_smiles: str) -> str:
    digest = hashlib.sha256(canonical_smiles.encode("utf-8")).hexdigest()
    return f"rdkit-smiles-sha256::{digest}"


def rdkit_identity_profile(
    smiles: Any,
    *,
    declared_formula: Any = "",
    declared_inchikey: Any = "",
    require_single_component: bool = True,
) -> dict[str, Any]:
    """Return a deterministic, auditable structure identity profile.

    Invalid structures and metadata conflicts return ``valid=False`` instead
    of being silently repaired.  This is intentional for authentic-standard
    truth and candidate-pair construction.
    """

    raw_smiles = str(smiles or "").strip()
    formula = normalize_formula_text(declared_formula)
    supplied_key = normalize_inchikey(declared_inchikey)
    base = {
        "input_smiles": raw_smiles,
        "declared_formula": formula,
        "declared_inchikey": supplied_key,
        "rdkit_version": rdkit.__version__,
        "identity_unit": "two_dimensional_connectivity",
        "stereochemistry_claimed_by_msms": False,
    }
    if not raw_smiles:
        return {**base, "valid": False, "reason": "missing_smiles"}
    try:
        mol = Chem.MolFromSmiles(raw_smiles, sanitize=True)
    except Exception:
        mol = None
    if mol is None:
        return {**base, "valid": False, "reason": "invalid_or_unsanitizable_smiles"}
    component_count = len(Chem.GetMolFrags(mol))
    if require_single_component and component_count != 1:
        return {
            **base,
            "valid": False,
            "reason": "multicomponent_structure_requires_explicit_parent_selection",
            "component_count": component_count,
        }
    try:
        canonical_isomeric = Chem.MolToSmiles(
            mol, canonical=True, isomericSmiles=True
        )
        canonical_connectivity = Chem.MolToSmiles(
            mol, canonical=True, isomericSmiles=False
        )
        structure_formula = normalize_formula_text(
            rdMolDescriptors.CalcMolFormula(mol)
        )
        derived_inchi = inchi.MolToInchi(mol)
        derived_key = normalize_inchikey(inchi.InchiToInchiKey(derived_inchi))
        scaffold_mol = MurckoScaffold.GetScaffoldForMol(mol)
        scaffold_smiles = (
            Chem.MolToSmiles(scaffold_mol, canonical=True, isomericSmiles=False)
            if scaffold_mol.GetNumAtoms()
            else ""
        )
        chiral_centers = Chem.FindMolChiralCenters(
            mol,
            includeUnassigned=True,
            useLegacyImplementation=False,
        )
    except Exception as exc:
        return {
            **base,
            "valid": False,
            "reason": "rdkit_identity_derivation_failed",
            "error": f"{type(exc).__name__}: {exc}"[:300],
        }
    if formula and structure_formula and formula != structure_formula:
        return {
            **base,
            "valid": False,
            "reason": "declared_formula_structure_conflict",
            "structure_formula": structure_formula,
            "canonical_isomeric_smiles": canonical_isomeric,
            "canonical_connectivity_smiles": canonical_connectivity,
        }
    if supplied_key and derived_key and supplied_key != derived_key:
        return {
            **base,
            "valid": False,
            "reason": "declared_inchikey_structure_conflict",
            "derived_inchikey": derived_key,
            "canonical_isomeric_smiles": canonical_isomeric,
            "canonical_connectivity_smiles": canonical_connectivity,
        }
    connectivity_key = (
        derived_key.split("-", 1)[0]
        if derived_key
        else _fallback_connectivity_key(canonical_connectivity)
    )
    scaffold_key = (
        f"murcko::{scaffold_smiles}"
        if scaffold_smiles
        else f"acyclic::{canonical_connectivity}"
    )
    return {
        **base,
        "valid": True,
        "reason": "valid_rdkit_identity",
        "component_count": component_count,
        "canonical_isomeric_smiles": canonical_isomeric,
        "canonical_connectivity_smiles": canonical_connectivity,
        "structure_formula": structure_formula,
        "normalized_formula": formula or structure_formula,
        "inchi": derived_inchi,
        "inchikey": supplied_key or derived_key,
        "connectivity_key": connectivity_key,
        "murcko_scaffold_smiles": scaffold_smiles,
        "scaffold_key": scaffold_key,
        "assigned_chiral_centers": sum(
            label != "?" for _, label in chiral_centers
        ),
        "unassigned_chiral_centers": sum(
            label == "?" for _, label in chiral_centers
        ),
    }


def same_2d_connectivity(left: dict[str, Any], right: dict[str, Any]) -> bool:
    """Return true only for two valid profiles sharing RDKit/InChI connectivity."""

    return bool(
        left.get("valid")
        and right.get("valid")
        and left.get("connectivity_key")
        and left.get("connectivity_key") == right.get("connectivity_key")
    )


__all__ = ["rdkit_identity_profile", "same_2d_connectivity"]
