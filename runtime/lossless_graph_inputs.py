"""Lossless molecular graph views and explicit, non-truncating input budgets.

This module is opt-in. It does not change the production ranker or its scores.
RDKit JSON provides the graph payload; the source isomeric identity remains a
separate, strictly checked field. Graph data are never reconstructed by an LLM.
"""
from __future__ import annotations

import hashlib
import json
from typing import Any

SCHEMA = "mssystem.lossless_structure_view.v1"


def _identity(mol: Any) -> str:
    from rdkit import Chem
    return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)


def _chemistry_snapshot(mol: Any) -> dict[str, Any]:
    """Atom order is fixed by the RDKit graph document; compare mapped chemistry."""
    from rdkit import Chem
    atoms = [{"index": a.GetIdx(), "element": a.GetAtomicNum(), "isotope": a.GetIsotope(),
              "charge": a.GetFormalCharge(), "radicals": a.GetNumRadicalElectrons(),
              "total_H": a.GetTotalNumHs(), "aromatic": a.GetIsAromatic(),
              "cip": a.GetProp("_CIPCode") if a.HasProp("_CIPCode") else None}
             for a in mol.GetAtoms()]
    bonds = [{"atoms": sorted([b.GetBeginAtomIdx(), b.GetEndAtomIdx()]),
              "type": str(b.GetBondType()), "aromatic": b.GetIsAromatic()}
             for b in mol.GetBonds()]
    return {"atoms": atoms, "bonds": sorted(bonds, key=lambda x: x["atoms"]),
            "canonical_cxsmiles": Chem.MolToCXSmiles(mol, isomericSmiles=True),
            "components": [list(x) for x in Chem.GetMolFrags(mol)]}


def build_structure_view(source_smiles: str) -> dict[str, Any]:
    """Return an explicitly validated full-identity RDKit graph and source text."""
    from rdkit import Chem, rdBase
    from rdkit.Chem import rdMolDescriptors

    if not isinstance(source_smiles, str) or not source_smiles.strip():
        raise ValueError("source_smiles must be a nonempty string")
    mol = Chem.MolFromSmiles(source_smiles)
    if mol is None:
        raise ValueError("source_smiles does not parse")
    graph_text = Chem.MolsToJSON([mol])
    graph = json.loads(graph_text)
    restored = Chem.JSONToMols(graph_text)
    if len(restored) != 1 or _identity(mol) != _identity(restored[0]):
        raise ValueError("RDKit graph changed full molecular identity")
    chemistry = _chemistry_snapshot(mol)
    if chemistry != _chemistry_snapshot(restored[0]):
        raise ValueError("RDKit graph changed mapped atoms, bonds or enhanced stereochemistry")
    canonical = _identity(mol)
    return {
        "schema": SCHEMA,
        "source_smiles": source_smiles,
        "source_smiles_sha256": hashlib.sha256(source_smiles.encode()).hexdigest(),
        "canonical_isomeric_smiles": canonical,
        "canonical_cxsmiles": chemistry["canonical_cxsmiles"],
        "source_atom_map": list(range(mol.GetNumAtoms())),
        "mapped_chemistry": chemistry,
        "full_inchikey": Chem.MolToInchiKey(mol),
        "identity_sha256": hashlib.sha256(canonical.encode()).hexdigest(),
        "rdkit_version": rdBase.rdkitVersion,
        "graph": graph,
        "formula": rdMolDescriptors.CalcMolFormula(mol),
        "exact_mass": rdMolDescriptors.CalcExactMolWt(mol),
        "formal_charge": Chem.GetFormalCharge(mol),
        "radical_electrons": sum(a.GetNumRadicalElectrons() for a in mol.GetAtoms()),
        "component_N": len(Chem.GetMolFrags(mol)),
        "atom_N": mol.GetNumAtoms(),
        "bond_N": mol.GetNumBonds(),
    }


def validate_structure_view(view: dict[str, Any]) -> None:
    """Reject a changed graph, identifier or source; do not silently normalize."""
    from rdkit import Chem
    if view.get("schema") != SCHEMA:
        raise ValueError("structure view schema mismatch")
    source = Chem.MolFromSmiles(view["source_smiles"])
    restored = Chem.JSONToMols(json.dumps(view["graph"], separators=(",", ":")))
    if source is None or len(restored) != 1:
        raise ValueError("structure view parse failed")
    canonical = _identity(source)
    chemistry = _chemistry_snapshot(source)
    if canonical != _identity(restored[0]) or canonical != view["canonical_isomeric_smiles"]:
        raise ValueError("structure graph/source identity mismatch")
    if hashlib.sha256(canonical.encode()).hexdigest() != view["identity_sha256"]:
        raise ValueError("structure identity hash mismatch")
    if hashlib.sha256(view["source_smiles"].encode()).hexdigest() != view["source_smiles_sha256"]:
        raise ValueError("source SMILES bytes changed")
    if chemistry != _chemistry_snapshot(restored[0]) or chemistry != view["mapped_chemistry"]:
        raise ValueError("mapped chemistry or enhanced stereochemistry changed")
    if view["source_atom_map"] != list(range(source.GetNumAtoms())):
        raise ValueError("source atom map mismatch")
    if chemistry["canonical_cxsmiles"] != view["canonical_cxsmiles"]:
        raise ValueError("canonical CXSMILES mismatch")
    if Chem.MolToInchiKey(source) != view["full_inchikey"]:
        raise ValueError("full InChIKey mismatch")


def pair_prompt(view_a: dict[str, Any], view_b: dict[str, Any], *,
                include_graph: bool = True, swapped_graph_control: bool = False) -> str:
    """Keep source identities in place; an explicit control can swap only graphs.

    The control is for frozen mechanism experiments, never production ranking.
    Both views must pass strict identity validation before any control is made.
    """
    validate_structure_view(view_a)
    validate_structure_view(view_b)
    if swapped_graph_control and not include_graph:
        raise ValueError("swapped graph control requires graph data")
    text = f"Candidate A: {view_a['source_smiles']}\nCandidate B: {view_b['source_smiles']}\n"
    if include_graph:
        ga, gb = (view_b, view_a) if swapped_graph_control else (view_a, view_b)
        text += (
            "RDKit graph A: " + json.dumps(ga["graph"], separators=(",", ":")) + "\n"
            "RDKit graph B: " + json.dumps(gb["graph"], separators=(",", ":")) + "\n"
        )
    return text + "Which candidate produced the spectrum? Answer:"


def tokenize_without_truncation(tokenizer: Any, text: str, *, prefix_rows: int,
                                max_total_tokens: int, native_max_positions: int) -> dict[str, Any]:
    """Fail before inference when any complete prompt exceeds a declared bound."""
    if not (isinstance(prefix_rows, int) and prefix_rows >= 0):
        raise ValueError("prefix_rows must be a nonnegative integer")
    if not (0 < max_total_tokens <= native_max_positions):
        raise ValueError("input budget exceeds native model positions")
    encoded = tokenizer(text, add_special_tokens=False, truncation=False,
                        return_attention_mask=False)
    ids = encoded["input_ids"]
    if not ids or len(ids) + prefix_rows > max_total_tokens:
        raise ValueError(f"complete input needs {len(ids) + prefix_rows} tokens; budget={max_total_tokens}")
    import struct
    raw = b"".join(struct.pack("<q", int(i)) for i in ids)
    return {"input_ids": ids, "prompt_token_N": len(ids),
            "prefix_rows": prefix_rows, "total_token_N": len(ids) + prefix_rows,
            "max_total_tokens": max_total_tokens,
            "final_hidden_absolute_index": len(ids) + prefix_rows - 1,
            "prompt_token_sha256": hashlib.sha256(raw).hexdigest(),
            "text_sha256": hashlib.sha256(text.encode()).hexdigest(),
            "truncated": False}
