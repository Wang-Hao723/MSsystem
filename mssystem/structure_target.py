from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, List, Optional

import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import rdMolDescriptors, AllChem


class RuleBasedTargetPredictor:
    def __init__(self, rule_map: Optional[Dict[str, str]] = None) -> None:
        self.rule_map = rule_map or {
            "c1cc2c(cc1)[nH]cc2": "Serotonin receptor",
            "NCCc1ccccc1": "Dopamine receptor",
            "c1ncc2n1ccc(c2)": "Kinase family",
            "c1ccc2c(c1)OCO2": "Nuclear receptor",
        }
        self.compiled = {
            smarts: Chem.MolFromSmarts(smarts) for smarts in self.rule_map.keys()
        }

    def predict(self, mol: Chem.Mol) -> List[str]:
        hits = []
        for smarts, patt in self.compiled.items():
            if patt is not None and mol.HasSubstructMatch(patt):
                hits.append(self.rule_map[smarts])
        return sorted(set(hits))


class StructureTargetPredictor:
    def __init__(self, known_target_df: Optional[pd.DataFrame] = None) -> None:
        self.known_target_df = known_target_df
        self.rule_predictor = RuleBasedTargetPredictor()

    def analyze_structure(self, smiles: str) -> Dict[str, Any]:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {"valid": False, "error": "Invalid SMILES"}

        logp, _ = rdMolDescriptors.CalcCrippenDescriptors(mol)
        desc = {
            "valid": True,
            "formula": rdMolDescriptors.CalcMolFormula(mol),
            "exact_mass": float(rdMolDescriptors.CalcExactMolWt(mol)),
            "logp": float(logp),
            "tpsa": float(rdMolDescriptors.CalcTPSA(mol)),
            "h_donors": int(rdMolDescriptors.CalcNumHBD(mol)),
            "h_acceptors": int(rdMolDescriptors.CalcNumHBA(mol)),
            "rotatable_bonds": int(rdMolDescriptors.CalcNumRotatableBonds(mol)),
            "lipinski_violations": self._lipinski_violations(mol),
        }
        return desc

    def predict_targets(self, smiles: str, top_k: int = 5) -> Dict[str, Any]:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {"similarity_targets": [], "rule_targets": []}

        sim_targets = self._similarity_targets(mol, top_k=top_k)
        rule_targets = self.rule_predictor.predict(mol)
        return {
            "similarity_targets": sim_targets,
            "rule_targets": rule_targets,
        }

    def _similarity_targets(
        self, query_mol: Chem.Mol, top_k: int = 5
    ) -> List[Dict[str, Any]]:
        if self.known_target_df is None or self.known_target_df.empty:
            return []

        import warnings

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            qfp = rdMolDescriptors.GetMorganFingerprintAsBitVect(
                query_mol, radius=2, nBits=2048
            )
        scored: List[Dict[str, Any]] = []

        for _, row in self.known_target_df.iterrows():
            smiles = row.get("smiles")
            target = row.get("target")
            if not smiles or not target:
                continue
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                continue

            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                rfp = rdMolDescriptors.GetMorganFingerprintAsBitVect(
                    mol, radius=2, nBits=2048
                )
            sim = DataStructs.TanimotoSimilarity(qfp, rfp)
            scored.append(
                {
                    "target": target,
                    "similarity": float(sim),
                    "reference_smiles": smiles,
                    "reference_name": row.get("compound_name"),
                }
            )

        scored.sort(key=lambda x: x["similarity"], reverse=True)
        return scored[:top_k]

    @staticmethod
    def _lipinski_violations(mol: Chem.Mol) -> int:
        mw = rdMolDescriptors.CalcExactMolWt(mol)
        logp, _ = rdMolDescriptors.CalcCrippenDescriptors(mol)
        hbd = rdMolDescriptors.CalcNumHBD(mol)
        hba = rdMolDescriptors.CalcNumHBA(mol)

        violations = 0
        if mw > 500:
            violations += 1
        if logp > 5:
            violations += 1
        if hbd > 5:
            violations += 1
        if hba > 10:
            violations += 1
        return violations

    def generate_selfies(self, smiles: str) -> Dict[str, Any]:
        """Convert SMILES to SELFIES representation.

        Returns a dict with 'valid', 'smiles', 'selfies', 'selfies_alphabet',
        'encoded_alphabet', and an optional 'accuracy' (fraction of atoms
        successfully encoded) if a ground-truth SMILES is provided.
        """
        import selfies as sf

        if not smiles:
            return {"valid": False, "error": "No SMILES provided"}

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {"valid": False, "error": "Invalid SMILES"}

        try:
            selfies_str = sf.encoder(smiles)
        except Exception as e:
            return {"valid": False, "error": f"SELFIES encoding failed: {e}"}

        if selfies_str is None or selfies_str == "":
            return {"valid": False, "error": "SELFIES encoder returned empty string"}

        # Alphabet tokens (each [X] chunk — split_selfies returns a generator)
        try:
            alphabet = list(sf.split_selfies(selfies_str))
        except Exception as e:
            return {"valid": False, "error": f"SELFIES tokenization failed: {e}"}

        # Encoded as integer labels using the default vocabulary
        try:
            encoded = list(sf.label_selfies(selfies_str))
        except Exception:
            encoded = []

        return {
            "valid": True,
            "smiles": smiles,
            "selfies": selfies_str,
            "selfies_length": len(alphabet),
            "selfies_tokens": alphabet,
            "num_atoms": mol.GetNumAtoms(),
            "encoded_indices": encoded,
        }
