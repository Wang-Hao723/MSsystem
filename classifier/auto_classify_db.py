"""Apply a structure classifier with explicit validation and abstention gates.

Unvalidated or legacy model bundles are rejected by default.  A researcher may
opt in with ``--allow-research-model``; such predictions are written with the
untrusted source ``ml_research_unvalidated`` and therefore cannot unlock
class-specific production scoring or confidence rules.
"""
from __future__ import annotations

import argparse
import math
import sqlite3
import sys
import time
from pathlib import Path

import joblib
import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
RDLogger.DisableLog("rdApp.*")

from matching.fragmentation_rules import normalize_classification
from mssystem.classification import RuleBasedClassifier


RC = RuleBasedClassifier()
MODEL_PATH = ROOT / "classifier" / "struct_clf.joblib"
CLASS_MERGE = {"fatty_acid": "lipid"}


def validate_model_bundle(bundle: dict, *, allow_research: bool = False) -> str:
    """Validate artifact semantics and return the ML provenance label."""
    if not isinstance(bundle, dict) or int(bundle.get("schema_version", 0)) < 2:
        raise ValueError(
            "Legacy structure model lacks dedup/scaffold/calibration provenance"
        )
    required = {
        "model", "classes", "fp_bits", "fp_radius", "temperature",
        "class_thresholds", "split_audit", "label_provenance",
    }
    missing = sorted(required.difference(bundle))
    if missing:
        raise ValueError(f"Structure model bundle is incomplete: {missing}")
    audit = bundle.get("split_audit") or {}
    overlap_fields = (
        "train_calibration_scaffold_overlap",
        "train_test_scaffold_overlap",
        "calibration_test_scaffold_overlap",
    )
    if any(int(audit.get(field, -1)) != 0 for field in overlap_fields):
        raise ValueError("Structure model scaffold split audit is not disjoint")
    if bundle.get("production_eligible") is True:
        return "ml_calibrated_qualified"
    if not allow_research:
        raise PermissionError(
            "Structure model is research-only; pass --allow-research-model "
            "to write explicitly untrusted research labels"
        )
    return "ml_research_unvalidated"


def calibrate_probabilities(raw_probabilities: np.ndarray, temperature: float) -> np.ndarray:
    clipped = np.clip(np.asarray(raw_probabilities, dtype=float), 1e-9, 1.0)
    logits = np.log(clipped) / max(float(temperature), 1e-6)
    logits -= logits.max(axis=1, keepdims=True)
    exp = np.exp(logits)
    return exp / exp.sum(axis=1, keepdims=True)


def selective_classifications(
    probabilities: np.ndarray,
    classes: list[str],
    class_thresholds: dict[str, float],
) -> list[tuple[str, float, bool]]:
    output = []
    for row in np.asarray(probabilities):
        class_index = int(np.argmax(row))
        label = classes[class_index]
        confidence = float(row[class_index])
        threshold = float(class_thresholds.get(label, 1.01))
        accepted = bool(
            math.isfinite(confidence)
            and confidence >= threshold
            and threshold <= 1.0
        )
        output.append((label if accepted else "other", confidence, accepted))
    return output


def _add_columns(cursor) -> None:
    existing = {row[1] for row in cursor.execute("PRAGMA table_info(compounds)")}
    for column, declaration in (
        ("classification_v3", "TEXT"),
        ("classification_source", "TEXT"),
        ("classification_confidence", "REAL"),
    ):
        if column not in existing:
            cursor.execute(
                f"ALTER TABLE compounds ADD COLUMN {column} {declaration}"
            )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--full", action="store_true")
    parser.add_argument("--model", type=Path, default=MODEL_PATH)
    parser.add_argument(
        "--allow-research-model",
        action="store_true",
        help="explicitly write untrusted research-only ML labels",
    )
    args = parser.parse_args()

    bundle = joblib.load(args.model)
    ml_source = validate_model_bundle(
        bundle, allow_research=args.allow_research_model
    )
    classifier = bundle["model"]
    classes = list(bundle["classes"])
    fp_bits = int(bundle["fp_bits"])
    fp_radius = int(bundle["fp_radius"])
    temperature = float(bundle["temperature"])
    thresholds = {
        str(label): float(value)
        for label, value in (bundle.get("class_thresholds") or {}).items()
    }

    database = ROOT / "data" / (
        "local_library.db" if args.full else "local_library_small.db"
    )
    connection = sqlite3.connect(str(database))
    cursor = connection.cursor()
    _add_columns(cursor)
    connection.commit()
    rows = cursor.execute(
        "SELECT id, smiles, classification FROM compounds"
    ).fetchall()

    from rdkit.DataStructs import ConvertToNumpyArray

    updates = []
    pending = []
    fingerprints = []
    counts = {
        "smarts_structure_rule": 0,
        "existing_unverified": 0,
        ml_source: 0,
        "other": 0,
    }
    started = time.time()
    for compound_id, smiles, existing_class in rows:
        detailed = RC.predict_detailed(smiles)
        if detailed["class"] != "Unknown":
            updates.append((
                detailed["class"], "smarts_structure_rule", 1.0, compound_id
            ))
            counts["smarts_structure_rule"] += 1
            continue
        normalized = normalize_classification(existing_class or "")
        if normalized and normalized != "Unknown":
            updates.append((
                normalized, "existing_unverified", None, compound_id
            ))
            counts["existing_unverified"] += 1
            continue
        molecule = Chem.MolFromSmiles(smiles) if smiles else None
        if molecule is None:
            updates.append(("other", "other", 0.0, compound_id))
            counts["other"] += 1
            continue
        bit_vector = AllChem.GetMorganFingerprintAsBitVect(
            molecule, fp_radius, nBits=fp_bits
        )
        array = np.zeros(fp_bits, dtype=np.uint8)
        ConvertToNumpyArray(bit_vector, array)
        pending.append(compound_id)
        fingerprints.append(array)

    if fingerprints:
        raw_probability = classifier.predict_proba(np.stack(fingerprints))
        calibrated = calibrate_probabilities(raw_probability, temperature)
        decisions = selective_classifications(calibrated, classes, thresholds)
        for compound_id, (label, confidence, accepted) in zip(pending, decisions):
            if accepted:
                updates.append((label, ml_source, confidence, compound_id))
                counts[ml_source] += 1
            else:
                updates.append(("other", "other", confidence, compound_id))
                counts["other"] += 1

    updates = [
        (CLASS_MERGE.get(label, label), source, confidence, compound_id)
        for label, source, confidence, compound_id in updates
    ]
    cursor.executemany(
        "UPDATE compounds SET classification_v3=?, "
        "classification_source=?, classification_confidence=? WHERE id=?",
        updates,
    )
    connection.commit()
    connection.close()
    print(
        f"Updated {len(updates)} compounds in {time.time()-started:.0f}s; "
        f"sources={counts}; ML source={ml_source}"
    )


if __name__ == "__main__":
    main()
