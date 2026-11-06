"""Train the research structure classifier without structure-family leakage.

The historical trainer randomly split rows carrying ECFP fingerprints.  Exact
duplicates and close scaffold analogues could therefore appear on both sides,
while SMARTS-generated labels made the reported validation score partly
circular.  This trainer instead:

* canonicalizes and deduplicates structures before any split;
* partitions by Bemis--Murcko scaffold into train/calibration/test;
* uses a separate scaffold holdout for early stopping;
* temperature-calibrates on calibration only;
* learns per-class selective thresholds with a Wilson precision lower bound;
* records weak-label provenance and always saves DB/SMARTS models as
  research-only until an external governance review promotes the artifact.

No model produced here is allowed to influence production merely because its
weak-label agreement is high.
"""
from __future__ import annotations

import argparse
import copy
import csv
import json
import math
import sqlite3
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

import joblib
import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn.metrics import classification_report, f1_score, log_loss
from sklearn.model_selection import GroupShuffleSplit, StratifiedGroupKFold
from sklearn.neural_network import MLPClassifier

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
RDLogger.DisableLog("rdApp.*")

from matching.fragmentation_rules import normalize_classification
from mssystem.classification import RuleBasedClassifier


RC = RuleBasedClassifier()
FP_BITS = 2048
FP_RADIUS = 2
MIN_CLASS_SAMPLES = 50
MIN_CLASS_SCAFFOLDS = 5
MODEL_PATH = ROOT / "classifier" / "struct_clf.joblib"
TARGET_SELECTIVE_PRECISION = 0.90
MIN_CALIBRATION_CALLS = 10
RANDOM_STATE = 42


def canonical_smiles(smiles: str) -> str | None:
    mol = Chem.MolFromSmiles(str(smiles or ""))
    if mol is None:
        return None
    return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)


def scaffold_key(smiles: str) -> str:
    """Return a stable structure-family group for leakage-safe splitting."""
    mol = Chem.MolFromSmiles(str(smiles or ""))
    if mol is None:
        return "invalid"
    scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    if scaffold is not None and scaffold.GetNumAtoms() > 0:
        generic = MurckoScaffold.MakeScaffoldGeneric(scaffold)
        return "murcko:" + Chem.MolToSmiles(generic, canonical=True)
    # Murcko is empty for acyclic molecules.  Formula grouping is deliberately
    # conservative: exact structures are already deduplicated, while common
    # acyclic homologues do not get scattered by arbitrary row order.
    return "acyclic_formula:" + rdMolDescriptors.CalcMolFormula(mol)


def smiles_to_fp(smiles: str) -> np.ndarray | None:
    mol = Chem.MolFromSmiles(str(smiles or ""))
    if mol is None:
        return None
    bit_vector = AllChem.GetMorganFingerprintAsBitVect(
        mol, FP_RADIUS, nBits=FP_BITS
    )
    array = np.zeros(FP_BITS, dtype=np.uint8)
    from rdkit.DataStructs import ConvertToNumpyArray

    ConvertToNumpyArray(bit_vector, array)
    return array


def deduplicate_labeled_records(records: list[dict]) -> tuple[list[dict], dict]:
    """Collapse exact structures and discard conflicting weak labels."""
    by_structure: dict[str, list[dict]] = defaultdict(list)
    invalid = 0
    for record in records:
        canonical = canonical_smiles(record.get("smiles", ""))
        if not canonical:
            invalid += 1
            continue
        normalized = dict(record)
        normalized["smiles"] = canonical
        by_structure[canonical].append(normalized)

    output: list[dict] = []
    conflicts = 0
    duplicate_rows = 0
    for canonical, group in by_structure.items():
        labels = {str(row.get("label") or "") for row in group}
        duplicate_rows += max(0, len(group) - 1)
        if len(labels) != 1:
            conflicts += 1
            continue
        representative = dict(group[0])
        representative["smiles"] = canonical
        representative["duplicate_count"] = len(group)
        output.append(representative)
    output.sort(key=lambda row: (row["smiles"], str(row.get("label") or "")))
    return output, {
        "input_rows": len(records),
        "unique_structures": len(output),
        "invalid_smiles_rows": invalid,
        "duplicate_rows_collapsed": duplicate_rows,
        "conflicting_structure_labels_discarded": conflicts,
    }


def scaffold_split_indices(
    labels: np.ndarray,
    groups: np.ndarray,
    *,
    random_state: int = RANDOM_STATE,
) -> dict[str, np.ndarray]:
    """Create mutually scaffold-disjoint train/calibration/test indices."""
    unique_groups = np.unique(groups)
    if len(unique_groups) < 5:
        raise ValueError("Need at least five distinct scaffolds for a 60/20/20 split")
    splitter = StratifiedGroupKFold(
        n_splits=5, shuffle=True, random_state=random_state
    )
    folds = [
        test_index
        for _train_index, test_index in splitter.split(
            np.zeros(len(labels)), labels, groups
        )
    ]
    all_classes = set(np.unique(labels).tolist())
    # Choose two distinct folds for untouched test/calibration while requiring
    # every class in all three partitions.  This avoids silently evaluating a
    # class only where it happened to occur.
    for test_fold in range(len(folds)):
        for calibration_fold in range(len(folds)):
            if calibration_fold == test_fold:
                continue
            test = folds[test_fold]
            calibration = folds[calibration_fold]
            reserved = np.concatenate([test, calibration])
            train = np.setdiff1d(np.arange(len(labels)), reserved)
            if all(
                set(np.unique(labels[index]).tolist()) == all_classes
                for index in (train, calibration, test)
            ):
                return {
                    "train": train,
                    "calibration": calibration,
                    "test": test,
                }
    raise ValueError(
        "Could not create scaffold-disjoint partitions containing every class; "
        "add independently labelled scaffolds or drop under-supported classes"
    )


def _softmax_from_probabilities(probabilities: np.ndarray, temperature: float) -> np.ndarray:
    clipped = np.clip(np.asarray(probabilities, dtype=float), 1e-9, 1.0)
    logits = np.log(clipped) / max(float(temperature), 1e-6)
    logits -= logits.max(axis=1, keepdims=True)
    exp = np.exp(logits)
    return exp / exp.sum(axis=1, keepdims=True)


def fit_temperature(probabilities: np.ndarray, labels: np.ndarray) -> float:
    candidates = np.geomspace(0.35, 5.0, 121)
    losses = [
        log_loss(
            labels,
            _softmax_from_probabilities(probabilities, temperature),
            labels=np.arange(probabilities.shape[1]),
        )
        for temperature in candidates
    ]
    return float(candidates[int(np.argmin(losses))])


def wilson_lower(successes: int, total: int, z: float = 1.959963984540054) -> float:
    if total <= 0:
        return 0.0
    proportion = successes / total
    denominator = 1.0 + z * z / total
    centre = proportion + z * z / (2.0 * total)
    radius = z * math.sqrt(
        (proportion * (1.0 - proportion) + z * z / (4.0 * total)) / total
    )
    return max(0.0, (centre - radius) / denominator)


def fit_class_thresholds(
    probabilities: np.ndarray,
    labels: np.ndarray,
    *,
    target_precision: float = TARGET_SELECTIVE_PRECISION,
    min_calls: int = MIN_CALIBRATION_CALLS,
) -> tuple[dict[int, float], dict[int, dict]]:
    predicted = probabilities.argmax(axis=1)
    confidence = probabilities.max(axis=1)
    thresholds: dict[int, float] = {}
    audit: dict[int, dict] = {}
    for class_index in range(probabilities.shape[1]):
        candidate_thresholds = sorted({
            float(value)
            for value in confidence[predicted == class_index]
        })
        selected = 1.01  # explicit abstention if no calibrated operating point
        selected_audit = {
            "calls": 0,
            "correct": 0,
            "precision": None,
            "wilson_lower": 0.0,
            "qualified": False,
        }
        for threshold in candidate_thresholds:
            mask = (predicted == class_index) & (confidence >= threshold)
            calls = int(mask.sum())
            if calls < min_calls:
                continue
            correct = int((labels[mask] == class_index).sum())
            lower = wilson_lower(correct, calls)
            if lower >= target_precision:
                selected = threshold
                selected_audit = {
                    "calls": calls,
                    "correct": correct,
                    "precision": correct / calls,
                    "wilson_lower": lower,
                    "qualified": True,
                }
                break
        thresholds[class_index] = float(selected)
        audit[class_index] = selected_audit
    return thresholds, audit


def selective_predictions(
    probabilities: np.ndarray,
    thresholds: dict[int, float],
) -> np.ndarray:
    predicted = probabilities.argmax(axis=1)
    confidence = probabilities.max(axis=1)
    return np.asarray([
        int(class_index)
        if confidence[row_index] >= thresholds.get(int(class_index), 1.01)
        else -1
        for row_index, class_index in enumerate(predicted)
    ], dtype=int)


def _fit_with_scaffold_early_stopping(
    X: np.ndarray,
    y: np.ndarray,
    groups: np.ndarray,
    class_count: int,
) -> tuple[MLPClassifier, int]:
    splitter = GroupShuffleSplit(
        n_splits=12, test_size=0.15, random_state=RANDOM_STATE
    )
    all_classes = set(range(class_count))
    inner_train = inner_validation = None
    for train_index, validation_index in splitter.split(X, y, groups):
        if set(np.unique(y[train_index]).tolist()) == all_classes:
            inner_train, inner_validation = train_index, validation_index
            break
    if inner_train is None or inner_validation is None:
        raise ValueError("Unable to create a class-complete scaffold early-stop split")

    classifier = MLPClassifier(
        hidden_layer_sizes=(512,),
        batch_size=min(256, max(16, len(inner_train))),
        alpha=1e-4,
        random_state=RANDOM_STATE,
    )
    best_model = None
    best_loss = float("inf")
    best_epoch = 0
    stale = 0
    for epoch in range(1, 201):
        classifier.partial_fit(
            X[inner_train], y[inner_train], classes=np.arange(class_count)
        )
        validation_probability = classifier.predict_proba(X[inner_validation])
        current_loss = log_loss(
            y[inner_validation], validation_probability,
            labels=np.arange(class_count),
        )
        if current_loss < best_loss - 1e-4:
            best_loss = current_loss
            best_model = copy.deepcopy(classifier)
            best_epoch = epoch
            stale = 0
        else:
            stale += 1
        if stale >= 12:
            break
    if best_model is None:
        raise RuntimeError("MLP failed to produce an early-stopped checkpoint")
    return best_model, best_epoch


def load_weak_db_labels(db_path: str) -> tuple[list[dict], dict]:
    connection = sqlite3.connect(db_path)
    rows = connection.execute(
        "SELECT id, smiles, classification FROM compounds "
        "WHERE smiles IS NOT NULL AND smiles != ''"
    ).fetchall()
    connection.close()
    records = []
    source_counts = Counter()
    for compound_id, smiles, existing_class in rows:
        detailed = RC.predict_detailed(smiles)
        if detailed["class"] != "Unknown":
            label = detailed["class"]
            source = "weak_smarts_generated"
        else:
            label = normalize_classification(existing_class or "")
            if not label or label == "Unknown":
                continue
            source = "weak_existing_database_label"
        source_counts[source] += 1
        records.append({
            "id": compound_id,
            "smiles": smiles,
            "label": label,
            "label_source": source,
        })
    return records, {
        "kind": "weak_structure_derived_and_existing_database_labels",
        "independent_ground_truth": False,
        "source_counts": dict(source_counts),
    }


def load_independent_csv(path: Path) -> tuple[list[dict], dict]:
    records = []
    with path.open(newline="", encoding="utf-8-sig") as handle:
        for row_number, row in enumerate(csv.DictReader(handle), start=2):
            label = normalize_classification(
                row.get("classification") or row.get("label") or ""
            )
            smiles = row.get("smiles") or ""
            if label and label != "Unknown" and smiles:
                records.append({
                    "id": row.get("id") or f"csv:{row_number}",
                    "smiles": smiles,
                    "label": label,
                    "label_source": row.get("label_source") or "external_csv",
                })
    return records, {
        "kind": "user_declared_independent_csv",
        "independent_ground_truth": True,
        "path": str(path.resolve()),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--full", action="store_true", help="use full DB")
    parser.add_argument(
        "--labels-csv",
        type=Path,
        help="independently curated CSV with smiles and classification columns",
    )
    parser.add_argument(
        "--labels-are-independent",
        action="store_true",
        help="required acknowledgement for --labels-csv",
    )
    parser.add_argument("--output", type=Path, default=MODEL_PATH)
    args = parser.parse_args()
    if args.labels_csv and not args.labels_are_independent:
        parser.error("--labels-csv requires --labels-are-independent")

    started = time.time()
    if args.labels_csv:
        raw_records, label_provenance = load_independent_csv(args.labels_csv)
    else:
        database = ROOT / "data" / (
            "local_library.db" if args.full else "local_library_small.db"
        )
        raw_records, label_provenance = load_weak_db_labels(str(database))

    records, dedup_audit = deduplicate_labeled_records(raw_records)
    sample_counts = Counter(row["label"] for row in records)
    scaffold_counts: dict[str, set[str]] = defaultdict(set)
    for row in records:
        scaffold_counts[row["label"]].add(scaffold_key(row["smiles"]))
    keep = {
        label for label, count in sample_counts.items()
        if count >= MIN_CLASS_SAMPLES
        and len(scaffold_counts[label]) >= MIN_CLASS_SCAFFOLDS
    }
    records = [row for row in records if row["label"] in keep]
    classes = sorted(keep)
    if len(classes) < 2:
        raise ValueError("Need at least two sufficiently supported classes")
    class_to_index = {label: index for index, label in enumerate(classes)}

    fingerprints = []
    labels = []
    groups = []
    retained = []
    for row in records:
        fingerprint = smiles_to_fp(row["smiles"])
        if fingerprint is None:
            continue
        fingerprints.append(fingerprint)
        labels.append(class_to_index[row["label"]])
        groups.append(scaffold_key(row["smiles"]))
        retained.append(row)
    X = np.asarray(fingerprints, dtype=np.uint8)
    y = np.asarray(labels, dtype=int)
    group_array = np.asarray(groups, dtype=object)
    splits = scaffold_split_indices(y, group_array)

    classifier, best_epoch = _fit_with_scaffold_early_stopping(
        X[splits["train"]],
        y[splits["train"]],
        group_array[splits["train"]],
        len(classes),
    )
    calibration_raw = classifier.predict_proba(X[splits["calibration"]])
    temperature = fit_temperature(calibration_raw, y[splits["calibration"]])
    calibration_probability = _softmax_from_probabilities(
        calibration_raw, temperature
    )
    thresholds, threshold_audit = fit_class_thresholds(
        calibration_probability, y[splits["calibration"]]
    )

    test_probability = _softmax_from_probabilities(
        classifier.predict_proba(X[splits["test"]]), temperature
    )
    test_selective = selective_predictions(test_probability, thresholds)
    accepted = test_selective >= 0
    test_coverage = float(accepted.mean()) if len(accepted) else 0.0
    selective_accuracy = (
        float((test_selective[accepted] == y[splits["test"]][accepted]).mean())
        if accepted.any() else None
    )
    raw_prediction = test_probability.argmax(axis=1)
    weak_macro_f1 = float(
        f1_score(y[splits["test"]], raw_prediction, average="macro")
    )
    report = classification_report(
        y[splits["test"]],
        raw_prediction,
        labels=np.arange(len(classes)),
        target_names=classes,
        output_dict=True,
        zero_division=0,
    )

    split_groups = {
        name: set(group_array[index].tolist()) for name, index in splits.items()
    }
    split_audit = {
        "policy": "canonical_smiles_dedup_then_bemis_murcko_group_split",
        "train_n": int(len(splits["train"])),
        "calibration_n": int(len(splits["calibration"])),
        "test_n": int(len(splits["test"])),
        "train_calibration_scaffold_overlap": len(
            split_groups["train"] & split_groups["calibration"]
        ),
        "train_test_scaffold_overlap": len(
            split_groups["train"] & split_groups["test"]
        ),
        "calibration_test_scaffold_overlap": len(
            split_groups["calibration"] & split_groups["test"]
        ),
    }
    production_eligible = False
    bundle = {
        "schema_version": 2,
        "model": classifier,
        "classes": classes,
        "fp_bits": FP_BITS,
        "fp_radius": FP_RADIUS,
        "temperature": temperature,
        "class_thresholds": {
            classes[index]: threshold for index, threshold in thresholds.items()
        },
        "threshold_audit": {
            classes[index]: values for index, values in threshold_audit.items()
        },
        "target_selective_precision": TARGET_SELECTIVE_PRECISION,
        "minimum_calibration_calls": MIN_CALIBRATION_CALLS,
        "label_provenance": label_provenance,
        "dedup_audit": dedup_audit,
        "split_audit": split_audit,
        "early_stopping": {
            "policy": "inner_scaffold_holdout",
            "best_epoch": best_epoch,
        },
        "test_evaluation": {
            "metric_name": (
                "external_label_macro_f1"
                if label_provenance["independent_ground_truth"]
                else "weak_label_agreement_macro_f1"
            ),
            "macro_f1": weak_macro_f1,
            "selective_coverage": test_coverage,
            "selective_accuracy": selective_accuracy,
            "classification_report": report,
        },
        "production_eligible": production_eligible,
        "production_eligibility_reason": (
            "requires_governance_promotion_after_independent_untouched_test_review"
            if label_provenance["independent_ground_truth"]
            else "training_and_evaluation_labels_are_not_independent_ground_truth"
        ),
        "default_prediction_policy": "abstain_unless_per_class_threshold_passes",
        "created_unix_time": time.time(),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    joblib.dump(bundle, args.output)
    audit_path = args.output.with_suffix(args.output.suffix + ".audit.json")
    audit_path.write_text(
        json.dumps(
            {key: value for key, value in bundle.items() if key != "model"},
            ensure_ascii=False,
            indent=2,
        ),
        encoding="utf-8",
    )
    print(
        f"Saved research-only model to {args.output}; "
        f"weak/external test macro-F1={weak_macro_f1:.3f}, "
        f"selective coverage={test_coverage:.3f}, elapsed={time.time()-started:.0f}s"
    )


if __name__ == "__main__":
    main()
