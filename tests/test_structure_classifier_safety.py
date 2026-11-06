"""Leakage, calibration, and deployment guards for the structure classifier."""

import numpy as np
import pytest

from classifier.auto_classify_db import (
    selective_classifications,
    validate_model_bundle,
)
from classifier.train_structure_clf import (
    deduplicate_labeled_records,
    fit_class_thresholds,
    scaffold_split_indices,
    selective_predictions,
)


def test_exact_structure_duplicates_and_conflicts_are_removed_before_split():
    records = [
        {"id": 1, "smiles": "CCO", "label": "alcohol"},
        {"id": 2, "smiles": "OCC", "label": "alcohol"},
        {"id": 3, "smiles": "c1ccccc1", "label": "aromatic"},
        {"id": 4, "smiles": "C1=CC=CC=C1", "label": "conflict"},
    ]

    unique, audit = deduplicate_labeled_records(records)

    assert [row["label"] for row in unique] == ["alcohol"]
    assert audit["duplicate_rows_collapsed"] == 2
    assert audit["conflicting_structure_labels_discarded"] == 1


def test_train_calibration_test_scaffolds_are_disjoint():
    labels = np.asarray([0, 1] * 10, dtype=int)
    groups = np.asarray([f"scaffold-{index}" for index in range(20)])

    splits = scaffold_split_indices(labels, groups, random_state=7)

    train = set(groups[splits["train"]])
    calibration = set(groups[splits["calibration"]])
    test = set(groups[splits["test"]])
    assert not train & calibration
    assert not train & test
    assert not calibration & test


def test_class_without_qualified_calibration_operating_point_abstains():
    probabilities = np.asarray([
        [0.60, 0.40],
        [0.61, 0.39],
        [0.62, 0.38],
        [0.63, 0.37],
    ])
    labels = np.asarray([0, 1, 1, 1])

    thresholds, audit = fit_class_thresholds(
        probabilities, labels, target_precision=0.90, min_calls=3
    )
    predictions = selective_predictions(probabilities, thresholds)

    assert thresholds[0] > 1.0
    assert audit[0]["qualified"] is False
    assert np.all(predictions == -1)


def test_runtime_refuses_unvalidated_bundle_by_default():
    bundle = {
        "schema_version": 2,
        "model": object(),
        "classes": ["a", "b"],
        "fp_bits": 8,
        "fp_radius": 2,
        "temperature": 1.0,
        "class_thresholds": {"a": 0.8, "b": 0.9},
        "split_audit": {
            "train_calibration_scaffold_overlap": 0,
            "train_test_scaffold_overlap": 0,
            "calibration_test_scaffold_overlap": 0,
        },
        "label_provenance": {"independent_ground_truth": False},
        "production_eligible": False,
    }

    with pytest.raises(PermissionError):
        validate_model_bundle(bundle)
    assert validate_model_bundle(bundle, allow_research=True) == (
        "ml_research_unvalidated"
    )


def test_runtime_uses_per_class_thresholds_instead_of_global_half():
    decisions = selective_classifications(
        np.asarray([[0.70, 0.30], [0.20, 0.80]]),
        ["a", "b"],
        {"a": 0.75, "b": 0.75},
    )

    assert decisions[0] == ("other", 0.70, False)
    assert decisions[1] == ("b", 0.80, True)
