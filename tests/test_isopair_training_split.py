"""Regression tests for IsoPair formula/connectivity split leakage."""

from scripts.train_isopair_ms import (
    SpectrumRow,
    enforce_formula_and_connectivity_disjoint_splits,
    formula_split,
    split_formula_key,
    split_overlap_audit,
)
from scripts.train_isopair_metric_ranker import (
    group_early_stopping_indices,
    split_validation_rows,
)
import numpy as np


def _row(index, formula, connectivity):
    return SpectrumRow(
        spectrum_id=index,
        accession=f"s{index}",
        formula=formula,
        connectivity=connectivity,
        adduct="[m+h]+",
        precursor_mz=300.0,
        collision_energy=None,
        collision_energy_unit="",
        instrument_type="qtof",
        dataset_source="source",
        peaks=((50.0, 1.0), (100.0, 0.5), (150.0, 0.2)),
    )


def test_bracketed_charged_formula_hashes_with_plain_formula():
    assert split_formula_key("[C27H31O16]+") == "C27H31O16"
    assert split_formula_key("C27H31O16") == "C27H31O16"
    assert formula_split("[C27H31O16]+", 3407) == formula_split(
        "C27H31O16", 3407
    )


def test_connectivity_bridge_forces_transitive_component_into_one_split():
    spectra = [
        _row(1, "[C27H31O16]+", "CONNECTIVITYAAA"),
        _row(2, "C27H31O16", "CONNECTIVITYBBB"),
        # Same connectivity under a second formula spelling/family bridges the
        # full component and must never land in another partition.
        _row(3, "C28H34O16", "CONNECTIVITYBBB"),
        _row(4, "C28H34O16", "CONNECTIVITYCCC"),
    ]
    deliberately_leaky = ["train", "validation", "test", "train"]

    repaired = enforce_formula_and_connectivity_disjoint_splits(
        spectra, deliberately_leaky, seed=3407
    )
    audit = split_overlap_audit(spectra, repaired)

    assert len(set(repaired)) == 1
    assert audit["zero_overlap"] is True
    assert all(
        values["formula_overlap"] == 0
        and values["connectivity_overlap"] == 0
        for values in audit["pairwise"].values()
    )


def test_metric_early_stopping_keeps_complete_families_together():
    groups = np.asarray([
        "f1", "f1", "f2", "f2", "f3", "f3", "f4", "f4",
        "f5", "f5", "f6", "f6", "f7", "f7", "f8", "f8",
    ])
    labels = np.asarray([0, 1] * 8)

    train, validation = group_early_stopping_indices(
        groups, labels, seed=3407, validation_fraction=0.25
    )

    assert not set(groups[train]) & set(groups[validation])
    assert set(labels[train]) == {0, 1}
    assert set(labels[validation]) == {0, 1}


def test_fusion_and_threshold_tuning_use_different_formula_families():
    rows = [
        {"formula": f"C{index}H{index * 2}", "candidates": []}
        for index in range(1, 13)
    ]

    fusion, threshold = split_validation_rows(rows, seed=3407)

    assert fusion and threshold
    assert not {row["formula"] for row in fusion} & {
        row["formula"] for row in threshold
    }
