from __future__ import annotations

import json

from eval.evograph.phase0_audit import (
    audit_phase0_manifests,
    load_manifest,
    write_audit_report,
)


def _row(key: str, scaffold: str, spectrum_hash: str, label: str, mz: float) -> dict:
    return {
        "inchikey": key,
        "canonical_smiles": f"SMILES-{key[:4]}",
        "formula": "C10H10O2",
        "murcko_scaffold": scaffold,
        "spectrum_hash": spectrum_hash,
        "spectrum_id": spectrum_hash,
        "label": label,
        "precursor_mz": mz,
        "source": "external-standards",
        "instrument": "QTOF",
    }


def test_phase0_clean_identity_and_scaffold_split_passes_with_negative_controls():
    train = [
        _row("AAAAAAAAAAAAAA-UHFFFAOYSA-N", "scaffold-a", "hash-a", "A", 101.0),
        _row("BBBBBBBBBBBBBB-UHFFFAOYSA-N", "scaffold-b", "hash-b", "B", 201.0),
    ]
    test = [
        _row("CCCCCCCCCCCCCC-UHFFFAOYSA-N", "scaffold-c", "hash-c", "A", 102.0),
        _row("DDDDDDDDDDDDDD-UHFFFAOYSA-N", "scaffold-d", "hash-d", "B", 202.0),
    ]

    report = audit_phase0_manifests(train, test)

    assert report["gate"]["passed"] is True
    assert report["gate"]["blocking_violations"] == []
    assert report["split_overlap"]["molecular_formula_diagnostic_only"]["overlap_unique"] == 1
    assert report["gate"]["formula_overlap_is_blocking"] is False
    controls = {item["name"]: item for item in report["negative_controls"]}
    assert controls["candidate_frequency_majority"]["status"] == "available"
    assert controls["mass_only_nearest_train"]["status"] == "available"
    assert all(item["uses_msms_peaks"] is False for item in controls.values())


def test_phase0_identity_spectrum_and_scaffold_overlap_blocks_evaluation():
    shared = _row(
        "AAAAAAAAAAAAAA-UHFFFAOYSA-N",
        "shared-scaffold",
        "shared-hash",
        "A",
        101.0,
    )

    report = audit_phase0_manifests([shared], [dict(shared)])

    assert report["gate"]["passed"] is False
    blockers = set(report["gate"]["blocking_violations"])
    assert "full_inchikey_overlap" in blockers
    assert "connectivity_inchikey_overlap" in blockers
    assert "spectrum_hash_overlap" in blockers
    assert "scaffold_overlap" in blockers


def test_phase0_manifest_loader_and_atomic_report(tmp_path):
    manifest = tmp_path / "rows.jsonl"
    manifest.write_text(
        json.dumps(_row("AAAAAAAAAAAAAA-UHFFFAOYSA-N", "s-a", "h-a", "A", 101.0))
        + "\n",
        encoding="utf-8",
    )
    rows = load_manifest(manifest)
    output = tmp_path / "audit.json"
    test_rows = [
        _row("BBBBBBBBBBBBBB-UHFFFAOYSA-N", "s-b", "h-b", "B", 201.0)
    ]
    report = audit_phase0_manifests(rows, test_rows, train_path=manifest)

    write_audit_report(report, output)

    assert len(rows) == 1
    assert json.loads(output.read_text(encoding="utf-8"))["train"]["sha256"]


def test_phase0_missing_split_provenance_cannot_pass_silently():
    report = audit_phase0_manifests(
        [{"label": "A", "precursor_mz": 100.0}],
        [{"label": "B", "precursor_mz": 200.0}],
    )

    assert report["gate"]["passed"] is False
    blockers = set(report["gate"]["blocking_violations"])
    assert "train_identity_coverage_incomplete" in blockers
    assert "test_spectrum_hash_coverage_incomplete" in blockers
