from __future__ import annotations

from copy import deepcopy

import pytest

from matching.reference_consensus import (
    CONSENSUS_SCHEMA_VERSION,
    replace_technical_replicates_with_consensus,
)


def _record(index: int, **overrides):
    record = {
        "spectrum_id": f"replicate-{index}",
        "source_file": f"replicate-{index}.mgf",
        "source_row_index": index,
        "source_asset_sha256": f"{index + 1:064x}",
        "name": "Ethanol reference",
        "formula": "C2H6O",
        "smiles": "CCO",
        "precursor_mz": 47.049141,
        "adduct": "[M+H]+",
        "ion_mode": "POS",
        "charge": 1,
        "instrument_type": "qtof",
        "instrument": "Agilent 6560 QTOF",
        "activation_method": "CID",
        "collision_cell_region": "collision cell",
        "instrument_domain_id": "agilent-qtof-cid",
        "acquisition_method": "method-20260826",
        "collision_energy": 20.0,
        "collision_energy_min": 20.0,
        "collision_energy_max": 20.0,
        "collision_energy_unit": "eV",
        "collision_energy_mode": "fixed",
        "reference_material_type": "authentic_standard",
        "standard_lot": "LOT-001",
        "standard_purity": 99.5,
        "source_provenance_verified": True,
        "source_tier": "experimental",
        "intensity_origin": "measured",
        "peaks": [
            [20.000, 10.0],
            [30.000, 100.0],
            [40.000, 50.0],
        ],
        "source_provenance": {"manifest_id": "unit-test"},
    }
    record.update(overrides)
    return record


def test_consensus_replaces_same_condition_replicates_without_mutating_inputs():
    first = _record(0, peaks=[
        [20.000, 10.0], [30.000, 100.0], [40.000, 50.0], [40.004, 20.0]
    ])
    second = _record(1, peaks=[
        [20.005, 20.0], [30.004, 80.0], [40.006, 40.0], [45.000, 10.0]
    ])
    original = deepcopy([first, second])

    output, summary = replace_technical_replicates_with_consensus([first, second])

    assert [first, second] == original
    assert summary["schema_version"] == CONSENSUS_SCHEMA_VERSION
    assert summary["consensus_groups_built"] == 1
    assert summary["member_records_replaced"] == 2
    assert len(output) == 1
    consensus = output[0]
    assert consensus["source_tier"] == "experimental"
    assert consensus["intensity_origin"] == "measured_consensus"
    assert consensus["consensus_member_count"] == 2
    assert consensus["technical_replicates_count_as_independent_evidence"] is False
    assert consensus["independent_evidence_count"] == 1
    assert [peak[0] for peak in consensus["peaks"]] == pytest.approx(
        [20.0025, 30.002, 40.003]
    )
    assert max(peak[1] for peak in consensus["peaks"]) == pytest.approx(100.0)
    provenance = consensus["source_provenance"]
    assert provenance["consensus_member_count"] == 2
    assert len(provenance["members"]) == 2
    assert all(member["source_asset_sha256"] for member in provenance["members"])


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("collision_energy", 30.0),
        ("collision_energy_unit", "NCE"),
        ("instrument", "Thermo Q Exactive"),
        ("adduct", "[M+Na]+"),
        ("standard_lot", "LOT-002"),
    ],
)
def test_cross_condition_records_are_never_combined(field, value):
    first = _record(0)
    second = _record(1, **{field: value})

    output, summary = replace_technical_replicates_with_consensus([first, second])

    assert len(output) == 2
    assert summary["consensus_groups_built"] == 0


def test_digitized_images_and_unverified_inchikeys_cannot_build_consensus():
    image = _record(
        0,
        intensity_origin="digitized_image",
        source_tier="user_unverified",
    )
    no_structure_a = _record(1, smiles="", reported_inchikey="LFQSCWFLJHTTHZ-UHFFFAOYSA-N")
    no_structure_b = _record(2, smiles="", reported_inchikey="LFQSCWFLJHTTHZ-UHFFFAOYSA-N")

    output, summary = replace_technical_replicates_with_consensus(
        [image, no_structure_a, no_structure_b]
    )

    assert len(output) == 3
    assert summary["eligible_records"] == 0
    assert summary["consensus_groups_built"] == 0


def test_same_uploaded_scan_is_not_counted_twice_as_two_replicates():
    first = _record(0)
    duplicate = deepcopy(first)

    output, summary = replace_technical_replicates_with_consensus([first, duplicate])

    assert len(output) == 2
    assert summary["consensus_groups_built"] == 0
    assert summary["technical_replicates_count_as_independent_evidence"] is False
