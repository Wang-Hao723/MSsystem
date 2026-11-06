import copy
import json
import random

import pytest

from eval.v16_evidence_validation import evaluate as audit
from eval.v16_evidence_validation import ms1_feature_challenger as challenger
from matching.v16_dual_polarity import semantic_sha256


@pytest.mark.parametrize("field", ["input", "output_features", "feature_count", "status"])
def test_ms1_manifest_must_match_exact_task_source(tmp_path, field):
    source = tmp_path / "source.mzML"
    xml = tmp_path / "features.featureXML"
    manifest = {"status": "ok", "feature_count": 0, "features": [],
        "input": str(source), "output_features": str(xml)}
    challenger.check_manifest_binding(manifest, {"feature_xml_path": str(xml)}, source)
    manifest[field] = 1 if field == "feature_count" else "wrong"
    with pytest.raises(ValueError, match="does not bind"):
        challenger.check_manifest_binding(manifest, {"feature_xml_path": str(xml)}, source)


def test_trace_control_is_seeded_preserves_masses_rt_values_and_inputs():
    traces = {"f": {"rt": [1, 2, 3, 4, 5, 6], "target_mz": 100.,
        "intensity": [0, 1, 2, 4, 2, 0]}}
    before = copy.deepcopy(traces)
    first = challenger.permute_trace_intensities(traces, random.Random(1))
    second = challenger.permute_trace_intensities(traces, random.Random(1))
    assert first == second
    assert traces == before
    assert first["f"]["rt"] == traces["f"]["rt"]
    assert first["f"]["target_mz"] == 100.
    assert sorted(first["f"]["intensity"]) == sorted(traces["f"]["intensity"])
    assert first["f"]["intensity"] != traces["f"]["intensity"]


@pytest.mark.parametrize("field", ["production_mutation_allowed", "automatic_promotion_allowed",
    "new_source_discovery_allowed", "labels_used_for_graph_construction", "rt_alignment_added"])
def test_protocol_cannot_enable_production_or_label_driven_graphs(tmp_path, monkeypatch, field):
    payload = json.loads(challenger.PROTOCOL.read_text())
    payload[field] = True
    path = tmp_path / "protocol.json"
    path.write_text(json.dumps(payload))
    monkeypatch.setattr(challenger, "PROTOCOL", path)
    with pytest.raises(ValueError, match="invariant changed"):
        challenger.load_protocol()


def test_bad_saved_graph_cannot_open_labels(tmp_path, monkeypatch):
    path = tmp_path / "controls.json"
    path.write_text(json.dumps({"schema_version": challenger.SCHEMA, "phase": "score"}))
    monkeypatch.setattr(audit, "bind_retrospective_labels",
        lambda *args: pytest.fail("Labels opened before controls were validated"))
    with pytest.raises(ValueError, match="binding mismatch"):
        challenger.score({}, {}, path)


@pytest.mark.parametrize("field", ["implementation_file_sha256", "evaluator_file_sha256",
    "protocol_semantic_sha256", "labels_opened", "variants"])
def test_valid_semantic_hash_does_not_override_code_or_phase_binding(field):
    protocol = {"protocol_id": "test"}
    frozen = {**challenger.base_receipt(protocol, "controls"),
        "labels_opened": False, "variants": {name: {} for name in challenger.VARIANTS}}
    frozen["semantic_sha256"] = audit.receipt_semantic_sha256(frozen)
    challenger.validate_controls_receipt(frozen, protocol)
    frozen[field] = True if field == "labels_opened" else {} if field == "variants" else "wrong"
    frozen["semantic_sha256"] = audit.receipt_semantic_sha256(frozen)
    with pytest.raises(ValueError, match="binding mismatch"):
        challenger.validate_controls_receipt(frozen, protocol)


def test_protocol_keeps_both_ablation_arms_and_all_original_gates():
    protocol, _, _ = challenger.load_protocol()
    assert protocol["variants"] == ["linked_ms1_coordinates", "all_ms1_partners"]
    assert protocol["locked_thresholds"]["minimum_shape_correlation"] == .85
    assert protocol["locked_thresholds"]["mass_ppm"] == 10.
    assert protocol["locked_thresholds"]["cross_polarity_rt_delta_seconds"] == 20.
    assert semantic_sha256(protocol)


def reconstruction_fixture(tmp_path, monkeypatch):
    result_path, source, xml, manifest_path = [tmp_path / name for name in
        ("results.json", "source.mzML", "features.featureXML", "features.json")]
    rows = [{"precursor_mz": 299., "retention_time": 1., "ion_mode": "NEG",
        "lcms_feature_id": "F1", "lcms_feature_mz": 299., "lcms_feature_rt_minutes": 1.,
        "lcms_feature_link_status": "linked_openms_feature", "match": {"name": "POISON"}}]
    manifest = {"status": "ok", "feature_count": 1, "features": [
        {"feature_id": "F1", "mz": 299., "rt_minutes": 1., "charge": -1}],
        "input": str(source), "output_features": str(xml)}
    result_path.write_text(json.dumps(rows))
    manifest_path.write_text(json.dumps(manifest))
    source.write_text("test mzML placeholder; extraction is injected")
    xml.write_text("test featureXML placeholder; hash verified")
    bound = {"results_path": str(result_path), "results_file_sha256": audit.file_sha256(result_path),
        "feature_manifest_path": str(manifest_path), "feature_manifest_sha256": audit.file_sha256(manifest_path),
        "feature_xml_path": str(xml), "feature_xml_sha256": audit.file_sha256(xml)}
    pair = {"task_pair": {"NEG": "NEG"}, "source_file_bindings": {"NEG": bound}}
    record = {"task_id": "NEG", "polarity": "NEG", "candidate_independent_result_rows": 1,
        "source_binding": {"results_file": {"path": str(result_path), "sha256": audit.file_sha256(result_path)},
            "source_data_file": {"path": str(source), "sha256": audit.file_sha256(source)},
            "candidate_independent_projection_sha256": semantic_sha256(challenger.project_candidate_independent_rows(rows))}}
    snapshots = {str(path): audit.file_sha256(path) for path in (result_path, source)}
    monkeypatch.setattr(challenger, "extract_targeted_ms1_traces",
        lambda *args: ({"F1": {"rt": [1.], "intensity": [2.]}}, {"status": "ok"}))
    return pair, record, snapshots, xml


def test_feature_hashes_can_be_bound_in_pair_without_prior_flat_hash_copy(tmp_path, monkeypatch):
    pair, record, snapshots, xml = reconstruction_fixture(tmp_path, monkeypatch)
    assert str(xml) not in snapshots
    prepared, traces, receipt = challenger.reconstruct_input(pair, record, snapshots)
    assert snapshots[str(xml)] == audit.file_sha256(xml)
    assert len(snapshots) == 4
    assert prepared["row_bindings"][0]["source_row_index"] == 0
    assert "POISON" not in json.dumps((prepared, receipt))


def test_changed_feature_xml_is_rejected_before_trace_extraction(tmp_path, monkeypatch):
    pair, record, snapshots, xml = reconstruction_fixture(tmp_path, monkeypatch)
    xml.write_text("changed")
    monkeypatch.setattr(challenger, "extract_targeted_ms1_traces",
        lambda *args: pytest.fail("unbound MS1 input was read"))
    with pytest.raises(RuntimeError, match="changed"):
        challenger.reconstruct_input(pair, record, snapshots)
