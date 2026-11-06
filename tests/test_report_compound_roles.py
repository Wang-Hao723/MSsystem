from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
import shutil
import subprocess

import pytest
from rdkit import Chem

from pipeline import compound_role_projection as roles
from pipeline.compound_conclusions import (
    build_compound_conclusions, filter_compound_conclusions,
    read_compound_conclusion_artifact, write_compound_conclusion_artifact,
    artifact_is_current,
)


KEY = "FKMHSNTVILORFA-UHFFFAOYSA-N"
SMILES = "CCCCCCCCCCCCOCCOCCOCCO"


def _candidate(smiles=SMILES, name="Triethylene Glycol Monododecyl Ether"):
    return {
        "compound_name": name,
        "inchikey": Chem.MolToInchiKey(Chem.MolFromSmiles(smiles)),
        "smiles": smiles, "formula": "C18H38O4", "has_structure": True,
        "confidence": "High", "identification_scope": "2d_connectivity",
        "cosine": 0.99, "final_score": 0.98, "rank_score": 0.98,
        "mass_ppm": 1.0, "compound_role": "unclassified",
        "potential_contaminant": False, "lib_spectrum_id": "reference-original",
        "dataset_source": "synthetic-test",
    }


def _row(match):
    return {
        "spectrum_id": "original-spectrum-id", "match": match,
        "precursor_mz": 319.2843, "retention_time": 3.1,
        "retention_time_unit": "minute", "ion_mode": "POS",
        "blank_used": True, "blank_evidence_category": "sample_enriched",
        "raw_peaks": [[70.0, 10.0], [120.0, 20.0]],
    }


def test_exact_structure_role_does_not_change_identity_scores_or_high():
    match = _candidate()
    original = deepcopy(match)
    roles.apply_report_compound_roles(match)
    assert match["compound_role"] == "external_contaminant"
    assert match["potential_contaminant"] is True
    assert roles.plant_report_selection(match)["excluded_from_plant_candidate_table"] is True
    for field in ("compound_name", "smiles", "inchikey", "formula", "confidence",
                  "cosine", "final_score", "rank_score", "mass_ppm", "lib_spectrum_id"):
        assert match[field] == original[field]
    assert match["report_role_projection"]["detection_confirmed_by_role"] is False
    assert "不证明样本中实际存在" in match["contaminant_note"]


def test_similar_alkyl_glycerol_structure_and_alias_only_do_not_trigger_rule():
    analogue = _candidate("CCCCCCCCCCCCCCCCCCOCC(O)CO", "Alkyl glycerol")
    roles.apply_report_compound_roles(analogue)
    assert analogue["potential_contaminant"] is False
    assert roles.plant_report_selection(analogue)["excluded_from_plant_candidate_table"] is False
    alias_only = {"compound_name": "Laureth-3"}
    roles.apply_report_compound_roles(alias_only)
    assert alias_only["potential_contaminant"] is False
    assert "report_role_projection" not in alias_only


def test_missing_key_can_bind_exact_structure_but_conflicting_key_cannot():
    match = _candidate()
    match.pop("inchikey")
    roles.apply_report_compound_roles(match)
    assert match["potential_contaminant"] is True
    assert "inchikey" not in match
    conflict = _candidate("CCCCCCCCCCCCCCCCCCOCC(O)CO")
    conflict["inchikey"] = KEY
    roles.apply_report_compound_roles(conflict)
    assert conflict["potential_contaminant"] is False
    assert conflict["report_role_projection"]["status"] == "identity_conflict_not_applied"
    assert roles.plant_report_selection(conflict)["review_category"] == "identity_conflict"


def test_nested_candidates_are_projected_and_input_evidence_is_preserved():
    match = _candidate("CCCCCCCCCCCCCCCCCCOCC(O)CO", "Alkyl glycerol")
    match["top_candidates"] = [_candidate()]
    original_score = match["top_candidates"][0]["cosine"]
    result = roles.annotate_result_compound_roles([_row(match)])
    assert result[0]["match"]["top_candidates"][0]["potential_contaminant"] is True
    assert result[0]["match"]["top_candidates"][0]["cosine"] == original_score
    assert result[0]["raw_peaks"] == [[70.0, 10.0], [120.0, 20.0]]


def test_high_external_candidate_is_filtered_from_plant_view_and_kept_in_qc():
    row = _row(_candidate())
    source = deepcopy(row)
    payload = build_compound_conclusions("role-filter", [row])
    compound = payload["compounds"][0]
    assert compound["representative"]["match"]["confidence"] == "High"
    assert compound["conclusion"]["counted_as_sample_compound"] is True
    assert compound["conclusion"]["counted_as_plant_candidate"] is False
    assert payload["metadata"]["main_result_count"] == 1
    assert payload["metadata"]["plant_candidate_count"] == 0
    assert payload["metadata"]["external_review_count"] == 1
    assert filter_compound_conclusions(payload, report_scope="plant_candidates")["compounds"] == []
    qc = filter_compound_conclusions(payload, report_scope="external_review")["compounds"]
    assert qc == payload["compounds"]
    assert qc[0]["representative"]["spectrum_id"] == row["spectrum_id"]
    assert qc[0]["representative"]["match"]["lib_spectrum_id"] == "reference-original"
    assert row == source


def test_atomic_primary_projection_cannot_restore_stale_external_role():
    candidate = _candidate()
    ranked = {
        "identity_key": KEY.split("-")[0], "entity_rank": 1,
        "entity_score": 0.98, "candidate": candidate,
    }
    primary = {
        "metadata": {"task_id": "primary-role", "production_identity_selector": True},
        "entities": [{
            "entity_id": "entity:primary-role", "member_result_indices": [0],
            "winner": ranked, "candidates": [ranked],
            "decision": {"status": "selected", "scope": "2d_connectivity_candidate", "blockers": []},
        }],
    }
    original = deepcopy(primary)
    payload = build_compound_conclusions("primary-role", [_row(candidate)], primary_identity_artifact=primary)
    match = payload["compounds"][0]["representative"]["match"]
    assert match["primary_identity_source"] == "v15_flagship"
    assert match["compound_role"] == "external_contaminant"
    assert match["potential_contaminant"] is True
    assert match["compound_name"] == candidate["compound_name"]
    assert payload["compounds"][0]["plant_report"]["excluded_from_plant_candidate_table"] is True
    assert primary == original


def test_mixed_isomer_group_is_retained_intact_in_review_without_new_winner():
    match = _candidate()
    match.update({"identification_scope": "isomer_group", "isomer_ambiguous": True,
                  "isomer_group_formula": "C18H38O4", "isomer_group_confidence": "High"})
    alternate = _candidate("CCCCCCCCCCCCCC(O)C(O)C(O)CCO", "Alternative polyol")
    alternate.update({"rank_score": 0.96, "final_score": 0.96, "spectrally_competitive": True})
    match["isomer_set"] = [match.copy(), alternate]
    payload = build_compound_conclusions("mixed-role", [_row(match)])
    compound = payload["compounds"][0]
    assert compound["plant_report"]["mixed_candidate_roles"] is True
    assert len(compound["candidate_members"]) == 2
    assert compound["lead_candidate"]["inchikey"] == KEY
    assert compound["conclusion"]["counted_as_plant_candidate"] is False
    assert filter_compound_conclusions(payload, report_scope="external_review")["compounds"] == [compound]


def test_report_reader_filters_before_pagination_and_catalog_binds_current_artifact(tmp_path):
    source = tmp_path / "saved.json"
    source.write_text("[]")
    from pipeline.compound_conclusions import compound_conclusion_artifact_path
    path = compound_conclusion_artifact_path(tmp_path, "reader")
    ordinary = _row(_candidate("CCCCCCCCCCCCCCCCCCOCC(O)CO", "Alkyl glycerol"))
    ordinary["spectrum_id"] = "ordinary"
    from pipeline.v15_flagship import write_v15_flagship_artifact, v15_flagship_artifact_path
    rows = [_row(_candidate()), ordinary]
    primary = v15_flagship_artifact_path(tmp_path, "reader")
    write_v15_flagship_artifact("reader", rows, primary)
    write_compound_conclusion_artifact("reader", rows, path, source_artifact=source,
                                       primary_identity_artifact=primary)
    assert artifact_is_current(path, source, primary)
    selected = read_compound_conclusion_artifact(path, limit=1, report_scope="plant_candidates")
    assert selected["compounds_total"] == 1
    assert selected["compounds"][0]["representative"]["match"]["inchikey"] == ordinary["match"]["inchikey"]
    payload = json.loads(path.read_text())
    payload["metadata"].pop("report_role_catalog")
    path.write_text(json.dumps(payload))
    assert not artifact_is_current(path, source, primary)
    with pytest.raises(ValueError, match="Unknown compound report scope"):
        read_compound_conclusion_artifact(path, report_scope="not-a-scope")


def test_catalog_rejects_identity_tampering(monkeypatch, tmp_path):
    payload = json.loads(roles._CATALOG_PATH.read_text())
    payload["annotations"][KEY]["smiles"] = "CCO"
    path = tmp_path / "invalid.json"
    path.write_text(json.dumps(payload))
    monkeypatch.setattr(roles, "_CATALOG_PATH", path)
    with pytest.raises(ValueError, match="structure does not match"):
        roles.report_role_catalog_summary()


def test_actual_overview_main_and_review_filters_keep_separate_evidence_views():
    node = shutil.which("node")
    if node is None:
        pytest.skip("Node.js is required to execute the browser selection logic")
    path = Path(__file__).resolve().parents[1] / "web/static/js/sample-overview.js"
    source = path.read_text()
    status_function = source.split("  function statusFromConclusion(conclusion) {", 1)[1].split(
        "\n  function groupCompounds(", 1
    )[0]
    selection = source.split("    const filtered = currentBundle.groups.filter(group => {", 1)[1].split(
        "    const active = document.getElementById('activeAnalysisFilter');", 1
    )[0]
    program = (
        "function statusFromConclusion(conclusion) {" + status_function + "\n"
        "function select(currentFilter, groups) {\n"
        " const currentBundle = {groups}; const search = '';\n"
        " const matchesAnalysisFilter = () => true;\n"
        " const filtered = currentBundle.groups.filter(group => {" + selection + "\n"
        " return filtered.map(group => group.id); }\n"
        "const group = (id, excluded, contaminant) => ({id,\n"
        " status:statusFromConclusion({key:'high',counted_as_sample_compound:true,counted_as_plant_candidate:!excluded}),\n"
        " compound:{plant_report:{excluded_from_plant_candidate_table:excluded}},\n"
        " rep:{match:{potential_contaminant:contaminant,confidence:'High'}}});\n"
        "const groups=[group('external-high',true,true),group('unknown-origin',false,false),group('mixed-isomers',true,false)];\n"
        "console.log(JSON.stringify({main:select('main',groups),qc:select('contaminant',groups),high:select('high',groups)}));\n"
    )
    result = subprocess.run([node, "-e", program], check=True, capture_output=True, text=True)
    assert json.loads(result.stdout) == {
        "main": ["unknown-origin"],
        "qc": ["external-high", "mixed-isomers"],
        "high": ["external-high", "unknown-origin", "mixed-isomers"],
    }


def test_saved_real_raw_through_server_loader_and_compound_projection_keeps_qc_high(monkeypatch, tmp_path):
    from web import server

    fixture_path = Path(__file__).with_name("fixtures") / "report_role_high_laur3.json"
    case = json.loads(fixture_path.read_text())
    task_id = case["task_id"]
    path = tmp_path / (task_id + "_results.json")
    path.write_text(json.dumps(case["rows"]))
    original_bytes = path.read_bytes()
    monkeypatch.setattr(server, "TASK_DIR", tmp_path)
    task = {"id": task_id, "status": "completed", "results": []}
    monkeypatch.setattr(server, "_tasks", {task_id: task})
    assert server._annotate_result_compound_roles.__module__ == roles.__name__
    # First load is from a saved raw artifact; repeats exercise the cached
    # in-memory path and both annotation passes in the real display pipeline.
    for _ in range(3):
        prepared = server._get_task_results(task_id)
        assert prepared[0]["match"]["confidence"] == "High"
        marker = prepared[0]["match"]["report_role_projection"]
        assert marker["pre_report_contaminant_veto"] is False
        payload = build_compound_conclusions(task_id, prepared, primary_identity_artifact=case["primary"])
        compound = payload["compounds"][0]
        match = compound["representative"]["match"]
        assert compound["conclusion"]["key"] == "high"
        assert compound["conclusion"]["high_confidence_2d_connectivity_detected"] is True
        assert compound["conclusion"]["counted_as_plant_candidate"] is False
        assert match["compound_role"] == "external_contaminant"
        assert match["potential_contaminant"] is True
        assert match["entity_high_evidence"]["supporting_spectra"][0]["report_role_excluded_from_plant_table_only"] is True
        assert filter_compound_conclusions(payload, report_scope="plant_candidates")["compounds"] == []
        assert filter_compound_conclusions(payload, report_scope="external_review")["compounds"] == [compound]
        for field, expected in case["expected_before"]["match"].items():
            assert match.get(field) == expected
        for field, expected in case["expected_before"]["conclusion"].items():
            assert compound["conclusion"].get(field) == expected
        for field in ("compound_id", "identity_key", "identity_scope"):
            assert compound[field] == case["expected_before"][field]
        for field in ("cosine", "mass_ppm", "final_score", "rank_score"):
            assert prepared[0]["match"].get(field) == case["rows"][0]["match"].get(field)
        assert match["rank_score"] == case["primary"]["entities"][0]["winner"]["entity_score"]
        task["results"] = prepared
    assert path.read_bytes() == original_bytes
