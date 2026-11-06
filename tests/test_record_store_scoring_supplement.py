import copy
import json

import pytest

from eval.v16_evidence_validation import record_store_scoring_supplement as supplement


def cohort(tmp_path):
    primary = {"scoring": {"details": []}, "input_file_sha256": {}}
    for index in range(46):
        task = f"sample{index:02d}"
        polarity = "POS" if index % 2 else "NEG"
        primary["scoring"]["details"].append({"task_id": task, "polarity": polarity})
        row = {"match": {"compound_name": "Never passed to the scorer", "final_score": 1.},
            "raw_peaks_complete": True, "raw_peaks": [[20., 100.], [35., 50.]],
            "precursor_mz": 100., "search_precursor_mz": 100.001, "ion_mode": polarity}
        rows = [{**row, "match": None}, {**row, "raw_peaks_complete": False}]
        rows += [copy.deepcopy(row) for _ in range(6)]
        path = tmp_path / f"{task}_results.json"
        path.write_text(json.dumps(rows))
        primary["input_file_sha256"][str(path)] = supplement.file_sha256(path)
    return primary


def test_supplement_freezes_first_five_eligible_rows_without_candidate_labels(tmp_path):
    primary = cohort(tmp_path)
    queries = supplement.select_queries(primary, 5)
    assert len(queries) == 230
    assert [item["row_index"] for item in queries[:5]] == [2, 3, 4, 5, 6]
    assert all(item["query_pmz"] == 100.001 for item in queries)
    assert all("match" not in item and "compound_name" not in item["kwargs"] for item in queries)
    assert {item["polarity"] for item in queries} == {"POS", "NEG"}


def test_supplement_rejects_incomplete_registered_task_cohort(tmp_path):
    primary = cohort(tmp_path)
    primary["scoring"]["details"].pop()
    with pytest.raises(ValueError, match="cohort"):
        supplement.select_queries(primary, 5)


def test_supplement_does_not_skip_a_registered_polarity_mismatch(tmp_path):
    primary = cohort(tmp_path)
    primary["scoring"]["details"][0]["polarity"] = "POS"
    with pytest.raises(ValueError, match="polarity"):
        supplement.select_queries(primary, 5)


def test_supplement_is_not_an_identity_accuracy_experiment():
    protocol = json.loads(supplement.PROTOCOL.read_text())
    assert protocol["identity_accuracy_claim_allowed"] is False
    assert protocol["minimum_nonempty_current_results"] == 100
    assert protocol["allowed_mismatches"] == 0
