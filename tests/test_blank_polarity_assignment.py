from pathlib import Path

from matching.blank_subtraction import (
    assign_blank_profiles_by_polarity,
    create_blank_profile_from_list,
)


def _profile(*polarities):
    profile = create_blank_profile_from_list([(100.0, 1000.0)])
    profile.available_polarities = frozenset(polarities)
    return profile


def test_single_polarity_sample_gets_only_its_matching_blank():
    profiles = [_profile("POS"), _profile("NEG")]

    positive = assign_blank_profiles_by_polarity(
        [{"polarity": "POS"}], profiles
    )
    negative = assign_blank_profiles_by_polarity(
        [{"polarity": "NEG"}], profiles
    )

    assert positive["assigned_indices"] == [0]
    assert positive["skipped_indices"] == [1]
    assert positive["assignment_mode"] == "automatic_by_sample_polarity"
    assert negative["assigned_indices"] == [1]
    assert negative["skipped_indices"] == [0]


def test_mixed_polarity_sample_gets_both_matching_blank_sets():
    assignment = assign_blank_profiles_by_polarity(
        [{"polarity": "POS"}, {"polarity": "NEG"}],
        [_profile("POS"), _profile("NEG")],
    )

    assert assignment["sample_polarities"] == ["NEG", "POS"]
    assert assignment["assigned_indices"] == [0, 1]
    assert assignment["skipped_indices"] == []


def test_mixed_polarity_blank_reports_only_the_sample_polarity_it_supplies():
    assignment = assign_blank_profiles_by_polarity(
        [{"polarity": "POS"}],
        [_profile("POS", "NEG")],
    )

    assert assignment["assigned_indices"] == [0]
    assert assignment["assigned_blank_polarities"] == ["POS"]


def test_unknown_blank_is_retained_but_audited_as_uncertain():
    assignment = assign_blank_profiles_by_polarity(
        [{"polarity": "POS"}],
        [_profile("POS"), _profile(), _profile("NEG")],
    )

    assert assignment["assigned_indices"] == [0, 1]
    assert assignment["unknown_assigned_indices"] == [1]
    assert assignment["skipped_indices"] == [2]
    assert assignment["assignment_mode"] == (
        "automatic_by_sample_polarity_with_unknown_fallback"
    )


def test_analysis_pipeline_passes_only_matching_blanks_to_assessment(
    monkeypatch, tmp_path,
):
    import matching.blank_subtraction as blank_module
    import web.server as server

    sample_path = tmp_path / "positive-sample.bin"
    positive_blank_path = tmp_path / "positive-blank.mzML"
    negative_blank_path = tmp_path / "negative-blank.mzML"
    for path in (sample_path, positive_blank_path, negative_blank_path):
        path.write_bytes(b"placeholder")

    positive_profile = _profile("POS")
    negative_profile = _profile("NEG")
    positive_profile.ms1_scan_count = 1
    negative_profile.ms1_scan_count = 1
    profiles_by_path = {
        positive_blank_path: positive_profile,
        negative_blank_path: negative_profile,
    }
    token_paths = {
        "positive-token": positive_blank_path,
        "negative-token": negative_blank_path,
    }
    captured = {}

    monkeypatch.setattr(server, "_tasks", {
        "auto-pair": {
            "id": "auto-pair",
            "filename": sample_path.name,
            "file_path": str(sample_path),
            "status": "queued",
            "params": {
                "ion_mode": "auto",
                "blank_files": ["positive-token", "negative-token"],
            },
            "results": [],
        }
    })
    monkeypatch.setattr(
        server,
        "_parse_spectrum_file",
        lambda *_args, **_kwargs: [{
            "spectrum_id": "scan=1",
            "precursor_mz": 279.159,
            "retention_time": 10.0,
            "retention_time_unit": "minutes",
            "polarity": "POS",
            "peaks": [[149.0, 1000.0], [205.0, 250.0]],
        }],
    )
    monkeypatch.setattr(
        server,
        "_resolve_blank_file_token",
        lambda token: token_paths[token],
    )
    monkeypatch.setattr(
        blank_module,
        "load_blank_file",
        lambda path: profiles_by_path[path],
    )

    def fake_assess(_spectra, profiles, roles, **_kwargs):
        captured["profiles"] = list(profiles)
        captured["roles"] = list(roles)
        return [{
            "blank_status": "no_blank_feature",
            "blank_feature_matched": False,
            "blank_excluded": False,
            "blank_clean_view_excluded": False,
            "blank_flagged": False,
            "blank_carryover_suspected": False,
            "blank_legacy_below_ratio": False,
        }]

    monkeypatch.setattr(blank_module, "assess_blank_role_groups", fake_assess)
    monkeypatch.setattr(server, "_save_tasks", lambda: None)
    monkeypatch.setattr(server, "_sse_broadcast", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(
        server,
        "search",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(
            RuntimeError("stop after blank assignment")
        ),
    )

    server._run_analysis("auto-pair")

    assert captured["profiles"] == [positive_profile]
    summary = server._tasks["auto-pair"]["blank_summary"]
    assert summary["sample_polarities"] == ["POS"]
    assert summary["assigned_blank_filenames"] == ["positive-blank.mzML"]
    assert summary["assigned_blank_polarities"] == ["POS"]
    assert summary["skipped_blank_filenames"] == ["negative-blank.mzML"]
    assert summary["skipped_blank_polarities"] == ["NEG"]
    assert summary["blank_pairing_status"] == (
        "polarity_only_unverified_method_pairing"
    )
    assert summary["method_comparability_verified"] is False
    assert summary["background_exclusion_eligible"] is False


def test_upload_ui_explains_per_sample_blank_polarity_assignment():
    project_root = Path(__file__).resolve().parents[1]
    upload_js = (
        project_root / "web" / "static" / "js" / "upload.js"
    ).read_text(encoding="utf-8")
    index_html = (
        project_root / "web" / "static" / "index.html"
    ).read_text(encoding="utf-8")

    assert "只把同极性空白分配给该样本" in upload_js
    assert "同极性不等于同批次／同方法配对" in upload_js
    assert "实际排除还要求确认同批次／同方法" in upload_js
    assert "不会退化成 0 扣除" in upload_js
    assert "逐样本分配同极性空白" in index_html
    assert "不证明同批次／同方法配对" in index_html


def test_legacy_unspecified_blank_summary_is_reprojected_as_warning_only():
    import web.server as server

    task = {
        "blank_summary": {
            "requested": True,
            "assigned_blank_replicate_count": 3,
            "blank_assignment_mode": "automatic_by_sample_polarity",
            "blank_roles": {"unspecified": 3},
        }
    }

    server._refresh_blank_pairing_semantics(task)

    summary = task["blank_summary"]
    assert summary["blank_pairing_status"] == (
        "polarity_only_unverified_method_pairing"
    )
    assert summary["background_exclusion_eligible"] is False
    assert summary["unverified_assigned_blank_count"] == 3
    assert "不据此剔除" in summary["blank_pairing_warning"]
