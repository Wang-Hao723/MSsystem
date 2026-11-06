from web.expected_compounds_export import expected_group_rows


def test_expected_group_projection_omits_raw_xic_rows():
    payload = {
        "features": [
            {
                "feature_id": "f1",
                "candidate_name": "Candidate A",
                "candidate_names": ["Candidate A"],
                "formula": "C10H12O2",
                "adduct": "[M+H]+",
                "polarity": "POS",
                "target_mz": 165.091,
                "observed_apex_mz": 165.0912,
                "mass_error_ppm": 1.2,
                "apex_rt_minutes": 2.5,
                "apex_intensity": 10000,
                "xic_area_intensity_seconds": 50000,
                "evidence_level": "MS1_ISOTOPE_SUPPORTED",
                "cd_pattern_accepted": True,
                "expected_ion_id": "ion-1",
                "dataset_sources": ["catalog"],
                "target_members": [
                    {
                        "compound_name": "Candidate A",
                        "identity_key": "ik::AAAAAAAAAAAAAA-BBBBBBBBBB-C",
                        "inchikey": "AAAAAAAAAAAAAA-BBBBBBBBBB-C",
                        "dataset_source": "catalog",
                    }
                ],
                "formula_isotope_evidence": {"status": "supporting"},
            }
        ]
    }

    groups = expected_group_rows(payload)

    assert len(groups) == 1
    assert groups[0]["candidate_display"] == "Candidate A"
    assert groups[0]["supported_feature_count"] == 1
    assert groups[0]["raw_peak_count"] == 1
    assert "_Raw_Rows" not in groups[0]
    assert set(groups[0]).isdisjoint({"target_members", "formula_isotope_evidence"})


def test_expected_group_uses_clean_primary_name_and_hides_catalogue_codes():
    payload = {
        "features": [{
            "feature_id": "f-mixed-names",
            "candidate_names": [
                "NP-000004", "NP-000062", "Flavone",
                "Schaftoside [M+H]+",
            ],
            "formula": "C26H28O14",
            "adduct": "[M-H]-",
            "polarity": "NEG",
            "target_mz": 563.141,
            "observed_apex_mz": 563.141,
            "evidence_level": "MS1_EXACT_MASS_CANDIDATE",
            "expected_ion_id": "ion-mixed",
            "dataset_sources": ["F001_catalogue"],
            "target_members": [
                {
                    "target_id": "NP-000004",
                    "compound_name": "NP-000004",
                    "identity_key": "target::NP-000004",
                },
                {
                    "target_id": "NP-000062",
                    "compound_name": "NP-000062",
                    "identity_key": "target::NP-000062",
                },
                {
                    "target_id": "schaft",
                    "compound_name": "Schaftoside [M+H]+",
                    "inchikey": "MMDUKUSNQNWVET-UHFFFAOYSA-N",
                    "identity_key": "ik::MMDUKUSNQNWVET-UHFFFAOYSA-N",
                },
                {
                    "target_id": "flavone",
                    "compound_name": "Flavone",
                    "identity_key": "target::flavone",
                },
            ],
        }]
    }

    group = expected_group_rows(payload)[0]
    assert group["candidate_display"] == "Schaftoside"
    assert "[M+H]+" not in group["candidate_display"]
    assert "NP-000004" not in group["candidate_names"]
    assert "NP-000062" not in group["candidate_names"]
    assert group["ambiguity_scope"] == (
        "same_formula_or_mass_database_candidates"
    )
