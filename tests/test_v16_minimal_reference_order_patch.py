from copy import deepcopy

import matching.flagship_ranker as ranker
from matching.v15_flagship_ranker import load_v15_flagship_policy
from eval.v16_fragment_discrimination_20260907.reference_order_patch import patched_candidate_rows
from eval.v16_fragment_discrimination_20260907.projection_repair import reverse_row
from tests.test_v16_reference_projection_repair import fixture


def scores(row):
    return {r["identity_key"]:r["spectral_score"] for r in ranker._candidate_rows(row, load_v15_flagship_policy())}


def test_current_and_legacy_variant_are_order_invariant_and_isolated():
    row = fixture()
    original = deepcopy(row)
    baseline = scores(row)
    reversed_baseline = scores(reverse_row(row))
    # The fix is already deployed: current production must not reproduce the
    # historical defect. Keep the fixture that previously exposed it.
    assert baseline == reversed_baseline
    with patched_candidate_rows():
        # This older compact-only experiment uses a different reference policy;
        # it need not equal the newer production score, but must be invariant.
        assert scores(row) == scores(reverse_row(row))
    assert row == original
    assert scores(reverse_row(row)) == reversed_baseline  # Isolated patch restores the module.
    assert scores(row) == baseline


def test_patch_context_restores_function_even_after_exception():
    original = ranker._candidate_rows
    try:
        with patched_candidate_rows():
            assert ranker._candidate_rows is not original
            raise RuntimeError("test cleanup")
    except RuntimeError:
        pass
    assert ranker._candidate_rows is original
