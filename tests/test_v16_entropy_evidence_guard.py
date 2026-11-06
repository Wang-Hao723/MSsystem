from eval.v16_entropy_evidence_guard_20260907.run import supported_maximum

def test_match_count_cannot_be_borrowed_from_lower_score_reference():
 assert not supported_maximum([(.9,2),(.5,10)])

def test_highest_reference_with_three_matches_is_accepted():
 assert supported_maximum([(.9,3),(.5,10)])

def test_tied_highest_reference_with_three_matches_is_accepted():
 assert supported_maximum([(.9,2),(.9,3)])

def test_empty_references_abstain():
 assert not supported_maximum([])

def test_reference_duplicates_and_order_do_not_change_decision():
 for scores in [[(.9,2),(.5,10)],[ (.9,3),(.5,10)]]:
  assert supported_maximum(scores)==supported_maximum(list(reversed(scores))*3)
