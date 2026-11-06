import numpy as np
from eval.v16_query_selected_reference_20260907.run import vector
from eval.v16_reference_cap_20260907.run import summarize

def test_reference_order_and_duplicates_are_invariant():
    query={'peaks':[[50,1],[70,.5],[90,.2],[110,.1]],'precursor_mz':200}
    rr=[{'peaks':[[50+d,1],[70,.5],[90,.2],[120,.1]],'precursor_mz':200} for d in [0,.001,1,2,3,4]]
    a=vector(query,rr);b=vector(query,rr[::-1]+rr)
    np.testing.assert_array_equal(a,b)

def test_extra_unmatched_references_do_not_dilute_best_three():
    query={'peaks':[[50,1],[70,.5],[90,.2],[110,.1]],'precursor_mz':200}
    good=[{'peaks':[[50+d,1],[70,.5],[90,.2],[120,.1]],'precursor_mz':200} for d in [.001,.002,.003]]
    bad=[{'peaks':[[30+d,1],[40,.5],[160,.2]],'precursor_mz':200} for d in [0,1,2,3,4]]
    np.testing.assert_array_equal(vector(query,good),vector(query,good+bad))

def test_best_reference_score_retained():
    query={'peaks':[[50,1],[70,.5],[90,.2],[110,.1]],'precursor_mz':200}
    rr=[{'peaks':[[50+d,1],[70,.5],[90,.2],[120,.1]],'precursor_mz':200} for d in [0,.001,1,2,3,4]]
    _,full,_=summarize(query,rr)
    assert vector(query,rr)[13]==full[13]

def test_no_reference_does_not_become_zero_score():
    assert vector({'peaks':[[50,1],[70,.5],[90,.2]],'precursor_mz':200},[]) is None
