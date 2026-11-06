import numpy as np
from eval.v16_reference_cap_20260907.run import summarize
from eval.v16_public_qtof_20260907.curate import fingerprint

def spectrum(masses,scale=1):return [[m,scale/(i+1)] for i,m in enumerate(masses)]

def test_uncapped_exact_duplicates_and_reverse_do_not_add_evidence():
    q={'peaks':spectrum([50,70,90,110]),'precursor_mz':200}
    refs=[{'peaks':spectrum([50+x,70,90,120]),'precursor_mz':200} for x in [.001,.005,1,2,3]]
    a,b,t=summarize(q,refs);aa,bb,tt=summarize(q,refs[::-1]+refs)
    np.testing.assert_array_equal(a,aa);np.testing.assert_array_equal(b,bb);assert t==tt and t['valid']==5

def test_processed_query_mirror_excluded_even_with_different_precursor_peak():
    q={'peaks':spectrum([50,70,90]),'precursor_mz':200}
    refs=[{'peaks':q['peaks']+[[200,.1]],'precursor_mz':200}]
    assert summarize(q,refs)==(None,None,{'valid':0})

def test_all_reference_summary_includes_nonfirst_three_information():
    q={'peaks':spectrum([50,70,90]),'precursor_mz':200}
    refs=[{'peaks':spectrum([30+i,60+i,95+i]),'precursor_mz':200} for i in range(12)]
    refs.sort(key=lambda r:fingerprint(r['peaks']))
    a,b,t=summarize(q,refs)
    assert t['valid']==12 and len(a)==len(b)==56
    assert t['max_sqrtcos_all']>=t['max_sqrtcos_first3']

def test_missing_reference_remains_unknown():
    assert summarize({'peaks':spectrum([50,70,90]),'precursor_mz':200},[])==(None,None,{'valid':0})
