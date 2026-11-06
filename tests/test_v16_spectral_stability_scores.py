import numpy as np
import pytest
from eval.v16_candidate_idf_20260907.run import mass_index,weights
from eval.v16_mass_precision_20260907.run import pair as mass_pair
from eval.v16_entropy_stability_20260907.run import pair as entropy_pair

def test_fragment_frequency_counts_candidates_not_references():
 a=np.array([[50.,1.],[70.,.5],[90.,.2]]);b=np.array([[50.,1.],[80.,.5],[100.,.2]])
 one=mass_index({'a':[(a,'a')],'b':[(b,'b')]});duplicates=mass_index({'a':[(a,'a')]*10,'b':[(b,'b')]})
 np.testing.assert_allclose(weights(np.array([50.,70.]),one,1),weights(np.array([50.,70.]),duplicates,1))
 assert weights(np.array([50.,70.]),one,1)[1]>weights(np.array([50.,70.]),one,1)[0]

@pytest.mark.parametrize('window,expected',[('5mDa',0.),('10mDa',1.),('20ppm_floor2mDa',0.)])
def test_mass_window_uses_requested_absolute_or_relative_tolerance(window,expected):
 q=np.array([[100.,1.],[200.,.5],[300.,.2]]);r=q.copy();r[:,0]+=.007
 score,_=mass_pair(q,r,window)
 assert score==pytest.approx(expected)

def test_ppm_window_scales_with_fragment_mass():
 q=np.array([[500.,1.],[600.,.5],[700.,.2]]);r=q.copy();r[:,0]+=.009
 assert mass_pair(q,r,'20ppm_floor2mDa')[0]==pytest.approx(1.)
 assert mass_pair(q,r,'5mDa')[0]==pytest.approx(0.)

@pytest.mark.parametrize('method',['entropy','unweighted_entropy'])
def test_entropy_respects_scaling_order_and_disjoint_spectra(method):
 q=np.array([[50.,1.],[70.,.5],[90.,.2]]);scaled=q.copy();scaled[:,1]*=100
 assert entropy_pair(q,scaled[::-1],method)[0]==pytest.approx(1,abs=1e-6)
 r=q.copy();r[:,0]+=5
 score,ratio=entropy_pair(q,r,method)
 assert score==pytest.approx(0.,abs=1e-6) and ratio==0
