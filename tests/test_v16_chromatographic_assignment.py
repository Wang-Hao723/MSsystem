import pytest
from eval.v16_chromatographic_assignment_20260907.run import solve,stereo_eligible

def nodes(overlap=False,mz=100.,polarity='POS'):
    return [dict(feature_id='a',polarity='POS',mz=100.,start=1.,end=2.,scores={'X':10.,'Y':8.}),
            dict(feature_id='b',polarity=polarity,mz=mz,start=1.5 if overlap else 3.,end=4.,scores={'X':9.,'Y':1.})]

def test_unique_global_solution_can_keep_stronger_pair_of_assignments():
    a,v,n=solve(nodes(),eligible={'X','Y'})
    assert a=={('POS','a'):'Y',('POS','b'):'X'} and v==17 and n==2
    _,alt,_=solve(nodes(),eligible={'X','Y'},forbidden=('POS','a','Y'))
    assert alt<v

@pytest.mark.parametrize('kwargs',[dict(overlap=True),dict(mz=100.1),dict(polarity='NEG')])
def test_no_exclusion_without_matching_ion_and_disjoint_interval(kwargs):
    a,v,n=solve(nodes(**kwargs),eligible={'X','Y'})
    assert set(a.values())=={'X'} and v==19 and n==0

def test_unknown_or_stereo_ambiguous_keys_not_exclusive():
    a,v,n=solve(nodes())
    assert set(a.values())=={'X'} and v==19 and n==0

def test_explicit_zero_scores_cannot_support_unique_replacement():
    _,v,_=solve(nodes(),eligible={'X','Y'},zero=True)
    _,v2,_=solve(nodes(),eligible={'X','Y'},zero=True,forbidden=('POS','a','Y'))
    assert v==v2==0

def test_input_reversal_invariant():
    assert solve(nodes(),eligible={'X','Y'})==solve(nodes()[::-1],eligible={'X','Y'})

@pytest.mark.parametrize('smiles',['CC(O)C(=O)O','C/C=C/C','FC(Cl)(Br)I','CC.O','invalid',None])
def test_exclude_potential_stereo_multicomponent_invalid(smiles):
    assert not stereo_eligible(smiles)

@pytest.mark.parametrize('smiles',['CCO','c1ccccc1','CC(=O)O'])
def test_achiral_single_molecule_eligible(smiles):
    assert stereo_eligible(smiles)
