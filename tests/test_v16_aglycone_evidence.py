from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from eval.v16_aglycone_evidence_20260907.chemistry import aglycones
from eval.v16_aglycone_evidence_20260907.run import veto

PARENT='COc1ccc(-c2coc3cc(OC4OC(CO)C(O)C(O)C4O)ccc3c2=O)cc1O'
OTHER='COc1ccc(-c2coc3cc(OC4OC(CO)C(O)C(O)C4O)cc(O)c3c2=O)cc1'

def test_neutral_glycosidic_loss_and_intact_aglycone():
    result=aglycones(PARENT);assert len(result)==1
    a=result[0];assert a['formula']=='C16H12O5'
    assert abs(rdMolDescriptors.CalcExactMolWt(Chem.MolFromSmiles(PARENT))-a['neutral_mass']-162.05282342)<1e-6
    assert not any(x.GetNumRadicalElectrons() for x in Chem.MolFromSmiles(a['smiles']).GetAtoms())

def test_distinct_aglycone_isomers_preserved():
    assert aglycones(PARENT)[0]['connectivity']!=aglycones(OTHER)[0]['connectivity']

def test_free_aglycone_or_benzene_not_removed():
    assert aglycones('c1ccccc1')==[]
    assert aglycones(aglycones(PARENT)[0]['smiles'])==[]

def test_same_aglycone_cannot_veto_linkage_isomers():
    class NoScoring:
        def score(self,*args):raise AssertionError('must abstain')
    after,evidence=veto(['A','B'],['B','A'],{'A':PARENT,'B':PARENT},'POS',[],NoScoring())
    assert after==['B','A'] and evidence['status']=='same_aglycone_cannot_distinguish_glycosides'

def test_missing_reference_not_negative_evidence():
    class Missing:
        def score(self,*args):return None,{'status':'missing'}
    after,evidence=veto(['A','B'],['B','A'],{'A':PARENT,'B':OTHER},'POS',[],Missing())
    assert after==['B','A'] and evidence['status']=='insufficient_paired_aglycone_evidence'
