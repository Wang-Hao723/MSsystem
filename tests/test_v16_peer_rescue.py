from eval.v16_peer_rescue_20260907.run import peer_choices

def row(qid,run,head='B',eligible=True):
 return {'query_id':qid,'run':run,'polarity':'NEG','mz':200.,'rt':2.,'base':[head],'before':[head],'eligible_targets':['B'],'vote_eligible':eligible}

def target():return row('q','r0','A',False)

def test_two_independent_runs_required():
 assert peer_choices([target(),row('p1','r1')])['q']['head'] is None
 assert peer_choices([target(),row('p1','r1'),row('p2','r2')])['q']['head']=='B'

def test_duplicate_scan_in_one_run_is_not_another_vote():
 assert peer_choices([target(),row('p1','r1'),row('p2','r1')])['q']['head'] is None

def test_current_run_cannot_vote_for_itself():
 assert peer_choices([target(),row('p1','r0'),row('p2','r1')])['q']['head'] is None

def test_conflicting_candidates_in_one_run_abstain():
 out=peer_choices([target(),row('p1','r1'),row('p2','r1','C'),row('p3','r2')])['q']
 assert out['head'] is None and out['ambiguous_runs']==['r1']

def test_peers_without_measured_signal_cannot_vote():
 assert peer_choices([target(),row('p1','r1',eligible=False),row('p2','r2')])['q']['head'] is None

def test_previously_accepted_replacement_is_not_overridden():
 q=target();q['base']=['D']
 assert peer_choices([q,row('p1','r1'),row('p2','r2')])['q']['head'] is None

def test_far_retention_time_does_not_support_target():
 r=row('p1','r1');r['rt']=3.
 assert peer_choices([target(),r,row('p2','r2')])['q']['head'] is None
