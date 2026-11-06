"""Retry only memory-limited read audits; never execute an analysis pipeline."""
from pathlib import Path
import json
import subprocess
import sys
import types

# This CLI imports the pure access bridge, never web.__init__/server or models.
_web = types.ModuleType('web')
_web.__path__ = [str(Path(__file__).resolve().parents[1] / 'web')]
sys.modules.setdefault('web', _web)

from web.completed_batches import CompletedBatchService
from web.frozen_batch_tasks import FrozenBatchTasks

ROOT = Path(__file__).resolve().parents[2]
WORKER = r'''
import collections,gc,json,resource,sys
resource.setrlimit(resource.RLIMIT_AS,(8*1024**3,8*1024**3))
request=json.load(sys.stdin)
namespace={'__name__':'saved_reader'}
exec(request.pop('reader_source'),namespace)
native=namespace['read_document'](request,'results')
rows=native['native_ms2_results']; layer=native['ms1_integration']; features=layer['features']
summary=request['summary']; ms1=request['ms1']
assert len(rows)==summary['parsed_searched_ms2_spectra']
assert sum(bool(row.get('match')) for row in rows)==summary['matched_ms2_spectra']
assert len(features)==summary['ms1_feature_records']==ms1['ms1_feature_records']
assert layer['static_artifact_sha256']==ms1['result_sha256']
assert layer['static_artifact']==ms1['artifact']
assert layer['static_protocol_sha256']==ms1['protocol_sha256']
assert not any(f.get('standard_confirmed') or f.get('exact_identity_confirmed') for f in features)
assert {b['path']:b['sha256'] for b in native['paired_blanks']}==request['paired_blanks']
assert all(b['role']=='blank' for b in native['paired_blanks'])
record={'sample_key':request['sample_key'],'filename':native['source']['filename'],
 'native_searched_msms':len(rows),'native_matched_msms':sum(bool(r.get('match')) for r in rows),
 'ms1_features':len(features),'ms1_only_features':sum(not f.get('associated_ms2') for f in features),
 'background_flagged':sum(bool(r.get('blank_flagged')) for r in rows),
 'clean_excluded':sum(bool(r.get('blank_clean_view_excluded')) for r in rows),
 'parameters':{k:native['parameters'].get(k) for k in ['enable_ms1_feature_detection','enable_expected_compounds','expected_compound_datasets','blank_filter_policy','blank_method_comparability_confirmed']}}
assert record['background_flagged']==summary['blank_flagged_ms2']
assert record['clean_excluded']==summary['blank_clean_view_excluded_ms2']
del native,rows,layer,features
gc.collect()
document=namespace['read_document'](request,'conclusions')
assert document['conclusions_protocol_sha256']==request['conclusion_protocol_sha256']
assert document['ms1_independent_layer']['static_artifact_sha256']==ms1['result_sha256']
assert document['ms1_independent_layer']['feature_count']==record['ms1_features']
full=document['compound_conclusions']; compounds=full['compounds']
assert len(compounds)==summary['compound_result_records']
# Saved identity IDs can repeat across feature/event rows; use array indices.
id_counts=collections.Counter(r.get('compound_id') for r in compounds)
subset=document['report_views']['xbt_named_candidates']['compounds']
assert len(subset)==summary['xbt_named_candidate_rows']
record.update({'native_conclusions':len(compounds),
 'native_record_types':dict(collections.Counter(r.get('record_type') for r in compounds)),
 'native_statuses':dict(collections.Counter((r.get('conclusion') or {}).get('key') for r in compounds)),
 'native_metadata':{k:v for k,v in full['metadata'].items() if k.endswith('_count') and isinstance(v,(int,float))},
 'xbt_joint_subset_rows':len(subset),
 'xbt_joint_record_types':dict(collections.Counter(r.get('record_type') for r in subset)),
 'xbt_joint_statuses':dict(collections.Counter((r.get('conclusion') or {}).get('key') for r in subset)),
 'matched_sha_verified':True,'ms1_binding_verified':True,
 'paired_blank_binding_verified':True,'native_full_conclusion_view_retains_all_records':True,
 'saved_compound_id_repeat_rows':sum(n-1 for n in id_counts.values()),
 'read_peak_rss_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss})
print(json.dumps(record,ensure_ascii=False))
'''

PARENT = r'''
import collections,json,os,shutil,subprocess,sys,time
from pathlib import Path
payload=json.load(sys.stdin)
root=Path('/data/Mssystem/xbt_real_batch_20261021_ms2/delivery_access_audit_20261022')
assert shutil.disk_usage(root).free>25*1024**3
old=json.loads((root/'RESULT.json').read_text())
failures={r['sample_key'] for r in old['errors']}
assert all(r['error'].startswith('MemoryError:') for r in old['errors'])
records=[];errors=[];start=time.monotonic()
def write(name,data):
 p=root/name;tmp=p.with_suffix(p.suffix+'.tmp');tmp.write_text(json.dumps(data,ensure_ascii=False,indent=2));tmp.replace(p)
for request in payload['requests']:
 if request['sample_key'] not in failures: continue
 request['reader_source']=payload['reader_source']
 child=subprocess.Popen(['python3','-c',payload['worker']],stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
 write('RETRY_STATE.json',{'state':'reading_saved_results','pid':os.getpid(),'child_pid':child.pid,'sample_key':request['sample_key'],'completed':len(records),'denominator':len(failures),'started_at':time.time()})
 out,err=child.communicate(json.dumps(request),timeout=240)
 if child.returncode: errors.append({'sample_key':request['sample_key'],'error':err[-2500:]})
 else: records.append(json.loads(out))
 write('RETRY_RESULT.json',{'state':'running','completed':len(records),'errors':errors,'samples':records})
 print(json.dumps({'retry_completed':len(records),'failed':len(errors),'denominator':len(failures)}),flush=True)
all_records=old['samples']+records
totals={}
for key in ['native_searched_msms','native_matched_msms','native_conclusions','ms1_features','ms1_only_features','xbt_joint_subset_rows','background_flagged','clean_excluded']:
 totals[key]=sum(r.get(key,0) for r in all_records)
for key in ['native_record_types','native_statuses','xbt_joint_record_types','xbt_joint_statuses']:
 counter=collections.Counter()
 for r in all_records:counter.update(r[key])
 totals[key]=dict(counter)
result={'state':'passed' if not errors and len(all_records)==149 else 'failed','completed':len(all_records),'failed':len(errors),'denominator':149,'totals':totals,'errors':errors,'samples':sorted(all_records,key=lambda r:r['sample_key']),'retry_elapsed_seconds':time.monotonic()-start,'original_failure_receipt':'RESULT.json','original_failure_preserved':True,'no_analysis_reexecution':True,'no_raw_files_read':True,'models_modified':False}
write('COMPLETE_RESULT.json',result);write('RETRY_STATE.json',{k:v for k,v in result.items() if k!='samples'})
print(json.dumps({k:v for k,v in result.items() if k!='samples'},ensure_ascii=False),flush=True)
'''

def main():
    batch = ROOT / 'research/20261021_xbt_ms2_integrated'
    service = CompletedBatchService(batch,Path('/Volumes/LaCie/Mssystem_V16/completed_batch_exports'),'msqc-server')
    bridge = FrozenBatchTasks(service)
    ms1 = json.loads((ROOT/'research/20261021_xbt_ms1_targeted/COMPLETE_MS1_RESULTS_SUMMARY.json').read_text())
    by_path = {r['file']:r for r in ms1['files']}
    manifest = json.loads((ROOT/'research/20261021_new_real_samples_preflight/INPUT_MANIFEST.json').read_text())
    paths = {r['path']:r for r in manifest}
    pairs = json.loads((ROOT/'research/20261021_new_real_samples_preflight/BLANK_PAIRING.json').read_text())
    pairs_by_sample = {s:p['blanks'] for p in pairs for s in p['samples']}
    conclusions = json.loads((batch/'CONCLUSIONS_STATE.json').read_text())
    requests = []
    for tid,entry in bridge.entries.items():
        source = bridge.tasks[tid]['file_path']
        requests.append({**entry,'summary':bridge.tasks[tid]['completed_batch']['summary'],
            'ms1':by_path[source],'paired_blanks':{p:paths[p]['sha256'] for p in pairs_by_sample[source]},
            'conclusion_protocol_sha256':conclusions['protocol_sha256']})
    reader = (ROOT/'mssystem_V16/web/frozen_batch_reader.py').read_text()
    program = ['ssh','-o','BatchMode=yes','msqc-server','python3 -c '+__import__('shlex').quote(PARENT)]
    run=subprocess.run(program,input=json.dumps({'reader_source':reader,'worker':WORKER,'requests':requests}),text=True)
    raise SystemExit(run.returncode)

if __name__=='__main__':
    main()
