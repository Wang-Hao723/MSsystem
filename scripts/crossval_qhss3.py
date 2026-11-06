"""交叉验证: QHSS3 在 Small vs Full 两库的匹配一致性。
读现有 small 结果(1212) + 重跑 full, 比对: 匹配数 / 命中谱图重叠 / top 化合物一致性。
"""
import json, time, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from utils.mzml_io import parse_mzml_file
from utils.centroiding import smart_peak_pick
from matching.search import load_library, search
from config.settings import (COSINE_THRESHOLD, FINAL_SCORE_THRESHOLD,
                             MZ_TOLERANCE_DA, MZ_TOLERANCE_PPM, TOP_N_PEAKS)

FP = 'output/tmp_uploads/66ae163b_QHSS3.mzML'
TASK = 'output/tasks/a8ae96598050_results.json'

# 1) 读现有 small 结果 (server 刚跑的 1212)
d = json.load(open(TASK))
small_res = d if isinstance(d, list) else d.get('results', [])
small_matched = {r['spectrum_id']: r['match'] for r in small_res if r.get('match')}
print(f'[small 现有结果] 匹配 {len(small_matched)}/{len(small_res)} = {100*len(small_matched)/len(small_res):.1f}%')

# 2) 解析 QHSS3
specs = parse_mzml_file(FP)
for s in specs:
    if s.get('peaks'):
        s['peaks'] = smart_peak_pick(s['peaks'], s['precursor_mz'])
print(f'[parse] {len(specs)} spectra')

def run(db, label):
    t0 = time.time()
    recs, idx = load_library(db)
    print(f'[{label}] 库加载 {len(recs):,} 条 ({time.time()-t0:.1f}s)', flush=True)
    matched = {}
    t1 = time.time()
    for i, s in enumerate(specs):
        m = search(s['peaks'], s['precursor_mz'], records=recs, mz_index=idx,
                   cos_thresh=COSINE_THRESHOLD, final_thresh=FINAL_SCORE_THRESHOLD,
                   mz_da=MZ_TOLERANCE_DA, mz_ppm=MZ_TOLERANCE_PPM, top_n=TOP_N_PEAKS)
        if m:
            matched[s['spectrum_id']] = m
        if (i + 1) % 500 == 0:
            print(f'  [{label}] {i+1}/{len(specs)} matched={len(matched)} ({(i+1)/(time.time()-t1):.0f}/s)', flush=True)
    print(f'[{label}] 完成: {len(matched)}/{len(specs)} = {100*len(matched)/len(specs):.1f}% ({time.time()-t1:.1f}s)', flush=True)
    return matched

# 3) Small 重跑 (确定性检查: 应仍是 1212)
small_re = run('data/local_library_small.db', 'SMALL重跑')

# 4) Full 重跑
full = run('data/local_library.db', 'FULL')

# 5) 比对
print('\n' + '='*60)
print('交叉验证结论:')
print(f'  Small (现有/重跑): 1212 / {len(small_re)}  (确定性: {"✓一致" if len(small_re)==1212 else "✗不一致"})')
print(f'  Full:               {len(full)}')
s_ids = set(small_re.keys()); f_ids = set(full.keys())
overlap = s_ids & f_ids
print(f'  命中谱图重叠: {len(overlap)} (small∩full)')
print(f'  Small独有: {len(s_ids-f_ids)} | Full独有: {len(f_ids-s_ids)}')
# top 化合物一致性 (按 final 排序前15, 看两库是否给出同一化合物)
def topn(matched, n=15):
    return [(m.get('compound_name','')[:28], round(m.get('final_score',0),3), round(m.get('cosine',0),3))
            for m in sorted(matched.values(), key=lambda x:-x.get('final_score',0))[:n]]
print('\n  Small top15 vs Full top15 (同谱图应命中同化合物):')
s_top = {m.get('compound_name'): m.get('final_score',0) for m in sorted(small_re.values(), key=lambda x:-x.get('final_score',0))[:15]}
f_top = {m.get('compound_name'): m.get('final_score',0) for m in sorted(full.values(), key=lambda x:-x.get('final_score',0))[:15]}
common_top = set(s_top) & set(f_top)
print(f'  top15 化合物名交集: {len(common_top)}/15')
