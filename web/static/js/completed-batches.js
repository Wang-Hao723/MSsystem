/** Browse frozen completed batches without entering the analysis task queue. */
const CompletedBatches = (() => {
  const base = '/api/completed-batches';
  const state = { batches: [], detail: null, batch: '', table: 'native_compounds', sample: '', offset: 0, limit: 50, request: 0, listPromise: null, tab: 'overview' };
  const fieldNames = {
    sample_key: '样本文件键', sample_id: '样本编号', polarity: '极性', replicate: '重复针',
    display_name: '原始候选／事件标题（QC）', candidate_name: 'MS1候选名称', candidate_names: 'MS1候选名称集合',
    frozen_compound_index: '冻结记录索引', record_type: '原生记录类型', conclusion_key: '原生判定状态',
    lead_compound_name: '谱库假设标题（QC）', representative_spectrum_id: '代表谱图编号',
    plant_report_review_category: '来源复核分类', plant_report_excluded: '植物正文表排除标记',
    formula: '分子式', inchikey: 'InChIKey', identity_scope: '身份限定', isomer_limited: '保留异构体限定',
    isomer_ambiguous: '异构体未区分', ms1_evidence_levels: 'MS1证据', msms_evidence: 'MS/MS证据',
    blank_statuses: '配对空白状态', blank_ratio: '空白比值', native_rank: '保存的排序',
    linked_feature_ids: '关联MS1特征', linked_spectrum_ids: '关联MS/MS扫描',
    spectrum_id: 'MS/MS扫描ID', precursor_mz: '前体 m/z', search_precursor_mz: '检索前体 m/z',
    retention_time: '保留时间', retention_time_unit: '时间单位', matched: '谱库匹配',
    match_name: '谱库候选名称', match_formula: '谱库分子式', match_inchikey: '谱库InChIKey',
    final_score: '原始综合分数', cosine: '余弦相似度', evidence_tier: '证据层级',
    identification_scope: '候选限定', confidence: '保存的置信层级', honesty_tier: '校准层级',
    potential_contaminant: '潜在污染物', exclude_from_natural_product_count: '排除植物内源计数',
    contaminant_category: '污染物类别', natural_product_candidate: '天然产物候选标记',
    name_resolution_status: '名称映射状态', isomer_group_size: '异构体候选数',
    isomer_resolution_status: '异构体解析状态', blank_status: '配对空白状态', blank_sample_ratio: '空白比值',
    native_lcms_feature_id: '原生LC-MS特征ID', native_precursor_feature_id: '原生前体特征ID',
    targeted_ms1_feature_ids_from_association: '靶向MS1关联特征ID', source_artifact_path: '完整来源文件',
    feature_id: 'MS1特征ID', observed_apex_mz: '实测峰顶 m/z', target_mz: '靶标 m/z',
    apex_rt_minutes: '峰顶RT（分钟）', apex_intensity: '峰顶强度',
    xic_area_intensity_seconds: '原始XIC面积（强度·秒）',
    baseline_corrected_area_intensity_seconds: '基线校正面积（强度·秒）', mass_error_ppm: '质量误差（ppm）',
    adduct: '加合物', isotope_envelope: '同位素包络', formula_isotope_evidence: '分子式同位素证据',
    blank_clean_view_excluded: '从空白净视图排除', associated_ms2_spectrum_ids: '关联MS/MS扫描',
    associated_ms2_count: '关联MS/MS扫描数', ms1_only: '仅MS1（无MS2关联）', evidence_level: '证据层级',
    qc_category: 'QC类别', feature_ids: 'MS1特征ID', spectrum_ids: 'MS/MS扫描ID',
    name: '候选名称', exclusion_reasons: 'QC复核原因', evidence: '保存的证据层级',
    source_sidecar_path: '完整联合结果文件', source_path: '原始样本路径', source_sha256: '原始样本SHA-256',
    declared_ms2_scans: '声明MS/MS扫描数', searched_ms2_spectra: '实际检索MS/MS事件',
    matched_ms2_spectra: 'MS/MS匹配记录', ms1_features: 'MS1特征数',
    main_path: '完整MS/MS文件', main_sha256: 'MS/MS文件SHA-256',
    conclusion_path: '完整联合结果文件', conclusion_sha256: '联合结果SHA-256',
  };
  const exportNames = { file_manifest: '149份文件索引', xbt_csv: 'XBT联合附加筛选（含未解析事件）', xbt_jsonl: 'XBT联合附加筛选完整字段', qc: '四类QC完整字段', ms1: '全部MS1特征完整字段', msms_table: '全部MS/MS前体事件表', ms1_table: '全部MS1特征表（含MS1-only）', qc_table: '四类QC记录表', summary: '既有MS1/MSMS/QC及XBT筛选摘要' };
  const priorityColumns = {
    native_compounds: ['sample_key', 'frozen_compound_index', 'record_type', 'conclusion_key', 'lead_compound_name', 'formula', 'identity_scope', 'representative_spectrum_id', 'blank_flagged', 'blank_clean_view_excluded', 'plant_report_review_category'],
    xbt_csv: ['display_name', 'formula', 'sample_id', 'polarity', 'replicate', 'isomer_limited', 'ms1_evidence_levels', 'msms_evidence', 'blank_statuses'],
    msms_table: ['sample_id', 'polarity', 'spectrum_id', 'precursor_mz', 'retention_time', 'matched', 'match_name', 'final_score', 'cosine', 'blank_status', 'blank_sample_ratio'],
    ms1_table: ['sample_id', 'polarity', 'feature_id', 'observed_apex_mz', 'apex_rt_minutes', 'candidate_name', 'evidence_level', 'blank_status', 'ms1_only', 'associated_ms2_count'],
    qc_table: ['qc_category', 'sample_id', 'polarity', 'name', 'formula', 'evidence', 'identity_scope', 'exclusion_reasons'],
    file_manifest: ['sample_id', 'polarity', 'replicate', 'searched_ms2_spectra', 'matched_ms2_spectra', 'ms1_features', 'sample_key'],
  };
  const $ = id => document.getElementById(id);
  const esc = value => String(value ?? '').replace(/[&<>"']/g, ch => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[ch]));
  const num = value => Number(value || 0).toLocaleString('zh-CN');
  const bytes = value => `${(Number(value || 0) / 1048576).toFixed(2)} MiB`;
  const sampleName = row => `${row.sample_id || row.sample_key || ''} · ${row.polarity || ''}${row.replicate ?? ''}`;

  async function json(url) {
    const response = await fetch(url, { cache: 'no-store' });
    if (!response.ok) {
      const error = await response.json().catch(() => ({}));
      throw new Error(typeof error.detail === 'string' ? error.detail : `读取失败（HTTP ${response.status}）`);
    }
    return response.json();
  }
  function status(id, message = '', error = false) {
    const node = $(id); if (!node) return;
    node.textContent = message; node.classList.toggle('is-error', error);
  }
  function batchHref(id) { return `/#completed-batches?batch=${encodeURIComponent(id)}`; }
  function entries() {
    document.querySelectorAll('[data-completed-batch-entry]').forEach(host => {
      host.innerHTML = state.batches.map(batch => `<section class="cb-entry"><div><strong>${esc(batch.title)}</strong><p>已完成 ${num(batch.completed_samples)} / ${num(batch.total_samples)} 份 · ${num(batch.paired_blank_files)} 份配对空白 · MS1、MS/MS 联合结果</p></div><div class="cb-action-row"><a class="btn btn-primary btn-sm" href="${batchHref(batch.id)}">查看整批结果</a><a class="btn btn-secondary btn-sm" href="${base}/${encodeURIComponent(batch.id)}/native-report" target="_blank" rel="noopener">完整原生结果报告</a></div></section>`).join('');
    });
  }
  async function loadList(force = false) {
    if (state.listPromise && !force) return state.listPromise;
    state.listPromise = (async () => {
      const payload = await json(base);
      state.batches = payload.batches || []; entries();
      const select = $('completedBatchSelect');
      if (select) select.innerHTML = state.batches.length ? state.batches.map(row => `<option value="${esc(row.id)}">${esc(row.title)}</option>`).join('') : '<option value="">暂无已完成的整批结果</option>';
      return state.batches;
    })().catch(error => { state.listPromise = null; throw error; });
    return state.listPromise;
  }
  async function enterPage(force = false) {
    status('completedBatchStatus', '读取已完成的联合结果…');
    try {
      const batches = await loadList(force);
      const requested = new URLSearchParams(window.location.hash.split('?')[1] || '').get('batch');
      const id = batches.find(row => row.id === requested)?.id || batches.find(row => row.id === state.batch)?.id || batches[0]?.id;
      if (!id) { status('completedBatchStatus', '暂无已完成的整批结果。'); return; }
      $('completedBatchSelect').value = id;
      if (state.detail && state.batch === id && !force) { status('completedBatchStatus'); return; }
      await selectBatch(id);
    } catch (error) { status('completedBatchStatus', error.message, true); }
  }
  async function selectBatch(id) {
    if (!id) return;
    status('completedBatchStatus', '读取批次索引…');
    try {
      const [batch, native] = await Promise.all([json(`${base}/${encodeURIComponent(id)}`), json(`${base}/${encodeURIComponent(id)}/native-delivery`)]);
      state.detail = {...batch, native_delivery: native}; state.batch = id;
      state.offset = 0; state.sample = '';
      renderDetail(); $('completedBatchContent').hidden = false;
      status('completedBatchStatus'); await loadPreview();
    } catch (error) { $('completedBatchContent').hidden = true; status('completedBatchStatus', error.message, true); }
  }
  function renderDetail() {
    const batch = state.detail, totals = batch.totals || {}, samples = batch.samples || [], native = batch.native_delivery;
    $('completedBatchTitle').textContent = batch.title;
    $('completedBatchSubtitle').textContent = `${batch.xbt_dataset || 'XBT'} · ${num(batch.completed_samples)} / ${num(batch.total_samples)} 份已完成 · ${num(batch.paired_groups)} 个同编号、同极性空白配对组`;
    const metrics = [ ['样本文件', `${num(batch.completed_samples)} / ${num(batch.total_samples)}`, '联合分析已完成'], ['配对空白文件', num(batch.paired_blank_files), `${num(batch.paired_groups)} 个配对组`], ['MS/MS前体事件', num(totals.searched_ms2_spectra), `有谱库候选 ${num(totals.matched_ms2_spectra)} 张谱`], ['MS1特征记录', num(totals.ms1_features), '包含无MS2关联的特征'], ['完整原生结论记录', num(native.native_conclusion_rows), '含候选、背景、未解析和未匹配记录'], ['XBT联合附加筛选', num(totals.xbt_named_candidate_rows), '按原生类型区分命名候选与未解析事件'] ];
    $('completedBatchMetrics').innerHTML = metrics.map(([label, value, note]) => `<div class="cb-metric"><span>${esc(label)}</span><strong>${esc(value)}</strong><small>${esc(note)}</small></div>`).join('');
    const levels = batch.ms1_evidence_counts || {};
    $('completedBatchEvidence').textContent = `按跨文件记录出现次数计数。MS1质量匹配候选 ${num(levels.MS1_EXACT_MASS_CANDIDATE)} 条，同位素支持 ${num(levels.MS1_ISOTOPE_SUPPORTED)} 条，MS1与MS/MS诊断碎片支持候选 ${num(levels.MS1_MS2_DIAGNOSTIC_PROBABLE)} 条。标准品确证 ${num(totals.standard_confirmed)} 条；MSClassify未运行。`;
    $('completedBatchDownloadAll').href = native.download_all_url;
    $('completedBatchOpenReport').href = native.report_url;
    $('completedBatchDownloadReport').href = `${native.report_url}?download=true`;
    $('completedBatchDownloadIndex').href = native.index_url;
    $('completedBatchReportFrame').removeAttribute('src');
    $('completedBatchSampleFilter').innerHTML = '<option value="">全部149份样本文件</option>' + samples.map(row => `<option value="${esc(row.sample_key)}">${esc(sampleName(row))}</option>`).join('');
    renderSamples();
    const qcTotal = Number(totals.background_qc_rows || 0) + Number(totals.external_qc_rows || 0) + Number(totals.name_qc_rows || 0) + Number(totals.unlinked_qc_rows || 0);
    const fullRecordCounts = { xbt_jsonl: totals.xbt_named_candidate_rows, ms1: totals.ms1_features, qc: qcTotal };
    const nativeExport = {name: 'native_compounds', label: '全部原生结论与QC记录', rows: native.native_conclusion_rows,
      ...native.artifacts.native_compounds, download_url: native.native_compounds_xz_url, csv_url: native.native_compounds_csv_url};
    $('completedBatchExports').innerHTML = [nativeExport, ...(batch.exports || [])].map(row => {
      const filename = String(row.path || '').split('/').pop() || row.name;
      const count = row.rows ?? fullRecordCounts[row.name];
      return `<tr><td><b>${esc(exportNames[row.name] || row.label || row.name)}</b><small>${esc(filename)}</small></td><td>${count == null ? '—' : num(count)}</td><td>${bytes(row.size)}</td><td><div class="cb-action-row">${row.csv_url ? `<a class="btn btn-primary btn-sm" href="${esc(row.csv_url)}" download>CSV</a>` : ''}<a class="btn btn-secondary btn-sm" href="${esc(row.download_url)}" download>${filename.endsWith('.xz') ? '完整XZ' : 'JSON'}</a></div></td></tr>`;
    }).join('');
    showTab('overview');
  }
  function renderSamples() {
    const all = state.detail?.samples || [], search = ($('completedBatchSampleSearch')?.value || '').trim().toLowerCase();
    const rows = all.filter(row => sampleName(row).toLowerCase().includes(search));
    $('completedBatchSampleCount').textContent = `显示 ${num(rows.length)} / ${num(all.length)} 份文件`;
    $('completedBatchSamples').innerHTML = rows.map(row => `<tr${state.sample === row.sample_key ? ' class="is-selected"' : ''}><td>${esc(row.sample_id || row.sample_key)}</td><td>${esc(row.polarity)}</td><td>${esc(row.replicate)}</td><td>${num(row.searched_ms2_spectra)}</td><td>${num(row.matched_ms2_spectra)}</td><td>${num(row.ms1_features)}</td><td>${num(row.native_conclusion_records)}</td><td>${num(row.xbt_named_candidate_rows)}</td><td>${num(row.blank_flagged_ms2)}</td><td>${num(row.blank_clean_view_excluded_ms2)}</td><td><button type="button" class="btn btn-secondary btn-sm" data-cb-sample="${esc(row.sample_key)}">查看该文件表格</button><button type="button" class="btn btn-secondary btn-sm" onclick="Tasks.viewResults('${esc(row.task_id)}')">鉴定结果</button><a class="btn btn-secondary btn-sm" href="/api/tasks/${encodeURIComponent(row.task_id)}/report" target="_blank" rel="noopener">本文件报告</a></td></tr>`).join('');
  }
  function showTab(tab) {
    state.tab = tab;
    document.querySelectorAll('[data-cb-panel]').forEach(node => { node.hidden = node.dataset.cbPanel !== tab; });
    document.querySelectorAll('.cb-tabs [data-cb-tab]').forEach(node => { node.classList.toggle('active', node.dataset.cbTab === tab); node.setAttribute('aria-selected', String(node.dataset.cbTab === tab)); });
    if (tab === 'report' && state.detail && !$('completedBatchReportFrame').getAttribute('src')) $('completedBatchReportFrame').src = state.detail.native_delivery.report_url;
  }
  async function loadPreview() {
    if (!state.batch) return;
    const request = ++state.request;
    const params = new URLSearchParams({ offset: state.offset, limit: state.limit });
    if (state.sample) params.set('sample_key', state.sample);
    status('completedBatchPreviewStatus', '读取结果表…');
    $('completedBatchPreviewHead').innerHTML = '';
    $('completedBatchPreviewRows').innerHTML = '';
    $('completedBatchPreviewCount').textContent = '';
    $('completedBatchPrevious').disabled = true; $('completedBatchNext').disabled = true;
    const source = (state.detail?.exports || []).find(row => row.name === state.table);
    $('completedBatchTableCsv').href = state.table === 'native_compounds' ? state.detail.native_delivery.native_compounds_csv_url : source?.csv_url || source?.download_url || '#';
    try {
      const tableUrl = state.table === 'native_compounds' ? state.detail.native_delivery.native_table_url : `${base}/${encodeURIComponent(state.batch)}/tables/${encodeURIComponent(state.table)}`;
      const payload = await json(`${tableUrl}?${params}`);
      if (request !== state.request) return;
      const savedColumns = payload.columns || [];
      const leading = (priorityColumns[state.table] || []).filter(column => savedColumns.includes(column));
      const columns = leading.concat(savedColumns.filter(column => !leading.includes(column)));
      const rows = payload.rows || [], total = Number(payload.total || 0);
      $('completedBatchPreviewHead').innerHTML = `<tr><th>序号</th>${columns.map(column => `<th title="${esc(column)}">${esc(fieldNames[column] || column)}</th>`).join('')}</tr>`;
      $('completedBatchPreviewRows').innerHTML = rows.map((row, index) => `<tr><td>${num(state.offset + index + 1)}</td>${columns.map(column => {
        let value = row[column]; if (typeof value === 'boolean') value = value ? '是' : '否';
        if (value == null || value === '') value = '—';
        const text = String(value);
        return `<td><span class="${text.length > 100 ? 'cb-cell-long' : ''}" title="${esc(text)}">${esc(text)}</span></td>`;
      }).join('')}</tr>`).join('');
      $('completedBatchPreviewCount').textContent = rows.length ? `第 ${num(state.offset + 1)}–${num(state.offset + rows.length)} 条，共 ${num(total)} 条` : `当前筛选共 ${num(total)} 条`;
      $('completedBatchPrevious').disabled = state.offset === 0;
      $('completedBatchNext').disabled = payload.has_more == null ? state.offset + rows.length >= total : !payload.has_more;
      status('completedBatchPreviewStatus', rows.length ? '' : '该文件在当前结果表中没有记录，可切换 MS1、MS/MS 或 QC 表查看。');
    } catch (error) { if (request === state.request) status('completedBatchPreviewStatus', error.message, true); }
  }
  function init() {
    $('completedBatchSelect')?.addEventListener('change', event => selectBatch(event.target.value));
    $('completedBatchRefresh')?.addEventListener('click', () => enterPage(true));
    $('completedBatchSampleSearch')?.addEventListener('input', renderSamples);
    $('completedBatchTableSelect')?.addEventListener('change', event => { state.table = event.target.value; state.offset = 0; loadPreview(); });
    $('completedBatchSampleFilter')?.addEventListener('change', event => { state.sample = event.target.value; state.offset = 0; renderSamples(); loadPreview(); });
    $('completedBatchPrevious')?.addEventListener('click', () => { state.offset = Math.max(0, state.offset - state.limit); loadPreview(); });
    $('completedBatchNext')?.addEventListener('click', () => { state.offset += state.limit; loadPreview(); });
    document.querySelectorAll('[data-cb-tab]').forEach(node => node.addEventListener('click', () => showTab(node.dataset.cbTab)));
    $('completedBatchSamples')?.addEventListener('click', event => {
      const button = event.target.closest('[data-cb-sample]'); if (!button) return;
      state.sample = button.dataset.cbSample; state.offset = 0; $('completedBatchSampleFilter').value = state.sample;
      renderSamples(); loadPreview(); $('completedBatchPreviewSection')?.scrollIntoView({ behavior: 'smooth', block: 'start' });
    });
    loadList().catch(() => {});
  }
  return { init, enterPage };
})();
