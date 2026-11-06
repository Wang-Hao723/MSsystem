/* Saved XBT results inside the existing task, identification and report pages. */
const FrozenTaskUI = (() => {
  const esc = value => Tasks.escapeHtml(String(value ?? ''));
  const num = value => Number(value || 0).toLocaleString();
  const tables = [
    ['native_compounds', '全部原生结论与QC'], ['xbt_csv', 'XBT联合筛选（附加）'], ['ms1_table', '全部 MS1 特征（含仅 MS1 特征）'],
    ['msms_table', '全部 MS/MS 前体事件'], ['qc_table', '背景及候选 QC'],
  ];
  function panel(task) {
    const batch = task?.completed_batch;
    if (!batch) return '';
    const summary = batch.summary || {};
    const taskId = encodeURIComponent(task.id);
    return `<section class="card" style="border-left:4px solid #246c71;padding:18px;overflow-wrap:anywhere" data-frozen-task="${esc(task.id)}">
      <div class="card-title">宣痹汤 / XBT · 完整原生MS/MS与独立MS1证据</div>
      <p>${esc(task.filename)}：${num(summary.parsed_searched_ms2_spectra)} 个 MS/MS 前体事件，${num(summary.matched_ms2_spectra)} 张MS/MS有谱库候选，${num(summary.ms1_feature_records)} 个 MS1 靶向特征。原生结论记录共 ${num(summary.compound_result_records ?? summary.parsed_searched_ms2_spectra)} 条，含候选、背景、未解析和未匹配记录。</p>
      <p>XBT与MS1联合筛选 ${num(summary.xbt_named_candidate_rows)} 条，作为附加视图；各条仍按原生记录类型区分候选与未解析事件。配对空白标记 ${num(summary.blank_flagged_ms2)} 个事件；净化视图排除 ${num(summary.blank_clean_view_excluded_ms2)} 个。标准品确证 0。</p>
      <div style="display:flex;gap:8px;flex-wrap:wrap;margin:10px 0">
        <a class="btn btn-primary btn-sm" href="/api/tasks/${taskId}/report" target="_blank" rel="noopener">当前样本全量证据报告</a>
        <a class="btn btn-primary btn-sm" href="/api/tasks/${taskId}/download/html" download>下载当前样本全量报告</a>
        <a class="btn btn-secondary btn-sm" href="/api/tasks/${taskId}/download/native-conclusions.json.xz" download>完整原生结论与证据JSON.XZ</a>
        <a class="btn btn-secondary btn-sm" href="/api/completed-batches/${encodeURIComponent(batch.id)}/native-download-all" download>149份完整原生结论、MS1/MSMS与QC结果包</a>
        ${tables.map(([key, label]) => `<a class="btn btn-secondary btn-sm" href="/api/tasks/${taskId}/frozen-tables/${key}">${esc(label)}：本文件 CSV</a>`).join('')}
      </div>
      <details><summary>查看本文件完整结果表</summary><div style="margin:12px 0;display:flex;gap:8px;flex-wrap:wrap">${tables.map(([key, label]) => `<button class="btn btn-secondary btn-sm" onclick="FrozenTaskUI.preview(this,'${esc(task.id)}','${key}',0)">${esc(label)}</button>`).join('')}</div><div data-frozen-preview></div></details>
      <small>全量原生结果直接读取已保存结论，不以XBT关联或MS1关联作为展示全部MS/MS的前提。记录出现次数、候选身份、未解析事件和标准品确证分别记录。<a href="${esc(batch.report_url)}" target="_blank" rel="noopener">既有XBT联合筛选摘要</a></small>
    </section>`;
  }
  async function preview(button, taskId, table, offset) {
    const root = button.closest('[data-frozen-task]');
    const area = root?.querySelector('[data-frozen-preview]');
    if (!area) return;
    area.textContent = '读取已保存结果…';
    try {
      const task = await API.get(`/tasks/${encodeURIComponent(taskId)}/status`);
      const batch = task.completed_batch;
      const params = new URLSearchParams({sample_key: batch.sample_key, limit: 50, offset});
      let data;
      if (table === 'native_compounds') {
        const result = await API.get(`/tasks/${encodeURIComponent(taskId)}/compounds?result_view=native&limit=50&offset=${offset}`);
        const columns = ['原生索引', '记录类型', '原生状态', '谱库假设标题（QC）', '分子式', '代表谱图', '前体m/z', '空白标记'];
        data = {columns, rows: result.compounds.map(row => Object.fromEntries(columns.map((key, index) => [key, [row.frozen_compound_index, row.record_type, row.conclusion?.key, row.lead_candidate?.compound_name || row.display_name, row.formula, row.representative?.spectrum_id, row.representative?.precursor_mz, row.representative?.blank_flagged][index]]))), total: result.filtered_total, has_more: result.has_more};
      } else data = await API.get(`/completed-batches/${encodeURIComponent(batch.id)}/tables/${table}?${params}`);
      area.innerHTML = `<p>第 ${num(offset + (data.rows.length ? 1 : 0))}–${num(offset + data.rows.length)} 条，共 ${num(data.total)} 条</p>
        <div class="table-wrap"><table><thead><tr>${data.columns.map(key => `<th>${esc(key)}</th>`).join('')}</tr></thead><tbody>${data.rows.map(row => `<tr>${data.columns.map(key => `<td>${esc(row[key] == null || row[key] === '' ? '—' : row[key])}</td>`).join('')}</tr>`).join('')}</tbody></table></div>
        <div class="btn-group" style="margin:10px 0"><button class="btn btn-secondary btn-sm" ${offset ? '' : 'disabled'} onclick="FrozenTaskUI.preview(this,'${esc(taskId)}','${table}',${Math.max(0, offset - 50)})">上一页</button><button class="btn btn-secondary btn-sm" ${data.has_more ? '' : 'disabled'} onclick="FrozenTaskUI.preview(this,'${esc(taskId)}','${table}',${offset + 50})">下一页</button></div>`;
    } catch (error) { area.textContent = error.message; }
  }
  function download(url) {
    const anchor = document.createElement('a');
    anchor.href = url; anchor.download = '';
    document.body.appendChild(anchor); anchor.click(); anchor.remove();
  }
  return {panel, preview, download};
})();
