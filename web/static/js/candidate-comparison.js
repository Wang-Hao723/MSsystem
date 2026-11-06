/** Complete-pool review. This component never mutates a saved candidate list. */
const CandidateComparison = (() => {
  const english = () => Boolean(window.MSSLanguage?.getLanguage?.()?.startsWith('en'));
  const text = (zh, en) => english() ? en : zh;
  const esc = value => String(value ?? '').replace(/[&<>"']/g, c => ({
    '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;',
  }[c]));
  const number = value => value == null || !Number.isFinite(Number(value)) ? '—' : Number(value).toFixed(4);

  function renderResult(payload) {
    const panel = payload.panel || {};
    const evidence = panel.comparison_evidence || {};
    const positions = new Map((payload.candidate_display_order || []).map((key, index) => [key, index]));
    const candidates = [...(evidence.candidates || [])].sort((a, b) =>
      (positions.get(a.connectivity) ?? Infinity) - (positions.get(b.connectivity) ?? Infinity)
      || a.connectivity.localeCompare(b.connectivity));
    const rows = candidates.map((c, index) => `<tr>
      <td>${index + 1}</td><td class="mono">${esc(c.connectivity)}</td>
      <td data-scientific-value>${esc((c.formulae || []).join(' / ') || '—')}</td>
      <td>${number(c.median_positive_fragment_contrast)}</td><td>${number(c.median_paired_cosine_delta)}</td>
      <td>${esc(c.compared_same_formula_candidates)}</td></tr>`).join('');
    const pairs = (evidence.pairs || []).filter(p => p.detailed_observations_included && p.status === 'compared');
    return `<p>${text(
      '按已保存的原生候选顺序、再按扩展召回顺序展示。差异碎片单独排序未通过回放检查，因此不替换正式 Top 1，不增加已鉴定化合物数量。碎片对比值是该候选相对同分子式竞争候选的实测差异离子支持中位数；共有离子不加分，未观察到的查询离子不扣分。数值不是概率。',
      'Saved native candidates are shown first, followed by production recall order. Fragment-only sorting failed replay checks and does not replace production Top 1 or add identified compounds. Fragment contrast is the median positive observed-ion support difference against same-formula competitors. Shared ions add zero; unobserved query ions incur no penalty. Scores are not probabilities.')}</p>
      <p>${text('固定候选', 'Frozen candidates')} ${esc(panel.candidate_count)} × ${text('查询谱', 'queries')} ${esc(panel.query_count)} × ${text('谱图视图', 'views')} ${esc(panel.view_count)} = ${esc(panel.expected_cells)} ${text('逐谱评分单元', 'scoring cells')}；
      ${text('已测量', 'measured')} ${esc(panel.cell_status_counts?.measured || 0)}；
      ${text('同条件候选对', 'condition-matched pairs')} ${esc(evidence.exact_condition_pair_count || 0)}；
      ${text('可比较 CE 轨迹的候选对', 'pairs with CE trajectories')} ${esc(evidence.trajectory_pair_count || 0)}。</p>
      <div class="table-wrap"><table><thead><tr>
      <th>${text('复核序号', 'Review order')}</th><th>${text('二维连接键', '2-D connectivity key')}</th>
      <th>${text('分子式', 'Formula')}</th><th>${text('碎片对比值', 'Fragment contrast')}</th>
      <th>${text('配对余弦分差', 'Paired cosine delta')}</th><th>${text('已比较竞争候选数', 'Compared competitors')}</th>
      </tr></thead><tbody>${rows}</tbody></table></div>
      <p>${text('— 表示未获得比较值，不是实测零分。同 CE 技术重复先取中位数；两个谱图视图不算两个独立实验。',
        '— means unavailable, not measured zero. Technical repeats are median-collapsed within CE; two spectral views are not two independent experiments.')}</p>
      ${pairs.length ? `<label>${text('逐对查看差异离子与 CE', 'Inspect paired ions and CE')} <select data-comparison-pair class="form-select">${pairs.map((p, i) => `<option value="${i}">${esc(p.a)} ↔ ${esc(p.b)}</option>`).join('')}</select></label><div data-comparison-detail></div>`
        : `<p>${text('没有可逐对比较的同分子式候选。缺分子式、参考谱或有效查询峰均不会被补成支持证据。',
          'No same-formula candidate pair is comparable. Missing formulae, references, or informative query peaks are not imputed as support.')}</p>`}
      <details><summary>${text('同条件参考谱排除记录', 'Conditioned-reference exclusion counts')}</summary>
      <pre style="white-space:pre-wrap;">${esc(JSON.stringify(evidence.strict_reference_rejection_tiers || {}, null, 2))}</pre></details>`;
  }

  function renderPair(pair, panel) {
    const conditions = new Map((panel.cells || []).map(c => [c.query_id, c.query_condition || {}]));
    const ions = values => (values || []).map(ion => `<details><summary>m/z ${Number(ion.query_mz).toFixed(5)} · ${number(ion.positive_support)}</summary>
      ${text('前体−产物 m/z 差（不额外计分）', 'Precursor–product m/z difference (not another vote)')} ${Number(ion.neutral_loss_da).toFixed(5)} Da
      <ul>${(ion.reference_observations || []).map(r => `<li class="mono">${esc(r.spectrum_id)} · Δ ${number(r.mass_error_da)} Da</li>`).join('')}</ul></details>`).join('') || '—';
    const observations = pair.observations || [];
    const rows = observations.slice(0, 40).map(o => {
      const condition = conditions.get(o.query_id) || {};
      return `<tr><td class="mono">${esc(o.query_id)}<br>${esc(o.view)}</td>
        <td>${esc(condition.collision_energy ?? '—')} ${esc(condition.collision_energy_unit || '')}</td>
        <td>${ions(o.a_observed_contrast_ions)}</td><td>${ions(o.b_observed_contrast_ions)}</td>
        <td>${esc(o.shared_observed_product_ions_excluded)}</td></tr>`;
    }).join('');
    const trajectories = pair.ce_trajectories || [];
    const ceRows = trajectories.flatMap(t => (t.edges || []).flatMap(edge => (edge.ions || []).map(ion => `<tr>
      <td>${esc(t.view)}</td><td>${esc(edge.low_ce)} → ${esc(edge.high_ce)} ${esc(t.query_context?.[4] || '')}</td>
      <td>${Number(ion.product_mz).toFixed(5)}</td><td>${number(ion.query_relative_intensity_change)}</td>
      <td>${number(ion.a_relative_intensity_change)}</td><td>${number(ion.b_relative_intensity_change)}</td></tr>`))).join('');
    return `<p><b>A:</b> <span class="mono">${esc(pair.a)}</span> · <b>B:</b> <span class="mono">${esc(pair.b)}</span><br>
      ${text('全条件碎片对比 A−B', 'All-condition fragment contrast A−B')}: ${number(pair.exploratory_fragment_delta_a_minus_b)}；
      ${text('全条件配对余弦 A−B', 'All-condition paired cosine A−B')}: ${number(pair.exploratory_cosine_delta_a_minus_b)}。
      ${text('正值表示 A 的实测支持更大，但全条件分差不能当作同条件证据。', 'Positive values indicate greater measured support for A; all-condition deltas are not condition-matched evidence.')}</p>
      <p>${text('逐谱明细', 'Per-query details')} ${Math.min(40, observations.length)} / ${observations.length}；
      ${text('明细行数不作为独立重复次数，完整明细可下载。某库谱未出现该离子不证明该化合物不会产生它。',
        'Detail rows are not independent repeats. Download JSON for all rows. An ion absent from a reference does not prove that its compound cannot produce it.')}</p>
      <div class="table-wrap"><table><thead><tr><th>${text('查询谱 / 视图', 'Query / view')}</th><th>CE</th>
      <th>${text('A 实测差异离子', 'A observed contrast ions')}</th><th>${text('B 实测差异离子', 'B observed contrast ions')}</th>
      <th>${text('不计分的共有离子', 'Shared ions, not scored')}</th></tr></thead><tbody>${rows}</tbody></table></div>
      <h5>${text('同条件 CE 强度轨迹', 'Condition-matched CE intensity trajectories')}</h5>
      ${ceRows ? `<p>${text('仅比较各能量下均实测到的离子；数值为该视图总离子强度归一化后的变化。',
        'Only ions observed at both energies are compared. Values are changes after normalization to total ion intensity within the view.')}</p>
        <div class="table-wrap"><table><thead><tr><th>${text('视图', 'View')}</th><th>CE</th><th>m/z</th><th>${text('查询 Δ', 'Query Δ')}</th><th>A Δ</th><th>B Δ</th></tr></thead><tbody>${ceRows}</tbody></table></div>`
        : `<p>${text('暂无可比较轨迹：必须有同仪器域、活化方式、离子形式、CE 单位和固定能量模式的参考谱，且至少两个相同 CE 点各有两条共同实测离子。缺字段、跨单位或单一 CE 不会被补成轨迹。',
          'No comparable trajectory. References must match platform, activation, ion form, CE unit and fixed-energy mode, with at least two identical CE points and two jointly observed ions. Missing metadata, mixed units or a single CE are not converted into a trajectory.')}</p>`}`;
  }

  function mount(container, taskId, entity) {
    if (!container || !entity?.entity_id) return;
    container.innerHTML = `<h4>${text('V16 · Top40 差异碎片与 CE 比较', 'V16 · Top40 fragment and CE comparison')}</h4>
      <p>${text('按需比较当前实体的完整候选池，不重新转换数据，不改写已有样本结果。',
        'Compare the complete candidate pool on demand, without reconversion or changes to saved sample results.')}</p>
      <button type="button" class="btn btn-primary" data-comparison-run>${text('计算完整候选比较', 'Compare complete candidate pool')}</button>
      <span role="status" aria-live="polite" data-comparison-status></span><div data-comparison-result></div>`;
    const button = container.querySelector('[data-comparison-run]');
    const status = container.querySelector('[data-comparison-status]');
    const result = container.querySelector('[data-comparison-result]');
    button.addEventListener('click', async () => {
      button.disabled = true;
      status.textContent = text(' 正在逐谱比较…', ' Comparing saved spectra…');
      try {
        const payload = await API.post(`/tasks/${encodeURIComponent(taskId)}/candidate-comparison`, {entity_id: entity.entity_id});
        if (!container.isConnected) return;
        result.innerHTML = renderResult(payload);
        const panel = payload.panel;
        const pairs = panel.comparison_evidence.pairs.filter(p => p.detailed_observations_included && p.status === 'compared');
        const select = result.querySelector('[data-comparison-pair]');
        const detail = result.querySelector('[data-comparison-detail]');
        if (select && detail) {
          const showPair = () => { detail.innerHTML = renderPair(pairs[Number(select.value)], panel); };
          select.addEventListener('change', showPair);
          showPair();
        }
        const download = document.createElement('a');
        download.className = 'btn btn-secondary';
        download.textContent = text('下载完整比较 JSON', 'Download complete comparison JSON');
        const params = new URLSearchParams({entity_id: entity.entity_id, comparison_id: payload.comparison_id});
        download.href = `${API.BASE}/tasks/${encodeURIComponent(taskId)}/candidate-comparison/export?${params.toString()}`;
        download.download = `${taskId}_v16_candidate_comparison.json`;
        result.appendChild(download);
        status.textContent = text(' 已完成；正式排序未改动。', ' Complete; production ranking unchanged.');
      } catch (error) {
        if (container.isConnected) status.textContent = `${text(' 比较未完成：', ' Comparison unavailable: ')}${error.message}`;
      } finally {
        button.disabled = false;
      }
    });
  }
  return {mount, renderResult, renderPair};
})();
