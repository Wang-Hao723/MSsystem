/** Read-only per-task module effect inspector. */
const ModuleEffects = (() => {
  const esc = value => Tasks.escapeHtml(String(value == null ? '' : value));
  const labels = {
    production: '当前鉴定流程', qualified: '条件性启用', shadow: '研究性建议', operational: '运行记录',
    operational_audit: '质量审计', independent_evidence_layer: '独立证据层',
  };
  const metricLabels = {
    matched_ms2: '获得候选的 MS2', total_ms2: 'MS2 总数', unique_library_candidates: '去重谱库候选',
    high_supported_identities: '较强证据身份', unmatched_with_fragments: '有碎片但未匹配',
    ms1_scans: 'MS1 扫描', lcms_features_with_ms2: '有 MS2 的 LC–MS 特征', linked_ms2_spectra: '已归属 MS2',
    accepted_spectra: '通过门控的谱图', high_supported_identity_delta: '较强证据身份净增',
    representative_high_delta: '代表谱 High 净增', high_spectra_delta: 'High 谱图净增', top1_changed: '首选候选改变',
    networks: '离子身份网络', networked_features: '联网特征', linked_annotations: '链接注释',
    detected_features: '检出的 XIC 峰', supported_features: '有正交支持的峰', mass_only_features: '仅精确质量的峰',
    ordinary_ms2_match_increment: '普通 MS2 命中增量', compiled: '已评估候选决策', released: '可作结论',
    calibrated: '有项目级校准', production_ranking_changed_records: '候选排序变化记录',
    identification_records: '鉴定记录', module_state_counts: '模块状态数', followups_linked: '已链接复采结果',
    plan_count: '确认计划', planned_actions: '备选采集方案',
    spectra_processed: '已构图谱图', formula_context_available: '可用分子式上下文', nodes: '图节点',
    truncated_by_task_limit: '达到影子任务上限',
    selected_edges: '选中边', heuristically_supported_edges: '启发式支持边', ambiguous_edges: '歧义边',
    unassigned_edges: '未归属边', cycle_count: '环路违规', high_promotions: 'High 提升',
  };
  const rankingLabels = {
    production_ranker: '产生候选排序', none: '不改排序', forbidden: '禁止改排序',
    selected_identity_preserved: '仅在已选身份内增分',
    separate_from_ordinary_ms2: '独立统计',
  };

  function userFacingEffect(effect) {
    if (effect.module_id === 'evigate_v2_claim_compiler') {
      const assessed = Number(effect.metrics?.compiled || 0);
      const eligible = Number(effect.metrics?.released || 0);
      return {
        ...effect,
        display_name: '鉴定结论条件评估',
        headline: `已评估 ${assessed} 条候选，${eligible} 条满足当前正式报告条件`,
        interpretation: '正式报告条件与候选鉴定分开评估；条件不足不会改写当前首选候选、原证据等级或候选鉴定范围。',
      };
    }
    if (effect.module_id === 'shadow_evidence_modules') {
      const states = Object.keys(effect.metrics?.module_state_counts || {}).length;
      return {
        ...effect,
        display_name: '补充证据与一致性检查',
        headline: `${states} 类补充检查留下逐记录状态`,
        interpretation: '碰撞能、规则证据与结构一致性等只用于显示和审计；按设计不改变当前候选排序。',
      };
    }
    return effect;
  }

  function metricRows(metrics) {
    return Object.entries(metrics || {}).filter(([, value]) => (
      ['string', 'number', 'boolean'].includes(typeof value)
    )).slice(0, 6).map(([key, value]) => {
      const label = metricLabels[key] || key.replaceAll('_', ' ');
      const display = typeof value === 'boolean' ? (value ? '是' : '否') : value;
      return `<span>${esc(label)} <b>${esc(display)}</b></span>`;
    }).join('');
  }

  function ratio(value, total) {
    const numerator = Number(value || 0);
    const denominator = Number(total || 0);
    return denominator > 0 ? Math.max(0, Math.min(100, numerator / denominator * 100)) : 0;
  }

  function meterRow(label, value, total, tone, note) {
    const percent = ratio(value, total);
    return `<div class="quality-meter-row meter-${esc(tone)}">
      <div><b>${esc(label)}</b><small>${esc(note || '')}</small></div>
      <div class="quality-meter-track" aria-label="${esc(label)} ${percent.toFixed(1)}%"><i style="width:${percent.toFixed(1)}%"></i></div>
      <strong>${Number(value || 0).toLocaleString()}<span> / ${Number(total || 0).toLocaleString()}</span></strong>
    </div>`;
  }

  function moduleTone(effect) {
    if (!effect.executed) return effect.applicable ? 'attention' : 'muted';
    if (effect.state === 'production') return 'core';
    if (effect.state === 'qualified') return 'qualified';
    if (effect.state === 'shadow') return 'research';
    return 'support';
  }

  function moduleRow(rawEffect, index) {
    const effect = userFacingEffect(rawEffect);
    const state = String(effect.state || 'unknown');
    const executed = !!effect.executed;
    const applicable = !!effect.applicable;
    const status = executed ? '已运行' : applicable ? '适用但无输出' : '本任务不适用';
    const metrics = metricRows(effect.metrics);
    return `<li class="quality-module-row quality-tone-${moduleTone(effect)}">
      <span class="quality-module-index">${String(index + 1).padStart(2, '0')}</span>
      <div class="quality-module-main">
        <div class="quality-module-title"><b>${esc(effect.display_name)}</b><span>${esc(labels[state] || state)}</span></div>
        <strong>${esc(status)} · ${esc(effect.headline)}</strong>
        <details class="quality-module-detail"><summary>查看指标与适用边界</summary><div><div class="quality-module-metrics">${metrics || '<span>无附加数值</span>'}</div><p>${esc(effect.interpretation)}</p></div></details>
      </div>
      <div class="quality-module-impact"><span>排序关系</span><b>${esc(rankingLabels[effect.ranking_effect] || effect.ranking_effect)}</b></div>
    </li>`;
  }

  async function render(taskId, containerId = 'moduleEffectPanel') {
    const container = document.getElementById(containerId);
    if (!container) return;
    if (!taskId) {
      container.innerHTML = '';
      return;
    }
    container.innerHTML = '<div class="analysis-panel-loading">正在核对各模块的执行状态与适用边界…</div>';
    try {
      const payload = await API.get(`/tasks/${encodeURIComponent(taskId)}/module-effects`);
      const changed = !!payload.core_ranking_changed;
      const effects = Array.isArray(payload.effects) ? payload.effects : [];
      const byId = moduleId => effects.find(effect => effect.module_id === moduleId) || {};
      const core = byId('core_spectral_ranker').metrics || {};
      const ms1 = byId('ms1_feature_detection').metrics || {};
      const expected = byId('expected_compounds').metrics || {};
      const deconvolution = byId('chimeric_deconvolution').metrics || {};
      const compiler = byId('evigate_v2_claim_compiler').metrics || {};
      const executedCount = effects.filter(effect => effect.executed).length;
      const visibleOutputCount = Number(payload.observable_nonranking_effects?.length || 0);
      const completion = effects.length ? executedCount / effects.length * 100 : 0;
      const moduleRows = effects.map(moduleRow).join('');
      const reasons = (payload.why_core_looks_unchanged || []).map(
        reason => `<li>${esc(reason)}</li>`
      ).join('');
      container.innerHTML = `<section class="analysis-quality-dashboard">
        <header class="analysis-dashboard-header"><div><span class="panel-eyebrow">QUALITY &amp; APPLICABILITY</span><h3>分析质量与适用性</h3><p>先看覆盖、风险和排序关系；每个模块的技术说明默认收起。</p></div><span class="quality-verdict ${changed ? 'verdict-alert' : 'verdict-stable'}">${changed ? '检测到排序变化' : '候选排序保持稳定'}</span></header>
        <div class="quality-overview">
          <div class="quality-completion-ring" style="--quality-progress:${(completion * 3.6).toFixed(1)}deg"><div><b>${executedCount}/${effects.length}</b><span>模块已运行</span></div></div>
          <div class="quality-overview-copy"><b>${visibleOutputCount} 个辅助模块形成可查看结果</b><p>模块运行不等于增加化合物身份；独立证据层和研究性建议按设计不直接改写候选排序。</p><div class="quality-overview-facts"><span><b>${Number(core.unique_library_candidates || 0).toLocaleString()}</b> 去重候选</span><span><b>${Number(core.high_supported_identities || 0).toLocaleString()}</b> 较强证据身份</span><span><b>${Number(compiler.calibrated || 0).toLocaleString()}</b> 项目级校准</span></div></div>
        </div>
        <section class="quality-signal-board"><div class="quality-section-heading"><div><span>KEY SIGNALS</span><h4>关键覆盖与证据边界</h4></div><small>横条只表示各自分母内的覆盖，不是鉴定概率</small></div>
          ${meterRow('获得谱库候选的 MS2', core.matched_ms2, core.total_ms2, 'teal', '当前谱库与检索阈值')}
          ${meterRow('归属到连续 MS1 特征的 MS2', ms1.linked_ms2_spectra, core.total_ms2, 'blue', '提供色谱上下文，不新增身份')}
          ${meterRow('有碎片但未获得谱库候选', core.unmatched_with_fragments, core.total_ms2, 'amber', '仍保留在未匹配谱图工作台')}
          ${meterRow('Expected Compounds 获得正交支持', expected.supported_features, expected.detected_features, 'violet', '其余仅为精确质量线索')}
        </section>
        <section class="quality-watchlist"><div><span class="watch-good">+${Number(deconvolution.high_supported_identity_delta || 0)}</span><p><b>去卷积带来的较强身份净增</b><small>Top-1 保持不变</small></p></div><div><span class="watch-attention">${Number(core.unmatched_with_fragments || 0).toLocaleString()}</span><p><b>有碎片但未匹配</b><small>优先查看未匹配谱图</small></p></div><div><span class="watch-attention">${Number(expected.mass_only_features || 0).toLocaleString()}</span><p><b>仅精确质量的 XIC 峰</b><small>不能当作 MS/MS 鉴定</small></p></div><div><span class="watch-good">${changed ? 1 : 0}</span><p><b>当前候选排序变化</b><small>${changed ? '需要复核变化来源' : '保持原排序'}</small></p></div></section>
        <section class="quality-module-section"><div class="quality-section-heading"><div><span>MODULE TRACE</span><h4>模块执行轨迹</h4></div><small>逐行查看，不再铺成并列卡片</small></div><ol class="quality-module-timeline">${moduleRows}</ol></section>
        <details class="quality-governance-note"><summary>为什么部分模块运行后候选名称仍不变？</summary><ol>${reasons}</ol></details>
      </section>`;
    } catch (error) {
      container.innerHTML = `<div class="analysis-panel-error">模块效果面板加载失败：${esc(error.message || error)}</div>`;
    }
  }

  return { render };
})();
