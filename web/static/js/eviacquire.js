/** EviAcquire-MS shadow active-evidence planning panel. */
const EviAcquire = (() => {
  let currentTaskId = null;
  let currentPayload = null;
  let currentContainerId = 'eviAcquirePanel';
  let showAllPlans = false;
  const INITIAL_PLAN_LIMIT = 8;
  const esc = value => Tasks.escapeHtml(String(value == null ? '' : value));
  const num = (value, digits = 2) => {
    const parsed = Number(value);
    return Number.isFinite(parsed) ? parsed.toFixed(digits) : '—';
  };
  const missingEvidenceLabels = {
    candidate_independent_isotope_support: '候选无关的同位素支持',
    continuous_ms1_feature_link: '连续 MS1 色谱特征链接',
    energy_resolved_fragmentation: '碰撞能分辨碎裂证据',
    independent_injection_replicate: '独立进样重复',
    narrow_isolation_confirmation: '窄隔离窗确认',
    same_method_retention_time_or_authentic_standard: '同方法 RT 证据缺失（未登记标品不生成动作）',
  };
  const feasibilityLabels = {
    explicit_authentic_standard_inventory_confirmed: '候选标品与本次使用预算均已登记',
    explicit_ion_mobility_and_ccs_reference_confirmed: '离子淌度能力与候选 CCS 参比均已登记',
    verify_precursor_density_and_instrument_limit: '需核对前体密度与仪器隔离能力',
    instrument_method_review_required: '需由仪器负责人审核方法',
    run_after_energy_scout_then_check_peak_width_and_cycle_time: '先完成能量侦察，再核对峰宽与循环时间',
    requires_opposite_polarity_detectability: '需确认反极性下具有可检测响应',
  };

  function energyDesignSummary(action) {
    const design = action.collision_energy_design || {};
    const stages = Array.isArray(design.stages) ? design.stages : [];
    const conditional = stages.slice(1).filter(stage =>
      Array.isArray(stage.collision_energies) && stage.collision_energies.length
    ).map(stage => `${stage.stage_id === 'labile_low_energy_refinement' ? '易碎补低能' : stage.stage_id === 'transition_midpoint_refinement' ? '转折补中能' : stage.stage_id === 'resistant_high_energy_extension' ? '耐碎补高能' : stage.stage_id} ${stage.collision_energies.map(value => num(value, 0)).join('/')} ${design.collision_energy_unit || ''}`);
    if (!conditional.length) return '';
    return `条件分支：${conditional.join('；')}（不是全部同时采集）`;
  }

  function actionSummary(action) {
    const ce = Array.isArray(action.collision_energies) && action.collision_energies.length
      ? `首轮 CE ${action.collision_energies.map(value => num(value, 0)).join('/')}${action.collision_energy_unit ? ' ' + action.collision_energy_unit : ''}`
      : action.action_kind === 'high_sensitivity_targeted_ms2' ? '能量侦察后选择 CE' : '无 CE 参数';
    const staged = energyDesignSummary(action);
    const isolation = action.isolation_width_da != null
      ? `隔离窗 ${num(action.isolation_width_da, 2)} Da` : '';
    const injection = action.max_injection_time_ms != null
      ? `最大注入 ${num(action.max_injection_time_ms, 0)} ms` : '';
    return [ce, staged, isolation, injection].filter(Boolean).join(' · ');
  }

  function candidatePair(plan) {
    const need = plan.need || {};
    const candidates = (need.candidates || []).map(item => item.compound_name).filter(Boolean);
    return candidates.slice(0, 3).join(' ↔ ') + (candidates.length > 3 ? ` (+${candidates.length - 3})` : '');
  }

  function actionParameterChips(action) {
    const chips = [];
    if (Array.isArray(action.collision_energies) && action.collision_energies.length) {
      chips.push(`CE ${action.collision_energies.map(value => num(value, 0)).join('/')} ${action.collision_energy_unit || ''}`.trim());
    }
    if (action.isolation_width_da != null) chips.push(`隔离窗 ${num(action.isolation_width_da, 2)} Da`);
    if (action.max_injection_time_ms != null) chips.push(`注入 ${num(action.max_injection_time_ms, 0)} ms`);
    return chips;
  }

  function planRow(plan, index, maximumUtility) {
    const need = plan.need || {};
    const actions = plan.ranked_actions || [];
    const top = actions[0] || {};
    const interval = top.utility_interval || [0, 0];
    const pair = candidatePair(plan);
    const priorityWidth = maximumUtility > 0
      ? Math.max(8, Math.min(100, Number(top.utility || 0) / maximumUtility * 100)) : 0;
    const parameters = actionParameterChips(top).map(item => `<span>${esc(item)}</span>`).join('');
    const missing = (need.missing_evidence || []).slice(0, 4).map(
      item => `<span>${esc(missingEvidenceLabels[item] || item.replaceAll('_', ' '))}</span>`
    ).join('');
    return `<article class="acquisition-plan-row">
      <div class="acquisition-plan-rank"><span>${String(index + 1).padStart(2, '0')}</span><i></i></div>
      <div class="acquisition-plan-main">
        <header><div><h4>${esc(pair || '未命名候选对')}</h4><p>m/z ${num(need.precursor_mz, 5)} · RT ${num(need.retention_time_min, 3)} min · ${esc(need.polarity || '极性未知')}</p></div><span>${Number((need.ambiguity_edges || []).length)} 条待区分关系</span></header>
        <div class="acquisition-primary-action"><b>${esc(top.label || '暂无可执行方案')}</b><div>${parameters || '<span>参数待方法负责人确认</span>'}</div></div>
        <div class="acquisition-priority"><span>计划优先级</span><div><i style="width:${priorityWidth.toFixed(1)}%"></i></div><b>${num(top.utility)}</b><small>预期区分 ${num(top.minimum_pair_separation)}</small></div>
        <details class="acquisition-plan-detail"><summary>查看依据、缺失证据与导出</summary><div>
          <div class="acquisition-missing-evidence">${missing || '<span>未记录额外缺失证据</span>'}</div>
          <p><b>参数摘要：</b>${esc(actionSummary(top))}</p>
          <p><b>依据：</b>${esc((top.rationale || []).slice(0, 3).join('；') || '未记录')}</p>
          <p class="acquisition-uncertainty">方案优先级 ${num(top.utility)} [${num(interval[0])}, ${num(interval[1])}]；该数值不是身份概率，预测信息不是已经观测到的证据。</p>
          <div class="acquisition-plan-actions"><button class="btn btn-primary btn-sm" onclick="EviAcquire.showPlan('${esc(plan.plan_id)}')">查看完整采集方案</button><button class="btn btn-secondary btn-sm" onclick="EviAcquire.download('${esc(plan.plan_id)}','vendor-neutral')">通用 CSV</button><button class="btn btn-secondary btn-sm" onclick="EviAcquire.download('${esc(plan.plan_id)}','thermo-qe')">Q Exactive</button><button class="btn btn-secondary btn-sm" onclick="EviAcquire.download('${esc(plan.plan_id)}','agilent-qtof')">Agilent Q-TOF</button></div>
        </div></details>
      </div>
    </article>`;
  }

  function renderAvailable(payload, container) {
    const plans = Array.isArray(payload.plans) ? payload.plans : [];
    const metadata = payload.metadata || {};
    const resources = metadata.resource_gate || {};
    const allActions = plans.flatMap(plan => plan.ranked_actions || []);
    const energyAction = allActions.find(
      action => action.collision_energy_design?.profile_id === 'agilent_qtof_ev_sequential_v1.0.0'
    );
    const agilentNarrowAction = allActions.find(
      action => action.action_kind === 'narrow_isolation_ms2'
        && action.isolation_width_profile_id === 'agilent_qtof_native_isolation_v1.0.0'
    );
    const maximumUtility = Math.max(0, ...plans.map(plan => Number(plan.ranked_actions?.[0]?.utility || 0)));
    const visiblePlans = showAllPlans ? plans : plans.slice(0, INITIAL_PLAN_LIMIT);
    const boundaryRows = [];
    if (resources.mode === 'sample_only_fail_closed') {
      boundaryRows.push('<div><span>资源</span><p><b>本次按无标品资源模式生成</b>：未登记的标品 RT 与离子淌度/CCS 已从动作空间、排序和导出中移除。只有任务明确登记候选标品或离子淌度/CCS 参比并给出预算后，相应动作才会进入。</p></div>');
    } else if (resources.sample_only_primary_workflow !== false) {
      boundaryRows.push(`<div><span>资源</span><p><b>已读取资源清单，但本批没有资源动作进入计划</b>：标品 ${Number(resources.authentic_standard_inventory_count || 0)} 个，CCS 参比 ${Number(resources.ccs_reference_count || 0)} 个；当前继续按原样本复采。</p></div>`);
    } else {
      boundaryRows.push(`<div><span>资源</span><p><b>本次使用显式资源清单</b>：标品动作 ${Number(resources.authentic_standard_actions_planned || 0)} 项，离子淌度动作 ${Number(resources.ion_mobility_actions_planned || 0)} 项；未登记身份仍按不可用处理。</p></div>`);
    }
    if (energyAction) {
      boundaryRows.push('<div><span>碰撞能</span><p><b>首轮 10/20/40 eV</b>；只有 40 eV 后仍欠碎裂才补 60 eV，易过碎裂或响应转折时分别补 5 或 30 eV。不声称所有化合物“完全碎裂”；RT 标品证据与碰撞能优化独立。</p></div>');
    }
    if (agilentNarrowAction) {
      boundaryRows.push('<div><span>隔离窗</span><p><b>原方法 Medium 约 4 amu</b>；窄窗复采使用仪器原生 Narrow 约 1.3 amu，不输出不属于该平台方法域的 0.4 Da。下发前仍需仪器负责人审核。</p></div>');
    }
    const scout = Array.isArray(energyAction?.collision_energies) && energyAction.collision_energies.length
      ? energyAction.collision_energies.map(value => num(value, 0)).join('/') + ' eV' : '按实测响应';
    const isolation = agilentNarrowAction?.isolation_width_da != null
      ? num(agilentNarrowAction.isolation_width_da, 2) + ' Da' : '按仪器域确认';
    container.innerHTML = `<section class="acquisition-dashboard">
      <header class="analysis-dashboard-header"><div><span class="panel-eyebrow">NEXT MEASUREMENT</span><h3>确认实验行动路线</h3><p>首屏只显示下一步做什么和关键参数；依据、边界与导出按需展开。</p></div><span class="acquisition-mode">原样本复采 · ${plans.length} 组</span></header>
      <div class="acquisition-overview"><div><b>${plans.length}</b><span>优先计划</span></div><div><b>${Number(metadata.ambiguity_edge_count || 0)}</b><span>待区分关系</span></div><div><b>${Number(metadata.planned_actions || 0)}</b><span>备选方案</span></div><p><b>当前鉴定排序未改变；预测信息不作为确认性证据。</b> 计划用于仪器负责人复核，不是可直接下发的方法文件。</p></div>
      <div class="acquisition-route" aria-label="确认实验执行路线"><div><span>01</span><b>锁定歧义事件</b><small>m/z + RT + 极性</small></div><i></i><div><span>02</span><b>降低共隔离</b><small>${esc(isolation)}</small></div><i></i><div><span>03</span><b>能量侦察</b><small>${esc(scout)}</small></div><i></i><div><span>04</span><b>按响应补采</b><small>条件 5 / 30 / 60 eV</small></div></div>
      <details class="acquisition-method-boundaries"><summary>方法边界与资源条件 <span>${boundaryRows.length}</span></summary><div>${boundaryRows.join('')}</div></details>
      <section class="acquisition-queue"><div class="quality-section-heading"><div><span>PRIORITY QUEUE</span><h4>优先行动队列</h4></div><small>当前显示 ${visiblePlans.length} / ${plans.length} 组</small></div><div class="acquisition-plan-list">${visiblePlans.map((plan, index) => planRow(plan, index, maximumUtility)).join('')}</div></section>
      <div class="acquisition-footer-actions">${plans.length > INITIAL_PLAN_LIMIT ? `<button class="btn btn-secondary" onclick="EviAcquire.toggleAll()">${showAllPlans ? `收起到前 ${INITIAL_PLAN_LIMIT} 组` : `显示全部 ${plans.length} 组`}</button>` : ''}<button class="btn btn-secondary" onclick="EviAcquire.generate(true)">按当前策略重新生成</button>${payload.has_more ? '<span>服务端仍有更多计划未加载</span>' : ''}</div>
    </section>`;
  }

  async function render(taskId, containerId = 'eviAcquirePanel') {
    if (currentTaskId !== (taskId || null)) showAllPlans = false;
    currentTaskId = taskId || null;
    currentContainerId = containerId;
    const container = document.getElementById(containerId);
    if (!container) return;
    if (!taskId) {
      container.innerHTML = '';
      return;
    }
    container.innerHTML = '<div class="analysis-panel-loading">正在读取确认实验计划…</div>';
    try {
      const payload = await API.get(`/tasks/${encodeURIComponent(taskId)}/acquisition-plans?limit=50`);
      currentPayload = payload;
      const plans = Array.isArray(payload.plans) ? payload.plans : [];
      if (!plans.length) {
        container.innerHTML = `<section class="acquisition-dashboard acquisition-empty"><header class="analysis-dashboard-header"><div><span class="panel-eyebrow">NEXT MEASUREMENT</span><h3>尚未生成确认实验路线</h3><p>默认按无标品、无离子淌度资源规划原样本复采；只有任务明确登记候选标品或离子淌度/CCS 参比后，相应动作才会进入。</p></div></header><button id="eviAcquireGenerate" class="btn btn-primary" onclick="EviAcquire.generate()">生成确认实验建议</button></section>`;
        return;
      }
      renderAvailable(payload, container);
    } catch (error) {
      container.innerHTML = `<div class="analysis-panel-error">确认实验建议加载失败：${esc(error.message || error)}</div>`;
    }
  }

  function toggleAll() {
    if (!currentPayload) return;
    const container = document.getElementById(currentContainerId);
    if (!container) return;
    showAllPlans = !showAllPlans;
    renderAvailable(currentPayload, container);
    container.scrollIntoView({ behavior: 'smooth', block: 'start' });
  }

  async function generate(refresh = false) {
    if (!currentTaskId) return;
    const button = document.getElementById('eviAcquireGenerate');
    if (button) {
      button.disabled = true;
      button.textContent = '正在构建歧义图并评分…';
    }
    try {
      await API.post(`/tasks/${encodeURIComponent(currentTaskId)}/acquisition-plans`, {
        max_plans: 25,
        max_actions_per_plan: 6,
        refresh: !!refresh,
      });
      Toast.show('确认实验建议已生成（不会改变当前鉴定排序）', 'success');
      await render(currentTaskId);
      if (typeof ModuleEffects !== 'undefined') ModuleEffects.render(currentTaskId);
    } catch (error) {
      Toast.show(`EviAcquire-MS 生成失败：${error.message || error}`, 'error');
      if (button) {
        button.disabled = false;
        button.textContent = '生成确认实验建议';
      }
    }
  }

  async function showPlan(planId) {
    try {
      const plan = await API.get(`/tasks/${encodeURIComponent(currentTaskId)}/acquisition-plans/${encodeURIComponent(planId)}`);
      const actions = (plan.ranked_actions || []).map((action, index) => `<tr><td>${index + 1}</td><td><b>${esc(action.label)}</b><br><span style="font-size:.64rem;color:#64748b;">${esc(feasibilityLabels[action.feasibility] || action.feasibility)}</span></td><td>${esc(actionSummary(action))}</td><td class="mono">${num(action.minimum_pair_separation)}</td><td class="mono">${num(action.utility)} [${num(action.utility_interval?.[0])}, ${num(action.utility_interval?.[1])}]</td><td style="font-size:.66rem;white-space:normal;">${esc((action.rationale || []).join('；'))}</td></tr>`).join('');
      App.showModal(`<div class="modal-overlay" onclick="if(event.target===this)App.closeModal()"><div class="modal modal-lg"><div class="flex-between"><h3>确认实验方案 ${esc(plan.plan_id)}</h3><span class="badge badge-medium">研究性建议</span></div><div class="table-wrap" style="margin-top:.7rem;max-height:520px;"><table><thead><tr><th>#</th><th>采集方案</th><th>建议参数</th><th>对最难候选的区分能力</th><th>方案优先级区间</th><th>依据与适用边界</th></tr></thead><tbody>${actions}</tbody></table></div><div style="margin-top:.55rem;font-size:.68rem;color:#64748b;">所有方案都需仪器负责人审核；未执行或未检出不能自动作为候选反证。</div><div class="modal-actions"><button class="btn btn-secondary" onclick="App.closeModal()">关闭</button></div></div></div>`);
    } catch (error) {
      Toast.show(`计划详情加载失败：${error.message || error}`, 'error');
    }
  }

  function download(planId, format) {
    const path = `/tasks/${encodeURIComponent(currentTaskId)}/acquisition-plans/${encodeURIComponent(planId)}/export?format=${encodeURIComponent(format)}`;
    API.download(path, `${currentTaskId}_${planId}_${format}_advisory.csv`);
  }

  return { render, generate, showPlan, download, toggleAll };
})();
