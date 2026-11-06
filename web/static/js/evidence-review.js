/**
 * Compound-level identification evidence workspace.
 *
 * The sample overview and this page share the same task/result indexes.  The
 * overview answers "what is in this sample"; this workspace explains the
 * deduplicated compound conclusion while retaining scans as raw evidence.
 */
const EvidenceReview = (() => {
  const TABS = new Set(['identification', 'unmatched', 'process', 'followup']);
  let activeTab = 'identification';
  let currentTaskId = null;
  let currentBundle = null;
  let selectedResultIndex = null;
  let candidateFilter = 'main';
  let candidateSearch = '';
  let unmatchedRows = [];
  let unmatchedTotal = 0;
  let unmatchedSearch = '';
  let selectedUnmatchedResultIndex = null;
  let unmatchedLoading = false;
  let restoringHash = false;
  let programmaticOpen = false;
  let candidateDetailRenderSequence = 0;
  let unmatchedDetailRenderSequence = 0;

  const esc = value => Tasks.escapeHtml(String(value == null ? '' : value));
  const integer = value => Number(value || 0).toLocaleString();

  function isEnglishInterface() {
    return window.MSSLanguage?.getLanguage?.() === 'en';
  }

  function spectrumPlotText() {
    return isEnglishInterface() ? {} : {
      query: '查询谱', library: '参考谱',
      matchedQuery: '查询谱 · 已记录匹配', otherQuery: '查询谱 · 其他峰',
      matchedLibrary: '参考谱 · 已记录匹配', otherLibrary: '参考谱 · 其他峰',
      locatedPairs: '已定位峰对／保存峰对', referenceMissing: '未保存参考谱',
      empty: '没有可绘制的峰', relativeIntensity: '相对强度（%）',
    };
  }

  function replicateMissingnessSummary(entity) {
    const english = window.MSSLanguage?.getLanguage?.()?.startsWith('en');
    const audit = entity?.replicate_missingness_audit;
    const bins = (entity?.candidates || []).flatMap(item => item.collision_energy_evidence || []);
    const prefix = english
      ? 'Technical replicates at one CE form one condition record, not independent samples. '
      : '同一 CE 的技术重复只形成一个条件记录，不等于独立样本。';
    if (audit?.status === 'incomplete_truncated_topk_vectors') {
      const change = audit.counterfactual_winner_differs
        ? (english ? 'The diagnostic changes the base-ranking leader.' : '诊断计算的基础排序首选候选发生变化。')
        : (english ? 'The diagnostic retains the base-ranking leader.' : '诊断计算的基础排序首选候选未变。');
      return prefix + (english
        ? 'The primary base ranking retains historical zero-filled within-CE medians. V16 separately records observed-only scores; they do not replace production ranking. '
        : '正式基础排序保留历史规则，以补零后的同 CE 分数计算中位数。V16 另存只使用已观测分数的诊断，不用于替换正式排序。')
        + change + (english
          ? ' Absence from Top-K is not a measured zero similarity; complete candidate-by-query rescoring is required to compare all candidates.'
          : '未进入 Top-K 不等于实测零匹配；须对固定候选集合逐谱评分后才能完成比较。');
    }
    if (bins.some(item => item.truncated_topk_absence_imputed_as_zero === true)) {
      return prefix + (english
        ? 'The base ranking retains the historical rule: missing within-CE candidates enter the median as zero. This does not establish a measured zero match.'
        : '基础排序沿用历史规则：同 CE 内缺失候选以零分参与中位数；这不代表实测匹配为零。');
    }
    if (bins.length && bins.every(item => item.truncated_topk_absence_imputed_as_zero === false)) {
      return prefix + (english
        ? 'Only observed scores enter the within-CE median; measured zero scores are retained.'
        : '同 CE 中位数只使用已观测分数，并保留实测零分。');
    }
    return prefix + (english
      ? 'This saved artifact does not declare how missing within-CE candidates were treated.'
      : '该旧工件未声明同 CE 内缺失候选的处理方式，不能推断是否补零。');
  }

  function interfaceCopy(value, fallback = '') {
    const source = String(value == null || value === '' ? fallback : value);
    if (!isEnglishInterface()) return source;
    const translated = window.MSSLanguage?.translateText?.(source) || source;
    return /[\u3400-\u4dbf\u4e00-\u9fff\uf900-\ufaff]/.test(translated)
      ? String(fallback || '')
      : translated;
  }

  function unresolvedEventKind(group) {
    const reason = String(group?.compound?.unresolved_event_reason || '');
    const statuses = Array.isArray(group?.compound?.precursor_event_resolution_statuses)
      ? group.compound.precursor_event_resolution_statuses : [];
    return reason === 'incomplete_candidate_formula_assignment'
      || statuses.includes('formula_assignment_incomplete')
      ? 'incomplete_formula'
      : 'formula_conflict';
  }

  function biosourceInterfaceCopy(value) {
    const source = String(value || '').trim();
    if (!source || !isEnglishInterface()) return source;
    let translated = window.MSSLanguage?.translateText?.(source) || source;
    const pathwayTerms = [
      ['莽草酸途径', 'Shikimate pathway'],
      ['莽草酸', 'shikimate'],
      ['醋酸-丙二酸途径', 'Acetate–malonate pathway'],
      ['乙酸-丙二酸途径', 'Acetate–malonate pathway'],
      ['氨基酸途径', 'Amino-acid pathway'],
      ['糖酵解相关途径', 'Glycolysis-related pathway'],
    ];
    for (const [sourceTerm, targetTerm] of pathwayTerms) {
      translated = translated.split(sourceTerm).join(targetTerm);
    }
    return /[\u3400-\u4dbf\u4e00-\u9fff\uf900-\ufaff]/.test(translated)
      ? 'Biosource pathway context is not available in the current English catalog'
      : translated;
  }

  function candidateName(match) {
    return match?.compound_display_name || match?.compound_name || 'Unnamed candidate';
  }

  function classificationText(match) {
    const display = String(match?.display_class || '').trim();
    const raw = String(match?.classification || '').trim();
    const displayUsable = display && !['unknown', 'other', 'unclassified'].includes(display.toLowerCase());
    const rawUsable = raw && !['unknown', 'other', 'unclassified'].includes(raw.toLowerCase());
    const label = displayUsable ? display : (rawUsable ? raw : '');
    if (!label) {
      const identityHigh = [
        match?.confidence,
        match?.exact_structure_confidence,
        match?.two_dimensional_connectivity_confidence,
      ].some(value => String(value || '').trim().toLowerCase() === 'high');
      return identityHigh
        ? '尚无可靠类别归属（谱库候选证据为 High）'
        : '尚无可靠类别归属';
    }
    const source = String(
      (displayUsable ? match?.display_classification_source : '')
      || match?.classification_source || ''
    ).trim().toLowerCase();
    const trusted = ['curated', 'manual', 'smarts', 'structure', 'rule', 'database_curated']
      .some(prefix => source === prefix || source.startsWith(`${prefix}_`))
      || source.startsWith('runtime_structure_rule');
    const untrusted = !source || ['unknown', 'other', 'existing', 'ml', 'msclass', 'msclass_spectrum', 'fragment', 'fragment_inference', 'spectrum'].includes(source)
      || ['ml_', 'model_', 'msclass_', 'fragment_'].some(prefix => source.startsWith(prefix));
    if (!untrusted && (trusted || match?.classification_reliable === true)) return label;
    return `待核查类别线索：${label}`;
  }

  function groupEvidence(group) {
    const matches = (group?.rows || []).map(row => row.match || {});
    const match = group?.rep?.match || matches[0] || {};
    const auxiliary = group?.compound?.auxiliary_interpretation
      || match.auxiliary_interpretation || {};
    const fragmentation = auxiliary.fragmentation || {};
    const expert = auxiliary.expert || {};
    const biosourceAux = auxiliary.biosource || {};
    const flags = matches.map(match => match.display_evidence_flags || {});
    const biosource = [...new Set([
      ...(Array.isArray(biosourceAux.pathways) ? biosourceAux.pathways : []),
      ...flags.flatMap(flag => Array.isArray(flag.biosynthetic_pathways)
        ? flag.biosynthetic_pathways : []),
    ].map(value => String(value || '').trim()
      .replace(/^生源\s*[:：]\s*/, '')).filter(Boolean))];
    const maximum = key => Math.max(0, ...flags.map(flag => Number(flag[key]) || 0));
    const fragmentationSteps = Number(fragmentation.step_count || 0)
      || maximum('fragmentation_step_count');
    const productStructures = Number(fragmentation.product_structure_count || 0)
      || maximum('product_structure_count');
    return {
      biosource,
      biosourceReasonZh: String(biosourceAux.reason_zh || ''),
      fragmentationSteps,
      productStructures,
      neutralLosses: Number(fragmentation.neutral_loss_count || 0)
        || maximum('neutral_loss_count'),
      diagnosticIons: Number(fragmentation.diagnostic_ion_count || 0)
        || maximum('diagnostic_ion_count'),
      hasFragmentation: fragmentationSteps > 0 || productStructures > 0
        || flags.some(flag => !!flag.has_fragmentation),
      hasLibraryPathwayReference:
        fragmentation.candidate_specific_library_pathway?.status === 'available'
        || flags.some(flag => !!flag.has_library_pathway_reference),
      hasExpert: ['evaluated', 'partially_evaluated'].includes(expert.status),
    };
  }

  function identificationScope(group) {
    if (group?.status?.key === 'isomer') {
      return group.status.hasIsomerHigh
        ? 'High 级同分子式异构体组'
        : '中等证据同分子式异构体组';
    }
    const labels = {
      confirmed: '单一化合物，已确证',
      high: '单一二维连接结构，High 级检出',
      medium: '中等证据结构候选',
      background: '空白／背景相关',
      unresolved_event: '未解析母离子事件；不计作独立化合物',
      low: '低证据谱库线索',
    };
    return labels[group?.status?.key] || '谱库候选';
  }

  function init() {
    const select = document.getElementById('reviewTaskSelect');
    if (select) {
      select.addEventListener('change', () => loadTask(select.value));
    }
    document.addEventListener('mssystem:languagechange', () => {
      if (!currentTaskId) return;
      SampleResultData.clear(currentTaskId);
      if (document.getElementById('page-review')?.classList.contains('active')) {
        loadTask(currentTaskId, selectedResultIndex);
      }
    });
    showTab(activeTab, false);
  }

  async function refreshSelect(preferredTaskId = null) {
    const select = document.getElementById('reviewTaskSelect');
    if (!select) return;
    try {
      const data = await API.get('/tasks?status=completed');
      const tasks = data.tasks || [];
      const target = preferredTaskId || select.value || currentTaskId
        || (typeof SampleOverview !== 'undefined' ? SampleOverview.currentTaskId : null);
      select.innerHTML = '<option value="">-- 选择已完成的任务 --</option>' + tasks.map(task => {
        const stamp = task.upload_time
          ? String(task.upload_time).replace('T', ' ').replace(/\.\d+.*$/, '').slice(0, 16)
          : '时间未记录';
        return `<option value="${esc(task.id)}">${esc(task.filename)} · ${esc(stamp)} · ${integer(task.matched_spectra)}/${integer(task.total_spectra)} MS2 有候选</option>`;
      }).join('');
      if (target && tasks.some(task => task.id === target)) select.value = target;
    } catch (error) {
      console.error('Failed to refresh identification workspace selector:', error);
    }
  }

  function showTab(tab, updateHash = true) {
    activeTab = TABS.has(tab) ? tab : 'identification';
    document.querySelectorAll('[data-review-tab]').forEach(button => {
      const selected = button.dataset.reviewTab === activeTab;
      button.classList.toggle('active', selected);
      button.setAttribute('aria-selected', selected ? 'true' : 'false');
    });
    document.querySelectorAll('[data-review-panel]').forEach(panel => {
      panel.classList.toggle('active', panel.dataset.reviewPanel === activeTab);
    });
    if (currentTaskId && activeTab === 'unmatched') ensureUnmatchedLoaded();
    if (currentTaskId && activeTab === 'process' && typeof ModuleEffects !== 'undefined') {
      ModuleEffects.render(currentTaskId, 'moduleEffectPanel');
    }
    if (currentTaskId && activeTab === 'followup' && typeof EviAcquire !== 'undefined') {
      EviAcquire.render(currentTaskId, 'eviAcquirePanel');
    }
    if (updateHash) updateHashState();
  }

  async function loadTask(taskId, preferredResultIndex = null) {
    const area = document.getElementById('resultsArea');
    if (!taskId) {
      currentTaskId = null;
      currentBundle = null;
      selectedResultIndex = null;
      if (area) area.innerHTML = '<div class="empty-state"><div class="empty-icon">⌕</div><div class="empty-title">选择样本以查看鉴定详情</div><div class="empty-desc">候选证据、分类、生源、规则证据和经验碎裂证据会在这里形成完整证据链</div></div>';
      return;
    }
    currentTaskId = taskId;
    candidateFilter = 'main';
    candidateSearch = '';
    unmatchedRows = [];
    unmatchedTotal = 0;
    unmatchedSearch = '';
    selectedUnmatchedResultIndex = null;
    const select = document.getElementById('reviewTaskSelect');
    if (select) select.value = taskId;
    if (area) area.innerHTML = '<div class="sample-loading"><span class="spinner"></span><b>正在建立候选证据工作台…</b><small>完整详情按候选懒加载，不改变候选排序</small></div>';
    try {
      currentBundle = await SampleResultData.load(taskId);
      if (currentTaskId !== taskId) return;
      if (currentBundle.task?.completed_batch) candidateFilter = 'all';
      if (currentBundle.savedOnly) candidateFilter = 'all';
      if (typeof Results !== 'undefined' && Results.setTaskContext) {
        Results.setTaskContext(taskId, currentBundle.matched, currentBundle.task);
      }
      renderCandidateWorkspace();
      const wanted = preferredResultIndex == null || preferredResultIndex === ''
        ? null : Number(preferredResultIndex);
      const first = currentBundle.groups[0]?.rep?.result_index;
      await selectCandidate(Number.isInteger(wanted) ? wanted : Number(first), false);
      const resultSelect = document.getElementById('resultTaskSelect');
      if (resultSelect && resultSelect.querySelector(`option[value="${CSS.escape(taskId)}"]`)) {
        resultSelect.value = taskId;
      }
      if (activeTab === 'process' && typeof ModuleEffects !== 'undefined') {
        ModuleEffects.render(taskId, 'moduleEffectPanel');
      }
      if (activeTab === 'followup' && typeof EviAcquire !== 'undefined') {
        EviAcquire.render(taskId, 'eviAcquirePanel');
      }
    } catch (error) {
      if (area) area.innerHTML = `<div class="empty-state"><div class="empty-icon">!</div><div class="empty-title">鉴定详情加载失败</div><div class="empty-desc">${esc(error.message || error)}</div></div>`;
    }
  }

  function renderCandidateWorkspace() {
    const area = document.getElementById('resultsArea');
    if (!area || !currentBundle) return;
    const task = currentBundle.task || {};
    const counts = currentBundle.compounds?.metadata || {};
    const exactHigh = Number(counts.high_confidence_2d_connectivity_count
      ?? counts.high_confidence_detected_count
      ?? counts.high_confidence_exact_structure_count ?? 0);
    const isomerHigh = Number(counts.high_confidence_isomer_group_count
      ?? counts.isomer_group_detected_count ?? 0);
    const highScopeGroups = Number(counts.high_confidence_scope_group_count
      ?? (exactHigh + isomerHigh));
    const unresolvedEvents = Number(counts.unresolved_precursor_event_count || 0);
    area.innerHTML = `
      ${typeof FrozenTaskUI !== 'undefined' ? FrozenTaskUI.panel(task) : ''}
      <section class="review-workspace-intro">
        <div><span class="panel-eyebrow">CANDIDATE EVIDENCE WORKSPACE</span><h3>候选鉴定证据工作台</h3>
          <p>${currentBundle.savedOnly ? '当前实体结论不可用，先查看已保存谱级候选与实测 MS/MS。' : '左侧选择候选，右侧对照结构与实测谱图。'}</p></div>
        <div class="review-workspace-kpis">${currentBundle.savedOnly ? `<span><b>${integer(currentBundle.groups.length)}</b> 保存候选身份</span><span><b>${integer(currentBundle.matched.length)}</b> 有候选 MS2</span>` : `<span><b>${integer(counts.main_result_count || 0)}</b> 主要结果</span><span><b>${integer(highScopeGroups)}</b> High 级鉴定组</span><span><b>${integer(exactHigh)}</b> 单一二维连接 / <b>${integer(isomerHigh)}</b> 异构体组</span><span><b>${integer(unresolvedEvents)}</b> 未解析事件</span>`}<span><b>${integer((task.identification_funnel || {}).unmatched_ms2_scans || currentBundle.unmatched?.length || 0)}</b> 未匹配 MS2</span></div>
      </section>
      <div class="identification-workspace">
        <aside class="candidate-rail" aria-label="候选列表">
          <div class="candidate-rail-tools">
            <label for="reviewCandidateSearch">搜索候选</label>
            <input id="reviewCandidateSearch" class="form-input" type="search" placeholder="名称、分子式、类别、InChIKey" oninput="EvidenceReview.setCandidateSearch(this.value)">
            <div class="candidate-rail-filters">
              ${[['main', '主要结果'], ['confirmed', '已确证'], ['high', 'High 单一二维连接'], ['isomer', 'High 异构体组'], ['medium', '中等证据候选'], ['background', '空白／背景'], ['unresolved_event', '未解析事件'], ['low', '证据不足'], ['all', '全部']].map(([key, label]) => `<button data-review-filter="${key}" class="${candidateFilter === key ? 'active' : ''}" onclick="EvidenceReview.setCandidateFilter('${key}')">${label}</button>`).join('')}
            </div>
            <div id="reviewCandidateCount" class="candidate-rail-count"></div>
          </div>
          <div id="reviewCandidateList" class="candidate-list"></div>
        </aside>
        <main id="candidateEvidencePane" class="candidate-evidence-pane">
          <div class="empty-state"><div class="empty-icon">⌕</div><div class="empty-title">请选择一个候选</div><div class="empty-desc">完整证据将按原始 result_index 加载</div></div>
        </main>
      </div>`;
    renderCandidateList();
  }

  function filteredCandidateGroups() {
    if (!currentBundle) return [];
    const query = candidateSearch.trim().toLowerCase();
    return currentBundle.groups.filter(group => {
      if (candidateFilter === 'main' && !group.status.countedAsSampleCompound) return false;
      if (candidateFilter === 'confirmed' && group.status.key !== 'confirmed') return false;
      if (candidateFilter === 'high' && group.status.key !== 'high') return false;
      if (candidateFilter === 'isomer' && !group.status.hasIsomerHigh) return false;
      if (['medium', 'background', 'unresolved_event', 'low'].includes(candidateFilter)
          && group.status.key !== candidateFilter) return false;
      if (!query) return true;
      const match = group.rep?.match || {};
      const evidence = groupEvidence(group);
      return [candidateName(match), ...(group.candidateNames || []), match.formula, match.inchikey, match.classification,
        match.display_class, ...evidence.biosource].filter(Boolean).join(' ').toLowerCase().includes(query);
    });
  }

  function renderCandidateList() {
    const list = document.getElementById('reviewCandidateList');
    const count = document.getElementById('reviewCandidateCount');
    if (!list || !currentBundle) return;
    const groups = filteredCandidateGroups();
    if (count) count.textContent = `显示 ${integer(groups.length)} / ${integer(currentBundle.groups.length)} 个候选`;
    list.innerHTML = groups.length ? groups.map(group => {
      const row = group.rep || {};
      const match = row.match || {};
      const evidence = groupEvidence(group);
      const resultIndex = Number(row.result_index);
      const fragments = [
        evidence.fragmentationSteps ? `${evidence.fragmentationSteps}步碎裂` : '',
        evidence.productStructures ? `${evidence.productStructures}个候选产物结构` : '',
        evidence.neutralLosses ? `${evidence.neutralLosses}项中性丢失` : '',
        evidence.hasLibraryPathwayReference ? '有谱库碎裂参考' : '',
      ].filter(Boolean).join(' · ') || '无可展示碎裂步骤';
      const leadName = String(group.leadCandidate?.compound_name
        || match.compound_name || 'Unnamed candidate');
      const candidateCount = Number(group.candidateMemberCount
        || match.isomer_group_size || 0);
      const candidateSummary = candidateName(match);
      const incompleteFormulaEvent = group.status.key === 'unresolved_event'
        && unresolvedEventKind(group) === 'incomplete_formula';
      const isomerSummary = group.status.key === 'isomer'
        ? `<span>Top1 排序先导 · ${candidateCount || 1} 个竞争候选</span>`
        : (group.status.key === 'unresolved_event'
          ? (isEnglishInterface()
            ? `<span>${candidateCount || 1} unresolved candidates; ${incompleteFormulaEvent ? 'at least one formula assignment is missing' : 'known formulae conflict'}; excluded from the compound count</span>`
            : `<span>${candidateCount || 1} 个未解析候选；${incompleteFormulaEvent ? '至少一个候选缺少分子式' : '已知分子式互相冲突'}；不计作化合物</span>`)
          : '');
      const formulaSummary = group.status.key === 'unresolved_event'
        ? ((group.compound?.candidate_formulas || []).join(' / ')
          || (isEnglishInterface() ? 'No usable candidate formula recorded' : '无可用候选分子式'))
        : (match.formula || '分子式未记录');
      return `<button class="candidate-list-item${resultIndex === selectedResultIndex ? ' selected' : ''}" data-review-result="${resultIndex}" onclick="EvidenceReview.selectCandidate(${resultIndex})">
        <span class="candidate-list-top"><b data-scientific-value>${esc(candidateSummary)}</b><em class="evidence-${esc(group.status.key)}">${esc(group.status.label)}</em></span>${isomerSummary}
        <span>${esc(formulaSummary)} · ${esc(identificationScope(group))}</span>
        <span>${esc(classificationText(match))}</span>
        <span>${esc(isEnglishInterface()
          ? (biosourceInterfaceCopy(evidence.biosource[0])
            || interfaceCopy(evidence.biosourceReasonZh, fragments))
          : (evidence.biosource[0] || evidence.biosourceReasonZh || fragments))}</span>
      </button>`;
    }).join('') : '<div class="candidate-list-empty">当前筛选条件下没有候选。</div>';
  }

  function setCandidateFilter(filter) {
    candidateFilter = ['main', 'confirmed', 'high', 'isomer', 'medium', 'background', 'unresolved_event', 'low', 'all'].includes(filter) ? filter : 'main';
    document.querySelectorAll('[data-review-filter]').forEach(button => {
      button.classList.toggle('active', button.dataset.reviewFilter === candidateFilter);
    });
    renderCandidateList();
  }

  function setCandidateSearch(value) {
    candidateSearch = String(value || '');
    renderCandidateList();
  }

  function decimal(value, digits = 4) {
    const parsed = Number(value);
    return Number.isFinite(parsed) ? parsed.toFixed(digits) : '—';
  }

  function flagshipDecisionText(entity) {
    const status = String(entity?.decision?.status || '');
    const labels = isEnglishInterface() ? {
      selected: 'Selected 2-D connectivity candidate',
      selected_isomer_group: 'Selected only to a same-formula isomer group',
      unresolved_formula_event: 'Unresolved precursor event; no compound identity released',
      selected_candidate_only: 'Limited-scope candidate; no complete connectivity is released',
      abstained: 'Flagship decision abstained',
      unmatched: 'No library candidate',
    } : {
      selected: '已选择二维连接结构候选',
      selected_isomer_group: '仅选择到同分子式异构体组',
      unresolved_formula_event: '未解析母离子事件，不释放化合物身份',
      selected_candidate_only: '有限范围候选，不释放为完整连接结构',
      abstained: '旗舰主动弃权',
      unmatched: '无谱库候选',
    };
    return labels[status] || status || (isEnglishInterface() ? 'Not evaluated' : '未判定');
  }

  async function renderFlagshipCandidateDetail(group) {
    const pane = document.getElementById('candidateEvidencePane');
    if (!pane || !group?.entity || !group?.flagship) return;
    const renderSequence = ++candidateDetailRenderSequence;
    const renderTaskId = currentTaskId;
    const renderResultIndex = Number(group.rep?.result_index);
    try {
      const resultIndex = Number(group.rep?.result_index);
      const detailPayload = Number.isInteger(resultIndex)
        ? await API.get(`/tasks/${encodeURIComponent(currentTaskId)}/results/item/${resultIndex}`)
        : { result: group.rep || {} };
      const raw = detailPayload.result || group.rep || {};
      if (renderSequence !== candidateDetailRenderSequence
          || renderTaskId !== currentTaskId
          || (Number.isInteger(renderResultIndex) && renderResultIndex !== selectedResultIndex)) return;
      if (group.entity._frozen_summary) {
        const payload = await API.get(`/tasks/${encodeURIComponent(currentTaskId)}/flagship?entity_id=${encodeURIComponent(group.entity.entity_id)}&materialize=false`);
        if (renderSequence !== candidateDetailRenderSequence || renderTaskId !== currentTaskId) return;
        if (!payload.entities?.length) throw new Error('保存的实体证据不存在');
        group.entity = payload.entities[0];
        group.flagship = group.entity.winner || group.entity.leading_candidate_when_abstained;
      }
      const entity = group.entity;
      const selected = group.flagship;
      const presentation = group.rep?.match || {};
      // Scientific identity is atomic: database name, formula, InChIKey,
      // SMILES and structure all come from the same V15 candidate.  The
      // compound sidecar contributes only event/sample context below.
      const candidate = { ...(selected.candidate || {}) };
      const decision = entity.decision || {};
      const chemistry = selected.chemistry_gate || {};
      const fragmentation = selected.fragmentation_explanation || {};
      const auxiliary = group.compound?.auxiliary_interpretation
        || presentation.auxiliary_interpretation || {};
      const adductEvidence = auxiliary.adduct_assignment || {};
      const expertEvidence = auxiliary.expert || {};
      const auxiliaryFragmentation = auxiliary.fragmentation || {};
      const productStructures = Array.isArray(auxiliaryFragmentation.plausible_product_structures)
        ? auxiliaryFragmentation.plausible_product_structures : [];
      const scoreComponents = selected.score_components || {};
      const isAbstained = !entity.winner;
      const ceGateRecord = entity.source_paired_ce_evidence || {};
      const ceGate = ceGateRecord.production_gate || {};
      const platformAudit = ceGate.reference_platform_audit || {};
      const exactCeAudit = ceGate.reference_exact_ce_audit || {};
      const featureFragments = selected.candidate?.feature_fragment_evidence || {};
      const featureRuntimeGate = featureFragments.runtime_export_gate || {};
      const featureReasonLabels = {
        no_feature_fragment_rule_passed_sealed_runtime_gate: '封存测试后没有规则达到运行时导出门槛',
        feature_fragment_runtime_gate_invalid: '运行时门文件缺失或校验失败，模块已默认关闭',
      };
      const featureReason = featureReasonLabels[featureFragments.reason]
        || featureFragments.reason || '没有可应用的已验证规则';
      const featureDomain = featureFragments.acquisition_domain || {};
      const featureScoreStandard = featureFragments.feature_fragment_score_standard || {};
      const featureScoreSteps = featureFragments.cumulative_evidence_by_independent_count || [];
      const featureRules = (featureFragments.rule_evidence || []).filter(item =>
        item && ['present', 'eligible_absence'].includes(item.status)
      );
      const featureRuleRows = featureRules.map(item => {
        const target = item.frequency_basis?.target || {};
        const background = item.frequency_basis?.background || {};
        const observation = item.observation || {};
        return `<tr>
          <td class="mono">${esc(item.rule_id || '—')}</td>
          <td>${item.status === 'present' ? `观测 ${decimal(observation.matched_mz, 5)}` : '满足缺失证据质量条件'}</td>
          <td>${integer(target.present)} / ${integer(target.total)} <span class="mono">${esc(target.counting_unit || '')}</span></td>
          <td>${integer(background.present)} / ${integer(background.total)} <span class="mono">${esc(background.counting_unit || '')}</span></td>
          <td class="mono">${decimal(item.likelihood_ratio, 4)}</td>
          <td class="mono">${decimal(item.score_nats ?? item.log_likelihood_ratio, 4)}</td>
          <td class="mono">${decimal(item.score_bits, 4)}</td>
          <td class="mono">${decimal(item.score_decibans, 3)}</td>
          <td>${esc(item.rule_standard?.collision_energy_domain?.raw || item.rule_standard?.collision_energy_domain?.unit || '—')}</td>
          <td class="mono">${esc(item.independent_fragmentation_group || '—')}</td>
          <td>${esc(item.evidence_level || '—')}</td>
        </tr>`;
      }).join('');
      const featureStepRows = featureScoreSteps.map(item => `<tr>
        <td>${integer(item.independent_evidence_count)}</td>
        <td class="mono">${esc(item.independent_fragmentation_group || '—')}</td>
        <td class="mono">${esc(item.selected_rule_id || '—')}</td>
        <td>${esc(item.observation_status || '—')}</td>
        <td class="mono">${decimal(item.increment_score_nats, 4)}</td>
        <td class="mono">${decimal(item.cumulative_score_nats, 4)}</td>
        <td class="mono">${decimal(item.cumulative_likelihood_ratio, 4)}</td>
      </tr>`).join('');
      const ceRows = (selected.collision_energy_evidence || []).map(item => `<tr>
        <td class="mono">${esc(item.collision_energy_bin || '—')}</td>
        <td class="mono">${decimal(item.support, 4)}</td>
        <td class="mono">${decimal(item.legacy_v14_support, 4)}</td>
        <td>${integer(item.replicates_with_candidate)} / ${integer(item.technical_replicates)}</td>
        <td>${decimal(100 * Number(item.scan_top1_fraction || 0), 1)}%</td>
        <td>${esc(item.conditioned_multilibrary_evidence?.selected_consensus_basis || '单一最佳参考谱')}</td>
        <td class="mono">${item.representative_result_index ?? '—'}</td>
      </tr>`).join('');
      const multilibraryEvidence = (selected.collision_energy_evidence || [])
        .map(item => ({ ce: item.collision_energy_bin || '—', evidence: item.conditioned_multilibrary_evidence || {} }))
        .filter(item => item.evidence && Object.keys(item.evidence).length);
      const multilibrarySummaryRows = multilibraryEvidence.map(item => `<tr>
        <td class="mono">${esc(item.ce)}</td>
        <td>${integer(item.evidence.reference_record_count)}</td>
        <td>${integer(item.evidence.unique_content_profile_count)}</td>
        <td>${integer(item.evidence.mirror_duplicate_record_count)}</td>
        <td>${integer(item.evidence.provenance_stratum_count)}</td>
        <td class="mono">${decimal(item.evidence.legacy_best_reference_score, 4)}</td>
        <td class="mono">${decimal(item.evidence.provenance_consensus_score, 4)}</td>
        <td class="mono">${decimal(item.evidence.conditioned_multilibrary_score, 4)}</td>
        <td>${esc(item.evidence.selected_consensus_basis || '—')}</td>
      </tr>`).join('');
      const libraryReferenceRows = multilibraryEvidence.flatMap(item =>
        (item.evidence.reference_records || []).map(reference => ({ ce: item.ce, reference }))
      ).slice(0, 120).map(item => `<tr>
        <td class="mono">${esc(item.ce)}</td>
        <td>${esc(item.reference.source_family || item.reference.dataset_source || '未记录')}</td>
        <td>${esc(item.reference.dataset_source || '—')}</td>
        <td class="mono">${esc(item.reference.source_accession || '—')}</td>
        <td class="mono">${esc(item.reference.spectrum_id || '—')}</td>
        <td class="mono">${decimal(item.reference.spectral_score, 4)}</td>
        <td>${esc(item.reference.condition?.tier || 'instrument_or_ce_metadata_incomplete')}</td>
        <td class="mono">${esc(String(item.reference.provenance_profile_id || '').slice(0, 14) || '—')}</td>
      </tr>`).join('');
      const libraryMatrixRows = multilibraryEvidence.flatMap(item =>
        (item.evidence.library_evidence_matrix || []).map(library => ({ ce: item.ce, evidence: item.evidence, library }))
      ).map(item => `<tr>
        <td class="mono">${esc(item.ce)}</td>
        <td>${esc(item.library.source_family || '—')}</td>
        <td>${esc(item.library.search_status || 'search_state_unknown')}</td>
        <td>${integer(item.library.reference_record_count)}</td>
        <td>${integer(item.library.condition_compatible_record_count)}</td>
        <td class="mono">${decimal(item.library.best_score, 4)}</td>
        <td class="mono">${decimal(item.library.best_condition_compatible_score, 4)}</td>
        <td class="mono">${esc(item.library.representative_source_accession || '—')}</td>
        <td class="mono">${esc(item.library.representative_spectrum_id || '—')}</td>
        <td>${esc(item.library.representative_condition?.tier || '—')}</td>
        <td>${item.library.rank_influence_eligible && item.evidence.rank_influence_allowed ? '是' : '否'}</td>
      </tr>`).join('');
      const edgeRows = (fragmentation.edges || []).slice(0, 40).map(edge => `<tr>
        <td class="mono">${esc(edge.parent_formula || '—')}</td>
        <td>→</td>
        <td class="mono">${esc(edge.child_formula || '—')}</td>
        <td class="mono">−${decimal(edge.observed_loss_da, 4)}</td>
        <td class="mono">${esc(edge.loss_formula || '—')}</td>
        <td class="mono">${decimal(edge.mass_error_da, 5)}</td>
      </tr>`).join('');
      const expertDimensionRows = (expertEvidence.dimensions || []).map(item => `<tr>
        <td>${esc(item.label || item.label_en || item.key || '—')}</td>
        <td>${esc(item.status || 'unavailable')}</td>
        <td class="mono">${item.score == null ? '—' : decimal(item.score, 3)}</td>
        <td>${integer(item.evidence_count || 0)}</td>
        <td>${esc(item.reason || '—')}</td>
      </tr>`).join('');
      const productStructureCards = productStructures.slice(0, 8).map((item, productIndex) => `<div style="width:190px;border:1px solid #bfdbfe;border-radius:8px;padding:.45rem;background:#fff;">
        <canvas id="flagshipProductStructure${productIndex}" width="180" height="120" style="width:180px;height:120px;"></canvas>
        <div style="font-size:.68rem;line-height:1.5;"><b>观测 m/z</b> ${decimal(item.observed_mz, 5)}<br><b>理论 m/z</b> ${decimal(item.expected_mz ?? item.fragment_mass, 5)}<br><b>Δ</b> ${decimal(item.delta, 5)} Da · ${esc(item.adduct || '')}<br><span class="mono">${esc(item.fragment_smiles || '')}</span></div>
      </div>`).join('');
      const candidateRows = (entity.candidates || []).slice(0, 5).map(item => `<tr${item.identity_key === selected.identity_key ? ' style="background:#f0fdf4;"' : ''}>
        <td>${integer(item.entity_rank)}</td>
        <td><b data-scientific-value>${esc(item.candidate?.compound_name || 'Unnamed candidate')}</b><br><span class="mono">${esc(item.identity_key || '—')}</span></td>
        <td class="mono">${decimal(item.entity_score, 4)}</td>
        <td>${integer(item.collision_energy_bins_supported)} / ${integer(item.collision_energy_bins_total)}</td>
        <td>${item.vetoed ? '化学矛盾否决' : '保留'}</td>
      </tr>`).join('');
      const recallSupport = entity.candidate_recall_support || {};
      const expandedCandidates = Array.isArray(recallSupport.expanded_candidates)
        ? recallSupport.expanded_candidates : [];
      const expandedCandidateRows = expandedCandidates.slice(0, 40).map(item => {
        const references = Object.values(item.best_reference_records || {})
          .filter(value => value && typeof value === 'object');
        const bestReference = references[0] || {};
        return `<tr>
          <td>${integer(item.retrieval_rank)}</td>
          <td><b data-scientific-value>${esc(item.candidate?.compound_name || 'Unnamed candidate')}</b><br><span class="mono">${esc(item.connectivity_identity || item.identity_key || '—')}</span></td>
          <td class="mono">${esc(item.candidate?.formula || '—')}</td>
          <td class="mono">${decimal(item.retrieval_rrf_utility, 6)}</td>
          <td>${integer(item.distinct_condition_support)}</td>
          <td>${esc((item.spectral_view_support || []).join(', ') || '—')}</td>
          <td>${esc(bestReference.dataset_source || '—')}<br><span class="mono">${esc(bestReference.source_accession || bestReference.spectrum_id || '—')}</span></td>
          <td>否</td>
        </tr>`;
      }).join('');
      const decisionAssurance = entity.decision_assurance || {};
      const queryQuality = decision.query_spectrum_release_gate
        || decisionAssurance.query_spectrum_quality || {};
      const abstentionGate = decision.cross_platform_abstention_gate
        || decisionAssurance.cross_platform_abstention || {};
      const calibration = decisionAssurance.project_target_decoy_calibration || {};
      const queryQualityRows = (queryQuality.collision_energy_evidence || []).map(item => {
        const scans = item.scans || [];
        const reasons = [...new Set(scans.flatMap(scan => scan.failed_checks || []))];
        const peakCounts = scans.map(scan => Number(scan.query_evidence_peak_count || 0));
        return `<tr>
          <td class="mono">${esc(item.collision_energy_bin || '—')}</td>
          <td>${item.qualified ? '通过' : '未通过'}</td>
          <td>${integer(item.technical_replicate_count)}</td>
          <td>${esc(peakCounts.join(', ') || '—')}</td>
          <td>${esc(reasons.join(', ') || '—')}</td>
        </tr>`;
      }).join('');
      const identityRepair = candidate.canonical_identity_resolution || {};
      const classification = classificationText(candidate);
      const savedScoringView = typeof SpectrumComparator !== 'undefined'
        ? SpectrumComparator.savedEvidence?.(raw.match) || null : null;
      const peaks = savedScoringView ? savedScoringView.query_peaks
        : Array.isArray(raw.query_peaks) ? raw.query_peaks
        : Array.isArray(raw.raw_peaks) ? raw.raw_peaks : [];
      const rt = SampleResultData.retentionTimeMinutes(raw);
      const blockers = Array.isArray(decision.blockers) ? decision.blockers : [];
      const libraryPeaks = savedScoringView ? savedScoringView.library_peaks
        : Array.isArray(raw.lib_peaks) ? raw.lib_peaks
        : Array.isArray(raw.library_peaks) ? raw.library_peaks
          : Array.isArray(raw.match?.lib_peaks) ? raw.match.lib_peaks : [];
      const matchedPeakPairs = savedScoringView ? savedScoringView.matched_pairs
        : Array.isArray(raw.matched_peaks) ? raw.matched_peaks
        : Array.isArray(raw.match?.matched_peaks) ? raw.match.matched_peaks : [];
      const hasMirrorEvidence = libraryPeaks.length > 0;
      const structureHtml = candidate.smiles
        ? `<img src="${API.BASE}/tasks/${encodeURIComponent(currentTaskId)}/flagship/structure?entity_id=${encodeURIComponent(entity.entity_id)}&format=svg" class="structure-img" alt="V15 flagship structure" onerror="this.outerHTML='<span class=\'text-muted\'>旗舰结构式渲染不可用</span>'">`
        : '<span class="text-muted">该旗舰候选没有可渲染结构</span>';
      const fragReason = String(fragmentation.reason || '当前没有可发布的分子式守恒边');
      const unresolvedEvent = group.status?.key === 'unresolved_event';
      const incompleteFormulaEvent = unresolvedEvent
        && unresolvedEventKind(group) === 'incomplete_formula';
      const unresolvedEventBanner = !unresolvedEvent ? '' : (isEnglishInterface()
        ? `<section class="pw-card" style="padding:.75rem .85rem;margin-bottom:1rem;border:1px solid #fdba74;border-left:5px solid #c2410c;background:#fff7ed;line-height:1.6;"><h4 style="margin:0 0 .35rem;">Unresolved precursor event — excluded from the independent compound count</h4><p style="margin:0;">${incompleteFormulaEvent ? 'Candidate structures disagree and at least one candidate lacks a usable molecular formula.' : 'Known candidate molecular formulae conflict.'} Database names shown below are ranking candidates retained for audit; none is released as the primary compound identity for this precursor event.</p></section>`
        : `<section class="pw-card" style="padding:.75rem .85rem;margin-bottom:1rem;border:1px solid #fdba74;border-left:5px solid #c2410c;background:#fff7ed;line-height:1.6;"><h4 style="margin:0 0 .35rem;">未解析母离子事件——不计作独立化合物</h4><p style="margin:0;">${incompleteFormulaEvent ? '候选结构不一致，且至少一个候选缺少可用分子式。' : '已知候选分子式互相冲突。'}下方数据库名称仅作为审计用排序候选保留；本母离子事件不释放主要化合物身份。</p></section>`);
      const isomerGroup = group.status?.key === 'isomer'
        || decision.scope === 'molecular_formula_isomer_group';
      const identityScopeLabel = unresolvedEvent
        ? (isEnglishInterface() ? 'Unresolved precursor event' : '未解析母离子事件')
        : isAbstained
          ? (isEnglishInterface() ? 'No identity released' : '未释放化合物身份')
        : isomerGroup
          ? (isEnglishInterface() ? 'Same-formula isomer group' : '同分子式异构体候选组')
          : (isEnglishInterface() ? '2D structure candidate' : '二维连接结构候选');
      const isomerGroupBanner = unresolvedEvent || !isomerGroup ? '' : (isEnglishInterface()
        ? '<section class="pw-card" style="padding:.75rem .85rem;margin-bottom:1rem;border:1px solid #93c5fd;border-left:5px solid #2563eb;background:#eff6ff;line-height:1.6;"><h4 style="margin:0 0 .35rem;">Unresolved same-formula isomer group</h4><p style="margin:0;">The structure and database name below belong to the ranking lead only. They do not mean that this individual connectivity has been resolved.</p></section>'
        : '<section class="pw-card" style="padding:.75rem .85rem;margin-bottom:1rem;border:1px solid #93c5fd;border-left:5px solid #2563eb;background:#eff6ff;line-height:1.6;"><h4 style="margin:0 0 .35rem;">同分子式异构体候选组尚未解析</h4><p style="margin:0;">下方结构与数据库名称仅属于排序先导，不表示该单一连接结构已经被解析。</p></section>');
      const releaseQualityCalibration = isEnglishInterface()
        ? `<section class="pw-card" style="padding:.75rem;margin-top:1rem;"><h4>Release, query-spectrum quality, and project calibration</h4><div style="font-size:.75rem;line-height:1.65;"><b>Cross-platform abstention gate:</b> ${esc(interfaceCopy(abstentionGate.status, 'not_applicable'))}; failed checks: ${esc(interfaceCopy((abstentionGate.failed_gates || []).join(', '), 'none'))}. This gate changes release scope only; it does not alter the Top 1, Top 5, or entity score.<br><b>Query-spectrum quality:</b> ${esc(interfaceCopy(queryQuality.status, 'not_applicable'))}; eligible CE conditions: ${integer(queryQuality.observed_quality_qualified_ce_bins ?? queryQuality.quality_qualified_ce_bins)} / ${integer(queryQuality.observed_total_ce_bins)}. Locked-transfer testing showed that this gate removed too many correct releases, so it is fixed as audit-only and does not change release scope.<br><b>Project target–decoy:</b> ${esc(interfaceCopy(calibration.status, 'unavailable'))}; q-value ${calibration.q_value == null ? 'not generated' : decimal(calibration.q_value, 4)}; reason: ${esc(interfaceCopy(calibration.reason, 'no eligible project-level calibration assignment is linked to this entity'))}. Neither raw similarity nor the entity score is a probability.</div>${queryQualityRows ? `<details style="margin-top:.55rem;"><summary>Per-CE query-spectrum quality audit</summary><div class="table-wrap"><table><thead><tr><th>CE bin</th><th>Status</th><th>Technical replicates</th><th>Valid peak count</th><th>Failed checks</th></tr></thead><tbody>${queryQualityRows}</tbody></table></div></details>` : ''}</section>`
        : `<section class="pw-card" style="padding:.75rem;margin-top:1rem;"><h4>发布、查询谱质量与项目校准</h4><div style="font-size:.75rem;line-height:1.65;"><b>跨平台弃权门：</b>${esc(abstentionGate.status || 'not_applicable')}；失败项 ${esc((abstentionGate.failed_gates || []).join(', ') || '无')}。该门只改变发布范围，不改变 Top 1、Top 5 或实体分。<br><b>查询谱质量：</b>${esc(queryQuality.status || 'not_applicable')}；合格 CE ${integer(queryQuality.observed_quality_qualified_ce_bins ?? queryQuality.quality_qualified_ce_bins)} / ${integer(queryQuality.observed_total_ce_bins)}。该门在锁定转移中误删过多正确发布，当前固定为 audit-only，不改变发布范围。<br><b>项目 target–decoy：</b>${esc(calibration.status || 'unavailable')}；q-value ${calibration.q_value == null ? '未生成' : decimal(calibration.q_value, 4)}；原因 ${esc(calibration.reason || '当前实体未绑定合格的项目级校准分配')}。原始相似度和实体分均不是概率。</div>${queryQualityRows ? `<details style="margin-top:.55rem;"><summary>逐 CE 查询谱质量审计</summary><div class="table-wrap"><table><thead><tr><th>CE 箱</th><th>状态</th><th>技术重复</th><th>有效峰数</th><th>失败检查</th></tr></thead><tbody>${queryQualityRows}</tbody></table></div></details>` : ''}</section>`;
      pane.innerHTML = `<article class="candidate-inline-detail">
        ${unresolvedEventBanner}
        ${isomerGroupBanner}
        <section style="padding:.75rem .85rem;border:1px solid ${isAbstained ? '#fcd34d' : '#86efac'};border-left:5px solid ${isAbstained ? '#d97706' : '#16a34a'};background:${isAbstained ? '#fffbeb' : '#f0fdf4'};border-radius:8px;line-height:1.55;">
          <div class="flex-between"><div><span class="panel-eyebrow">V15 CONDITIONED MULTI-LIBRARY + CE DECISION</span><h3 data-scientific-value>${esc(candidateName(candidate))}</h3></div><span class="sample-evidence-badge evidence-${esc(group.status.key)}">${esc(flagshipDecisionText(entity))}</span></div>
          <div>实体 <span class="mono">${esc(entity.entity_source_value || entity.entity_id)}</span> · ${integer(entity.member_scan_count)} 张 MS2 · ${integer(entity.distinct_collision_energy_bins)} 个碰撞能量箱 · RT ${rt == null ? '—' : rt.toFixed(3) + ' min'} · 前体 m/z ${decimal(entity.representative_precursor_mz, 5)}</div>
          <div><b>实体排序分 ${decimal(selected.entity_score, 4)}</b>（不是概率） · 胜者间隔 ${decimal(decision.winner_margin, 4)} · 声明范围 <span class="mono">${esc(decision.scope || 'none')}</span>${blockers.length ? (isEnglishInterface() ? ` · limitations: ${esc(blockers.join(', '))}` : ` · 限制：${esc(blockers.join(', '))}`) : ''}</div>
        </section>

        <section class="result-evidence-hero" aria-label="${esc(isEnglishInterface() ? 'Candidate visual evidence' : '候选可视证据')}">
          <div class="result-identity-strip ${unresolvedEvent || isAbstained ? 'scope-unresolved' : isomerGroup ? 'scope-isomer' : 'scope-candidate'}">
            <span class="result-identity-kicker">${esc(isEnglishInterface() ? 'IDENTITY SCOPE' : '身份范围')}</span>
            <strong>${esc(identityScopeLabel)}</strong>
            <span>${esc(isEnglishInterface()
              ? (isAbstained ? 'The drawing is retained as an audit candidate and is not a released identity.' : isomerGroup ? 'The drawing is the ranking lead within the group.' : 'Library-supported candidate; no standard confirmation is claimed.')
              : (isAbstained ? '结构图仅为审计用排序候选，不是已发布身份。' : isomerGroup ? '结构图为组内排序先导，不表示单一异构体已解析。' : '谱库支持的结构候选，未声称标准品确证。'))}</span>
          </div>
          <div class="result-visual-grid">
            <section class="result-structure-panel">
              <div class="result-visual-heading"><span>${esc(isEnglishInterface() ? 'Candidate structure' : '候选二维结构')}</span><small>${esc(candidate.formula || '—')}</small></div>
              <div class="result-structure-stage">${structureHtml}</div>
              <div class="result-structure-caption"><b data-scientific-value>${esc(candidateName(candidate))}</b><span class="mono">${esc(candidate.inchikey || '—')}</span></div>
            </section>
            <section class="result-spectrum-panel">
              <div class="result-visual-heading"><span>${esc(hasMirrorEvidence ? (isEnglishInterface() ? 'Measured vs library MS/MS' : '实测谱 / 馆藏谱镜像对照') : (isEnglishInterface() ? 'Measured MS/MS' : '实测 MS/MS 棒谱'))}</span><small>${integer(peaks.length)} ${esc(isEnglishInterface() ? 'query peaks' : '个查询峰')}</small></div>
              <div id="flagshipEvidenceSpectrum" class="result-spectrum-stage" aria-live="polite"><div class="result-visual-empty">${esc(isEnglishInterface() ? 'Preparing spectrum…' : '正在准备谱图…')}</div></div>
            </section>
          </div>
          <div class="result-metric-ribbon">
            <div><span>${esc(isEnglishInterface() ? 'Entity rank score' : '实体排序分')}</span><b>${decimal(selected.entity_score, 4)}</b><small>${esc(isEnglishInterface() ? 'ranking score' : '排序分，不是概率')}</small></div>
            <div><span>${esc(isEnglishInterface() ? 'Winner margin' : '胜者间隔')}</span><b>${decimal(decision.winner_margin, 4)}</b><small>${esc(isEnglishInterface() ? 'score difference' : '分数差')}</small></div>
            <div><span>${esc(isEnglishInterface() ? 'Collision energy' : '碰撞能量')}</span><b>${integer(entity.distinct_collision_energy_bins)}</b><small>${esc(isEnglishInterface() ? 'CE bins' : '个 CE 条件箱')}</small></div>
            <div><span>${esc(isEnglishInterface() ? 'Measured evidence' : '实测证据')}</span><b>${integer(entity.member_scan_count)}</b><small>MS/MS scans</small></div>
            <div><span>RT</span><b>${rt == null ? '—' : rt.toFixed(3)}</b><small>min</small></div>
            <div><span>${esc(isEnglishInterface() ? 'Precursor' : '前体')}</span><b>${decimal(entity.representative_precursor_mz, 5)}</b><small>m/z</small></div>
          </div>
        </section>

        <details class="result-evidence-details">
          <summary><span>${esc(isEnglishInterface() ? 'Detailed evidence' : '详细证据')}</span><small>${esc(isEnglishInterface() ? 'methods, candidate tables, raw metadata and expert interpretation' : '方法、候选表、原始元数据与完整专家解释')}</small></summary>
          <div class="result-evidence-details-body">

        <div style="display:flex;flex-wrap:wrap;gap:1rem;margin-top:1rem;">
          <section class="pw-card" style="padding:.75rem;text-align:center;flex:1 1 250px;"><h4>旗舰候选结构</h4>${structureHtml}<div style="text-align:left;font-size:.74rem;line-height:1.55;margin-top:.5rem;"><b>Formula</b> ${esc(candidate.formula || '—')}<br><b>InChIKey</b> <span class="mono">${esc(candidate.inchikey || '—')}</span><br><b>类别</b> ${esc(classification)}<br><span style="color:var(--text-secondary);">分类只作辅助显示，身份排序影响固定为 0。</span>${identityRepair.status ? `<br><b>身份修复</b> ${esc(identityRepair.status)} · ${esc(identityRepair.source || '')}` : ''}</div></section>
          <section class="pw-card" style="padding:.75rem;flex:2 1 380px;"><h4>V15 实体基础分与最终选择</h4><div class="score-row">
            ${[['平均 CE 支持', scoreComponents.mean_collision_energy_support], ['中位 CE 支持', scoreComponents.median_collision_energy_support], ['CE 覆盖', scoreComponents.collision_energy_coverage], ['CE Top-1 比例', scoreComponents.collision_energy_top1_fraction]].map(([label, value]) => `<div class="score-chip"><div class="score-val">${decimal(value, 4)}</div><div class="score-lbl">${label}</div></div>`).join('')}
          </div><div style="margin-top:.65rem;padding:.55rem .65rem;background:#eff6ff;border-left:4px solid #2563eb;font-size:.74rem;line-height:1.55;"><b>最终选择依据：</b>${esc(decision.selection_basis || 'conditioned_multilibrary_entity_score')}。实体分不是概率；若 CE 配对门通过，页面保留原实体分与原排名，并单独记录选择覆盖，不把 CE 差值伪装成概率。<br><b>被排除的身份输入：</b>规则聚合值 0 · 分类 0 · 生源通路 0 · 候选通路 0。</div></section>
        </div>

        <section class="pw-card" style="padding:.75rem;margin-top:1rem;"><h4>跨碰撞能量证据</h4><div class="table-wrap"><table><thead><tr><th>CE 箱</th><th>V15 支持分</th><th>V14 审计分</th><th>有候选/技术重复</th><th>扫描 Top-1 比例</th><th>多谱库聚合依据</th><th>代表 result_index</th></tr></thead><tbody>${ceRows || '<tr><td colspan="7">没有可用碰撞能量证据</td></tr>'}</tbody></table></div><p style="font-size:.7rem;color:var(--text-secondary);">${esc(replicateMissingnessSummary(entity))}</p></section>

        <section class="pw-card" style="padding:.75rem;margin-top:1rem;"><h4>多谱库证据矩阵</h4><div style="font-size:.74rem;line-height:1.55;margin-bottom:.55rem;">同一峰内容的镜像记录先按归一化峰表哈希合并；accession、谱图编号和来源名称不会各自增加一票。只有仪器域/机型或供应商+分析器族锦标匹配、碎裂方式一致、CE 均为 fixed 模式、单位相同且差值在门限内的谱图才能改分。全采集条件中位数只展示，条件元数据缺失时保留 V14 最佳参考谱分。</div><div class="table-wrap"><table><thead><tr><th>查询 CE</th><th>参考记录</th><th>唯一峰谱</th><th>镜像重复</th><th>来源层</th><th>单谱最佳</th><th>条件共识</th><th>V15 谱图分</th><th>采用依据</th></tr></thead><tbody>${multilibrarySummaryRows || '<tr><td colspan="9">本任务由旧版结果生成，未保存逐参考谱峰哈希；V15 不会据此伪造多库增益。</td></tr>'}</tbody></table></div><details open style="margin-top:.55rem;"><summary>逐谱库命中/未命中、条件相容性和排序影响</summary><div class="table-wrap"><table><thead><tr><th>查询 CE</th><th>谱库</th><th>检索状态</th><th>参考谱</th><th>条件相容</th><th>单库最佳</th><th>相容最佳</th><th>accession</th><th>谱图 ID</th><th>条件层级</th><th>改分</th></tr></thead><tbody>${libraryMatrixRows || '<tr><td colspan="11">旧任务未保存谱库参与清单，无法区分“未命中”和“未检索”</td></tr>'}</tbody></table></div></details><details style="margin-top:.55rem;"><summary>逐参考谱峰哈希追溯（最多 120 条）</summary><div class="table-wrap"><table><thead><tr><th>查询 CE</th><th>来源族</th><th>原始库名</th><th>accession</th><th>谱图 ID</th><th>单谱分</th><th>条件关系</th><th>峰谱哈希组</th></tr></thead><tbody>${libraryReferenceRows || '<tr><td colspan="8">没有逐参考谱记录</td></tr>'}</tbody></table></div></details></section>

        <section class="pw-card" style="padding:.75rem;margin-top:1rem;"><h4>同分子式 CE 配对选择门</h4><div style="font-size:.75rem;line-height:1.6;"><b>结果：</b>${esc(ceGate.status || '未评价')} · ${esc(ceGate.reason || '该实体不满足近分同分子式候选条件')}<br><b>基础候选 → 最终候选：</b><span class="mono">${esc(ceGate.baseline_identity_key || selected.identity_key || '—')}</span> → <span class="mono">${esc(ceGate.selected_identity_key || selected.identity_key || '—')}</span><br><b>同平台参考子集：</b>${integer(platformAudit.paired_observation_count)} 个配对观测，${integer(platformAudit.distinct_query_collision_energy_condition_count)} 个查询 CE 条件，${integer(platformAudit.spectral_view_count)} 个谱图视图，中位候选差 ${decimal(platformAudit.median_score_delta, 4)}，胜出比例 ${decimal(100 * Number(platformAudit.observation_win_fraction || 0), 1)}%。<br><b>参考谱 CE 明确且近邻：</b>${integer(exactCeAudit.paired_observation_count)} 个配对观测；为 0 时表示谱库参考记录没有足够 CE 元数据，系统不会把缺失当作同 CE。<br><span style="color:var(--text-secondary);">该门只允许 V15 前五中的同分子式连接体切换；扩展召回候选不能直接成为生产 Top‑1。通过后仍标记为 MS/MS 谱库连接体候选，不能写成标准品确证。</span></div></section>

        <section class="pw-card" style="padding:.75rem;margin-top:1rem;"><h4>Top40 扩展候选支持</h4><div style="font-size:.75rem;line-height:1.6;margin-bottom:.55rem;"><b>状态：</b>${esc(recallSupport.status || (Object.keys(recallSupport).length ? 'evaluated' : 'unavailable'))} · 新增 ${integer(expandedCandidates.length)} 个二维候选假说；目标候选池 ${integer(recallSupport.target_total_connectivities || 40)}。扩展项只能用于后续复核，不能替换当前赢家、不能取得二维结构发布、不能进入残差模型评分；RRF 数值是检索效用，不是概率。</div>${expandedCandidateRows ? `<div class="table-wrap"><table><thead><tr><th>召回名次</th><th>候选</th><th>分子式</th><th>RRF 效用</th><th>CE 条件</th><th>谱图视图</th><th>最佳参考</th><th>可替换赢家</th></tr></thead><tbody>${expandedCandidateRows}</tbody></table></div>` : '<div class="sample-no-results">该实体没有新增扩展候选，或旧任务尚未生成 V15 Top40 支持工件。</div>'}</section>

        <section id="v16CandidateComparison" class="pw-card" style="padding:.75rem;margin-top:1rem;"></section>

        ${releaseQualityCalibration}

        <section class="pw-card" style="padding:.75rem;margin-top:1rem;"><h4>候选无关化学门</h4><div style="font-size:.76rem;line-height:1.6;"><b>${esc(chemistry.status || 'unavailable')}</b> · ${esc(chemistry.reason || '未获得独立分子式/加合物证据')}<br>Formula ${esc(chemistry.formula || candidate.formula || '—')} · adduct ${esc(chemistry.adduct || '—')} · 理论/实测 m/z ${decimal(chemistry.expected_precursor_mz, 5)} / ${decimal(chemistry.observed_precursor_mz, 5)} · ${decimal(chemistry.mass_error_ppm, 2)} ppm<br><span style="color:var(--text-secondary);">馆藏谱加合物仅审计，不能单独否决中性化学身份；缺失化学元数据保持中性。</span></div></section>

        <section class="pw-card" style="padding:.75rem;margin-top:1rem;"><h4>样品离子形式与规则证据覆盖</h4><div style="font-size:.75rem;line-height:1.6;margin-bottom:.55rem;"><b>样品离子形式</b> ${esc(adductEvidence.sample_ion_form || '未可靠判定')} · ${esc(adductEvidence.sample_ion_form_source || 'unresolved')} · ${adductEvidence.candidate_independent ? '候选无关样品证据' : '候选依赖假设'}<br><b>谱库加合物</b> ${esc(adductEvidence.library_adduct || '—')} · 一致性 ${esc(adductEvidence.sample_library_consistency || 'not_tested')}<br><b>规则证据覆盖</b> 有支持证据 ${integer(expertEvidence.coverage?.scored || 0)} 维；已评价但未观察 ${integer(expertEvidence.coverage?.no_evidence || 0)} 维；不适用 ${integer(expertEvidence.coverage?.not_applicable || 0)} 维；未建立 ${integer(expertEvidence.coverage?.unavailable || 0)} 维。规则证据不参与身份排序。</div><details style="font-size:.7rem;color:var(--text-secondary);margin-bottom:.55rem;"><summary>展开旧版数值审计字段</summary><div style="padding:.35rem 0 0 .6rem;">未校准规则聚合 ${decimal(expertEvidence.score, 3)}；不是概率，不得用于身份等级或排序。</div></details>${expertDimensionRows ? `<div class="table-wrap"><table><thead><tr><th>维度</th><th>状态</th><th>维度分</th><th>证据数</th><th>原因</th></tr></thead><tbody>${expertDimensionRows}</tbody></table></div>` : '<div class="sample-no-results">没有保存可逐维审计的规则证据结果</div>'}</section>

        <section class="pw-card" style="padding:.75rem;margin-top:1rem;"><h4>V15 特征碎片计分与审计</h4><div style="font-size:.75rem;line-height:1.6;margin-bottom:.55rem;"><b>状态：</b>${esc(featureFragments.status || 'unavailable')} · ${esc(featureReason)}<br><b>运行时导出：</b>${integer(featureRuntimeGate.exported_rule_count)} 条；30 谱试运行规则 ${featureRuntimeGate.pilot_v1_runtime_allowed ? '被错误启用' : '未启用'}；研究 catalog ${featureRuntimeGate.research_catalog_evaluated ? '已计算' : '未计算'}。<br><b>一个独立碎片：</b>ln[P(该观测|目标类,采集条件)/P(该观测|背景类,采集条件)]；<b>两个独立碎片：</b>两个 lnLR 相加。同一独立碎裂组的相关规则只取 |lnLR| 最大的一条，所以“碎片个数”本身没有固定加分。<br><b>CE 适用标准：</b>${esc(featureScoreStandard.collision_energy_requirement || '单位、fixed 模式、仪器域和规则 CE 范围必须一致')}${isEnglishInterface() ? '.' : '。'}<br><b>独立类别前提：</b>${esc(featureFragments.classification_premise?.chemical_class || '未建立')} · ${featureFragments.classification_premise?.reliable ? '结构或人工整理来源可核查' : '未通过来源检查'}<br><b>采集域：</b><span class="mono">${esc(featureDomain.instrument_domain_id || '未解析')}</span> · ${esc(featureDomain.activation || '活化方式缺失')} · ${esc(featureDomain.collision_energy_unit || 'CE 单位缺失')} · ${esc(featureDomain.collision_energy_mode || '采集模式缺失')}；跨单位换算 ${featureDomain.cross_unit_conversion_applied ? '发生（结果无效）' : '未发生'}。<br><b>当前证据分：</b>${decimal(featureFragments.combined_score_nats ?? featureFragments.combined_log_likelihood_ratio, 4)} nats · ${decimal(featureFragments.combined_score_bits, 4)} bits · ${decimal(featureFragments.combined_score_decibans, 3)} decibans${isEnglishInterface() ? ';' : '；'}<b>生产排序贡献：</b>0，不改变 Top 1、Top 5、最终分、证据等级或标准品确证状态。</div>${featureRuleRows ? `<div class="table-wrap"><table><thead><tr><th>规则</th><th>查询观测</th><th>目标身份</th><th>背景身份</th><th>LR</th><th>lnLR (nats)</th><th>bits</th><th>decibans</th><th>规则 CE</th><th>独立碎裂组</th><th>验证层级</th></tr></thead><tbody>${featureRuleRows}</tbody></table></div>` : '<div class="sample-no-results">当前运行时没有通过封存测试门槛的特征碎片规则；系统仍显示计分公式和 CE 适用标准，但不伪造规则命中。</div>'}${featureStepRows ? `<details open style="margin-top:.55rem;"><summary>按独立碎片个数展开累积分</summary><div class="table-wrap"><table><thead><tr><th>独立证据数</th><th>碎裂组</th><th>选中规则</th><th>观测</th><th>本条 lnLR</th><th>累积 lnLR</th><th>累积 LR</th></tr></thead><tbody>${featureStepRows}</tbody></table></div></details>` : ''}</section>

        <section class="pw-card" style="padding:.75rem;margin-top:1rem;"><h4>观测峰支持的候选产物结构</h4>${productStructureCards ? `<div style="display:flex;gap:.65rem;flex-wrap:wrap;">${productStructureCards}</div>` : '<div class="sample-no-results">当前没有可展示的候选产物结构；这不等于没有 MS/MS，质量差与分子式候选仍在碎裂证据中保留。</div>'}<p style="font-size:.7rem;color:var(--text-secondary);">结构由候选前体断键枚举与精确质量匹配得到，不等于已确认气相机理或唯一产物结构。</p></section>

        <section class="pw-card" style="padding:.75rem;margin-top:1rem;"><h4>分子式守恒碎裂解释</h4>${edgeRows ? `<div style="font-size:.74rem;margin-bottom:.45rem;">状态 <b>${esc(fragmentation.status)}</b> · ${integer(fragmentation.assigned_fragment_nodes)}/${integer(fragmentation.fragment_nodes)} 个碎片节点获分子式 · ${integer(fragmentation.formula_conserving_edge_count)} 条守恒边。仅显示前 40 条。</div><div class="table-wrap"><table><thead><tr><th>母离子式</th><th></th><th>子离子式</th><th>实测 Δm</th><th>元素差</th><th>误差 Da</th></tr></thead><tbody>${edgeRows}</tbody></table></div>` : `<div style="padding:.6rem;background:#f8fafc;border-left:4px solid #94a3b8;font-size:.75rem;line-height:1.55;"><b>${esc(fragmentation.status || 'unavailable')}：</b>${esc(fragReason)}。系统不会用纯质量差编造碎裂机理。</div>`}<p style="font-size:.7rem;color:var(--text-secondary);">这些边只证明元素组成可守恒，不等同于已确认的气相反应机理，也不支持立体化学声明。</p></section>

        <section class="pw-card" style="padding:.75rem;margin-top:1rem;"><h4>实体候选短名单</h4><div class="table-wrap"><table><thead><tr><th>#</th><th>二维身份</th><th>实体分</th><th>CE 支持</th><th>化学门</th></tr></thead><tbody>${candidateRows || '<tr><td colspan="5">无候选</td></tr>'}</tbody></table></div></section>

        <section class="pw-card" style="padding:.75rem;margin-top:1rem;"><h4>代表原始 MS/MS 定位</h4><div style="font-size:.72rem;color:var(--text-secondary);margin-bottom:.45rem;">${esc(currentBundle?.task?.filename || '未记录源文件')} · spectrum ${esc(raw.spectrum_id || '—')} · result_index ${Number.isInteger(resultIndex) ? resultIndex : '—'}。首屏谱图使用该扫描的真实查询峰；此扫描的旧 Top‑1、规则证据和经验通路不能覆盖旗舰实体决策。</div></section>
          </div>
        </details>
      </article>`;
      if (typeof CandidateComparison !== 'undefined') {
        CandidateComparison.mount(document.getElementById('v16CandidateComparison'), currentTaskId, entity);
      }
      setTimeout(() => {
        if (renderSequence !== candidateDetailRenderSequence || renderTaskId !== currentTaskId) return;
        const element = document.getElementById('flagshipEvidenceSpectrum');
        if (!element) return;
        element.innerHTML = '';
        if (peaks.length && hasMirrorEvidence && typeof SpectrumComparator !== 'undefined') {
          new SpectrumComparator().render('flagshipEvidenceSpectrum',
            { peaks, name: savedScoringView
                ? (isEnglishInterface() ? 'Scoring query' : '评分查询谱')
                : (isEnglishInterface() ? 'Measured query' : '实测查询谱') },
            { peaks: libraryPeaks, name: candidate.compound_name || (isEnglishInterface() ? 'Library' : '馆藏谱') },
            matchedPeakPairs, {
              height: 330,
              scorePrecision: 4,
              score: savedScoringView ? savedScoringView.comparison_cosine : raw.match?.cosine,
              scoreLabel: savedScoringView
                ? (isEnglishInterface() ? 'Saved comparison cosine' : '已保存候选比较余弦')
                : (isEnglishInterface() ? 'Saved cosine' : '已保存余弦'),
              title: isEnglishInterface() ? 'MS/MS spectrum comparison' : 'MS/MS 谱图对照',
              text: spectrumPlotText(),
            });
        } else if (peaks.length && typeof SpectrumViewer !== 'undefined') {
          new SpectrumViewer('flagshipEvidenceSpectrum').render(peaks, { title: 'MS/MS', height: 330 });
        } else {
          element.innerHTML = `<div class="result-visual-empty">${esc(peaks.length
            ? (isEnglishInterface() ? 'Spectrum renderer is unavailable; the measured peaks remain in the detailed evidence.' : '谱图组件未加载；实测峰仍保留在详细证据中。')
            : (isEnglishInterface() ? 'No measured fragment peaks were saved for this scan.' : '该代表扫描没有保存可绘制的实测碎片峰。'))}</div>`;
        }
        if (typeof MoleculeViewer !== 'undefined') {
          const viewer = new MoleculeViewer({ width: 180, height: 120 });
          productStructures.slice(0, 8).forEach((item, productIndex) => {
            viewer.render(`flagshipProductStructure${productIndex}`, item.fragment_smiles || '');
          });
        }
      }, 80);
    } catch (error) {
      if (renderSequence !== candidateDetailRenderSequence || renderTaskId !== currentTaskId) return;
      pane.innerHTML = `<div class="empty-state"><div class="empty-icon">!</div><div class="empty-title">旗舰实体详情加载失败</div><div class="empty-desc">${esc(error.message || error)}</div></div>`;
    }
  }

  async function selectCandidate(resultIndex, updateHash = true) {
    const wanted = Number(resultIndex);
    if (!currentBundle || !Number.isInteger(wanted)) return;
    const group = currentBundle.groups.find(item => Number(item.rep?.result_index) === wanted);
    if (!group) {
      Toast.show('未能在当前样本定位该候选', 'error');
      return;
    }
    selectedResultIndex = wanted;
    document.querySelectorAll('[data-review-result]').forEach(button => {
      button.classList.toggle('selected', Number(button.dataset.reviewResult) === wanted);
    });
    const pane = document.getElementById('candidateEvidencePane');
    if (pane) pane.innerHTML = '<div class="sample-loading"><span class="spinner"></span><b>正在加载完整候选证据…</b><small>分类、生源、规则证据、碎裂和原始谱图均来自该原始结果</small></div>';
    if (updateHash) updateHashState();
    if (group.entity && group.flagship) {
      await renderFlagshipCandidateDetail(group);
    } else if (typeof Results !== 'undefined' && Results.showDetailFromSummary) {
      await Results.showDetailFromSummary(
        currentTaskId, group.rep, 'candidateEvidencePane', currentBundle.matched, currentBundle.task
      );
    }
  }

  async function ensureUnmatchedLoaded() {
    if (unmatchedRows.length || unmatchedLoading || !currentTaskId) {
      if (unmatchedRows.length) renderUnmatchedWorkspace();
      return;
    }
    await loadMoreUnmatched(true);
  }

  async function loadMoreUnmatched(reset = false) {
    if (!currentTaskId || unmatchedLoading) return;
    unmatchedLoading = true;
    const panel = document.getElementById('unmatchedWorkbenchPanel');
    if (reset && panel) panel.innerHTML = '<div class="sample-loading"><span class="spinner"></span><b>正在加载未匹配谱图…</b><small>默认每次加载 50 张，避免一次创建数千行</small></div>';
    try {
      if (reset) unmatchedRows = [];
      const offset = unmatchedRows.length;
      const payload = await API.get(`/tasks/${encodeURIComponent(currentTaskId)}/results?projection=summary&match_state=unmatched&limit=50&offset=${offset}`);
      unmatchedRows = unmatchedRows.concat(Array.isArray(payload.results) ? payload.results : []);
      unmatchedTotal = Number(payload.total || (payload.view_counts || {}).unmatched || unmatchedRows.length);
      renderUnmatchedWorkspace();
      if (selectedUnmatchedResultIndex == null && unmatchedRows.length) {
        await selectUnmatched(Number(unmatchedRows[0].result_index), true);
      }
    } catch (error) {
      if (panel) panel.innerHTML = `<div class="empty-state"><div class="empty-icon">!</div><div class="empty-title">未匹配谱图加载失败</div><div class="empty-desc">${esc(error.message || error)}</div></div>`;
    } finally {
      unmatchedLoading = false;
    }
  }

  function renderUnmatchedWorkspace() {
    const panel = document.getElementById('unmatchedWorkbenchPanel');
    if (!panel) return;
    const funnel = currentBundle?.task?.identification_funnel || {};
    const withPeaks = Number(funnel.unmatched_with_nonzero_query_peaks || 0);
    panel.innerHTML = `
      <section class="review-workspace-intro unmatched-intro">
        <div><span class="panel-eyebrow">UNMATCHED SPECTRA WORKSPACE</span><h3>未匹配谱图分析工作台</h3>
          <p>未匹配只表示当前谱库与阈值下没有足够相似的参比谱，不表示谱图或样本成分被删除。MSClassify Top-3 是未校准的宽类别线索。</p></div>
        <div class="review-workspace-kpis"><span><b>${integer(unmatchedTotal)}</b> 未匹配 MS2</span><span><b>${integer(withPeaks)}</b> 仍有可用碎片</span><span><b>${integer(unmatchedRows.length)}</b> 已加载</span></div>
      </section>
      <div class="identification-workspace unmatched-workspace">
        <aside class="candidate-rail" aria-label="未匹配谱图列表">
          <div class="candidate-rail-tools"><label for="unmatchedSpectrumSearch">搜索已加载谱图</label>
            <input id="unmatchedSpectrumSearch" class="form-input" type="search" value="${esc(unmatchedSearch)}" placeholder="Scan、m/z、RT 或类别" oninput="EvidenceReview.setUnmatchedSearch(this.value)">
            <div id="unmatchedLoadedCount" class="candidate-rail-count"></div>
          </div>
          <div id="unmatchedSpectrumList" class="candidate-list unmatched-list"></div>
          ${unmatchedRows.length < unmatchedTotal ? `<button class="btn btn-secondary unmatched-load-more" onclick="EvidenceReview.loadMoreUnmatched(false)">再加载 50 张（${integer(unmatchedRows.length)}/${integer(unmatchedTotal)}）</button>` : ''}
        </aside>
        <main id="unmatchedEvidencePane" class="candidate-evidence-pane">
          <div class="empty-state"><div class="empty-icon">≋</div><div class="empty-title">请选择一张未匹配谱图</div><div class="empty-desc">这里显示原始碎片、Top-3 类别和逐谱定位</div></div>
        </main>
      </div>`;
    renderUnmatchedList();
  }

  function renderUnmatchedList() {
    const list = document.getElementById('unmatchedSpectrumList');
    const count = document.getElementById('unmatchedLoadedCount');
    if (!list) return;
    const query = unmatchedSearch.trim().toLowerCase();
    const rows = unmatchedRows.filter(row => {
      if (!query) return true;
      return [row.spectrum_id, row.precursor_mz, row.retention_time, row.display_class,
        ...(row.dl_top3 || []).flat()].filter(value => value != null).join(' ').toLowerCase().includes(query);
    });
    if (count) count.textContent = `显示 ${integer(rows.length)} / ${integer(unmatchedRows.length)} 张已加载谱图`;
    list.innerHTML = rows.map(row => {
      const resultIndex = Number(row.result_index);
      const probability = Number(row.dl_confidence || 0);
      return `<button class="candidate-list-item unmatched-list-item${resultIndex === selectedUnmatchedResultIndex ? ' selected' : ''}" data-unmatched-result="${resultIndex}" onclick="EvidenceReview.selectUnmatched(${resultIndex})">
        <span class="candidate-list-top"><b>MS2 ${esc(row.spectrum_id || '—')}</b><em>${integer(row.peaks_count)} 原始峰</em></span>
        <span>m/z ${Number(row.precursor_mz || 0).toFixed(5)} · RT ${row.retention_time == null ? '—' : Number(row.retention_time).toFixed(3)} min</span>
        <span>${esc(row.display_class || '未能判类')}${probability ? ` · 模型概率 ${(probability * 100).toFixed(1)}%` : ''}</span>
      </button>`;
    }).join('') || '<div class="candidate-list-empty">当前搜索条件下没有已加载谱图。</div>';
  }

  function setUnmatchedSearch(value) {
    unmatchedSearch = String(value || '');
    renderUnmatchedList();
  }

  function flagshipEntityForResultIndex(resultIndex) {
    const wanted = Number(resultIndex);
    if (!Number.isInteger(wanted)) return null;
    return (currentBundle?.flagship?.entities || []).find(entity =>
      (entity?.member_result_indices || []).some(index => Number(index) === wanted)
    ) || null;
  }

  function retrievalOnlyCandidateSupport(row, resultIndex) {
    const entity = flagshipEntityForResultIndex(resultIndex ?? row?.result_index);
    const nativeCandidates = Array.isArray(entity?.candidates) ? entity.candidates : [];
    const support = entity?.candidate_recall_support || {};
    const candidates = Array.isArray(support.expanded_candidates)
      ? support.expanded_candidates : [];
    if (nativeCandidates.length || !candidates.length) return '';

    const candidateRows = candidates.slice(0, 40).map(item => {
      const references = Object.values(item.best_reference_records || {})
        .filter(value => value && typeof value === 'object');
      const bestReference = references[0] || {};
      const databaseName = String(item.candidate?.compound_name || 'Unnamed database candidate');
      return `<tr>
        <td>${integer(item.retrieval_rank)}</td>
        <td><b data-scientific-value>${esc(databaseName)}</b><br><span class="mono" data-scientific-value>${esc(item.connectivity_identity || item.identity_key || '—')}</span></td>
        <td class="mono" data-scientific-value>${esc(item.candidate?.formula || '—')}</td>
        <td class="mono">${decimal(item.retrieval_rrf_utility, 6)}</td>
        <td>${esc((item.spectral_view_support || []).join(', ') || '—')}</td>
        <td data-scientific-value>${esc(bestReference.dataset_source || '—')}<br><span class="mono" data-scientific-value>${esc(bestReference.source_accession || bestReference.spectrum_id || '—')}</span></td>
      </tr>`;
    }).join('');
    const count = candidates.length;
    if (isEnglishInterface()) {
      return `<section class="pw-card unmatched-retrieval-support" style="padding:.75rem;margin-top:1rem;border-left:5px solid #2563eb;">
        <h4>Spectral-retrieval candidate hypotheses — no primary result</h4>
        <p>V15 has no native candidate for this precursor event. The independent Top-40 retrieval lane retained <b>${integer(count)}</b> database candidates for review. They do not replace a winner, do not constitute a unique 2-D structure assignment, and add zero to the compound count. Database names below are displayed exactly as stored and are never translated by the interface language. RRF utility is a retrieval score, not a probability.</p>
        <details><summary>Review ${integer(count)} retrieval candidates</summary><div class="table-wrap"><table><thead><tr><th>Retrieval rank</th><th>Database candidate</th><th>Formula</th><th>RRF utility</th><th>Spectral views</th><th>Best reference</th></tr></thead><tbody>${candidateRows}</tbody></table></div></details>
      </section>`;
    }
    return `<section class="pw-card unmatched-retrieval-support" style="padding:.75rem;margin-top:1rem;border-left:5px solid #2563eb;">
      <h4>谱库检索候选假设——暂无主要结果</h4>
      <p>V15 对这个母离子事件没有原生候选；独立 Top‑40 检索保留了 <b>${integer(count)}</b> 个数据库候选供复核。它们不能替换赢家，不构成唯一二维结构指认，化合物计数增加 0。下表名称严格使用数据库原名，不随界面语言翻译。RRF 效用是检索分，不是概率。</p>
      <details><summary>查看 ${integer(count)} 个检索候选</summary><div class="table-wrap"><table><thead><tr><th>检索名次</th><th>数据库候选</th><th>分子式</th><th>RRF 效用</th><th>谱图视图</th><th>最佳参考</th></tr></thead><tbody>${candidateRows}</tbody></table></div></details>
    </section>`;
  }

  async function selectUnmatched(resultIndex, updateHash = true) {
    const wanted = Number(resultIndex);
    if (!currentTaskId || !Number.isInteger(wanted)) return;
    const renderSequence = ++unmatchedDetailRenderSequence;
    const renderTaskId = currentTaskId;
    selectedUnmatchedResultIndex = wanted;
    document.querySelectorAll('[data-unmatched-result]').forEach(button => {
      button.classList.toggle('selected', Number(button.dataset.unmatchedResult) === wanted);
    });
    const pane = document.getElementById('unmatchedEvidencePane');
    if (pane) pane.innerHTML = '<div class="sample-loading"><span class="spinner"></span><b>正在加载原始未匹配谱图…</b></div>';
    if (updateHash) updateHashState();
    try {
      const payload = await API.get(`/tasks/${encodeURIComponent(renderTaskId)}/results/item/${wanted}`);
      if (renderSequence !== unmatchedDetailRenderSequence || renderTaskId !== currentTaskId) return;
      const row = payload.result || {};
      renderUnmatchedDetail(row, wanted);
    } catch (error) {
      if (renderSequence !== unmatchedDetailRenderSequence || renderTaskId !== currentTaskId) return;
      if (pane) pane.innerHTML = `<div class="empty-state"><div class="empty-icon">!</div><div class="empty-title">谱图加载失败</div><div class="empty-desc">${esc(error.message || error)}</div></div>`;
    }
  }

  function renderUnmatchedDetail(row, resultIndex = null) {
    const pane = document.getElementById('unmatchedEvidencePane');
    if (!pane) return;
    const renderSequence = unmatchedDetailRenderSequence;
    const renderTaskId = currentTaskId;
    const english = isEnglishInterface();
    const copy = (zh, en) => english ? en : zh;
    const finite = value => !['number', 'string'].includes(typeof value)
      || (typeof value === 'string' && !value.trim()) || !Number.isFinite(Number(value)) ? null : Number(value);
    const savedPeaks = Array.isArray(row.query_peaks) ? row.query_peaks
      : Array.isArray(row.peaks) ? row.peaks : [];
    const peakRecords = savedPeaks.map((peak, index) => ({
      sourceIndex: index,
      mz: Array.isArray(peak) ? finite(peak[0]) : null,
      intensity: Array.isArray(peak) ? finite(peak[1]) : null,
    })).filter(peak => peak.mz != null && peak.mz > 0 && peak.intensity != null && peak.intensity >= 0);
    const peaks = peakRecords.map(peak => [peak.mz, peak.intensity]);
    const top3 = Array.isArray(row.dl_top3) ? row.dl_top3.slice(0, 3) : [];
    const probability = finite(row.dl_confidence);
    const probabilityText = probability != null && probability >= 0 && probability <= 1
      ? `${(probability * 100).toFixed(1)}%` : '—';
    const precursor = finite(row.precursor_mz);
    const precursorText = precursor != null && precursor > 0 ? precursor.toFixed(5) : copy('未记录', 'not recorded');
    const rt = SampleResultData.retentionTimeMinutes(row);
    const topRows = top3.map((item, index) => {
      const value = finite(item?.[1]);
      const available = value != null && value >= 0 && value <= 1;
      const label = english ? interfaceCopy(item?.[0], 'Unclassified category') : (item?.[0] || '未能判类');
      return `<div class="unmatched-class-option"><div><span>${index + 1}. ${esc(label)}</span><b>${available ? `${(value * 100).toFixed(2)}%` : '—'}</b></div><div class="unmatched-class-track" aria-hidden="true"><span style="width:${available ? value * 100 : 0}%"></span></div></div>`;
    }).join('');
    const peakRows = peakRecords.slice(0, 300).map(peak => `<tr><td>${peak.sourceIndex + 1}</td><td class="mono">${peak.mz.toFixed(5)}</td><td class="mono">${peak.intensity.toFixed(2)}</td></tr>`).join('');
    const rawReason = String(currentBundle?.task?.identification_funnel?.unmatched_reason_scope || '');
    const unmatchedReason = rawReason.startsWith('No library candidate passed configured retrieval and score thresholds')
      ? '在当前同极性谱库、检索范围和评分阈值下，没有候选通过；没有谱库匹配不等于没有 MS1 色谱特征。'
      : (rawReason || '当前同极性谱库和配置阈值下没有候选通过检索与评分门。');
    const retrievalSupport = retrievalOnlyCandidateSupport(row, resultIndex);
    const standardReason = 'No candidate passed the same-polarity library retrieval and configured score thresholds.';
    const resultJsonLink = Number.isInteger(resultIndex)
      ? `<a class="btn btn-secondary btn-sm" href="${API.BASE}/tasks/${encodeURIComponent(renderTaskId)}/results/item/${resultIndex}" target="_blank" rel="noopener">${copy('查看该谱全部保存数据（JSON）', 'View all saved spectrum data (JSON)')}</a>` : '';
    pane.innerHTML = `<article class="unmatched-evidence-detail">
      <header><div><span class="panel-eyebrow">MEASURED MS/MS</span><h3>${copy('未匹配谱图', 'Unmatched spectrum')} · MS2 <span data-scientific-value>${esc(row.spectrum_id || '—')}</span></h3><p><span data-scientific-value>${esc(currentBundle?.task?.filename || copy('未记录源文件', 'Source file not recorded'))}</span> · RT ${rt == null ? copy('未记录', 'not recorded') : rt.toFixed(3) + ' min'} · m/z ${precursorText} · ${integer(peaks.length)} / ${integer(savedPeaks.length)} ${copy('个保存峰可绘制', 'saved peaks drawable')}</p></div><span class="sample-evidence-badge evidence-low">${copy('类别模型', 'Class model')} · ${probabilityText}</span></header>
      <div class="unmatched-visual-grid">
        <section class="result-spectrum-panel"><div class="result-visual-heading"><span>${copy('实测 MS/MS 谱图', 'Measured MS/MS spectrum')}</span><small>${copy('峰顶悬停查看 m/z 与强度；拖动放大', 'Hover for m/z and intensity; drag to zoom')}</small></div><div id="reviewUnmatchedSpectrum" class="result-spectrum-stage"></div></section>
        <section class="unmatched-class-card"><div class="result-visual-heading"><span>${copy('类别模型 Top-3', 'Class model Top-3')}</span></div><div class="unmatched-class-bars">${topRows || `<div class="sample-no-results">${copy('未保存类别概率', 'No class probabilities saved')}</div>`}</div><p class="unmatched-class-note">${copy('分类模型输出，未经化合物鉴定置信度校准。', 'Class-model outputs; uncalibrated for compound identification.')}</p></section>
      </div>
      <div class="unmatched-reason-note"><b>${copy('谱库检索状态：', 'Library retrieval: ')}</b>${esc(english ? interfaceCopy(rawReason, standardReason) : unmatchedReason)}</div>
      ${retrievalSupport}
      <details class="result-evidence-details unmatched-peak-details"><summary><span>${copy('查看保存峰数据', 'Inspect saved peak data')}</span><small>${copy('当前显示', 'Showing')} ${Math.min(300, peakRecords.length)} / ${peakRecords.length}</small></summary><div class="result-evidence-details-body"><div class="table-wrap unmatched-peak-table"><table><thead><tr><th>${copy('原峰序号', 'Source peak')}</th><th>m/z</th><th>${copy('保存强度', 'Stored intensity')}</th></tr></thead><tbody>${peakRows || `<tr><td colspan="3">${copy('没有可绘制峰', 'No drawable peaks')}</td></tr>`}</tbody></table></div>${resultJsonLink}</div></details>
    </article>`;
    setTimeout(() => {
      if (renderSequence !== unmatchedDetailRenderSequence || renderTaskId !== currentTaskId) return;
      const element = document.getElementById('reviewUnmatchedSpectrum');
      if (!element) return;
      if (peaks.length && typeof SpectrumComparator !== 'undefined') {
        new SpectrumComparator().render('reviewUnmatchedSpectrum',
          { peaks, name: copy('实测查询谱', 'Measured query') }, null, [],
          { title: copy('实测 MS/MS', 'Measured MS/MS'), height: 330, text: spectrumPlotText() });
      } else if (peaks.length && typeof SpectrumViewer !== 'undefined') {
        new SpectrumViewer('reviewUnmatchedSpectrum').render(peaks, { title: copy('实测 MS/MS', 'Measured MS/MS'), height: 330 });
      } else {
        element.innerHTML = `<div class="result-visual-empty">${copy('未保存可绘制的 MS/MS 碎片峰', 'No drawable MS/MS fragment peaks saved')}</div>`;
      }
    }, 80);
  }

  async function open(taskId, tab = 'identification', compoundSearch = '', resultIndex = null, updateHash = true) {
    if (!taskId) return;
    programmaticOpen = true;
    try {
      App.navigateTo('review');
      await refreshSelect(taskId);
      const select = document.getElementById('reviewTaskSelect');
      if (select) select.value = taskId;
      showTab(tab, false);
      await loadTask(taskId, resultIndex);
      if (compoundSearch) {
        candidateSearch = String(compoundSearch);
        const input = document.getElementById('reviewCandidateSearch');
        if (input) input.value = candidateSearch;
        renderCandidateList();
      }
      showTab(tab, false);
      if (updateHash) updateHashState();
    } finally {
      programmaticOpen = false;
    }
  }

  function updateHashState() {
    if (!currentTaskId || restoringHash) return;
    const params = new URLSearchParams({ task: currentTaskId, tab: activeTab });
    if (activeTab === 'identification' && Number.isInteger(selectedResultIndex)) {
      params.set('result', String(selectedResultIndex));
    }
    if (activeTab === 'unmatched' && Number.isInteger(selectedUnmatchedResultIndex)) {
      params.set('result', String(selectedUnmatchedResultIndex));
    }
    history.replaceState(null, '', `#review?${params.toString()}`);
  }

  async function restoreHashState() {
    if (restoringHash) return true;
    const raw = window.location.hash.replace(/^#/, '');
    const [page, query = ''] = raw.split('?');
    if (page !== 'review' || !query) return false;
    const params = new URLSearchParams(query);
    const taskId = params.get('task');
    const tab = params.get('tab') || 'identification';
    const rawResultIndex = params.get('result');
    const resultIndex = rawResultIndex == null || rawResultIndex === '' ? null : Number(rawResultIndex);
    if (!taskId) return false;
    restoringHash = true;
    try {
      if (currentTaskId !== taskId || !currentBundle) {
        await open(taskId, tab, '', Number.isInteger(resultIndex) ? resultIndex : null, false);
      } else {
        showTab(tab, false);
        if (tab === 'identification' && Number.isInteger(resultIndex)) {
          await selectCandidate(resultIndex, false);
        } else if (tab === 'unmatched') {
          await ensureUnmatchedLoaded();
          if (Number.isInteger(resultIndex)) await selectUnmatched(resultIndex, false);
        }
      }
    } finally {
      restoringHash = false;
      if (currentTaskId) updateHashState();
    }
    return true;
  }

  async function enterPage() {
    if (programmaticOpen) return;
    const [, query = ''] = window.location.hash.replace(/^#/, '').split('?');
    if (!query) activeTab = 'identification';
    const hashTaskId = query ? new URLSearchParams(query).get('task') : null;
    const overviewTaskId = typeof SampleOverview !== 'undefined' ? SampleOverview.currentTaskId : null;
    await refreshSelect(hashTaskId || overviewTaskId || null);
    if (await restoreHashState()) return;
    const select = document.getElementById('reviewTaskSelect');
    const taskId = overviewTaskId || select?.value || currentTaskId;
    if (!taskId) return;
    if (currentTaskId !== taskId || !currentBundle) {
      await loadTask(taskId);
    } else {
      showTab(activeTab, false);
    }
    updateHashState();
  }

  async function backToResults() {
    if (currentTaskId && typeof SampleOverview !== 'undefined') {
      await SampleOverview.openTask(currentTaskId);
      return;
    }
    App.navigateTo('results');
  }

  return {
    init, refreshSelect, showTab, loadTask, open, restoreHashState, enterPage, backToResults,
    setCandidateFilter, setCandidateSearch, selectCandidate,
    loadMoreUnmatched, setUnmatchedSearch, selectUnmatched, replicateMissingnessSummary,
    get currentTaskId() { return currentTaskId; },
  };
})();
