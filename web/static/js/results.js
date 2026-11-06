/**
 * MSsystem — Results Viewer Module
 * Task result viewing with summary, match table, detail modal
 */
const Results = (() => {
  let _currentTaskId = null;
  let _currentMatched = [];   // 缓存当前任务的已排序matched列表, showDetail 直接索引，避免与API顺序错位 (修复列表vs详情化合物不一致bug)
  let _currentUnmatched = []; // 未匹配谱图也必须可审查：缓存用于类别Top-3/碎片详情弹窗
  let _currentUnmatchedTotal = 0;
  let _unmatchedDisplayLimit = 50;
  let _currentIonMode = null;
  let _currentTaskFilename = '';
  let _blankView = 'all';
  let _detailRenderSequence = 0;

  function _retentionTimeMinutes(row) {
    const parsed = Number(row?.retention_time);
    if (!Number.isFinite(parsed)) return null;
    const unit = String(row?.retention_time_unit || 'minutes').trim().toLowerCase();
    if (unit.includes('millisecond') || ['ms', 'msec', 'msecs'].includes(unit)) return parsed / 60000;
    if (unit.includes('second') || ['s', 'sec', 'secs'].includes(unit)) return parsed / 60;
    if (unit.includes('hour') || ['h', 'hr', 'hrs'].includes(unit)) return parsed * 60;
    return parsed;
  }

  function _precursorMassAxis(row) {
    const selectedIonMz = Number(row?.precursor_mz);
    const searchPrecursorMz = Number(row?.search_precursor_mz);
    const searchSource = String(row?.search_precursor_mz_source || '');
    const featureUsed = row?.lcms_feature_used_for_search === true
      || searchSource === 'measured_continuous_ms1_feature_centroid';
    const deltaValue = row?.search_precursor_delta_ppm;
    return {
      selectedIonMz: Number.isFinite(selectedIonMz) ? selectedIonMz : null,
      searchPrecursorMz: Number.isFinite(searchPrecursorMz)
        ? searchPrecursorMz
        : (Number.isFinite(selectedIonMz) ? selectedIonMz : null),
      searchSource,
      featureUsed,
      deltaPpm: deltaValue == null || deltaValue === ''
        ? null : Number(deltaValue),
    };
  }

  function _precursorMassCell(row, digits = 4) {
    const mass = _precursorMassAxis(row);
    const selected = mass.selectedIonMz == null
      ? '—' : mass.selectedIonMz.toFixed(digits);
    if (!mass.featureUsed || mass.searchPrecursorMz == null) return selected;
    const delta = Number.isFinite(mass.deltaPpm)
      ? ` · Δ ${mass.deltaPpm.toFixed(2)} ppm` : '';
    return `<span title="mzML selected-ion / 四极杆隔离目标">${selected}</span>`
      + `<small style="display:block;color:#166534;white-space:nowrap;" title="连续 MS1 特征实测质心；该质量用于 HRMS 候选检索">MS1检索 ${mass.searchPrecursorMz.toFixed(digits)}${delta}</small>`;
  }

  function _precursorMassInline(row, digits = 5) {
    const mass = _precursorMassAxis(row);
    const selected = mass.selectedIonMz == null
      ? '隔离目标未记录' : `MS2 隔离目标 m/z ${mass.selectedIonMz.toFixed(digits)}`;
    if (!mass.featureUsed || mass.searchPrecursorMz == null) {
      return `${selected}（检索质量轴回退）`;
    }
    const delta = Number.isFinite(mass.deltaPpm)
      ? `，Δ ${mass.deltaPpm.toFixed(2)} ppm` : '';
    return `${selected} · 检索用实测 MS1 m/z ${mass.searchPrecursorMz.toFixed(digits)}${delta}`;
  }

  function _evographFragmentationHtml(graph) {
    if (!graph || typeof graph !== 'object') return '';
    const esc = value => Tasks.escapeHtml(String(value == null ? '—' : value));
    const number = (value, digits = 4) => {
      const parsed = Number(value);
      return Number.isFinite(parsed) ? parsed.toFixed(digits) : '—';
    };
    const context = graph.context || {};
    const audit = graph.graph_audit || {};
    const nodes = Array.isArray(graph.nodes) ? graph.nodes : [];
    const nodeById = new Map(nodes.map(node => [String(node.node_id || ''), node]));
    const fragmentRows = nodes
      .filter(node => node.node_type === 'fragment')
      .slice(0, 16)
      .map(node => {
        const candidates = Array.isArray(node.formula_candidates)
          ? node.formula_candidates : [];
        const best = candidates[0] || null;
        const status = String(node.incoming_edge_status || 'unassigned');
        const statusLabel = status === 'heuristically_supported' ? '启发支持'
          : status === 'ambiguous' ? '多解'
            : '未归属';
        return `<tr>
          <td class="mono">${number(node.mz, 5)}</td>
          <td class="mono">${best ? esc(best.formula) : '—'}</td>
          <td class="mono">${best ? number(best.mass_error_ppm, 2) : '—'}</td>
          <td>${esc(statusLabel)}</td>
          <td>${candidates.length > 1 ? `${candidates.length} 个候选组成` : (best ? '唯一保留组成' : '质量/母式约束下无解')}</td>
        </tr>`;
      }).join('');
    const edgeRows = (Array.isArray(graph.edges) ? graph.edges : [])
      .filter(edge => edge.selected && edge.kind === 'formula_conserving')
      .slice(0, 12)
      .map(edge => {
        const parent = nodeById.get(String(edge.parent_id || '')) || {};
        const child = nodeById.get(String(edge.child_id || '')) || {};
        const prior = edge.common_loss_prior || {};
        const selection = edge.selection_status === 'heuristically_supported'
          ? '启发支持' : '多解/保留';
        return `<tr>
          <td class="mono">${number(parent.mz, 5)} → ${number(child.mz, 5)}</td>
          <td class="mono">${esc(edge.parent_formula)} → ${esc(edge.child_formula)}</td>
          <td class="mono">−${number(edge.observed_loss_da, 5)}</td>
          <td class="mono">${esc(edge.loss_formula)}</td>
          <td class="mono">${number(edge.mass_error_ppm, 2)}</td>
          <td>${esc(prior.name || selection)}</td>
        </tr>`;
      }).join('');
    const formulaAvailable = graph.formula_graph_available === true;
    const statusMessage = formulaAvailable
      ? `已为 ${Number(audit.fragment_node_count || 0)} 个观测碎片检查母式子集组成；${Number(audit.selected_edge_count || 0)} 条边被保留。`
      : `未生成分子式图：${esc((context.reasons || []).join('、') || graph.reason || '候选分子式/明确加合物/电荷条件不完整')}。`;
    const selectedMz = number(graph.selected_ion_mz, 6);
    const searchMz = number(graph.search_precursor_mz ?? graph.precursor_mz, 6);
    const parentError = number(context.mass_consistency_error_ppm, 2);
    return `<details class="pw-card" open>
      <summary>🧬 观测碎片组成图（EvoGraph-F）</summary>
      <div style="padding:.5rem .65rem;font-size:.72rem;line-height:1.5;">
        <div style="padding:.42rem .55rem;border:1px solid ${formulaAvailable ? '#a7f3d0' : '#fed7aa'};background:${formulaAvailable ? '#ecfdf5' : '#fff7ed'};border-radius:5px;">
          <b>${formulaAvailable ? '已完成候选条件下的分子式守恒检查' : '当前组成解释不可用'}</b> · ${statusMessage}<br>
          隔离目标 m/z ${selectedMz} · 检索/计算质量 m/z ${searchMz} · 母离子式质量误差 ${parentError} ppm
        </div>
        ${fragmentRows ? `<details class="pw-det" open><summary>观测产物离子与候选组成</summary><table class="pw-tbl"><thead><tr><th>观测 m/z</th><th>候选中性组成</th><th>误差 ppm</th><th>入边状态</th><th>组成歧义</th></tr></thead><tbody>${fragmentRows}</tbody></table></details>` : ''}
        ${edgeRows ? `<details class="pw-det" open><summary>分子式守恒的碎裂边</summary><table class="pw-tbl"><thead><tr><th>观测路径</th><th>母式 → 子式</th><th>观测 Δm</th><th>丢失组成</th><th>误差 ppm</th><th>注释</th></tr></thead><tbody>${edgeRows}</tbody></table></details>` : ''}
        <div style="margin-top:.4rem;padding:.42rem .55rem;background:#f8fafc;border:1px solid #cbd5e1;border-radius:5px;color:#475569;">
          边界：这些分子式是以当前谱库候选、明确加合物和电荷保持假设为条件的<b>组成假设</b>，不是已确认的产物结构或键断裂机理；启发权重不是概率。本模块不参与候选检索、排序、置信度或 High 升档。若有 CFM/RDKit 二维产物结构，将在独立附加区展示。
        </div>
      </div>
    </details>`;
  }

  function _hasStructure(match) {
    if (match.has_structure !== undefined) return !!match.has_structure;
    const value = String(match.smiles || '').trim().toLowerCase();
    return !!value && !['nan', 'none', 'null', 'n/a'].includes(value);
  }

  function _normalizeInChIKey(value) {
    const key = String(value || '').trim().toUpperCase().replace(/^INCHIKEY=/, '').replace(/\s+/g, '');
    return /^[A-Z]{14}-[A-Z]{10}-[A-Z]$/.test(key) ? key : '';
  }

  function _compoundIdentityKey(match) {
    const ik = _normalizeInChIKey(match.inchikey);
    if (ik) return `ik::${ik}`;
    const smiles = String(match.smiles || '').trim();
    if (smiles) return `sm::${smiles}`;
    const name = String(match.compound_name || '').toLowerCase().replace(/[^a-z0-9]+/g, '');
    const formula = String(match.formula || '').trim();
    return `nf::${name}|${formula}`;
  }

  function _compareMatched(a, b) {
    const am = a.match || {}, bm = b.match || {};
    const num = (m, key, fallback = 0) => Number(m[key] ?? fallback) || 0;
    const ppm = (m) => m.mass_ppm == null ? 999999 : Math.abs(Number(m.mass_ppm) || 0);
    const numeric = [
      [num(bm, 'rank_score', num(bm, 'final_score')), num(am, 'rank_score', num(am, 'final_score'))],
      [num(bm, 'final_score'), num(am, 'final_score')],
      [num(bm, 'cosine'), num(am, 'cosine')],
      [num(bm, 'reverse_match'), num(am, 'reverse_match')],
      [ppm(am), ppm(bm)],
    ];
    for (const [left, right] of numeric) {
      if (left !== right) return left - right;
    }
    const aKey = `${_compoundIdentityKey(am)}|${am.dataset_source || ''}|${am.lib_spectrum_id || ''}|${a.spectrum_id || ''}`;
    const bKey = `${_compoundIdentityKey(bm)}|${bm.dataset_source || ''}|${bm.lib_spectrum_id || ''}|${b.spectrum_id || ''}`;
    return aKey.localeCompare(bKey);
  }

  function _classificationPresentation(match) {
    const display = String(match.display_class || '').trim();
    const raw = String(match.classification || '').trim();
    const displayUsable = display && !['unknown', 'other', 'unclassified'].includes(display.toLowerCase());
    const rawUsable = raw && !['unknown', 'other', 'unclassified'].includes(raw.toLowerCase());
    const label = displayUsable ? display : (rawUsable ? raw : '');
    if (!label) {
      const identityHigh = [
        match?.confidence,
        match?.exact_structure_confidence,
        match?.two_dimensional_connectivity_confidence,
      ].some(value => String(value || '').trim().toLowerCase() === 'high');
      return {
        label: identityHigh
          ? '尚无可靠类别归属（谱库候选证据为 High）'
          : (_hasStructure(match) ? '未能可靠判类' : '不可分类（无结构）'),
        reliable: false,
        hasRaw: false,
      };
    }
    const source = String(
      (displayUsable ? match.display_classification_source : '')
      || match.classification_source || ''
    ).trim().toLowerCase();
    const provenance = match.expert_details?.classification_provenance || {};
    const provenanceReliable = match.expert_details?.classification_provenance?.reliable;
    const provenanceMatches = provenanceReliable === true
      && String(provenance.analysis_class || '').trim().toLowerCase() === label.toLowerCase();
    const trustedSource = ['curated', 'manual', 'smarts', 'structure', 'rule', 'database_curated']
      .some(prefix => source === prefix || source.startsWith(`${prefix}_`))
      || source.startsWith('runtime_structure_rule');
    const untrustedSource = !source || ['unknown', 'other', 'existing', 'ml', 'msclass', 'msclass_spectrum', 'fragment', 'fragment_inference', 'spectrum'].includes(source)
      || ['ml_', 'model_', 'msclass_', 'fragment_'].some(prefix => source.startsWith(prefix));
    if (!untrustedSource && (match.classification_reliable === true || provenanceMatches || trustedSource)) {
      return { label, reliable: true, hasRaw: true, source };
    }
    const modelSource = source === 'ml' || source === 'msclass_spectrum'
      || source.startsWith('ml_') || source.startsWith('model_') || source.startsWith('msclass_');
    if (modelSource) {
      return { label: `类别待核查（模型提示：${label}）`, reliable: false, hasRaw: true, source };
    }
    return { label: `类别线索（待核查）：${label}`, reliable: false, hasRaw: true, source };
  }

  function _classificationLabel(match) {
    return _classificationPresentation(match).label;
  }

  function _confidenceLabel(value, match) {
    const raw = String(value || '').trim();
    if (raw && raw.toLowerCase() !== 'unknown') return raw;
    return _hasStructure(match) ? '证据不足' : '仅谱图线索';
  }

  const _governanceBlockerLabels = {
    cosine_below_locked_release_threshold: 'MS/MS 相似度未达到正式报告要求',
    final_score_below_locked_release_threshold: '综合匹配分未达到正式报告要求',
    expert_not_applicable: '当前结构缺少可用的特征碎片规则',
    scoped_expert_support_below_threshold: '可解释的诊断碎片或中性丢失支持不足',
    no_validated_structure: '参比记录缺少已验证结构',
    library_source_not_experimental: '参比来源不是合格实验谱',
    precursor_formula_adduct_conflict: '前体质量、分子式或加合物不一致',
    mass_error_above_locked_release_limit: '前体质量误差偏大',
    severe_library_quality_flag: '参比谱质量存在重要警告',
    blank_background_conflict: '候选与空白背景存在冲突',
    potential_external_contaminant: '候选具有潜在外源用途',
    declared_blank_query: '当前输入被标记为空白样',
    scope_downgraded_due_to_unresolved_connectivity_ambiguity: '近邻结构尚未排除，只能报告到候选组',
    project_level_qvalue_unavailable: '项目级 q 值校准尚未建立',
  };

  const _reportScopeLabels = {
    authentic_standard: '同方法标准品确认',
    exact_structure: '具体结构候选',
    '2d_connectivity': '具体化合物候选（连接结构）',
    isomer_group: '分子式异构体组',
    molecular_formula_isomer_group: '分子式异构体组',
    molecular_formula: '分子式',
    chemical_class: '化学类别',
    unresolved: '谱库候选',
  };

  const _projectLevelReasonCodes = new Set([
    'project_level_qvalue_unavailable',
    'project_level_target_decoy_qvalue_not_yet_available',
    'insufficient_independent_injections',
  ]);

  function _reportScopeLabel(value) {
    const key = String(value || 'unresolved');
    return _reportScopeLabels[key] || '谱库候选';
  }

  function _candidateReasonCodes(match) {
    const decision = match.final_claim_decision || {};
    const passport = match.evidence_governance || {};
    return [...new Set([
      ...(decision.reason_codes || []),
      ...(passport.release_blockers || []),
    ].map(String))].filter(code => !_projectLevelReasonCodes.has(code));
  }

  function _candidateNextStep(match, reasonCodes) {
    const scope = String((match.final_claim_decision || {}).scope
      || (match.evidence_governance || {}).claim_scope || '');
    const reasons = new Set(reasonCodes);
    if (scope.includes('isomer') || match.isomer_ambiguous
      || reasons.has('scope_downgraded_due_to_unresolved_connectivity_ambiguity')) {
      return '优先使用同方法标准品 RT；必要时结合差异碎片、CCS 或离子迁移区分异构体。';
    }
    if (reasons.has('precursor_formula_adduct_conflict')
      || reasons.has('mass_error_above_locked_release_limit')) {
      return '先核对前体离子、加合物和质量校准，再复查或复采目标 MS/MS。';
    }
    if (reasons.has('cosine_below_locked_release_threshold')
      || reasons.has('final_score_below_locked_release_threshold')
      || reasons.has('scoped_expert_support_below_threshold')) {
      return '建议用多碰撞能、较窄隔离窗的靶向 MS/MS 补充可区分碎片。';
    }
    if (reasons.has('library_source_not_experimental')
      || reasons.has('no_validated_structure')) {
      return '补充合格实验参比谱，或使用同方法标准品进行确认。';
    }
    if (reasons.has('blank_background_conflict')
      || reasons.has('potential_external_contaminant')) {
      return '加入同极性空白并比较样品/空白强度与批次复现情况。';
    }
    if ((match.final_claim_decision || {}).released) {
      return '按项目 SOP 归档原始谱图定位、校准记录与确认依据。';
    }
    return '按页面顶部说明补充独立进样，并完成项目级正式报告校准。';
  }

  const _evidenceNodeLabels = {
    spectral_similarity: 'MS/MS 谱库相似度',
    precursor_mass: '前体质量与加合物',
    library_provenance: '参比谱来源',
    scoped_expert_reasoning: '适用范围内的碎裂解释',
    candidate_competition: '候选与异构体竞争',
    blank_contrast: '空白对照',
    local_rt_standard: '同方法标准品保留时间',
    within_injection_stability: '单次进样内多能量一致性',
  };

  const _evidenceNodeStatusLabels = {
    supporting: '支持', contradictory: '存在冲突', insufficient: '支持不足',
    not_applicable: '不适用', unavailable: '不可用', not_run: '未执行',
    scope_limited: '报告层级受限', unevaluable: '不可评价',
    context_warning: '背景提示', single_observation: '仅一次观测',
    conflicting: '存在冲突', role_veto: '样本角色不支持',
  };

  function _evidenceNodeExplanation(node) {
    const purpose = {
      spectral_similarity: '核对查询谱与参比谱的碎片相似度和覆盖度。',
      precursor_mass: '核对实测前体、候选分子式与加合物是否一致。',
      library_provenance: '核对参比谱是否来自可追溯的实验来源。',
      scoped_expert_reasoning: '检查诊断离子和中性丢失是否在适用范围内支持该候选。',
      candidate_competition: '检查是否仍有分数接近或同分子式的竞争候选。',
      blank_contrast: '检查候选是否在同极性空白中重复出现。',
      local_rt_standard: '检查是否有同方法标准品保留时间支持。',
      within_injection_stability: '检查同一进样内不同碰撞能结果是否一致；这些扫描不算独立重复。',
    }[String(node.key || '')] || '记录该证据项目的适用性与结论。';
    return `${purpose} 当前结论：${_evidenceNodeStatusLabels[node.status] || node.status || '未记录'}。`;
  }

  function _finalClaimBadge(match) {
    const decision = match.final_claim_decision || {};
    if (!decision.policy_id) {
      return '<span class="formal-report-badge unassessed" title="尚未评估正式报告条件">报告条件未评估</span>';
    }
    const released = !!decision.released;
    const title = `正式报告条件：${released ? '已满足' : '未满足'}；当前候选可报告范围：${_reportScopeLabel(decision.scope)}；不改变候选排序和证据等级`;
    return `<span class="formal-report-badge ${released ? 'eligible' : 'candidate'}" title="${Tasks.escapeHtml(title)}">${released ? '报告条件已满足' : '报告条件未满足'}</span>`;
  }

  function _evidencePassportHtml(match) {
    const passport = match.evidence_governance || {};
    const finalDecision = match.final_claim_decision || {};
    if (!passport.policy_id && !finalDecision.policy_id) return '';
    const status = String(passport.claim_status || 'insufficient_for_release');
    const reasonCodes = _candidateReasonCodes(match);
    const mainReasons = reasonCodes.slice(0, 2).map(code =>
      _governanceBlockerLabels[code] || '存在需要进一步核查的候选级条件'
    );
    const released = !!finalDecision.released;
    const assessed = !!finalDecision.policy_id;
    const reportState = assessed
      ? (released ? '正式报告条件已满足' : '正式报告条件未满足')
      : '尚未评估正式报告条件';
    const stateClass = assessed ? (released ? 'eligible' : 'candidate') : 'unassessed';
    const scope = _reportScopeLabel(finalDecision.scope || passport.claim_scope);
    const nodeRows = (passport.nodes || []).map(node => {
      const good = node.status === 'supporting';
      const bad = ['contradictory', 'conflicting', 'role_veto'].includes(node.status);
      const stateColor = good ? '#166534' : bad ? '#991b1b' : '#64748b';
      const independent = node.independent_measurement
        ? '<span style="color:#166534;">独立观测</span>'
        : '<span style="color:#64748b;">同一数据派生或背景信息</span>';
      return `<tr>
        <td>${Tasks.escapeHtml(_evidenceNodeLabels[node.key] || node.label || node.key || '—')}</td>
        <td style="color:${stateColor};font-weight:700;">${Tasks.escapeHtml(_evidenceNodeStatusLabels[node.status] || node.status || '—')}</td>
        <td>${independent}</td>
        <td style="font-size:.67rem;color:var(--text-secondary);">${Tasks.escapeHtml(_evidenceNodeExplanation(node))}</td>
      </tr>`;
    }).join('');
    const ambiguity = (passport.ambiguity_set || []).map((candidate, index) =>
      `${index + 1}. ${candidate.compound_name || candidate.inchikey || 'candidate'}${candidate.formula ? ' [' + candidate.formula + ']' : ''}`
    ).join('；');
    const missingText = mainReasons.length
      ? mainReasons.join('；')
      : '未记录额外的候选级缺失条件；项目级要求见页面顶部。';
    const auditHtml = nodeRows || ambiguity ? `<details class="candidate-report-audit">
      <summary>技术审计信息（默认收起）</summary>
      ${nodeRows ? `<div class="table-wrap"><table><thead><tr><th>证据项目</th><th>结论</th><th>证据关系</th><th>说明</th></tr></thead><tbody>${nodeRows}</tbody></table></div>` : ''}
      ${ambiguity ? `<div class="candidate-report-audit-note"><b>当前近似候选集合：</b>${Tasks.escapeHtml(ambiguity)}<br>该集合用于提示近邻候选尚未完全排除，不代表统计覆盖率保证。</div>` : ''}
      <div class="candidate-report-audit-note">同一前体事件内的阶梯碰撞能扫描只算一次进样内观测，不能充当独立重复；当前候选排序和原置信字段均未被此检查改写。</div>
    </details>` : '';
    return `<section class="candidate-report-check">
      <header><b>正式报告条件（不改变候选鉴定）</b><span class="formal-report-badge ${stateClass}">${Tasks.escapeHtml(reportState)}</span></header>
      <div class="candidate-report-grid">
        <div><span>正式报告条件</span><b>${Tasks.escapeHtml(reportState)}</b></div>
        <div><span>正式报告允许范围</span><b>${Tasks.escapeHtml(scope)}</b></div>
        <div><span>主要缺失条件（最多两项）</span><b>${Tasks.escapeHtml(missingText)}</b></div>
        <div><span>建议下一步</span><b>${Tasks.escapeHtml(_candidateNextStep(match, reasonCodes))}</b></div>
      </div>
      ${auditHtml}
    </section>`;
  }

  function _pathwayData(value) {
    if (!value) return null;
    if (typeof value === 'object') return value;
    try { return JSON.parse(value); } catch (e) { return null; }
  }

  function _pathwayStatus(match, ionMode) {
    if (match.pathway_polarity_status === 'mismatch') return { state: 'mismatch', data: null };
    const data = _pathwayData(match.pathway_json);
    if (!data) return { state: 'missing', data: null };
    const wanted = String(ionMode || '').toLowerCase().startsWith('neg') ? 'neg'
      : String(ionMode || '').toLowerCase().startsWith('pos') ? 'pos' : '';
    const polarities = (data.polarities || []).map(p => String(p).toLowerCase());
    const normalized = new Set(polarities.map(p => p.startsWith('neg') ? 'neg' : p.startsWith('pos') ? 'pos' : p));
    if (wanted && normalized.size && !normalized.has(wanted)) return { state: 'mismatch', data };
    return { state: normalized.size ? 'available' : 'unknown_polarity', data };
  }

  function _analysisBadges(match, ionMode) {
    const ed = match.expert_details || {};
    const nl = ed.neutral_loss_analysis || {};
    const compactFlags = match.display_evidence_flags || null;
    const path = compactFlags
      ? { state: compactFlags.pathway_status || 'missing', data: null }
      : _pathwayStatus(match, ionMode);
    const hasNl = compactFlags
      ? !!compactFlags.has_neutral_loss
      : !!((nl.neutral_losses || []).length || (nl.sequential_losses || []).length);
    const hasFrag = compactFlags
      ? !!compactFlags.has_fragmentation
      : !!((nl.best_pathway || []).length
        || (nl.sequential_losses || []).length
        || (ed.structure_fragment_analysis?.plausible_fragments || []).length);
    const hasLibraryPathwayReference = compactFlags
      ? !!compactFlags.has_library_pathway_reference
      : path.state === 'available' || path.state === 'unknown_polarity';
    const v3 = ed.v3 || {};
    const dimensionAssessments = Object.values(ed.dimension_assessments || {})
      .filter(value => value && typeof value === 'object');
    const ruleEvaluated = dimensionAssessments.length > 0;
    const hasRuleSupport = dimensionAssessments.some(value => (
      value.status === 'scored' && Number(value.evidence_count || 0) > 0
    ));
    const badge = (label, ok, title) => `<span title="${Tasks.escapeHtml(title)}" style="display:inline-block;margin:1px 2px;padding:1px 5px;border-radius:8px;background:${ok ? '#dcfce7' : '#f1f5f9'};color:${ok ? '#166534' : '#64748b'};font-size:0.6rem;font-weight:700;white-space:nowrap;">${label}${ok ? '✓' : '—'}</span>`;
    const classPresentation = _classificationPresentation(match);
    const hasClass = classPresentation.reliable;
    let html = badge('分类', hasClass, hasClass ? '经策展或结构规则支持的类别结果' : (classPresentation.hasRaw ? `存在尚未验证的类别线索：${classPresentation.label}` : (_hasStructure(match) ? '已分析，未得到可靠类别' : '无结构谱图线索不能进行结构分类'))) +
      badge('中性丢失', hasNl, hasNl ? '已完成中性丢失分析' : '已分析，未提取到可解释中性丢失') +
      badge('样品碎裂', hasFrag, hasFrag ? '已从本次样品 MS/MS 得到碎裂步骤或候选产物结构' : '本次样品未重建出可审计的碎裂步骤') +
      badge('库碎裂参考', hasLibraryPathwayReference, hasLibraryPathwayReference ? '有同候选谱库重复碎裂特征；仅作参考，不代表本样品已重建途径' : '无可用的候选谱库碎裂参考') +
      badge('规则证据', hasRuleSupport, hasRuleSupport
        ? '逐维规则评价中观察到正向支持；详情可审计'
        : (ruleEvaluated
          ? '逐维规则评价已完成，本次未观察到正向支持'
          : '没有保存可审计的逐维规则记录'));
    if (path.state === 'mismatch') {
      html += '<span title="馆藏预存路径来自另一离子极性，本次未用于评分或展示" style="display:inline-block;margin:1px 2px;padding:1px 5px;border-radius:8px;background:#fef3c7;color:#92400e;font-size:0.6rem;font-weight:700;white-space:nowrap;">库路径极性不符</span>';
    }
    const deconv = match.chimeric_deconvolution || {};
    if (deconv.status === 'accepted_score_applied') {
      html += badge('去卷积', true, `嵌合谱 L1 去卷积通过；target fraction ${deconv.target_fraction ?? '—'}；分数增益 ${deconv.score_contribution ?? 0}`);
    } else if (['mixed', 'severe'].includes(String(deconv.isolation_status || ''))) {
      html += badge('去卷积', false, `隔离窗 ${deconv.isolation_status}；${deconv.status || '未通过去卷积守卫'}`);
    }
    const iin = match.ion_identity_evidence || {};
    if (iin.network_id) {
      html += badge('IIMN', !iin.identity_conflict && !match.ion_identity_adduct_conflict,
        `${iin.network_id} · ${iin.ion_identity || 'ion'} · neutral mass ${iin.neutral_mass ?? '—'}${iin.adduct_consistency ? ' · ' + iin.adduct_consistency : ''}`);
    }
    return html;
  }

  // ── Initialization ────────────────────────────────────────────────
  function init() {
    const sel = document.getElementById('reviewTaskSelect');
    if (sel && typeof EvidenceReview === 'undefined') {
      sel.addEventListener('change', () => loadForTask(sel.value));
    }
  }

  function setTaskContext(taskId, matchedRows = [], task = {}) {
    _currentTaskId = taskId || null;
    _currentMatched = Array.isArray(matchedRows) ? matchedRows : [];
    _currentIonMode = String(task.ion_mode || task.params?.ion_mode || '').toUpperCase() || null;
    _currentTaskFilename = String(task.filename || '');
  }

  // ── Refresh Task Select ───────────────────────────────────────────
  async function refreshSelect() {
    try {
      const data = await API.get('/tasks?status=completed');
      const sel = document.getElementById('reviewTaskSelect');
      if (!sel) return;
      const currentVal = sel.value;
      const latestVal = sel.value;
      sel.innerHTML = '<option value="">-- 选择已完成的任务 --</option>' +
        (data.tasks || []).map(t => {
          const stamp = t.upload_time
            ? String(t.upload_time).replace('T', ' ').replace(/\.\d+.*$/, '').slice(0, 16)
            : '时间未记录';
          const ms1 = Number(t.ms1_scans || 0);
          return `<option value="${t.id}">${Tasks.escapeHtml(t.filename)} · ${Tasks.escapeHtml(stamp)} · #${Tasks.escapeHtml(String(t.id).slice(0, 8))} (${ms1} MS1, ${t.matched_spectra}/${t.total_spectra} MS2)</option>`;
        }).join('');
      const targetVal = latestVal || currentVal || _currentTaskId;
      if (targetVal && data.tasks?.find(t => t.id === targetVal)) sel.value = targetVal;
    } catch (e) {
      console.error('Failed to refresh result select:', e);
    }
  }

  async function _loadAllResultSummaries(taskId) {
    const pageSize = 5000;
    let offset = 0;
    let rows = [];
    let governanceBlockedCandidates = [];
    let firstPayload = null;
    while (true) {
      const payload = await API.get(
        `/tasks/${encodeURIComponent(taskId)}/results?projection=summary&limit=${pageSize}&offset=${offset}`
      );
      if (!firstPayload) firstPayload = payload;
      const page = Array.isArray(payload.results) ? payload.results : [];
      rows = rows.concat(page);
      if (Array.isArray(payload.governance_blocked_candidates)) {
        governanceBlockedCandidates = governanceBlockedCandidates.concat(
          payload.governance_blocked_candidates
        );
      }
      if (!payload.has_more || !page.length) break;
      offset += page.length;
    }
    return {
      results: rows,
      view_counts: firstPayload?.view_counts || {},
      claim_projection: firstPayload?.claim_projection || {},
      governance_blocked_candidates: governanceBlockedCandidates.length
        ? governanceBlockedCandidates : null,
    };
  }

  function _projectBlockedFromSummaries(results) {
    const blocked = new Set(['insufficient_for_release', 'reported_but_governance_blocked']);
    return (results || []).flatMap((row, index) => {
      const match = row?.match || {};
      const passport = match.evidence_governance || {};
      const claimStatus = String(passport.claim_status || '');
      if (!blocked.has(claimStatus)) return [];
      const aliases = Array.isArray(match.compound_name_cn)
        ? match.compound_name_cn : (match.compound_name_cn ? [match.compound_name_cn] : []);
      const blockers = Array.isArray(passport.release_blockers) ? passport.release_blockers : [];
      const featureId = String(row.lcms_feature_id || row.precursor_feature_id || '');
      const searchText = [match.compound_name, ...aliases, match.formula, featureId,
        row.spectrum_id, claimStatus, ...blockers].filter(Boolean).join(' ').toLowerCase();
      return [{
        source_result_index: Number(row.result_index ?? index),
        source_feature_id: featureId,
        spectrum_id: String(row.spectrum_id || ''),
        compound_name: String(match.compound_name || ''),
        compound_name_cn: aliases,
        formula: String(match.formula || ''),
        retention_time: row.retention_time,
        retention_time_unit: String(row.retention_time_unit || ''),
        precursor_mz: row.precursor_mz,
        cosine: match.cosine,
        claim_status: claimStatus,
        claim_scope: String(passport.claim_scope || ''),
        policy_id: String(passport.policy_id || ''),
        release_eligible_before_replication: !!passport.release_eligible_before_replication,
        release_blockers: blockers,
        reason_text: blockers.join('; '),
        search_text: searchText,
      }];
    });
  }

  function _governanceBlockedCard(rows, taskId, totalCount = null, hasMore = false) {
    if (!Array.isArray(rows) || !rows.length) return '';
    const total = Number.isFinite(Number(totalCount)) && Number(totalCount) >= rows.length
      ? Number(totalCount) : rows.length;
    const body = rows.map((row) => {
      const aliases = Array.isArray(row.compound_name_cn) ? row.compound_name_cn : [];
      const aliasChip = aliases.length
        ? ` <span class="badge" style="background:#ecfeff;color:#0e7490;">${Tasks.escapeHtml(aliases.join('·'))}</span>` : '';
      const mz = Number(row.precursor_mz);
      const rt = Number(row.retention_time);
      const cosine = Number(row.cosine);
      const searchText = String(row.search_text || [row.compound_name, ...aliases, row.formula,
        row.source_feature_id, row.spectrum_id, row.claim_status, ...(row.release_blockers || [])]
        .filter(Boolean).join(' ')).toLowerCase();
      return `<tr class="governance-blocked-row" data-search="${Tasks.escapeHtml(searchText)}"
          data-mz="${Number.isFinite(mz) ? mz : ''}" style="cursor:pointer;"
          onclick="Results.showDetailByResultIndex('${taskId}', ${Number(row.source_result_index || 0)})">
        <td>${Number(row.source_result_index || 0)}</td><td>${Tasks.escapeHtml(row.spectrum_id || '—')}</td>
        <td>${Tasks.escapeHtml(row.source_feature_id || '—')}</td>
        <td><b>${Tasks.escapeHtml(row.compound_name || 'Unnamed candidate')}</b>${aliasChip}</td>
        <td>${Tasks.escapeHtml(row.formula || '—')}</td>
        <td class="mono">${Number.isFinite(rt) ? rt.toFixed(3) : '—'}</td>
        <td class="mono">${Number.isFinite(mz) ? mz.toFixed(5) : '—'}</td>
        <td class="mono">${Number.isFinite(cosine) ? cosine.toFixed(4) : '—'}</td>
        <td>${Tasks.escapeHtml(row.claim_status || '—')}</td>
        <td style="max-width:360px;white-space:normal;">${Tasks.escapeHtml(row.reason_text || (row.release_blockers || []).join('; ') || '—')}</td>
      </tr>`;
    }).join('');
    return `<div class="card">
      <details class="governance-blocked-panel">
        <summary style="cursor:pointer;font-weight:700;">治理未放行候选 ${total} 条——按原判定保留，仅供排查（扫描候选记录数）</summary>
        ${hasMore ? `<div style="padding:.45rem 0;color:#92400e;font-size:.72rem;">当前页面显示前 ${rows.length}/${total} 条；完整记录请使用结果页分页或 CSV/Excel 导出。</div>` : ''}
        <div style="padding:.65rem 0;display:flex;gap:.5rem;align-items:center;">
          <input id="governanceBlockedSearch" class="form-input" type="search" placeholder="筛选名称、中文名、分子式、m/z 或原因" oninput="Results.filterGovernanceBlockedRows(this.value)" style="max-width:620px;">
          <span id="governanceBlockedSearchCount" style="font-size:.75rem;color:var(--text-secondary);"></span>
        </div>
        <div class="table-wrap"><table><thead><tr><th>索引</th><th>谱图</th><th>特征</th><th>候选名称</th><th>分子式</th><th>RT</th><th>前体 m/z</th><th>余弦</th><th>治理状态</th><th>原始拦截原因</th></tr></thead>
        <tbody id="governanceBlockedBody">${body}</tbody></table></div>
      </details></div>`;
  }

  function _expectedFeatureEvidenceRank(feature) {
    const explicit = Number(feature?.evidence_rank);
    if (Number.isFinite(explicit)) return explicit;
    return {
      MS1_EXACT_MASS_CANDIDATE: 1,
      MS1_ISOTOPE_SUPPORTED: 2,
      MS1_MS2_DIAGNOSTIC_PROBABLE: 3,
      EXPERIMENTAL_LIBRARY_CONFIRMED: 4,
    }[String(feature?.evidence_level || '')] || 0;
  }

  function _expectedFeatureGroupKey(feature) {
    const formula = String(feature?.formula || '').trim();
    const members = Array.isArray(feature?.target_members)
      ? feature.target_members : [];
    const identities = Array.from(new Set(
      members.map(member => String(
        member?.identity_key
        || (member?.inchikey ? `ik::${member.inchikey}` : '')
        || ''
      ).trim().toLowerCase()).filter(Boolean)
    )).sort();
    const names = Array.from(new Set(
      (Array.isArray(feature?.candidate_names)
        ? feature.candidate_names
        : (feature?.candidate_name ? [feature.candidate_name] : []))
        .map(name => String(name || '').trim().toLowerCase())
        .filter(Boolean)
    )).sort();
    return JSON.stringify([formula, identities.length ? identities : names]);
  }

  function _expectedRepresentativeSort(a, b) {
    const blankRank = {
      sample_enriched: 3,
      not_detected: 2,
      blank_associated: 1,
      blank_dominant: 0,
    };
    const isotopeRank = feature =>
      feature?.formula_isotope_evidence?.status === 'supporting' ? 1 : 0;
    const arrayLength = value => Array.isArray(value) ? value.length : 0;
    return (
      _expectedFeatureEvidenceRank(b) - _expectedFeatureEvidenceRank(a)
      || Number(!!b?.cd_pattern_accepted) - Number(!!a?.cd_pattern_accepted)
      || arrayLength(b?.experimental_library_support)
        - arrayLength(a?.experimental_library_support)
      || arrayLength(b?.ms2_diagnostic_support)
        - arrayLength(a?.ms2_diagnostic_support)
      || isotopeRank(b) - isotopeRank(a)
      || Number(!!b?.openms_feature_supported)
        - Number(!!a?.openms_feature_supported)
      || Number(blankRank[String(b?.blank_evidence_category || '')] || 0)
        - Number(blankRank[String(a?.blank_evidence_category || '')] || 0)
      || Number(b?.peak_quality_score || 0)
        - Number(a?.peak_quality_score || 0)
      || Number(b?.baseline_corrected_area_intensity_seconds || 0)
        - Number(a?.baseline_corrected_area_intensity_seconds || 0)
      || Number(b?.apex_intensity || 0) - Number(a?.apex_intensity || 0)
      || Math.abs(Number(a?.mass_error_ppm || 0))
        - Math.abs(Number(b?.mass_error_ppm || 0))
      || Number(a?.apex_rt_minutes || 0) - Number(b?.apex_rt_minutes || 0)
    );
  }

  function _expectedGroupDisplayName(names, isomerAmbiguous) {
    const label = names.slice().sort((left, right) => (
      String(left).localeCompare(String(right))
    )).join(' / ')
      || '未命名靶标';
    return isomerAmbiguous ? `${label}（同质量候选组）` : label;
  }

  function _groupExpectedFeatures(features) {
    const grouped = new Map();
    (features || []).forEach(feature => {
      const key = _expectedFeatureGroupKey(feature);
      if (!grouped.has(key)) grouped.set(key, []);
      grouped.get(key).push(feature);
    });

    const groups = Array.from(grouped.values()).map(rawFeatures => {
      const rankedFeatures = rawFeatures.slice().sort(_expectedRepresentativeSort);
      const representative = rankedFeatures[0] || {};
      const names = Array.from(new Set(
        rankedFeatures.flatMap(feature =>
          Array.isArray(feature.candidate_names)
            ? feature.candidate_names
            : (feature.candidate_name ? [feature.candidate_name] : [])
        ).map(name => String(name || '')
          .replace(/\s+(?:\[\d*M[^\]]*\]\d*[+-]?|\d*M[+-][^\s]+)\s*$/i, '')
          .trim())
          .filter(name => name && !/^(?:NP[-_ ]?\d+|[A-Z]{1,8}[-_]\d+)$/i.test(name))
      )).sort((a, b) => a.localeCompare(b));
      const formulas = Array.from(new Set(
        rankedFeatures.map(feature => String(feature.formula || '').trim())
          .filter(Boolean)
      )).sort();
      const adducts = Array.from(new Set(
        rankedFeatures.map(feature => String(feature.adduct || '').trim())
          .filter(Boolean)
      )).sort();
      const ionHypotheses = new Set(
        rankedFeatures.map(feature => String(feature.expected_ion_id || ''))
          .filter(Boolean)
      );
      const datasets = Array.from(new Set(
        rankedFeatures.flatMap(feature =>
          Array.isArray(feature.dataset_sources)
            ? feature.dataset_sources : []
        ).map(item => String(item || '').trim()).filter(Boolean)
      )).sort();
      const identityKeys = new Set(
        rankedFeatures.flatMap(feature =>
          Array.isArray(feature.target_members)
            ? feature.target_members : []
        ).map(member => String(
          member?.identity_key || member?.inchikey || ''
        ).trim()).filter(Boolean)
      );
      const rtValues = rankedFeatures
        .map(feature => Number(feature.apex_rt_minutes))
        .filter(Number.isFinite);
      const supportedFeatureCount = rankedFeatures.filter(
        feature => !!feature.cd_pattern_accepted
      ).length;
      const isomerAmbiguous = identityKeys.size > 1
        || (identityKeys.size === 0 && names.length > 1);
      return {
        representative,
        features: rankedFeatures,
        candidateNames: names,
        displayName: _expectedGroupDisplayName(names, isomerAmbiguous),
        formulas,
        adducts,
        datasets,
        identityCount: identityKeys.size || names.length || 1,
        isomerAmbiguous,
        multipleNamesSameIdentity: identityKeys.size === 1
          && names.length > 1,
        rawFeatureCount: rankedFeatures.length,
        supportedFeatureCount,
        ionHypothesisCount: ionHypotheses.size,
        rtMinimum: rtValues.length ? Math.min(...rtValues) : null,
        rtMaximum: rtValues.length ? Math.max(...rtValues) : null,
      };
    });

    return groups.sort((a, b) => (
      _expectedRepresentativeSort(a.representative, b.representative)
      || b.supportedFeatureCount - a.supportedFeatureCount
      || b.rawFeatureCount - a.rawFeatureCount
      || a.displayName.localeCompare(b.displayName)
    ));
  }

  function _expectedCatalogAudit(payload) {
    const summary = payload?.summary || {};
    const audit = summary.target_catalog_quality_audit
      || payload?.catalog?.quality_audit || {};
    const total = Number(
      summary.target_catalog_rows ?? payload?.catalog?.target_rows ?? 0
    );
    const quarantined = Number(
      summary.quarantined_target_rows
      ?? audit.quarantined_target_rows
      ?? 0
    );
    const active = Number(
      summary.active_target_rows
      ?? audit.active_target_rows
      ?? Math.max(0, total - quarantined)
    );
    const reasonLabels = {
      'identity_status:rejected_identity_formula_mismatch': '已拒绝的名称/分子式身份冲突',
      invalid_nonempty_molecular_formula: '无效分子式',
      no_usable_formula_exact_mass_or_precursor: '缺少可用分子式、精确质量和前体',
    };
    const reasons = Object.entries(audit.quarantine_reason_counts || {})
      .filter(([, count]) => Number(count) > 0)
      .sort((left, right) => Number(right[1]) - Number(left[1]))
      .map(([reason, count]) => `${reasonLabels[reason] || reason} ${Number(count)}`);
    return {
      total,
      active,
      quarantined,
      reasonText: reasons.join('；'),
      note: quarantined
        ? `源库 ${total.toLocaleString()} 行中，仅 ${active.toLocaleString()} 行进入计算；${quarantined.toLocaleString()} 行因身份或化学元数据不合格被运行时隔离，原始数据库记录未删除${reasons.length ? `（${reasons.join('；')}）` : ''}。`
        : `源库 ${total.toLocaleString()} 行全部通过运行时化学质量门控。`,
    };
  }

  function _expectedGroupedCard(payload) {
    const groups = Array.isArray(payload.groups) ? payload.groups : [];
    const summary = payload.summary || {};
    const supportedGroups = groups.filter(
      group => Number(group.supported_feature_count || 0) > 0
    );
    const massOnlyGroups = groups.filter(
      group => Number(group.supported_feature_count || 0) === 0
    );
    const evidenceLabels = {
      MS1_EXACT_MASS_CANDIDATE: 'MS1 精确质量候选',
      MS1_ISOTOPE_SUPPORTED: 'MS1 + 同位素支持',
      MS1_MS2_DIAGNOSTIC_PROBABLE: 'MS1 + 诊断离子筛查支持',
      EXPERIMENTAL_LIBRARY_CONFIRMED: 'MS/MS 实验谱库候选支持',
    };
    const evidenceColors = {
      MS1_EXACT_MASS_CANDIDATE: ['#e2e8f0', '#475569'],
      MS1_ISOTOPE_SUPPORTED: ['#dbeafe', '#1d4ed8'],
      MS1_MS2_DIAGNOSTIC_PROBABLE: ['#fef3c7', '#92400e'],
      EXPERIMENTAL_LIBRARY_CONFIRMED: ['#dcfce7', '#166534'],
    };
    const num = (value, digits = 3) => {
      const parsed = Number(value);
      return Number.isFinite(parsed) ? parsed.toFixed(digits) : '—';
    };
    const scientific = value => {
      const parsed = Number(value);
      return Number.isFinite(parsed) ? parsed.toExponential(3) : '—';
    };
    const badge = group => {
      const level = String(group.evidence_level || 'MS1_EXACT_MASS_CANDIDATE');
      const colors = evidenceColors[level] || evidenceColors.MS1_EXACT_MASS_CANDIDATE;
      return `<span style="display:inline-block;padding:2px 6px;border-radius:9px;background:${colors[0]};color:${colors[1]};font-size:.66rem;font-weight:700;">${Tasks.escapeHtml(evidenceLabels[level] || level)}</span>`;
    };
    const rows = candidateGroups => candidateGroups.map((group, index) => {
      const ambiguity = group.isomer_ambiguous
        ? `<span style="display:inline-block;margin-top:2px;padding:1px 5px;border-radius:8px;background:#fef3c7;color:#92400e;font-size:.62rem;font-weight:700;">同质量数据库候选 ${Number(group.candidate_count || 2)}</span>`
        : '';
      const aliases = group.multiple_names_same_identity
        ? '<span style="display:inline-block;margin-top:2px;padding:1px 5px;border-radius:8px;background:#e0f2fe;color:#075985;font-size:.62rem;font-weight:700;">同一结构多名称</span>'
        : '';
      const rtRange = Number.isFinite(Number(group.rt_range_min))
        && Number.isFinite(Number(group.rt_range_max))
        ? `${num(group.rt_range_min)}–${num(group.rt_range_max)}` : '—';
      return `<tr>
        <td>${index + 1}</td>
        <td style="max-width:320px;white-space:normal;"><b>${Tasks.escapeHtml(group.candidate_display || '未命名靶标')}</b> ${ambiguity} ${aliases}<div style="font-size:.64rem;color:#64748b;">${Tasks.escapeHtml(group.target_database || '')}</div></td>
        <td class="mono">${Tasks.escapeHtml(group.formula || '—')}<br><span style="color:#64748b;">${Tasks.escapeHtml(group.related_adducts || group.adduct || '—')}</span></td>
        <td class="mono"><b>${Number(group.ion_hypothesis_count || 0)}</b> 离子 / <b>${Number(group.raw_peak_count || 0)}</b> 峰<br><span style="font-size:.64rem;color:${Number(group.supported_feature_count || 0) ? '#166534' : '#64748b'};">${Number(group.supported_feature_count || 0)} 峰有正交支持</span></td>
        <td class="mono">${num(group.target_mz, 6)}<br><span style="color:#64748b;">→ ${num(group.observed_mz, 6)} · ${num(group.mass_error_ppm, 2)} ppm</span></td>
        <td class="mono">${num(group.rt_min)}<br><span style="font-size:.64rem;color:#64748b;">全组 ${rtRange}</span></td>
        <td class="mono">${scientific(group.peak_height)}<br><span style="font-size:.64rem;color:#64748b;">A=${scientific(group.peak_area)}</span></td>
        <td>${badge(group)}</td>
      </tr>`;
    }).join('');
    const table = candidateGroups => `<div class="table-wrap" style="margin-top:.5rem;max-height:620px;overflow:auto;"><table><thead><tr><th>#</th><th>去重候选/歧义组</th><th>分子式/加合物</th><th>离子/峰</th><th>代表 m/z</th><th>代表 RT/范围</th><th>强度/面积</th><th>最高证据</th></tr></thead><tbody>${rows(candidateGroups)}</tbody></table></div>`;
    const totalGroups = Number(payload.group_total ?? groups.length);
    const supportedTotal = Number(payload.supported_group_total ?? supportedGroups.length);
    const massOnlyTotal = Number(payload.mass_only_group_total ?? massOnlyGroups.length);
    const rawFeatures = Number(summary.detected_feature_count || summary.features_found || 0);
    const supportedFeatures = Number(summary.supported_features || 0);
    const selectedDatasets = Array.isArray(summary.selected_dataset_sources)
      ? summary.selected_dataset_sources.filter(Boolean).join('、') : '';
    const catalogAudit = _expectedCatalogAudit(payload);
    return `<div class="card" style="border-left:4px solid #2563eb;background:#f8fbff;">
      <div class="card-header"><div><div class="card-title">🎯 靶向数据库结果（Expected Compounds，独立连续 MS1 层）</div><div class="card-subtitle">已应用数据集：${Tasks.escapeHtml(selectedDatasets || '未记录')}。服务端先按化合物/歧义组去重，不与普通 MS2 谱库结果混计。</div></div><span class="badge badge-high">已运行</span></div>
      <div class="grid-4" style="margin-top:.5rem;">
        <div class="stat-card"><div class="stat-value">${catalogAudit.total.toLocaleString()}</div><div class="stat-label">源库目标行</div></div>
        <div class="stat-card"><div class="stat-value" style="color:#166534;">${catalogAudit.active.toLocaleString()}</div><div class="stat-label">进入计算的有效行</div></div>
        <div class="stat-card"><div class="stat-value" style="color:${catalogAudit.quarantined ? '#b45309' : '#166534'};">${catalogAudit.quarantined.toLocaleString()}</div><div class="stat-label">运行时隔离行</div></div>
        <div class="stat-card"><div class="stat-value">${Number(summary.expected_ion_hypotheses || 0).toLocaleString()}</div><div class="stat-label">离子假设</div></div>
        <div class="stat-card" style="border:2px solid #2563eb;background:#eff6ff;"><div class="stat-value" style="color:#1d4ed8;">${totalGroups.toLocaleString()}</div><div class="stat-label"><b>候选组（已去重）</b></div></div>
        <div class="stat-card"><div class="stat-value">${supportedTotal.toLocaleString()}</div><div class="stat-label">有正交支持的组</div></div>
        <div class="stat-card"><div class="stat-value">${rawFeatures.toLocaleString()}</div><div class="stat-label">原始 XIC 峰（审计）</div></div>
        <div class="stat-card"><div class="stat-value">${supportedFeatures.toLocaleString()}</div><div class="stat-label">有支持原始峰</div></div>
      </div>
      <div style="margin-top:.55rem;padding:.5rem .7rem;background:${catalogAudit.quarantined ? '#fffbeb' : '#f0fdf4'};border-left:3px solid ${catalogAudit.quarantined ? '#d97706' : '#16a34a'};border-radius:6px;color:${catalogAudit.quarantined ? '#78350f' : '#166534'};font-size:.72rem;line-height:1.55;"><b>靶标库质量门控：</b>${Tasks.escapeHtml(catalogAudit.note)}</div>
      <div style="margin-top:.55rem;padding:.55rem .7rem;background:#eff6ff;border-radius:7px;color:#1e3a8a;font-size:.73rem;line-height:1.55;"><b>轻量结果视图：</b>这里保留候选、歧义、m/z、RT、峰强和证据等级；逐 XIC 峰的完整审计仍保存在任务工件和 Excel/CSV 导出中。精确质量不是异构体特异证据，<b>普通实验 MS2 谱库匹配增量 = 0。</b>${payload.has_more ? ' 当前仅显示首批候选组。' : ''}</div>
      ${supportedGroups.length ? `<details open style="margin-top:.65rem;"><summary style="cursor:pointer;font-size:.78rem;font-weight:700;">有正交证据支持的候选组（${supportedTotal}）</summary>${table(supportedGroups)}</details>` : '<div style="margin-top:.65rem;color:#64748b;font-size:.77rem;">本次没有获得正交支持的候选组。</div>'}
      ${massOnlyGroups.length ? `<details style="margin-top:.65rem;padding:.55rem .65rem;border:1px solid #cbd5e1;border-radius:7px;background:#fff;"><summary style="cursor:pointer;font-size:.78rem;font-weight:700;color:#475569;">仅精确质量的候选组（${massOnlyTotal}）— 点击展开</summary><div style="margin-top:.45rem;font-size:.71rem;color:#92400e;">不能据此区分同分异构体或确认具体化合物。</div>${table(massOnlyGroups)}</details>` : ''}
    </div>`;
  }

  function _expectedCompoundsCard(payload) {
    if (!payload || payload._request_error) {
      const reason = Tasks.escapeHtml(String(payload?._request_error || '结果接口不可用'));
      return `
        <div class="card" style="border-left:4px solid #94a3b8;background:#f8fafc;">
          <div class="card-title">🎯 靶向数据库结果（Expected Compounds，独立 MS1 层）</div>
          <div style="margin-top:.4rem;font-size:.77rem;color:#64748b;">
            本任务没有可读取的独立结果：${reason}。普通 MS2 谱库结果不受影响。
          </div>
        </div>`;
    }

    const status = String(payload.status || payload.summary?.status || 'not_run');
    const summary = payload.summary || {};
    const warning = String(summary.warning || payload.warning || '');
    const unavailable = status !== 'ok';
    const statusLabels = {
      disabled: '本次分析已关闭',
      not_run: '旧任务未运行',
      queued: '等待运行',
      failed: '运行失败',
      unavailable_requires_mzml_with_ms1: '不可用：缺少连续 MS1',
      unavailable_unknown_or_mixed_polarity: '不可用：极性未解析',
      unavailable_no_selected_target_dataset: '不可用：未绑定样本专属靶向库',
      unavailable_no_presence_only_targets: '不可用：没有文献靶标',
      unavailable_no_valid_expected_ions: '不可用：没有可搜索离子',
      unavailable_no_ms1_scans: '不可用：mzML 中没有 survey MS1',
      unavailable_no_usable_ms1_scans: '不可用：mzML 中没有可用 survey MS1',
    };
    if (unavailable) {
      const explicitInputWarning = status === 'unavailable_requires_mzml_with_ms1'
        ? 'MGF/HTML 只有离散 MS2 峰表，不能重建连续色谱 XIC；请使用保留 MS1+MS2 的 mzML。'
        : '';
      const catalogAudit = _expectedCatalogAudit(payload);
      const auditNote = catalogAudit.total
        ? `<div style="margin-top:.35rem;"><b>靶标库质量门控：</b>${Tasks.escapeHtml(catalogAudit.note)}</div>`
        : '';
      return `
        <div class="card" style="border-left:4px solid #d97706;background:#fffbeb;">
          <div class="card-header">
            <span class="card-title">🎯 靶向数据库结果（Expected Compounds，独立 MS1 层）</span>
            <span class="badge badge-medium">${Tasks.escapeHtml(statusLabels[status] || status)}</span>
          </div>
          <div style="font-size:.77rem;line-height:1.6;color:#78350f;">
            ${warning ? Tasks.escapeHtml(warning) : Tasks.escapeHtml(explicitInputWarning || '本任务未形成可用的 Expected Compounds 结果。')}
            ${warning && explicitInputWarning ? `<div>${Tasks.escapeHtml(explicitInputWarning)}</div>` : ''}
            ${auditNote}
            <div style="margin-top:.35rem;"><b>普通实验 MS2 谱库匹配保持原样，增量 = 0。</b></div>
          </div>
        </div>`;
    }

    if (Array.isArray(payload.groups)) {
      return _expectedGroupedCard(payload);
    }

    const features = Array.isArray(payload.features) ? payload.features : [];
    const evidenceCounts = summary.evidence_counts || {};
    const found = Number(
      payload.features_found ?? summary.detected_feature_count ?? features.length
    );
    const supported = Number(
      payload.supported_features ?? summary.supported_features
      ?? features.filter(feature => feature.cd_pattern_accepted).length
    );
    const evidenceLabels = {
      MS1_EXACT_MASS_CANDIDATE: 'MS1 精确质量候选',
      MS1_ISOTOPE_SUPPORTED: 'MS1 + 同位素支持',
      MS1_MS2_DIAGNOSTIC_PROBABLE: 'MS1 + 诊断离子筛查支持',
      EXPERIMENTAL_LIBRARY_CONFIRMED: 'MS/MS 实验谱库候选支持',
    };
    const evidenceColors = {
      MS1_EXACT_MASS_CANDIDATE: ['#e2e8f0', '#475569'],
      MS1_ISOTOPE_SUPPORTED: ['#dbeafe', '#1d4ed8'],
      MS1_MS2_DIAGNOSTIC_PROBABLE: ['#fef3c7', '#92400e'],
      EXPERIMENTAL_LIBRARY_CONFIRMED: ['#dcfce7', '#166534'],
    };
    const num = (value, digits = 3) => {
      const parsed = Number(value);
      return Number.isFinite(parsed) ? parsed.toFixed(digits) : '—';
    };
    const scientific = (value) => {
      const parsed = Number(value);
      return Number.isFinite(parsed) ? parsed.toExponential(3) : '—';
    };
    const evidenceBadge = feature => {
      const level = String(
        feature?.evidence_level || 'MS1_EXACT_MASS_CANDIDATE'
      );
      const colors = evidenceColors[level]
        || evidenceColors.MS1_EXACT_MASS_CANDIDATE;
      return `<span style="display:inline-block;padding:2px 6px;border-radius:9px;background:${colors[0]};color:${colors[1]};font-size:.66rem;font-weight:700;">${Tasks.escapeHtml(evidenceLabels[level] || level)}</span>`;
    };
    const groups = _groupExpectedFeatures(features);
    const supportedGroups = groups.filter(
      group => group.supportedFeatureCount > 0
    );
    const massOnlyGroups = groups.filter(
      group => group.supportedFeatureCount === 0
    );
    const groupRows = candidateGroups => candidateGroups.map(
      (group, index) => {
        const feature = group.representative;
        const ambiguity = group.isomerAmbiguous
          ? `<span style="display:inline-block;margin-top:2px;padding:1px 5px;border-radius:8px;background:#fef3c7;color:#92400e;font-size:.62rem;font-weight:700;">同质量数据库候选 ${Number(group.identityCount || 2)}</span>`
          : '';
        const aliases = group.multipleNamesSameIdentity
          ? `<span style="display:inline-block;margin-top:2px;padding:1px 5px;border-radius:8px;background:#e0f2fe;color:#075985;font-size:.62rem;font-weight:700;">同一结构多名称 ${group.candidateNames.length}</span>`
          : '';
        const rtRange = (
          Number.isFinite(group.rtMinimum)
          && Number.isFinite(group.rtMaximum)
        )
          ? `${num(group.rtMinimum, 3)}–${num(group.rtMaximum, 3)}`
          : '—';
        const rawFeatureAudit = group.features.map((rawFeature, rawIndex) => {
          const rawLevel = String(
            rawFeature.evidence_level || 'MS1_EXACT_MASS_CANDIDATE'
          );
          return `
            <div style="display:grid;grid-template-columns:34px minmax(72px,.65fr) minmax(160px,1.5fr) minmax(82px,.7fr) minmax(94px,.8fr) minmax(112px,1fr);gap:6px;padding:4px 6px;border-top:${rawIndex ? '1px solid #e2e8f0' : '0'};align-items:center;">
              <span>${rawIndex + 1}</span>
              <span class="mono">${Tasks.escapeHtml(String(rawFeature.adduct || '—'))}</span>
              <span class="mono">${num(rawFeature.target_mz, 6)} → ${num(rawFeature.observed_apex_mz, 6)} (${num(rawFeature.mass_error_ppm, 2)} ppm)</span>
              <span class="mono">RT ${num(rawFeature.apex_rt_minutes, 3)}</span>
              <span class="mono">I ${scientific(rawFeature.apex_intensity)}</span>
              <span>${Tasks.escapeHtml(evidenceLabels[rawLevel] || rawLevel)}</span>
            </div>`;
        }).join('');
        return `
          <tr>
            <td>${index + 1}</td>
            <td style="max-width:300px;white-space:normal;">
              <b>${Tasks.escapeHtml(group.displayName)}</b>${ambiguity}${aliases}
              <div style="font-size:.64rem;color:#64748b;">${Tasks.escapeHtml(group.datasets.join(' · '))}</div>
            </td>
            <td class="mono">${Tasks.escapeHtml(group.formulas.join(' / ') || '—')}<br><span style="color:#64748b;">${Tasks.escapeHtml(group.adducts.join(' · ') || '—')}</span></td>
            <td style="min-width:155px;"><span class="mono"><b>${group.ionHypothesisCount}</b> 离子 / <b>${group.rawFeatureCount}</b> 峰</span><br><span style="font-size:.64rem;color:${group.supportedFeatureCount ? '#166534' : '#64748b'};">${group.supportedFeatureCount} 峰有正交支持</span>
              <details style="margin-top:3px;">
                <summary style="cursor:pointer;color:#1d4ed8;font-size:.66rem;font-weight:700;">展开 ${group.rawFeatureCount} 个原始 XIC 峰</summary>
                <div style="margin-top:5px;min-width:690px;max-height:280px;overflow:auto;border:1px solid #dbeafe;border-radius:6px;background:#fff;font-size:.64rem;color:#475569;">
                  ${rawFeatureAudit}
                </div>
              </details>
            </td>
            <td class="mono">${num(feature.target_mz, 6)}<br><span style="color:#64748b;">→ ${num(feature.observed_apex_mz, 6)} · ${num(feature.mass_error_ppm, 2)} ppm</span></td>
            <td class="mono">${num(feature.apex_rt_minutes, 3)}<br><span style="font-size:.64rem;color:#64748b;">全组 ${rtRange}</span></td>
            <td class="mono">${scientific(feature.apex_intensity)}<br><span style="font-size:.64rem;color:#64748b;">A=${scientific(feature.xic_area_intensity_seconds)}</span></td>
            <td>${evidenceBadge(feature)}</td>
          </tr>`;
      }
    ).join('');
    const groupTable = candidateGroups => `
      <div class="table-wrap" style="margin-top:.5rem;max-height:620px;overflow:auto;">
        <table>
          <thead><tr><th>#</th><th>去重候选/歧义组</th><th>分子式/检测加合物</th><th>离子/原始峰审计</th><th>代表 m/z</th><th>代表 RT/全组范围</th><th>代表强度/面积</th><th>最高证据</th></tr></thead>
          <tbody>${groupRows(candidateGroups)}</tbody>
        </table>
      </div>`;
    const evidenceSummary = Object.entries(evidenceLabels).map(([key, label]) =>
      `${label} ${Number(evidenceCounts[key] || 0)}`
    ).join(' · ');
    const pagingNote = payload.has_more
      ? `当前载入 ${features.length}/${Number(payload.filtered_total || found)} 个原始峰；本页去重组只基于已载入部分。`
      : `已将 ${features.length} 个原始 XIC 峰折叠为 ${groups.length} 个化合物/歧义候选组。`;
    const selectedDatasets = Array.isArray(summary.selected_dataset_sources)
      ? summary.selected_dataset_sources.filter(Boolean).join('、') : '';
    const catalogAudit = _expectedCatalogAudit(payload);

    return `
      <div class="card" style="border-left:4px solid #2563eb;background:#f8fbff;">
        <div class="card-header">
          <div>
            <div class="card-title">🎯 靶向数据库结果（Expected Compounds，独立连续 MS1 层）</div>
            <div class="card-subtitle">已应用数据集：${Tasks.escapeHtml(selectedDatasets || '未记录')}。按文献靶标提取 XIC；不混入普通 MS2 谱库命中数、High 数或去重候选数。</div>
          </div>
          <span class="badge badge-high">已运行</span>
        </div>
        <div class="grid-4" style="margin-top:.5rem;">
          <div class="stat-card"><div class="stat-value">${catalogAudit.total.toLocaleString()}</div><div class="stat-label">源库目标行</div></div>
          <div class="stat-card"><div class="stat-value" style="color:#166534;">${catalogAudit.active.toLocaleString()}</div><div class="stat-label">进入计算的有效行</div></div>
          <div class="stat-card"><div class="stat-value" style="color:${catalogAudit.quarantined ? '#b45309' : '#166534'};">${catalogAudit.quarantined.toLocaleString()}</div><div class="stat-label">运行时隔离行</div></div>
          <div class="stat-card"><div class="stat-value">${Number(summary.expected_ion_hypotheses || 0).toLocaleString()}</div><div class="stat-label">离子假设</div></div>
          <div class="stat-card" style="border:2px solid #2563eb;background:#eff6ff;"><div class="stat-value" style="color:#1d4ed8;">${groups.length.toLocaleString()}</div><div class="stat-label"><b>候选组（已去重）</b></div></div>
          <div class="stat-card"><div class="stat-value">${supportedGroups.length.toLocaleString()}</div><div class="stat-label">有正交支持的候选组</div></div>
          <div class="stat-card"><div class="stat-value">${found.toLocaleString()}</div><div class="stat-label">原始 XIC 峰（审计）</div></div>
          <div class="stat-card"><div class="stat-value">${supported.toLocaleString()}</div><div class="stat-label">有正交支持的原始峰</div></div>
        </div>
        <div style="margin-top:.55rem;padding:.5rem .7rem;background:${catalogAudit.quarantined ? '#fffbeb' : '#f0fdf4'};border-left:3px solid ${catalogAudit.quarantined ? '#d97706' : '#16a34a'};border-radius:6px;color:${catalogAudit.quarantined ? '#78350f' : '#166534'};font-size:.72rem;line-height:1.55;"><b>靶标库质量门控：</b>${Tasks.escapeHtml(catalogAudit.note)}</div>
        <div style="margin-top:.55rem;padding:.55rem .7rem;background:#eff6ff;border-radius:7px;color:#1e3a8a;font-size:.73rem;line-height:1.55;">
          <b>已去重显示：</b>同一化合物或同质量歧义组的不同加合物、不同 RT 色谱峰不再重复铺成多行；主表每个候选组只显示一次，全部原始峰仍可在行末点击展开审计。
          原始峰证据分布：${Tasks.escapeHtml(evidenceSummary)}。精确质量单独命中不是异构体特异鉴定。
          <b>普通实验 MS2 谱库匹配增量 = 0。</b> ${Tasks.escapeHtml(pagingNote)}
        </div>
        ${supportedGroups.length ? `
          <details open style="margin-top:.65rem;">
            <summary style="cursor:pointer;font-size:.78rem;font-weight:700;">有正交证据支持的去重候选组（${supportedGroups.length}）</summary>
            ${groupTable(supportedGroups)}
          </details>` : '<div style="margin-top:.65rem;color:#64748b;font-size:.77rem;">本次没有获得同位素、诊断离子或严格实验谱库正交支持的候选组。</div>'}
        ${massOnlyGroups.length ? `
          <details style="margin-top:.65rem;padding:.55rem .65rem;border:1px solid #cbd5e1;border-radius:7px;background:#fff;">
            <summary style="cursor:pointer;font-size:.78rem;font-weight:700;color:#475569;">仅精确质量的去重候选组（${massOnlyGroups.length}）— 点击展开</summary>
            <div style="margin-top:.45rem;font-size:.71rem;color:#92400e;">这部分只有目标 m/z 与色谱峰形，默认收起；不能据此区分同分异构体或确认具体化合物。</div>
            ${groupTable(massOnlyGroups)}
          </details>` : ''}
      </div>`;
  }

  // ── Load Results for Selected Task ────────────────────────────────
  async function loadForTask(taskIdOverride = null) {
    const sel = document.getElementById('reviewTaskSelect');
    const taskId = taskIdOverride || sel?.value;
    const area = document.getElementById('resultsArea');

    if (taskIdOverride && sel) sel.value = taskIdOverride;

    if (!taskId) {
      _currentTaskFilename = '';
      if (typeof ModuleEffects !== 'undefined') ModuleEffects.render(null);
      if (typeof EviAcquire !== 'undefined') EviAcquire.render(null);
      if (area) area.innerHTML = `
        <div class="empty-state"><div class="empty-icon">🔬</div>
        <div class="empty-title">选择任务以查看结果</div>
        <div class="empty-desc">请先选择一个已完成的分析任务</div></div>`;
      return;
    }
    if (_currentTaskId !== taskId) {
      _unmatchedDisplayLimit = 50;
      _currentMatched = [];
      _currentUnmatched = [];
      _currentUnmatchedTotal = 0;
    }
    _currentTaskId = taskId;
    if (typeof ModuleEffects !== 'undefined') {
      ModuleEffects.render(taskId, 'moduleEffectPanel');
    }
    if (typeof EviAcquire !== 'undefined') {
      EviAcquire.render(taskId, 'eviAcquirePanel');
    }

    try {
      const [task, resultPayload, expectedCompounds] = await Promise.all([
        API.get(`/tasks/${encodeURIComponent(taskId)}/status`),
        _loadAllResultSummaries(taskId),
        API.get(`/tasks/${encodeURIComponent(taskId)}/expected-compounds/groups?limit=1000`)
          .catch(error => ({ _request_error: error.message })),
      ]);
      _currentIonMode = task.ion_mode || task.params?.ion_mode || null;
      _currentTaskFilename = String(task.filename || '');
      const results = resultPayload.results || [];
      const governanceBlockedCandidates = Array.isArray(
        resultPayload.governance_blocked_candidates
      ) ? resultPayload.governance_blocked_candidates : _projectBlockedFromSummaries(results);
      // 仅保留真实库匹配: 排除历史 ICEBERG 结构推理兜底结果 (已下线; 旧任务 results 里残留的
      // putative 结果 cosine/final=0.70 为占位假值非真实相似度, 且每条都是假 InChIKey 无法去重)。
      const viewResults = results.filter(r => {
        if (_blankView === 'clean') return !r.blank_clean_view_excluded;
        if (_blankView === 'flagged') return !!r.blank_flagged;
        return true;
      });
      const matched = viewResults
        .filter(r => r.match && (r.match.dataset_source !== 'iceberg_fallback'))
        .sort(_compareMatched);
      _currentMatched = matched;   // 缓存排序后的列表(未折叠), 供 showDetail 按 index 取

      // 按化合物折叠去重（完整 InChIKey；缺失时才按结构/name+formula 降级）。同物多谱合并为一行，
      // 代表=首次出现(matched 已按 raw rank_score 优先排序). 解决"同一化合物反复列出上百次"
      // (样本内同一化合物的多张 DDA 谱图各自成行). 异构体(不同连接层)不合并。
      const _gmap = new Map();
      matched.forEach((r, i) => {
        const m = r.match;
        const k = _compoundIdentityKey(m);
        if (!_gmap.has(k)) _gmap.set(k, { repIdx: i, rep: r, nSpectra: 1 });
        else _gmap.get(k).nSpectra++;
      });
      const groups = [..._gmap.values()];
      const contaminantGroups = groups.filter(g => !!g.rep.match.potential_contaminant);
      const blankExcluded = results.filter(r => !!r.blank_excluded);
      const blankFlagged = results.filter(r => !!r.blank_flagged);
      const blankCleanExcluded = results.filter(r => !!r.blank_clean_view_excluded);
      const unmatchedAll = viewResults.filter(r => !r.match && !r.blank_excluded);
      const unmatched = unmatchedAll.slice(0, _unmatchedDisplayLimit);
      _currentUnmatched = unmatched;
      _currentUnmatchedTotal = unmatchedAll.length;
      const total = results.length;
      const viewTotal = viewResults.length;
      const rate = viewTotal ? (matched.length / viewTotal * 100).toFixed(1) : 0;
      const blankViewLabel = _blankView === 'clean'
        ? '空白清洁视图'
        : (_blankView === 'flagged' ? '仅空白相关' : '完整鉴定');
      // High 必须同时展示两个不可互换的口径：
      //   1) 至少一张独立 MS2 被校准为 High 的去重身份（High-supported identity）；
      //   2) 该身份按原始 rank_score 选出的最高分代表谱本身仍为 High。
      // 后者可因同分异构体近邻而是 Medium，即使同一身份另有一张 High 支持谱。
      // 这里只统计，不重排候选、不提升任何置信度。
      const frontendHighSupported = new Set(
        matched
          .filter(r => (r.match.confidence || 'Unknown') === 'High')
          .map(r => _compoundIdentityKey(r.match))
      ).size;
      const frontendRepresentativeHigh = groups.filter(
        g => (g.rep.match.confidence || 'Unknown') === 'High'
      ).length;
      const confidenceSummary = task.compound_confidence_summary || {};
      const summaryCount = (key, fallback) => {
        const value = Number(confidenceSummary[key]);
        return Number.isFinite(value) && value >= 0 ? value : fallback;
      };
      const nHighSupportedIdentities = summaryCount(
        'full_identity_high_supported_unique_compounds',
        summaryCount('high_supported_unique_compounds', frontendHighSupported)
      );
      const nRepresentativeHigh = summaryCount(
        'representative_full_identity_high_unique_compounds',
        summaryCount('representative_high_unique_compounds', frontendRepresentativeHigh)
      );
      const nHighSpectra = matched.filter(r => (r.match.confidence || 'Unknown') === 'High').length;
      const frontendIsomerGroupHigh = groups.filter(g => {
        const m = g.rep.match;
        return String(m.identification_scope || '') === 'isomer_group'
          && (m.isomer_group_confidence || 'Unknown') === 'High';
      }).length;
      // Backward-compatible total class-High audit count.  It is hierarchical
      // and therefore not displayed as the disjoint "class-only" card below.
      const nClassHigh = groups.filter(
        g => (g.rep.match.confidence_class || 'Unknown') === 'High'
      ).length;
      const frontendClassOnlyHigh = groups.filter(g => {
        const m = g.rep.match;
        return (m.confidence_class || 'Unknown') === 'High'
          && (m.confidence || 'Unknown') !== 'High'
          && !(String(m.identification_scope || '') === 'isomer_group'
            && (m.isomer_group_confidence || 'Unknown') === 'High')
          && !(
            String(m.two_dimensional_connectivity_confidence || '') === 'High'
            && String(m.identification_scope || '') !== 'isomer_group'
          );
      }).length;
      const nIsomerGroupHigh = summaryCount(
        'representative_isomer_group_high_unique_compounds',
        frontendIsomerGroupHigh
      );
      const compoundConclusionSummary = task.compound_conclusion_summary || {};
      const nHighExactStructures = Number(
        compoundConclusionSummary.high_confidence_2d_connectivity_count
        ?? compoundConclusionSummary.high_confidence_detected_count
        ?? nHighSupportedIdentities
      );
      const nHighScopeIsomerGroups = Number(
        compoundConclusionSummary.high_confidence_isomer_group_count
        ?? nIsomerGroupHigh
      );
      const nHighScopeGroups = Number(
        compoundConclusionSummary.high_confidence_scope_group_count
        ?? confidenceSummary.high_reportable_scope_unique_compounds
        ?? (nHighExactStructures + nHighScopeIsomerGroups)
      );
      const nClassOnlyHigh = summaryCount(
        'representative_class_only_high_unique_compounds',
        frontendClassOnlyHigh
      );
      const nStrong = groups.filter(g => (g.rep.match.evidence_tier || '') === 'Strong Candidate').length;
      const nInsufficient = groups.filter(g => (g.rep.match.confidence || 'Unknown') === 'Unknown').length;
      // Backend PF IDs use deterministic m/z + RT + polarity acquisition-event
      // clustering.  The fallback keeps historical cached tasks readable.
      const _fallbackPrecursorKey = (r) => {
        const mz = Number(r.precursor_mz || 0);
        const rt = Number(r.retention_time);
        const mode = r.ion_mode || '';
        const mzBin = Math.round(mz * 1000);       // about 0.001 Da
        const rtBin = Number.isFinite(rt) ? Math.round(rt / 0.05) : String(r.spectrum_id || '');
        return `${mode}|${mzBin}|${rtBin}`;
      };
      const _precursorKey = (r) => r.precursor_feature_id || _fallbackPrecursorKey(r);
      const precursorClusters = new Set(results.map(_precursorKey)).size;
      const matchedPrecursorClusters = new Set(matched.map(_precursorKey)).size;
      const highPrecursorClusters = new Set(
        matched.filter(r => (r.match.confidence || 'Unknown') === 'High').map(_precursorKey)
      ).size;
      const limitedInput = results.filter(r => r.limited_fragment_list);
      const ionNetworks = new Set(results.map(r => r.ion_identity_network_id).filter(Boolean));
      const ionNetworkedFeatures = new Set(results.filter(r => r.ion_identity_network_id).map(_precursorKey));
      const ionLinked = results.filter(r => r.ion_identity_linked_annotation).length;
      const chimericRows = matched.filter(r => ['mixed', 'severe'].includes(String((r.match.chimeric_deconvolution || {}).isolation_status || '')));
      const deconvolvedRows = matched.filter(r => (r.match.chimeric_deconvolution || {}).status === 'accepted_score_applied');
      const moduleAvailability = task.module_availability || {};
      const iimnAvailability = moduleAvailability.iimn || {};
      const deconvolutionAvailability = moduleAvailability.chimeric_deconvolution || {};
      const identificationFunnel = task.identification_funnel || {};
      const deconvolutionImpact = task.deconvolution_impact || {};
      const evidenceGovernance = task.evidence_governance_summary || {};
      const claimContract = task.claim_contract_summary || {};
      const compiledFinalClaims = groups.filter(
        g => !!(g.rep.match.final_claim_decision || {}).policy_id
      );
      const releasedFinalClaims = compiledFinalClaims.filter(
        g => !!(g.rep.match.final_claim_decision || {}).released
      );
      const governanceIdentityCounts = evidenceGovernance.claim_status_unique_identities || {};
      const governancePreeligible = Number(
        evidenceGovernance.release_eligible_before_replication_unique_identities
        ?? groups.filter(g => !!(g.rep.match.evidence_governance || {}).release_eligible_before_replication).length
      );
      const governanceNeedsReplicate = Number(
        governanceIdentityCounts.provisional_needs_independent_injection || 0
      );
      const governanceBlocked = Number(
        governanceIdentityCounts.reported_but_governance_blocked || 0
      ) + Number(governanceIdentityCounts.insufficient_for_release || 0)
        + Number(governanceIdentityCounts.blank_conflicted_candidate || 0);
      const governanceIsomerOnly = Number(
        governanceIdentityCounts.isomer_group_only || 0
      );
      const untargetedWorkflow = task.untargeted_workflow || {};
      const featureAttachment = untargetedWorkflow.ms2_attachment || {};
      const openmsDiagnostics = Array.isArray(untargetedWorkflow.algorithm_diagnostics)
        ? untargetedWorkflow.algorithm_diagnostics
        : [];
      const lcmsFeatureIds = new Set(
        results.map(r => r.lcms_feature_id).filter(Boolean)
      );
      const lcmsFeatureAvailable = lcmsFeatureIds.size > 0;
      const matchedLcmsFeatures = new Set(
        matched.map(r => r.lcms_feature_id).filter(Boolean)
      ).size;
      const highLcmsFeatures = new Set(
        matched
          .filter(r => (r.match.confidence || 'Unknown') === 'High')
          .map(r => r.lcms_feature_id)
          .filter(Boolean)
      ).size;
      const lcmsLinkedMs2 = Number(
        featureAttachment.linked_ms2_spectra
        || results.filter(r => r.lcms_feature_id).length
      );
      const peakProcessing = task.peak_processing_summary || {};
      const representation = task.spectrum_representation_processing || {};
      const consensusSummary = task.consensus_spectrum_summary
        || untargetedWorkflow.consensus_spectra || {};
      const ms1ScanCount = Number(moduleAvailability.ms1_scan_count || 0);
      const iimnUnavailable = !!iimnAvailability.requested
        && String(iimnAvailability.status || '').startsWith('unavailable');
      const deconvolutionUnavailable = !!deconvolutionAvailability.requested
        && String(deconvolutionAvailability.status || '').startsWith('unavailable');
      const analysisWarnings = Array.isArray(task.analysis_warnings)
        ? task.analysis_warnings
        : (Array.isArray(moduleAvailability.warnings) ? moduleAvailability.warnings : []);
      const moduleWarningHtml = analysisWarnings.map(warning => {
        const message = Tasks.escapeHtml(String(warning?.message || warning || ''));
        const action = warning?.action
          ? `<div style="margin-top:.2rem;">建议：${Tasks.escapeHtml(String(warning.action))}</div>`
          : '';
        return `<div style="margin-top:.35rem;">⚠ ${message}${action}</div>`;
      }).join('');
      const moduleWarningCount = analysisWarnings.length;
      const workflowDescription = moduleAvailability.continuous_ms1_available === false
        ? `本文件没有连续 MS1，IIMN 与依赖 MS1 隔离纯度的嵌合谱去卷积均未实际执行；当前 ${groups.length} 个去重候选只来自原始 MS2 谱库检索。`
        : (lcmsFeatureAvailable
          ? `OpenMS 检出 ${Number(untargetedWorkflow.feature_count || 0)} 个连续 MS1 色谱特征；${lcmsLinkedMs2}/${total} 张 MS2 成功归属到 ${lcmsFeatureIds.size} 个有 MS2 的特征，并生成 ${Number(consensusSummary.consensus_spectra || 0)} 张特征级共识谱。`
          : `连续 MS1 色谱特征未形成可用 MS2 链接；当前网络必须把谱图标为前体事件回退节点，不能把这些节点称为定量 LC–MS feature。`);
      const peakStageLabels = {
        valid_numeric_positive: '有效正强度峰',
        relative_intensity: '相对强度过滤',
        precursor_exclusion: '前体邻域排除',
        upper_mass_bound: '质量上限',
        top_n: 'Top-N',
      };
      const peakStageHtml = (peakProcessing.stages || []).map(stage =>
        `<div style="padding:.45rem .55rem;border:1px solid #e2e8f0;border-radius:7px;background:#fff;"><b>${Tasks.escapeHtml(peakStageLabels[stage.stage] || stage.stage)}</b><br><span>${Number(stage.input_peaks || 0).toLocaleString()} → ${Number(stage.retained_peaks || 0).toLocaleString()}（${stage.retained_from_previous_pct == null ? '—' : Number(stage.retained_from_previous_pct).toFixed(1) + '%'}）</span></div>`
      ).join('');
      const representationAudit = representation.representation_audit_before
        || representation.representation_audit || {};
      const representationCounts = representationAudit.counts || {};
      const representationText = ['ms1', 'ms2'].map(level => {
        const counts = representationCounts[level] || {};
        return `${level.toUpperCase()}：centroid ${Number(counts.centroid || 0)} / profile ${Number(counts.profile || 0)} / 未知 ${Number(counts.unknown || 0)}`;
      }).join('；');
      const openmsDiagnosticHtml = openmsDiagnostics.map(item =>
        `<div style="margin-top:.35rem;color:${item.severity === 'warning' ? '#92400e' : '#526579'};">${item.severity === 'warning' ? '⚠' : 'ℹ'} ${Tasks.escapeHtml(String(item.interpretation || item.code || 'OpenMS 诊断'))}</div>`
      ).join('');
      const blockerLabels = {
        no_validated_structure: '无验证结构',
        cosine_below_0_70: 'MS/MS 相似度 < 0.70',
        cosine_below_0_80: 'MS/MS 相似度 0.70–0.80',
        cosine_below_high_threshold: 'MS/MS 相似度 0.80–0.90',
        final_score_below_high_threshold: '综合匹配分 < 0.50',
        precursor_formula_adduct_conflict: '前体/分子式/加合物冲突',
        mass_error_above_50_ppm: '|质量误差| > 50 ppm',
        reverse_coverage_below_0_60: 'Reverse coverage < 0.60',
        competitive_same_formula_isomer: '同分子式竞争异构体',
        class_consensus_below_0_50: '类别共识 < 0.50',
        ion_identity_conflict: 'IIMN 身份/加合物冲突',
        library_source_confidence_cap: '谱库来源等级封顶',
        other_calibration_guard: '其他校准门控',
      };
      const primaryBlockers = Object.entries(
        identificationFunnel.non_high_primary_blockers || {}
      ).sort((a, b) => Number(b[1]) - Number(a[1]));
      const primaryBlockerHtml = primaryBlockers.length
        ? primaryBlockers.slice(0, 6).map(([key, value]) =>
            `<span style="display:inline-block;margin:.18rem .28rem .18rem 0;padding:.16rem .42rem;border-radius:10px;background:#fff7ed;color:#9a3412;">${Tasks.escapeHtml(blockerLabels[key] || key)} <b>${Number(value)}</b></span>`
          ).join('')
        : '<span>无</span>';
      const impactBaseline = deconvolutionImpact.baseline || {};
      const impactCurrent = deconvolutionImpact.current || {};
      const impactDelta = deconvolutionImpact.delta || {};
      const signed = (value) => `${Number(value) >= 0 ? '+' : ''}${Number(value || 0)}`;
      const deconvolutionInvariantOk = (
        deconvolutionImpact.top1_changed !== true
        && deconvolutionImpact.candidate_margin_changed !== true
        && deconvolutionImpact.raw_rank_score_changed !== true
      );
      const deconvolutionInvariantText = deconvolutionInvariantOk
        ? '逐谱反事实审计通过：Top-1、原始 rank 与候选 margin 均未改变。'
        : '<b style="color:#b91c1c;">不变量审计失败：检测到 Top-1、原始 rank 或候选 margin 被改变；该结果不得作为去卷积分数增益使用。</b>';
      const deconvolutionImpactHtml = deconvolutionImpact.method ? `
        <div style="margin-top:.65rem;padding:.65rem .8rem;border-left:4px solid #0f766e;background:#f0fdfa;color:#134e4a;font-size:.76rem;line-height:1.6;">
          <b>严格去卷积净增（冻结原始谱反事实）</b>：
          通过全部门控 ${Number(deconvolutionImpact.accepted_score_increase_spectra || 0)} 张；
          High-supported 身份 ${Number(impactBaseline.high_supported_unique_compounds || 0)}→${Number(impactCurrent.high_supported_unique_compounds || 0)}
         （${signed(impactDelta.high_supported_unique_compounds)}），
          最高分代表 High ${Number(impactBaseline.representative_high_unique_compounds || 0)}→${Number(impactCurrent.representative_high_unique_compounds || 0)}
         （${signed(impactDelta.representative_high_unique_compounds)}），
          High 前体事件 ${Number(impactBaseline.high_precursor_features || 0)}→${Number(impactCurrent.high_precursor_features || 0)}
         （${signed(impactDelta.high_precursor_features)}），
          High 谱 ${Number(impactBaseline.high_spectra || 0)}→${Number(impactCurrent.high_spectra || 0)}
         （${signed(impactDelta.high_spectra)}）。${deconvolutionInvariantText}
          </div>` : '';
      const reportContext = claimContract.context || {};
      const independentInjectionCount = new Set(
        reportContext.independent_injection_groups || []
      ).size;
      const minimumIndependentInjections = Number(
        reportContext.min_independent_injections || 2
      );
      const globalMissingConditions = [];
      if (Number(claimContract.calibrated_count || 0) < Number(claimContract.record_count || 0)) {
        globalMissingConditions.push('项目级 target–decoy q 值校准尚未建立');
      }
      if (independentInjectionCount < minimumIndependentInjections) {
        globalMissingConditions.push(`独立进样 ${independentInjectionCount}/${minimumIndependentInjections}`);
      }
      const formalReportConditionHtml = (claimContract.policy_id || evidenceGovernance.policy_id) ? `
        <section class="formal-report-overview">
          <div class="formal-report-overview-head">
            <b>正式报告条件（不改变候选鉴定）</b>
            <span class="formal-report-badge ${releasedFinalClaims.length ? 'eligible' : 'candidate'}">${releasedFinalClaims.length ? `${releasedFinalClaims.length} 个代表候选满足报告条件` : '当前正式报告条件未满足'}</span>
          </div>
          <div style="margin-top:.35rem;">已评估 ${compiledFinalClaims.length || groups.length} 个去重代表候选；${releasedFinalClaims.length ? `其中 ${releasedFinalClaims.length} 个满足当前正式报告条件。` : '尚无代表候选满足全部正式报告条件，但候选评分、证据等级与鉴定范围仍可正常查看。'}</div>
          <div class="formal-report-overview-reasons">${(globalMissingConditions.length ? globalMissingConditions : ['本任务未记录全项目共同缺失条件']).map(item => `<span>${Tasks.escapeHtml(item)}</span>`).join('')}</div>
          <div style="margin-top:.35rem;">以上是本任务全部候选共同适用的项目级条件，仅在此说明一次；逐候选详情只列候选自身最重要的缺失条件。综合匹配分不是鉴定概率。</div>
          <details><summary>查看适用性与独立性检查摘要</summary><div style="margin-top:.3rem;">${Number(evidenceGovernance.unique_identities_audited || groups.length)} 个去重身份完成检查：${governancePreeligible} 个通过来源/质量前置条件，${governanceNeedsReplicate} 个仍需独立进样，${governanceIsomerOnly} 个当前只到异构体组，${governanceBlocked} 个还需补充候选级证据。同一进样内的多碰撞能扫描不能充当独立重复，且这些检查不改变当前候选排序。</div></details>
          <button class="btn btn-secondary" style="margin-top:.48rem;font-size:.7rem;padding:.28rem .65rem;" onclick="Results.openGovernanceCohort()">检查跨进样一致性</button>
        </section>` : '';

      let html = `
      <div class="card">
        <div class="card-header">
          <span class="card-title">${Tasks.escapeHtml(task.filename)} — 完整技术摘要</span>
          <span style="display:flex;gap:.55rem;align-items:center;flex-wrap:wrap;">
            <label style="font-size:.72rem;color:var(--text-secondary);">结果视图
              <select class="form-select" onchange="Results.setBlankView(this.value)"
                      style="width:auto;padding:2px 6px;margin-left:4px;font-size:.72rem;">
                <option value="all"${_blankView === 'all' ? ' selected' : ''}>完整鉴定（推荐）</option>
                <option value="clean"${_blankView === 'clean' ? ' selected' : ''}>空白清洁视图</option>
                <option value="flagged"${_blankView === 'flagged' ? ' selected' : ''}>仅看空白相关</option>
              </select>
            </label>
            ${Tasks.statusBadge(task.status)}
          </span>
        </div>
        <div class="grid-4">
          <div class="stat-card" title="连续 survey scans；用于 XIC、同位素、隔离纯度和 IIMN，不直接等同于化合物鉴定"><div class="stat-value">${ms1ScanCount}</div><div class="stat-label">MS1 survey 扫描</div></div>
          <div class="stat-card"><div class="stat-value">${total}</div><div class="stat-label">MS2 碎裂扫描</div></div>
          <div class="stat-card" title="${lcmsFeatureAvailable ? '由连续 MS1 色谱质量轨迹检测并与 MS2 链接' : '按前体 m/z、RT和极性聚合同一 DDA/阶梯碰撞能事件；不是基于连续MS1的定量 feature'}"><div class="stat-value">${lcmsFeatureAvailable ? lcmsFeatureIds.size : precursorClusters}</div><div class="stat-label">${lcmsFeatureAvailable ? '有 MS2 的 LC–MS 特征' : 'MS2 前体事件（回退）'}</div></div>
          <div class="stat-card"><div class="stat-value">${matched.length}</div><div class="stat-label">有候选谱图</div></div>
          <div class="stat-card" title="至少一张归属 MS2 获得谱库候选的节点"><div class="stat-value">${lcmsFeatureAvailable ? matchedLcmsFeatures : matchedPrecursorClusters}</div><div class="stat-label">${lcmsFeatureAvailable ? '有候选 LC–MS 特征' : '有候选回退事件'}</div></div>
          <div class="stat-card"><div class="stat-value">${groups.length}</div><div class="stat-label">谱库匹配候选（去重）</div></div>
          <div class="stat-card" title="按当前证据可报告范围统计：单一二维连接 High 与未解异构体组 High 相加，不把异构体组误写成唯一结构"><div class="stat-value">${nHighScopeGroups}<span style="font-size:.45em;color:var(--text-secondary);font-weight:400">（${nHighExactStructures} 单一二维连接 + ${nHighScopeIsomerGroups} 异构体组）</span></div><div class="stat-label">High 级鉴定组</div></div>
          <div class="stat-card" title="完整结构身份（含立体／定义层）达到 High 的旧严格口径；二维连接 High 另由前一张卡按可报告范围统计"><div class="stat-value">${nHighSupportedIdentities}<span style="font-size:.45em;color:var(--text-secondary);font-weight:400">（${lcmsFeatureAvailable ? highLcmsFeatures : highPrecursorClusters}${lcmsFeatureAvailable ? ' LC–MS特征' : '回退事件'}/${nHighSpectra}谱）</span></div><div class="stat-label">完整结构身份 High 支持</div></div>
          <div class="stat-card" title="每个身份按原始 rank_score 选择最高分代表；该数仅表示完整结构身份 High，不替代二维连接层结论"><div class="stat-value">${nRepresentativeHigh}</div><div class="stat-label">完整身份代表 High</div></div>
          <div class="stat-card"><div class="stat-value">${nStrong}</div><div class="stat-label">Strong Candidate</div></div>
          <div class="stat-card" title="同分子式竞争候选仍无法点名，但原始谱图、质量和来源证据足以支持该明确分子式异构体组"><div class="stat-value">${nIsomerGroupHigh}</div><div class="stat-label">异构体组 High</div></div>
          <div class="stat-card" title="最高可支持到化学类别，且尚未达到二维连接 High 或异构体组 High；与两者不混计。层级式类别 High 总数为 ${nClassHigh}"><div class="stat-value">${nClassOnlyHigh}</div><div class="stat-label">仅类别 High</div></div>
          <div class="stat-card"><div class="stat-value">${nInsufficient}</div><div class="stat-label">证据不足候选</div></div>
          <div class="stat-card" title="鉴定结论与候选证据等级分开统计"><div class="stat-value">${releasedFinalClaims.length}<span style="font-size:.45em;color:var(--text-secondary);font-weight:400">/${compiledFinalClaims.length}</span></div><div class="stat-label">可作结论 / 已评估</div></div>
          <div class="stat-card" title="连续 MS1 峰形相关 + 中性质量一致的加合物/源内碎片网络"><div class="stat-value">${iimnUnavailable ? '<span style="font-size:.62em;color:#b45309">未运行</span>' : `${ionNetworks.size}<span style="font-size:.45em;color:var(--text-secondary);font-weight:400">（${ionNetworkedFeatures.size}特征）</span>`}</div><div class="stat-label">IIMN 离子身份网络</div></div>
          <div class="stat-card" title="链接注释不等同于独立 MS/MS High"><div class="stat-value">${ionLinked}</div><div class="stat-label">IIMN 链接离子</div></div>
          <div class="stat-card" title="只有稳定的非负 L1 多组分拟合且谱图证据改善才计为成功"><div class="stat-value">${deconvolutionUnavailable ? '<span style="font-size:.62em;color:#b45309">未运行</span>' : `${deconvolvedRows.length}<span style="font-size:.45em;color:var(--text-secondary);font-weight:400">/${chimericRows.length}</span>`}</div><div class="stat-label">嵌合谱去卷积通过</div></div>
        </div>
        <div style="margin-top:.55rem;font-size:.75rem;line-height:1.55;color:var(--text-secondary);">
          当前为<b>${blankViewLabel}</b>：${rate}% 的可见 MS2 有谱库候选。${ms1ScanCount} 张 MS1 与 ${total} 张 MS2 都是扫描次数，不是化合物数；同一色谱峰可被重复 DDA，并在不同碰撞能量下再次采集。${workflowDescription} 当前视图得到 ${groups.length} 个去重谱库候选。“单一二维连接 High”“异构体组 High”“仅类别 High”是不同层级；立体化学／完整结构定义另行报告，IIMN 链接离子不冒充独立 High，同分异构体仍需 RT/标品/CCS 等正交证据解析。
        </div>
        <div style="margin-top:.65rem;padding:.65rem .8rem;border-left:4px solid #d97706;background:#fffbeb;color:#78350f;font-size:.76rem;line-height:1.55;">
          <b>未达到 High 的首要原因（按去重代表身份）</b>：${primaryBlockerHtml}
          <div style="margin-top:.25rem;">未匹配 MS2 ${Number(identificationFunnel.unmatched_ms2_scans || _currentUnmatchedTotal)} 张，其中 ${Number(identificationFunnel.unmatched_with_nonzero_query_peaks || 0)} 张仍有非零清洗后碎片峰；这表示没有候选通过检索/分数门槛，不表示 MS1 中没有色谱特征。</div>
        </div>
        ${deconvolutionImpactHtml}
        ${formalReportConditionHtml}
        ${moduleWarningHtml ? `<details style="margin-top:.65rem;padding:.55rem .75rem;border-left:4px solid #d97706;background:#fffbeb;color:#92400e;font-size:.78rem;line-height:1.55;"><summary style="cursor:pointer;font-weight:700;">⚠ MS1 模块适用性提示（${moduleWarningCount}）<span style="margin-left:.35rem;font-weight:400;color:#a16207;">点击展开</span></summary><div style="margin-top:.35rem;">${moduleWarningHtml}</div></details>` : ''}
        <div id="ticOverviewArea" style="margin-top:0.75rem;width:100%;height:280px;"></div>
      </div>`;

      html += _expectedCompoundsCard(expectedCompounds);

      html += `
      <div class="card">
        <div class="card-header"><div><div class="card-title">处理链与峰保留审计</div><div class="card-subtitle">centroiding、MS1 feature finding 和 MS/MS 检索过滤分别记录，互不混称。</div></div></div>
        <div style="font-size:.77rem;line-height:1.65;color:var(--text-secondary);">
          <div><b>逐谱表征：</b>${Tasks.escapeHtml(String(representation.status || '未记录'))}${representationText ? `；${Tasks.escapeHtml(representationText)}` : ''}</div>
          <div><b>MS1 特征：</b>${Tasks.escapeHtml(String(untargetedWorkflow.status || '未记录'))}；总特征 ${Number(untargetedWorkflow.feature_count || 0).toLocaleString()}；同位素特征 ${Number(untargetedWorkflow.isotope_grouped_feature_count || 0).toLocaleString()}；加合物赋值 ${Number(untargetedWorkflow.adduct_assignment_count || 0).toLocaleString()}（${Number(untargetedWorkflow.adduct_assignment_coverage_pct || 0).toFixed(1)}%），多成员离子身份组 ${Number(untargetedWorkflow.adduct_group_count || 0).toLocaleString()}；MS2 链接 ${lcmsLinkedMs2.toLocaleString()}/${total.toLocaleString()}；共识谱 ${Number(consensusSummary.consensus_spectra || 0).toLocaleString()}。加合物/中性质量均为算法假设，不等同于化合物鉴定。</div>
          <div><b>MS/MS 比较视图：</b>${Number(peakProcessing.input_peaks || 0).toLocaleString()} 峰 → ${Number(peakProcessing.output_peaks || 0).toLocaleString()} 峰（总保留 ${peakProcessing.retained_pct == null ? '—' : Number(peakProcessing.retained_pct).toFixed(1) + '%'}）；查询谱与库谱采用相同参数。</div>
        </div>
        ${peakStageHtml ? `<div style="display:grid;grid-template-columns:repeat(auto-fit,minmax(150px,1fr));gap:7px;margin-top:.65rem;font-size:.72rem;color:#526579;">${peakStageHtml}</div>` : ''}
        ${untargetedWorkflow.warning ? `<div style="margin-top:.6rem;color:#92400e;font-size:.75rem;">⚠ ${Tasks.escapeHtml(String(untargetedWorkflow.warning))}</div>` : ''}
        ${openmsDiagnosticHtml ? `<div style="margin-top:.5rem;font-size:.75rem;line-height:1.5;">${openmsDiagnosticHtml}</div>` : ''}
      </div>`;

      const blankSummary = task.blank_summary || {};
      if (blankSummary.requested) {
        const polarityMismatchCount = Number(
          blankSummary.status_counts?.blank_polarity_not_applicable || 0
        );
        const blankNotApplicable = blankSummary.applicable === false
          || (total > 0 && polarityMismatchCount >= total);
        const samplePolarity = Array.isArray(blankSummary.sample_polarities)
          ? blankSummary.sample_polarities.join('/') : (task.ion_mode || '当前样本');
        const blankPolarity = Array.isArray(blankSummary.blank_polarities)
          ? blankSummary.blank_polarities.join('/') : '所上传空白';
        const mismatchText = `空白过滤未执行：${samplePolarity} 样本与 ${blankPolarity} 空白极性不匹配；“0”不代表已完成扣除。`;
        const warningText = blankNotApplicable
          ? (blankSummary.applicable === false && blankSummary.warning
            ? blankSummary.warning : mismatchText)
          : blankSummary.warning;
        const warning = warningText
          ? `<div style="margin-top:.5rem;color:#991b1b;font-weight:600;">⚠ ${Tasks.escapeHtml(warningText)}</div>`
          : '';
        const blankEvidenceRows = blankFlagged.length ? blankFlagged : blankExcluded;
        const blankRows = blankEvidenceRows.slice(0, 30).map((r, i) => `
          <tr>
            <td>${i + 1}</td><td class="mono">${Tasks.escapeHtml(String(r.spectrum_id || '—'))}</td>
            <td class="mono">${Number(r.precursor_mz || 0).toFixed(4)}</td>
            <td class="mono">${r.retention_time != null ? Number(r.retention_time).toFixed(3) : '—'}</td>
            <td class="mono">${r.blank_sample_ratio != null ? Number(r.blank_sample_ratio).toFixed(3) : '—'}</td>
            <td>${Tasks.escapeHtml(String(r.blank_evidence_category || r.blank_status || '—'))}</td>
            <td>${Tasks.escapeHtml(String(r.blank_evidence_role || 'unspecified'))}</td>
            <td>${r.match?.compound_name ? Tasks.escapeHtml(String(r.match.compound_name)) : '—'}</td>
          </tr>`).join('');
        const softPolicy = blankSummary.filter_policy === 'soft_evidence';
        html += `
        <div class="card" style="border-left:4px solid ${blankNotApplicable ? '#b91c1c' : (blankFlagged.length || blankExcluded.length ? '#b45309' : '#15803d')};background:${blankNotApplicable ? '#fef2f2' : (blankFlagged.length || blankExcluded.length ? '#fffbeb' : '#f0fdf4')};">
          <div class="card-title">🧪 ${blankNotApplicable ? '空白评估未执行' : (softPolicy ? '空白背景证据（完整鉴定后分层）' : '旧版硬过滤')}</div>
          <div style="display:flex;gap:1.2rem;flex-wrap:wrap;font-size:.8rem;line-height:1.7;">
            <span>匹配空白特征 <b>${blankNotApplicable ? '不可计算' : Number(blankSummary.feature_matches || 0)}</b></span>
            <span>完成谱库检索 <b>${blankNotApplicable ? '未执行' : Number(blankSummary.searched_spectra ?? (total - Number(blankSummary.excluded_spectra || 0)))}</b></span>
            <span>空白相关标记 <b>${blankNotApplicable ? '未执行' : Number(blankSummary.flagged_spectra ?? blankFlagged.length)}</b></span>
            <span>清洁视图隐藏 <b>${blankNotApplicable ? '未执行' : Number(blankSummary.clean_view_excluded_spectra ?? blankCleanExcluded.length)}</b></span>
            <span>检索前排除 <b>${blankNotApplicable ? '未执行' : Number(blankSummary.excluded_spectra || 0)}</b></span>
            <span>carryover 证据 <b>${blankNotApplicable ? '未执行' : Number(blankSummary.carryover_suspected_spectra || 0)}</b></span>
            <span>不可评价/未检出 <b>${Number(blankSummary.unevaluable_spectra || 0)}</b></span>
            <span>样品富集：样品/空白 ≥ <b>${Number(blankSummary.sample_to_blank_ratio_threshold || 3).toFixed(1)}</b>，RT ±${Number(blankSummary.rt_tolerance_minutes || .3).toFixed(2)} min</span>
          </div>${warning}
          ${softPolicy ? '<div style="margin-top:.45rem;font-size:.75rem;color:#78350f;">完整鉴定视图不会因空白证据删除化合物；清洁视图只隐藏在明确背景空白中重复出现且样品响应不高于空白的特征。carryover 与未注明类型空白只作标记。</div>' : ''}
          ${blankRows ? `<details style="margin-top:.65rem;"><summary style="cursor:pointer;font-size:.78rem;font-weight:600;">查看空白证据条目（最多显示 30 条）</summary><div class="table-wrap" style="margin-top:.5rem;"><table><thead><tr><th>#</th><th>Scan</th><th>m/z</th><th>RT</th><th>样品/空白</th><th>证据类别</th><th>空白类型</th><th>鉴定候选</th></tr></thead><tbody>${blankRows}</tbody></table></div></details>` : ''}
        </div>`;
      }

      if (contaminantGroups.length > 0) {
        const contaminantNames = contaminantGroups.slice(0, 8)
          .map(g => g.rep.match.compound_name).filter(Boolean).join('、');
        html += `
        <div class="card" style="border-left:4px solid #c2410c;background:#fff7ed;">
          <div class="card-title" style="color:#9a3412;">⚠ 检出 ${contaminantGroups.length} 个潜在外源污染物候选</div>
          <p style="margin:.35rem 0 0;font-size:.8rem;line-height:1.55;color:#7c2d12;">
            ${Tasks.escapeHtml(contaminantNames || '外源污染物候选')}。这些条目可以来自真实实验谱命中，但化合物角色属于耗材、塑料添加剂或环境背景；请结合空白样、前体 XIC 与保留时间复核。
          </p>
        </div>`;
      }

      if (limitedInput.length > 0) {
        const minPeaks = Math.min(...limitedInput.map(r => Number(r.peaks_count || 0)));
        const maxPeaks = Math.max(...limitedInput.map(r => Number(r.peaks_count || 0)));
        const peakRange = minPeaks === maxPeaks ? `${minPeaks}` : `${minPeaks}–${maxPeaks}`;
        html += `
        <div class="card" style="border-left:4px solid #d97706;background:#fffbeb;">
          <div class="card-title" style="color:#92400e;">⚠ 输入是有限 MS2 摘要，不是完整原始谱图</div>
          <p style="margin:.35rem 0 0;font-size:.8rem;line-height:1.55;color:#78350f;">
            此 HTML 报告的 ${limitedInput.length} 张谱图每张仅提供 ${peakRange} 个碎片峰。系统已完成谱库检索，
            系统已自动启用严格证据门槛，避免少峰偶然重叠被列为化合物候选。缺失的低丰度诊断离子会显著降低化合物名级鉴定率；“未匹配”不代表样本中没有化合物。
            下方类别预测只能作为筛查线索。若要可靠点名化合物，请导入报告标题所指原始 <b>.d</b> 转出的完整 mzML/MGF。
          </p>
        </div>`;
      }

      if (groups.length > 0) {
        const dedupHint = groups.length < matched.length
          ? ` · 已按化合物去重 (${groups.length} 个化合物 / 共 ${matched.length} 张匹配谱图)`
          : ' · 点击行查看详情';
        html += `
        <div class="card">
          <div class="card-header">
            <span class="card-title">逐候选鉴定依据表 (${groups.length})</span>
            <span style="font-size:0.78rem;color:var(--text-secondary);">候选不等同于确认鉴定 · 按当前候选排序展示${dedupHint}</span>
          </div>
          <div style="padding:0 0 0.7rem;display:flex;gap:0.5rem;align-items:center;">
            <input id="compoundResultSearch" class="form-input" type="search"
              placeholder="搜索数据库名称、分子式或完整 InChIKey"
              oninput="Results.filterCompoundRows(this.value)" style="max-width:620px;">
            <span id="compoundResultSearchCount" style="font-size:0.75rem;color:var(--text-secondary);"></span>
          </div>
          <div class="table-wrap">
          <table>
            <thead><tr>
              <th>#</th><th>隔离 / 检索 m/z</th><th>化合物名称</th><th>分子式 / 离子形式</th>
              <th title="查询谱与参比谱的碎片相似度">MS/MS 相似度</th><th title="多项匹配证据的综合分，不是鉴定概率">综合匹配分</th><th title="越接近 0 越好">质量误差 (ppm)</th><th>分类</th><th>分析证据</th><th>证据等级</th><th>鉴定结论</th>
            </tr></thead>
            <tbody id="compoundResultsBody">
              ${groups.map((g, i) => {
                const r = g.rep;
                const m = r.match;
                const confCls = m.confidence === 'High' ? 'badge-high' : m.confidence === 'Medium' ? 'badge-medium' : 'badge-low';
                // P4 诚实性校准: 类级置信 + 异构警示 + 诚实标题
                const _cc = m.confidence_class || m.confidence;
                const _gc = m.isomer_group_confidence || 'Unknown';
                const _confLabel = _confidenceLabel(m.confidence, m);
                const _classConfLabel = _confidenceLabel(_cc, m);
                const _groupConfLabel = _confidenceLabel(_gc, m);
                const classConfCls = _cc === 'High' ? 'badge-high' : _cc === 'Medium' ? 'badge-medium' : 'badge-low';
                const groupConfCls = _gc === 'High' ? 'badge-high' : _gc === 'Medium' ? 'badge-medium' : 'badge-low';
                const _tier = m.honesty_tier || '';
                const _eTier = m.evidence_tier || '';
                const _isoN = m.isomer_n || 0;
                const _isoBadge = (_tier === 'isomer_ambiguous' || _isoN >= 1)
                  ? ` <span style="display:inline-block;margin-left:4px;padding:0 6px;border-radius:9px;background:#fef3c7;color:#92400e;font-size:0.62rem;font-weight:700;vertical-align:middle;" title="${Tasks.escapeHtml((m.honest_headline||'').substring(0,200))}">⚠异构×${_isoN+1}</span>`
                  : '';
                const _honest = (_tier === 'isomer_ambiguous' || _tier === 'class_level' || _tier === 'putative');
                // The database identity is a scientific value, not interface
                // copy.  Never decorate it with a translated alias.
                const _databaseName = String(m.compound_name || 'Unnamed candidate');
                const _nameDisp = `<span data-scientific-value style="font-weight:600;">${Tasks.escapeHtml(_databaseName.substring(0, 80))}</span>`;
                const _honestSub = (_honest && m.honest_headline)
                  ? `<div style="font-size:0.62rem;color:#92400e;line-height:1.25;margin-top:2px;overflow:hidden;text-overflow:ellipsis;white-space:nowrap;" title="${Tasks.escapeHtml(m.honest_headline)}">${Tasks.escapeHtml(m.honest_headline.substring(0, 58))}</div>`
                  : '';
                const nSpec = g.nSpectra;
                // 无结构兼容判断；通过严格门控的 name-only 命中是谱图线索，不是匿名鉴定。
                const _hasStruct = _hasStructure(m);
                const noStructBadge = !_hasStruct
                  ? (m.name_only_lead
                    ? ` <span style="display:inline-block;margin-left:4px;padding:0 6px;border-radius:9px;background:#dbeafe;color:#1e40af;font-size:0.62rem;font-weight:700;" title="高相似馆藏谱图线索；无结构，不能作为结构鉴定">谱图线索·无结构</span>`
                    : ` <span style="display:inline-block;margin-left:4px;padding:0 6px;border-radius:9px;background:#fee2e2;color:#991b1b;font-size:0.62rem;font-weight:700;" title="匿名谱图/无结构，无法进行结构鉴定">⚠匿名</span>`)
                  : '';
                const contaminantBadge = m.potential_contaminant
                  ? ` <span style="display:inline-block;margin-left:4px;padding:0 6px;border-radius:9px;background:#ffedd5;color:#9a3412;font-size:0.62rem;font-weight:700;" title="${Tasks.escapeHtml(m.contaminant_note || '潜在外源污染物；不计入天然产物')}">⚠外源污染物</span>`
                  : '';
                const multiBadge = nSpec > 1
                  ? ` <span style="display:inline-block;margin-left:4px;padding:0 6px;border-radius:9px;background:#dcfce7;color:#166534;font-size:0.62rem;font-weight:700;vertical-align:middle;" title="该化合物共命中 ${nSpec} 张谱图(已合并展示)">×${nSpec}谱</span>`
                  : '';
                const evidenceBadge = (_eTier && !['High', 'Putative'].includes(_eTier))
                  ? ` <span style="display:inline-block;margin-left:4px;padding:0 6px;border-radius:9px;background:${_eTier === 'Strong Candidate' ? '#fef3c7' : '#dbeafe'};color:${_eTier === 'Strong Candidate' ? '#92400e' : '#1e40af'};font-size:0.62rem;font-weight:700;vertical-align:middle;" title="${Tasks.escapeHtml(m.annotation_level || _eTier)}">${Tasks.escapeHtml(_eTier)}</span>`
                  : '';
                const chimericWarningBadge = m.chimeric_isolation_warning
                  ? ` <span style="display:inline-block;margin-left:4px;padding:0 6px;border-radius:9px;background:#fff7ed;color:#9a3412;font-size:0.62rem;font-weight:700;vertical-align:middle;" title="${Tasks.escapeHtml(m.chimeric_isolation_warning)}">⚠嵌合·原谱置信</span>`
                  : '';
                const _ionFlags = m.display_evidence_flags || {};
                const _sampleIon = String(_ionFlags.sample_ion_form || '').trim();
                const _libraryIon = String(_ionFlags.library_adduct || m.adduct || '').trim();
                const _candidateIon = String(
                  _ionFlags.candidate_ion_form_hypothesis
                  || m.candidate_adduct_hypothesis || ''
                ).trim();
                const _ionCell = `<b class="mono">${Tasks.escapeHtml(String(m.formula || '—'))}</b>`
                  + `<div style="font-size:.62rem;line-height:1.35;margin-top:2px;color:${_sampleIon ? '#166534' : '#92400e'};">样品 ${Tasks.escapeHtml(_sampleIon || '未解析')}</div>`
                  + `<div style="font-size:.62rem;line-height:1.35;color:#64748b;">谱库 ${Tasks.escapeHtml(_libraryIon || '未记录')}</div>`
                  + (_candidateIon
                    ? `<div style="font-size:.62rem;line-height:1.35;color:#1d4ed8;" title="由当前候选结构/分子式与谱库前体反推；不是样品独立 adduct 证据">候选假设 ${Tasks.escapeHtml(_candidateIon)}</div>`
                    : '');
                const _searchText = [
                  m.compound_name, m.source_compound_name,
                  m.formula, m.inchikey, m.classification, _sampleIon, _libraryIon,
                  _candidateIon
                ].filter(Boolean).join(' ').toLowerCase();
                return `
                <tr class="compound-result-row" data-search="${Tasks.escapeHtml(_searchText)}" style="cursor:pointer;" onclick="Results.showDetail('${taskId}', ${g.repIdx})">
                  <td>${i + 1}</td>
                  <td class="mono">${_precursorMassCell(r, 4)}${r.blank_feature_matched ? ` <span style="font-size:0.62rem;color:${r.blank_flagged ? '#b45309' : '#15803d'};" title="空白证据：${Tasks.escapeHtml(String(r.blank_evidence_category || r.blank_status || '已评价'))}；类型 ${Tasks.escapeHtml(String(r.blank_evidence_role || 'unspecified'))}；样品/空白 ${r.blank_sample_ratio != null ? Number(r.blank_sample_ratio).toFixed(2) : '不可计算'}${r.blank_clean_view_excluded ? '；仅从清洁视图隐藏' : ''}">🧪${r.blank_sample_ratio != null ? Number(r.blank_sample_ratio).toFixed(1) : '?'}</span>` : ''}</td>
                  <td style="font-weight:500;max-width:240px;overflow:hidden;text-overflow:ellipsis;white-space:nowrap;"
                      data-scientific-name="${Tasks.escapeHtml(_databaseName)}" title="${Tasks.escapeHtml(_databaseName)}${nSpec > 1 ? ` (共 ${nSpec} 张谱图)` : ''}${!_hasStruct ? (m.name_only_lead ? ' [谱图线索/无结构]' : ' [匿名谱图/无结构-不可信]') : ''}${_isoN>=1 ? ` [${Tasks.escapeHtml((m.honest_headline||'').substring(0,150))}]` : ''}">${_nameDisp}${contaminantBadge}${noStructBadge}${multiBadge}${_isoBadge}${evidenceBadge}${chimericWarningBadge}${_honestSub}</td>
                  <td>${_ionCell}</td>
                  <td class="mono" style="color:#276749;font-weight:700;">${m.cosine.toFixed(4)}</td>
                  <td class="mono" style="color:#9C4221;font-weight:700;">${m.final_score.toFixed(4)}</td>
                  <td class="mono">${m.mass_ppm != null ? m.mass_ppm.toFixed(2) : '—'}</td>
                  <td>${Tasks.escapeHtml(_classificationLabel(m))}</td>
                  <td style="min-width:130px;line-height:1.35;">${_analysisBadges(m, _currentIonMode)}</td>
                  <td><span class="badge ${confCls}" title="具体结构置信度">具体:${Tasks.escapeHtml(_confLabel)}</span>${m.isomer_ambiguous ? `<br><span class="badge ${groupConfCls}" style="font-size:0.6rem;" title="同分子式竞争候选作为一个异构体组的置信度">异构组:${Tasks.escapeHtml(_groupConfLabel)}</span>` : ''}<br><span class="badge ${classConfCls}" style="font-size:0.6rem;opacity:0.8;" title="化学类别置信度">类别:${Tasks.escapeHtml(_classConfLabel)}</span></td>
                  <td>${_finalClaimBadge(m)}</td>
                </tr>`;
              }).join('')}
            </tbody>
          </table>
          </div>
        </div>`;
      }

      html += _governanceBlockedCard(governanceBlockedCandidates, taskId);

      if (_currentUnmatchedTotal > 0) {
        // MSClassify 的 softmax 最大概率未经样本域校准，不能把 0.30 冒充“高置信”。
        // 只把 >=0.70 称为模型强信号；0.45~0.70 中等；其余明确为低。
        const _tallyStrong = {}, _tallyMedium = {}, _tallyLow = {};
        let _wc = 0, _nStrong = 0, _nMedium = 0, _nLow = 0;
        unmatched.forEach(r => {
          const c = r.display_class;
          if (!c) return;
          _wc++;
          const p = Number(r.dl_confidence || 0);
          if (p >= 0.70) { _tallyStrong[c] = (_tallyStrong[c] || 0) + 1; _nStrong++; }
          else if (p >= 0.45) { _tallyMedium[c] = (_tallyMedium[c] || 0) + 1; _nMedium++; }
          else { _tallyLow[c] = (_tallyLow[c] || 0) + 1; _nLow++; }
        });
        const _keys = [...new Set([...Object.keys(_tallyStrong), ...Object.keys(_tallyMedium), ...Object.keys(_tallyLow)])];
        const _badges = _keys.sort((a, b) => ((_tallyStrong[b]||0)+(_tallyMedium[b]||0)+(_tallyLow[b]||0)) - ((_tallyStrong[a]||0)+(_tallyMedium[a]||0)+(_tallyLow[a]||0)))
          .map(c => {
            const strong = _tallyStrong[c] || 0, medium = _tallyMedium[c] || 0, low = _tallyLow[c] || 0;
            const total = strong + medium + low;
            const detail = [strong ? `${strong}强` : '', medium ? `${medium}中` : '', low ? `${low}低` : ''].filter(Boolean).join('/');
            return `<span style="display:inline-block;margin:2px 4px;padding:2px 9px;border-radius:10px;background:#e8f0fe;color:#1a56c4;font-size:0.74rem;" title="模型概率分层：>=0.70强，0.45–0.70中，<0.45低">${Tasks.escapeHtml(c)} · ${total} <span style="font-size:.64rem;color:#64748b;">(${detail})</span></span>`;
          }).join('');
        const _unmatchedRows = unmatched.map((r, i) => {
          const cls = r.display_class || '未能判类';
          const prob = Number(r.dl_confidence || 0);
          const confText = prob > 0 ? `${(prob * 100).toFixed(1)}%` : '—';
          const confColor = prob >= 0.70 ? '#166534' : prob >= 0.45 ? '#1d4ed8' : prob >= 0.15 ? '#92400e' : '#64748b';
          const top3 = (r.dl_top3 || []).map(x => `${x[0]} ${(Number(x[1] || 0) * 100).toFixed(1)}%`).join(' · ');
          return `<tr style="cursor:pointer;" onclick="Results.showUnmatchedDetail(${i})" title="点击查看碎片谱图和 Top-3 类别概率">
            <td>${i + 1}</td>
            <td class="mono">${Tasks.escapeHtml(String(r.spectrum_id || '—'))}</td>
            <td class="mono">${Number(r.precursor_mz || 0).toFixed(4)}</td>
            <td class="mono">${r.retention_time != null ? Number(r.retention_time).toFixed(3) : '—'}</td>
            <td>${Number(r.peaks_count || 0)}</td>
            <td><b>${Tasks.escapeHtml(cls)}</b></td>
            <td style="color:${confColor};font-weight:700;">${confText}</td>
            <td style="max-width:360px;white-space:normal;font-size:.7rem;color:var(--text-secondary);">${Tasks.escapeHtml(top3 || '—')}</td>
          </tr>`;
        }).join('');
        html += `
        <div class="card">
          <div class="card-title" style="margin-bottom:0.4rem;">未获得谱库候选的 MS2 (${_currentUnmatchedTotal})</div>
          <p style="font-size:0.78rem;color:var(--text-secondary);">以下先显示 ${unmatched.length}/${_currentUnmatchedTotal} 张，避免一次创建数千个 DOM 行。下方<b>分类来自 MS2 谱图推断(MSClassify)</b>，性质不同于上方“匹配化合物”的库分类，仅供参考。</p>
          ${_wc > 0 ? `<div style="margin-top:0.5rem;font-size:0.8rem;"><b>🧪 当前已显示谱图判类</b> <span style="color:var(--text-secondary);">(MSClassify 未校准模型概率, ${_wc}/${unmatched.length} 有输出 · ${_nStrong} 强 / ${_nMedium} 中 / <span style="color:#b9770e;">${_nLow} 低</span>)</span>:<br>${_badges}</div>` : `<div style="margin-top:0.5rem;font-size:0.8rem;color:var(--text-secondary);">当前页无谱图可判类(峰数过少或信号不足)。</div>`}
          <div class="table-wrap" style="margin-top:.65rem;">
            <table>
              <thead><tr><th>#</th><th>Scan</th><th>前体 m/z</th><th>RT(min)</th><th>碎片数</th><th>预测类别</th><th>概率</th><th>Top-3 类别</th></tr></thead>
              <tbody>${_unmatchedRows}</tbody>
            </table>
          </div>
          ${unmatched.length < _currentUnmatchedTotal ? `<div style="margin-top:.55rem;text-align:center;"><button class="btn btn-secondary btn-sm" onclick="Results.showMoreUnmatched()">再显示 50 张（已显示 ${unmatched.length}/${_currentUnmatchedTotal}）</button></div>` : ''}
          <div style="margin-top:.4rem;font-size:.7rem;color:var(--text-secondary);">点击任一行可查看该谱图的全部可用碎片和分类概率；概率不是化合物鉴定置信度。</div>
        </div>`;
      }

      if (area) area.innerHTML = html;

      // 渲染 TIC 全扫描概览 (优先用后端 tic_overview.json: MS1 连续色谱图 + MS2 红色匹配点;
      // 旧任务无该文件则 fallback 用 results 里的 MS2 stem 兜底)
      setTimeout(async function () {
        try {
          const ticEl = document.getElementById("ticOverviewArea");
          if (!ticEl || typeof SpectrumViewer === "undefined") return;

          // 1) 先尝试新 API
          let overview = null;
          try {
            overview = await API.get(`/tasks/${taskId}/tic_overview`);
          } catch (e) { /* 旧任务 404, 继续走 fallback */ }

          if (overview && (overview.ms1?.length || overview.ms2?.length)) {
            // 用 results 里的 spectrum_id 标记哪些 MS2 已匹配 + 化合物名
            const matchedIds = new Set();
            const nameMap = new Map();
            results.forEach(r => {
              if (r.match && r.match.compound_name && r.spectrum_id != null) {
                const sidRaw = String(r.spectrum_id);
                const sidNum = Number(r.spectrum_id);
                matchedIds.add(sidRaw);
                nameMap.set(sidRaw, r.match.compound_name);
                if (Number.isFinite(sidNum)) {
                  matchedIds.add(sidNum);
                  nameMap.set(sidNum, r.match.compound_name);
                }
              }
            });
            new SpectrumViewer("ticOverviewArea").renderTIC({
              ms1: overview.ms1 || [],
              ms2: overview.ms2 || [],
              matched_scan_ids: matchedIds,
              scan_name_map: nameMap,
            });
            return;
          }

          // 2) Fallback: 用 results 里的 retention_time + tic 画 stem
          const points = results
            .filter(r => r.retention_time != null && r.tic != null)
            .map(r => ({
              rt: Number(r.retention_time) || 0,
              tic: Number(r.tic) || 0,
              matched: !!(r.match && r.match.compound_name),
              name: r.match ? r.match.compound_name : "",
            }))
            .filter(p => p.tic > 0)
            .sort((a, b) => a.rt - b.rt);
          if (points.length) {
            new SpectrumViewer("ticOverviewArea").renderTIC(points);
          } else {
            ticEl.innerHTML = '<div style="padding:2rem;text-align:center;color:#999;font-size:0.85rem;">无 RT/TIC 数据 (旧分析任务请重新分析以填充)</div>';
          }
        } catch (e) {
          console.warn('TIC overview render failed:', e);
        }
      }, 60);
    } catch (e) {
      if (area) area.innerHTML = `
        <div class="empty-state"><div class="empty-icon">❌</div>
        <div class="empty-title">加载失败</div>
        <div class="empty-desc">${Tasks.escapeHtml(e.message)}</div></div>`;
    }
  }

  function showMoreUnmatched() {
    _unmatchedDisplayLimit += 50;
    loadForTask();
  }

  async function showUnmatchedDetail(index) {
    const summaryRow = _currentUnmatched[index];
    if (!summaryRow) {
      Toast.show('未匹配谱图数据未加载，请重新选择任务', 'error');
      return;
    }
    let r = summaryRow;
    if (Number.isInteger(Number(summaryRow.result_index))) {
      try {
        const payload = await API.get(
          `/tasks/${encodeURIComponent(_currentTaskId)}/results/item/${Number(summaryRow.result_index)}`
        );
        r = payload.result || summaryRow;
      } catch (error) {
        Toast.show(`完整谱图加载失败：${error.message || 'unknown error'}`, 'error');
        return;
      }
    }
    const cls = r.display_class || '未能判类';
    const prob = Number(r.dl_confidence || 0);
    const top3 = r.dl_top3 || [];
    const topRows = top3.length ? top3.map((x, i) => {
      const p = Number(x[1] || 0);
      return `<tr><td>${i + 1}</td><td><b>${Tasks.escapeHtml(String(x[0] || ''))}</b></td><td class="mono">${(p * 100).toFixed(2)}%</td></tr>`;
    }).join('') : '<tr><td colspan="3">无可用类别概率</td></tr>';
    const peaks = r.query_peaks || [];
    const rtMinutes = _retentionTimeMinutes(r);
    const peakRows = peaks.map((p, i) => `<tr><td>${i + 1}</td><td class="mono">${Number(p[0]).toFixed(4)}</td><td class="mono">${Number(p[1]).toFixed(2)}</td></tr>`).join('');
    const html = `
      <div class="modal-overlay" onclick="if(event.target===this)App.closeModal()">
      <div class="modal modal-lg">
        <div class="flex-between" style="margin-bottom:.65rem;">
          <div><h3 style="font-size:1.05rem;">未匹配谱图 · Scan ${Tasks.escapeHtml(String(r.spectrum_id || '—'))}</h3>
          <div style="font-size:.8rem;color:var(--text-secondary);">源文件 ${Tasks.escapeHtml(_currentTaskFilename || '未记录')} · RT ${rtMinutes == null ? '未记录' : rtMinutes.toFixed(3) + ' min'} · ${Tasks.escapeHtml(_precursorMassInline(r, 5))} · ${peaks.length} 个可用碎片</div></div>
          <span class="badge ${prob >= .70 ? 'badge-high' : prob >= .45 ? 'badge-medium' : 'badge-low'}">${Tasks.escapeHtml(cls)} · 模型概率 ${(prob * 100).toFixed(1)}%</span>
        </div>
        ${r.input_evidence_note ? `<div style="padding:.45rem .6rem;background:#fffbeb;border:1px solid #fcd34d;border-radius:5px;color:#92400e;font-size:.76rem;margin-bottom:.65rem;">⚠ ${Tasks.escapeHtml(r.input_evidence_note)}</div>` : ''}
        <div style="display:grid;grid-template-columns:minmax(240px,.8fr) minmax(300px,1.2fr);gap:1rem;align-items:start;">
          <div><div style="font-weight:700;font-size:.8rem;margin-bottom:.35rem;">MSClassify Top-3（类别线索）</div>
            <table><thead><tr><th>#</th><th>类别</th><th>概率</th></tr></thead><tbody>${topRows}</tbody></table>
          </div>
          <div><div style="font-weight:700;font-size:.8rem;margin-bottom:.35rem;">HTML 中实际可用碎片</div>
            <div class="table-wrap" style="max-height:210px;"><table><thead><tr><th>#</th><th>m/z</th><th>相对强度</th></tr></thead><tbody>${peakRows || '<tr><td colspan="3">无碎片</td></tr>'}</tbody></table></div>
          </div>
        </div>
        <details class="pw-card" open style="margin-top:.75rem;"><summary>📊 MS/MS 谱图</summary><div id="unmatchedSpectrum" style="width:100%;height:280px;"></div></details>
        <div style="margin-top:.55rem;font-size:.72rem;color:var(--text-secondary);">未匹配表示当前谱库和阈值下没有足够相似的同极性参比谱；Top-3 仅为宽类别模型输出，不能替代化合物名称鉴定。</div>
        <div class="modal-actions"><button class="btn btn-secondary" onclick="App.closeModal()">关闭</button></div>
      </div></div>`;
    App.showModal(html);
    setTimeout(function () {
      const el = document.getElementById('unmatchedSpectrum');
      if (!el) return;
      if (peaks.length && typeof SpectrumViewer !== 'undefined') {
        new SpectrumViewer('unmatchedSpectrum').render(peaks, { title: `Unmatched m/z ${Number(r.precursor_mz || 0).toFixed(4)}` });
      } else {
        el.innerHTML = '<div style="padding:1rem;text-align:center;color:#999;">无 MS/MS 碎片数据</div>';
      }
    }, 80);
  }

  async function _loadRemoteIsomerLiterature(request) {
    const el = document.getElementById(request.containerId);
    if (!el) return;
    const directUrl = String(request.searchUrl || '');
    const directLink = directUrl.startsWith('https://pubmed.ncbi.nlm.nih.gov/')
      ? `<a href="${Tasks.escapeHtml(directUrl)}" target="_blank" rel="noopener noreferrer">在 PubMed 打开检索式</a>`
      : '';
    el.innerHTML = '<span style="color:var(--text-secondary);">正在通过 NCBI PubMed 检索候选对；此步骤不阻塞分析队列…</span>';
    try {
      const response = await fetch(`${API.BASE}/isomer-literature/search`, {
        method: 'POST',
        headers: {'Content-Type': 'application/json'},
        body: JSON.stringify({
          candidate_names: request.candidateNames,
          formula: request.formula || '',
          max_results: 5,
        }),
      });
      if (!response.ok) {
        throw new Error(`HTTP ${response.status}`);
      }
      const payload = await response.json();
      const records = Array.isArray(payload.records) ? payload.records : [];
      if (!records.length) {
        const statusText = payload.status === 'search_terms_insufficient'
          ? '候选名称不足以形成可靠的化合物检索式'
          : payload.status === 'lookup_temporarily_unavailable'
            ? 'PubMed 暂时不可用，本次未把网络故障误写成“无文献”'
            : '已检索，但未发现可直接对应这组候选的文献线索';
        el.innerHTML = `<span style="color:#92400e;">${Tasks.escapeHtml(statusText)}。</span> ${directLink}`;
        return;
      }
      const rows = records.map(record => {
        const url = String(record.url || '');
        const title = Tasks.escapeHtml(record.title || `PMID ${record.pmid || ''}`);
        const linkedTitle = url.startsWith('https://pubmed.ncbi.nlm.nih.gov/')
          ? `<a href="${Tasks.escapeHtml(url)}" target="_blank" rel="noopener noreferrer">${title}</a>`
          : title;
        const ids = [
          record.pmid ? `PMID ${record.pmid}` : '',
          record.doi ? `DOI ${record.doi}` : '',
        ].filter(Boolean).join(' · ');
        const tags = (record.evidence_tags || []).join(', ');
        const pairSpecific = record.evidence_scope
          === 'pair_specific_bibliographic_method_lead';
        const scopeLabel = pairSpecific
          ? '题录/摘要同时命中至少两个不同候选'
          : '仅命中一个候选或一般方法线索';
        const candidateHits = Array.isArray(record.candidate_name_hits)
          ? record.candidate_name_hits.filter(Boolean).join(' / ') : '';
        return `<li style="margin:.25rem 0;">${linkedTitle}
          <div style="color:var(--text-secondary);font-size:.66rem;">${Tasks.escapeHtml([record.year, record.journal, ids].filter(Boolean).join(' · '))}${tags ? ` · 范围标签: ${Tasks.escapeHtml(tags)}` : ''}</div>
          <div style="font-size:.66rem;color:${pairSpecific ? '#166534' : '#92400e'};">${Tasks.escapeHtml(scopeLabel)}${candidateHits ? ` · 命中：${Tasks.escapeHtml(candidateHits)}` : ''}</div>
        </li>`;
      }).join('');
      const distinctNames = Array.from(new Set(
        (request.candidateNames || [])
          .map(name => String(name || '').trim().toLocaleLowerCase())
          .filter(Boolean)
      ));
      const nameScopeWarning = distinctNames.length < 2
        ? '<div style="color:#92400e;margin-bottom:.2rem;">⚠ 当前歧义组只有一个可区分的候选名称；检索结果只能作为单化合物/一般方法线索，不能证明候选对已被文献区分。</div>'
        : '';
      el.innerHTML = `<div style="color:#166534;margin-bottom:.2rem;">PubMed 返回 ${records.length} 条可追溯文献线索${payload.cache_hit ? '（缓存）' : ''}：</div>
        ${nameScopeWarning}
        <ol style="margin:.15rem 0 .25rem 1rem;padding:0;">${rows}</ol>
        <div style="color:#92400e;">⚠ 只完成题录/摘要范围筛查；尚未从全文提取可比 CE、仪器、加合物或实测诊断离子，因而不参与评分或精确身份升级。${directLink ? ` ${directLink}` : ''}</div>`;
    } catch (error) {
      el.innerHTML = `<span style="color:#b91c1c;">PubMed 联机检索失败（${Tasks.escapeHtml(error.message || 'network error')}），未把失败解释成“没有文献”。</span> ${directLink}`;
    }
  }

  function _mountInlineDetail(html, targetId) {
    const target = document.getElementById(targetId);
    if (!target) return false;
    const template = document.createElement('template');
    template.innerHTML = String(html || '').trim();
    const detail = template.content.querySelector('.modal');
    if (!detail) return false;
    detail.querySelector('.modal-actions')?.remove();
    detail.classList.remove('modal', 'modal-lg');
    detail.classList.add('candidate-inline-detail');
    target.replaceChildren(detail);
    target.scrollTop = 0;
    return true;
  }

  // ── Match Detail (modal compatibility + inline workspace) ─────────
  async function showDetail(taskId, index, inlineTargetId = null, summaryRowOverride = null) {
    const renderSequence = ++_detailRenderSequence;
    try {
      const summaryRow = summaryRowOverride || _currentMatched[index];
      if (!summaryRow) {
        Toast.show("结果数据未加载, 请重新选择任务", "error");
        return;
      }
      const resultIndex = Number(summaryRow.result_index);
      const detailPayload = Number.isInteger(resultIndex)
        ? await API.get(`/tasks/${encodeURIComponent(taskId)}/results/item/${resultIndex}`)
        : { result: summaryRow };
      if (renderSequence !== _detailRenderSequence) return;
      const r = detailPayload.result || summaryRow;
      const m = r.match;
      if (!m) return;

      // 该化合物在本任务内命中的总谱图数(折叠视图的补充, 让用户知道主表合并了哪些谱图)
      const _cmpKey = (mm) => _compoundIdentityKey(mm);
      const _nSpecThis = _currentMatched.filter(rr => rr.match && _cmpKey(rr.match) === _cmpKey(m)).length;

      const confCls = m.confidence === 'High' ? 'badge-high' : m.confidence === 'Medium' ? 'badge-medium' : 'badge-low';
      const confLabel = _confidenceLabel(m.confidence, m);
      const isIsomerGroupScope = String(m.identification_scope || '') === 'isomer_group'
        || !!m.isomer_ambiguous;
      const isomerGroupConfidence = String(m.isomer_group_confidence || 'Unknown');
      const exactStructureConfidence = String(
        m.exact_structure_confidence || m.confidence || 'Unknown'
      );
      const isomerGroupConfCls = isomerGroupConfidence === 'High'
        ? 'badge-high' : isomerGroupConfidence === 'Medium' ? 'badge-medium' : 'badge-low';
      const exactStructureConfCls = exactStructureConfidence === 'High'
        ? 'badge-high' : exactStructureConfidence === 'Medium' ? 'badge-medium' : 'badge-low';
      const headerConfidenceHtml = isIsomerGroupScope
        ? `<div style="display:flex;gap:.35rem;align-items:center;flex-wrap:wrap;justify-content:flex-end;">
            <span class="badge ${isomerGroupConfCls}" style="font-size:0.72rem;padding:3px 10px;" title="同分子式竞争候选作为一个异构体组的证据等级">异构体组：${Tasks.escapeHtml(_confidenceLabel(isomerGroupConfidence, m))}</span>
            <span class="badge ${exactStructureConfCls}" style="font-size:0.72rem;padding:3px 10px;" title="当前 Top1 候选的完整结构身份证据等级；不等于异构体组等级">Top1 完整结构：${Tasks.escapeHtml(_confidenceLabel(exactStructureConfidence, m))}</span>
          </div>`
        : `<span class="badge ${confCls}" style="font-size:0.78rem;padding:3px 12px;">${Tasks.escapeHtml(confLabel)}</span>`;
      const classLabel = _classificationLabel(m);
      const smiles = m.smiles || '';
      const rtMinutes = _retentionTimeMinutes(r);
      const rtText = rtMinutes == null ? '未记录' : `${rtMinutes.toFixed(3)} min`;
      const spectrumId = String(r.spectrum_id || '').trim() || '未记录';
      const savedScoringView = typeof SpectrumComparator !== 'undefined'
        ? SpectrumComparator.savedEvidence?.(m) || null : null;
      const detailQueryPeaks = savedScoringView ? savedScoringView.query_peaks
        : Array.isArray(r.query_peaks) ? r.query_peaks
        : Array.isArray(r.peaks) ? r.peaks : [];
      const detailLibraryPeaks = savedScoringView ? savedScoringView.library_peaks
        : Array.isArray(r.lib_peaks) ? r.lib_peaks
        : Array.isArray(m.lib_peaks) ? m.lib_peaks
          : Array.isArray(m.library_peaks) ? m.library_peaks : [];
      const detailMatchedPairs = savedScoringView ? savedScoringView.matched_pairs
        : Array.isArray(r.matched_peaks) ? r.matched_peaks
        : Array.isArray(m.matched_peaks) ? m.matched_peaks : [];
      const hasDetailMirrorEvidence = detailLibraryPeaks.length > 0;
      const exactRtStandard = m.rt_evidence?.eligible_for_exact_identity === true;
      const detailIdentityLabel = exactRtStandard
        ? '标准品确证（同方法 RT）'
        : isIsomerGroupScope
          ? '同分子式异构体候选组'
          : 'MS/MS 谱库二维结构候选';
      const detailIdentityNote = exactRtStandard
        ? `标准品 ${m.rt_evidence?.standard_name || '未记录名称'} · ΔRT ${m.rt_evidence?.rt_delta_sec != null ? m.rt_evidence.rt_delta_sec + ' s' : '未记录'}`
        : isIsomerGroupScope
          ? '结构图为当前 Top 1 排序候选，未解析到单一异构体。'
          : '该名称与结构是谱库候选；相似度和综合分均不是身份概率。';

      // Structure image
      let structHtml = '<span class="text-muted" style="font-size:0.85rem;">无结构式</span>';
      if (smiles) {
        structHtml = `<img src="${API.BASE}/tasks/${taskId}/download/structure/result/${Number.isInteger(resultIndex) ? resultIndex : index}?format=svg"
          class="structure-img" alt="Structure"
          onerror="this.outerHTML='<span class=\\'text-muted\\'>结构式不可用</span>'">`;
      }

      // Expert details
      let expertHtml = '';
      let isomerLiteratureRequest = null;
      const ed = m.expert_details || {};
      const expertScore = m.expert_score || 0;
      const auxiliaryInterpretation = m.auxiliary_interpretation || {};
      const auxiliaryExpert = auxiliaryInterpretation.expert || {};
      const adductEvidence = auxiliaryInterpretation.adduct_assignment || {};
      const expertCoverage = auxiliaryExpert.coverage || {};
      const evographHtml = _evographFragmentationHtml(
        r.evograph_fragmentation
      );
      // v3 Phase 4: 结构推理链 (来自 v3 裂解规则库)
      let v3ReasoningHtml = '';
      const v3 = (ed && ed.v3) || {};
      if (v3.reasoning && v3.reasoning.length) {
        const ions = (v3.matched_ions || []).map(i => i.mz.toFixed(2)).join(', ');
        const losses = (v3.matched_losses || []).map(l => `${l.loss_da.toFixed(1)}(${l.annotation || '?'})`).join(', ');
        v3ReasoningHtml = `<details class="pw-card" open><summary>🧪 结构推理 (v3 规则库 · 一致性 ${v3.v3_score != null ? v3.v3_score.toFixed(2) : '—'})</summary>
          <div style="padding:0.5rem 0.6rem;font-size:0.8rem;">
            ${ions ? `<div style="margin-bottom:4px;">🎯 命中诊断离子: <span class="mono">${Tasks.escapeHtml(ions)}</span></div>` : ''}
            ${losses ? `<div style="margin-bottom:4px;">➖ 命中中性丢失: <span class="mono">${Tasks.escapeHtml(losses)}</span></div>` : ''}
            <div style="margin-top:6px;color:var(--text-secondary);">${v3.reasoning.map(r => '• ' + Tasks.escapeHtml(r)).join('<br>')}</div>
          </div></details>`;
      }
      // v3 Phase 5 / Path1: Top-K 近邻候选 (按化合物去重, 含异构体, 首位=当前 best) + 异构体歧义警示
      let topCandHtml = '';
      let isopairHtml = '';
      const isopair = m.isopair_evidence || {};
      if (Object.keys(isopair).length) {
        const statusLabels = {
          applied: '已进入近并列适用域',
          clear_baseline_winner_no_rerank: '已评估·原始胜者清晰，不重排',
          not_applicable: '不在适用域',
          evaluation_failed: '模型评估失败',
        };
        const gateLabels = {
          supported: 'High 精度门通过',
          vetoed_to_medium: '精度门未通过，已限制为 Medium',
          evaluated_not_high_candidate: '已评估，本条不是待发布 High',
          rt_standard_supersedes_model_gate: '同方法 RT 标准品优先于模型门',
          not_applicable: '不适用',
        };
        const numberText = (value, digits = 3) => {
          const parsed = Number(value);
          return Number.isFinite(parsed) ? parsed.toFixed(digits) : '—';
        };
        const scoreRows = (isopair.candidate_scores || []).slice(0, 12).map((candidate, rank) => `
          <tr${candidate.connectivity === isopair.selected_connectivity ? ' style="background:#f0fdf4;"' : ''}>
            <td>${rank + 1}</td>
            <td class="mono">${Tasks.escapeHtml(candidate.connectivity || '—')}</td>
            <td class="mono">${numberText(candidate.baseline_modified_cosine)}</td>
            <td class="mono">${numberText(candidate.same_identity_probability)}</td>
            <td class="mono"><b>${numberText(candidate.fused_score)}</b></td>
            <td>${Number(candidate.reference_spectra_scored || 0)}</td>
          </tr>`).join('');
        const changed = !!isopair.ranking_changed;
        const gateResult = m.isopair_high_precision_gate_result
          || (isopair.high_precision_support ? 'supported' : 'not_applicable');
        const border = changed ? '#7c3aed'
          : gateResult === 'vetoed_to_medium' ? '#dc2626' : '#0f766e';
        isopairHtml = `<details class="pw-card" open>
          <summary>🧩 IsoPair-MS 同分子式 2-D 连接异构体判别${changed ? ' · 已重排' : ''}</summary>
          <div style="padding:.5rem .65rem;border-left:4px solid ${border};font-size:.72rem;line-height:1.5;">
            <div><b>${Tasks.escapeHtml(statusLabels[isopair.status] || isopair.status || '未评估')}</b>${isopair.reason ? ` · ${Tasks.escapeHtml(isopair.reason)}` : ''}</div>
            <div>适用组：${Tasks.escapeHtml(isopair.formula || '—')} · ${Tasks.escapeHtml(isopair.adduct || '—')} · ${Number(isopair.connectivity_candidate_count || 0)} 个 connectivity</div>
            <div>原始 cosine 分差 ${numberText(isopair.baseline_margin)} / 重排上限 ${numberText(isopair.maximum_influence_margin)} · 融合胜者 ${numberText(isopair.fused_winner_score)} · 融合分差 ${numberText(isopair.fused_winner_margin)}</div>
            <div><b>${Tasks.escapeHtml(gateLabels[gateResult] || gateResult)}</b> · 排名改变 ${changed ? '是' : '否'} · 立体化学声明：否</div>
            ${scoreRows ? `<div class="table-wrap" style="margin-top:.4rem;"><table><thead><tr><th>#</th><th>InChIKey-14</th><th>Modified Cos</th><th>同连接概率</th><th>融合分</th><th>参比谱</th></tr></thead><tbody>${scoreRows}</tbody></table></div>` : ''}
            <div style="margin-top:.35rem;color:#92400e;">该模块只在同分子式、同加合物、2–12 个实验谱候选的适用域内工作；训练标签为 InChIKey 前 14 位的 2-D connectivity。重排后的近并列候选不会被自动升为 High，也不声称区分一般立体异构体。</div>
          </div></details>`;
      }
      const topCRaw = m.top_candidates || [];
      // 兜底去重: 按完整 InChIKey（或结构/name+formula）去重, 防止旧任务 results.json 里
      // 同物多谱生成的重复候选(当前 search.py 已在生成时去重, 此处兼容旧任务陈旧数据)。
      // 例: 旧数据里 34_dihydroxyCinnamicAcid 重复5次 → 去重后只显示1次。
      const _seenK = new Set();
      const topC = topCRaw.filter(c => {
        const k = _compoundIdentityKey(c);
        if (_seenK.has(k)) return false;
        _seenK.add(k); return true;
      });
      if (topC.length > 1) {
        // 异构体歧义: 同分子式竞争结构尚未被当前证据解析。
        const f0 = (topC[0].formula || '').trim();
        const f1 = (topC[1].formula || '').trim();
        const hasExplicitIsomerState = (
          typeof m.isomer_ambiguous === 'boolean'
        );
        const legacyFormulaFallback = !!(
          f0 && f0 === f1
          && _compoundIdentityKey(topC[0]) !== _compoundIdentityKey(topC[1])
        );
        const isoAmbig = hasExplicitIsomerState
          ? m.isomer_ambiguous
          : legacyFormulaFallback;
        const _shadowMetric = (c, key) => {
          const metric = ((c.multimetric_shadow || {})[key] || {});
          const value = Number(metric.score);
          return metric.status === 'available' && Number.isFinite(value)
            ? value.toFixed(3) : '—';
        };
        const ambiguityLabels = {
          same_connectivity_variant_group: '同连接骨架的完整结构定义差异',
          constitutional_or_positional_isomer_group: '连接/位置异构体',
          mixed_isomer_group: '混合异构体组',
          same_formula_identity_unresolved_group: '同分子式身份未充分解析',
          none: '无竞争异构体',
        };
        const resolutionLabels = {
          same_connectivity_variant_not_resolved_by_msms: '同连接骨架变体不能由当前 MS/MS 点名',
          same_connectivity_structure_definition_unresolved: '当前候选结构定义不完整，不能点名到完整立体/同位素身份',
          same_connectivity_structure_definition_specificity_ambiguity: '候选间结构定义精度不同；先补全未指定结构，再用完整身份标准品确认',
          same_connectivity_protonation_variant_requires_standardization: '先统一母体、盐型/反离子与质子化状态；这不等同于手性异构',
          mixed_isomer_group_requires_multi_axis_resolution: '混合歧义同时含连接层与完整结构定义差异，需逐轴补全并用标准品/RT/CCS验证',
          candidate_metadata_not_comparable: '候选实验条件不可比',
          conflicting_diagnostic_evidence: 'query 同时支持多个候选，证据冲突',
          conflicting_features_observed_but_unreplicated: '观察到双向候选特征，但馆藏谱重复不足，不能判定冲突方向',
          conflicting_features_observed_but_metadata_incomplete: '观察到双向候选特征，但加合物/CE/仪器信息不完整',
          diagnostic_evidence_favors_competitor: '诊断证据更支持竞争候选',
          competitor_features_observed_but_unreplicated: '观察到竞争候选特征，但馆藏谱未重复复现，不能判为更支持竞争者',
          competitor_features_observed_but_metadata_incomplete: '观察到竞争候选特征，但加合物/CE/仪器信息不完整',
          competitor_supported_pending_validation: '重复且条件可比的特征偏向竞争候选，仍待独立真值集验证',
          selected_candidate_supported_but_unreplicated: '支持当前候选，但馆藏谱未重复复现',
          selected_candidate_supported_but_metadata_incomplete: '支持当前候选，但加合物/CE/仪器信息不完整',
          selected_candidate_supported_pending_validation: '支持当前候选，仍待独立真值校准',
          single_selected_feature_insufficient: '仅一个候选特异特征，不足以点名',
          diagnostically_uninformative: '未获得可区分诊断证据',
          no_spectrally_competitive_isomer: '同分子式近邻已由生产谱图分差排除，不构成当前异构体歧义',
          same_connectivity_requires_chiral_or_orthogonal_standard: '需具立体分离能力的标准品方法或正交证据',
          resolved_by_same_method_exact_rt_standard: '由同方法、完整身份标准品 RT 解析',
        };
        const rows = topC.slice(0, 50).map((c, i) => {
          const isIso = i > 0 && f0 && (c.formula || '').trim() === f0;  // 与 Top-1 同分子式 = 异构体
          const isoMark = isIso ? '<span style="color:#b9770e;" title="与Top-1同分子式异构体">Ⓘ</span> ' : '';
          const detail = c.isomer_fragment_detail || {};
          const explicitlyNotAmbiguous = (
            i === 0 && m.isomer_ambiguous === false
          );
          const productIons = (detail.candidate_specific_product_ions || []).filter(x => x.observed);
          const neutralLosses = (detail.candidate_specific_neutral_losses || []).filter(x => x.observed);
          const opposingCount = explicitlyNotAmbiguous
            ? 0 : Number(detail.opposing_supported_features || 0);
          const backgroundOpposingCount = explicitlyNotAmbiguous
            ? Number(detail.opposing_supported_features || 0)
            : Number(detail.all_same_formula_opposing_supported_features || 0);
          const replication = detail.reference_replication || {};
          const productText = productIons.length
            ? productIons.slice(0, 4).map(x => Number(x.query_mz).toFixed(4)).join(', ')
            : '';
          const lossText = neutralLosses.length
            ? neutralLosses.slice(0, 4).map(x => `−${Number(x.query_loss_da).toFixed(4)}`).join(', ')
            : '';
          const rowResolutionStatus = explicitlyNotAmbiguous
            ? 'no_spectrally_competitive_isomer'
            : detail.resolution_status;
          const diagnosticText = Object.keys(detail).length
            ? `<div title="${Tasks.escapeHtml(detail.diagnostic_status || '')}">
                产物 ${Number(detail.unique_supported || 0)}/${Number(detail.unique_total || 0)}
                · 丢失 ${Number(detail.neutral_loss_supported || 0)}/${Number(detail.neutral_loss_total || 0)}
                · 反向支持 ${opposingCount}
                ${backgroundOpposingCount ? ` · 非竞争背景 ${backgroundOpposingCount}` : ''}
                ${productText ? `<br><span class="mono">m/z ${Tasks.escapeHtml(productText)}</span>` : ''}
                ${lossText ? `<br><span class="mono">${Tasks.escapeHtml(lossText)}</span>` : ''}
                <br><span style="color:var(--text-secondary);">${Tasks.escapeHtml(resolutionLabels[rowResolutionStatus] || rowResolutionStatus || '未判定')} · 馆藏谱 ${Number(replication.selected_reference_spectra || 0)} 张</span>
                ${detail.metadata_comparability !== 'complete'
                  ? '<br><span style="color:#b45309;">⚠ 加合物/CE/仪器元数据不完整或不可比</span>' : ''}
              </div>`
            : '—';
          const finalValue = Number(c.final_score || 0);
          const cosineValue = Number(c.cosine || 0);
          return `<tr${isIso ? ' style="background:#FEF5E7;"' : ''}>
            <td style="color:${i === 0 ? '#276749' : 'var(--text-secondary)'};font-weight:${i === 0 ? 700 : 400};">${i + 1}</td>
            <td style="max-width:220px;overflow:hidden;text-overflow:ellipsis;white-space:nowrap;" title="${Tasks.escapeHtml(c.compound_name)}">${isoMark}${Tasks.escapeHtml((c.compound_name || '').substring(0, 36))}</td>
            <td class="mono">${finalValue.toFixed(3)}</td>
            <td class="mono">${cosineValue.toFixed(3)}</td>
            <td class="mono" title="matchms direct Cosine Hungarian">${_shadowMetric(c, 'direct_cosine_hungarian')}</td>
            <td class="mono" title="matchms Modified Cosine Hungarian">${_shadowMetric(c, 'modified_cosine_hungarian')}</td>
            <td class="mono" title="matchms neutral-loss cosine">${_shadowMetric(c, 'neutral_loss_cosine')}</td>
            <td style="min-width:205px;font-size:.68rem;line-height:1.35;">${diagnosticText}</td>
          </tr>`;
        }).join('');
        const ambiguityType = m.isomer_ambiguity_type
          || (m.isomer_group_evidence || {}).ambiguity_type
          || 'same_formula_identity_unresolved_group';
        const resolutionStatus = m.isomer_resolution_status
          || (m.isomer_group_evidence || {}).resolution_status
          || 'unresolved';
        const variantSubtypeCounts = m.isomer_variant_subtype_counts
          || (m.isomer_group_evidence || {}).same_connectivity_variant_subtype_counts
          || {};
        let sameConnectivityNote = '';
        if (ambiguityType === 'mixed_isomer_group') {
          sameConnectivityNote = '该组同时含不同连接结构与同连接骨架的完整结构定义差异；必须分别补全结构定义，并用同方法标准品 RT、可比 CCS/离子迁移或经验证诊断证据逐轴排除，不能只解决其中一类歧义。';
        } else if (Number(variantSubtypeCounts.structure_definition_specificity_variant || 0) > 0) {
          sameConnectivityNote = '其中至少一条库记录只给出较少的立体/同位素结构定义；这是“结构定义精度不同”，不能简单当成两个已证实的异构体，也不能靠普通 MS/MS 给未指定结构补成立体身份。';
        } else if (Number(variantSubtypeCounts.protonation_or_charge_layer_variant || 0) > 0) {
          sameConnectivityNote = '差异位于质子化/电荷层，应先统一中性母体、盐型与离子化状态；不能误提示必须使用手性柱。';
        } else if (ambiguityType === 'same_connectivity_variant_group' || ambiguityType === 'mixed_isomer_group') {
          sameConnectivityNote = '若差异确属立体层，普通 achiral LC-MS/MS 不能证明完整立体身份，需具分离能力的标准品方法或经验证的正交证据。';
        }
        const literatureEvidence = m.isomer_literature_evidence
          || (m.isomer_group_evidence || {}).literature_evidence
          || {};
        const literatureSources = Array.isArray(literatureEvidence.sources)
          ? literatureEvidence.sources : [];
        const literatureRows = literatureSources.map(source => {
          const url = String(source.url || '');
          const label = Tasks.escapeHtml(source.title || source.source_id || '文献来源');
          const linked = url.startsWith('https://pubmed.ncbi.nlm.nih.gov/')
            ? `<a href="${Tasks.escapeHtml(url)}" target="_blank" rel="noopener noreferrer">${label}</a>`
            : label;
          const provenance = [
            source.pmid ? `PMID ${source.pmid}` : '',
            source.doi ? `DOI ${source.doi}` : '',
            source.evidence_scope || '',
          ].filter(Boolean).join(' · ');
          return `<li style="margin:.22rem 0;">${linked}
            <div style="font-size:.66rem;color:var(--text-secondary);">${Tasks.escapeHtml(provenance)}</div>
            <div style="font-size:.68rem;">${Tasks.escapeHtml(source.finding || '')}</div>
          </li>`;
        }).join('');
        const literatureContainerId = `isomerLiteratureRemote${index}`;
        const directPubMedUrl = String(literatureEvidence.pubmed_query_url || '');
        const literatureHtml = isoAmbig ? `<div style="margin-top:.45rem;padding:.45rem .55rem;border:1px solid #bfdbfe;background:#eff6ff;border-radius:4px;font-size:.7rem;line-height:1.45;">
          <b>📚 文献与正交证据</b>
          ${literatureRows ? `<ol style="margin:.25rem 0 .3rem 1rem;padding:0;">${literatureRows}</ol>` : '<div style="margin:.2rem 0;">当前仅有通用方法学边界，尚无经策展的候选对专属证据。</div>'}
          <div style="color:#92400e;">文献报道 ≠ 本系统实测馆藏谱，也不等于同方法标准品；这里的来源不改变综合匹配分、Top-1 或 High。</div>
          ${literatureEvidence.pubmed_lookup_required !== false
            ? `<div id="${literatureContainerId}" style="margin-top:.35rem;">等待 PubMed 自动补查…</div>` : ''}
        </div>` : '';
        if (isoAmbig && literatureEvidence.pubmed_lookup_required !== false) {
          isomerLiteratureRequest = {
            containerId: literatureContainerId,
            candidateNames: (
              literatureEvidence.candidate_names
              || [m.compound_name, ...(m.isomer_set || []).map(item => item.compound_name)]
            ).slice(0, 8),
            formula: literatureEvidence.formula || m.isomer_group_formula || f0 || '',
            searchUrl: directPubMedUrl,
          };
        }
        const isoBanner = isoAmbig
          ? `<div style="padding:0.4rem 0.6rem;background:#FEF5E7;border:1px solid #F5B041;border-radius:4px;font-size:0.72rem;color:#b9770e;margin-bottom:0.3rem;">⚠ <b>${Tasks.escapeHtml(ambiguityLabels[ambiguityType] || ambiguityType)}</b>：具体结构 ${Tasks.escapeHtml(m.exact_structure_confidence || m.confidence || 'Unknown')}；${Tasks.escapeHtml(m.isomer_group_formula || f0 || '同分子式')} 异构体组 ${Tasks.escapeHtml(m.isomer_group_confidence || 'Unknown')}。当前状态：${Tasks.escapeHtml(resolutionLabels[resolutionStatus] || resolutionStatus)}。产物离子与中性丢失按 0.02 Da、馆藏重复复现、竞争候选反向证据和实验条件共同审计；未经独立真值校准不改变 Top-1。${sameConnectivityNote ? ` ${Tasks.escapeHtml(sameConnectivityNote)}` : ''}</div>`
          : '';
        topCandHtml = `<details class="pw-card"${isoAmbig ? ' open' : ''}><summary>🏆 Top-K 候选 (${topC.length}；第 1 名为本针排序先导)${isoAmbig ? ' ⚠' : ''}</summary>
          <div style="padding:0.4rem 0.6rem;font-size:0.78rem;">${isoBanner}${literatureHtml}
          <div style="margin-bottom:.35rem;color:var(--text-secondary);font-size:.68rem;">Direct、Mod-H 与 NL 为候选选定后的 matchms 只读对照；均不参与综合匹配分、排序或置信度。</div>
          <table style="width:100%;border-collapse:collapse;">
          <thead><tr style="color:var(--text-muted);font-size:0.7rem;"><th>#</th><th style="text-align:left;">化合物</th><th>综合匹配分</th><th>MS/MS 相似度</th><th>Direct</th><th>Mod-H</th><th>NL</th><th style="text-align:left;">候选专属诊断证据</th></tr></thead>
          <tbody>${rows}</tbody></table></div></details>`;
      }
      if (ed && Object.keys(ed).length) {
        const dimSpecs = [
          ['neutral_loss', 'neutral_loss_score', '中性丢失'],
          ['rule_coverage', 'rule_coverage_score', '规则覆盖'],
          ['fg_signature', 'fg_signature_score', '官能团特征'],
          ['feature_library', 'feature_library_score', '大类特征库'],
          ['structure_fragment', 'structure_fragment_score', '结构-碎片'],
          ['diagnostic_ion', 'diagnostic_ion_score', '诊断离子'],
          ['bond_energy', 'bond_energy_score', '键能支持'],
          ['isomer_differentiation', 'isomer_differentiation_score', '异构区分'],
        ];
        const assessments = ed.dimension_assessments || {};
        const dims = dimSpecs.map(([key, scoreKey, fallbackLabel]) => {
          const assessment = assessments[key] || {};
          const rawValue = assessment.score != null ? assessment.score : (ed[scoreKey] || 0);
          return {
            key,
            name: assessment.label || fallbackLabel,
            value: Number(rawValue) || 0,
            status: assessment.status || ((Number(rawValue) || 0) > 0 ? 'scored' : 'no_evidence'),
            reason: assessment.reason || '',
          };
        });
        const reasonLabels = {
          classification_unconfirmed: '尚无可靠类别归属，未启用按大类查找的规则',
          no_supported_rule_family_for_class: '该分类暂不在运行时规则族覆盖范围内',
          all_rules_outside_polarity_adduct_structure_or_energy_scope: '规则与当前极性、加合物、结构或碰撞能适用域不符',
          no_executable_rules_for_class: '该类当前没有可执行规则',
          no_functional_group_signature_for_class: '该类暂无官能团签名档案',
          no_class_feature_profile: '该类暂无大类特征档案',
          no_scope_applicable_class_features: '当前极性/结构下没有适用的大类特征',
          no_isomer_rules_for_class: '该类暂无可用的异构诊断离子对',
          feature_library_unavailable: '特征库不可用',
          class_signature_not_observed: '规则适用，但本谱未观察到类别签名',
          no_curated_neutral_loss_match: '规则适用，但本谱未命中策展中性丢失',
          no_scope_guarded_diagnostic_ion_match: '规则适用，但本谱未命中范围受控诊断离子',
          executable_rules_present_but_not_covered: '存在适用规则，但观测碎片未覆盖',
          class_features_present_but_not_observed: '存在适用大类特征，但本谱未命中',
          isomer_rules_present_but_no_diagnostic_pair: '存在异构规则，但未观察到诊断离子对',
          no_predicted_single_bond_fragment_match: '候选结构可分析，但未命中当前单键断裂模型',
        };
        const fusion = ed.score_fusion || {};
        const fusionComponents = fusion.components || m.expert_score_components || {};
        const selectedSource = fusion.selected_source || m.expert_score_source || 'none';
        const sourceOrder = ['base_8d', 'v3_kb', 'compound_pathway'];
        const sourceBadges = sourceOrder
          .filter(key => fusionComponents[key] && fusionComponents[key].available)
          .map(key => {
            const item = fusionComponents[key];
            const selected = key === selectedSource;
            return `<span style="display:inline-flex;gap:4px;align-items:center;padding:3px 7px;border-radius:999px;border:1px solid ${selected ? '#805AD5' : '#B8C5CA'};background:${selected ? '#FAF5FF' : '#F7FAFC'};color:${selected ? '#553C9A' : '#43585F'};font-size:.68rem;">
              ${Tasks.escapeHtml(item.label || key)}${selected ? ' · 记录来源' : ''}
            </span>`;
          }).join(' ');
        const classProvenance = ed.classification_provenance || {};
        const classReliable = m.classification_reliable != null
          ? !!m.classification_reliable
          : classProvenance.reliable != null
            ? !!classProvenance.reliable : !!ed.classification_reliable;
        const className = m.display_class || m.classification
          || classProvenance.analysis_class || ed.analysis_class || 'Unknown';
        const classSource = m.display_classification_source || m.classification_source
          || classProvenance.input_source || 'unknown';
        const classConfidence = classProvenance.input_confidence;
        const classAuditHtml = `<div style="margin:6px 0 8px;padding:6px 8px;border-radius:5px;background:${classReliable ? '#F0FFF4' : '#FFFAF0'};color:${classReliable ? '#276749' : '#805014'};font-size:.7rem;line-height:1.45;">
          分类门控：<b>${Tasks.escapeHtml(className)}</b> · 来源 ${Tasks.escapeHtml(classSource)}
          ${classConfidence != null ? ` · 置信 ${Number(classConfidence).toFixed(2)}` : ''}
          · ${classReliable ? '结构规则已形成可靠类别归属，可用大类规则' : '尚无可靠类别归属；大类规则显示 N/A，避免跨类误报'}
        </div>`;
        expertHtml = `
          <div style="margin-top:0.75rem;padding:0.65rem;background:#FAFBFC;border:1px solid var(--border);border-radius:6px;">
            <div style="font-size:0.78rem;font-weight:600;color:#553C9A;margin-bottom:0.3rem;">🧪 结构规则逐维证据（已评价 ${Number(expertCoverage.scored || 0)}/${Number(expertCoverage.total || 8)} 维）</div>
            <div style="display:flex;flex-wrap:wrap;gap:4px;margin-bottom:4px;">${sourceBadges}</div>
            <div style="font-size:.7rem;color:#4F636A;margin-bottom:6px;">下方逐维报告“有支持证据／已评价但未观察／不适用／未建立”；这才是可解释结论。旧式聚合值不是概率，也不参与身份排序。</div>
            <details style="margin:0 0 6px;font-size:.68rem;color:#64748B;"><summary style="cursor:pointer;">展开旧版数值审计字段</summary><div style="padding:4px 0 0 8px;">未校准规则聚合 ${expertScore.toFixed(3)}；记录来源 ${Tasks.escapeHtml(selectedSource)}。该数值不得用于身份概率、证据等级或候选排序。</div></details>
            ${classAuditHtml}
            ${dims.map(dim => {
              const pct = Math.round(dim.value * 100);
              const isNA = dim.status === 'not_applicable' || dim.status === 'unavailable';
              const color = isNA ? '#7A8990' : (
                dim.status === 'no_evidence' ? '#718086' : (
                  dim.value >= 0.6 ? '#38A169' : dim.value >= 0.3 ? '#D69E2E' : '#3182CE'
                )
              );
              const statusText = dim.status === 'not_applicable' ? 'N/A'
                : dim.status === 'unavailable' ? '无档案'
                : dim.status === 'no_evidence' ? '未观察'
                : dim.value.toFixed(2);
              const reason = reasonLabels[dim.reason] || dim.reason || '';
              return `<div style="display:flex;align-items:center;gap:6px;margin:4px 0;">
                <span style="font-size:.68rem;color:#52666C;width:92px;text-align:right;" title="${Tasks.escapeHtml(reason)}">${Tasks.escapeHtml(dim.name)}</span>
                <div style="flex:1;background:#E7ECEE;border-radius:3px;height:8px;">
                  <div style="width:${pct}%;background:${color};height:100%;border-radius:3px;"></div>
                </div>
                <span style="font-size:.68rem;font-weight:700;width:48px;color:${color};" title="${Tasks.escapeHtml(reason)}">${statusText}</span>
              </div>`;
            }).join('')}
            ${ed.chemical_reasoning && !(classReliable && classProvenance.reliable === false) ? `<div style="margin-top:7px;padding:7px 9px;background:#FFFBEB;border:1px solid #E8C95A;border-radius:4px;font-size:.72rem;color:#694116;max-height:110px;overflow:auto;line-height:1.5;">${Tasks.escapeHtml(ed.chemical_reasoning)}</div>` : (classReliable && classProvenance.reliable === false ? '<div style="margin-top:7px;padding:7px 9px;background:#eff6ff;border:1px solid #93c5fd;border-radius:4px;font-size:.72rem;color:#1e40af;">当前结构分类已更新；旧任务中基于“类别证据不足”的压缩推理文字已停止展示。逐维原始记录仍保留供审计，重新分析可刷新规则评价。</div>' : '')}
          </div>`;
      }

      // ── 碎裂途径面板 (级联/步骤/诊断离子/中性丢失) ──
      // 数据已在 expert_details 中 (后端零改动), 此前仅压缩 reasoning 呈现
      let pathwayHtml = '';
      if (ed && Object.keys(ed).length) {
        const nl = ed.neutral_loss_analysis || {};
        const currentFragmentation = (m.auxiliary_interpretation || {}).fragmentation || {};
        const fmz = v => (v == null || isNaN(v)) ? '—' : Number(v).toFixed(2);
        const esc = s => Tasks.escapeHtml(String(s == null ? '—' : s));
        const sections = [];
        const dbPath = _pathwayStatus(m, _currentIonMode);

        // 数据库预存通路必须与本次离子极性一致；不一致时明确提示，避免把正离子
        // 碎裂特征套到负离子样本上。样本自身的中性丢失/碎裂步骤仍在下方独立展示。
        if (dbPath.state === 'mismatch') {
          sections.push('<div style="padding:0.45rem 0.6rem;background:#fffbeb;border:1px solid #fcd34d;border-radius:5px;color:#92400e;font-size:0.73rem;">⚠ 馆藏预存碎裂路径仅包含另一离子极性，本次未采用；下方仍展示由当前负离子谱图现场计算的中性丢失与结构推理。</div>');
        } else if (currentFragmentation.candidate_specific_library_pathway?.status === 'available' || dbPath.data) {
          const currentLibraryPath = currentFragmentation.candidate_specific_library_pathway || {};
          const dbIons = (currentLibraryPath.recurring_ions || dbPath.data?.diagnostic_ions || []).slice(0, 8);
          const dbLosses = (currentLibraryPath.recurring_losses || dbPath.data?.neutral_losses || []).slice(0, 8);
          if (dbIons.length || dbLosses.length) {
            const ionText = dbIons.map(x => `${fmz(x.mz)}${x.annotation ? ' (' + esc(x.annotation) + ')' : ''}`).join(' · ');
            const lossText = dbLosses.map(x => `−${fmz(x.loss_da)}${x.annotation ? ` (${esc(x.annotation)}${x.annotation_status === 'legacy_annotation_rejected_by_current_policy' ? '；旧归属已驳回' : ''})` : '（仅质量差）'}`).join(' · ');
            sections.push(`<details class="pw-det"><summary>📚 候选相关参考谱复现特征（${Number(currentLibraryPath.n_reference_spectra || dbPath.data?.n_spectra || 0)} 张库谱）</summary><div style="padding:0.45rem 0.55rem;font-size:0.73rem;line-height:1.5;">${ionText ? `<div><b>复现碎片离子:</b> ${ionText}</div>` : ''}${lossText ? `<div><b>前体–碎片质量差:</b> ${lossText}</div>` : ''}<div style="color:var(--text-secondary);margin-top:4px;">候选依赖的参考谱统计，不作为独立身份分数，也不等于已确认机理。</div></div></details>`);
          }
        }

        // 糖苷糖基级联
        const sugar = nl.sugar_loss_cascade;
        if (sugar && sugar.cascade && sugar.cascade.length) {
          const chain = sugar.cascade.map(s =>
            `<span class="pw-node">${fmz(s.from_mz)}</span>` +
            `<span class="pw-arrow">−${fmz(s.mass_difference)} (${esc(s.sugar_name_cn)})→</span>`
          ).join('') + `<span class="pw-node pw-agly">${fmz(sugar.aglycone_mz)}<span class="pw-tag">苷元</span></span>`;
          const comp = sugar.sugar_composition ? Object.entries(sugar.sugar_composition).map(([n, c]) => `${c} ${n}`).join(' + ') : '';
          sections.push(
            `<div class="pw-cascade pw-sugar"><b>🍬 糖基顺序丢失级联</b> <span class="pw-sub">糖苷/皂苷</span>` +
            `<div class="pw-chain">${chain}</div>` +
            `<div class="pw-meta">糖组成: <b>${esc(comp)}</b> · 共 ${sugar.n_sugars} 个糖 · 苷元 m/z ${fmz(sugar.aglycone_mz)}</div></div>`
          );
        }

        // 乌头碱酯基级联
        const ester = nl.alkaloid_ester_cascade;
        if (ester && ester.cascade && ester.cascade.length) {
          const chain = ester.cascade.map(s =>
            `<span class="pw-node">${fmz(s.from_mz)}</span>` +
            `<span class="pw-arrow">−${fmz(s.mass_difference)} (${esc(s.loss_name_cn)})→</span>`
          ).join('') + `<span class="pw-node pw-agly">${fmz(ester.aglycone_mz)}<span class="pw-tag">乌头原碱</span></span>`;
          const comp = ester.loss_composition ? Object.entries(ester.loss_composition).map(([n, c]) => `${c}×${n}`).join(' + ') : '';
          sections.push(
            `<div class="pw-cascade pw-ester"><b>⚗️ 乌头碱型酯基丢失级联</b> <span class="pw-subtype">${esc(ester.aconitine_subtype)}</span>` +
            `<div class="pw-chain">${chain}</div>` +
            `<div class="pw-meta">丢失组成: <b>${esc(comp)}</b> · 共 ${ester.n_losses} 步</div></div>`
          );
        }

        // 碎裂步骤表
        const best = currentFragmentation.steps || [];
        if (best.length) {
          const rows = best.slice(0, 6).map(s =>
            `<tr><td>${fmz(s.from_mz)}</td><td>${fmz(s.to_mz)}</td><td>${fmz(s.mass_difference)}</td>` +
            `<td>${esc(s.neutral_loss)}</td><td>${esc(s.mechanism)}</td><td>${esc(s.priority)}</td></tr>`
          ).join('');
          sections.push(
            `<details class="pw-det" open><summary>🧬 碎裂步骤 <span class="pw-sub">fragment→fragment</span></summary>` +
            `<table class="pw-tbl"><thead><tr><th>From</th><th>To</th><th>Δm</th><th>中性丢失</th><th>机理</th><th>优先级</th></tr></thead><tbody>${rows}</tbody></table></details>`
          );
        }

        // 诊断离子表
        const di = ed.diagnostic_ion_analysis || {};
        const ions = (di.matched_ions || []).slice(0, 6);
        if (ions.length) {
          const rows = ions.map(io => {
            const sp = io.intensity_characteristic || '';
            const cls = sp === 'high' ? 'pw-sp-high' : (sp === 'medium' ? 'pw-sp-med' : '');
            return `<tr><td>${fmz(io.observed_mz)}</td><td>${esc(io.ion_name)}</td><td>${esc(io.structural_origin)}</td><td><span class="pw-sp ${cls}">${esc(sp)}</span></td></tr>`;
          }).join('');
          sections.push(
            `<details class="pw-det"><summary>🔬 诊断离子 (${ions.length})</summary>` +
            `<table class="pw-tbl"><thead><tr><th>m/z</th><th>名称</th><th>结构来源</th><th>特异性</th></tr></thead><tbody>${rows}</tbody></table></details>`
          );
        }

        // 中性丢失表
        const nls = currentFragmentation.neutral_loss_assignments || nl.neutral_losses || [];
        if (nls.length) {
          const rows = [];
          const explained = nls.filter(it => it && it.neutral_loss);
          const unassigned = nls.filter(it => it && !it.neutral_loss);
          for (const it of explained.slice(0, 8)) {
            const tentative = it.assignment_level === 'tentative_composite'
              ? ' <span class="pw-sub">（累积候选）</span>' : '';
            rows.push(`<tr><td>${fmz(it.mass_difference)}</td><td class="pw-expl">${esc(it.neutral_loss)}${tentative}</td><td>${esc(it.formula)}</td><td>${esc(it.mechanism)}</td></tr>`);
          }
          if (rows.length) {
            sections.push(
              `<details class="pw-det"><summary>📉 已解释中性丢失 (${explained.length}) <span class="pw-sub">precursor→fragment</span></summary>` +
              `<table class="pw-tbl"><thead><tr><th>Δm</th><th>中性丢失</th><th>分子式</th><th>经验规则类型</th></tr></thead><tbody>${rows.join('')}</tbody></table></details>`
            );
          }
          if (unassigned.length) {
            const values = unassigned.slice(0, 8).map(it => fmz(it.mass_difference)).join('、');
            sections.push(
              `<details class="pw-det pw-unassigned"><summary>🧩 未归属前体–碎片差值 (${unassigned.length})</summary>` +
              `<div class="pw-unexpl">这些是尚未匹配可靠化学规则的碎片质量差，不等同于已确认中性丢失。` +
              `${values ? ` 代表性 Δm：${esc(values)}。` : ''}</div></details>`
            );
          }
        }
        const rejectedAssignments = currentFragmentation.rejected_legacy_assignments || [];
        if (rejectedAssignments.length) {
          const rejectedRows = rejectedAssignments.slice(0, 12).map(it => `<tr><td>${fmz(it.mass_difference)}</td><td>${esc(it.neutral_loss || it.neutral_loss_en)}</td><td>${esc(it.formula)}</td><td>${esc(it.reason || 'current_policy_rejected')}</td></tr>`).join('');
          sections.push(`<details class="pw-det pw-unassigned"><summary>⚠ 旧版归属的当前策略复核（已驳回 ${rejectedAssignments.length} 条）</summary><table class="pw-tbl"><thead><tr><th>Δm</th><th>旧归属</th><th>分子式</th><th>驳回原因</th></tr></thead><tbody>${rejectedRows}</tbody></table></details>`);
        }

        if (sections.length) {
          const _structInPathway = smiles ? `<div id="cfmCascade${index}" style="margin-bottom:0.5rem;padding:0.4rem;background:#FAFBFC;border-radius:4px;">
            <div style="font-size:0.72rem;color:var(--text-secondary);margin-bottom:3px;font-weight:600;">🧬 候选产物结构解释（观测峰支持优先，模型预测兜底）</div>
            <div class="cfm-loading" style="padding:0.5rem;text-align:center;color:#999;font-size:0.7rem;">⏳ 正在读取候选碎片结构证据…</div></div>` : '';
          pathwayHtml = `<details class="pw-card" open><summary>🧬 经验碎裂证据与中性丢失序列 (Empirical Fragment Evidence)</summary><div style="padding:0.5rem 0.6rem;">${_structInPathway}${sections.join('')}</div></details>`;
        }
      }

	      const governancePassportHtml = _evidencePassportHtml(m);
	      const precursorHighGateHtml = m.precursor_high_evidence_gate_required
	        ? `<div style="margin-top:0.55rem;padding:0.55rem 0.7rem;border:1px solid #f0b44c;border-left:4px solid #c47b00;background:#fff8e8;border-radius:6px;font-size:0.75rem;color:#7a4b00;line-height:1.55;">
	          <b>High 证据门控：该候选已降为 Medium。</b><br>${Tasks.escapeHtml(String(m.precursor_high_evidence_gate_reason_zh || m.precursor_high_evidence_gate_reason || '碎片配对稀疏且缺少独立前体支持'))}
	        </div>` : '';
	      const rtEv = m.rt_evidence || null;
	      const rtEvidenceHtml = rtEv ? `<div style="margin-top:0.55rem;padding:0.5rem 0.65rem;border:1px solid ${rtEv.eligible_for_exact_identity ? '#bbf7d0' : '#fed7aa'};background:${rtEv.eligible_for_exact_identity ? '#f0fdf4' : '#fff7ed'};border-radius:6px;font-size:0.75rem;color:${rtEv.eligible_for_exact_identity ? '#166534' : '#9a3412'};">
	        <b>RT evidence:</b> ${Tasks.escapeHtml(rtEv.status || '—')} · ΔRT ${rtEv.rt_delta_sec != null ? rtEv.rt_delta_sec + 's' : '—'} · tol ${rtEv.tolerance_sec != null ? rtEv.tolerance_sec + 's' : '—'}<br>
	        <span style="color:var(--text-secondary);">Standard: ${Tasks.escapeHtml(rtEv.standard_name || '—')} · identity ${Tasks.escapeHtml(rtEv.identity_scope || 'unverified')} · method ${Tasks.escapeHtml(rtEv.method_status || 'unverified')} · ${rtEv.eligible_for_exact_identity ? '可用于精确身份' : '不可用于精确身份升档'}</span>
	      </div>` : '';
	      const sc = m.sample_consensus_support || null;
		      const sampleSupportHtml = sc ? `<div style="margin-top:0.55rem;padding:0.5rem 0.65rem;border:1px solid #bfdbfe;background:#eff6ff;border-radius:6px;font-size:0.75rem;color:#1e40af;">
		        <b>Sample consensus:</b> ${sc.n_spectra || 0} MS2 spectra · RT span ${sc.rt_span_min != null ? sc.rt_span_min + ' min' : '—'} · best cosine ${sc.best_cosine != null ? sc.best_cosine : '—'} · best final ${sc.best_final != null ? sc.best_final : '—'}${sc.n_collision_energies ? ` · ${sc.n_collision_energies} collision energies (${(sc.collision_energies || []).join('/') } eV)` : ''}
		      </div>` : '';
		      const iso = m.formula_isotope_evidence || null;
		      const isoStatus = r.isotope_evidence_status || 'not_available';
		      const isotopeEvidenceHtml = iso ? `<div style="margin-top:0.55rem;padding:0.5rem 0.65rem;border:1px solid ${iso.status === 'supporting' ? '#bbf7d0' : iso.status === 'consistent' ? '#bfdbfe' : '#fed7aa'};background:${iso.status === 'supporting' ? '#f0fdf4' : iso.status === 'consistent' ? '#eff6ff' : '#fff7ed'};border-radius:6px;font-size:0.75rem;line-height:1.5;">
		        <b>MS1 isotope/formula evidence:</b> ${Tasks.escapeHtml(iso.status)} · score ${Number(iso.score || 0).toFixed(3)}<br><span style="color:var(--text-secondary);">${Tasks.escapeHtml(iso.interpretation || '')}</span>
		      </div>` : `<div style="margin-top:0.55rem;padding:0.45rem 0.65rem;border:1px solid var(--border);background:#f8fafc;border-radius:6px;font-size:0.72rem;color:var(--text-secondary);">
		        <b>MS1 isotope evidence:</b> ${isoStatus === 'referenced_ms1_not_in_file' ? '输入文件仅含 MS2，引用的 MS1 未包含，不能计算同位素包络' : Tasks.escapeHtml(isoStatus)}
		      </div>`;
		      const contaminantHtml = m.potential_contaminant ? `<div style="margin-top:0.55rem;padding:0.55rem 0.7rem;border:1px solid #fdba74;background:#fff7ed;border-radius:6px;font-size:0.76rem;color:#9a3412;line-height:1.5;">
		        <b>⚠ 潜在外源污染物（不计入天然产物）</b><br>${Tasks.escapeHtml(m.contaminant_note || '请结合空白样、前体 XIC 与保留时间复核。')}
		      </div>` : '';
	      const deconv = m.chimeric_deconvolution || null;
	      const deconvOk = deconv && deconv.status === 'accepted_score_applied';
	      const deconvEvidenceHtml = deconv ? `<div style="margin-top:0.55rem;padding:0.55rem 0.7rem;border:1px solid ${deconvOk ? '#86efac' : '#cbd5e1'};background:${deconvOk ? '#f0fdf4' : '#f8fafc'};border-radius:6px;font-size:0.75rem;color:${deconvOk ? '#166534' : '#475569'};line-height:1.55;">
	        <b>嵌合 MS/MS 稀疏去卷积:</b> ${Tasks.escapeHtml(deconv.status || 'not evaluated')} · isolation ${Tasks.escapeHtml(deconv.isolation_status || 'unavailable')} · purity ${deconv.isolation_purity != null ? Number(deconv.isolation_purity).toFixed(3) : '—'}<br>
		        <span>重建 cosine ${deconv.reconstruction_cosine != null ? Number(deconv.reconstruction_cosine).toFixed(3) : '—'} · target fraction ${deconv.target_fraction != null ? Number(deconv.target_fraction).toFixed(3) : '—'} · active components ${deconv.active_component_count ?? '—'}${deconvOk ? ` · raw→deconvolved cosine ${Number(m.raw_cosine || 0).toFixed(3)}→${Number(m.cosine || 0).toFixed(3)} · final ${Number(m.raw_final_score || 0).toFixed(3)}→${Number(m.final_score || 0).toFixed(3)}` : ''} · raw candidate cosine margin ${m.candidate_margin != null ? Number(m.candidate_margin).toFixed(3) : '—'}</span><br>
		        <span style="color:var(--text-secondary);">非负 L1 系数仅用于谱图分解，不解释为定量丰度；Top-1、Top-K、异构体竞争和候选 margin 均锁定原始实测谱。</span>
	      </div>` : '';
	      const iin = r.ion_identity_evidence || m.ion_identity_evidence || null;
	      const iinConflict = !!(m.ion_identity_conflict || m.ion_identity_adduct_conflict);
	      const ionIdentityHtml = iin ? `<div style="margin-top:0.55rem;padding:0.55rem 0.7rem;border:1px solid ${iinConflict ? '#fca5a5' : '#93c5fd'};background:${iinConflict ? '#fef2f2' : '#eff6ff'};border-radius:6px;font-size:0.75rem;color:${iinConflict ? '#991b1b' : '#1e40af'};line-height:1.55;">
	        <b>Ion Identity Molecular Network:</b> ${Tasks.escapeHtml(iin.network_id || '—')} · ${Tasks.escapeHtml(iin.ion_identity || r.ion_identity || '—')} · neutral mass ${iin.neutral_mass != null ? Number(iin.neutral_mass).toFixed(5) : '—'}<br>
	        <span>峰形 Pearson r ${iin.best_shape_correlation != null ? Number(iin.best_shape_correlation).toFixed(3) : '—'} · overlap ${iin.best_overlap_fraction != null ? Number(iin.best_overlap_fraction).toFixed(3) : '—'} · representative ${Tasks.escapeHtml(iin.representative_feature_id || '—')}${iin.adduct_consistency ? ` · adduct check ${Tasks.escapeHtml(iin.adduct_consistency)}` : ''}</span>${iinConflict ? '<br><b>⚠ 网络身份或加合物冲突：High 已被限制。</b>' : ''}
	      </div>` : '';

	      const html = `
      <div class="modal-overlay" onclick="if(event.target===this)App.closeModal()">
      <div class="modal modal-lg">
        <div class="flex-between" style="margin-bottom:0.65rem;">
          <div>
            <h3 style="font-size:1.05rem;">${Tasks.escapeHtml((m.compound_name || '').substring(0, 80))}
              ${_nSpecThis > 1 ? `<span style="font-size:0.72rem;font-weight:600;color:#166534;background:#dcfce7;padding:2px 8px;border-radius:9px;margin-left:6px;vertical-align:middle;" title="该化合物共命中 ${_nSpecThis} 张谱图(主表已合并)">命中 ${_nSpecThis} 谱</span>` : ''}
            </h3>
            <div style="font-size:0.8rem;color:var(--text-secondary);">
              RT ${Tasks.escapeHtml(rtText)} · MS2 ${Tasks.escapeHtml(spectrumId)} · ${Tasks.escapeHtml(_precursorMassInline(r, 5))} · ${Tasks.escapeHtml(String(m.formula || '—'))} · ${Tasks.escapeHtml(classLabel)}
            </div>
          </div>
          ${headerConfidenceHtml}
        </div>

        <section class="result-evidence-hero">
          <div class="result-identity-strip ${exactRtStandard ? 'scope-standard' : isIsomerGroupScope ? 'scope-isomer' : 'scope-candidate'}">
            <span class="result-identity-kicker">身份范围</span>
            <strong>${Tasks.escapeHtml(detailIdentityLabel)}</strong>
            <span>${Tasks.escapeHtml(detailIdentityNote)}</span>
          </div>
          <div class="result-visual-grid">
            <section class="result-structure-panel">
              <div class="result-visual-heading"><span>候选二维结构</span><small>${Tasks.escapeHtml(String(m.formula || '—'))}</small></div>
              <div class="result-structure-stage">${structHtml}</div>
              <div class="result-structure-caption"><b>${Tasks.escapeHtml((m.compound_name || '未命名候选').substring(0, 80))}</b><span class="mono">${Tasks.escapeHtml(String(m.inchikey?.substring(0, 27) || '—'))}</span></div>
            </section>
            <section class="result-spectrum-panel">
              <div class="result-visual-heading"><span>${hasDetailMirrorEvidence ? '实测谱 / 馆藏谱镜像对照' : '实测 MS/MS 棒谱'}</span><small>${detailQueryPeaks.length} 个查询峰${hasDetailMirrorEvidence ? ` · ${detailLibraryPeaks.length} 个馆藏峰` : ''}</small></div>
              <div id="detailPrimarySpectrum" class="result-spectrum-stage" aria-live="polite"><div class="result-visual-empty">正在准备谱图…</div></div>
            </section>
          </div>
          <div class="result-metric-ribbon">
            <div><span>MS/MS 相似度</span><b>${Number(m.cosine || 0).toFixed(4)}</b><small>无量纲相似度</small></div>
            <div><span>综合匹配分</span><b>${Number(m.final_score || 0).toFixed(4)}</b><small>排序分，不是概率</small></div>
            <div><span>质量误差</span><b>${m.mass_ppm != null ? Number(m.mass_ppm).toFixed(2) : '—'}</b><small>ppm</small></div>
            <div><span>保留时间</span><b>${rtMinutes == null ? '—' : rtMinutes.toFixed(3)}</b><small>min</small></div>
            <div><span>前体</span><b>${Number(r.precursor_mz || 0).toFixed(5)}</b><small>m/z</small></div>
            <div><span>实测碎片</span><b>${detailQueryPeaks.length}</b><small>peaks</small></div>
          </div>
        </section>

        <details class="result-evidence-details">
          <summary><span>详细证据</span><small>方法、原始元数据、完整专家解释与候选表</small></summary>
          <div class="result-evidence-details-body">
        <div style="display:flex;gap:1.25rem;flex-wrap:wrap;">
          <div style="flex:1;min-width:280px;">
            <div style="margin-bottom:.65rem;padding:.55rem .7rem;border:1px solid #9cc9c6;border-left:4px solid #176b68;background:#f0f8f7;border-radius:6px;font-size:.76rem;line-height:1.55;color:#214c50;"><b>原始 MS2 定位与检索质量轴</b><br>源文件 ${Tasks.escapeHtml(_currentTaskFilename || '未记录')} · RT ${Tasks.escapeHtml(rtText)} · ${Tasks.escapeHtml(_precursorMassInline(r, 5))} · MS2 谱图编号/native ID ${Tasks.escapeHtml(spectrumId)}${r.precursor_scan_ref ? ` · 关联前级扫描 ${Tasks.escapeHtml(String(r.precursor_scan_ref))}` : ''}<br><span style="color:var(--text-secondary);">优先按谱图编号/native ID 回查；无法检索编号时，用 RT 与隔离目标定位。只有通过 RT/电荷/唯一性检查的连续 MS1 特征质心才用于 HRMS 检索。</span></div>
            <table style="font-size:0.8rem;width:100%;">
              <tr><td style="color:var(--text-secondary);padding:3px 6px;">SMILES</td>
                  <td class="mono" style="word-break:break-all;font-size:0.72rem;">${Tasks.escapeHtml(smiles.substring(0, 150)) || '—'}</td></tr>
              <tr><td style="color:var(--text-secondary);padding:3px 6px;">InChIKey</td>
                  <td class="mono" style="font-size:0.72rem;">${Tasks.escapeHtml(String(m.inchikey?.substring(0, 27) || '—'))}</td></tr>
              <tr><td style="color:var(--text-secondary);padding:3px 6px;">样品离子形式</td>
                  <td>${Tasks.escapeHtml(String(adductEvidence.sample_ion_form || '未可靠判定'))}<span style="font-size:.68rem;color:var(--text-secondary);"> · ${Tasks.escapeHtml(String(adductEvidence.sample_ion_form_source || 'unresolved'))}${adductEvidence.candidate_independent ? ' · 候选无关' : ' · 候选依赖'}</span></td></tr>
              <tr><td style="color:var(--text-secondary);padding:3px 6px;">谱库加合物</td>
                  <td>${Tasks.escapeHtml(String(adductEvidence.library_adduct || m.adduct || '—'))}<span style="font-size:.68rem;color:var(--text-secondary);"> · 与样品 ${Tasks.escapeHtml(String(adductEvidence.sample_library_consistency || 'not_tested'))}</span></td></tr>
              ${(adductEvidence.candidate_ion_form_hypothesis || m.candidate_adduct_hypothesis) ? `<tr><td style="color:var(--text-secondary);padding:3px 6px;">候选离子假设</td>
                  <td>${Tasks.escapeHtml(String(adductEvidence.candidate_ion_form_hypothesis || m.candidate_adduct_hypothesis))}<span style="font-size:.68rem;color:#1d4ed8;"> · 由候选结构/分子式与谱库前体反推，非样品独立证据</span></td></tr>` : ''}
	              <tr><td style="color:var(--text-secondary);padding:3px 6px;">前体电荷 / 碰撞能量</td>
	                  <td>${r.precursor_charge != null ? Math.abs(Number(r.precursor_charge)) + (String(r.ion_mode || '').toUpperCase() === 'NEG' ? '−' : '+') : '—'} / ${r.collision_energy != null ? Number(r.collision_energy).toFixed(1) + (r.collision_energy_unit ? ' ' + Tasks.escapeHtml(String(r.collision_energy_unit)) : ' (unit unknown)') : '—'}</td></tr>
	              <tr><td style="color:var(--text-secondary);padding:3px 6px;">类别 (Class)</td>
	                  <td>${Tasks.escapeHtml(classLabel)}<span style="font-size:0.7rem;color:var(--text-secondary);"> · 来源: ${Tasks.escapeHtml(m.dataset_source || '—')}</span></td></tr>
	              <tr><td style="color:var(--text-secondary);padding:3px 6px;">谱库证据层级</td>
	                  <td>${Tasks.escapeHtml(m.source_tier || 'unknown')}</td></tr>
              ${m.source_compound_name ? `<tr><td style="color:var(--text-secondary);padding:3px 6px;">原始库名称</td><td class="mono" style="word-break:break-all;font-size:0.68rem;color:var(--text-secondary);">${Tasks.escapeHtml(m.source_compound_name)}</td></tr>` : ''}
              <tr><td style="color:var(--text-secondary);padding:3px 6px;">MS2 谱图编号/native ID</td>
                  <td class="mono" style="font-size:0.72rem;">${Tasks.escapeHtml(spectrumId)}</td></tr>
              <tr><td style="color:var(--text-secondary);padding:3px 6px;">样本实测保留时间</td>
                  <td class="mono" style="font-size:0.72rem;">${Tasks.escapeHtml(rtText)}</td></tr>
              <tr><td style="color:var(--text-secondary);padding:3px 6px;">检索质量轴</td>
                  <td class="mono" style="font-size:0.72rem;">${Tasks.escapeHtml(_precursorMassInline(r, 7))}</td></tr>
            </table>

            <div class="score-row" style="margin-top:0.65rem;">
              ${[
                ['MS/MS 相似度', m.cosine], ['综合匹配分', m.final_score], ['质量匹配分', m.mass_score || 0],
                ['碎片匹配分', m.fragment_score || 0], ['强度匹配分', m.intensity_score || 0],
              ].map(([l, v]) => `<div class="score-chip"><div class="score-val">${v.toFixed(4)}</div><div class="score-lbl">${l}</div></div>`).join('')}
	              <div class="score-chip"><div class="score-val">${m.mass_ppm != null ? m.mass_ppm.toFixed(2) : '—'}</div><div class="score-lbl">质量误差 ppm</div></div>
	            </div>
	            ${precursorHighGateHtml}
		            ${contaminantHtml}
		            ${rtEvidenceHtml}
	            ${sampleSupportHtml}
	            ${isotopeEvidenceHtml}
	            ${deconvEvidenceHtml}
	            ${ionIdentityHtml}
	            ${expertHtml}
            ${evographHtml}
            ${pathwayHtml}
            ${v3ReasoningHtml}
            ${isopairHtml}
            ${topCandHtml}
            ${governancePassportHtml}
          </div>
          <div style="text-align:center;">
            ${structHtml}
            <div style="font-size:0.68rem;color:var(--text-muted);margin-top:3px;">2D Structure</div>
          </div>
        </div>

        ${(r.deconvolved_query_peaks || []).length ? `<details class="pw-card" style="margin-top:0.75rem;">
          <summary>🧩 去卷积目标组分谱（仅分配自实测强度）</summary>
          <div id="deconvolvedDetailSpectrum" style="width:100%;height:260px;"><div class="result-visual-empty">展开后绘制去卷积目标组分谱</div></div>
        </details>` : ''}
          </div>
        </details>

        <div class="modal-actions">
          <button class="btn btn-secondary" onclick="App.closeModal()">关闭</button>
          <button class="btn btn-primary" onclick="App.navigateTo('export');document.getElementById('exportTaskSelect').value='${taskId}';ExportMod.updateOptions();App.closeModal();">📥 导出</button>
        </div>
      </div></div>`;

      if (inlineTargetId) {
        if (!_mountInlineDetail(html, inlineTargetId)) {
          Toast.show('候选详情工作台容器不可用', 'error');
          return;
        }
      } else {
        App.showModal(html);
      }
      if (isomerLiteratureRequest) {
        setTimeout(function () {
          _loadRemoteIsomerLiterature(isomerLiteratureRequest);
        }, 120);
      }
      // 渲染 MS/MS 棒状图 (SpectrumViewer, GNPS风格 m/z vs intensity)
      setTimeout(function () {
        try {
          if (renderSequence !== _detailRenderSequence) return;
          var qp = detailQueryPeaks;
          var el = document.getElementById("detailPrimarySpectrum");
          if (!el) return;
          el.innerHTML = '';
          const plotEnglish = window.MSSLanguage?.getLanguage?.()?.startsWith('en');
          if (qp.length && hasDetailMirrorEvidence && typeof SpectrumComparator !== "undefined") {
            new SpectrumComparator().render("detailPrimarySpectrum",
              { peaks: qp, name: savedScoringView
                  ? (plotEnglish ? "Scoring query" : "评分查询谱")
                  : (plotEnglish ? "Measured query" : "实测查询谱") },
              { peaks: detailLibraryPeaks, name: m.compound_name || "馆藏谱", smiles: m.smiles || "" },
              detailMatchedPairs, {
                height: 330, scorePrecision: 4,
                score: savedScoringView ? savedScoringView.comparison_cosine : m.cosine,
                scoreLabel: savedScoringView
                  ? (plotEnglish ? "Saved comparison cosine" : "已保存候选比较余弦")
                  : (plotEnglish ? "Saved cosine" : "已保存余弦"),
                title: plotEnglish ? "MS/MS spectrum comparison" : "MS/MS 谱图对照",
                text: plotEnglish ? {} : {
                  query: "查询谱", library: "参考谱",
                  matchedQuery: "查询谱 · 已记录匹配", otherQuery: "查询谱 · 其他峰",
                  matchedLibrary: "参考谱 · 已记录匹配", otherLibrary: "参考谱 · 其他峰",
                  locatedPairs: "已定位峰对／保存峰对", referenceMissing: "未保存参考谱",
                  empty: "没有可绘制的峰", relativeIntensity: "相对强度（%）",
                },
              });
          } else if (qp.length && typeof SpectrumViewer !== "undefined") {
            new SpectrumViewer("detailPrimarySpectrum").render(qp, { title: "MS/MS", height: 330 });
          } else {
            el.innerHTML = `<div class="result-visual-empty">${qp.length ? '谱图组件未加载；峰数据仍保留在结果记录中。' : '该结果没有保存可绘制的 MS/MS 碎片峰。'}</div>`;
          }
          var dp = r.deconvolved_query_peaks || [];
          var de = document.getElementById("deconvolvedDetailSpectrum");
          if (de && dp.length && typeof SpectrumViewer !== "undefined") {
            var deDetails = de.closest("details");
            var renderDeconvolved = function () {
              if (renderSequence !== _detailRenderSequence || de.dataset.rendered === "true") return;
              de.dataset.rendered = "true";
              de.innerHTML = "";
              new SpectrumViewer("deconvolvedDetailSpectrum").render(dp, { title: (m.compound_name || "") + " deconvolved target" });
            };
            if (deDetails && !deDetails.open) {
              deDetails.addEventListener("toggle", function () {
                if (deDetails.open) renderDeconvolved();
              }, { once: true });
            } else {
              renderDeconvolved();
            }
          }
        } catch (e) {
          var failedEl = document.getElementById("detailPrimarySpectrum");
          if (failedEl) failedEl.innerHTML = '<div class="result-visual-empty">谱图渲染失败；原始峰数据仍保留在结果记录中。</div>';
        }
      }, 80);

      // 异步加载 CFM-ID 碎裂产物结构级联 (前体→碎片 2D)
      if (smiles) {
        setTimeout(function () {
          fetch(`${API.BASE}/tasks/${taskId}/download/fragmentation/${index}?result_index=${Number.isInteger(resultIndex) ? resultIndex : index}`)
            .then(function (resp) { return resp.ok ? resp.text() : null; })
            .then(function (txt) {
              if (renderSequence !== _detailRenderSequence) return;
              const el = document.getElementById('cfmCascade' + index);
              if (!el) return;
              if (txt && txt.trim() && txt.indexOf('<div') >= 0) {
                el.innerHTML = '<div style="font-size:0.72rem;color:var(--text-secondary);margin-bottom:4px;font-weight:600;">🧬 候选产物结构证据</div>' + txt;
              } else {
                el.innerHTML = '<div style="font-size:0.7rem;color:#64748b;padding:0.3rem 0;">当前没有可展示的候选产物结构；这不等于没有 MS/MS，观测碎片与未归属质量差仍保留在下方证据表。</div>';
              }
            })
            .catch(function () {
              if (renderSequence !== _detailRenderSequence) return;
              const el = document.getElementById('cfmCascade' + index);
              if (el) el.innerHTML = '<div style="font-size:0.7rem;color:#999;padding:0.3rem 0;">碎裂结构加载失败</div>';
            });
        }, 120);
      }
    } catch (e) {
      if (inlineTargetId) {
        const target = document.getElementById(inlineTargetId);
        if (target) target.innerHTML = `<div class="empty-state"><div class="empty-icon">!</div><div class="empty-title">候选详情加载失败</div><div class="empty-desc">${Tasks.escapeHtml(e.message || String(e))}</div></div>`;
      }
      Toast.show(`加载详情失败: ${e.message}`, 'error');
    }
  }

  async function showDetailByResultIndex(taskId, resultIndex) {
    const wanted = Number(resultIndex);
    if (!taskId || !Number.isInteger(wanted)) {
      Toast.show('候选原始结果定位不可用', 'error');
      return;
    }
    if (_currentTaskId !== taskId || !_currentMatched.length) {
      await loadForTask(taskId);
    }
    const index = _currentMatched.findIndex(row => Number(row.result_index) === wanted);
    if (index < 0) {
      Toast.show('未能在当前任务中定位该候选，请重新选择样本', 'error');
      return;
    }
    await showDetail(taskId, index);
  }

  async function showDetailFromSummary(taskId, summaryRow, targetId, matchedRows = [], task = {}) {
    const resultIndex = Number(summaryRow?.result_index);
    if (!taskId || !summaryRow || !Number.isInteger(resultIndex)) {
      Toast.show('候选原始结果定位不可用', 'error');
      return;
    }
    setTaskContext(taskId, matchedRows, task);
    await showDetail(taskId, resultIndex, targetId, summaryRow);
  }

  function filterCompoundRows(value) {
    const query = String(value || '').trim().toLowerCase();
    const rows = Array.from(document.querySelectorAll('#compoundResultsBody .compound-result-row'));
    let visible = 0;
    rows.forEach((row) => {
      const show = !query || String(row.dataset.search || '').includes(query);
      row.style.display = show ? '' : 'none';
      if (show) visible++;
    });
    const count = document.getElementById('compoundResultSearchCount');
    if (count) count.textContent = query ? `显示 ${visible}/${rows.length}` : '';
  }

  function filterGovernanceBlockedRows(value) {
    const query = String(value || '').trim().toLowerCase();
    const numeric = /^[0-9]+([.][0-9]+)?$/.test(query);
    const queryMz = numeric ? Number(query) : NaN;
    const rows = Array.from(document.querySelectorAll(
      '#governanceBlockedBody .governance-blocked-row'
    ));
    let visible = 0;
    rows.forEach((row) => {
      const mz = Number(row.dataset.mz);
      const show = !query || (numeric
        ? Number.isFinite(mz) && Math.abs(mz - queryMz) <= 0.5
        : String(row.dataset.search || '').includes(query));
      row.style.display = show ? '' : 'none';
      if (show) visible++;
    });
    const count = document.getElementById('governanceBlockedSearchCount');
    if (count) count.textContent = query ? `显示 ${visible}/${rows.length}` : '';
  }

  async function openGovernanceCohort() {
    try {
      const data = await API.get('/tasks?status=completed');
      const tasks = (data.tasks || []).slice(0, 100);
      if (tasks.length < 2) {
        Toast.show('至少需要两个已完成的独立进样任务', 'warning');
        return;
      }
      const rows = tasks.map(task => {
        const checked = task.id === _currentTaskId ? ' checked' : '';
        return `<tr class="evigate-cohort-row" data-task-id="${Tasks.escapeHtml(task.id)}">
          <td><input class="evigate-task-check" type="checkbox"${checked}></td>
          <td style="font-size:.72rem;max-width:270px;overflow:hidden;text-overflow:ellipsis;white-space:nowrap;" title="${Tasks.escapeHtml(task.filename || task.id)}">${Tasks.escapeHtml(task.filename || task.id)}</td>
          <td><select class="form-select evigate-task-role" style="min-width:95px;padding:2px 5px;font-size:.69rem;">
            <option value="sample">样品进样</option><option value="qc">QC</option><option value="blank">空白</option>
          </select></td>
          <td><input class="form-input evigate-independence" value="${Tasks.escapeHtml(task.id)}" style="min-width:150px;padding:3px 6px;font-size:.68rem;" title="相同标签会被折叠为一个独立投票"></td>
        </tr>`;
      }).join('');
      App.showModal(`<div class="modal-overlay" onclick="if(event.target===this)App.closeModal()">
        <div class="modal modal-lg">
          <div class="flex-between"><div><h3>跨进样一致性评估</h3><div style="font-size:.72rem;color:var(--text-secondary);margin-top:.2rem;">先折叠同一前体事件的多碰撞能扫描，再按 m/z/RT 对齐独立进样；不会合并谱图或改变当前首选候选。</div></div><button class="btn btn-secondary" onclick="App.closeModal()">关闭</button></div>
          <div style="margin:.6rem 0;padding:.5rem .65rem;background:#fffbeb;border-left:4px solid #d97706;color:#78350f;font-size:.72rem;line-height:1.5;">请只勾选应直接比较的重复进样。默认每个任务是一个独立<b>技术进样</b>，不是生物学重复；若两个任务不独立，请给它们填写相同“独立组”标签。把空白任务标记为“空白”后，系统会把重复出现的身份记为背景冲突，但不会删除 QC 背景记录。</div>
          <div class="table-wrap" style="max-height:360px;overflow:auto;"><table><thead><tr><th>选</th><th>任务</th><th>角色</th><th>独立组</th></tr></thead><tbody>${rows}</tbody></table></div>
          <div id="evigateCohortOutput" style="margin-top:.65rem;"></div>
          <div class="modal-actions"><button class="btn btn-primary" onclick="Results.runGovernanceCohort()">运行跨进样一致性评估</button><button class="btn btn-secondary" onclick="App.closeModal()">关闭</button></div>
        </div></div>`);
    } catch (error) {
      Toast.show(`无法加载进样任务: ${error.message}`, 'error');
    }
  }

  async function runGovernanceCohort() {
    const selected = Array.from(document.querySelectorAll('.evigate-cohort-row'))
      .filter(row => row.querySelector('.evigate-task-check')?.checked);
    if (selected.length < 2) {
      Toast.show('请至少选择两个已完成的进样任务', 'warning');
      return;
    }
    const taskIds = selected.map(row => row.dataset.taskId);
    const taskRoles = {};
    const independenceGroups = {};
    selected.forEach(row => {
      const taskId = row.dataset.taskId;
      taskRoles[taskId] = row.querySelector('.evigate-task-role')?.value || 'sample';
      independenceGroups[taskId] = String(
        row.querySelector('.evigate-independence')?.value || taskId
      ).trim() || taskId;
    });
    const output = document.getElementById('evigateCohortOutput');
    if (output) output.innerHTML = '<div style="padding:.8rem;text-align:center;color:#1e40af;">正在折叠前体事件并进行跨进样证据审计…</div>';
    try {
      const payload = await API.post('/evidence-governance/cohort', {
        task_ids: taskIds,
        task_roles: taskRoles,
        independence_groups: independenceGroups,
        limit: 500,
      });
      const summary = payload.summary || {};
      const actions = summary.action_counts || {};
      const actionLabels = {
        release_replicate_stable_candidate: '跨进样稳定候选',
        release_isomer_group_candidate: '一致到异构体组',
        abstain_blank_recurrent: '空白中重复出现',
        abstain_identity_conflict: '进样间候选不一致',
        abstain_scope_or_quality_gate: '证据条件不足',
        needs_more_independent_injections: '仍需更多独立进样',
        no_candidate: '无候选',
      };
      const actionColor = action => action === 'release_replicate_stable_candidate'
        ? '#166534' : action === 'release_isomer_group_candidate' ? '#92400e'
          : action.includes('abstain') ? '#991b1b' : '#475569';
      const rows = (payload.clusters || []).slice(0, 150).map(cluster => {
        const candidate = cluster.dominant_candidate || {};
        return `<tr>
          <td>${Tasks.escapeHtml(cluster.cluster_id || '—')}</td>
          <td>${Tasks.escapeHtml(cluster.polarity || '—')}</td>
          <td class="mono">${cluster.precursor_mz == null ? '—' : Number(cluster.precursor_mz).toFixed(5)}</td>
          <td class="mono">${cluster.retention_time_min == null ? '—' : Number(cluster.retention_time_min).toFixed(3)}</td>
          <td>${Tasks.escapeHtml(candidate.compound_name || '—')}</td>
          <td>${Number(cluster.independent_injection_groups || 0)}/${Number(cluster.evaluable_injection_groups || 0)}</td>
          <td style="color:${actionColor(String(cluster.action || ''))};font-weight:700;">${Tasks.escapeHtml(actionLabels[cluster.action] || cluster.action || '—')}</td>
        </tr>`;
      }).join('');
      const actionBadges = Object.entries(actions).map(([action, count]) =>
        `<span style="display:inline-block;margin:2px 4px;padding:2px 8px;border-radius:10px;background:#f1f5f9;color:${actionColor(action)};">${Tasks.escapeHtml(actionLabels[action] || action)} <b>${Number(count)}</b></span>`
      ).join('');
      if (output) output.innerHTML = `<div style="padding:.6rem .7rem;border-left:4px solid #2563eb;background:#eff6ff;font-size:.74rem;line-height:1.55;">
          <b>评估完成：</b>${Number(summary.input_features || 0)} 个进样内特征折叠为 ${Number(summary.aligned_feature_clusters || 0)} 个对齐簇；形成 ${Number(summary.released_scope_unique_units ?? summary.released_unique_identities ?? 0)} 个跨进样一致的去重报告单元（其中 ${Number(summary.released_unique_identities || 0)} 个一致到连接结构，其余一致到异构体组）；折叠 ${Number(summary.pseudoreplicate_scans_collapsed || 0)} 个相关扫描投票。<b>当前候选排序未改变。</b><br>${actionBadges}
        </div>
        <div class="table-wrap" style="max-height:330px;overflow:auto;margin-top:.45rem;"><table><thead><tr><th>簇</th><th>极性</th><th>m/z</th><th>RT</th><th>主候选</th><th>独立支持</th><th>跨进样结论</th></tr></thead><tbody>${rows}</tbody></table></div>
        ${Number(payload.clusters_total || 0) > 150 ? '<div style="font-size:.68rem;color:#64748b;margin-top:.25rem;">界面仅显示前 150 个簇；完整结果可由 API 分页读取。</div>' : ''}`;
    } catch (error) {
      if (output) output.innerHTML = `<div style="color:#991b1b;padding:.6rem;">审计失败：${Tasks.escapeHtml(error.message || 'unknown error')}</div>`;
    }
  }

  function setBlankView(value) {
    const requested = String(value || 'all');
    _blankView = ['all', 'clean', 'flagged'].includes(requested)
      ? requested : 'all';
    void loadForTask();
  }

  return {
    init, refreshSelect, loadForTask, setTaskContext,
    showDetail, showDetailByResultIndex, showDetailFromSummary,
    showUnmatchedDetail,
    showMoreUnmatched,
    filterCompoundRows, filterGovernanceBlockedRows, setBlankView, openGovernanceCohort,
    runGovernanceCohort, expectedCompoundsCard: _expectedCompoundsCard,
    renderGovernanceBlockedCard: _governanceBlockedCard,
    get currentTaskId() { return _currentTaskId; },
  };
})();
