/**
 * MSsystem — Result Export Module
 * Single-task and batch export with search, multiple formats
 */
const ExportMod = (() => {
  let _selectedBatchTasks = new Set();
  let _searchQuery = '';
  let _batchTasksById = new Map();
  let _supportsTargetedFeatureAudit = null;
  let _supportsCompleteMs2Csv = null;
  let _supportsIsomerCandidateCsv = null;
  let _csvCapabilityPromise = null;
  let _batchExportInFlight = false;
  let _singleExportInFlight = 0;
  const _frozenTasks = new Map();

  const exportPhaseNames = {
    verifying_library: '核验谱库文件', preparing: '准备导出', waiting_export_slot: '等待前一个导出释放资源',
    reading_results: '读取样本结果', checking_evidence: '核验已保存证据',
    validating_sources: '核验报告数据', projecting_results: '整理报告结果',
    rendering_worksheets: '写入工作表', styling_workbook: '设置表格格式',
    rendering_raw_peaks: '写入完整原始峰', finalizing_xlsx: '压缩工作簿', xlsx_ready: '工作簿已生成',
    writing_csv: '写入 CSV', writing_report: '生成报告', packing_raw_peaks: '打包原始峰附件',
    saving: '保存并校验文件', ready: '文件可下载', failed: '导出失败',
  };

  function showExportProgress(job) {
    const single = job.request?.mode === 'single';
    const host = document.getElementById(single ? 'singleExportProgress' : 'batchExportProgress');
    if (!host) return;
    host.hidden = false;
    const cardId = `exportProgress-${job.id}`;
    let card = document.getElementById(cardId);
    if (!card) {
      card = document.createElement('section');
      card.id = cardId;
      card.className = 'card export-progress-card';
      card.style.cssText = 'margin:12px 0;padding:16px 20px;border:1px solid var(--border-color,#dce4ee);border-left:4px solid #3972bd;overflow-wrap:anywhere';
      card.innerHTML = '<div style="display:flex;align-items:center;justify-content:space-between;gap:12px;flex-wrap:wrap"><strong data-title></strong><div style="display:flex;gap:8px;align-items:center"><a data-download class="btn btn-primary btn-sm" hidden>下载文件</a><button data-toggle class="btn btn-secondary btn-sm" type="button" aria-expanded="true">收起详情</button></div></div><div data-file style="margin-top:8px;font-size:.85rem;color:var(--text-secondary)"></div><div data-details style="margin-top:10px"><div data-count role="status" aria-live="polite"></div><progress aria-label="导出生成进度" style="width:100%;height:8px;margin:8px 0"></progress><div data-phase style="font-size:.85rem"></div><a data-status target="_blank" rel="noopener" style="display:inline-block;margin-top:8px;font-size:.85rem">在独立页面查看进度</a></div>';
      card.querySelector('[data-toggle]').addEventListener('click', () => {
        const details = card.querySelector('[data-details]');
        details.hidden = !details.hidden;
        const toggle = card.querySelector('[data-toggle]');
        toggle.textContent = details.hidden ? '展开详情' : '收起详情';
        toggle.setAttribute('aria-expanded', String(!details.hidden));
      });
      host.prepend(card);
    }
    const previous = card.dataset.status;
    card.dataset.status = job.status;
    const done = job.status === 'completed';
    const failed = job.status === 'failed';
    const format = String(job.request?.format || '').toUpperCase();
    card.querySelector('[data-title]').textContent = `${single ? '单样本' : '批量'} ${format} · ${done ? (job.stale ? '历史文件' : '导出完成') : failed ? '导出失败' : job.status === 'queued' ? '排队中' : '正在生成'}`;
    card.style.borderLeftColor = failed ? '#b34135' : done ? '#2b8060' : '#3972bd';
    const taskId = job.current_task || job.request?.task_ids?.[0];
    const option = Array.from(document.getElementById('exportTaskSelect')?.options || []).find(item => item.value === taskId);
    const sampleName = _batchTasksById.get(taskId)?.filename || option?.textContent?.replace(/ \(.*$/, '') || taskId || '';
    const fileName = done ? _responseDownloadName({headers:new Headers(job.headers || {})}, sampleName) : sampleName;
    card.querySelector('[data-file]').textContent = (single ? fileName : `${job.total} 个样本${done ? ' · ' + fileName : ' · 当前：' + sampleName}`)
      + (done && job.size_bytes ? ` · ${(job.size_bytes / 1048576).toFixed(1)} MB` : '');
    card.querySelector('[data-count]').textContent = `已完成 ${job.completed || 0} / ${job.total || 0} 个样本；失败 ${job.failed || 0} 个`;
    const progress = card.querySelector('progress');
    progress.max = job.total || 1;
    if (job.status === 'queued' || job.phase === 'preparing' || (single && job.status === 'running')) progress.removeAttribute('value');
    else progress.value = job.completed || 0;
    const elapsed = Math.max(0, Math.floor((done || failed ? job.updated_at || job.created_at : Date.now()/1000) - job.created_at));
    card.querySelector('[data-phase]').textContent = job.error || job.notice || `${exportPhaseNames[job.phase] || job.status} · 已用时 ${elapsed} 秒`;
    const download = card.querySelector('[data-download]');
    download.hidden = !done;
    download.style.display = done ? 'inline-flex' : 'none';
    download.href = `${API.BASE}/export/jobs/${job.id}/download`;
    download.setAttribute('download', '');
    card.querySelector('[data-status]').href = `/static/export-progress.html?job=${job.id}`;
    if (done && previous !== 'completed') {
      card.querySelector('[data-details]').hidden = true;
      card.querySelector('[data-toggle]').textContent = '展开详情';
      card.querySelector('[data-toggle]').setAttribute('aria-expanded', 'false');
    }
    // Keep recent files nearby without turning the format picker into a long log.
    Array.from(host.querySelectorAll('[data-status="completed"], [data-status="failed"]')).slice(3).forEach(old => old.remove());
  }

  async function watchActiveExports() {
    const pending = new Set();
    try {
      JSON.parse(localStorage.getItem('mssystem.exportWatchIds') || '[]').forEach(id => pending.add(id));
      const single = localStorage.getItem('mssystem.lastSingleExport');
      const batch = JSON.parse(localStorage.getItem('mssystem.pendingBatchExport') || 'null');
      if (single) pending.add(single);
      if (batch?.id) pending.add(batch.id);
    } catch {}
    for (;;) {
      try {
        if (!_batchExportInFlight && !_singleExportInFlight) {
          const response = await fetch(`${API.BASE}/export/jobs`, {cache:'no-store'});
          if (response.ok) {
            const jobs = (await response.json()).jobs || [];
            const activeIds = new Set(jobs.map(job => job.id));
            const finished = await Promise.all(Array.from(pending).filter(id => !activeIds.has(id)).map(async id => {
              const result = await fetch(`${API.BASE}/export/jobs/${id}`, {cache:'no-store'});
              if (!result.ok) { pending.delete(id); return null; }
              return result.json();
            }));
            [...jobs, ...finished.filter(Boolean)].forEach(job => {
              showExportProgress(job);
              if (['completed', 'failed'].includes(job.status)) pending.delete(job.id);
              else pending.add(job.id);
            });
            localStorage.setItem('mssystem.exportWatchIds', JSON.stringify(Array.from(pending)));
          }
        }
      } catch {}
      await new Promise(resolve => setTimeout(resolve, 3000));
    }
  }

  function interfaceLanguage(requested) {
    if (requested === 'en' || requested === 'zh') return requested;
    return window.MSSLanguage?.getLanguage() || 'zh';
  }

  // ── Initialization ────────────────────────────────────────────────
  function init() {
    document.getElementById('exportLiveProgress')?.remove();
    const sel = document.getElementById('exportTaskSelect');
    if (sel) sel.addEventListener('change', () => updateOptions());

    const searchInput = document.getElementById('exportSearch');
    if (searchInput) {
      searchInput.addEventListener('input', (e) => {
        _searchQuery = e.target.value;
        filterExportSelect();
      });
    }
    _csvCapabilityPromise = refreshCsvLayerCapabilities();
    watchActiveExports();
  }

  async function refreshCsvLayerCapabilities() {
    const auditCard = document.getElementById('targetedFeatureAuditCard');
    try {
      const response = await fetch(
        `/openapi.json?capability_nonce=${Date.now()}`,
        { cache: 'no-store' }
      );
      if (!response.ok) {
        _supportsTargetedFeatureAudit = false;
        _supportsCompleteMs2Csv = false;
        _supportsIsomerCandidateCsv = false;
        return;
      }
      const schema = await response.json();
      const parameters = schema?.paths?.[
        '/api/tasks/{task_id}/download/csv'
      ]?.get?.parameters || [];
      const layerParameter = parameters.find(
        parameter => parameter.name === 'result_layer'
      );
      const supported = layerParameter?.schema?.enum || [];
      _supportsTargetedFeatureAudit = supported.includes(
        'targeted_features'
      );
      _supportsIsomerCandidateCsv = supported.includes(
        'isomer_candidates'
      );
      _supportsCompleteMs2Csv = parameters.some(
        parameter => parameter.name === 'ms2_schema'
      );
      if (auditCard && _supportsTargetedFeatureAudit) {
        auditCard.style.display = '';
      }
      const title = document.getElementById('targetedCsvTitle');
      const description = document.getElementById(
        'targetedCsvDescription'
      );
      const contract = document.getElementById('csvLayerContract');
      if (title) title.textContent = 'CSV 靶向化合物汇总';
      if (description) {
        description.textContent =
          '每个去重靶向化合物一行；候选数单列，代表峰与相关加合物已归并';
      }
      if (contract) {
        contract.innerHTML =
          'Excel 已按“结果导航 → 靶向高证据 → MS2 谱库鉴定 → 异构体竞争候选明细 → 同位素支持 → 仅质量候选”分层；' +
          'MS2 主 CSV 已将最佳候选与全部竞争异构体逐行合并，并保留碎片、分类和谱库追溯证据；候选专表仅作可选审计。' +
          '靶向层继续分为“化合物汇总”和“逐峰审计”。' +
          '<b>该层不计入普通 MS2 的匹配数、High 数或去重候选数。</b>';
      }
    } catch (error) {
      _supportsTargetedFeatureAudit = false;
      _supportsCompleteMs2Csv = false;
      _supportsIsomerCandidateCsv = false;
      console.debug('CSV layer capability check unavailable:', error);
    }
  }

  // ── Refresh Task Select ───────────────────────────────────────────
  async function refreshSelect() {
    try {
      const data = await API.get('/tasks?status=completed');
      _frozenTasks.clear();
      (data.tasks || []).filter(task => task.completed_batch).forEach(task => _frozenTasks.set(task.id, task));
      const sel = document.getElementById('exportTaskSelect');
      if (!sel) return;
      const currentVal = sel.value;
      sel.innerHTML = '<option value="">-- 选择已完成的任务 --</option>' +
        (data.tasks || []).map(t =>
          `<option value="${t.id}" data-name="${Tasks.escapeHtml(t.filename).toLowerCase()}" data-filename="${Tasks.escapeHtml(t.filename)}" data-id="${t.id}" data-compounds="${Tasks.escapeHtml((t.compound_names || []).join('; ').toLowerCase())}">
            ${Tasks.escapeHtml(t.filename)} (${t.matched_spectra}/${t.total_spectra} 匹配, ${t.match_rate}%)
          </option>`
        ).join('');
      if (currentVal && data.tasks?.find(t => t.id === currentVal)) sel.value = currentVal;
      updateOptions();
      loadBatchList();
    } catch (e) {
      console.error('Failed to refresh export select:', e);
    }
  }

  function filterExportSelect() {
    const sel = document.getElementById('exportTaskSelect');
    if (!sel) return;
    const q = _searchQuery.toLowerCase();
    Array.from(sel.options).forEach(opt => {
      if (!opt.value) return; // skip placeholder
      const name = (opt.dataset.name || '').toLowerCase();
      const id = (opt.dataset.id || '').toLowerCase();
      const compounds = (opt.dataset.compounds || '').toLowerCase();
      opt.style.display = (!q || name.includes(q) || id.includes(q) || compounds.includes(q)) ? '' : 'none';
    });
  }

  function updateOptions() {
    const taskId = document.getElementById('exportTaskSelect')?.value;
    const opts = document.getElementById('exportOptions');
    const frozen = _frozenTasks.get(taskId);
    if (opts) opts.classList.toggle('hidden', !taskId || !!frozen);
    const saved = document.getElementById('frozenSingleExport');
    if (saved) saved.innerHTML = frozen && typeof FrozenTaskUI !== 'undefined' ? FrozenTaskUI.panel(frozen) : '';
  }

  function getSelectedTaskId() {
    return document.getElementById('exportTaskSelect')?.value || '';
  }

  const TARGET_EVIDENCE_RANK = {
    MS1_EXACT_MASS_CANDIDATE: 1,
    MS1_ISOTOPE_SUPPORTED: 2,
    MS1_MS2_DIAGNOSTIC_PROBABLE: 3,
    EXPERIMENTAL_LIBRARY_CONFIRMED: 4,
  };

  const TARGET_EVIDENCE_LABEL = {
    MS1_EXACT_MASS_CANDIDATE: 'MS1 精确质量候选',
    MS1_ISOTOPE_SUPPORTED: 'MS1 + 同位素支持',
    MS1_MS2_DIAGNOSTIC_PROBABLE: 'MS1 + 诊断离子筛查支持',
    EXPERIMENTAL_LIBRARY_CONFIRMED: 'MS/MS 实验谱库候选支持',
  };

  function _uniqueText(values) {
    const seen = new Set();
    const output = [];
    (values || []).forEach(value => {
      const text = String(value || '').trim();
      if (text && !seen.has(text)) {
        seen.add(text);
        output.push(text);
      }
    });
    return output;
  }

  function _targetIdentity(member) {
    const inchikey = String(member?.inchikey || '').trim().toUpperCase();
    if (inchikey) return `inchikey::${inchikey}`;
    const stored = String(member?.identity_key || '').trim();
    if (stored) return stored;
    const targetId = String(
      member?.target_id || member?.catalog_row_id || ''
    ).trim();
    if (targetId) return `target::${targetId}`;
    return `name::${String(member?.compound_name || '').trim().toLowerCase()}|${String(member?.formula || '').trim()}`;
  }

  function _selectedTargetMembers(feature) {
    const uniqueMembers = new Map();
    (feature.target_members || []).forEach(member => {
      if (member && typeof member === 'object') {
        const key = _targetIdentity(member);
        if (key && !uniqueMembers.has(key)) uniqueMembers.set(key, member);
      }
    });
    const members = [...uniqueMembers.values()];
    const experimental = (feature.experimental_library_support || []).filter(
      item => item && typeof item === 'object'
    );
    const resolving = experimental.length
      ? experimental
      : (feature.ms2_diagnostic_support || []).filter(
          item => item && typeof item === 'object'
        );
    if (!resolving.length) return members;
    const targetIds = new Set(
      resolving.map(item => String(item.target_id || '').trim()).filter(Boolean)
    );
    const inchikeys = new Set(
      experimental
        .map(item => String(item.target_inchikey || '').trim().toUpperCase())
        .filter(Boolean)
    );
    const selected = members.filter(member =>
      targetIds.has(String(member.target_id || '').trim()) ||
      inchikeys.has(String(member.inchikey || '').trim().toUpperCase())
    );
    return selected.length ? selected : members;
  }

  function _targetDatabaseLabel(values) {
    return _uniqueText(values).map(value =>
      value.replace(/^(user_contrib:|literature:|database:)/, '')
    ).join(' | ');
  }

  function _openmsLabel(feature) {
    return ({
      linked_openms_feature: '已由独立 OpenMS 峰验证',
      ambiguous_openms_features: '已关联多个 OpenMS 峰',
      no_openms_feature_in_mz_rt_window: '未关联到 OpenMS 峰',
      openms_features_unavailable: 'OpenMS 峰表不可用',
      target_feature_missing_mz_or_rt: '缺少 m/z 或 RT',
    })[String(feature.openms_feature_link_status || '')] || '旧任务未计算';
  }

  function _backgroundLabel(feature) {
    const key = String(
      feature.blank_evidence_category || feature.blank_status || ''
    );
    return ({
      not_requested: '',
      no_blank_feature: '空白中未检出',
      sample_enriched: '样本富集',
      blank_associated: '空白相关',
      blank_dominant: '背景主导',
      low_blank_prevalence: '空白中低频检出',
      unevaluable: '空白证据不可判定',
    })[key] || '';
  }

  function _supportCount(items, library = false) {
    return new Set((items || []).map(item =>
      `${String(item?.spectrum_id || '')}|${
        library ? String(item?.library_spectrum_id || '') : ''
      }`
    )).size;
  }

  function _targetRawRow(feature) {
    const members = _selectedTargetMembers(feature);
    let names = _uniqueText(members.map(member => member.compound_name));
    if (!names.length) names = _uniqueText(feature.candidate_names || []);
    if (!names.length && feature.candidate_name) {
      names = [String(feature.candidate_name)];
    }
    const identityKeys = _uniqueText(members.map(_targetIdentity));
    const diagnostic = feature.ms2_diagnostic_support || [];
    const library = feature.experimental_library_support || [];
    const resolved = diagnostic.length || library.length;
    const candidateCount = resolved
      ? identityKeys.length
      : Number(
          feature.identity_count || identityKeys.length || names.length || 0
        );
    const databases = _uniqueText(
      (feature.dataset_sources || []).length
        ? feature.dataset_sources
        : members.map(member => member.dataset_source)
    );
    return {
      feature,
      names,
      identityKeys,
      candidateCount,
      databases,
      diagnostic,
      library,
      evidenceLevel: String(
        feature.evidence_level || 'MS1_EXACT_MASS_CANDIDATE'
      ),
      formula: String(feature.formula || ''),
      polarity: String(feature.polarity || ''),
      adduct: String(feature.adduct || ''),
      targetMz: Number(feature.target_mz),
      observedMz: Number(feature.observed_apex_mz),
      massError: Number(feature.mass_error_ppm),
      rt: Number(feature.apex_rt_minutes),
      area: Number(
        feature.baseline_corrected_area_intensity_seconds ||
        feature.xic_area_intensity_seconds || 0
      ),
      height: Number(feature.apex_intensity || 0),
      snr: Number(feature.signal_to_noise || 0),
      peakQuality: Number(feature.peak_quality_score),
      openms: _openmsLabel(feature),
      background: _backgroundLabel(feature),
    };
  }

  function _targetRowPriority(row) {
    return [
      TARGET_EVIDENCE_RANK[row.evidenceLevel] || 0,
      row.feature.blank_flagged ? 0 : 1,
      ['linked_openms_feature', 'ambiguous_openms_features'].includes(
        String(row.feature.openms_feature_link_status || '')
      ) ? 1 : 0,
      Number.isFinite(row.peakQuality) ? row.peakQuality : 0,
      row.height,
      row.area,
    ];
  }

  function _comparePriority(left, right) {
    const a = _targetRowPriority(left);
    const b = _targetRowPriority(right);
    for (let index = 0; index < a.length; index += 1) {
      if (a[index] !== b[index]) return a[index] - b[index];
    }
    return 0;
  }

  function _groupTargetedFeatures(features) {
    const groups = new Map();
    (features || []).map(_targetRawRow).forEach(row => {
      const identities = row.identityKeys.length
        ? [...row.identityKeys].sort()
        : row.names.map(name => `name::${name.toLowerCase()}`).sort();
      const key = JSON.stringify([identities, row.formula, row.polarity]);
      if (!groups.has(key)) groups.set(key, []);
      groups.get(key).push(row);
    });

    const output = [];
    groups.forEach(rows => {
      const representative = rows.reduce(
        (best, row) => _comparePriority(row, best) > 0 ? row : best,
        rows[0]
      );
      const component = rows.filter(
        row => Number.isFinite(row.rt) &&
          Number.isFinite(representative.rt) &&
          Math.abs(row.rt - representative.rt) <= 0.05
      );
      const componentRows = component.length ? component : [representative];
      const adductAreas = new Map();
      componentRows.forEach(row => {
        const adduct = row.adduct || 'unassigned';
        adductAreas.set(
          adduct,
          Math.max(adductAreas.get(adduct) || 0, row.area)
        );
      });
      const names = _uniqueText(rows.flatMap(row => row.names));
      const identityCount = Math.max(
        ...rows.map(row => row.candidateCount),
        ...rows.map(row => row.identityKeys.length),
        0
      );
      const primary = names[0] || '';
      const displayedName = identityCount <= 1
        ? primary
        : `${primary} 等 ${identityCount} 个候选`;
      const databases = _uniqueText(rows.flatMap(row => row.databases));
      const diagnostic = componentRows.flatMap(row => row.diagnostic);
      const library = componentRows.flatMap(row => row.library);
      const openmsRow = componentRows.find(row =>
        ['linked_openms_feature', 'ambiguous_openms_features'].includes(
          String(row.feature.openms_feature_link_status || '')
        )
      ) || representative;
      const backgroundRank = {
        blank_dominant: 4,
        blank_associated: 3,
        low_blank_prevalence: 2,
        sample_enriched: 1,
      };
      const backgroundRow = componentRows.reduce((worst, row) => {
        const rowKey = String(
          row.feature.blank_evidence_category ||
          row.feature.blank_status || ''
        );
        const worstKey = String(
          worst.feature.blank_evidence_category ||
          worst.feature.blank_status || ''
        );
        return (backgroundRank[rowKey] || 0) > (
          backgroundRank[worstKey] || 0
        ) ? row : worst;
      }, representative);
      output.push({
        evidenceLevel: representative.evidenceLevel,
        values: [
          TARGET_EVIDENCE_LABEL[representative.evidenceLevel] ||
            representative.evidenceLevel,
          displayedName,
          identityCount,
          _targetDatabaseLabel(databases),
          representative.formula,
          representative.adduct,
          [...adductAreas.keys()].sort().join(' | '),
          representative.polarity,
          Number.isFinite(representative.targetMz)
            ? representative.targetMz : '',
          Number.isFinite(representative.observedMz)
            ? representative.observedMz : '',
          Number.isFinite(representative.massError)
            ? representative.massError : '',
          Number.isFinite(representative.rt) ? representative.rt : '',
          [...adductAreas.values()].reduce((sum, value) => sum + value, 0),
          Number.isFinite(representative.peakQuality)
            ? representative.peakQuality : '',
          openmsRow.openms,
          backgroundRow.background,
          rows.length,
          representative.snr,
          _supportCount(diagnostic),
          _supportCount(library, true),
        ],
      });
    });
    output.sort((left, right) => {
      const evidence = (
        TARGET_EVIDENCE_RANK[right.evidenceLevel] || 0
      ) - (TARGET_EVIDENCE_RANK[left.evidenceLevel] || 0);
      if (evidence) return evidence;
      return Number(right.values[12] || 0) - Number(left.values[12] || 0);
    });
    return output;
  }

  async function _fetchAllExpectedFeatures(taskId) {
    const features = [];
    let offset = 0;
    const limit = 5000;
    while (true) {
      const response = await fetch(
        `${API.BASE}/tasks/${encodeURIComponent(taskId)}/expected-compounds?limit=${limit}&offset=${offset}`,
        { cache: 'no-store' }
      );
      if (!response.ok) {
        throw new Error(`靶向结果读取失败 (HTTP ${response.status})`);
      }
      const payload = await response.json();
      const page = Array.isArray(payload.features) ? payload.features : [];
      features.push(...page);
      if (!payload.has_more || !page.length) break;
      offset += page.length;
    }
    return features;
  }

  function _csvSafeText(value) {
    let text = String(value ?? '')
      .replace(/[\r\n\u2028\u2029]+/g, ' ')
      .trim();
    const limit = 32000;
    const marker = ' …[单元格内容过长，标准表已截断；完整数据请查阅详细审计或任务原始数据]';
    if (text.length > limit) {
      text = text.slice(0, limit - marker.length) + marker;
    }
    return text;
  }

  function _csvCell(value) {
    const text = _csvSafeText(value);
    return /[",\r\n]/.test(text)
      ? `"${text.replaceAll('"', '""')}"`
      : text;
  }

  function _csvText(headers, rows) {
    if (!Array.isArray(headers) || !headers.length) {
      throw new Error('CSV 字段定义为空');
    }
    rows.forEach((row, index) => {
      if (!Array.isArray(row) || row.length !== headers.length) {
        throw new Error(
          `CSV 第 ${index + 1} 行字段数不一致：` +
          `${Array.isArray(row) ? row.length : '非数组'} != ${headers.length}`
        );
      }
    });
    return '\uFEFF' + [
      headers.map(_csvCell).join(','),
      ...rows.map(row => row.map(_csvCell).join(',')),
    ].join('\r\n') + '\r\n';
  }

  function _saveBlob(filename, blob) {
    const url = URL.createObjectURL(
      blob
    );
    const anchor = document.createElement('a');
    anchor.href = url;
    anchor.download = filename;
    anchor.style.display = 'none';
    document.body.appendChild(anchor);
    anchor.click();
    anchor.remove();
    // Large batch archives are not always opened by the browser immediately.
    // Revoking after one second made otherwise successful exports appear to do
    // nothing on slower disks/browsers.  Keep the URL alive long enough for
    // the download manager to take ownership without retaining it forever.
    setTimeout(() => URL.revokeObjectURL(url), 60_000);
  }

  function _saveCsv(filename, headers, rows) {
    _saveBlob(
      filename,
      new Blob(
        [_csvText(headers, rows)],
        { type: 'text/csv;charset=utf-8' }
      )
    );
  }

  const MS2_READER_HEADERS = [
    'Feature ID', 'Result ID', '结果状态', '鉴定组',
    '候选化合物', '组内排名', '候选角色', '组内候选数',
    '分子式', '鉴定结论', '鉴定组置信度', '候选置信度', 'RT (min)',
    'MS2 隔离目标 m/z', '检索用实测/回退 m/z', '检索质量轴来源',
    '离子模式', '样本 CE', '样本 CE 单位',
    'CE 可比性状态', 'CE 是否影响评分排序',
    '碎片信息状态', '有效碎片数',
    '背景 B1', '背景 B2', '背景 B3', '背景处置',
    '样本仪器类型', '样品离子形式', '样品离子形式依据',
    '候选谱库加合物', '代表馆藏 CE', '代表馆藏 CE 单位',
    '代表馆藏仪器类型', '候选 MS/MS 相似度', '候选谱熵相似度',
    '候选综合匹配分', '候选质量误差 (ppm)', '代表谱匹配碎片数',
    '代表谱匹配碎片 m/z', '代表谱匹配碎片 (m/z:强度)',
    '样本 MS2 碎片离子（实际检索谱，m/z:强度）', '实际检索 MS2 峰数',
    '原始 MS2 碎片离子（已保存输入峰，m/z:强度）', '原始 MS2 输入峰数',
    '候选化学类别', '鉴定组化学类别', '分类置信度', '分类依据',
    '分类证据状态', '分类可靠性说明', '谱库来源', '馆藏谱图编号',
    'InChIKey', 'SMILES', '样本谱图编号', '备注',
  ];

  const FORMULA_MONOISOTOPIC_MASS = {
    H: 1.00782503223, B: 11.00930536, C: 12.0,
    N: 14.00307400443, O: 15.99491461957,
    F: 18.99840316273, Na: 22.9897692820,
    Mg: 23.985041697, Al: 26.98153853,
    Si: 27.97692653465, P: 30.97376199842,
    S: 31.9720711744, Cl: 34.968852682,
    K: 38.9637064864, Ca: 39.962590863,
    Fe: 55.93493633, Br: 78.9183376, I: 126.904473,
  };

  const EXPORT_ADDUCTS = {
    POS: [
      ['[M+H]+', 1.007276],
      ['[M+Na]+', 22.989218],
      ['[M+K]+', 38.963158],
      ['[M+NH4]+', 18.033822],
      ['[M-H2O+H]+', -17.003288],
    ],
    NEG: [
      ['[M-H]-', -1.007276],
      ['[M+Cl]-', 34.968853],
      ['[M+FA-H]-', 44.998201],
      ['[M+Hac-H]-', 59.013851],
      ['[M-H-H2O]-', -19.018391],
    ],
  };

  function _numberOrBlank(value) {
    if (value === null || value === undefined || value === '' ||
        typeof value === 'boolean') return '';
    const numeric = Number(value);
    return Number.isFinite(numeric) ? numeric : value;
  }

  function _booleanLabel(value) {
    return value === true ? '是' : value === false ? '否' : '';
  }

  function _formulaExactMass(formula) {
    let text = String(formula || '').trim().replaceAll(' ', '');
    if (!text) return null;
    if (text.startsWith('[') && text.includes(']')) {
      const closing = text.lastIndexOf(']');
      text = text.slice(1, closing) + text.slice(closing + 1);
    }
    text = text.replace(/(?:\^\d*[+-]|[+-]\d*)$/, '');
    const token = /([A-Z][a-z]?)(\d*)/g;
    let cursor = 0;
    let total = 0;
    let match;
    while ((match = token.exec(text)) !== null) {
      if (match.index !== cursor ||
          FORMULA_MONOISOTOPIC_MASS[match[1]] === undefined) return null;
      total += FORMULA_MONOISOTOPIC_MASS[match[1]] *
        Number(match[2] || 1);
      cursor = token.lastIndex;
    }
    return cursor === text.length && total > 0 ? total : null;
  }

  function _specificAdduct(value) {
    const text = String(value || '').trim();
    if ([
      '', 'nan', 'none', 'null', 'n/a', 'unknown',
      'pos', 'positive', '+', 'neg', 'negative', '-',
    ].includes(text.toLowerCase())) return '';
    return text;
  }

  function _resolvedMs2Adduct(row, match) {
    const reported = _specificAdduct(match.adduct);
    if (reported) {
      return { value: reported, source: '谱库原始记录', errorPpm: '' };
    }
    const query = _specificAdduct(
      match.query_adduct || row.query_adduct
    );
    if (query && Boolean(
      match.query_adduct_reliable || row.query_adduct_reliable
    )) {
      return { value: query, source: '样本元数据', errorPpm: '' };
    }
    const charge = _numberOrBlank(row.precursor_charge);
    if (charge !== '' && Math.abs(Number(charge)) !== 1) {
      return {
        value: '未能可靠判定',
        source: '当前电荷态不支持安全推断',
        errorPpm: '',
      };
    }
    const precursor = Number(row.precursor_mz);
    const neutralMass = _formulaExactMass(match.formula);
    const polarity = String(row.ion_mode || match.query_ion_mode || '')
      .trim().toUpperCase();
    if (!Number.isFinite(precursor) || neutralMass === null ||
        !EXPORT_ADDUCTS[polarity]) {
      return {
        value: '未能可靠判定',
        source: !Number.isFinite(precursor)
          ? '缺少前体 m/z，未猜测'
          : '缺少分子式或极性，未猜测',
        errorPpm: '',
      };
    }
    const tolerance = Math.max(0.001, neutralMass * 10e-6);
    const matches = EXPORT_ADDUCTS[polarity].map(([name, shift]) => ({
      name,
      errorDa: Math.abs(precursor - shift - neutralMass),
    })).filter(candidate => candidate.errorDa <= tolerance);
    if (new Set(matches.map(candidate => candidate.name)).size !== 1) {
      return {
        value: '未能可靠判定',
        source: '无法唯一判定（未猜测）',
        errorPpm: '',
      };
    }
    const best = matches.sort(
      (left, right) => left.errorDa - right.errorDa
    )[0];
    return {
      value: best.name,
      source: '分子式 + 前体 m/z 唯一支持（10 ppm）',
      errorPpm: Number((best.errorDa / neutralMass * 1e6).toFixed(4)),
    };
  }

  function _cleanMs2CompoundName(match) {
    const key = String(match.inchikey || '').trim().toUpperCase();
    const exactAliases = {
      'DLIKSSGEMUFQOK-UHFFFAOYSA-N': 'Naringenin 7-O-glucoside',
      'YNWXJFQOCHMPCK-LXGDFETPSA-N': 'Isoliquiritin',
      'HKQYGTCOTHHOMP-UHFFFAOYSA-N': 'Formononetin',
      'MPYSZNHNOUOLSB-UHFFFAOYSA-N': '5-Hydroxybornyl 6-O-apiosyl-beta-D-glucoside',
    };
    if (exactAliases[key]) return exactAliases[key];
    let text = String(match.compound_name || '').trim()
      .replace(/^(?:MassBank|GNPS|ReSpect)\s*:\s*\S+\s+/i, '')
      .split('|', 1)[0].trim();
    text = text.replace(
      /\s*\(\s*known\s+isomers\s*:\s*\d+\s*;\s*isobaric\s+peaks\s*:\s*\d+\s*\)\s*$/i,
      ''
    ).trim();
    if (/^Prunin,\s*Naringenin-7-O-glucoside$/i.test(text)) {
      return 'Naringenin 7-O-glucoside';
    }
    return text;
  }

  function _collisionEnergyValue(value) {
    if (value && typeof value === 'object') {
      return value.value ?? value.collision_energy ?? '';
    }
    return value ?? '';
  }

  function _collisionEnergyUnit(value, explicitUnit) {
    const explicit = String(explicitUnit || '').trim();
    if (explicit) return explicit;
    return value && typeof value === 'object'
      ? String(value.unit || '').trim()
      : '';
  }

  function _libraryCollisionEnergy(match) {
    const raw = String(match.lib_collision_energy_raw ?? '').trim();
    return raw || _numberOrBlank(match.lib_collision_energy);
  }

  function _libraryCeStatus(match) {
    const value = _libraryCollisionEnergy(match);
    if (value === '') return '';
    return String(match.lib_collision_energy_unit || '').trim()
      ? '已报告（含单位）'
      : '已报告但单位缺失';
  }

  function _ms2PeakText(row) {
    return (row.query_peaks || row.raw_peaks || []).flatMap(peak => {
      const mz = Number(peak?.[0]);
      const intensity = Number(peak?.[1]);
      return Number.isFinite(mz) && Number.isFinite(intensity)
        ? [`${mz.toFixed(5)}:${intensity.toFixed(0)}`]
        : [];
    }).join(' ');
  }

  function _rawMs2PeakText(row) {
    return (row.raw_peaks || []).flatMap(peak => {
      const mz = Number(peak?.[0]);
      const intensity = Number(peak?.[1]);
      return Number.isFinite(mz) && Number.isFinite(intensity)
        ? [`${mz.toFixed(5)}:${intensity.toFixed(0)}`]
        : [];
    }).join(' ');
  }

  function _rawMs2PeakCount(row) {
    const explicit = Number(row.raw_peaks_input_count);
    if (Number.isFinite(explicit)) return explicit;
    if (Array.isArray(row.raw_peaks)) return row.raw_peaks.length;
    return _numberOrBlank(row.peaks_count);
  }

  function _matchedFragmentDetails(row, match) {
    const targets = String(match.ms_ms_fragments || '')
      .split(/[;,\s]+/)
      .map(Number)
      .filter(Number.isFinite);
    const peaks = (row.query_peaks || row.raw_peaks || []).flatMap(peak => {
      const mz = Number(peak?.[0]);
      const intensity = Number(peak?.[1]);
      return Number.isFinite(mz) && Number.isFinite(intensity)
        ? [[mz, intensity]]
        : [];
    });
    return targets.map(target => {
      const nearest = peaks.reduce((best, peak) =>
        !best || Math.abs(peak[0] - target) < Math.abs(best[0] - target)
          ? peak : best, null);
      return nearest && Math.abs(nearest[0] - target) <= 0.02
        ? `${nearest[0].toFixed(5)}:${Number(nearest[1].toPrecision(6))}`
        : target.toFixed(5);
    }).join('; ');
  }

  function _ms2IdentityKey(match) {
    const inchikey = String(match.inchikey || '').trim().toUpperCase();
    if (/^[A-Z]{14}-[A-Z]{10}-[A-Z]$/.test(inchikey)) {
      return `ik::${inchikey}`;
    }
    const smiles = String(match.smiles || '').trim();
    if (smiles) return `sm::${smiles}`;
    const name = String(match.compound_name || '')
      .toLowerCase().replace(/[^a-z0-9]+/g, '');
    const formula = String(match.formula || '').replaceAll(' ', '');
    if (name || formula) return `nf::${name}|${formula}`;
    return `sid::${String(match.lib_spectrum_id || '')}`;
  }

  function _ms2Priority(row) {
    const match = row.match || {};
    const ppm = Number(match.mass_ppm);
    return [
      match.potential_contaminant === true ? 1 : 0,
      -Number(match.rank_score ?? match.final_score ?? 0),
      -Number(match.final_score || 0),
      -Number(match.cosine || 0),
      Number.isFinite(ppm) ? Math.abs(ppm) : 999999,
      _ms2IdentityKey(match),
      String(match.dataset_source || ''),
      String(match.lib_spectrum_id || ''),
      String(row.spectrum_id || ''),
    ];
  }

  function _compareMs2Priority(left, right) {
    const a = _ms2Priority(left);
    const b = _ms2Priority(right);
    for (let index = 0; index < a.length; index += 1) {
      if (a[index] === b[index]) continue;
      return a[index] < b[index] ? -1 : 1;
    }
    return 0;
  }

  function _dedupMs2Results(results) {
    const best = new Map();
    (results || []).forEach(row => {
      const match = row?.match;
      if (!match || match.dataset_source === 'iceberg_fallback') return;
      const smiles = String(match.smiles || '').trim().toLowerCase();
      if (!smiles || ['nan', 'none', 'null'].includes(smiles)) return;
      const key = _ms2IdentityKey(match);
      if (!key) return;
      const previous = best.get(key);
      if (!previous || _compareMs2Priority(row, previous) < 0) {
        best.set(key, row);
      }
    });
    return [...best.values()].sort(_compareMs2Priority);
  }

  function _isomerMemberIdentity(member) {
    const inchikey = String(member?.inchikey || '').trim().toUpperCase();
    if (inchikey) return `inchikey::${inchikey}`;
    return `name::${String(member?.compound_name || '').trim().toLowerCase()}|${String(member?.formula || '').trim().toUpperCase()}`;
  }

  function _dedupIsomerMembers(group) {
    const supplied = Array.isArray(group?.candidate_members)
      ? group.candidate_members : [];
    const fallbackNames = Array.isArray(group?.candidate_names)
      ? group.candidate_names : [];
    const candidates = supplied.length ? supplied : fallbackNames.map(
      (name, index) => ({
        rank: index + 1,
        compound_name: name,
        formula: group?.formula || '',
        is_ranking_lead: index === 0,
        candidate_record_status: '旧任务仅保存候选名称',
      })
    );
    const seen = new Set();
    return candidates.filter(member => {
      const key = _isomerMemberIdentity(member);
      if (seen.has(key)) return false;
      seen.add(key);
      return true;
    }).map((member, index) => ({
      ...member,
      rank: index + 1,
      is_ranking_lead: index === 0,
      role: index === 0 ? 'ranking_lead' : 'competitive_isomer',
    }));
  }

  function _isomerGroupsByFormula(compounds) {
    const groups = new Map();
    (compounds || []).forEach(group => {
      const scope = String(group?.identity_scope || '').toLowerCase();
      const unresolved = scope === 'isomer_group'
        || scope === 'molecular_formula_isomer_group'
        || scope === 'stereoisomer_group'
        || group?.conclusion?.isomer_group_only === true;
      const formula = String(group?.formula || '').trim().toUpperCase();
      if (unresolved && formula && !groups.has(formula)) {
        groups.set(formula, group);
      }
    });
    return groups;
  }

  function _isUnresolvedCompound(group) {
    const scope = String(group?.identity_scope || '').toLowerCase();
    return scope === 'isomer_group'
      || scope === 'molecular_formula_isomer_group'
      || scope === 'stereoisomer_group'
      || group?.conclusion?.isomer_group_only === true;
  }

  function _classificationStatusLabel(match) {
    const status = String(match.classification_status || '').toLowerCase();
    if (status === 'reliable') return '可靠';
    if (status === 'tentative') return '待核查';
    if (status === 'unresolved') return '未可靠确定';
    const label = String(
      match.display_class || match.classification || ''
    ).trim().toLowerCase();
    if (!label || ['unknown', 'other', 'unclassified'].includes(label)) {
      return '未可靠确定';
    }
    return match.classification_reliable === true ? '可靠' : '待核查';
  }

  function _searchPrecursorSourceLabel(value) {
    const labels = {
      measured_continuous_ms1_feature_centroid: '连续 MS1 特征实测质心',
      ms2_selected_ion_mz_fallback: 'mzML selected-ion / 隔离目标回退',
    };
    const key = String(value || '').trim().toLowerCase();
    return labels[key] || String(value || '');
  }

  function _classificationSourceLabel(value) {
    const labels = {
      'runtime_structure_rule:stored_class_agrees_with_smiles_structure_rule': 'SMILES 结构规则与谱库类别一致',
      'runtime_structure_rule:stored_and_structure_classes_compatible': 'SMILES 结构规则与谱库类别相容',
      'runtime_structure_rule:trusted_source:curated_ciead_npclassifier': '结构规则 + 策展 NPClassifier 类别',
      'runtime_structure_rule:smiles_structure_rule_overrode_provisional_stored_class': 'SMILES 结构规则校正暂定谱库类别',
      ml: 'MS/MS 谱图分类模型',
      library: '谱库类别',
      missing: '未提供可靠分类依据',
    };
    const key = String(value || '').trim().toLowerCase();
    if (!key) return '';
    return labels[key] || '其他分类依据（详细代码见审计 CSV）';
  }

  function _classificationReasonLabel(value) {
    const labels = {
      stored_class_agrees_with_smiles_structure_rule: '谱库类别与 SMILES 结构规则一致',
      stored_and_structure_classes_compatible: '谱库类别与 SMILES 结构规则相容',
      'trusted_source:curated_ciead_npclassifier': '由策展 NPClassifier 来源支持',
      smiles_structure_rule_overrode_provisional_stored_class: 'SMILES 结构规则已校正暂定谱库类别',
      'untrusted_source:ml': '仅有模型分类，缺少独立结构证据',
      'untrusted_source:library': '仅有谱库类别标签，缺少独立结构证据',
      'untrusted_source:missing': '未记录可靠的分类来源',
      no_nontrivial_class: '未得到可解释的化学类别',
      independent_classification_evidence_unavailable: '缺少独立分类证据',
      missing_match: '缺少可用于分类的谱库匹配',
    };
    const key = String(value || '').trim().toLowerCase();
    if (!key) return '';
    if (labels[key]) return labels[key];
    if (key.startsWith('untrusted_source:')) {
      return `分类来源未达到独立可靠性门槛（${key.split(':', 2)[1].replaceAll('_', ' ')}）`;
    }
    return '其他分类可靠性原因（详细代码见审计 CSV）';
  }

  function _ms2ReaderRow(
    row, group = null, candidate = null, groupNumber = ''
  ) {
    const match = row.match || {};
    const unannotated = String(group?.conclusion?.key || '') === 'unannotated'
      || String(group?.record_type || '') === 'unannotated_ion_feature';
    const adduct = _resolvedMs2Adduct(row, match);
    const unresolved = _isUnresolvedCompound(group);
    const member = candidate || {};
    const lead = unresolved
      ? (member.is_ranking_lead === true
        || String(member.role || '') === 'ranking_lead'
        || Number(member.rank || 0) === 1)
      : true;
    const candidateName = unannotated
      ? ''
      : unresolved
      ? String(member.compound_name || '')
      : _cleanMs2CompoundName(match);
    const candidateFormula = unresolved
      ? String(member.formula || group?.formula || match.formula || '')
      : String(match.formula || '');
    const members = unresolved ? _dedupIsomerMembers(group) : [];
    const conclusion = String(
      group?.conclusion?.label
      || match.reportable_identity_label
      || match.evidence_tier
      || match.confidence
      || ''
    );
    const groupConfidence = String(
      group?.conclusion?.isomer_group_confidence
      || match.isomer_group_confidence
      || match.confidence
      || ''
    );
    const rawSpectrum = _rawMs2PeakText(row);
    const ce = match.ce_comparability || {};
    const background = row.background_evidence || {};
    const notes = [];
    const recordStatus = String(member.candidate_record_status || '').trim();
    if (recordStatus && recordStatus.toLowerCase() !== 'available') {
      notes.push(recordStatus);
    }
    if (row.raw_peaks_complete === false) {
      notes.push('原始峰按保存策略未完整保留');
    }
    if (rawSpectrum.length > 32000) {
      notes.push('原始峰列表超过单个表格单元格容量；该单元格将明确截断');
    }
    return [
      String(group?.feature_id || row.lcms_feature_id
        || row.precursor_feature_id || ''),
      String(group?.result_id || ''),
      unannotated ? '未注释离子特征' : '谱库候选结果',
      groupNumber,
      candidateName,
      unannotated ? '' : (unresolved ? _numberOrBlank(member.rank) : 1),
      unannotated
        ? '无谱库候选'
        : unresolved
        ? (lead ? '最佳候选（异构体未唯一解析）' : '竞争异构体候选')
        : '单一候选',
      unannotated ? 0 : (unresolved ? Math.max(1, members.length) : 1),
      candidateFormula,
      conclusion,
      groupConfidence,
      String(unresolved ? (member.confidence || '') : (match.confidence || '')),
      _numberOrBlank(row.retention_time),
      _numberOrBlank(row.precursor_mz),
      _numberOrBlank(row.search_precursor_mz ?? row.precursor_mz),
      _searchPrecursorSourceLabel(row.search_precursor_mz_source),
      String(row.ion_mode || ''),
      _numberOrBlank(_collisionEnergyValue(row.collision_energy)),
      _collisionEnergyUnit(row.collision_energy, row.collision_energy_unit),
      String(ce.status || (unannotated ? 'no_candidate' : '')),
      ce.ranking_influence === true ? '是' : '否',
      String(row.fragmentation_status || row.fragmentation_annotation?.status || ''),
      _numberOrBlank(row.effective_fragment_count
        ?? row.fragmentation_annotation?.effective_fragment_count),
      String(background.B1_cross_precursor || 'unknown'),
      String(background.B2_incoherent_chromatography || 'unknown'),
      String(background.B3_blank_confirmed || 'unknown'),
      String(background.action || 'retain'),
      String(row.instrument_type || match.query_instrument_type || ''),
      adduct.value,
      adduct.source,
      String(unresolved ? (member.adduct || '') : (match.adduct || '')),
      _libraryCollisionEnergy(match),
      String(match.lib_collision_energy_unit || ''),
      String(match.lib_instrument_type || ''),
      _numberOrBlank(unresolved ? member.cosine : match.cosine),
      _numberOrBlank(unresolved ? member.entropy : match.entropy),
      _numberOrBlank(unresolved ? member.final_score : match.final_score),
      _numberOrBlank(unresolved ? member.mass_ppm : match.mass_ppm),
      _numberOrBlank(match.fragment_count),
      String(match.ms_ms_fragments || ''),
      _matchedFragmentDetails(row, match),
      _ms2PeakText(row),
      _numberOrBlank(
        Array.isArray(row.query_peaks) && row.query_peaks.length
          ? row.query_peaks.length
          : row.processed_peaks_count ?? row.peaks_count
      ),
      rawSpectrum,
      _rawMs2PeakCount(row),
      String(unresolved ? (member.classification || '') : (match.classification || '')),
      String(match.display_class || match.classification || ''),
      _numberOrBlank(match.classification_confidence),
      _classificationSourceLabel(
        match.display_classification_source || match.classification_source
      ),
      _classificationStatusLabel(match),
      _classificationReasonLabel(match.classification_reliability_reason),
      String(unresolved ? (member.dataset_source || '') : (match.dataset_source || '')),
      String(unresolved ? (member.lib_spectrum_id || '') : (match.lib_spectrum_id || '')),
      unresolved ? String(member.inchikey || '') : String(match.inchikey || ''),
      unresolved ? String(member.smiles || '') : String(match.smiles || ''),
      String(row.spectrum_id || ''),
      [...new Set(notes)].join('; '),
    ];
  }

  async function _fetchAllFeatureResults(taskId, blankView) {
    const results = [];
    let offset = 0;
    const limit = 5000;
    while (true) {
      const query = new URLSearchParams({
        matched_only: 'false',
        blank_view: blankView,
        limit: String(limit),
        offset: String(offset),
      });
      const response = await fetch(
        `${API.BASE}/tasks/${encodeURIComponent(taskId)}/results?${query}`,
        { cache: 'no-store' }
      );
      if (!response.ok) {
        throw new Error(`MS2 结果读取失败 (HTTP ${response.status})`);
      }
      const payload = await response.json();
      const page = Array.isArray(payload.results) ? payload.results : [];
      results.push(...page);
      if (!payload.has_more || !page.length) break;
      offset += page.length;
    }
    return results;
  }

  async function _ms2CompatTable(taskId, blankView) {
    const [featureResults, compounds] = await Promise.all([
      _fetchAllFeatureResults(taskId, blankView),
      _fetchAllCompoundConclusions(taskId),
    ]);
    const matched = featureResults.filter(row => row?.match);
    const deduplicated = _dedupMs2Results(matched);
    const finalGroups = (compounds || []).filter(group => {
      const representative = group?.representative || {};
      if (blankView === 'clean') {
        return representative.blank_clean_view_excluded !== true;
      }
      if (blankView === 'flagged') {
        return representative.blank_flagged === true;
      }
      return true;
    });
    const rows = [];
    let isomerGroupCount = 0;
    let isomerCandidateCount = 0;
    const resultBySpectrum = new Map(
      featureResults
        .filter(row => String(row?.spectrum_id || ''))
        .map(row => [String(row.spectrum_id), row])
    );
    if (finalGroups.length) {
      finalGroups.forEach((group, groupIndex) => {
        const representative = group.representative || {};
        const source = resultBySpectrum.get(
          String(representative.spectrum_id || '')
        ) || {};
        const enrichedRepresentative = {
          ...source,
          ...representative,
          match: {
            ...(source.match || {}),
            ...(representative.match || {}),
          },
        };
        const groupNumber = `G${String(groupIndex + 1).padStart(3, '0')}`;
        if (_isUnresolvedCompound(group)) {
          const members = _dedupIsomerMembers(group);
          isomerGroupCount += 1;
          isomerCandidateCount += members.length;
          members.forEach(member => {
            rows.push(_ms2ReaderRow(
              enrichedRepresentative, group, member, groupNumber
            ));
          });
        } else {
          rows.push(_ms2ReaderRow(
            enrichedRepresentative, group, null, groupNumber
          ));
        }
      });
    } else {
      deduplicated.forEach((row, groupIndex) => rows.push(
        _ms2ReaderRow(
          row,
          null,
          null,
          `G${String(groupIndex + 1).padStart(3, '0')}`
        )
      ));
    }
    return {
      headers: MS2_READER_HEADERS,
      rows,
      matchedCount: matched.length,
      deduplicatedCount: finalGroups.length || deduplicated.length,
      exportedRowCount: rows.length,
      isomerGroupCount,
      isomerCandidateCount,
      compounds,
    };
  }

  async function _downloadMs2Compat(taskId, blankView) {
    Toast.show('正在生成可直接阅读的 MS2 鉴定结果表…', 'info', 5000);
    const table = await _ms2CompatTable(taskId, blankView);
    _saveCsv(
      `${_selectedSampleStem(taskId)}_ms2_${blankView}.csv`,
      table.headers,
      table.rows
    );
    Toast.show(
      `MS2 结果表已生成：${table.deduplicatedCount} 个鉴定组，` +
      `${table.isomerGroupCount} 个异构体组已展开为 ` +
      `${table.isomerCandidateCount} 条候选记录（共 ${table.exportedRowCount} 行）`,
      'success',
      7000
    );
  }

  async function _targetedCompatTable(taskId) {
    const features = await _fetchAllExpectedFeatures(taskId);
    const grouped = _groupTargetedFeatures(features);
    const headers = [
      '证据等级', '主要候选', '候选数', '靶向数据库', '分子式',
      '代表加合物', '共洗脱加合物', '离子模式', '理论 m/z', '实测 m/z',
      '误差 (ppm)', '主峰 RT (min)', '组分峰面积', '峰质量 (0–10)',
      'OpenMS 峰验证', '背景状态', '原始峰段数', '信噪比',
      '诊断 MS2 谱图数', '实验谱库谱图数',
    ];
    return {
      headers,
      rows: grouped.map(row => row.values),
      featureCount: features.length,
      compoundCount: grouped.length,
    };
  }

  async function _downloadTargetedSummaryCompat(taskId) {
    Toast.show('正在归并靶向峰段并生成化合物汇总…', 'info', 5000);
    const table = await _targetedCompatTable(taskId);
    _saveCsv(
      `${taskId}_targeted_compounds.csv`,
      table.headers,
      table.rows
    );
    Toast.show(
      `靶向汇总已生成：${table.featureCount} 个峰段 → ${table.compoundCount} 个去重靶向化合物`,
      'success',
      7000
    );
  }

  async function _fetchAllCompoundConclusions(taskId) {
    const compounds = [];
    let offset = 0;
    const limit = 2000;
    while (true) {
      const query = new URLSearchParams({
        limit: String(limit),
        offset: String(offset),
        materialize: 'true',
      });
      const response = await fetch(
        `${API.BASE}/tasks/${encodeURIComponent(taskId)}/compounds?${query}`,
        { cache: 'no-store' }
      );
      if (!response.ok) {
        throw new Error(`异构体候选读取失败 (HTTP ${response.status})`);
      }
      const payload = await response.json();
      const page = Array.isArray(payload.compounds) ? payload.compounds : [];
      compounds.push(...page);
      if (!payload.has_more || !page.length) break;
      offset += page.length;
    }
    return compounds;
  }

  function _selectedSampleStem(taskId) {
    const select = document.getElementById('exportTaskSelect');
    const selected = select?.selectedOptions?.[0];
    const selectedFilename = String(selected?.value || '') === String(taskId)
      ? selected?.dataset?.filename
      : '';
    const filename = String(
      selectedFilename
      || _batchTasksById.get(taskId)?.filename
      || taskId
    ).split(/[\\/]/).pop();
    return filename.replace(/\.(?:mzML|mzXML|mgf)$/i, '')
      .replace(/[<>:"/\\|?*\x00-\x1f]+/g, '_')
      .replace(/\s+/g, '_')
      .replace(/^[._-]+|[._-]+$/g, '') || taskId;
  }

  const ISOMER_CANDIDATE_HEADERS = [
    'Feature ID', 'Result ID', '代表谱图编号', '前体 m/z', 'RT (min)', '鉴定层级',
    '异构体组置信度', '异构体解析状态', '异构体歧义类型',
    '组内候选总数', '候选记录状态', '候选排名', '候选角色',
    '是否当前最佳候选', '候选化合物名称', '分子式', 'InChIKey',
    'SMILES', '候选谱库加合物', 'MS/MS 相似度', '谱熵相似度',
    '反向匹配覆盖', '质量误差 (ppm)', '综合匹配分', '候选排序分',
    '候选置信度', '化学类别', '与最佳候选关系', '谱库来源',
    '馆藏谱图编号',
  ];

  function _isomerCandidateCompatTable(compounds) {
    const groups = compounds.filter(compound => {
      const scope = String(compound?.identity_scope || '').toLowerCase();
      return scope === 'isomer_group'
        || scope === 'molecular_formula_isomer_group'
        || scope === 'stereoisomer_group'
        || compound?.conclusion?.isomer_group_only === true;
    });
    const rows = [];
    groups.forEach(group => {
      const representative = group.representative || {};
      const match = representative.match || {};
      const members = _dedupIsomerMembers(group);
      members.forEach((member, index) => {
        const lead = member.is_ranking_lead === true
          || String(member.role || '') === 'ranking_lead'
          || index === 0;
        rows.push([
          group.feature_id || '',
          group.result_id || '',
          representative.spectrum_id || '',
          _numberOrBlank(representative.precursor_mz),
          _numberOrBlank(representative.retention_time),
          group.identity_scope || '',
          group.conclusion?.isomer_group_confidence
            || match.isomer_group_confidence || '',
          match.isomer_resolution_status || '',
          match.isomer_ambiguity_type || '',
          Math.max(
            Number(group.candidate_member_count || 0),
            members.length
          ),
          member.candidate_record_status || 'available',
          member.rank || index + 1,
          lead
            ? '当前最佳候选（排名首位，未代表唯一解析）'
            : '竞争异构体候选',
          lead ? '是' : '否',
          member.compound_name || '',
          member.formula || group.formula || '',
          member.inchikey || '',
          member.smiles || '',
          member.adduct || '',
          _numberOrBlank(member.cosine),
          _numberOrBlank(member.entropy),
          _numberOrBlank(member.reverse_match),
          _numberOrBlank(member.mass_ppm),
          _numberOrBlank(member.final_score),
          _numberOrBlank(member.rank_score),
          member.confidence || '',
          member.classification || '',
          member.relationship || '',
          member.dataset_source || '',
          member.lib_spectrum_id || '',
        ]);
      });
    });
    return {
      headers: ISOMER_CANDIDATE_HEADERS,
      rows,
      groupCount: groups.length,
    };
  }

  async function _downloadIsomerCandidatesCompat(taskId) {
    Toast.show('正在整理完整异构体竞争候选清单…', 'info', 5000);
    const compounds = await _fetchAllCompoundConclusions(taskId);
    const table = _isomerCandidateCompatTable(compounds);
    _saveCsv(
      `${_selectedSampleStem(taskId)}_isomer_competing_candidates.csv`,
      table.headers,
      table.rows
    );
    Toast.show(
      `异构体候选表已生成：${table.groupCount} 个未解析组，${table.rows.length} 个逐候选条目`,
      'success',
      7000
    );
  }

  async function queueSingleExport(format, params = {}, requestedTaskId = null) {
    const taskId = requestedTaskId || getSelectedTaskId();
    if (!taskId) { Toast.show('请先选择任务', 'warning'); return; }
    const savedTask = _frozenTasks.get(taskId) || await API.get(`/tasks/${encodeURIComponent(taskId)}/status`);
    if (savedTask.completed_batch) {
      if (format === 'csv') FrozenTaskUI.download(`${API.BASE}/tasks/${encodeURIComponent(taskId)}/download/csv?${new URLSearchParams(params)}`);
      else if (format === 'html') FrozenTaskUI.download(`${API.BASE}/tasks/${encodeURIComponent(taskId)}/download/html`);
      else if (format === 'json') FrozenTaskUI.download(`${API.BASE}/tasks/${encodeURIComponent(taskId)}/download/native-conclusions.json.xz`);
      else Toast.show('本样本已保存全量HTML、CSV和完整JSON.XZ，请使用上方对应下载按钮。', 'info', 6000);
      return;
    }
    const query = new URLSearchParams(params);
    _singleExportInFlight += 1;
    Toast.show(`正在提交 ${format.toUpperCase()} 导出，随后显示生成进度…`, 'info', 5000);
    try {
      const response = await fetch(`${API.BASE}/tasks/${encodeURIComponent(taskId)}/export/jobs/${format}?${query}`, {method:'POST'});
      const data = await response.json();
      if (!response.ok) throw new Error(typeof data.detail === 'string' ? data.detail : '导出参数错误');
      let job = data;
      localStorage.setItem('mssystem.lastSingleExport', job.id);
      showExportProgress(job);
      while (!['completed', 'failed'].includes(job.status)) {
        await new Promise(resolve => setTimeout(resolve, 2000));
        const status = await fetch(`${API.BASE}/export/jobs/${job.id}`, {cache:'no-store'});
        if (!status.ok) throw new Error('进度连接暂时中断，可通过进度链接继续查看');
        job = await status.json();
        showExportProgress(job);
      }
      if (job.status === 'failed') throw new Error(job.error || '文件生成失败');
      startArtifactDownload(job.id);
      Toast.show('文件已生成，已交给浏览器下载；进度卡中可再次下载', 'success', 7000);
    } catch (error) {
      Toast.show(`导出失败：${error.message}`, 'error', 10000);
    } finally {
      _singleExportInFlight -= 1;
    }
  }

  function startArtifactDownload(jobId) {
    const anchor = document.createElement('a');
    anchor.href = `${API.BASE}/export/jobs/${jobId}/download`;
    anchor.download = '';
    document.body.appendChild(anchor); anchor.click(); anchor.remove();
  }

  // ── Download Functions ────────────────────────────────────────────
  async function downloadCSV(formatType, resultLayer = 'ms2') {
    const taskId = getSelectedTaskId();
    if (!taskId) { Toast.show('请先选择任务', 'warning'); return; }
    if (
      _supportsTargetedFeatureAudit === null &&
      _csvCapabilityPromise
    ) {
      await _csvCapabilityPromise;
    }
    if (
      formatType === 'results' &&
      resultLayer === 'targeted' &&
      _supportsTargetedFeatureAudit === false
    ) {
      try {
        await _downloadTargetedSummaryCompat(taskId);
      } catch (error) {
        Toast.show(`靶向汇总导出失败：${error.message}`, 'error');
      }
      return;
    }
    if (
      formatType === 'results' &&
      resultLayer === 'isomer_candidates' &&
      _supportsIsomerCandidateCsv === false
    ) {
      try {
        await _downloadIsomerCandidatesCompat(taskId);
      } catch (error) {
        Toast.show(`异构体候选导出失败：${error.message}`, 'error');
      }
      return;
    }
    const blankView = document.getElementById('exportBlankView')?.value || 'all';
    if (
      formatType === 'results' &&
      resultLayer === 'ms2' &&
      (
        _supportsCompleteMs2Csv === false
        || _supportsIsomerCandidateCsv === false
      )
    ) {
      try {
        await _downloadMs2Compat(taskId, blankView);
      } catch (error) {
        Toast.show(`MS2 鉴定表导出失败：${error.message}`, 'error');
      }
      return;
    }
    const query = new URLSearchParams({
      format_type: formatType,
      result_layer: resultLayer,
      ms2_schema: 'reader_v3',
      blank_view: blankView,
      download_nonce: String(Date.now()),
    });
    return queueSingleExport('csv', Object.fromEntries(query));
  }

  async function downloadMs2Csv() {
    const taskId = getSelectedTaskId();
    if (!taskId) {
      Toast.show('请先选择任务', 'warning');
      return;
    }
    const button = document.getElementById('ms2CsvPackageButton');
    const previousText = button?.textContent || '';
    if (button) {
      button.disabled = true;
      button.textContent = '正在整理结果 CSV…';
    }
    try {
      Toast.show(
        '正在生成 MS2 鉴定结果表（含全部异构体候选）…',
        'info',
        6000
      );
      await downloadCSV('results', 'ms2');

    } catch (error) {
      Toast.show(`MS2 结果 CSV 导出失败：${error.message}`, 'error');
    } finally {
      if (button) {
        button.disabled = false;
        button.textContent = previousText || '下载 MS2 结果 CSV';
      }
    }
  }

  // Compatibility alias for a browser tab that retained the previous button
  // handler.  The output is now one CSV, not a two-file package.
  const downloadMs2CsvPackage = downloadMs2Csv;

  function downloadAllMs2Events() {
    return queueSingleExport('csv', {
      format_type: 'results', result_layer: 'ms2',
      ms2_schema: 'complete_v2', blank_view: 'all',
    });
  }

  function downloadExcel() { return queueSingleExport('excel'); }
  function downloadRawPeaks() { return queueSingleExport('raw-peaks'); }
  function downloadMultilibraryJson() { return queueSingleExport('multilibrary'); }
  function downloadHTML(formatType, language = null) {
    return queueSingleExport('html', {format_type: formatType, language: interfaceLanguage(language)});
  }
  function downloadPDF(language = null) {
    return queueSingleExport('pdf', {language: interfaceLanguage(language)});
  }

  function openStructureDownload() {
    const taskId = getSelectedTaskId();
    if (!taskId) { Toast.show('请先选择任务', 'warning'); return; }
    const index = prompt('输入要下载的匹配序号（从1开始，对应结果表格中的 # 列）:');
    if (!index) return;
    const i = parseInt(index) - 1;
    if (i < 0 || isNaN(i)) { Toast.show('无效序号', 'warning'); return; }
    const format = confirm('点击"确定"下载 PNG 格式，点击"取消"下载 SVG 格式') ? 'png' : 'svg';
    return queueSingleExport('structure', {spectrum_index: i, image_format: format});
  }

  // ── Batch Export ──────────────────────────────────────────────────
  let _groupBy = 'project'; // 'project' | 'date' | 'flat'

  /** 从文件名提取项目/样本分组键。
   *  兼容多种命名约定:
   *    F编号:      20260528174445-F004021.mzML  → "F004021"
   *    短F编号:    F002_NEG_C1_2.mzML / F002-NEG-3-7.mzML → "F002"  (F+1~6位数字, 下划线/连字符均兼容)
   *    时间戳+样本: 20260520223333-DHJS1.mzML    → "DHJS"   (取样本前缀)
   *    纯样本名:   QHSS1.mzML / BK_NEG1.mzML     → "QHSS" / "BK"  (取开头字母, 忽略极性/重复后缀)
   *  规则: 去扩展名 → 去时间戳前缀 → 优先 F编号(开头 F+1~6位数字; 6位长编号兼容出现在中间),
   *        否则取开头字母样本前缀。
   */
  function _extractProject(filename) {
    if (!filename) return '其他';
    let name = filename.replace(/\.[^.]+$/, '');          // 去扩展名
    name = name.replace(/^\d{8,14}[-_]/, '');              // 去时间戳前缀 (YYYYMMDD[-HHMMSS]-)
    const fcode = name.match(/^(F\d{1,6})/i);              // F编号(F+1~6位数字) → 整个编号作为项目
    if (fcode) return fcode[1].toUpperCase();
    const fcodeMid = name.match(/(F\d{6})/i);              // 兼容: 6位长F编号出现在名称中间
    if (fcodeMid) return fcodeMid[1].toUpperCase();
    const prefix = name.match(/^([A-Za-z]+)/);             // 否则取开头字母 (样本/项目前缀, 忽略 NEG/POS/数字后缀)
    if (prefix) return prefix[1].toUpperCase();
    return name || '其他';
  }

  /** 从文件名提取月份 (YYYY-MM) */
  function _extractMonth(filename) {
    const m = filename.match(/(\d{4})(\d{2})/);
    if (m) return `${m[1]}-${m[2]}`;
    return '其他';
  }

  /** 对任务列表分组 */
  function _groupTasks(tasks, groupBy) {
    if (groupBy === 'flat') {
      return { '全部任务': tasks };
    }
    const groups = {};
    for (const t of tasks) {
      let key;
      if (groupBy === 'project') {
        key = _extractProject(t.filename);
      } else {
        key = _extractMonth(t.filename);
      }
      if (!groups[key]) groups[key] = [];
      groups[key].push(t);
    }
    // 按组内文件数降序排
    const sorted = {};
    for (const k of Object.keys(groups).sort((a, b) => groups[b].length - groups[a].length)) {
      sorted[k] = groups[k];
    }
    return sorted;
  }

  async function loadBatchList() {
    try {
      const data = await API.get('/tasks?status=completed');
      _batchTasksById = new Map((data.tasks || []).map(t => [t.id, t]));
      for (const id of Array.from(_selectedBatchTasks)) {
        if (!_batchTasksById.has(id)) _selectedBatchTasks.delete(id);
      }
      const container = document.getElementById('batchTaskList');
      if (!container) return;
      if (!data.tasks?.length) {
        container.innerHTML = '<p class="text-muted" style="font-size:0.82rem;">暂无已完成任务</p>';
        return;
      }

      // 分组
      const groups = _groupTasks(data.tasks, _groupBy);
      const groupKeys = Object.keys(groups);

      let html = '';

      // 分组切换按钮
      html += `<div class="batch-group-bar">
        <span style="font-size:0.78rem;color:#666;margin-right:6px;">分组:</span>
        <button class="batch-grp-btn${_groupBy === 'project' ? ' active' : ''}" onclick="ExportMod.setGroupBy('project')">按项目</button>
        <button class="batch-grp-btn${_groupBy === 'date' ? ' active' : ''}" onclick="ExportMod.setGroupBy('date')">按月份</button>
        <button class="batch-grp-btn${_groupBy === 'flat' ? ' active' : ''}" onclick="ExportMod.setGroupBy('flat')">不分组</button>
      </div>`;

      // 全选/全不选
      html += `<div class="batch-select-bar">
        <button class="batch-grp-btn" onclick="ExportMod.selectAll(true)">全选</button>
        <button class="batch-grp-btn" onclick="ExportMod.selectAll(false)">全不选</button>
        <span id="batchSelectedCount" style="font-size:0.78rem;color:#666;margin-left:8px;">已选 0 项</span>
      </div>`;

      // 渲染各组
      for (const key of groupKeys) {
        const tasks = groups[key];
        const isSingle = groupKeys.length === 1;

        // 分组头 (可折叠, 默认折叠除非是唯一组)
        html += `<div class="batch-group${isSingle ? '' : ' collapsed'}">`;
        if (!isSingle) {
          html += `<div class="batch-group-header" onclick="this.parentElement.classList.toggle('collapsed')">
            <span class="batch-group-arrow">▼</span>
            <span class="batch-group-name">${Tasks.escapeHtml(key)}</span>
            <span class="batch-group-count">${tasks.length} 个</span>
            <button class="batch-grp-btn" style="margin-left:auto;font-size:0.72rem;padding:1px 6px;" onclick="event.stopPropagation();ExportMod.selectGroup('${Tasks.escapeHtml(key)}', true)">全选</button>
            <button class="batch-grp-btn" style="font-size:0.72rem;padding:1px 6px;" onclick="event.stopPropagation();ExportMod.selectGroup('${Tasks.escapeHtml(key)}', false)">取消</button>
          </div>`;
        }
        html += `<div class="batch-group-items">`;
        for (const t of tasks) {
          html += `<label class="batch-item" data-group="${Tasks.escapeHtml(key)}">
            <input type="checkbox" value="${t.id}"
              ${_selectedBatchTasks.has(t.id) ? 'checked' : ''}
              onchange="ExportMod.toggleBatch('${t.id}', this.checked)">
            ${Tasks.escapeHtml(t.filename)} (${t.matched_spectra}/${t.total_spectra} 匹配)
          </label>`;
        }
        html += `</div></div>`;
      }

      container.innerHTML = html;
      _updateSelectedCount();
    } catch (e) {
      console.error('Failed to load batch list:', e);
    }
  }

  function setGroupBy(mode) {
    _groupBy = mode;
    loadBatchList();
  }

  function toggleBatch(taskId, checked) {
    if (checked) _selectedBatchTasks.add(taskId);
    else _selectedBatchTasks.delete(taskId);
    _updateSelectedCount();
  }

  function selectAll(checked) {
    document.querySelectorAll('#batchTaskList .batch-item input[type="checkbox"]').forEach(cb => {
      cb.checked = checked;
      if (checked) _selectedBatchTasks.add(cb.value);
      else _selectedBatchTasks.delete(cb.value);
    });
    _updateSelectedCount();
  }

  function selectGroup(groupName, checked) {
    document.querySelectorAll(`#batchTaskList .batch-item[data-group="${groupName}"] input[type="checkbox"]`).forEach(cb => {
      cb.checked = checked;
      if (checked) _selectedBatchTasks.add(cb.value);
      else _selectedBatchTasks.delete(cb.value);
    });
    _updateSelectedCount();
  }

  function _updateSelectedCount() {
    const n = _selectedBatchTasks.size;
    const projects = [...new Set(Array.from(_selectedBatchTasks)
      .map(id => _batchTasksById.get(id))
      .filter(Boolean)
      .map(t => _extractProject(t.filename)))].sort();
    const projectText = projects.length === 1
      ? projects[0]
      : projects.length > 1 ? `${projects.length} 个项目` : '';
    const summary = `已选 ${n} 项${projectText ? ` · ${projectText}` : ''}`;
    // 更新内联计数
    const el = document.getElementById('batchSelectedCount');
    if (el) el.textContent = summary;
    // 更新悬浮条
    const floatBar = document.getElementById('batchFloatBar');
    const floatCount = document.getElementById('batchFloatCount');
    if (floatBar) floatBar.classList.toggle('hidden', n === 0);
    if (floatCount) floatCount.textContent = summary;
  }

  function _taskFilename(taskId) {
    const task = _batchTasksById.get(taskId);
    if (task?.filename) return String(task.filename);
    const selected = document.getElementById(
      'exportTaskSelect'
    )?.selectedOptions?.[0];
    if (String(selected?.value || '') === String(taskId)) {
      return String(selected?.dataset?.filename || taskId);
    }
    return String(taskId);
  }

  function _uniqueCompatArchiveName(name, taskId, usedNames) {
    let candidate = name;
    const lower = () => candidate.toLocaleLowerCase();
    if (!usedNames.has(lower())) {
      usedNames.add(lower());
      return candidate;
    }
    const dot = name.lastIndexOf('.');
    const suffix = `__${String(taskId).slice(0, 8)}`;
    candidate = dot > 0
      ? `${name.slice(0, dot)}${suffix}${name.slice(dot)}`
      : `${name}${suffix}`;
    let serial = 2;
    while (usedNames.has(lower())) {
      candidate = dot > 0
        ? `${name.slice(0, dot)}${suffix}_${serial}${name.slice(dot)}`
        : `${name}${suffix}_${serial}`;
      serial += 1;
    }
    usedNames.add(lower());
    return candidate;
  }

  async function _fetchCsvLayerText(taskId, resultLayer, optional = false) {
    const query = new URLSearchParams({
      format_type: 'results',
      result_layer: resultLayer,
      ms2_schema: 'reader_v3',
      blank_view: 'all',
      download_nonce: String(Date.now()),
    });
    const response = await fetch(
      `${API.BASE}/tasks/${encodeURIComponent(taskId)}/download/csv?${query}`,
      { cache: 'no-store' }
    );
    if (optional && response.status === 404) return '';
    if (!response.ok) {
      throw new Error(
        `${resultLayer} CSV 读取失败 (HTTP ${response.status})`
      );
    }
    return response.text();
  }

  async function _buildCompatCsvBatch(
    format, taskIds, onProgress = null
  ) {
    if (
      typeof SampleReport === 'undefined'
      || typeof SampleReport.buildZipEntries !== 'function'
    ) {
      throw new Error('完整 CSV 打包模块未加载，请刷新页面后重试');
    }
    const taskRows = taskIds.map(taskId => ({
      taskId,
      filename: _taskFilename(taskId),
    }));
    const projects = [...new Set(
      taskRows.map(row => _extractProject(row.filename))
    )];
    const multiProject = projects.length > 1;
    const entries = [];
    const manifestRows = [[
      'Project', 'Source_File', 'Task_ID', 'Output_File', 'Format',
      'Matched_Spectra', 'Unique_Compounds', 'Exported_At',
    ]];
    const usedNames = new Set();
    const failed = [];
    let exported = 0;
    const exportedAt = new Date().toISOString();

    for (let index = 0; index < taskRows.length; index += 1) {
      const { taskId, filename } = taskRows[index];
      const project = _extractProject(filename);
      const folder = multiProject ? `${project}/` : '';
      const stem = _selectedSampleStem(taskId);
      const outputNames = [];
      let matchedCount = '';
      let uniqueCount = '';
      try {
        let ms2Table = null;
        let compounds = null;
        if (format === 'csv' || format === 'csv_ms2') {
          ms2Table = await _ms2CompatTable(taskId, 'all');
          compounds = ms2Table.compounds;
          matchedCount = ms2Table.matchedCount;
          uniqueCount = ms2Table.deduplicatedCount;
          const ms2Name = _uniqueCompatArchiveName(
            `${folder}${stem}_${format === 'csv_ms2' ? 'ms2_identifications' : 'ms2'}.csv`,
            taskId,
            usedNames
          );
          entries.push({
            name: ms2Name,
            content: _csvText(ms2Table.headers, ms2Table.rows),
          });
          outputNames.push(ms2Name);
        }

        if (format === 'csv_isomer_candidates') {
          compounds = compounds || await _fetchAllCompoundConclusions(taskId);
          const candidates = _isomerCandidateCompatTable(compounds);
          const candidateName = _uniqueCompatArchiveName(
            `${folder}${stem}_isomer_competing_candidates.csv`,
            taskId,
            usedNames
          );
          entries.push({
            name: candidateName,
            content: _csvText(candidates.headers, candidates.rows),
          });
          outputNames.push(candidateName);
        }

        if (format === 'csv') {
          const targeted = await _targetedCompatTable(taskId);
          const targetedName = _uniqueCompatArchiveName(
            `${folder}${stem}_targeted_ms1.csv`,
            taskId,
            usedNames
          );
          entries.push({
            name: targetedName,
            content: _csvText(targeted.headers, targeted.rows),
          });
          outputNames.push(targetedName);

          const auditName = _uniqueCompatArchiveName(
            `${folder}${stem}_targeted_peak_audit.csv`,
            taskId,
            usedNames
          );
          entries.push({
            name: auditName,
            content: await _fetchCsvLayerText(
              taskId, 'targeted_features'
            ),
          });
          outputNames.push(auditName);

          const referenceQc = await _fetchCsvLayerText(
            taskId, 'reference_qc', true
          );
          if (referenceQc) {
            const qcName = _uniqueCompatArchiveName(
              `${folder}${stem}_reference_standard_qc.csv`,
              taskId,
              usedNames
            );
            entries.push({ name: qcName, content: referenceQc });
            outputNames.push(qcName);
          }
        }

        exported += 1;
        manifestRows.push([
          project,
          filename,
          taskId,
          outputNames.join(' | '),
          format === 'csv_ms2'
            ? 'CSV MS2 鉴定结果表（含异构体候选）'
            : format === 'csv_isomer_candidates'
            ? 'CSV 异构体竞争候选明细'
            : 'CSV 完整数据包',
          matchedCount,
          uniqueCount,
          exportedAt,
        ]);
      } catch (error) {
        failed.push({ taskId, filename, error: String(error?.message || error) });
      }
      onProgress?.({
        phase: 'csv',
        completed: index + 1,
        total: taskRows.length,
        exported,
        failed: failed.length,
      });
    }
    if (!exported) {
      throw new Error('所选样本均未能生成 CSV 数据包');
    }
    failed.forEach(item => manifestRows.push([
      '', item.filename, item.taskId, '', `导出失败：${item.error}`,
      '', '', exportedAt,
    ]));
    entries.unshift({
      name: '00_export_manifest.csv',
      content: _csvText(manifestRows[0], manifestRows.slice(1)),
    });
    const blob = await SampleReport.buildZipEntries(entries, onProgress);
    return {
      blob,
      exported,
      failed: failed.length,
      downloadName: _fallbackBatchFilename(format, taskIds),
    };
  }

  function _fallbackBatchFilename(format, taskIds) {
    const samples = [...new Set(taskIds
      .map(id => _taskFilename(id))
      .map(filename => String(filename || 'sample')
        .replace(/\.[^.]+$/, '')
        .replace(/[<>:"/\\|?*\x00-\x1f]/g, '_')
        .replace(/\s+/g, '_')))].filter(Boolean);
    let samplePart = 'MSsystem';
    if (samples.length <= 4 && samples.length) samplePart = samples.join('__');
    else if (samples.length > 4) {
      samplePart = `${samples[0]}__${samples[1]}__等${samples.length}个样本`;
    }
    samplePart = samplePart.slice(0, 170).replace(/[<>:"/\\|?*\x00-\x1f]/g, '_');
    const token = ({
      csv: 'CSV',
      csv_ms2: 'CSV_MS2',
      csv_isomer_candidates: 'CSV_ISOMER_CANDIDATES',
      csv_targeted: 'CSV_TARGETED',
      csv_targeted_features: 'CSV_TARGETED_AUDIT',
      csv_audit: 'CSV_MS2_AUDIT',
      excel: 'EXCEL_LAYERED',
      report: 'HTML_DETAILED',
      html_standard: 'HTML_STANDARD',
      html_detailed: 'HTML_DETAILED',
      html_integrated: 'HTML_INTEGRATED',
      html_targeted: 'HTML_TARGETED',
      pdf: 'PDF',
    })[format] || format.toUpperCase();
    const now = new Date();
    const stamp = `${now.getFullYear()}${String(now.getMonth() + 1).padStart(2, '0')}${String(now.getDate()).padStart(2, '0')}`;
    const countToken = taskIds.length > 1 ? `_${taskIds.length}samples` : '';
    return `${samplePart}_${token}${countToken}_${stamp}.zip`;
  }

  function _responseDownloadName(resp, fallback) {
    const disposition = resp.headers.get('Content-Disposition') || '';
    const utf8 = disposition.match(/filename\*\s*=\s*UTF-8''([^;]+)/i);
    if (utf8) {
      try { return decodeURIComponent(utf8[1].trim().replace(/^"|"$/g, '')); }
      catch { /* fall through to the ASCII filename */ }
    }
    const quoted = disposition.match(/filename\s*=\s*"([^"]+)"/i);
    if (quoted) return quoted[1];
    const plain = disposition.match(/filename\s*=\s*([^;\s]+)/i);
    return plain ? plain[1] : fallback;
  }

  async function batchExport(format = 'csv', language = null) {
    language = interfaceLanguage(language);
    if (_batchExportInFlight) {
      Toast.show('批量导出正在生成，请等待当前下载完成', 'info');
      return;
    }
    if (_selectedBatchTasks.size === 0) {
      Toast.show('请至少选择一个任务', 'warning');
      return;
    }
    const frozen = Array.from(_selectedBatchTasks).map(id => _batchTasksById.get(id)).filter(task => task?.completed_batch);
    if (frozen.length) {
      if (frozen.length !== _selectedBatchTasks.size) {
        Toast.show('请将宣痹汤完整批次与其他任务分别导出', 'warning');
        return;
      }
      FrozenTaskUI.download(`${API.BASE}/completed-batches/${encodeURIComponent(frozen[0].completed_batch.id)}/native-download-all`);
      Toast.show('开始下载完整宣痹汤批次：149份样本的全部结果表与报告；单文件CSV可在上方选择样本下载', 'info', 8000);
      return;
    }
    const languageLabel = language === 'en' && (
      ['report', 'pdf'].includes(format) || format.startsWith('html_')
    ) ? '英文' : '';
    const formatNames = {
      csv: 'CSV 完整数据包（每样本：MS2 鉴定结果表 + 靶向化合物汇总 + 原始峰审计）',
      csv_ms2: 'CSV MS2 鉴定结果表（含异构体候选；每样本 1 份）',
      csv_isomer_candidates: 'CSV 异构体竞争候选明细',
      csv_targeted: 'CSV 靶向化合物汇总',
      csv_targeted_features: 'CSV 靶向离子 / 原始峰审计（逐峰）',
      csv_audit: 'CSV MS2 详细审计',
      excel: 'Excel 分层工作簿',
      report: `${languageLabel} HTML 详细综合报告`,
      html_standard: `${languageLabel} HTML 标准综合报告`,
      html_detailed: `${languageLabel} HTML 详细综合报告`,
      html_integrated: 'HTML 新版样本综合报告',
      html_targeted: `${languageLabel} HTML 靶向库专用报告`,
      pdf: `${languageLabel} PDF 综合报告`,
    };
    const exportButtons = [
      document.getElementById('batchCsvBundleButton'),
      document.getElementById('batchExportButton'),
    ].filter(Boolean);
    const originalButtonStates = exportButtons.map(button => ({
      button,
      disabled: button.disabled,
      text: button.textContent,
    }));
    _batchExportInFlight = true;
    try {
      const taskIds = Array.from(_selectedBatchTasks);
      Toast.show(`正在生成 ${formatNames[format] || format}，请稍候...`, 'info', 8000);
      exportButtons.forEach(button => {
        button.disabled = true;
        button.textContent = `正在打包 ${taskIds.length} 个样本…`;
      });
      if (format === 'html_integrated') {
        if (typeof SampleReport === 'undefined' || typeof SampleReport.generateBatch !== 'function') {
          throw new Error('新版样本综合报告模块未加载，请刷新页面后重试');
        }
        const outcome = await SampleReport.generateBatch(taskIds, progress => {
          exportButtons.forEach(button => {
            button.textContent = progress.phase === 'zip'
              ? `正在压缩 ${progress.completed}/${progress.total}…`
              : `正在生成报告 ${progress.completed}/${progress.total}…`;
          });
        });
        const failedText = outcome.failed ? `，${outcome.failed} 个失败已写入清单` : '';
        Toast.show(`已打包 ${outcome.exported} 个样本综合报告${failedText}：${outcome.downloadName}`, outcome.failed ? 'warning' : 'success', 7000);
        return;
      }
      // Submit quickly; generation continues even when the browser disconnects.
      const requestBody = { format, taskIds, language };
      const resumeKey = 'mssystem.pendingBatchExport';
      let saved = null;
      try { saved = JSON.parse(localStorage.getItem(resumeKey) || 'null'); } catch {}
      let job;
      if (saved && JSON.stringify(saved.request) === JSON.stringify(requestBody)) {
        const existing = await fetch(`${API.BASE}/export/jobs/${saved.id}`, { cache: 'no-store' });
        if (existing.ok) job = await existing.json();
        if (job?.status === 'failed' || job?.stale) job = null;
      }
      if (!job) {
        const started = await fetch(`${API.BASE}/export/jobs/${encodeURIComponent(format)}?language=${encodeURIComponent(language)}`, {
          method: 'POST', headers: { 'Content-Type': 'application/json' }, body: JSON.stringify(taskIds),
        });
        if (!started.ok) {
          const error = await started.json();
          throw new Error(error.detail || `HTTP ${started.status}`);
        }
        job = await started.json();
        localStorage.setItem(resumeKey, JSON.stringify({ id: job.id, request: requestBody }));
      }
      let progressLink = document.getElementById('batchExportProgressLink');
      if (!progressLink && exportButtons[0]) {
        progressLink = document.createElement('a');
        progressLink.id = 'batchExportProgressLink';
        progressLink.target = '_blank';
        progressLink.rel = 'noopener';
        progressLink.style.marginLeft = '12px';
        exportButtons[0].insertAdjacentElement('afterend', progressLink);
      }
      if (progressLink) {
        progressLink.href = `/static/export-progress.html?job=${job.id}`;
        progressLink.textContent = '查看进度 / 下载已生成文件';
      }
      while (job.status !== 'completed' && job.status !== 'failed') {
        showExportProgress(job);
        exportButtons.forEach(button => {
          button.textContent = job.status === 'queued' ? '导出排队中…'
            : `已处理 ${job.completed}/${job.total}，失败 ${job.failed || 0}…`;
          button.title = job.current_task ? `当前任务：${job.current_task}` : '';
        });
        await new Promise(resolve => setTimeout(resolve, 2000));
        const statusResponse = await fetch(`${API.BASE}/export/jobs/${job.id}`, { cache: 'no-store' });
        if (!statusResponse.ok) throw new Error('无法获取导出进度；再次点击可恢复查询');
        job = await statusResponse.json();
      }
      showExportProgress(job);
      if (job.status === 'failed') {
        localStorage.removeItem(resumeKey);
        throw new Error(job.error || '导出生成失败');
      }
      startArtifactDownload(job.id);
      Toast.show(`文件已生成：${job.exported || 0} 个样本成功，${job.failed || 0} 个失败。已交给浏览器下载，可在进度卡中再次下载。`, job.failed ? 'warning' : 'success', 7000);

    } catch (e) {
      Toast.show(`批量导出失败: ${e.message}`, 'error');
    } finally {
      _batchExportInFlight = false;
      originalButtonStates.forEach(({ button, disabled, text }) => {
        button.disabled = disabled;
        button.textContent = text;
      });
    }
  }

  function batchExportSelected() {
    const selector = document.getElementById('batchExportType');
    const [format, configuredLanguage = 'auto'] = String(
      selector?.value || 'csv_ms2|auto'
    ).split('|');
    batchExport(format, interfaceLanguage(configuredLanguage));
  }

  return {
    queueSingleExport,
    init, refreshSelect, updateOptions, filterExportSelect,
    downloadCSV, downloadMs2Csv, downloadMs2CsvPackage, downloadAllMs2Events,
    downloadExcel, downloadRawPeaks, downloadMultilibraryJson, downloadHTML, downloadPDF,
    openStructureDownload,
    loadBatchList, toggleBatch, batchExport, batchExportSelected,
    setGroupBy, selectAll, selectGroup,
    get selectedBatchTasks() { return _selectedBatchTasks; },
  };
})();
