/**
 * MSsystem — Task Management Module
 * Task list, search, filter, progress, actions (delete/retry/view)
 */
const Tasks = (() => {
  let _searchQuery = '';
  let _statusFilter = '';
  let _pollingTimers = {};
  let _selected = new Set();        // 批量操作选中的 task_id 集合
  let _renderedTasks = [];          // 当前渲染的任务列表 (用于全选)
  let _recycleVisible = false;
  let _recycleOffset = 0;
  const _recycleLimit = 50;
  let _recyclePayload = null;
  let _recyclePurgePreview = null;

  // ── Initialization ────────────────────────────────────────────────
  function init() {
    // SSE listener for real-time updates
    API.onTaskUpdate(onSSEUpdate);
  }

  // ── SSE Handler ───────────────────────────────────────────────────
  function onSSEUpdate(event) {
    if (event.type === 'queue_reordered') {
      loadTasks();
      return;
    }
    if (event.type === 'progress' || event.type === 'status_change') {
      updateTaskRow(event.task_id, event);
      const terminal = ['completed', 'failed', 'cancelled'].includes(event.status);
      // Status transitions also change which actions are safe for the row.
      if (event.type === 'status_change' && !terminal) loadTasks();
      if (terminal) {
        loadTasks();
        Results.refreshSelect();
        ExportMod.refreshSelect();
        if (typeof populateTaskSelects === 'function') populateTaskSelects();
        if (event.status === 'completed') {
          Toast.show(`分析完成: ${event.matched_spectra || 0}/${event.total_spectra || 0} 匹配`, 'success');
          sendBrowserNotification('MSsystem 分析完成',
            `${event.filename || '任务'}: ${event.matched_spectra || 0}/${event.total_spectra || 0} 匹配`);
        } else if (event.status === 'failed') {
          Toast.show(`分析失败: ${event.error || '未知错误'}`, 'error');
        } else {
          Toast.show(`分析已安全停止: ${event.filename || event.task_id}`, 'warning');
        }
      }
      updateRunningBadge();
    }
  }

  // ── Task List ─────────────────────────────────────────────────────
  async function loadTasks() {
    try {
      let path = '/tasks';
      const params = [];
      if (_statusFilter) params.push(`status=${_statusFilter}`);
      if (params.length) path += '?' + params.join('&');

      const data = await API.get(path);
      // If search query, filter client-side
      let tasks = data.tasks || [];
      if (_searchQuery) {
        const q = _searchQuery.toLowerCase();
        tasks = tasks.filter(t =>
          t.id.toLowerCase().includes(q) ||
          t.filename.toLowerCase().includes(q) ||
          (t.completed_batch && ('宣痹汤 xbt ' + t.completed_batch.id).includes(q)) ||
          (t.compound_names || []).some(c => c.toLowerCase().includes(q))
        );
      }
      renderTaskTable(tasks);
      updateRunningBadge(tasks);
    } catch (e) {
      console.error('Failed to load tasks:', e);
    }
  }

  function renderTaskTable(tasks) {
    const area = document.getElementById('taskTableArea');
    if (!area) return;

    if (tasks.length === 0) {
      area.innerHTML = `
        <div class="empty-state">
          <div class="empty-icon">📋</div>
          <div class="empty-title">暂无分析任务</div>
          <div class="empty-desc">上传文件并开始分析后，任务会显示在这里</div>
        </div>`;
      return;
    }

    area.innerHTML = `
      <div class="table-wrap">
      <table>
        <thead><tr>
          <th style="width:34px;text-align:center;"><input type="checkbox" id="taskSelectAll" onchange="Tasks.toggleSelectAll(this)"></th>
          <th>任务 ID</th><th>文件名</th><th>离子模式</th><th>上传／批处理时间</th><th>状态</th>
          <th>进度</th><th>匹配数</th><th>匹配率</th><th>分析参数</th><th>操作</th>
        </tr></thead>
        <tbody>
          ${tasks.map(t => renderTaskRow(t)).join('')}
        </tbody>
      </table>
      </div>`;
    _renderedTasks = tasks;
    _syncSelectionUI();

    // Resume polling for live tasks (fallback if SSE fails)
    for (const t of tasks) {
      if (['running', 'cancelling'].includes(t.status) && !_pollingTimers[t.id]) {
        startPolling(t.id);
      }
    }
  }

  function renderTaskRow(t) {
    const usesFinalCount = t.compound_count_semantics
      === 'sample_supported_main_results_excluding_unresolved_events';
    const unresolvedCount = Number(t.unresolved_precursor_event_count || 0);
    const uniqueInfo = t.unique_compounds == null ? '' : (usesFinalCount
      ? `<br><span style="font-size:0.68rem;color:var(--text-secondary);font-weight:650;">${t.unique_compounds} 个主要化合物结果</span>${unresolvedCount
        ? `<br><span style="font-size:0.68rem;color:#9a3412;">${unresolvedCount} 个未解析事件（未计数）</span>` : ''}`
      : `<br><span style="font-size:0.68rem;color:var(--text-secondary);">${t.unique_compounds} 个谱库候选身份</span>`);
    const nativeInfo = t.completed_batch
      ? `<br><span style="font-size:0.68rem;">${Number(t.completed_batch.summary.compound_result_records || 0).toLocaleString()} 条原生结论 · ${Number(t.completed_batch.summary.ms1_feature_records || 0).toLocaleString()} 个MS1特征</span>` : '';
    const blankMismatchCount = Number(
      t.blank_summary?.status_counts?.blank_polarity_not_applicable || 0
    );
    const blankNotApplicable = t.blank_summary?.applicable === false
      || (t.total_spectra > 0 && blankMismatchCount >= t.total_spectra);
    const blankMismatchTitle = t.blank_summary?.applicable === false
      ? (t.blank_summary.warning || '样本与空白极性不匹配')
      : '样本与空白极性不匹配；空白过滤未执行，“0”不代表已完成扣除';
    const blankInfo = t.completed_batch
      ? `<br><span style="font-size:0.68rem;" title="保留所有原始MS/MS事件及其配对空白证据">🧪 配对空白：标记 ${t.completed_batch.summary.blank_flagged_ms2 || 0} · 净化视图排除 ${t.completed_batch.summary.blank_clean_view_excluded_ms2 || 0}</span>`
      : t.blank_summary?.requested
      ? (blankNotApplicable
        ? `<br><span style="font-size:0.68rem;color:#b91c1c;font-weight:700;" title="${escapeHtml(blankMismatchTitle)}">🧪 空白未执行（极性不匹配）</span>`
        : (t.blank_summary?.blank_operation === 'comparison_only'
          || t.blank_summary?.filter_policy === 'soft_evidence'
          ? `<br><span style="font-size:0.68rem;color:#92400e;" title="匹配空白特征 ${t.blank_summary.feature_matches || 0} 条；未跨针相减 MS/MS">🧪 空白仅比较 · ${t.blank_summary.excluded_spectra || 0} 条被扣除</span>`
          : `<br><span style="font-size:0.68rem;color:#92400e;" title="仅在方法可比性与排除授权均已确认时允许">🧪 已确认背景排除 ${t.blank_summary.excluded_spectra || 0} 条</span>`))
      : '';
    const assignedBlankCount = Number(
      t.blank_summary?.assigned_blank_replicate_count || 0
    );
    const skippedBlankCount = Number(
      t.blank_summary?.skipped_blank_replicate_count || 0
    );
    const assignedBlankPolarity = Array.isArray(
      t.blank_summary?.assigned_blank_polarities
    ) ? t.blank_summary.assigned_blank_polarities.join('/') : '';
    const blankPairingConfirmed = t.blank_summary?.method_comparability_verified === true
      && t.blank_summary?.background_exclusion_authorized === true;
    const blankAssignmentBoundary = blankPairingConfirmed
      ? '（方法可比性与背景排除已确认）'
      : '（不等于方法配对确认）';
    const blankAssignmentInfo = (
      t.blank_summary?.blank_assignment_mode
      && t.blank_summary?.applicable !== false
    )
      ? `<br><span style="font-size:0.68rem;color:#276749;" title="${escapeHtml(
          `使用：${(t.blank_summary.assigned_blank_filenames || []).join('、') || '极性未知空白'}`
          + `${skippedBlankCount ? `；跳过：${(t.blank_summary.skipped_blank_filenames || []).join('、')}` : ''}`
        )}">↔ 按极性分配 ${assignedBlankPolarity || '未知极性'} 空白 ${assignedBlankCount} 个${blankAssignmentBoundary}${skippedBlankCount ? `，跳过异极性 ${skippedBlankCount} 个` : ''}</span>`
      : '';
    return `
      <tr id="task-row-${t.id}" class="${_selected.has(t.id) ? 'is-selected' : ''}">
        <td style="text-align:center;"><input type="checkbox" class="task-row-check" data-id="${t.id}" onchange="Tasks.toggleSelect('${t.id}', this.checked)" ${_selected.has(t.id) ? 'checked' : ''}></td>
        <td class="mono" style="font-size:0.7rem;">${t.id}</td>
        <td data-i18n-skip style="font-weight:500;max-width:180px;overflow:hidden;text-overflow:ellipsis;white-space:nowrap;"
            title="${escapeHtml(t.filename)}">${escapeHtml(t.filename)}</td>
        <td>${ionModeBadge(t)}</td>
        <td style="font-size:0.75rem;color:var(--text-secondary);">${t.completed_batch ? '批处理开始<br>' : ''}${formatTime(t.upload_time)}</td>
        <td class="task-status">${statusBadge(t.status)}</td>
        <td class="task-progress">${progressHtml(t)}</td>
        <td class="task-match" title="匹配谱图 / 总谱图">${t.total_spectra ? `${t.matched_spectra}/${t.total_spectra}${uniqueInfo}${nativeInfo}${blankInfo}${blankAssignmentInfo}` : '—'}</td>
        <td>${t.total_spectra ? `${t.match_rate}%` : '—'}</td>
        <td>${paramsHtml(t)}</td>
        <td class="task-actions">${renderTaskActions(t)}</td>
      </tr>`;
  }

  function renderTaskActions(t) {
    if (t.completed_batch) return `<div class="btn-group"><button class="btn btn-secondary btn-sm" onclick="Tasks.viewResults('${t.id}')">🔍 查看</button><a class="btn btn-secondary btn-sm" href="/api/tasks/${encodeURIComponent(t.id)}/report" target="_blank" rel="noopener">本文件全量报告</a></div>`;
    const canDelete = ['pending', 'queued', 'completed', 'failed', 'cancelled'].includes(t.status);
    const canRetry = ['completed', 'failed', 'cancelled'].includes(t.status);
    const deleteTitle = t.status === 'queued' ? '移出分析队列并放入回收站' : '将任务及相关数据移入回收站';
    const queuePosition = Number(t.queue_position || 0);
    const queueLength = Number(t.queue_length || 0);
    return `
      <div class="btn-group">
        ${t.status === 'pending' ? `<button class="btn btn-primary btn-sm" onclick="Tasks.startSingle('${t.id}')">▶ 开始</button>` : ''}
        ${canRetry ? `<button class="btn btn-secondary btn-sm" onclick="Tasks.retry('${t.id}')">🔄 重试</button>` : ''}
        ${t.status === 'completed' ? `<button class="btn btn-secondary btn-sm" onclick="Tasks.viewResults('${t.id}')">🔍 查看</button>` : ''}
        ${t.status === 'running' ? `<button class="btn btn-danger btn-sm" onclick="Tasks.cancelTask('${t.id}')" title="安全停止分析">⏹ 停止</button>` : ''}
        ${t.status === 'cancelling' ? '<button class="btn btn-secondary btn-sm" disabled title="正在等待安全停止点">⏳ 停止中</button>' : ''}
        ${t.status === 'finalizing' ? '<button class="btn btn-secondary btn-sm" disabled title="正在生成最终排序和结果文件，暂不可停止或删除">⚙️ 最终处理中</button>' : ''}
        ${t.queue_reorderable ? `<button class="btn btn-secondary btn-sm" onclick="Tasks.reorderTask('${t.id}', 'front')" title="移到等待队列第一位" ${queuePosition <= 1 ? 'disabled' : ''}>⏫ 置顶</button>` : ''}
        ${t.queue_reorderable ? `<button class="btn btn-secondary btn-sm btn-icon" onclick="Tasks.reorderTask('${t.id}', 'up')" title="在等待队列中上移一位" ${queuePosition <= 1 ? 'disabled' : ''}>↑</button>` : ''}
        ${t.queue_reorderable ? `<button class="btn btn-secondary btn-sm btn-icon" onclick="Tasks.reorderTask('${t.id}', 'down')" title="在等待队列中下移一位" ${queuePosition >= queueLength ? 'disabled' : ''}>↓</button>` : ''}
        ${canDelete ? `<button class="btn btn-danger btn-sm btn-icon" onclick="Tasks.deleteTask('${t.id}')" title="${deleteTitle}">🗑</button>` : ''}
      </div>`;
  }

  function progressHtml(t) {
    if (t.status === 'running') {
      return `
        <div style="font-size:0.72rem;color:var(--text-secondary);">${t.progress}%</div>
        <div class="progress-wrap"><div class="progress-bar running" style="width:${t.progress}%"></div></div>
        <div class="progress-text">${escapeHtml(t.message || '')}</div>`;
    } else if (t.status === 'completed') {
      return `
        <div style="font-size:0.72rem;color:var(--text-secondary);">100%</div>
        <div class="progress-wrap"><div class="progress-bar" style="width:100%"></div></div>`;
    } else if (t.status === 'queued') {
      const position = Number(t.queue_position || 0);
      return `<div style="font-size:0.72rem;color:var(--text-muted);">排队中${position ? ` · 第 ${position} 位` : ''}</div>`;
    } else if (t.status === 'cancelling') {
      return `<div style="font-size:0.72rem;color:var(--danger);">正在安全停止…</div>`;
    } else if (t.status === 'cancelled') {
      return `<div style="font-size:0.72rem;color:var(--text-muted);">已停止</div>`;
    } else if (t.status === 'finalizing') {
      return `
        <div style="font-size:0.72rem;color:var(--text-secondary);">${Number(t.progress || 96)}%</div>
        <div class="progress-wrap"><div class="progress-bar running" style="width:${Number(t.progress || 96)}%"></div></div>
        <div class="progress-text">${escapeHtml(t.message || '正在生成最终排序和结果文件…')}</div>`;
    } else if (t.status === 'failed') {
      return `<div style="font-size:0.72rem;color:var(--danger);">失败</div>`;
    }
    return `<div style="font-size:0.72rem;color:var(--text-muted);">等待中</div>`;
  }

  // ── Inline Row Update ─────────────────────────────────────────────
  function updateTaskRow(taskId, data) {
    const row = document.getElementById(`task-row-${taskId}`);
    if (!row) return;

    const statusCell = row.querySelector('.task-status');
    const progressCell = row.querySelector('.task-progress');
    const matchCell = row.querySelector('.task-match');

    if (statusCell && data.status) statusCell.innerHTML = statusBadge(data.status);
    if (progressCell) {
      const t = { status: data.status, progress: data.progress || 0, message: data.message || '' };
      progressCell.innerHTML = progressHtml(t);
    }
    if (matchCell && data.matched_spectra !== undefined) {
      matchCell.textContent = `${data.matched_spectra}/${data.total_spectra || 0}`;
    }
  }

  // ── Task Actions ──────────────────────────────────────────────────
  async function startSingle(taskId) {
    try {
      const params = Upload.analysisParams('ionMode');
      await API.post(`/tasks/${taskId}/analyze`, params);
      startPolling(taskId);
      loadTasks();
      Toast.show('分析已启动', 'success');
    } catch (e) {
      Toast.show(`启动失败: ${e.message}`, 'error');
    }
  }

  async function retry(taskId) {
    if (!confirm('确认重新分析此任务？之前的结果将被覆盖。')) return;
    try {
      const task = _renderedTasks.find(t => t.id === taskId);
      const params = Object.assign({
        db_mode: document.getElementById('dbMode')?.value || 'full',
        cosine_threshold: 0.3, final_score_threshold: 0.25,
        mz_tolerance_da: 0.02, mz_tolerance_ppm: 10.0,
        top_n_peaks: 15, strict_evidence_gate: false,
        enable_chimeric_deconvolution: true,
        enable_ion_identity_networking: true,
      }, (task && task.params) || {});
      if (params.strict_evidence_gate == null) params.strict_evidence_gate = false;
      await API.post(`/tasks/${taskId}/retry`, params);
      startPolling(taskId);
      loadTasks();
      Toast.show('重新分析已启动', 'success');
    } catch (e) {
      Toast.show(`重新分析失败: ${e.message}`, 'error');
    }
  }

  async function deleteTask(taskId) {
    const task = _renderedTasks.find(t => t.id === taskId);
    try {
      await API.del(`/tasks/${taskId}`);
      stopPolling(taskId);
      _selected.delete(taskId);
      await loadTasks();
      if (_recycleVisible) await loadRecycleBin();
      Results.refreshSelect();
      ExportMod.refreshSelect();
      if (typeof populateTaskSelects === 'function') populateTaskSelects();
      Toast.show(task?.status === 'queued' ? '样本已移出队列并放入回收站' : '任务已移入回收站', 'success');
    } catch (e) {
      Toast.show(`删除失败: ${e.message}`, 'error');
    }
  }

  async function cancelTask(taskId) {
    const task = _renderedTasks.find(t => t.id === taskId);
    const label = task?.filename || taskId;
    if (!confirm(`确认安全停止“${label}”的分析？停止完成后可删除该任务。`)) return;
    try {
      const res = await API.post(`/tasks/${taskId}/cancel`, {});
      startPolling(taskId);
      await loadTasks();
      if (_recycleVisible) await loadRecycleBin();
      Toast.show(
        res.status === 'cancelled' ? '任务已停止' : '已请求安全停止，完成后即可删除',
        'warning'
      );
    } catch (e) {
      Toast.show(`停止失败: ${e.message}`, 'error');
    }
  }

  async function reorderTask(taskId, action) {
    const labels = { front: '置顶', up: '上移', down: '下移' };
    try {
      const result = await API.post('/tasks/queue/reorder', {
        task_id: taskId,
        action,
      });
      await loadTasks();
      const position = (result.queued || []).find(item => item.id === taskId)?.position;
      Toast.show(
        `队列已${labels[action] || '调整'}${position ? `，当前第 ${position} 位` : ''}`,
        'success'
      );
    } catch (e) {
      await loadTasks();
      Toast.show(`调整队列失败: ${e.message}`, 'error');
    }
  }

  // ── Batch Selection / Delete ─────────────────────────────────────
  function toggleSelect(taskId, checked) {
    if (checked) _selected.add(taskId);
    else _selected.delete(taskId);
    document.getElementById(`task-row-${taskId}`)?.classList.toggle('is-selected', checked);
    _syncSelectionUI();
  }

  function toggleSelectAll(checkbox) {
    if (checkbox.checked) {
      _renderedTasks.forEach(t => _selected.add(t.id));
    } else {
      _renderedTasks.forEach(t => _selected.delete(t.id));
    }
    // 更新每行的勾选状态
    document.querySelectorAll('.task-row-check').forEach(cb => { cb.checked = checkbox.checked; });
    document.querySelectorAll('#taskTableArea tbody tr').forEach(row => {
      row.classList.toggle('is-selected', checkbox.checked);
    });
    _syncSelectionUI();
  }

  function _syncSelectionUI() {
    const btn = document.getElementById('btnBatchDelete');
    const runBtn = document.getElementById('btnBatchRun');
    const consensusBtn = document.getElementById('btnCrossOccurrence');
    const selAll = document.getElementById('taskSelectAll');
    if (selAll && _renderedTasks.length) {
      selAll.checked = _renderedTasks.every(t => _selected.has(t.id));
      selAll.indeterminate = _selected.size > 0 && !selAll.checked;
    }
    if (btn) {
      const n = _selected.size;
      btn.disabled = n === 0;
      btn.textContent = n > 0 ? `🗑 批量移除 (${n})` : '🗑 批量移除';
    }
    if (runBtn) {
      const runnable = Array.from(_selected).filter(id => {
        const task = _renderedTasks.find(t => t.id === id);
        return task && ['pending', 'completed', 'failed', 'cancelled'].includes(task.status);
      }).length;
      runBtn.disabled = runnable === 0;
      runBtn.textContent = runnable > 0 ? `▶ 批量运行 (${runnable})` : '▶ 批量运行';
    }
    if (consensusBtn) {
      const completed = Array.from(_selected).filter(id => {
        const task = _renderedTasks.find(t => t.id === id);
        return task && task.status === 'completed';
      }).length;
      consensusBtn.disabled = completed < 3;
      consensusBtn.textContent = completed >= 3
        ? crossOccurrenceCopy(
          `⇄ 跨进样共识排序 (${completed})`,
          `⇄ Cross-run Consensus Ranking (${completed})`
        )
        : crossOccurrenceCopy('⇄ 跨进样共识排序', '⇄ Cross-run Consensus Ranking');
    }
  }

  function currentAnalysisParams() {
    return Upload.analysisParamsForRun('batchIonMode');
  }

  function localizedCopy(zh, en) {
    return globalThis.I18n?.getLanguage?.() === 'en' ? en : zh;
  }

  function blankRunLabel(params) {
    const blankFiles = Array.isArray(params?.blank_files)
      ? params.blank_files : [];
    if (!blankFiles.length) {
      return localizedCopy('空白：未启用', 'Blanks: not enabled');
    }
    const roleCounts = Object.values(params.blank_roles || {}).reduce(
      (counts, role) => {
        const key = String(role || 'unspecified');
        counts[key] = (counts[key] || 0) + 1;
        return counts;
      },
      {}
    );
    const roleLabels = {
      background_blank: localizedCopy('背景空白', 'background blank'),
      solvent_blank: localizedCopy('溶剂空白', 'solvent blank'),
      method_blank: localizedCopy('方法空白', 'method blank'),
      extraction_blank: localizedCopy('提取空白', 'extraction blank'),
      carryover_blank: localizedCopy('残留空白', 'carryover blank'),
      unspecified: localizedCopy('类型未注明', 'unspecified'),
    };
    const roles = Object.entries(roleCounts).map(([role, count]) =>
      `${roleLabels[role] || role} ${count}`
    ).join(localizedCopy('、', ', '));
    const operation = params.blank_filter_policy === 'hard_legacy'
      ? localizedCopy('执行特征级背景排除', 'feature-level background exclusion')
      : localizedCopy('仅作背景比较', 'comparison only');
    return localizedCopy(
      `空白：${blankFiles.length} 个（${roles}；${operation}）`,
      `Blanks: ${blankFiles.length} (${roles}; ${operation})`
    );
  }

  async function batchRun() {
    const ids = Array.from(_selected).filter(id => {
      const task = _renderedTasks.find(t => t.id === id);
      return task && ['pending', 'completed', 'failed', 'cancelled'].includes(task.status);
    });
    if (!ids.length) return;
    let params;
    try {
      params = currentAnalysisParams();
    } catch (error) {
      Toast.show(
        localizedCopy(
          `批量运行参数无效：${error.message}`,
          `Invalid batch-run parameters: ${error.message}`
        ),
        'warning',
        9000
      );
      return;
    }
    const mode = params.ion_mode;
    const modeLabel = mode === 'POS' ? '正离子' : mode === 'NEG' ? '负离子' : '自动识别极性';
    const hasExisting = ids.some(id => {
      const task = _renderedTasks.find(t => t.id === id);
      return task && (task.status === 'completed' || task.status === 'failed');
    });
    const overwrite = hasExisting ? '；已完成/失败任务会重新分析并覆盖旧结果' : '';
    const blankLabel = blankRunLabel(params);
    const confirmation = localizedCopy(
      `确认批量运行 ${ids.length} 个样本？离子模式：${modeLabel}；${blankLabel}${overwrite}。`,
      `Run ${ids.length} samples? Ion mode: ${mode}; ${blankLabel}${hasExisting ? '; completed/failed runs will be overwritten' : ''}.`
    );
    if (!confirm(confirmation)) return;

    const runBtn = document.getElementById('btnBatchRun');
    if (runBtn) { runBtn.disabled = true; runBtn.innerHTML = '<span class="spinner spinner-sm"></span> 正在排队...'; }
    try {
      const res = await API.post('/tasks/batch/analyze', {
        task_ids: ids,
        params,
      });
      (res.queued || []).forEach(startPolling);
      _selected.clear();
      await loadTasks();
      const failed = (res.failed || []).length;
      const skipped = (res.skipped || []).length;
      Toast.show(
        `已将 ${res.queued_count || 0} 个样本加入队列${skipped ? `，跳过 ${skipped} 个` : ''}${failed ? `，失败 ${failed} 个` : ''}`,
        failed ? 'warning' : 'success'
      );
    } catch (e) {
      Toast.show(`批量运行失败: ${e.message}`, 'error');
    } finally {
      _syncSelectionUI();
    }
  }

  function crossOccurrenceCopy(zh, en) {
    return window.MSSLanguage?.getLanguage?.() === 'en' ? en : zh;
  }

  async function crossOccurrenceRank() {
    const ids = Array.from(_selected).filter(id => {
      const task = _renderedTasks.find(t => t.id === id);
      return task && task.status === 'completed';
    });
    if (ids.length < 3) {
      Toast.show(
        crossOccurrenceCopy(
          '请至少勾选 3 个同一 LC 方法的已完成独立进样',
          'Select at least three completed independent injections from one LC method.'
        ),
        'warning'
      );
      return;
    }
    const confirmed = confirm(crossOccurrenceCopy(
      `确认对 ${ids.length} 个已完成任务执行跨进样共识排序？\n\n`
        + '只有当这些任务来自同一 LC 方法、RT 可直接比较且每个任务代表独立进样时才能继续。'
        + '系统仅在原有同分子式 Top-5 内换序；不会把扩展候选升为首位，也不会把结果标为标准品确证。',
      `Apply cross-run consensus ranking to ${ids.length} completed runs?\n\n`
        + 'Continue only if all runs use one LC method, their retention times are directly comparable, '
        + 'and every run is an independent injection. The system only reorders the existing same-formula '
        + 'Top-5; expanded candidates cannot become the lead and no result is marked as reference-standard confirmed.'
    ));
    if (!confirmed) return;
    const button = document.getElementById('btnCrossOccurrence');
    if (button) {
      button.disabled = true;
      button.innerHTML = `<span class="spinner spinner-sm"></span> ${crossOccurrenceCopy(
        '正在计算共识…', 'Computing consensus…'
      )}`;
    }
    try {
      const result = await API.post('/workflows/cross-occurrence-ranking', {
        task_ids: ids,
        same_lc_method_confirmed: true,
        activate: true,
      });
      await loadTasks();
      Results.refreshSelect();
      if (typeof SampleOverview !== 'undefined') SampleOverview.refreshSelect();
      if (typeof SampleReport !== 'undefined') SampleReport.refreshSelect();
      const changed = Number(result.applied_count || 0);
      if (result.status === 'activated') {
        Toast.show(
          crossOccurrenceCopy(
            `V16 跨进样共识已启用：${changed} 个前体事件换序；回滚编号 ${result.cohort_id}`,
            `V16 cross-run consensus activated: ${changed} precursor events reordered; rollback ID ${result.cohort_id}`
          ),
          'success',
          12000
        );
      } else {
        Toast.show(
          crossOccurrenceCopy(
            '所选任务没有满足严格条件的候选换序，原排序保持不变',
            'No candidate reorder met the strict criteria; the original ranking was retained.'
          ),
          'warning',
          9000
        );
      }
    } catch (error) {
      Toast.show(
        crossOccurrenceCopy(
          `跨进样共识排序失败：${error.message}`,
          `Cross-run consensus ranking failed: ${error.message}`
        ),
        'error',
        12000
      );
    } finally {
      _syncSelectionUI();
    }
  }

  async function batchDelete() {
    const ids = Array.from(_selected);
    if (ids.length === 0) return;
    try {
      // Clicking “批量移除” is the user's confirmation.  Keep the backend's
      // short-lived, state-bound preview token so a changed task set still
      // fails closed, but do not interrupt deletion with another modal.
      const preview = await API.post('/tasks/batch/delete/preview', {
        task_ids: ids,
      });
      const res = await API.post('/tasks/batch/delete', {
        task_ids: ids,
        confirmation_token: preview.confirmation_token,
        confirmation_phrase: preview.confirmation_phrase,
      });
      (res.deleted || []).forEach(id => stopPolling(id));
      (res.stop_requested || []).forEach(id => startPolling(id));
      _selected.clear();
      await loadTasks();
      if (_recycleVisible) await loadRecycleBin();
      Results.refreshSelect();
      if (typeof SampleOverview !== 'undefined') SampleOverview.refreshSelect();
      if (typeof SampleReport !== 'undefined') SampleReport.refreshSelect();
      ExportMod.refreshSelect();
      if (typeof populateTaskSelects === 'function') populateTaskSelects();
      const ok = res.deleted_count || 0;
      const stopping = res.stop_requested_count || 0;
      const fail = (res.failed || []).length;
      Toast.show(
        `已将 ${ok} 个任务移入回收站${stopping ? `，已请求安全停止 ${stopping} 个（停止后可再次移除）` : ''}${fail ? `，${fail} 个未处理` : ''}`,
        fail || stopping ? 'warning' : 'success'
      );
    } catch (e) {
      Toast.show(`批量移除失败: ${e.message}`, 'error');
    }
  }

  function viewResults(taskId) {
    if (typeof SampleOverview !== 'undefined' && SampleOverview.openTask) {
      SampleOverview.openTask(taskId);
    } else {
      App.navigateTo('results');
    }
  }

  async function toggleRecycleBin() {
    const area = document.getElementById('taskRecycleBinArea');
    const button = document.getElementById('btnTaskRecycleBin');
    if (!area) return;
    _recycleVisible = !_recycleVisible;
    area.hidden = !_recycleVisible;
    if (button) {
      button.textContent = _recycleVisible
        ? uiText('♻ 收起回收站') : uiText('♻ 回收站');
    }
    if (_recycleVisible) {
      _recyclePurgePreview = null;
      await loadRecycleBin(0);
    }
  }

  async function loadRecycleBin(offset = _recycleOffset) {
    const area = document.getElementById('taskRecycleBinArea');
    if (!area || !_recycleVisible) return;
    _recycleOffset = Math.max(0, Number(offset) || 0);
    area.innerHTML = '<div class="empty-state"><span class="spinner spinner-sm"></span> 正在读取回收站…</div>';
    try {
      const payload = await API.get(
        `/task-recycle-bin?limit=${_recycleLimit}&offset=${_recycleOffset}`
      );
      _recyclePayload = payload;
      renderRecycleBin(payload);
    } catch (error) {
      area.innerHTML = `<div class="empty-state"><div class="empty-title">回收站读取失败</div><div class="empty-desc">${escapeHtml(error.message)}</div></div>`;
    }
  }

  function renderRecycleBin(payload) {
    const area = document.getElementById('taskRecycleBinArea');
    if (!area) return;
    const entries = payload?.entries || [];
    const summary = payload?.summary || {};
    const total = Number(payload?.total || 0);
    const totalSize = Number(summary.size_bytes || 0);
    const totalFiles = Number(summary.file_count || 0);
    const blocked = Number(summary.purge_blocked_count || 0);
    const currentStart = total ? Number(payload.offset || 0) + 1 : 0;
    const currentEnd = Math.min(
      Number(payload.offset || 0) + entries.length,
      total
    );
    const blockedMessage = uiText(
      '有 {count} 个项目未通过安全检查；清空功能已锁定，请查看审计记录。',
      { count: blocked }
    );
    const pageMessage = uiText(
      '显示 {start}–{end}，共 {total} 个项目',
      { start: currentStart, end: currentEnd, total }
    );
    const confirmation = _recyclePurgePreview
      ? `
        <div class="recycle-purge-confirm" role="alert">
          <div>
            <strong>${escapeHtml(uiText('永久清空 {count} 个项目？', {
              count: _recyclePurgePreview.task_count || 0,
            }))}</strong>
            <p>${escapeHtml(uiText('将释放 {size}；原始文件和分析结果将无法恢复。', {
              size: formatBytes(_recyclePurgePreview.size_bytes || 0),
            }))}</p>
          </div>
          <div class="btn-group">
            <button type="button" class="btn btn-secondary btn-sm" onclick="Tasks.cancelEmptyRecycleBin()">取消</button>
            <button type="button" class="btn btn-danger btn-sm" id="btnConfirmRecyclePurge" onclick="Tasks.confirmEmptyRecycleBin()">永久清空</button>
          </div>
        </div>` : '';

    area.innerHTML = `
      <section class="recycle-panel" aria-label="任务回收站">
        <div class="recycle-panel-header">
          <div class="recycle-title-wrap">
            <span class="recycle-icon" aria-hidden="true">♻</span>
            <div>
              <h3>任务回收站</h3>
              <p>查看、恢复或永久清理已移除的分析项目。</p>
            </div>
          </div>
          <div class="btn-group">
            <button type="button" class="btn btn-secondary btn-sm" onclick="Tasks.loadRecycleBin()">🔄 刷新</button>
            <button type="button" class="btn btn-danger btn-sm" onclick="Tasks.beginEmptyRecycleBin()"
              ${!total || blocked ? 'disabled' : ''} title="永久删除回收站内全部项目">🗑 清空回收站</button>
          </div>
        </div>
        <div class="recycle-summary-grid">
          <div class="recycle-summary-card"><span>可恢复项目</span><strong>${escapeHtml(String(total))}</strong></div>
          <div class="recycle-summary-card"><span>占用空间</span><strong>${escapeHtml(formatBytes(totalSize))}</strong></div>
          <div class="recycle-summary-card"><span>归档文件</span><strong>${escapeHtml(formatCount(totalFiles))}</strong></div>
          <div class="recycle-summary-card"><span>到期自动清理</span><strong>${payload.automatic_permanent_deletion ? `保留 ${escapeHtml(String(payload.retention_days))} 天` : '未启用'}</strong><button type="button" class="btn btn-secondary btn-sm" onclick="MSsystemStorage.open()">管理</button></div>
        </div>
        ${blocked ? `<div class="recycle-warning">${escapeHtml(blockedMessage)}</div>` : ''}
        ${confirmation}
        ${!entries.length ? `
          <div class="recycle-empty">
            <span aria-hidden="true">✓</span>
            <div><strong>回收站为空</strong><p>当前没有可恢复的已删除项目。</p></div>
          </div>` : `
          <div class="table-wrap recycle-table-wrap"><table class="recycle-table">
            <thead><tr><th>样本 / 项目</th><th>任务 ID</th><th>移除时间</th><th>删除前状态</th><th>空间 / 文件</th><th>归档状态</th><th>操作</th></tr></thead>
            <tbody>${entries.map(item => `
              <tr>
                <td data-scientific-value title="${escapeHtml(item.filename || '')}"><strong>${escapeHtml(item.filename || item.task_id || '')}</strong></td>
                <td data-scientific-value><code>${escapeHtml(item.task_id || '')}</code></td>
                <td>${escapeHtml(formatTime(item.deleted_at))}</td>
                <td>${statusBadge(item.previous_status || 'unknown')}</td>
                <td><strong>${escapeHtml(formatBytes(item.size_bytes || 0))}</strong><br><span class="recycle-muted">${escapeHtml(uiText('{count} 个文件', { count: formatCount(item.file_count || 0) }))}</span></td>
                <td>${item.recoverable && item.safe_to_purge
                  ? '<span class="recycle-state-ok">可恢复</span>'
                  : '<span class="recycle-state-warning">需审计</span>'}</td>
                <td><button type="button" class="btn btn-secondary btn-sm task-restore-button"
                  data-deletion-id="${escapeHtml(item.deletion_id || '')}"
                  data-task-id="${escapeHtml(item.task_id || '')}">↩ 恢复</button></td>
              </tr>`).join('')}</tbody>
          </table></div>
          <div class="recycle-pagination">
            <span>${escapeHtml(pageMessage)}</span>
            <div class="btn-group">
              <button type="button" class="btn btn-secondary btn-sm" onclick="Tasks.loadRecycleBin(${Math.max(0, Number(payload.offset || 0) - _recycleLimit)})" ${Number(payload.offset || 0) <= 0 ? 'disabled' : ''}>上一页</button>
              <button type="button" class="btn btn-secondary btn-sm" onclick="Tasks.loadRecycleBin(${Number(payload.offset || 0) + _recycleLimit})" ${payload.has_more ? '' : 'disabled'}>下一页</button>
            </div>
          </div>`}
      </section>`;
    area.querySelectorAll('.task-restore-button').forEach(button => {
      button.addEventListener('click', () => restoreTask(
        button.dataset.deletionId || '', button.dataset.taskId || ''
      ));
    });
    window.MSSLanguage?.translate?.(area);
  }

  async function beginEmptyRecycleBin() {
    try {
      const preview = await API.post('/task-recycle-bin/purge/preview', {});
      if (!Number(preview.task_count || 0)) {
        _recyclePurgePreview = null;
        await loadRecycleBin(0);
        return;
      }
      _recyclePurgePreview = preview;
      renderRecycleBin(_recyclePayload || { entries: [], total: 0 });
    } catch (error) {
      Toast.show(uiText('无法清空回收站: {error}', {
        error: error.message,
      }), 'error');
      await loadRecycleBin();
    }
  }

  function cancelEmptyRecycleBin() {
    _recyclePurgePreview = null;
    if (_recyclePayload) renderRecycleBin(_recyclePayload);
  }

  async function confirmEmptyRecycleBin() {
    const preview = _recyclePurgePreview;
    if (!preview?.confirmation_token) return;
    const button = document.getElementById('btnConfirmRecyclePurge');
    if (button) {
      button.disabled = true;
      button.textContent = uiText('正在永久清理…');
    }
    try {
      const result = await API.del('/task-recycle-bin', {
        confirmation_token: preview.confirmation_token,
      });
      _recyclePurgePreview = null;
      await loadRecycleBin(0);
      Toast.show(
        uiText('已永久清理 {count} 个项目，释放 {size}', {
          count: result.removed_task_count || 0,
          size: formatBytes(result.removed_size_bytes || 0),
        }),
        result.failure_count ? 'warning' : 'success'
      );
    } catch (error) {
      _recyclePurgePreview = null;
      Toast.show(uiText('永久清理失败: {error}', {
        error: error.message,
      }), 'error');
      await loadRecycleBin();
    }
  }

  async function restoreTask(deletionId, taskId) {
    if (!deletionId || !taskId) return;
    if (!confirm(uiText('确认恢复任务 {task_id} 及其原始数据？排队中的任务会恢复为“待处理”，不会自动运行。', { task_id: taskId }))) return;
    try {
      await API.post(`/task-recycle-bin/${encodeURIComponent(deletionId)}/${encodeURIComponent(taskId)}/restore`, {});
      await loadTasks();
      await loadRecycleBin();
      Results.refreshSelect();
      ExportMod.refreshSelect();
      Toast.show(uiText('任务及其原始数据已恢复'), 'success');
    } catch (error) {
      Toast.show(uiText('恢复失败: {error}', {
        error: error.message,
      }), 'error');
    }
  }

  // ── Search & Filter ───────────────────────────────────────────────
  function setSearch(query) {
    _searchQuery = query;
    loadTasks();
  }

  function setStatusFilter(status) {
    _statusFilter = status;
    loadTasks();
  }

  // ── Polling (fallback for SSE) ────────────────────────────────────
  function startPolling(taskId) {
    if (_pollingTimers[taskId]) clearInterval(_pollingTimers[taskId]);
    _pollingTimers[taskId] = setInterval(async () => {
      try {
        const status = await API.get(`/tasks/${taskId}/status`);
        updateTaskRow(taskId, status);
        if (['completed', 'failed', 'cancelled'].includes(status.status)) {
          stopPolling(taskId);
          loadTasks();
          Results.refreshSelect();
          ExportMod.refreshSelect();
          if (typeof populateTaskSelects === 'function') populateTaskSelects();
          if (status.status === 'completed') {
            Toast.show(`分析完成: ${status.matched_spectra}/${status.total_spectra} 匹配`, 'success');
          }
        }
      } catch (e) {
        stopPolling(taskId);
      }
    }, 2000);
  }

  function stopPolling(taskId) {
    if (_pollingTimers[taskId]) {
      clearInterval(_pollingTimers[taskId]);
      delete _pollingTimers[taskId];
    }
  }

  // ── Badge ─────────────────────────────────────────────────────────
  function updateRunningBadge(tasks) {
    const badge = document.getElementById('taskCount');
    if (!badge) return;
    // If we don't have tasks list, fetch count
    if (!tasks) {
      API.get('/tasks?status=running').then(d => {
        const c = d.tasks?.length || 0;
        badge.style.display = c > 0 ? 'inline-flex' : 'none';
        badge.textContent = c;
      }).catch(() => {});
      return;
    }
    const runningCount = tasks.filter(t => t.status === 'running').length;
    badge.style.display = runningCount > 0 ? 'inline-flex' : 'none';
    badge.textContent = runningCount;
  }

  // ── Helpers ───────────────────────────────────────────────────────
  function uiText(source, params = null) {
    if (window.MSSLanguage?.t) return window.MSSLanguage.t(source, params);
    return String(source).replace(/\{([A-Za-z0-9_]+)\}/g, (match, key) => (
      params && Object.prototype.hasOwnProperty.call(params, key)
        ? String(params[key]) : match
    ));
  }

  function formatCount(value) {
    const numeric = Number(value || 0);
    if (window.MSSLanguage?.formatNumber) {
      return window.MSSLanguage.formatNumber(numeric);
    }
    return numeric.toLocaleString();
  }

  function formatBytes(value) {
    let amount = Math.max(0, Number(value || 0));
    const units = ['B', 'KB', 'MB', 'GB', 'TB'];
    let unit = 0;
    while (amount >= 1024 && unit < units.length - 1) {
      amount /= 1024;
      unit += 1;
    }
    const digits = unit === 0 ? 0 : amount >= 100 ? 0 : amount >= 10 ? 1 : 2;
    const locale = window.MSSLanguage?.getLocale?.() || undefined;
    return `${amount.toLocaleString(locale, {
      minimumFractionDigits: digits,
      maximumFractionDigits: digits,
    })} ${units[unit]}`;
  }

  function statusBadge(status) {
    const map = {
      pending: '<span class="badge badge-pending"><span class="status-dot pending"></span>待处理</span>',
      queued: '<span class="badge badge-pending"><span class="status-dot pending"></span>排队中</span>',
      running: '<span class="badge badge-running"><span class="status-dot running"></span>运行中</span>',
      finalizing: '<span class="badge badge-running"><span class="status-dot running"></span>写入中</span>',
      cancelling: '<span class="badge badge-running"><span class="status-dot running"></span>停止中</span>',
      cancelled: '<span class="badge badge-failed"><span class="status-dot failed"></span>已停止</span>',
      completed: '<span class="badge badge-completed"><span class="status-dot completed"></span>已完成</span>',
      failed: '<span class="badge badge-failed"><span class="status-dot failed"></span>失败</span>',
    };
    return map[status] || status;
  }

  function ionModeBadge(task) {
    const mode = task.ion_mode || task.params?.ion_mode;
    const counts = task.ion_mode_counts || {};
    if (mode === 'POS') return `<span class="polarity-badge pos" title="${counts.POS || 0} 条正离子谱图">⊕ POS</span>`;
    if (mode === 'NEG') return `<span class="polarity-badge neg" title="${counts.NEG || 0} 条负离子谱图">⊖ NEG</span>`;
    if (mode === 'MIXED') return `<span class="polarity-badge mixed" title="POS ${counts.POS || 0} / NEG ${counts.NEG || 0}">◐ MIXED</span>`;
    if (mode === 'auto') return '<span class="polarity-badge pending">待识别</span>';
    return '<span class="polarity-badge unknown">未记录</span>';
  }

  function paramsHtml(t) {
    const p = (t && t.params) || {};
    if (!Object.keys(p).length) return '<span style="color:var(--text-muted);">—</span>';
    const db = p.db_mode === 'full' ? 'Full' : p.db_mode === 'small' ? 'Small' : (p.db_mode || '—');
    const cos = p.cosine_threshold != null ? Number(p.cosine_threshold).toFixed(2) : '—';
    const fin = p.final_score_threshold != null ? Number(p.final_score_threshold).toFixed(2) : '—';
    const mzDa = p.mz_tolerance_da != null ? p.mz_tolerance_da : '—';
    const mzPpm = p.mz_tolerance_ppm != null ? p.mz_tolerance_ppm : '—';
    const topN = p.top_n_peaks != null ? p.top_n_peaks : '—';
    const blankCount = Array.isArray(p.blank_files) && p.blank_files.length
      ? p.blank_files.length : (p.blank_file ? 1 : 0);
    const blank = blankCount
      ? `\n空白特征过滤: ${blankCount} 个空白, 样品/空白≥${p.blank_sample_ratio || 3}, RT±${p.blank_rt_tolerance_minutes || .3} min`
      : '';
    const gate = p.strict_evidence_gate === true
      ? '严格门槛'
      : p.strict_evidence_gate === false ? '常规门槛' : '历史模式（未记录）';
    const deconv = p.enable_chimeric_deconvolution === false ? '关' : '开';
    const iimn = p.enable_ion_identity_networking === false ? '关' : '开';
    const isolation = p.isolation_window_width_da != null
      ? `\n隔离窗全宽: ${p.isolation_window_width_da} Da（显式覆盖缺失元数据）` : '';
    const compact = `${db} · ${cos}/${fin} · ${gate}`;
    const full = `分析参数:\n库: ${db}\nCosine 阈值: ${cos}\nFinal 阈值: ${fin}\n证据模式: ${gate}\nm/z 容差: ±${mzDa} Da (${mzPpm} ppm)\nTopN 峰: ${topN}\n嵌合去卷积: ${deconv}\nIIMN: ${iimn}${isolation}${blank}`;
    return `<span style="font-size:0.7rem;color:var(--text-secondary);cursor:help;white-space:nowrap;" title="${escapeHtml(full)}">${escapeHtml(compact)}</span>`;
  }

  function formatTime(isoStr) {
    if (!isoStr) return '—';
    try {
      const d = new Date(isoStr);
      const locale = window.MSSLanguage?.getLocale?.() || undefined;
      return d.toLocaleString(locale, {
        year: 'numeric', month: '2-digit', day: '2-digit',
        hour: '2-digit', minute: '2-digit',
      });
    } catch { return isoStr; }
  }

  function escapeHtml(str) {
    const div = document.createElement('div');
    div.textContent = str || '';
    return div.innerHTML;
  }

  function sendBrowserNotification(title, body) {
    if (!('Notification' in window)) return;
    if (Notification.permission === 'granted') {
      new Notification(title, { body, icon: '🔬' });
    } else if (Notification.permission !== 'denied') {
      Notification.requestPermission().then(perm => {
        if (perm === 'granted') new Notification(title, { body, icon: '🔬' });
      });
    }
  }

  return {
    init, loadTasks, setSearch, setStatusFilter,
    startSingle, retry, deleteTask, cancelTask, reorderTask, viewResults,
    toggleSelect, toggleSelectAll, batchRun, crossOccurrenceRank, batchDelete,
    toggleRecycleBin, loadRecycleBin, restoreTask,
    beginEmptyRecycleBin, cancelEmptyRecycleBin, confirmEmptyRecycleBin,
    updateRunningBadge, formatTime, statusBadge, escapeHtml,
  };
})();
