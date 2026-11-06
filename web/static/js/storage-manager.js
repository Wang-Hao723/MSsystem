/* Recoverable tasks and reproducible export scratch have separate policies. */
(function () {
  'use strict';
  let dialog, state, preview, policyDraft, busy = false, previousFocus;
  const esc = value => String(value ?? '').replace(/[&<>"']/g, ch => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[ch]));
  const bytes = value => {
    let number = Number(value) || 0, unit = 0;
    const units = ['B', 'KB', 'MB', 'GB', 'TB'];
    while (number >= 1024 && unit < units.length - 1) { number /= 1024; unit++; }
    return `${number.toLocaleString('zh-CN', {maximumFractionDigits:unit ? 1 : 0})} ${units[unit]}`;
  };
  const time = value => value ? new Date(value).toLocaleString('zh-CN', {hour12:false}) : '尚无清理记录';
  async function request(path = '', options = {}) {
    const response = await fetch(`/api/storage-management${path}`, {
      cache:'no-store', ...options, headers:{'Content-Type':'application/json', ...(options.headers || {})}
    });
    let payload;
    try { payload = await response.json(); } catch (_) { payload = {}; }
    if (!response.ok) {
      if (response.status === 404) throw new Error('后台尚未加载存储管理模块。请保存工作，通过“完全退出并停止服务”退出后重新打开 MSsystem。');
      throw new Error(typeof payload.detail === 'string' ? payload.detail : `操作未完成（${response.status}），请刷新后重试。`);
    }
    return payload;
  }
  function notice(text, kind = '') {
    const element = dialog.querySelector('[data-notice]');
    element.textContent = text;
    element.dataset.kind = kind;
  }
  function setBusy(value) {
    busy = value;
    dialog.querySelectorAll('button, input, select').forEach(element => {
      if (element.dataset.action !== 'close') element.disabled = value;
    });
    if (!value && preview) {
      const execute = dialog.querySelector('[data-action="execute"]');
      if (execute) execute.disabled = !preview.token || !preview.item_count;
    }
    if (!state) dialog.querySelectorAll('[data-action="save"], [data-action="history"]').forEach(element => { element.disabled = true; });
  }
  function card(target, title, description, choices) {
    const policy = state.policy[target], summary = state[target];
    return `<section class="ms-storage-card" aria-labelledby="ms-storage-${target}-title">
      <h3 id="ms-storage-${target}-title">${title}</h3><p class="ms-storage-description">${description}</p>
      <div class="ms-storage-stats"><div><span>当前占用</span><strong>${bytes(summary.total_bytes)}</strong><small>${esc(summary.total_count)} 项</small></div>
        <div><span>已到保留期限</span><strong>${bytes(summary.expired_bytes)}</strong><small>${esc(summary.expired_count)} 项可供预览核查</small></div></div>
      <div class="ms-storage-setting"><label><input type="checkbox" data-enabled="${target}" ${policy.automatic_cleanup_enabled ? 'checked' : ''}>到期自动清理</label><span>${policy.automatic_cleanup_enabled ? '已开启' : '已关闭'}</span></div>
      <div class="ms-storage-setting"><label for="ms-storage-days-${target}">保留时间</label><select id="ms-storage-days-${target}" data-days="${target}">${choices.map(days => `<option value="${days}" ${days === policy.retention_days ? 'selected' : ''}>${days} 天</option>`).join('')}</select></div>
      <p class="ms-storage-fine">${target === 'recycle' ? '从移入回收站时计时；到期永久清理后无法恢复。' : '仅清理过期的导出临时文件；以后导出会重新生成。'}<br>首次开启或重新开启后，已有项目也会完整保留所选天数。</p>
      <button type="button" data-action="preview" data-target="${target}">预览到期清理</button>
    </section>`;
  }
  function render() {
    preview = null; policyDraft = null;
    const active = state.active_tasks.count || 0;
    dialog.querySelector('[data-body]').innerHTML = `<div class="ms-storage-grid">
      ${card('recycle', '任务回收站', '移除的分析项目先放在这里，永久清理前可到运行队列的“回收站”恢复。', [7,30,90])}
      ${card('cache', '导出临时缓存', '生成报告时留下的临时文件，与分析结果、正式报告和谱库分别管理。', [1,7,30])}
      </div><div class="ms-storage-protected"><strong>保留你的科研记录</strong>原始数据、当前分析结果、正式报告、日志、谱库和模型均不在清理范围内。${active ? `当前有 ${esc(active)} 个待处理或活动任务，临时缓存清理会暂停。` : '分析、排队或报告导出期间，临时缓存清理会暂停。'}<br>自动检查仅在软件后台运行时进行，最多每天清理一次；关闭自动清理后仍可手动预览。</div><section data-preview hidden class="ms-storage-preview" aria-live="polite"></section>`;
    const last = state.last_cleanup;
    dialog.querySelector('[data-last]').textContent = last ? `最近清理：${time(last.at)} · 已清理 ${last.removed_count || 0} 项 / ${bytes(last.removed_bytes)}${last.failure_count ? `，${last.failure_count} 项未完成` : ''}` : '尚无清理记录';
    notice('回收站和临时缓存分别设置；更改后点击“保存设置”。清理磁盘文件与释放后台内存是两个操作。');
    setBusy(false);
  }
  async function load() {
    if (busy) return;
    setBusy(true);
    notice('正在读取存储摘要…');
    try { state = await request(); render(); }
    catch (error) { state = null; dialog.querySelector('[data-body]').innerHTML = ''; notice(error.message, 'error'); setBusy(false); }
  }
  function clearPreview() {
    preview = null; policyDraft = null;
    const section = dialog.querySelector('[data-preview]');
    if (section) { section.hidden = true; section.innerHTML = ''; }
  }
  function draft() {
    return {
      recycle_automatic_cleanup_enabled:dialog.querySelector('[data-enabled="recycle"]').checked,
      recycle_retention_days:Number(dialog.querySelector('[data-days="recycle"]').value),
      cache_automatic_cleanup_enabled:dialog.querySelector('[data-enabled="cache"]').checked,
      cache_retention_days:Number(dialog.querySelector('[data-days="cache"]').value)
    };
  }
  async function save(confirmed = false) {
    if (busy || !state) return;
    const values = confirmed ? policyDraft : draft();
    if (!values) return;
    const enabling = ['recycle', 'cache'].some(key => values[`${key}_automatic_cleanup_enabled`] && (!state.policy[key].automatic_cleanup_enabled || values[`${key}_retention_days`] < state.policy[key].retention_days));
    if (enabling && !confirmed) {
      clearPreview(); policyDraft = values;
      const section = dialog.querySelector('[data-preview]');
      section.hidden = false;
      section.innerHTML = `<h3>确认到期自动清理设置？</h3><p>任务回收站：${values.recycle_automatic_cleanup_enabled ? `保留 ${values.recycle_retention_days} 天后永久清理` : '关闭'}；导出临时缓存：${values.cache_automatic_cleanup_enabled ? `保留 ${values.cache_retention_days} 天后清理` : '关闭'}。从关闭切换为开启时，已有项目也会获得完整保留期；缩短已开启的保留时间，会使已有项目提前到期。到期清理会自动执行并留下记录，回收站项目永久清理后无法恢复。</p><button type="button" data-action="confirm-policy">确认并保存</button><button type="button" data-action="cancel-preview">取消</button>`;
      return;
    }
    setBusy(true);
    try {
      await request('/policy', {method:'PUT', body:JSON.stringify({...values, confirm_auto_delete:confirmed})});
      state = await request(); render(); notice('设置已保存。自动清理仅处理达到保留期限且通过安全检查的项目。', 'success');
      if (window.Tasks?.loadRecycleBin && !document.getElementById('taskRecycleBinArea')?.hidden) window.Tasks.loadRecycleBin();
    } catch (error) { notice(error.message, 'error'); }
    finally { setBusy(false); }
  }
  function blockedMessage(reasons) {
    if (reasons.some(value => value.includes('active_tasks'))) return '当前有待处理或活动任务，请在任务结束后重新预览。';
    if (reasons.some(value => /lock|open_file/.test(value))) return '有导出任务正在使用临时文件，或暂时无法确认文件使用状态；本次未清理。';
    return '部分项目未通过完整性检查，本次未清理。请保留现有文件，核查后再操作。';
  }
  async function showPreview(target) {
    if (busy || !state) return;
    const current = draft();
    if (['recycle','cache'].some(key => current[`${key}_retention_days`] !== state.policy[key].retention_days || current[`${key}_automatic_cleanup_enabled`] !== state.policy[key].automatic_cleanup_enabled)) {
      notice('有尚未保存的设置。请先保存设置，再按新保留时间预览。', 'error'); return;
    }
    clearPreview(); setBusy(true);
    try {
      preview = await request('/preview', {method:'POST', body:JSON.stringify({target})});
      const section = dialog.querySelector('[data-preview]');
      const blocked = preview.blocked_reasons || [];
      section.hidden = false;
      section.innerHTML = `<h3>${target === 'recycle' ? '回收站' : '临时缓存'} · 到期清理预览</h3>
        <p>${blocked.length ? esc(blockedMessage(blocked)) : `本次将${target === 'recycle' ? '永久清理' : '清理'} ${esc(preview.item_count)} 项，预计释放 ${bytes(preview.size_bytes)}。${target === 'recycle' ? '清理后无法恢复。' : ''}预览有效期为 2 分钟，执行前会再次检查文件和任务状态。`}</p>
        ${preview.items?.length ? `<ul>${preview.items.slice(0,100).map(item => `<li>${esc(item.filename || item.name || item.task_id)} · ${bytes(item.size_bytes)}</li>`).join('')}${preview.items.length > 100 ? '<li>更多项目已列入本次预览。</li>' : ''}</ul>` : ''}
        <button type="button" data-action="execute" ${!preview.token || !preview.item_count ? 'disabled' : ''}>确认${target === 'recycle' ? '永久' : ''}清理 ${esc(preview.item_count)} 项</button><button type="button" data-action="cancel-preview">取消</button>`;
    } catch (error) { notice(error.message, 'error'); }
    finally { setBusy(false); }
  }
  async function execute() {
    if (busy || !preview?.token || !preview.item_count) return;
    const token = preview.token; preview = null; setBusy(true);
    try {
      const result = await request('/clean', {method:'POST', body:JSON.stringify({token, confirmed:true})});
      state = await request(); render();
      notice(result.failure_count ? `已清理 ${result.removed_count || 0} 项，${result.failure_count} 项未完成，逐项记录已保存。` : `已清理 ${result.removed_count || 0} 项，释放 ${bytes(result.removed_bytes)}，逐项记录已保存。`, result.failure_count ? 'error' : 'success');
      if (window.Tasks?.loadRecycleBin && !document.getElementById('taskRecycleBinArea')?.hidden) window.Tasks.loadRecycleBin();
    } catch (error) { clearPreview(); notice(`${error.message} 未继续清理，请重新预览。`, 'error'); }
    finally { setBusy(false); }
  }
  async function history() {
    if (busy || !state) return;
    clearPreview(); setBusy(true);
    try {
      const data = await request('/history?limit=10');
      const section = dialog.querySelector('[data-preview]');
      section.hidden = false;
      section.innerHTML = `<h3>最近清理记录</h3>${data.records.length ? data.records.map(record => `<p><strong>${esc(time(record.finished_at))} · ${record.target === 'recycle' ? '任务回收站' : '导出临时缓存'}</strong><br>已清理 ${esc(record.removed_count)} 项 / ${bytes(record.removed_bytes)}，未完成 ${esc(record.failure_count)} 项。</p>${record.removed_items?.length ? `<details><summary>查看已清理项目</summary><ul>${record.removed_items.map(item => `<li>${esc(item.name || item.filename || item.task_id || '项目')} · ${bytes(item.size_bytes)}</li>`).join('')}</ul></details>` : ''}${record.failed_items?.length ? `<ul>${record.failed_items.map(item => `<li>${esc(item.name || item.filename || item.task_id || '项目')}：${esc(item.error)}</li>`).join('')}</ul>` : ''}`).join('') : '<p>尚无清理记录。</p>'}${data.rejected_count ? '<p>部分记录无法完整读取，已保留原文件。</p>' : ''}<button data-action="cancel-preview">收起记录</button>`;
    } catch (error) { notice(error.message, 'error'); }
    finally { setBusy(false); }
  }
  function create() {
    dialog = document.createElement('dialog'); dialog.className = 'ms-storage-dialog';
    dialog.setAttribute('aria-labelledby', 'ms-storage-title');
    dialog.innerHTML = `<header class="ms-storage-head"><div><h2 id="ms-storage-title">存储与回收站管理</h2><p>设置保留时间，预览到期项目，按需释放磁盘空间。</p></div><button class="ms-storage-close" data-action="close" aria-label="关闭存储管理">×</button></header><div class="ms-storage-content"><div class="ms-storage-notice" data-notice role="status"></div><div data-body></div></div><footer class="ms-storage-footer"><p data-last>尚无清理记录</p><div><button data-action="history">清理记录</button><button data-action="refresh">刷新</button><button data-action="save">保存设置</button></div></footer>`;
    document.body.appendChild(dialog);
    dialog.addEventListener('click', event => {
      const button = event.target.closest('button[data-action]');
      if (!button) return;
      switch (button.dataset.action) {
        case 'close': dialog.close(); break;
        case 'refresh': load(); break;
        case 'save': save(); break;
        case 'confirm-policy': save(true); break;
        case 'preview': showPreview(button.dataset.target); break;
        case 'execute': execute(); break;
        case 'history': history(); break;
        case 'cancel-preview': clearPreview(); break;
      }
    });
    dialog.addEventListener('change', () => {
      clearPreview(); notice('设置尚未保存，点击“保存设置”后生效。');
    });
    dialog.addEventListener('close', () => { clearPreview(); previousFocus?.focus(); });
  }
  function open() {
    if (!dialog) create();
    if (!dialog.open) { previousFocus = document.activeElement; dialog.showModal(); }
    load();
  }
  window.MSsystemStorage = {open};
})();
