/**
 * MSsystem analysis workbench.
 * Presents one polarity-neutral workflow over runtime, library and task facts.
 */
var Dashboard = (function () {
  'use strict';

  var refreshTimer = null;

  function tr(source) {
    return window.MSSLanguage?.t(source) || source;
  }

  function isEnglish() {
    return window.MSSLanguage?.getLanguage() === 'en';
  }

  function esc(value) {
    return String(value == null ? '' : value)
      .replace(/&/g, '&amp;')
      .replace(/</g, '&lt;')
      .replace(/>/g, '&gt;')
      .replace(/"/g, '&quot;')
      .replace(/'/g, '&#039;');
  }

  function setText(id, value) {
    var element = document.getElementById(id);
    if (element) element.textContent = value;
  }

  function number(value) {
    var numeric = Number(value);
    return Number.isFinite(numeric)
      ? (window.MSSLanguage?.formatNumber(numeric) || numeric.toLocaleString())
      : '—';
  }

  function taskState(status) {
    var labels = {
      pending: '待处理',
      queued: '排队中',
      running: '运行中',
      completed: '已完成',
      failed: '失败'
    };
    return tr(labels[status] || status || '未知');
  }

  function taskPolarity(task) {
    var mode = String(task.ion_mode || '').toUpperCase();
    if (mode === 'POS') return 'POS';
    if (mode === 'NEG') return 'NEG';
    var counts = task.ion_mode_counts || {};
    if (Number(counts.POS || 0) && Number(counts.NEG || 0)) return 'POS + NEG';
    return tr('自动极性');
  }

  function renderRuntime(info) {
    setText('dashInstanceName', 'MSsystem · ' + info.port);
    setText('dashRuntimePath', info.project_root);
    setText('dashServiceStatus', tr('服务在线 · UI/API 同源'));
    setText('sidebarRuntimeState', 'I.D.E.A platform · ' + info.port);
    var runtimeChip = document.getElementById('topbarRuntimeState');
    if (runtimeChip) {
      runtimeChip.innerHTML = '<i></i>' + esc(info.host || tr('本机')) + ' · ' +
        esc(info.port);
    }
    setText(
      'dashStorageValue',
      info.storage && info.storage.free_gib != null
        ? Number(info.storage.free_gib).toFixed(1) + ' GiB'
        : '—'
    );
    var dot = document.getElementById('dashServiceDot');
    if (dot) dot.classList.add('online');
  }

  function renderLibrary(library) {
    var loaded = !!library.loaded;
    setText('dashLibraryValue', loaded ? (library.name || tr('已加载')) : tr('未加载'));
    setText(
      'dashLibraryMeta',
      loaded
        ? (isEnglish()
          ? number(library.spectra_count) + ' searchable spectra'
          : number(library.spectra_count) + ' 张可检索谱图')
        : tr('分析前选择 Full 或 Small')
    );
    setText(
      'dashLibraryReadiness',
      loaded
        ? (library.name || tr('参考库')) + (isEnglish()
          ? ' is loaded. Library analysis can start.'
          : ' 已加载，可开始谱库分析。')
        : tr('需要在分析前明确选择 Full 或 Small；系统不会静默加载大库。')
    );
    var badge = document.getElementById('dashLibraryBadge');
    if (badge) {
      badge.textContent = loaded ? tr('已加载') : tr('待加载');
      badge.classList.toggle('ready', loaded);
    }
    var topbar = document.getElementById('topbarLibraryState');
    if (topbar) {
      topbar.classList.toggle('is-ready', loaded);
      topbar.innerHTML = '<i></i>' + (
        loaded
          ? esc(library.name || tr('谱库已加载')) + ' · ' + esc(number(library.spectra_count))
          : tr('谱库未加载')
      );
    }
  }

  function renderTasks(payload) {
    var tasks = (payload && payload.tasks) || [];
    var active = tasks.filter(function (task) {
      return task.status === 'running' || task.status === 'queued' || task.status === 'pending';
    });
    var completed = tasks.filter(function (task) {
      return task.status === 'completed';
    });
    setText('dashTaskValue', active.length
      ? (isEnglish() ? number(active.length) + ' active' : number(active.length) + ' 个活动')
      : tr('空闲'));
    setText('dashTaskMeta', isEnglish()
      ? number(tasks.length) + ' run records'
      : number(tasks.length) + ' 个任务记录');
    setText('dashCompletedValue', number(completed.length));

    var target = document.getElementById('dashRecentTasks');
    if (!target) return;
    if (!tasks.length) {
      target.innerHTML =
        '<div class="dashboard-empty"><div>' + esc(tr('尚无分析任务')) + '<br>' +
        '<button class="text-action" style="margin-top:.55rem" ' +
        'onclick="App.navigateTo(\'upload\')">' + esc(tr('导入第一批数据 →')) + '</button></div></div>';
      return;
    }
    target.innerHTML = tasks.slice(0, 5).map(function (task) {
      return '<div class="recent-task"><div><strong data-i18n-skip>' + esc(task.filename || task.id) +
        '</strong><small>' + esc(taskPolarity(task)) + ' · ' +
        esc(taskState(task.status)) + ' · ' +
        esc(number(task.matched_spectra || 0)) + '/' +
        esc(number(task.total_spectra || 0)) + (isEnglish() ? ' matched' : ' 匹配') + '</small></div>' +
        '<button data-task-id="' + esc(task.id) + '" data-task-status="' +
        esc(task.status) + '" onclick="Dashboard.openTask(this.dataset.taskId,' +
        'this.dataset.taskStatus)">' + esc(tr('打开 →')) + '</button></div>';
    }).join('');
  }

  function renderUnavailable() {
    setText('dashServiceStatus', tr('服务状态读取失败'));
    var dot = document.getElementById('dashServiceDot');
    if (dot) dot.classList.remove('online');
  }

  async function load(force) {
    if (force) {
      setText('dashServiceStatus', tr('正在刷新系统状态…'));
    }
    var responses = await Promise.allSettled([
      API.get('/system/info'),
      API.get('/library/status'),
      API.get('/tasks')
    ]);
    if (responses[0].status === 'fulfilled') renderRuntime(responses[0].value);
    else renderUnavailable();
    if (responses[1].status === 'fulfilled') renderLibrary(responses[1].value);
    if (responses[2].status === 'fulfilled') renderTasks(responses[2].value);
  }

  function scheduleRefresh() {
    clearTimeout(refreshTimer);
    refreshTimer = setTimeout(function () { load(false); }, 450);
  }

  async function openTask(taskId, status) {
    if (status === 'completed') {
      if (typeof SampleOverview !== 'undefined' && SampleOverview.openTask) {
        await SampleOverview.openTask(taskId);
      } else {
        App.navigateTo('results');
      }
      return;
    }
    App.navigateTo('tasks');
  }

  function init() {
    load(false);
    document.addEventListener('mssystem:languagechange', function () {
      load(false);
    });
    if (typeof API !== 'undefined' && API.onTaskUpdate) {
      API.onTaskUpdate(scheduleRefresh);
    }
  }

  return {
    init: init,
    load: load,
    openTask: openTask
  };
})();
