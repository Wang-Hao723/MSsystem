/**
 * MSsystem — Main Application Controller
 * Navigation, initialization, onboarding, toast, modal management
 */
const App = (() => {
  let currentPage = 'dashboard';

  function tr(keyOrSource, params) {
    return window.MSSLanguage?.t(keyOrSource, params) ?? keyOrSource;
  }

  // ── Initialization ────────────────────────────────────────────────
  function init() {
    // Init modules
    Upload.init();
    Tasks.init();
    Results.init();
    if (typeof SampleOverview !== 'undefined') SampleOverview.init();
    if (typeof EvidenceReview !== 'undefined') EvidenceReview.init();
    if (typeof SampleReport !== 'undefined') SampleReport.init();
    ExportMod.init();
    if (typeof CompletedBatches !== 'undefined') CompletedBatches.init();
    if (typeof Dashboard !== "undefined") Dashboard.init();
    if (typeof LibAugmentPage !== "undefined") LibAugmentPage.init();
    if (typeof LibraryCenterPage !== "undefined") LibraryCenterPage.init();

    // Navigation
    document.querySelectorAll('.nav-item').forEach(btn => {
      btn.addEventListener('click', () => {
        navigateTo(btn.dataset.page);
        if (window.innerWidth < 1024) toggleSidebar(false);
      });
    });

    // Enter through the analysis workbench. Hashes keep refresh/deep-link behavior.
    const requestedPage = window.location.hash.replace(/^#/, '').split('?')[0];
    const initialPage = document.getElementById(`page-${requestedPage}`)
      ? requestedPage
      : 'dashboard';
    navigateTo(initialPage, false);
    window.addEventListener('hashchange', () => {
      const page = window.location.hash.replace(/^#/, '').split('?')[0];
      if (document.getElementById(`page-${page}`)) {
        navigateTo(page, false);
      }
    });
    document.addEventListener('keydown', event => {
      if (event.key === 'Escape' && document.querySelector('#modalContainer .modal-overlay')) {
        closeModal();
      }
    });
    document.addEventListener('mssystem:languagechange', () => {
      updateCurrentPageLabel(currentPage);
    });

    // Check library status
    Upload.checkLibraryStatus();

    // Connect SSE for real-time updates
    API.connectSSE();

    // The dashboard itself is the onboarding surface; no blocking first-run modal.
  }

  // ── Navigation ────────────────────────────────────────────────────
  function updateCurrentPageLabel(page) {
    const parentPage = page === 'libaugment'
      ? 'library'
      : (page === 'report' ? 'export' : page);
    const activeNav = document.querySelector(`.nav-item[data-page="${parentPage}"]`);
    const pageLabel = document.getElementById('currentPageLabel');
    if (!pageLabel) return;
    const isEnglish = window.MSSLanguage?.getLanguage() === 'en';
    if (page === 'libaugment') {
      pageLabel.textContent = tr('page.libaugment');
    } else if (page === 'report') {
      pageLabel.textContent = tr('page.report');
    } else if (isEnglish) {
      pageLabel.textContent = activeNav?.dataset.labelEn || 'I.D.E.A platform';
    } else {
      const primaryLabel = activeNav?.dataset.label || 'MSsystem';
      const englishLabel = activeNav?.dataset.labelEn || 'I.D.E.A platform';
      pageLabel.textContent = `${primaryLabel} · ${englishLabel}`;
    }
  }

  function navigateTo(page, updateHash = true) {
    currentPage = page;
    const parentPage = page === 'libaugment'
      ? 'library'
      : (page === 'report' ? 'export' : page);
    document.querySelectorAll('.nav-item').forEach(b => b.classList.remove('active'));
    const activeNav = document.querySelector(`.nav-item[data-page="${parentPage}"]`);
    activeNav?.classList.add('active');
    document.querySelectorAll('.page').forEach(p => p.classList.remove('active'));
    document.getElementById(`page-${page}`)?.classList.add('active');
    updateCurrentPageLabel(page);
    if (updateHash && window.location.hash !== `#${page}`) {
      const url = new URL(window.location.href);
      url.hash = page;
      history.replaceState(null, '', `${url.pathname}${url.search}${url.hash}`);
    }

    // Load data for the target page
    if (page === 'dashboard' && typeof Dashboard !== 'undefined') Dashboard.load(false);
    if (page === 'tasks') Tasks.loadTasks();
    if (page === 'results' && typeof SampleOverview !== 'undefined') {
      SampleOverview.refreshSelect();
    }
    if (page === 'review' && typeof EvidenceReview !== 'undefined') {
      EvidenceReview.enterPage();
    }
    if (page === 'export') { ExportMod.refreshSelect(); }
    if (page === 'completed-batches' && typeof CompletedBatches !== 'undefined') {
      CompletedBatches.enterPage();
    }
    if (page === 'library' && typeof LibraryCenterPage !== 'undefined') {
      LibraryCenterPage.ensureLoaded();
    }
    // 谱图/分类/分子网络/报告 页面: 导航时刷新任务下拉框,
    // 确保新增/删除任务后这些下拉里的列表是最新 (修复下拉不自动更新)
    if (page === 'spectrum' || page === 'classification' ||
        page === 'network') {
      if (typeof populateTaskSelects === 'function') populateTaskSelects();
    }
    if (page === 'report' && typeof SampleReport !== 'undefined') {
      SampleReport.refreshSelect();
    }
  }

  // ── Sidebar Toggle ────────────────────────────────────────────────
  function toggleSidebar(show) {
    const sidebar = document.getElementById('sidebar');
    const overlay = document.getElementById('sidebarOverlay');
    if (show === undefined) show = !sidebar.classList.contains('open');
    sidebar.classList.toggle('open', show);
    overlay.classList.toggle('show', show);
  }

  // ── Toast Notifications ───────────────────────────────────────────
  function showToast(message, type = 'info', duration = 4500) {
    const container = document.getElementById('toastContainer');
    if (!container) return;
    const icons = { success: '✅', error: '❌', warning: '⚠️', info: 'ℹ️' };
    const toast = document.createElement('div');
    toast.className = `toast ${type}`;
    const icon = document.createElement('span');
    icon.textContent = icons[type] || icons.info;
    icon.setAttribute('aria-hidden', 'true');
    const copy = document.createElement('span');
    copy.style.flex = '1';
    copy.textContent = tr(message);
    const close = document.createElement('button');
    close.type = 'button';
    close.className = 'toast-close';
    close.setAttribute('aria-label', tr('toast.close'));
    close.textContent = '×';
    close.addEventListener('click', () => toast.remove());
    toast.append(icon, copy, close);
    container.appendChild(toast);
    if (duration > 0) setTimeout(() => toast.remove(), duration);
  }

  // ── Modal ─────────────────────────────────────────────────────────
  function showModal(html) {
    const container = document.getElementById('modalContainer');
    if (!container) return;
    container.innerHTML = html;

    // Detail dialogs can be much taller than the viewport (spectra, evidence
    // passports and literature audit).  Give every modal a real event-listener
    // close control that remains visible while the modal body scrolls.  This is
    // deliberately not an inline onclick handler: it stays usable under a
    // stricter Content-Security-Policy and avoids trapping users at the bottom
    // of a long result view.
    const modal = container.querySelector('.modal');
    if (modal && !modal.querySelector('.modal-close-floating')) {
      const closeButton = document.createElement('button');
      closeButton.type = 'button';
      closeButton.className = 'modal-close-floating';
      closeButton.setAttribute('aria-label', tr('modal.close'));
      closeButton.setAttribute('title', tr('modal.closeWithEsc'));
      closeButton.textContent = '×';
      closeButton.addEventListener('click', closeModal);
      modal.prepend(closeButton);
    }
    window.MSSLanguage?.translate(container);
  }

  function closeModal() {
    const container = document.getElementById('modalContainer');
    if (container) container.innerHTML = '';
  }

  // ── Onboarding Guide ──────────────────────────────────────────────
  function showOnboarding() {
    const html = `
    <div class="guide-overlay" id="guideOverlay">
      <div class="guide-modal">
        <h2>欢迎使用 MSsystem</h2>
        <p><strong>I.D.E.A platform</strong><br>正、负离子统一的 LC–MS/MS 非靶向分析与分子注释平台</p>
        <div class="guide-steps">
          <div class="guide-step">
            <div class="guide-step-num">1</div>
            <h4>数据与质控<br><small>Data &amp; Quality Control</small></h4>
            <p>核对 MS1/MS2、极性和空白</p>
          </div>
          <div class="guide-step">
            <div class="guide-step-num">2</div>
            <h4>分析与解释<br><small>Analysis &amp; Interpretation</small></h4>
            <p>谱库鉴定、谱图和分子网络</p>
          </div>
          <div class="guide-step">
            <div class="guide-step-num">3</div>
            <h4>报告交付<br><small>Report Delivery</small></h4>
            <p>统一生成数据表与完整报告</p>
          </div>
        </div>
        <div class="guide-actions">
          <button class="btn btn-secondary" onclick="App.dismissGuide()">跳过</button>
          <button class="btn btn-primary" onclick="App.dismissGuide()">开始使用</button>
        </div>
      </div>
    </div>`;
    document.body.insertAdjacentHTML('beforeend', html);
    window.MSSLanguage?.translate(document.getElementById('guideOverlay'));
  }

  function dismissGuide() {
    const overlay = document.getElementById('guideOverlay');
    if (overlay) overlay.remove();
    localStorage.setItem('mssystem_v7_workspace_guide_seen', '1');
  }

  // ── DOMContentLoaded ──────────────────────────────────────────────
  document.addEventListener('DOMContentLoaded', init);

  return {
    navigateTo, toggleSidebar, showToast, showModal, closeModal,
    showOnboarding, dismissGuide,
  };
})();

// ── Global shortcuts for inline event handlers ──
const Toast = { show: App.showToast };

// ═══════════════════════════════════════════════════════════
// Evidence modules — task-select population and deferred initialization.
// ═══════════════════════════════════════════════════════════

// Initialize evidence viewers after their scripts and the main controller exist.
var origInit = App._origInit || null;
document.addEventListener("DOMContentLoaded", function() {
    // Wait for original init to complete
    setTimeout(function() {
        if (typeof SpectrumPage !== "undefined") SpectrumPage.init();
        if (typeof ClassPage !== "undefined") ClassPage.init();
        if (typeof NetworkPage !== "undefined") NetworkPage.init();
        if (typeof RAScreeningPage !== "undefined") RAScreeningPage.init();
        populateTaskSelects();
        console.log("MSsystem modules initialized");
    }, 500);
});

function populateTaskSelects() {
    fetch("/api/tasks?status=completed").then(function(r) { return r.json(); }).then(function(data) {
        if (!data || !data.tasks) return;
        var taskList = data.tasks;
        var selects = ["specTaskSelect", "classTaskSelect", "networkTaskSelect", "raTaskSelect"];
        selects.forEach(function(id) {
            var sel = document.getElementById(id);
            if (!sel) return;
            var val = sel.value;
            const placeholder = window.MSSLanguage?.t('-- 选择已完成的任务 --')
                || '-- 选择已完成的任务 --';
            sel.innerHTML = `<option value="">${placeholder}</option>`;
            taskList.forEach(function(t) {
                var opt = document.createElement("option");
                opt.value = t.id;
                opt.textContent = window.MSSLanguage?.getLanguage() === 'en'
                    ? t.filename + " (" + t.matched_spectra + "/" + t.total_spectra + " matched)"
                    : t.filename + " (" + t.matched_spectra + "/" + t.total_spectra + " 匹配)";
                sel.appendChild(opt);
            });
            sel.value = val;
        });
        var alignSelect = document.getElementById("rtAlignmentTasks");
        if (alignSelect) {
            var selected = Array.from(alignSelect.selectedOptions || []).map(function(option) { return option.value; });
            alignSelect.innerHTML = "";
            taskList.forEach(function(t) {
                var workflow = t.untargeted_workflow || {};
                var option = document.createElement("option");
                option.value = t.id;
                option.textContent = t.filename + " · " +
                    (workflow.feature_count || 0) + " MS1 features · " +
                    (workflow.status || window.MSSLanguage?.t("未运行") || "未运行");
                option.disabled = !(workflow.feature_xml || workflow.output_features);
                option.selected = selected.indexOf(t.id) >= 0;
                alignSelect.appendChild(option);
            });
        }
    }).catch(function() {});
}
