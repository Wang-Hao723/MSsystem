/** Polarity-neutral library governance center. */
var LibraryCenterPage = (function () {
  'use strict';

  var loaded = false;

  function esc(value) {
    return String(value == null ? '' : value)
      .replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;')
      .replace(/"/g, '&quot;').replace(/'/g, '&#039;');
  }

  function number(value) {
    var numeric = Number(value);
    return Number.isFinite(numeric) ? numeric.toLocaleString('zh-CN') : '动态/待审计';
  }

  function init() {
    // Styling lives in workspace-v7.css so the library page shares one design system.
  }

  async function loadRuntime() {
    var target = document.getElementById('v7RuntimeInfo');
    if (!target) return;
    try {
      var info = await API.get('/system/info');
      target.innerHTML =
        '<strong>MSsystem · I.D.E.A platform</strong>　同源 UI/API：<code>' +
        esc(info.base_url) + '</code><br>' +
        '当前运行库：<strong>' +
        esc(info.database && info.database.loaded
          ? (info.database.loaded_name || '已加载')
          : '未加载') + '</strong>　可检索谱图：' +
        esc(number(info.database && info.database.loaded_spectra)) + '<br>' +
        '项目路径：<code>' + esc(info.project_root) + '</code><br>' +
        '共享谱库：<code>' + esc(info.shared_data_root) + '</code>';
    } catch (error) {
      target.innerHTML = '<span style="color:#b42318">无法读取当前 MSsystem 实例：' +
        esc(error.message || error) + '</span>';
    }
  }

  function renderCoverage(payload) {
    var snapshot = payload.database_snapshot || {};
    var summary = payload.coverage_summary || {};
    var structure = payload.structure_quality_snapshot || {};
    var aliases = payload.compound_name_alias_audit || {};
    var cards = [
      ['全部物理谱图', number(snapshot.total_spectra), 'POS、NEG 与待核验极性'],
      ['正离子覆盖', number(snapshot.compounds_with_positive),
        (summary.positive_compound_coverage_pct == null ? '—' : summary.positive_compound_coverage_pct + '%') + ' mapped compounds'],
      ['负离子覆盖', number(snapshot.compounds_with_negative),
        (summary.negative_compound_coverage_pct == null ? '—' : summary.negative_compound_coverage_pct + '%') + ' mapped compounds'],
      ['双极性覆盖', number(snapshot.both_polarities),
        (summary.dual_polarity_coverage_pct == null ? '—' : summary.dual_polarity_coverage_pct + '%') + ' mapped compounds'],
      ['极性未知谱图', number(snapshot.unknown_polarity_spectra),
        '元数据不足，不从 precursor m/z 猜测'],
      ['极性元数据冲突', number(snapshot.polarity_metadata_conflicts),
        snapshot.polarity_rule_version || '统一规则待记录'],
      ['结构来源已记录', number(structure.with_structure_identity_source),
        '结构派生验证不冒充外部身份认证'],
      ['危险别名已隔离', number(aliases.quarantined_aliases),
        '碎片词、内部编号和 putative 标签不再作名称']
    ];
    var target = document.getElementById('negativeCoverageCards');
    if (!target) return;
    target.innerHTML = cards.map(function (card) {
      return '<div class="library-coverage-card"><span>' + esc(card[0]) +
        '</span><strong>' + esc(card[1]) + '</strong><small>' +
        esc(card[2]) + '</small></div>';
    }).join('');
  }

  function sourceLink(source) {
    var url = source.download_url || source.landing_page ||
      (source.landing_pages || [])[0] || source.download_docs;
    if (!url) return '';
    return '<a href="' + esc(url) + '" target="_blank" rel="noopener">官方来源 ↗</a>';
  }

  function renderSources(payload) {
    var target = document.getElementById('negativeSourceList');
    if (!target) return;
    var sources = payload.sources || [];
    target.innerHTML = sources.map(function (source) {
      var count = source.official_relevant_records;
      if (count == null) count = source.official_negative_records;
      if (count == null) count = source.official_negative_search_hits;
      if (count == null) count = source.official_lcmsms_records;
      var current = source.current_exact_source_records;
      return '<div class="registry-source">' +
        '<div><h4>' + esc(source.name) + '</h4><span class="registry-pill">' +
        esc(source.experimental_status || source.kind) + '</span>　' +
        (source.polarity_scope ? '<span class="registry-pill">' +
          esc(source.polarity_scope) + '</span>　' : '') +
        sourceLink(source) + '<p>' + esc(source.natural_product_relevance || '') + '</p></div>' +
        '<div><span style="color:#8190a0">公开记录（未跨库去重）</span><br><strong>' +
        esc(number(count)) + '</strong></div>' +
        '<div><span style="color:#8190a0">当前精确来源</span><br><strong>' +
        esc(current == null ? '待去重核验' : number(current)) +
        '</strong><br><span class="registry-priority">' + esc(source.priority) + '</span></div>' +
        '<div class="registry-source-action"><strong>' + esc(source.status) +
        '</strong><p>' + esc(source.action || '') + '</p></div></div>';
    }).join('') || '<div style="padding:1rem;color:#8190a0">暂无数据源记录。</div>';
  }

  async function load(force) {
    if (loaded && !force) return;
    loaded = true;
    loadRuntime();
    var target = document.getElementById('negativeSourceList');
    if (target) target.innerHTML = '<div style="color:#8190a0;padding:1rem 0;">正在加载数据源清单…</div>';
    try {
      var payload = await API.get('/library/source-registry');
      renderCoverage(payload);
      renderSources(payload);
    } catch (error) {
      loaded = false;
      if (target) {
        target.innerHTML = '<div style="color:#b42318;padding:1rem 0;">加载失败：' +
          esc(error.message || error) + '</div>';
      }
    }
  }

  return { init: init, load: function () { return load(true); }, ensureLoaded: load };
})();
