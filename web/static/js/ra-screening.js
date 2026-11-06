/**
 * MSsystem — Rheumatoid Arthritis (RA) Targeted Screening Page
 * 对已完成分析结果做抗 RA 活性化合物解读 (后处理视图, 不改核心检索):
 *   - direct     命中化合物在抗 RA 参考库内 → 已证实抗 RA 活性
 *   - propagated 节点所在分子家族含 RA 直接命中 → 候选结构类似物
 *   - class      命中分类属 RA 高相关类 → RA 相关化学类别
 * 依赖: plotly (已加载), API 层 (fetch /api/tasks/{id}/ra-screening)
 */
var RAScreeningPage = (function () {
  var currentTask = null;
  var lastReport = null;

  var MECH_COLORS = {
    "免疫调节": "#22c55e", "免疫调节/抗炎": "#16a34a", "免疫调节/镇痛": "#65a30d",
    "NF-κB抑制": "#ef4444", "COX-2抑制": "#f59e0b", "JAK抑制": "#3b82f6",
    "抗炎": "#84cc16", "抗炎/抗氧化": "#a3e635", "镇痛": "#8b5cf6", "镇痛/抗炎": "#7c3aed",
    "叶酸代谢拮抗": "#0ea5e9", "嘧啶合成抑制": "#06b6d4",
  };
  var TOX_CLASS = { "高": "ra-tox-high", "中": "ra-tox-mid", "低": "ra-tox-low" };

  function mechColor(m) { return MECH_COLORS[m] || "#94a3b8"; }
  function badge(flag) {
    var map = { direct: ["ra-badge-direct", "直接命中"], propagated: ["ra-badge-propagated", "候选类似物"], class: ["ra-badge-class", "相关类别"] };
    var m = map[flag] || ["ra-badge-class", flag];
    return '<span class="ra-badge ' + m[0] + '">' + m[1] + '</span>';
  }
  function toxBadge(t) {
    if (!t) return "";
    return '<span class="ra-badge ' + (TOX_CLASS[t] || "ra-tox-low") + '">毒性 ' + t + '</span>';
  }
  function esc(s) {
    if (s === null || s === undefined) return "";
    return String(s).replace(/[&<>"]/g, function (c) { return { "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c]; });
  }
  function fmt(v, d) { d = d === undefined ? 4 : d; return (v === null || v === undefined || v === "") ? "—" : (typeof v === "number" ? Number(v).toFixed(d) : v); }

  function init() {
    var sel = document.getElementById('raTaskSelect');
    if (sel) sel.addEventListener('change', function () { loadTask(this.value); });
  }

  function loadTask(taskId) {
    currentTask = taskId || null;
    if (!taskId) { showEmpty("请在上方选择一个已完成分析的任务，查看抗 RA 成分筛查结果。"); return; }
    showLoading();
    fetch("/api/tasks/" + encodeURIComponent(taskId) + "/ra-screening")
      .then(function (r) {
        if (!r.ok) throw new Error("HTTP " + r.status);
        return r.json();
      })
      .then(function (rep) {
        lastReport = rep;
        render(rep);
      })
      .catch(function (err) {
        showEmpty("RA 筛查加载失败：" + err.message + "（请确认任务已完成分析）");
        if (typeof Toast !== "undefined" && Toast.show) Toast.show("RA 筛查加载失败: " + err.message, "error");
      });
  }

  // ── 渲染 ────────────────────────────────────────────────────────────
  function render(rep) {
    var sum = rep.summary || {};
    var nDirect = sum.n_direct || 0;
    var nProp = sum.n_propagated || 0;
    var nClass = sum.n_class_relevant || 0;

    if (nDirect === 0 && nProp === 0 && nClass === 0) {
      showEmpty("该任务未检出与抗 RA 参考库相关的化合物。");
      return;
    }

    show("raEmpty", false); show("raCharts", true); show("raDirectCard", true);
    show("raPropCard", nProp > 0); show("raClassCard", nClass > 0);

    // 汇总卡片
    var s = document.getElementById("raSummary");
    s.style.display = "flex";
    s.innerHTML = statCard("ra-stat-direct", nDirect, "抗 RA 直接命中") +
      statCard("ra-stat-prop", nProp, "候选类似物(家族传播)") +
      statCard("ra-stat-class", nClass, "RA 相关化学类别") +
      statCard("", (rep.summary || {}).panel_size || "—", "参考库化合物数");

    renderCharts(sum);
    renderDirectTable(rep.direct_hits || []);
    renderPropTable(rep.propagated || [], sum);
    renderClassTable(rep.class_relevant || []);
  }

  function statCard(cls, num, lbl) {
    return '<div class="ra-stat ' + cls + '"><span class="ra-stat-num">' + num + '</span><span class="ra-stat-lbl">' + lbl + '</span></div>';
  }

  function renderCharts(sum) {
    var mech = sum.by_mechanism || {};
    var herb = sum.by_herb || {};
    barChart("raMechChart", "按抗 RA 机制", mech, function (k) { return mechColor(k); });
    barChart("raHerbChart", "按来源药材", herb, "#0ea5e9");
  }

  function barChart(elId, title, obj, colorFn) {
    var el = document.getElementById(elId);
    if (!el) return;
    var entries = Object.keys(obj).map(function (k) { return [k, obj[k]]; }).sort(function (a, b) { return b[1] - a[1]; });
    if (!entries.length) { el.innerHTML = '<div class="ra-empty">无数据</div>'; return; }
    var colorFn2 = typeof colorFn === "function" ? colorFn : function () { return colorFn; };
    Plotly.newPlot(el, [{
      type: "bar", orientation: "h",
      x: entries.map(function (e) { return e[1]; }),
      y: entries.map(function (e) { return e[0]; }),
      marker: { color: entries.map(function (e) { return colorFn2(e[0]); }), line: { width: 1, color: "#fff" } },
      text: entries.map(function (e) { return e[1]; }), textposition: "outside",
      hovertemplate: "%{y}: %{x}<extra></extra>",
    }], {
      title: { text: title, font: { size: 13 } },
      margin: { t: 30, b: 20, l: 110, r: 20 },
      height: Math.max(160, entries.length * 28 + 40),
      yaxis: { autorange: "reversed", tickfont: { size: 11 } },
      xaxis: { tickfont: { size: 10 }, dtick: 1 },
      paper_bgcolor: "white",
    }, { responsive: true, displayModeBar: false });
  }

  function renderDirectTable(rows) {
    var html = '<thead><tr>' +
      th("化合物") + th("抗 RA 活性成分") + th("来源药材") + th("机制") +
      th("m/z") + th("RT") + th("匹配分") + th("毒性") + th("活性说明") +
      '</tr></thead><tbody>';
    rows.forEach(function (r) {
      var e = r.ra_entry || {};
      html += '<tr>' +
        td(esc(r.compound_name), "font-weight:600;") +
        td(esc(e.name_en) + (e.name_cn ? '<br><span style="color:#718096;font-size:0.74rem;">' + esc(e.name_cn) + '</span>' : '')) +
        td(esc(e.source_herb)) +
        td('<span class="ra-badge" style="background:' + mechColor(e.mechanism) + '22;color:#1f2937;">' + esc(e.mechanism) + '</span>') +
        td(fmt(r.precursor_mz, 2)) + td(fmt(r.retention_time, 2)) +
        td(scoreCell(r.final_score)) + td(toxBadge(e.toxicity_flag)) +
        td('<span style="font-size:0.74rem;color:#64748b;">' + esc(e.activity_note) + '</span>') +
        '</tr>';
    });
    html += '</tbody>';
    document.getElementById("raDirectTable").innerHTML = html;
  }

  function renderPropTable(rows, sum) {
    var html = '<thead><tr>' +
      th("节点#") + th("谱图ID") + th("m/z") + th("RT") + th("峰数") +
      th("注释来源") + th("传播化合物") + th("传播分类") + '</tr></thead><tbody>';
    rows.forEach(function (r) {
      html += '<tr>' +
        td(esc(r.node_idx)) +
        td('<span style="font-family:monospace;font-size:0.74rem;">' + esc(r.spectrum_id) + '</span>') +
        td(fmt(r.precursor_mz, 2)) + td(fmt(r.retention_time, 2)) +
        td(esc(r.peaks_count)) +
        td(annotBadge(r.annotation_source)) +
        td(esc(r.propagated_compound)) + td(esc(r.propagated_class)) +
        '</tr>';
    });
    html += '</tbody>';
    var card = document.getElementById("raPropTable").closest(".card");
    if (sum && sum.n_propagated > (rows.length)) {
      html = '<caption style="caption-side:top;text-align:left;padding:4px 8px;font-size:0.76rem;color:#92400e;">显示前 ' + rows.length + ' 条，共 ' + sum.n_propagated + ' 个候选节点</caption>' + html;
    }
    document.getElementById("raPropTable").innerHTML = html;
  }

  function renderClassTable(rows) {
    var html = '<thead><tr>' +
      th("化合物") + th("分类") + th("m/z") + th("RT") + th("匹配分") + th("InChIKey") +
      '</tr></thead><tbody>';
    rows.forEach(function (r) {
      html += '<tr>' +
        td(esc(r.compound_name), "font-weight:600;") +
        td(esc(r.classification)) +
        td(fmt(r.precursor_mz, 2)) + td(fmt(r.retention_time, 2)) +
        td(scoreCell(r.final_score)) +
        td('<span style="font-family:monospace;font-size:0.72rem;color:#94a3b8;">' + esc((r.inchikey || "").slice(0, 14)) + '</span>') +
        '</tr>';
    });
    html += '</tbody>';
    document.getElementById("raClassTable").innerHTML = html;
  }

  function annotBadge(src) {
    var m = { library: ["ra-badge-class", "库匹配"], propagated: ["ra-badge-propagated", "传播"], unknown: ["ra-badge-class", "未知"] };
    var v = m[src] || ["ra-badge-class", esc(src)];
    return '<span class="ra-badge ' + v[0] + '">' + v[1] + '</span>';
  }
  function scoreCell(s) {
    if (s === null || s === undefined) return td("—");
    var v = Number(s); var col = v >= 0.7 ? "#16a34a" : (v >= 0.4 ? "#d97706" : "#94a3b8");
    return '<td style="font-weight:600;color:' + col + ';">' + v.toFixed(3) + '</td>';
  }
  function th(t) { return '<th>' + t + '</th>'; }
  function td(t, style) { return '<td' + (style ? ' style="' + style + '"' : '') + '>' + t + '</td>'; }

  // ── 状态切换 ────────────────────────────────────────────────────────
  function show(id, on) { var el = document.getElementById(id); if (el) el.style.display = on ? "" : "none"; }
  function showEmpty(msg) {
    show("raSummary", false); show("raCharts", false);
    show("raDirectCard", false); show("raPropCard", false); show("raClassCard", false);
    var e = document.getElementById("raEmpty"); e.style.display = ""; e.textContent = msg;
  }
  function showLoading() {
    show("raSummary", false); show("raCharts", false);
    show("raDirectCard", false); show("raPropCard", false); show("raClassCard", false);
    var e = document.getElementById("raEmpty"); e.style.display = ""; e.textContent = "⏳ 正在加载 RA 筛查结果…";
  }

  return { init: init, loadTask: loadTask, refresh: function () { if (currentTask) loadTask(currentTask); } };
})();
