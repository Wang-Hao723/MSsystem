/**
 * mssystem molecular-network workspace.
 * Cytoscape is used only for presentation; network construction stays in the
 * backend so the same nodes/edges are exported to GraphML.
 */
var NetworkPage = (function () {
  'use strict';

  var cy = null;
  var currentTask = null;
  var lastNetwork = null;
  var labelsVisible = false;
  var singletonsVisible = false;
  var activeRequest = 0;

  function byId(id) { return document.getElementById(id); }
  function esc(value) {
    return String(value == null ? '' : value)
      .replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;')
      .replace(/"/g, '&quot;').replace(/'/g, '&#039;');
  }
  function number(value, digits, fallback) {
    var n = Number(value);
    return Number.isFinite(n) ? n.toFixed(digits) : (fallback || '—');
  }

  function init() {
    [
      ['netCosThresh', 'netCosThreshVal', 2],
      ['netMinPeaks', 'netMinPeaksVal', 0],
      ['netTopk', 'netTopkVal', 0],
    ].forEach(function (item) {
      var input = byId(item[0]);
      var output = byId(item[1]);
      if (!input || !output) return;
      input.addEventListener('input', function () {
        output.textContent = Number(input.value).toFixed(item[2]);
      });
    });
  }

  function refreshSelect() {
    if (typeof populateTaskSelects === 'function') populateTaskSelects();
  }

  function loadTask(taskId) {
    currentTask = taskId || null;
    lastNetwork = null;
    activeRequest += 1;
    var button = byId('btnBuildNetwork');
    if (button) button.disabled = !currentTask;
    if (cy) { cy.destroy(); cy = null; }
    if (byId('networkArea')) byId('networkArea').classList.add('hidden');
    if (byId('networkStats')) byId('networkStats').innerHTML = '';
    if (byId('networkFamilies')) byId('networkFamilies').innerHTML = '';
    if (byId('networkNodeInfo')) byId('networkNodeInfo').innerHTML = '可查看前体 m/z、注释来源、匹配分数与原始 MS/MS 谱图。';
    if (byId('networkNodeSpectrum')) byId('networkNodeSpectrum').innerHTML = '';
    if (byId('networkWorkflowAudit')) byId('networkWorkflowAudit').innerHTML = '';
    if (currentTask) {
      API.get('/tasks/' + encodeURIComponent(currentTask) + '/status')
        .then(function (task) { renderWorkflowAudit(task); })
        .catch(function () {});
    }
  }

  function renderWorkflowAudit(task) {
    var area = byId('networkWorkflowAudit');
    if (!area) return;
    var workflow = task.untargeted_workflow || {};
    var attachment = workflow.ms2_attachment || {};
    var consensus = task.consensus_spectrum_summary || workflow.consensus_spectra || {};
    var status = workflow.status || '未记录';
    var linked = Number(attachment.linked_ms2_spectra || 0);
    var total = Number(attachment.ms2_spectra || task.total_spectra || 0);
    var featureCount = Number(workflow.feature_count || 0);
    var fallback = Math.max(0, total - linked);
    area.innerHTML = '<strong>节点依据：</strong>' +
      'OpenMS 状态 ' + esc(status) + '；检测 ' + featureCount.toLocaleString() +
      ' 个 MS1 色谱特征，链接 MS2 ' + linked.toLocaleString() + '/' +
      total.toLocaleString() + '，生成 ' +
      Number(consensus.consensus_spectra || 0).toLocaleString() +
      ' 张特征级共识谱。' +
      (fallback ? ' 其余 ' + fallback.toLocaleString() + ' 张谱图会明确标为前体事件回退，不与真实 LC–MS feature 混称。' : '');
  }

  function params() {
    return {
      cos_thresh: byId('netCosThresh').value,
      min_peaks: byId('netMinPeaks').value,
      topk: byId('netTopk').value,
      topk_policy: byId('netTopkPolicy').value,
      mz_tol: byId('netMzTol').value,
      max_component_size: byId('netMaxComponent').value,
    };
  }

  function paramStr() {
    var p = params();
    return '?' + Object.keys(p).map(function (key) {
      return encodeURIComponent(key) + '=' + encodeURIComponent(p[key]);
    }).join('&');
  }

  function build() {
    if (!currentTask) { Toast.show('请先选择任务', 'warning'); return; }
    if (typeof cytoscape === 'undefined') {
      Toast.show('网络绘图库未加载，请检查网络连接后刷新页面', 'error');
      return;
    }

    var requestId = ++activeRequest;
    var button = byId('btnBuildNetwork');
    button.disabled = true;
    button.textContent = '正在构建…';
    byId('networkArea').classList.remove('hidden');
    byId('networkEmpty').style.display = 'grid';
    byId('networkEmpty').innerHTML = '<div><strong>正在计算谱图相似性</strong>大任务可能需要一些时间，请稍候…</div>';

    API.get('/tasks/' + encodeURIComponent(currentTask) + '/network' + paramStr())
      .then(function (network) {
        if (requestId !== activeRequest) return;
        lastNetwork = network;
        labelsVisible = false;
        singletonsVisible = false;
        renderStats(network.stats || {});
        renderFamilies(network.families || []);
        renderGraph(network);
        button.disabled = false;
        button.textContent = '重新构建';
        var familyCount = network.stats.n_molecular_families || 0;
        var singletonCount = network.stats.singleton_count || 0;
        Toast.show('网络构建完成：' + familyCount + ' 个连接家族，' +
          singletonCount + ' 个孤立点', 'success');
        (network.stats.warnings || []).forEach(function (warning) {
          if (warning && warning.message) {
            Toast.show(warning.message, 'warning', 7000);
          }
        });
      })
      .catch(function (error) {
        if (requestId !== activeRequest) return;
        button.disabled = false;
        button.textContent = '重新构建';
        byId('networkEmpty').innerHTML = '<div><strong style="color:#b42318">构建失败</strong>' + esc(error.message || error) + '</div>';
        Toast.show('分子网络构建失败：' + (error.message || error), 'error');
      });
  }

  function renderStats(stats) {
    var cards = [
      ['#4f7ea8', stats.n_nodes, '共识节点'],
      ['#57738c', stats.n_edges, '证据边'],
      ['#6f65a8', stats.n_molecular_families, '连接家族'],
      ['#7b8794', stats.singleton_count, '孤立点'],
      ['#009e73', stats.n_matched, '谱库匹配'],
      ['#e69f00', stats.n_propagated, '注释传播'],
      ['#40a070', stats.n_ion_identity_linked, '离子身份链接'],
      ['#2f6f9f', stats.largest_family, '最大家族'],
      ['#446b82', (stats.input_node_coverage_pct == null ? 0 : stats.input_node_coverage_pct) + '%', '输入节点覆盖'],
      ['#3978a8', (stats.full_peak_coverage_pct == null ? (stats.raw_peak_coverage_pct || 0) : stats.full_peak_coverage_pct) + '%', '完整/共识峰覆盖'],
      ['#2c7a7b', stats.lcms_feature_nodes || 0, 'LC–MS 特征节点'],
      ['#b7791f', stats.fallback_precursor_event_nodes || 0, '回退节点'],
      ['#9a6b35', stats.insufficient_fragment_nodes || 0, '低信息孤立节点'],
    ];
    byId('networkStats').innerHTML = cards.map(function (card) {
      return '<div class="net-stat" style="--accent:' + card[0] + '"><span class="net-stat-num">' +
        esc(card[1] == null ? 0 : card[1]) + '</span><span class="net-stat-lbl">' + card[2] + '</span></div>';
    }).join('') + '<div class="net-stat" style="--accent:#64748b" title="阈值建边 ' +
      esc(stats.n_edges_threshold == null ? stats.n_edges : stats.n_edges_threshold) +
      '；' + esc(stats.topk_policy || 'mutual') + ' Top-K 后 ' + esc(stats.n_edges_after_topk == null ? stats.n_edges : stats.n_edges_after_topk) +
      '；家族上限额外移除 ' + esc(stats.n_edges_removed_component_cap || 0) + '">' +
      '<span class="net-stat-num">' + esc(stats.n_edges_removed_component_cap || 0) +
      '</span><span class="net-stat-lbl">家族上限删边</span></div>';
  }

  function renderFamilies(families) {
    var connected = families.filter(function (family) { return family.size > 1; }).slice(0, 14);
    if (!connected.length) {
      byId('networkFamilies').innerHTML = '<div style="color:#82909e;font-size:.73rem;line-height:1.55">当前阈值下没有多节点家族。可适当降低 Cosine 阈值或最少匹配峰数后重建。</div>';
      return;
    }
    byId('networkFamilies').innerHTML = connected.map(function (family) {
      var title = family.top_compound || family.consensus_class || '未注释家族';
      var subtitle = family.consensus_class || (family.matched_count ? family.matched_count + ' 个库匹配' : '暂无库匹配');
      return '<button class="net-family" data-family="' + Number(family.id) + '" onclick="NetworkPage.focusFamily(' + Number(family.id) + ')">' +
        '<span class="net-family-id">F' + Number(family.id) + '</span>' +
        '<span class="net-family-size">' + Number(family.size) + '</span>' +
        '<span class="net-family-name"><strong title="' + esc(title) + '">' + esc(title) + '</strong>' + esc(subtitle) + '</span></button>';
    }).join('');
  }

  function makeLayout(count) {
    if (count <= 650) {
      return {
        name: 'cose', animate: false, randomize: true, fit: true, padding: 44,
        idealEdgeLength: 92, nodeRepulsion: 12000, edgeElasticity: 90,
        nestingFactor: 1.15, gravity: 0.22, componentSpacing: 120, numIter: 1800,
      };
    }
    if (count <= 2200) {
      return {
        name: 'cose', animate: false, randomize: true, fit: true, padding: 36,
        idealEdgeLength: 70, nodeRepulsion: 9000, gravity: 0.28,
        componentSpacing: 75, numIter: 420,
      };
    }
    return {
      name: 'concentric', animate: false, fit: true, padding: 35, minNodeSpacing: 8,
      concentric: function (node) { return node.data('degree') || 0; },
      levelWidth: function () { return 2; },
    };
  }

  function assignFamilyPositions(nodes) {
    var groups = {};
    nodes.forEach(function (element) {
      if (element.data.degree === 0) return;
      var key = String(element.data.family_id);
      if (!groups[key]) groups[key] = [];
      groups[key].push(element);
    });
    var families = Object.keys(groups).map(function (key) { return groups[key]; })
      .sort(function (a, b) { return b.length - a.length; });
    var columns = Math.max(1, Math.ceil(Math.sqrt(families.length)));
    var cell = 285;
    families.forEach(function (family, familyIndex) {
      family.sort(function (a, b) { return b.data.degree - a.data.degree; });
      var cx = (familyIndex % columns) * cell;
      var cyPos = Math.floor(familyIndex / columns) * cell;
      family.forEach(function (element, index) {
        // Golden-angle spiral keeps each capped family compact without overlap.
        var angle = index * 2.3999632297;
        var radius = index === 0 ? 0 : 13.5 * Math.sqrt(index);
        element.position = { x: cx + Math.cos(angle) * radius, y: cyPos + Math.sin(angle) * radius };
        element.data.familyX = element.position.x;
        element.data.familyY = element.position.y;
      });
    });
    nodes.forEach(function (element) {
      if (!element.position) element.position = { x: -500, y: -500 };
    });
    return families.length;
  }

  function renderGraph(network) {
    var nodes = (network.nodes || []).map(function (node) {
      var fullLabel = node.compound_name || node.propagated_compound || ('m/z ' + number(node.precursor_mz, 2));
      var shortLabel = fullLabel.length > 28 ? fullLabel.slice(0, 27) + '…' : fullLabel;
      return { data: {
        id: 'n' + node.id,
        rawId: node.id,
        label: shortLabel,
        fullLabel: fullLabel,
        ann: node.annotation_source || 'unknown',
        degree: Number(node.degree || 0),
        precursor_mz: Number(node.precursor_mz || 0),
        family_id: node.family_id == null ? -1 : Number(node.family_id),
        classification: node.classification || node.propagated_class || '',
        score: node.final_score == null ? -1 : Number(node.final_score),
        feature_id: node.feature_id || '',
        spectrum_count: Number(node.spectrum_count || 1),
        peak_source: node.peak_source || 'unknown',
        network_eligibility: node.network_eligibility || 'similarity_eligible',
        ion_identity_network_id: node.ion_identity_network_id || '',
      }, classes: node.annotation_source || 'unknown' };
    });
    var edges = (network.edges || []).map(function (edge, index) {
      return { data: {
        id: 'e' + index, source: 'n' + edge.source, target: 'n' + edge.target,
        cosine: Number(edge.cosine || 0), matched_peaks: Number(edge.matched_peaks || 0),
        edge_type: edge.edge_type || 'spectral',
        evidence_layers: (edge.evidence_layers || []).join('|'),
        modified_cosine: Number(edge.modified_cosine || 0),
        neutral_loss_cosine: Number(edge.neutral_loss_cosine || 0),
      }, classes: edge.edge_type || 'spectral' };
    });

    if (cy) cy.destroy();
    byId('networkEmpty').style.display = nodes.length ? 'none' : 'grid';
    if (!nodes.length) {
      byId('networkEmpty').innerHTML = '<div><strong>没有可显示的网络节点</strong>请检查前体 m/z、峰表和空白过滤状态。</div>';
      return;
    }

    var connectedCount = nodes.filter(function (node) { return node.data.degree > 0; }).length;
    var useFamilyOverview = connectedCount > 1200;
    if (useFamilyOverview) assignFamilyPositions(nodes);

    cy = cytoscape({
      container: byId('networkGraph'),
      elements: { nodes: nodes, edges: edges },
      style: [
        { selector: 'node', style: {
          'width': 'mapData(degree, 0, 20, 24, 54)', 'height': 'mapData(degree, 0, 20, 24, 54)',
          'background-color': '#7b8794', 'border-width': 2, 'border-color': '#55616d',
          'label': '', 'font-family': 'Inter, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif',
          'font-size': 10, 'font-weight': 600, 'color': '#26394c', 'text-valign': 'bottom',
          'text-margin-y': 7, 'text-wrap': 'ellipsis', 'text-max-width': 110,
          'text-background-color': '#ffffff', 'text-background-opacity': 0.9,
          'text-background-shape': 'roundrectangle', 'text-background-padding': 3,
          'min-zoomed-font-size': 8,
        }},
        { selector: 'node.library', style: { 'background-color': '#009e73', 'border-color': '#006b4f' } },
        { selector: 'node.propagated', style: { 'background-color': '#e69f00', 'border-color': '#9f6d00' } },
        { selector: 'node.ion_identity', style: { 'background-color': '#40a070', 'border-color': '#226748' } },
        { selector: 'node.unknown', style: { 'background-color': '#7b8794', 'border-color': '#55616d' } },
        { selector: 'node.labels-on', style: { 'label': 'data(label)' } },
        { selector: 'edge', style: {
          'width': 'mapData(cosine, 0.5, 1, 0.7, 4.5)', 'line-color': '#8fa2b3',
          'line-opacity': 0.34, 'curve-style': 'haystack', 'haystack-radius': 0.7,
        }},
        { selector: 'edge.analog, edge.neutral_loss', style: {
          'line-color': '#7762a8', 'line-style': 'dashed', 'line-opacity': 0.52,
          'curve-style': 'bezier',
        }},
        { selector: 'edge.ion_identity', style: {
          'line-color': '#2f9b68', 'line-style': 'dotted', 'line-opacity': 0.72,
          'width': 2.2, 'curve-style': 'bezier',
        }},
        { selector: 'node:selected', style: {
          'border-width': 5, 'border-color': '#2563eb', 'overlay-color': '#60a5fa', 'overlay-opacity': 0.12,
        }},
        { selector: '.dim', style: { 'opacity': 0.075 } },
        { selector: 'node.focused', style: { 'opacity': 1, 'border-width': 4, 'border-color': '#2563eb', 'z-index': 20 } },
        { selector: 'edge.focused', style: { 'opacity': 0.88, 'line-color': '#4b7fa8', 'z-index': 12 } },
        { selector: 'node.hovered', style: { 'border-width': 5, 'border-color': '#1d4ed8', 'z-index': 30 } },
        { selector: 'edge.hovered', style: { 'line-opacity': 0.95, 'line-color': '#326b98', 'z-index': 25 } },
      ],
      layout: useFamilyOverview ? { name: 'preset', fit: true, padding: 44 } : makeLayout(connectedCount),
      wheelSensitivity: 0.18,
      minZoom: 0.03,
      maxZoom: 5,
      boxSelectionEnabled: true,
    });

    cy.nodes().filter(function (node) { return node.data('degree') === 0; }).style('display', 'none');
    bindGraphEvents();
    setTimeout(function () { if (cy) cy.fit(cy.elements(':visible'), 42); }, 0);

    if (useFamilyOverview) {
      Toast.show('大型网络已按分子家族分区排列；点击右侧家族可聚焦查看', 'info');
    }
  }

  function bindGraphEvents() {
    cy.on('tap', 'node', function (event) {
      showNode(event.target.data('rawId'));
    });
    cy.on('mouseover', 'node', function (event) {
      var node = event.target;
      node.addClass('hovered');
      node.connectedEdges().addClass('hovered');
      showTooltip(node, event.renderedPosition || node.renderedPosition());
    });
    cy.on('mouseout', 'node', function () {
      cy.elements().removeClass('hovered');
      hideTooltip();
    });
    cy.on('tap', function (event) {
      if (event.target === cy) hideTooltip();
    });
  }

  function showTooltip(node, position) {
    var tip = byId('networkTooltip');
    var familyLabel = Number(node.data('family_id')) >= 0
      ? '家族 F' + esc(node.data('family_id')) : '孤立点';
    tip.innerHTML = '<strong>' + esc(node.data('fullLabel')) + '</strong><br>' +
      'm/z ' + number(node.data('precursor_mz'), 4) + ' · ' + familyLabel +
      '<br>连接 ' + esc(node.data('degree')) + ' · ' + esc(node.data('classification') || '未分类');
    tip.style.left = Math.min(position.x + 14, byId('networkGraph').clientWidth - 270) + 'px';
    tip.style.top = Math.max(8, position.y - 12) + 'px';
    tip.style.display = 'block';
  }

  function hideTooltip() {
    var tip = byId('networkTooltip');
    if (tip) tip.style.display = 'none';
  }

  function nodeMeta(rawId) {
    if (!lastNetwork) return null;
    return (lastNetwork.nodes || []).find(function (node) { return Number(node.id) === Number(rawId); }) || null;
  }

  function showNode(rawId) {
    if (!currentTask) return;
    var meta = nodeMeta(rawId) || {};
    byId('networkSelectionHint').textContent = '#' + rawId;
    byId('networkNodeInfo').innerHTML = '<div style="color:#82909e;font-size:.74rem">正在加载谱图…</div>';

    API.get('/tasks/' + encodeURIComponent(currentTask) + '/network/node/' + encodeURIComponent(rawId) + paramStr())
      .then(function (detail) {
        var match = detail.match || {};
        var source = meta.annotation_source || (match.compound_name ? 'library' : 'unknown');
        var badge = source === 'library' ? '谱库匹配' :
          source === 'propagated' ? '谱图传播' :
          source === 'ion_identity' ? '离子身份链接' : '未知';
        var name = match.compound_name || meta.propagated_compound || ('m/z ' + number(detail.precursor_mz, 4));
        var peakEvidence = detail.peak_source === 'lcms_feature_representative'
          ? '实测单扫描代表谱'
          : (detail.peak_source === 'lcms_feature_consensus'
          ? 'LC–MS 特征级 MS2 共识谱'
          : (detail.peak_source === 'raw'
            ? '完整 centroid MS/MS 碎片峰'
            : 'Top-N 搜索峰回退'));
        var peakMeta = detail.peak_evidence || {};
        var ce = peakMeta.collision_energy == null
          ? '—'
          : number(peakMeta.collision_energy, 1) +
            (peakMeta.collision_energy_unit ? ' ' + peakMeta.collision_energy_unit : '');
        var rows = [
          ['Precursor m/z', number(detail.precursor_mz, 4)],
          ['峰谱图编号', peakMeta.spectrum_id || detail.spectrum_id || '—'],
          ['峰谱产物编号', peakMeta.artifact_id || detail.lcms_consensus_spectrum_id || '—'],
          ['峰谱结果索引', peakMeta.result_index == null ? '—' : peakMeta.result_index],
          ['峰谱碰撞能', ce],
          ['峰谱活化方式', peakMeta.activation_method || detail.activation_method || '—'],
          ['候选注释谱图', detail.annotation_spectrum_id || '—'],
          ['候选注释结果索引', detail.annotation_result_index == null ? '—' : detail.annotation_result_index],
          ['特征节点', meta.feature_id || '—'],
          ['节点来源', detail.feature_source === 'openms_lcms_feature' ? 'OpenMS LC–MS 色谱特征' : 'MS2 前体事件回退'],
          ['合并 MS2 数', meta.spectrum_count || 1],
          ['峰证据', peakEvidence],
          ['RT', detail.retention_time == null ? '—' : detail.retention_time],
          ['分子家族', meta.family_id == null ? '孤立点' : 'F' + meta.family_id],
          ['连接数', meta.degree == null ? '—' : meta.degree],
          ['分类', match.classification || meta.classification || meta.propagated_class || '—'],
        ];
        if (match.confidence) rows.push(['匹配置信度', match.confidence]);
        if (match.final_score != null) rows.push(['Final Score', number(match.final_score, 3)]);
        if (match.cosine != null) rows.push(['Cosine', number(match.cosine, 3)]);
        if (source === 'propagated') rows.push(['证据等级', '家族传播候选（非确认）']);
        if (source === 'ion_identity') rows.push(['证据等级', 'MS1 共洗脱/离子形式链接（非独立鉴定）']);
        if (meta.ion_identity_network_id) rows.push(['IIMN', meta.ion_identity_network_id]);

        byId('networkNodeInfo').innerHTML = '<div class="net-node-head"><span class="net-badge net-badge-' + source + '">' + badge + '</span>' +
          '<span class="net-node-name">' + esc(name) + '</span></div><div class="net-node-grid">' +
          rows.map(function (row) { return '<div class="net-kv"><span class="net-kv-k">' + esc(row[0]) + '</span><span class="net-kv-v">' + esc(row[1]) + '</span></div>'; }).join('') + '</div>';
        renderSpectrum(detail.peaks || [], name);
      })
      .catch(function (error) {
        byId('networkNodeInfo').innerHTML = '<div style="color:#b42318;font-size:.74rem">节点加载失败：' + esc(error.message || error) + '</div>';
      });
  }

  function renderSpectrum(peaks, name) {
    if (!peaks.length) {
      byId('networkNodeSpectrum').innerHTML = '<div style="padding:25px 0;color:#82909e;font-size:.72rem;text-align:center">无可显示碎片峰</div>';
      return;
    }
    var maxIntensity = Math.max.apply(null, peaks.map(function (peak) { return Number(peak[1]) || 0; })) || 1;
    var x = [], y = [], hover = [];
    peaks.forEach(function (peak) {
      var mz = Number(peak[0]);
      var relative = Number(peak[1]) / maxIntensity * 100;
      x.push(mz, mz, null); y.push(0, relative, null);
      hover.push('', 'm/z ' + mz.toFixed(4) + '<br>' + relative.toFixed(1) + '%', '');
    });
    Plotly.react('networkNodeSpectrum', [{
      type: 'scatter', mode: 'lines', x: x, y: y, text: hover,
      hovertemplate: '%{text}<extra></extra>', line: { color: '#3978a8', width: 1.4 },
    }], {
      title: { text: 'MS/MS · ' + String(name).slice(0, 32), font: { size: 11, color: '#41576b' }, x: 0.03 },
      xaxis: { title: 'm/z', titlefont: { size: 10 }, tickfont: { size: 9 }, gridcolor: '#edf1f5', zeroline: false },
      yaxis: { title: '%', range: [0, 108], titlefont: { size: 10 }, tickfont: { size: 9 }, gridcolor: '#edf1f5', zeroline: false },
      margin: { t: 38, b: 38, l: 38, r: 8 }, height: 235, paper_bgcolor: 'rgba(0,0,0,0)', plot_bgcolor: '#fff',
    }, { responsive: true, displayModeBar: false });
  }

  function focusElements(nodes) {
    if (!cy || !nodes || !nodes.length) return;
    var edges = nodes.edgesWith(nodes);
    cy.batch(function () {
      cy.elements().removeClass('focused dim');
      cy.elements().addClass('dim');
      nodes.removeClass('dim').addClass('focused');
      edges.removeClass('dim').addClass('focused');
    });
    cy.animate({ fit: { eles: nodes.union(edges), padding: 75 }, duration: 320 });
  }

  function focusFamily(familyId) {
    if (!cy) return;
    var nodes = cy.nodes().filter(function (node) { return Number(node.data('family_id')) === Number(familyId); });
    if (!nodes.length) return;
    focusElements(nodes);
    document.querySelectorAll('.net-family').forEach(function (button) {
      button.classList.toggle('active', Number(button.dataset.family) === Number(familyId));
    });
  }

  function search() {
    if (!cy) return;
    var query = (byId('netSearch').value || '').trim().toLowerCase();
    if (!query) { resetFocus(); return; }
    var nodes = cy.nodes().filter(function (node) {
      return String(node.data('fullLabel') || '').toLowerCase().indexOf(query) >= 0 ||
        String(node.data('classification') || '').toLowerCase().indexOf(query) >= 0 ||
        String(node.data('precursor_mz') || '').indexOf(query) >= 0;
    });
    if (!nodes.length) { Toast.show('没有找到匹配节点', 'warning'); return; }
    focusElements(nodes);
    nodes.select();
    if (nodes.length === 1) showNode(nodes[0].data('rawId'));
    Toast.show('找到 ' + nodes.length + ' 个节点', 'success');
  }

  function resetFocus() {
    if (!cy) return;
    cy.elements().removeClass('focused dim');
    cy.nodes().unselect();
    document.querySelectorAll('.net-family').forEach(function (button) { button.classList.remove('active'); });
    fit();
  }

  function toggleLabels() {
    if (!cy) return;
    labelsVisible = !labelsVisible;
    cy.nodes().toggleClass('labels-on', labelsVisible);
    byId('netLabelBtn').textContent = labelsVisible ? '隐藏标签' : '显示标签';
  }

  function toggleSingletons() {
    if (!cy) return;
    singletonsVisible = !singletonsVisible;
    cy.nodes().filter(function (node) { return node.data('degree') === 0; })
      .style('display', singletonsVisible ? 'element' : 'none');
    var button = byId('netSingletonBtn');
    if (button) button.textContent = singletonsVisible ? '隐藏孤立点' : '显示孤立点';
    fit();
  }

  function fit() {
    if (cy) cy.animate({ fit: { eles: cy.elements(':visible'), padding: 44 }, duration: 260 });
  }

  function relayout() {
    if (!cy) return;
    var visibleCount = cy.nodes(':visible').length;
    var hasFamilyPositions = cy.nodes(':visible').filter(function (node) {
      return Number.isFinite(Number(node.data('familyX')));
    }).length > 0;
    var layout = hasFamilyPositions && visibleCount > 1200
      ? { name: 'preset', fit: true, padding: 44, positions: function (node) {
          return { x: Number(node.data('familyX')), y: Number(node.data('familyY')) };
        }}
      : makeLayout(visibleCount);
    cy.layout(layout).run();
  }

  function downloadPng() {
    if (!cy) { Toast.show('请先构建分子网络', 'warning'); return; }
    var link = document.createElement('a');
    link.href = cy.png({ full: true, scale: 3, bg: '#ffffff', maxWidth: 12000, maxHeight: 12000 });
    link.download = (currentTask || 'network') + '_molecular_network.png';
    document.body.appendChild(link); link.click(); link.remove();
  }

  function downloadGraphml() {
    if (!currentTask) { Toast.show('请先选择任务', 'warning'); return; }
    API.download('/tasks/' + encodeURIComponent(currentTask) + '/network/graphml' + paramStr(), currentTask + '_network.graphml');
  }

  function downloadNodesCsv() {
    if (!currentTask) { Toast.show('请先选择任务', 'warning'); return; }
    API.download('/tasks/' + encodeURIComponent(currentTask) + '/network/nodes.csv' + paramStr(), currentTask + '_network_nodes.csv');
  }

  function downloadEdgesCsv() {
    if (!currentTask) { Toast.show('请先选择任务', 'warning'); return; }
    API.download('/tasks/' + encodeURIComponent(currentTask) + '/network/edges.csv' + paramStr(), currentTask + '_network_edges.csv');
  }

  async function alignSelectedTasks() {
    var select = byId('rtAlignmentTasks');
    var status = byId('rtAlignmentStatus');
    var taskIds = Array.from(select?.selectedOptions || []).map(function (option) {
      return option.value;
    });
    if (taskIds.length < 2) {
      Toast.show('请至少选择两个包含 OpenMS featureXML 的任务', 'warning');
      return;
    }
    if (status) status.textContent = '正在执行 RT 对齐和跨样本 feature grouping…';
    try {
      var result = await API.post('/workflows/rt-align', {
        task_ids: taskIds,
        rt_tolerance_seconds: Number(byId('rtAlignmentTolerance')?.value || 20),
        mz_tolerance: Number(byId('rtAlignmentMzTolerance')?.value || 10),
        mz_unit: 'ppm',
      });
      if (status) {
        status.textContent = '完成：' + Number(result.consensus_count || 0).toLocaleString() +
          ' 个跨样本共识特征；参考样本 ' + (result.reference || '—') +
          '；ID ' + (result.alignment_id || '—');
      }
      Toast.show('跨样本 RT 对齐完成', 'success');
    } catch (error) {
      if (status) status.textContent = '对齐失败：' + (error.message || error);
      Toast.show('RT 对齐失败：' + (error.message || error), 'error');
    }
  }

  return {
    init: init,
    refreshSelect: refreshSelect,
    loadTask: loadTask,
    build: build,
    showNode: showNode,
    search: search,
    focusFamily: focusFamily,
    resetFocus: resetFocus,
    toggleLabels: toggleLabels,
    toggleSingletons: toggleSingletons,
    fit: fit,
    relayout: relayout,
    downloadPng: downloadPng,
    downloadGraphml: downloadGraphml,
    downloadNodesCsv: downloadNodesCsv,
    downloadEdgesCsv: downloadEdgesCsv,
    alignSelectedTasks: alignSelectedTasks,
  };
})();
