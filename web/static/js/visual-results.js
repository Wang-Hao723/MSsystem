/** Visual projection of saved results; never changes candidate identities or scores. */
const VisualResults = (() => {
  const esc = value => String(value ?? '').replace(/[&<>"']/g, c => ({
    '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;',
  })[c]);
  const finite = value => value == null || String(value).trim() === ''
    ? null : (Number.isFinite(Number(value)) ? Number(value) : null);
  const fmt = (value, digits = 3) => finite(value) == null ? '—' : Number(value).toFixed(digits);
  const n = value => Number(value || 0).toLocaleString();
  const states = {
    confirmed: ['已确证', 'Confirmed', '#17605b'],
    high: ['High · 二维连接', 'High · 2-D connectivity', '#248c75'],
    isomer: ['异构体组', 'Isomer group', '#7b68ad'],
    medium: ['Medium 候选', 'Medium candidate', '#bd8c38'],
    background: ['空白／背景', 'Blank / background', '#ad6477'],
    unresolved_event: ['未解析事件', 'Unresolved event', '#b16c42'],
    low: ['证据不足', 'Insufficient evidence', '#8093a6'],
    unannotated: ['未注释特征', 'Unannotated feature', '#bec9d2'],
    saved_candidate: ['已保存谱级候选', 'Saved spectral candidate', '#628da2'],
  };
  const info = group => states[group?.status?.key] || states.low;
  function rt(row) {
    if (finite(row?.retention_time) == null) return null;
    return SampleResultData.retentionTimeMinutes(row);
  }
  function name(group, english = false) {
    if (group.status.key === 'unresolved_event') return english ? 'Unresolved precursor event' : '未解析母离子事件';
    if (group.status.key === 'unannotated') return english ? 'Unannotated ion feature' : '未注释离子特征';
    return String(group.compound?.display_name || group.rep?.match?.compound_name || (english ? 'Unnamed candidate' : '未命名候选'));
  }
  function structure(group, taskId, english) {
    const match = group.rep?.match || {};
    const smiles = String(match.smiles || '').trim();
    const index = finite(group.rep?.result_index);
    const selectedSmiles = String(group.flagship?.candidate?.smiles || '').trim();
    const entityId = group.entity?.entity_id;
    let path = '';
    const drawable = smiles && match.has_structure !== false && !['none', 'null', 'nan', 'n/a'].includes(smiles.toLowerCase());
    if (drawable) {
      if (entityId && selectedSmiles === smiles) {
        path = `/api/tasks/${encodeURIComponent(taskId)}/flagship/structure?entity_id=${encodeURIComponent(entityId)}&format=svg`;
      } else if (index != null && Number.isInteger(index) && index >= 0) {
        path = `/api/tasks/${encodeURIComponent(taskId)}/download/structure/result/${index}?format=svg`;
      }
    }
    const empty = english ? 'Structure unavailable' : '无可用结构图';
    const label = group.status.key === 'isomer'
      ? (english ? 'One candidate in this isomer group' : '异构体组中的一个候选')
      : group.status.key === 'unresolved_event'
      ? (english ? 'Ranking lead for review' : '供复核的排序候选')
      : (english ? 'Candidate 2-D structure' : '候选二维结构');
    return `<div class="vr-structure ${drawable ? '' : 'vr-structure-empty'}">
      ${drawable ? `<canvas id="vr-molecule-${Number(group.index)}" width="220" height="145" data-vr-smiles="${esc(smiles)}" role="img" aria-label="${esc(label)}"></canvas>` : ''}
      ${path ? `<img loading="lazy" decoding="async" data-vr-fallback-src="${esc(path)}" alt="${esc(label)}" style="display:none" onerror="this.parentElement.classList.add('vr-structure-empty');this.remove()">` : ''}
      <span class="vr-structure-placeholder">${empty}</span>
      ${drawable ? `<small>${label}</small>` : ''}
    </div>`;
  }
  function mountStructures(containerId) {
    const container = document.getElementById(containerId);
    if (!container) return;
    const viewer = typeof MoleculeViewer !== 'undefined' ? new MoleculeViewer({ width: 220, height: 145 }) : null;
    container.querySelectorAll('canvas[data-vr-smiles]').forEach(canvas => {
      if (viewer?.drawer) {
        viewer.render(canvas.id, canvas.dataset.vrSmiles);
        return;
      }
      const stage = canvas.parentElement;
      const fallback = stage.querySelector('img[data-vr-fallback-src]');
      canvas.remove();
      if (fallback) {
        fallback.style.display = '';
        fallback.src = fallback.dataset.vrFallbackSrc;
      } else stage.classList.add('vr-structure-empty');
    });
  }
  function badge(group, english) {
    const state = info(group);
    const isomerLevel = group.status.key === 'isomer'
      ? ` · ${group.status.hasIsomerHigh ? 'High' : 'Medium'}` : '';
    return `<span class="vr-badge" style="--vr-tone:${state[2]}">${esc(state[english ? 1 : 0])}${isomerLevel}</span>`;
  }
  function similarity(value, english) {
    const parsed = finite(value);
    return `<div class="vr-similarity"><span>${english ? 'Cosine' : '碎片余弦'} <b>${fmt(parsed)}</b></span><span class="vr-meter" aria-hidden="true"><i style="width:${parsed == null ? 0 : Math.max(0, Math.min(1, parsed)) * 100}%"></i></span></div>`;
  }
  function scatter(groups, english) {
    const points = groups.map(group => ({ group, x: rt(group.rep), y: finite(group.rep?.precursor_mz) }))
      .filter(p => p.x != null && p.y != null && p.x >= 0 && p.y > 0);
    if (!points.length) return `<div class="vr-chart-empty">${english ? 'No recorded RT–m/z coordinates' : '没有可绘制的 RT–m/z 记录'}</div>`;
    const maxX = Math.max(...points.map(p => p.x), .1) * 1.04;
    const minY = Math.max(0, Math.min(...points.map(p => p.y)) * .9);
    const maxY = Math.max(...points.map(p => p.y)) * 1.05;
    const px = x => 48 + x / maxX * 488;
    const py = y => 174 - (y - minY) / (maxY - minY || 1) * 152;
    const ticks = Array.from({ length: 5 }, (_, i) => {
      const x = maxX * i / 4, y = minY + (maxY - minY) * i / 4;
      return `<line x1="48" x2="536" y1="${py(y)}" y2="${py(y)}" class="vr-grid"/><text x="40" y="${py(y) + 4}" text-anchor="end">${y.toFixed(0)}</text><text x="${px(x)}" y="195" text-anchor="middle">${x.toFixed(1)}</text>`;
    }).join('');
    return `<svg class="vr-scatter" viewBox="0 0 560 220" role="group" aria-label="${english ? 'Representative RT–m/z points, select a point to view evidence' : '代表谱RT与m/z分布，选点查看证据'}">
      ${ticks}<text x="7" y="12">m/z</text><text x="535" y="216" text-anchor="end">RT / min</text>
      ${points.map(p => `<circle cx="${px(p.x).toFixed(2)}" cy="${py(p.y).toFixed(2)}" r="4" fill="${info(p.group)[2]}" tabindex="0" role="button" aria-label="${esc(name(p.group, english))}, RT ${fmt(p.x)} min, m/z ${fmt(p.y, 5)}" onclick="SampleOverview.openEvidence(${Number(p.group.index)})" onkeydown="if(event.key==='Enter'||event.key===' '){event.preventDefault();SampleOverview.openEvidence(${Number(p.group.index)})}"><title>${esc(name(p.group, english))}\nRT ${fmt(p.x)} min · m/z ${fmt(p.y, 5)}</title></circle>`).join('')}
      </svg><small class="vr-plot-note">${n(points.length)} / ${n(groups.length)} ${english ? 'groups/events with recorded coordinates · one representative point per row' : '个身份／事件有坐标 · 每个点为一条结果的代表谱'}</small>`;
  }
  function summary(bundle, counts, english = false) {
    const groups = bundle.groups || [];
    const byState = (bundle.savedOnly ? ['saved_candidate'] : Object.keys(states).filter(key => key !== 'saved_candidate'))
      .map(key => ({ key, count: groups.filter(g => (states[g.status.key] ? g.status.key : 'low') === key).length }));
    const bars = byState.map(({ key, count }) => `<span style="width:${groups.length ? count / groups.length * 100 : 0}%;background:${states[key][2]}" title="${esc(states[key][english ? 1 : 0])}: ${count}"></span>`).join('');
    const legend = byState.map(({ key, count }) => {
      const filter = key === 'isomer' ? 'isomer_all' : key === 'medium' ? 'medium' : key;
      return `<button class="vr-legend-item" data-overview-filter="${filter}" onclick="SampleOverview.setFilter('${filter}')"><i style="background:${states[key][2]}"></i><span>${esc(states[key][english ? 1 : 0])}</span><b>${n(count)}</b></button>`;
    }).join('');
    const metric = (value, label, filter) => `<button class="vr-stat" onclick="SampleOverview.setFilter('${filter}')"><b>${value == null ? '—' : n(value)}</b><span>${label}</span></button>`;
    const entityCount = value => bundle.savedOnly ? null : value;
    const sourceNote = bundle.savedOnly ? `<div class="vr-source-note" role="status"><b>${english ? 'Saved spectral results' : '已保存谱级结果'}</b><span>${english ? 'Current entity decisions are unavailable. Original candidates and spectra remain available below.' : '当前实体结论不可用；下方可查看原候选与实测谱图。'}</span><details><summary>${english ? 'View reason' : '查看原因'}</summary>${esc(bundle.sourceFailure?.message || bundle.sourceWarning || '')}</details></div>` : '';
    return `<section class="vr-summary" aria-label="${english ? 'Visual result overview' : '结果图形概览'}">
      <div class="vr-summary-top"><div><span class="vr-eyebrow">RESULT EXPLORER</span><h3>${esc(bundle.task?.filename || (english ? 'Sample results' : '样本结果'))}</h3></div><div class="vr-top-actions"><button class="btn btn-secondary btn-sm" onclick="SampleOverview.openReport()">${english ? 'Report' : '综合报告'}</button><button class="btn btn-primary btn-sm" onclick="SampleOverview.openReview()">${english ? 'Spectral evidence →' : '谱图与证据 →'}</button></div></div>
      ${sourceNote}
      <div class="vr-stats">${metric(groups.length, english ? 'Identity / event rows' : '去重身份／事件', 'all')}${metric(entityCount(counts.high), english ? 'High · 2-D connectivity' : 'High · 单一二维连接', 'high')}${metric(entityCount(counts.isomerHigh), english ? 'High · isomer groups' : 'High · 异构体组', 'isomer')}${metric(entityCount(counts.confirmed), english ? 'Confirmed identities' : '已确证身份', 'confirmed')}</div>
      <div class="vr-chart-grid"><section class="vr-chart"><header><b>${english ? 'RT–m/z distribution' : 'RT–m/z 分布'}</b><small>${english ? 'Select a point for evidence' : '选点查看谱图证据'}</small></header>${scatter(groups, english)}</section><section class="vr-chart"><header><b>${english ? 'Evidence categories' : '证据分层'}</b><small>${n(groups.length)} ${english ? 'identity / event rows' : '个身份／事件'}</small></header><div class="vr-stacked-bar" aria-hidden="true">${bars}</div><div class="vr-legend">${legend}</div><small class="vr-plot-note">${english ? 'Category counts; library support and standard confirmation are separate.' : '按结果类别计数；谱库支持与标准品确证分别列出。'}</small></section></div>
    </section>`;
  }
  function catalog(groups, { taskId, mode = 'cards', page = 0, pageSize = 24, english = false } = {}) {
    const pages = Math.max(1, Math.ceil(groups.length / pageSize));
    page = Math.max(0, Math.min(pages - 1, page));
    const visible = groups.slice(page * pageSize, (page + 1) * pageSize);
    const footer = `<div class="vr-pagination"><span>${groups.length ? page * pageSize + 1 : 0}–${Math.min((page + 1) * pageSize, groups.length)} / ${n(groups.length)}</span><div><button class="btn btn-secondary btn-sm" ${page ? '' : 'disabled'} onclick="SampleOverview.changeVisualPage(-1)">${english ? 'Previous' : '上一页'}</button><span>${page + 1} / ${pages}</span><button class="btn btn-secondary btn-sm" ${page + 1 < pages ? '' : 'disabled'} onclick="SampleOverview.changeVisualPage(1)">${english ? 'Next' : '下一页'}</button></div></div>`;
    if (!visible.length) return `<div class="vr-chart-empty">${english ? 'No results match these filters. Select another category or clear the search.' : '当前筛选没有结果，请切换类别或清空搜索。'}</div>${footer}`;
    const items = visible.map(group => {
      const match = group.rep?.match || {};
      const title = name(group, english);
      const compoundName = group.status.key === 'isomer'
        ? `${String(group.leadCandidate?.compound_name || match.compound_name || title)} ${english ? 'and isomer candidates' : '及异构体候选'}` : title;
      const formula = group.compound?.formula || match.formula || '—';
      const coord = `<span>m/z <b>${fmt(group.rep?.precursor_mz, 5)}</b></span><span>RT <b>${fmt(rt(group.rep))}</b> min</span>`;
      const candidateN = Math.max(Number(group.candidateMemberCount || 0), Number(match.isomer_group_size || 0));
      const flags = [group.status.key === 'isomer' ? `${n(candidateN)} ${english ? 'competing candidates' : '个竞争候选'}` : '', group.compound?.plant_report?.excluded_from_plant_candidate_table || match.potential_contaminant ? (english ? 'External / role review' : '外源／身份复核') : ''].filter(Boolean);
      const open = `<button class="vr-open" onclick="SampleOverview.openEvidence(${Number(group.index)})">${english ? 'View spectrum & evidence' : '查看谱图与证据'} <span aria-hidden="true">↗</span></button>`;
      if (mode === 'table') return `<tr><td>${structure(group, taskId, english)}</td><td><b class="vr-table-name" data-scientific-value title="${esc(title)}">${esc(compoundName)}</b><small>${esc(formula)}</small>${badge(group, english)}${flags.length ? `<small class="vr-caution">${esc(flags.join(' · '))}</small>` : ''}</td><td>${similarity(match.cosine, english)}<small>${english ? 'Mass error' : '质量误差'} ${fmt(match.mass_ppm, 2)} ppm</small></td><td class="vr-coordinates">${coord}<small>${n(group.spectra)} MS2</small></td><td>${open}</td></tr>`;
      return `<article class="vr-compound" style="--vr-tone:${info(group)[2]}"><div class="vr-compound-top">${badge(group, english)}<span>${n(group.spectra)} MS2</span></div>${structure(group, taskId, english)}<div class="vr-compound-copy"><h4 data-scientific-value title="${esc(title)}">${esc(compoundName)}</h4><small class="vr-formula">${esc(formula)}</small>${similarity(match.cosine, english)}<div class="vr-coordinates">${coord}</div><div class="vr-card-foot"><span>Δ ${fmt(match.mass_ppm, 2)} ppm</span><span>${esc(flags.join(' · '))}</span></div>${open}</div></article>`;
    }).join('');
    return mode === 'table' ? `<div class="table-wrap vr-table-wrap"><table class="vr-table"><thead><tr><th>${english ? 'Structure' : '候选结构'}</th><th>${english ? 'Identity & evidence' : '身份与证据'}</th><th>${english ? 'Spectrum match' : '谱图匹配'}</th><th>${english ? 'Measured location' : '实测位置'}</th><th></th></tr></thead><tbody>${items}</tbody></table></div>${footer}` : `<div class="vr-compound-grid">${items}</div>${footer}`;
  }
  return { summary, catalog, mountStructures, finite, rt };
})();
