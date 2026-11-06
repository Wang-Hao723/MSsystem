/**
 * LibAugmentPage — 谱库扩增模块
 * 上传表格、原始 MS2 或经复核图片峰表 → 预览 → 确认追加到双库
 *
 * 后端端点：
 *   POST /api/library/import/preview  — 预览（不写库）
 *   POST /api/library/import?db_mode=both|full|small — 确认导入
 */
var LibAugmentPage = (function () {
  var stagedFiles = [];
  var lastPreview = null;

  function init() {
    var zone = document.getElementById('libaugmentZone');
    var input = document.getElementById('libaugmentFileInput');
    if (input) {
      input.addEventListener('change', function () { handleFiles(this.files); this.value = ''; });
    }
    if (zone) {
      zone.addEventListener('click', function () { if (input) input.click(); });
      ['dragenter', 'dragover'].forEach(function (ev) {
        zone.addEventListener(ev, function (e) { e.preventDefault(); zone.classList.add('dragover'); });
      });
      ['dragleave', 'drop'].forEach(function (ev) {
        zone.addEventListener(ev, function (e) { e.preventDefault(); zone.classList.remove('dragover'); });
      });
      zone.addEventListener('drop', function (e) {
        if (e.dataTransfer && e.dataTransfer.files) handleFiles(e.dataTransfer.files);
      });
    }
    var previewBtn = document.getElementById('libaugmentPreviewBtn');
    if (previewBtn) previewBtn.addEventListener('click', doPreview);
    var importBtn = document.getElementById('libaugmentImportBtn');
    if (importBtn) importBtn.addEventListener('click', doImport);
    renderFileList();
  }

  function handleFiles(fileList) {
    var allowed = ['.csv', '.xlsx', '.xls', '.mgf', '.msp', '.mzml', '.mzxml', '.xml',
      '.json', '.png', '.jpg', '.jpeg', '.tif', '.tiff'];
    for (var i = 0; i < fileList.length; i++) {
      var f = fileList[i];
      var ext = '.' + (f.name.split('.').pop() || '').toLowerCase();
      if (allowed.indexOf(ext) === -1) {
        if (typeof Toast !== 'undefined') Toast.show('跳过不支持的文件: ' + f.name, 'warning');
        continue;
      }
      // 去重（同名）
      var dup = stagedFiles.some(function (s) { return s.name === f.name && s.size === f.size; });
      if (!dup) stagedFiles.push(f);
    }
    renderFileList();
  }

  function sourceFileExtension(name) {
    return '.' + String(name || '').split('.').pop().toLowerCase();
  }

  function renderReferenceSourceOptions() {
    var select = document.getElementById('laRefSource');
    if (!select) return;
    var previous = select.value;
    var sourceExtensions = ['.mgf', '.msp', '.mzml', '.mzxml', '.xml', '.png', '.jpg', '.jpeg', '.tif', '.tiff'];
    var options = stagedFiles.filter(function (file) {
      return sourceExtensions.indexOf(sourceFileExtension(file.name)) !== -1;
    });
    select.innerHTML = '<option value="">请选择已加入的 MGF/MSP/mzML 或图片</option>' + options.map(function (file) {
      return '<option value="' + esc(file.name) + '">' + esc(file.name) + '</option>';
    }).join('');
    if (options.some(function (file) { return file.name === previous; })) select.value = previous;
  }

  function removeFile(i) { stagedFiles.splice(i, 1); renderFileList(); }

  function renderFileList() {
    var box = document.getElementById('libaugmentFileList');
    if (!box) return;
    if (!stagedFiles.length) {
      box.innerHTML = '<div class="la-muted">尚未选择文件</div>';
      toggleButtons(false);
      renderReferenceSourceOptions();
      return;
    }
    box.innerHTML = stagedFiles.map(function (f, i) {
      return '<div class="la-file-row"><span>📄 ' + esc(f.name) + '</span>' +
        '<span class="la-muted">' + (f.size / 1024).toFixed(1) + ' KB</span>' +
        '<button class="btn btn-sm" onclick="LibAugmentPage.removeFile(' + i + ')">✕</button></div>';
    }).join('');
    toggleButtons(true);
    renderReferenceSourceOptions();
  }

  function toggleButtons(enabled) {
    var p = document.getElementById('libaugmentPreviewBtn');
    if (p) p.disabled = !enabled;
  }

  function setBusy(busy, msg) {
    var p = document.getElementById('libaugmentPreviewBtn');
    var im = document.getElementById('libaugmentImportBtn');
    if (p) p.disabled = busy;
    if (im) im.disabled = busy;
    var area = document.getElementById('libaugmentResult');
    if (area && msg) area.innerHTML = '<div class="la-muted">' + msg + '</div>';
  }

  async function doPreview() {
    if (!stagedFiles.length) { if (typeof Toast !== 'undefined') Toast.show('请先选择文件', 'warning'); return; }
    setBusy(true, '🔄 转换并与现有库比对中...');
    var fd = new FormData();
    stagedFiles.forEach(function (f) { fd.append('files', f); });
    try {
      var res = await API.post('/library/import/preview', fd);
      lastPreview = res;
      renderPreview(res);
      var im = document.getElementById('libaugmentImportBtn');
      if (im) {
        im.classList.remove('hidden');
        im.disabled = !(Number(res.importable_records || 0) > 0);
      }
    } catch (e) {
      if (typeof Toast !== 'undefined') Toast.show('预览失败: ' + e.message, 'error');
      var area = document.getElementById('libaugmentResult');
      if (area) area.innerHTML = '<div class="la-error">❌ 预览失败: ' + esc(e.message) + '</div>';
    } finally {
      var p = document.getElementById('libaugmentPreviewBtn');
      if (p) p.disabled = false;
    }
  }

  function renderPreview(res) {
    var area = document.getElementById('libaugmentResult');
    var quality = res.quality || {};
    var tierCounts = quality.source_tier_counts || {};
    var conv = (res.convert || []).map(function (c) {
      if (c.error) return '<div class="la-muted">📄 ' + esc(c.file) + ': <span class="la-error">转换错误 - ' + esc(c.error) + '</span></div>';
      return '<div class="la-muted">📄 ' + esc(c.file) + ': 转换 ' + (c.written || 0) + ' 条' +
        (c.skipped_no_frag ? '（跳过 ' + c.skipped_no_frag + ' 无碎片）' : '') + ' [' + (c.format || 'standard') + ']</div>';
    }).join('');

    var summary = '<div class="la-card">' +
      '<h4>📊 预览结果</h4>' +
      '<div class="la-stat-grid">' +
      statBox(res.importable_records || 0, '实际可导入', 'la-new') +
      statBox(res.spectral_records || 0, '可检索实验谱', 'la-new') +
      statBox(res.literature_target_records || 0, '文献靶标记录', 'la-warn') +
      statBox(res.provenance_held_records || 0, '来源待核验·暂存', 'la-warn') +
      statBox(quality.accepted || 0, '内部一致性通过', '') +
      statBox(quality.rejected || 0, '隔离/拒绝', 'la-warn') +
      statBox(res.new_compounds || 0, '🆕 新化合物', 'la-new') +
      statBox(res.existing_compounds || 0, '✅ 库中已有', '') +
      statBox(res.existing_with_structure || 0, '✅已有·含结构', '') +
      statBox(res.new_no_structure || 0, '🆕新·需补结构', 'la-warn') +
      '</div>' +
      (res.dedup_removed ? '<div class="la-muted">去重合并: ' + res.dedup_removed + ' 条重复</div>' : '') +
      ((tierCounts.literature_presence_only || 0) ? '<div class="la-warn-card"><b>无实测碎片强度：</b>' + tierCounts.literature_presence_only + ' 条仅作为文献靶标/诊断离子线索，不能获得实验谱 Level 2a/High，且不会计作已确认化合物。</div>' : '') +
      ((res.provenance_held_records || 0) ? '<div class="la-warn-card"><b>数值强度不等于实验证明：</b>' + res.provenance_held_records + ' 条虽含强度，但缺少可核验的仪器/来源/DOI/accession，已标记为 user_unverified 并暂存，不进入检索 Top-K。</div>' : '') +
      ((quality.rejected || 0) ? '<div class="la-error">已隔离 ' + quality.rejected + ' 条内部不一致记录。主要错误：' + esc(JSON.stringify(quality.error_counts || {})) + '</div>' : '') +
      conv + '</div>';

    var samples = res.samples || [];
    var rows = samples.slice(0, 50).map(function (s) {
      var status;
      if (s.is_new) {
        status = s.has_smiles
          ? '<span class="la-badge la-badge-new">🆕新·有结构</span>'
          : '<span class="la-badge la-badge-warn">🆕新·无结构(不自动借用)</span>';
      } else {
        status = (s.db_has_smiles === true)
          ? '<span class="la-badge la-badge-exist">✅已有·库中含结构</span>'
          : (s.db_has_smiles === false)
            ? '<span class="la-badge la-badge-warn">✅已有·库中也无结构</span>'
            : '<span class="la-badge la-badge-exist">✅已有</span>';
      }
      return '<tr><td>' + esc(s.name) + '</td><td><code>' + esc(s.formula) + '</code></td>' +
        '<td>' + (s.precursor_mz ? (+s.precursor_mz).toFixed(4) : '-') + '</td>' +
        '<td>' + s.n_peaks + '</td><td>' + status + '</td></tr>';
    }).join('');
    var table = '<div class="la-card">' +
      '<h4>样本明细（前 ' + Math.min(50, samples.length) + ' / ' + samples.length + '）</h4>' +
      '<div class="la-table-wrap"><table class="la-table"><thead><tr>' +
      '<th>化合物名</th><th>Formula</th><th>Precursor m/z</th><th>碎片数</th><th>状态</th>' +
      '</tr></thead><tbody>' + rows + '</tbody></table></div></div>';

    area.innerHTML = summary + table;
  }

  function statBox(num, lbl, cls) {
    return '<div class="la-stat ' + (cls || '') + '"><div class="la-stat-num">' + num + '</div><div class="la-stat-lbl">' + lbl + '</div></div>';
  }

  async function doImport() {
    if (!stagedFiles.length) { if (typeof Toast !== 'undefined') Toast.show('请先选择文件', 'warning'); return; }
    var q = (lastPreview && lastPreview.quality) || {};
    var importable = Number((lastPreview && lastPreview.importable_records) || 0);
    var spectra = Number((lastPreview && lastPreview.spectral_records) || 0);
    var targets = Number((lastPreview && lastPreview.literature_target_records) || 0);
    var held = Number((lastPreview && lastPreview.provenance_held_records) || 0);
    if (!importable) { if (typeof Toast !== 'undefined') Toast.show('没有可导入的实验谱或文献靶标，数据库不会修改', 'warning'); return; }
    if (!confirm('确认导入 ' + importable + ' 条记录到 Full + Small 双库？\n' +
      '其中可检索实验谱 ' + spectra + ' 条、文献靶标 ' + targets + ' 条；' +
      held + ' 条来源待核验记录不会进入检索，' + (q.rejected || 0) + ' 条不一致记录将被隔离。\n' +
      '图片数字化峰表不会直接进入普通实验谱检索。')) return;
    setBusy(true, '⏳ 导入中（追加到双库，可能需要数秒到数十秒）...');
    var fd = new FormData();
    stagedFiles.forEach(function (f) { fd.append('files', f); });
    try {
      var res = await API.post('/library/import?db_mode=both', fd);
      renderImportResult(res);
    } catch (e) {
      if (typeof Toast !== 'undefined') Toast.show('导入失败: ' + e.message, 'error');
      var area = document.getElementById('libaugmentResult');
      if (area) area.innerHTML = '<div class="la-error">❌ 导入失败: ' + esc(e.message) + '</div>';
    } finally {
      var p = document.getElementById('libaugmentPreviewBtn');
      var im = document.getElementById('libaugmentImportBtn');
      if (p) p.disabled = false;
      if (im) im.disabled = false;
    }
  }

  function renderImportResult(res) {
    var area = document.getElementById('libaugmentResult');
    var stats = '';
    if (res.full) stats += statBox(res.full.new_compounds || 0, 'Full 新化合物', 'la-new') + statBox(res.full.new_spectra || 0, 'Full 新谱图', '');
    if (res.small) stats += statBox(res.small.new_compounds || 0, 'Small 新化合物<br>(&lt;800Da)', 'la-new') + statBox(res.small.new_spectra || 0, 'Small 新谱图', '');
    var html = '<div class="la-card la-success">' +
      '<h4>✅ 导入完成</h4>' +
      '<div class="la-stat-grid">' + stats + '</div>' +
      '<div class="la-muted">实际导入 ' + (res.importable_records || 0) + ' 条（实验谱 ' + (res.spectral_records || 0) + ' · 文献靶标 ' + (res.literature_target_records || 0) + '）· 来源待核验暂存 ' + (res.provenance_held_records || 0) + ' · 去重 ' + (res.dedup_removed || 0) + ' 条</div>' +
      (res.quality ? '<div class="la-muted">质量门控：通过 ' + (res.quality.accepted || 0) + ' · 隔离 ' + (res.quality.rejected || 0) + '</div>' : '') +
      '<div class="la-warn-card">⚠️ 检索缓存已失效。请前往「数据上传」页重新加载数据库（Full 或 Small），新导入的化合物才会在检索中命中。</div>' +
      '<button class="btn btn-primary mt-1" onclick="LibAugmentPage.goReloadLib()">🔄 去重新加载数据库</button>' +
      '</div>';
    area.innerHTML = html;
    if (typeof Toast !== 'undefined') Toast.show('导入成功！', 'success');
    // 库已变化 → 刷新数据上传页下拉的真实规模数字
    if (typeof Upload !== 'undefined' && Upload.refreshDbCounts) Upload.refreshDbCounts();
  }

  function goReloadLib() {
    if (typeof App !== 'undefined' && App.navigateTo) App.navigateTo('upload');
  }

  function value(id) {
    var element = document.getElementById(id);
    return element ? String(element.value || '').trim() : '';
  }

  function numericValue(id) {
    var raw = value(id);
    if (!raw) return null;
    var number = Number(raw);
    return Number.isFinite(number) ? number : null;
  }

  function parseDigitizedPeaks(text) {
    return String(text || '').replace(/;/g, '\n').split(/\r?\n/).map(function (line) {
      var parts = line.trim().split(/[\s,\t]+/);
      return parts.length >= 2 ? [Number(parts[0]), Number(parts[1])] : null;
    }).filter(function (peak) {
      return peak && Number.isFinite(peak[0]) && Number.isFinite(peak[1]) && peak[0] > 0 && peak[1] > 0;
    }).sort(function (left, right) { return left[0] - right[0]; });
  }

  function referenceManifestPayload(sourceName, includeExample) {
    var extension = sourceFileExtension(sourceName);
    var isImage = ['.png', '.jpg', '.jpeg', '.tif', '.tiff'].indexOf(extension) !== -1;
    var spectrum = {
      spectrum_id: value('laRefSpectrumId') || undefined,
      compound: {
        name: value('laRefName'),
        formula: value('laRefFormula'),
        smiles: value('laRefSmiles'),
        inchikey: value('laRefInchikey'),
      },
      precursor_mz: numericValue('laRefPrecursor'),
      adduct: value('laRefAdduct'),
      ion_mode: value('laRefPolarity'),
      charge: numericValue('laRefCharge'),
      instrument_metadata: {
        model: value('laRefInstrument'),
        analyzer: value('laRefAnalyzer'),
        activation: value('laRefActivation'),
        collision_cell_region: value('laRefCollisionRegion'),
        acquisition_method: value('laRefMethod'),
      },
      collision_energy: {
        value: numericValue('laRefCe'),
        minimum: numericValue('laRefCeMin'),
        maximum: numericValue('laRefCeMax'),
        unit: value('laRefCeUnit'),
        mode: value('laRefCeMode'),
      },
      reference_standard: {
        material_type: value('laRefMaterial'),
        lot: value('laRefLot'),
        purity_percent: numericValue('laRefPurity'),
      },
      provenance: {
        verified: !!document.getElementById('laRefVerified')?.checked,
        source_accession: value('laRefAccession'),
        doi: value('laRefDoi'),
      },
      source_provenance_verified: !!document.getElementById('laRefVerified')?.checked,
    };
    if (isImage) {
      spectrum.image_file = sourceName;
      spectrum.peaks = includeExample
        ? [[184.0733, 100], [104.1070, 35], [86.0964, 18]]
        : parseDigitizedPeaks(value('laRefPeaks'));
      spectrum.digitization = {
        method: 'manual_axis_calibrated_peak_pick',
        reviewer_confirmed: !!document.getElementById('laRefVerified')?.checked,
      };
      spectrum.intensity_origin = 'digitized_image';
    } else {
      spectrum.source_file = sourceName || 'reference_standard.mzML';
    }
    return {
      schema_version: 'mssystem-reference-intake/1.0',
      manifest_id: 'v15-reference-' + new Date().toISOString().replace(/[^0-9]/g, '').slice(0, 14),
      spectra: [spectrum],
    };
  }

  function addReferenceManifest() {
    var sourceName = value('laRefSource');
    if (!sourceName) {
      if (typeof Toast !== 'undefined') Toast.show('请先在上方加入并选择 MGF、MSP、mzML 或谱图图片', 'warning');
      return;
    }
    if (!value('laRefName') || !value('laRefFormula') || !value('laRefAdduct')) {
      if (typeof Toast !== 'undefined') Toast.show('化合物名称、分子式和加合物不能为空', 'warning');
      return;
    }
    if (!value('laRefInstrument') || !value('laRefActivation') || numericValue('laRefCe') == null || !value('laRefCeUnit') || !value('laRefCeMode')) {
      if (typeof Toast !== 'undefined') Toast.show('仪器型号、活化方式、CE 数值、CE 单位和 CE 模式必须明确', 'warning');
      return;
    }
    if (numericValue('laRefPrecursor') == null || numericValue('laRefPrecursor') <= 0) {
      if (typeof Toast !== 'undefined') Toast.show('前体 m/z 必须是正数', 'warning');
      return;
    }
    var ceMode = value('laRefCeMode');
    if (ceMode !== 'fixed' && (numericValue('laRefCeMin') == null || numericValue('laRefCeMax') == null)) {
      if (typeof Toast !== 'undefined') Toast.show('stepped/ramp 谱必须保留 CE 最小值和最大值；系统不会用中点代替整张复合谱', 'warning');
      return;
    }
    var material = value('laRefMaterial');
    if ((material === 'authentic_standard' || material === 'certified_reference_material') &&
        (!value('laRefLot') || numericValue('laRefPurity') == null || !value('laRefMethod'))) {
      if (typeof Toast !== 'undefined') Toast.show('真实对照品记录必须填写批号、纯度和采集方法 ID', 'warning');
      return;
    }
    if (material === 'public_experimental_spectrum' && !value('laRefAccession') && !value('laRefDoi')) {
      if (typeof Toast !== 'undefined') Toast.show('公共实验谱必须填写 accession 或 DOI', 'warning');
      return;
    }
    if (!document.getElementById('laRefVerified')?.checked) {
      if (typeof Toast !== 'undefined') Toast.show('请先完成来源与采集条件人工核对', 'warning');
      return;
    }
    var extension = sourceFileExtension(sourceName);
    if (['.png', '.jpg', '.jpeg', '.tif', '.tiff'].indexOf(extension) !== -1 && parseDigitizedPeaks(value('laRefPeaks')).length < 3) {
      if (typeof Toast !== 'undefined') Toast.show('图片记录至少需要 3 个经复核的 m/z–intensity 数值', 'warning');
      return;
    }
    var payload = referenceManifestPayload(sourceName, false);
    var blob = new Blob([JSON.stringify(payload, null, 2)], { type: 'application/json' });
    var filename = 'v15_reference_manifest_' + Date.now() + '.json';
    stagedFiles.push(new File([blob], filename, { type: 'application/json' }));
    renderFileList();
    if (typeof Toast !== 'undefined') Toast.show('V15 参考谱清单已加入待预览文件', 'success');
  }

  function downloadManifestTemplate() {
    var payload = referenceManifestPayload(value('laRefSource') || 'reference_standard.mzML', true);
    var blob = new Blob([JSON.stringify(payload, null, 2)], { type: 'application/json' });
    var link = document.createElement('a');
    link.href = URL.createObjectURL(blob);
    link.download = 'mssystem_v15_reference_intake_template.json';
    document.body.appendChild(link);
    link.click();
    link.remove();
    URL.revokeObjectURL(link.href);
  }

  function esc(s) {
    return String(s == null ? '' : s).replace(/[&<>"]/g, function (c) {
      return { '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c];
    });
  }

  return {
    init: init,
    removeFile: removeFile,
    goReloadLib: goReloadLib,
    addReferenceManifest: addReferenceManifest,
    downloadManifestTemplate: downloadManifestTemplate,
  };
})();
