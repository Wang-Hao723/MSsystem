/**
 * MSsystem global language runtime.
 *
 * Locale resolution order:
 *   1. ?lang=en / ?lang=zh URL override
 *   2. the user's saved choice
 *   3. the browser language
 *   4. English fallback
 *
 * The runtime intentionally translates interface copy only. Analytical values,
 * filenames, compound names, formulae and identifiers are never sent to a
 * translation service and are not modified by heuristic machine translation.
 * Renderers mark scientific values with ``data-scientific-value``; that marker
 * is a hard translation firewall even when the source value contains Han text.
 */
(function installMssystemLanguage(global) {
  'use strict';

  const STORAGE_KEY = 'mssystem.locale.v1';
  const DEFAULT_LOCALE = 'en-US';
  const SUPPORTED_LOCALES = Object.freeze(['zh-CN', 'en-US']);
  const HAN_RE = /[\u3400-\u4dbf\u4e00-\u9fff\uf900-\ufaff]/;
  const catalog = global.MSS_I18N_CATALOG || { messages: {}, legacy: {} };
  const textRecords = new WeakMap();
  const attributeRecords = new WeakMap();
  const translatedElements = new WeakSet();
  const missing = new Set();
  const translatableAttributes = Object.freeze([
    'placeholder', 'title', 'aria-label', 'alt', 'label', 'data-empty-text',
  ]);
  const skipSelector = [
    'script', 'style', 'code', 'pre', 'textarea', 'template',
    '[contenteditable="true"]', '[data-i18n-skip]', '[data-scientific-value]',
  ].join(',');
  const skipElementSelector = [
    'script', 'style', 'template', '[data-i18n-skip]', '[data-scientific-value]',
  ].join(',');

  let locale = DEFAULT_LOCALE;
  let localeSource = 'fallback';
  let observer = null;
  let observerRoot = null;
  let observerOptions = null;
  let translationQueued = false;
  const pendingRoots = new Set();

  function normalizeLocale(value) {
    const raw = String(value || '').trim().replace('_', '-').toLowerCase();
    if (!raw || raw === 'auto' || raw === 'system' || raw === 'browser') return null;
    if (raw === 'zh' || raw.startsWith('zh-')) return 'zh-CN';
    if (raw === 'en' || raw.startsWith('en-')) return 'en-US';
    return null;
  }

  function languageCode(value = locale) {
    return normalizeLocale(value) === 'zh-CN' ? 'zh' : 'en';
  }

  function safeStorageGet() {
    try { return global.localStorage?.getItem(STORAGE_KEY) || ''; }
    catch { return ''; }
  }

  function browserLocale() {
    const candidates = [];
    if (Array.isArray(global.navigator?.languages)) {
      candidates.push(...global.navigator.languages);
    }
    candidates.push(global.navigator?.language);
    for (const candidate of candidates) {
      const normalized = normalizeLocale(candidate);
      if (normalized) return normalized;
    }
    return DEFAULT_LOCALE;
  }

  function urlLocale() {
    try {
      return normalizeLocale(new URL(global.location.href).searchParams.get('lang'));
    } catch {
      return null;
    }
  }

  function resolveLocale() {
    const fromUrl = urlLocale();
    if (fromUrl) return { locale: fromUrl, source: 'url' };
    const fromStorage = normalizeLocale(safeStorageGet());
    if (fromStorage) return { locale: fromStorage, source: 'saved' };
    const fromBrowser = browserLocale();
    if (fromBrowser) return { locale: fromBrowser, source: 'browser' };
    return { locale: DEFAULT_LOCALE, source: 'fallback' };
  }

  function interpolate(message, params) {
    if (!params) return String(message);
    return String(message).replace(/\{([A-Za-z0-9_]+)\}/g, (match, key) => (
      Object.prototype.hasOwnProperty.call(params, key) ? String(params[key]) : match
    ));
  }

  function preserveWhitespace(source, translated) {
    const leading = String(source).match(/^\s*/)?.[0] || '';
    const trailing = String(source).match(/\s*$/)?.[0] || '';
    return `${leading}${translated}${trailing}`;
  }

  function normalizedText(value) {
    return String(value || '').replace(/\s+/g, ' ').trim();
  }

  function translatePattern(source, targetLocale) {
    if (targetLocale !== 'en-US') return null;
    if (source.includes('，推测丢失／转化 ')) {
      const fragmentation = source
        .replace(/，推测丢失／转化 /g, ', proposed loss/transformation: ')
        .replace(/；/g, '; ')
        .replace(/甲基自由基/g, 'methyl radical')
        .replace(/二氧化碳/g, 'carbon dioxide')
        .replace(/一氧化碳/g, 'carbon monoxide')
        .replace(/水/g, 'water')
        .replace(/氨/g, 'ammonia')
        .replace(/甲醇/g, 'methanol')
        .replace(/甲醛/g, 'formaldehyde')
        .replace(/乙烯/g, 'ethylene')
        .replace(/戊糖残基/g, 'pentose residue')
        .replace(/己糖残基/g, 'hexose residue')
        .replace(/苯氧基/g, 'phenoxy group');
      if (!HAN_RE.test(fragmentation)) return fragmentation;
    }
    if (source.startsWith('分析参数:')) {
      const structured = source
        .replace(/^分析参数:/, 'Analysis parameters:')
        .replace(/(^|\s)库:/g, '$1Library:')
        .replace(/Cosine 阈值:/g, 'Cosine threshold:')
        .replace(/Final 阈值:/g, 'Final threshold:')
        .replace(/证据模式:/g, 'Evidence mode:')
        .replace(/严格门槛/g, 'Strict gate')
        .replace(/常规门槛/g, 'Routine gate')
        .replace(/历史模式（未记录）/g, 'Legacy mode (not recorded)')
        .replace(/m\/z 容差:/g, 'm/z tolerance:')
        .replace(/TopN 峰:/g, 'Top-N peaks:')
        .replace(/嵌合去卷积: 开/g, 'Chimeric deconvolution: On')
        .replace(/嵌合去卷积: 关/g, 'Chimeric deconvolution: Off')
        .replace(/IIMN: 开/g, 'IIMN: On')
        .replace(/IIMN: 关/g, 'IIMN: Off')
        .replace(/隔离窗全宽:/g, 'Full isolation-window width:')
        .replace(/（显式覆盖缺失元数据）/g, ' (explicit override for missing metadata)')
        .replace(/空白特征过滤:/g, 'Blank-feature filtering:')
        .replace(/([\d,]+) 个空白/g, '$1 blanks')
        .replace(/样品\/空白≥/g, 'sample/blank ≥ ');
      if (!HAN_RE.test(structured)) return structured;
    }
    const registeredPatterns = catalog.patterns?.[targetLocale] || [];
    for (const [pattern, replacement] of registeredPatterns) {
      if (!(pattern instanceof RegExp)) continue;
      pattern.lastIndex = 0;
      if (!pattern.test(source)) continue;
      pattern.lastIndex = 0;
      return source.replace(pattern, replacement);
    }
    const patterns = [
      [/^--\s*选择已完成的任务\s*--$/, '-- Select a completed run --'],
      [/^--\s*选择已完成任务\s*--$/, '-- Select a completed run --'],
      [/^--\s*选择任务\s*--$/, '-- Select a run --'],
      [/^(\d+)\s*个样本$/, '$1 samples'],
      [/^(\d+)\s*个任务$/, '$1 runs'],
      [/^(\d+)\s*个候选$/, '$1 candidates'],
      [/^(\d+)\s*个谱图$/, '$1 spectra'],
      [/^(\d+)\s*个$/, '$1'],
      [/^(\d+)\/(\d+)\s*匹配$/, '$1/$2 matched'],
      [/^已选\s*(\d+)\s*项$/, '$1 selected'],
      [/^上传进度:\s*(\d+)%$/, 'Upload progress: $1%'],
      [/^第\s*(\d+)\s*页$/, 'Page $1'],
      [/^共\s*(\d+)\s*页$/, '$1 pages'],
      [/^共\s*(\d+)\s*条$/, '$1 records'],
      [/^正在打包\s*(\d+)\s*个样本…?$/, 'Packaging $1 samples…'],
      [/^正在生成报告\s*(\d+)\/(\d+)…?$/, 'Generating reports $1/$2…'],
      [/^正在压缩\s*(\d+)\/(\d+)…?$/, 'Compressing $1/$2…'],
      [/^([\d,]+)\s*个任务记录$/, '$1 run records'],
      [/^([\d,]+)\s*张可检索谱图$/, '$1 searchable spectra'],
      [/^([\d,]+)\s*行\s*·\s*([\d,]+)\s*个名称\/分子式实体$/, '$1 rows · $2 name/formula entities'],
      [/^(POS|NEG|POS \+ NEG)\s*·\s*已完成\s*·\s*([\d,]+)\/([\d,]+)\s*匹配$/, '$1 · Completed · $2/$3 matched'],
      [/^—\s*([\d,]+)\s*张谱图已加载\s*\(([\d.]+)s\)$/, '— $1 spectra loaded ($2s)'],
      [/^Full\s*\(([\d,]+)\)\s*已加载，可开始谱库分析。$/, 'Full ($1) is loaded. Library analysis can start.'],
      [/^(Full|Small <800 Da)\s*\(([^/]+)\s*当前已载入\s*\/\s*([^/]+)\s*来源合格\s*\/\s*([^；]+)\s*物理记录；([^\)]+)\s*化合物\)$/, '$1 ($2 currently loaded / $3 source-qualified / $4 physical records; $5 compounds)'],
      [/^(Full|Small <800 Da)\s*\(([^/]+)\s*运行时可检索\s*\/\s*([^/]+)\s*来源合格\s*\/\s*([^；]+)\s*物理记录；([^\)]+)\s*化合物\)$/, '$1 ($2 runtime-searchable / $3 source-qualified / $4 physical records; $5 compounds)'],
      [/^全选可删除（([\d,]+)）$/, 'Select all removable ($1)'],
      [/^已选\s*([\d,]+)$/, '$1 selected'],
      [/^🗑\s*删除所选\s*\(([\d,]+)\)$/, '🗑 Remove Selected ($1)'],
      [/^尚未选择目标库\s*·\s*靶向层不会运行（可用\s*([\d,]+)\s*个）$/, 'No target library selected · the targeted layer will not run ($1 available)'],
      [/^选择删除\s+(.+)$/, 'Select for deletion: $1'],
      [/^([\d,]+)\s*个去重化合物$/, '$1 deduplicated compounds'],
      [/^([\d,]+)\s*个主要化合物结果$/, '$1 main compound results'],
      [/^([\d,]+)\s*个未解析事件（未计数）$/, (_, count) => `${count} unresolved event${count === '1' ? '' : 's'} (excluded from count)`],
      [/^([\d,]+)\s*个谱库候选身份$/, '$1 library candidate identities'],
      [/^🧪\s*已确认背景排除\s*([\d,]+)\s*条$/, '🧪 Confirmed background exclusions: $1'],
      [/^🧪\s*空白仅比较\s*·\s*([\d,]+)\s*条被扣除$/, '🧪 Blank comparison only · $1 excluded'],
      [/^↔\s*按极性分配\s*(POS|NEG|POS\/NEG|未知极性)\s*空白\s*([\d,]+)\s*个（不等于方法配对确认）(?:，跳过异极性\s*([\d,]+)\s*个)?$/, (_, polarity, assigned, skipped) => `↔ ${assigned} ${polarity} blank${assigned === '1' ? '' : 's'} assigned by polarity (method pairing not confirmed)${skipped ? `; ${skipped} opposite-polarity blank${skipped === '1' ? '' : 's'} skipped` : ''}`],
      [/^↔\s*按极性分配\s*(POS|NEG|POS\/NEG|未知极性)\s*空白\s*([\d,]+)\s*个（方法可比性与背景排除已确认）(?:，跳过异极性\s*([\d,]+)\s*个)?$/, (_, polarity, assigned, skipped) => `↔ ${assigned} ${polarity} blank${assigned === '1' ? '' : 's'} assigned; method comparability and background exclusion confirmed${skipped ? `; ${skipped} opposite-polarity blank${skipped === '1' ? '' : 's'} skipped` : ''}`],
      [/^(Full|Small)\s*·\s*([\d.]+)\/([\d.]+)\s*·\s*严格门槛$/, '$1 · $2/$3 · Strict gate'],
      [/^(Full|Small)\s*·\s*([\d.]+)\/([\d.]+)\s*·\s*常规门槛$/, '$1 · $2/$3 · Routine gate'],
      [/^可检索谱图：([\d,]+)$/, 'Searchable spectra: $1'],
      [/^([\d,]+)\s*条正离子谱图$/, '$1 positive-ion spectra'],
      [/^([\d,]+)\s*条负离子谱图$/, '$1 negative-ion spectra'],
      [/^匹配空白特征\s*([\d,]+)\s*条；未跨针相减\s*MS\/MS$/, '$1 matched blank features; no cross-injection MS/MS subtraction'],
      [/^使用：(.+?)；跳过：(.+)$/, 'Used: $1; skipped: $2'],
      [/^使用：(.+)$/, 'Used: $1'],
      [/^(.+)\s+\((\d+)\/(\d+)\s*匹配,\s*([\d.]+)%\)$/, '$1 ($2/$3 matched, $4%)'],
      [/^(.+)\s+\((\d+)\/(\d+)\s*匹配\)$/, '$1 ($2/$3 matched)'],
      [/^(.+)\s+·\s+(\d+)\s+MS1 features\s+·\s+未运行$/, '$1 · $2 MS1 features · Not run'],
      [/^(.+?)\s+·\s+([0-9-]+\s+[0-9:]+)\s+·\s+无空白证据\s+·\s+未绑定靶向库\s+·\s+([\d,]+)\/([\d,]+)\s+MS2\s+有候选$/, '$1 · $2 · No blank evidence · No targeted library · $3/$4 MS2 with candidates'],
      [/^(.+?)\s+·\s+([0-9-]+\s+[0-9:]+)\s+·\s+已确认背景空白\s+([\d,]+)（扣除\s+([\d,]+)\s+条）\s+·\s+未绑定靶向库\s+·\s+([\d,]+)\/([\d,]+)\s+MS2\s+有候选$/, '$1 · $2 · $3 confirmed background blanks ($4 excluded) · No targeted library · $5/$6 MS2 with candidates'],
      [/^(.+?)\s+·\s+([0-9-]+\s+[0-9:]+)\s+·\s+(.+?)\s+·\s+([\d,]+)\/([\d,]+)\s+MS2\s+有候选$/, '$1 · $2 · $3 · $4/$5 MS2 with candidates'],
      [/^(.+?)\s+·\s+([\d,]+)\/([\d,]+)\s+MS2\s+有候选$/, '$1 · $2/$3 MS2 with candidates'],
      [/^(POS|NEG)\s+离子模式$/, '$1 ion mode'],
      [/^(TIC|BPC)\s+·\s+([\d,]+)\s+张连续\s+MS1$/, '$1 · $2 consecutive MS1 scans'],
      [/^([\d,]+)\s+张连续\s+MS1$/, '$1 consecutive MS1 scans'],
      [/^主要\s+(TIC|BPC)\s+峰\s+([\d,]+)\s+·\s+RT\s+([\d.]+)\s+min\s+·\s+附近采集\s+([\d,]+)\s+张\s+MS2，但没有候选通过当前谱库与分数门槛；不代表该峰没有化学信号$/, 'Major $1 peak $2 · RT $3 min · $4 nearby MS2 scans; no candidate passed the current library and score thresholds. This does not imply absence of a chemical signal.'],
      [/^主要\s+(TIC|BPC)\s+峰\s+([\d,]+)\s+·\s+RT\s+([\d.]+)\s+min\s+·\s+附近采集\s+([\d,]+)\s+张\s+MS2，其中\s+([\d,]+)\s+张获得候选，但最高仅为\s+Low\/Unknown\s+线索；不代表该峰没有化学信号$/, 'Major $1 peak $2 · RT $3 min · $4 nearby MS2 scans; $5 produced a candidate, but the strongest result is only a Low/Unknown lead. This does not imply absence of a chemical signal.'],
      [/^(中等证据候选|高可信异构体组检出|中等证据异构体组|高可信单一二维连接结构检出)：(.+?)\s+·\s+RT\s+([\d.]+)\s+min\s+·\s+m\/z\s+([\d.]+)\s+·\s+([\d,]+)\s+张局部\s+MS2$/, (_, label, name, rt, mz, scans) => {
        const labels = {
          '中等证据候选': 'Moderate-evidence candidate',
          '高可信异构体组检出': 'High-confidence isomer-group detection',
          '中等证据异构体组': 'Moderate-evidence isomer group',
          '高可信单一二维连接结构检出': 'High-confidence single 2-D connectivity detection',
        };
        return `${labels[label]}: ${name} · RT ${rt} min · m/z ${mz} · ${scans} local MS2 scan${scans === '1' ? '' : 's'}`;
      }],
      [/^主要结果检测位置\s+([\d,]+)\s+·\s+暂无主要结果邻近的主要峰\s+([\d,]+)$/, '$1 primary-result acquisition locations · $2 major peaks without a nearby primary result'],
      [/^([\d,]+)\s+个\s+OpenMS\s+特征$/, '$1 OpenMS features'],
      [/^([\d,]+)\s+个特征关联\s+([\d,]+)\s+张\s+MS2$/, '$1 features linked to $2 MS2 scans'],
      [/^([\d,]+)\s+个主要结果$/, '$1 primary results'],
      [/^([\d,]+)\s+个化合物／异构体组线索，含背景与低证据项$/, '$1 compound/isomer-group leads, including background-related and low-evidence entries'],
      [/^([\d,]+)\s+可靠分类\s+·\s+([\d,]+)\s+类别级生源背景$/, '$1 reliable classifications · $2 class-level biosource contexts'],
      [/^([\d,]+)\s+个结果重建了样品碎裂；([\d,]+)\s+个有谱库参考$/, '$1 results with sample fragmentation reconstruction; $2 with library references'],
      [/^([\d,]+)\s+个\s+IIMN$/, '$1 IIMN relationships'],
      [/^([\d,]+)\s+张\s+MS2$/, '$1 MS2 scans'],
      [/^([\d,]+)\s+张仍有可用碎片$/, '$1 still contain usable fragments'],
      [/^([\d,]+)\s+可靠\s+·\s+([\d,]+)\s+待核查$/, '$1 reliable · $2 requiring review'],
      [/^([\d,]+)\s*\/\s*([\d,]+)\s+候选$/, '$1 / $2 candidates'],
      [/^([\d,]+)\s+个结果有正向规则支持$/, '$1 results with positive rule support'],
      [/^([\d,]+)\s+个结果有样品步骤\/产物$/, '$1 results with sample-derived steps/products'],
      [/^([\d,]+)\s+张当前没有谱库候选，其中\s+([\d,]+)\s+张仍有可用碎片。未匹配只表示当前谱库与阈值下没有足够相似的参比谱。$/, '$1 scans currently have no library candidate; $2 still contain usable fragments. Unmatched means only that no sufficiently similar reference spectrum was found under the current library and thresholds.'],
      [/^显示\s+([\d,]+)\s*\/\s*([\d,]+)\s+个去重结果$/, 'Showing $1 / $2 deduplicated results'],
      [/^张\s+MS2，已合并为\s+([\d,]+)\s+个色谱峰区$/, ' MS2 scans, consolidated into $1 chromatographic peak region(s)'],
      [/^查看全部\s+([\d,]+)\s+个竞争候选及证据\s+→$/, 'View all $1 competing candidates and evidence →'],
      [/^本组\s+([\d,]+)\s+个谱学竞争候选；其余名称不再挤在概览中$/, '$1 spectrally competing candidates in this group; remaining names are available in the detailed view'],
      [/^([\d,]+)\s+个观测峰支持的候选产物结构；([\d,]+)\s+个质量差步骤$/, '$1 observed peaks support candidate product structures; $2 mass-difference steps'],
      [/^([\d,]+)\s+个观测峰支持的候选产物结构$/, '$1 observed peaks support candidate product structures'],
      [/^有支持证据\s+([\d,]+)\s+维；已评价未观察\s+([\d,]+)\s+维；不适用\s+([\d,]+)\s+维；未建立\s+([\d,]+)\s+维$/, '$1 supported dimensions; $2 evaluated but not observed; $3 not applicable; $4 not established'],
      [/^综合排序分\s+([\d.]+)；规则证据、分类、生源和碎裂解释不参与身份排序。$/, 'Composite ranking score $1; rule evidence, classification, biosource context and fragmentation interpretation do not contribute to identity ranking.'],
      [/^([\d,]+)\s+个色谱峰区\s+·\s+([\d,]+)\s+张\s+MS2$/, '$1 chromatographic peak region(s) · $2 MS2 scans'],
      [/^代表谱\s+(.+)$/, 'Representative spectrum $1'],
      [/^候选离子假设\s+(.+)（非独立证据）$/, 'Candidate ion-form hypothesis $1 (not independent evidence)'],
      [/^类别待核查（模型提示：(.+)）$/, 'Class requires review (model suggestion: $1)'],
      [/^分类来源：(.+)；分类可靠性与谱库身份匹配等级独立评估$/, 'Classification source: $1; classification reliability is evaluated independently of library identity-match level'],
      [/^(总离子流|基峰)色谱图（(TIC|BPC)），保留时间\s+([\d.]+)\s+至\s+([\d.]+)\s+分钟；标出\s+([\d,]+)\s+个已确证或高可信检测位置、([\d,]+)\s+个异构体组检测位置、([\d,]+)\s+个中等证据候选检测位置和\s+([\d,]+)\s+个(暂无重点结果邻近的主要峰|未建立该基峰离子 MS2 结果的主要峰)$/, (_, kind, abbr, start, end, high, isomer, medium, unassigned, unresolvedLabel) => `${kind === '总离子流' ? 'Total ion' : 'Base-peak'} chromatogram (${abbr}), retention time ${start}–${end} min; ${high} confirmed or high-confidence locations, ${isomer} isomer-group locations, ${medium} moderate-evidence candidate locations and ${unassigned} ${unresolvedLabel.startsWith('未建立') ? 'major peaks without MS2 evidence for the base-peak ion' : 'major peaks without a nearby priority result'}`],
    ];
    for (const [pattern, replacement] of patterns) {
      if (pattern.test(source)) return source.replace(pattern, replacement);
    }
    return null;
  }

  function translateDualLabel(source, targetLocale) {
    if (targetLocale !== 'en-US' || !HAN_RE.test(source)) return null;
    for (const delimiter of [' · ', ' / ']) {
      const parts = source.split(delimiter);
      if (parts.length === 2 && HAN_RE.test(parts[0]) && !HAN_RE.test(parts[1])) {
        return parts[1].trim();
      }
    }
    return null;
  }

  function lookupMessage(key, targetLocale = locale, params) {
    const localeMessages = catalog.messages?.[targetLocale] || {};
    const fallbackMessages = catalog.messages?.[DEFAULT_LOCALE] || {};
    const value = localeMessages[key] ?? fallbackMessages[key];
    return value == null ? null : interpolate(value, params);
  }

  function translateText(value, params, targetLocale = locale) {
    const source = normalizedText(value);
    if (!source) return String(value || '');

    const semantic = lookupMessage(source, targetLocale, params);
    if (semantic != null) return semantic;

    if (targetLocale === 'zh-CN') return interpolate(source, params);

    const exact = catalog.legacy?.[targetLocale]?.[source];
    if (exact != null) return interpolate(exact, params);

    const pattern = translatePattern(source, targetLocale);
    if (pattern != null) return interpolate(pattern, params);

    const dual = translateDualLabel(source, targetLocale);
    if (dual != null) return interpolate(dual, params);

    if (HAN_RE.test(source)) missing.add(source);
    return interpolate(source, params);
  }

  function t(keyOrSource, params) {
    return translateText(keyOrSource, params, locale);
  }

  function shouldSkip(element) {
    if (!element || element.nodeType !== 1) return false;
    return Boolean(element.closest?.(skipSelector));
  }

  function renderTextNode(node) {
    if (!node || node.nodeType !== 3 || shouldSkip(node.parentElement)) return;
    const current = node.nodeValue || '';
    const compact = normalizedText(current);
    if (!compact || (!HAN_RE.test(compact) && !textRecords.has(node))) return;

    let record = textRecords.get(node);
    if (!record || current !== record.rendered) {
      record = { source: current, rendered: current };
      textRecords.set(node, record);
    }
    const sourceCompact = normalizedText(record.source);
    const translated = translateText(sourceCompact, null, locale);
    const next = preserveWhitespace(record.source, translated);
    record.rendered = next;
    if (node.nodeValue !== next) node.nodeValue = next;
  }

  function attributeRecord(element) {
    let record = attributeRecords.get(element);
    if (!record) {
      record = new Map();
      attributeRecords.set(element, record);
    }
    return record;
  }

  function renderAttribute(element, name, semanticKey = null) {
    if (!element?.hasAttribute?.(name) || shouldSkip(element)) return;
    const current = element.getAttribute(name) || '';
    const records = attributeRecord(element);
    let record = records.get(name);
    if (!record || current !== record.rendered) {
      record = { source: current, rendered: current };
      records.set(name, record);
    }
    const translated = semanticKey
      ? (lookupMessage(semanticKey, locale) ?? record.source)
      : translateText(record.source, null, locale);
    record.rendered = translated;
    if (current !== translated) element.setAttribute(name, translated);
  }

  function renderElement(element) {
    if (!element || element.nodeType !== 1 || element.closest?.(skipElementSelector)) return;
    const semanticKey = element.getAttribute('data-i18n');
    if (semanticKey) {
      const translated = lookupMessage(semanticKey, locale);
      if (translated != null) {
        let record = textRecords.get(element);
        if (!record) {
          record = { source: element.textContent, rendered: element.textContent };
          textRecords.set(element, record);
        }
        record.rendered = translated;
        if (element.textContent !== translated) element.textContent = translated;
        translatedElements.add(element);
      }
    }
    translatableAttributes.forEach(name => {
      const key = element.getAttribute(`data-i18n-${name}`);
      if (element.hasAttribute(name)) renderAttribute(element, name, key);
    });
  }

  function walk(root) {
    if (!root) return;
    if (root.nodeType === 3) {
      renderTextNode(root);
      return;
    }
    if (root.nodeType !== 1 && root.nodeType !== 9 && root.nodeType !== 11) return;
    if (root.nodeType === 1) renderElement(root);
    const documentRef = root.ownerDocument || root;
    if (!documentRef.createTreeWalker || typeof NodeFilter === 'undefined') return;

    const elements = documentRef.createTreeWalker(root, NodeFilter.SHOW_ELEMENT);
    let element = elements.nextNode();
    while (element) {
      renderElement(element);
      element = elements.nextNode();
    }
    const textNodes = documentRef.createTreeWalker(root, NodeFilter.SHOW_TEXT);
    let textNode = textNodes.nextNode();
    while (textNode) {
      renderTextNode(textNode);
      textNode = textNodes.nextNode();
    }
  }

  function pauseObserver(callback) {
    const active = Boolean(observer && observerRoot);
    if (active) observer.disconnect();
    try { callback(); }
    finally {
      if (active) observer.observe(observerRoot, observerOptions);
    }
  }

  function translate(root = global.document) {
    if (!root) return;
    pauseObserver(() => walk(root));
    syncControls();
  }

  function queueTranslation(root) {
    if (root) pendingRoots.add(root);
    if (translationQueued) return;
    translationQueued = true;
    const schedule = global.queueMicrotask || (callback => Promise.resolve().then(callback));
    schedule(() => {
      translationQueued = false;
      const roots = Array.from(pendingRoots);
      pendingRoots.clear();
      pauseObserver(() => roots.forEach(walk));
      syncControls();
    });
  }

  function startObserver() {
    if (!global.document?.body || typeof global.MutationObserver !== 'function') return;
    observer?.disconnect();
    observerRoot = global.document.body;
    observerOptions = {
      childList: true,
      subtree: true,
      characterData: true,
      attributes: true,
      attributeFilter: [...translatableAttributes],
    };
    observer = new global.MutationObserver(mutations => {
      mutations.forEach(mutation => {
        if (mutation.type === 'characterData') {
          queueTranslation(mutation.target);
        } else if (mutation.type === 'attributes') {
          queueTranslation(mutation.target);
        } else {
          mutation.addedNodes.forEach(node => queueTranslation(node));
        }
      });
    });
    observer.observe(observerRoot, observerOptions);
  }

  function updateUrl(targetLocale, mode) {
    try {
      const url = new URL(global.location.href);
      if (mode === 'auto') url.searchParams.delete('lang');
      else url.searchParams.set('lang', languageCode(targetLocale));
      global.history?.replaceState?.(global.history.state, '', `${url.pathname}${url.search}${url.hash}`);
    } catch { /* Non-browser test environment. */ }
  }

  function applyDocumentLocale() {
    const root = global.document?.documentElement;
    if (!root) return;
    root.lang = locale;
    root.dataset.locale = locale;
    root.dataset.localeSource = localeSource;
  }

  function syncControls() {
    const controls = global.document?.querySelectorAll?.('[data-language-control]') || [];
    controls.forEach(control => {
      const requested = control.value;
      const isAuto = !normalizeLocale(safeStorageGet()) && !urlLocale();
      const desired = isAuto ? 'auto' : locale;
      if (requested !== desired) control.value = desired;
      const ariaLabel = t('language.control.label');
      const title = t('language.control.title');
      if (control.getAttribute('aria-label') !== ariaLabel) {
        control.setAttribute('aria-label', ariaLabel);
      }
      if (control.getAttribute('title') !== title) {
        control.setAttribute('title', title);
      }
    });
    const options = global.document?.querySelectorAll?.('[data-language-option]') || [];
    options.forEach(option => {
      const selected = normalizeLocale(option.dataset.languageOption) === locale;
      option.classList?.toggle?.('active', selected);
      const pressed = selected ? 'true' : 'false';
      if (option.getAttribute('aria-pressed') !== pressed) {
        option.setAttribute('aria-pressed', pressed);
      }
    });
    const switcher = global.document?.querySelector?.('[data-language-switch]');
    if (switcher) {
      const ariaLabel = t('language.control.label');
      const title = t('language.control.title');
      if (switcher.getAttribute('aria-label') !== ariaLabel) {
        switcher.setAttribute('aria-label', ariaLabel);
      }
      if (switcher.getAttribute('title') !== title) {
        switcher.setAttribute('title', title);
      }
    }
    const current = global.document?.querySelector?.('[data-current-language]');
    const currentLabel = languageCode(locale).toUpperCase();
    if (current && current.textContent !== currentLabel) current.textContent = currentLabel;
  }

  function dispatchChange(previousLocale, source) {
    if (!global.document?.dispatchEvent || typeof global.CustomEvent !== 'function') return;
    global.document.dispatchEvent(new global.CustomEvent('mssystem:languagechange', {
      detail: {
        locale,
        language: languageCode(locale),
        previousLocale,
        source,
      },
    }));
  }

  function setLocale(requestedLocale, options = {}) {
    const mode = String(requestedLocale || '').toLowerCase() === 'auto' ? 'auto' : 'manual';
    const normalized = mode === 'auto' ? browserLocale() : normalizeLocale(requestedLocale);
    if (!normalized) return locale;
    const previousLocale = locale;
    locale = normalized;
    localeSource = options.source || (mode === 'auto' ? 'browser' : 'user');

    if (options.persist !== false) {
      try {
        if (mode === 'auto') global.localStorage?.removeItem(STORAGE_KEY);
        else global.localStorage?.setItem(STORAGE_KEY, locale);
      } catch { /* Storage can be disabled by browser policy. */ }
    }
    if (options.updateUrl !== false) updateUrl(locale, mode);
    applyDocumentLocale();
    translate(global.document);
    dispatchChange(previousLocale, localeSource);
    return locale;
  }

  function toggle() {
    return setLocale(locale === 'en-US' ? 'zh-CN' : 'en-US');
  }

  function formatNumber(value, options) {
    return new Intl.NumberFormat(locale, options).format(value);
  }

  function formatDateTime(value, options = {}) {
    const date = value instanceof Date ? value : new Date(value);
    if (Number.isNaN(date.getTime())) return String(value ?? '');
    return new Intl.DateTimeFormat(locale, options).format(date);
  }

  function formatList(values, options) {
    if (typeof Intl.ListFormat !== 'function') return (values || []).join(', ');
    return new Intl.ListFormat(locale, options).format(values || []);
  }

  function registerCatalog(bundle = {}) {
    Object.entries(bundle.messages || {}).forEach(([targetLocale, entries]) => {
      catalog.messages[targetLocale] = Object.assign(catalog.messages[targetLocale] || {}, entries);
    });
    Object.entries(bundle.legacy || {}).forEach(([targetLocale, entries]) => {
      catalog.legacy[targetLocale] = Object.assign(catalog.legacy[targetLocale] || {}, entries);
    });
    Object.entries(bundle.patterns || {}).forEach(([targetLocale, entries]) => {
      catalog.patterns[targetLocale] = [
        ...(catalog.patterns[targetLocale] || []),
        ...(Array.isArray(entries) ? entries : []),
      ];
    });
    translate(global.document);
  }

  function audit(root = global.document) {
    const residual = new Set();
    if (!root?.createTreeWalker || typeof NodeFilter === 'undefined') return [];
    const scope = root.body || root;
    const nodes = root.createTreeWalker(scope, NodeFilter.SHOW_TEXT);
    let node = nodes.nextNode();
    while (node) {
      const value = normalizedText(node.nodeValue);
      if (value && HAN_RE.test(value) && !shouldSkip(node.parentElement)) residual.add(value);
      node = nodes.nextNode();
    }
    const attributeSelector = translatableAttributes.map(name => `[${name}]`).join(',');
    scope.querySelectorAll?.(attributeSelector).forEach(element => {
      if (element.closest?.(skipElementSelector)) return;
      translatableAttributes.forEach(name => {
        const value = normalizedText(element.getAttribute(name));
        if (value && HAN_RE.test(value)) residual.add(value);
      });
    });
    return Array.from(residual).sort();
  }

  function initControl(control) {
    if (!control || control.dataset.languageBound === 'true') return;
    control.dataset.languageBound = 'true';
    control.addEventListener('change', event => setLocale(event.target.value));
  }

  function initLanguageOption(option) {
    if (!option || option.dataset.languageBound === 'true') return;
    option.dataset.languageBound = 'true';
    option.addEventListener('click', () => setLocale(option.dataset.languageOption));
  }

  function initializeDom() {
    global.document?.querySelectorAll?.('[data-language-control]').forEach(initControl);
    global.document?.querySelectorAll?.('[data-language-option]').forEach(initLanguageOption);
    translate(global.document);
    startObserver();
    syncControls();
    if (new URL(global.location.href).searchParams.get('i18nDebug') === '1') {
      global.setTimeout?.(() => {
        const residual = audit();
        if (residual.length) console.warn('[MSsystem i18n] Untranslated UI strings:', residual);
      }, 1200);
    }
  }

  const initial = resolveLocale();
  locale = initial.locale;
  localeSource = initial.source;
  applyDocumentLocale();

  const api = Object.freeze({
    STORAGE_KEY,
    supportedLocales: SUPPORTED_LOCALES,
    normalizeLocale,
    resolveLocale,
    getLocale: () => locale,
    getLanguage: () => languageCode(locale),
    getSource: () => localeSource,
    setLocale,
    toggle,
    t,
    translateText,
    translate,
    registerCatalog,
    formatNumber,
    formatDateTime,
    formatList,
    audit,
    getMissingTranslations: () => Array.from(missing).sort(),
  });

  global.MSSLanguage = api;
  global.I18n = api;
  if (global.document?.readyState === 'loading') {
    global.document.addEventListener('DOMContentLoaded', initializeDom, { once: true });
  } else {
    initializeDom();
  }
})(typeof window !== 'undefined' ? window : globalThis);
