/**
 * Shared confidence/score interpretation renderer.
 *
 * All substantive wording comes from the server's versioned policy payload so
 * the result page and integrated report cannot silently drift from exports.
 */
const ConfidenceGuide = (() => {
  const cache = new Map();
  const esc = value => {
    if (window.Tasks?.escapeHtml) return Tasks.escapeHtml(String(value == null ? '' : value));
    const node = document.createElement('div');
    node.textContent = String(value == null ? '' : value);
    return node.innerHTML;
  };

  function fallbackPayload(language = 'zh') {
    if (String(language).toLowerCase().startsWith('en')) {
      return {
        language: 'en', schema_version: 'fallback', policy_id: 'fallback-boundary-only',
        title: 'How to interpret confidence',
        subtitle: 'The complete guide could not be loaded; the safety boundary remains available.',
        read_first: 'High is not confirmation, and confidence or any displayed score is not the probability that the compound is present.',
        threshold_note: 'Use definitive wording only when the formal conclusion is explicitly released/confirmed.',
        levels: [], identity_scopes: [], scores: [],
      };
    }
    return {
      language: 'zh', schema_version: 'fallback', policy_id: 'fallback-boundary-only',
      title: '怎样理解置信度',
      subtitle: '完整说明暂时未加载，但以下结论边界仍然有效。',
      read_first: 'High 不等于已经确证该化合物；置信度和各项分数也不是“该化合物存在的概率”。',
      threshold_note: '只有正式结论状态显示“已确证／已发布”时，才可使用确定性表述。',
      levels: [], identity_scopes: [], scores: [],
    };
  }

  async function load(language = 'zh', refresh = false) {
    const normalized = String(language).toLowerCase().startsWith('en') ? 'en' : 'zh';
    if (refresh) cache.delete(normalized);
    if (!cache.has(normalized)) {
      cache.set(normalized, fetch('/static/policies/confidence-guide.json?v=1.0.0', {
        cache: 'no-store',
      }).then(response => {
        if (!response.ok) throw new Error(`Static confidence guide HTTP ${response.status}`);
        return response.json();
      }).then(payloads => {
        const payload = payloads?.[normalized];
        if (!payload || !Array.isArray(payload.levels)) {
          throw new Error('Static confidence guide payload is invalid');
        }
        return payload;
      }).catch(staticError => API.get(
        `/interpretation/confidence-guide?language=${normalized}`
      ).catch(apiError => {
        cache.delete(normalized);
        console.warn('Confidence guide unavailable:', staticError, apiError);
        return fallbackPayload(normalized);
      })));
    }
    return cache.get(normalized);
  }

  const token = value => String(value || 'unknown')
    .toLowerCase().replace(/[^a-z0-9_-]/g, '') || 'unknown';

  function shortLevelLabel(item, zh) {
    const labels = zh
      ? { confirmed: '已确证', high: 'High', medium: 'Medium', low: 'Low', unknown: 'Unknown', background: '背景相关' }
      : { confirmed: 'Confirmed', high: 'High', medium: 'Medium', low: 'Low', unknown: 'Unknown', background: 'Background' };
    return labels[token(item.key || item.tone)] || item.label;
  }

  function styles(id) {
    return `<style>
      #${id}{--guide-navy:#102a43;--guide-teal:#16746f;--guide-line:#d9e3eb;--guide-muted:#64748b;--guide-soft:#f5f8fb;margin:1.05rem 0;color:#172b3a;font:14px/1.65 Inter,-apple-system,BlinkMacSystemFont,"Segoe UI","PingFang SC","Microsoft YaHei",sans-serif;box-sizing:border-box}
      #${id},#${id} *{box-sizing:border-box} #${id} p,#${id} h3,#${id} h4,#${id} ol{margin:0} #${id} summary{list-style:none} #${id} summary::-webkit-details-marker{display:none}
      #${id} .confidence-guide-shell{overflow:hidden;border:1px solid #d5e0e8;border-radius:18px;background:#fff;box-shadow:0 14px 38px rgba(22,41,60,.09)}
      #${id} .confidence-guide-hero{position:relative;display:grid;grid-template-columns:54px minmax(0,1fr) auto;gap:1rem;align-items:center;min-height:146px;padding:1.35rem 1.45rem;overflow:hidden;background:linear-gradient(124deg,#10273f 0%,#123c50 56%,#176e70 100%);color:#fff;cursor:pointer;isolation:isolate}
      #${id} .confidence-guide-hero:before,#${id} .confidence-guide-hero:after{position:absolute;border:1px solid rgba(255,255,255,.1);border-radius:50%;content:"";z-index:-1} #${id} .confidence-guide-hero:before{width:260px;height:260px;right:9%;top:-194px} #${id} .confidence-guide-hero:after{width:170px;height:170px;right:-58px;bottom:-112px;background:rgba(255,255,255,.035)}
      #${id} .confidence-guide-mark{display:grid;width:52px;height:52px;place-items:center;border:1px solid rgba(255,255,255,.28);border-radius:16px;background:rgba(255,255,255,.11);box-shadow:inset 0 1px 0 rgba(255,255,255,.16);font:800 26px/1 Georgia,serif}
      #${id} .confidence-guide-eyebrow{display:block;margin-bottom:.18rem;color:#9fe4de;font-size:.61rem;font-weight:800;letter-spacing:.16em;text-transform:uppercase} #${id} .confidence-guide-heading h3{color:#fff;font-size:clamp(1.12rem,2vw,1.45rem);font-weight:760;letter-spacing:-.02em;line-height:1.3} #${id} .confidence-guide-subtitle{margin-top:.2rem!important;color:rgba(255,255,255,.72);font-size:.73rem}
      #${id} .confidence-guide-boundaries{display:flex;flex-wrap:wrap;gap:.42rem;margin-top:.72rem} #${id} .confidence-guide-boundaries span{display:inline-flex;align-items:center;gap:.35rem;padding:.25rem .52rem;border:1px solid rgba(255,255,255,.17);border-radius:999px;background:rgba(255,255,255,.09);color:#e7f6f4;font-size:.65rem;font-weight:650} #${id} .confidence-guide-boundaries b{color:#fff;font-weight:800}
      #${id} .confidence-guide-toggle{display:inline-flex;align-items:center;gap:.5rem;justify-self:end;padding:.48rem .72rem;border:1px solid rgba(255,255,255,.22);border-radius:999px;background:rgba(255,255,255,.1);color:#fff;font-size:.66rem;font-weight:750;white-space:nowrap;transition:background .18s ease,transform .18s ease} #${id} .confidence-guide-toggle:after{width:7px;height:7px;border-right:2px solid currentColor;border-bottom:2px solid currentColor;content:"";transform:rotate(45deg) translateY(-2px);transition:transform .18s ease} #${id} .confidence-guide-shell[open] .confidence-guide-toggle{background:rgba(255,255,255,.17)} #${id} .confidence-guide-shell[open] .confidence-guide-toggle:after{transform:rotate(225deg) translate(-1px,-1px)} #${id} .confidence-guide-shell[open] .confidence-guide-toggle .when-closed,#${id} .confidence-guide-shell:not([open]) .confidence-guide-toggle .when-open{display:none}
      #${id} .confidence-guide-body{padding:1.25rem 1.35rem 1rem;background:linear-gradient(180deg,#f8fafc 0%,#fff 48%)}
      #${id} .confidence-guide-core{display:grid;grid-template-columns:minmax(0,1.25fr) minmax(270px,.75fr);gap:.8rem;margin-bottom:1rem} #${id} .confidence-guide-principle{display:grid;grid-template-columns:auto minmax(0,1fr);gap:.75rem;align-items:start;padding:.9rem 1rem;border:1px solid #eed8a9;border-radius:12px;background:#fffaf0;color:#60491d} #${id} .confidence-guide-principle>span{display:grid;width:28px;height:28px;place-items:center;border-radius:9px;background:#f2c45e;color:#563c08;font-weight:850} #${id} .confidence-guide-principle b{display:block;margin-bottom:.15rem;color:#654914;font-size:.72rem} #${id} .confidence-guide-principle p{font-size:.72rem;line-height:1.65}
      #${id} .confidence-guide-rule-note{padding:.9rem 1rem;border:1px solid #d6e3eb;border-radius:12px;background:#fff;color:#536777;font-size:.69rem;line-height:1.62} #${id} .confidence-guide-rule-note b{display:block;margin-bottom:.16rem;color:#28485e;font-size:.72rem}
      #${id} .confidence-guide-route{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:.55rem;margin-bottom:1.15rem} #${id} .confidence-guide-route>div{position:relative;min-width:0;padding:.72rem .78rem .7rem 2.65rem;border:1px solid #dbe5ec;border-radius:11px;background:#fff} #${id} .confidence-guide-route>div:not(:last-child):after{position:absolute;right:-.38rem;top:50%;display:grid;width:.75rem;height:.75rem;place-items:center;border-radius:50%;background:#d9ebe8;color:#26736d;font-size:.68rem;line-height:1;content:"›";transform:translateY(-50%);z-index:1} #${id} .confidence-guide-route span{position:absolute;left:.72rem;top:.7rem;color:#2a7974;font-size:.61rem;font-weight:850;letter-spacing:.08em} #${id} .confidence-guide-route b{display:block;color:#243d50;font-size:.71rem} #${id} .confidence-guide-route small{display:block;margin-top:.14rem;color:#718196;font-size:.62rem;line-height:1.45}
      #${id} .confidence-guide-section-head{display:flex;align-items:flex-end;justify-content:space-between;gap:1rem;margin:0 0 .6rem} #${id} .confidence-guide-section-kicker{display:block;color:#258078;font-size:.6rem;font-weight:850;letter-spacing:.13em;text-transform:uppercase} #${id} .confidence-guide-section-head h4{margin-top:.06rem;color:#1d3649;font-size:.9rem} #${id} .confidence-guide-section-head p{color:#718196;font-size:.64rem;text-align:right}
      #${id} .confidence-guide-level-switcher{position:relative;margin-bottom:.85rem} #${id} .confidence-guide-level-radio{position:absolute;width:1px;height:1px;overflow:hidden;clip:rect(0 0 0 0);clip-path:inset(50%);white-space:nowrap} #${id} .confidence-guide-level-tablist{display:grid;grid-template-columns:repeat(6,minmax(0,1fr));gap:.35rem;padding:.34rem;border:1px solid #d9e3eb;border-radius:12px;background:#edf2f6} #${id} .confidence-guide-level-tab{display:flex;min-width:0;align-items:center;justify-content:center;gap:.35rem;padding:.48rem .3rem;border:1px solid transparent;border-radius:9px;color:#66788a;font-size:.64rem;font-weight:760;line-height:1.25;text-align:center;cursor:pointer;transition:all .16s ease} #${id} .confidence-guide-level-tab:before{width:7px;height:7px;flex:0 0 auto;border-radius:50%;background:var(--level-color,#7b8794);content:""} #${id} .confidence-guide-level-tab:hover{background:rgba(255,255,255,.72);color:#294255}
      #${id} .confidence-guide-level-panels{margin-top:.55rem} #${id} .confidence-guide-level-panel{display:none;position:relative;overflow:hidden;padding:1rem 1.05rem;border:1px solid color-mix(in srgb,var(--level-color) 28%,#dce5eb);border-radius:13px;background:linear-gradient(135deg,var(--level-soft) 0%,#fff 58%)} #${id} .confidence-guide-level-panel:before{position:absolute;inset:0 auto 0 0;width:4px;background:var(--level-color);content:""} #${id} .confidence-guide-level-panel-head{display:flex;align-items:center;gap:.65rem;margin-bottom:.55rem} #${id} .confidence-guide-level-index{display:grid;width:30px;height:30px;flex:0 0 auto;place-items:center;border-radius:9px;background:var(--level-color);color:#fff;font-size:.62rem;font-weight:850} #${id} .confidence-guide-level-panel h4{color:var(--level-ink);font-size:.92rem;line-height:1.35} #${id} .confidence-guide-level-meaning{max-width:880px;color:#405768;font-size:.72rem;line-height:1.65}
      #${id} .confidence-guide-actions{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:.55rem;margin-top:.72rem} #${id} .confidence-guide-action{padding:.66rem .72rem;border:1px solid #dce5ea;border-radius:9px;background:rgba(255,255,255,.82)} #${id} .confidence-guide-action b{display:block;margin-bottom:.16rem;color:#2e4a5e;font-size:.64rem} #${id} .confidence-guide-action p{color:#627483;font-size:.64rem;line-height:1.55} #${id} .confidence-guide-action-use{border-color:#cce5d8;background:#f4fbf7} #${id} .confidence-guide-action-use b{color:#207052} #${id} .confidence-guide-action-avoid{border-color:#efd5cb;background:#fff8f5} #${id} .confidence-guide-action-avoid b{color:#a44c2b} #${id} .confidence-guide-action-next{border-color:#cfdeee;background:#f6f9fd} #${id} .confidence-guide-action-next b{color:#35699a}
      #${id} .guide-tone-confirmed{--level-color:#147d64;--level-soft:#edf9f4;--level-ink:#125b4a} #${id} .guide-tone-high{--level-color:#2474b8;--level-soft:#eef6fd;--level-ink:#1f5e91} #${id} .guide-tone-medium{--level-color:#c28a20;--level-soft:#fff9eb;--level-ink:#865f14} #${id} .guide-tone-low{--level-color:#bc5a31;--level-soft:#fff5f0;--level-ink:#8c4022} #${id} .guide-tone-unknown{--level-color:#748191;--level-soft:#f3f5f7;--level-ink:#53606e} #${id} .guide-tone-background{--level-color:#7456a6;--level-soft:#f6f2fb;--level-ink:#5e438d}
      #${id} .confidence-guide-level-radio[data-tone="confirmed"]:checked~.confidence-guide-level-tablist [data-tone="confirmed"],#${id} .confidence-guide-level-radio[data-tone="high"]:checked~.confidence-guide-level-tablist [data-tone="high"],#${id} .confidence-guide-level-radio[data-tone="medium"]:checked~.confidence-guide-level-tablist [data-tone="medium"],#${id} .confidence-guide-level-radio[data-tone="low"]:checked~.confidence-guide-level-tablist [data-tone="low"],#${id} .confidence-guide-level-radio[data-tone="unknown"]:checked~.confidence-guide-level-tablist [data-tone="unknown"],#${id} .confidence-guide-level-radio[data-tone="background"]:checked~.confidence-guide-level-tablist [data-tone="background"]{border-color:color-mix(in srgb,var(--level-color) 34%,#fff);background:#fff;box-shadow:0 3px 10px rgba(30,51,70,.08);color:var(--level-ink)}
      #${id} .confidence-guide-level-radio[data-tone="confirmed"]:checked~.confidence-guide-level-panels [data-tone="confirmed"],#${id} .confidence-guide-level-radio[data-tone="high"]:checked~.confidence-guide-level-panels [data-tone="high"],#${id} .confidence-guide-level-radio[data-tone="medium"]:checked~.confidence-guide-level-panels [data-tone="medium"],#${id} .confidence-guide-level-radio[data-tone="low"]:checked~.confidence-guide-level-panels [data-tone="low"],#${id} .confidence-guide-level-radio[data-tone="unknown"]:checked~.confidence-guide-level-panels [data-tone="unknown"],#${id} .confidence-guide-level-radio[data-tone="background"]:checked~.confidence-guide-level-panels [data-tone="background"]{display:block}
      #${id} .confidence-guide-disclosure{margin-top:.55rem;border:1px solid #d8e3ea;border-radius:12px;background:#fff} #${id} .confidence-guide-disclosure>summary{display:flex;align-items:center;justify-content:space-between;gap:1rem;padding:.78rem .9rem;cursor:pointer} #${id} .confidence-guide-disclosure>summary>span:first-child{display:flex;min-width:0;align-items:center;gap:.65rem} #${id} .confidence-guide-disclosure-index{display:grid;width:29px;height:29px;flex:0 0 auto;place-items:center;border-radius:9px;background:#e7f3f1;color:#28736d;font-size:.61rem;font-weight:850} #${id} .confidence-guide-disclosure-title b{display:block;color:#294457;font-size:.72rem} #${id} .confidence-guide-disclosure-title small{display:block;margin-top:.05rem;color:#7a8999;font-size:.61rem} #${id} .confidence-guide-disclosure-meta{display:flex;align-items:center;gap:.55rem;color:#6d7e8e;font-size:.61rem;font-weight:700;white-space:nowrap} #${id} .confidence-guide-disclosure-meta:after{width:7px;height:7px;border-right:2px solid #728495;border-bottom:2px solid #728495;content:"";transform:rotate(45deg) translateY(-2px);transition:transform .16s ease} #${id} .confidence-guide-disclosure[open]>.confidence-guide-disclosure-meta,#${id} .confidence-guide-disclosure[open]>summary .confidence-guide-disclosure-meta:after{transform:rotate(225deg) translate(-1px,-1px)} #${id} .confidence-guide-disclosure[open]>summary{border-bottom:1px solid #e3eaef;background:#f8fbfc}
      #${id} .confidence-guide-scope-list{display:grid;gap:0;padding:.25rem .9rem .7rem;list-style:none;counter-reset:guide-scope} #${id} .confidence-guide-scope-list li{position:relative;display:grid;grid-template-columns:32px minmax(0,1fr);gap:.65rem;padding:.7rem 0;counter-increment:guide-scope} #${id} .confidence-guide-scope-list li+li{border-top:1px solid #e8edf1} #${id} .confidence-guide-scope-list li:before{display:grid;width:28px;height:28px;place-items:center;border:1px solid #cddce5;border-radius:50%;background:#f6fafb;color:#557487;font-size:.58rem;font-weight:800;content:counter(guide-scope,decimal-leading-zero)} #${id} .confidence-guide-scope-list b{display:block;color:#29485a;font-size:.7rem} #${id} .confidence-guide-scope-list p{margin-top:.1rem;color:#617483;font-size:.64rem;line-height:1.55} #${id} .confidence-guide-scope-list em{display:block;margin-top:.22rem;color:#2b746d;font-size:.62rem;font-style:normal;font-weight:650}
      #${id} .confidence-guide-score-grid{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:.5rem;padding:.7rem} #${id} .confidence-guide-score{margin:0;border:1px solid #dfe7ec;border-radius:9px;background:#fbfcfd} #${id} .confidence-guide-score>summary{display:grid;grid-template-columns:26px minmax(0,1fr) auto;gap:.55rem;align-items:center;padding:.62rem .68rem;cursor:pointer} #${id} .confidence-guide-score-number{display:grid;width:25px;height:25px;place-items:center;border-radius:7px;background:#edf3f7;color:#567083;font:750 .56rem/1 ui-monospace,SFMono-Regular,Menlo,monospace} #${id} .confidence-guide-score>summary b{min-width:0;color:#324c60;font-size:.65rem;line-height:1.35} #${id} .confidence-guide-score>summary:after{color:#6f8495;font-size:.9rem;content:"+"} #${id} .confidence-guide-score[open]>summary:after{content:"−"} #${id} .confidence-guide-score[open]{border-color:#bfd3df;background:#fff} #${id} .confidence-guide-score-body{padding:.02rem .72rem .68rem 3.92rem} #${id} .confidence-guide-score-body p{color:#607484;font-size:.63rem;line-height:1.55} #${id} .confidence-guide-score-boundary{margin-top:.4rem!important;padding:.42rem .5rem;border-radius:7px;background:#fff5f1;color:#91462c!important} #${id} .confidence-guide-score-boundary b{font-size:.6rem}
      #${id} .confidence-guide-policy{display:flex;justify-content:flex-end;gap:.35rem;margin-top:.75rem;color:#96a2af;font-size:.54rem;letter-spacing:.01em} #${id} .confidence-guide-policy span{font-family:ui-monospace,SFMono-Regular,Menlo,monospace}
      #${id}.confidence-guide-report{margin:.9rem 0 1.2rem} #${id} .confidence-guide-report-shell{overflow:hidden;border:1px solid #d5e0e8;border-radius:15px;background:#fff;box-shadow:0 9px 28px rgba(22,41,60,.07)}
      #${id} .confidence-guide-report-head{display:grid;grid-template-columns:minmax(0,1fr) auto;gap:1.15rem;align-items:center;padding:1rem 1.15rem;border-left:6px solid #16746f;background:linear-gradient(108deg,#f3f9f8 0%,#fff 72%)} #${id} .confidence-guide-report-kicker{display:block;color:#178078;font-size:.58rem;font-weight:850;letter-spacing:.16em} #${id} .confidence-guide-report-head h3{margin:.08rem 0 .12rem;color:#17364b;font-size:1.12rem;line-height:1.3} #${id} .confidence-guide-report-head p{color:#697c89;font-size:.68rem} #${id} .confidence-guide-report-boundaries{display:flex;flex-wrap:wrap;justify-content:flex-end;gap:.38rem} #${id} .confidence-guide-report-boundaries span{padding:.3rem .52rem;border:1px solid #cedde4;border-radius:999px;background:#fff;color:#526875;font-size:.61rem;white-space:nowrap} #${id} .confidence-guide-report-boundaries b{color:#1d475e}
      #${id} .confidence-guide-report-body{padding:.9rem 1.05rem .8rem;background:#fbfcfd} #${id} .confidence-guide-report-body .confidence-guide-core{margin-bottom:.7rem} #${id} .confidence-guide-report-body .confidence-guide-route{margin-bottom:.8rem}
      #${id} .confidence-guide-report-tier-map{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:.45rem;margin-bottom:.55rem} #${id} .confidence-guide-report-tier{display:grid;grid-template-columns:25px minmax(0,1fr);gap:.52rem;min-width:0;padding:.58rem .62rem;border:1px solid #dce5ea;border-radius:9px;background:#fff} #${id} .confidence-guide-report-tier>span{display:grid;width:24px;height:24px;place-items:center;border-radius:7px;background:var(--level-color);color:#fff;font-size:.52rem;font-weight:850} #${id} .confidence-guide-report-tier b{display:block;color:var(--level-ink);font-size:.65rem;line-height:1.35} #${id} .confidence-guide-report-tier p{margin-top:.08rem;color:#6d7e89;font-size:.59rem;line-height:1.45} #${id} .confidence-guide-report-full .confidence-guide-level-switcher{padding:.72rem;margin:0} #${id} .confidence-guide-report-full .confidence-guide-level-panels{margin-bottom:0}
      @media(max-width:820px){#${id} .confidence-guide-hero{grid-template-columns:46px minmax(0,1fr);padding:1.15rem} #${id} .confidence-guide-mark{width:44px;height:44px;border-radius:13px;font-size:22px} #${id} .confidence-guide-toggle{grid-column:2;justify-self:start} #${id} .confidence-guide-core{grid-template-columns:1fr} #${id} .confidence-guide-level-tablist{grid-template-columns:repeat(3,minmax(0,1fr))} #${id} .confidence-guide-actions{grid-template-columns:1fr} #${id} .confidence-guide-report-head{grid-template-columns:1fr} #${id} .confidence-guide-report-boundaries{justify-content:flex-start} #${id} .confidence-guide-report-tier-map{grid-template-columns:repeat(2,minmax(0,1fr))}}
      @media(max-width:580px){#${id}{margin:.75rem 0} #${id} .confidence-guide-shell,#${id} .confidence-guide-report-shell{border-radius:14px} #${id} .confidence-guide-hero{grid-template-columns:1fr;gap:.7rem} #${id} .confidence-guide-mark{display:none} #${id} .confidence-guide-toggle{grid-column:1} #${id} .confidence-guide-body,#${id} .confidence-guide-report-body{padding:1rem .8rem .8rem} #${id} .confidence-guide-route{grid-template-columns:1fr} #${id} .confidence-guide-route>div:not(:last-child):after{display:none} #${id} .confidence-guide-score-grid,#${id} .confidence-guide-report-tier-map{grid-template-columns:1fr} #${id} .confidence-guide-section-head{align-items:flex-start;flex-direction:column;gap:.2rem} #${id} .confidence-guide-section-head p{text-align:left} #${id} .confidence-guide-disclosure-meta{font-size:0} #${id} .confidence-guide-score-body{padding-left:.72rem}}
      @media print{#${id}{margin:0 0 12px;break-inside:auto} #${id} .confidence-guide-shell,#${id} .confidence-guide-report-shell{box-shadow:none} #${id} .confidence-guide-hero{min-height:0;padding:16px 18px;-webkit-print-color-adjust:exact;print-color-adjust:exact} #${id} .confidence-guide-toggle{display:none} #${id} details:not([open])>*:not(summary){display:block!important} #${id} .confidence-guide-report-tier-map{display:none} #${id} .confidence-guide-level-tablist{display:none} #${id} .confidence-guide-level-panels{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:7px} #${id} .confidence-guide-level-panel{display:block!important;break-inside:avoid} #${id} .confidence-guide-actions{grid-template-columns:1fr} #${id} .confidence-guide-disclosure>summary{break-after:avoid} #${id} .confidence-guide-score-grid{grid-template-columns:repeat(2,minmax(0,1fr))} #${id} .confidence-guide-score-body{display:block!important} #${id} .confidence-guide-policy{display:none}}
    </style>`;
  }

  function render(payload, options = {}) {
    const guide = payload || fallbackPayload(options.language || 'zh');
    const zh = guide.language !== 'en';
    const id = String(options.id || 'mssystem-confidence-guide').replace(/[^A-Za-z0-9_-]/g, '-');
    const reportMode = options.mode === 'report';
    const levelItems = guide.levels || [];
    const levelRadios = levelItems.map((item, index) => {
      const tone = token(item.tone || item.key);
      const checked = tone === 'high' || (!levelItems.some(level => token(level.tone || level.key) === 'high') && index === 0);
      return `<input class="confidence-guide-level-radio" type="radio" name="${id}-evidence-level" id="${id}-level-${index}" data-tone="${esc(tone)}" ${checked ? 'checked' : ''}>`;
    }).join('');
    const levelTabs = levelItems.map((item, index) => {
      const tone = token(item.tone || item.key);
      return `<label class="confidence-guide-level-tab guide-tone-${esc(tone)}" for="${id}-level-${index}" data-tone="${esc(tone)}">${esc(shortLevelLabel(item, zh))}</label>`;
    }).join('');
    const levelPanels = levelItems.map((item, index) => {
      const tone = token(item.tone || item.key);
      return `<article class="confidence-guide-level-panel guide-tone-${esc(tone)}" data-tone="${esc(tone)}">
        <div class="confidence-guide-level-panel-head"><span class="confidence-guide-level-index">${String(index + 1).padStart(2, '0')}</span><h4>${esc(item.label)}</h4></div>
        <p class="confidence-guide-level-meaning">${esc(item.meaning)}</p>
        <div class="confidence-guide-actions">
          <div class="confidence-guide-action confidence-guide-action-use"><b>${zh ? '✓ 可以这样使用' : '✓ Appropriate use'}</b><p>${esc(item.use)}</p></div>
          <div class="confidence-guide-action confidence-guide-action-avoid"><b>${zh ? '× 不要这样理解' : '× Do not infer'}</b><p>${esc(item.avoid)}</p></div>
          <div class="confidence-guide-action confidence-guide-action-next"><b>${zh ? '→ 建议下一步' : '→ Recommended next step'}</b><p>${esc(item.next_step)}</p></div>
        </div>
      </article>`;
    }).join('');
    const scopes = (guide.identity_scopes || []).map(item => `<li><div><b>${esc(item.label)}</b><p>${esc(item.meaning)}</p><em>${zh ? '建议表述：' : 'Recommended wording: '}${esc(item.report)}</em></div></li>`).join('');
    const scores = (guide.scores || []).map((item, index) => `<details class="confidence-guide-score">
      <summary><span class="confidence-guide-score-number">${String(index + 1).padStart(2, '0')}</span><b>${esc(item.label)}</b></summary>
      <div class="confidence-guide-score-body"><p>${esc(item.meaning)}</p><p class="confidence-guide-score-boundary"><b>${zh ? '注意：' : 'Boundary: '}</b>${esc(item.boundary)}</p></div>
    </details>`).join('');
    const route = zh
      ? [
        ['01', '先看证据等级', 'High／Medium／Low 回答证据强弱'],
        ['02', '再看身份范围', '完整结构、二维连接与异构体组不可混用'],
        ['03', '最后看正式状态', '只有已确证／已发布才可确定性表述'],
      ]
      : [
        ['01', 'Read the evidence tier', 'High / Medium / Low describes evidence strength'],
        ['02', 'Check identity scope', 'Full identity, 2-D structure and isomer group differ'],
        ['03', 'Check formal status', 'Use definitive wording only after formal release'],
      ];
    const routeHtml = route.map(item => `<div><span>${item[0]}</span><b>${item[1]}</b><small>${item[2]}</small></div>`).join('');
    const reportTierMap = levelItems.map((item, index) => {
      const tone = token(item.tone || item.key);
      return `<article class="confidence-guide-report-tier guide-tone-${esc(tone)}"><span>${String(index + 1).padStart(2, '0')}</span><div><b>${esc(item.label)}</b><p>${esc(item.meaning)}</p></div></article>`;
    }).join('');
    if (reportMode) {
      return `${styles(id)}<section id="${id}" class="confidence-guide confidence-guide-report" data-policy-id="${esc(guide.policy_id || '')}" data-schema-version="${esc(guide.schema_version || '')}">
        <div class="confidence-guide-report-shell">
          <header class="confidence-guide-report-head"><div><span class="confidence-guide-report-kicker">REPORT READING GUIDE</span><h3>${esc(guide.title)}</h3><p>${esc(guide.subtitle)}</p></div><div class="confidence-guide-report-boundaries"><span><b>High</b> ≠ ${zh ? '已确证' : 'confirmed'}</span><span><b>${zh ? '分数' : 'Score'}</b> ≠ ${zh ? '存在概率' : 'presence probability'}</span></div></header>
          <div class="confidence-guide-report-body">
            <div class="confidence-guide-core"><div class="confidence-guide-principle"><span>!</span><div><b>${zh ? '阅读报告前先记住' : 'Remember before reading'}</b><p>${esc(guide.read_first)}</p></div></div><div class="confidence-guide-rule-note"><b>${zh ? '为什么不能只看一个分数？' : 'Why is one score not enough?'}</b>${esc(guide.threshold_note)}</div></div>
            <div class="confidence-guide-route">${routeHtml}</div>
            ${reportTierMap ? `<div class="confidence-guide-section-head"><div><span class="confidence-guide-section-kicker">EVIDENCE TIER · QUICK MAP</span><h4>${zh ? '证据等级速查' : 'Evidence-tier quick map'}</h4></div><p>${zh ? '先定位等级，再按需展开完整定义' : 'Locate the tier first; open definitions when needed'}</p></div><div class="confidence-guide-report-tier-map">${reportTierMap}</div><details class="confidence-guide-disclosure confidence-guide-report-full"><summary><span><span class="confidence-guide-disclosure-index">01</span><span class="confidence-guide-disclosure-title"><b>${zh ? '完整证据等级说明' : 'Full evidence-tier definitions'}</b><small>${zh ? '使用边界、误读风险与建议下一步' : 'Use, misinterpretation boundary and next step'}</small></span></span><span class="confidence-guide-disclosure-meta">${levelItems.length} ${zh ? '级' : 'tiers'}</span></summary><div class="confidence-guide-level-switcher">${levelRadios}<div class="confidence-guide-level-tablist">${levelTabs}</div><div class="confidence-guide-level-panels">${levelPanels}</div></div></details>` : ''}
            ${scopes ? `<details class="confidence-guide-disclosure"><summary><span><span class="confidence-guide-disclosure-index">02</span><span class="confidence-guide-disclosure-title"><b>${zh ? '身份范围：最多能把结构说到哪一步' : 'Identity scope: how specifically the structure may be named'}</b><small>${zh ? '完整身份、二维连接、异构体组与类别是不同结论层级' : 'Full identity, 2-D connectivity, isomer group and class are different claim levels'}</small></span></span><span class="confidence-guide-disclosure-meta">${(guide.identity_scopes || []).length} ${zh ? '层' : 'levels'}</span></summary><ol class="confidence-guide-scope-list">${scopes}</ol></details>` : ''}
            ${scores ? `<details class="confidence-guide-disclosure"><summary><span><span class="confidence-guide-disclosure-index">03</span><span class="confidence-guide-disclosure-title"><b>${zh ? '各项评分与辅助字段' : 'Scores and auxiliary fields'}</b><small>${zh ? '按名称找到指标，点开查看含义与误读边界' : 'Find a metric by name, then open its meaning and boundary'}</small></span></span><span class="confidence-guide-disclosure-meta">${(guide.scores || []).length} ${zh ? '项' : 'fields'}</span></summary><div class="confidence-guide-score-grid">${scores}</div></details>` : ''}
            <div class="confidence-guide-policy">${zh ? '解释规则' : 'Interpretation policy'} <span>${esc(guide.policy_id || '—')} · ${esc(guide.schema_version || '—')}</span></div>
          </div>
        </div>
      </section>`;
    }
    return `${styles(id)}<section id="${id}" class="confidence-guide" data-policy-id="${esc(guide.policy_id || '')}" data-schema-version="${esc(guide.schema_version || '')}">
      <details class="confidence-guide-shell">
        <summary class="confidence-guide-hero">
          <span class="confidence-guide-mark" aria-hidden="true">i</span>
          <span class="confidence-guide-heading"><span class="confidence-guide-eyebrow">${zh ? 'RESULT READING GUIDE' : 'RESULT READING GUIDE'}</span><h3>${esc(guide.title)}</h3><span class="confidence-guide-subtitle">${esc(guide.subtitle)}</span><span class="confidence-guide-boundaries"><span><b>High</b> ≠ ${zh ? '已确证' : 'confirmed'}</span><span><b>${zh ? '分数' : 'Score'}</b> ≠ ${zh ? '存在概率' : 'presence probability'}</span></span></span>
          <span class="confidence-guide-toggle"><span class="when-closed">${zh ? '展开解读' : 'Open guide'}</span><span class="when-open">${zh ? '收起说明' : 'Close guide'}</span></span>
        </summary>
        <div class="confidence-guide-body">
          <div class="confidence-guide-core"><div class="confidence-guide-principle"><span>!</span><div><b>${zh ? '先记住这条边界' : 'The boundary to remember'}</b><p>${esc(guide.read_first)}</p></div></div><div class="confidence-guide-rule-note"><b>${zh ? '为什么不能只看一个分数？' : 'Why is one score not enough?'}</b>${esc(guide.threshold_note)}</div></div>
          <div class="confidence-guide-route">${routeHtml}</div>
          ${levelPanels ? `<div class="confidence-guide-section-head"><div><span class="confidence-guide-section-kicker">01 · EVIDENCE TIER</span><h4>${zh ? '证据等级怎样用' : 'How to use evidence tiers'}</h4></div><p>${zh ? '点击等级，查看使用边界和建议下一步' : 'Select a tier to see its use, boundary and next step'}</p></div><div class="confidence-guide-level-switcher">${levelRadios}<div class="confidence-guide-level-tablist">${levelTabs}</div><div class="confidence-guide-level-panels">${levelPanels}</div></div>` : ''}
          ${scopes ? `<details class="confidence-guide-disclosure"><summary><span><span class="confidence-guide-disclosure-index">02</span><span class="confidence-guide-disclosure-title"><b>${zh ? '身份范围：最多能把结构说到哪一步' : 'Identity scope: how specifically the structure may be named'}</b><small>${zh ? '完整身份、二维连接、异构体组与类别是不同结论层级' : 'Full identity, 2-D connectivity, isomer group and class are different claim levels'}</small></span></span><span class="confidence-guide-disclosure-meta">${(guide.identity_scopes || []).length} ${zh ? '层' : 'levels'}</span></summary><ol class="confidence-guide-scope-list">${scopes}</ol></details>` : ''}
          ${scores ? `<details class="confidence-guide-disclosure"><summary><span><span class="confidence-guide-disclosure-index">03</span><span class="confidence-guide-disclosure-title"><b>${zh ? '各项评分与辅助字段' : 'Scores and auxiliary fields'}</b><small>${zh ? '按名称找到指标，点开查看含义与误读边界' : 'Find a metric by name, then open its meaning and boundary'}</small></span></span><span class="confidence-guide-disclosure-meta">${(guide.scores || []).length} ${zh ? '项' : 'fields'}</span></summary><div class="confidence-guide-score-grid">${scores}</div></details>` : ''}
          <div class="confidence-guide-policy">${zh ? '解释规则' : 'Interpretation policy'} <span>${esc(guide.policy_id || '—')} · ${esc(guide.schema_version || '—')}</span></div>
        </div>
      </details>
    </section>`;
  }

  return { load, render, fallbackPayload };
})();

window.ConfidenceGuide = ConfidenceGuide;
