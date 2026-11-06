"""Read-only projection for candidates retained below the governance release gate.

The projection is deliberately separate from identification tables and counts.
It never mutates saved rows, candidate order, confidence, or release decisions.
One saved scan candidate produces one diagnostic row; no identity deduplication is
performed because scan-specific governance reasons must remain inspectable.
"""

from __future__ import annotations

import csv
import io
import math
import re
from html import escape
from typing import Any, Iterable, Mapping


BLOCKED_CLAIM_STATUSES = frozenset({
    "insufficient_for_release",
    "reported_but_governance_blocked",
})

# Shared presentation-only aliases. They never participate in ranking, identity
# grouping, evidence assessment, or release governance.
ZH_COMPOUND_ALIASES = {
    "sinomenine": ["青藤碱", "清风藤碱"],
    "magnoflorine": ["木兰碱", "广玉兰碱"],
    "berberine": ["小檗碱"], "palmatine": ["巴马汀"],
    "jatrorrhizine": ["药根碱"], "coptisine": ["黄连碱"],
    "epiberberine": ["表小檗碱"],
    "tetrandrine": ["粉防己碱", "汉防己甲素"],
    "fangchinoline": ["防己诺林碱"], "dauricine": ["蝙蝠葛碱"],
    "puerarin": ["葛根素"], "daidzein": ["大豆素"],
    "daidzin": ["大豆苷"], "baicalin": ["黄芩苷"],
    "baicalein": ["黄芩素"], "wogonin": ["汉黄芩素"],
    "wogonoside": ["汉黄芩苷"], "astragaloside": ["黄芪甲苷"],
    "glycyrrhizic acid": ["甘草酸"], "glycyrrhizin": ["甘草酸"],
    "liquiritin": ["甘草苷"], "paeoniflorin": ["芍药苷"],
    "albiflorin": ["白芍苷"], "aconitine": ["乌头碱"],
    "mesaconitine": ["新乌头碱", "中乌头碱"],
    "hypaconitine": ["次乌头碱"], "benzoylaconine": ["苯甲酰乌头原碱"],
    "benzoylmesaconine": ["苯甲酰新乌头原碱"],
    "matrine": ["苦参碱"], "oxymatrine": ["氧化苦参碱"],
    "sophocarpine": ["槐果碱"], "ginsenoside": ["人参皂苷"],
    "notoginsenoside": ["三七皂苷"], "salidroside": ["红景天苷"],
    "tetramethylpyrazine": ["川芎嗪"], "ligustrazine": ["川芎嗪"],
    "ferulic acid": ["阿魏酸"], "tetrahydropalmatine": ["延胡索乙素"],
    "corydaline": ["延胡索甲素", "紫堇碱"],
    "cryptotanshinone": ["隐丹参酮"], "tanshinone": ["丹参酮"],
    "danshensu": ["丹参素"], "salvianolic acid": ["丹酚酸"],
    "protocatechuic acid": ["原儿茶酸"], "hesperidin": ["橙皮苷"],
    "hesperetin": ["橙皮素"], "naringin": ["柚皮苷"],
    "neohesperidin": ["新橙皮苷"], "synephrine": ["辛弗林"],
    "ephedrine": ["麻黄碱"], "pseudoephedrine": ["伪麻黄碱"],
    "quercetin": ["槲皮素"], "kaempferol": ["山奈酚"],
    "rutin": ["芦丁"], "hyperoside": ["金丝桃苷"],
    "isoquercitrin": ["异槲皮苷"], "vitexin": ["牡荆素"],
    "emodin": ["大黄素"], "aloe-emodin": ["芦荟大黄素"],
    "rhein": ["大黄酸"], "chrysophanol": ["大黄酚"],
    "physcion": ["大黄素甲醚"], "ligustilide": ["藁本内酯"],
    "catalpol": ["梓醇"], "acteoside": ["毛蕊花糖苷"],
    "verbascoside": ["毛蕊花糖苷"], "echinacoside": ["松果菊苷"],
    "forsythiaside": ["连翘酯苷"], "phillyrin": ["连翘苷"],
    "arctiin": ["牛蒡苷"], "arctigenin": ["牛蒡苷元"],
    "chlorogenic acid": ["绿原酸"], "caffeic acid": ["咖啡酸"],
    "gallic acid": ["没食子酸"], "ellagic acid": ["鞣花酸"],
    "bergenin": ["岩白菜素"], "schisandrin": ["五味子素"],
    "deoxyschisandrin": ["五味子甲素"], "curcumin": ["姜黄素"],
    "gingerol": ["姜辣素"], "aucubin": ["桃叶珊瑚苷"],
    "geniposide": ["栀子苷"], "genipin": ["栀子苷元"],
    "ursolic acid": ["熊果酸"], "oleanolic acid": ["齐墩果酸"],
    "betulinic acid": ["白桦脂酸"], "sitosterol": ["谷甾醇"],
    "stigmasterol": ["豆甾醇"],
}

GOVERNANCE_BLOCKED_EXPORT_COLUMNS = [
    "Source_Result_Index", "Feature_ID", "Spectrum_ID", "Identity_Key", "Compound_Name",
    "Compound_Name_CN", "Formula", "Retention_Time", "Retention_Time_Unit",
    "Precursor_mz", "Cosine", "Claim_Status", "Claim_Scope", "Policy_ID",
    "Release_Eligible_Before_Replication", "Release_Blockers",
]


def zh_aliases_for(name: Any) -> list[str]:
    """Return exact-token aliases for presentation and search only."""
    text = str(name or "").lower()
    if not text:
        return []
    tokens = {token for token in re.split(r"[^a-z0-9]+", text) if token}
    tokens.add(text.strip())
    aliases: list[str] = []
    for key in sorted(ZH_COMPOUND_ALIASES, key=len, reverse=True):
        if key not in tokens:
            continue
        for alias in ZH_COMPOUND_ALIASES[key]:
            if alias not in aliases:
                aliases.append(alias)
        tokens.discard(key)
    return aliases


def _finite(value: Any) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def project_governance_blocked_candidates(
    results: Iterable[Mapping[str, Any]],
) -> list[dict[str, Any]]:
    """Project the two locked blocked statuses without changing source rows."""
    projected: list[dict[str, Any]] = []
    for fallback_index, row in enumerate(results or []):
        if not isinstance(row, Mapping):
            continue
        explicit_index = row.get("_result_index", row.get("result_index"))
        try:
            result_index = int(explicit_index)
        except (TypeError, ValueError):
            result_index = fallback_index
        match = row.get("match")
        if not isinstance(match, Mapping):
            continue
        passport = match.get("evidence_governance")
        if not isinstance(passport, Mapping):
            continue
        claim_status = str(passport.get("claim_status") or "")
        if claim_status not in BLOCKED_CLAIM_STATUSES:
            continue
        name = str(match.get("compound_name") or "")
        aliases = zh_aliases_for(name)
        blockers = list(passport.get("release_blockers") or [])
        mz = _finite(row.get("precursor_mz"))
        cosine = _finite(match.get("cosine"))
        feature_id = str(
            row.get("lcms_feature_id") or row.get("precursor_feature_id") or ""
        )
        identity_key = str(
            match.get("identity_key") or match.get("inchikey")
            or match.get("smiles") or f"name-formula::{name}|{match.get('formula') or ''}"
        )
        search_terms = [
            name, *aliases, str(match.get("formula") or ""), feature_id,
            str(row.get("spectrum_id") or ""), claim_status, *map(str, blockers),
        ]
        projected.append({
            "source_result_index": result_index,
            "source_feature_id": feature_id,
            "spectrum_id": str(row.get("spectrum_id") or ""),
            "identity_key": identity_key,
            "compound_name": name,
            "compound_name_cn": aliases,
            "formula": str(match.get("formula") or ""),
            "retention_time": _finite(row.get("retention_time")),
            "retention_time_unit": str(row.get("retention_time_unit") or ""),
            "precursor_mz": mz,
            "cosine": cosine,
            "claim_status": claim_status,
            "claim_scope": str(passport.get("claim_scope") or ""),
            "policy_id": str(passport.get("policy_id") or ""),
            "release_eligible_before_replication": bool(
                passport.get("release_eligible_before_replication")
            ),
            "release_blockers": blockers,
            "reason_text": "; ".join(map(str, blockers)),
            "search_text": " ".join(
                dict.fromkeys(str(value).strip().lower() for value in search_terms if value)
            ),
        })
    return projected


def governance_blocked_export_rows(
    rows: Iterable[Mapping[str, Any]],
) -> list[dict[str, Any]]:
    """Return stable flat rows for a CSV section or independent Excel sheet."""
    output = []
    for row in rows or []:
        output.append({
            "Source_Result_Index": row.get("source_result_index", ""),
            "Feature_ID": row.get("source_feature_id", ""),
            "Spectrum_ID": row.get("spectrum_id", ""),
            "Identity_Key": row.get("identity_key", ""),
            "Compound_Name": row.get("compound_name", ""),
            "Compound_Name_CN": "; ".join(row.get("compound_name_cn") or []),
            "Formula": row.get("formula", ""),
            "Retention_Time": row.get("retention_time", ""),
            "Retention_Time_Unit": row.get("retention_time_unit", ""),
            "Precursor_mz": row.get("precursor_mz", ""),
            "Cosine": row.get("cosine", ""),
            "Claim_Status": row.get("claim_status", ""),
            "Claim_Scope": row.get("claim_scope", ""),
            "Policy_ID": row.get("policy_id", ""),
            "Release_Eligible_Before_Replication": row.get(
                "release_eligible_before_replication", False
            ),
            "Release_Blockers": "; ".join(map(str, row.get("release_blockers") or [])),
        })
    return output


def append_governance_blocked_csv_section(
    csv_text: str,
    rows: Iterable[Mapping[str, Any]],
    *,
    language: str = "zh",
) -> str:
    """Append an isolated diagnostic section without changing existing rows."""
    flat_rows = governance_blocked_export_rows(rows)
    if not flat_rows:
        return csv_text
    buffer = io.StringIO()
    if csv_text:
        buffer.write(csv_text.rstrip("\r\n"))
        buffer.write("\n\n")
    writer = csv.writer(buffer, lineterminator="\n")
    title = (
        f"=== Governance not released: {len(flat_rows)} saved scan candidate records ==="
        if str(language).lower().startswith("en")
        else f"=== 治理未放行候选：{len(flat_rows)} 条保存的扫描候选记录 ==="
    )
    writer.writerow([title])
    writer.writerow(GOVERNANCE_BLOCKED_EXPORT_COLUMNS)
    for row in flat_rows:
        writer.writerow([row.get(column, "") for column in GOVERNANCE_BLOCKED_EXPORT_COLUMNS])
    return buffer.getvalue()


def render_governance_blocked_html(
    rows: Iterable[Mapping[str, Any]],
    *,
    language: str = "zh",
) -> str:
    """Render a closed, independently filtered diagnostic table."""
    records = list(rows or [])
    if not records:
        return ""
    english = str(language).lower().startswith("en")
    table_rows = []
    for row in records:
        aliases = "·".join(map(str, row.get("compound_name_cn") or []))
        alias_html = (
            f'<span class="meta-chip zh-alias">{escape(aliases)}</span>' if aliases else ""
        )
        mz = row.get("precursor_mz")
        rt = row.get("retention_time")
        cosine = row.get("cosine")
        search_text = escape(str(row.get("search_text") or ""), quote=True)
        table_rows.append(
            f'<tr class="governance-blocked-row" data-search="{search_text}" '
            f'data-mz="{escape(str(mz if mz is not None else ""), quote=True)}">'
            f'<td>{int(row.get("source_result_index") or 0)}</td>'
            f'<td>{escape(str(row.get("spectrum_id") or ""))}</td>'
            f'<td>{escape(str(row.get("source_feature_id") or ""))}</td>'
            f'<td><b>{escape(str(row.get("compound_name") or ""))}</b>{alias_html}</td>'
            f'<td>{escape(str(row.get("formula") or ""))}</td>'
            f'<td>{"—" if rt is None else escape(str(rt))}</td>'
            f'<td>{"—" if mz is None else escape(str(mz))}</td>'
            f'<td>{"—" if cosine is None else escape(f"{float(cosine):.4f}")}</td>'
            f'<td>{escape(str(row.get("claim_status") or ""))}</td>'
            f'<td>{escape(str(row.get("reason_text") or ""))}</td></tr>'
        )
    heading = (
        f"Governance not released: {len(records)} saved scan candidate records"
        if english else f"治理未放行候选 {len(records)} 条——按原判定保留，仅供排查（扫描候选记录数）"
    )
    placeholder = "Filter name, formula, m/z or reason" if english else "筛选名称、中文名、分子式、m/z 或原因"
    headers = (
        ["Index", "Spectrum", "Feature", "Candidate", "Formula", "RT", "Precursor m/z", "Cosine", "Claim status", "Original blockers"]
        if english else ["索引", "谱图", "特征", "候选名称", "分子式", "RT", "前体 m/z", "余弦", "治理状态", "原始拦截原因"]
    )
    return f'''<details class="governance-blocked-panel" style="margin-top:16px">
<summary>{escape(heading)}</summary>
<div style="margin:8px 0"><input id="governanceBlockedSearch" type="search" placeholder="{escape(placeholder, quote=True)}" oninput="filterGovernanceBlockedRows(this.value)"><span id="governanceBlockedSearchCount"></span></div>
<div class="table-wrap"><table><thead><tr>{''.join(f'<th>{escape(label)}</th>' for label in headers)}</tr></thead>
<tbody id="governanceBlockedBody">{''.join(table_rows)}</tbody></table></div>
</details>
<script>
function filterGovernanceBlockedRows(value) {{
  var query = String(value || '').trim().toLowerCase();
  var numeric = /^[0-9]+([.][0-9]+)?$/.test(query);
  var queryMz = numeric ? parseFloat(query) : NaN;
  var rows = Array.from(document.querySelectorAll('#governanceBlockedBody .governance-blocked-row'));
  var visible = 0;
  rows.forEach(function(row) {{
    var mz = parseFloat(row.dataset.mz || '');
    var show = !query || (numeric ? Number.isFinite(mz) && Math.abs(mz - queryMz) <= 0.5 : String(row.dataset.search || '').includes(query));
    row.style.display = show ? '' : 'none';
    if (show) visible += 1;
  }});
  var count = document.getElementById('governanceBlockedSearchCount');
  if (count) count.textContent = query ? ' ' + visible + '/' + rows.length : '';
}}
</script>'''


def inject_governance_blocked_html(
    html_text: str,
    rows: Iterable[Mapping[str, Any]],
    *,
    language: str = "zh",
) -> str:
    """Insert the diagnostic projection before the document body closes."""
    section = render_governance_blocked_html(rows, language=language)
    if not section:
        return html_text
    if "</body>" in html_text:
        return html_text.replace("</body>", section + "</body>", 1)
    return html_text + section


__all__ = [
    "BLOCKED_CLAIM_STATUSES", "ZH_COMPOUND_ALIASES",
    "GOVERNANCE_BLOCKED_EXPORT_COLUMNS", "zh_aliases_for",
    "project_governance_blocked_candidates", "governance_blocked_export_rows",
    "append_governance_blocked_csv_section", "render_governance_blocked_html",
    "inject_governance_blocked_html",
]
