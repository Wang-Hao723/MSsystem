"""Separate reviewed external candidates from a plant candidate reading table."""
from pipeline.compound_role_projection import apply_report_compound_roles, plant_report_selection


def split_plant_reader_rows(rows, *, language="zh"):
    """Legacy rows: one excluded member places its complete original group in QC."""
    plant, review = [], []
    grouped = {}
    for original in rows:
        identity = {
            "compound_name": original.get("Candidate_Compound_Name", ""),
            "smiles": original.get("Candidate_SMILES", ""),
            "inchikey": original.get("Candidate_InChIKey", ""),
        }
        apply_report_compound_roles(identity)
        selection = plant_report_selection(identity)
        key = original.get("Group_Number") or original.get("Spectrum")
        grouped.setdefault(key, []).append((original, selection))
    for members in grouped.values():
        excluded = any(selection["excluded_from_plant_candidate_table"] for _, selection in members)
        for original, selection in members:
            row = _annotate_row(original, selection, excluded, language)
            (review if excluded else plant).append(row)
    return plant, review


def _annotate_row(original, selection, excluded, language):
    row = dict(original)
    if not selection:
        identity = {"compound_name": row.get("Candidate_Compound_Name", ""),
                    "smiles": row.get("Candidate_SMILES", ""),
                    "inchikey": row.get("Candidate_InChIKey", "")}
        apply_report_compound_roles(identity)
        selection = plant_report_selection(identity)
    if excluded:
        note = ("Complete group excluded from plant candidate table; retained for review"
                if language == "en" else "原组排除于植物候选表；完整保留供外源或身份冲突复核")
    else:
        note = ("Plant candidate reading table; passing this filter does not verify plant origin"
                if language == "en" else "植物候选阅读表；通过此筛选不等于已核验植物来源")
    note += f"; role={selection.get('role', 'group')}; reason={selection.get('reason', 'group_source_selection')}"
    if selection.get("role_annotation_source"):
        note += f"; source={selection['role_annotation_source']}"
    note += ("; original candidate rank, group size and identity uncertainty retained"
             if language == "en" else "；保留原候选排名、原组内候选数及异构体未解析限定")
    row["Export_Note"] = "; ".join(filter(None, [str(row.get("Export_Note") or ""), note]))
    return row


def plant_reader_views(matched, *, compounds=None, language="zh"):
    """Use the same persisted compound-group selection as the main API."""
    from pipeline.compound_conclusions import filter_compound_conclusions
    from web.report_common import ms2_reader_rows
    if compounds is None:
        return split_plant_reader_rows(ms2_reader_rows(matched, language=language), language=language)
    views = []
    for scope in ("plant_candidates", "external_review"):
        selected = filter_compound_conclusions({"metadata": {}, "compounds": compounds}, report_scope=scope)
        rows = ms2_reader_rows(matched, compounds=selected["compounds"], language=language)
        views.append([_annotate_row(row, {}, scope == "external_review", language) for row in rows])
    return tuple(views)
