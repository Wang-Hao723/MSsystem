from pathlib import Path
import json

import matplotlib.pyplot as plt
from matplotlib import font_manager
from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER, TA_LEFT
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import mm
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.cidfonts import UnicodeCIDFont
from reportlab.platypus import (
    SimpleDocTemplate, Paragraph, Spacer, PageBreak, Table, TableStyle,
    Image, KeepTogether, LongTable, HRFlowable
)

ROOT = Path(__file__).resolve().parent.parent
OUT_DIR = ROOT.parent / "output" / "pdf"
TMP = ROOT.parent / "tmp" / "pdfs" / "latest_report"
OUT_DIR.mkdir(parents=True, exist_ok=True)
TMP.mkdir(parents=True, exist_ok=True)
OUT = OUT_DIR / "民族药LCMS_AI结构解析最新研究报告_20260713.pdf"

BLUE = colors.HexColor("#123F92")
BLUE2 = colors.HexColor("#1A5BC4")
PALE = colors.HexColor("#EAF2FF")
CYAN = colors.HexColor("#E6F6FA")
YELLOW = colors.HexColor("#FFF3C4")
RED = colors.HexColor("#C62828")
GREEN = colors.HexColor("#168A55")
GRAY = colors.HexColor("#5E6B7A")
LIGHT = colors.HexColor("#F6F8FB")

pdfmetrics.registerFont(UnicodeCIDFont("STSong-Light"))
FONT = "STSong-Light"

def set_mpl_font():
    candidates = [
        "/System/Library/Fonts/STHeiti Light.ttc",
        "/System/Library/AssetsV2/com_apple_MobileAsset_Font7/3419f2a427639ad8c8e139149a287865a90fa17e.asset/AssetData/PingFang.ttc",
        "/System/Library/Fonts/ヒラギノ角ゴシック W3.ttc",
    ]
    for p in candidates:
        if Path(p).exists():
            font_manager.fontManager.addfont(p)
            plt.rcParams["font.family"] = font_manager.FontProperties(fname=p).get_name()
            break
    plt.rcParams["axes.unicode_minus"] = False

set_mpl_font()

styles = getSampleStyleSheet()
styles.add(ParagraphStyle(name="CN", fontName=FONT, fontSize=9.2, leading=15, textColor=colors.HexColor("#26364A"), spaceAfter=5))
styles.add(ParagraphStyle(name="SmallCN", fontName=FONT, fontSize=7.5, leading=11, textColor=GRAY))
styles.add(ParagraphStyle(name="TitleCN", fontName=FONT, fontSize=24, leading=34, textColor=BLUE, alignment=TA_CENTER, spaceAfter=12))
styles.add(ParagraphStyle(name="SubTitleCN", fontName=FONT, fontSize=12, leading=19, textColor=GRAY, alignment=TA_CENTER))
styles.add(ParagraphStyle(name="H1CN", fontName=FONT, fontSize=16, leading=22, textColor=BLUE, spaceBefore=4, spaceAfter=9))
styles.add(ParagraphStyle(name="H2CN", fontName=FONT, fontSize=12, leading=18, textColor=BLUE2, spaceBefore=5, spaceAfter=5))
styles.add(ParagraphStyle(name="Callout", fontName=FONT, fontSize=10, leading=16, textColor=colors.HexColor("#5B3B00")))
styles.add(ParagraphStyle(name="FootCN", fontName=FONT, fontSize=7, leading=9, textColor=GRAY))

def P(text, style="CN"):
    return Paragraph(text, styles[style])

def section(title):
    return [Spacer(1, 2*mm), HRFlowable(width="100%", thickness=2, color=BLUE), Spacer(1, 2*mm), P(title, "H1CN")]

def table(data, widths, header=True, font=7.3, repeat=1):
    t = LongTable([[P(str(c), "SmallCN") for c in row] for row in data], colWidths=widths, repeatRows=repeat if header else 0, hAlign="LEFT")
    cmd = [
        ("FONTNAME", (0,0), (-1,-1), FONT), ("FONTSIZE", (0,0), (-1,-1), font),
        ("LEADING", (0,0), (-1,-1), font+3), ("VALIGN", (0,0), (-1,-1), "MIDDLE"),
        ("GRID", (0,0), (-1,-1), 0.35, colors.HexColor("#B9C7DA")),
        ("LEFTPADDING", (0,0), (-1,-1), 4), ("RIGHTPADDING", (0,0), (-1,-1), 4),
        ("TOPPADDING", (0,0), (-1,-1), 4), ("BOTTOMPADDING", (0,0), (-1,-1), 4),
    ]
    if header:
        cmd += [("BACKGROUND", (0,0), (-1,0), BLUE), ("TEXTCOLOR", (0,0), (-1,0), colors.white)]
    for r in range(1 if header else 0, len(data)):
        if r % 2 == 0: cmd.append(("BACKGROUND", (0,r), (-1,r), LIGHT))
    t.setStyle(TableStyle(cmd))
    return t

def callout(title, body, color=YELLOW):
    t=Table([[P(f"<b>{title}</b><br/>{body}", "Callout")]], colWidths=[174*mm])
    t.setStyle(TableStyle([("BACKGROUND",(0,0),(-1,-1),color),("BOX",(0,0),(-1,-1),0.8,colors.HexColor("#D5B957")),("LEFTPADDING",(0,0),(-1,-1),9),("RIGHTPADDING",(0,0),(-1,-1),9),("TOPPADDING",(0,0),(-1,-1),7),("BOTTOMPADDING",(0,0),(-1,-1),7)]))
    return t

def metric_cards(cards):
    cells=[]
    for value,label,color in cards:
        cells.append(P(f'<font size="17" color="{color}"><b>{value}</b></font><br/><font size="8">{label}</font>', "CN"))
    t=Table([cells], colWidths=[174*mm/len(cells)]*len(cells), rowHeights=[24*mm])
    t.setStyle(TableStyle([("BACKGROUND",(0,0),(-1,-1),PALE),("GRID",(0,0),(-1,-1),0.6,colors.HexColor("#B7C9E8")),("VALIGN",(0,0),(-1,-1),"MIDDLE"),("ALIGN",(0,0),(-1,-1),"CENTER")]))
    return t

def make_charts():
    fig, ax=plt.subplots(figsize=(7.2,3.5))
    labels=["优化前 Top-1","优化后 Top-1","优化前 Top-5","优化后 Top-5"]
    vals=[69.7,87.9,81.8,93.9]
    bars=ax.bar(labels,vals,color=["#9DB7DF","#1647D8","#F3C969","#168A55"],width=.62)
    ax.set_ylim(0,105); ax.set_ylabel("准确率 / 召回率（%）"); ax.grid(axis="y",alpha=.25)
    for b,v in zip(bars,vals): ax.text(b.get_x()+b.get_width()/2,v+2,f"{v:.1f}%",ha="center",fontweight="bold")
    ax.set_title("mix33 严格评估：优化前后对比")
    fig.tight_layout(); fig.savefig(TMP/"retrieval_improvement.png",dpi=190,bbox_inches="tight"); plt.close(fig)

    fig, ax=plt.subplots(figsize=(7.2,3.5))
    labels=["有匹配MS/MS谱","去重候选化合物","Top-1正确标品","Top-5召回标品"]
    vals=[3123,862,29,31]
    bars=ax.bar(labels,vals,color=["#1647D8","#55A4D9","#F3A13B","#168A55"])
    ax.set_yscale("log"); ax.set_ylabel("数量（对数坐标）"); ax.grid(axis="y",alpha=.25)
    for b,v in zip(bars,vals): ax.text(b.get_x()+b.get_width()/2,v*1.16,str(v),ha="center",fontweight="bold")
    ax.set_title("“总匹配量”与“标准品准确度”是不同口径")
    fig.tight_layout(); fig.savefig(TMP/"counts_scope.png",dpi=190,bbox_inches="tight"); plt.close(fig)

    fig, ax=plt.subplots(figsize=(7.2,3.5))
    labels=["Tanimoto@1","Exact@1","Tan≥0.6","Valid Ratio","|Δmass|≤5 Da","Oracle Recall@10"]
    vals=[37.2,15.0,26.5,100,41,18]
    bars=ax.bar(labels,vals,color=["#1647D8","#F28C28","#55A4D9","#168A55","#D8A72B","#C94C4C"])
    ax.set_ylim(0,108); ax.set_ylabel("比例（%）"); ax.grid(axis="y",alpha=.25); ax.tick_params(axis="x",rotation=18)
    for b,v in zip(bars,vals): ax.text(b.get_x()+b.get_width()/2,v+2,f"{v:g}%",ha="center",fontsize=9)
    ax.set_title("生成式模型独立200谱图评估（Beam=10）")
    fig.tight_layout(); fig.savefig(TMP/"generative_metrics.png",dpi=190,bbox_inches="tight"); plt.close(fig)

    fig, ax=plt.subplots(figsize=(7.2,3.5))
    names=["水晶兰苷","牛膝皂苷I","β-细辛醚","细叶远志皂苷"]
    rank=[6,3,2,6]
    colors_=["#C94C4C","#F3A13B","#F3A13B","#C94C4C"]
    ax.barh(names,rank,color=colors_); ax.invert_yaxis(); ax.set_xlim(0,7)
    ax.axvline(1,color="#168A55",ls="--",label="Top-1"); ax.axvline(5,color="#1647D8",ls=":",label="Top-5边界")
    ax.set_xlabel("正确候选名次（6表示未进入Top-5）"); ax.legend(loc="lower right"); ax.grid(axis="x",alpha=.2)
    ax.set_title("4个非Top-1标品的候选位置")
    fig.tight_layout(); fig.savefig(TMP/"failure_ranks.png",dpi=190,bbox_inches="tight"); plt.close(fig)

make_charts()

eval_data=json.load(open(ROOT/"eval/v5output_33std_beta_massfix_strict.json"))
train_curve=Path("/Users/zizi/Desktop/MsQCModel_v2/logs_zmem6base_s777/training_curves.png")

def header_footer(canvas, doc):
    canvas.saveState()
    canvas.setStrokeColor(BLUE); canvas.setLineWidth(0.7); canvas.line(20*mm, 18*mm, 190*mm, 18*mm)
    canvas.setFont(FONT, 7); canvas.setFillColor(GRAY)
    canvas.drawString(20*mm, 11*mm, "民族药资源库建设与应用 | LC-MS/MS AI结构解析最新研究报告")
    canvas.drawRightString(190*mm, 11*mm, f"第 {doc.page} 页")
    canvas.restoreState()

doc=SimpleDocTemplate(str(OUT),pagesize=A4,rightMargin=18*mm,leftMargin=18*mm,topMargin=17*mm,bottomMargin=23*mm,title="LC-MS/MS AI结构解析最新研究报告")
story=[]

# Cover
story += [Spacer(1,27*mm), P("最新研究报告", "TitleCN"), P("民族药资源库建设与应用：LC-MS/MS检索式与生成式AI结构解析", "SubTitleCN"), Spacer(1,15*mm)]
story.append(callout("报告口径", "FULL谱库生产重跑 + 33标准品严格评估 + 生成式独立测试。检索结果分别报告总匹配谱图、去重候选、Top-1准确度和Top-5召回率，不再用“只识别29个”概括系统表现。", PALE))
story += [Spacer(1,20*mm), metric_cards([("3,123","有匹配的MS/MS谱","#1647D8"),("862","去重候选化合物","#1647D8"),("29/33","Top-1正确","#C62828"),("31/33","Top-5召回","#168A55")]), Spacer(1,18*mm)]
story += [P("报告日期：2026-07-13", "SubTitleCN"), P("系统：mssystem_v5（FULL library 1,330,239 spectra）", "SubTitleCN"), P("生成模型：MsQCModel_v2 阶段检查点与独立生成评估", "SubTitleCN"), PageBreak()]

# Executive summary
story += section("一、执行摘要与本次更正")
story.append(P("本次复核的核心更正是：<b>29/33是严格Top-1准确度，不是系统总共只检索出29个化合物</b>。FULL库生产任务在正、负离子模式共输出3,123条有匹配的MS/MS谱，去重后形成862个候选化合物。33个已知标品中，29个正确结构排在第一，31个在前五名内；真正未进入Top-5的只有水晶兰苷和细叶远志皂苷。"))
story.append(Image(str(TMP/"counts_scope.png"),width=166*mm,height=80*mm))
story.append(callout("一句话结论", "mssystem已达到较强的已知物谱库检索水平；当前瓶颈集中于同分异构体排序、跨仪器/碰撞能域偏移，以及少数目标在对应离子模式下缺少高质量参考谱。"))
story += [Spacer(1,4*mm), table([
    ["问题","审计结论","是否已修复/下一步"],
    ["“只分析出27/29个”","混淆了Top-k标准品准确度与总匹配量","已更正报告口径"],
    ["β-细辛醚质量信息","库中前体值仅写209.11且加合物不明确","已修为209.11722、[M+H]+"],
    ["水晶兰苷","负离子实测谱被同分子式候选压制","需同仪器多能量单标谱+RT"],
    ["细叶远志皂苷","正离子[M+NH4]+峰已采集，但仅返回同分子式异构体","补同仪器正/负离子多能量实验谱与RT"],
], [38*mm,76*mm,60*mm]), PageBreak()]

# Background
story += section("二、项目背景与两条技术路线")
story.append(P("项目面向民族药复杂体系中的化学成分发现与结构注释，形成两条互补路线：检索式用于已知结构的高通量去重复和可信注释；生成式用于谱库未覆盖结构的候选生成。两条路线不应使用同一评价口径，也不能把候选生成等同于结构确证。"))
route_data=[
    ["维度","检索式 Database Search","生成式 Generative AI"],
    ["输入","实验MS/MS、前体m/z、极性、加合物、RT","实验MS/MS，可选分子式/质量/类别先验"],
    ["输出","谱库中的Top-k候选与匹配证据","新结构候选、SELFIES/SMILES/分子图"],
    ["优势","速度快、可追溯、已知物准确度高","可提出库外候选，覆盖未知化学空间"],
    ["主要风险","库覆盖、异构体、跨仪器域偏移","化学有效但结构错误、质量漂移、置信度失准"],
    ["本项目定位","当前生产主线","研究型补充路线"],
]
story.append(table(route_data,[34*mm,70*mm,70*mm]))
story += [Spacer(1,6*mm), P("推荐的混合工作流", "H2CN")]
flow=[P("谱图质控<br/>去噪/空白/极性","CN"),P("检索式Top-k<br/>cosine+entropy","CN"),P("生成式补充<br/>质量/分子式约束","CN"),P("预测谱回打分<br/>RT/分类/诊断离子","CN"),P("专家复核<br/>标准品确认","CN")]
ft=Table([flow],colWidths=[34.8*mm]*5,rowHeights=[24*mm])
ft.setStyle(TableStyle([("BACKGROUND",(0,0),(-1,-1),PALE),("BOX",(0,0),(-1,-1),0.7,BLUE),("INNERGRID",(0,0),(-1,-1),0.5,colors.white),("ALIGN",(0,0),(-1,-1),"CENTER"),("VALIGN",(0,0),(-1,-1),"MIDDLE")]))
story.append(ft)
story += [Spacer(1,7*mm), callout("证据分级", "标准品+RT一致可作为Level 1；实验谱库匹配通常为Level 2；异构体组或分类注释应降至Level 3；生成式候选在正交验证前只属于结构假设。", CYAN), PageBreak()]

# Data and method
story += section("三、数据、系统与严格评价方法")
story.append(table([
    ["项目","本次使用数据"],
    ["真实样本","MIX_33_POS2.mzML、MIX_33_NEG2.mzML"],
    ["目标清单","33个标准品，包含β-细辛醚"],
    ["谱图库","FULL：加载后1,330,239条可检索谱"],
    ["检索阈值","cosine≥0.30、final≥0.25、Top15峰、前体窗口0.5 Da；严格证据门关闭"],
    ["评价门控","目标前体窗口内，名称/公认同义名与分子式同时一致"],
    ["生产任务","POS 376bd31ccdec；NEG 3af46d212955"],
], [48*mm,126*mm]))
story += [Spacer(1,5*mm), P("严格评价为什么重要", "H2CN"), P("旧的名称子串匹配会把同名片段、泛化皂苷名或错误分子式计为正确。新版评估要求名称/同义名和分子式一致，并只在目标前体窗口内寻找；Top-1和Top-5分别报告。该口径更严格，也更接近可复核的化合物级鉴定。"), Spacer(1,5*mm)]
story.append(callout("边界", "关闭严格证据门是为了恢复用户设置阈值的真实含义，增加的低分结果仍只能作为Putative候选；正式报告需结合空白、RT和标准品验证。"))
story += [Spacer(1,8*mm), P("数据库规模", "H2CN"), table([
    ["数据库层","谱图/化合物规模","用途"],
    ["Small","约62.7万条可检索谱","交互探索、快速试跑"],
    ["Full","133.0万条可检索谱","正式分析与最大覆盖"],
    ["目标化合物元数据","名称、分子式、SMILES/InChIKey、分类、裂解通路","过滤、重排和解释"],
], [42*mm,62*mm,70*mm]), PageBreak()]

# Retrieval results
story += section("四、检索式结果：Top-1 29/33，Top-5 31/33")
story.append(metric_cards([("2,228","POS有匹配谱","#1647D8"),("895","NEG有匹配谱","#1647D8"),("87.9%","Top-1准确度","#C62828"),("93.9%","Top-5召回率","#168A55")]))
story += [Spacer(1,5*mm), Image(str(TMP/"retrieval_improvement.png"),width=166*mm,height=80*mm)]
story.append(P("优化前严格重算为Top-1 23/33（69.7%）、Top-5 27/33（81.8%）；补库、名称/精确质量修正、生产重跑和候选重排后，Top-1提高18.2个百分点，Top-5提高12.1个百分点。"))
story += [Spacer(1,5*mm), table([
    ["指标","优化前","优化后","变化"],
    ["Top-1正确","23/33（69.7%）","29/33（87.9%）","+6个；+18.2个百分点"],
    ["Top-5召回","27/33（81.8%）","31/33（93.9%）","+4个；+12.1个百分点"],
    ["有匹配MS/MS谱","未统一FULL重跑口径","3,123","POS 2,228；NEG 895"],
    ["去重候选","未统一口径","862","InChIKey优先、名称兜底"],
], [40*mm,43*mm,43*mm,48*mm]), PageBreak()]

# Per standard
story += section("五、33个标准品逐项结果")
rows=[["序号","标品","英文名","Top-1","Top-5","系统结果/名次"]]
for i,r in enumerate(eval_data["per_standard"],1):
    top1="正确" if r.get("hit") else "否"
    top5="召回" if r.get("hit_top5") else "未召回"
    detail=r.get("v5_name") or (f"{r.get('top5_name')}（Top-{r.get('top5_rank')}）" if r.get("top5_name") else "未进入Top-5")
    rows.append([i,r.get("zh"),r.get("name"),top1,top5,detail])
story.append(table(rows,[9*mm,27*mm,42*mm,18*mm,20*mm,58*mm],font=6.7))
story += [PageBreak()]

# Failures
story += section("六、非Top-1标品的逐项诊断")
story.append(Image(str(TMP/"failure_ranks.png"),width=166*mm,height=80*mm))
fail=[
    ["标品","当前候选位置","可核查根因","解决方案"],
    ["牛膝皂苷I（Cyaonoside B）","Top-3","正确[M+NH4]+参考谱已存在；同分子式三萜皂苷碎裂接近，另一异构体获得更高候选特异碎片分","补同仪器多碰撞能谱；建立糖链中性丢失和皂苷元诊断离子；加入单标RT"],
    ["β-细辛醚","Top-2","仅1条参考谱、4个有效碎片；库谱与实测碰撞能覆盖不一致；旧前体值仅209.11","已修前体209.11722/[M+H]+；仍需补3-5档碰撞能标准谱"],
    ["水晶兰苷（Monotropein）","未进Top-5","有负离子参考谱，但实测谱更接近同分子式去乙酰车叶草酸；属于真实位置/构型异构体竞争","同仪器负离子标准谱、RT、能量梯度或MS3；无法区分时报告候选组"],
    ["细叶远志皂苷（Tenuifolin）","未进Top-5","正离子文件已采到[M+NH4]+ m/z 698.414734（RT约12.13 min），并返回3个C36H56O12候选，但InChIKey均非Tenuifolin；正确结构因同分异构体竞争及跨加合物/碰撞能差异未进入候选。负离子覆盖不足是附加因素","同时使用正、负离子窗口；补同仪器[M+NH4]+/[M+H]+/[M-H]-多碰撞能单标谱和RT；增加糖链/皂苷元诊断规则"],
]
story.append(table(fail,[29*mm,25*mm,62*mm,58*mm],font=6.8))
story += [Spacer(1,6*mm), callout("不能采用的“优化”", "不能把当前mix33谱反灌数据库后再用同一批数据计算准确率，也不能针对33个名称设置硬编码奖励。这两种做法会造成训练-测试泄漏或面板过拟合。"), PageBreak()]

# DB quality
story += section("七、数据库信息完整性审计与已实施修复")
story.append(P("数据库大并不自动等于目标可识别。真正决定召回和排序的是目标结构在正确离子模式、加合物、仪器类型和碰撞能下是否有高质量参考谱，以及名称、分子式、结构标识是否一致。"))
story.append(table([
    ["字段/证据","作用","当前问题","处理状态"],
    ["compound_name / synonyms","名称归一化和结果展示","历史存在拼写、泛名和长IUPAC名","已修Rehmannioside D；严格同义名表"],
    ["formula / exact_mass","质量门控和错误候选排除","Aucubin曾有错误精确质量；部分记录缺失","已修目标错误；需全库一致性扫描"],
    ["adduct / polarity","选择正确前体和离子模式","β-细辛醚原标记Positive，前体仅两位小数","已改[M+H]+、209.11722"],
    ["InChIKey / SMILES","结构去重、异构体识别","部分公共谱缺结构或立体层不一致","加载时过滤无结构记录；需结构标准化"],
    ["collision energy / instrument","处理跨仪器域偏移","目标谱覆盖极不均匀","优先补本仪器多能量标准谱"],
    ["fragmentation_pathway","诊断离子和裂解重排","皂苷/环烯醚萜规则仍不足","建立类别专属规则库"],
    ["classification","类别先验和冲突降权","公共库标签不一致","采用v3分类与可靠类白名单"],
], [32*mm,42*mm,57*mm,43*mm],font=6.7))
story += [Spacer(1,6*mm), P("建议的数据库准入规则", "H2CN"), P("每条实验参考谱至少应具备：来源、许可证、仪器类型、碰撞能、极性、加合物、前体m/z、分子式、精确质量、SMILES、InChIKey和分类标签。缺少关键元数据的谱可以保留用于探索，但应降低可信度，不能无条件压过元数据完整的实验标准谱。"), PageBreak()]

# Generative training
story += section("八、生成式模型：训练收敛与独立生成能力")
story.append(P("阶段训练日志显示模型能够稳定学习序列，但teacher-forcing下的token/sequence准确率明显高于自由生成Exact，说明曝光偏差和质量约束仍是主要瓶颈。"))
if train_curve.exists():
    story.append(Image(str(train_curve),width=166*mm,height=89*mm))
story.append(table([
    ["Epoch 4指标","数值","解释"],
    ["Validation loss","0.536","验证损失已明显下降"],
    ["Token accuracy","82.1%","teacher-forcing逐token正确率"],
    ["Sequence accuracy","81.2%","训练日志定义的序列指标"],
    ["Teacher-forced Exact","7.3%","仍不能代表自由生成"],
    ["Generated Exact","0%","阶段检查点自由生成尚未形成精确匹配"],
    ["Generated Tanimoto","0.257","能学习部分结构片段"],
    ["Molecule valid","100%","语法/价态有效，不等于结构正确"],
], [45*mm,35*mm,94*mm]))
story += [Spacer(1,4*mm), P("数据来源：MsQCModel_v2/logs_zmem6base_s777/train_summary.json、train_metrics.json、training_curves.png。该页为阶段快照。", "FootCN"), PageBreak()]

# Generative metrics
story += section("九、生成式模型：独立200谱图评估")
story.append(Image(str(TMP/"generative_metrics.png"),width=166*mm,height=80*mm))
story.append(metric_cards([("0.372","Tanimoto@1","#1647D8"),("15.0%","Exact@1","#C62828"),("26.5%","Tan≥0.6","#F28C28"),("100%","Valid Ratio","#168A55")]))
story += [Spacer(1,5*mm), table([
    ["质量/召回指标","结果","判断"],
    ["质量误差均值 / 中位数","34.2 / 12.0 Da","质量约束仍弱"],
    ["|Δmass|≤5 Da / ≤0.1 Da","41% / 19%","高精度质量一致率偏低"],
    ["Oracle Recall@10","18%","候选池本身仍常缺正确结构"],
    ["ECE","0.675","置信度尚未校准"],
], [55*mm,45*mm,74*mm])]
story.append(callout("mix33生成式验证状态", "当前没有可复核的mix33生成推理输出，因此正确生成数量应标记为N/A，而不是0或任意推测值。必须按固定谱图、加合物和分子式可见性协议运行Top-1/5/10后再报告。"))
story += [Spacer(1,4*mm), P("数据来源：MsQCModel_v2/benchmark_id_gen_eval_zmem2.json。", "FootCN"), PageBreak()]

# Ablation compare
story += section("十、消融实验与生成式前沿对比")
story.append(table([
    ["实验/模块","变化","结果","可支持结论"],
    ["内部历史v2→v7（300查询）","模型和解码策略演进","Tan@1 0.295→0.325；Exact@1 0.7%→3.7%；Exact@10 4.3%→10.3%","生成质量有持续提升，但与最新200谱独立口径不可混算"],
    ["Mayer键级特征（QC代理任务）","GIN/GINE加入Mayer信息","验证MAE 0.1145→0.1056，降低约7.8%","键级信息有助于性质预测；不是结构Exact消融"],
    ["MARLIN conditioning diversity","移除条件多样性","Top-1 16.94%→12.83%；Top-10 23.54%→14.94%","多条件训练显著贡献Top-k"],
], [45*mm,43*mm,50*mm,36*mm],font=6.8))
story += [Spacer(1,6*mm), table([
    ["模型","任务/先验","公开结果或定位","与本项目关系"],
    ["本模型","内部独立200谱","Exact@1 15.0%；Tan@1 0.372；Valid 100%","质量误差和Top-k召回仍需提升"],
    ["MARLIN（2026预印本）","formula-free / NPLIB1","Top-1 16.94%；Top-10 23.54%","Top-1接近；需同基准复现"],
    ["DiffMS（ICML 2025）","已知分子式图扩散","论文报告SOTA","分子式先验更强，不能直接横比"],
    ["MSGo（2025）","特定含氟化合物域","报告48%结构验证","领域与评价口径不同"],
    ["CFM-ID","谱图预测/候选排序","不是无约束de novo生成器","适合作为预测谱回打分模块"],
], [37*mm,45*mm,48*mm,44*mm],font=6.8)]
story += [Spacer(1,4*mm), P("跨数据集数字只用于技术定位，不作为排行榜。外部来源：MassSpecGym（NeurIPS 2024）、DiffMS（ICML 2025）、MSGo（Nature Machine Intelligence 2025）、MARLIN arXiv:2607.04774。", "FootCN"), PageBreak()]

# Retrieval compare
story += section("十一、检索式外部工具对照：当前证据与公平实验设计")
story.append(P("用户要求在相同数据集上与GNPS、NIST、MassBank比较。当前只有mssystem在锁定的mix33生产数据上完成了可复核运行；其他工具没有对应运行日志，因此本报告不填虚构数值。"))
story.append(table([
    ["工具","mix33同数据实跑","Top-1","Top-5","当前结论"],
    ["mssystem_v5 / FULL","已完成","29/33（87.9%）","31/33（93.9%）","强已知库检索；异构体和缺谱仍受限"],
    ["GNPS Library Search","未完成","-","-","需统一MGF、容差、库版本和Top-k"],
    ["MassBank/MoNA","未完成","-","-","需固定实验谱快照、去重口径"],
    ["NIST MS Search","未完成","-","-","需确认商业库版本和授权环境"],
], [42*mm,30*mm,31*mm,31*mm,40*mm]))
story += [Spacer(1,7*mm), P("公平比较协议", "H2CN"), P("锁定每个标准品的代表性MS/MS谱并导出统一MGF；统一前体/碎片容差、最小峰数、强度阈值、离子模式、Top-k和同义名/分子式门控；所有工具盲评；报告Top-1、Top-5、未返回率、运行时间和bootstrap 95%置信区间。"), Spacer(1,5*mm)]
story.append(callout("水平判断", "mssystem在已知库覆盖场景已达到较强去重复和注释水平；内部检索基准已知结构Exact@1为95.1%，但novel Exact@1=0，说明纯检索无法解决库外新结构。"))
story += [Spacer(1,5*mm), P("MassBank是开放、FAIR的质谱资源，公开数据以实验谱为主并包含正负离子和多种仪器元数据；本项目应优先补充可追溯实验谱，而不是用预测谱冒充Level 1/2证据。", "FootCN"), PageBreak()]

# Plan
story += section("十二、针对性优化计划与验收标准")
story.append(table([
    ["优先级","任务","具体动作","验收指标"],
    ["P0","补齐4个困难标品证据","单标同仪器正负模式；10/20/40/60 eV；记录RT；必要时MS3","4个目标均有匹配模式≥4档碰撞能实验谱"],
    ["P0","提升β-细辛醚","补[M+H]+多能量谱，复核209.117与209.154两个特征的归属","独立重复样本Top-1稳定"],
    ["P0","解决两类异构体","水晶兰苷/去乙酰车叶草酸、远志皂苷异构体建立诊断离子和中性丢失","正确候选进入Top-5；不能区分时输出候选组"],
    ["P1","数据库全库质控","formula-exact mass-adduct-InChIKey一致性；重复谱按能量/仪器分层","关键字段完整率>95%；错误记录有审计日志"],
    ["P1","多证据重排","entropy、modified cosine、reverse、RT、诊断离子、类别冲突和校准概率","留出标准品Top-1≥90%，Top-5≥95%，且空白假阳性不升高"],
    ["P2","外部工具同数据比较","GNPS/MassBank/NIST统一协议盲评","形成同数据、同口径对照表和置信区间"],
    ["P2","生成式mix33验证","固定公式可见性，Top-1/5/10；Exact/Tan/Valid/质量误差","首次形成可复核mix33生成基线"],
], [15*mm,34*mm,79*mm,46*mm],font=6.7))
story += [Spacer(1,6*mm), P("预期结果与科学边界", "H2CN"), P("在补齐同仪器单标谱和RT后，Top-5达到33/33具有现实可行性；Top-1预计可提升至31-33/33。但如果异构体在当前MS/MS条件下不可区分，系统应诚实输出候选组，而不是强制选择一个名称。准确度提升必须在独立重复或留一碰撞能数据上验证。"), PageBreak()]

# Conclusions/references
story += section("十三、结论")
story.append(metric_cards([("862","去重候选化合物","#1647D8"),("29/33","Top-1正确","#C62828"),("31/33","Top-5召回","#168A55"),("2个","真正未进Top-5","#C62828")]))
story += [Spacer(1,7*mm), P("1. mssystem并非只检索出29个化合物；29/33是严格Top-1标准品准确度。生产任务实际得到3,123条匹配谱和862个去重候选。"), P("2. 当前Top-5召回率93.9%，真正未召回的是水晶兰苷和细叶远志皂苷；牛膝皂苷I和β-细辛醚已分别在Top-3和Top-2。"), P("3. 细叶远志皂苷在正离子模式已采到[M+NH4]+目标峰，但系统仅返回同分子式异构体；因此瓶颈是异构体排序与参考谱域偏移，而不是没有分析正离子。"), P("4. 生成式模型Exact@1为15.0%、Tanimoto@1为0.372、Valid为100%，但质量误差和Oracle Recall@10仍限制真实应用；mix33生成验证尚待完成。"), P("5. 最优落地方案是检索式生产主线与生成式候选补充的混合架构，并通过预测谱、RT、分类和标准品做逐级验证。")]
story += [Spacer(1,7*mm), P("主要参考文献与数据来源", "H2CN")]
refs=[
    "1. Li et al. Spectral entropy for mass spectra similarity. Nature Methods, 2021.",
    "2. Huber et al. MS2Query. Nature Communications, 2023.",
    "3. Litsa et al. MassSpecGym. NeurIPS, 2024.",
    "4. DreaMS. Nature Biotechnology, 2025.",
    "5. Denoising Search. Nature Methods, 2025.",
    "6. DiffMS: molecular structure elucidation with graph diffusion. ICML, 2025.",
    "7. MSGo. Nature Machine Intelligence, 2025.",
    "8. MARLIN. arXiv:2607.04774, 2026 preprint.",
    "9. MassBank open and FAIR mass spectral data resource, 2025/2026 release documentation.",
    "10. 本项目：v5output_33std_beta_massfix_strict.json、benchmark_id_gen_eval_zmem2.json、训练日志与生产任务结果。",
]
for r in refs: story.append(P(r,"SmallCN"))

doc.build(story,onFirstPage=header_footer,onLaterPages=header_footer)
print(OUT)
