"""Build a 細明體 (PMingLiU) DOCX diagnostic report from sample data.

Layout follows docs/reference (報告模板20260506.docx) — five top-level
sections (一/二/三/四/五) with sub-blocks for SNV/Indel, Mitochondrial,
and CNV variants. Reviewer state (status=1 → Causative, status=2 →
Other; status=C/Candidate is NOT included in the formal report) drives
which variants appear; gene-list mode controls how the "本次檢測基因
包括" section in §五 is rendered.

The font is forced to PMingLiU on every run + on Normal style; Word
will substitute when opening on a machine without 細明體 installed.
"""
from __future__ import annotations

import io
import re
from typing import Iterable

from docx import Document
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Pt

from ..config import TERTIARY_OUTPUT_ROOT  # noqa: F401 (kept for callers)
from . import hpo_ontology, phenotype_scorer, report_store, sample_loader

REPORT_FONT     = "PMingLiU"   # 細明體
TITLE_FONT_SIZE = Pt(13)
BODY_FONT_SIZE  = Pt(11)

# ClinVar release these reports reference. Hardcoded — bump when the
# clinvar.vcf.gz used by annotate_clinvar.py is refreshed. Format
# matches the template: ISO-like "YYYYMMDD" and natural "YYYY 年 MM 月
# DD 日".
CLINVAR_DATE      = "20260510"
CLINVAR_DATE_HUMN = "2026 年 05 月 10 日"


# ── Lookup tables ─────────────────────────────────────────────────

# VEP SO term → glossary key. Severity order is the standard VEP
# severity ranking; we pick the first match among the variant's
# `Consequence` field (which VEP joins with "&" or ",").
_CONSEQUENCE_SEVERITY: list[tuple[str, str]] = [
    ("stop_gained",                       "nonsense"),
    ("frameshift_variant",                "frameshift"),
    ("splice_acceptor_variant",           "splice"),
    ("splice_donor_variant",              "splice"),
    ("splice_region_variant",             "splice"),
    ("start_lost",                        "start_loss"),
    ("stop_lost",                         "stop_loss"),
    ("inframe_insertion",                 "in_frame"),
    ("inframe_deletion",                  "in_frame"),
    ("missense_variant",                  "missense"),
    ("synonymous_variant",                "synonymous"),
    ("5_prime_UTR_variant",               "noncoding"),
    ("3_prime_UTR_variant",               "noncoding"),
    ("intron_variant",                    "noncoding"),
    ("non_coding_transcript_exon_variant","noncoding"),
    ("non_coding_transcript_variant",     "noncoding"),
    ("upstream_gene_variant",             "noncoding"),
    ("downstream_gene_variant",           "noncoding"),
    ("regulatory_region_variant",         "noncoding"),
    ("intergenic_variant",                "noncoding"),
]

# Glossary key → 中文制式句 (from 報告模板20260506.docx lines 81-89)
_CONSEQUENCE_ZH: dict[str, str] = {
    "nonsense":   "無義突變 (Nonsense mutation)，此變異會形成過早的終止密碼子，使蛋白質轉譯提前終止，通常導致蛋白質功能喪失。",
    "missense":   "誤義突變 (Missense mutation)，此變異使密碼子改變，造成氨基酸改變，可能影響蛋白質的結構穩定性、功能或相互作用。",
    "synonymous": "同義突變 (Synonymous mutation)，此變異雖不改變蛋白質胺基酸序列，但可能影響 mRNA 剪接效率、穩定性或轉譯效率，因此仍可能具有功能影響。",
    "frameshift": "移碼突變 (Frameshift mutation)，此變異會改變原有閱讀框架，造成後續胺基酸順序大幅改變，通常導致不完整或功能喪失的蛋白質。",
    "in_frame":   "非移碼突變 (In-frame mutation)，此變異不改變閱讀框架，但會增加或刪除完整密碼子，可能影響蛋白質結構、穩定性或功能性區域。",
    "splice":     "剪接位點突變 (Splice site mutation)，此變異可能破壞正常 RNA 剪接，造成外顯子缺失、內含子保留或異常剪接，進而改變蛋白質長度或功能。",
    "start_loss": "起始密碼子缺失 (Start-loss mutation)，此變異可能影響正常翻譯啟動，使蛋白質無法正常生成，或自替代啟動位置開始產生異常蛋白。",
    "stop_loss":  "終止密碼子缺失 (Stop-loss mutation)，此變異會破壞終止密碼子，使轉譯延伸至下游序列生成額外胺基酸，可能改變蛋白質長度與功能。",
    "noncoding":  "非編碼區變異 (Noncoding region mutation)，此變異不直接改變蛋白序列，但可能影響基因調控、剪接、mRNA 穩定性或表達量。",
}

# OMIM inheritance code → 中文
_INHERITANCE_ZH: dict[str, str] = {
    "AD":   "體染色體顯性遺傳",
    "AR":   "體染色體隱性遺傳",
    "XL":   "X 染色體連鎖遺傳",
    "XLD":  "X 染色體連鎖顯性遺傳",
    "XLR":  "X 染色體連鎖隱性遺傳",
    "YL":   "Y 染色體連鎖遺傳",
    "Mi":   "粒線體遺傳",
    "DR":   "雙基因隱性遺傳",
    "DD":   "雙基因顯性遺傳",
    "Smu":  "體細胞突變",
    "Mu":   "多因子遺傳",
    "Isol": "散發性",
}

_ZYGOSITY_ZH: dict[str, str] = {
    "Heterozygous":    "異合子",
    "Homozygous":      "同合子",
    "Hemizygous":      "半合子",
    "het":             "異合子",
    "hom":             "同合子",
    "Heteroplasmic":   "異質性",
    "Homoplasmic":     "同質性",
}

_CNV_KIND_ZH: dict[str, str] = {
    "DEL": "缺失",
    "DUP": "重複",
    "INS": "插入",
    "INV": "倒轉",
}


def _consequence_zh(vep_csq: str) -> str:
    """Pick the most severe glossary entry for a VEP Consequence string.

    The Consequence field is `&`- or `,`-joined VEP SO terms. We walk
    `_CONSEQUENCE_SEVERITY` (already in severity order) and return the
    first 中文 句 whose term appears in the input.
    """
    if not vep_csq:
        return _CONSEQUENCE_ZH["noncoding"]
    tokens = set(re.split(r"[&,;|/\s]+", vep_csq))
    for so_term, key in _CONSEQUENCE_SEVERITY:
        if so_term in tokens:
            return _CONSEQUENCE_ZH[key]
    return _CONSEQUENCE_ZH["noncoding"]


def _inheritance_zh(code: str) -> str:
    """OMIM inheritance code → 中文. Unknown codes pass through verbatim
    so the reviewer can spot weird values (rather than silently mapping
    to a wrong term).
    """
    if not code:
        return ""
    # OMIM inheritance is sometimes a multi-code "AD; AR" — split and
    # translate each, join with 「、」.
    out = []
    for part in re.split(r"[;,/\s]+", code.strip()):
        out.append(_INHERITANCE_ZH.get(part, part))
    return "、".join(p for p in out if p)


def _zygosity_zh(z: str) -> str:
    return _ZYGOSITY_ZH.get((z or "").strip(), z or "")


# ── Font + paragraph helpers ──────────────────────────────────────

def _set_run_font(run, name: str = REPORT_FONT) -> None:
    """Set both Western and East-Asian fonts on a run."""
    run.font.name = name
    rPr = run._element.get_or_add_rPr()
    rFonts = rPr.find(qn("w:rFonts"))
    if rFonts is None:
        rFonts = OxmlElement("w:rFonts")
        rPr.append(rFonts)
    rFonts.set(qn("w:eastAsia"), name)
    rFonts.set(qn("w:ascii"),    name)
    rFonts.set(qn("w:hAnsi"),    name)


def _add_paragraph(doc, text: str, bold: bool = False,
                   size: Pt = BODY_FONT_SIZE, align: str | None = None):
    p = doc.add_paragraph()
    if align == "center":
        p.alignment = 1
    run = p.add_run(text)
    run.bold = bold
    run.font.size = size
    _set_run_font(run)
    return p


def _blank(doc):
    p = doc.add_paragraph()
    run = p.add_run("")
    run.font.size = BODY_FONT_SIZE
    _set_run_font(run)


def _apply_normal_font(doc) -> None:
    normal = doc.styles["Normal"]
    normal.font.name = REPORT_FONT
    normal.font.size = BODY_FONT_SIZE
    rPr = normal.element.get_or_add_rPr()
    rFonts = rPr.find(qn("w:rFonts"))
    if rFonts is None:
        rFonts = OxmlElement("w:rFonts")
        rPr.append(rFonts)
    rFonts.set(qn("w:eastAsia"), REPORT_FONT)
    rFonts.set(qn("w:ascii"),    REPORT_FONT)
    rFonts.set(qn("w:hAnsi"),    REPORT_FONT)


# ── Variant-bucketing helpers ─────────────────────────────────────

def _buckets_for_type(variants: dict, statuses: dict, status: str) -> list[dict]:
    """Pull variants whose reviewer-set status matches; preserve order
    from the variants dict (tier-sorted by the adapter).
    """
    out = []
    for vid in variants:
        s = (statuses.get(vid) or "").strip()
        if s == status:
            out.append(variants[vid])
    return out


def _manual_for(report: dict, status: str) -> list[dict]:
    out = []
    for m in (report.get("manual_variants") or []):
        if (m or {}).get("status") == status and (m.get("position") or "").strip():
            out.append(m)
    return out


# ── Sections ──────────────────────────────────────────────────────

def _section_test_info(doc, test_type: str) -> None:
    """一、檢驗項目"""
    is_wgs = (test_type or "").upper() == "WGS"
    label = "次世代定序全基因組定序檢測-單基因遺傳疾病" if is_wgs \
        else "次世代定序全外顯子定序檢測-單基因遺傳疾病"
    _add_paragraph(doc, f"一、檢驗項目: {label}")
    _blank(doc)


def _section_panel_set(doc) -> None:
    """二、檢驗套組  (固定「非特定」for now — reviewer 套組概念之後再加)"""
    _add_paragraph(doc, "二、檢驗套組: 非特定")
    _blank(doc)


def _summary_line(name: str, hpo_labels: list[str]) -> str:
    """『在非特定 (term, term, ...) 檢驗套組中未找到已知致病性位點。』"""
    if hpo_labels:
        return f"    在{name} ({', '.join(hpo_labels)}) 檢驗套組中未找到已知致病性位點。"
    return f"    在{name}檢驗套組中未找到已知致病性位點。"


def _section_results(doc, sample: dict, report: dict, test_type: str) -> None:
    """三、檢測結果 — the meat of the report."""
    _add_paragraph(doc, "三、檢測結果")
    _add_paragraph(doc, "  檢體說明:")
    _add_paragraph(doc, "    檢體類別：血液")

    # — 綜合說明
    _add_paragraph(doc, "  綜合說明:")
    statuses = report.get("status", {}) or {}
    variants = sample.get("variants", {}) or {}
    cnv_vars = (sample.get("cnv") or {}).get("variants", {}) or {}
    sv_vars  = (sample.get("sv")  or {}).get("variants", {}) or {}
    mito_vars = (sample.get("mito") or {}).get("variants", {}) or {}

    snv_t1  = _buckets_for_type(variants,  statuses, "1")
    snv_t2  = _buckets_for_type(variants,  statuses, "2")
    cnv_t1  = _buckets_for_type(cnv_vars,  statuses, "1")
    cnv_t2  = _buckets_for_type(cnv_vars,  statuses, "2")
    sv_t1   = _buckets_for_type(sv_vars,   statuses, "1")
    sv_t2   = _buckets_for_type(sv_vars,   statuses, "2")
    mito_t1 = _buckets_for_type(mito_vars, statuses, "1")
    mito_t2 = _buckets_for_type(mito_vars, statuses, "2")

    has_t1 = bool(snv_t1 or cnv_t1 or sv_t1 or mito_t1
                  or _manual_for(report, "1"))
    has_t2 = bool(snv_t2 or cnv_t2 or sv_t2 or mito_t2
                  or _manual_for(report, "2"))

    _add_paragraph(doc, "    第一類：與臨床症狀相關基因之已知致病性變異位點")
    if not has_t1:
        # Mirror the template wording when the bucket is empty.
        hpo_labels = [r.get("label", "") or r.get("phenotype", "")
                      for r in (sample.get("patient_phenotype") or [])]
        hpo_labels = [x for x in hpo_labels if x][:5]
        _add_paragraph(doc, _summary_line("非特定", hpo_labels))

    _blank(doc)
    _add_paragraph(doc, "    第二類：其他變異位點")
    if not has_t2:
        _add_paragraph(doc, "    未找到其他變異位點。")
    _blank(doc)

    # — Variant-type sections (only render headings/contents when present)
    if snv_t1 or snv_t2 or _manual_for(report, "1") or _manual_for(report, "2"):
        _subsection_snv(doc, snv_t1, snv_t2, report)
    if mito_t1 or mito_t2:
        _subsection_mito(doc, mito_t1, mito_t2, report)
    # CNV + SV share the same template (染色體 / 變異位置 / 拷貝數 / ACMG).
    if cnv_t1 or cnv_t2 or sv_t1 or sv_t2:
        _subsection_cnv(doc, cnv_t1 + sv_t1, cnv_t2 + sv_t2, test_type, report)

    # — Footer recommendation (always shown)
    _add_paragraph(doc, "    建議比對臨床表徵並進行父母親與家族成員之變異位點檢測，"
                        "以釐清上述變異致病之可能性；根據家族成員變異位點檢測報告或"
                        "相關資料庫更新，可能影響變異位點ACMG判讀結果。")
    _blank(doc)


# ── Subsections: per variant type ─────────────────────────────────

def _afs_str(variant: dict) -> str:
    """Return e.g. '0.01%' or '未報導過發生率'. Prefer the most-relevant
    AF column (gnomAD genome > exome).
    """
    for key in ("AF", "AF_exome"):
        v = variant.get(key)
        if v in (None, "", "."):
            continue
        try:
            f = float(v)
            return f"{f * 100:.4f}%".rstrip("0").rstrip(".")
        except ValueError:
            continue
    return ""


def _clinvar_label(variant: dict) -> str:
    sig = (variant.get("CLNSIG") or "").strip()
    return sig.replace("_", " ") if sig else "—"


def _acmg_label(variant: dict, edits_for_v: dict) -> str:
    cls = edits_for_v.get("ACMG_classification") \
        or variant.get("ACMG_classification") or ""
    return (cls or "—").strip()


def _omim_block_for_snv(v: dict) -> str:
    """『GENE 為 DISEASE 的致病基因之一，其遺傳模式屬於 X
    (Phenotype MIM number: M)』 — fall back to whatever is present."""
    gene = v.get("gene_symbol") or v.get("GENE") or "?"
    disease = (v.get("Disease1") or v.get("OMIM_disease") or "").splitlines()
    disease_str = disease[0].strip() if disease else ""
    inh_zh = _inheritance_zh(v.get("Inheritance", "") or "")
    mim = (v.get("OMIM_id") or "").strip()

    parts = [f"{gene}為"]
    if disease_str:
        # Disease strings sometimes carry a trailing inheritance tag in
        # parens — strip "(AD)" etc. so we don't double-print.
        disease_str = re.sub(r"\s*\([A-Za-z?;,/ ]+\)\s*$", "", disease_str)
        parts.append(f"{disease_str}的致病基因之一")
    else:
        parts.append("此疾病的致病基因之一")
    if inh_zh:
        parts.append(f"，其遺傳模式屬於{inh_zh}")
    if mim:
        parts.append(f" (Phenotype MIM number: {mim})")
    return "".join(parts) + "。"


def _subsection_snv(doc, t1: list[dict], t2: list[dict], report: dict) -> None:
    _add_paragraph(doc, "## SNV/indel", bold=True)
    edits = report.get("edits", {}) or {}
    for v in t1:
        _snv_variant_block(doc, v, tier="1", edits=edits.get(v["id"], {}))
    for v in t2:
        _snv_variant_block(doc, v, tier="2", edits=edits.get(v["id"], {}))
    # manual SNV entries (reviewer-typed "＋ 新增 variant")
    for m in _manual_for(report, "1") + _manual_for(report, "2"):
        _add_paragraph(doc, f"    {m.get('position', '')}", bold=True)
        if m.get("disease"):
            _add_paragraph(doc, f"    {m['disease']}")
        if m.get("comment"):
            _add_paragraph(doc, f"    {m['comment']}")
        _blank(doc)


def _snv_variant_block(doc, v: dict, *, tier: str, edits: dict) -> None:
    gene = v.get("gene_symbol") or "?"
    tx   = v.get("transcript")  or v.get("MANE_SELECT") or ""
    _add_paragraph(doc, f"    {gene} ({tx})", bold=True)
    # ASCII-table-style block matching the template
    rs    = ""   # pipeline doesn't carry rsID yet — leave blank
    struc = (v.get("exon") or v.get("intron") or "").strip()
    hgvs_c = v.get("hgvs_c") or v.get("HGVS_C") or ""
    hgvs_p = v.get("hgvs_p") or v.get("HGVS_P") or ""
    nuc    = f"{hgvs_c}" + (f"({hgvs_p})" if hgvs_p else "")
    zyg    = _zygosity_zh(v.get("zygosity", ""))
    clnsig = _clinvar_label(v)
    acmg   = _acmg_label(v, edits)

    sep = "=" * 78
    _add_paragraph(doc, f"    {sep}")
    _add_paragraph(doc, f"     類別  基因          RS ID  結構     核苷酸                       基因型     ClinVar              ACMG&AMP指引")
    _add_paragraph(doc, f"    {sep}")
    _add_paragraph(doc, f"     {tier:<5}{gene:<14}{rs:<7}{struc:<9}{nuc:<29}{zyg:<11}{clnsig:<21}{acmg}")
    _add_paragraph(doc, f"    {sep}")
    # Numbered notes
    _add_paragraph(doc, f"    1. {_omim_block_for_snv(v)}")
    if tier == "1":
        _add_paragraph(doc, "    2. 此為致病性之變異位點，與臨床症狀相關。")
    else:
        _add_paragraph(doc, "    2. 此變異位點之臨床意義須由醫師配合其他相關資料進行最佳綜合判斷。")
    # Reviewer comment
    cmt = (edits.get("comment") or "").strip()
    if cmt:
        _add_paragraph(doc, f"    3. {cmt}")

    # — 參考資料 (free-text per variant)
    _blank(doc)
    _add_paragraph(doc, "  參考資料:")
    _add_paragraph(doc, _snv_reference_text(v, edits))
    _blank(doc)


def _snv_reference_text(v: dict, edits: dict) -> str:
    gene = v.get("gene_symbol") or "?"
    hgvs_c = v.get("hgvs_c") or v.get("HGVS_C") or ""
    hgvs_p = v.get("hgvs_p") or v.get("HGVS_P") or ""
    nuc    = hgvs_c + (f" ({hgvs_p})" if hgvs_p else "")
    cq_zh  = _consequence_zh(v.get("Consequence", ""))
    af     = _afs_str(v)
    clinvar = _clinvar_label(v)
    acmg    = _acmg_label(v, edits)

    af_sent = (f"該變異位點在族群資料庫 gnomAD 中發生率為 {af}，"
               "顯示其為罕見變異位點") if af \
        else "該變異位點在族群資料庫 gnomAD 中未報導過發生率，顯示其為罕見變異位點"
    clnv_sent = (f"在疾病資料庫 (ClinVar) 中此變異位點被報導為「{clinvar}」"
                 if clinvar and clinvar != "—" else
                 "在疾病資料庫 (ClinVar) 中尚未有此變異位點之紀錄")

    # PubMed sentence intentionally omitted — the new pipeline carries
    # no HGMD reported-cases column.
    return (
        f"    在個案之檢體中，檢測到1個位於基因{gene}的變異位點。"
        f"變異位點 {nuc} 為{cq_zh}"
        f"{af_sent}。{clnv_sent}。"
        "根據美國醫學遺傳學暨基因體學學會 (American College of Medical "
        "Genetics and Genomics) 與分子病理學學會 (Association for "
        f"Molecular Pathology) 於2015年發表之準則，評測此變異位點為「{acmg}」。"
        "此報告僅供參考，臨床判斷仍應以病患的實際狀況為主。"
    )


# ── Mitochondrial ─────────────────────────────────────────────────

def _subsection_mito(doc, t1: list[dict], t2: list[dict], report: dict) -> None:
    _add_paragraph(doc, "## Mitochondrial variants", bold=True)
    edits = report.get("edits", {}) or {}
    for v in t1:
        _mito_variant_block(doc, v, tier="1", edits=edits.get(v["id"], {}))
    for v in t2:
        _mito_variant_block(doc, v, tier="2", edits=edits.get(v["id"], {}))


def _mito_variant_block(doc, v: dict, *, tier: str, edits: dict) -> None:
    gene  = v.get("gene_symbol") or "?"
    _add_paragraph(doc, f"    {gene} (NC_012920.1)", bold=True)
    hgvs_m = v.get("HGVS_M") or v.get("id") or ""
    aa     = (v.get("aa_change") or "").strip()
    nuc    = hgvs_m + (f"({aa})" if aa else "")
    rs     = ""
    het    = v.get("heteroplasmy")
    het_s  = f"{float(het) * 100:.1f}%" if het not in (None, "", ".") else "—"
    clnsig = _clinvar_label(v)
    acmg   = _acmg_label(v, edits)

    sep = "=" * 78
    _add_paragraph(doc, f"    {sep}")
    _add_paragraph(doc, f"     類別  基因         RS ID  核苷酸                       異質性比例   ClinVar              ACMG&AMP指引")
    _add_paragraph(doc, f"    {sep}")
    _add_paragraph(doc, f"     {tier:<5}{gene:<13}{rs:<7}{nuc:<29}{het_s:<13}{clnsig:<21}{acmg}")
    _add_paragraph(doc, f"    {sep}")

    mt_disease = (v.get("mitomap_disease") or "").strip()
    # 1. 致病基因之一 — MITOMAP disease name; 遺傳模式 fixed to 粒線體遺傳; no MIM
    if mt_disease:
        _add_paragraph(doc, f"    1. {gene}為{mt_disease}的致病基因之一，其遺傳模式屬於粒線體遺傳。")
    else:
        _add_paragraph(doc, f"    1. {gene}為粒線體疾病之相關基因，其遺傳模式屬於粒線體遺傳。")
    if tier == "1":
        _add_paragraph(doc, "    2. 此為致病性之變異位點，與臨床症狀相關。")
    else:
        _add_paragraph(doc, "    2. 此變異位點臨床意義須由醫師配合其他相關資料進行最佳綜合判斷。")
    _add_paragraph(doc, "    3. 此檢驗結果須由臨床醫師來分析粒線體DNA基因檢驗結果與受檢者臨床症狀的相關性"
                        "並考量家族史的關聯性。此外，粒線體DNA基因變異有組織的特異性，粒線體DNA基因致病變異的"
                        "異質性（heteroplasmy）在不同組織間會有差異。")
    _add_paragraph(doc, "    4. 對於所有粒線體DNA基因變異所造成的臨床症狀表現的差異，可能取決於下列三項因素："
                        "（1）異質性，即致病變異及正常粒線體 DNA 存在量的比例（2）變異的粒線體 DNA 於組織的"
                        "分布情形（3）限界效應（Threshold effect），即每種組織對於氧化壓力代謝影響的易受性。"
                        "由於患者情況各異，此檢驗結果須由臨床醫師判讀檢驗結果與受檢者臨床症狀的相關性。")
    cmt = (edits.get("comment") or "").strip()
    if cmt:
        _add_paragraph(doc, f"    5. {cmt}")

    _blank(doc)
    _add_paragraph(doc, "  參考資料:")
    _add_paragraph(doc, _mito_reference_text(v, edits))
    _blank(doc)


def _mito_reference_text(v: dict, edits: dict) -> str:
    gene = v.get("gene_symbol") or "?"
    hgvs_m = v.get("HGVS_M") or v.get("id") or ""
    aa = (v.get("aa_change") or "").strip()
    nuc = hgvs_m + (f"({aa})" if aa else "")
    cq_zh = _consequence_zh(v.get("consequence", ""))
    clinvar = _clinvar_label(v)
    acmg    = _acmg_label(v, edits)
    clnv_sent = (f"在疾病資料庫 (ClinVar) 中此變異位點被報導為「{clinvar}」"
                 if clinvar and clinvar != "—"
                 else "在疾病資料庫 (ClinVar) 中尚未有此變異位點之紀錄")
    return (
        f"    在個案之檢體中，檢測到1個位於基因{gene}的變異位點。"
        f"變異位點 {nuc} 為{cq_zh}"
        f"{clnv_sent}。"
        "根據美國醫學遺傳學暨基因體學學會 (American College of Medical "
        "Genetics and Genomics) 與分子病理學學會 (Association for "
        "Molecular Pathology) 於2015年發表之準則，並參照ClinGen Mitochondrial "
        "Disease Expert Panel於 2020 年發布之粒線體DNA變異判讀專用ACMG/AMP標準"
        f"進行評估，評測此變異位點為「{acmg}」。"
        "此報告僅供參考，臨床判斷仍應以病患的實際狀況為主。"
    )


# ── CNV / SV ──────────────────────────────────────────────────────

def _subsection_cnv(doc, t1: list[dict], t2: list[dict], test_type: str,
                    report: dict) -> None:
    _add_paragraph(doc, "## CNV", bold=True)
    edits = report.get("edits", {}) or {}
    is_wgs = (test_type or "").upper() == "WGS"
    for v in t1:
        _cnv_variant_block(doc, v, tier="1", is_wgs=is_wgs,
                           edits=edits.get(v["id"], {}))
    for v in t2:
        _cnv_variant_block(doc, v, tier="2", is_wgs=is_wgs,
                           edits=edits.get(v["id"], {}))


def _coords_str(v: dict, is_wgs: bool) -> str:
    """WES uses imprecise breakpoint notation `(?_start)_(end_?)`;
    WGS uses the precise start_end form. The adapter's `coords` is the
    raw AnnotSV CHROM:POS-END so we annotate accordingly.
    """
    coords = v.get("coords") or ""
    if is_wgs or not coords:
        return coords
    # WES form
    m = re.match(r"(chr[\dXYM]+):(\d+)-(\d+)", coords)
    if not m:
        return coords
    ch, s, e = m.group(1), m.group(2), m.group(3)
    return f"{ch}:(?_{s})_({e}_?)"


def _sv_kind_zh(v: dict) -> str:
    return _CNV_KIND_ZH.get((v.get("sv_type") or "").upper(), "變異")


def _cnv_variant_block(doc, v: dict, *, tier: str, is_wgs: bool,
                       edits: dict) -> None:
    coords    = _coords_str(v, is_wgs)
    sv_type   = (v.get("sv_type") or "").lower()  # del/dup/...
    cn        = v.get("copy_number")
    cn_s      = "—" if cn in (None, "", ".") else str(cn)
    zyg       = _zygosity_zh(v.get("zygosity", ""))
    acmg      = _acmg_label(v, edits)
    chrom     = (v.get("CHROM") or "").replace("chr", "")
    coords_disp = coords.split(":", 1)[1] if ":" in coords else coords

    _add_paragraph(doc, f"    [GRCh38] {coords}{'' if is_wgs else 'del/dup'}",
                   bold=True)
    sep = "=" * 78
    _add_paragraph(doc, f"    {sep}")
    _add_paragraph(doc, f"     類別  染色體   變異位置                                  拷貝數   基因型     ACMG&AMP指引")
    _add_paragraph(doc, f"    {sep}")
    _add_paragraph(doc, f"     {tier:<5}{chrom:<9}{coords_disp:<42}{cn_s:<9}{zyg:<11}{acmg}")
    _add_paragraph(doc, f"    {sep}")

    # 1. 片段位置描述 — single-gene vs multi-gene region
    gene_count = int(v.get("gene_count") or 0)
    genes      = v.get("genes") or []   # may be top-10 only
    chrom_num  = (v.get("CHROM") or "").replace("chr", "")

    # Reviewer can prune which genes go on the report via edits.report_genes
    report_genes = edits.get("report_genes") or {}
    if isinstance(report_genes, dict):
        kept = [g for g in genes if report_genes.get(g.get("gene") if isinstance(g, dict) else g)]
        if kept:
            genes = kept

    if gene_count <= 1 and genes:
        g = genes[0]
        gname = g["gene"] if isinstance(g, dict) else g
        # AnnotSV's location string is e.g. "exon1-exon8" or "intron3-exon5"
        loc = (v.get("location") or "").strip() or "整個區域"
        _add_paragraph(doc, f"    1. 此片段位於第 {chrom_num} 號染色體上 {gname} 基因 {loc}。")
        # 2. OMIM + inheritance
        ph  = (v.get("omim_phenotype") or "").strip()
        inh = _inheritance_zh(v.get("omim_inheritance", "") or "")
        mim = (v.get("omim_id") or "").strip()
        bits = [f"{gname}為"]
        bits.append(f"{ph}的致病基因之一" if ph else "此疾病的致病基因之一")
        if inh: bits.append(f"，其遺傳模式屬於{inh}")
        if mim: bits.append(f" (Phenotype MIM number: {mim})")
        _add_paragraph(doc, f"    2. {''.join(bits)}。")
    else:
        # Multi-gene region template
        names = [g["gene"] if isinstance(g, dict) else g for g in genes[:10]]
        more  = "等" if gene_count > len(names) else ""
        gene_list_str = ", ".join(names)
        _add_paragraph(doc, f"    1. 此片段位於第 {chrom_num} 號染色體上，"
                            f"包含 {gene_list_str}{more} OMIM 疾病基因。")
        # 2. 區域對應疾病 — left to reviewer comment (no clean source).
        cmt = (edits.get("comment") or "").strip()
        if cmt:
            _add_paragraph(doc, f"    2. {cmt}")

    note_idx = 3
    if tier == "1":
        _add_paragraph(doc, f"    {note_idx}. 此為致病性之變異位點，與臨床症狀相關。")
    else:
        _add_paragraph(doc, f"    {note_idx}. 此變異位點之臨床意義須由醫師配合其他相關資料進行最佳綜合判斷。")
    note_idx += 1
    if not is_wgs:
        _add_paragraph(doc, f"    {note_idx}. 由於此檢驗技術為全外顯子定序，"
                            "若缺失片段之斷點(Breakpoints)發生於內含子(Intron) ，"
                            "則無法明確判別起始及末端位置。")

    _blank(doc)
    _add_paragraph(doc, "  參考資料:")
    _add_paragraph(doc, _cnv_reference_text(v, edits, gene_count, genes))
    _blank(doc)


def _cnv_reference_text(v: dict, edits: dict, gene_count: int,
                        genes: Iterable[dict]) -> str:
    coords = v.get("coords") or ""
    sv_type = (v.get("sv_type") or "").upper()
    kind_zh = _CNV_KIND_ZH.get(sv_type, "變異")
    zyg     = _zygosity_zh(v.get("zygosity", "")) or "異合子"
    acmg    = _acmg_label(v, edits)

    if gene_count <= 1 and genes:
        g = next(iter(genes))
        gname = g["gene"] if isinstance(g, dict) else g
        loc = (v.get("location") or "").strip() or "整個區域"
        span_desc = f"此段{kind_zh}涵蓋 {gname} 基因 {loc}"
    else:
        names = [g["gene"] if isinstance(g, dict) else g for g in list(genes)[:10]]
        more = "等" if gene_count > len(names) else ""
        span_desc = f"此段{kind_zh}涵蓋 {', '.join(names)}{more} OMIM 疾病基因"

    return (
        f"    在個案之檢體中，檢測到位於 {coords} 之{zyg}片段{kind_zh}變異，"
        f"{span_desc}。"
        "根據美國醫學遺傳學暨基因體學學會 (American College of Medical "
        "Genetics and Genomics) 與分子病理學學會 (Association for Molecular "
        "Pathology) 於2015年發表之準則，並參照ClinGen 及Riggs等人於 2020 年發布之"
        f"拷貝數變異判讀專用ACMG/ClinGen技術標準進行評估，評測此變異位點為「{acmg}」。"
        "此報告僅供參考，臨床判斷仍應以病患的實際狀況為主。"
    )


# ── §四 方法、§五 注釋 ───────────────────────────────────────────

def _section_methods(doc, test_type: str) -> None:
    is_wgs = (test_type or "").upper() == "WGS"
    seq    = "Illumina NovaSeq X Plus" if is_wgs else "Illumina NextSeq 2000"
    depth  = "30X" if is_wgs else "50X"   # WGS standard ~30X mean coverage
    _add_paragraph(doc, "四、檢測方法說明")
    _add_paragraph(doc, f"  1. 本次檢測使用次世代定序儀分析 ({seq})。")
    _add_paragraph(doc, "  2. 本次檢測變異位點的錯誤率 ≦ 0.1% (Phred-scaled Q score ≧ 30)。")
    _add_paragraph(doc, f"  3. 本次檢測平均定序深度 ≧ {depth}。")
    _add_paragraph(doc, "  4. 本檢測僅能檢測出基因內單一核苷酸變異 (single nucleotide variant) 、"
                        "小片段的缺失或插入 (small indel)及部分拷貝數變異 (copy number variant)，"
                        "無法檢測出轉位 (translocation)、倒轉 (inversion) 或其他複雜性結構變異 "
                        "(complex structural variant)、組織特異性的鑲嵌 (tissue-specific mosaicism) "
                        "以及未包含在本次定序範圍之區域。")
    _add_paragraph(doc, "  5. 本實驗方法以次世代方法定序粒線體DNA基因序列，"
                        "變異點位判讀之cut-off值定為5%異質性（heteroplasmy）。")
    _add_paragraph(doc, "  6. 本檢測報告僅供醫療專業人員參考，需配合其他相關臨床資料與家族成員之相關檢驗。"
                        "依衛福部規定，目前次世代定序分子遺傳診斷皆屬研究性質。")
    _blank(doc)


def _section_annotations(doc, sample: dict, gene_list_mode: str) -> None:
    _add_paragraph(doc, "五、檢測結果注釋")
    _add_paragraph(doc, "  1. 本檢測結果比對參考序列為人類hg38版本。")
    _add_paragraph(doc, f"  2. ClinVar及ACMG&AMP指引：引用ClinVar資料庫截至{CLINVAR_DATE_HUMN}更新的註解，"
                        "及美國醫學遺傳學暨基因體學學會 (ACMG) 與分子病理學學會 (AMP) 2015年頒佈的指引，"
                        "並且主要列入致病(Pathogenic) 及可能致病 (Likely pathogenic) 變異；"
                        "其他類別變異經醫師判斷認為與疾病相關時亦可列入。")
    _add_paragraph(doc, "  3. 參考資料:")
    _add_paragraph(doc, f"     a. 疾病資料庫: OMIM、ClinVar ({CLINVAR_DATE})")
    _add_paragraph(doc, "     b. 族群資料庫: gnomAD (v4.1 genome)")
    _add_paragraph(doc, "     c. 序列資料庫: RefSeqGene (105.20220307)")
    _add_paragraph(doc, "  4. 本次檢測基因包括")
    _render_gene_list(doc, sample, gene_list_mode)


def _hpo_label_for(hid: str) -> str:
    term = hpo_ontology.get(hid)
    return term.name if term else hid


def _genes_for_term_or_panel(key: str) -> list[str]:
    """phenotype_scorer._HPO_TO_GENES is the union of HPO→gene map and
    every panel's gene set (HPO ids + panel names live in the same dict)."""
    return sorted(phenotype_scorer._HPO_TO_GENES.get(key, set()))


def _render_gene_list(doc, sample: dict, mode: str) -> None:
    """mode = 'grouped' → one paragraph per HPO term / panel,
       mode = 'merged'  → single deduped list.
    """
    hpo_rows: list = sample.get("patient_phenotype") or []
    panel_names: list[str] = sample.get("selected_panels") or []

    # Build [(display_name, [genes...])] preserving order.
    sections: list[tuple[str, list[str]]] = []
    for r in hpo_rows:
        hid = r.get("phenotype") or ""
        label = r.get("label") or _hpo_label_for(hid)
        if not hid: continue
        genes = _genes_for_term_or_panel(hid)
        sections.append((f"{label} ({hid})", genes))
    for pname in panel_names:
        genes = _genes_for_term_or_panel(pname)
        sections.append((pname, genes))

    if not sections:
        _add_paragraph(doc, "    （未設定 HPO / panel — 無檢測基因清單）")
        return

    if mode == "merged":
        merged: set[str] = set()
        for _, gs in sections:
            merged |= set(gs)
        gene_str = ", ".join(sorted(merged))
        _add_paragraph(doc, f"    {gene_str}")
        return

    # grouped (default)
    for name, gs in sections:
        _add_paragraph(doc, f"    {name}:")
        if gs:
            _add_paragraph(doc, f"    {', '.join(gs)}")
        else:
            _add_paragraph(doc, "    （無對應基因）")


# ── Top-level entrypoint ──────────────────────────────────────────

def build_diagnosis_docx(sample_id: str, *, gene_list_mode: str = "grouped") -> bytes:
    """Render the full diagnostic report as DOCX bytes.

    gene_list_mode = "grouped" (default) or "merged" controls how §五.4
    "本次檢測基因包括" is rendered.
    """
    sample = sample_loader.load_sample(sample_id, include_aux=True)
    if sample is None:
        raise FileNotFoundError(f"sample not found: {sample_id}")

    report = report_store.load(sample_id)
    meta   = sample.get("meta") or {}
    test_type = meta.get("Test", "") or "WES"

    doc = Document()
    _apply_normal_font(doc)

    _section_test_info(doc, test_type)
    _section_panel_set(doc)
    _section_results(doc, sample, report, test_type)
    _section_methods(doc, test_type)
    _section_annotations(doc, sample, gene_list_mode)

    buf = io.BytesIO()
    doc.save(buf)
    return buf.getvalue()
