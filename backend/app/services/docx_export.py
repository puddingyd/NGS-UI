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
import unicodedata
from typing import Iterable

from docx import Document
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Cm, Pt

from ..config import TERTIARY_OUTPUT_ROOT  # noqa: F401 (kept for callers)
from . import hpo_ontology, phenotype_scorer, report_store, sample_loader

# 細明體 (MingLiU) — the *monospace* CJK family. (新細明體 = PMingLiU
# is proportional, which would break the ASCII tables below.) Word
# falls back to whatever the system has if MingLiU isn't installed.
REPORT_FONT     = "MingLiU"
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
    ("missense",                          "missense"),    # mito script short form
    ("synonymous_variant",                "synonymous"),
    ("synonymous",                        "synonymous"),  # mito script short form
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
    "XL":   "性聯遺傳",
    "XLD":  "性染色體顯性遺傳",
    "XLR":  "性染色體隱性遺傳",
    "YL":   "Y 染色體連鎖遺傳",
    "MT":   "粒線體遺傳",
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

# Amino-acid single-letter → three-letter mapping. Used to normalise
# MITOMAP's compact `A64V` style to HGVS `p.Ala64Val` so the report
# matches the SNV nucleotide column format.
_AA_1_TO_3: dict[str, str] = {
    "A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys",
    "E": "Glu", "Q": "Gln", "G": "Gly", "H": "His", "I": "Ile",
    "L": "Leu", "K": "Lys", "M": "Met", "F": "Phe", "P": "Pro",
    "S": "Ser", "T": "Thr", "W": "Trp", "Y": "Tyr", "V": "Val",
    "*": "Ter", "X": "Ter", "U": "Sec", "O": "Pyl",
}

# AnnotSV's ACMG/ClinGen classification is an int 1-5 on full rows.
_CNV_ACMG_INT_TO_LABEL: dict[int, str] = {
    1: "Benign",
    2: "Likely benign",
    3: "Uncertain significance",
    4: "Likely pathogenic",
    5: "Pathogenic",
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


_ZYG_LONG: dict[str, str] = {
    "het":          "Heterozygous",
    "hom":          "Homozygous",
    "hemi":         "Hemizygous",
    "homo":         "Homozygous",
    "heterozygous": "Heterozygous",
    "homozygous":   "Homozygous",
    "hemizygous":   "Hemizygous",
}


def _zygosity_long(z: str) -> str:
    """Expand short forms (het/hom/hemi) to the canonical English long
    form used in the diagnostic report tables."""
    s = (z or "").strip()
    return _ZYG_LONG.get(s.lower(), s)


def _strip_tx_prefix(s: str) -> str:
    """Pipeline HGVS values sometimes come prefixed with the transcript
    ID (e.g. "NM_000295.5:c.1075A>G"); the report's 核苷酸 column shows
    transcript-only on the gene header line, so strip the prefix here.
    """
    return re.sub(r"^[A-Z]+_\d+(\.\d+)?:", "", s or "")


def _mito_aa_to_hgvsp(aa: str) -> str:
    """MITOMAP gives "A64V" — convert to "p.Ala64Val" so the report's
    nucleotide cell matches the SNV HGVS format. Already-3-letter input
    (with or without `p.` prefix) is normalised to carry the `p.`.
    Empty / unparseable input passes through.
    """
    if not aa:
        return ""
    s = aa.strip()
    if s.startswith("p."):
        s = s[2:]
    # 1-letter form: A64V or *64L (stop)
    m1 = re.match(r"^([A-Z*])(\d+)([A-Z*])$", s)
    if m1:
        a1, pos, a2 = m1.groups()
        return f"p.{_AA_1_TO_3.get(a1, a1)}{pos}{_AA_1_TO_3.get(a2, a2)}"
    # 3-letter form: Ala64Val
    m3 = re.match(r"^([A-Z][a-z]{2}|Ter)(\d+)([A-Z][a-z]{2}|Ter)$", s)
    if m3:
        return f"p.{s}"
    return aa


def _structure_label(v: dict) -> str:
    """VEP's EXON / INTRON fields look like "5/22" (rank/total). The
    report writes them as `exon5` or `intron6` — pick whichever has
    a value and prefix with the type."""
    exon = (v.get("exon") or "").strip()
    intron = (v.get("intron") or "").strip()
    if exon:
        return "exon" + exon.split("/")[0].strip()
    if intron:
        return "intron" + intron.split("/")[0].strip()
    return ""


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
    # Word's default Normal style adds 8pt space after each paragraph;
    # the reviewer was manually clearing that on every export. Set it
    # to 0 so the report opens already tight.
    pf = normal.paragraph_format
    pf.space_after  = Pt(0)
    pf.space_before = Pt(0)
    # Tighten line spacing too — single (1.0) instead of 1.15 multiple.
    from docx.shared import Pt as _Pt  # noqa: F401 (silence linter)
    pf.line_spacing = 1.0


def _apply_page_margins(doc) -> None:
    """A4 portrait with 1.5cm L/R margins so the 89-char-wide variant
    tables (4-space indent + 85-char ===== box) fit on one line in
    細明體 11pt."""
    for section in doc.sections:
        section.left_margin  = Cm(1.5)
        section.right_margin = Cm(1.5)
        section.top_margin    = Cm(2.0)
        section.bottom_margin = Cm(2.0)


# ── ASCII-table helpers ───────────────────────────────────────────
# 細明體 renders CJK chars at 2x ASCII width and ASCII chars at uniform
# half-width, so a fixed-width column layout drawn with spaces lines up
# visually. We treat East-Asian Wide/Fullwidth (W/F) as 2 columns and
# everything else as 1.

def _ea_width(c: str) -> int:
    return 2 if unicodedata.east_asian_width(c) in ("W", "F") else 1


def _str_width(s: str) -> int:
    return sum(_ea_width(c) for c in s)


def _pad_right(s: str, w: int) -> str:
    """Right-pad with spaces up to display width w."""
    return s + " " * max(0, w - _str_width(s))


def _wrap_to_cols(text: str, width: int, mode: str = "char") -> list[str]:
    """Wrap `text` so each chunk fits within `width` display columns.

    mode="char" → greedy char-by-char wrap respecting CJK double width
                  (default; mid-token break, matches template's HGVS
                  wrap behaviour).
    mode="token" → one whitespace-delimited token per line; also splits
                   after "/" so "Pathogenic/Likely pathogenic" becomes
                   3 lines (Pathogenic/ · Likely · pathogenic).
                   Tokens longer than `width` fall back to char wrap.
    """
    if text is None:
        return [""]
    s = str(text)
    if not s:
        return [""]

    if mode == "token":
        # Split on whitespace; also break after "/" so combined sigs
        # like "Pathogenic/Likely pathogenic" land on three lines.
        tokens = [t for t in re.split(r"\s+|(?<=/)", s) if t]
        out: list[str] = []
        for tok in tokens:
            if _str_width(tok) <= width:
                out.append(tok)
            else:
                out.extend(_wrap_to_cols(tok, width, mode="char"))
        return out or [""]

    if mode == "hgvs":
        # Keep one blank display column between HGVS content and the next
        # table cell. Protein notation always starts on its own line.
        content_width = max(1, width - 1)
        m = re.search(r"(\(p\.[^)]*\)?)$", s)
        if m and m.start() > 0:
            prefix = s[: m.start()]
            suffix = m.group(1)
            return (
                _wrap_to_cols(prefix, content_width, mode="char")
                + _wrap_to_cols(suffix, content_width, mode="char")
            )
        return _wrap_to_cols(s, content_width, mode="char")

    out, cur, cur_w = [], "", 0
    for c in s:
        cw = _ea_width(c)
        if cur and cur_w + cw > width:
            out.append(cur)
            cur, cur_w = c, cw
        else:
            cur += c
            cur_w += cw
    if cur:
        out.append(cur)
    return out


def _ascii_table(doc,
                 columns: list[tuple],  # (header, width) or (header, width, mode)
                 rows: list[list[str]],
                 indent: str = "    ") -> None:
    """Render a ==== bounded table where each column has a fixed display
    width; over-long cells wrap per the column's `mode` ("char" default
    / "token" for clinical sigs). Columns are packed without inter-
    column whitespace — pad each column to its full width on its own;
    the data row gets a single space of side-padding inside the ====
    boundary.

    With column widths summing to N and (k) columns, the data line is
    N+2 chars wide (1-space pad each side) and the sep is `=`*(N+2).
    """
    # Normalize each column spec to (header, width, mode).
    cols = [(c[0], c[1], (c[2] if len(c) > 2 else "char")) for c in columns]
    widths = [w for _, w, _ in cols]
    modes  = [m for _, _, m in cols]
    total  = sum(widths) + 2   # +1 pad on each side
    sep    = "=" * total

    def emit_row(cells: list[str], header: bool = False) -> None:
        # Headers don't wrap (kept short by design).
        if header:
            parts = [_pad_right(str(c or ""), w) for c, w in zip(cells, widths)]
            _add_paragraph(doc, f"{indent} {''.join(parts)} ")
            return
        wrapped = [
            _wrap_to_cols(str(c or ""), w, mode=m)
            for c, w, m in zip(cells, widths, modes)
        ]
        n = max((len(w) for w in wrapped), default=1)
        for i in range(n):
            parts = [
                _pad_right(cell_lines[i] if i < len(cell_lines) else "", w)
                for cell_lines, w in zip(wrapped, widths)
            ]
            _add_paragraph(doc, f"{indent} {''.join(parts)} ")

    _add_paragraph(doc, f"{indent}{sep}")
    emit_row([h for h, _, _ in cols], header=True)
    _add_paragraph(doc, f"{indent}{sep}")
    for row in rows:
        emit_row(row)
    _add_paragraph(doc, f"{indent}{sep}")


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
    """三、檢測結果 — the meat of the report.

    第一類 (status=1) and 第二類 (status=2) sections each contain every
    variant of that status with the matching per-type block (SNV /
    Mito / CNV+SV) rendered inline. No "## SNV/indel" labels — those
    were template structure markers, not actual report headings.
    """
    _add_paragraph(doc, "三、檢測結果")
    _add_paragraph(doc, "  檢體說明:")
    _add_paragraph(doc, "    檢體類別：血液")
    _add_paragraph(doc, "  綜合說明:")

    statuses = report.get("status", {}) or {}
    edits    = report.get("edits", {})  or {}
    is_wgs   = (test_type or "").upper() == "WGS"

    snv_vars  = sample.get("variants", {})       or {}
    cnv_vars  = sample.get("cnv_variants", {})   or {}
    sv_vars   = sample.get("sv_variants",  {})   or {}
    mito_vars = sample.get("mito_variants", {})  or {}

    # Group by reviewer status. Each entry is ("kind", variant_dict).
    # Insertion order = (snv → mito → cnv → sv), so 第一類/第二類 list
    # SNVs first (the most common case in past reports).
    def _collect(status: str) -> list[tuple[str, dict]]:
        out: list[tuple[str, dict]] = []
        for src_kind, src in (("snv",  snv_vars),
                              ("mito", mito_vars),
                              ("cnv",  cnv_vars),
                              ("sv",   sv_vars)):
            for vid, v in src.items():
                if (statuses.get(vid) or "").strip() == status:
                    out.append((src_kind, v))
        return out

    bucket1 = _collect("1")
    bucket2 = _collect("2")
    man1 = _manual_for(report, "1")
    man2 = _manual_for(report, "2")

    # — 第一類
    _add_paragraph(doc, "    第一類：與臨床症狀相關基因之已知致病性變異位點")
    if bucket1 or man1:
        for kind, v in bucket1:
            _render_variant(doc, kind, v, tier="1",
                            edits=edits.get(v.get("id", ""), {}),
                            is_wgs=is_wgs)
        for m in man1:
            _render_manual_variant(doc, m)
    else:
        # Empty bucket — match the template's wording. We mirror past
        # reports that surface the HPO labels in the empty notice.
        hpo_labels = [r.get("label", "") or r.get("phenotype", "")
                      for r in (sample.get("patient_phenotype") or [])]
        hpo_labels = [x for x in hpo_labels if x][:8]
        _add_paragraph(doc, _summary_line("非特定", hpo_labels))
    _blank(doc)

    # — 第二類
    _add_paragraph(doc, "    第二類：其他變異位點")
    if bucket2 or man2:
        for kind, v in bucket2:
            _render_variant(doc, kind, v, tier="2",
                            edits=edits.get(v.get("id", ""), {}),
                            is_wgs=is_wgs)
        for m in man2:
            _render_manual_variant(doc, m)
    else:
        _add_paragraph(doc, "    未找到其他變異位點。")
    _blank(doc)

    # Footer recommendation (always shown)
    _add_paragraph(doc, "    建議比對臨床表徵並進行父母親與家族成員之變異位點檢測，"
                        "以釐清上述變異致病之可能性；根據家族成員變異位點檢測報告或"
                        "相關資料庫更新，可能影響變異位點ACMG判讀結果。")
    _blank(doc)
    referenced = bucket1 + bucket2
    if referenced:
        _add_paragraph(doc, "  參考資料:")
        for kind, v in referenced:
            _add_paragraph(doc, _variant_reference_text(
                kind, v, edits=edits.get(v.get("id", ""), {}), is_wgs=is_wgs,
            ))
            _blank(doc)


def _render_variant(doc, kind: str, v: dict, *, tier: str, edits: dict,
                    is_wgs: bool) -> None:
    if kind == "snv":
        _snv_variant_block(doc, v, tier=tier, edits=edits)
    elif kind == "mito":
        _mito_variant_block(doc, v, tier=tier, edits=edits)
    else:  # cnv / sv share the same template
        _cnv_variant_block(doc, v, tier=tier, is_wgs=is_wgs, edits=edits)
    _blank(doc)


def _variant_reference_text(kind: str, v: dict, *, edits: dict,
                            is_wgs: bool) -> str:
    if kind == "snv":
        return _snv_reference_text(v, edits)
    if kind == "mito":
        return _mito_reference_text(v, edits)
    omim_genes = _omim_genes(v)
    report_genes = edits.get("report_genes") or {}
    if isinstance(report_genes, dict):
        kept = [g for g in omim_genes if report_genes.get(g.get("gene"))]
        if kept:
            omim_genes = kept
    return _cnv_reference_text(v, edits, omim_genes, _sv_kind_zh(v), is_wgs)


def _render_manual_variant(doc, m: dict) -> None:
    _add_paragraph(doc, f"    {m.get('position', '')}", bold=True)
    if m.get("disease"):
        _add_paragraph(doc, f"    {m['disease']}")
    if m.get("comment"):
        _add_paragraph(doc, f"    {m['comment']}")
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
    return sig.replace("_", " ") if sig and sig != "." else "—"


def _acmg_label(variant: dict, edits_for_v: dict) -> str:
    """Return the ACMG/AMP classification label for any variant kind.

    Source priority (most-specific first):
      1. Mito reviewer manual entry   (edits.ACMG_classification_mito)
      2. SNV reviewer override        (edits.ACMG_classification)
      3. CNV/SV reviewer override     (edits.ACMG_class_sv → int 1-5)
      4. SNV pipeline value           (variant.ACMG_classification)
      5. CNV/SV pipeline value, int   (variant.acmg_class 1-5)
    Mito has no automated ACMG source — when no manual entry is set
    the function returns an empty string so the report's table cell
    renders blank, per spec.
    """
    cls = (edits_for_v.get("ACMG_classification_mito") or "").strip()
    if cls:
        return cls
    cls = (edits_for_v.get("ACMG_classification") or "").strip()
    if cls:
        return cls
    raw_cnv = (edits_for_v.get("ACMG_class_sv") or "").strip()
    if raw_cnv:
        try:
            n = int(raw_cnv)
            if n in _CNV_ACMG_INT_TO_LABEL:
                return _CNV_ACMG_INT_TO_LABEL[n]
        except ValueError:
            return raw_cnv
    cls = (variant.get("ACMG_classification") or "").strip()
    if cls:
        return cls
    n = variant.get("acmg_class")
    if isinstance(n, int) and n in _CNV_ACMG_INT_TO_LABEL:
        return _CNV_ACMG_INT_TO_LABEL[n]
    return ""


def _patho_sentence(acmg_class: str) -> str:
    """Note-2 sentence in §三 driven by the actual ACMG class — so a
    VUS variant reads "此為不確定意義之變異位點..." not "此為致病性...".
    """
    cls = (acmg_class or "").strip().lower().replace("_", " ")
    if cls == "pathogenic":
        return "此為致病性之變異位點，與臨床症狀相關。"
    if cls in ("likely pathogenic",):
        return "此為疑似致病性之變異位點，與臨床症狀相關。"
    if "pathogenic/likely" in cls or "likely pathogenic/pathogenic" in cls:
        return "此為致病性 / 疑似致病性之變異位點，與臨床症狀相關。"
    if cls in ("uncertain significance", "vus"):
        return ("此為不確定意義之變異位點，其臨床意義須由醫師配合其他"
                "相關資料進行最佳綜合判斷。")
    if cls in ("likely benign",):
        return "此為可能良性之變異位點。"
    if cls == "benign":
        return "此為良性之變異位點。"
    # Mito MITOMAP statuses
    if any(k in cls for k in ("cfrm", "confirmed", "[p]", "[lp]")):
        return "此為致病性之變異位點，與臨床症狀相關。"
    if "reported" in cls:
        return ("此為與疾病相關之變異位點，其臨床意義須由醫師配合其他"
                "相關資料進行最佳綜合判斷。")
    # Unknown / empty
    return "此變異位點之臨床意義須由醫師配合其他相關資料進行最佳綜合判斷。"


def _picked_disease_for_snv(v: dict, edits: dict) -> str:
    """Return the reviewer-picked Disease cell, falling back to the first."""
    picked = edits.get("report_diseases") or {}
    if isinstance(picked, dict):
        for idx in range(1, 6):
            if picked.get(str(idx)) or picked.get(idx):
                disease = (v.get(f"Disease{idx}") or "").strip()
                if disease and disease != "NA":
                    return disease
    for idx in range(1, 6):
        disease = (v.get(f"Disease{idx}") or "").strip()
        if disease and disease != "NA":
            return disease
    return (v.get("OMIM_disease") or "").strip()


def _disease_info(disease: str) -> tuple[str, str, str]:
    """Parse disease name, disease-specific inheritance, and phenotype MIM."""
    first_line = (disease or "").splitlines()[0].strip()
    inheritance_match = re.search(
        r"\((AD|AR|XLD|XLR|XL|YL|MT|Mi|DR|DD|Smu|Mu|Isol)"
        r"(?:\s*[/,;]\s*(?:AD|AR|XLD|XLR|XL|YL|MT|Mi|DR|DD|Smu|Mu|Isol))*\)",
        first_line,
        re.IGNORECASE,
    )
    phenotype_mim_match = re.search(r"\((\d{6})\)", first_line)
    inheritance = inheritance_match.group(0)[1:-1] if inheritance_match else ""
    phenotype_mim = phenotype_mim_match.group(1) if phenotype_mim_match else ""
    metadata_starts = [
        match.start() for match in (phenotype_mim_match, inheritance_match)
        if match
    ]
    name = first_line[:min(metadata_starts)] if metadata_starts else first_line
    return name.rstrip(" :,;").strip(), inheritance, phenotype_mim


def _omim_block_for_snv(v: dict, edits: dict) -> str:
    """『GENE 為 DISEASE 的致病基因之一，其遺傳模式屬於 X
    (Phenotype MIM number: M)』 — fall back to whatever is present."""
    gene = v.get("gene_symbol") or v.get("GENE") or "?"
    disease_str, inheritance, phenotype_mim = _disease_info(
        _picked_disease_for_snv(v, edits)
    )
    inh_zh = _inheritance_zh(inheritance)
    mim = phenotype_mim or (v.get("OMIM_id") or "").strip()

    parts = [f"{gene}為"]
    if disease_str:
        parts.append(f"{disease_str}的致病基因之一")
    else:
        parts.append("此疾病的致病基因之一")
    if inh_zh:
        parts.append(f"，其遺傳模式屬於{inh_zh}")
    if mim:
        parts.append(f" (Phenotype MIM number: {mim})")
    return "".join(parts) + "。"


def _snv_variant_block(doc, v: dict, *, tier: str, edits: dict) -> None:
    gene = v.get("gene_symbol") or "?"
    tx   = v.get("transcript")  or v.get("MANE_SELECT") or ""
    _add_paragraph(doc, f"    {gene} ({tx})", bold=True)
    rs    = v.get("rs_id") or v.get("RS_ID") or ""
    struc = _structure_label(v)
    hgvs_c = _strip_tx_prefix(v.get("hgvs_c") or v.get("HGVS_C") or "")
    hgvs_p = _strip_tx_prefix(v.get("hgvs_p") or v.get("HGVS_P") or "")
    nuc    = hgvs_c + (f"({hgvs_p})" if hgvs_p else "")
    # 基因型 column stays in English per spec (Heterozygous / Homozygous)
    zyg    = _zygosity_long(v.get("zygosity", ""))
    clnsig = _clinvar_label(v)
    acmg   = _acmg_label(v, edits)

    _ascii_table(doc, columns=[
        ("類別",          5),
        ("基因",          9),
        ("RS ID",         8),
        ("結構",          9),
        ("核苷酸",       14, "hgvs"),
        ("基因型",       13),
        ("ClinVar",      13, "token"),
        ("ACMG&AMP指引", 13, "token"),
    ], rows=[[tier, gene, rs, struc, nuc, zyg, clnsig, acmg]])

    _add_paragraph(doc, f"    1. {_omim_block_for_snv(v, edits)}")
    _add_paragraph(doc, f"    2. {_patho_sentence(acmg)}")

def _snv_reference_text(v: dict, edits: dict) -> str:
    gene = v.get("gene_symbol") or "?"
    # Drop the "NM_xxx:" transcript prefix from HGVS values here too —
    # the transcript already appears beside the gene name on the block
    # header, no need to repeat it in the ref text.
    hgvs_c = _strip_tx_prefix(v.get("hgvs_c") or v.get("HGVS_C") or "")
    hgvs_p = _strip_tx_prefix(v.get("hgvs_p") or v.get("HGVS_P") or "")
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
                 "在疾病資料庫 (ClinVar) 中未被報導過")

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

def _mito_variant_block(doc, v: dict, *, tier: str, edits: dict) -> None:
    gene  = v.get("gene_symbol") or "?"
    _add_paragraph(doc, f"    {gene} (NC_012920.1)", bold=True)
    hgvs_m = v.get("HGVS_M") or v.get("id") or ""
    aa     = _mito_aa_to_hgvsp(v.get("aa_change") or "")
    nuc    = hgvs_m + (f"({aa})" if aa else "")
    het    = v.get("heteroplasmy")
    try:
        het_s = f"{float(het) * 100:.1f}%" if het not in (None, "", ".") else "—"
    except (TypeError, ValueError):
        het_s = "—"
    clnsig = _clinvar_label(v)
    acmg   = _acmg_label(v, edits)

    # RS ID column dropped — MITOMAP doesn't expose rsIDs in our
    # current adapter. If the pipeline starts emitting one, conditionally
    # add the column back.
    _ascii_table(doc, columns=[
        ("類別",          5),
        ("基因",          9),
        ("核苷酸",       14, "hgvs"),
        ("異質性比例",   13),
        ("ClinVar",      13, "token"),
        ("ACMG&AMP指引", 13, "token"),
    ], rows=[[tier, gene, nuc, het_s, clnsig, acmg]])

    mt_disease = (v.get("mitomap_disease") or "").strip()
    # 1. 致病基因之一 — MITOMAP disease name; 遺傳模式 fixed to 粒線體遺傳; no MIM
    if mt_disease:
        _add_paragraph(doc, f"    1. {gene}為{mt_disease}的致病基因之一，其遺傳模式屬於粒線體遺傳。")
    else:
        _add_paragraph(doc, f"    1. {gene}為粒線體疾病之相關基因，其遺傳模式屬於粒線體遺傳。")
    _add_paragraph(doc, f"    2. {_patho_sentence(acmg)}")
    _add_paragraph(doc, "    3. 此檢驗結果須由臨床醫師來分析粒線體DNA基因檢驗結果與受檢者臨床症狀的相關性"
                        "並考量家族史的關聯性。此外，粒線體DNA基因變異有組織的特異性，粒線體DNA基因致病變異的"
                        "異質性（heteroplasmy）在不同組織間會有差異。")
    _add_paragraph(doc, "    4. 對於所有粒線體DNA基因變異所造成的臨床症狀表現的差異，可能取決於下列三項因素："
                        "（1）異質性，即致病變異及正常粒線體 DNA 存在量的比例（2）變異的粒線體 DNA 於組織的"
                        "分布情形（3）限界效應（Threshold effect），即每種組織對於氧化壓力代謝影響的易受性。"
                        "由於患者情況各異，此檢驗結果須由臨床醫師判讀檢驗結果與受檢者臨床症狀的相關性。")
def _mito_reference_text(v: dict, edits: dict) -> str:
    gene = v.get("gene_symbol") or "?"
    hgvs_m = v.get("HGVS_M") or v.get("id") or ""
    aa = _mito_aa_to_hgvsp(v.get("aa_change") or "")
    nuc = hgvs_m + (f"({aa})" if aa else "")
    cq_zh = _consequence_zh(v.get("consequence", ""))
    clinvar = _clinvar_label(v)
    acmg    = _acmg_label(v, edits)
    clnv_sent = (f"在疾病資料庫 (ClinVar) 中此變異位點被報導為「{clinvar}」"
                 if clinvar and clinvar != "—"
                 else "在疾病資料庫 (ClinVar) 中未被報導過")
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

def _build_coords(v: dict) -> str:
    """Prefer the adapter's `coords` field; fall back to constructing
    `chrN:POS-END` from CHROM / POS / END when it's missing or blank.
    Always include the `chr` prefix on the chromosome name (UCSC
    convention; matches the UI header `[GRCh38] chr1:...`).
    """
    coords = (v.get("coords") or "").strip()
    if not coords:
        chrom = (v.get("CHROM") or "").strip()
        pos   = v.get("POS")
        end   = v.get("END")
        if chrom and pos is not None and end is not None:
            coords = f"{chrom}:{pos}-{end}"
    if coords and not coords.startswith("chr"):
        coords = "chr" + coords
    return coords


def _location_zh(g: dict) -> str:
    """AnnotSV's per-gene `Location` field (e.g. 'txStart-txEnd',
    'txStart-intron3', 'exon2-exon8') → 中文 phrase for the report.

    Whole-gene coverage ('txStart-txEnd') always returns 「整個區域」.
    Partial coverage parses each endpoint:
      txStart → 基因起始, txEnd → 基因末端,
      exonN   → "Exon N 區域", intronN → "Intron N 區域"
    Joined with 「至」. A leading ASCII letter on the second piece gets
    a leading space ("基因起始至 Intron 3 區域") so Latin tokens don't
    sit right against the Chinese 至 separator.
    Unrecognised endpoints pass through verbatim so weird AnnotSV
    values stay visible for review.
    """
    loc = (g.get("location") or "").strip() if isinstance(g, dict) else ""
    if not loc or loc == "txStart-txEnd":
        return "整個區域"
    parts = loc.split("-", 1)
    if len(parts) != 2:
        return loc

    def piece(p: str) -> str:
        p = p.strip()
        if p == "txStart": return "基因起始"
        if p == "txEnd":   return "基因末端"
        m = re.match(r"^exon(\d+)$", p, re.IGNORECASE)
        if m: return f"Exon {m.group(1)} 區域"
        m = re.match(r"^intron(\d+)$", p, re.IGNORECASE)
        if m: return f"Intron {m.group(1)} 區域"
        return p

    p1 = piece(parts[0])
    p2 = piece(parts[1])
    sep = "至" + (" " if p2 and "A" <= p2[0] <= "Z" else "")
    return f"{p1}{sep}{p2}"


def _gene_loc_phrase(gname: str, loc_zh: str) -> str:
    """Build「GENE 基因之 …」 / 「GENE 基因之 Exon X 區域」 — always
    with 之, with a space inserted when the location phrase starts
    with a Latin letter so "基因之 Exon" reads cleanly."""
    sep = "之" + (" " if loc_zh and "A" <= loc_zh[0] <= "Z" else "")
    return f"{gname} 基因{sep}{loc_zh}"


def _omim_genes(v: dict) -> list[dict]:
    """Filter the variant's `genes` list to those carrying an OMIM_ID —
    the only ones worth surfacing in the diagnostic report."""
    out = []
    for g in (v.get("genes") or []):
        if not isinstance(g, dict):
            continue
        if (g.get("omim_id") or "").strip():
            out.append(g)
    return out


def _coords_str(v: dict, is_wgs: bool) -> str:
    """WES uses imprecise breakpoint notation `(?_start)_(end_?)`;
    WGS uses the precise start_end form. The adapter's `coords` is the
    raw AnnotSV CHROM:POS-END so we annotate accordingly.
    """
    coords = _build_coords(v)
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
    sv_upper  = (v.get("sv_type") or "").upper()
    cn        = v.get("copy_number")
    cn_s      = "—" if cn in (None, "", ".") else str(cn)
    # 基因型 only meaningful for DEL — DUP / INS / INV reuse the column
    # for "—" so the reviewer doesn't read het/hom on a copy gain.
    raw_zyg   = (v.get("zygosity", "") or "").strip()
    if sv_upper == "DEL":
        zyg = _zygosity_long(raw_zyg) or "—"
    else:
        zyg = "—"
    acmg      = _acmg_label(v, edits)
    chrom     = (v.get("CHROM") or "").replace("chr", "")
    coords_disp = coords.split(":", 1)[1] if ":" in coords else coords

    sv_tag = "del" if sv_upper == "DEL" else ("dup" if sv_upper == "DUP" else "")
    _add_paragraph(doc, f"    [GRCh38] {coords}{sv_tag}", bold=True)
    _ascii_table(doc, columns=[
        ("類別",          5),
        ("染色體",        7),
        ("變異位置",     26),
        ("拷貝數",        7),
        ("基因型",       13),
        ("ACMG&AMP指引", 13, "token"),
    ], rows=[[tier, chrom, coords_disp, cn_s, zyg, acmg]])

    # 1. 片段位置描述 — only surface OMIM-tagged genes (per spec);
    # decide single-gene vs multi-gene template by how many are left.
    chrom_num  = (v.get("CHROM") or "").replace("chr", "")
    omim_genes = _omim_genes(v)

    # Reviewer can prune which genes appear via edits.report_genes
    # (frontend toggle on the CNV/SV card). Honour the prune only when
    # at least one OMIM gene survives.
    report_genes = edits.get("report_genes") or {}
    if isinstance(report_genes, dict):
        kept = [g for g in omim_genes if report_genes.get(g.get("gene"))]
        if kept:
            omim_genes = kept

    kind_zh = _CNV_KIND_ZH.get((v.get("sv_type") or "").upper(), "變異")

    if len(omim_genes) == 1:
        g = omim_genes[0]
        gname = g.get("gene") or "?"
        loc_zh = _location_zh(g)
        _add_paragraph(doc, f"    1. 此片段位於第 {chrom_num} 號染色體上 {_gene_loc_phrase(gname, loc_zh)}。")
        # 2. OMIM phenotype + inheritance, per-gene
        ph, ph_inheritance, phenotype_mim = _disease_info(
            (g.get("omim_phenotype") or "").strip()
        )
        inh = _inheritance_zh(
            ph_inheritance or g.get("omim_inheritance", "") or ""
        )
        mim = phenotype_mim or (g.get("omim_id") or "").strip()
        bits = [f"{gname}為"]
        bits.append(f"{ph}的致病基因之一" if ph else "此疾病的致病基因之一")
        if inh: bits.append(f"，其遺傳模式屬於{inh}")
        if mim: bits.append(f" (Phenotype MIM number: {mim})")
        _add_paragraph(doc, f"    2. {''.join(bits)}。")
    elif len(omim_genes) > 1:
        names = [g.get("gene", "") for g in omim_genes[:10] if g.get("gene")]
        _add_paragraph(doc, f"    1. 此片段位於第 {chrom_num} 號染色體上，"
                            f"包含 {', '.join(names)} 等 OMIM 疾病基因。")
    else:
        # No OMIM-tagged gene in this CNV — list whatever genes are
        # present so the reviewer has context, but skip the OMIM line.
        all_names = [
            g.get("gene", "") for g in (v.get("genes") or [])
            if isinstance(g, dict) and g.get("gene")
        ][:10]
        if all_names:
            _add_paragraph(doc, f"    1. 此片段位於第 {chrom_num} 號染色體上，"
                                f"涵蓋 {', '.join(all_names)} 等基因區域（無 OMIM 疾病基因紀錄）。")
        else:
            _add_paragraph(doc, f"    1. 此片段位於第 {chrom_num} 號染色體上。")

    note_idx = 3
    _add_paragraph(doc, f"    {note_idx}. {_patho_sentence(acmg)}")
    note_idx += 1
    if not is_wgs:
        _add_paragraph(doc, f"    {note_idx}. 由於此檢驗技術為全外顯子定序，"
                            "若缺失片段之斷點(Breakpoints)發生於內含子(Intron) ，"
                            "則無法明確判別起始及末端位置。")

def _cnv_reference_text(v: dict, edits: dict, omim_genes: list[dict],
                        kind_zh: str, is_wgs: bool) -> str:
    coords = _coords_str(v, is_wgs) or _build_coords(v)
    # Zygosity in the sentence only when DEL — DUP / INS / INV don't
    # carry meaningful zygosity (a triploid duplication isn't "het").
    sv_upper = (v.get("sv_type") or "").upper()
    if sv_upper == "DEL":
        zyg_text = _zygosity_zh(v.get("zygosity", "")) or "異合子"
        zyg_phrase = f"之{zyg_text}"
    else:
        zyg_phrase = "之"
    acmg    = _acmg_label(v, edits)

    if len(omim_genes) == 1:
        g = omim_genes[0]
        gname = g.get("gene") or "?"
        loc_zh = _location_zh(g)
        span_desc = f"此段{kind_zh}涵蓋 {_gene_loc_phrase(gname, loc_zh)}"
    elif len(omim_genes) > 1:
        names = [g.get("gene", "") for g in omim_genes[:10] if g.get("gene")]
        span_desc = f"此段{kind_zh}涵蓋 {', '.join(names)} 等 OMIM 疾病基因"
    else:
        all_names = [
            g.get("gene", "") for g in (v.get("genes") or [])
            if isinstance(g, dict) and g.get("gene")
        ][:10]
        span_desc = (f"此段{kind_zh}涵蓋 {', '.join(all_names)} 等基因區域（無 OMIM 疾病基因紀錄）"
                     if all_names else "此段未涵蓋已知疾病基因")

    return (
        f"    在個案之檢體中，檢測到位於 {coords} {zyg_phrase}片段{kind_zh}變異，"
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
                        "並且主要列入致病(Pathogenic) 及疑似致病 (Likely pathogenic) 變異；"
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
    # selected_panels is a list of {name, weight} dicts (see
    # phenotype_io.parse) — pull the name field.
    panel_entries: list = sample.get("selected_panels") or []

    # Build [(display_name, [genes...])] preserving order.
    sections: list[tuple[str, list[str]]] = []
    for r in hpo_rows:
        hid = r.get("phenotype") or ""
        label = r.get("label") or _hpo_label_for(hid)
        if not hid: continue
        genes = _genes_for_term_or_panel(hid)
        sections.append((label, genes))
    for entry in panel_entries:
        pname = entry.get("name") if isinstance(entry, dict) else str(entry)
        if not pname: continue
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
        _add_paragraph(doc, gene_str)
        return

    # grouped (default)
    for name, gs in sections:
        _add_paragraph(doc, f"{name}:")
        if gs:
            _add_paragraph(doc, ", ".join(gs))
        else:
            _add_paragraph(doc, "（無對應基因）")


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
    _apply_page_margins(doc)

    _section_test_info(doc, test_type)
    _section_panel_set(doc)
    _section_results(doc, sample, report, test_type)
    _section_methods(doc, test_type)
    _section_annotations(doc, sample, gene_list_mode)

    buf = io.BytesIO()
    doc.save(buf)
    return buf.getvalue()
