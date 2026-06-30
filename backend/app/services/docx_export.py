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
from . import (
    cnv_sv_merge,
    hpo_ontology,
    panel_deadzone,
    phenotype_scorer,
    report_store,
    sample_loader,
)

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
    mode="buffered" → greedy char wrap with one trailing display column
                      reserved as whitespace before the next table cell.
    mode="token-buffered" → token wrap with one trailing display column
                            reserved before the next table cell.
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

    if mode == "token-buffered":
        return _wrap_to_cols(s, max(1, width - 1), mode="token")

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

    if mode == "buffered":
        return _wrap_to_cols(s, max(1, width - 1), mode="char")

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

def _section_test_info(doc, test_type: str, *, health: bool = False) -> None:
    """一、檢驗項目"""
    is_wgs = (test_type or "").upper() == "WGS"
    if health:
        label = "次世代定序全基因體定序檢測" if is_wgs else "次世代定序全外顯子定序檢測"
    else:
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
    cnv_vars, sv_vars = cnv_sv_merge.apply_confirmed_merges(
        sample.get("cnv_variants", {}) or {},
        sample.get("sv_variants", {}) or {},
        report.get("cnv_sv_merges") or [],
    )
    mito_vars = sample.get("mito_variants", {})  or {}

    # Group by reviewer status. Each entry is ("kind", variant_dict).
    # Insertion order = (snv → mito → cnv → sv), so 第一類/第二類 list
    # SNVs first (the most common case in past reports). CNV/SV entries
    # are sorted inside their source by combined phenotype + AnnotSV score.
    def _ranked_items(src: dict, kind: str):
        items = list(src.items())
        if kind in {"cnv", "sv"}:
            items.sort(key=lambda item: (
                -float(item[1].get("cnv_sv_sort_score") or -999),
                str(item[0]),
            ))
        return items

    def _collect(status: str) -> list[tuple[str, dict]]:
        out: list[tuple[str, dict]] = []
        for src_kind, src in (("snv",  snv_vars),
                              ("mito", mito_vars),
                              ("cnv",  cnv_vars),
                              ("sv",   sv_vars)):
            for vid, v in _ranked_items(src, src_kind):
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
      4. SNV GeneBe value             (variant.genebe_acmg_class)
      5. SNV pipeline value           (variant.ACMG_classification)
      6. CNV/SV pipeline value, int   (variant.acmg_class 1-5)
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
    cls = (variant.get("genebe_acmg_class") or "").strip()
    if cls:
        return cls
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
    lines = (disease or "").splitlines()
    first_line = lines[0].strip() if lines else ""
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


def _picked_clinvar_diseases_for_mito(v: dict, edits: dict) -> list[str]:
    diseases = v.get("clinvar_diseases") or []
    if not isinstance(diseases, list):
        diseases = []
    picked = edits.get("report_diseases_clinvar") or {}
    out: list[str] = []
    if isinstance(picked, dict):
        for idx, disease in enumerate(diseases):
            if picked.get(str(idx)) or picked.get(idx):
                d = str(disease or "").strip()
                if d:
                    out.append(d)
    return out


def _mito_disease_text(v: dict, edits: dict) -> str:
    picked = _picked_clinvar_diseases_for_mito(v, edits)
    if picked:
        return "、".join(picked)
    return ""


def _selected_snv_transcript(v: dict, edits: dict) -> dict:
    options = v.get("transcript_options") or []
    selected = str(edits.get("selected_transcript_key") or "").strip()
    if selected:
        for opt in options:
            if str(opt.get("key") or "") == selected:
                return opt
    default_key = str(v.get("default_transcript_key") or "").strip()
    if default_key:
        for opt in options:
            if str(opt.get("key") or "") == default_key:
                return opt
    return {}


def _snv_tx_field(v: dict, edits: dict, field: str, *fallbacks: str) -> str:
    selected = _selected_snv_transcript(v, edits)
    value = selected.get(field)
    if value not in (None, ""):
        return str(value)
    for key in fallbacks:
        value = v.get(key)
        if value not in (None, ""):
            return str(value)
    value = v.get(field)
    return "" if value in (None, "") else str(value)


def _strip_hgvs_prefix(value: str) -> str:
    text = str(value or "").strip()
    if not text or text in {".", "-", "NA", "N/A"}:
        return ""
    for marker in ("c.", "n.", "m.", "g.", "p."):
        idx = text.find(marker)
        if idx >= 0:
            return text[idx:]
    return text.split(":", 1)[1] if ":" in text else text


def _snv_transcript_label(v: dict, edits: dict) -> str:
    enst = _snv_tx_field(v, edits, "ensembl_transcript")
    refseq = _snv_tx_field(v, edits, "refseq_transcript")
    if enst and refseq and enst != refseq:
        return f"{enst}; {refseq}"
    return refseq or enst or _snv_tx_field(v, edits, "transcript", "MANE_SELECT")


def _snv_variant_block(doc, v: dict, *, tier: str, edits: dict) -> None:
    gene = _snv_tx_field(v, edits, "gene_symbol") or "?"
    tx   = _snv_transcript_label(v, edits)
    _add_paragraph(doc, f"    {gene} ({tx})", bold=True)
    rs    = v.get("rs_id") or v.get("RS_ID") or ""
    struc = _structure_label({
        **v,
        "exon": _snv_tx_field(v, edits, "exon"),
        "intron": _snv_tx_field(v, edits, "intron"),
    })
    hgvs_c = _strip_hgvs_prefix(_snv_tx_field(v, edits, "HGVS_C", "hgvs_c"))
    hgvs_p = _strip_hgvs_prefix(_snv_tx_field(v, edits, "HGVS_P", "hgvs_p"))
    nuc    = hgvs_c + (f"({hgvs_p})" if hgvs_p else "")
    # 基因型 column stays in English per spec (Heterozygous / Homozygous)
    zyg    = _zygosity_long(v.get("zygosity", ""))
    clnsig = _clinvar_label(v)
    acmg   = _acmg_label(v, edits)

    _ascii_table(doc, columns=[
        ("類別",          5),
        ("基因",          9),
        ("RS ID",         7, "buffered"),
        ("結構",          9),
        ("核苷酸",       14, "hgvs"),
        ("基因型",       13),
        ("ClinVar",      16, "token-buffered"),
        ("ACMG&AMP指引", 12, "token"),
    ], rows=[[tier, gene, rs, struc, nuc, zyg, clnsig, acmg]])

    _add_paragraph(doc, f"    1. {_omim_block_for_snv(v, edits)}")
    _add_paragraph(doc, f"    2. {_patho_sentence(acmg)}")

def _snv_reference_text(v: dict, edits: dict) -> str:
    gene = _snv_tx_field(v, edits, "gene_symbol") or "?"
    # Drop the "NM_xxx:" transcript prefix from HGVS values here too —
    # the transcript already appears beside the gene name on the block
    # header, no need to repeat it in the ref text.
    hgvs_c = _strip_hgvs_prefix(_snv_tx_field(v, edits, "HGVS_C", "hgvs_c"))
    hgvs_p = _strip_hgvs_prefix(_snv_tx_field(v, edits, "HGVS_P", "hgvs_p"))
    nuc    = hgvs_c + (f" ({hgvs_p})" if hgvs_p else "")
    cq_zh  = _consequence_zh(_snv_tx_field(v, edits, "Consequence"))
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

    mt_disease = _mito_disease_text(v, edits)
    # 1. 致病基因之一 — ClinVar disease name; 遺傳模式 fixed to 粒線體遺傳; no MIM
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


def _cnv_report_disease(edits: dict) -> str:
    """Reviewer-entered CNV/SV disease label for the formal report."""
    return str(edits.get("disease") or "").strip()


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
        ("變異位置",     26, "buffered"),
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

    next_idx = 2
    disease_override = _cnv_report_disease(edits)
    if len(omim_genes) == 1:
        g = omim_genes[0]
        gname = g.get("gene") or "?"
        loc_zh = _location_zh(g)
        _add_paragraph(doc, f"    1. 此片段位於第 {chrom_num} 號染色體上 {_gene_loc_phrase(gname, loc_zh)}。")
        # 2. OMIM phenotype + inheritance, per-gene
        ph, ph_inheritance, phenotype_mim = _disease_info(
            disease_override or (g.get("omim_phenotype") or "").strip()
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
        next_idx = 3
    elif len(omim_genes) > 1:
        names = [g.get("gene", "") for g in omim_genes[:10] if g.get("gene")]
        disease_suffix = f"，與 {disease_override} 相關" if disease_override else ""
        _add_paragraph(doc, f"    1. 此片段位於第 {chrom_num} 號染色體上，"
                            f"包含 {', '.join(names)} 等 OMIM 疾病基因"
                            f"{disease_suffix}。")
    else:
        disease_suffix = f"，與 {disease_override} 相關" if disease_override else ""
        _add_paragraph(doc, f"    1. 此片段位於第 {chrom_num} 號染色體上，"
                            "未涵蓋 OMIM 疾病相關基因"
                            f"{disease_suffix}。")

    _add_paragraph(doc, f"    {next_idx}. {_patho_sentence(acmg)}")
    next_idx += 1
    if not is_wgs:
        _add_paragraph(doc, f"    {next_idx}. 由於此檢驗技術為全外顯子定序，"
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
        span_desc = f"此段{kind_zh}未涵蓋 OMIM 疾病相關基因"

    disease = _cnv_report_disease(edits)
    disease_sent = f"此變異與「{disease}」相關。" if disease else ""
    return (
        f"    在個案之檢體中，檢測到位於 {coords} {zyg_phrase}片段{kind_zh}變異，"
        f"{span_desc}。"
        f"{disease_sent}"
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
    final_note_number = 6 if is_wgs else 5
    if is_wgs:
        _add_paragraph(doc, "  5. 本實驗方法以次世代方法定序粒線體DNA基因序列，"
                            "變異點位判讀之cut-off值定為5%異質性（heteroplasmy）。")
    _add_paragraph(doc, f"  {final_note_number}. 本檢測報告僅供醫療專業人員參考，需配合其他相關臨床資料與家族成員之相關檢驗。"
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
    genes = []
    for raw in phenotype_scorer._HPO_TO_GENES.get(key, set()):
        gene, _hid = panel_deadzone.canonical_gene_symbol(raw)
        if panel_deadzone.is_disease_associated_gene(gene):
            genes.append(gene)
    return sorted(set(genes))


def _gene_list_label(gene: str, test_type: str) -> str:
    return gene


def _add_dead_zone_gene_list_note(doc, threshold: int) -> None:
    _blank(doc)
    _add_paragraph(doc, f"註：括號中標示之 exon 為 cohort dead-zone，代表該 exon coverage 低於本檢測判讀門檻（<{threshold}X）。")


def _render_gene_list(doc, sample: dict, mode: str) -> None:
    """mode = 'grouped' → one paragraph per HPO term / panel,
       mode = 'merged'  → single deduped list.
    """
    hpo_rows: list = sample.get("patient_phenotype") or []
    # selected_panels is a list of {name, weight} dicts (see
    # phenotype_io.parse) — pull the name field.
    panel_entries: list = sample.get("selected_panels") or []
    test_type = ((sample.get("meta") or {}).get("Test") or "WES").upper()
    include_docx_dead_zone = False

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
        gene_str = ", ".join(_gene_list_label(g, test_type) for g in sorted(merged))
        _add_paragraph(doc, gene_str)
        if include_docx_dead_zone:
            _add_dead_zone_gene_list_note(doc, threshold)
        return

    # grouped (default)
    for idx, (name, gs) in enumerate(sections):
        if idx:
            _blank(doc)
        _add_paragraph(doc, f"{name}:")
        if gs:
            _add_paragraph(doc, ", ".join(_gene_list_label(g, test_type) for g in gs))
        else:
            _add_paragraph(doc, "（無對應基因）")
    if include_docx_dead_zone:
        _add_dead_zone_gene_list_note(doc, threshold)


# ── Health screening DOCX ─────────────────────────────────────────

_NO_HEALTH_VARIANT_TEXT = "於本次檢測涵蓋之基因中未檢出致病或疑似致病之基因變異。"

_HEALTH_SECTION_ORDER = [
    ("acmg_sf", "第一類：重大可預防疾病風險基因"),
    ("lipid_fh", "第二類：血脂相關基因"),
    ("hereditary_cancer", "第三類：腫瘤相關基因"),
    ("stroke", "第四類：中風相關基因"),
    ("carrier", "第五類：帶因者篩查"),
    ("proactive", "第六類：主動篩查"),
    ("pgx", "藥物基因體學"),
]

_HEALTH_ACMG_GENE_LIST_TITLE = (
    "第一類：重大可預防疾病風險基因（參考美國遺傳醫學會 "
    "(American College of Medical Genetics and Genomics) 於 2025 年所公告之"
    "次發現 (Secondary findings) 基因清單 v3.3 版；PMID: 40568962）"
)


def _meta_value(meta: dict, key: str) -> str:
    value = meta.get(key)
    return "" if value in (None, "") else str(value)


def _meta_pair(label: str, value: str, width: int = 39) -> str:
    return _pad_right(f"{label}: {value or '—'}", width)


def _section_health_patient_header(doc, meta: dict) -> None:
    """Basic case metadata above the health-report sections."""
    rows = [
        (("檢體編號", _meta_value(meta, "LIS_ID")), ("姓名", _meta_value(meta, "Name"))),
        (("病歷號", _meta_value(meta, "MRN")), ("性別", _meta_value(meta, "Sex"))),
        (("出生日期", _meta_value(meta, "DOB")), ("簽收日期", _meta_value(meta, "SignReceivedAt"))),
        (("科別", _meta_value(meta, "Department")), ("開單醫師", _meta_value(meta, "Physician"))),
    ]
    _add_paragraph(doc, "受檢者資料", bold=True)
    for left, right in rows:
        _add_paragraph(doc, f"{_meta_pair(*left)}{_meta_pair(*right)}")
    _blank(doc)

_ACMG_SF_GROUPS = [
    {
        "key": "lipid",
        "title": "血脂相關基因",
        "genes": ["APOB", "LDLR", "PCSK9"],
    },
    {
        "key": "tumor",
        "title": "腫瘤相關基因",
        "genes": [
            "APC", "BMPR1A", "BRCA1", "BRCA2", "MAX", "MEN1", "MLH1",
            "MSH2", "MSH6", "MUTYH", "NF2", "PALB2", "PMS2", "PTEN",
            "RB1", "RET", "SDHAF2", "SDHB", "SDHC", "SDHD", "SMAD4",
            "STK11", "TMEM127", "TP53", "TSC1", "TSC2", "VHL", "WT1",
        ],
    },
    {
        "key": "cardiovascular",
        "title": (
            "心血管疾病相關基因，包含心肌病變相關基因（ACTC1, BAG3, DES, DSC2, "
            "DSG2, DSP, FLNC, LMNA, MYBPC3, MYH7, MYL2, MYL3, PKP2, PLN, "
            "PRKAG2, RBM20, TMEM43, TNNC1, TNNI3, TNNT2, TPM1, TTN）、"
            "心律不整相關基因（CALM1, CALM2, CALM3, CASQ2, KCNH2, KCNQ1, "
            "RYR2, SCN5A, TRDN）、主動脈及血管疾病相關基因（ACTA2, COL3A1, "
            "FBN1, MYH11, SMAD3, TGFBR1, TGFBR2, ACVRL1, ENG）"
        ),
        "genes": [
            "ACTC1", "BAG3", "DES", "DSC2", "DSG2", "DSP", "FLNC",
            "LMNA", "MYBPC3", "MYH7", "MYL2", "MYL3", "PKP2", "PLN",
            "PRKAG2", "RBM20", "TMEM43", "TNNC1", "TNNI3", "TNNT2",
            "TPM1", "TTN", "CALM1", "CALM2", "CALM3", "CASQ2", "KCNH2",
            "KCNQ1", "RYR2", "SCN5A", "TRDN", "ACTA2", "COL3A1",
            "FBN1", "MYH11", "SMAD3", "TGFBR1", "TGFBR2", "ACVRL1",
            "ENG",
        ],
    },
    {
        "key": "metabolic_endocrine",
        "title": "代謝與內分泌疾病相關基因",
        "genes": ["ABCD1", "ATP7B", "BTD", "CYP27A1", "GAA", "GLA", "HFE", "HNF1A", "OTC"],
    },
    {
        "key": "anesthesia",
        "title": "麻醉用藥風險相關基因",
        "genes": ["CACNA1S", "RYR1"],
    },
    {
        "key": "other",
        "title": "其它基因",
        "genes": ["RPE65", "TTR"],
    },
]

_ACMG_SF_DISEASES = {
    "ABCD1": [{"disease": "X-linked adrenoleukodystrophy", "inheritance": "XL"}],
    "ACTA2": [{"disease": "Familial thoracic aortic aneurysm", "inheritance": "AD"}],
    "ACTC1": [{"disease": "Hypertrophic cardiomyopathy", "inheritance": "AD"}],
    "ACVRL1": [{"disease": "Hereditary hemorrhagic telangiectasia", "inheritance": "AD"}],
    "APC": [{"disease": "Familial adenomatous polyposis", "inheritance": "AD"}],
    "APOB": [{"disease": "Familial hypercholesterolemia", "inheritance": "AD"}],
    "ATP7B": [{"disease": "Wilson disease", "inheritance": "AR"}],
    "BAG3": [{"disease": "Dilated cardiomyopathy", "inheritance": "AD"}, {"disease": "Myofibrillar myopathy", "inheritance": "AD"}],
    "BMPR1A": [{"disease": "Juvenile polyposis syndrome", "inheritance": "AD"}],
    "BRCA1": [{"disease": "Hereditary breast and ovarian cancer", "inheritance": "AD"}],
    "BRCA2": [{"disease": "Hereditary breast and ovarian cancer", "inheritance": "AD"}],
    "BTD": [{"disease": "Biotinidase deficiency", "inheritance": "AR"}],
    "CACNA1S": [{"disease": "Malignant hyperthermia", "inheritance": "AD"}],
    "CALM1": [{"disease": "Long-QT syndrome type 14", "inheritance": "AD"}, {"disease": "Catecholaminergic polymorphic ventricular tachycardia", "inheritance": "AD"}],
    "CALM2": [{"disease": "Long-QT syndrome type 15", "inheritance": "AD"}, {"disease": "Catecholaminergic polymorphic ventricular tachycardia", "inheritance": "AD"}],
    "CALM3": [{"disease": "Long-QT syndrome type 16", "inheritance": "AD"}, {"disease": "Catecholaminergic polymorphic ventricular tachycardia", "inheritance": "AD"}],
    "CASQ2": [{"disease": "Catecholaminergic polymorphic ventricular tachycardia", "inheritance": "AR"}],
    "COL3A1": [{"disease": "Ehlers-Danlos syndrome, vascular type", "inheritance": "AD"}],
    "CYP27A1": [{"disease": "Cerebrotendinous xanthomatosis", "inheritance": "AR"}],
    "DES": [{"disease": "Dilated cardiomyopathy", "inheritance": "AD"}, {"disease": "Myofibrillar myopathy", "inheritance": "AD"}],
    "DSC2": [{"disease": "Arrhythmogenic right ventricular cardiomyopathy", "inheritance": "AD"}],
    "DSG2": [{"disease": "Arrhythmogenic right ventricular cardiomyopathy", "inheritance": "AD"}],
    "DSP": [{"disease": "Arrhythmogenic right ventricular cardiomyopathy", "inheritance": "AD"}, {"disease": "Dilated cardiomyopathy", "inheritance": "AD"}],
    "ENG": [{"disease": "Hereditary hemorrhagic telangiectasia", "inheritance": "AD"}],
    "FBN1": [{"disease": "Marfan syndrome", "inheritance": "AD"}],
    "FLNC": [{"disease": "Dilated cardiomyopathy", "inheritance": "AD"}, {"disease": "Hypertrophic cardiomyopathy", "inheritance": "AD"}, {"disease": "Myofibrillar myopathy", "inheritance": "AD"}],
    "GAA": [{"disease": "Pompe disease", "inheritance": "AR"}],
    "GLA": [{"disease": "Fabry disease", "inheritance": "XL"}],
    "HFE": [{"disease": "Hereditary hemochromatosis (c.845G>A; p.C282Y homozygotes only)", "inheritance": "AR"}],
    "HNF1A": [{"disease": "Maturity-Onset of Diabetes of the Young", "inheritance": "AD"}],
    "KCNH2": [{"disease": "Long-QT syndrome type 2", "inheritance": "AD"}],
    "KCNQ1": [{"disease": "Long-QT syndrome type 1", "inheritance": "AD"}],
    "LDLR": [{"disease": "Familial hypercholesterolemia", "inheritance": "AD"}],
    "LMNA": [{"disease": "Dilated cardiomyopathy", "inheritance": "AD"}],
    "MAX": [{"disease": "Hereditary paraganglioma-pheochromocytoma syndrome", "inheritance": "AD"}],
    "MEN1": [{"disease": "Multiple endocrine neoplasia type 1", "inheritance": "AD"}],
    "MLH1": [{"disease": "Lynch syndrome", "inheritance": "AD"}],
    "MSH2": [{"disease": "Lynch syndrome", "inheritance": "AD"}],
    "MSH6": [{"disease": "Lynch syndrome", "inheritance": "AD"}],
    "MUTYH": [{"disease": "MUTYH-associated polyposis", "inheritance": "AR"}],
    "MYBPC3": [{"disease": "Hypertrophic cardiomyopathy", "inheritance": "AD"}],
    "MYH11": [{"disease": "Familial thoracic aortic aneurysm", "inheritance": "AD"}],
    "MYH7": [{"disease": "Hypertrophic cardiomyopathy", "inheritance": "AD"}, {"disease": "Dilated cardiomyopathy", "inheritance": "AD"}],
    "MYL2": [{"disease": "Hypertrophic cardiomyopathy", "inheritance": "AD"}],
    "MYL3": [{"disease": "Hypertrophic cardiomyopathy", "inheritance": "AD"}],
    "NF2": [{"disease": "NF2-related schwannomatosis", "inheritance": "AD"}],
    "OTC": [{"disease": "Ornithine transcarbamylase deficiency", "inheritance": "XL"}],
    "PALB2": [{"disease": "Hereditary breast cancer", "inheritance": "AD"}],
    "PCSK9": [{"disease": "Familial hypercholesterolemia", "inheritance": "AD"}],
    "PKP2": [{"disease": "Arrhythmogenic right ventricular cardiomyopathy", "inheritance": "AD"}],
    "PLN": [{"disease": "Dilated cardiomyopathy", "inheritance": "AD"}],
    "PMS2": [{"disease": "Lynch syndrome", "inheritance": "AD"}],
    "PRKAG2": [{"disease": "Hypertrophic cardiomyopathy", "inheritance": "AD"}],
    "PTEN": [{"disease": "PTEN hamartoma tumor syndrome", "inheritance": "AD"}],
    "RB1": [{"disease": "Retinoblastoma", "inheritance": "AD"}],
    "RBM20": [{"disease": "Dilated cardiomyopathy", "inheritance": "AD"}],
    "RET": [{"disease": "Familial medullary thyroid cancer", "inheritance": "AD"}, {"disease": "Multiple endocrine neoplasia type 2A", "inheritance": "AD"}, {"disease": "Multiple endocrine neoplasia type 2B", "inheritance": "AD"}],
    "RPE65": [{"disease": "RPE65-related retinopathy", "inheritance": "AR"}],
    "RYR1": [{"disease": "Malignant hyperthermia", "inheritance": "AD"}],
    "RYR2": [{"disease": "Catecholaminergic polymorphic ventricular tachycardia", "inheritance": "AD"}],
    "SCN5A": [{"disease": "Long QT syndrome type 3", "inheritance": "AD"}, {"disease": "Brugada syndrome", "inheritance": "AD"}, {"disease": "Dilated cardiomyopathy", "inheritance": "AD"}],
    "SDHAF2": [{"disease": "Hereditary paraganglioma-pheochromocytoma syndrome", "inheritance": "AD"}],
    "SDHB": [{"disease": "Hereditary paraganglioma-pheochromocytoma syndrome", "inheritance": "AD"}],
    "SDHC": [{"disease": "Hereditary paraganglioma-pheochromocytoma syndrome", "inheritance": "AD"}],
    "SDHD": [{"disease": "Hereditary paraganglioma-pheochromocytoma syndrome", "inheritance": "AD"}],
    "SMAD3": [{"disease": "Loeys-Dietz syndrome", "inheritance": "AD"}],
    "SMAD4": [{"disease": "Juvenile polyposis syndrome", "inheritance": "AD"}, {"disease": "Hereditary hemorrhagic telangiectasia", "inheritance": "AD"}],
    "STK11": [{"disease": "Peutz-Jeghers syndrome", "inheritance": "AD"}],
    "TGFBR1": [{"disease": "Loeys-Dietz syndrome", "inheritance": "AD"}],
    "TGFBR2": [{"disease": "Loeys-Dietz syndrome", "inheritance": "AD"}],
    "TMEM127": [{"disease": "Hereditary paraganglioma-pheochromocytoma syndrome", "inheritance": "AD"}],
    "TMEM43": [{"disease": "Arrhythmogenic right ventricular cardiomyopathy", "inheritance": "AD"}],
    "TNNC1": [{"disease": "Dilated cardiomyopathy", "inheritance": "AD"}],
    "TNNI3": [{"disease": "Hypertrophic cardiomyopathy", "inheritance": "AD"}],
    "TNNT2": [{"disease": "Dilated cardiomyopathy", "inheritance": "AD"}, {"disease": "Hypertrophic cardiomyopathy", "inheritance": "AD"}],
    "TP53": [{"disease": "Li-Fraumeni syndrome", "inheritance": "AD"}],
    "TPM1": [{"disease": "Hypertrophic cardiomyopathy", "inheritance": "AD"}],
    "TRDN": [{"disease": "Catecholaminergic polymorphic ventricular tachycardia", "inheritance": "AR"}, {"disease": "Long QT syndrome", "inheritance": "AR"}],
    "TSC1": [{"disease": "Tuberous sclerosis complex", "inheritance": "AD"}],
    "TSC2": [{"disease": "Tuberous sclerosis complex", "inheritance": "AD"}],
    "TTN": [{"disease": "Dilated cardiomyopathy (truncating variants only)", "inheritance": "AD"}],
    "TTR": [{"disease": "Hereditary transthyretin-related amyloidosis", "inheritance": "AD"}],
    "VHL": [{"disease": "Von Hippel-Lindau syndrome", "inheritance": "AD"}],
    "WT1": [{"disease": "WT1-related Wilms tumor", "inheritance": "AD"}],
}

_PGX_REPORTABLE_STRENGTHS = {"strong", "moderate", "optional"}

_PGX_LEVEL_A_GENE_INFO = {
    "CYP2D6": ("StellarPGx（WGS）", "抗憂鬱劑、止痛藥、抗精神病藥"),
    "CYP2C19": ("PharmCAT（VCF）", "clopidogrel、PPI、抗憂鬱劑"),
    "CYP2C9": ("PharmCAT（VCF）", "warfarin、NSAIDs、phenytoin"),
    "DPYD": ("PharmCAT（VCF）", "fluoropyrimidine（5-FU）毒性"),
    "TPMT": ("PharmCAT（VCF）", "thiopurine（azathioprine）毒性"),
    "NUDT15": ("PharmCAT（VCF）", "thiopurine（azathioprine）毒性（亞洲族群重要）"),
    "SLCO1B1": ("PharmCAT（VCF）", "simvastatin 肌肉毒性"),
    "HLA-A": ("OptiType（WGS）", "abacavir 過敏（*31:01 positive）"),
    "HLA-B": ("OptiType（WGS）", "abacavir（*57:01）、carbamazepine（*15:02）、allopurinol（*58:01）過敏"),
    "UGT1A1": ("PharmCAT（VCF）", "irinotecan 毒性"),
    "G6PD": ("PharmCAT（VCF）", "多種藥物溶血風險"),
    "MT-RNR1": ("mito pipeline + BAM mpileup", "aminoglycoside 致聾風險（CPIC Level A）"),
}


def _health_variant_gene(v: dict, edits: dict) -> str:
    return _snv_tx_field(v, edits, "gene_symbol") or v.get("gene_symbol") or v.get("GENE") or ""


def _health_variant_hgvs(v: dict, edits: dict) -> str:
    hgvs_c = _strip_hgvs_prefix(_snv_tx_field(v, edits, "HGVS_C", "hgvs_c"))
    hgvs_p = _strip_hgvs_prefix(_snv_tx_field(v, edits, "HGVS_P", "hgvs_p"))
    return hgvs_c + (f" ({hgvs_p})" if hgvs_p else "")


def _is_health_clinvar_plp(v: dict) -> bool:
    return _clinvar_label(v).lower().replace("_", " ") in {
        "pathogenic",
        "likely pathogenic",
        "pathogenic/likely pathogenic",
        "likely pathogenic/pathogenic",
    }


def _health_selected_ids(report: dict, category: str, candidate_ids: list[str], variants: dict) -> list[str]:
    section = (report.get("secondary_findings") or {}).get(category) or {}
    selected = {str(x) for x in section.get("selected") or []}
    dismissed = {str(x) for x in section.get("dismissed") or []}
    legacy = report.get("panels") or {}
    out = []
    for vid in candidate_ids:
        if vid in dismissed or (legacy.get(vid, {}) or {}).get(category) == "0":
            continue
        if vid in selected or (legacy.get(vid, {}) or {}).get(category) == "V" or _is_health_clinvar_plp(variants.get(vid, {})):
            out.append(vid)
    return out


def _acmg_disease_text(gene: str) -> str:
    rows = _ACMG_SF_DISEASES.get(gene) or []
    parts: list[str] = []
    for row in rows:
        disease = row.get("disease") or ""
        if disease:
            parts.append(disease)
    return "；".join(parts)


def _acmg_inheritance_text(gene: str) -> str:
    rows = _ACMG_SF_DISEASES.get(gene) or []
    seen: set[str] = set()
    out: list[str] = []
    for row in rows:
        inh = row.get("inheritance") or ""
        zh = _inheritance_zh(inh)
        if zh and zh not in seen:
            seen.add(zh)
            out.append(zh)
    return "、".join(out)


def _health_snv_variant_block(doc, v: dict, *, tier: str, edits: dict,
                              disease_text: str = "",
                              inheritance_text: str = "") -> None:
    gene = _health_variant_gene(v, edits) or "?"
    tx = _snv_transcript_label(v, edits)
    _add_paragraph(doc, f"    {gene} ({tx})", bold=True)
    rs = v.get("rs_id") or v.get("RS_ID") or ""
    struc = _structure_label({
        **v,
        "exon": _snv_tx_field(v, edits, "exon"),
        "intron": _snv_tx_field(v, edits, "intron"),
    })
    hgvs_c = _strip_hgvs_prefix(_snv_tx_field(v, edits, "HGVS_C", "hgvs_c"))
    hgvs_p = _strip_hgvs_prefix(_snv_tx_field(v, edits, "HGVS_P", "hgvs_p"))
    nuc = hgvs_c + (f"({hgvs_p})" if hgvs_p else "")
    zyg = _zygosity_long(v.get("zygosity", ""))
    clnsig = _clinvar_label(v)
    acmg = _acmg_label(v, edits)

    _ascii_table(doc, columns=[
        ("類別",          5),
        ("基因",          9),
        ("RS ID",         7, "buffered"),
        ("結構",          9),
        ("核苷酸",       14, "hgvs"),
        ("基因型",       13),
        ("ClinVar",      16, "token-buffered"),
        ("ACMG&AMP指引", 12, "token"),
    ], rows=[[tier, gene, rs, struc, nuc, zyg, clnsig, acmg]])

    if disease_text:
        inh_clause = f"，其遺傳模式屬於{inheritance_text}" if inheritance_text else ""
        _add_paragraph(doc, f"    1. {gene} 為 {disease_text} 之相關基因{inh_clause}。")
    else:
        _add_paragraph(doc, f"    1. {_omim_block_for_snv(v, edits)}")
    _add_paragraph(doc, f"    2. {_patho_sentence(acmg)}")


def _render_health_secondary_section(doc, title: str, ids: list[str], variants: dict, report: dict) -> None:
    _add_paragraph(doc, title, bold=True)
    if not ids:
        _add_paragraph(doc, f"  {_NO_HEALTH_VARIANT_TEXT}")
    else:
        edits = report.get("edits") or {}
        for idx, vid in enumerate(ids, start=1):
            v = variants.get(vid) or {}
            _health_snv_variant_block(doc, v, tier=str(idx), edits=edits.get(vid) or {})
            _blank(doc)
    _blank(doc)


def _render_health_acmg_section(doc, title: str, ids: list[str], variants: dict, report: dict) -> None:
    _add_paragraph(doc, title, bold=True)
    remaining = set(ids)
    for idx, group in enumerate(_ACMG_SF_GROUPS, start=1):
        group_genes = set(group["genes"])
        group_ids = [
            vid for vid in ids
            if _health_variant_gene(variants.get(vid, {}), (report.get("edits") or {}).get(vid) or {}) in group_genes
        ]
        remaining -= set(group_ids)
        genes_label = ", ".join(group["genes"])
        if group["key"] == "cardiovascular":
            _add_paragraph(doc, f"  {idx}. {group['title']}")
        else:
            _add_paragraph(doc, f"  {idx}. {group['title']}，包含 {genes_label}")
        if group_ids:
            edits = report.get("edits") or {}
            for vid_idx, vid in enumerate(group_ids, start=1):
                v = variants.get(vid) or {}
                gene = _health_variant_gene(v, edits.get(vid) or {})
                _health_snv_variant_block(
                    doc,
                    v,
                    tier=str(vid_idx),
                    edits=edits.get(vid) or {},
                    disease_text=_acmg_disease_text(gene),
                    inheritance_text=_acmg_inheritance_text(gene),
                )
                _blank(doc)
        else:
            _add_paragraph(doc, f"    {_NO_HEALTH_VARIANT_TEXT}")
    if remaining:
        _add_paragraph(doc, "  其它未分類 ACMG SF 基因")
        edits = report.get("edits") or {}
        for vid_idx, vid in enumerate(sorted(remaining), start=1):
            v = variants.get(vid) or {}
            gene = _health_variant_gene(v, edits.get(vid) or {})
            _health_snv_variant_block(
                doc,
                v,
                tier=str(vid_idx),
                edits=edits.get(vid) or {},
                disease_text=_acmg_disease_text(gene),
                inheritance_text=_acmg_inheritance_text(gene),
            )
            _blank(doc)
    _blank(doc)


def _pgx_actionable_groups(pgx: dict) -> list[dict]:
    groups = []
    genes = pgx.get("genes") or {}
    for gene in pgx.get("actionable") or []:
        g = genes.get(gene) or {}
        diplotype = g.get("diplotype") or (g.get("details") or {}).get("label") or ""
        phenotype = g.get("phenotype") or g.get("mtrn1_risk") or ""
        recs = []
        for drug in g.get("drugs") or []:
            for rec in drug.get("recommendations") or []:
                if "CPIC" not in (rec.get("source") or "").upper():
                    continue
                level = rec.get("cpic_level") or rec.get("evidence") or ""
                if level.strip().lower() not in _PGX_REPORTABLE_STRENGTHS:
                    continue
                recs.append({
                    "drug": drug.get("drug") or "",
                    "level": level,
                    "recommendation": rec.get("recommendation") or "",
                })
        if recs:
            groups.append({
                "gene": gene,
                "diplotype": diplotype,
                "phenotype": phenotype,
                "recommendations": recs,
            })
    return groups


def _render_health_pgx_section(doc, title: str, pgx: dict) -> None:
    _add_paragraph(doc, title, bold=True)
    groups = _pgx_actionable_groups(pgx or {})
    if not groups:
        _add_paragraph(doc, "  本次藥物基因體學分析未檢出具 CPIC 用藥建議之 actionable 結果。")
    else:
        for group in groups:
            _ascii_table(doc, columns=[
                ("基因", 12),
                ("基因型", 28, "buffered"),
                ("表型", 44, "buffered"),
            ], rows=[[group["gene"], group["diplotype"], group["phenotype"]]])
            _add_paragraph(doc, "    相關用藥建議：")
            for idx, rec in enumerate(group["recommendations"], start=1):
                level = f" (CPIC Recommendation strength: {rec['level']})" if rec.get("level") else ""
                _add_paragraph(doc, f"    {idx}. {rec['drug']}: {rec['recommendation']}{level}")
            _blank(doc)
    _blank(doc)


def _health_panel_gene_sections(requested_set: set[str]) -> list[tuple[str, list[str]]]:
    out: list[tuple[str, list[str]]] = []
    panel_keys = {
        "acmg_sf": "ACMG_SF_v3.3",
        "lipid_fh": "lipid_fh",
        "hereditary_cancer": "WES-I__腫瘤醫學__遺傳癌症",
        "stroke": "WGS__神經科__Stroke",
        "carrier": "carrier_mackenzie_1300+",
        "proactive": "proactive",
    }
    titles = dict(_HEALTH_SECTION_ORDER)
    titles["acmg_sf"] = _HEALTH_ACMG_GENE_LIST_TITLE
    for key, panel_name in panel_keys.items():
        if key not in requested_set:
            continue
        payload = phenotype_scorer.genes_for_key(panel_name, kind="panel")
        genes = sorted({
            panel_deadzone.canonical_panel_gene_symbol(g)
            for g in payload.get("genes", [])
            if g
        })
        out.append((titles.get(key, key), [g for g in genes if g]))
    return out


def _pgx_level_a_gene_rows(pgx: dict) -> list[list[str]]:
    genes = pgx.get("genes") or {}
    tsv_genes = [g for g in (pgx.get("gene_order") or []) if g in genes and not (genes.get(g) or {}).get("additional")]
    rows: list[list[str]] = []
    for gene in _PGX_LEVEL_A_GENE_INFO:
        if gene not in tsv_genes:
            continue
        caller, meaning = _PGX_LEVEL_A_GENE_INFO[gene]
        rows.append([gene, caller, meaning])
    return rows


def _section_health_annotations(doc, requested_set: set[str], pgx: dict | None = None) -> None:
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
    sections = _health_panel_gene_sections(requested_set)
    if not sections:
        _add_paragraph(doc, "    （未選擇檢測基因項目）")
        return
    for idx, (name, genes) in enumerate(sections):
        if idx:
            _blank(doc)
        _add_paragraph(doc, f"{name}:")
        if genes:
            _add_paragraph(doc, ", ".join(genes))
        else:
            _add_paragraph(doc, "依藥物基因體學分析模組可判讀基因輸出。")
    if "pgx" in requested_set:
        if sections:
            _blank(doc)
        _add_paragraph(doc, "藥物基因體學:")
        rows = _pgx_level_a_gene_rows(pgx or {})
        if rows:
            _add_paragraph(doc, "CPIC Level A 基因清單")
            _ascii_table(doc, columns=[
                ("基因", 10),
                ("Outside caller", 24, "buffered"),
                ("主要臨床意義", 46, "buffered"),
            ], rows=rows)
        else:
            _add_paragraph(doc, "本次藥物基因體學分析未輸出 CPIC Level A 基因。")


def build_health_docx(sample_id: str, *, sections: Iterable[str] | None = None) -> bytes:
    sample = sample_loader.load_sample(sample_id, include_aux=False)
    if sample is None:
        raise FileNotFoundError(f"sample not found: {sample_id}")
    secondary = sample_loader.load_sample_secondary_snv(sample_id) or {}
    pgx_payload = sample_loader.load_sample_pgx(sample_id) or {}

    requested = [str(s).strip() for s in (sections or []) if str(s).strip()]
    if not requested:
        requested = ["acmg_sf", "pgx"]
    requested_set = set(requested)
    report = report_store.load(sample_id)
    variants = secondary.get("variants") or {}
    categories = secondary.get("categories") or {}

    doc = Document()
    _apply_normal_font(doc)
    _apply_page_margins(doc)

    meta = sample.get("meta") or {}
    test_type = meta.get("Test", "") or "WES"

    _section_health_patient_header(doc, meta)
    _section_test_info(doc, test_type, health=True)
    _add_paragraph(doc, "二、檢驗套組: 重大疾病基因篩檢")
    _blank(doc)

    _add_paragraph(doc, "三、檢測結果")
    _add_paragraph(doc, "  檢體說明:")
    _add_paragraph(doc, "    檢體類別：血液")
    _add_paragraph(doc, "  綜合說明:")

    referenced: list[dict] = []
    referenced_ids: set[str] = set()
    for key, title in _HEALTH_SECTION_ORDER:
        if key not in requested_set:
            continue
        if key == "pgx":
            _render_health_pgx_section(doc, title, pgx_payload.get("pgx") or pgx_payload.get("pharmcat") or {})
            continue
        ids = _health_selected_ids(report, key, categories.get(key) or [], variants)
        for vid in ids:
            if vid in referenced_ids:
                continue
            referenced_ids.add(vid)
            referenced.append(variants.get(vid) or {})
        if key == "acmg_sf":
            _render_health_acmg_section(doc, title, ids, variants, report)
        else:
            _render_health_secondary_section(doc, title, ids, variants, report)

    _add_paragraph(
        doc,
        "    建議比對個人臨床資料與家族病史，並由醫師或遺傳諮詢人員進行綜合判斷；"
        "根據家族成員變異位點檢測報告或相關資料庫更新，可能影響變異位點ACMG判讀結果。",
    )
    _blank(doc)
    if referenced:
        _add_paragraph(doc, "  參考資料:")
        edits = report.get("edits") or {}
        for v in referenced:
            if not v:
                continue
            _add_paragraph(doc, _snv_reference_text(v, edits.get(v.get("id", ""), {})))
            _blank(doc)

    _section_methods(doc, test_type)
    _section_health_annotations(doc, requested_set, pgx_payload.get("pgx") or pgx_payload.get("pharmcat") or {})

    buf = io.BytesIO()
    doc.save(buf)
    return buf.getvalue()


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
