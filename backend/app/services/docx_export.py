"""Build a 細明體 (PMingLiU) DOCX diagnostic report from sample data.

Layout (v1, intentionally minimal — iterate as the form takes shape):

  Title              "NGS 三級分析 — 診斷報告"
  Patient block      LIS_ID / Name / MRN / Test / Build / Run date
  Phenotype          ・bullet list of HPO label (HP:id) [w]
  Causative variants Section heading + per-variant box for status='1'
  Other variants     Same template for status='2'

The font name is set to PMingLiU on every run; East-Asian text gets
the same via w:eastAsia. Word will substitute when opening on a
machine without 細明體 installed.
"""
from __future__ import annotations

import io
from datetime import datetime
from typing import Iterable

from docx import Document
from docx.oxml.ns import qn
from docx.shared import Pt

from ..config import TERTIARY_OUTPUT_ROOT
from . import report_store, sample_loader

REPORT_FONT     = "PMingLiU"   # 細明體
TITLE_FONT_SIZE = Pt(16)
H2_FONT_SIZE    = Pt(13)
BODY_FONT_SIZE  = Pt(11)


def _set_run_font(run, name: str = REPORT_FONT) -> None:
    """Set both Western and East-Asian fonts on a run.

    python-docx's `run.font.name` only sets the Western face. East-Asian
    Chinese/Japanese/Korean glyphs need <w:eastAsia> set explicitly via
    raw XML, otherwise Word picks the default theme font.
    """
    run.font.name = name
    rPr = run._element.get_or_add_rPr()
    rFonts = rPr.find(qn("w:rFonts"))
    if rFonts is None:
        rFonts = rPr.makeelement(qn("w:rFonts"), {})
        rPr.append(rFonts)
    rFonts.set(qn("w:eastAsia"), name)
    rFonts.set(qn("w:ascii"),    name)
    rFonts.set(qn("w:hAnsi"),    name)


def _add_paragraph(doc, text: str, bold: bool = False, size: Pt = BODY_FONT_SIZE,
                   align: str | None = None):
    p = doc.add_paragraph()
    if align == "center":
        p.alignment = 1
    run = p.add_run(text)
    run.bold = bold
    run.font.size = size
    _set_run_font(run)
    return p


def _add_heading(doc, text: str, size: Pt = H2_FONT_SIZE):
    return _add_paragraph(doc, text, bold=True, size=size)


def _kv_table(doc, rows: list[tuple[str, str]]):
    """Two-column key/value table with bordered cells."""
    if not rows:
        return
    table = doc.add_table(rows=len(rows), cols=2)
    table.autofit = True
    for i, (k, v) in enumerate(rows):
        c0 = table.cell(i, 0); c1 = table.cell(i, 1)
        for cell, txt, bold in ((c0, k, True), (c1, v or "—", False)):
            cell.text = ""
            run = cell.paragraphs[0].add_run(txt)
            run.bold = bold
            run.font.size = BODY_FONT_SIZE
            _set_run_font(run)


def _variant_block(doc, v: dict, edits: dict) -> None:
    """One variant per box. Pulls user edits over the upstream values."""
    gene = v.get("gene_symbol") or "?"
    hgvs = v.get("HGVS") or v.get("id") or ""
    _add_paragraph(doc, f"{gene}  {hgvs}", bold=True, size=Pt(12))

    acmg_class = edits.get("ACMG_classification") or v.get("ACMG_classification") or ""
    acmg_score = edits.get("ACMG_score") if "ACMG_score" in edits else v.get("ACMG_score")
    acmg_crit  = edits.get("ACMG_criteria") or v.get("ACMG_criteria") or ""
    comment    = edits.get("comment") or ""
    clnsig     = v.get("CLNSIG") or ""
    stars      = v.get("clinvar_stars")
    clinvar    = f"{clnsig}{f' ({stars}★)' if stars not in (None, '') else ''}" if clnsig else "—"
    af         = v.get("AF")
    af_eas     = v.get("AF_eas")
    twb        = v.get("TaiwanBioBank")

    rows = [
        ("Zygosity",       str(v.get("zygosity") or "—")),
        ("Consequence",    str(v.get("Consequence") or "—")),
        ("Phase",          str(v.get("phase_result") or "—")),
        ("ClinVar",        clinvar),
        ("ACMG class",     str(acmg_class or "—")),
        ("ACMG score",     "—" if acmg_score in (None, "") else str(acmg_score)),
        ("ACMG criteria",  str(acmg_crit or "—")),
        ("AF (gnomAD G)",  "—" if af in (None, "") else f"{af}"),
        ("AF EAS",         "—" if af_eas in (None, "") else f"{af_eas}"),
        ("TWB AF",         "—" if twb in (None, "") else f"{twb}"),
        ("In panel",       "Yes" if v.get("in_panel") else "No"),
        ("Pheno score",    "—" if v.get("pheno_score") in (None, "") else f"{v.get('pheno_score')}"),
    ]
    _kv_table(doc, rows)
    if comment:
        _add_paragraph(doc, f"Comment: {comment}")
    # Diseases the user picked for the report
    picked = (edits.get("report_diseases") or {}) if isinstance(edits.get("report_diseases"), dict) else {}
    for i in range(1, 6):
        if not picked.get(str(i)) and not picked.get(i):
            continue
        d = v.get(f"Disease{i}")
        if d:
            _add_paragraph(doc, f"Disease {i}: {d}")
    doc.add_paragraph()  # spacer


def build_diagnosis_docx(sample_id: str) -> bytes:
    sample = sample_loader.load_sample(sample_id)
    if sample is None:
        raise FileNotFoundError(f"sample not found: {sample_id}")
    report = report_store.load(sample_id)

    doc = Document()
    # Make Normal style use 細明體 too so anything we forget to wrap still works.
    normal = doc.styles["Normal"]
    normal.font.name = REPORT_FONT
    normal.font.size = BODY_FONT_SIZE
    rPr = normal.element.get_or_add_rPr()
    rFonts = rPr.find(qn("w:rFonts"))
    if rFonts is None:
        from docx.oxml import OxmlElement
        rFonts = OxmlElement("w:rFonts")
        rPr.append(rFonts)
    rFonts.set(qn("w:eastAsia"), REPORT_FONT)
    rFonts.set(qn("w:ascii"),    REPORT_FONT)
    rFonts.set(qn("w:hAnsi"),    REPORT_FONT)

    # Title
    _add_paragraph(doc, "NGS 三級分析 — 診斷報告", bold=True, size=TITLE_FONT_SIZE,
                   align="center")
    doc.add_paragraph()

    # Patient block
    m = sample["meta"]
    _add_heading(doc, "病人資訊")
    _kv_table(doc, [
        ("LIS_ID",    m.get("LIS_ID") or "—"),
        ("Name",      m.get("Name")   or "—"),
        ("MRN",       m.get("MRN")    or "—"),
        ("Test",      m.get("Test")   or "—"),
        ("Build",     sample.get("genome_build") or "—"),
        ("Category",  m.get("Category") or "—"),
        ("Generated", sample.get("generated_at") or "—"),
        ("Reported",  datetime.now().strftime("%Y-%m-%d %H:%M")),
    ])
    doc.add_paragraph()

    # Phenotype
    pheno = sample.get("patient_phenotype") or []
    if pheno:
        _add_heading(doc, "臨床表型 (HPO)")
        for r in pheno:
            label = r.get("label") or r.get("phenotype") or ""
            hid   = r.get("phenotype") or ""
            w     = r.get("weight")
            tag   = f"{label}  ({hid})" + (f"  [w={w}]" if w not in (None, "") else "")
            p = doc.add_paragraph(style="List Bullet")
            run = p.add_run(tag)
            run.font.size = BODY_FONT_SIZE
            _set_run_font(run)
        doc.add_paragraph()

    variants = sample["variants"]
    statuses = report.get("status", {}) or {}
    edits    = report.get("edits", {})  or {}
    manuals  = report.get("manual_variants", []) or []

    def _section(title: str, target_status: str, manual_status: str):
        ids = [vid for vid, s in statuses.items() if s == target_status and vid in variants]
        ms  = [m for m in manuals if (m or {}).get("status") == manual_status
                                  and (m.get("position") or "").strip()]
        if not ids and not ms:
            _add_heading(doc, title)
            _add_paragraph(doc, "（無）")
            doc.add_paragraph()
            return
        _add_heading(doc, title)
        for vid in ids:
            _variant_block(doc, variants[vid], edits.get(vid, {}))
        for m in ms:
            _add_paragraph(doc, f"{m.get('position', '')}", bold=True, size=Pt(12))
            if m.get("disease"):
                _add_paragraph(doc, f"Disease: {m['disease']}")
            if m.get("comment"):
                _add_paragraph(doc, f"Comment: {m['comment']}")
            doc.add_paragraph()

    _section("Causative variants",  "1", "1")
    _section("Other variants",      "2", "2")

    # Optional free-text comment
    comment = (report.get("comment") or "").strip()
    if comment:
        _add_heading(doc, "備註")
        _add_paragraph(doc, comment)

    buf = io.BytesIO()
    doc.save(buf)
    return buf.getvalue()
