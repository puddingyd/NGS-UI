#!/usr/bin/env python3
"""Render the health report's ACMG education without loading patient data."""
from __future__ import annotations

import argparse
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "backend"))

from docx import Document
from app.services import docx_export
from app.services.acmg_sf_education import load_catalogue, render_acmg_sf_education


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    doc = Document()
    docx_export._apply_normal_font(doc)
    docx_export._apply_page_margins(doc)
    docx_export._add_paragraph(doc, "附錄", bold=True, align="center")
    docx_export._blank(doc)
    render_acmg_sf_education(doc, font_name=docx_export.REPORT_FONT)
    doc.core_properties.title = load_catalogue()["title"]
    doc.core_properties.author = ""
    doc.core_properties.last_modified_by = ""
    args.output.parent.mkdir(parents=True, exist_ok=True)
    doc.save(args.output)
    print(args.output)


if __name__ == "__main__":
    main()
