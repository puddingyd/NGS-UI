"""Shared CNV/SV report-disease selection helpers.

The reviewer can select OMIM workbook diseases and pathogenic-region overlap
diseases, then optionally add the legacy free-text Disease value.  Every
consumer uses the same ordered, de-duplicated union.
"""
from __future__ import annotations

from typing import Any


def selected_items(edits: dict | None) -> list[dict[str, str]]:
    raw = (edits or {}).get("report_disease_items") or {}
    values: list[Any]
    if isinstance(raw, dict):
        values = list(raw.values())
    elif isinstance(raw, list):
        values = raw
    else:
        values = []

    out: list[dict[str, str]] = []
    for value in values:
        if isinstance(value, str):
            label = value.strip()
            item = {"label": label}
        elif isinstance(value, dict):
            label = str(
                value.get("label")
                or value.get("display_name")
                or value.get("disease")
                or ""
            ).strip()
            item = {
                "label": label,
                "source": str(value.get("source") or "").strip(),
                "gene": str(value.get("gene") or "").strip(),
                "phenotype_mim": str(value.get("phenotype_mim") or "").strip(),
                "inheritance": str(value.get("inheritance") or "").strip(),
                "overlap_type": str(value.get("overlap_type") or "").strip(),
            }
        else:
            continue
        if label:
            out.append(item)
    return out


def disease_labels(edits: dict | None) -> list[str]:
    """Return selected labels followed by legacy free text, de-duplicated."""
    candidates = [item["label"] for item in selected_items(edits)]
    manual = str((edits or {}).get("disease") or "").strip()
    if manual:
        # The report itself uses the ideographic comma as its disease
        # separator, so treat a legacy value already written in that form as
        # individual entries.  Other punctuation remains untouched because it
        # may be part of the disease name.
        candidates.extend(part.strip() for part in manual.split("、"))

    out: list[str] = []
    seen: set[str] = set()
    for label in candidates:
        clean = " ".join(str(label or "").split())
        if not clean:
            continue
        key = clean.casefold()
        if key in seen:
            continue
        seen.add(key)
        out.append(clean)
    return out


def disease_text(edits: dict | None) -> str:
    return "、".join(disease_labels(edits))
