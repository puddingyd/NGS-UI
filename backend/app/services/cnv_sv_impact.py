"""Display-only impact summaries; never change CNV/SV eligibility or ACMG.

Use individual AnnotSV split rows before trimming/merging. Location2=CDS
includes introns and must never, by itself, imply coding-exon overlap.
"""
from __future__ import annotations

import math
import re


def number(value):
    try:
        result = float(value)
        return result if math.isfinite(result) else None
    except (ValueError, TypeError):
        return None


def gene_impact(gene: dict, sv_type: str) -> tuple[str, str]:
    location = str(gene.get("location") or "").lower()
    coding = str(gene.get("location2") or "").lower()
    cds = number(gene.get("overlap_cds_len"))
    distance = number(gene.get("splice_distance"))
    # Balanced/complex SV spans do not necessarily disrupt the enclosed genes.
    # Retain them as unknown until breakpoint-specific effects are available.
    if sv_type.upper() not in {"DEL", "DUP"}:
        return "unknown", "斷點影響待釐清"
    if cds is not None and cds > 0:
        if location in {"txstart-txend", "txend-txstart"}:
            return "functional", "完整基因"
        return "functional", "編碼外顯子"
    splice_type = str(gene.get("splice_type") or "").strip().lower()
    if distance is not None and 0 <= distance <= 2 and splice_type not in {"", ".", "na", "nan"}:
        return "functional", "剪接位置（距離 ≤2 bp）"
    # Only call noncoding when the annotation explicitly establishes it.
    # Missing splice distance cannot establish that an intronic call is remote.
    if cds == 0 and distance is not None and distance > 2:
        introns = re.fullmatch(r"intron(\d+)-intron(\d+)", location)
        if introns and introns[1] == introns[2]:
            return "noncoding", "純內含子"
        if "utr" in coding and "cds" not in coding:
            return "noncoding", "僅 UTR／非編碼轉錄本"
    return "unknown", "位置註解不足"


def summarize(genes: list[dict], sv_type: str) -> dict:
    if not genes:
        return {"category": "unknown", "reasons": [], "hpo_score": None, "mechanism": 0}
    classified = [(gene, *gene_impact(gene, sv_type)) for gene in genes]
    # Unknown alongside a known noncoding hit must not cause automatic hiding.
    category = next(c for c in ("functional", "unknown", "noncoding")
                    if any(item[1] == c for item in classified))
    relevant = [(g, reason) for g, c, reason in classified if c == category]
    scores = [number(g.get("hpo_score")) for g, _ in relevant]
    mechanism = 0
    for g, reason in relevant:
        if category != "functional":
            continue
        if sv_type.upper() == "DEL" and number(g.get("hi")) == 3:
            mechanism = 1
        elif sv_type.upper() == "DUP" and reason == "完整基因" and number(g.get("ts")) == 3:
            mechanism = 1
    return {
        "category": category,
        "reasons": [{"gene": g.get("gene", ""), "impact": reason,
                     "hpo_score": g.get("hpo_score")} for g, reason in relevant[:3]],
        "hpo_score": max((s for s in scores if s is not None), default=None),
        "mechanism": mechanism,
    }


def attach(variant: dict) -> None:
    genes = variant.get("genes") or []
    matched = [g for g in genes if (number(g.get("pheno_score")) or 0) > 0]
    variant["impact_clinical"] = summarize(matched, variant.get("sv_type", ""))
    variant["impact_all"] = summarize(genes, variant.get("sv_type", ""))


def merge_summaries(summaries: list[dict]) -> dict:
    """Combine actual segments, never annotate the union's intervening gaps."""
    if not summaries:
        return summarize([], "")
    category = next(c for c in ("functional", "unknown", "noncoding")
                    if any(s.get("category") == c for s in summaries))
    chosen = [s for s in summaries if s.get("category") == category]
    scores = [number(s.get("hpo_score")) for s in chosen]
    return {"category": category,
            "reasons": [r for s in chosen for r in s.get("reasons", [])][:3],
            "hpo_score": max((s for s in scores if s is not None), default=None),
            "mechanism": max(s.get("mechanism", 0) for s in chosen)}
