"""Rebuild reviewer-confirmed CNV/SV parent records from original segments."""
from __future__ import annotations

from copy import deepcopy


def _merge_id(source: str, chrom: str, start: int, end: int, sv_type: str) -> str:
    return f"MERGED-{source.upper()}-{chrom}-{start}-{end}-{sv_type.upper()}"


def _dedupe_genes(segments: list[dict]) -> list[dict]:
    out: list[dict] = []
    seen: set[str] = set()
    for seg in segments:
        for gene in (seg.get("genes") or []) + (seg.get("genes_overflow") or []):
            name = (gene.get("gene") or "").strip() if isinstance(gene, dict) else ""
            if not name or name in seen:
                continue
            seen.add(name)
            out.append(gene)
    return out


def build_parent(merge: dict, variants: dict[str, dict]) -> dict | None:
    member_ids = merge.get("member_ids") or []
    segments = [variants[mid] for mid in member_ids if mid in variants]
    if len(segments) < 2:
        return None
    source = str(merge.get("source") or segments[0].get("source") or "cnv").lower()
    chrom = str(segments[0].get("CHROM") or "")
    sv_type = str(segments[0].get("sv_type") or "").upper()
    if not chrom or sv_type not in ("DEL", "DUP"):
        return None
    if any(str(v.get("CHROM") or "") != chrom or
           str(v.get("sv_type") or "").upper() != sv_type for v in segments):
        return None
    starts = [int(v["POS"]) for v in segments if v.get("POS") is not None]
    ends = [int(v["END"]) for v in segments if v.get("END") is not None]
    if len(starts) != len(segments) or len(ends) != len(segments):
        return None
    start, end = min(starts), max(ends)
    parent = deepcopy(max(segments, key=lambda v: v.get("ranking_score") or -999))
    genes = _dedupe_genes(segments)
    parent.update({
        "id": str(merge.get("id") or _merge_id(source, chrom, start, end, sv_type)),
        "source": source,
        "CHROM": chrom,
        "POS": start,
        "END": end,
        "length": end - start,
        "gene_count": len(genes),
        "genes": genes,
        "genes_overflow": [],
        "genes_compact": [],
        "genes_total": len(genes),
        "merged_segment_ids": list(member_ids),
        "is_merged_parent": True,
    })
    return parent


def apply_confirmed_merges(
    cnv_variants: dict[str, dict],
    sv_variants: dict[str, dict],
    merges: list[dict] | None,
) -> tuple[dict[str, dict], dict[str, dict]]:
    """Return copies where confirmed parents replace their member segments."""
    cnv_out = dict(cnv_variants or {})
    sv_out = dict(sv_variants or {})
    for merge in merges or []:
        if not isinstance(merge, dict):
            continue
        source = str(merge.get("source") or "cnv").lower()
        target = cnv_out if source == "cnv" else sv_out
        parent = build_parent(merge, target)
        if parent:
            for member_id in merge.get("member_ids") or []:
                target.pop(member_id, None)
            target[parent["id"]] = parent
    return cnv_out, sv_out
