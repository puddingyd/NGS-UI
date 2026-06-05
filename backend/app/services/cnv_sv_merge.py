"""Build parent CNV/SV records by joining compatible nearby segments."""
from __future__ import annotations

from copy import deepcopy

MERGE_GAP_THRESHOLD = 250_000


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


def _compatible(a: dict, b: dict) -> bool:
    if str(a.get("CHROM") or "") != str(b.get("CHROM") or ""):
        return False
    if str(a.get("sv_type") or "").upper() != str(b.get("sv_type") or "").upper():
        return False
    if str(a.get("sv_type") or "").upper() not in ("DEL", "DUP"):
        return False
    acn, bcn = a.get("copy_number"), b.get("copy_number")
    return acn is None or bcn is None or str(acn) == str(bcn)


def automatic_merges(variants: dict[str, dict], source: str) -> list[dict]:
    """Group same-type adjacent DEL/DUP segments using a 250 kb max gap."""
    sorted_vars = sorted(
        variants.values(),
        key=lambda v: (
            str(v.get("CHROM") or ""),
            str(v.get("sv_type") or "").upper(),
            int(v.get("POS") or 0),
        ),
    )
    groups: list[list[dict]] = []
    current: list[dict] = []
    for variant in sorted_vars:
        previous = current[-1] if current else None
        gap = (
            int(variant.get("POS") or 0) - int(previous.get("END") or 0)
            if previous else MERGE_GAP_THRESHOLD + 1
        )
        if previous and _compatible(previous, variant) and 0 <= gap <= MERGE_GAP_THRESHOLD:
            current.append(variant)
        else:
            if len(current) >= 2:
                groups.append(current)
            current = [variant]
    if len(current) >= 2:
        groups.append(current)
    return [{
        "id": _merge_id(
            source,
            str(group[0].get("CHROM") or ""),
            min(int(v["POS"]) for v in group),
            max(int(v["END"]) for v in group),
            str(group[0].get("sv_type") or ""),
        ),
        "source": source,
        "member_ids": [str(v["id"]) for v in group],
    } for group in groups]


def apply_confirmed_merges(
    cnv_variants: dict[str, dict],
    sv_variants: dict[str, dict],
    merges: list[dict] | None,
) -> tuple[dict[str, dict], dict[str, dict]]:
    """Return copies where automatic parents replace nearby member segments."""
    cnv_out = dict(cnv_variants or {})
    sv_out = dict(sv_variants or {})
    automatic = automatic_merges(cnv_out, "cnv") + automatic_merges(sv_out, "sv")
    # Keep saved definitions readable for data created by the earlier
    # reviewer-confirmation version, while automatic grouping is now default.
    all_merges = automatic + list(merges or [])
    accepted: dict[str, list[tuple[dict, dict]]] = {"cnv": [], "sv": []}
    consumed: set[str] = set()
    for merge in all_merges:
        if not isinstance(merge, dict):
            continue
        member_ids = merge.get("member_ids") or []
        if any(member_id in consumed for member_id in member_ids):
            continue
        source = str(merge.get("source") or "cnv").lower()
        target = cnv_out if source == "cnv" else sv_out
        parent = build_parent(merge, target)
        if parent:
            for member_id in member_ids:
                consumed.add(member_id)
            accepted[source].append((merge, parent))

    def replace_in_order(source: str, variants: dict[str, dict]) -> dict[str, dict]:
        parent_by_member: dict[str, dict] = {}
        first_members: set[str] = set()
        positions = {variant_id: index for index, variant_id in enumerate(variants)}
        for merge, parent in accepted[source]:
            members = [mid for mid in (merge.get("member_ids") or []) if mid in variants]
            if not members:
                continue
            first = min(members, key=positions.__getitem__)
            first_members.add(first)
            for member_id in members:
                parent_by_member[member_id] = parent
        out: dict[str, dict] = {}
        for vid, variant in variants.items():
            if vid in first_members:
                parent = parent_by_member[vid]
                out[parent["id"]] = parent
            elif vid not in parent_by_member:
                out[vid] = variant
        return out

    return replace_in_order("cnv", cnv_out), replace_in_order("sv", sv_out)
