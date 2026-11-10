"""Gene lookup for the /omim-capture/ page.

Shows exactly what variant cards and reports use for a gene (the rows
``omim_store`` serves, from the platform OMIM database or OMIM.xlsx), plus
database-only context: where each disease text came from, adoption and OMIM
edit dates, other genes sharing a disease, removed diseases, the gene's work
queue steps, open review items and recent changes.
"""
from __future__ import annotations

import re

from . import disease_picks, omim_db, omim_queue, omim_store

_MIM_RE = re.compile(r"^\d{6}$")


def _rows() -> list[dict]:
    omim_store.ensure_loaded()
    seen, rows = set(), []
    for row in list(omim_store._state.get("by_omim_id", {}).values()) + list(
        omim_store._state.get("by_gene", {}).values()
    ):
        if id(row) not in seen:
            seen.add(id(row))
            rows.append(row)
    return rows


def search(q: str, limit: int = 20) -> list[dict]:
    """Gene symbols (incl. aliases) starting with / containing ``q``."""
    q = str(q or "").strip().upper()
    if not q:
        return []
    omim_store.ensure_loaded()
    by_gene = omim_store._state.get("by_gene", {})
    starts, contains = [], []
    for symbol, row in by_gene.items():
        upper = symbol.upper()
        if upper.startswith(q):
            starts.append((upper != q, len(symbol), symbol, row))
        elif q in upper:
            contains.append((True, len(symbol), symbol, row))
    out = []
    for _exact, _n, symbol, row in sorted(starts) + sorted(contains):
        out.append({
            "symbol": symbol,
            "gene": row.get("gene_symbol") or symbol,
            "gene_mim": str(row.get("OMIM_id") or ""),
        })
        if len(out) >= limit:
            break
    return out


def resolve(q: str) -> list[dict]:
    """Rows matching a gene symbol/alias, a gene MIM or a phenotype MIM."""
    q = str(q or "").strip()
    if not q:
        return []
    omim_store.ensure_loaded()
    if _MIM_RE.match(q):
        row = omim_store._state.get("by_omim_id", {}).get(int(q))
        if row:
            return [row]
        tag = f"({q})"
        return [
            row for row in _rows()
            if tag in (row.get("OMIM_disease") or "")
            or any(tag in disease_picks.first_line(t) for _s, t in disease_picks.disease_cells(row))
        ]
    by_gene = omim_store._state.get("by_gene", {})
    row = by_gene.get(q)
    if row is None:
        upper = q.upper()
        row = next((r for sym, r in by_gene.items() if sym.upper() == upper), None)
    return [row] if row else []


def summary(row: dict) -> dict:
    return {
        "gene": row.get("gene_symbol") or "",
        "gene_mim": str(row.get("OMIM_id") or ""),
        "diseases": len(disease_picks.disease_cells(row)),
    }


def detail(row: dict) -> dict:
    gene_mim = str(row.get("OMIM_id") or "")
    extras = omim_db.gene_extras(gene_mim) if omim_db.exists() and gene_mim else {}
    texts = extras.get("texts", {})
    linked = extras.get("linked_genes", {})
    usage = omim_queue.reported_usage()
    meta = {m.get("slot"): m for m in row.get("Disease_meta") or []}
    diseases = []
    for slot, text in disease_picks.disease_cells(row):
        label = disease_picks.first_line(text)
        pm = disease_picks.phenotype_mim(label)
        lines = len([ln for ln in text.splitlines() if ln.strip()])
        source = (meta.get(slot) or {}).get("source") or ("xlsx" if lines > 1 else "label")
        info = texts.get(pm, {}) if source == "capture" else {}
        diseases.append({
            "slot": slot,
            "label": label,
            "phenotype_mim": pm,
            "url": f"https://omim.org/entry/{pm}" if pm else "",
            "text": text,
            "source": source,
            "adopted_at": info.get("adopted_at", ""),
            "adopted_by": info.get("adopted_by", ""),
            "omim_edit_date": info.get("omim_edit_date", ""),
            "description_source": info.get("description_source", ""),
            "synopsis_status": info.get("synopsis_status", ""),
            "linked_genes": sorted(set(linked.get(pm, []))),
            "reported": usage.get(pm, 0) if pm else 0,
        })
    queue_item = next(
        (g for g in omim_queue.build_genes() if g["gene_mim"] == gene_mim and gene_mim), None
    )
    return {
        "gene": row.get("gene_symbol") or "",
        "gene_mim": gene_mim,
        "url": f"https://omim.org/entry/{gene_mim}" if _MIM_RE.match(gene_mim) else "",
        "source": omim_store.source(),
        "aliases": extras.get("aliases", []),
        "inheritance": row.get("Inheritance") or "",
        "done": row.get("Done") or "",
        "legacy_done": extras.get("legacy_done", row.get("Legacy_done", "")),
        "synced_at": extras.get("synced_at", ""),
        "omim_disease": row.get("OMIM_disease") or "",
        "diseases": diseases,
        "removed": extras.get("removed", []),
        "changes": extras.get("changes", []),
        "queue": {
            "status": queue_item["status"],
            "steps": queue_item["steps"],
            "reviews": queue_item["reviews"],
        } if queue_item else None,
    }
