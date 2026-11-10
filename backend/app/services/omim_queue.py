"""OMIM work queue, organised by gene.

Each queue item is one gene with ordered steps:

1. the OMIM *gene page* (``*``/``+`` entry) when it has never been adopted —
   this syncs the gene's disease list, so diseases OMIM lists but the
   platform lacks get added (and then show up as missing);
2. one *phenotype page* per disease of that gene that still needs text:
   ``missing`` (card shows only the name, the ``*``) or ``stale`` (curated
   text older than ``STALE_DAYS`` that was never refreshed from a capture).

Diseases without a phenotype MIM are handled by the gene step. Open review
items (``omim_db.review_items``) are attached to their gene.

Gene order: genes in progress (gene page adopted, but diseases still to
fill or items to confirm) → disease-associated (panel_loose) genes whose
gene page has not been synced yet → everything else; inside each group,
genes whose diseases were ticked in case reports come first, then genes
with review items, then genes with diseases to fill. Within a gene,
the gene page comes first, then diseases (missing before stale; entries in
{}, [] or with "?" last).

``next_step`` keeps the reviewer on the same gene until its steps are done,
skipping pages captured but not adopted yet, pages adopted in the last 7
days, and anything the extension asked to skip (``MIM`` or ``g<geneMIM>``).
"""
from __future__ import annotations

import json
import threading
import time
from datetime import date, timedelta

from . import disease_picks, omim_db, omim_store, panel_deadzone

STALE_DAYS = 365
KIND_RANK = {"missing": 2, "stale": 1}
STATUS_LABELS = {
    "in_progress": "進行中（基因頁已同步，尚有疾病待補或待確認）",
    "gene_sync": "基因頁待同步",
    "phenotype": "只有疾病待補",
}
_REPORT_TTL = 300

_lock = threading.Lock()
_cache: dict = {"sig": None, "genes": []}
_report_cache: dict = {"at": 0.0, "usage": {}}


def reported_usage(ttl: int = _REPORT_TTL) -> dict[str, int]:
    """phenotype MIM -> number of variant reports that ticked it."""
    now = time.monotonic()
    if _report_cache["usage"] and now - _report_cache["at"] < ttl:
        return _report_cache["usage"]
    from . import sample_layout

    usage: dict[str, int] = {}
    for sample_id in sample_layout.iter_sample_ids():
        path = sample_layout.state_file(sample_id, "sample_metadata.json")
        try:
            meta = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            continue
        edits = meta.get("edits") if isinstance(meta, dict) else None
        if not isinstance(edits, dict):
            continue
        for edit in edits.values():
            picks = edit.get(disease_picks.PICKS_FIELD) if isinstance(edit, dict) else None
            if not isinstance(picks, dict):
                continue
            for key, snapshot in picks.items():
                if not snapshot:
                    continue
                mim = str(key).split("|", 1)[0]
                if mim:
                    usage[mim] = usage.get(mim, 0) + 1
    _report_cache.update(at=now, usage=usage)
    return usage


def _rows() -> list[dict]:
    omim_store.ensure_loaded()
    seen = set()
    rows = []
    for row in list(omim_store._state.get("by_omim_id", {}).values()) + list(
        omim_store._state.get("by_gene", {}).values()
    ):
        if id(row) in seen:
            continue
        seen.add(id(row))
        rows.append(row)
    return rows


def _stale(done: str, today: date) -> bool:
    if not done:
        return True
    try:
        return date.fromisoformat(done[:10]) < today - timedelta(days=STALE_DAYS)
    except ValueError:
        return True


def _is_mim(value: str) -> bool:
    return len(value) == 6 and value.isdigit()


def _signature() -> tuple:
    return (omim_store.cache_signature(), omim_db.signature() if omim_db.exists() else None)


def build_genes() -> list[dict]:
    """Every gene with work to do (cached until OMIM data or the DB changes)."""
    sig = _signature()
    with _lock:
        if _cache["sig"] == sig:
            return _cache["genes"]
    today = date.today()
    panel = panel_deadzone.disease_associated_genes()
    synced = omim_db.gene_sync_map() if omim_db.exists() else {}
    reviews = omim_db.open_reviews_by_gene() if omim_db.exists() else {}
    genes = []
    for row in _rows():
        gene = row.get("gene_symbol") or ""
        gene_mim = str(row.get("OMIM_id") or "").split(".")[0]
        done = row.get("Done") or ""
        # Curated workbook text is judged by the workbook's own Done date.
        curated_done = row.get("Legacy_done", done) or ""
        meta = {m.get("slot"): m for m in row.get("Disease_meta") or []}
        phenotype_steps: dict[str, dict] = {}
        disease_mims = []
        no_mim = False
        for slot, text in disease_picks.disease_cells(row):
            label = disease_picks.first_line(text)
            pm = disease_picks.phenotype_mim(label)
            if not pm:
                no_mim = True
                continue
            disease_mims.append(pm)
            lines = len([ln for ln in text.splitlines() if ln.strip()])
            source = (meta.get(slot) or {}).get("source", "")
            kinds = []
            if lines <= 1:
                kinds.append("missing")
            elif source != "capture" and _stale(curated_done, today):
                kinds.append("stale")
            if not kinds:
                continue
            step = phenotype_steps.setdefault(pm, {
                "type": "phenotype", "mim": pm, "url": f"https://omim.org/entry/{pm}",
                "labels": [], "kinds": [],
            })
            if label not in step["labels"]:
                step["labels"].append(label)
            for kind in kinds:
                if kind not in step["kinds"]:
                    step["kinds"].append(kind)
        gene_reviews = reviews.get(gene_mim, [])
        in_panel = gene in panel
        openable = _is_mim(gene_mim)
        is_synced = gene_mim in synced
        need_gene_step = openable and not is_synced and bool(
            in_panel or no_mim or gene_reviews or phenotype_steps
        )
        if not need_gene_step and not phenotype_steps and not gene_reviews:
            continue
        ph_steps = list(phenotype_steps.values())
        for step in ph_steps:
            step["minor"] = all(l[:1] in "{[?" for l in step["labels"])
            step["rank"] = max(KIND_RANK[k] for k in step["kinds"])
        ph_steps.sort(key=lambda s: (s["minor"], -s["rank"], s["mim"]))
        steps = ([{
            "type": "gene", "mim": gene_mim, "url": f"https://omim.org/entry/{gene_mim}",
            "labels": [f"{gene} 基因頁"], "kinds": ["gene_page"],
        }] if need_gene_step else []) + ph_steps
        if is_synced and (ph_steps or gene_reviews):
            status = "in_progress"
        elif need_gene_step:
            status = "gene_sync"
        else:
            status = "phenotype"
        genes.append({
            "gene": gene,
            "gene_mim": gene_mim,
            "url": f"https://omim.org/entry/{gene_mim}" if openable else "",
            "panel": in_panel,
            "synced_at": synced.get(gene_mim, ""),
            "done": done,
            "status": status,
            "steps": steps,
            "reviews": gene_reviews,
            "disease_mims": disease_mims,
            "minor": bool(ph_steps) and all(s["minor"] for s in ph_steps) and not need_gene_step,
        })
    with _lock:
        _cache.update(sig=sig, genes=genes)
    return genes


_STATUS_ORDER = {"in_progress": 0, "gene_sync": 1, "phenotype": 2}


def _decorate(genes: list[dict], scope: str) -> list[dict]:
    usage = reported_usage()
    out = []
    for base in genes:
        item = dict(base)
        item["reported"] = sum(usage.get(m, 0) for m in item["disease_mims"])
        if scope == "panel" and not item["panel"]:
            continue
        if scope == "reported" and not item["reported"]:
            continue
        out.append(item)

    def key(i):
        status = _STATUS_ORDER[i["status"]]
        # Panel genes waiting for their first gene-page sync come before
        # other genes that only have diseases to fill.
        group = status if not (status == 1 and not i["panel"]) else 3
        has_work = bool(i["reviews"]) + any(s["type"] == "phenotype" for s in i["steps"])
        return (group, -bool(i["reported"]), -i["panel"], -bool(i["reviews"]), -has_work, i["minor"],
                -i["reported"], i["gene"], i["gene_mim"])

    out.sort(key=key)
    return out


def queue(*, status: str = "", scope: str = "all", pending_mims: set[str] | None = None,
          recent_mims: set[str] | None = None, offset: int = 0, limit: int = 100) -> dict:
    pending_mims = pending_mims or set()
    recent_mims = recent_mims or set()
    genes = _decorate(build_genes(), scope)
    counts = {k: 0 for k in STATUS_LABELS}
    counts["panel_gene_sync"] = 0
    counts["phenotype_steps"] = 0
    counts["reviews"] = 0
    items = []
    for gene in genes:
        counts[gene["status"]] += 1
        if gene["status"] == "gene_sync" and gene["panel"]:
            counts["panel_gene_sync"] += 1
        counts["phenotype_steps"] += sum(1 for s in gene["steps"] if s["type"] == "phenotype")
        counts["reviews"] += len(gene["reviews"])
        if status == "review" and not gene["reviews"]:
            continue
        if status and status != "review" and gene["status"] != status:
            continue
        item = dict(gene)
        item["steps"] = [
            dict(s, captured=s["mim"] in pending_mims, recently_adopted=s["mim"] in recent_mims)
            for s in gene["steps"]
        ]
        items.append(item)
    return {
        "total": len(items),
        "counts": counts,
        "items": items[offset:offset + limit],
        "statuses": STATUS_LABELS,
    }


def next_step(*, current: str = "", gene: str = "", skip: set[str] | None = None,
              scope: str = "all", pending_mims: set[str] | None = None,
              recent_mims: set[str] | None = None) -> dict:
    """Next page to open, staying on the current gene until it is done."""
    skip = set(skip or ())
    pending_mims = pending_mims or set()
    recent_mims = recent_mims or set()
    genes = _decorate(build_genes(), scope)

    def available(g):
        if f"g{g['gene_mim']}" in skip:
            return []
        return [
            s for s in g["steps"]
            if s["mim"] != current and s["mim"] not in skip
            and s["mim"] not in pending_mims and s["mim"] not in recent_mims
        ]

    def answer(g, steps, continuing):
        total = len(g["steps"])
        return {
            "step": steps[0],
            "gene": {
                "gene": g["gene"], "gene_mim": g["gene_mim"], "url": g["url"], "panel": g["panel"],
                "status": g["status"], "reported": g["reported"], "reviews": len(g["reviews"]),
                "total_steps": total, "remaining_steps": len(steps),
            },
            "continuing": continuing,
            "remaining_genes": sum(1 for x in genes if available(x)),
        }

    # Stay on the gene being worked on.
    candidates = []
    if gene:
        candidates += [g for g in genes if g["gene_mim"] == gene]
    if current:
        candidates += [g for g in genes if g["gene_mim"] == current]
        candidates += [g for g in genes if any(s["mim"] == current for s in g["steps"])
                       or current in g["disease_mims"]]
    for g in candidates:
        steps = available(g)
        if steps:
            return answer(g, steps, True)
    for g in genes:
        steps = available(g)
        if steps:
            return answer(g, steps, False)
    return {"step": None, "gene": None, "continuing": False, "remaining_genes": 0}
