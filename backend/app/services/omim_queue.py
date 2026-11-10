"""OMIM work queue: which OMIM pages to capture next.

Built from the OMIM rows ``omim_store`` currently serves (platform database
or OMIM.xlsx), one item per OMIM page to open:

* ``missing`` — the disease shows only its name on the card (the ``*``).
* ``no_mim``  — the disease line has no phenotype MIM; open the gene page.
* ``stale``   — curated text older than ``STALE_DAYS`` (gene Done date) and
                not refreshed from a capture since.

Diseases already ticked in any case report come first, then genes in the
disease-associated list (panel_loose), then the rest; within each group
{susceptibility} / [non-disease] / ?provisional entries go last. Items whose newest
capture is still waiting for adoption are marked ``captured`` and skipped
by ``next_item`` so the reviewer is not sent back to a page just captured;
pages adopted within the last 7 days are skipped the same way (e.g. a gene
page whose old MIM-less label still needs a review decision).
"""
from __future__ import annotations

import json
import threading
import time
from datetime import date, timedelta

from . import disease_picks, omim_store, panel_deadzone

STALE_DAYS = 365
KIND_RANK = {"missing": 3, "no_mim": 2, "stale": 1}
KIND_LABELS = {
    "missing": "缺說明（*）",
    "no_mim": "沒有 phenotype MIM",
    "stale": f"超過 {STALE_DAYS} 天未更新",
}
_REPORT_TTL = 300

_lock = threading.Lock()
_items_cache: dict = {"sig": None, "items": []}
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


def build_items() -> list[dict]:
    """All queue items (cached until the OMIM source changes)."""
    sig = omim_store.cache_signature()
    with _lock:
        if _items_cache["sig"] == sig:
            return _items_cache["items"]
    today = date.today()
    panel = panel_deadzone.disease_associated_genes()
    targets: dict[str, dict] = {}
    for row in _rows():
        gene = row.get("gene_symbol") or ""
        gene_mim = str(row.get("OMIM_id") or "").split(".")[0]
        done = row.get("Done") or ""
        meta = {m.get("slot"): m for m in row.get("Disease_meta") or []}
        for slot, text in disease_picks.disease_cells(row):
            label = disease_picks.first_line(text)
            pm = disease_picks.phenotype_mim(label)
            lines = len([ln for ln in text.splitlines() if ln.strip()])
            source = (meta.get(slot) or {}).get("source", "")
            kinds = []
            if lines <= 1:
                kinds.append("missing")
            if not pm:
                kinds.append("no_mim")
            if lines > 1 and source != "capture" and _stale(done, today):
                kinds.append("stale")
            if not kinds:
                continue
            target = pm or (gene_mim if gene_mim.isdigit() and len(gene_mim) == 6 else "")
            if not target:
                continue
            item = targets.get(target)
            if item is None:
                item = targets[target] = {
                    "target_mim": target,
                    "target_type": "phenotype" if pm else "gene",
                    "url": f"https://omim.org/entry/{target}",
                    "genes": [],
                    "labels": [],
                    "kinds": [],
                    "panel": False,
                    "done": done,
                }
            if gene and gene not in item["genes"]:
                item["genes"].append(gene)
            if label not in item["labels"]:
                item["labels"].append(label)
            for kind in kinds:
                if kind not in item["kinds"]:
                    item["kinds"].append(kind)
            item["panel"] = item["panel"] or gene in panel
            if done and (not item["done"] or done < item["done"]):
                item["done"] = done
    items = list(targets.values())
    for item in items:
        item["kinds"].sort(key=lambda k: -KIND_RANK[k])
        item["rank"] = max(KIND_RANK[k] for k in item["kinds"])
        # {susceptibility}, [non-disease] and ?provisional entries go last.
        item["minor"] = all(l[:1] in "{[?" for l in item["labels"])
    with _lock:
        _items_cache.update(sig=sig, items=items)
    return items


def queue(*, kind: str = "", scope: str = "all", pending_mims: set[str] | None = None,
          recent_mims: set[str] | None = None, offset: int = 0, limit: int = 200) -> dict:
    usage = reported_usage()
    pending_mims = pending_mims or set()
    recent_mims = recent_mims or set()
    items = []
    counts = {k: 0 for k in KIND_LABELS}
    for base in build_items():
        item = dict(base)
        item["reported"] = usage.get(item["target_mim"], 0)
        item["captured"] = item["target_mim"] in pending_mims
        item["recently_adopted"] = item["target_mim"] in recent_mims
        if scope == "panel" and not item["panel"]:
            continue
        if scope == "reported" and not item["reported"]:
            continue
        for k in item["kinds"]:
            counts[k] += 1
        if kind and kind not in item["kinds"]:
            continue
        items.append(item)
    items.sort(key=lambda i: (
        i["captured"], -bool(i["reported"]), -i["panel"], i["minor"],
        -(i["rank"] if not kind else KIND_RANK[kind]),
        -i["reported"], (i["genes"] or [""])[0], i["target_mim"],
    ))
    return {
        "total": len(items),
        "counts": counts,
        "items": items[offset:offset + limit],
        "kinds": KIND_LABELS,
    }


def next_item(*, current: str = "", skip: set[str] | None = None, kind: str = "",
              scope: str = "all", pending_mims: set[str] | None = None,
              recent_mims: set[str] | None = None) -> dict:
    """First item not just captured, recently adopted, current or skipped."""
    skip = set(skip or ())
    if current:
        skip.add(current)
    data = queue(kind=kind, scope=scope, pending_mims=pending_mims, recent_mims=recent_mims,
                 limit=100000)
    remaining = [
        i for i in data["items"]
        if not i["captured"] and not i["recently_adopted"] and i["target_mim"] not in skip
    ]
    return {"item": remaining[0] if remaining else None, "remaining": len(remaining)}
