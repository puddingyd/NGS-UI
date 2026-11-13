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

With the OMIM Update List stored (``omim_updates``), every OMIM update to
content the platform holds is queued when it is newer than the platform's
version: ``updated`` (phenotype page) and ``gene_updated`` (gene page —
the adopted gene page, else the workbook Done date; a gene never synced
and without a Done date counts as older than any update). Update List
entries the platform does not carry at all become their own items:
``new_disease`` (``#``/``%`` no platform gene has) and ``new_gene``. Workbook text whose Done date is covered by
the stored history with no update since is ``confirmed_current``; with the
``stale_mode`` setting ``updates_only`` such text is no longer queued as
"older than a year".

Gene order: genes in progress (gene page adopted, or some disease already
updated from a capture — possibly via another gene sharing the phenotype —
but first-pass steps or items to confirm remain) → genes with OMIM updates
to content already done → disease-associated (panel_loose) genes whose
gene page has not been synced yet → everything else; inside each group,
genes whose diseases were ticked in case reports come first, then genes
with review items, then genes with diseases to fill. Within a gene,
the gene page comes first, then diseases (missing before stale; entries in
{}, [] or with "?" last).

``next_step`` keeps the reviewer on the same gene until its steps are done,
then moves to in-progress genes sharing a disease with it,
skipping pages captured but not adopted yet, pages adopted in the last 7
days, and anything the extension asked to skip (``MIM`` or ``g<geneMIM>``).
"""
from __future__ import annotations

import json
import threading
import time
from datetime import date, timedelta

from . import disease_picks, omim_capture, omim_db, omim_store, omim_updates, panel_deadzone

STALE_DAYS = 365
KIND_RANK = {"missing": 3, "updated": 2, "new_disease": 2, "stale": 1}
UPDATE_KINDS = {"gene_updated", "updated", "new_disease", "new_gene"}
FIRST_PASS_KINDS = {"gene_page", "gene_moved", "missing", "stale"}
STATUS_LABELS = {
    "in_progress": "進行中（已同步基因頁或已有疾病被更新，但尚未完成）",
    "updated": "OMIM 有更新（做過的基因或疾病，OMIM 之後又改了）",
    "gene_sync": "基因頁待同步",
    "phenotype": "只有疾病待補",
}
_REPORT_TTL = 300

_lock = threading.Lock()
_cache: dict = {"sig": None, "genes": [], "stats": {}}
_report_cache: dict = {"at": 0.0, "usage": {}, "samples": {}}


def reported_usage(ttl: int = _REPORT_TTL) -> dict[str, int]:
    """phenotype MIM -> number of variant reports that ticked it."""
    _scan_reports(ttl)
    return _report_cache["usage"]


def reported_samples(ttl: int = _REPORT_TTL) -> dict[str, list[str]]:
    """phenotype MIM -> sample IDs (檢體編號) whose report ticked it."""
    _scan_reports(ttl)
    return _report_cache["samples"]


def _scan_reports(ttl: int) -> None:
    now = time.monotonic()
    if _report_cache["at"] and now - _report_cache["at"] < ttl:
        return
    from . import sample_layout

    usage: dict[str, int] = {}
    samples: dict[str, set[str]] = {}
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
                    samples.setdefault(mim, set()).add(str(sample_id))
    _report_cache.update(at=now, usage=usage, samples={m: sorted(s) for m, s in samples.items()})


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
    return (omim_store.cache_signature(), omim_db.signature() if omim_db.exists() else None,
            omim_capture.index_signature(), omim_updates.us_today())


def build_genes() -> list[dict]:
    """Every gene with work to do (cached until OMIM data or the DB changes)."""
    sig = _signature()
    with _lock:
        if _cache["sig"] == sig:
            return _cache["genes"]
    today = date.today()
    panel = panel_deadzone.disease_associated_genes()
    db = omim_db.exists()
    synced = omim_db.gene_sync_map() if db else {}
    reviews = omim_db.open_reviews_by_gene() if db else {}
    adopted = omim_db.latest_adoptions() if db else {}
    redirects = omim_db.mim_redirects() if db else {}
    texts = omim_db.phenotype_text_sources() if db else {}
    entry_latest, any_latest = omim_updates.latest_dates() if db else ({}, {})
    coverage_start = omim_updates.coverage().get("coverage_start", "") if db else ""
    mode = omim_updates.stale_mode() if db else "both"
    times = omim_capture.capture_times() if (adopted or texts) else {}
    stats = {"stale_hidden": 0}

    def capture_day(capture_id: str, adopted_at: str) -> str:
        return omim_updates.us_date(times.get(capture_id) or adopted_at)

    genes = []
    for row in _rows():
        gene = row.get("gene_symbol") or ""
        gene_mim = str(row.get("OMIM_id") or "").split(".")[0]
        done = row.get("Done") or ""
        # Curated workbook text is judged by the workbook's own Done date.
        curated_done = (row.get("Legacy_done", done) or "")[:10]
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
            listed = any_latest.get(pm, "")
            kinds = []
            confirmed = False
            if pm == gene_mim and source != "gene_page" and lines <= 1:
                # The phenotype is listed under the gene's own MIM: capturing
                # the gene page is the step (no separate disease page).
                no_mim = True
                continue
            if source == "gene_page":
                continue  # follows the gene page (gene step covers updates)
            if lines <= 1:
                kinds.append("missing")
            elif source == "capture":
                src = texts.get(pm)
                if src and listed and listed > max(capture_day(src[0], src[1]), src[2] or ""):
                    kinds.append("updated")
            elif listed and curated_done and listed > curated_done:
                kinds.append("updated")
            elif _stale(curated_done, today):
                # With the Update List read back past the Done date and no
                # update since, the curated text is known to be current.
                confirmed = bool(coverage_start and curated_done and curated_done[:7] >= coverage_start)
                if confirmed and mode == "updates_only":
                    stats["stale_hidden"] += 1
                else:
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
            if "updated" in kinds:
                step["update_date"] = listed
            if confirmed:
                step["confirmed_current"] = True
        gene_reviews = reviews.get(gene_mim, [])
        in_panel = gene in panel
        # OMIM moved the gene's entry (e.g. 610815 SLC25A32 -> 138480): the
        # gene page to capture is the new MIM; a removed entry has no page.
        moved = redirects.get(gene_mim)
        page_mim = moved[0] if moved else gene_mim
        openable = _is_mim(page_mim) and moved != []
        # Any adopted capture of the gene's page MIM counts as its page synced
        # — also when OMIM lists that MIM as a locus (%) or phenotype entry
        # (e.g. HYT1, SCZD1, BMND7 rows in the workbook).
        is_synced = page_mim in synced or page_mim in adopted if moved else (
            gene_mim in synced or gene_mim in adopted)
        # Gene page OMIM updated after the version the platform holds: the
        # latest adopted capture of that MIM (whatever its entry type — old
        # or pseudogene entries have no "*" prefix), else the workbook Done
        # date (none = unknown, so any listed update counts). Re-capturing
        # unchanged content refreshes the capture time.
        listed = entry_latest.get(page_mim, "")
        content = capture_day(*adopted[page_mim]) if page_mim in adopted else (
            "" if moved else curated_done)
        gene_update = listed if listed and listed > content else ""
        need_gene_step = openable and (bool(gene_update) or (not is_synced and bool(
            in_panel or no_mim or gene_reviews or phenotype_steps or moved
        )))
        if not need_gene_step and not phenotype_steps and not gene_reviews:
            continue
        ph_steps = list(phenotype_steps.values())
        for step in ph_steps:
            step["minor"] = all(l[:1] in "{[?" for l in step["labels"])
            step["rank"] = max(KIND_RANK[k] for k in step["kinds"])
        ph_steps.sort(key=lambda s: (s["minor"], -s["rank"], s["mim"]))
        gene_step = []
        if need_gene_step:
            kind = "gene_updated" if gene_update and is_synced else ("gene_moved" if moved else (
                "gene_updated" if gene_update else "gene_page"))
            gene_step = [{
                "type": "gene", "mim": page_mim, "url": f"https://omim.org/entry/{page_mim}",
                "labels": [f"{gene} 基因頁" + (f"（OMIM 已從 {gene_mim} 搬到 {page_mim}）" if moved else "")],
                "kinds": [kind],
            }]
            if gene_update:
                gene_step[0]["update_date"] = gene_update
        steps = gene_step + ph_steps
        # Some disease of this gene already carries adopted capture text —
        # possibly via another gene sharing the phenotype (e.g. 167800 for
        # CFTR/CTRC/SPINK1/PRSS1/PRSS2) — so finish it soon.
        touched = any(m.get("source") == "capture" for m in meta.values())
        kinds_all = {k for st in steps for k in st["kinds"]}
        first_pass = bool(kinds_all & FIRST_PASS_KINDS) or bool(gene_reviews)
        has_update = bool(kinds_all & UPDATE_KINDS)
        if (is_synced or touched) and first_pass:
            status = "in_progress"
        elif has_update and (is_synced or touched or not need_gene_step):
            status = "updated"
        elif need_gene_step:
            status = "gene_sync"
        else:
            status = "phenotype"
        genes.append({
            "gene": gene,
            "gene_mim": gene_mim,
            "url": f"https://omim.org/entry/{page_mim}" if openable else "",
            "moved_to": moved,
            "panel": in_panel,
            "synced_at": synced.get(gene_mim, "") or (adopted[gene_mim][1] if gene_mim in adopted else ""),
            "done": done,
            "status": status,
            "has_update": has_update,
            "partially_updated": touched and not is_synced,
            "steps": steps,
            "reviews": gene_reviews,
            "disease_mims": disease_mims,
            "minor": bool(ph_steps) and all(s["minor"] for s in ph_steps) and not need_gene_step,
        })
    # Update List entries the platform does not carry: genes missing from
    # the platform, and diseases no platform gene has. Done once a capture
    # of that page newer than OMIM's update was adopted.
    for item in omim_updates.missing_entries() if db else []:
        mim = item["mim"]
        if mim in adopted and capture_day(*adopted[mim]) >= item["listed_date"]:
            continue
        is_gene = item["cls"] == "gene"
        title = item["title"] or ""
        symbol = title.rsplit("; ", 1)[1].strip() if is_gene and "; " in title else ""
        genes.append({
            "gene": symbol or ("平台沒有的基因" if is_gene else "平台沒有的疾病"),
            "gene_mim": f"new{mim}", "url": "", "panel": False,
            "synced_at": "", "done": "", "status": "updated", "has_update": True,
            "partially_updated": False, "new_entry": True, "prefix": item["prefix"],
            "entry": {k: item[k] for k in ("mim", "prefix", "cls", "listed_date", "first_date", "title", "is_new")},
            # Diseases with a known gene come with the other OMIM updates;
            # missing genes and %-entries (no known gene) after everything.
            "late": is_gene or item["prefix"] != "#",
            "steps": [{
                "type": "gene" if is_gene else "phenotype", "mim": mim, "url": f"https://omim.org/entry/{mim}",
                "labels": [f"{item['prefix']}{mim} {title}".strip()],
                "kinds": ["new_gene" if is_gene else "new_disease"],
                "update_date": item["listed_date"], "minor": False, "rank": KIND_RANK["new_disease"],
            }],
            "reviews": [], "disease_mims": [] if is_gene else [mim], "minor": False,
        })
    with _lock:
        _cache.update(sig=sig, genes=genes, stats=stats)
    return genes


_STATUS_ORDER = {"in_progress": 0, "updated": 1, "gene_sync": 2, "phenotype": 3}


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
        group = status if not (i["status"] == "gene_sync" and not i["panel"]) else 4
        if i.get("late"):
            group = 5
        has_work = bool(i["reviews"]) + any(s["type"] == "phenotype" for s in i["steps"])
        return (group, -bool(i["reported"]), -i["panel"], -i.get("has_update", False), -bool(i["reviews"]),
                -has_work, i["minor"], -i["reported"], i["gene"], i["gene_mim"])

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
    counts["updated_steps"] = 0
    counts["new_disease_steps"] = 0
    counts["new_gene_steps"] = 0
    counts["gene_update_steps"] = 0
    counts["stale_steps"] = 0
    counts["stale_confirmed"] = 0
    counts["stale_hidden"] = (_cache.get("stats") or {}).get("stale_hidden", 0)
    items = []
    for gene in genes:
        counts[gene["status"]] += 1
        if gene["status"] == "gene_sync" and gene["panel"]:
            counts["panel_gene_sync"] += 1
        counts["phenotype_steps"] += sum(1 for s in gene["steps"] if s["type"] == "phenotype")
        counts["reviews"] += len(gene["reviews"])
        for st in gene["steps"]:
            kinds = set(st["kinds"])
            counts["updated_steps"] += bool(kinds & {"updated", "gene_updated"})
            counts["new_disease_steps"] += "new_disease" in kinds
            counts["new_gene_steps"] += "new_gene" in kinds
            counts["gene_update_steps"] += "gene_updated" in kinds
            counts["stale_steps"] += kinds == {"stale"}
            counts["stale_confirmed"] += kinds == {"stale"} and bool(st.get("confirmed_current"))
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
              recent_mims: set[str] | None = None, captured_days: dict[str, str] | None = None) -> dict:
    """Next page to open, staying on the current gene until it is done."""
    skip = set(skip or ())
    pending_mims = pending_mims or set()
    recent_mims = recent_mims or set()
    captured_days = captured_days or {}
    genes = _decorate(build_genes(), scope)

    def available(g):
        if f"g{g['gene_mim']}" in skip:
            return []
        # A page adopted recently is not suggested again, unless OMIM
        # updated it after the adopted version — but a page already captured
        # on or after that update never comes back (no capture loops).
        return [
            s for s in g["steps"]
            if s["mim"] != current and s["mim"] not in skip and s["mim"] not in pending_mims
            and (s["mim"] not in recent_mims or set(s["kinds"]) & UPDATE_KINDS)
            and not (s.get("update_date") and captured_days.get(s["mim"], "") >= s["update_date"])
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
    # Then genes sharing a disease with the one just finished (they were
    # partly updated by the same adoption).
    shared = {m for g in candidates for m in g["disease_mims"]}
    if shared:
        for g in genes:
            if g["status"] == "in_progress" and shared.intersection(g["disease_mims"]):
                steps = available(g)
                if steps:
                    result = answer(g, steps, False)
                    result["related"] = True
                    return result
    for g in genes:
        steps = available(g)
        if steps:
            return answer(g, steps, False)
    return {"step": None, "gene": None, "continuing": False, "remaining_genes": 0}
