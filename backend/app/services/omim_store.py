"""OMIM annotation lookup keyed by OMIM_id (primary) + gene_symbol.

Reads OMIM.xlsx (path from config.OMIM_XLSX) into two in-memory
dicts on first use, then serves O(1) lookups. The file mtime is
checked on every lookup so the operator can drop in a refreshed
xlsx without a server restart — a stale cache would silently serve
old disease text otherwise.

Workbook layout (one row per OMIM phenotype-bearing gene):
    OMIM_id | gene_symbol | OMIM_disease | Inheritance |
    Disease1 .. Disease5 | Done

`Disease1..5` carry the curator-written rich text (~9% of rows have
Disease1; rarer for the rest). For the ~28% of genes that only have
`OMIM_disease` (multi-line, each line `<phenotype name> (<inh>)`),
we synthesise Disease1..N from those lines so the variant card has
something to show even when the curator hasn't filled detail yet.

A single gene_symbol can map to multiple OMIM_id rows (~248 cases);
we pick by OMIM_id first when the variant TSV's OMIM_LINK supplies
one, falling back to the first gene_symbol hit.
"""
from __future__ import annotations

import re
import threading
from pathlib import Path
from typing import Optional

from ..config import OMIM_XLSX

_DISEASE_FIELDS = ("Disease1", "Disease2", "Disease3", "Disease4", "Disease5")
_OMIM_URL_RE = re.compile(r"/entry/(\d+)\b")

_lock = threading.Lock()
_state: dict = {
    "mtime": None,
    "by_omim_id": {},   # int → row dict
    "by_gene":    {},   # str → row dict (first hit wins; see _ingest_row)
}


def _empty_row() -> dict:
    return {
        "OMIM_id": "",
        "OMIM_disease": "",
        "Inheritance": "",
        "Disease1": "", "Disease2": "", "Disease3": "",
        "Disease4": "", "Disease5": "",
    }


def _row_to_dict(headers: tuple, row: tuple) -> dict:
    """Project an xlsx row to the fields we expose. Missing columns
    just stay empty so future column re-orderings don't crash."""
    out = _empty_row()
    h2i = {h: i for i, h in enumerate(headers) if h}
    def get(name):
        i = h2i.get(name)
        if i is None or i >= len(row):
            return ""
        v = row[i]
        return "" if v is None else str(v).strip()
    out["OMIM_id"]      = get("OMIM_id")
    out["OMIM_disease"] = get("OMIM_disease")
    out["Inheritance"]  = get("Inheritance")
    for f in _DISEASE_FIELDS:
        out[f] = get(f)
    return out


def _synthesize_diseases(row: dict) -> None:
    """When Disease1..5 are empty but OMIM_disease has lines, fill
    Disease1..N with those lines so the UI renders summaries."""
    if any(row[f] for f in _DISEASE_FIELDS):
        return
    od = row.get("OMIM_disease") or ""
    if not od.strip():
        return
    lines = [ln.strip() for ln in od.splitlines() if ln.strip()]
    for i, ln in enumerate(lines[:5]):
        row[_DISEASE_FIELDS[i]] = ln


def _load(path: Path) -> tuple[dict, dict]:
    """Read the xlsx into (by_omim_id, by_gene)."""
    import openpyxl  # imported lazily so a missing dep doesn't break boot
    wb = openpyxl.load_workbook(str(path), read_only=True, data_only=True)
    ws = wb.worksheets[0]
    rows = ws.iter_rows(values_only=True)
    headers = next(rows, None)
    if not headers:
        wb.close()
        return {}, {}
    by_omim_id: dict[int, dict] = {}
    by_gene: dict[str, dict] = {}
    for raw in rows:
        if not raw:
            continue
        rec = _row_to_dict(headers, raw)
        gene = ""
        # gene_symbol isn't a key we expose downstream but we need it
        # to fan out the by_gene index. Re-read it directly here.
        h2i = {h: i for i, h in enumerate(headers) if h}
        gi = h2i.get("gene_symbol")
        if gi is not None and gi < len(raw) and raw[gi] is not None:
            gene = str(raw[gi]).strip()
        if not rec["OMIM_id"] and not gene:
            continue
        _synthesize_diseases(rec)
        try:
            oid = int(rec["OMIM_id"]) if rec["OMIM_id"] else None
        except ValueError:
            oid = None
        if oid is not None:
            by_omim_id.setdefault(oid, rec)
        if gene:
            by_gene.setdefault(gene, rec)
    wb.close()
    return by_omim_id, by_gene


def _ensure_loaded() -> None:
    """Load (or reload) the xlsx if it changed on disk. Silent on
    errors — variant payloads just lose their OMIM annotation."""
    path = OMIM_XLSX
    if not path or not path.is_file():
        # Wipe the cache so a previously-loaded xlsx isn't served
        # after the operator removed it.
        if _state["mtime"] is not None:
            with _lock:
                _state["mtime"] = None
                _state["by_omim_id"] = {}
                _state["by_gene"] = {}
        return
    try:
        mtime = path.stat().st_mtime
    except OSError:
        return
    if _state["mtime"] == mtime:
        return
    with _lock:
        if _state["mtime"] == mtime:
            return
        try:
            by_id, by_gene = _load(path)
        except Exception:
            # Don't poison the cache on parse error — keep whatever
            # we had so the UI degrades gracefully.
            return
        _state["mtime"] = mtime
        _state["by_omim_id"] = by_id
        _state["by_gene"] = by_gene


def parse_omim_id_from_link(link: str) -> Optional[int]:
    """Extract the trailing integer from an omim.org/entry/<id> URL."""
    if not link:
        return None
    m = _OMIM_URL_RE.search(link)
    if not m:
        return None
    try:
        return int(m.group(1))
    except ValueError:
        return None


def lookup(*, omim_id: Optional[int] = None, gene: str = "") -> Optional[dict]:
    """Return the joined OMIM record (OMIM_id + Disease1..5 + ...).

    Prefers OMIM_id (precise; the xlsx is keyed there). Falls back to
    gene_symbol (first hit) when OMIM_id was missing or unmatched.
    Returns None when neither key resolves.
    """
    _ensure_loaded()
    by_id = _state["by_omim_id"]
    by_gene = _state["by_gene"]
    if omim_id is not None:
        rec = by_id.get(omim_id)
        if rec:
            return rec
    if gene:
        return by_gene.get(gene)
    return None


def is_loaded() -> bool:
    """True when the xlsx was successfully ingested at least once."""
    _ensure_loaded()
    return bool(_state["by_omim_id"] or _state["by_gene"])
