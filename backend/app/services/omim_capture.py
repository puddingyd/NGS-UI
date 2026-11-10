"""OMIM page captures sent by the browser extension (proof of concept).

The extension (``browser_extension/omim_capture``) runs only on
``omim.org/entry/*`` and only when the reviewer presses Capture in their own
browser session. It sends the rendered page text (what a manual Ctrl+A /
Ctrl+C would copy) plus a few DOM-extracted extras. This module

* parses that text into title / relationships / Description / Clinical
  Synopsis / OMIM edit dates (``parse_entry_text``);
* stores every capture append-only under ``OMIM_CAPTURE_DIR/{MIM}/`` as
  gzip JSON with an ``index.jsonl`` log (``store_capture``);
* compares the capture with the curator workbook currently loaded by
  ``omim_store`` so the extension can say whether the platform copy is
  missing, outdated or current (``compare_with_workbook``).

Nothing here writes to OMIM.xlsx yet; applying captures is a later phase.
"""
from __future__ import annotations

import gzip
import hashlib
import json
import re
import threading
from datetime import date, datetime, timezone
from pathlib import Path

from ..config import OMIM_CAPTURE_DIR
from . import disease_picks, omim_store

PARSER_VERSION = 1
MAX_TEXT_CHARS = 3_000_000

ENTRY_TYPES = {
    "*": "gene",
    "+": "gene_phenotype",
    "#": "phenotype",
    "%": "phenotype_locus",
    "^": "moved_removed",
    "": "other",
}

_HEADER_RE = re.compile(r"^\s*([#%*+^])?\s*(\d{6})\s*$")
_DATE_RE = re.compile(r"(\d{1,2})/(\d{1,2})/(\d{4})")
_MIM6_RE = re.compile(r"^\d{6}$")
_KEY_RE = re.compile(r"^[1-4]$")
_SECTION_RE = re.compile(r"^\s*[▼▲►]\s*(.*?)\s*:?\s*$")
_RELATIONSHIP_HEADINGS = {
    "phenotype-gene relationships",
    "gene-phenotype relationships",
}
_FOOTER_STARTS = (
    "Contributors:",
    "Creation Date:",
    "Edit History:",
    "NOTE: OMIM is intended",
)
_SYNOPSIS_TRAILERS = (
    "▲ Close",
    "▲Close",
)

_lock = threading.Lock()
_link_cache: dict = {"sig": None, "by_mim": {}, "by_gene_mim": {}}


# --------------------------------------------------------------------------
# Parsing
# --------------------------------------------------------------------------

def _lines(text: str) -> list[str]:
    text = str(text or "").replace("\r\n", "\n").replace("\r", "\n")
    return [line.replace(" ", " ").rstrip() for line in text.split("\n")]


def _strip_block(lines: list[str]) -> str:
    out = [line.rstrip() for line in lines]
    while out and not out[0].strip():
        out.pop(0)
    while out and not out[-1].strip():
        out.pop()
    # Collapse runs of blank lines to one so stored text stays compact.
    compact: list[str] = []
    for line in out:
        if not line.strip() and compact and not compact[-1].strip():
            continue
        compact.append(line.strip() if not line.strip() else line)
    return "\n".join(compact)


def _iso(month: str, day: str, year: str) -> str:
    try:
        return date(int(year), int(month), int(day)).isoformat()
    except ValueError:
        return ""


def _max_date(text: str) -> str:
    dates = [_iso(m, d, y) for m, d, y in _DATE_RE.findall(text or "")]
    dates = [d for d in dates if d]
    return max(dates) if dates else ""


def _find_header(lines: list[str], url_mim: str) -> tuple[int, str, str]:
    """Return (index, prefix, mim) of the entry header line.

    The table of contents repeats ``#232300`` right before "Table of
    Contents"; the real header (``# 232300``) is followed by the title.
    """
    candidates = []
    for idx, line in enumerate(lines):
        match = _HEADER_RE.match(line)
        if not match:
            continue
        if url_mim and match.group(2) != url_mim:
            continue
        nxt = next((ln.strip() for ln in lines[idx + 1:] if ln.strip()), "")
        if nxt.lower() == "table of contents":
            continue
        candidates.append((idx, match.group(1) or "", match.group(2)))
    if candidates:
        return candidates[0]
    return -1, "", url_mim


def _split_relationship_row(cells: list[str]) -> dict | None:
    cells = [c.strip() for c in cells]
    gene = gene_mim = ""
    rest = cells
    if (
        len(cells) >= 2
        and _MIM6_RE.match(cells[-1])
        and cells[-2]
        and not _MIM6_RE.match(cells[-2])
        and not _KEY_RE.match(cells[-2])
    ):
        gene, gene_mim = cells[-2], cells[-1]
        rest = cells[:-2]
    key_idx = None
    for idx in range(len(rest) - 1, -1, -1):
        if _KEY_RE.match(rest[idx]):
            key_idx = idx
            break
    if key_idx is None or key_idx < 2:
        return None
    phenotype = rest[key_idx - 3] if key_idx >= 3 else ""
    mim = rest[key_idx - 2]
    if mim and not _MIM6_RE.match(mim):
        # Shifted row (no inheritance cell): name sits right before the MIM.
        phenotype, mim, inheritance = rest[key_idx - 2], rest[key_idx - 1], ""
        if not _MIM6_RE.match(mim):
            return None
    else:
        inheritance = rest[key_idx - 1]
    if not phenotype and key_idx >= 3:
        return None
    location = rest[key_idx - 4] if key_idx >= 4 else ""
    if not phenotype:
        return None
    return {
        "location": location,
        "phenotype": phenotype,
        "phenotype_mim": mim if _MIM6_RE.match(mim or "") else "",
        "inheritance": inheritance,
        "mapping_key": rest[key_idx],
        "gene": gene,
        "gene_mim": gene_mim,
    }


def _parse_rows_from_cells(rows: list[list[str]]) -> list[dict]:
    out = []
    for cells in rows or []:
        if not isinstance(cells, list):
            continue
        row = _split_relationship_row([str(c or "") for c in cells])
        if row:
            out.append(row)
    return out


SYNOPSIS_CATEGORIES = {
    "INHERITANCE", "GROWTH", "HEAD & NECK", "CARDIOVASCULAR", "RESPIRATORY",
    "CHEST", "ABDOMEN", "GENITOURINARY", "SKELETAL", "SKIN, NAILS, & HAIR",
    "MUSCLE, SOFT TISSUES", "NEUROLOGIC", "VOICE", "METABOLIC FEATURES",
    "ENDOCRINE FEATURES", "HEMATOLOGY", "IMMUNOLOGY", "NEOPLASIA",
    "PRENATAL MANIFESTATIONS", "LABORATORY ABNORMALITIES", "MISCELLANEOUS",
    "MOLECULAR BASIS",
}


def _synopsis_from_page(text: str) -> str:
    """Clinical Synopsis block from a /clinicalSynopsis/{mim} page's text."""
    lines = _lines(text)
    start = next(
        (i for i, line in enumerate(lines) if line.strip() in SYNOPSIS_CATEGORIES),
        None,
    )
    if start is None:
        return ""
    block = []
    for line in lines[start:]:
        s = line.strip()
        if s in _SYNOPSIS_TRAILERS or s.startswith(_FOOTER_STARTS):
            break
        block.append(line)
    return _strip_block(block)


def parse_entry_text(text: str, *, url_mim: str = "", extracted: dict | None = None) -> dict:
    """Parse rendered OMIM entry text into structured fields."""
    extracted = extracted or {}
    lines = _lines(text)
    header_idx, prefix, mim = _find_header(lines, url_mim)
    warnings: list[str] = []
    if header_idx < 0:
        warnings.append("找不到 entry 標題列（# 123456）")
    if url_mim and mim and url_mim != mim:
        warnings.append(f"網址 MIM {url_mim} 與頁面 MIM {mim} 不一致")

    body = lines[header_idx + 1:] if header_idx >= 0 else lines
    stripped = [ln.strip() for ln in body]

    title = symbol = ""
    for line in stripped:
        if line:
            title, _, symbol = line.partition("; ")
            break

    toc: list[str] = []
    if header_idx >= 0:
        try:
            toc_start = [ln.strip() for ln in lines[:header_idx]].index("Table of Contents")
            for line in lines[toc_start + 1:header_idx]:
                s = line.strip()
                if s.startswith(("▼", "►", "▲")):
                    break
                if s:
                    toc.append(s)
        except ValueError:
            pass

    def index_of(pred, start=0):
        for idx in range(start, len(stripped)):
            if pred(stripped[idx]):
                return idx
        return -1

    first_section = index_of(lambda s: bool(_SECTION_RE.match(s)) and s.startswith("▼"))
    rel_idx = index_of(lambda s: s.lower() in _RELATIONSHIP_HEADINGS)
    alt_idx = index_of(lambda s: s.lower() == "alternative titles; symbols")
    incl_idx = index_of(lambda s: s.lower().startswith("other entities represented in this entry"))

    def block_until(start: int, stops: list[int]) -> list[str]:
        if start < 0:
            return []
        ends = [s for s in stops if s > start] or [len(stripped)]
        return [s for s in stripped[start + 1:min(ends)] if s]

    stops = [i for i in (rel_idx, incl_idx, first_section) if i >= 0]
    alt_titles = block_until(alt_idx, stops)
    included = block_until(incl_idx, [i for i in (rel_idx, first_section) if i >= 0])

    relationships: list[dict] = []
    synopsis = ""
    synopsis_source = ""
    if rel_idx >= 0:
        last_row = -1
        end = first_section if first_section > rel_idx else len(body)
        for idx in range(rel_idx + 1, end):
            raw = body[idx]
            if "\t" not in raw:
                if last_row >= 0:
                    break
                continue
            row = _split_relationship_row(raw.split("\t"))
            if row:
                relationships.append(row)
                last_row = idx
        if last_row >= 0:
            close_idx = index_of(lambda s: s in _SYNOPSIS_TRAILERS, last_row + 1)
            if 0 <= close_idx <= end:
                synopsis = _strip_block(body[last_row + 1:close_idx])
                synopsis_source = "entry_page"
    if not relationships and extracted.get("relationship_rows"):
        relationships = _parse_rows_from_cells(extracted.get("relationship_rows"))
        if relationships:
            warnings.append("關係表改用頁面表格結構解析")
    if not synopsis and extracted.get("synopsis_text"):
        synopsis = _synopsis_from_page(extracted["synopsis_text"])
        if synopsis:
            synopsis_source = "clinical_synopsis_page"

    entry_type = ENTRY_TYPES.get(prefix, "other")
    if entry_type in {"gene", "gene_phenotype"}:
        own_symbol = symbol.split(",")[0].strip()
        for row in relationships:
            if not row["gene"]:
                row["gene"], row["gene_mim"] = own_symbol, mim

    has_synopsis_in_toc = any(t.lower() == "clinical synopsis" for t in toc)
    if synopsis:
        synopsis_status = "present"
    elif entry_type in {"gene"}:
        synopsis_status = "not_applicable"
    elif has_synopsis_in_toc or not toc:
        synopsis_status = "not_loaded"
        warnings.append("Clinical Synopsis 尚未展開或沒有載入")
    else:
        synopsis_status = "absent"

    sections: dict[str, str] = {}
    order: list[str] = []
    current = None
    buf: list[str] = []
    footer_idx = len(body)
    for idx, raw in enumerate(body):
        s = raw.strip()
        if s.startswith(_FOOTER_STARTS):
            footer_idx = idx
            break
        match = _SECTION_RE.match(s) if s[:1] in "▼▲►" else None
        if match:
            if current and current not in sections:
                sections[current] = _strip_block(buf)
                order.append(current)
            name = match.group(1).strip().rstrip(":").strip()
            current = None if name.lower() == "close" else name
            buf = []
            continue
        if current:
            buf.append(raw)
    if current and current not in sections:
        sections[current] = _strip_block(buf)
        order.append(current)

    def section(name: str) -> str:
        for key in order:
            if key.lower() == name.lower():
                return sections[key]
        return ""

    description = section("Description")
    description_source = "Description" if description else ""
    if not description:
        description = section("Clinical Features")
        description_source = "Clinical Features" if description else ""
    if entry_type not in {"gene"} and not description:
        warnings.append("沒有 Description，也沒有 Clinical Features")
    elif description_source == "Clinical Features":
        warnings.append("沒有 Description，改用 Clinical Features")

    footer = body[footer_idx:]
    footer_text = {"contributors": "", "creation": "", "edit_history": ""}
    current_key = None
    for raw in footer:
        s = raw.strip()
        if s.startswith("Contributors:"):
            current_key = "contributors"
        elif s.startswith("Creation Date:"):
            current_key = "creation"
        elif s.startswith("Edit History:"):
            current_key = "edit_history"
        elif s.startswith("NOTE:") or not s:
            current_key = None if s.startswith("NOTE:") else current_key
            continue
        if current_key:
            footer_text[current_key] += s + "\n"

    own_rows = [r for r in relationships if r["phenotype_mim"] == mim]
    if entry_type in {"phenotype", "gene_phenotype", "phenotype_locus"} and relationships and not own_rows:
        warnings.append("關係表裡沒有這個 phenotype MIM 的列")
    genes = sorted({r["gene"] for r in own_rows if r["gene"]})
    if len(genes) > 1:
        warnings.append(f"此疾病連到 {len(genes)} 個基因：{', '.join(genes)}")

    parsed = {
        "parser_version": PARSER_VERSION,
        "mim": mim,
        "prefix": prefix,
        "entry_type": entry_type,
        "title": title,
        "symbol": symbol,
        "alternative_titles": alt_titles,
        "included_entities": included,
        "toc": toc,
        "relationships": relationships,
        "description": description,
        "description_source": description_source,
        "synopsis": synopsis,
        "synopsis_status": synopsis_status,
        "synopsis_source": synopsis_source,
        "sections": order,
        "omim_edit_date": _max_date(footer_text["edit_history"]),
        "omim_contributor_date": _max_date(footer_text["contributors"]),
        "omim_creation_date": _max_date(footer_text["creation"]),
        "warnings": warnings,
    }
    parsed["content_hash"] = content_hash(parsed)
    return parsed


def content_hash(parsed: dict) -> str:
    material = {
        k: parsed.get(k)
        for k in ("mim", "prefix", "title", "relationships", "description", "synopsis")
    }
    raw = json.dumps(material, ensure_ascii=False, sort_keys=True)
    return hashlib.sha256(raw.encode("utf-8")).hexdigest()


def disease_cell(parsed: dict, row: dict) -> str:
    """Compose the workbook Disease{N} cell exactly like process_genes.py."""
    parts = [f"{row['phenotype']} ({row['phenotype_mim']})({row['inheritance']})"]
    if parsed.get("description"):
        parts.append(parsed["description"])
    if parsed.get("synopsis"):
        parts.append(parsed["synopsis"])
    return "\n\n".join(parts)


# --------------------------------------------------------------------------
# Workbook comparison
# --------------------------------------------------------------------------

def _workbook_index() -> tuple[dict, dict]:
    omim_store.ensure_loaded()
    sig = omim_store.cache_signature()
    with _lock:
        if _link_cache["sig"] == sig:
            return _link_cache["by_mim"], _link_cache["by_gene_mim"]
        by_mim: dict[str, list[dict]] = {}
        by_gene_mim: dict[str, dict] = {}
        seen = set()
        rows = list(omim_store._state.get("by_gene", {}).values()) + list(
            omim_store._state.get("by_omim_id", {}).values()
        )
        for row in rows:
            if id(row) in seen:
                continue
            seen.add(id(row))
            gene_mim = str(row.get("OMIM_id") or "").split(".")[0]
            gene = row.get("gene_symbol") or ""
            if gene_mim:
                by_gene_mim.setdefault(gene_mim, row)
            lines = [ln.strip() for ln in str(row.get("OMIM_disease") or "").splitlines() if ln.strip()]
            cells = disease_picks.disease_cells(row)
            mims = {disease_picks.phenotype_mim(ln) for ln in lines}
            mims |= {disease_picks.phenotype_mim(t) for _s, t in cells}
            mims.discard("")
            for pm in mims:
                slot_hits = [
                    (slot, text) for slot, text in cells
                    if disease_picks.phenotype_mim(text) == pm
                ]
                rich = any(len([l for l in t.splitlines() if l.strip()]) > 1 for _s, t in slot_hits)
                by_mim.setdefault(pm, []).append({
                    "gene": gene,
                    "gene_mim": gene_mim,
                    "done": row.get("Done") or "",
                    "in_omim_disease": any(disease_picks.phenotype_mim(ln) == pm for ln in lines),
                    "slots": [s for s, _t in slot_hits],
                    "labels": [disease_picks.first_line(t) for _s, t in slot_hits]
                    or [ln for ln in lines if disease_picks.phenotype_mim(ln) == pm],
                    "has_detail": rich,
                })
        _link_cache.update(sig=sig, by_mim=by_mim, by_gene_mim=by_gene_mim)
        return by_mim, by_gene_mim


def _status_for_link(link: dict, edit_date: str) -> str:
    if not link["has_detail"]:
        return "missing_detail"
    if edit_date and link["done"] and edit_date > link["done"][:10]:
        return "outdated"
    if not link["done"]:
        return "no_done_date"
    return "up_to_date"


_STATUS_RANK = {
    "not_in_workbook": 5,
    "missing_detail": 4,
    "outdated": 3,
    "no_done_date": 2,
    "up_to_date": 1,
}


def compare_with_workbook(mim: str, *, entry_type: str = "", edit_date: str = "",
                          relationships: list[dict] | None = None) -> dict:
    """Describe how the curator workbook currently covers this MIM."""
    by_mim, by_gene_mim = _workbook_index()
    out: dict = {"mim": mim, "workbook_loaded": bool(by_gene_mim)}
    links = []
    for link in by_mim.get(mim, []):
        item = dict(link)
        item["status"] = _status_for_link(link, edit_date)
        links.append(item)
    out["phenotype_links"] = links

    if entry_type in {"gene", "gene_phenotype"} or mim in by_gene_mim:
        row = by_gene_mim.get(mim)
        gene_view = {"in_workbook": bool(row)}
        if row:
            lines = [ln.strip() for ln in str(row.get("OMIM_disease") or "").splitlines() if ln.strip()]
            wb = {}
            for ln in lines:
                pm = disease_picks.phenotype_mim(ln)
                wb.setdefault(pm or f"name:{disease_picks.name_key(ln)}", ln)
            omim_rows = [r for r in (relationships or []) if r.get("phenotype_mim")]
            omim_mims = {r["phenotype_mim"] for r in omim_rows}
            gene_view.update({
                "gene": row.get("gene_symbol") or "",
                "done": row.get("Done") or "",
                "workbook_lines": lines,
                "added_in_omim": [
                    r for r in omim_rows if r["phenotype_mim"] not in wb
                ],
                "missing_from_omim": [
                    ln for key, ln in wb.items()
                    if not key.startswith("name:") and key not in omim_mims
                ] if relationships else [],
                "workbook_lines_without_mim": [
                    ln for key, ln in wb.items() if key.startswith("name:")
                ],
                "renamed": [
                    {"phenotype_mim": r["phenotype_mim"], "omim": r["phenotype"],
                     "workbook": wb[r["phenotype_mim"]]}
                    for r in omim_rows
                    if r["phenotype_mim"] in wb
                    and disease_picks.name_key(r["phenotype"])
                    != disease_picks.name_key(wb[r["phenotype_mim"]])
                ],
            })
        out["gene"] = gene_view

    if links:
        out["status"] = max((l["status"] for l in links), key=lambda s: _STATUS_RANK[s])
    elif entry_type in {"gene", "gene_phenotype"} and out.get("gene", {}).get("in_workbook"):
        out["status"] = "gene_in_workbook"
    else:
        out["status"] = "not_in_workbook"
    return out


# --------------------------------------------------------------------------
# Storage
# --------------------------------------------------------------------------

def _capture_dir(mim: str) -> Path:
    if not _MIM6_RE.match(mim or ""):
        raise ValueError("invalid MIM")
    return OMIM_CAPTURE_DIR / mim


def store_capture(payload: dict, parsed: dict, *, user: dict) -> dict:
    """Write one capture; identical content to the latest one is not duplicated."""
    mim = parsed.get("mim") or ""
    directory = _capture_dir(mim)
    directory.mkdir(parents=True, exist_ok=True)
    digest = parsed["content_hash"]
    existing = sorted(directory.glob("*.json.gz"))
    duplicate = bool(existing) and existing[-1].name.endswith(f"_{digest[:12]}.json.gz")
    now = datetime.now(timezone.utc)
    capture_id = f"{mim}/{now.strftime('%Y%m%dT%H%M%S%fZ')}_{digest[:12]}"
    record = {
        "id": capture_id,
        "captured_at": now.isoformat(timespec="seconds"),
        "user": user.get("username") or "",
        "user_id": user.get("id"),
        "url": str(payload.get("url") or "")[:500],
        "extension_version": str(payload.get("extension_version") or "")[:40],
        "duplicate_of": existing[-1].name if duplicate else "",
    }
    if not duplicate:
        with gzip.open(directory / f"{capture_id.split('/')[1]}.json.gz", "wt", encoding="utf-8") as fh:
            json.dump({"meta": record, "payload": payload, "parsed": parsed}, fh, ensure_ascii=False)
    index_line = {
        **record,
        "mim": mim,
        "prefix": parsed.get("prefix"),
        "title": parsed.get("title"),
        "content_hash": digest,
        "omim_edit_date": parsed.get("omim_edit_date"),
        "synopsis_status": parsed.get("synopsis_status"),
        "description_source": parsed.get("description_source"),
        "relationship_count": len(parsed.get("relationships") or []),
        "warnings": parsed.get("warnings") or [],
    }
    OMIM_CAPTURE_DIR.mkdir(parents=True, exist_ok=True)
    with _lock, (OMIM_CAPTURE_DIR / "index.jsonl").open("a", encoding="utf-8") as fh:
        fh.write(json.dumps(index_line, ensure_ascii=False) + "\n")
    return {**record, "duplicate": duplicate}


def recent_captures(limit: int = 50) -> list[dict]:
    path = OMIM_CAPTURE_DIR / "index.jsonl"
    if not path.is_file():
        return []
    lines = path.read_text(encoding="utf-8").splitlines()
    out = []
    for line in reversed(lines):
        try:
            out.append(json.loads(line))
        except json.JSONDecodeError:
            continue
        if len(out) >= limit:
            break
    return out


def load_capture(capture_id: str) -> dict | None:
    mim, _, name = str(capture_id or "").partition("/")
    if not _MIM6_RE.match(mim) or not re.match(r"^[0-9TZ]+_[0-9a-f]{12}$", name):
        return None
    path = OMIM_CAPTURE_DIR / mim / f"{name}.json.gz"
    if not path.is_file():
        return None
    with gzip.open(path, "rt", encoding="utf-8") as fh:
        return json.load(fh)


def summarize(parsed: dict, comparison: dict) -> dict:
    """Compact response for the extension panel / capture page."""
    own = [r for r in parsed.get("relationships") or [] if r.get("phenotype_mim") == parsed.get("mim")]
    return {
        "mim": parsed.get("mim"),
        "prefix": parsed.get("prefix"),
        "entry_type": parsed.get("entry_type"),
        "title": parsed.get("title"),
        "symbol": parsed.get("symbol"),
        "relationships": parsed.get("relationships") or [],
        "description_chars": len(parsed.get("description") or ""),
        "description_source": parsed.get("description_source"),
        "description_preview": (parsed.get("description") or "")[:400],
        "synopsis_lines": len([l for l in (parsed.get("synopsis") or "").splitlines() if l.strip()]),
        "synopsis_status": parsed.get("synopsis_status"),
        "synopsis_source": parsed.get("synopsis_source"),
        "omim_edit_date": parsed.get("omim_edit_date"),
        "omim_contributor_date": parsed.get("omim_contributor_date"),
        "warnings": parsed.get("warnings") or [],
        "content_hash": parsed.get("content_hash"),
        "disease_cells": [
            {"gene": r.get("gene"), "cell": disease_cell(parsed, r)} for r in own
        ],
        "comparison": comparison,
    }
