"""Supplemental gene-disease associations for variant-card display.

OMIM.xlsx remains the source of curator-maintained rich text. This module
loads an optional local TSV of supplemental GenCC / ClinGen / MONDO-aligned
associations and merges them with the current OMIM row at runtime, so edits
to OMIM Disease1..5 keep showing up without rebuilding this file.
"""
from __future__ import annotations

import re
import sqlite3
import threading
from pathlib import Path

from ..config import GENE_DISEASE_DB, GENE_DISEASE_TSV
from . import panel_deadzone

_DISEASE_FIELDS = ("Disease1", "Disease2", "Disease3", "Disease4", "Disease5")
_INHERITANCE_RE = re.compile(
    r"\((AD|AR|XLD|XLR|XL|YL|MT|Mi|DR|DD|Smu|Mu|Isol)"
    r"(?:\s*[/,;]\s*(?:AD|AR|XLD|XLR|XL|YL|MT|Mi|DR|DD|Smu|Mu|Isol))*\)",
    re.IGNORECASE,
)
_MIM_RE = re.compile(r"\((\d{6})\)")
_WORD_RE = re.compile(r"[^a-z0-9]+")

_EVIDENCE_RANK = {
    "definitive": 50,
    "strong": 40,
    "moderate": 30,
    "limited": 20,
    "supportive": 10,
    "animal model only": 5,
}
_EXCLUDED_EVIDENCE = {
    "disputed",
    "refuted",
    "no known disease relationship",
    "no reported evidence",
}

_lock = threading.Lock()
_state = {
    "path": "",
    "mtime": None,
    "size": 0,
    "count": 0,
    "by_gene": {},
}


def _norm_header(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", str(value or "").strip().lower()).strip("_")


def _first(row: dict, *names: str) -> str:
    for name in names:
        value = row.get(name)
        if value not in (None, ""):
            return str(value).strip()
    return ""


def _evidence_rank(value: str) -> int:
    text = str(value or "").strip().lower()
    return _EVIDENCE_RANK.get(text, 0)


def _canonical_gene(value: str) -> str:
    gene = panel_deadzone.canonical_panel_gene_symbol(value or "")
    return gene or str(value or "").strip()


def _name_key(value: str) -> str:
    text = str(value or "").lower()
    text = re.sub(r"\([^)]*\)", " ", text)
    text = text.replace("&", " and ")
    return _WORD_RE.sub(" ", text).strip()


def _parse_omim_disease(text: str) -> tuple[str, str, str]:
    first_line = str(text or "").splitlines()[0].strip()
    mim_match = _MIM_RE.search(first_line)
    inh_match = _INHERITANCE_RE.search(first_line)
    metadata_starts = [m.start() for m in (mim_match, inh_match) if m]
    name = first_line[:min(metadata_starts)] if metadata_starts else first_line
    return (
        name.rstrip(" :,;").strip(),
        mim_match.group(1) if mim_match else "",
        inh_match.group(0)[1:-1] if inh_match else "",
    )


def _association_key(item: dict) -> str:
    mim = str(item.get("phenotype_mim") or "").strip()
    if mim:
        return f"mim:{mim}"
    mondo = str(item.get("mondo_id") or "").strip()
    if mondo:
        return f"mondo:{mondo.lower()}"
    source_id = str(item.get("source_disease_id") or "").strip()
    if source_id:
        return f"source:{source_id.lower()}"
    return f"name:{_name_key(item.get('display_name') or '')}"


def _source_label(source: str, evidence: str = "") -> str:
    source = str(source or "").strip()
    evidence = str(evidence or "").strip()
    if source and evidence:
        return f"{source} {evidence}"
    return source or evidence


def _clean_sources(values: list[str]) -> list[str]:
    out = []
    seen = set()
    for value in values:
        text = str(value or "").strip()
        if not text:
            continue
        key = text.lower()
        if key in seen:
            continue
        seen.add(key)
        out.append(text)
    return out


def _empty_indexes() -> dict:
    return {}


def _row_to_item(row: dict) -> dict:
    source = _first(row, "source", "database", "submitter")
    evidence = _first(
        row, "classification", "evidence_level", "validity",
        "gene_disease_validity", "assertion",
    )
    rank = _evidence_rank(evidence)
    item = {
        "id": "",
        "source_kind": "supplemental",
        "display_name": _first(
            row, "disease_name", "disease", "disease_label",
            "condition", "condition_name", "mondo_label",
        ),
        "phenotype_mim": _first(
            row, "phenotype_mim", "phenotype_mim_id", "omim_id",
            "omim_phenotype_id", "mim", "mim_number",
        ),
        "inheritance": _first(row, "inheritance", "mode_of_inheritance", "moi"),
        "detail": "",
        "mondo_id": _first(row, "mondo_id", "mondo", "disease_mondo"),
        "source_disease_id": _first(
            row, "source_disease_id", "disease_id", "condition_id",
            "curie", "iri", "online_report",
        ),
        "sources": _clean_sources([source]),
        "evidence": _clean_sources([_source_label(source, evidence)]),
        "evidence_rank": rank,
        "is_omim_rich": False,
        "needs_description": False,
    }
    item["id"] = _association_key(item)
    return item


def _load_tsv(path: Path) -> dict[str, list[dict]]:
    import csv

    by_gene: dict[str, list[dict]] = {}
    with path.open("r", encoding="utf-8-sig", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not reader.fieldnames:
            return by_gene
        for raw in reader:
            row = {_norm_header(k): (v or "").strip() for k, v in raw.items()}
            gene = _canonical_gene(_first(
                row, "gene_symbol", "gene", "hgnc_symbol", "symbol"
            ))
            if not gene:
                continue
            evidence = _first(
                row, "classification", "evidence_level", "validity",
                "gene_disease_validity", "assertion",
            )
            evidence_key = evidence.lower()
            if evidence_key in _EXCLUDED_EVIDENCE:
                continue
            rank = _evidence_rank(evidence)
            if evidence_key in _EVIDENCE_RANK and rank < _EVIDENCE_RANK["limited"]:
                continue
            item = _row_to_item(row)
            if not item["display_name"]:
                continue
            by_gene.setdefault(gene, []).append(item)
    return by_gene


def _load_sqlite(path: Path) -> tuple[dict[str, list[dict]], int]:
    by_gene: dict[str, list[dict]] = {}
    count = 0
    con = sqlite3.connect(str(path))
    con.row_factory = sqlite3.Row
    try:
        rows = con.execute(
            """
            SELECT gene, disease_name, source, classification, mondo_id,
                   phenotype_mim, source_disease_id, inheritance, evidence_rank
            FROM gene_disease_associations
            WHERE include_in_ui = 1
            ORDER BY gene, evidence_rank DESC, source, disease_name
            """
        )
        for row in rows:
            gene = row["gene"]
            item = _row_to_item(dict(row))
            if not gene or not item["display_name"]:
                continue
            by_gene.setdefault(gene, []).append(item)
            count += 1
    finally:
        con.close()
    return by_gene, count


def _active_path() -> Path:
    return GENE_DISEASE_DB if GENE_DISEASE_DB.is_file() else GENE_DISEASE_TSV


def _ensure_loaded() -> None:
    path = _active_path()
    if not path or not path.is_file():
        if _state["mtime"] is not None:
            with _lock:
                _state["path"] = ""
                _state["mtime"] = None
                _state["size"] = 0
                _state["count"] = 0
                _state["by_gene"] = _empty_indexes()
        return
    try:
        st = path.stat()
    except OSError:
        return
    path_key = str(path)
    if (
        _state["path"] == path_key
        and _state["mtime"] == st.st_mtime_ns
        and _state["size"] == st.st_size
    ):
        return
    with _lock:
        if (
            _state["path"] == path_key
            and _state["mtime"] == st.st_mtime_ns
            and _state["size"] == st.st_size
        ):
            return
        try:
            if path.suffix.lower() in {".sqlite", ".db", ".sqlite3"}:
                by_gene, count = _load_sqlite(path)
            else:
                by_gene = _load_tsv(path)
                count = sum(len(v) for v in by_gene.values())
        except Exception:
            return
        _state["path"] = path_key
        _state["mtime"] = st.st_mtime_ns
        _state["size"] = st.st_size
        _state["count"] = count
        _state["by_gene"] = by_gene


def ensure_loaded() -> None:
    _ensure_loaded()


def cache_signature() -> tuple:
    _ensure_loaded()
    return (
        _state["path"] or str(_active_path()),
        _state["mtime"] or 0,
        _state["size"] or 0,
        _state["count"] or 0,
        len(_state["by_gene"]),
    )


def lookup_cached(gene: str) -> list[dict]:
    if not gene:
        return []
    canonical = _canonical_gene(gene)
    return list((_state["by_gene"] or {}).get(canonical, []))


def _omim_associations(omim_row: dict | None) -> list[dict]:
    if not omim_row:
        return []
    out = []
    for idx, field in enumerate(_DISEASE_FIELDS, start=1):
        detail = str((omim_row or {}).get(field) or "").strip()
        if not detail or detail == "NA":
            continue
        name, mim, inheritance = _parse_omim_disease(detail)
        if not name and not mim:
            continue
        line_count = len([ln for ln in detail.splitlines() if ln.strip()])
        item = {
            "id": f"omim:{mim}" if mim else f"omim-slot:{idx}",
            "source_kind": "omim",
            "display_name": name or detail.splitlines()[0].strip(),
            "phenotype_mim": mim,
            "inheritance": inheritance or str(omim_row.get("Inheritance") or "").strip(),
            "detail": detail,
            "mondo_id": "",
            "source_disease_id": f"OMIM:{mim}" if mim else "",
            "sources": ["OMIM"],
            "evidence": ["OMIM"],
            "evidence_rank": 60,
            "is_omim_rich": True,
            "needs_description": line_count <= 1,
            "omim_slot": idx,
        }
        out.append(item)
    return out


def merged_associations(
    gene: str,
    omim_row: dict | None = None,
    *,
    refresh: bool = True,
) -> list[dict]:
    """Return OMIM rich rows plus non-duplicate supplemental associations."""
    if refresh:
        _ensure_loaded()
    clusters: dict[str, dict] = {}
    order: list[str] = []

    for item in _omim_associations(omim_row):
        key = _association_key(item)
        clusters[key] = item
        order.append(key)

    for item in lookup_cached(gene):
        key = _association_key(item)
        existing = clusters.get(key)
        if existing:
            existing["sources"] = _clean_sources(existing.get("sources", []) + item.get("sources", []))
            existing["evidence"] = _clean_sources(existing.get("evidence", []) + item.get("evidence", []))
            existing["evidence_rank"] = max(existing.get("evidence_rank") or 0, item.get("evidence_rank") or 0)
            if not existing.get("mondo_id") and item.get("mondo_id"):
                existing["mondo_id"] = item.get("mondo_id")
            continue
        clusters[key] = dict(item)
        order.append(key)

    return [clusters[key] for key in order]
