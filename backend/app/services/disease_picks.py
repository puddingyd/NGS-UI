"""Stable reviewer selection of SNV OMIM diseases for reports.

Historically the reviewer's disease ticks were stored as slot numbers
(``edits["report_diseases"] = {"2": True}``) and every consumer re-read
``Disease{N}`` at export time. Any change to the OMIM workbook order then
silently pointed an already-reviewed case at a different disease.

The current format stores a content key plus a label snapshot::

    edits["report_disease_picks"] = {
        "610125|microphthalmia syndromic 5": {
            "label": "Microphthalmia, syndromic 5 (610125)(AD)",
            "slot": 1,          # informational only: slot at pick time
        },
    }

Resolution against the *current* Disease1..N cells, in order:

1. exact key (phenotype MIM + normalized name);
2. same phenotype MIM when exactly one current cell carries it (rename);
3. same normalized name when exactly one current cell has it (MIM added);
4. otherwise the snapshot label is used, so the report still names the
   disease the reviewer actually ticked.

When ``report_disease_picks`` is absent the legacy slot behaviour is kept
unchanged (old slot numbers are translated through ``Disease_legacy_slots``
when the platform OMIM database renders diseases in a different order). The frontend (``frontend/app.js`` ``diseasePickKey``) must
compute the same key; ``tests/test_disease_picks.py`` pins both.
"""
from __future__ import annotations

import re
from typing import Iterable

from . import omim_store

PICKS_FIELD = "report_disease_picks"
LEGACY_FIELD = "report_diseases"

_MIM_RE = re.compile(r"\((\d{6})\)")
_TRAILING_PAREN_RE = re.compile(r"\s*\([^()]*\)\s*$")
_NON_ALNUM_RE = re.compile(r"[^a-z0-9]+")


def first_line(text) -> str:
    value = str(text or "").strip()
    return value.splitlines()[0].strip() if value else ""


def name_key(text) -> str:
    """Normalized disease name: first line without trailing (...) groups."""
    line = first_line(text)
    previous = None
    while previous != line:
        previous = line
        line = _TRAILING_PAREN_RE.sub("", line)
    line = line.lower().replace("&", " and ")
    return _NON_ALNUM_RE.sub(" ", line).strip()


def phenotype_mim(text) -> str:
    match = _MIM_RE.search(first_line(text))
    return match.group(1) if match else ""


def pick_key(text) -> str:
    return f"{phenotype_mim(text)}|{name_key(text)}"


def disease_cells(record: dict | None) -> list[tuple[int, str]]:
    """(slot, text) for every non-empty Disease1..N cell, in slot order."""
    out: list[tuple[int, str]] = []
    for idx, field in enumerate(omim_store.DISEASE_FIELDS, start=1):
        text = str((record or {}).get(field) or "").strip()
        if text and text != "NA":
            out.append((idx, text))
    return out


def has_picks(edits: dict | None) -> bool:
    return isinstance((edits or {}).get(PICKS_FIELD), dict)


def _legacy_slots(edits: dict) -> list[int]:
    picked = edits.get(LEGACY_FIELD) or {}
    if not isinstance(picked, dict):
        return []
    slots = []
    for raw_key, value in picked.items():
        if not value:
            continue
        try:
            slot = int(raw_key)
        except (TypeError, ValueError):
            continue
        if 1 <= slot <= omim_store.DISEASE_SLOT_COUNT:
            slots.append(slot)
    return sorted(set(slots))


def resolve(record: dict | None, edits: dict | None) -> dict:
    """Resolve reviewer picks against the current Disease cells.

    Returns ``{"mode", "matched": [(slot, text, key)], "unmatched":
    [(key, label)]}``; ``matched`` is in current slot order.
    """
    edits = edits or {}
    cells = disease_cells(record)
    if not has_picks(edits):
        wanted = set(_legacy_slots(edits))
        # The platform OMIM database may render diseases in a different
        # order than the old workbook; translate old slot numbers.
        legacy_map = (record or {}).get("Disease_legacy_slots") or (record or {}).get(
            "OMIM_legacy_slots"
        )
        if isinstance(legacy_map, dict) and legacy_map:
            wanted = {int(legacy_map[str(s)]) for s in wanted if str(s) in legacy_map}
        matched = [
            (slot, text, pick_key(text)) for slot, text in cells if slot in wanted
        ]
        return {"mode": "legacy", "matched": matched, "unmatched": []}

    picks = edits.get(PICKS_FIELD) or {}
    cell_keys = [(slot, text, pick_key(text)) for slot, text in cells]
    used: set[int] = set()
    matched: list[tuple[int, str, str]] = []
    unmatched: list[tuple[str, str]] = []

    def unique(predicate) -> tuple | None:
        hits = [c for c in cell_keys if c[0] not in used and predicate(c)]
        return hits[0] if len(hits) == 1 else None

    for key, snapshot in picks.items():
        if not snapshot:
            continue
        label = ""
        if isinstance(snapshot, dict):
            label = str(snapshot.get("label") or "").strip()
        key = str(key or "")
        mim, _, nkey = key.partition("|")
        hit = unique(lambda c: c[2] == key)
        if hit is None and mim:
            hit = unique(lambda c: phenotype_mim(c[1]) == mim)
        if hit is None and nkey:
            hit = unique(lambda c: name_key(c[1]) == nkey)
        if hit is None:
            unmatched.append((key, label))
            continue
        used.add(hit[0])
        matched.append((hit[0], hit[1], key))
    matched.sort(key=lambda item: item[0])
    return {"mode": "keys", "matched": matched, "unmatched": unmatched}


def selected_texts(
    record: dict | None,
    edits: dict | None,
    *,
    fallback_first: bool = True,
) -> list[str]:
    """Disease texts to report: matched cells, then unmatched snapshots.

    With nothing ticked, ``fallback_first`` returns the first disease cell
    (the long-standing DOCX behaviour).
    """
    result = resolve(record, edits)
    texts = [text for _slot, text, _key in result["matched"]]
    texts.extend(label for _key, label in result["unmatched"] if label)
    if texts or not fallback_first:
        return texts
    cells = disease_cells(record)
    return [cells[0][1]] if cells else []


def picks_from_slots(record: dict | None, slots: Iterable[int]) -> dict:
    """Build the new-format picks from legacy slot numbers (migration)."""
    wanted = {int(s) for s in slots}
    legacy_map = (record or {}).get("Disease_legacy_slots") or (record or {}).get("OMIM_legacy_slots")
    if isinstance(legacy_map, dict) and legacy_map:
        wanted = {int(legacy_map[str(s)]) for s in wanted if str(s) in legacy_map}
    picks: dict[str, dict] = {}
    for slot, text in disease_cells(record):
        if slot in wanted:
            picks[pick_key(text)] = {"label": first_line(text), "slot": slot}
    return picks
