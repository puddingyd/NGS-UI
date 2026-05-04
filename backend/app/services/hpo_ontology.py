"""hp.obo parser + in-memory HPO term search.

Loaded once at startup (~17 000 terms, ~10 MB on disk → a few MB in
memory). Search ranks results by:
  1. Exact ID match  (HP:0001250 → that term first)
  2. Exact name match
  3. Name starts-with the query
  4. Name contains the query
  5. Any synonym contains the query
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable

from ..config import REPO_ROOT

HP_OBO_PATH = REPO_ROOT / "phenotype_data" / "hp.obo"

_ID_RE = re.compile(r"^id:\s*(HP:\d+)\s*$")
_NAME_RE = re.compile(r"^name:\s*(.+?)\s*$")
_DEF_RE = re.compile(r'^def:\s*"([^"]+)"')
_SYN_RE = re.compile(r'^synonym:\s*"([^"]+)"')
_OBSOLETE_RE = re.compile(r"^is_obsolete:\s*true\s*$")


@dataclass
class HpoTerm:
    id: str
    name: str
    synonyms: list[str] = field(default_factory=list)
    definition: str = ""

    def to_dict(self) -> dict:
        return {
            "hpo_id": self.id,
            "name": self.name,
            "synonyms": self.synonyms,
            "definition": self.definition,
        }


_TERMS: dict[str, HpoTerm] = {}
_NAME_INDEX: list[tuple[str, str]] = []  # (lowercased_name, hpo_id)
_SYN_INDEX:  list[tuple[str, str]] = []  # (lowercased_synonym, hpo_id)


def _parse_obo(path: Path) -> dict[str, HpoTerm]:
    terms: dict[str, HpoTerm] = {}
    cur: HpoTerm | None = None
    in_term = False
    obsolete = False
    with path.open("r", encoding="utf-8") as f:
        for line in f:
            line = line.rstrip("\n")
            if line == "[Term]":
                if cur and not obsolete and cur.id and cur.name:
                    terms[cur.id] = cur
                cur = HpoTerm(id="", name="")
                in_term = True
                obsolete = False
                continue
            if line.startswith("[") and line != "[Term]":
                if cur and not obsolete and cur.id and cur.name:
                    terms[cur.id] = cur
                cur = None
                in_term = False
                continue
            if not in_term or cur is None:
                continue
            if _OBSOLETE_RE.match(line):
                obsolete = True
                continue
            m = _ID_RE.match(line)
            if m:
                cur.id = m.group(1); continue
            m = _NAME_RE.match(line)
            if m:
                cur.name = m.group(1); continue
            m = _DEF_RE.match(line)
            if m:
                cur.definition = m.group(1); continue
            m = _SYN_RE.match(line)
            if m:
                cur.synonyms.append(m.group(1)); continue
    if cur and not obsolete and cur.id and cur.name:
        terms[cur.id] = cur
    return terms


def _build_indexes(terms: dict[str, HpoTerm]) -> None:
    _NAME_INDEX.clear()
    _SYN_INDEX.clear()
    for t in terms.values():
        _NAME_INDEX.append((t.name.lower(), t.id))
        for syn in t.synonyms:
            _SYN_INDEX.append((syn.lower(), t.id))


def load(path: Path = HP_OBO_PATH) -> int:
    """Idempotent loader. Returns the term count."""
    if _TERMS:
        return len(_TERMS)
    if not path.exists():
        return 0
    parsed = _parse_obo(path)
    _TERMS.update(parsed)
    _build_indexes(_TERMS)
    return len(_TERMS)


def get(hpo_id: str) -> HpoTerm | None:
    if not _TERMS:
        load()
    return _TERMS.get(hpo_id)


def search(query: str, limit: int = 20) -> list[dict]:
    """Rank-order search; see module docstring for ranking."""
    if not _TERMS:
        load()
    q = (query or "").strip()
    if not q:
        return []
    q_lower = q.lower()
    seen: set[str] = set()
    out: list[HpoTerm] = []

    # 1. Exact ID
    if q.upper().startswith("HP:") and q.upper() in _TERMS:
        out.append(_TERMS[q.upper()]); seen.add(q.upper())

    # 2. Exact name
    for name_lc, hid in _NAME_INDEX:
        if hid in seen: continue
        if name_lc == q_lower:
            out.append(_TERMS[hid]); seen.add(hid)

    # 3. Name starts-with
    for name_lc, hid in _NAME_INDEX:
        if hid in seen: continue
        if name_lc.startswith(q_lower):
            out.append(_TERMS[hid]); seen.add(hid)
            if len(out) >= limit: break

    # 4. Name contains
    if len(out) < limit:
        for name_lc, hid in _NAME_INDEX:
            if hid in seen: continue
            if q_lower in name_lc:
                out.append(_TERMS[hid]); seen.add(hid)
                if len(out) >= limit: break

    # 5. Synonym contains
    if len(out) < limit:
        for syn_lc, hid in _SYN_INDEX:
            if hid in seen: continue
            if q_lower in syn_lc:
                out.append(_TERMS[hid]); seen.add(hid)
                if len(out) >= limit: break

    return [t.to_dict() for t in out[:limit]]
