"""hp.obo parser + in-memory HPO term search.

Loaded once at startup (~17 000 terms, ~10 MB on disk → a few MB in
memory). Search ranks results by:
  1. Exact ID match  (HP:0001250 → that term first)
  2. Exact name match
  3. Name starts-with the query
  4. Name contains the query
  5. Any synonym contains the query
  6. Fuzzy name / synonym match (small spelling differences)

Each result also carries its immediate ``is_a`` parent terms so every HPO
picker can offer the same parent shortcut as the standalone phenotype tool.
"""
from __future__ import annotations

import re
import threading
import unicodedata
from dataclasses import dataclass, field
from difflib import SequenceMatcher
from pathlib import Path

from ..config import PHENO_DATA_DIR

HP_OBO_PATH = PHENO_DATA_DIR / "hp.obo"

_ID_RE = re.compile(r"^id:\s*(HP:\d+)\s*$")
_NAME_RE = re.compile(r"^name:\s*(.+?)\s*$")
_DEF_RE = re.compile(r'^def:\s*"([^"]+)"')
_SYN_RE = re.compile(r'^synonym:\s*"([^"]+)"')
_IS_A_RE = re.compile(r"^is_a:\s*(HP:\d+)\b")
_OBSOLETE_RE = re.compile(r"^is_obsolete:\s*true\s*$")
_SEARCH_SEP_RE = re.compile(r"[^a-z0-9]+")


@dataclass
class HpoTerm:
    id: str
    name: str
    synonyms: list[str] = field(default_factory=list)
    definition: str = ""
    parent_ids: list[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        parents = []
        for parent_id in self.parent_ids:
            parent = _TERMS.get(parent_id)
            if parent is not None:
                parents.append({"hpo_id": parent.id, "name": parent.name})
        return {
            "hpo_id": self.id,
            "name": self.name,
            "synonyms": self.synonyms,
            "definition": self.definition,
            "parents": parents,
        }


_TERMS: dict[str, HpoTerm] = {}
_NAME_INDEX: list[tuple[str, str]] = []  # (lowercased_name, hpo_id)
_SYN_INDEX:  list[tuple[str, str]] = []  # (lowercased_synonym, hpo_id)
_SEARCH_TEXTS: dict[str, tuple[str, ...]] = {}
_GRAM_INDEX: dict[str, list[str]] = {}
_LOAD_LOCK = threading.Lock()


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
            m = _IS_A_RE.match(line)
            if m and m.group(1) not in cur.parent_ids:
                cur.parent_ids.append(m.group(1)); continue
    if cur and not obsolete and cur.id and cur.name:
        terms[cur.id] = cur
    return terms


def _build_indexes(terms: dict[str, HpoTerm]) -> None:
    _NAME_INDEX.clear()
    _SYN_INDEX.clear()
    _SEARCH_TEXTS.clear()
    _GRAM_INDEX.clear()
    for t in terms.values():
        _NAME_INDEX.append((t.name.lower(), t.id))
        for syn in t.synonyms:
            _SYN_INDEX.append((syn.lower(), t.id))
        texts = tuple(dict.fromkeys(
            normalized
            for raw in (t.name, *t.synonyms)
            if (normalized := _normalize_search_text(raw))
        ))
        _SEARCH_TEXTS[t.id] = texts
        term_grams: set[str] = set()
        for text in texts:
            term_grams.update(_search_grams(text))
        for gram in term_grams:
            _GRAM_INDEX.setdefault(gram, []).append(t.id)


def _normalize_search_text(value: str) -> str:
    """Normalize punctuation/spacing while preserving fuzzy-search words."""
    folded = unicodedata.normalize("NFKD", value or "").casefold()
    ascii_text = "".join(ch for ch in folded if not unicodedata.combining(ch))
    return " ".join(part for part in _SEARCH_SEP_RE.split(ascii_text) if part)


def _search_grams(value: str) -> set[str]:
    compact = (value or "").replace(" ", "")
    if len(compact) < 3:
        return {compact} if compact else set()
    return {compact[i:i + 3] for i in range(len(compact) - 2)}


def _fuzzy_similarity(query: str, candidate: str) -> float:
    """Approximate Fuse.js-style matching for names and synonyms.

    Full-string similarity handles misspelled phrases; token similarity lets a
    misspelled word still find a longer HPO label (for example ``seizuer`` →
    ``Febrile seizure``).
    """
    if not query or not candidate:
        return 0.0
    full = SequenceMatcher(None, query, candidate).ratio()
    query_tokens = query.split()
    candidate_tokens = candidate.split()
    if not query_tokens or not candidate_tokens:
        return full
    token_scores = [
        max(SequenceMatcher(None, q_token, c_token).ratio() for c_token in candidate_tokens)
        for q_token in query_tokens
    ]
    token_score = sum(token_scores) / len(token_scores)
    # Prefer the concise direct term when several longer labels contain the
    # same nearly-matching word, while still allowing a typo to find a longer
    # phrase such as "Febrile seizure".
    specificity = min(1.0, len(query) / len(candidate))
    return max(full, (0.85 * token_score) + (0.15 * specificity))


def load(path: Path = HP_OBO_PATH) -> int:
    """Idempotent loader. Returns the term count."""
    if _TERMS:
        return len(_TERMS)
    with _LOAD_LOCK:
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


def _normalize_hpo_id(q: str) -> str | None:
    """If the query looks like an HPO ID fragment, return canonical
    'HP:NNNNNNN'. Handles bare digits and HP:-prefixed forms with or
    without leading zeros: '1250' / 'HP:1250' / 'hp:0001250' /
    '0001250' all → 'HP:0001250'. Returns None when the query has any
    non-digit content after stripping 'HP:'.
    """
    s = (q or "").strip().upper()
    if s.startswith("HP:"):
        s = s[3:]
    if s and s.isdigit() and len(s) <= 7:
        return "HP:" + s.zfill(7)
    return None


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

    # 1. Exact ID (canonical form 'HP:NNNNNNN' or any fragment that
    #    normalises to one — bare digits, missing leading zeros, …).
    norm_id = _normalize_hpo_id(q)
    if norm_id and norm_id in _TERMS:
        out.append(_TERMS[norm_id]); seen.add(norm_id)
    if q.upper().startswith("HP:") and q.upper() in _TERMS and q.upper() not in seen:
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

    # 6. ID substring — for digit-only queries that didn't resolve via
    #    zero-pad (e.g. '125' should also surface HP:0012500 / 0125000).
    if q.isdigit() and len(out) < limit:
        for hid, t in _TERMS.items():
            if hid in seen: continue
            if q in hid:
                out.append(t); seen.add(hid)
                if len(out) >= limit: break

    # 7. Fuzzy name / synonym fallback. Exact/substring matches stay ahead of
    # fuzzy matches, while misspellings such as "seizuer" can still surface
    # "Seizure". A trigram index keeps the candidate set small enough for an
    # interactive typeahead without adding a third-party Python dependency.
    normalized_query = _normalize_search_text(q)
    if len(out) < limit and len(normalized_query.replace(" ", "")) >= 3:
        candidates: set[str] = set()
        for gram in _search_grams(normalized_query):
            candidates.update(_GRAM_INDEX.get(gram, ()))
        ranked: list[tuple[float, str, str]] = []
        for hid in candidates:
            if hid in seen:
                continue
            score = max(
                (_fuzzy_similarity(normalized_query, text) for text in _SEARCH_TEXTS.get(hid, ())),
                default=0.0,
            )
            if score >= 0.70:
                ranked.append((-score, _TERMS[hid].name.casefold(), hid))
        ranked.sort()
        for _negative_score, _name, hid in ranked:
            out.append(_TERMS[hid]); seen.add(hid)
            if len(out) >= limit:
                break

    return [t.to_dict() for t in out[:limit]]
