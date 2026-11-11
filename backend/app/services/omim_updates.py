"""OMIM Update List (omim.org/statistics/updates) as a source of change dates.

The extension reads the monthly Update List pages (New / Updated Entries,
New / Updated Clinical Synopses) and sends the rows here. With those dates
the work queue can tell which platform content OMIM changed afterwards:

* content date — the OMIM version a piece of platform content reflects:
  an adopted capture → the capture day (US Eastern) or the page's Edit
  History date, whichever is later; OMIM.xlsx text → the workbook Done date;
* a MIM listed in the Update List after its content date needs a new capture.

Months are recorded so the page can show which part of the history is
covered; a month counts as complete once it was read after it ended.
"""
from __future__ import annotations

import calendar
import re
from datetime import date, datetime, timedelta, timezone

from . import omim_db

KINDS = {
    ("new", "entries"): "new_entry",
    ("updated", "entries"): "updated_entry",
    ("new", "clinical synopses"): "new_synopsis",
    ("updated", "clinical synopses"): "updated_synopsis",
}
KIND_LABELS = {
    "new_entry": "新 entry",
    "updated_entry": "更新 entry",
    "new_synopsis": "新 Clinical Synopsis",
    "updated_synopsis": "更新 Clinical Synopsis",
}
ENTRY_KINDS = {"new_entry", "updated_entry"}
FIRST_MONTH = "1995-12"
STALE_MODES = {"both", "updates_only"}
MAX_ROWS = 5000

_MONTHS = {name.lower(): i for i, name in enumerate(calendar.month_name) if name}
_DATE_RE = re.compile(
    r"^\s*(January|February|March|April|May|June|July|August|September|October|November|December)"
    r"\s+(\d{1,2})(?:st|nd|rd|th)?,?\s+(\d{4})\s*$",
    re.IGNORECASE,
)
_SECTION_RE = re.compile(r"^\s*(New|Updated)\s+(Entries|Clinical Synopses)\s*:?\s*$", re.IGNORECASE)
# "*108731", "[*108731](https://omim.org/entry/108731)", "*108731\tTITLE", "619701"
_ROW_RE = re.compile(
    r"^\s*\[?\s*([#%*+^])?\s*(\d{6})\s*\]?(?:\((?:[^)]*/(entry|clinicalSynopsis)/\d{6}[^)]*)\))?\s*(.*)$"
)
_MONTH_RE = re.compile(r"^\d{4}-\d{2}$")
_DAY_RE = re.compile(r"^\d{4}-\d{2}-\d{2}$")


# --------------------------------------------------------------------------
# Dates
# --------------------------------------------------------------------------

def us_date(iso: str) -> str:
    """Calendar date in US Eastern time (OMIM's) for a UTC ISO timestamp."""
    if not iso:
        return ""
    try:
        dt = datetime.fromisoformat(str(iso).replace("Z", "+00:00"))
    except ValueError:
        return str(iso)[:10]
    if dt.tzinfo is None:
        dt = dt.replace(tzinfo=timezone.utc)
    return (dt.astimezone(timezone.utc) - timedelta(hours=5)).date().isoformat()


def us_today() -> str:
    return us_date(datetime.now(timezone.utc).isoformat())


def parse_date(text: str) -> str:
    match = _DATE_RE.match(str(text or ""))
    if not match:
        return ""
    month = _MONTHS[match.group(1).lower()]
    try:
        return date(int(match.group(3)), month, int(match.group(2))).isoformat()
    except ValueError:
        return ""


def month_of(day: str) -> str:
    return str(day or "")[:7]


def months_between(start: str, end: str) -> list[str]:
    """Inclusive list of YYYY-MM from ``start`` to ``end``."""
    if not (_MONTH_RE.match(start or "") and _MONTH_RE.match(end or "")) or start > end:
        return []
    year, month = int(start[:4]), int(start[5:])
    out = []
    while f"{year:04d}-{month:02d}" <= end:
        out.append(f"{year:04d}-{month:02d}")
        month += 1
        if month > 12:
            year, month = year + 1, 1
    return out


def _previous_month(month: str) -> str:
    year, mon = int(month[:4]), int(month[5:])
    return f"{year - 1:04d}-12" if mon == 1 else f"{year:04d}-{mon - 1:02d}"


def _month_end(month: str) -> str:
    year, mon = int(month[:4]), int(month[5:])
    return date(year, mon, calendar.monthrange(year, mon)[1]).isoformat()


# --------------------------------------------------------------------------
# Parsing
# --------------------------------------------------------------------------

def _kind(section: str, link_kind: str = "") -> str:
    match = _SECTION_RE.match(str(section or ""))
    if not match:
        return ""
    return KINDS.get((match.group(1).lower(), match.group(2).lower()), "")


def normalize_rows(rows: list) -> list[dict]:
    """Rows extracted by the extension from the page DOM."""
    out = []
    for row in rows[:MAX_ROWS] if isinstance(rows, list) else []:
        if not isinstance(row, dict):
            continue
        mim = str(row.get("mim") or "").strip()
        listed = parse_date(row.get("date") or "")
        kind = _kind(row.get("section") or "")
        if not (re.match(r"^\d{6}$", mim) and listed and kind):
            continue
        prefix = str(row.get("prefix") or "").strip()
        out.append({
            "mim": mim,
            "kind": kind,
            "listed_date": listed,
            "prefix": prefix if prefix in {"#", "%", "*", "+", "^"} else "",
            "title": " ".join(str(row.get("title") or "").split())[:300],
        })
    return out


def parse_text(text: str) -> list[dict]:
    """Fallback: rows from the page's rendered text (or a copied page)."""
    lines = [ln.strip() for ln in str(text or "").replace("\r", "").split("\n")]
    out: list[dict] = []
    listed = kind = ""
    pending: dict | None = None
    for line in lines:
        if not line:
            continue
        day = parse_date(line)
        if day:
            listed, kind, pending = day, "", None
            continue
        section = _kind(line)
        if section:
            kind, pending = section, None
            continue
        if not (listed and kind):
            continue
        match = _ROW_RE.match(line)
        if match:
            pending = {
                "mim": match.group(2), "kind": kind, "listed_date": listed,
                "prefix": match.group(1) or "", "title": match.group(4).strip(),
            }
            out.append(pending)
            if pending["title"]:
                pending = None
            continue
        if pending is not None and not pending["title"]:
            pending["title"] = line[:300]
            pending = None
        elif line.startswith(("NOTE:", "OMIM®", "Copyright")):
            listed = kind = ""
    return out


# --------------------------------------------------------------------------
# Storage
# --------------------------------------------------------------------------

def ingest(year: int, month: int, rows: list | None = None, text: str = "", *, user: str = "",
           path=None) -> dict:
    """Store one month's Update List; returns counts per kind."""
    ym = f"{int(year):04d}-{int(month):02d}"
    if not _MONTH_RE.match(ym) or not (FIRST_MONTH <= ym <= us_today()[:7]) or not 1 <= int(month) <= 12:
        raise ValueError("年月不在 OMIM Update List 的範圍內")
    items = normalize_rows(rows or [])
    if not items and text:
        items = parse_text(text)
    items = [i for i in items if month_of(i["listed_date"]) == ym]
    if not items:
        raise ValueError(f"{ym} 的頁面讀不到任何更新項目")
    counts: dict[str, int] = {}
    now = datetime.now(timezone.utc).isoformat(timespec="seconds")
    complete = 1 if us_today() > _month_end(ym) else 0
    with omim_db._write_lock, omim_db.connect(path, write=True) as con:
        for item in items:
            con.execute(
                "INSERT INTO omim_updates(mim, kind, listed_date, prefix, title) VALUES (?,?,?,?,?) "
                "ON CONFLICT(mim, kind, listed_date) DO UPDATE SET prefix=excluded.prefix, title=excluded.title",
                (item["mim"], item["kind"], item["listed_date"], item["prefix"], item["title"]),
            )
            counts[item["kind"]] = counts.get(item["kind"], 0) + 1
        con.execute(
            "INSERT INTO omim_update_months(month, captured_at, captured_by, rows, complete) VALUES (?,?,?,?,?) "
            "ON CONFLICT(month) DO UPDATE SET captured_at=excluded.captured_at, captured_by=excluded.captured_by, "
            "rows=excluded.rows, complete=MAX(omim_update_months.complete, excluded.complete)",
            (ym, now, user, len(items), complete),
        )
    return {"month": ym, "rows": len(items), "counts": counts, "complete": bool(complete)}


def latest_dates(path=None) -> tuple[dict[str, str], dict[str, str]]:
    """(entry updates, any update incl. synopses): MIM -> latest listed date."""
    entry: dict[str, str] = {}
    anything: dict[str, str] = {}
    if not omim_db.exists(path):
        return entry, anything
    with omim_db.connect(path) as con:
        for r in con.execute("SELECT mim, kind, MAX(listed_date) AS d FROM omim_updates GROUP BY mim, kind"):
            if r["d"] > anything.get(r["mim"], ""):
                anything[r["mim"]] = r["d"]
            if r["kind"] in ENTRY_KINDS and r["d"] > entry.get(r["mim"], ""):
                entry[r["mim"]] = r["d"]
    return entry, anything


def get_setting(key: str, default: str = "", path=None) -> str:
    if not omim_db.exists(path):
        return default
    with omim_db.connect(path) as con:
        row = con.execute("SELECT value FROM meta WHERE key = ?", (f"setting:{key}",)).fetchone()
    return row["value"] if row else default


def set_setting(key: str, value: str, path=None) -> None:
    with omim_db._write_lock, omim_db.connect(path, write=True) as con:
        con.execute(
            "INSERT INTO meta(key, value) VALUES (?, ?) ON CONFLICT(key) DO UPDATE SET value = excluded.value",
            (f"setting:{key}", value),
        )


def stale_mode(path=None) -> str:
    mode = get_setting("stale_mode", "both", path)
    return mode if mode in STALE_MODES else "both"


def coverage(path=None) -> dict:
    """Which Update List months are stored, and which are still needed."""
    current = us_today()[:7]
    out = {"needed_from": "", "current_month": current, "months": [], "missing": [], "untitled": [],
           "coverage_start": "", "rows": 0, "last_listed_date": "", "last_captured_at": ""}
    if not omim_db.exists(path):
        return out
    with omim_db.connect(path) as con:
        oldest = con.execute(
            "SELECT MIN(legacy_done) FROM genes WHERE legacy_done GLOB '[12][0-9][0-9][0-9]-[01][0-9]*'"
        ).fetchone()[0]
        months = {r["month"]: dict(r) for r in con.execute("SELECT * FROM omim_update_months")}
        out["rows"] = con.execute("SELECT COUNT(*) FROM omim_updates").fetchone()[0]
        out["last_listed_date"] = con.execute("SELECT MAX(listed_date) FROM omim_updates").fetchone()[0] or ""
        # Months stored before titles were read (extension < 0.6.1).
        untitled = [r[0] for r in con.execute(
            "SELECT DISTINCT substr(listed_date, 1, 7) FROM omim_updates WHERE title = '' "
            "AND prefix != '' ORDER BY 1 DESC")]
    needed = max(FIRST_MONTH, (oldest or "")[:7]) if oldest else ""
    if not needed:
        year, mon = int(current[:4]), int(current[5:])
        needed = f"{year - 1:04d}-{mon:02d}"
    out["needed_from"] = needed
    out["months"] = sorted(months.values(), key=lambda m: m["month"], reverse=True)
    out["last_captured_at"] = max((m["captured_at"] for m in months.values()), default="")
    out["missing"] = [m for m in reversed(months_between(needed, current))
                      if not (months.get(m) or {}).get("complete")]
    # Earliest month from which every month up to the last full one is complete.
    start = ""
    month = _previous_month(current)
    while (months.get(month) or {}).get("complete"):
        start = month
        if month <= FIRST_MONTH:
            break
        month = _previous_month(month)
    out["coverage_start"] = start
    out["untitled"] = untitled
    out["history_complete"] = bool(start) and start <= needed
    return out


GENE_PREFIXES = {"*", "+"}
DISEASE_PREFIXES = {"#", "%"}


def missing_entries(path=None) -> list[dict]:
    """Update List MIMs the platform does not carry at all.

    * gene entries (``*``/``+``) not among the platform genes;
    * disease entries (``#``/``%``) no platform gene carries (e.g. a disease
      newly linked to a gene whose workbook row has no disease).

    The Update List does not name a disease's gene, so each becomes its own
    queue item; adopting the page links it from the page's relationships
    table. MIMs later moved or removed (``^``) and rows whose entry type is
    unknown are left out. Each row: mim, prefix, cls (gene/disease),
    listed_date (latest), first_date, title.
    """
    if not omim_db.exists(path):
        return []
    with omim_db.connect(path) as con:
        rows = con.execute(
            "SELECT mim, MAX(prefix) AS prefix, MAX(listed_date) AS listed_date, MIN(listed_date) AS first_date, "
            "MAX(title) AS title, MAX(kind = 'new_entry') AS is_new FROM omim_updates "
            "WHERE mim NOT IN (SELECT mim FROM omim_updates WHERE prefix = '^') "
            "GROUP BY mim"
        ).fetchall()
        genes = {r[0] for r in con.execute("SELECT gene_mim FROM genes")}
        diseases = {r[0] for r in con.execute(
            "SELECT DISTINCT phenotype_mim FROM gene_diseases WHERE active = 1 AND phenotype_mim != ''")}
    out = []
    for r in rows:
        prefix = r["prefix"] or ""
        if prefix in GENE_PREFIXES and r["mim"] not in genes:
            cls = "gene"
        elif prefix in DISEASE_PREFIXES and r["mim"] not in diseases:
            cls = "disease"
        else:
            continue
        out.append({**dict(r), "cls": cls, "is_new": bool(r["is_new"])})
    out.sort(key=lambda e: (e["listed_date"], e["mim"]), reverse=True)
    return out


def is_stale(captured_at: str, mim: str, latest_any: dict[str, str]) -> str:
    """The OMIM update date when a capture predates it, else ''."""
    listed = latest_any.get(mim or "", "")
    return listed if listed and listed > us_date(captured_at) else ""
