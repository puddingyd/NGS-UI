"""Patient roster built from the hospital's 未完成報告清單 xlsx.

The lab exports a weekly "未完成報告明細" spreadsheet listing the
specimen ID (檢體編號), MRN (病歷號), patient name (姓名), test name
(檢驗名稱) and ordering department (科別) for every pending case.

Uploading it here:
  1. archives the raw xlsx under patient_list/{ts}_{name}.xlsx
  2. merges each row into patient_list/roster.json keyed by LIS_ID

LIS_ID is the specimen ID with its fixed "8BB1" prefix stripped
(清單 "8BB126WE0092" → LIS_ID "26WE0092").

The 載入新個案 modal then reads the roster to auto-fill MRN / 姓名 /
Test type when the reviewer picks a pipeline-dropped TSV — replacing
the old "parse the MRN out of the phenotype.txt filename" trick.
Merge semantics are additive: a fresh weekly upload adds new LIS_IDs
and updates changed fields, but never removes a LIS_ID that's no
longer in the latest list (so already-handled cases keep their
roster entry).
"""
from __future__ import annotations

import json
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from ..config import PATIENT_LIST_DIR

_ROSTER_NAME = "roster.json"
_UPLOADS_NAME = "uploads.json"
_SPECIMEN_PREFIX = "8BB1"
_HEADER_KEY = "檢體編號"     # the cell that marks the header row in col 0
_UI_ALIAS_SUFFIXES = ("-dragen", "-nckuh", "-inhouse")


def _roster_path() -> Path:
    return PATIENT_LIST_DIR / _ROSTER_NAME


def _uploads_path() -> Path:
    return PATIENT_LIST_DIR / _UPLOADS_NAME


def list_uploads() -> list[dict]:
    """Read patient_list/uploads.json — append-only log of every
    successful ingest. Latest first. Empty when nothing uploaded yet.
    """
    p = _uploads_path()
    if not p.is_file():
        return []
    try:
        data = json.loads(p.read_text(encoding="utf-8"))
        if not isinstance(data, list):
            return []
        return sorted(data, key=lambda r: r.get("uploaded_at", ""), reverse=True)
    except (json.JSONDecodeError, OSError):
        return []


def _append_upload(record: dict) -> None:
    p = _uploads_path()
    cur: list[dict] = []
    if p.is_file():
        try:
            data = json.loads(p.read_text(encoding="utf-8"))
            if isinstance(data, list):
                cur = data
        except (json.JSONDecodeError, OSError):
            cur = []
    cur.append(record)
    p.write_text(json.dumps(cur, ensure_ascii=False, indent=2),
                 encoding="utf-8")


def _now_iso() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def _strip(v: Any) -> str:
    return str(v).strip() if v is not None else ""


def strip_ui_alias_suffix(sample_id: str) -> str:
    """Return the LIS_ID candidate behind a caller-disambiguated UI ID.

    Tertiary jobs append -dragen / -nckuh / -inhouse to avoid output
    collisions. The clinic roster is keyed by the original LIS_ID, so
    roster lookups may need this caller suffix removed. Exact roster
    matches still win before this fallback is used.
    """
    sid = _strip(sample_id)
    sid_l = sid.lower()
    for suffix in _UI_ALIAS_SUFFIXES:
        if sid_l.endswith(suffix):
            return sid[: -len(suffix)]
    return sid


def lookup_candidates(sample_id: str, *aliases: str) -> list[str]:
    """Ordered roster keys to try for a UI sample ID and source aliases."""
    out: list[str] = []
    for raw in (sample_id, *aliases):
        sid = _strip(raw)
        if not sid:
            continue
        for candidate in (sid, strip_ui_alias_suffix(sid)):
            if candidate and candidate not in out:
                out.append(candidate)
    return out


def _lis_id_from_specimen(specimen: str) -> str:
    """8BB126WE0092 → 26WE0092. Leaves non-prefixed values as-is."""
    s = _strip(specimen)
    if s.startswith(_SPECIMEN_PREFIX):
        return s[len(_SPECIMEN_PREFIX):]
    return s


def _test_type_from_name(test_name: str) -> str:
    """檢驗名稱 → WES / WGS (best-effort; defaults WES)."""
    n = (test_name or "").upper()
    if "WGS" in n:
        return "WGS"
    return "WES"


def parse_xlsx(path: Path) -> list[dict]:
    """Pull the data rows out of the 未完成報告清單 xlsx.

    Returns a list of {lis_id, mrn, name, test_name, test_type,
    department, specimen}. Robust to the variable header-row position
    (the file has a few title rows before the real header): we find
    the row whose first cell is 檢體編號, then take subsequent rows
    whose first cell looks like a real specimen ID.
    """
    import openpyxl
    wb = openpyxl.load_workbook(str(path), read_only=True, data_only=True)
    ws = wb.worksheets[0]
    rows = list(ws.iter_rows(values_only=True))
    wb.close()

    # Locate the header row + the column indices we care about.
    header_idx = None
    col = {}
    for i, r in enumerate(rows):
        cells = [_strip(c) for c in (r or [])]
        if cells and cells[0] == _HEADER_KEY:
            header_idx = i
            for j, name in enumerate(cells):
                col[name] = j
            break
    if header_idx is None:
        return []

    def get(r: tuple, name: str) -> str:
        j = col.get(name)
        if j is None or j >= len(r):
            return ""
        return _strip(r[j])

    out: list[dict] = []
    seen: set[str] = set()
    for r in rows[header_idx + 1:]:
        if not r:
            continue
        specimen = get(r, _HEADER_KEY)
        if not specimen:
            continue
        lis_id = _lis_id_from_specimen(specimen)
        if not lis_id or lis_id in seen:
            # Dedupe by LIS_ID — a patient can appear on multiple rows
            # (multiple 醫令). First occurrence wins.
            continue
        seen.add(lis_id)
        test_name = get(r, "檢驗名稱")
        out.append({
            "lis_id":     lis_id,
            "specimen":   specimen,
            "mrn":        get(r, "病歷號"),
            "name":       get(r, "姓名"),
            "test_name":  test_name,
            "test_type":  _test_type_from_name(test_name),
            "department": get(r, "科別"),
            "physician":  get(r, "開單醫師"),
            # 簽收時間 — when LIS booked the sample in. Used as the
            # "日期" shown on the sample card and on the diagnostic
            # report header. Format in the source xlsx is
            # "YYYY-MM-DD HH:MM" (text or datetime cell — _strip handles
            # both).
            "sign_received_at": get(r, "簽收時間"),
        })
    return out


def list_options() -> dict[str, list[str]]:
    """Distinct department / physician values seen across the roster.

    Feeds the editable datalists on the sample card so reviewers can
    pick a known value or type a new one. Sorted for stable UI.
    """
    roster = load_roster()
    deps:  set[str] = set()
    docs:  set[str] = set()
    for entry in roster.values():
        d = _strip(entry.get("department"))
        if d: deps.add(d)
        p = _strip(entry.get("physician"))
        if p: docs.add(p)
    return {
        "departments": sorted(deps),
        "physicians":  sorted(docs),
    }


def load_roster() -> dict[str, dict]:
    """Return the merged roster ({lis_id: {...}}). Empty when none yet."""
    p = _roster_path()
    if not p.is_file():
        return {}
    try:
        data = json.loads(p.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, OSError):
        return {}
    return data if isinstance(data, dict) else {}


def lookup(lis_id: str) -> dict | None:
    """Roster entry for one LIS_ID, or None."""
    roster = load_roster()
    for key in lookup_candidates(lis_id):
        hit = roster.get(key)
        if hit:
            return hit
    return None


def lookup_with_key(lis_id: str, *aliases: str, roster: dict[str, dict] | None = None) -> tuple[dict | None, str]:
    """Return the first roster entry and key matching a sample ID.

    `aliases` is intended for original source sample IDs from
    pipeline_source.json. Returns ({...}, key) or (None, "").
    """
    roster = roster if roster is not None else load_roster()
    for key in lookup_candidates(lis_id, *aliases):
        hit = roster.get(key)
        if hit:
            return hit, key
    return None, ""


def ingest_xlsx(content: bytes, original_filename: str) -> dict:
    """Archive the upload + merge its rows into roster.json.

    Returns {added, updated, total, parsed, archive_name}.
    Raises ValueError when the xlsx has no recognisable data rows.
    """
    PATIENT_LIST_DIR.mkdir(parents=True, exist_ok=True)

    # 1. Archive the raw upload (timestamped so re-uploads don't clash).
    ts = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    safe_name = re.sub(r"[^A-Za-z0-9._-]+", "_", original_filename or "upload.xlsx")
    if not safe_name.lower().endswith((".xlsx", ".xlsm")):
        safe_name += ".xlsx"
    archive = PATIENT_LIST_DIR / f"{ts}_{safe_name}"
    archive.write_bytes(content)

    # 2. Parse it.
    try:
        parsed = parse_xlsx(archive)
    except Exception as exc:  # noqa: BLE001 — surface a clean message
        raise ValueError(f"無法解析 xlsx：{exc}") from exc
    if not parsed:
        raise ValueError("xlsx 裡找不到可辨識的資料列（缺『檢體編號』標題列？）")

    # 3. Merge into roster.json (additive — never drop existing keys).
    roster = load_roster()
    added = updated = 0
    now = _now_iso()
    for rec in parsed:
        lid = rec["lis_id"]
        prev = roster.get(lid)
        entry = {
            "lis_id":           lid,
            "specimen":         rec["specimen"],
            "mrn":              rec["mrn"],
            "name":             rec["name"],
            "test_name":        rec["test_name"],
            "test_type":        rec["test_type"],
            "department":       rec["department"],
            "physician":        rec["physician"],
            "sign_received_at": rec["sign_received_at"],
            "updated_at":       now,
        }
        if prev is None:
            entry["created_at"] = now
            roster[lid] = entry
            added += 1
        else:
            # Keep the original created_at; bump updated_at; overwrite
            # the data fields with the freshest values.
            entry["created_at"] = prev.get("created_at") or now
            # Only count as "updated" if something actually changed.
            changed = any(
                prev.get(k) != entry.get(k)
                for k in ("specimen", "mrn", "name", "test_name", "test_type",
                          "department", "physician", "sign_received_at")
            )
            roster[lid] = entry
            if changed:
                updated += 1

    _roster_path().write_text(
        json.dumps(roster, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )

    # Append to the upload-history log so reviewers can audit which xlsx
    # was ingested when and what each batch contributed.
    upload_record = {
        "uploaded_at":       now,
        "original_filename": original_filename or "",
        "archive_name":      archive.name,
        "parsed":            len(parsed),
        "added":             added,
        "updated":           updated,
        "total_after":       len(roster),
    }
    _append_upload(upload_record)

    return {
        "added":        added,
        "updated":      updated,
        "total":        len(roster),
        "parsed":       len(parsed),
        "archive_name": archive.name,
    }
