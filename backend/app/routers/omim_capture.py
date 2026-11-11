"""OMIM capture endpoints for the browser extension (proof of concept).

Extension calls authenticate with ``Authorization: Bearer ngsc_…`` (a
personal capture token); the platform's own pages may use the normal
session cookie. Captures are stored and compared only — nothing is
written to OMIM.xlsx yet.
"""
from __future__ import annotations

import io
import json
import re
import zipfile

from fastapi import APIRouter, Depends, HTTPException, Request
from fastapi.responses import Response

from ..auth import current_user
from ..config import REPO_ROOT
from ..services import capture_tokens, omim_capture, omim_db, omim_lookup, omim_queue, omim_store, omim_updates

router = APIRouter(prefix="/api/omim", tags=["omim-capture"])

EXTENSION_DIR = REPO_ROOT / "reference" / "omim_capture"
MAX_BODY_BYTES = 8 * 1024 * 1024
_MIM_RE = re.compile(r"^\d{6}$")


def capture_user(request: Request) -> dict:
    """Token (extension) or session (platform page) authentication."""
    header = request.headers.get("authorization") or ""
    if header.lower().startswith("bearer "):
        user = capture_tokens.verify(header[7:])
        if not user:
            raise HTTPException(401, "invalid capture token")
        return user
    return current_user(request)


@router.get("/ping")
def ping(user: dict = Depends(capture_user)):
    return {"ok": True, "username": user.get("username")}


@router.get("/status")
def status(mim: str, edit_date: str = "", user: dict = Depends(capture_user)):
    """Platform view of one MIM; ``edit_date`` is the page's OMIM edit date."""
    if not _MIM_RE.match(mim or ""):
        raise HTTPException(400, "mim must be 6 digits")
    if edit_date and not re.match(r"^\d{4}-\d{2}-\d{2}$", edit_date):
        raise HTTPException(400, "edit_date must be YYYY-MM-DD")
    latest = next((c for c in omim_capture.recent_captures(500) if c.get("mim") == mim), None)
    comparison = omim_capture.compare_with_workbook(
        mim,
        entry_type="",
        edit_date=edit_date or (latest or {}).get("omim_edit_date") or "",
    )
    return {"mim": mim, "comparison": comparison, "latest_capture": latest}


@router.post("/captures")
async def create_capture(request: Request, user: dict = Depends(capture_user)):
    body = await request.body()
    if len(body) > MAX_BODY_BYTES:
        raise HTTPException(413, "capture too large")
    try:
        payload = json.loads(body.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError):
        raise HTTPException(400, "body must be JSON")
    if not isinstance(payload, dict):
        raise HTTPException(400, "body must be a JSON object")
    text = str(payload.get("text") or "")
    if not text.strip():
        raise HTTPException(400, "empty page text")
    if len(text) > omim_capture.MAX_TEXT_CHARS:
        raise HTTPException(413, "page text too large")
    url_mim = str(payload.get("mim") or "")
    if url_mim and not _MIM_RE.match(url_mim):
        raise HTTPException(400, "mim must be 6 digits")
    extracted = payload.get("extracted") if isinstance(payload.get("extracted"), dict) else {}
    parsed = omim_capture.parse_entry_text(text, url_mim=url_mim, extracted=extracted)
    if not parsed.get("mim"):
        raise HTTPException(422, "could not identify the OMIM entry on this page")
    stored = omim_capture.store_capture(payload, parsed, user=user)
    comparison = omim_capture.compare_with_workbook(
        parsed["mim"],
        entry_type=parsed.get("entry_type") or "",
        edit_date=parsed.get("omim_edit_date") or "",
        relationships=parsed.get("relationships") or [],
    )
    summary = omim_capture.summarize(parsed, comparison)
    summary["db_enabled"] = omim_db.exists()
    summary["adopt_preview"] = omim_db.preview_adoption(parsed) if summary["db_enabled"] else None
    summary["blockers"] = omim_db.adoption_blockers(parsed)
    return {"capture": stored, "summary": summary}


def _stale_blocker(mim: str, captured_at: str, latest_any: dict[str, str]) -> list[str]:
    listed = omim_updates.is_stale(captured_at, mim, latest_any)
    return [f"OMIM 之後又更新（{listed}），請重新擷取"] if listed else []


def _index_blockers(event: dict, latest_any: dict[str, str] | None = None) -> list[str]:
    """Bulk-adopt blockers estimated from an index.jsonl event."""
    stale = _stale_blocker(event.get("mim") or "", event.get("captured_at") or "", latest_any or {})
    return stale + omim_db.adoption_blockers({
        "mim": event.get("mim"),
        "entry_type": omim_capture.ENTRY_TYPES.get(event.get("prefix") or "", "other"),
        "warnings": event.get("warnings") or [],
        "synopsis_status": event.get("synopsis_status") or "",
        "description": event.get("description_source") or "",
        "relationships": [None] * int(event.get("relationship_count") or 0),
        "relationship_section": event.get("relationship_section") or "",
    })


@router.get("/captures")
def list_captures(limit: int = 100, view: str = "pending", user: dict = Depends(current_user)):
    """``view=pending``: newest capture per MIM not yet adopted; ``all``: raw log."""
    limit = max(1, min(int(limit), 1000))
    if view == "all":
        return {"captures": omim_capture.recent_captures(limit)}
    adopted = omim_db.adopted_capture_ids()
    data = omim_capture.pending_captures(adopted, limit)
    latest_any = omim_updates.latest_dates()[1]
    for event in data["captures"]:
        if event.get("parser_version") != omim_capture.PARSER_VERSION:
            # Captured with older parser rules: re-read so the list shows
            # what adopting would actually use (e.g. old-style synopses).
            mim, _, name = str(event.get("id") or "").partition("/")
            try:
                _data, parsed = _load_parsed(mim, name)
            except HTTPException:
                parsed = None
            if parsed:
                event.update(
                    synopsis_status=parsed.get("synopsis_status"),
                    description_source=parsed.get("description_source"),
                    warnings=parsed.get("warnings") or [],
                    relationship_count=len(parsed.get("relationships") or []),
                    relationship_section=parsed.get("relationship_section") or "",
                    title=parsed.get("title") or event.get("title"),
                )
        event["blockers"] = _index_blockers(event, latest_any) if omim_db.exists() else []
    data["db_enabled"] = omim_db.exists()
    return data


_QUEUE_STATUSES = {"", "in_progress", "updated", "gene_sync", "phenotype", "review"}
_QUEUE_SCOPES = {"all", "panel", "reported"}
_SKIP_RE = re.compile(r"^g?(new)?\d{6}$")


def _pending_mims() -> set[str]:
    """MIMs captured but not adopted yet — unless OMIM updated them after
    the capture, so the queue sends the reviewer back for the new version."""
    pending = omim_capture.pending_captures(omim_db.adopted_capture_ids(), 100000)["captures"]
    latest_any = omim_updates.latest_dates()[1]
    return {e.get("mim") for e in pending if e.get("mim")
            and not omim_updates.is_stale(e.get("captured_at") or "", e.get("mim"), latest_any)}


@router.get("/queue")
def work_queue(status: str = "", scope: str = "all", offset: int = 0, limit: int = 100,
               user: dict = Depends(current_user)):
    """Genes to work through, with their gene-page / disease-page steps."""
    if status not in _QUEUE_STATUSES or scope not in _QUEUE_SCOPES:
        raise HTTPException(400, "invalid status or scope")
    return omim_queue.queue(status=status, scope=scope, pending_mims=_pending_mims(),
                            recent_mims=omim_db.recently_adopted_mims(),
                            offset=max(0, int(offset)), limit=max(1, min(int(limit), 1000)))


@router.get("/queue/next")
def work_queue_next(current: str = "", gene: str = "", skip: str = "", scope: str = "all",
                    user: dict = Depends(capture_user)):
    """Next page for the extension; stays on ``gene`` until its steps are done.

    ``skip``: comma-separated MIMs, or ``g<geneMIM>`` to skip a whole gene.
    """
    if scope not in _QUEUE_SCOPES:
        raise HTTPException(400, "invalid scope")
    skipped = {m for m in skip.split(",") if _SKIP_RE.match(m)}
    return omim_queue.next_step(
        current=current if _MIM_RE.match(current or "") else "",
        gene=gene if _MIM_RE.match(gene or "") else "",
        skip=skipped, scope=scope, pending_mims=_pending_mims(),
        recent_mims=omim_db.recently_adopted_mims(),
    )


@router.post("/updates/month")
async def updates_month(request: Request, user: dict = Depends(capture_user)):
    """One month of the OMIM Update List, read by the extension."""
    if not omim_db.exists():
        raise HTTPException(409, "平台 OMIM 資料庫尚未建立（請先執行 scripts/omim_db.py init）")
    body = await request.body()
    if len(body) > MAX_BODY_BYTES:
        raise HTTPException(413, "page too large")
    try:
        payload = json.loads(body.decode("utf-8"))
        year, month = int(payload.get("year")), int(payload.get("month"))
    except (UnicodeDecodeError, json.JSONDecodeError, AttributeError, TypeError, ValueError):
        raise HTTPException(400, "body must be JSON with year and month")
    rows = payload.get("rows") if isinstance(payload.get("rows"), list) else []
    text = str(payload.get("text") or "")[:omim_capture.MAX_TEXT_CHARS]
    try:
        result = omim_updates.ingest(year, month, rows, text, user=user.get("username") or "")
    except ValueError as exc:
        raise HTTPException(422, str(exc))
    return {**result, "coverage": _coverage_summary()}


def _coverage_summary() -> dict:
    cov = omim_updates.coverage()
    return {k: cov[k] for k in ("needed_from", "current_month", "missing", "coverage_start",
                                "history_complete", "rows", "last_listed_date", "last_captured_at")}


@router.get("/updates/status")
def updates_status(user: dict = Depends(capture_user)):
    """Update List coverage, plus what it means for the work queue."""
    if not omim_db.exists():
        return {"db_enabled": False}
    cov = omim_updates.coverage()
    counts = omim_queue.queue(limit=1)["counts"]
    diseases = omim_updates.new_diseases()
    return {
        "db_enabled": True,
        **cov,
        "months": cov["months"][:240],
        "stale_mode": omim_updates.stale_mode(),
        "queue": {k: counts.get(k, 0) for k in ("updated", "updated_steps", "new_disease_steps", "stale_steps",
                                                 "stale_confirmed", "stale_hidden")},
        "new_diseases": diseases[:500],
        "new_diseases_total": len(diseases),
    }


@router.post("/settings")
def update_settings(payload: dict, user: dict = Depends(current_user)):
    mode = str((payload or {}).get("stale_mode") or "")
    if mode not in omim_updates.STALE_MODES:
        raise HTTPException(400, "stale_mode must be both or updates_only")
    if not omim_db.exists():
        raise HTTPException(409, "平台 OMIM 資料庫尚未建立")
    omim_updates.set_setting("stale_mode", mode)
    return {"stale_mode": mode}


@router.get("/genes/search")
def gene_search(q: str = "", limit: int = 20, user: dict = Depends(current_user)):
    """Gene symbol typeahead for the lookup box."""
    return {"items": omim_lookup.search(q[:40], max(1, min(int(limit), 50)))}


@router.get("/gene")
def gene_lookup(q: str = "", user: dict = Depends(current_user)):
    """A gene's OMIM data as the platform uses it.

    ``q``: gene symbol or alias, gene MIM, or phenotype MIM (may match
    several genes; then only a short list is returned to pick from).
    """
    q = q.strip()[:40]
    if not q:
        raise HTTPException(400, "請輸入基因名稱或 MIM")
    rows = omim_lookup.resolve(q)
    if not rows:
        return {"query": q, "matches": [], "gene": None,
                "suggestions": omim_lookup.search(q, 10)}
    if len(rows) > 1:
        return {"query": q, "matches": [omim_lookup.summary(r) for r in rows], "gene": None}
    return {"query": q, "matches": [omim_lookup.summary(rows[0])], "gene": omim_lookup.detail(rows[0])}


@router.get("/adoptions")
def list_adoptions(limit: int = 200, user: dict = Depends(current_user)):
    return {"adoptions": omim_db.list_adoptions(max(1, min(int(limit), 2000)))}


@router.post("/captures/adopt-pending")
def adopt_pending(payload: dict | None = None, user: dict = Depends(current_user)):
    """Adopt many pending captures; ones with blockers are skipped, not forced."""
    if not omim_db.exists():
        raise HTTPException(409, "平台 OMIM 資料庫尚未建立（請先執行 scripts/omim_db.py init）")
    ids = [str(i) for i in (payload or {}).get("ids") or []]
    pending_events = omim_capture.pending_captures(omim_db.adopted_capture_ids(), 100000)["captures"]
    pending = {e["id"]: e for e in pending_events}
    latest_any = omim_updates.latest_dates()[1]
    adopted, skipped = [], []
    for capture_id in ids[:500]:
        mim, _, name = capture_id.partition("/")
        if capture_id not in pending:
            skipped.append({"id": capture_id, "reasons": ["已採用或已有較新的擷取"]})
            continue
        try:
            _data, parsed = _load_parsed(mim, name)
        except HTTPException:
            skipped.append({"id": capture_id, "reasons": ["找不到擷取檔"]})
            continue
        reasons = _stale_blocker(mim, pending[capture_id].get("captured_at") or "", latest_any) \
            + omim_db.adoption_blockers(parsed)
        if reasons:
            skipped.append({"id": capture_id, "mim": mim, "title": parsed.get("title"), "reasons": reasons})
            continue
        try:
            result = omim_db.adopt_capture(parsed, user=user.get("username") or "", capture_id=capture_id)
        except ValueError as exc:
            skipped.append({"id": capture_id, "mim": mim, "reasons": [str(exc)]})
            continue
        adopted.append({"id": capture_id, "mim": mim, "title": parsed.get("title"),
                        "genes": [g.get("gene") for g in result.get("genes") or []],
                        "reviews": result.get("reviews") or 0})
    omim_store.ensure_loaded()
    return {"adopted": adopted, "skipped": skipped}


def _load_parsed(mim: str, name: str) -> tuple[dict, dict]:
    data = omim_capture.load_capture(f"{mim}/{name}")
    if not data:
        raise HTTPException(404, "capture not found")
    parsed = data.get("parsed") or {}
    payload = data.get("payload") or {}
    if parsed.get("parser_version") != omim_capture.PARSER_VERSION and payload.get("text"):
        # Older captures are re-read with the current parser rules.
        parsed = omim_capture.parse_entry_text(
            payload["text"],
            url_mim=str(payload.get("mim") or mim),
            extracted=payload.get("extracted") if isinstance(payload.get("extracted"), dict) else {},
        )
    return data, parsed


def _summary(parsed: dict) -> dict:
    comparison = omim_capture.compare_with_workbook(
        parsed.get("mim") or "",
        entry_type=parsed.get("entry_type") or "",
        edit_date=parsed.get("omim_edit_date") or "",
        relationships=parsed.get("relationships") or [],
    )
    summary = omim_capture.summarize(parsed, comparison)
    summary["db_enabled"] = omim_db.exists()
    summary["adopt_preview"] = omim_db.preview_adoption(parsed) if summary["db_enabled"] else None
    summary["blockers"] = omim_db.adoption_blockers(parsed)
    return summary


@router.get("/captures/{mim}/{name}")
def get_capture(mim: str, name: str, user: dict = Depends(current_user)):
    data, parsed = _load_parsed(mim, name)
    return {"meta": data.get("meta"), "summary": _summary(parsed)}


@router.post("/captures/{mim}/{name}/adopt")
def adopt_capture(mim: str, name: str, user: dict = Depends(capture_user)):
    """Write a capture into the platform OMIM database."""
    if not omim_db.exists():
        raise HTTPException(409, "平台 OMIM 資料庫尚未建立（請先執行 scripts/omim_db.py init）")
    _data, parsed = _load_parsed(mim, name)
    try:
        result = omim_db.adopt_capture(
            parsed, user=user.get("username") or "", capture_id=f"{mim}/{name}",
        )
    except ValueError as exc:
        raise HTTPException(422, str(exc))
    omim_store.ensure_loaded()
    by_gene = omim_db.open_reviews_by_gene()
    result["gene_reviews"] = [
        item for g in result.get("genes") or [] for item in by_gene.get(g.get("gene_mim"), [])
    ]
    return {"result": result, "summary": _summary(parsed)}


@router.get("/db/status")
def db_status(user: dict = Depends(current_user)):
    return {**omim_db.status(), "source": omim_store.source()}


@router.get("/reviews")
def list_reviews(status: str = "open", kind: str = "", limit: int = 500,
                 user: dict = Depends(current_user)):
    if not omim_db.exists():
        return {"items": [], "open_counts": {}, "db_enabled": False}
    if status not in {"open", "resolved"}:
        raise HTTPException(400, "status must be open or resolved")
    return {**omim_db.list_reviews(status, kind, max(1, min(int(limit), 2000))), "db_enabled": True,
            "kinds": omim_db.REVIEW_KINDS}


@router.post("/reviews/{item_id}")
def resolve_review(item_id: int, payload: dict, user: dict = Depends(capture_user)):
    if not omim_db.exists():
        raise HTTPException(409, "平台 OMIM 資料庫尚未建立")
    try:
        result = omim_db.resolve_review(item_id, str((payload or {}).get("action") or ""),
                                        user=user.get("username") or "")
    except KeyError:
        raise HTTPException(404, "review item not found")
    except ValueError as exc:
        raise HTTPException(400, str(exc))
    omim_store.ensure_loaded()
    return result


@router.get("/export.xlsx")
def export_xlsx(user: dict = Depends(current_user)):
    if not omim_db.exists():
        raise HTTPException(404, "平台 OMIM 資料庫尚未建立")
    from datetime import date as _date
    return Response(
        omim_db.export_xlsx(),
        media_type="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
        headers={"Content-Disposition": f'attachment; filename="OMIM_{_date.today():%Y%m%d}.xlsx"'},
    )


@router.get("/capture-tokens")
def list_tokens(user: dict = Depends(current_user)):
    return {"tokens": capture_tokens.list_for_user(user)}


@router.post("/capture-tokens")
def create_token(payload: dict | None = None, user: dict = Depends(current_user)):
    return capture_tokens.create(user, (payload or {}).get("label") or "")


@router.delete("/capture-tokens/{token_id}")
def delete_token(token_id: int, user: dict = Depends(current_user)):
    if not capture_tokens.revoke(user, token_id):
        raise HTTPException(404, "token not found")
    return {"ok": True}


@router.get("/extension.zip")
def extension_zip(user: dict = Depends(current_user)):
    if not EXTENSION_DIR.is_dir():
        raise HTTPException(404, "extension sources not found")
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
        for path in sorted(EXTENSION_DIR.rglob("*")):
            if path.is_file() and not path.name.startswith("."):
                zf.write(path, f"omim_capture/{path.relative_to(EXTENSION_DIR)}")
    return Response(
        buf.getvalue(),
        media_type="application/zip",
        headers={"Content-Disposition": 'attachment; filename="ngs-omim-capture.zip"'},
    )
