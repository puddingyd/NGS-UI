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
from ..services import capture_tokens, omim_capture, omim_db, omim_store

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
    return {"capture": stored, "summary": summary}


def _index_blockers(event: dict) -> list[str]:
    """Bulk-adopt blockers estimated from an index.jsonl event."""
    return omim_db.adoption_blockers({
        "mim": event.get("mim"),
        "entry_type": omim_capture.ENTRY_TYPES.get(event.get("prefix") or "", "other"),
        "warnings": event.get("warnings") or [],
        "synopsis_status": event.get("synopsis_status") or "",
        "description": event.get("description_source") or "",
        "relationships": [None] * int(event.get("relationship_count") or 0),
    })


@router.get("/captures")
def list_captures(limit: int = 100, view: str = "pending", user: dict = Depends(current_user)):
    """``view=pending``: newest capture per MIM not yet adopted; ``all``: raw log."""
    limit = max(1, min(int(limit), 1000))
    if view == "all":
        return {"captures": omim_capture.recent_captures(limit)}
    adopted = omim_db.adopted_capture_ids()
    data = omim_capture.pending_captures(adopted, limit)
    for event in data["captures"]:
        event["blockers"] = _index_blockers(event) if omim_db.exists() else []
    data["db_enabled"] = omim_db.exists()
    return data


@router.get("/adoptions")
def list_adoptions(limit: int = 200, user: dict = Depends(current_user)):
    return {"adoptions": omim_db.list_adoptions(max(1, min(int(limit), 2000)))}


@router.post("/captures/adopt-pending")
def adopt_pending(payload: dict | None = None, user: dict = Depends(current_user)):
    """Adopt many pending captures; ones with blockers are skipped, not forced."""
    if not omim_db.exists():
        raise HTTPException(409, "平台 OMIM 資料庫尚未建立（請先執行 scripts/omim_db.py init）")
    ids = [str(i) for i in (payload or {}).get("ids") or []]
    pending = {e["id"] for e in omim_capture.pending_captures(omim_db.adopted_capture_ids(), 100000)["captures"]}
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
        reasons = omim_db.adoption_blockers(parsed)
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
def resolve_review(item_id: int, payload: dict, user: dict = Depends(current_user)):
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
