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
from ..services import capture_tokens, omim_capture

router = APIRouter(prefix="/api/omim", tags=["omim-capture"])

EXTENSION_DIR = REPO_ROOT / "browser_extension" / "omim_capture"
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
    return {"capture": stored, "summary": omim_capture.summarize(parsed, comparison)}


@router.get("/captures")
def list_captures(limit: int = 50, user: dict = Depends(current_user)):
    return {"captures": omim_capture.recent_captures(max(1, min(int(limit), 500)))}


@router.get("/captures/{mim}/{name}")
def get_capture(mim: str, name: str, user: dict = Depends(current_user)):
    data = omim_capture.load_capture(f"{mim}/{name}")
    if not data:
        raise HTTPException(404, "capture not found")
    parsed = data.get("parsed") or {}
    comparison = omim_capture.compare_with_workbook(
        parsed.get("mim") or mim,
        entry_type=parsed.get("entry_type") or "",
        edit_date=parsed.get("omim_edit_date") or "",
        relationships=parsed.get("relationships") or [],
    )
    return {"meta": data.get("meta"), "summary": omim_capture.summarize(parsed, comparison)}


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
