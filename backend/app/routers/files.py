"""Authenticated file inventory, bounded tertiary browsing and downloads."""
from __future__ import annotations

import asyncio
from html import escape
import os
from pathlib import Path, PurePosixPath
import re
import sqlite3
import time
from urllib.parse import quote, urlencode

from fastapi import APIRouter, BackgroundTasks, Depends, HTTPException, Query, Request
from fastapi.responses import HTMLResponse, StreamingResponse
from starlette.background import BackgroundTask

from ..auth import current_user
from ..services import file_index

router = APIRouter(prefix="/api/files", tags=["files"], dependencies=[Depends(current_user)])


async def ensure_index(background):
    meta = file_index.metadata()
    if meta is None:
        try:
            await asyncio.to_thread(file_index.refresh)
        except (OSError, ValueError, sqlite3.Error) as exc:
            raise HTTPException(503, str(exc))
    elif file_index.stale(meta):
        background.add_task(file_index.refresh)


@router.get("")
async def inventory(background: BackgroundTasks, q: str = Query("", max_length=128),
                    offset: int = Query(0, ge=0), limit: int = Query(50, ge=1, le=100)):
    await ensure_index(background)
    return await asyncio.to_thread(file_index.list_samples, q.strip(), offset, limit)


@router.post("/refresh")
async def refresh():
    requested = time.time_ns()
    try:
        views = await asyncio.to_thread(file_index.refresh, force=True, requested_at_ns=requested)
        return views["meta"]
    except (OSError, ValueError, sqlite3.Error) as exc:
        raise HTTPException(503, str(exc))


def _open(path, roots, directory=False):
    try:
        return file_index.open_bounded(path, roots, directory=directory)
    except FileNotFoundError:
        raise HTTPException(404, "找不到檔案或資料夾")
    except (OSError, ValueError):
        raise HTTPException(403, "無法讀取此路徑")


def _download(request, path, roots):
    fd, st = _open(path, roots)
    stream = os.fdopen(fd, "rb")
    size, start, end = st.st_size, 0, st.st_size-1
    headers = {"content-disposition": "attachment; filename*=UTF-8''" + quote(path.name, safe=""),
               "accept-ranges": "bytes", "cache-control": "private, no-store",
               "x-content-type-options": "nosniff", "content-encoding": "identity"}
    status = 200
    rng = request.headers.get("range")
    if rng:
        match = re.fullmatch(r"bytes=(\d*)-(\d*)", rng.strip()) if len(rng) <= 128 else None
        if not match or not any(match.groups()) or not size:
            stream.close()
            raise HTTPException(416, "無效的下載範圍", headers={"Content-Range": f"bytes */{size}"})
        first, last = match.groups()
        if first:
            start, end = int(first), min(int(last), end) if last else end
        else:
            start = max(size-int(last), 0)
        if start > end or start >= size:
            stream.close()
            raise HTTPException(416, "無效的下載範圍", headers={"Content-Range": f"bytes */{size}"})
        headers["content-range"] = f"bytes {start}-{end}/{size}"
        status = 206
    length = max(0, end-start+1)
    headers["content-length"] = str(length)

    def chunks():
        try:
            if request.method == "HEAD":
                return
            stream.seek(start)
            remaining = length
            while remaining:
                block = stream.read(min(1024*1024, remaining))
                if not block:
                    break
                remaining -= len(block)
                yield block
        finally:
            stream.close()
    return StreamingResponse(chunks(), status_code=status, media_type="application/octet-stream",
                             headers=headers, background=BackgroundTask(stream.close))


@router.api_route("/download/{file_id}", methods=["GET", "HEAD"])
def download(request: Request, file_id: str):
    if not re.fullmatch(r"[0-9a-f]{64}", file_id):
        raise HTTPException(404, "找不到檔案")
    try:
        path = file_index.download_path(file_id)
    except (FileNotFoundError, sqlite3.Error):
        raise HTTPException(404, "找不到檔案")
    return _download(request, path, [p for roots in file_index.data_roots().values() for p in roots])


def _tertiary_path(root_id, relative):
    root = next((r for r in file_index.tertiary_roots() if file_index.identity(r) == root_id), None)
    parts = PurePosixPath(relative).parts
    if root is None or relative.startswith("/") or "\\" in relative or "\x00" in relative or any(p.startswith(".") for p in parts):
        raise HTTPException(403, "路徑超出三級資料範圍")
    return root, root.joinpath(*parts)


@router.api_route("/tertiary-download", methods=["GET", "HEAD"])
def tertiary_download(request: Request, root: str, path: str):
    base, target = _tertiary_path(root, path)
    return _download(request, target, [base])


@router.get("/browse", response_class=HTMLResponse)
def browse(root: str, path: str = ""):
    base, target = _tertiary_path(root, path)
    fd, _ = _open(target, [base], directory=True)
    entries = []
    try:
        with os.scandir(fd) as directory:
            for entry in directory:
                if entry.name.startswith("."):
                    continue
                candidate = target / entry.name
                try:
                    resolved = candidate.resolve(strict=True)
                    if not resolved.is_relative_to(base):
                        continue
                    is_dir = entry.is_dir()
                    if not is_dir and not entry.is_file():
                        continue
                    st = entry.stat()
                    entries.append((not is_dir, entry.name, st.st_size, st.st_mtime))
                except OSError:
                    continue
    finally:
        os.close(fd)

    def url(endpoint, relative):
        return "/api/files/" + endpoint + "?" + urlencode(dict(root=root, path=relative))

    crumbs = [f'<a href="{escape(url("browse", ""), quote=True)}">三級資料</a>']
    parts = PurePosixPath(path).parts
    for i, part in enumerate(parts):
        crumbs.append(f'<a href="{escape(url("browse", "/".join(parts[:i+1])), quote=True)}">{escape(part)}</a>')
    rows = []
    if parts:
        parent = "/".join(parts[:-1])
        rows.append(f'<tr><td><a href="{escape(url("browse", parent), quote=True)}">↑ 上一層</a></td><td></td><td></td></tr>')
    for is_file, name, size, mtime in sorted(entries, key=lambda e: (e[0], e[1].lower())):
        relative = "/".join((*parts, name))
        href = url("tertiary-download" if is_file else "browse", relative)
        stamp = time.strftime("%Y-%m-%d %H:%M", time.gmtime(mtime + 8*3600))
        rows.append(f'<tr><td><a href="{escape(href, quote=True)}">{escape(name)}{"" if is_file else "/"}</a></td>'
                    f'<td>{format(size, ",") if is_file else ""}</td><td>{stamp}</td></tr>')
    title = escape(target.name if parts else "三級資料")
    return HTMLResponse(f'''<!doctype html><html lang="zh-Hant"><head><meta charset="utf-8">
      <meta name="viewport" content="width=device-width,initial-scale=1"><title>{title} · 三級資料</title>
      <link rel="stylesheet" href="/style.css"><link rel="stylesheet" href="/files.css"></head>
      <body class="file-browser"><main><h1>三級資料</h1><nav aria-label="資料夾路徑">{' / '.join(crumbs)}</nav>
      <div class="file-table-scroll"><table class="file-browser-table"><thead><tr><th>檔名</th><th>大小（bytes）</th><th>修改時間（台北）</th></tr></thead>
      <tbody>{''.join(rows)}</tbody></table></div></main></body></html>''',
      headers={"Cache-Control": "private, no-store", "X-Content-Type-Options": "nosniff",
               "Content-Security-Policy": "default-src 'none'; style-src 'self'; base-uri 'none'; frame-ancestors 'self'"})
