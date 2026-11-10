"""Personal API tokens for the OMIM capture browser extension.

The extension posts from the omim.org tab's background worker, where the
platform's SameSite session cookie is not available, so each reviewer
creates a token on the platform and pastes it into the extension options.
Only a SHA-256 hash is stored (in ``users.db``); the plain token is shown
once at creation.
"""
from __future__ import annotations

import hashlib
import secrets
from datetime import datetime, timezone

from . import users

TOKEN_PREFIX = "ngsc_"


def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def _hash(token: str) -> str:
    return hashlib.sha256(token.encode("utf-8")).hexdigest()


def _conn():
    c = users._conn()
    c.execute("""
        CREATE TABLE IF NOT EXISTS capture_tokens (
            id           INTEGER PRIMARY KEY,
            token_hash   TEXT UNIQUE NOT NULL,
            user_id      INTEGER NOT NULL,
            label        TEXT NOT NULL DEFAULT '',
            created_at   TEXT NOT NULL,
            last_used_at TEXT
        )
    """)
    return c


def create(user: dict, label: str = "") -> dict:
    token = TOKEN_PREFIX + secrets.token_urlsafe(32)
    label = str(label or "").strip()[:60]
    with _conn() as c:
        cur = c.execute(
            "INSERT INTO capture_tokens(token_hash, user_id, label, created_at) VALUES (?,?,?,?)",
            (_hash(token), int(user["id"]), label, _now()),
        )
        token_id = cur.lastrowid
    return {"id": token_id, "token": token, "label": label}


def list_for_user(user: dict) -> list[dict]:
    with _conn() as c:
        rows = c.execute(
            "SELECT id, label, created_at, last_used_at FROM capture_tokens "
            "WHERE user_id = ? ORDER BY id DESC",
            (int(user["id"]),),
        ).fetchall()
    return [dict(r) for r in rows]


def revoke(user: dict, token_id: int) -> bool:
    with _conn() as c:
        cur = c.execute(
            "DELETE FROM capture_tokens WHERE id = ? AND user_id = ?",
            (int(token_id), int(user["id"])),
        )
    return cur.rowcount > 0


def verify(token: str) -> dict | None:
    token = str(token or "").strip()
    if not token.startswith(TOKEN_PREFIX):
        return None
    with _conn() as c:
        row = c.execute(
            "SELECT id, user_id FROM capture_tokens WHERE token_hash = ?",
            (_hash(token),),
        ).fetchone()
        if not row:
            return None
        c.execute(
            "UPDATE capture_tokens SET last_used_at = ? WHERE id = ?",
            (_now(), row["id"]),
        )
    return users.get_by_id(row["user_id"])
