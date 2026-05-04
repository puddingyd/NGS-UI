"""SQLite-backed user store with bcrypt hashing.

Single-table schema: id / username / password_hash / created_at.
DB lives at `data/users.db`. Bootstrap accounts are created via
`python -m app.cli create-user`.
"""
from __future__ import annotations

import os
import secrets
import sqlite3
from datetime import datetime, timezone
from pathlib import Path

import bcrypt

from ..config import DATA_ROOT


_DB_PATH = DATA_ROOT / "users.db"
_BCRYPT_MAX_BYTES = 72  # bcrypt's hard cap; we surface a clear error instead


def _hash(password: str) -> str:
    if len(password.encode("utf-8")) > _BCRYPT_MAX_BYTES:
        raise ValueError(f"password must be ≤ {_BCRYPT_MAX_BYTES} bytes (bcrypt limit)")
    return bcrypt.hashpw(password.encode("utf-8"), bcrypt.gensalt()).decode("utf-8")


def _verify(password: str, hashed: str) -> bool:
    try:
        return bcrypt.checkpw(password.encode("utf-8"), hashed.encode("utf-8"))
    except (ValueError, TypeError):
        return False


def _conn() -> sqlite3.Connection:
    DATA_ROOT.mkdir(parents=True, exist_ok=True)
    c = sqlite3.connect(_DB_PATH)
    c.row_factory = sqlite3.Row
    c.execute("""
        CREATE TABLE IF NOT EXISTS users (
            id            INTEGER PRIMARY KEY,
            username      TEXT UNIQUE NOT NULL,
            password_hash TEXT NOT NULL,
            created_at    TEXT NOT NULL
        )
    """)
    return c


def create_user(username: str, password: str) -> dict:
    if not username or not password:
        raise ValueError("username + password required")
    if len(password) < 8:
        raise ValueError("password must be at least 8 characters")
    h = _hash(password)
    now = datetime.now(timezone.utc).isoformat(timespec="seconds")
    with _conn() as c:
        c.execute(
            "INSERT INTO users(username, password_hash, created_at) VALUES (?,?,?)",
            (username, h, now),
        )
        row = c.execute(
            "SELECT id, username, created_at FROM users WHERE username=?",
            (username,),
        ).fetchone()
    return dict(row)


def verify_user(username: str, password: str) -> dict | None:
    with _conn() as c:
        row = c.execute(
            "SELECT id, username, password_hash, created_at FROM users WHERE username=?",
            (username,),
        ).fetchone()
    if not row or not _verify(password, row["password_hash"]):
        return None
    return {"id": row["id"], "username": row["username"], "created_at": row["created_at"]}


def get_by_id(user_id: int) -> dict | None:
    with _conn() as c:
        row = c.execute(
            "SELECT id, username, created_at FROM users WHERE id=?",
            (user_id,),
        ).fetchone()
    return dict(row) if row else None


def list_users() -> list[dict]:
    with _conn() as c:
        rows = c.execute(
            "SELECT id, username, created_at FROM users ORDER BY id"
        ).fetchall()
    return [dict(r) for r in rows]


def session_secret() -> str:
    """Read the cookie-signing secret. Lazily generated on first run."""
    env = os.environ.get("SESSION_SECRET")
    if env:
        return env
    p = DATA_ROOT / "session_secret"
    if p.exists():
        return p.read_text(encoding="utf-8").strip()
    DATA_ROOT.mkdir(parents=True, exist_ok=True)
    s = secrets.token_urlsafe(48)
    p.write_text(s, encoding="utf-8")
    p.chmod(0o600)
    return s
