"""Auth dependencies.

Uses Starlette SessionMiddleware (signed cookie) so we don't need a
session table. Session contents: {"user_id": int, "username": str}.
"""
from __future__ import annotations

from fastapi import HTTPException, Request

from .services import users


def current_user(request: Request) -> dict:
    sess = getattr(request, "session", None) or {}
    uid = sess.get("user_id")
    if not uid:
        raise HTTPException(401, "not authenticated")
    user = users.get_by_id(uid)
    if not user:
        # The session points at a deleted user; force re-login.
        sess.clear()
        raise HTTPException(401, "user no longer exists")
    return user
