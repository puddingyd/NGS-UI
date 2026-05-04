from __future__ import annotations

from fastapi import APIRouter, Depends, HTTPException, Request

from ..auth import current_user
from ..services import users

router = APIRouter(prefix="/api/auth", tags=["auth"])


@router.post("/login")
def login(request: Request, payload: dict):
    username = (payload or {}).get("username", "").strip()
    password = (payload or {}).get("password", "")
    user = users.verify_user(username, password)
    if not user:
        raise HTTPException(401, "incorrect username or password")
    request.session["user_id"]  = user["id"]
    request.session["username"] = user["username"]
    return {"username": user["username"]}


@router.post("/logout")
def logout(request: Request):
    request.session.clear()
    return {"ok": True}


@router.get("/me")
def me(user: dict = Depends(current_user)):
    return {"username": user["username"], "id": user["id"]}
