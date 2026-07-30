"""Structured manual ACMG editing and observed-case lookup."""
from __future__ import annotations

import json

from fastapi import APIRouter, Depends, HTTPException

from ..auth import current_user
from ..services import manual_acmg, report_store, sample_layout, sample_loader


router = APIRouter(prefix="/api", tags=["acmg"], dependencies=[Depends(current_user)])


def _sample_meta(sample_id: str) -> dict:
    path = sample_layout.state_file(sample_id, "sample_metadata.json")
    if not path.is_file():
        raise HTTPException(404, f"sample not found: {sample_id}")
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        payload = {}
    return payload if isinstance(payload, dict) else {}


@router.get("/acmg/catalog")
def get_acmg_catalog():
    return manual_acmg.catalog()


@router.get("/samples/{sample_id}/variants/{variant_id}/acmg")
def get_variant_acmg(sample_id: str, variant_id: str):
    meta = _sample_meta(sample_id)
    try:
        normalized = manual_acmg.normalize_variant_id(variant_id)
        current = manual_acmg.current_assertion(
            meta.get("genome_build") or "hg38", normalized
        )
    except ValueError as exc:
        raise HTTPException(400, str(exc)) from exc
    report = report_store.load(sample_id)
    edit = (report.get("edits") or {}).get(variant_id) or {}
    return {
        "genome_build": manual_acmg.normalize_build(meta.get("genome_build") or "hg38"),
        "variant_id": normalized,
        "sample_snapshot": edit.get("manual_acmg"),
        "manual_current": current,
    }


@router.put("/samples/{sample_id}/variants/{variant_id}/acmg")
def put_variant_acmg(
    sample_id: str,
    variant_id: str,
    payload: dict,
    user: dict = Depends(current_user),
):
    meta = _sample_meta(sample_id)
    try:
        assertion = manual_acmg.save_assertion(
            meta.get("genome_build") or "hg38",
            variant_id,
            (payload or {}).get("criteria") or {},
            reviewer_user_id=user.get("id"),
            reviewer_username=user.get("username") or "",
            source_sample_id=sample_id,
        )
    except ValueError as exc:
        raise HTTPException(400, str(exc)) from exc
    report_store.save_manual_acmg(
        sample_id, variant_id, assertion, user=user
    )
    loaded = sample_loader.load_sample(sample_id, include_aux=False)
    variant = (loaded or {}).get("variants", {}).get(variant_id)
    if variant is None:
        # Coordinate normalization can add ``chr`` to old identifiers.
        variant = (loaded or {}).get("variants", {}).get(assertion["variant_id"])
    loaded_variants = (loaded or {}).get("variants") or {}
    categories = {
        key: value
        for key, value in ((loaded or {}).get("categories") or {}).items()
        if key in {"1A", "1B", "1C", "2"}
    }
    categories.update(
        sample_loader._build_secondary_snv_categories(loaded_variants)
    )
    return {
        "saved": assertion,
        "variant": variant,
        "categories": categories,
    }


@router.get("/samples/{sample_id}/variants/{variant_id}/observed")
def get_observed_cases(sample_id: str, variant_id: str):
    meta = _sample_meta(sample_id)
    try:
        cases = manual_acmg.observed_cases(
            meta.get("genome_build") or "hg38",
            variant_id,
            exclude_sample_id=sample_id,
        )
    except ValueError as exc:
        raise HTTPException(400, str(exc)) from exc
    return {
        "variant_id": manual_acmg.normalize_variant_id(variant_id),
        "cases": cases,
        "count": len(cases),
    }
