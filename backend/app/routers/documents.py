from __future__ import annotations

from fastapi import APIRouter, Depends, File, Form, HTTPException, Query, UploadFile
from fastapi.responses import FileResponse, Response, StreamingResponse

from ..auth import current_user
from ..services import patient_documents


router = APIRouter(
    prefix="/api/documents",
    tags=["documents"],
    dependencies=[Depends(current_user)],
)


def _raise(exc: Exception) -> None:
    if isinstance(exc, patient_documents.InvalidDocument):
        raise HTTPException(400, str(exc)) from exc
    if isinstance(exc, patient_documents.DocumentConflict):
        raise HTTPException(409, str(exc)) from exc
    if isinstance(exc, patient_documents.DocumentNotFound):
        raise HTTPException(404, str(exc)) from exc
    if isinstance(exc, patient_documents.DocumentStorageFull):
        raise HTTPException(507, str(exc)) from exc
    if isinstance(exc, patient_documents.PreviewUnavailable):
        raise HTTPException(422, str(exc)) from exc
    if isinstance(exc, OSError):
        raise HTTPException(500, f"文件儲存失敗：{exc}") from exc
    raise exc


@router.get("")
def list_patient_documents(mrn: str = Query(...)):
    try:
        return patient_documents.list_documents(mrn)
    except Exception as exc:
        _raise(exc)


@router.post("")
async def upload_patient_document(
    mrn: str = Form(...),
    display_name: str = Form(""),
    source_sample_id: str = Form(""),
    file: UploadFile = File(...),
    user: dict = Depends(current_user),
):
    try:
        return await patient_documents.save_upload(
            file,
            mrn=mrn,
            display_name=display_name,
            source_sample_id=source_sample_id,
            user=user,
        )
    except Exception as exc:
        _raise(exc)


@router.get("/archive.zip")
def download_patient_document_archive(
    mrn: str = Query(...),
    user: dict = Depends(current_user),
):
    try:
        safe_mrn = patient_documents.validate_mrn(mrn)
        content, count = patient_documents.stream_archive(safe_mrn, user=user)
        return StreamingResponse(
            content,
            media_type="application/zip",
            headers={
                "Content-Disposition": f'attachment; filename="{safe_mrn}_documents.zip"',
                "Cache-Control": "private, no-store",
                "X-Content-Type-Options": "nosniff",
                "X-Document-Count": str(count),
            },
        )
    except Exception as exc:
        _raise(exc)


@router.patch("/{document_id}")
def rename_patient_document(
    document_id: str,
    payload: dict,
    user: dict = Depends(current_user),
):
    try:
        return patient_documents.rename_document(
            document_id,
            str((payload or {}).get("display_name") or ""),
            user=user,
        )
    except Exception as exc:
        _raise(exc)


@router.delete("/{document_id}")
def delete_patient_document(
    document_id: str,
    user: dict = Depends(current_user),
):
    try:
        return patient_documents.delete_document(document_id, user=user)
    except Exception as exc:
        _raise(exc)


@router.get("/{document_id}/download")
def download_patient_document(
    document_id: str,
    user: dict = Depends(current_user),
):
    try:
        path, info = patient_documents.document_file(
            document_id,
            user=user,
            action="download",
        )
        return FileResponse(
            path,
            filename=info["display_name"],
            media_type=info["content_type"],
            headers={
                "Cache-Control": "private, no-store",
                "X-Content-Type-Options": "nosniff",
            },
        )
    except Exception as exc:
        _raise(exc)


@router.get("/{document_id}/preview")
def preview_patient_document(
    document_id: str,
    page: int = Query(0, ge=0),
    user: dict = Depends(current_user),
):
    try:
        content, info = patient_documents.render_preview(
            document_id,
            page=page,
            user=user,
        )
        return Response(
            content=content,
            media_type="image/png",
            headers={
                "Cache-Control": "private, no-store",
                "X-Content-Type-Options": "nosniff",
                "X-Image-Pages": str(info.get("image_pages") or 1),
            },
        )
    except Exception as exc:
        _raise(exc)
