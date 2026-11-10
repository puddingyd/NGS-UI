"""Authenticated PRS estimates and shared methodology documentation."""
from fastapi import APIRouter, Depends, HTTPException
from fastapi.responses import Response

from ..auth import current_user
from .. import config
from ..services import prs

router = APIRouter(prefix='/api/samples/{sample_id}/prs', tags=['prs'],
                   dependencies=[Depends(current_user)])
documentation_router = APIRouter(prefix='/api/prs', tags=['prs'],
                                dependencies=[Depends(current_user)])


@documentation_router.get('/documentation')
def documentation():
    path = config.REPO_ROOT / 'docs/research/prs/PRS_WORKFLOW.md'
    try:
        content = path.read_bytes()
    except OSError as exc:
        raise HTTPException(503, 'PRS 說明文件暫時無法讀取') from exc
    return Response(content, media_type='text/markdown; charset=utf-8', headers={
        'Content-Disposition': 'attachment; filename="PRS_WORKFLOW.md"',
        'Cache-Control': 'no-store',
    })


def checked(call, *args):
    try:
        return call(*args)
    except FileNotFoundError as exc:
        raise HTTPException(404, '找不到個案或必要資料') from exc
    except prs.NotReady as exc:
        raise HTTPException(409, str(exc)) from exc
    except (ValueError, KeyError, TypeError) as exc:
        raise HTTPException(400, str(exc) if isinstance(exc, ValueError) else 'PRS 輸入格式不正確') from exc
    except OSError as exc:
        raise HTTPException(503, 'PRS 資料暫時無法讀寫') from exc


@router.get('')
def status(sample_id: str):
    return checked(prs.status, sample_id)


@router.get('/report-readiness')
def report_readiness(sample_id: str):
    checked(prs.case_metadata, sample_id)
    return prs.report_readiness(sample_id)


@router.post('/assessments')
def assess(sample_id: str, payload: dict, user=Depends(current_user)):
    return checked(prs.assess, sample_id, payload, user.get('username', ''))


@router.post('/jobs')
def start_score(sample_id: str, payload: dict, user=Depends(current_user)):
    from ..services import prs_jobs
    checked(prs.case_metadata, sample_id)
    if set(payload) - {'candidate_id', 'new_run'}:
        raise HTTPException(400, 'PRS 只接受伺服器列出的 gVCF 選項')
    if type(payload.get('new_run', False)) is not bool:
        raise HTTPException(400, 'new_run 必須是布林值')
    return checked(prs_jobs.start, sample_id, payload.get('candidate_id'), user.get('username', ''), None, payload.get('new_run', False))


@router.post('/jobs/{job_id}/cancel')
def cancel_score(sample_id: str, job_id: str):
    from ..services import prs_jobs
    checked(prs.case_metadata, sample_id)
    return checked(prs_jobs.cancel, sample_id, job_id)


@router.get('/jobs/{job_id}/log')
def score_log(sample_id: str, job_id: str):
    from ..services import prs_jobs
    checked(prs.case_metadata, sample_id)
    return checked(prs_jobs.log, sample_id, job_id)


@router.get('/jobs/{job_id}/raw-log')
def raw_score_log(sample_id: str, job_id: str):
    from ..services import prs_jobs
    checked(prs.case_metadata, sample_id)
    path = checked(prs_jobs.raw_log_path, sample_id, job_id)
    return Response(checked(path.read_bytes), media_type='text/plain; charset=utf-8',
                    headers={'Content-Disposition': f'attachment; filename="prs-{job_id}.log"'})
