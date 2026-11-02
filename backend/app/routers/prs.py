"""Authenticated CHD research estimates for an already registered sample."""
from fastapi import APIRouter, Depends, HTTPException

from ..auth import current_user
from ..services import prs

router = APIRouter(prefix='/api/samples/{sample_id}/prs', tags=['prs'],
                   dependencies=[Depends(current_user)])


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


@router.post('/assessments')
def assess(sample_id: str, payload: dict, user=Depends(current_user)):
    return checked(prs.assess, sample_id, payload, user.get('username', ''))
