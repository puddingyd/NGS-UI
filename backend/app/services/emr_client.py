"""NCKU intranet EMR fetch + parse.

Two upstream APIs (mirrors VIP_API.sh from the legacy system):

  GetPhenotypeList — http://hisweb.hosp.ncku/...
    POST form-urlencoded with JasonInputValue={"ChartNo": <mrn>}.
    Returns broken JSON (response.text mixes a leading payload with
    a trailing ASP.NET viewstate HTML form). Strip the HTML and
    json.loads on the leading chunk; fall back to the heuristic
    bracket/comma repair if that fails. Empty-result sentinel: [{}].

  easyform/getdata — https://apigw-i.apim.hosp.ncku.edu.tw/...
    POST JSON {"chartNo": <mrn>, "tcode": "EMR-3-GC-002"} with
    X-IBM-Client-Id header. Returns clean JSON. Empty-result: [].

The intranet-only nature plus the Client-Id quirks means everything
is best-effort: any fetch returns a dict with `found: bool`, the UI
shows a hint when the call yields nothing rather than blocking.
"""
from __future__ import annotations

import json
import re
import ssl
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from typing import Any

from ..config import EMR_CLIENT_ID


PHENO_URL   = "http://hisweb.hosp.ncku/hisservice/opd/nckuhisweb/aspx/DelegateExamServiceGate.aspx/GetPhenotypeList"
CONSULT_URL = "https://apigw-i.apim.hosp.ncku.edu.tw/rd/prod-i/easyform/getdata"
CONSULT_TCODE = "EMR-3-GC-002"

_HP_ID_RE = re.compile(r"\bHP:\d{7}\b")
_PANEL_HINT_RE = re.compile(r"panel|panelapp|\(Version\s+\d", re.IGNORECASE)


def is_enabled() -> bool:
    return bool(EMR_CLIENT_ID)


def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def _post(url: str, body: bytes, headers: dict, timeout: int = 15) -> tuple[int, str]:
    """Minimal urllib POST that surfaces 4xx/5xx bodies as data."""
    req = urllib.request.Request(url, data=body, headers=headers, method="POST")
    ctx = ssl.create_default_context()
    ctx.check_hostname = False
    ctx.verify_mode = ssl.CERT_NONE
    try:
        with urllib.request.urlopen(req, timeout=timeout, context=ctx) as resp:
            return resp.status, resp.read().decode("utf-8", errors="replace")
    except urllib.error.HTTPError as e:
        return e.code, e.read().decode("utf-8", errors="replace")


# ---- broken-JSON repair (ported from VIP_API.sh::try_fix_json) ------

def _try_fix_json(raw: str) -> Any:
    cleaned = re.split(r"</html>|<html", raw, flags=re.IGNORECASE)[0].strip()
    cleaned = re.sub(r'([}\]])\s*"(?=\w+"\s*:)', r'\1, "', cleaned)
    cleaned = re.sub(r"([}\]])\s*(?=[{\[])", r"\1, ", cleaned)
    cleaned = re.sub(r",\s*([}\]])", r"\1", cleaned)

    stack: list[str] = []
    out: list[str] = []
    quote = False
    esc = False
    for ch in cleaned:
        if ch == "\\" and not esc:
            esc = True
            out.append(ch); continue
        if ch == '"' and not esc:
            quote = not quote
        if not quote:
            if ch in "{[":
                stack.append(ch)
            elif ch in "}]":
                if not stack:
                    esc = False; continue
                last = stack.pop()
                if (last == "{" and ch != "}") or (last == "[" and ch != "]"):
                    stack.append(last)
                    esc = False; continue
        esc = False
        out.append(ch)
    while stack:
        out.append("}" if stack.pop() == "{" else "]")
    try:
        return json.loads("".join(out))
    except json.JSONDecodeError:
        return None


# ---- phenotype content parser ---------------------------------------

def parse_phenotype_content(text: str) -> tuple[list[dict], list[str]]:
    """Split the multi-line `phenotypes[].content` into HPO terms + notes.

    Returns (hpo_list, notes). hpo_list entries are
    {phenotype: HP:..., label: human-name, weight: 1}. Lines without a
    HP:nnnnnnn pattern (panel references, free-text annotations like
    'Growth panel' or '無其他pheotype') land in notes verbatim so the
    UI can surface them without losing information.
    """
    hpo: list[dict] = []
    notes: list[str] = []
    for line in (text or "").splitlines():
        line = line.strip()
        if not line:
            continue
        m = _HP_ID_RE.search(line)
        if m:
            hp_id = m.group(0)
            # Strip the HP: token to leave the human-readable label.
            label = _HP_ID_RE.sub("", line).strip(" \t·-:：")
            # Some entries put HP:ID at the beginning with extra spaces;
            # collapse any internal whitespace runs.
            label = re.sub(r"\s+", " ", label).strip()
            hpo.append({
                "phenotype": hp_id,
                "label":     label or hp_id,
                "weight":    1,
            })
        else:
            notes.append(line)
    return hpo, notes


# ---- fetchers --------------------------------------------------------

def fetch_phenotype(mrn: str) -> dict:
    """Hit GetPhenotypeList for one MRN and shape the result.

    Returns:
        {found: bool, hpo: [...], notes: [...], date: "YYYY/MM/DD",
         raw_content: str, error: str?}
    """
    if not mrn:
        return {"found": False, "hpo": [], "notes": [], "date": "", "raw_content": ""}
    form = urllib.parse.urlencode({
        "JasonInputValue": json.dumps({"ChartNo": mrn}),
    }).encode("utf-8")
    headers = {"Content-Type": "application/x-www-form-urlencoded"}
    try:
        status, raw = _post(PHENO_URL, form, headers)
    except Exception as exc:
        return {"found": False, "hpo": [], "notes": [], "date": "",
                "raw_content": "", "error": f"fetch failed: {exc}"}
    json_part = raw.replace("\r", "").split("\n\n")[0].replace("\n", "\\n")
    try:
        parsed = json.loads(json_part, strict=False)
    except json.JSONDecodeError:
        parsed = _try_fix_json(json_part)
        if not isinstance(parsed, list):
            return {"found": False, "hpo": [], "notes": [], "date": "",
                    "raw_content": "", "error": "phenotype JSON unparseable"}
    # Empty-result sentinel: [{}]
    if not parsed or not isinstance(parsed[0], dict) or not parsed[0]:
        return {"found": False, "hpo": [], "notes": [], "date": "",
                "raw_content": ""}
    entry = parsed[0]
    content = ""
    for p in entry.get("phenotypes", []) or []:
        if isinstance(p, dict) and p.get("content"):
            content = p["content"]
            break
    hpo, notes = parse_phenotype_content(content)
    return {
        "found":       bool(content.strip()),
        "hpo":         hpo,
        "notes":       notes,
        "date":        entry.get("date", ""),
        "raw_content": content,
    }


def fetch_consultation(mrn: str) -> dict:
    """Hit easyform/getdata for one MRN and shape the result.

    Returns:
        {found, gender_raw ('男'/'女'), sex ('M'/'F'/''),
         date_of_birth, records: [{date_of_consult, reason, diagnosis, record}],
         text: combined dedupe'd record text, error?}
    """
    if not mrn:
        return {"found": False}
    if not is_enabled():
        return {"found": False, "error": "EMR client_id not configured"}
    body = json.dumps({"chartNo": mrn, "tcode": CONSULT_TCODE}).encode("utf-8")
    headers = {
        "Content-Type":    "application/json",
        "Accept":          "application/json",
        "X-IBM-Client-Id": EMR_CLIENT_ID,
    }
    try:
        status, raw = _post(CONSULT_URL, body, headers)
    except Exception as exc:
        return {"found": False, "error": f"fetch failed: {exc}"}
    if status != 200:
        return {"found": False, "error": f"HTTP {status}", "raw": raw[:200]}
    try:
        parsed = json.loads(raw)
    except json.JSONDecodeError as e:
        return {"found": False, "error": f"JSON parse: {e}"}
    if not parsed or not isinstance(parsed[0], dict):
        return {"found": False}
    entry = parsed[0]
    gender_raw = (entry.get("gender") or "").strip()
    sex = {"男": "M", "女": "F"}.get(gender_raw, "")
    consult = entry.get("consult") or []
    # dedupe by `record` content; ignore empty records.
    seen: set[str] = set()
    records: list[dict] = []
    for r in consult:
        if not isinstance(r, dict):
            continue
        rec = (r.get("record") or "").strip()
        if not rec or rec in seen:
            continue
        seen.add(rec)
        records.append({
            "date_of_consult": r.get("date_of_consult", ""),
            "reason":          r.get("reason", ""),
            "diagnosis":       r.get("diagnosis", ""),
            "record":          rec,
        })
    # Combined text: each record on its own block, joined by separator.
    chunks = []
    for r in records:
        chunks.append(r["record"])
        diag = (r.get("diagnosis") or "").strip()
        if diag and diag not in r["record"]:
            chunks.append(f"--- 診斷 ---\n{diag}")
    text = "\n\n=====\n\n".join(chunks)
    return {
        "found":         bool(records),
        "gender_raw":    gender_raw,
        "sex":           sex,
        "date_of_birth": entry.get("date_of_birth", ""),
        "records":       records,
        "text":          text,
    }


def fetch(mrn: str) -> dict:
    """Combined helper. Returns both API results in one dict."""
    return {
        "mrn":          mrn,
        "phenotype":    fetch_phenotype(mrn),
        "consultation": fetch_consultation(mrn),
        "fetched_at":   _now(),
    }
