#!/usr/bin/env python3
"""Probe the NCKU EMR APIs from VIP_API.sh and dump shapes for design.

Reads two internal endpoints — GetPhenotypeList (phenotype API,
HTTP, returns broken JSON) and easyform/getdata (consultation
records, HTTPS, IBM APIM gateway) — for a list of MRNs. For each
MRN it:

  1. saves the raw HTTP response body to /tmp/emr_probe/<MRN>_*.{txt,json}
  2. attempts JSON parsing (with the .sh's heuristic repair as a
     fallback for the phenotype API)
  3. prints a one-line summary to stdout: HTTP status / size / parse
     result / sample fields seen at the top level

Run on the hospital intranet (192.168.84.91) with the same Python
that has `requests` available — i.e. inside the project venv:

    /home/n102968/NGS_UI/NGS-UI/.venv/bin/python3 \\
        scripts/probe_emr_api.py

Edit MRNS below if you want to add or skip any. The IBM APIM client
id is the same one the legacy VIP_API.sh has hard-coded; if it
expires, replace it inline for the probe.

The script writes everything under /tmp/emr_probe so it's easy to
clean up afterwards (`rm -rf /tmp/emr_probe`).
"""
from __future__ import annotations

import json
import re
import sys
import time
from pathlib import Path

import requests


PHENO_URL = "http://hisweb.hosp.ncku/hisservice/opd/nckuhisweb/aspx/DelegateExamServiceGate.aspx/GetPhenotypeList"
CONSULT_URL = "https://apigw-i.apim.hosp.ncku.edu.tw/rd/prod-i/easyform/getdata"
CONSULT_TCODE = "EMR-3-GC-002"
APIM_CLIENT_ID = "9c03b0c83c562ffa22d1b4ff0e54d41d"

MRNS = [
    "18281656", "23197691", "18061494", "22433814", "20124265",
    "23163111", "23188065", "23127742", "23051665", "22263217",
    "21843765", "18233518", "23243811", "15986525",
]

OUT_DIR = Path("/tmp/emr_probe")


# ---- ported from VIP_API.sh's try_fix_json --------------------------
def try_fix_json(raw: str):
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
    except json.JSONDecodeError as e:
        return f"<repair failed: {e}>"


def fetch_phenotype(mrn: str) -> dict:
    """Returns {http_status, raw_size, parse, parsed} for one MRN."""
    payload = {"JasonInputValue": json.dumps({"ChartNo": mrn})}
    try:
        r = requests.post(PHENO_URL, data=payload, timeout=15)
    except Exception as exc:
        return {"http_status": "EXC", "raw_size": 0, "parse": "exc",
                "parsed": str(exc)}
    raw = r.text
    (OUT_DIR / f"{mrn}_phenotype.txt").write_text(raw, encoding="utf-8")
    # The .sh splits on \n\n then escapes inner newlines. Replicate.
    json_part = raw.replace("\r", "").split("\n\n")[0].replace("\n", "\\n")
    try:
        parsed = json.loads(json_part, strict=False)
        return {"http_status": r.status_code, "raw_size": len(raw),
                "parse": "ok", "parsed": parsed}
    except json.JSONDecodeError as e:
        repaired = try_fix_json(json_part)
        return {"http_status": r.status_code, "raw_size": len(raw),
                "parse": "repaired" if isinstance(repaired, (list, dict))
                                     else "broken",
                "parsed": repaired, "json_error": str(e)}


def fetch_consultation(mrn: str) -> dict:
    body = json.dumps({"chartNo": mrn, "tcode": CONSULT_TCODE})
    headers = {
        "Content-Type":   "application/json",
        "Accept":         "application/json",
        "X-IBM-Client-Id": APIM_CLIENT_ID,
    }
    try:
        r = requests.post(CONSULT_URL, data=body, headers=headers, timeout=15)
    except Exception as exc:
        return {"http_status": "EXC", "raw_size": 0, "parse": "exc",
                "parsed": str(exc)}
    (OUT_DIR / f"{mrn}_consultation.json").write_text(r.text, encoding="utf-8")
    if r.status_code != 200:
        return {"http_status": r.status_code, "raw_size": len(r.text),
                "parse": "http", "parsed": r.text[:500]}
    try:
        return {"http_status": r.status_code, "raw_size": len(r.text),
                "parse": "ok", "parsed": r.json()}
    except json.JSONDecodeError as e:
        return {"http_status": r.status_code, "raw_size": len(r.text),
                "parse": "broken", "parsed": r.text[:500],
                "json_error": str(e)}


def _shape(v) -> str:
    """Compact 'shape' summary for sanity-checking responses at a glance."""
    if isinstance(v, dict):
        keys = list(v.keys())
        return "{" + ", ".join(keys[:8]) + ("…" if len(keys) > 8 else "") + "}"
    if isinstance(v, list):
        if not v:
            return "[]"
        return f"[{len(v)}× {_shape(v[0])}]"
    if isinstance(v, str):
        return f"<str len={len(v)}>"
    return type(v).__name__


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    print(f"writing raw bodies to {OUT_DIR}/")
    print()
    print(f"{'MRN':12}  {'pheno status':14}  {'pheno parse':12}  {'pheno shape':40}  "
          f"{'consult status':14}  {'consult parse':12}  consult shape")
    print("-" * 160)
    for mrn in MRNS:
        p = fetch_phenotype(mrn)
        c = fetch_consultation(mrn)
        p_shape = _shape(p.get("parsed"))[:38]
        c_shape = _shape(c.get("parsed"))[:60]
        print(f"{mrn:12}  {str(p['http_status']):14}  {p['parse']:12}  "
              f"{p_shape:40}  "
              f"{str(c['http_status']):14}  {c['parse']:12}  {c_shape}")
        time.sleep(0.3)  # tiny gap; legacy notes the API misbehaves on bursts
    print()
    print(f"Inspect a single sample:")
    print(f"  cat {OUT_DIR}/<MRN>_phenotype.txt   | head -40")
    print(f"  cat {OUT_DIR}/<MRN>_consultation.json | python3 -m json.tool | head -80")
    print()
    print("When you paste back, redact name/dob/etc as needed; we mainly")
    print("need to see the field names + types, not the values.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
