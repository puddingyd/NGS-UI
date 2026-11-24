"""ClinVar affected-case counts for PS4 B (proband counting, rulebook ``ps4``).

Reads the ClinVar VCV XML release (``ClinVarVCVRelease_00-latest.xml.gz``) and
keeps, per Variation ID, the submitters that reported the variant in an
affected individual:

* the submission (ClinicalAssertion) is current and not classified benign /
  likely benign by that submitter;
* an ObservedIn sample is germline (not somatic) with ``AffectedStatus`` yes;
* its method is one of ``methods`` (clinical testing, research, provider
  interpretation, case-control by default): literature-only and curation
  submissions repeat published cases and are left out.

Output ``clinvar_cases.sqlite`` (``cases``: variation_id, submitters,
individuals, submitter_ids; ``meta``).  The engine caps this evidence at
Supporting because case details and independence cannot be checked.
"""
from __future__ import annotations

import gzip
import os
import shutil
import sqlite3
import subprocess
import time
import xml.etree.ElementTree as ET
from pathlib import Path

CASES_DB_NAME = "clinvar_cases.sqlite"
XML_NAME = "ClinVarVCVRelease_00-latest.xml.gz"
XML_URL = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/" + XML_NAME
DEFAULT_METHODS = ("clinical testing", "research", "provider interpretation", "case-control")


def _local(tag) -> str:
    tag = tag if isinstance(tag, str) else ""
    return tag.rsplit("}", 1)[-1]


def _children(elem, name: str):
    return [c for c in elem if _local(c.tag) == name]


def _descendants(elem, name: str):
    return [c for c in elem.iter() if _local(c.tag) == name]


def _text(elem, name: str) -> str:
    for c in _descendants(elem, name):
        if c.text and c.text.strip():
            return c.text.strip()
    return ""


def _classification(assertion) -> str:
    """Submitter's germline classification (schema 2: Classification/GermlineClassification;
    schema 1: Interpretation/Description)."""
    for c in _children(assertion, "Classification"):
        value = _text(c, "GermlineClassification")
        if value:
            return value
    for c in _children(assertion, "Interpretation"):
        value = _text(c, "Description")
        if value:
            return value
    return ""


def _assertion_cases(assertion, methods: set) -> int | None:
    """Affected individuals (≥1) reported by one submission, or None if it does not qualify."""
    status = ""
    for c in _children(assertion, "RecordStatus"):
        status = (c.text or "").strip().lower()
    if status and status != "current":
        return None
    cls = _classification(assertion).lower()
    if "benign" in cls and "pathogenic" not in cls:
        return None
    best = None
    for observed in _descendants(assertion, "ObservedIn"):
        samples = _children(observed, "Sample")
        sample = samples[0] if samples else None
        if sample is None:
            continue
        if _text(sample, "AffectedStatus").lower() != "yes":
            continue
        origin = _text(sample, "Origin").lower()
        if "somatic" in origin:
            continue
        found = {(_text(m, "MethodType") or "").lower() for m in _descendants(observed, "Method")}
        if not found & methods:
            continue
        number = _text(sample, "NumberOfIndividuals")
        count = int(number) if number.isdigit() and int(number) > 0 else 1
        best = max(best or 0, count)
    return best


def iter_cases(stream, methods=DEFAULT_METHODS):
    """Yield (variation_id, {submitter: individuals}) per VariationArchive with affected cases."""
    methods = {m.lower() for m in methods}
    variation = ""
    found: dict[str, int] = {}
    root = None
    for event, elem in ET.iterparse(stream, events=("start", "end")):
        name = _local(elem.tag)
        if event == "start":
            if root is None:
                root = elem
            if name == "VariationArchive":
                variation = elem.get("VariationID") or ""
                found = {}
            continue
        if name == "ClinicalAssertion":
            cases = _assertion_cases(elem, methods)
            if cases:
                acc = _children(elem, "ClinVarAccession")
                who = (acc[0].get("OrgID") or acc[0].get("SubmitterName") or acc[0].get("Accession")) if acc else ""
                who = who or f"scv{len(found)}"
                found[who] = max(found.get(who, 0), cases)
            elem.clear()
        elif name == "VariationArchive":
            if variation and found:
                yield variation, found
            elem.clear()
            if root is not None:
                root.clear()


def open_xml(path: Path):
    """Decompressed stream; uses an external gzip when available (faster)."""
    path = Path(path)
    if str(path).endswith(".gz"):
        gz = shutil.which("pigz") or shutil.which("gzip")
        if gz:
            proc = subprocess.Popen([gz, "-dc", str(path)], stdout=subprocess.PIPE)
            return proc.stdout, proc
        return gzip.open(path, "rb"), None
    return open(path, "rb"), None


def build(xml_path: Path, out: Path, *, release: str = "", methods=DEFAULT_METHODS, progress=print) -> dict:
    """Parse ``xml_path`` and atomically write ``out``; returns stats."""
    started = time.time()
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = out.with_name(f".{out.name}.{os.getpid()}.tmp")
    tmp.unlink(missing_ok=True)
    conn = sqlite3.connect(tmp)
    stream, proc = open_xml(xml_path)
    variants = 0
    try:
        conn.executescript("""
            PRAGMA journal_mode=OFF; PRAGMA synchronous=OFF;
            CREATE TABLE meta (key TEXT PRIMARY KEY, value TEXT NOT NULL);
            CREATE TABLE cases (variation_id TEXT PRIMARY KEY, submitters INTEGER NOT NULL,
                                individuals INTEGER NOT NULL, submitter_ids TEXT NOT NULL) WITHOUT ROWID;
        """)
        batch = []
        for variation, found in iter_cases(stream, methods):
            batch.append((variation, len(found), sum(found.values()), ",".join(sorted(found))))
            variants += 1
            if len(batch) >= 5000:
                conn.executemany("INSERT OR REPLACE INTO cases VALUES (?,?,?,?)", batch)
                batch.clear()
                progress(f"  ClinVar 受影響個案：{variants:,} 個變異")
        conn.executemany("INSERT OR REPLACE INTO cases VALUES (?,?,?,?)", batch)
        meta = {"schema_version": "1", "release": release, "source": str(xml_path),
                "methods": ",".join(methods), "variants": str(variants),
                "built_at": time.strftime("%Y-%m-%dT%H:%M:%S")}
        conn.executemany("INSERT INTO meta VALUES (?,?)", sorted(meta.items()))
        conn.commit()
    except BaseException:
        conn.close()
        tmp.unlink(missing_ok=True)
        raise
    finally:
        if proc is not None:
            stream.close()
            proc.wait()
    conn.close()
    if proc is not None and proc.returncode not in (0, None):
        tmp.unlink(missing_ok=True)
        raise RuntimeError(f"解壓縮 {xml_path} 失敗（exit {proc.returncode}）")
    os.replace(tmp, out)
    return {"variants": variants, "seconds": round(time.time() - started, 1), "path": str(out)}
