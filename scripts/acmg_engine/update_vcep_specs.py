#!/usr/bin/env python3
"""Rebuild the rule engine's VCEP specification and collagen tables.

Downloads ClinGen CSpec criteria specifications (Released / Approved For
Release) and UniProt sequences of the collagen genes, then writes, inside the
engine package (so the result signature changes and stored results become
stale)::

    backend/app/services/acmg_engine/vcep_specs.json
    backend/app/services/acmg_engine/collagen_gly.json

and a review table for people to check what was read from the free text::

    docs/research/acmg_vcep_specs_review.tsv

    .venv/bin/python scripts/acmg_engine/update_vcep_specs.py
    .venv/bin/python scripts/acmg_engine/update_vcep_specs.py --cache-dir /tmp/cspec   # reuse downloads

Needs network access to cspec.genome.network and rest.uniprot.org.  After an
update, recompute engine results (build_engine_results.py --all --force).
"""
from __future__ import annotations

import argparse
import concurrent.futures as cf
import csv
import json
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "backend"))

from app.services.acmg_engine import cspec  # noqa: E402

PACKAGE = REPO / "backend/app/services/acmg_engine"
CSPEC_LIST = "https://cspec.genome.network/cspec/api/svis?detail=high"
UNIPROT = ("https://rest.uniprot.org/uniprotkb/search?query={q}"
           "&fields=accession,gene_primary,length,ft_region,ft_chain,sequence&format=json")
COLLAGEN_GENES = ("COL1A1 COL1A2 COL2A1 COL3A1 COL4A1 COL4A2 COL4A3 COL4A4 COL4A5 COL4A6 COL5A1 COL5A2 "
                  "COL6A1 COL6A2 COL6A3 COL7A1 COL9A1 COL9A2 COL9A3 COL10A1 COL11A1 COL11A2 COL17A1").split()
USER_AGENT = "ngs-ui-acmg-engine/1.0"


def _get(url: str, cache: Path | None, name: str) -> bytes:
    path = cache / name if cache else None
    if path and path.is_file():
        return path.read_bytes()
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT, "Accept": "application/json"})
    for attempt in range(4):
        try:
            with urllib.request.urlopen(request, timeout=90) as response:
                data = response.read()
            break
        except OSError:
            if attempt == 3:
                raise
            time.sleep(2 ** (attempt + 1))
    if path:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(data)
    return data


def fetch_cspec(cache: Path | None) -> list[dict]:
    listing = json.loads(_get(CSPEC_LIST, cache, "svis.json"))["data"]
    wanted = [x for x in listing if x.get("status") in cspec.USABLE_STATUS]

    def one(item):
        sid = item["@id"].rsplit("/", 1)[-1]
        return json.loads(_get(item["@id"], cache, f"svi/{sid}.json"))

    with cf.ThreadPoolExecutor(6) as pool:
        return list(pool.map(one, wanted))


def fetch_collagen(cache: Path | None) -> dict:
    out = {}
    for gene in COLLAGEN_GENES:
        q = urllib.parse.quote(f"gene_exact:{gene} AND organism_id:9606 AND reviewed:true")
        data = json.loads(_get(UNIPROT.format(q=q), cache, f"uniprot/{gene}.json"))
        results = data.get("results") or []
        if not results:
            print(f"[warn] UniProt: {gene} not found", file=sys.stderr)
            continue
        entry = cspec.uniprot_collagen_entry(results[0])
        if entry["positions"]:
            out[gene] = entry
    return out


def write_review(path: Path, specs: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["gene", "svi", "version", "vcep", "not_applicable", "BA1", "BS1", "PM2",
                         "PM1_rules", "PM1_unparsed", "threshold_notes", "pm1_text"])
        for gene, s in sorted(specs["genes"].items()):
            na = ",".join(code for code, allowed in s["criteria"].items() if not allowed)
            thr = s["thresholds"]
            pm1 = "; ".join(f"{r['strength']}:" + ",".join(f"{a}-{b}" if a != b else str(a) for a, b in r["ranges"])
                            for r in s["pm1"])
            writer.writerow([gene, s["svi"], s["version"], s["vcep"], na,
                             thr.get("ba1", ""), thr.get("bs1", ""), thr.get("pm2", ""), pm1,
                             "yes" if s["pm1_unparsed"] else "",
                             " | ".join(f"{k}: {v}" for k, v in s["threshold_notes"].items()),
                             " | ".join(f"{k}: {v}" for k, v in s["pm1_text"].items())])
        for gene, reason in sorted(specs["skipped"].items()):
            writer.writerow([gene, "", "", f"略過：{reason}", "", "", "", "", "", "", "", ""])


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--cache-dir", type=Path, default=None, help="keep/reuse raw downloads here")
    ap.add_argument("--review", type=Path, default=REPO / "docs/research/acmg_vcep_specs_review.tsv")
    args = ap.parse_args(argv)
    stamp = time.strftime("%Y-%m-%d")

    documents = fetch_cspec(args.cache_dir)
    specs = cspec.build_vcep_specs(documents)
    specs = {"source": {"name": "ClinGen Criteria Specification Registry", "url": CSPEC_LIST,
                        "retrieved": stamp, "documents": len(documents)}, **specs}
    (PACKAGE / "vcep_specs.json").write_text(json.dumps(specs, ensure_ascii=False, indent=1, sort_keys=True) + "\n",
                                            encoding="utf-8")
    collagen = {"source": {"name": "UniProtKB (reviewed, human)", "retrieved": stamp,
                           "min_gly_xy_run": cspec.MIN_GXY_RUN}, "genes": fetch_collagen(args.cache_dir)}
    (PACKAGE / "collagen_gly.json").write_text(json.dumps(collagen, sort_keys=True) + "\n", encoding="utf-8")
    write_review(args.review, specs)
    print(f"VCEP specs: {len(specs['genes'])} genes from {len(documents)} documents "
          f"(skipped {len(specs['skipped'])}); collagen genes: {len(collagen['genes'])}")
    print(f"review table: {args.review}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
