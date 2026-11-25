"""Disease → gene links for the phenotype engine (plan §4.2.6).

* Platform OMIM: phenotype MIM ↔ gene from the OMIM rows ``omim_store`` serves
  (DB or xlsx).  ``?`` (provisional) and ``{}`` (susceptibility) are used with a
  flag; ``[]`` (non-disease traits) and somatic entries are kept but not used.
* Orphanet product6: germline disease-causing associations are used; somatic,
  biomarker, candidate-gene and fusion-gene associations are kept but not
  used; other types (role in phenotype, susceptibility factor, modifier) are
  flagged ``review`` and not used until validated.
* HPO ``genes_to_disease.txt``: cross-check only; disagreements go to QC.

Gene symbols are brought to the platform's HGNC-current names.
"""
from __future__ import annotations

import re
import xml.etree.ElementTree as ET
from pathlib import Path

_MIM = re.compile(r"\((\d{6})\)")
_INH = re.compile(r"\(([^()]*)\)\s*$")

ORPHA_USE = {
    "Disease-causing germline mutation(s) in",
    "Disease-causing germline mutation(s) (loss of function) in",
    "Disease-causing germline mutation(s) (gain of function) in",
}
ORPHA_EXCLUDE = {
    "Disease-causing somatic mutation(s) in",
    "Biomarker tested in",
    "Candidate gene tested in",
    "Part of a fusion gene in",
}

COLUMNS = ("disease_id", "gene_symbol", "source", "association_type", "use", "flags", "status", "source_gene")


def _canonical(gene: str) -> str:
    gene = (gene or "").strip()
    try:
        from .. import panel_deadzone
        return panel_deadzone.canonical_panel_gene_symbol(gene) or gene
    except Exception:  # noqa: BLE001
        return gene


def omim_links(rows, canonical=_canonical) -> list[dict]:
    """Links from OMIM gene rows (``gene_symbol`` + multi-line ``OMIM_disease``)."""
    out = {}
    for row in rows:
        symbol = str(row.get("gene_symbol") or "").strip()
        if not symbol:
            continue
        gene = canonical(symbol)
        for line in str(row.get("OMIM_disease") or "").splitlines():
            text = line.strip()
            m = _MIM.search(text)
            if not text or not m:
                continue
            flags, use, kind = [], 1, "phenotype"
            if text.startswith("["):
                use, kind = 0, "non_disease"
                flags.append("non_disease")
            if text.startswith("{") or "{" in text[:2]:
                flags.append("susceptibility")
            if text.lstrip("{[").startswith("?"):
                flags.append("provisional")
            inh = _INH.search(text)
            inheritance = inh.group(1) if inh and not inh.group(1).isdigit() else ""
            if "SMu" in inheritance.split(",") or "SMu" in inheritance or "somatic" in text.lower():
                use, kind = 0, "somatic"
                flags.append("somatic")
            key = (f"OMIM:{m.group(1)}", gene)
            prev = out.get(key)
            rec = {"disease_id": key[0], "gene_symbol": gene, "source": "OMIM", "association_type": kind,
                   "use": use, "flags": ",".join(sorted(set(flags))), "status": "", "source_gene": symbol}
            if prev is None or (rec["use"] and not prev["use"]):
                out[key] = rec
    return [out[k] for k in sorted(out)]


def orphanet_links(path: Path, canonical=_canonical) -> tuple[list[dict], dict]:
    """Links and ``{"version", "date"}`` from Orphanet product6 XML (``en_product6.xml``)."""
    info = {"version": "", "date": ""}
    out = {}
    for event, elem in ET.iterparse(str(path), events=("start", "end")):
        if event == "start" and elem.tag == "JDBOR":
            info["version"] = elem.attrib.get("version", "")
            info["date"] = elem.attrib.get("date", "")
            continue
        if event != "end" or elem.tag != "Disorder":
            continue
        code = (elem.findtext("OrphaCode") or "").strip()
        if not code:
            elem.clear()
            continue
        for assoc in elem.iter("DisorderGeneAssociation"):
            symbol = (assoc.findtext("Gene/Symbol") or "").strip()
            kind = (assoc.findtext("DisorderGeneAssociationType/Name") or "").strip()
            status = (assoc.findtext("DisorderGeneAssociationStatus/Name") or "").strip()
            if not symbol:
                continue
            if kind in ORPHA_USE:
                use, flags = 1, ""
            elif kind in ORPHA_EXCLUDE:
                use, flags = 0, "excluded_type"
            else:
                use, flags = 0, "review"
            if status and status.lower() != "assessed":
                flags = ",".join(x for x in (flags, "not_assessed") if x)
            gene = canonical(symbol)
            key = (f"ORPHA:{code}", gene)
            prev = out.get(key)
            rec = {"disease_id": key[0], "gene_symbol": gene, "source": "Orphanet", "association_type": kind,
                   "use": use, "flags": flags, "status": status, "source_gene": symbol}
            if prev is None or (rec["use"] and not prev["use"]):
                out[key] = rec
        elem.clear()
    return [out[k] for k in sorted(out)], info


def orphanet_links_from_g2d(g2d: dict[str, dict[str, str]], canonical=_canonical,
                            skip: set[str] | None = None) -> list[dict]:
    """Orphanet links taken from HPO genes_to_disease.txt.

    Used for every ORPHA disease when product6 is unavailable, and (``skip`` =
    the diseases product6 lists) for ORPHA diseases newer than the local
    product6 — HPO reads a more recent product6 than the copy the server can
    download.  HPO drops the association type, so every row is used and
    flagged ``via_genes_to_disease`` (plus ``not_in_product6`` in the second case)."""
    flags = "via_genes_to_disease" if skip is None else "via_genes_to_disease,not_in_product6"
    skip = skip or set()
    out = {}
    for disease, genes in g2d.items():
        if not disease.startswith("ORPHA:") or disease in skip:
            continue
        for symbol in genes:
            gene = canonical(symbol)
            out[(disease, gene)] = {"disease_id": disease, "gene_symbol": gene, "source": "Orphanet",
                                    "association_type": "unknown (HPO genes_to_disease)", "use": 1,
                                    "flags": flags, "status": "", "source_gene": symbol}
    return [out[k] for k in sorted(out)]


def cross_check(links: list[dict], g2d: dict[str, dict[str, str]], universe: set[str],
                canonical=_canonical) -> tuple[list[tuple], dict]:
    """Disagreements between the used links and HPO genes_to_disease, per disease prefix."""
    used = {(r["disease_id"], r["gene_symbol"]) for r in links if r["use"]}
    ref = {(d, canonical(g)) for d, genes in g2d.items() for g in genes}
    rows = []
    stats = {}
    for prefix in ("OMIM", "ORPHA"):
        a = {p for p in used if p[0].startswith(prefix + ":") and p[0] in universe}
        b = {p for p in ref if p[0].startswith(prefix + ":") and p[0] in universe}
        stats[prefix] = {"kb": len(a), "genes_to_disease": len(b), "both": len(a & b),
                         "kb_only": len(a - b), "genes_to_disease_only": len(b - a)}
        rows += [(d, g, 1, 0) for d, g in sorted(a - b)] + [(d, g, 0, 1) for d, g in sorted(b - a)]
    return sorted(rows), stats
