"""PGx TSV + PharmCAT report adapter."""
from __future__ import annotations

import csv
import html
import json
import re
from collections import Counter
from pathlib import Path

EVIDENCE_ORDER = {"Strong": 0, "Moderate": 1, "Optional": 2, "Unspecified": 3, ".": 4, "": 5}


def _clean(value) -> str:
    s = "" if value is None else str(value).strip()
    return "" if not s or s.upper() in {"NA", "N/A", "."} else s


def _is_no_action(text: str) -> bool:
    t = (text or "").strip().lower()
    return (
        "no action is required" in t
        or "no action is needed" in t
        or "no reason to avoid" in t
        or t == "no recommendation"
    )


def _is_actionable(gene: str, phenotype: str, rows: list[dict]) -> bool:
    ph = (phenotype or "").strip().lower()
    sym = (gene or "").strip().upper()
    if sym.startswith("HLA"):
        return "positive" in ph
    if sym == "MT-RNR1":
        return any((r.get("MTRN1_RISK") or "").upper() == "HIGH" for r in rows)
    if "uncertain susceptibility" in ph:
        return False
    if ph in {"normal", "normal metabolizer", "normal function"}:
        return False
    if "normal metabolizer" in ph or "normal function" in ph:
        return False
    if "unknown" in ph or "indeterminate" in ph:
        return False
    if any(
        (r.get("EVIDENCE_STRENGTH") in {"Strong", "Moderate"}
         or r.get("CPIC_LEVEL") in {"Strong", "Moderate"}
         or r.get("DPWG_LEVEL") in {"Strong", "Moderate"})
        and not _is_no_action(r.get("RECOMMENDATION") or "")
        for r in rows
    ):
        return True
    return any(word in ph for word in ("poor", "intermediate", "rapid", "ultrarapid", "increased risk", "decreased"))


def _evidence_rank(row: dict) -> int:
    vals = [
        row.get("EVIDENCE_STRENGTH") or row.get("evidence") or "",
        row.get("CPIC_LEVEL") or row.get("cpic_level") or "",
        row.get("DPWG_LEVEL") or row.get("dpwg_level") or "",
    ]
    return min(EVIDENCE_ORDER.get(v, 5) for v in vals)


def _plain_text(value) -> str:
    text = html.unescape(str(value or ""))
    text = re.sub(r"<li[^>]*>", "- ", text, flags=re.I)
    text = re.sub(r"</?(?:ul|ol|p|div|br)[^>]*>", "\n", text, flags=re.I)
    text = re.sub(r"<[^>]+>", "", text)
    return re.sub(r"\s+", " ", text).strip()


def _annotation_genes(annotation: dict) -> list[str]:
    genes: list[str] = []
    for genotype_group in annotation.get("genotypes") or []:
        for diplotype in genotype_group.get("diplotypes") or []:
            gene = _clean(diplotype.get("gene"))
            if gene:
                genes.append(gene)
    for item in annotation.get("lookupKey") or []:
        if isinstance(item, dict):
            genes.extend(_clean(gene) for gene in item if _clean(gene))
    for implication in annotation.get("implications") or []:
        match = re.match(r"\s*([A-Za-z0-9-]+)\s*:", str(implication))
        if match:
            genes.append(match.group(1))
    return list(dict.fromkeys(genes))


def _fda_category(section: str, guideline_name: str) -> str:
    if section != "FDA PGx Association":
        return ""
    text = guideline_name.lower()
    if "therapeutic management" in text:
        return "therapeutic_management"
    if "potential impact" in text:
        return "potential_impact"
    if "pharmacokinetic properties only" in text:
        return "pharmacokinetic_only"
    return "unspecified"


def _compact_drug_annotations(data: dict) -> list[dict]:
    out: list[dict] = []
    for section, drugs in (data.get("drugs") or {}).items():
        if not isinstance(drugs, dict):
            continue
        for drug_key, payload in drugs.items():
            if not isinstance(payload, dict):
                continue
            drug = _clean(payload.get("name")) or _clean(drug_key)
            for guideline in payload.get("guidelines") or []:
                guideline_name = _clean(guideline.get("name"))
                source = _clean(guideline.get("source")) or _clean(payload.get("source"))
                for annotation in guideline.get("annotations") or []:
                    recommendation = _plain_text(annotation.get("drugRecommendation"))
                    if not recommendation:
                        continue
                    out.append({
                        "section": _clean(section),
                        "source": source,
                        "drug": drug,
                        "guideline": guideline_name,
                        "url": _clean(guideline.get("url")),
                        "classification": _clean(annotation.get("classification")) or "Unspecified",
                        "recommendation": recommendation,
                        "implications": [
                            _plain_text(value)
                            for value in annotation.get("implications") or []
                            if _plain_text(value)
                        ],
                        "genes": _annotation_genes(annotation),
                        "fda_category": _fda_category(section, guideline_name),
                        "dosing_information": bool(annotation.get("dosingInformation")),
                        "alternate_drug_available": bool(annotation.get("alternateDrugAvailable")),
                        "other_prescribing_guidance": bool(annotation.get("otherPrescribingGuidance")),
                    })
    return out


def _compact_report_json(path: Path) -> dict:
    if not path.exists():
        return {}
    try:
        with path.open("r", encoding="utf-8") as fh:
            data = json.load(fh)
    except (OSError, json.JSONDecodeError):
        return {}

    genes = {}
    for gene, payload in (data.get("genes") or {}).items():
        source = (payload.get("sourceDiplotypes") or [{}])[0] or {}
        allele1 = source.get("allele1") or {}
        allele2 = source.get("allele2") or {}
        variants = []
        for variant in payload.get("variants") or []:
            call = _clean(variant.get("call"))
            ref = _clean(variant.get("referenceAllele"))
            if ref and call == f"{ref}/{ref}":
                continue
            variants.append({
                "rsid": _clean(variant.get("dbSnpId")),
                "chr": _clean(variant.get("chromosome")),
                "pos": variant.get("position"),
                "call": call,
                "alleles": ", ".join(variant.get("alleles") or []),
            })
        genes[gene] = {
            "label": _clean(source.get("label")),
            "phenotypes": [
                _clean(value)
                for value in source.get("phenotypes") or []
                if _clean(value)
            ],
            "allele1_name": _clean(allele1.get("name")),
            "allele1_function": _clean(allele1.get("function")),
            "allele2_name": _clean(allele2.get("name")),
            "allele2_function": _clean(allele2.get("function")),
            "uncalled": payload.get("uncalledHaplotypes") or [],
            "messages": payload.get("messages") or [],
            "variants": variants[:80],
        }
    return {
        "pharmcat_version": _clean(data.get("pharmcatVersion")),
        "data_version": _clean(data.get("dataVersion")),
        "timestamp": _clean(data.get("timestamp")),
        "messages": data.get("messages") or [],
        "genes": genes,
        "guideline_annotations": _compact_drug_annotations(data),
    }


def load_pgx(pgx_tsv: Path, pharmcat_json: Path | None = None) -> dict:
    report = _compact_report_json(pharmcat_json) if pharmcat_json else {}
    result = {
        "genes": {},
        "gene_order": [],
        "actionable": [],
        "routine": [],
        "additional_genes": [],
        "summary": {
            "called_genes": 0,
            "actionable_genes": 0,
            "recommendations": 0,
            "strong_recommendations": 0,
        },
        "pharmcat_version": report.get("pharmcat_version", ""),
        "data_version": report.get("data_version", ""),
        "timestamp": report.get("timestamp", ""),
        "pharmcat_available": bool(report),
        "messages": report.get("messages", []),
        "guideline_annotations": report.get("guideline_annotations", []),
    }
    rows_by_gene: dict[str, list[dict]] = {}
    if pgx_tsv.exists():
        with pgx_tsv.open("r", encoding="utf-8-sig", newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for raw in reader:
                row = {key: _clean(value) for key, value in raw.items()}
                gene = row.get("GENE")
                if not gene:
                    continue
                rows_by_gene.setdefault(gene, []).append(row)

    details_by_gene = report.get("genes") or {}
    for gene, rows in rows_by_gene.items():
        first = rows[0]
        drug_map: dict[str, list[dict]] = {}
        evidence_counts = Counter()
        for row in rows:
            drug = row.get("DRUG") or "general"
            evidence = row.get("EVIDENCE_STRENGTH") or row.get("CPIC_LEVEL") or row.get("DPWG_LEVEL") or "Unspecified"
            evidence_counts[evidence] += 1
            drug_map.setdefault(drug, []).append({
                "source": row.get("GUIDELINE_SOURCE"),
                "recommendation": row.get("RECOMMENDATION"),
                "implication": row.get("IMPLICATION"),
                "cpic_level": row.get("CPIC_LEVEL") or ".",
                "dpwg_level": row.get("DPWG_LEVEL") or ".",
                "evidence": row.get("EVIDENCE_STRENGTH") or "Unspecified",
            })
        drugs = []
        for drug, recs in drug_map.items():
            recs.sort(key=_evidence_rank)
            drugs.append({
                "drug": drug,
                "recommendations": recs,
                "best_rank": min((_evidence_rank(r) for r in recs), default=9),
            })
        drugs.sort(key=lambda d: (d["best_rank"], d["drug"]))

        actionable = _is_actionable(gene, first.get("PHENOTYPE") or "", rows)
        gene_payload = {
            "gene": gene,
            "pipeline": first.get("PIPELINE"),
            "diplotype": first.get("DIPLOTYPE"),
            "activity_score": first.get("ACTIVITY_SCORE"),
            "phenotype": first.get("PHENOTYPE"),
            "outside_caller": first.get("OUTSIDE_CALLER"),
            "mtrn1_risk": first.get("MTRN1_RISK"),
            "notes": first.get("NOTES"),
            "evidence_counts": dict(evidence_counts),
            "drugs": drugs,
            "actionable": actionable,
            "details": details_by_gene.get(gene, {}),
        }
        result["genes"][gene] = gene_payload
        result["gene_order"].append(gene)

    def gene_sort_key(gene: str):
        g = result["genes"][gene]
        best = min((d["best_rank"] for d in g.get("drugs") or []), default=9)
        return (0 if g.get("actionable") else 1, best, gene)

    result["gene_order"].sort(key=gene_sort_key)
    result["actionable"] = [g for g in result["gene_order"] if result["genes"][g].get("actionable")]
    result["routine"] = [g for g in result["gene_order"] if not result["genes"][g].get("actionable")]
    for gene in sorted(set(details_by_gene) - set(result["genes"])):
        details = details_by_gene.get(gene) or {}
        result["genes"][gene] = {
            "gene": gene,
            "pipeline": "",
            "diplotype": details.get("label") or "",
            "activity_score": "",
            "phenotype": "",
            "outside_caller": "PharmCAT",
            "mtrn1_risk": "",
            "notes": "",
            "evidence_counts": {},
            "drugs": [],
            "actionable": False,
            "details": details,
            "additional": True,
        }
        result["additional_genes"].append(gene)
    result["summary"]["called_genes"] = len(result["gene_order"])
    result["summary"]["actionable_genes"] = len(result["actionable"])
    result["summary"]["recommendations"] = sum(len(rows) for rows in rows_by_gene.values())
    result["summary"]["strong_recommendations"] = sum(
        1
        for rows in rows_by_gene.values()
        for row in rows
        if "Strong" in {row.get("EVIDENCE_STRENGTH"), row.get("CPIC_LEVEL"), row.get("DPWG_LEVEL")}
    )
    return result
