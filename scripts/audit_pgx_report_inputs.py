#!/usr/bin/env python3
"""Audit raw PharmCAT JSON/PGx TSV files against the health-report table.

The script is intentionally standalone (Python standard library only) so it can
be copied to, or executed directly on, the analysis server.  It never modifies
sample outputs.  When --output-dir is provided it writes audit files only to
that directory.
"""
from __future__ import annotations

import argparse
import csv
import json
import re
import sys
import unicodedata
from collections import Counter, defaultdict
from pathlib import Path


REPORT_GENES = (
    "ABCG2", "CACNA1S", "CFTR", "CYP2B6", "CYP2C19", "CYP2C9",
    "CYP2D6", "CYP3A5", "CYP4F2", "DPYD", "G6PD", "HLA-A", "HLA-B",
    "MT-RNR1", "NAT2", "NUDT15", "RYR1", "SLCO1B1", "TPMT",
    "UGT1A1", "VKORC1",
)
HLA_TARGETS = (
    ("HLA-A", "*31:01"),
    ("HLA-B", "*15:02"),
    ("HLA-B", "*57:01"),
    ("HLA-B", "*58:01"),
)
MARKERS = {
    "ABCG2": {
        "test": "ABCG2 rs2231142 (c.421C>A)",
        "reference": "C",
        "variant": "A",
        "base_map": {"G": "C", "T": "A", "C": "C", "A": "A"},
    },
    "VKORC1": {
        "test": "VKORC1 rs9923231 (c.-1639G>A)",
        "reference": "G",
        "variant": "A",
        "base_map": {"C": "G", "T": "A", "G": "G", "A": "A"},
    },
}
NO_RESULT_VALUES = {
    "", "-", "—", ".", "n/a", "na", "unknown", "no result",
    "not available",
}
CASE_FIELDS = (
    "sample", "row_kind", "gene", "test", "source_rule", "json_present",
    "source_diplotype_count", "json_label", "json_allele1", "json_allele2",
    "json_phenotypes", "tsv_present", "tsv_diplotypes", "tsv_phenotypes",
    "allele1", "allele2", "result_span", "compatibility", "issue_codes",
)
OBSERVED_FIELDS = (
    "gene", "test", "source_rule", "json_label", "json_allele1",
    "json_allele2", "json_phenotypes", "tsv_diplotypes", "tsv_phenotypes",
    "allele1", "allele2", "result_span", "count", "sample_examples",
)
ISSUE_FIELDS = ("sample", "gene", "severity", "code", "detail")
SOURCE_FIELDS = (
    "sample", "gene", "source_index", "label", "allele1", "allele2",
    "phenotypes", "combination", "inferred", "match_score",
)
VARIANT_FIELDS = (
    "sample", "gene", "rsid", "chromosome", "position", "call",
    "reference", "alleles", "phased", "phase_set",
)


def _clean(value) -> str:
    return "" if value is None else str(value).strip()


def _is_no_result(value) -> bool:
    return _clean(value).lower() in NO_RESULT_VALUES


def _normalize_allele(value) -> str:
    text = _clean(value)
    if _is_no_result(text):
        return "No Result"
    if text.lower() == "reference":
        return "Reference"
    if text.lower() == "variant":
        return "Variant"
    text = re.sub(r"\breference\b", "Reference", text, flags=re.I)
    text = re.sub(r"\bvariant\b", "Variant", text, flags=re.I)
    if re.match(r"^[cmng]\.[^\s()]+\(", text, flags=re.I):
        text = re.sub(r"\(", " (", text, count=1)
    return text


def _split_label(label: str) -> tuple[str, str]:
    text = _clean(label)
    if "/" not in text:
        return text, ""
    allele1, allele2 = text.split("/", 1)
    return allele1.strip(), allele2.strip()


def _source_details(payload: dict) -> tuple[list[dict], dict]:
    sources = payload.get("sourceDiplotypes") or []
    sources = [source for source in sources if isinstance(source, dict)]
    return sources, (sources[0] if sources else {})


def _source_alleles(
    gene: str,
    payload: dict,
    tsv_diplotypes: list[str],
) -> tuple[str, str, str]:
    _sources, source = _source_details(payload)
    label1, label2 = _split_label(source.get("label") or "")
    raw1 = _clean((source.get("allele1") or {}).get("name"))
    raw2 = _clean((source.get("allele2") or {}).get("name"))
    if (
        gene in {"HLA-A", "HLA-B"}
        and not raw1
        and not raw2
        and re.search(r"\b(?:positive|negative)\b", source.get("label") or "", re.I)
    ):
        label1, label2 = "", ""
    if gene == "CFTR" and label1 and label2:
        allele1, allele2 = label1, label2
    else:
        allele1 = raw1 if not _is_no_result(raw1) else label1
        allele2 = raw2 if not _is_no_result(raw2) else label2
    source_rule = "JSON"
    if gene == "MT-RNR1" and _is_no_result(allele1) and _is_no_result(allele2):
        if len(tsv_diplotypes) == 1 and not _is_no_result(tsv_diplotypes[0]):
            allele1, allele2 = tsv_diplotypes[0], ""
            source_rule = "TSV fallback"
        else:
            source_rule = "No usable result"
    elif not payload:
        source_rule = "No JSON gene"
    elif _is_no_result(allele1) and _is_no_result(allele2):
        source_rule = "JSON no result"
    return _normalize_allele(allele1), _normalize_allele(allele2), source_rule


def _hla_display(gene: str, allele: str) -> str:
    value = _normalize_allele(allele)
    if value == "No Result":
        return value
    value = re.sub(r"^HLA-", "", value, flags=re.I)
    if value.startswith("*"):
        return f"{gene.rsplit('-', 1)[-1]}{value}"
    return value


def _hla_matches(gene: str, allele: str, target: str) -> bool:
    match = re.search(
        r"\*([0-9]+(?::[0-9A-Z]+)+)",
        _hla_display(gene, allele),
        flags=re.I,
    )
    if not match:
        return False
    called = match.group(1).upper()
    expected = target.lstrip("*").upper()
    return called == expected or called.startswith(f"{expected}:")


def _hla_explicit_status(
    payload: dict,
    target: str,
    tsv_phenotypes: list[str],
) -> str:
    _sources, source = _source_details(payload)
    candidates = [
        *(source.get("phenotypes") or []),
        source.get("label") or "",
        *tsv_phenotypes,
    ]
    token = re.escape(target.lstrip("*")).replace(":", r"\s*:\s*")
    pattern = re.compile(
        rf"(?<!\d){token}(?!\d)\s*(?:[:：=]\s*)?(positive|negative)\b",
        re.I,
    )
    for candidate in candidates:
        match = pattern.search(_clean(candidate))
        if match:
            return match.group(1).title()
    return ""


def _hla_statuses(
    gene: str,
    target: str,
    called_alleles: tuple[str, str],
    payload: dict,
    tsv_phenotypes: list[str],
) -> tuple[str, str]:
    statuses = [
        "" if _is_no_result(called) else (
            "Positive" if _hla_matches(gene, called, target) else "Negative"
        )
        for called in called_alleles
    ]
    explicit = _hla_explicit_status(payload, target, tsv_phenotypes)
    missing = [index for index, status in enumerate(statuses) if not status]
    if explicit == "Negative":
        for index in missing:
            statuses[index] = "Negative"
    elif explicit == "Positive" and missing and "Positive" not in statuses:
        statuses[missing[0]] = "Positive"
    return tuple(status or "No Result" for status in statuses)


def _marker_allele(gene: str, allele: str) -> tuple[str, bool]:
    value = _normalize_allele(allele)
    if value == "No Result":
        return value, True
    marker = MARKERS[gene]
    lower = value.lower()
    status = (
        "Reference" if "reference" in lower
        else "Variant" if "variant" in lower
        else ""
    )
    match = re.search(r"\(([ACGT])\)", value, flags=re.I)
    if match:
        base = match.group(1)
    else:
        bases = re.findall(r"(?<![A-Z])([ACGT])(?![A-Z])", value, flags=re.I)
        base = bases[-1] if bases else ""
    transcript_base = marker["base_map"].get(base.upper(), "")
    if not status and transcript_base:
        if transcript_base == marker["reference"]:
            status = "Reference"
        elif transcript_base == marker["variant"]:
            status = "Variant"
    if transcript_base and status:
        return f"{transcript_base} ({status})", True
    return status or transcript_base or value, False


def _display_width(value: str) -> int:
    return sum(
        2 if unicodedata.east_asian_width(char) in {"W", "F"} else 1
        for char in _clean(value)
    )


def _read_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        data = json.load(handle)
    if not isinstance(data, dict):
        raise ValueError("top-level JSON is not an object")
    return data


def _read_tsv(path: Path) -> dict[str, list[dict[str, str]]]:
    rows_by_gene: dict[str, list[dict[str, str]]] = defaultdict(list)
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or "GENE" not in reader.fieldnames:
            raise ValueError("TSV has no GENE column")
        for row in reader:
            gene = _clean(row.get("GENE"))
            if gene:
                rows_by_gene[gene].append({key: _clean(value) for key, value in row.items()})
    return rows_by_gene


def _unique(rows: list[dict[str, str]], key: str) -> list[str]:
    return sorted({row.get(key, "") for row in rows if not _is_no_result(row.get(key, ""))})


def _is_reference_call(call: str, reference: str) -> bool:
    ref = _clean(reference)
    values = [value for value in re.split(r"[/|]", _clean(call)) if value]
    return bool(ref and values and all(value == ref for value in values))


def _variant_zygosity(variant: dict) -> str:
    reference = _clean(variant.get("referenceAllele"))
    values = [
        value for value in re.split(r"[/|]", _clean(variant.get("call")))
        if value
    ]
    if not values:
        return ""
    unique = set(values)
    if len(unique) == 1:
        return "ref" if reference and unique == {reference} else "hom"
    if reference and reference in unique:
        return "het"
    return "non-reference"


def _dpyd_unphased_variants(payload: dict) -> list[str]:
    sources, _source = _source_details(payload)
    if len(sources) <= 1 or payload.get("effectivelyPhased"):
        return []

    def key(value: str) -> str:
        return re.sub(r"\s+", "", _clean(value)).lower()

    variants = [value for value in payload.get("variants") or [] if isinstance(value, dict)]
    observed: list[str] = []
    seen: set[str] = set()
    for source in sources:
        for raw_name in (
            (source.get("allele1") or {}).get("name"),
            (source.get("allele2") or {}).get("name"),
        ):
            name = _normalize_allele(raw_name)
            name_key = key(name)
            if _is_no_result(name) or name == "Reference" or name_key in seen:
                continue
            seen.add(name_key)
            zygosity = ""
            for variant in variants:
                if any(key(candidate) == name_key for candidate in variant.get("alleles") or []):
                    zygosity = _variant_zygosity(variant)
                    break
            observed.append(f"{name} ({zygosity})" if zygosity else name)
    return observed


def _add_issue(
    issues: list[dict[str, str]],
    sample: str,
    gene: str,
    severity: str,
    code: str,
    detail: str,
) -> None:
    issues.append({
        "sample": sample,
        "gene": gene,
        "severity": severity,
        "code": code,
        "detail": detail,
    })


def _case_row(
    sample: str,
    gene: str,
    json_payload: dict,
    tsv_rows: list[dict[str, str]],
    issues: list[dict[str, str]],
) -> tuple[dict[str, str], tuple[str, str], list[str]]:
    sources, source = _source_details(json_payload)
    json_label = _clean(source.get("label"))
    json_a1 = _clean((source.get("allele1") or {}).get("name"))
    json_a2 = _clean((source.get("allele2") or {}).get("name"))
    json_phenotypes = sorted({_clean(value) for value in source.get("phenotypes") or [] if _clean(value)})
    tsv_diplotypes = _unique(tsv_rows, "DIPLOTYPE")
    tsv_phenotypes = _unique(tsv_rows, "PHENOTYPE")
    allele1, allele2, source_rule = _source_alleles(gene, json_payload, tsv_diplotypes)
    test = MARKERS.get(gene, {}).get("test", gene)
    result_span = ""
    issue_codes: list[str] = []

    def issue(severity: str, code: str, detail: str) -> None:
        issue_codes.append(code)
        _add_issue(issues, sample, gene, severity, code, detail)

    if not json_payload:
        issue("WARN", "JSON_GENE_MISSING", "fixed report gene is absent from JSON")
    elif not sources:
        issue("WARN", "SOURCE_DIPLOTYPE_MISSING", "JSON gene has no sourceDiplotypes")
    elif len(sources) > 1 and gene == "DPYD" and _dpyd_unphased_variants(json_payload):
        issue(
            "INFO", "DPYD_MULTIPLE_UNPHASED_VARIANTS",
            f"health report preserves all {len(sources)} sourceDiplotypes in a spanning result",
        )
    elif len(sources) > 1:
        issue(
            "WARN", "MULTIPLE_SOURCE_DIPLOTYPES",
            f"unexpected {len(sources)} sourceDiplotypes require review",
        )
    if len(tsv_diplotypes) > 1:
        issue("WARN", "MULTIPLE_TSV_DIPLOTYPES", " | ".join(tsv_diplotypes))
    if json_label and len(tsv_diplotypes) == 1 and json_label != tsv_diplotypes[0]:
        severity = "INFO" if gene == "MT-RNR1" else "WARN"
        issue(
            severity, "JSON_TSV_DIPLOTYPE_DIFFER",
            f"JSON={json_label!r}; TSV={tsv_diplotypes[0]!r}",
        )

    marker_ok = True
    if gene in MARKERS:
        allele1, ok1 = _marker_allele(gene, allele1)
        allele2, ok2 = _marker_allele(gene, allele2)
        marker_ok = ok1 and ok2
        if not marker_ok:
            issue(
                "WARN", "MARKER_FORMAT_UNRECOGNIZED",
                f"raw alleles={json_a1!r}, {json_a2!r}",
            )
    elif gene in {"HLA-A", "HLA-B"}:
        allele1 = _hla_display(gene, allele1)
        allele2 = _hla_display(gene, allele2)
        for position, allele in (("Allele 1", allele1), ("Allele 2", allele2)):
            if allele != "No Result" and not re.search(r"\*[0-9]+(?::[0-9A-Z]+)+", allele, re.I):
                issue("WARN", "HLA_TYPING_UNRECOGNIZED", f"{position}={allele!r}")
    elif gene == "MT-RNR1":
        allele2 = "N/A"
    elif gene == "G6PD" and allele1 != "No Result" and allele2 == "No Result":
        allele2 = "N/A"

    if gene == "DPYD":
        unphased_variants = _dpyd_unphased_variants(json_payload)
        if unphased_variants:
            test = "DPYD（相位未定）"
            allele1 = ""
            allele2 = ""
            result_span = "檢出變異：" + "；".join(unphased_variants)

    if not result_span and allele1 == "No Result" and allele2 in {"No Result", "N/A"}:
        issue("WARN", "NO_EFFECTIVE_RESULT", f"source rule={source_rule}")
    elif not result_span and gene not in {"MT-RNR1", "G6PD"} and allele2 == "No Result":
        issue("WARN", "SECOND_ALLELE_MISSING", f"Allele 1={allele1!r}")
    if _display_width(allele1) > 22:
        issue("INFO", "ALLELE1_WRAPS", f"display width={_display_width(allele1)}")
    if _display_width(allele2) > 22:
        issue("INFO", "ALLELE2_WRAPS", f"display width={_display_width(allele2)}")
    if _display_width(result_span) > 45:
        issue("INFO", "RESULT_SPAN_WRAPS", f"display width={_display_width(result_span)}")

    compatibility = "OK"
    if any(code in {
        "JSON_GENE_MISSING", "SOURCE_DIPLOTYPE_MISSING", "MULTIPLE_SOURCE_DIPLOTYPES",
        "MULTIPLE_TSV_DIPLOTYPES", "MARKER_FORMAT_UNRECOGNIZED",
        "HLA_TYPING_UNRECOGNIZED", "NO_EFFECTIVE_RESULT", "SECOND_ALLELE_MISSING",
    } for code in issue_codes):
        compatibility = "REVIEW"
    elif any(code.endswith("WRAPS") for code in issue_codes):
        compatibility = "WRAPS"

    row = {
        "sample": sample,
        "row_kind": "gene",
        "gene": gene,
        "test": test,
        "source_rule": source_rule,
        "json_present": "yes" if json_payload else "no",
        "source_diplotype_count": str(len(sources)),
        "json_label": json_label,
        "json_allele1": json_a1,
        "json_allele2": json_a2,
        "json_phenotypes": " | ".join(json_phenotypes),
        "tsv_present": "yes" if tsv_rows else "no",
        "tsv_diplotypes": " | ".join(tsv_diplotypes),
        "tsv_phenotypes": " | ".join(tsv_phenotypes),
        "allele1": allele1,
        "allele2": allele2,
        "result_span": result_span,
        "compatibility": compatibility,
        "issue_codes": " | ".join(dict.fromkeys(issue_codes)),
    }
    return row, (allele1, allele2), tsv_phenotypes


def _hla_rows(
    base_row: dict[str, str],
    called_alleles: tuple[str, str],
    payload: dict,
    tsv_phenotypes: list[str],
) -> list[dict[str, str]]:
    rows = []
    gene = base_row["gene"]
    for parent, target in HLA_TARGETS:
        if parent != gene:
            continue
        allele1, allele2 = _hla_statuses(
            gene, target, called_alleles, payload, tsv_phenotypes,
        )
        rows.append({
            **base_row,
            "row_kind": "hla_screen",
            "gene": f"{gene}{target}",
            "test": f"{gene}{target}",
            "source_rule": "JSON HLA typing/status",
            "allele1": allele1,
            "allele2": allele2,
            "result_span": "",
            "compatibility": (
                "REVIEW" if "No Result" in {allele1, allele2} else "OK"
            ),
            "issue_codes": (
                "HLA_SCREEN_PARTIAL" if "No Result" in {allele1, allele2} else ""
            ),
        })
    return rows


def audit(root: Path, prefix: str) -> dict:
    sample_dirs = sorted(
        path for path in root.glob(f"{prefix}*")
        if path.is_dir()
    )
    cases: list[dict[str, str]] = []
    issues: list[dict[str, str]] = []
    sample_files: list[dict[str, str]] = []
    source_diplotypes: list[dict[str, str]] = []
    nonreference_variants: list[dict[str, str]] = []

    for sample_dir in sample_dirs:
        sample = sample_dir.name
        pgx_dir = sample_dir / "07_pgx"
        json_paths = sorted(pgx_dir.glob("*.pharmcat.report.json")) if pgx_dir.is_dir() else []
        tsv_paths = sorted(pgx_dir.glob("*.pgx.tsv")) if pgx_dir.is_dir() else []
        json_path = json_paths[0] if len(json_paths) == 1 else None
        tsv_path = tsv_paths[0] if len(tsv_paths) == 1 else None
        sample_files.append({
            "sample": sample,
            "json_file": str(json_path or ""),
            "tsv_file": str(tsv_path or ""),
            "json_file_count": str(len(json_paths)),
            "tsv_file_count": str(len(tsv_paths)),
        })
        if len(json_paths) != 1:
            _add_issue(
                issues, sample, "", "ERROR", "JSON_FILE_COUNT",
                f"expected 1, found {len(json_paths)} in {pgx_dir}",
            )
        if len(tsv_paths) != 1:
            _add_issue(
                issues, sample, "", "ERROR", "TSV_FILE_COUNT",
                f"expected 1, found {len(tsv_paths)} in {pgx_dir}",
            )
        try:
            json_data = _read_json(json_path) if json_path else {}
        except (OSError, ValueError, json.JSONDecodeError) as exc:
            json_data = {}
            _add_issue(issues, sample, "", "ERROR", "JSON_READ_ERROR", str(exc))
        try:
            rows_by_gene = _read_tsv(tsv_path) if tsv_path else {}
        except (OSError, ValueError, csv.Error) as exc:
            rows_by_gene = {}
            _add_issue(issues, sample, "", "ERROR", "TSV_READ_ERROR", str(exc))

        genes = json_data.get("genes") or {}
        if not isinstance(genes, dict):
            genes = {}
            _add_issue(issues, sample, "", "ERROR", "JSON_GENES_INVALID", "genes is not an object")
        for gene in REPORT_GENES:
            payload = genes.get(gene) if isinstance(genes.get(gene), dict) else {}
            sources, _source = _source_details(payload)
            for source_index, source in enumerate(sources, start=1):
                source_diplotypes.append({
                    "sample": sample,
                    "gene": gene,
                    "source_index": str(source_index),
                    "label": _clean(source.get("label")),
                    "allele1": _clean((source.get("allele1") or {}).get("name")),
                    "allele2": _clean((source.get("allele2") or {}).get("name")),
                    "phenotypes": " | ".join(
                        _clean(value)
                        for value in source.get("phenotypes") or []
                        if _clean(value)
                    ),
                    "combination": str(bool(source.get("combination"))).lower(),
                    "inferred": str(bool(source.get("inferred"))).lower(),
                    "match_score": _clean(source.get("matchScore")),
                })
            for variant in payload.get("variants") or []:
                if not isinstance(variant, dict):
                    continue
                call = _clean(variant.get("call"))
                reference = _clean(variant.get("referenceAllele"))
                if _is_reference_call(call, reference):
                    continue
                nonreference_variants.append({
                    "sample": sample,
                    "gene": gene,
                    "rsid": _clean(variant.get("dbSnpId")),
                    "chromosome": _clean(variant.get("chromosome")),
                    "position": _clean(variant.get("position")),
                    "call": call,
                    "reference": reference,
                    "alleles": " | ".join(
                        _clean(value)
                        for value in variant.get("alleles") or []
                        if _clean(value)
                    ),
                    "phased": str(bool(variant.get("phased"))).lower(),
                    "phase_set": _clean(variant.get("phaseSet")),
                })
            base_row, called_alleles, tsv_phenotypes = _case_row(
                sample,
                gene,
                payload,
                rows_by_gene.get(gene, []),
                issues,
            )
            cases.append(base_row)
            cases.extend(_hla_rows(
                base_row, called_alleles, payload, tsv_phenotypes,
            ))

    observed_groups: dict[tuple[str, ...], list[str]] = defaultdict(list)
    for row in cases:
        key = tuple(row.get(field, "") for field in OBSERVED_FIELDS[:-2])
        observed_groups[key].append(row["sample"])
    observed = []
    for key, samples in observed_groups.items():
        row = dict(zip(OBSERVED_FIELDS[:-2], key))
        row["count"] = str(len(samples))
        row["sample_examples"] = ", ".join(samples[:10])
        observed.append(row)
    observed.sort(key=lambda row: (
        REPORT_GENES.index(row["gene"]) if row["gene"] in REPORT_GENES else 999,
        row["gene"], row["allele1"], row["allele2"], row["json_label"],
    ))

    summary = []
    rows_by_report_item: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in cases:
        rows_by_report_item[row["gene"]].append(row)
    report_order = []
    for report_gene in REPORT_GENES:
        report_order.append(report_gene)
        report_order.extend(
            f"{parent}{target}"
            for parent, target in HLA_TARGETS
            if parent == report_gene
        )
    for gene in report_order:
        gene_rows = rows_by_report_item.get(gene, [])
        summary.append({
            "gene": gene,
            "samples": str(len(gene_rows)),
            "json_present": str(sum(row["json_present"] == "yes" for row in gene_rows)),
            "tsv_present": str(sum(row["tsv_present"] == "yes" for row in gene_rows)),
            "json_source": str(sum(row["source_rule"].startswith("JSON") for row in gene_rows)),
            "tsv_fallback": str(sum(row["source_rule"] == "TSV fallback" for row in gene_rows)),
            "ok": str(sum(row["compatibility"] == "OK" for row in gene_rows)),
            "wraps": str(sum(row["compatibility"] == "WRAPS" for row in gene_rows)),
            "review": str(sum(row["compatibility"] == "REVIEW" for row in gene_rows)),
            "distinct_outputs": str(len({
                (row["allele1"], row["allele2"], row["result_span"])
                for row in gene_rows
            })),
        })

    severity_counts = Counter(issue["severity"] for issue in issues)
    code_counts = Counter(issue["code"] for issue in issues)
    return {
        "root": str(root),
        "prefix": prefix,
        "sample_count": len(sample_dirs),
        "complete_file_pairs": sum(
            row["json_file_count"] == "1" and row["tsv_file_count"] == "1"
            for row in sample_files
        ),
        "summary": summary,
        "observed": observed,
        "cases": cases,
        "issues": issues,
        "sample_files": sample_files,
        "source_diplotypes": source_diplotypes,
        "nonreference_variants": nonreference_variants,
        "severity_counts": dict(sorted(severity_counts.items())),
        "issue_code_counts": dict(sorted(code_counts.items())),
    }


def _write_tsv(path: Path, rows: list[dict], fields: tuple[str, ...]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_outputs(report: dict, output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    summary_fields = (
        "gene", "samples", "json_present", "tsv_present", "json_source",
        "tsv_fallback", "ok", "wraps", "review", "distinct_outputs",
    )
    sample_fields = (
        "sample", "json_file", "tsv_file", "json_file_count", "tsv_file_count",
    )
    _write_tsv(output_dir / "pgx_audit_summary.tsv", report["summary"], summary_fields)
    _write_tsv(output_dir / "pgx_audit_observed.tsv", report["observed"], OBSERVED_FIELDS)
    _write_tsv(output_dir / "pgx_audit_cases.tsv", report["cases"], CASE_FIELDS)
    _write_tsv(output_dir / "pgx_audit_issues.tsv", report["issues"], ISSUE_FIELDS)
    _write_tsv(output_dir / "pgx_audit_sample_files.tsv", report["sample_files"], sample_fields)
    _write_tsv(
        output_dir / "pgx_audit_source_diplotypes.tsv",
        report["source_diplotypes"],
        SOURCE_FIELDS,
    )
    _write_tsv(
        output_dir / "pgx_audit_nonreference_variants.tsv",
        report["nonreference_variants"],
        VARIANT_FIELDS,
    )
    with (output_dir / "pgx_audit_report.json").open("w", encoding="utf-8") as handle:
        json.dump(report, handle, ensure_ascii=False, indent=2)
        handle.write("\n")


def print_report(report: dict, output_dir: Path | None) -> None:
    print(f"Root: {report['root']}")
    print(f"Prefix: {report['prefix']}")
    print(
        f"Samples: {report['sample_count']} "
        f"(complete JSON+TSV pairs: {report['complete_file_pairs']})"
    )
    severities = report["severity_counts"]
    print(
        "Issues: "
        f"ERROR={severities.get('ERROR', 0)} "
        f"WARN={severities.get('WARN', 0)} "
        f"INFO={severities.get('INFO', 0)}"
    )
    print()
    print("Report row                 JSON  TSVfb   OK Wrap Review Distinct")
    for row in report["summary"]:
        print(
            f"{row['gene']:<26}"
            f"{int(row['json_source']):>5}"
            f"{int(row['tsv_fallback']):>7}"
            f"{int(row['ok']):>5}"
            f"{int(row['wraps']):>5}"
            f"{int(row['review']):>7}"
            f"{int(row['distinct_outputs']):>9}"
        )
    if report["issue_code_counts"]:
        print("\nIssue codes:")
        for code, count in report["issue_code_counts"].items():
            print(f"  {code}: {count}")
    if output_dir:
        print(f"\nAudit files: {output_dir}")


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--root",
        type=Path,
        default=Path("/home/datalake_Intermediate/pipeline/tertiary_output"),
        help="tertiary_output root",
    )
    parser.add_argument("--prefix", default="26T", help="sample directory prefix")
    parser.add_argument("--output-dir", type=Path, help="optional audit output directory")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    if not args.root.is_dir():
        print(f"error: root does not exist: {args.root}", file=sys.stderr)
        return 2
    report = audit(args.root, args.prefix)
    if args.output_dir:
        write_outputs(report, args.output_dir)
    print_report(report, args.output_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
