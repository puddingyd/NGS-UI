"""Self-contained DRAGEN targeted-caller bundle. Only prepare_outputs reads inputs.

Clinical observations, locus dosage and molecular alleles remain distinct. Raw
records (including reference/no-call records) are retained for local reprocessing.
"""
from __future__ import annotations

import copy
import csv
import gzip
import hashlib
import json
import re
import shutil
import subprocess
import uuid
from collections import Counter
from datetime import datetime, timezone
from functools import lru_cache
from pathlib import Path

from . import sample_layout

SCHEMA_VERSION = 1
MAX_BYTES = 100_000_000
LOCUS_GENES = {"gba": ["GBA1"], "cyp21a2": ["CYP21A2"], "smn": ["SMN1", "SMN2"], "hba": ["HBA1", "HBA2"]}
NAMES = {part: f"specialized.{part}.{ext}" for part, ext in
         (("manifest", "json"), ("source", "json"), ("evidence", "jsonl.gz"), ("calls", "json"))}
MANAGED_NAMES = set(NAMES.values())
BASE_URL = "https://help.dragen.illumina.com/dragen-v4.3/product-guide/dragen-v4.3/dragen-dna-pipeline/targeted-caller/"
REGIONS = {
    "hg38": {"gba": ("1", 155180000, 155250000), "cyp21a2": ("6", 31990000, 32060000),
             "smn": ("5", 69200000, 71000000), "hba": ("16", 150000, 190000)},
    "hg19": {"gba": ("1", 155170000, 155230000), "cyp21a2": ("6", 31980000, 32070000),
             "smn": ("5", 69200000, 70400000), "hba": ("16", 200000, 240000)},
}
DISEASES = {"CYP21A2": "21-羥化酶缺乏相關先天性腎上腺增生症", "SMN1": "脊髓性肌肉萎縮症",
            "HBA1": "α 型地中海型貧血", "HBA2": "α 型地中海型貧血", "HBA1/HBA2": "α 型地中海型貧血"}


def _bytes(data):
    return json.dumps(data, ensure_ascii=False, sort_keys=True, separators=(",", ":")).encode()


def _sha(data):
    return hashlib.sha256(data).hexdigest()


def _read_small(path):
    with Path(path).open("rb") as f:
        data = f.read(MAX_BYTES)
    if len(data) >= MAX_BYTES:
        raise ValueError("specialized 資料達 100 MB 上限；拒絕複製或截斷")
    return data


def _write_small(path, data):
    if len(data) >= MAX_BYTES:
        raise ValueError("specialized 擷取資料達 100 MB 上限")
    Path(path).write_bytes(data)


def _id(locus, build, kind, identity):
    return f"{locus.upper()}-" + _sha(_bytes([build, kind, identity]))[:24]


def _number(value):
    try:
        n = float(value)
        return n if n >= 0 and n < float("inf") else None
    except (ValueError, TypeError):
        return None


def _cn(value):
    return value if type(value) is int and value >= 0 else None


def _locus(chrom, pos, build):
    return next((key for key, (c, lo, hi) in REGIONS[build].items()
                 if chrom.removeprefix("chr") == c and lo <= pos <= hi), None)


def _parse_record(line, source, build):
    cols = line.rstrip().split("\t")
    if len(cols) != 10:
        raise ValueError("specialized VCF 必須是單樣本且具完整欄位")
    info = dict(x.split("=", 1) if "=" in x else (x, True) for x in cols[7].split(";"))
    loc = _locus(cols[0], int(cols[1]), build)
    if not loc:
        return None
    # A region slice of hard-filtered must not turn ordinary calls into targeted calls.
    if source == "hard_filtered" and "TARGETED" not in info and not any(k in info for k in ("ALLELE_ID", "JIDS", "EVENT")):
        return None
    return {"record_id": _sha(_bytes([source, line.rstrip()]))[:24], "source": source, "locus": loc,
            "raw": line.rstrip(), "CHROM": "chr" + cols[0].removeprefix("chr"), "POS": int(cols[1]),
            "vcf_id": cols[2], "REF": cols[3], "ALTS": cols[4].split(","), "qual": cols[5],
            "filter": cols[6], "info": info, "format": dict(zip(cols[8].split(":"), cols[9].split(":")))}


def _extract(path, source, sample, build):
    """Bounded streaming; indexed WGS extraction. No unbounded capture_output."""
    executable = shutil.which("bcftools")
    indexed = executable and any(Path(str(path) + ext).is_file() for ext in (".tbi", ".csi"))
    if source == "hard_filtered" and path.stat().st_size >= MAX_BYTES and not indexed:
        raise ValueError("大型 hard-filtered VCF 需 bcftools 與索引，請補齊後重做 specialized post-processing")
    header, rows, size = [], [], 0
    def consume(lines, header_only=False):
        nonlocal size
        checked = False
        for line in lines:
            size += len(line.encode())
            if size >= MAX_BYTES:
                raise ValueError("specialized VCF 擷取內容達 100 MB 上限")
            if line.startswith("#"):
                header.append(line.rstrip())
                if line.startswith("#CHROM"):
                    if line.rstrip().split("\t")[9:] != [sample]:
                        raise ValueError("specialized VCF sample 與來源樣本不符")
                    checked = True
            elif not header_only:
                record = _parse_record(line, source, build)
                if record:
                    rows.append(record)
        return checked
    if indexed:
        def run(args, header_only=False):
            with subprocess.Popen([executable, "view", *args, str(path)], stdout=subprocess.PIPE,
                                  stderr=subprocess.DEVNULL, text=True) as proc:
                try:
                    checked = consume(proc.stdout, header_only)
                except BaseException:
                    proc.kill()
                    raise
                if proc.wait() != 0:
                    raise ValueError("specialized indexed VCF 查詢失敗")
                return checked
        if not run(["-h"], True):
            raise ValueError("specialized VCF 缺 sample header")
        prefix = "chr" if any(x.startswith("##contig=<ID=chr") for x in header) else ""
        regions = ",".join(f"{prefix}{c}:{lo}-{hi}" for c, lo, hi in REGIONS[build].values())
        run(["-H", "-r", regions])
    else:
        opener = gzip.open if path.suffix == ".gz" else open
        with opener(path, "rt") as f:
            if not consume(f):
                raise ValueError("specialized VCF 缺 sample header")
    return header, rows


def allele_evidence(records):
    """One projection per ALT; raw multiallelic record remains lossless in bundle."""
    out = []
    for r in records:
        gt = re.split(r"[/|]", r["format"].get("GT", "."))
        aids = str(r["info"].get("ALLELE_ID", "")).split(",")
        for i, alt in enumerate(r["ALTS"], 1):
            if str(i) not in gt or alt in {"*", "<NON_REF>"}:
                continue
            out.append({**r, "id": f"{r['CHROM']}-{r['POS']}-{r['REF']}-{alt}", "ALT": alt,
                        "allele_index": i, "allele_id": aids[i] if len(aids) == len(r["ALTS"]) + 1 and aids[i] != "." else "",
                        "gt": r["format"].get("GT"), "copies": gt.count(str(i))})
    return out


def _best(rows):
    grouped = {}
    for r in rows:
        grouped.setdefault(r["id"], []).append(r)
    out, conflicts = [], []
    for key, values in grouped.items():
        if len({r["copies"] for r in values}) > 1:
            conflicts.append(f"{key} 來源 genotype 不一致")
        out.append(max(values, key=lambda r: (sum(r["format"].get(k) not in (None, "", ".") for k in ("JDP", "JAD", "JAF")), bool(r["allele_id"]))))
    return out, conflicts


def _annotations(path, evidence):
    from ..adapters.snv_tsv import _row_to_variant, merge_snv_variant_row
    wanted = {r["id"] for r in evidence}
    out = {}
    if path.is_file() and wanted:
        with path.open() as f:
            for row in csv.DictReader(f, delimiter="\t"):
                vid = "-".join(row.get(k, "") for k in ("CHROM", "POS", "REF", "ALT"))
                if vid in wanted:
                    merge_snv_variant_row(out, _row_to_variant(row))
    return out


def _finding(locus, build, kind, name, *, identity=None, gene=None, evidence=(), annotations=None,
             problems=(), composition="", zygosity="", category=None, fields=None):
    from ..adapters.snv_tsv import _row_to_variant
    gene = gene or LOCUS_GENES[locus][0]
    identity = identity if identity is not None else name
    vid = _id(locus, build, kind, identity)
    base = {}
    if kind == "variant" and len(evidence) == 1:
        vid = evidence[0]["id"]
        candidate = (annotations or {}).get(vid, {})
        if candidate.get("gene_symbol") == gene:
            base = candidate
    v = copy.deepcopy(base) if base else _row_to_variant({"GENE": gene, "CHROM": "", "POS": "0", "REF": "", "ALT": ""})
    d = {"locus": locus, "kind": kind, "name": name, "genes": LOCUS_GENES[locus],
         "evidence": list(evidence), "evidence_ids": [r["record_id"] for r in evidence],
         "review_required": list(problems), "health_category": category,
         "reference": BASE_URL + locus + "-calling", **(fields or {})}
    v.update(id=vid, gene_symbol=gene, callers="DRAGEN " + locus.upper(), HGVS_C=name,
             HGVS=(f"{gene}:{name}" if name.startswith("NM_") else f"{gene}::{name}"), zygosity=zygosity, low_alt_support=False, low_depth=False,
             depth=None, alt_depth=None, alt_af=None, AD="", disease_associated=True, specialized=d, specialized_composition=composition)
    if kind != "variant" or not base:
        v.update(REF="", ALT="", rs_id="", RS_ID="", HGVS_P="", refseq_transcript="", transcript_options=[], MANE_ALL=[],
                 CLNSIG="", CLNSIG_old="", clinvar_stars=None, clinvar_stars_old=None,
                 ACMG_classification="Uncertain significance", ACMG_score=None, tier="2",
                 AF=None, AF_eas=None)
    if name.startswith("NM_") and ":" in name and not base:
        v["refseq_transcript"] = name.split(":", 1)[0]
    if evidence:
        v.update(CHROM=evidence[0]["CHROM"], POS=min(r["POS"] for r in evidence), END=max(r["POS"] for r in evidence))
    else:
        v.update(CHROM="", POS=None, END=None)
    if kind == "variant" and len(evidence) == 1:
        v.update({k: evidence[0][k] for k in ("CHROM", "POS", "REF", "ALT")})
    v["Disease1"] = DISEASES.get(gene, "") + "(AR)"
    v["OMIM_disease"] = v["Disease1"]
    v["Inheritance"] = "AR"
    v["specialized"]["review_token"] = _sha(_bytes([build, identity, d]))[:24]
    return v


def _cyp(node, build, rows, annotations):
    findings, replaced = {}, set()
    warnings = []
    cn = _cn(node.get("totalCopyNumber"))
    haps = node.get("recombinantHaplotypes")
    if not isinstance(haps, list) or len(haps) != 2 or not all(isinstance(x, str) for x in haps):
        warnings.append("CYP21A2 recombinant haplotypes 未完整判定")
        haps = []
    if cn is None:
        warnings.append("CYP21A2 combined CN 未判定")
    elif cn != 4:
        warnings.append("CYP21A2 combined CN 異常，需核對 allele 劑量與基因型")
    if not isinstance(node.get("variants"), list):
        warnings.append("CYP21A2 variants 未完整判定")
    calls = [("haplotype", name, count, None, None) for name, count in Counter("+".join(sorted(re.split(r"\+(?=NM_)", x))) for x in haps if x).items()]
    calls += [("variant", x.get("alleleId", ""), x.get("alleleCopyNumber"), x.get("filter"), x.get("genotypeQuality"))
              for x in node.get("variants") or []]
    for kind, name, copies, filt, gq in calls:
        if copies == 0:
            continue
        parts = re.split(r"\+(?=NM_)", name)
        transcript_ids = {p.split(":")[0] for p in parts}
        matched = [r for r in rows if r["allele_id"] in parts]
        matched, conflicts = _best(matched)
        problems = list(warnings) + conflicts
        if not name or len(transcript_ids) != 1 or not all(re.fullmatch(r"NM_000500\.\d+:c\..+", p) for p in parts):
            problems.append("未知 CYP21A2 allele 定義，需核對組成")
        if copies not in (1, 2) or (kind == "variant" and copies == 2 and cn != 4):
            problems.append("CYP21A2 allele copy number 需複核")
        if kind == "variant" and filt != "PASS":
            problems.append("CYP21A2 FILTER 未通過")
        if set(parts) - {r["allele_id"] for r in matched} or any(r["filter"] != "PASS" for r in matched):
            problems.append("缺少完整 PASS allele 證據")
        if kind == "variant" and any("JGT" in r["format"] for r in matched):
            problems.append("同源區小變異的 CYP21A2 定位需確認")
        c = "c.[" + ";".join(p.split(":c.")[-1] for p in parts) + "]" if kind == "haplotype" else name
        composition = (f"變異 {c} 的組成由 specialized caller 判定在同一等位基因上。" if kind == "haplotype"
                       else f"變異 {name} 的 allele copies 為 {copies}；同源區證據須配合基因定位判讀。")
        v = _finding("cyp21a2", build, kind, c, identity=[sorted(parts), node.get("deletionBreakpointInGene")],
                     evidence=matched, annotations=annotations, problems=problems, composition=composition,
                     zygosity=("het" if copies == 1 else "hom") if copies in (1, 2) and (kind == "haplotype" or cn == 4) else "",
                     fields={"components": parts, "allele_copy_number": copies, "total_copy_number": cn,
                             "deletion_breakpoint_in_gene": node.get("deletionBreakpointInGene"), "filter": filt,
                             "genotype_quality": gq, "hard_blocked": not name or len(transcript_ids) != 1 or not all(re.fullmatch(r"NM_000500\.\d+:c\..+", p) for p in parts) or copies not in (1, 2) or cn is None, "phase": "caller_haplotype" if kind == "haplotype" else "unknown"})
        if len(transcript_ids) == 1 and name.startswith("NM_"):
            transcript = next(iter(transcript_ids))
            v["refseq_transcript"] = transcript
            if kind == "haplotype":
                v["HGVS"] = f"CYP21A2:{transcript}:{c}"
        findings[v["id"]] = v
        if kind == "haplotype" and not problems:
            replaced.update(r["id"] for r in matched)
    return findings, replaced, warnings


def _smn(node, build, rows):
    cn1, cn2, delta = (_cn(node.get(k)) for k in ("smn1CopyNumber", "smn2CopyNumber", "smn2Delta78CopyNumber"))
    warnings = [] if cn1 is not None else ["SMN1 copy number 未判定"]
    marker = next((x for x in node.get("variants") or [] if x.get("hgvs") == "NM_000344.4:c.*3+80T>G"), {})
    mcn = _cn(marker.get("altCopyNumber"))
    fields = {"smn1_copy_number": cn1, "smn2_copy_number": cn2, "smn2_delta78_copy_number": delta,
              "total_copy_number_float": _number(node.get("totalCopyNumberFloat", node.get("totalCopyNumber"))),
              "full_length_copy_number_float": _number(node.get("fullLengthCopyNumberFloat", node.get("fullLengthCopyNumber"))),
              "marker": marker, "phase": "unknown", "filter": None}
    summary = f"SMN1 CN: {cn1 if cn1 is not None else '未知'} · SMN2 CN: {cn2 if cn2 is not None else '未知'} · SMN2Δ7–8 CN: {delta if delta is not None else '未知'}"
    summary += f" · c.*3+80T>G ALT CN: {mcn if mcn is not None else '未知'}"
    findings = {}
    if cn1 in (0, 1):
        text = f"SMN1 完整拷貝數為 {cn1}；SMN2 完整拷貝數為 {cn2 if cn2 is not None else '未知'}。此結果不提供缺失的精確斷點。"
        v = _finding("smn", build, "copy_number", "SMN1 拷貝數減少", identity=["SMN1", cn1],
                     evidence=rows, composition=text, zygosity=f"CN = {cn1}", category="carrier" if cn1 == 1 else "risk", fields=fields)
        findings[v["id"]] = v
    elif mcn is not None and mcn > 0:
        text = "檢出 SMN1 c.*3+80T>G 相關 marker；此 marker 不能單獨確診 2+0 帶因者，須結合 SMN1 拷貝數判讀。"
        v = _finding("smn", build, "risk_marker", "c.*3+80T>G", evidence=rows, composition=text, fields=fields,
                     problems=["silent-carrier associated marker 不等同確診帶因者"])
        findings[v["id"]] = v
    return findings, set(), warnings, {**fields, "summary": summary,
        "note": "SMN1 CN ≥2 仍不能完全排除 2+0 或其他未涵蓋變異；SMN2 CN 不是 SMN1 疾病 allele。"}


def _hba(node, build, rows, annotations):
    raw_gt = node.get("genotype")
    gt = raw_gt if isinstance(raw_gt, str) else ""
    parts = gt.split("/")
    copies = {"aa": 2, "-a3.7": 1, "-a4.2": 1, "--": 0, "aaa3.7": 3, "aaa4.2": 3}
    supported = {tuple(sorted(x.split("/"))) for x in (
        "aa/aa", "aaa3.7/aa", "aaa4.2/aa", "-a3.7/aa", "-a4.2/aa", "--/aaa3.7", "--/aaa4.2",
        "-a3.7/-a3.7", "-a4.2/-a4.2", "-a3.7/-a4.2", "--/aa", "--/-a3.7", "--/-a4.2", "--/--")}
    known = tuple(sorted(parts)) in supported
    warnings = []
    if not known:
        warnings.append("HBA genotype 未完整判定或不在已支援定義內")
    if node.get("genotypeFilter") != "PASS":
        warnings.append("HBA genotype FILTER 未通過")
    findings, replaced = {}, set()
    total = sum(copies[p] for p in parts) if known else None
    if known and _cn(node.get("totalCopyNumber")) is not None and total != node["totalCopyNumber"]:
        warnings.append("HBA genotype 與 total CN 不一致")
    fields = {"genotype": gt or None, "total_copy_number": node.get("totalCopyNumber"),
              "genotype_quality": node.get("genotypeQuality", node.get("genotypeQual")),
              "filter": node.get("genotypeFilter"), "hard_blocked": not known, "phase": "caller_genotype"}
    if gt not in ("", "None", "aa/aa"):
        category = "risk" if total in (0, 1) else "carrier" if total in (2, 3) and known else None
        display = gt.replace("aaa", "ααα").replace("aa", "αα").replace("a", "α")
        text = f"HBA1／HBA2 locus genotype 為 {display}；保留兩條 haplotype 的組合。未提供的缺失名稱與精確斷點不予推定。"
        if any(p.startswith("aaa") for p in parts) and "--" not in parts:
            category = None
        v = _finding("hba", build, "locus_genotype", gt, identity=sorted(parts), gene="HBA1/HBA2",
                     evidence=rows, problems=warnings, composition=text, zygosity=gt, category=category, fields=fields)
        findings[v["id"]] = v
    for call in node.get("variants") or []:
        if call.get("alleleCopyNumber") == 0:
            continue
        name = call.get("alleleId", "")
        matched, conflicts = _best([r for r in rows if r["allele_id"] == name])
        gene = "HBA2" if name.startswith("NM_000517.") else "HBA1" if name.startswith("NM_000558.") else "HBA1/HBA2"
        problems = list(warnings) + conflicts
        if not name or not re.fullmatch(r"NM_000(?:517|558)\.\d+:c\..+", name) or _cn(call.get("alleleCopyNumber")) is None:
            problems.append("HBA allele 定義或 copy number 未完整判定")
        if call.get("filter") != "PASS":
            problems.append("HBA 小變異 FILTER 未通過")
        # DRAGEN HBA small variants use both homologous loci; a single VCF
        # coordinate (even without JGT) does not establish locus specificity.
        problems.append("HBA1／HBA2 同源區小變異定位與相位需確認")
        v = _finding("hba", build, "variant", name, gene=gene, evidence=matched, annotations=annotations,
                     problems=problems, composition=f"{name} 的 allele copies 為 {call.get('alleleCopyNumber')}；須與 locus genotype {gt or '未知'} 合併判讀。",
                     fields={**fields, "hard_blocked": not known or not name or gene == "HBA1/HBA2" or _cn(call.get("alleleCopyNumber")) is None,
                             "allele_copy_number": call.get("alleleCopyNumber"), "filter": call.get("filter"), "phase": "unknown"})
        findings[v["id"]] = v
    return findings, replaced, warnings


def project(payload, records, annotations, build):
    from . import gba
    variants, replaced, loci = {}, set(), {}
    alleles = allele_evidence(records)
    signatures = {}
    for r in records:
        gt = re.split(r"[/|]", r["format"].get("GT", "."))
        for index, alt in enumerate(r["ALTS"], 1):
            if alt in {"*", "<NON_REF>"}:
                continue
            key = f"{r['CHROM']}-{r['POS']}-{r['REF']}-{alt}"
            signatures.setdefault(key, set()).add(None if "." in gt else gt.count(str(index)))
    conflicting = {key for key, values in signatures.items() if len(values) > 1}
    for locus, genes in LOCUS_GENES.items():
        node = payload.get(locus)
        if not isinstance(node, dict):
            loci[locus] = {"genes": genes, "status": "not_available", "warnings": ["來源未提供此 caller 結果"]}
            continue
        rows = [r for r in alleles if r["locus"] == locus]
        summary = {}
        if locus == "gba":
            # Preserve existing biological IDs and whole-allele curation.
            lo, hi = (155234000, 155245000) if build == "hg38" else (155204000, 155216000)
            best, conflicts = _best([r for r in rows if lo <= r["POS"] <= hi])
            for r in best:
                r["info"] = {**r["info"], "TARGETED": True}
            vs, rep, warn = gba.project(node, best, annotations, build)
            warn += conflicts
            for v in vs.values():
                v["gba"]["review_required"].extend(conflicts)
        elif locus == "cyp21a2":
            vs, rep, warn = _cyp(node, build, rows, annotations)
        elif locus == "smn":
            vs, rep, warn, summary = _smn(node, build, rows)
        else:
            vs, rep, warn = _hba(node, build, rows, annotations)
        for v in vs.values():
            d = v.get("specialized") or v.get("gba") or {}
            if any(r["id"] in conflicting for r in d.get("evidence", [])):
                d.setdefault("review_required", []).append("targeted／hard-filtered genotype 不一致")
            if d.get("locus"):
                d["review_token"] = _sha(_bytes([node, d.get("evidence", []), d.get("review_required", [])]))[:24]
        variants.update(vs)
        replaced.update(rep)
        loci[locus] = {"genes": genes, "status": "review_required" if warn or any((v.get("specialized") or v.get("gba") or {}).get("review_required") for v in vs.values()) else "complete", "warnings": warn,
                       "evidence_status": "available" if any(r["locus"] == locus for r in records) else "not_provided",
                       "raw": node, "finding_ids": list(vs), "replaced_ids": sorted(rep), **summary}
    replaced.difference_update(variants)
    return {"variants": variants, "replaced_ids": sorted(replaced), "loci": loci}


def prepare_outputs(*, mode, source_vcf, source_sample_id, sample_id, post_dir, raw_tsv):
    from . import gba
    source_vcf, raw_tsv, post_dir = Path(source_vcf), Path(raw_tsv), Path(post_dir)
    post_dir.mkdir(parents=True, exist_ok=True)
    path = gba.discover_json(source_vcf, source_sample_id) if mode == "dragen" else None
    payload, records, evidence, sources = {}, [], [], {}
    build = "hg38"
    if path:
        content = _read_small(path)
        payload = json.loads(content)
        if payload.get("sampleId") != source_sample_id:
            raise ValueError("specialized JSON sampleId 與來源樣本不符")
        build = {"GRCh38": "hg38", "GRCh37": "hg19"}.get(payload.get("genomeBuild"), payload.get("genomeBuild"))
        if build not in REGIONS:
            raise ValueError("不支援的 specialized genome build")
        sources["json"] = {"path": str(path), "sha256": _sha(content), "size": len(content)}
        targeted = source_vcf.parent / f"{source_sample_id}.targeted.vcf.gz"
        paths = [("targeted", targeted), ("hard_filtered", source_vcf)]
        for source, vcf in paths:
            if not vcf.is_file():
                continue
            header, rows = _extract(vcf, source, source_sample_id, build)
            stat = vcf.stat()
            sources[source] = {"path": str(vcf), "size": stat.st_size, "mtime_ns": stat.st_mtime_ns,
                               "header_sha256": _sha(_bytes(header)), "records_sha256": _sha(_bytes(rows))}
            evidence.append({"type": "header", "source": source, "lines": header})
            evidence.extend({"type": "record", **r} for r in rows)
            records.extend(rows)
    result = project(payload, records, _annotations(raw_tsv, allele_evidence(records)), build)
    return write_bundle(result=result, payload=payload, evidence=evidence, sources=sources, build=build,
                        sample_id=sample_id, source_sample_id=source_sample_id, post_dir=post_dir,
                        status="complete" if path else "not_available")


def write_bundle(*, result, payload, evidence, sources, build, sample_id, source_sample_id, post_dir, status):
    generation = uuid.uuid4().hex
    common = {"schema_version": SCHEMA_VERSION, "sample_id": sample_id, "source_sample_id": source_sample_id,
              "generation_id": generation, "genome_build": build, "software_version": payload.get("softwareVersion"),
              "status": status}
    result.update(common)
    files = {"source": _bytes(payload), "calls": _bytes(result),
             "evidence": gzip.compress(b"\n".join(_bytes(r) for r in evidence), mtime=0)}
    if len(_bytes(evidence)) >= MAX_BYTES or sum(map(len, files.values())) >= MAX_BYTES:
        raise ValueError("specialized 資料包達 100 MB 上限")
    manifest = {**common, "created_at": datetime.now(timezone.utc).isoformat(), "sources": sources,
                "parser_version": "1", "files": {}, "loci": {k: v["status"] for k, v in result["loci"].items()}}
    for part, data in files.items():
        dest = post_dir / f"{sample_id}.{NAMES[part]}"
        _write_small(dest, data)
        manifest["files"][part] = {"name": dest.name, "sha256": _sha(data), "size": len(data)}
    # Caller publishes staging transactionally. Manifest is written last.
    _write_small(post_dir / f"{sample_id}.{NAMES['manifest']}", _bytes(manifest))
    return result


def rebuild_local(*, sample_id, post_dir, raw_tsv):
    """Re-parse an existing snapshot without resolving any original source path."""
    calls = load(sample_id)  # validate all local checksums first
    if calls["status"] == "not_available":
        raise ValueError("沒有可重新解析的 specialized 快照")
    manifest_path = sample_layout.state_file(sample_id, NAMES["manifest"])
    manifest = json.loads(_read_small(manifest_path))
    payload = json.loads(_read_small(manifest_path.parent / manifest["files"]["source"]["name"]))
    evidence = []
    with gzip.open(manifest_path.parent / manifest["files"]["evidence"]["name"], "rt") as f:
        size = 0
        for line in f:
            size += len(line.encode())
            if size >= MAX_BYTES:
                raise ValueError("specialized evidence 解壓超過上限")
            if line.strip():
                evidence.append(json.loads(line))
    records = [r for r in evidence if r.get("type") == "record"]
    build = calls["genome_build"]
    result = project(payload, records, _annotations(Path(raw_tsv), allele_evidence(records)), build)
    return write_bundle(result=result, payload=payload, evidence=evidence, sources=manifest["sources"], build=build,
                        sample_id=sample_id, source_sample_id=calls["source_sample_id"], post_dir=Path(post_dir), status="complete")


@lru_cache(maxsize=8)
def _load_cached(path_string, signatures):
    path = Path(path_string)
    manifest = json.loads(_read_small(path))
    if manifest.get("schema_version") != SCHEMA_VERSION:
        raise ValueError("specialized schema 不支援，請重新匯入")
    if set(manifest.get("files", {})) != {"source", "calls", "evidence"}:
        raise ValueError("specialized manifest 缺少必要檔案")
    contents = {}
    for key, item in manifest["files"].items():
        name = item["name"]
        if Path(name).name != name:
            raise ValueError("specialized bundle 路徑必須位於同一資料夾")
        content = _read_small(path.parent / name)
        if _sha(content) != item["sha256"]:
            raise ValueError("specialized 資料包 checksum 不符，請重新匯入")
        contents[key] = content
    calls = json.loads(contents["calls"])
    if any(calls.get(key) != manifest.get(key) for key in
           ("generation_id", "schema_version", "sample_id", "source_sample_id", "genome_build")):
        raise ValueError("specialized generation／sample／schema 不一致")
    return calls


def load(sample_id):
    path = sample_layout.state_file(sample_id, NAMES["manifest"])
    if not path.is_file():
        return {"status": "not_available", "variants": {}, "replaced_ids": [], "loci": {}}
    # Signatures are strictly local; provenance paths are never dereferenced.
    signatures = tuple((p.name, p.stat().st_size, p.stat().st_mtime_ns)
                       for p in [path] + [path.parent / f"{sample_id}.{NAMES[k]}" for k in ("source", "calls", "evidence")])
    result = copy.deepcopy(_load_cached(str(path), signatures))
    if result.get("sample_id") != sample_id:
        raise ValueError("specialized 資料包樣本不符")
    return result


def merge(variants, sample_id, *, genes=None):
    result = load(sample_id)
    additions = {k: v for k, v in result["variants"].items()
                 if not v.get("gba") and (genes is None or set((v.get("specialized") or {}).get("genes", [])).intersection(genes))}
    for locus, loc in result.get("loci", {}).items():
        if locus != "gba" and (genes is None or set(LOCUS_GENES[locus]).intersection(genes)):
            for vid in loc.get("replaced_ids", []):
                variants.pop(vid, None)
    for vid, v in additions.items():
        if vid in variants and v["specialized"]["kind"] == "variant" and variants[vid].get("gene_symbol") == v.get("gene_symbol"):
            current = variants[vid]
            preserve = {k: v[k] for k in ("specialized", "specialized_composition", "zygosity", "depth", "alt_depth", "alt_af", "AD", "callers")}
            additions[vid] = {**v, **current, **preserve}
    variants.update(additions)
    return additions


def composition_text(v):
    from . import gba
    return v.get("specialized_composition") or gba.composition_text(v)


def summary(sample_id):
    pipeline_type = str(sample_layout._pipeline_source(sample_id).get("pipeline_type") or "").strip().lower()
    if not pipeline_type:
        name = sample_id.lower()
        pipeline_type = "dragen" if name.endswith("-dragen") else (
            "nckuh" if name.endswith(("-nckuh", "-inhouse")) else ""
        )
    r = load(sample_id)
    result = {k: r.get(k) for k in ("status", "loci", "replaced_ids", "software_version")}
    result["pipeline_type"] = pipeline_type
    return result


def validate_selection(sample, report, *, panels=None, selected_ids=None):
    selected = set(selected_ids or [])
    if panels is None:
        selected.update(k for k, v in (report.get("status") or {}).items() if str(v) in {"1", "2"})
    else:
        for panel in panels:
            selected.update((report.get("secondary_findings", {}).get(panel) or {}).get("selected") or [])
        selected.update(k for k, v in (report.get("panels") or {}).items() if any(v.get(p) == "V" for p in panels))
    replaced = set((sample.get("specialized") or {}).get("replaced_ids") or [])
    if selected & replaced:
        raise ValueError("specialized 組成位點已整合，請重新選取完整 allele")
    selected_hba = [vid for vid in selected if ((sample.get("variants") or {}).get(vid, {}).get("specialized") or {}).get("locus") == "hba"]
    hba_categories = set()
    for vid in selected:
        v = (sample.get("variants") or {}).get(vid)
        if not v and re.match(r"(?:CYP21A2|SMN|HBA)-", vid):
            raise ValueError("已選取的 specialized 結果不在目前資料包，請重新核對")
        d = (v or {}).get("specialized") or {}
        if d.get("kind") == "haplotype" and selected.intersection(r["id"] for r in d.get("evidence", [])):
            raise ValueError("不可同時回報 specialized allele 與其組成位點，請核對選取")
        if d.get("review_required"):
            edits = (report.get("edits") or {}).get(vid) or {}
            confirmed = edits.get("specialized_review") == d.get("review_token") and bool(str(edits.get("comment") or "").strip())
            if d.get("hard_blocked") or d.get("kind") == "risk_marker" or not confirmed:
                raise ValueError("specialized 結果尚待釐清：" + "；".join(d["review_required"]))
        if panels is not None and d.get("locus") == "hba":
            edits = (report.get("edits") or {}).get(vid) or {}
            category = d.get("health_category")
            if d.get("kind") == "variant" or category is None or len(selected_hba) > 1 or edits.get("specialized_health_category"):
                category = edits.get("specialized_health_category")
                if (category not in {"carrier", "risk"} or edits.get("specialized_health_token") != d.get("review_token")
                        or not str(edits.get("comment") or "").strip()):
                    raise ValueError("HBA 複合結果需核對完整 locus genotype，指定健檢第一／第二類並於 Comment 記錄依據")
            hba_categories.add(category)
        if panels is not None and (d.get("kind") in {"copy_number", "locus_genotype"} or d.get("locus") == "hba"):
            edits = (report.get("edits") or {}).get(vid) or {}
            manual = edits.get("manual_acmg") or {}
            classification = manual.get("classification") or (v or {}).get("effective_acmg_class") or (v or {}).get("ACMG_classification")
            if classification not in {"Pathogenic", "Likely pathogenic"}:
                raise ValueError("specialized CN／locus 結果須完成 P/LP 人工分類後才能納入健檢疾病報告")
    if len(hba_categories) > 1:
        raise ValueError("同一 HBA locus 的健檢分類不一致，請合併核對")


def report_limitations(summary, genes):
    messages = []
    for locus, loc in (summary.get("loci") or {}).items():
        if set(loc.get("genes", [])).intersection(genes) and loc.get("status") != "complete":
            messages.append(f"{'／'.join(loc['genes'])} specialized caller 結果未完整判定，不能據此視為陰性。")
    return messages
