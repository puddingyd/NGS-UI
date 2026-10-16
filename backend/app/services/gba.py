"""DRAGEN GBA evidence and allele-level projection; never mutate 03 annotations.

The named haplotype is the unit of review. Joint GBA/GBAP1 depths are retained
as evidence, never converted to gene-specific DP/AD/VAF or ACMG points.
"""
from __future__ import annotations

import csv
import gzip
import hashlib
import json
import re
import shutil
import subprocess
from collections import Counter
from pathlib import Path

from . import sample_layout

MANAGED_NAMES = {"gba.json"}
SCHEMA_VERSION = 1
TRANSCRIPT = "NM_000157.4"
CALLER_REFERENCE = "https://help.dragen.illumina.com/dragen-v4.3/product-guide/dragen-v4.3/dragen-dna-pipeline/targeted-caller/gba-calling"
# DRAGEN v4.3 table of associated HGVS alleles. These are not a blanket
# pathogenicity classification and do not describe exact structural breakpoints.
HAPLOTYPES = {
    "A495P": ["c.1483G>C"],
    "L483P": ["c.1448T>C"],
    "D448H": ["c.1342G>C"],
    "c.1263del": ["c.1265_1319del"],
    "RecNciI": ["c.1448T>C", "c.1483G>C", "c.1497G>C"],
    "RecTL": ["c.1342G>C", "c.1448T>C", "c.1483G>C", "c.1497G>C"],
    "c.1263del+RecTL": ["c.1265_1319del", "c.1342G>C", "c.1448T>C", "c.1483G>C", "c.1497G>C"],
}
PROTEINS = {"c.1342G>C": "Asp448His", "c.1448T>C": "Leu483Pro",
            "c.1483G>C": "Ala495Pro", "c.1497G>C": "Val499="}
ALIASES = {"1263del": "c.1263del", "1263del+RecTL": "c.1263del+RecTL"}
# Whole-allele assertion, NOT the classification of a constituent SNV.
RECNCII_CLINVAR = {"id": "4297", "classification": "Pathogenic",
                  "accession": "VCV000004297.9", "checked": "2026-10-15",
                  "url": "https://www.ncbi.nlm.nih.gov/clinvar/variation/4297/"}


def _digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def _identity(build: str, kind: str, name: str, deletion=None) -> str:
    key = json.dumps([build, kind, name, deletion], separators=(",", ":"))
    return "GBA1-" + hashlib.sha256(key.encode()).hexdigest()[:24]


def _components(name: str) -> list[str]:
    if name in HAPLOTYPES:
        return HAPLOTYPES[name]
    # New caller versions can use HGVS allele IDs joined by '+'. Unknown
    # names stay visible but are never guessed into a known clinical allele.
    values = re.split(r"\+(?=NM_000157\.4:)", name)
    if all(value.startswith(TRANSCRIPT + ":c.") for value in values):
        return [value.split(":", 1)[1] for value in values]
    return []


def _canonical_haplotype(name: str) -> str:
    name = ALIASES.get(name, name)
    if name in HAPLOTYPES:
        return name
    components = set(_components(name))
    for known, values in HAPLOTYPES.items():
        # The caller table lists only associated diagnostic HGVS IDs; the
        # classical synonymous marker is part of the named Rec allele.
        if components and components in (set(values), set(values) - {"c.1497G>C"}):
            return known
    return name


def composition_text(variant: dict) -> str:
    data = variant.get("gba") or {}
    if data.get("kind") != "haplotype":
        return ""
    name = data["name"]
    components = data.get("components") or []
    if not components:
        return f"變異 {name} 為 GBA1 specialized caller 偵測的重組等位基因，組成尚待釐清。"
    c = "c.[" + ";".join(x.removeprefix("c.") for x in components) + "]"
    p = ""
    if all(x in PROTEINS for x in components):
        p = "（p.[" + ";".join(PROTEINS[x] for x in components) + "]）"
    return f"變異 {name} 包含位於同一等位基因上的 {c}{p}。"


def discover_json(source_vcf: Path, source_id: str) -> Path | None:
    candidates = [source_vcf.parent / f"{source_id}.targeted.json",
                  source_vcf.parent.parent / "other" / source_id / "germline_seq" / f"{source_id}.targeted.json"]
    found = list(dict.fromkeys(p.resolve() for p in candidates if p.is_file()))
    if len(found) > 1 and len({_digest(p) for p in found}) > 1:
        raise ValueError("GBA1: 同一來源有不同 targeted.json，請核對來源")
    return found[0] if found else None


def _vcf_lines(path: Path, build: str):
    # Indexed reading avoids scanning a WGS VCF. The fallback is streaming,
    # used for a small standalone targeted VCF or installations without tabix.
    region = "1:155234000-155245000" if build == "hg38" else "1:155204000-155216000"
    executable = shutil.which("bcftools")
    if executable and Path(str(path) + ".tbi").exists():
        header = subprocess.run([executable, "view", "-h", str(path)], check=True, capture_output=True, text=True).stdout
        prefix = "chr" if "##contig=<ID=chr1," in header else ""
        result = subprocess.run([executable, "view", "-H", "-r", prefix + region, str(path)], check=True, capture_output=True, text=True)
        yield from header.splitlines()
        yield from result.stdout.splitlines()
    else:
        opener = gzip.open if path.suffix == ".gz" else open
        with opener(path, "rt") as handle:
            yield from handle


def read_evidence(path: Path, source_id: str, build: str) -> list[dict]:
    records = []
    checked = False
    for line in _vcf_lines(path, build):
        if line.startswith("#CHROM"):
            if line.strip().split("\t")[9:] != [source_id]:
                raise ValueError("GBA1: VCF sample 與 targeted.json/source sample 不一致")
            checked = True
        if line.startswith("#"):
            continue
        cols = line.rstrip().split("\t")
        if len(cols) < 10:
            continue
        info = dict(item.split("=", 1) if "=" in item else (item, True) for item in cols[7].split(";"))
        lo, hi = (155234000, 155245000) if build == "hg38" else (155204000, 155216000)
        if cols[0].removeprefix("chr") != "1" or not lo <= int(cols[1]) <= hi:
            continue
        if "TARGETED" not in info and "TargetedConflict" not in cols[6]:
            continue
        fmt = dict(zip(cols[8].split(":"), cols[9].split(":")))
        alleles = re.split(r"[/|]", fmt.get("GT", "."))
        ids = str(info.get("ALLELE_ID", "")).split(",")
        alts = cols[4].split(",")
        for index, alt in enumerate(alts, 1):
            if str(index) not in alleles:
                continue
            allele_id = ids[index] if len(ids) == len(alts) + 1 else ""
            records.append({"id": f"chr1-{cols[1]}-{cols[3]}-{alt}",
                            "CHROM": "chr1", "POS": int(cols[1]), "REF": cols[3], "ALT": alt,
                            "allele_id": allele_id if allele_id != "." else "",
                            "copies": alleles.count(str(index)), "gt": fmt.get("GT"),
                            "filter": cols[6], "info": info, "format": fmt})
    if not checked:
        raise ValueError("GBA1: VCF 缺少 sample header")
    return records


def _annotations(raw_tsv: Path, evidence: list[dict]) -> dict:
    from ..adapters.snv_tsv import _row_to_variant, merge_snv_variant_row
    wanted = {row["id"] for row in evidence}
    out = {}
    if not wanted or not raw_tsv.is_file():
        return out
    with raw_tsv.open(encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            vid = "-".join(row.get(k, "") for k in ("CHROM", "POS", "REF", "ALT"))
            if vid in wanted and row.get("GENE") in {"GBA", "GBA1"}:
                merge_snv_variant_row(out, _row_to_variant(row))
    return out


def project(data: dict, evidence: list[dict], annotations: dict, build: str) -> tuple[dict, list[str], list[str]]:
    """Normalize all caller alleles; return cards, replaced ids, warnings."""
    from ..adapters.snv_tsv import _row_to_variant
    variants, replaced, warnings = {}, set(), []
    haplotypes = data.get("recombinantHaplotypes")
    if not isinstance(haplotypes, list) or len(haplotypes) != 2 or not all(isinstance(x, str) for x in haplotypes):
        warnings.append("GBA1 recombinant haplotypes 未完整判定")
        haplotypes = []
    counts = Counter(_canonical_haplotype(x) for x in haplotypes if x)
    cn = data.get("totalCopyNumber")
    deletion = data.get("deletionBreakpointInGene")
    if not isinstance(data.get("variants"), list):
        warnings.append("GBA1 nonrecombinant variants 未完整判定")
    if type(cn) is not int or cn < 0:
        warnings.append("GBA1 combined copy number 未判定")
    by_allele = {}
    for record in evidence:
        if record["allele_id"]:
            by_allele.setdefault(record["allele_id"], []).append(record)
    calls = [("haplotype", name, copies, "PASS", None) for name, copies in counts.items()]
    for call in data.get("variants") or []:
        if not isinstance(call, dict) or not isinstance(call.get("alleleId"), str):
            raise ValueError("GBA1: invalid nonrecombinant variant")
        copies = call.get("alleleCopyNumber")
        if type(copies) is int and copies == 0:
            continue
        calls.append(("variant", call["alleleId"], copies, call.get("filter"), call.get("genotypeQuality")))
    for kind, name, copies, filt, gq in calls:
        components = _components(name) if kind == "haplotype" else [name.split(":", 1)[-1]]
        matched = [r for c in components for r in by_allele.get(TRANSCRIPT + ":" + c, []) if "TARGETED" in r["info"]]
        # Include the synonymous RecNciI marker omitted from ALLELE_ID in v4.3.
        if kind == "haplotype" and "c.1497G>C" in components and build == "hg38":
            matched += [r for r in evidence if r["id"] == "chr1-155235203-C-G" and r["info"].get("Recombinant")]
        base = annotations.get(matched[0]["id"], {}) if len(components) == 1 and matched else {}
        v = dict(base) if base else _row_to_variant({"CHROM": "chr1", "POS": "0", "REF": "", "ALT": "", "GENE": "GBA1"})
        vid = matched[0]["id"] if kind == "variant" and len(matched) == 1 else _identity(build, kind, name, deletion if kind == "haplotype" else None)
        v.update(id=vid, gene_symbol="GBA1", refseq_transcript=TRANSCRIPT, callers="DRAGEN GBA",
                 zygosity="het" if copies == 1 else "hom" if copies == 2 and (kind == "haplotype" or cn == 4) else "",
                 HGVS_C=name if kind == "haplotype" else components[0], HGVS_P="" if kind == "haplotype" else v.get("HGVS_P", ""),
                 HGVS=f"GBA1:{TRANSCRIPT}:" + (name if kind == "haplotype" else components[0]), low_alt_support=False, low_depth=False)
        # Synthetic IDs must not masquerade as a single genomic allele.
        if kind == "haplotype" or not matched:
            v.update(REF="", ALT="", rs_id="", RS_ID="", transcript_options=[], MANE_ALL=[])
        if matched:
            v.update(CHROM="chr1", POS=min(r["POS"] for r in matched), END=max(r["POS"] + len(r["REF"]) - 1 for r in matched))
        else:
            v.update(POS=None, END=None)
        # DP/AD/AF from the conventional caller must not replace joint metrics.
        v.update(depth=None, alt_depth=None, alt_af=None, AD="")
        if kind == "variant" and len(matched) == 1:
            v.update({key: matched[0][key] for key in ("CHROM", "POS", "REF", "ALT")})
        problems = list(warnings)
        if not components:
            problems.append("未知 haplotype，需人工確認組成與分類")
        if type(copies) is not int or copies not in (1, 2) or (kind == "variant" and copies == 2 and cn != 4):
            problems.append("allele copy number／基因型需複核")
        if filt != "PASS" or any(r["filter"] != "PASS" for r in matched):
            problems.append("specialized caller FILTER 未通過")
        if kind == "variant" and matched and any(r["copies"] != copies for r in matched):
            problems.append("JSON／VCF allele copy number 不一致")
        if kind == "variant" and len(matched) != 1:
            problems.append("缺少唯一對應的 GBA1 VCF allele")
        if kind == "haplotype" and cn is None and not warnings:
            problems.append("combined copy number 未判定")
        v["gba"] = {"kind": kind, "name": name, "components": components,
                    "allele_copy_number": copies, "filter": filt, "genotype_quality": gq,
                    "total_copy_number": cn, "deletion_breakpoint_in_gene": deletion,
                    "evidence": matched, "review_required": problems,
                    "phase": "caller_haplotype" if kind == "haplotype" else "unknown"}
        if name == "RecNciI" and kind == "haplotype":
            v.update(CLNSIG="Pathogenic", CLNSIG_old="Pathogenic", clinvar_stars=2,
                     clinvar_stars_old=2, clinvar_variation_id_old="4297", tier="1A",
                     ACMG_classification="Pathogenic", ACMG_score=None)
            v["gba"]["clinical_assertion"] = RECNCII_CLINVAR
        v["gba_composition"] = composition_text(v)
        variants[vid] = v
        replaced.update(r["id"] for r in matched if r["id"] != vid)
    if counts:
        # Recombinant differentiating markers (including noncoding ones) are
        # evidence for the parent allele, not additional reportable variants.
        replaced.update(r["id"] for r in evidence if r["info"].get("Recombinant"))
    replaced.difference_update(variants)
    return variants, sorted(replaced), warnings


def _joint_support_count(rows: list[dict]) -> int:
    return sum(any(row.get("format", {}).get(key) not in (None, "", ".")
                   for key in ("JDP", "JAD", "JAF")) for row in rows)


def _select_evidence(source_vcf: Path, source_sample_id: str,
                     build: str, *, indexed_only: bool = False) -> tuple[Path, list[dict]]:
    targeted = source_vcf.parent / f"{source_sample_id}.targeted.vcf.gz"
    paths = list(dict.fromkeys(path for path in (targeted, source_vcf) if path.is_file()))
    if indexed_only:
        if not shutil.which("bcftools"):
            return source_vcf, []
        paths = [path for path in paths if Path(str(path) + ".tbi").is_file()]
    # The standalone targeted VCF can omit records retained in hard-filtered.
    # Select the richer record for each genomic ALT, retaining other GBA1
    # alleles from either file instead of choosing one file wholesale.
    selected: dict[str, tuple[tuple[int, int, int], Path, list[dict]]] = {}
    for path in paths:
        rows = read_evidence(path, source_sample_id, build)
        by_id: dict[str, list[dict]] = {}
        for row in rows:
            by_id.setdefault(row["id"], []).append(row)
        for vid, group in by_id.items():
            score = (
                _joint_support_count(group),
                sum(bool(row["allele_id"]) for row in group),
                len(group),
            )
            if vid not in selected or score > selected[vid][0]:
                selected[vid] = score, path, group
    evidence = [dict(row, source_vcf=str(path)) for _, path, group in selected.values()
                for row in group]
    chosen = next((path for path in paths if any(item[1] == path for item in selected.values())),
                  source_vcf)
    return chosen, evidence


def prepare_outputs(*, mode: str, source_vcf: str | Path, source_sample_id: str,
                    sample_id: str, post_dir: Path, raw_tsv: Path) -> dict:
    source_vcf = Path(source_vcf)
    result = {"schema_version": SCHEMA_VERSION, "sample_id": sample_id,
              "source_sample_id": source_sample_id, "status": "not_available",
              "variants": {}, "replaced_ids": [], "warnings": []}
    path = discover_json(source_vcf, source_sample_id) if mode == "dragen" else None
    if path:
        payload = json.loads(path.read_text(encoding="utf-8"))
        if payload.get("sampleId") != source_sample_id:
            raise ValueError("GBA1: targeted.json sampleId 與 source sample 不一致")
        build = {"GRCh38": "hg38", "GRCh37": "hg19"}.get(payload.get("genomeBuild"), payload.get("genomeBuild"))
        if build not in {"hg38", "hg19"}:
            raise ValueError("GBA1: unsupported genome build")
        if isinstance(payload.get("gba"), dict):
            vcf, evidence = _select_evidence(source_vcf, source_sample_id, build)
            variants, replaced, warnings = project(payload["gba"], evidence, _annotations(raw_tsv, evidence), build)
            result.update(status="complete", genome_build=build, variants=variants,
                          replaced_ids=replaced, warnings=warnings, raw_gba=payload["gba"], evidence=evidence,
                          software_version=payload.get("softwareVersion"),
                          sources={"json": str(path), "json_sha256": _digest(path), "vcf": str(vcf),
                                   "evidence_vcfs": sorted({r["source_vcf"] for r in evidence}),
                                   "evidence_vcf_signatures": {
                                       str(source): [source.stat().st_size, source.stat().st_mtime_ns]
                                       for source in {Path(r["source_vcf"]) for r in evidence}
                                   }},
                          caller_reference=CALLER_REFERENCE)
            if raw_tsv.is_file():
                stat = raw_tsv.stat()
                result["raw_signature"] = [stat.st_size, stat.st_mtime_ns]
            stat = vcf.stat()
            result["sources"]["vcf_signature"] = [stat.st_size, stat.st_mtime_ns]
    dest = sample_layout.scoped_file(post_dir, sample_id, "gba.json", for_write=True, force_prefixed=True)
    temp = dest.with_suffix(".json.tmp")
    temp.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    temp.replace(dest)
    return result


def load(sample_id: str) -> dict:
    from . import specialized
    bundle_path = sample_layout.state_file(sample_id, specialized.NAMES["manifest"])
    if bundle_path.is_file():
        bundle = specialized.load(sample_id)
        locus = bundle.get("loci", {}).get("gba", {})
        return {"status": locus.get("status", "not_available"),
                "variants": {k: v for k, v in bundle.get("variants", {}).items() if v.get("gba")},
                "replaced_ids": locus.get("replaced_ids", []), "warnings": locus.get("warnings", []),
                "raw_gba": locus.get("raw", {}), "software_version": bundle.get("software_version")}

    path = sample_layout.state_file(sample_id, "gba.json")
    if not path.is_file():
        return {"status": "not_available", "variants": {}, "replaced_ids": []}
    result = json.loads(path.read_text(encoding="utf-8"))
    if result.get("schema_version") != SCHEMA_VERSION or result.get("sample_id") != sample_id:
        raise ValueError("GBA1: 衍生結果版本／樣本不符，請重新分析")
    raw = sample_layout.snv_raw_tsv(sample_id)
    if result.get("raw_signature") and raw.is_file():
        stat = raw.stat()
        if result["raw_signature"] != [stat.st_size, stat.st_mtime_ns]:
            raise ValueError("GBA1: 原始三級結果已改變，請重新產生 GBA1 衍生結果")
    # Legacy files are self-contained. Missing evidence is never recovered at
    # request time: use the explicit post-processing backfill command instead.
    if any((v.get("gba") or {}).get("kind") == "haplotype"
           and not _joint_support_count((v.get("gba") or {}).get("evidence") or [])
           for v in (result.get("variants") or {}).values()):
        result.setdefault("warnings", []).append("舊 GBA1 結果缺 joint evidence，請補建 specialized 資料包")
    return result


def merge(variants: dict, sample_id: str, *, genes: set[str] | None = None) -> dict:
    if genes is not None and not genes.intersection({"GBA", "GBA1"}):
        return {}
    data = load(sample_id)
    for vid in data.get("replaced_ids", []):
        variants.pop(vid, None)
    additions = data.get("variants") or {}
    for vid, variant in additions.items():
        if vid in variants and variant["gba"]["kind"] == "variant":
            current = variants[vid]
            # Keep current transcript/predictor/ClinVar overlays, while the
            # specialized caller remains the sole source of genotype/support.
            caller_fields = {key: variant.get(key) for key in (
                "gba", "gba_composition", "zygosity", "depth", "alt_depth", "alt_af", "AD",
                "callers", "low_depth", "low_alt_support")}
            additions[vid] = {**variant, **current, **caller_fields}
    variants.update(additions)
    return additions


def validate_report(variants: dict, selected_ids) -> None:
    for vid in selected_ids:
        variant = variants.get(vid) or {}
        if vid.startswith("GBA1-") and not variant:
            raise ValueError("已選取的 GBA1 allele 已不在目前結果，請重新核對報告選取")
        problems = (variant.get("gba") or {}).get("review_required") or []
        if problems:
            raise ValueError("GBA1 結果尚待釐清：" + "；".join(problems))


def validate_saved_selection(sample: dict, report: dict, *, panels=None) -> None:
    selected = set()
    if panels is None:
        selected.update(vid for vid, status in (report.get("status") or {}).items() if str(status) in {"1", "2"})
    else:
        for panel in panels:
            selected.update((report.get("secondary_findings", {}).get(panel) or {}).get("selected") or [])
        selected.update(vid for vid, values in (report.get("panels") or {}).items()
                        if any(values.get(p) == "V" for p in panels))
    if selected.intersection((sample.get("gba") or {}).get("replaced_ids") or []):
        raise ValueError("已選取的 GBA1 組成位點已整合為 allele，請重新選取完整 allele 後匯出")
    validate_report(sample.get("variants") or {}, selected)
