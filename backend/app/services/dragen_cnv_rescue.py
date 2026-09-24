"""DRAGEN CNV review from integrated CNV/SV PASS calls plus Rule B.

No calling, ACMG override, or modification of the 00-07 pipeline artifacts.
Coordinates used by the matching rule are VCF (POS, END], i.e. END-POS bp.
"""
from __future__ import annotations

import csv
import gzip
import hashlib
import json
import os
import shutil
import tempfile
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

POLICY_VERSION = 2
ALLOWED_FILTERS = frozenset({"cnvLength", "cnvQual"})
REVIEW_NAME = "cnv.review.tsv"
RESCUED_NAME = "cnv.rescued.annotated.tsv"
MANIFEST_NAME = "cnv_rescue.json"
EVIDENCE_COLUMN = "NGS_UI_CNV_RESCUE"
MANAGED_NAMES = {REVIEW_NAME, RESCUED_NAME, MANIFEST_NAME}
REVERSE_LINKS = {
    "LEFT_BND_OF": "LEFT_BND", "END_LEFT_BND_OF": "LEFT_BND",
    "RIGHT_BND_OF": "RIGHT_BND", "END_RIGHT_BND_OF": "RIGHT_BND",
}


def _chrom(value: str) -> str:
    return value.removeprefix("chr")


@dataclass
class Record:
    fields: list[str]
    info: dict[str, str]

    @property
    def id(self) -> str:
        return self.fields[2]

    @property
    def call(self) -> dict[str, str]:
        return dict(zip(self.fields[8].split(":"), self.fields[9].split(":")))

    @property
    def kind(self) -> str:
        return self.fields[4].strip("<>")

    @property
    def key(self) -> tuple[str, int, int, str]:
        return (_chrom(self.fields[0]), int(self.fields[1]),
                int(self.info.get("END", self.fields[1])), self.kind)

    def evidence(self) -> dict:
        return {"id": self.id, "chrom": self.fields[0], "pos": int(self.fields[1]),
                "end": self.key[2], "type": self.kind, "qual": self.fields[5],
                "filter": self.fields[6], "format": self.call}


def read_vcf(path: Path, source_sample: str) -> tuple[list[str], list[Record]]:
    """Require the exact single sample; never accidentally use the first proband."""
    headers, records = [], []
    opener = gzip.open if path.suffix == ".gz" else open
    sample_checked = False
    with opener(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                headers.append(line.rstrip("\n"))
                if line.startswith("#CHROM\t"):
                    if line.rstrip().split("\t")[9:] != [source_sample]:
                        raise ValueError(f"CNV rescue sample mismatch: {path}")
                    sample_checked = True
                continue
            if not sample_checked:
                raise ValueError(f"Missing VCF sample header: {path}")
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 10:
                raise ValueError(f"Invalid single-sample VCF record: {path}")
            info = dict(item.split("=", 1) if "=" in item else (item, "")
                        for item in fields[7].split(";"))
            record = Record(fields, info)
            # Validate interval values even when the record will be rejected.
            record.key
            records.append(record)
    if not sample_checked:
        raise ValueError(f"Missing VCF sample header: {path}")
    return headers, records


def _has_alt(gt: str) -> bool:
    return any(a.isdigit() and int(a) > 0 for a in gt.replace("|", "/").split("/"))


def _genotype_class(gt: str) -> str | None:
    alleles = gt.replace("|", "/").split("/")
    if not all(a.isdigit() for a in alleles):
        return None
    if all(a == "0" for a in alleles):
        return "ref"
    if len(alleles) == 1:
        return "hemi"
    return "hom" if len(set(alleles)) == 1 else "het"


def _supported(cnv: Record, sv: Record) -> tuple[float, float] | None:
    cchrom, start, end, kind = cnv.key
    schrom, left, right, skind = sv.key
    if (cchrom != schrom or kind != skind or end <= start or right <= left
            or sv.fields[6] != "PASS" or sv.call.get("FT", "PASS") != "PASS"):
        return None
    cgt, sgt = cnv.call.get("GT", ""), sv.call.get("GT", "")
    if not _has_alt(cgt) or not _has_alt(sgt):
        return None
    cclass, sclass = _genotype_class(cgt), _genotype_class(sgt)
    if cclass and sclass and cclass != sclass:
        return None
    overlap = max(0, min(end, right) - max(start, left))
    # Integer comparison includes exactly 50%; no rounding at the boundary.
    if 2 * overlap < end - start or 2 * overlap < right - left:
        return None
    return overlap / (end - start), overlap / (right - left)


def select_integrated_cnvs(integrated: list[Record],
                           identity_ids: dict | None = None) -> tuple[list[dict], dict]:
    """Select using the integrated file alone; legacy IDs never affect eligibility."""
    identity_ids = identity_ids or {}
    links = defaultdict(list)
    for record in integrated:
        if record.info.get("SVCLAIM") != "J" or record.kind not in ("DEL", "DUP"):
            continue
        for reverse, forward in REVERSE_LINKS.items():
            for cnv_id in record.info.get(reverse, "").split(","):
                if cnv_id:
                    links[cnv_id].append((record, forward))
    selected, counts, seen = [], Counter(), set()
    for cnv in integrated:
        if cnv.info.get("SVTYPE") != "CNV" or cnv.fields[4] not in ("<DEL>", "<DUP>"):
            continue
        if cnv.info.get("SVCLAIM") not in ("D", "DJ"):
            raise ValueError(f"Missing or unsupported SVCLAIM for CNV: {cnv.id}")
        if cnv.key[0] not in {str(n) for n in range(1, 23)} | {"X", "Y"}:
            counts["non_primary_contig"] += 1
            continue
        key = (cnv.key[0], int(cnv.info.get("OrigCnvPos", cnv.fields[1])),
               int(cnv.info.get("OrigCnvEnd", cnv.info["END"])), cnv.kind)
        rule, support = "", []
        if cnv.fields[6] == "PASS":
            rule = "PASS"
            if cnv.info.get("MatchSv") not in (None, "", "."):
                support = [{"original_sv_id": cnv.info["MatchSv"]}]
        elif (cnv.info.get("SVCLAIM") == "DJ"
              and set(cnv.fields[6].split(";")) <= ALLOWED_FILTERS):
            for sv, forward in links[cnv.id]:
                reference = cnv.info.get(forward, "")
                if not reference or reference == ".":
                    continue
                overlap = _supported(cnv, sv)
                if overlap is not None:
                    support.append({**sv.evidence(), "link": forward,
                                    "original_sv_id": reference.removesuffix(".end"),
                                    "cnv_overlap": overlap[0], "sv_overlap": overlap[1]})
            if support:
                rule = "B"
        if not rule:
            counts["no_qualifying_sv_support"] += 1
            continue
        if key in seen:
            raise ValueError(f"Multiple integrated records for CNV identity: {key}")
        seen.add(key)
        stable_id = identity_ids.get(key) or f"CNVRESCUE-{key[0]}-{key[1]}-{key[2]}-{key[3]}"
        evidence = {"policy_version": POLICY_VERSION, "rule": rule,
                    "svclaim": cnv.info["SVCLAIM"], "integrated": cnv.evidence(),
                    "sv_support": support}
        selected.append({"id": stable_id, "record": cnv, "evidence": evidence})
        counts["integrated_pass" if rule == "PASS" else "rescued_B"] += 1
    return selected, dict(counts)


def _signature(path: Path) -> dict:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    stat = path.stat()
    return {"path": str(path.resolve()), "size": stat.st_size, "sha256": digest.hexdigest()}


def _read_tsv(path: Path) -> tuple[list[str], list[dict]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        headers = reader.fieldnames or []
        if not {"AnnotSV_ID", "Annotation_mode", "SV_chrom", "SV_start", "SV_end", "SV_type"} <= set(headers):
            raise ValueError(f"Incomplete AnnotSV header: {path}")
        rows = list(reader)
        if any(None in row for row in rows):
            raise ValueError(f"Malformed AnnotSV rows: {path}")
        return headers, rows


def _write_tsv(path: Path, headers: list[str], rows: list[dict]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=headers, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def legacy_identity_ids(base: Path | None) -> dict:
    """Read only the old event identities, never reuse annotation or filter values."""
    if base is None or not base.is_file():
        return {}
    ids = {}
    with base.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row.get("Annotation_mode") != "full":
                continue
            key = (_chrom(row["SV_chrom"]), int(row["SV_start"]),
                   int(row["SV_end"]), row["SV_type"])
            if key in ids:
                raise ValueError(f"Ambiguous legacy CNV identity: {key}")
            ids[key] = row["AnnotSV_ID"]
    return ids


def validate_annotations(annotated: Path, selected: list[dict],
                         review_out: Path, source_sample: str) -> dict:
    """Publish exclusively fresh annotations, with one full row per selected event."""
    new_header, new_rows = _read_tsv(annotated)
    if ("FORMAT" not in new_header or new_header.index("FORMAT") + 1 >= len(new_header)
            or new_header[new_header.index("FORMAT") + 1] != source_sample):
        raise ValueError("AnnotSV sample column does not match source sample")
    by_id = {item["id"]: item for item in selected}
    if len(by_id) != len(selected):
        raise ValueError("Duplicate stable CNV IDs")
    full_seen = set()
    for row in new_rows:
        rescue_id = row.get("ID", "")
        if rescue_id not in by_id:
            raise ValueError("AnnotSV lost rescue VCF ID; cannot safely join evidence")
        row["AnnotSV_ID"] = rescue_id
        if row["Annotation_mode"] == "full":
            if rescue_id in full_seen:
                raise ValueError(f"Duplicate full annotation for {rescue_id}")
            full_seen.add(rescue_id)
            record = by_id[rescue_id]["record"]
            key = (_chrom(row["SV_chrom"]), int(row["SV_start"]),
                   int(row["SV_end"]), row["SV_type"])
            if key != record.key or row.get("FILTER") != record.fields[6]:
                raise ValueError(f"AnnotSV changed integrated CNV coordinates/type/FILTER: {rescue_id}")
            row[EVIDENCE_COLUMN] = json.dumps(by_id[rescue_id]["evidence"], separators=(",", ":"))
        else:
            row[EVIDENCE_COLUMN] = ""
    if full_seen != set(by_id):
        raise ValueError(f"AnnotSV omitted {len(set(by_id) - full_seen)} selected CNVs")

    headers = list(new_header)
    if EVIDENCE_COLUMN not in headers:
        headers.append(EVIDENCE_COLUMN)
    _write_tsv(review_out, headers, new_rows)
    return {"annotated_events": len(full_seen)}


def annotsv_command(vcf: Path, output: Path, ngs_home: Path) -> list[str]:
    """Use an explicit/local AnnotSV or the same 3.5.10 image as tertiary."""
    binary = os.environ.get("ANNOTSV_BIN")
    local = ngs_home / "biotools/AnnotSV/bin/AnnotSV"
    if not binary and local.is_file():
        binary = str(local)
    annotation_env = os.environ.get("ANNOTSV_ANNOTATIONS")
    candidates = [Path(annotation_env)] if annotation_env else [
        Path("/home/pipeline/reference/hg38/tertiary/annotsv_annotations/share/AnnotSV"),
        ngs_home / "biotools/AnnotSV/share/AnnotSV",
    ]
    annotations = next((p for p in candidates if p.is_dir()), None)
    if annotations is None:
        raise FileNotFoundError("AnnotSV annotations missing; set ANNOTSV_ANNOTATIONS")
    if binary:
        command = [binary]
    else:
        image_env = os.environ.get("NGS_UI_ANNOTSV_SIF")
        images = [Path(image_env)] if image_env else [
            Path("/home/datalake_Intermediate/pipeline/nextflow_containers/annotsv_3.5.10.sif"),
            Path("/home/pipeline/nextflow_containers/annotsv_3.5.10.sif"),
        ]
        image = next((p for p in images if p.is_file()), None)
        runtime = shutil.which("apptainer")
        if image is None or runtime is None:
            raise FileNotFoundError("AnnotSV runtime missing; set ANNOTSV_BIN or NGS_UI_ANNOTSV_SIF")
        command = [runtime, "exec", "--bind", f"{vcf.parent}:{vcf.parent}",
                   "--bind", f"{annotations}:{annotations}", "--bind", "/tmp", str(image), "AnnotSV"]
    return command + ["-SVinputFile", str(vcf), "-outputDir", str(output.parent),
                      "-outputFile", output.name, "-genomeBuild", "GRCh38",
                      "-annotationsDir", str(annotations), "-SVinputInfo", "1",
                      "-annotationMode", "both", "-SVminSize", "1", "-includeCI", "0"]


def build_rescue(*, joint_cnv: Path, base_tsv: Path | None = None,
                 post_dir: Path, sample_id: str, source_sample: str,
                 annotate: Callable[[Path, Path], None]) -> dict:
    """Called inside the worker's private staging tree; marker/promotion is external."""
    post_dir.mkdir(parents=True, exist_ok=True)
    outputs = {name: post_dir / f"{sample_id}.{name}" for name in MANAGED_NAMES}
    # A rerun, including a missing-input run, must not revive stale rescues.
    for path in outputs.values():
        path.unlink(missing_ok=True)
    manifest = {"policy_version": POLICY_VERSION, "sample_id": sample_id,
                "source_sample_id": source_sample, "pipeline": "dragen",
                "allowed_rule_b_filters": sorted(ALLOWED_FILTERS),
                "annotation_source": "integrated_cnv_only",
                "reciprocal_overlap": 0.5, "coordinate_basis": "integrated VCF (POS,END]",
                "status": "skipped", "records": []}
    missing = [str(joint_cnv)] if not joint_cnv.is_file() else []
    if missing:
        manifest["reason"] = "missing_input"
        manifest["missing"] = missing
    else:
        headers, joint_records = read_vcf(joint_cnv, source_sample)
        selected, counts = select_integrated_cnvs(joint_records, legacy_identity_ids(base_tsv))
        manifest.update(status="complete", counts=counts,
                        inputs={"cnv_sv": _signature(joint_cnv)},
                        records=[{"id": item["id"], **item["evidence"]} for item in selected],
                        annotated_events=0)
        if base_tsv is not None and base_tsv.is_file():
            manifest["inputs"]["base_tsv"] = {**_signature(base_tsv), "purpose": "legacy_ids_only"}
        if selected:
            with tempfile.TemporaryDirectory(prefix=".cnv-rescue-", dir=post_dir) as tmp:
                work = Path(tmp)
                vcf, annotation = work / "integrated_cnv.vcf", work / "annotated.tsv"
                with vcf.open("w", encoding="utf-8") as handle:
                    handle.write("\n".join(headers) + "\n")
                    for item in selected:
                        record = item["record"]
                        fields = list(record.fields)
                        fields[2] = item["id"]
                        # Explicit DEL/DUP for AnnotSV; input files remain immutable.
                        info = {**record.info, "SVTYPE": record.kind}
                        fields[7] = ";".join(f"{k}={v}" if v else k for k, v in info.items())
                        handle.write("\t".join(fields) + "\n")
                annotate(vcf, annotation)
                manifest.update(validate_annotations(annotation, selected, work / REVIEW_NAME, source_sample))
                os.replace(work / REVIEW_NAME, outputs[REVIEW_NAME])
        else:
            # An authoritative empty review must not fall back to old 06 events.
            _write_tsv(outputs[REVIEW_NAME], ["AnnotSV_ID", "SV_chrom", "SV_start", "SV_end",
                       "SV_type", "Annotation_mode", "FORMAT", source_sample, EVIDENCE_COLUMN], [])
    pending = outputs[MANIFEST_NAME].with_suffix(".json.tmp")
    pending.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    os.replace(pending, outputs[MANIFEST_NAME])
    return manifest
