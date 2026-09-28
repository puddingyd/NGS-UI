"""BAM → targeted Mutect2 → full-germline subtraction → VEP → UI rows.

Run on the UI worker host (native binaries or configured container prefixes).
All execution uses argument arrays, with cooperative process-group cancellation.
"""
from __future__ import annotations

import csv
import gzip
import json
import os
import signal
import sqlite3
import subprocess
import sys
import time
import traceback
from urllib.parse import quote, unquote
from pathlib import Path

from .. import config
from ..services import somatic as store, sample_layout
from ..services.snv_rows import is_reportable_raw_row
from .dragen_run import _acquire_sample_locks, _release_sample_locks


class Cancelled(Exception):
    pass


def in_targets(chrom: str, pos: int, ref: str, intervals: list) -> bool:
    end = pos + len(ref) - 1
    return any(c == chrom and pos <= e and end >= s for c, s, e in intervals)


def variant_id(parts: list[str]) -> str:
    return "-".join(parts[i] for i in (0, 1, 3, 4))


def prepare_allele_filters(source: Path, output: Path, intervals: list | None = None):
    """GATK uses pipe-separated AS_FilterStatus; make it safely splittable.

    Keep the caller's original compressed VCF untouched. Encode commas within
    an allele's filter list before bcftools splits Number=A fields.
    """
    opener = gzip.open if source.suffix == ".gz" else open
    with opener(source, "rt") as inp, output.open("w") as out:
        for line in inp:
            if line.startswith("#CHROM"):
                out.write('##INFO=<ID=NGS_UI_AS_FILTER,Number=A,Type=String,Description="URL-encoded per-allele Mutect2 filters">\n')
                out.write('##INFO=<ID=NGS_UI_TARGET,Number=1,Type=Integer,Description="Overlaps requested intervals before normalization">\n')
            if line.startswith("#"):
                out.write(line)
                continue
            p = line.rstrip().split("\t")
            fields = [x for x in p[7].split(";") if x not in {".", ""}]
            if intervals is not None:
                fields.append("NGS_UI_TARGET=" + str(int(in_targets(p[0], int(p[1]), p[3], intervals))))
            allele_filter = next((x.split("=", 1)[1] for x in fields if x.startswith("AS_FilterStatus=")), "")
            if allele_filter:
                filters = allele_filter.split("|")
                count = len(p[4].split(","))
                if len(filters) != count:
                    # Some GATK versions use standard comma-separated Number=A.
                    filters = allele_filter.split(",")
                if len(filters) != count:
                    raise ValueError("無法對應 AS_FilterStatus 與 ALT，拒絕發布")
                fields = [x for x in fields if not x.startswith("AS_FilterStatus=")]
                fields.append("NGS_UI_AS_FILTER=" + ",".join(quote(x, safe="") for x in filters))
            p[7] = ";".join(fields) or "."
            out.write("\t".join(p) + "\n")


def subtract_vcf(source: Path, output: Path, db: sqlite3.Connection, intervals: list,
                 requested_ids: set[str] | None = None) -> dict:
    counts = {"germline_excluded": 0, "outside_targets": 0, "new_candidates": 0, "pass": 0}
    seen = set()
    with source.open() as inp, output.open("w") as out:
        for line in inp:
            if line.startswith("#"):
                out.write(line)
                continue
            p = line.rstrip().split("\t")
            info = dict(x.split("=", 1) for x in p[7].split(";") if "=" in x)
            allele_filters = unquote(info.get("NGS_UI_AS_FILTER", ""))
            if allele_filters not in {"", ".", "SITE", "PASS"}:
                filters = set(p[6].split(";")) - {".", "PASS"}
                filters.update(allele_filters.replace(",", ";").split(";"))
                p[6] = ";".join(sorted(filters))
            vid = variant_id(p)
            if db.execute("SELECT 1 FROM variants WHERE id=?", (vid,)).fetchone():
                counts["germline_excluded"] += 1
                continue
            if not (info.get("NGS_UI_TARGET") == "1" or vid in (requested_ids or set())
                    or in_targets(p[0], int(p[1]), p[3], intervals)):
                counts["outside_targets"] += 1
                continue
            if vid in seen:
                continue
            seen.add(vid)
            counts["new_candidates"] += 1
            counts["pass"] += p[6] == "PASS"
            out.write("\t".join(p) + "\n")
    return counts


def vep_rows(path: Path) -> list[dict]:
    """VEP JSON retains the original VCF input including Mutect2 FORMAT."""
    rows = []
    with path.open() as handle:
        for line in handle:
            obj = json.loads(line)
            p = obj["input"].split()
            if len(p) != 10:
                raise ValueError("Somatic VCF 必須只有一個檢體")
            fmt = dict(zip(p[8].split(":"), p[9].split(":")))
            info = dict(x.split("=", 1) for x in p[7].split(";") if "=" in x)
            base = {"CHROM": p[0], "POS": p[1], "REF": p[3], "ALT": p[4],
                    "RS_ID": ",".join(c["id"] for c in obj.get("colocated_variants", []) if str(c.get("id", "")).startswith("rs")),
                    "CALLERS": "Mutect2", "DP": fmt.get("DP", ""), "AD": fmt.get("AD", ""),
                    "VAF": fmt.get("AF", ""), "SOMATIC_FILTER": p[6],
                    "SOMATIC_TLOD": info.get("TLOD", ""), "SOMATIC_GT": fmt.get("GT", ""),
                    "SOMATIC_MMQ": info.get("MMQ", ""), "SOMATIC_MBQ": info.get("MBQ", ""),
                    "SOMATIC_MPOS": info.get("MPOS", ""), "SOMATIC_F1R2": fmt.get("F1R2", ""),
                    "SOMATIC_F2R1": fmt.get("F2R1", ""),
                    # Somatic GT=0/1 is not a constitutional heterozygosity assertion.
                    "ZYGOSITY": "", "ACMG_CRITERIA": "", "ACMG_CLASS": "",
                    "ACMG_SCORE": "", "STRAND_BIAS": "", "CLINVAR_SIG": "",
                    "CLINVAR_STARS": "", "CLINVAR_DN": "", "GNOMAD_G_AF": ""}
            afs = []
            for coloc in obj.get("colocated_variants", []):
                for value in coloc.get("frequencies", {}).get(p[4], {}).items():
                    if value[0] in {"gnomadg", "gnomadg_af"}:
                        afs.append(float(value[1]))
            if afs:
                base["GNOMAD_G_AF"] = str(max(afs))
            for tx in obj.get("transcript_consequences") or [{}]:
                rows.append(dict(base, GENE=tx.get("gene_symbol", ""), HGNC_ID=tx.get("hgnc_id", ""),
                                 TRANSCRIPT=tx.get("transcript_id", ""),
                                 TRANSCRIPT_TYPE="MANE_SELECT" if tx.get("mane_select") else "",
                                 HGVS_C=tx.get("hgvsc", "").split(":", 1)[-1],
                                 HGVS_P=tx.get("hgvsp", "").split(":", 1)[-1],
                                 REFSEQ_NUC=tx.get("mane_select", ""),
                                 EXON=tx.get("exon", ""), INTRON=tx.get("intron", ""),
                                 LOFTEE_HC=tx.get("lof", ""),
                                 MANE_STATUS="MANE_SELECT" if tx.get("mane_select") else "",
                                 CONSEQUENCE="&".join(tx.get("consequence_terms", [obj.get("most_severe_consequence", "")])),
                                 IMPACT=tx.get("impact", ""),
                                 SIFT=tx.get("sift_score", ""), POLYPHEN=tx.get("polyphen_score", "")))
    return rows


def write_rows(path: Path, rows: list[dict]) -> None:
    fields = list(dict.fromkeys(k for row in rows for k in row)) or [
        "CHROM", "POS", "REF", "ALT", "GENE", "ACMG_CRITERIA", "SOMATIC_FILTER"]
    with path.open("w") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


class Worker:
    def __init__(self, run_id: str):
        self.run_id = run_id
        self.directory = store.job_dir(run_id)
        self.job = store.read_json(self.directory / "state.json")
        self.cfg = store.read_json(self.directory / "config.json")
        self.stage = self.directory / "staging"
        self.stage.mkdir(exist_ok=True)

    def check_cancel(self):
        if (self.directory / "cancel").exists():
            raise Cancelled("使用者已取消分析")

    def update(self, step: str, **extra):
        self.job.update(step=step, updated=time.time(), **extra)
        store.atomic_json(self.directory / "state.json", self.job)
        print(step, flush=True)

    def run(self, args: list, output: Path | None = None):
        self.check_cancel()
        print("RUN", json.dumps([str(a) for a in args]), flush=True)
        stream = output.open("w") if output else None
        try:
            proc = subprocess.Popen([str(a) for a in args], stdout=stream, start_new_session=True)
            try:
                while proc.poll() is None:
                    self.check_cancel()
                    time.sleep(0.25)
                if proc.returncode:
                    raise RuntimeError(f"{args[0]} 執行失敗 (exit {proc.returncode})，請查看 log")
            except BaseException:
                if proc.poll() is None:
                    os.killpg(proc.pid, signal.SIGTERM)
                    try:
                        proc.wait(timeout=5)
                    except subprocess.TimeoutExpired:
                        os.killpg(proc.pid, signal.SIGKILL)
                        proc.wait()
                raise
        finally:
            if stream:
                stream.close()

    def execute(self):
        sid = self.job["sample_id"]
        cfg, stage = self.cfg, self.stage
        ref, bam = cfg["reference"], self.job["bam_path"]
        gatk, bcftools, samtools, vep = [cfg[t + "_command"] for t in ("gatk", "bcftools", "samtools", "vep")]
        self.update("preflight", status="running", pid=os.getpid())
        if not sample_layout.state_file(sid, "sample_metadata.json").is_file():
            raise ValueError("個案已取消登錄")
        raw = sample_layout.snv_raw_tsv(sid)
        if store.signature(raw) != self.job["raw_signature"] or store.signature(Path(bam)) != self.job["bam_signature"]:
            raise ValueError("輸入資料已變更，請重新送出分析")
        for tool, command in (("gatk", gatk), ("bcftools", bcftools), ("samtools", samtools), ("vep", vep)):
            self.run(command + ["--help" if tool == "vep" else "--version"], stage / f"{tool}.version.txt")
        self.run(samtools + ["quickcheck", bam])
        self.run(samtools + ["view", "-H", bam], stage / "bam.header.sam")
        expected = store.contigs(cfg)
        sq, samples = {}, set()
        for line in (stage / "bam.header.sam").read_text().splitlines():
            fields = dict(p.split(":", 1) for p in line.split("\t")[1:] if ":" in p)
            if line.startswith("@SQ"):
                sq[fields["SN"]] = int(fields["LN"])
            elif line.startswith("@RG") and fields.get("SM"):
                samples.add(fields["SM"])
        if len(samples) != 1:
            raise ValueError("BAM 必須有且只有一個 RG sample name")
        for chrom, length in expected.items():
            if chrom in {"chr" + str(i) for i in range(1, 23)} | {"chrX", "chrY"} and sq.get(chrom) != length:
                raise ValueError(f"BAM 與 GRCh38 reference 不相容：{chrom}")
        self.job["bam_sample"] = next(iter(samples))
        self.job["resource_signatures"] = {key: store.signature(Path(cfg[key])) for key in
            ("reference", "gene_regions", "germline_resource", "clinvar_vcf", "pon", "contamination_sites") if cfg.get(key)}
        target_bed = stage / "targets.bed"
        with target_bed.open("w") as handle:
            for chrom, start, end in self.job["targets"]["intervals"]:
                handle.write(f"{chrom}\t{start - 1}\t{end}\n")
        alleles = [p for p in self.job["targets"]["positions"] if p["ref"]]
        force_args = []
        requested_ids = set()
        if alleles:
            force = stage / "requested.vcf"
            with force.open("w") as out:
                out.write("##fileformat=VCFv4.2\n")
                for chrom, length in expected.items():
                    out.write(f"##contig=<ID={chrom},length={length}>\n")
                out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
                for p in alleles:
                    out.write(f'{p["chrom"]}\t{p["start"]}\t.\t{p["ref"]}\t{p["alt"]}\t.\t.\t.\n')
            self.run(bcftools + ["norm", "-f", ref, "-c", "e", "-Ov", force], stage / "requested.norm.vcf")
            with (stage / "requested.norm.vcf").open() as handle:
                requested_ids = {variant_id(line.rstrip().split("\t")) for line in handle if not line.startswith("#")}
            self.run(bcftools + ["sort", "-Oz", "-o", stage / "requested.vcf.gz", stage / "requested.norm.vcf"])
            self.run(bcftools + ["index", "-t", stage / "requested.vcf.gz"])
            force_args = ["--alleles", stage / "requested.vcf.gz"]
        self.update("mutect2")
        raw_vcf = stage / "mutect2.vcf.gz"
        cmd = gatk + ["Mutect2", "-R", ref, "-I", bam, "-L", target_bed, "--interval-padding", "100",
                      "--germline-resource", cfg["germline_resource"],
                      "--f1r2-tar-gz", stage / "f1r2.tar.gz", "-O", raw_vcf] + force_args
        if cfg.get("pon"):
            cmd += ["--panel-of-normals", cfg["pon"]]
        self.run(cmd)
        self.update("filtering")
        orientation = stage / "orientation.tar.gz"
        self.run(gatk + ["LearnReadOrientationModel", "-I", stage / "f1r2.tar.gz", "-O", orientation])
        contamination_args = []
        warnings = []
        if cfg.get("contamination_sites"):
            self.run(gatk + ["GetPileupSummaries", "-I", bam, "-V", cfg["contamination_sites"],
                             "-L", cfg["contamination_sites"], "-O", stage / "pileups.table"])
            self.run(gatk + ["CalculateContamination", "-I", stage / "pileups.table",
                             "-O", stage / "contamination.table"])
            contamination_args = ["--contamination-table", stage / "contamination.table"]
        else:
            warnings.append("未設定 contamination_sites；本次未估計污染比例")
        if not cfg.get("pon"):
            warnings.append("未設定相容的 panel of normals")
        filtered = stage / "filtered.vcf.gz"
        self.run(gatk + ["FilterMutectCalls", "-R", ref, "-V", raw_vcf,
                         "--stats", str(raw_vcf) + ".stats", "--ob-priors", orientation,
                         "-O", filtered] + contamination_args)
        self.update("subtract-germline")
        normalized = stage / "filtered.normalized.vcf"
        prepared = stage / "filtered.allele_filters.vcf"
        prepare_allele_filters(filtered, prepared, self.job["targets"]["intervals"])
        self.run(bcftools + ["norm", "-f", ref, "-c", "e", "-m", "-any", "-Ov", prepared], normalized)
        germline = stage / ".germline.vcf"
        with germline.open("w") as out, raw.open() as inp:
            out.write("##fileformat=VCFv4.2\n")
            for chrom, length in expected.items():
                out.write(f"##contig=<ID={chrom},length={length}>\n")
            out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            for n, row in enumerate(csv.DictReader(inp, delimiter="\t")):
                if n % 10000 == 0:
                    self.check_cancel()
                if not is_reportable_raw_row(row):
                    continue
                chrom = "chr" + re_strip_chr(row["CHROM"])
                if chrom not in expected or chrom in {"chrM", "chrMT"}:
                    continue
                if not row["REF"] or any(a not in "ACGTNacgtn," for a in row["REF"] + row["ALT"]):
                    continue
                out.write(f'{chrom}\t{row["POS"]}\t.\t{row["REF"]}\t{row["ALT"]}\t.\tPASS\t.\n')
        norm_germline = stage / ".germline.normalized.vcf"
        self.run(bcftools + ["norm", "-f", ref, "-c", "e", "-m", "-any", "-Ov", germline], norm_germline)
        # Only candidate keys are needed in memory, even for whole-genome raw TSVs.
        with normalized.open() as handle:
            candidate_ids = {variant_id(line.rstrip().split("\t")) for line in handle if not line.startswith("#")}
        with sqlite3.connect(stage / ".germline.sqlite") as db:
            db.execute("CREATE TABLE variants (id TEXT PRIMARY KEY)")
            with norm_germline.open() as handle:
                for n, line in enumerate(handle):
                    if n % 10000 == 0:
                        self.check_cancel()
                    if not line.startswith("#"):
                        vid = variant_id(line.rstrip().split("\t"))
                        if vid in candidate_ids:
                            db.execute("INSERT OR IGNORE INTO variants VALUES (?)", (vid,))
            counts = subtract_vcf(normalized, stage / "novel.vcf", db, self.job["targets"]["intervals"], requested_ids)
        germline.unlink()
        norm_germline.unlink()
        (stage / ".germline.sqlite").unlink()
        self.update("annotation", counts=counts)
        annotation = stage / "annotations.tsv"
        if counts["new_candidates"]:
            self.run(vep + ["--offline", "--cache", "--dir_cache", cfg["vep_cache"], "--assembly", "GRCh38",
                            "--cache_version", str(cfg["vep_cache_version"]),
                            "--fasta", ref, "--format", "vcf", "--json", "--everything", "--no_stats",
                            "--force_overwrite", "--input_file", stage / "novel.vcf",
                            "--output_file", stage / "vep.json"])
            rows = vep_rows(stage / "vep.json")
            ids = {f'{r["CHROM"]}-{r["POS"]}-{r["REF"]}-{r["ALT"]}' for r in rows}
            if len(ids) != counts["new_candidates"]:
                raise ValueError("VEP 未完整保留所有新增點位，拒絕發布")
            write_rows(annotation, rows)
            self.run([sys.executable, config.REPO_ROOT / "scripts/annotate_clinvar.py",
                      "--tsv", annotation, "--clinvar", cfg["clinvar_vcf"]])
            if config.GENEBE_DB.is_file():
                self.run([sys.executable, config.REPO_ROOT / "scripts/annotate_acmg_genebe.py",
                          "--tsv", annotation, "--genebe-db", config.GENEBE_DB, "--skip-api", "--test-type", "WGS"])
            else:
                warnings.append("GeneBe 本地資料庫不存在；ACMG 保留未分類")
        else:
            write_rows(annotation, [])
        self.update("coverage")
        # Quality-filtered, non-overlapping read depth; no LOD claim is derived from it.
        depth_file = stage / "depth.tsv"
        self.run(samtools + ["depth", "-aa", "-s", "-q", "20", "-Q", "20", "-b", target_bed, bam], depth_file)
        positions = [dict(p, covered_bases=0, depth_sum=0, min_depth=None, max_depth=0) for p in self.job["targets"]["positions"]]
        bases = covered = depth_sum = 0
        with depth_file.open() as handle:
            for line in handle:
                c, p, d = line.split()[:3]
                p, d = int(p), int(d)
                bases += 1
                covered += d > 0
                depth_sum += d
                for pos in positions:
                    if pos["chrom"] == c and pos["start"] <= p <= pos["end"]:
                        pos["covered_bases"] += d > 0
                        pos["depth_sum"] += d
                        pos["min_depth"] = min(pos["min_depth"], d) if pos["min_depth"] is not None else d
                        pos["max_depth"] = max(pos["max_depth"], d)
        for pos in positions:
            pos["mean_depth"] = pos.pop("depth_sum") / (pos["end"] - pos["start"] + 1)
            pos["min_depth"] = pos["min_depth"] or 0
        depth_file.unlink()
        coverage = {"target_bases": self.job["targets"]["total_bases"], "covered_bases": covered,
                    "mean_depth": depth_sum / self.job["targets"]["total_bases"], "positions": positions,
                    "method": "samtools depth -aa -s -q 20 -Q 20; 未檢出不等於排除低比例變異"}
        store.atomic_json(stage / "coverage.json", coverage)
        candidates = {}
        with annotation.open() as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                vid = f'{row["CHROM"]}-{row["POS"]}-{row["REF"]}-{row["ALT"]}'
                candidates.setdefault(vid, {"id": vid, "gene": row.get("GENE", ""), "filter": row["SOMATIC_FILTER"],
                                           "dp": row.get("DP", ""), "ad": row.get("AD", ""), "vaf": row.get("VAF", ""),
                                           "qc": {key: row.get("SOMATIC_" + key, "") for key in
                                                  ("TLOD", "MMQ", "MBQ", "MPOS", "F1R2", "F2R1")}})
        store.atomic_json(stage / "candidates.json", list(candidates.values()))
        self.check_cancel()
        if store.signature(raw) != self.job["raw_signature"] or store.signature(Path(bam)) != self.job["bam_signature"]:
            raise ValueError("分析期間輸入已變更；拒絕發布")
        for key, sig in self.job["resource_signatures"].items():
            if store.signature(Path(cfg[key])) != sig:
                raise ValueError(f"分析期間 reference 資源已變更：{key}")
        if not sample_layout.state_file(sid, "sample_metadata.json").is_file():
            raise ValueError("個案已取消登錄；拒絕發布")
        record = {"run_id": self.run_id, "created": self.job["created"], "raw_signature": self.job["raw_signature"],
                  "bam_path": bam, "clinvar_release": cfg["clinvar_release"], "counts": counts}
        store.atomic_json(stage / "manifest.json", dict(self.job, warnings=warnings, coverage=coverage,
                                                        config=cfg, clinvar_release=cfg["clinvar_release"]))
        self.update("publishing")
        target = store.result_dir(sid, self.run_id)
        target.parent.mkdir(parents=True, exist_ok=True)
        # Copy staging to a hidden directory on the destination filesystem first.
        import shutil
        temporary = target.with_name("." + target.name)
        shutil.copytree(stage, temporary)
        self.check_cancel()
        os.replace(temporary, target)
        with store.submission_lock():
            self.check_cancel()
            data = store.manifest(sid)
            data.setdefault("runs", [])
            for prior in data["runs"]:
                if prior["run_id"] == self.job.get("replace_run_id"):
                    prior["archived"] = True
            data["runs"].append(record)
            store.atomic_json(store.index_path(sid), data)
            self.update("completed", status="completed", counts=counts, warnings=warnings, coverage=coverage)


def re_strip_chr(chrom: str) -> str:
    return chrom[3:] if chrom.lower().startswith("chr") else chrom


def main(run_id: str) -> int:
    worker = Worker(run_id)
    def interrupted(signum, frame):
        raise Cancelled("Somatic worker 收到停止訊號")
    signal.signal(signal.SIGTERM, interrupted)
    signal.signal(signal.SIGINT, interrupted)
    locks = []
    try:
        locks = _acquire_sample_locks([worker.job["sample_id"]])
        worker.execute()
        return 0
    except Cancelled as exc:
        worker.update("cancelled", status="cancelled", error=str(exc))
        return 1
    except Exception as exc:
        traceback.print_exc()
        worker.update("failed", status="failed", error=str(exc))
        return 1
    finally:
        import shutil
        _release_sample_locks(locks)
        shutil.rmtree(worker.stage, ignore_errors=True)
        target = store.result_dir(worker.job["sample_id"], run_id)
        shutil.rmtree(target.with_name("." + target.name), ignore_errors=True)


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1]))
