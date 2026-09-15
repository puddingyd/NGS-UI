#!/usr/bin/env python3
"""DGX-side WES QC. Standalone Python 3.8+ / Samtools 1.13+; no UI imports.

Install at ${PIPELINE_CODE}/scripts/secondary_qc_report.py. Run --check-only
before Nextflow, then run without it after successful completion. See
docs/ops/SECONDARY_QC_REPORT.md for definitions and deployment instructions.
"""
from __future__ import annotations

import argparse
import ast
import bisect
import csv
import fcntl
import hashlib
import io
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from collections import Counter
from contextlib import contextmanager
from datetime import datetime, timezone
from fractions import Fraction
from functools import lru_cache
from pathlib import Path

VERSION = "1.1.0"
FIELDS = ["Sample ID", "Total reads", "Duplicated rate", "Mapping rate",
          "On target rate", "Mean depth", "Uniformity", "QC"]
PRIMARY_EXCLUDE = 0x100 | 0x800
# Depth follows the legacy Samtools defaults, including supplementary alignments.
# Total/mapped/on-target read counts retain their separate primary-only definition.
DEPTH_EXCLUDE = 0x4 | 0x100 | 0x200 | 0x400
DEPTH_MIN_MQ = 0
DEPTH_MIN_BQ = 0
METHOD = {
    "version": VERSION, "seq_type": "WES", "read_unit": "read end (R1/R2 counted separately)",
    "primary_exclude_flags": PRIMARY_EXCLUDE, "target_denominator": "mapped primary reads",
    "target_overlap": "at least one M/= /X reference base; each read counted once",
    "depth_exclude_flags": DEPTH_EXCLUDE,
    "min_mapping_quality": DEPTH_MIN_MQ, "min_base_quality": DEPTH_MIN_BQ,
    "overlap_removal": "none; both read ends contribute", "zero_depth_targets": "included",
    "uniformity": "bases with DP >= ceil(unrounded mean / 5) / all target bases",
    "thresholds": {"total_reads": 30000000, "mapping_rate": 0.95,
                   "on_target_rate": 0.40, "mean_depth": 50, "uniformity": 0.90},
}
SAMPLE_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]{0,127}$")


def now():
    return datetime.now(timezone.utc).isoformat()


def file_signature(path):
    p = Path(path).resolve(strict=True)
    st = p.stat()
    return {"path": str(p), "size": st.st_size, "mtime_ns": st.st_mtime_ns,
            "ctime_ns": st.st_ctime_ns}


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def atomic_text(path, text):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(prefix="." + path.name + ".", dir=str(path.parent))
    try:
        with os.fdopen(fd, "w", encoding="utf-8", newline="") as handle:
            handle.write(text)
            handle.flush()
            os.fsync(handle.fileno())
            os.fchmod(handle.fileno(), 0o664)
        os.replace(tmp, path)
    finally:
        if os.path.exists(tmp):
            os.unlink(tmp)


def atomic_json(path, value):
    atomic_text(path, json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n")


def samples_from_sheet(path):
    samples = []
    with Path(path).open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if "sample" not in (reader.fieldnames or []):
            raise ValueError("Samplesheet requires a sample column")
        for row in reader:
            sid = (row.get("sample") or "").strip()
            if not SAMPLE_RE.fullmatch(sid) or sid in {".", ".."}:
                raise ValueError("Invalid samplesheet sample ID: {!r}".format(sid))
            if sid not in samples:
                samples.append(sid)  # Multi-lane samples produce one row, in sheet order.
    if not samples:
        raise ValueError("Samplesheet has no samples")
    return samples


def resolved_settings(path, runtime_override=False):
    """Read literal paths from `nextflow config -flat`; never execute config text."""
    values = {}
    if path:
        for line in Path(path).read_text(encoding="utf-8").splitlines():
            key, sep, value = line.partition(" = ")
            if not sep:
                continue
            try:
                parsed = ast.literal_eval(value.strip())
            except (ValueError, SyntaxError):
                continue
            if isinstance(parsed, str) and "${" not in parsed:
                values[key.strip()] = parsed
    target = values.get("params.wes_targets")
    containers = {v for k, v in values.items()
                  if "SAMTOOLS" in k and k.endswith(".container")}
    if len(containers) > 1 and not runtime_override:
        raise ValueError("Multiple SAMTOOLS containers in config; specify --samtools-sif")
    return target, next(iter(containers), None) if len(containers) == 1 else None


class Targets:
    def __init__(self, path):
        intervals = {}
        with Path(path).open(encoding="utf-8-sig") as handle:
            for n, line in enumerate(handle, 1):
                if not line.strip() or line.startswith(("#", "track ", "browser ")):
                    continue
                cols = line.split()
                if len(cols) < 3:
                    raise ValueError("BED line {} needs at least 3 columns".format(n))
                chrom, start, end = cols[0], int(cols[1]), int(cols[2])
                if start < 0 or end <= start:
                    raise ValueError("Invalid BED interval at line {}".format(n))
                intervals.setdefault(chrom, []).append((start, end))
        self.regions = {}
        for chrom, pairs in intervals.items():
            merged = []
            for start, end in sorted(pairs):
                if merged and start <= merged[-1][1]:
                    merged[-1] = (merged[-1][0], max(end, merged[-1][1]))
                else:
                    merged.append((start, end))
            self.regions[chrom] = merged
        self.starts = {c: [s for s, _ in pairs] for c, pairs in self.regions.items()}
        self.ends = {c: [e for _, e in pairs] for c, pairs in self.regions.items()}
        self.length = sum(e - s for pairs in self.regions.values() for s, e in pairs)
        if not self.length:
            raise ValueError("Target BED is empty")

    def overlap(self, chrom, start, end):
        ends = self.ends.get(chrom, [])
        i = bisect.bisect_right(ends, start)
        return i < len(ends) and self.starts[chrom][i] < end

    def validate_header(self, header, sample):
        refs, sample_names = {}, set()
        sorted_bam = False
        for line in header.splitlines():
            fields = dict(c.split(":", 1) for c in line.split("\t")[1:] if ":" in c)
            if line.startswith("@HD\t"):
                sorted_bam = fields.get("SO") == "coordinate"
            elif line.startswith("@SQ\t"):
                refs[fields["SN"]] = int(fields["LN"])
            elif line.startswith("@RG\t") and fields.get("SM"):
                sample_names.add(fields["SM"])
        if not sorted_bam:
            raise ValueError("BAM must declare coordinate sort order")
        if sample_names and sample_names != {sample}:
            raise ValueError("BAM read-group sample does not match samplesheet: {}".format(sample_names))
        for chrom, regions in self.regions.items():
            if chrom not in refs or regions[-1][1] > refs[chrom]:
                raise ValueError("BED contig/coordinates do not match BAM: {}".format(chrom))

    def write(self, path):
        atomic_text(path, "".join("{}\t{}\t{}\n".format(c, s, e)
                                  for c in sorted(self.regions) for s, e in self.regions[c]))


class Samtools:
    def __init__(self, executable=None, sif=None, bind_dirs=()):
        self.image_signature = None
        if executable:
            self.prefix = [str(executable)]
        elif sif:
            image = Path(sif).resolve(strict=True)
            self.image_signature = file_signature(image)
            engine = shutil.which("apptainer") or shutil.which("singularity")
            if not engine:
                raise ValueError("Apptainer/Singularity is required for --samtools-sif")
            self.prefix = [engine, "exec", "--cleanenv"]
            roots = {str(Path(p).resolve()) for p in bind_dirs}
            roots.update(p for p in ("/datalake_Intermediate", "/datalake_Raw", "/raid") if Path(p).is_dir())
            for root in sorted(roots):
                if "," in root or ":" in root:
                    raise ValueError("Container bind path cannot contain ':' or ',': " + root)
                self.prefix.extend(["--bind", root + ":" + root])
            self.prefix.extend([str(image), "samtools"])
        else:
            self.prefix = [shutil.which("samtools") or "samtools"]
        self.version = self.capture(["--version"]).splitlines()[0]
        match = re.match(r"samtools (\d+)\.(\d+)", self.version)
        if not match or tuple(map(int, match.groups())) < (1, 13):
            raise ValueError("Samtools >=1.13 is required")

    def capture(self, args):
        proc = subprocess.run(self.prefix + list(map(str, args)), stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE, text=True)
        if proc.returncode:
            raise RuntimeError("Samtools {} failed ({}): {}".format(args[0], proc.returncode, proc.stderr[-4000:]))
        return proc.stdout

    @contextmanager
    def stream(self, args):
        # A temporary stderr file avoids pipe deadlock; large SAM/depth stdout is streamed.
        with tempfile.TemporaryFile(mode="w+t") as errors:
            proc = subprocess.Popen(self.prefix + list(map(str, args)), stdout=subprocess.PIPE,
                                    stderr=errors, text=True, bufsize=1024 * 1024)
            try:
                yield proc.stdout
                code = proc.wait()
                if code:
                    errors.seek(0)
                    raise RuntimeError("Samtools {} failed ({}): {}".format(args[0], code, errors.read()[-4000:]))
            finally:
                proc.stdout.close()
                if proc.poll() is None:
                    proc.terminate()
                    try:
                        proc.wait(timeout=10)
                    except subprocess.TimeoutExpired:
                        proc.kill()
                        proc.wait()


@lru_cache(maxsize=4096)
def aligned_blocks(cigar):
    """M/= /X consume actual aligned bases; D/N alone cannot make a read on-target."""
    offset, blocks, consumed = 0, [], 0
    for match in re.finditer(r"(\d+)([MIDNSHP=X])", cigar):
        size, op = int(match[1]), match[2]
        consumed += len(match[0])
        if op in "M=X":
            blocks.append((offset, offset + size))
        if op in "MDN=X":
            offset += size
    if consumed != len(cigar) or not blocks:
        if cigar == "*":
            return ()
        raise ValueError("Invalid mapped CIGAR: " + cigar)
    return tuple(blocks)


def count_target_reads(lines, targets):
    count = 0
    for line in lines:
        cols = line.split("\t", 6)
        if len(cols) < 6:
            raise ValueError("Malformed Samtools SAM output")
        chrom, start, cigar = cols[2], int(cols[3]) - 1, cols[5]
        if any(targets.overlap(chrom, start + s, start + e) for s, e in aligned_blocks(cigar)):
            count += 1
    return count


def depth_histogram(lines, targets):
    hist, last_position = Counter(), {}
    observed = 0
    for line in lines:
        chrom, position, value = line.rstrip().split("\t")
        pos, depth = int(position) - 1, int(value)
        if depth < 0 or pos <= last_position.get(chrom, -1) or not targets.overlap(chrom, pos, pos + 1):
            raise ValueError("Invalid, duplicate or out-of-target depth position: " + line.strip())
        last_position[chrom] = pos
        hist[depth] += 1
        observed += 1
    if observed > targets.length:
        raise ValueError("Depth output exceeds target territory")
    hist[0] += targets.length - observed  # Includes contigs/intervals with no alignments at all.
    return hist


def duplication_fraction(path):
    with Path(path).open(encoding="utf-8") as handle:
        rows = iter(handle)
        for line in rows:
            if line.startswith("LIBRARY\t"):
                fields = line.rstrip().split("\t")
                numerator = denominator = libraries = 0
                for row in rows:
                    if not row.strip() or row.startswith("#"):
                        break
                    values = dict(zip(fields, row.rstrip().split("\t")))
                    single = int(values["UNPAIRED_READS_EXAMINED"])
                    paired = int(values["READ_PAIRS_EXAMINED"])
                    dup_single = int(values["UNPAIRED_READ_DUPLICATES"])
                    dup_paired = int(values["READ_PAIR_DUPLICATES"])
                    if min(single, paired, dup_single, dup_paired) < 0 or dup_single > single or dup_paired > paired:
                        raise ValueError("Invalid duplication counts")
                    numerator += dup_single + 2 * dup_paired
                    denominator += single + 2 * paired
                    libraries += 1
                if libraries:
                    return numerator, denominator
    raise ValueError("Missing DuplicationMetrics table: " + str(path))


def percent(value):
    return "NA" if value is None else "{:.2f}%".format(float(value * 100))


def summarize(sample, total, mapped, on_target, dup_counts, hist, territory):
    if not 0 <= on_target <= mapped <= total or sum(hist.values()) != territory:
        raise ValueError("Inconsistent read/depth counts")
    depth_sum = sum(d * n for d, n in hist.items())
    mean = Fraction(depth_sum, territory)
    cutoff = (depth_sum + 5 * territory - 1) // (5 * territory)
    uniform_bases = sum(n for d, n in hist.items() if d >= cutoff) if depth_sum else 0
    mapping = Fraction(mapped, total) if total else None
    target = Fraction(on_target, mapped) if mapped else None
    uniformity = Fraction(uniform_bases, territory) if depth_sum else None
    dup_num, dup_den = dup_counts
    duplicate = Fraction(dup_num, dup_den) if dup_den else None
    checks = {
        "Total reads >= 30000000": total >= 30000000,
        "Mapping rate >= 95%": mapping is not None and mapping >= Fraction(95, 100),
        "On target rate >= 40%": target is not None and target >= Fraction(40, 100),
        "Mean depth >= 50X": mean >= 50,
        "Uniformity >= 90%": uniformity is not None and uniformity >= Fraction(90, 100),
    }
    row = dict(zip(FIELDS, [sample, total, percent(duplicate), percent(mapping), percent(target),
                            "{:.2f}".format(float(mean)), percent(uniformity),
                            "PASS" if all(checks.values()) else "FAIL"]))
    return {"row": row, "failed_checks": [k for k, v in checks.items() if not v],
            "counts": {"total_reads": total, "mapped_reads": mapped, "on_target_reads": on_target,
                       "duplicate_reads": dup_num, "duplicate_denominator": dup_den,
                       "target_bases": territory, "depth_sum": depth_sum,
                       "zero_depth_bases": hist[0], "uniformity_depth_cutoff": cutoff,
                       "uniformity_bases": uniform_bases},
            "unrounded": {"mean_depth": float(mean),
                          "uniformity": float(uniformity) if uniformity is not None else None}}


def analyze_sample(sample, out_dir, targets, merged_bed, tool, threads, method_signature, force=False):
    alignment = out_dir / sample / "02_alignment"
    bam = alignment / (sample + ".aligned.sorted.bam")
    dup = alignment / (sample + ".duplicate_metrics.txt")
    indexes = [Path(str(bam) + suffix) for suffix in (".bai", ".csi")] + [bam.with_suffix(".bai")]
    index = next((p for p in indexes if p.is_file()), None)
    if index is None:
        raise ValueError("Missing BAM index: " + str(bam))
    signature = {"method": method_signature, "bam": file_signature(bam),
                 "index": file_signature(index), "duplicates": file_signature(dup)}
    cache = out_dir / sample / "03_alignment_qc" / (sample + ".report_qc.json")
    if not force and cache.is_file():
        try:
            previous = json.loads(cache.read_text(encoding="utf-8"))
            if (isinstance(previous, dict) and previous.get("signature") == signature
                    and isinstance(previous.get("row"), dict) and set(previous["row"]) == set(FIELDS)
                    and isinstance(previous.get("failed_checks"), list)
                    and previous["row"].get("Sample ID") == sample
                    and previous["row"].get("QC") in {"PASS", "FAIL"}):
                return previous, True
        except (ValueError, OSError):
            pass
    tool.capture(["quickcheck", "-v", bam])
    targets.validate_header(tool.capture(["view", "-H", bam]), sample)
    stats = json.loads(tool.capture(["flagstat", "-@", threads, "-O", "json", bam]))
    total = sum(stats[k]["primary"] for k in ("QC-passed reads", "QC-failed reads"))
    mapped = sum(stats[k]["primary mapped"] for k in ("QC-passed reads", "QC-failed reads"))
    with tool.stream(["view", "-@", threads, "-F", PRIMARY_EXCLUDE | 0x4,
                      "-M", "-L", merged_bed, bam]) as lines:
        on_target = count_target_reads(lines, targets)
    with tool.stream(["depth", "-b", merged_bed, "-q", DEPTH_MIN_BQ, "-Q", DEPTH_MIN_MQ,
                      "-G", DEPTH_EXCLUDE, bam]) as lines:
        hist = depth_histogram(lines, targets)
    result = summarize(sample, total, mapped, on_target, duplication_fraction(dup), hist, targets.length)
    if any(file_signature(p) != signature[k] for k, p in (("bam", bam), ("index", index), ("duplicates", dup))):
        raise ValueError("BAM/index/metrics changed during QC; retry after analysis completes")
    result.update({"signature": signature, "computed_at": now(), "method": METHOD})
    atomic_json(cache, result)
    return result, False


def run(args):
    out_dir = args.out_dir.resolve()
    samples = samples_from_sheet(args.samplesheet)
    config_target, config_sif = resolved_settings(args.nextflow_config, bool(args.samtools or args.samtools_sif))
    target_path = args.target_bed or (Path(config_target) if config_target else None)
    if target_path is None or not target_path.is_absolute():
        raise ValueError("Provide absolute --target-bed or resolved params.wes_targets in --nextflow-config")
    target_path = target_path.resolve(strict=True)
    targets = Targets(target_path)
    sif = args.samtools_sif or config_sif
    if args.nextflow_config and not args.samtools and not sif:
        raise ValueError("No SAMTOOLS container in resolved config; supply --samtools-sif")
    info_dir = out_dir / "pipeline_info"
    info_dir.mkdir(parents=True, exist_ok=True)
    tool = Samtools(args.samtools, sif, [out_dir, target_path.parent])
    method_signature = {"script_sha256": sha256(__file__), "method": METHOD,
                        "target": file_signature(target_path), "target_sha256": sha256(target_path),
                        "samtools": tool.version, "samtools_command": tool.prefix,
                        "samtools_image": tool.image_signature}
    if args.check_only:
        print("[QC] Ready: {} samples; {} target bases; {}".format(len(samples), targets.length, tool.version))
        print("[QC] Target BED: " + str(target_path))
        return 0
    with (info_dir / ".report_summary.lock").open("a") as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise ValueError("Another QC summary is running for this batch")
        results = []
        details_path = info_dir / "report_summary.details.json"
        atomic_json(details_path, {"state": "RUNNING", "started_at": now(), "samples": samples})
        with tempfile.TemporaryDirectory(prefix=".report_qc_", dir=str(info_dir)) as tmp:
            merged_bed = Path(tmp) / "targets.bed"
            targets.write(merged_bed)
            for i, sample in enumerate(samples, 1):
                print("[QC] {}/{} {}".format(i, len(samples), sample), flush=True)
                try:
                    result, cached = analyze_sample(sample, out_dir, targets, merged_bed, tool,
                                                    args.threads, method_signature, args.force)
                    print("[QC] {}: {}{}{}".format(sample, result["row"]["QC"],
                          " (cached)" if cached else "",
                          "; " + ", ".join(result["failed_checks"]) if result["failed_checks"] else ""), flush=True)
                except (OSError, ValueError, RuntimeError, KeyError) as error:
                    result = {"row": {k: "NA" for k in FIELDS}, "error": str(error)}
                    result["row"].update({"Sample ID": sample, "QC": "ERROR"})
                    print("[QC] {}: ERROR: {}".format(sample, error), file=sys.stderr, flush=True)
                results.append(result)
        if sha256(target_path) != method_signature["target_sha256"]:
            atomic_json(details_path, {"state": "ERROR", "completed_at": now(),
                                      "error": "Target BED changed during QC; rerun the batch"})
            raise ValueError("Target BED changed during QC; rerun the batch")
        buffer = io.StringIO(newline="")
        writer = csv.DictWriter(buffer, fieldnames=FIELDS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(r["row"] for r in results)
        output = info_dir / "report_summary.csv"
        atomic_text(output, buffer.getvalue())
        states = Counter(r["row"]["QC"] for r in results)
        atomic_json(details_path, {"state": "ERROR" if states["ERROR"] else "COMPLETE",
                                  "completed_at": now(), "method": method_signature, "samples": results})
        print("[QC] {}: PASS={} FAIL={} ERROR={}".format(output, states["PASS"], states["FAIL"], states["ERROR"]), flush=True)
        return 2 if states["ERROR"] else 0  # A measured QC FAIL is a completed report, not a tool error.


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--samplesheet", type=Path, required=True)
    parser.add_argument("--target-bed", type=Path)
    parser.add_argument("--nextflow-config", type=Path, help="Resolved nextflow config -flat output")
    runtime = parser.add_mutually_exclusive_group()
    runtime.add_argument("--samtools", help="Use a native Samtools executable instead of a container")
    runtime.add_argument("--samtools-sif", type=Path)
    parser.add_argument("--threads", type=int, default=2, help="Samtools additional decompression threads")
    parser.add_argument("--check-only", action="store_true", help="Check deployment/BED/runtime before Nextflow; no BAMs needed")
    parser.add_argument("--force", action="store_true", help="Recalculate even when a verified sample cache exists")
    args = parser.parse_args(argv)
    if args.threads < 1 or args.threads > 16:
        parser.error("--threads must be between 1 and 16")
    try:
        return run(args)
    except (OSError, ValueError, RuntimeError) as error:
        print("[QC] ERROR: " + str(error), file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
