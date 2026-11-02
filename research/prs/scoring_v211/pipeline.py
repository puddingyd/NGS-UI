#!/usr/bin/env python3
"""PGS004879 integrated runner, v1.0; Python stdlib only.

Never mutates legacy runs. No background execution. Serial foreground by design.
Raw stage delegates all genotype science to the validated FINAL gVCF script.
"""
import argparse
import contextlib
import csv
import datetime as dt
import fcntl
import gzip
import hashlib
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
SID = re.compile(r"^[A-Za-z0-9_.-]+$")


def failure(message):
    raise RuntimeError(message)


def now():
    return dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")


def read_json(p):
    return json.loads(Path(p).read_text())


def atomic_json(path, obj):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(json.dumps(obj, indent=2, ensure_ascii=False) + "\n")
    os.replace(tmp, path)


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1024*1024), b""):
            h.update(b)
    return h.hexdigest()


def stat_stamp(path):
    p = Path(path).resolve()
    if not p.is_file():
        failure(f"Missing file: {p}")
    st = p.stat()
    return {"path": str(p), "size": st.st_size, "mtime_ns": st.st_mtime_ns}


def cfg_read(path):
    p = Path(path).resolve()
    c = read_json(p)
    if not isinstance(c, dict): failure("Configuration must be a JSON object")
    for k in ("project", "raw_script", "raw_work", "original_scorefile", "fasta", "dbsnp", "plink_img", "raw_nextflow_config", "ancestry_nextflow_config", "reference_panel"):
        if k not in c or not c[k]: failure(f"Missing config field: {k}")
    c = {k: (str(Path(v).expanduser().resolve()) if k in ("project", "raw_script", "raw_work", "original_scorefile", "fasta", "dbsnp", "plink_img", "raw_nextflow_config", "ancestry_nextflow_config", "reference_panel") else v) for k, v in c.items()}
    c.setdefault("min_dp",10);c.setdefault("min_gq",20);c.setdefault("min_overlap",0.75)
    c.setdefault("palindromic_mode","assume-direct")
    c.setdefault("raw_threads",2);c.setdefault("raw_memory_mb",8192)
    c.setdefault("max_cpus",8);c.setdefault("max_memory","96.GB")
    c.setdefault("reference_label","SuperPop")
    c.setdefault("projection_method","oadp")
    c.setdefault("pca_maf_target",0);c.setdefault("pca_geno_miss_target",0.1)
    c.setdefault("boundary_review_pp",2.0)
    c.setdefault("high_weight_abs",0.001)
    c.setdefault("min_match_fraction",0.995)
    return c


def check_paths(c, mode, gvcf, raw_source, sample):
    for k in ("raw_script", "original_scorefile", "fasta", "dbsnp", "plink_img", "raw_nextflow_config", "ancestry_nextflow_config", "reference_panel"):
        if not Path(c[k]).is_file(): failure(f"Missing configured file {k}: {c[k]}")
    for k in ("fasta",):
        if not Path(c[k] + ".fai").is_file(): failure(f"Missing FASTA index: {c[k]}.fai")
    for k in ("dbsnp",):
        if not (Path(c[k] + ".tbi").exists() or Path(c[k] + ".csi").exists()):
            failure(f"Missing dbSNP index for {c[k]}")
    for tool in ("nextflow", "apptainer", "bcftools", "samtools"):
        if not shutil.which(tool): failure(f"Tool not in PATH: {tool}")
    with gzip.open(c["original_scorefile"], "rt") as f:
        rows=sum(1 for line in f if line.strip() and not line.startswith("#"))
    if rows != 610678:failure(f"Original scorefile must have header + 610677 loci; observed {rows} non-comment rows")
    if mode == "fresh":
        if not gvcf: failure("--gvcf is required with --mode fresh")
        if not Path(gvcf).is_file(): failure(f"Missing gVCF: {gvcf}")
        if not (Path(gvcf + ".tbi").exists() or Path(gvcf + ".csi").exists()):
            failure(f"Missing gVCF index: {gvcf}")
        found=subprocess.check_output(["bcftools","query","-l",gvcf],text=True).strip().splitlines()
        if found != [sample]:failure(f"gVCF sample header {found!r} differs from requested {sample}")
    if mode == "reuse" and not raw_source: failure("--raw-source required with --mode reuse")


def fingerprint(c, sample, mode, gvcf, raw_source):
    # Hash small configs/script/PGS; preserve stat identity for huge reference/gVCF.
    paths = ["raw_script", "original_scorefile", "raw_nextflow_config", "ancestry_nextflow_config"]
    d = {"sample":sample,"mode":mode,"parameters":{k:v for k,v in c.items() if k not in ("project",)},
         "sha256":{k:sha(c[k]) for k in paths},
         "source_metadata": {k:stat_stamp(c[k]) for k in ("fasta","dbsnp","plink_img","reference_panel")}}
    if gvcf: d["gvcf_metadata"] = stat_stamp(gvcf)
    if raw_source: d["raw_source"] = str(Path(raw_source).resolve())
    return d


@contextlib.contextmanager
def sample_lock(project, sample):
    lockdir = Path(project) / "integrated_runs" / ".locks"
    lockdir.mkdir(parents=True, exist_ok=True)
    fp = open(lockdir / f"{sample}.lock", "a+")
    try:
        try:
            fcntl.flock(fp, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            failure(f"Sample {sample} already has an integrated run in progress")
        fp.seek(0);fp.truncate(0)
        fp.write(f"pid={os.getpid()}\nstart={now()}\n");fp.flush()
        yield
    finally:
        fp.close()


def log_run(cmd, logfile, cwd=None):
    logfile = Path(logfile)
    logfile.parent.mkdir(parents=True, exist_ok=True)
    print("COMMAND:", shlex.join(str(z) for z in cmd), flush=True)
    with logfile.open("a") as out:
        out.write(f"\n===== {now()} =====\n{shlex.join(str(z) for z in cmd)}\n")
        out.flush()
        with subprocess.Popen([str(x) for x in cmd], cwd=cwd,
                              stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                              text=True, errors="replace", bufsize=1) as p:
            for line in p.stdout:
                print(line, end="", flush=True)
                out.write(line)
            ret = p.wait()
    if ret: failure(f"Command returned {ret}; see {logfile}")


def raw_integrity(rawdir, sample):
    rawdir = Path(rawdir)
    f = rawdir / "audit/final_reconciliation.json"
    if not f.is_file():failure(f"Raw run incomplete (missing reconciliation): {f}")
    d = read_json(f)
    if str(d.get("iid")) != sample or str(d.get("status")) != "PASS":
        failure(f"Raw reconciliation sample/status mismatch: {f}")
    corrected = rawdir / "scorefiles/PGS004879_hmPOS_GRCh38.txt.gz"
    if not corrected.is_file(): failure(f"Missing raw-stage corrected scorefile: {corrected}")
    sheet = rawdir / "samplesheet.csv"
    if not sheet.is_file():failure(f"Missing raw-stage samplesheet: {sheet}")
    with sheet.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:failure("Empty raw-stage samplesheet")
    if any(r.get("sampleset") != f"case{sample}" for r in rows):
        failure(f"Raw sampleset mismatch in {sheet}")
    for r in rows:
        prefix = r.get("path_prefix", "")
        if not prefix or r.get("format") != "pfile":failure(f"Invalid row in raw samplesheet: {r}")
        for ext in (".pgen", ".psam", ".pvar"):
            if not Path(prefix + ext).is_file():failure(f"Missing raw-stage genotype: {prefix}{ext}")
    if not (rawdir / f"results/case{sample}/score/aggregated_scores.txt.gz").is_file():
        failure(f"Missing original raw PLINK score: {rawdir}")
    return d


def build_samplesheet(rawdir, ancestry, sample):
    source = Path(rawdir) / "samplesheet.csv"
    dest = Path(ancestry) / "inputs/samplesheet.csv"
    dest.parent.mkdir(parents=True, exist_ok=True)
    # Keep the same absolute path prefixes, do not copy multi-GB PGEN data.
    rows = list(csv.DictReader(source.open(newline="")))
    if any(not Path(r["path_prefix"] + ".pgen").is_file() for r in rows):
        failure("Raw samplesheet genotype paths are stale")
    with dest.open("w", newline="") as o:
        writer=csv.DictWriter(o,fieldnames=["sampleset","path_prefix","chrom","format"])
        writer.writeheader(); writer.writerows(rows)
    return dest


def cmd_run(a):
    c=cfg_read(a.config)
    sample=a.sample
    if (not SID.fullmatch(sample) or sample in (".", "..") or
        not SID.fullmatch(a.run_id) or a.run_id in (".", "..")):
        failure("Unsafe sample or run ID")
    if a.mode == "fresh" and a.raw_source:failure("--raw-source applies only to reuse mode")
    raw_source = a.raw_source or (str(Path(c["project"])/"results"/f"{sample}_PGS004879_FINAL_gVCF_to_raw_score") if a.mode=="reuse" else None)
    gvcf = str(Path(a.gvcf).resolve()) if a.gvcf else None
    check_paths(c, a.mode, gvcf, raw_source, sample)
    dest=Path(c["project"])/"integrated_runs"/a.run_id/sample
    ancestry=dest/"ancestry";qc=dest/"qc"; rawdest=dest/"raw"
    with sample_lock(c["project"], sample):
        dest.mkdir(parents=True,exist_ok=True)
        fingerprint_file=dest/"run_inputs.json"
        fp=fingerprint(c,sample,a.mode,gvcf,raw_source)
        if fingerprint_file.is_file():
            old=read_json(fingerprint_file)
            if old != fp:failure(f"Run inputs changed, use a NEW --run-id: {fingerprint_file}")
        else:atomic_json(fingerprint_file,fp)
        status=dest/"status.json"
        atomic_json(status,{"stage":"STARTED","updated_at":now()})
        try:
            if a.mode == "fresh":
                rawdir=rawdest
                try:
                    raw_integrity(rawdir, sample)
                    print("RAW_ALREADY_VALID; reusing completed raw stage", flush=True)
                except RuntimeError:
                    if rawdir.exists():
                        failure("Raw directory exists but is incomplete. The validated raw script refuses overwrite; use a new run ID (do not delete original audit).")
                    rawcmd=["bash",c["raw_script"],"--work",c["raw_work"],"--gvcf",gvcf,
                            "--out",str(rawdir),"--scorefile",c["original_scorefile"],"--ref",c["fasta"],
                            "--dbsnp",c["dbsnp"],"--plink-img",c["plink_img"],"--container-config",c["raw_nextflow_config"],
                            "--min-dp",str(c["min_dp"]),"--min-gq",str(c["min_gq"]),
                            "--min-overlap",str(c["min_overlap"]),"--palindromic-mode",str(c["palindromic_mode"]),
                            "--threads",str(c["raw_threads"]),"--memory-mb",str(c["raw_memory_mb"])]
                    log_run(rawcmd, dest/"logs/raw_stage.log")
            else:rawdir=Path(raw_source).resolve()
            raw=raw_integrity(rawdir,sample)
            # Guard against editing/replacing a completed raw reference after ancestry was run.
            raw_manifest={
                "raw_dir":str(rawdir.resolve()),
                "sha256_reconciliation":sha(rawdir/"audit/final_reconciliation.json"),
                "sha256_corrected_scorefile":sha(rawdir/"scorefiles/PGS004879_hmPOS_GRCh38.txt.gz"),
                "sha256_samplesheet":sha(rawdir/"samplesheet.csv"),
                "pgen_metadata":[stat_stamp(row["path_prefix"]+".pgen") for row in csv.DictReader((rawdir/"samplesheet.csv").open(newline=""))]
            }
            raw_manifest_path=dest/"raw_inputs.json"
            if raw_manifest_path.is_file() and read_json(raw_manifest_path)!=raw_manifest:
                failure(f"Completed raw inputs have changed; use NEW --run-id: {raw_manifest_path}")
            if not raw_manifest_path.is_file():atomic_json(raw_manifest_path,raw_manifest)
            atomic_json(status,{"stage":"RAW_VALID","updated_at":now(),"raw_dir":str(rawdir)})
            inp=build_samplesheet(rawdir,ancestry,sample)
            correction=str(rawdir/"scorefiles/PGS004879_hmPOS_GRCh38.txt.gz")
            normdir=ancestry/"results"/f"case{sample}"
            pgsout=normdir/"score"/f"case{sample}_pgs.txt.gz"
            popout=normdir/"score"/f"case{sample}_popsimilarity.txt.gz"
            matchout=normdir/"match"/f"case{sample}_summary.csv"
            if pgsout.is_file() and popout.is_file() and matchout.is_file() and not a.recheck_ancestry:
                print("ANCESTRY_OUTPUTS_EXIST; validating outputs in QC report",flush=True)
            else:
                cmd=["nextflow","run","pgscatalog/pgsc_calc","-r","v2.3.0","-profile","singularity",
                     "-c",c["ancestry_nextflow_config"],"-resume","-w",str(ancestry/"work"),
                     "--input",str(inp),"--target_build","GRCh38","--scorefile",correction,
                     "--run_ancestry",c["reference_panel"],"--projection_method",c["projection_method"],
                     "--ref_label",c["reference_label"],"--pca_maf_target",str(c["pca_maf_target"]),
                     "--pca_geno_miss_target",str(c["pca_geno_miss_target"]),"--max_cpus",str(c["max_cpus"]),
                     "--max_memory",c["max_memory"],"--outdir",str(ancestry/"results")]
                log_run(cmd,dest/"logs/ancestry_stage.log",cwd=ancestry)
            atomic_json(status,{"stage":"ANCESTRY_FINISHED","updated_at":now(),"raw_dir":str(rawdir)})
            qcmd=[sys.executable,str(HERE/"qc_report.py"),"--sample",sample,"--raw",str(rawdir),
                  "--ancestry",str(ancestry),"--out",str(qc),"--score-sha256",sha(correction),
                  "--ref-panel",c["reference_panel"],"--min-match",str(c["min_match_fraction"]),
                  "--high-weight",str(c["high_weight_abs"]),"--boundary-pp",str(c["boundary_review_pp"])]
            log_run(qcmd,dest/"logs/qc_stage.log")
            result=read_json(qc/"qc.json")
            if result["qc_status"] == "FAIL":failure(f"QC hard gates failed; see {qc/'report.md'}")
            atomic_json(status,{"stage":"COMPLETE","updated_at":now(),"qc_status":result["qc_status"],
                                "raw_dir":str(rawdir),"qc_report":str(qc/"report.md")})
            print("COMPLETED:", sample, "SUM=", result["normalized_sum"], "PCT=", result["percentile"],
                  "QC=", result["qc_status"],flush=True)
        except Exception as exc:
            atomic_json(status,{"stage":"FAILED","updated_at":now(),"message":str(exc)})
            raise


def cmd_regression(a):
    c=cfg_read(a.config)
    manifest=Path(a.manifest)
    if not manifest.is_file():failure(f"Missing manifest: {manifest}")
    rows=list(csv.DictReader(manifest.open(newline=""),delimiter="\t"))
    if len(rows)!=3 or {r.get("sample") for r in rows}!={"26T00970","26T00969","26T00968"}:
        failure("Regression manifest must contain exactly the three validated sample IDs")
    results=[]
    for row in rows:
        arg=argparse.Namespace(config=a.config,sample=row["sample"],run_id=a.run_id,
          mode=a.mode,gvcf=row.get("gvcf","") or None,
          raw_source=row.get("raw_source","") or None,
          recheck_ancestry=a.recheck_ancestry)
        try:
            cmd_run(arg)
            dest=Path(c["project"])/"integrated_runs"/a.run_id/arg.sample
            results.append((arg.sample,"PASS",str(dest/"qc/qc.json")))
        except Exception as e:
            results.append((arg.sample,"FAIL",str(e)))
            if not a.keep_going:break
    report_root=Path(c["project"])/"integrated_runs"/a.run_id
    report_root.mkdir(parents=True,exist_ok=True)
    with (report_root/"regression_execution.tsv").open("w",newline="") as f:
        wr=csv.writer(f,delimiter="\t");wr.writerow(["sample","execution_status","qc_or_error"]);wr.writerows(results)
    if len(results)!=3 or any(x[1] != "PASS" for x in results):
        failure(f"Regression execution incomplete: {results}")
    command=[sys.executable,str(HERE/"compare_benchmarks.py"),"--run-dir",str(report_root),
             "--expected",str(Path(a.expected).resolve())]
    log_run(command,report_root/"regression_comparison.log")


def cli():
    parser=argparse.ArgumentParser(description=__doc__)
    sub=parser.add_subparsers(dest="command",required=True)
    run=sub.add_parser("run",help="Run one sample; fresh gVCF or previously completed raw score")
    run.add_argument("--config",required=True)
    run.add_argument("--sample",required=True)
    run.add_argument("--run-id",required=True)
    run.add_argument("--mode",choices=["fresh","reuse"],required=True)
    run.add_argument("--gvcf")
    run.add_argument("--raw-source")
    run.add_argument("--recheck-ancestry",action="store_true",help="rerun nextflow -resume even if outputs already exist")
    reg=sub.add_parser("regression",help="Run exact three pilot samples and compare to golden results")
    reg.add_argument("--config",required=True)
    reg.add_argument("--manifest",required=True)
    reg.add_argument("--expected",default=str(HERE/"benchmarks.json"))
    reg.add_argument("--run-id",required=True)
    reg.add_argument("--mode",choices=["fresh","reuse"],required=True)
    reg.add_argument("--keep-going",action="store_true")
    reg.add_argument("--recheck-ancestry",action="store_true")
    a=parser.parse_args()
    try:
        if a.command=="run":cmd_run(a)
        else:cmd_regression(a)
    except Exception as e:
        print("ERROR:",e,file=sys.stderr)
        sys.exit(1)

if __name__=="__main__":cli()
