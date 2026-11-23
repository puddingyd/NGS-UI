#!/usr/bin/env python3
"""Tertiary post-processing for one or many samples, step by step.

Every annotation step runs once for the whole batch: each resource (GeneBe
SQLite, SpliceAI/VEP, in-house AF sites VCF, ClinVar / LitVar2 SQLite, GIAB
BEDs, MANE summary, GPN score tables) is opened and queried once for the
union of all samples' variants, and the results are merged back into each
sample's disposable working TSV.  All lookups are per variant, so each
sample's output is identical to running the samples one at a time; only the
repeated per-sample database work disappears.

Order (unchanged from the per-sample chain):
    1. clinvar-latest   weekly ClinVar comparison (UI only)
    2. genebe           GeneBe ACMG second opinion (local DB → cache → API)
    3. spliceai         SpliceAI via VEP (Research-only)
    4. giab-strata      GIAB difficult-region labels
    5. inhouse-af       in-house allele frequency
    6. mane-refseq      Ensembl → MANE RefSeq IDs
    7. litvar2          local LitVar2 literature lookup
    8. gpn-prefetch     GPN-MSA / GPN-Star scores for all review keys
    then per sample, in parallel: snv-overlay → review-tsv → gene-index → acmg-engine
    (acmg-engine is best-effort: a failure is logged and does not fail the batch)

Single-sample use (``run_stopgaps.sh`` is a thin wrapper around this)::

    scripts/postprocess_batch.py --work-tsv W --raw-tsv R --post-dir D \\
        --sample SID --seq-type WES [--research-only] [...]

Batch use (the tertiary worker)::

    scripts/postprocess_batch.py --manifest batch.json [--research-only] [...]

with ``batch.json`` = ``{"samples": [{"sample", "work_tsv", "raw_tsv",
"post_dir", "seq_type"}, ...]}``.

Markers on stdout drive the worker's progress/log summary:
    [post-processing-step] <step> start|done elapsed=Ns       (batch-wide)
    [sample-step] <step> start|done sample=<SID> elapsed=Ns   (per sample)
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import os
import subprocess
import sys
import tempfile
import threading
import time
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timedelta, timezone
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent
sys.path.insert(0, str(REPO_ROOT / "backend"))

from app.services.parallel_tsv import fork_map  # noqa: E402

TAIPEI = timezone(timedelta(hours=8))
_PRINT_LOCK = threading.Lock()


def emit(line: str = "") -> None:
    with _PRINT_LOCK:
        print(line, flush=True)


def _stamp() -> str:
    return datetime.now(TAIPEI).strftime("%Y-%m-%dT%H:%M:%S+08:00")


class Step:
    """Context manager printing start/done markers with elapsed seconds."""

    def __init__(self, name: str, marker: str = "post-processing-step", sample: str = ""):
        self.name, self.marker, self.sample = name, marker, sample

    def __enter__(self):
        self.started = time.monotonic()
        suffix = f" sample={self.sample}" if self.sample else ""
        emit(f"[{_stamp()}] [{self.marker}] {self.name} start{suffix}")
        return self

    def __exit__(self, exc_type, *_rest):
        if exc_type is None:
            suffix = f" sample={self.sample}" if self.sample else ""
            elapsed = int(round(time.monotonic() - self.started))
            emit(f"[{_stamp()}] [{self.marker}] {self.name} done{suffix} elapsed={elapsed}s")
        return False


def _load_script(name: str):
    spec = importlib.util.spec_from_file_location(name, SCRIPT_DIR / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    assert spec and spec.loader
    spec.loader.exec_module(module)
    return module


def parse_args(argv=None) -> argparse.Namespace:
    home = Path(os.environ.get("NGS_UI_HOME") or Path.home() / "NGS_UI")
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--manifest", type=Path, help="batch JSON (see module docstring)")
    ap.add_argument("--tsv", "--work-tsv", dest="work_tsv", type=Path)
    ap.add_argument("--raw-tsv", type=Path)
    ap.add_argument("--post-dir", type=Path)
    ap.add_argument("--sample", default="")
    ap.add_argument("--seq-type", default=os.environ.get("SEQ_TYPE", "WES"))
    ap.add_argument("--research-only", action="store_true")
    from app.services.parallel_tsv import default_workers
    ap.add_argument("--workers", type=int, default=default_workers(),
                    help="parallel per-sample TSV passes and review/overlay/index builds "
                         "(env NGS_UI_POSTPROC_WORKERS; default min(8, CPUs))")
    ap.add_argument("--raw-log", type=Path, help="append verbose per-sample tool output here")
    ap.add_argument("--candidate-bed", type=Path,
                    default=Path(os.environ.get("NGS_UI_CDS_CANDIDATE_BED") or home / "biotools/cds_combined.bed"))
    ap.add_argument("--skip-candidate-bed", action="store_true")
    ap.add_argument("--clinvar-latest-db", type=Path,
                    default=Path(os.environ.get("NGS_UI_CLINVAR_LATEST_DB")
                                 or home / "biotools/clinvar_latest/clinvar_latest.sqlite"))
    ap.add_argument("--skip-clinvar-latest", action="store_true")
    ap.add_argument("--genebe-db", type=Path,
                    default=Path(os.environ.get("NGS_UI_GENEBE_DB") or home / "biotools/genebe/genebe_hg38.tsv.gz"))
    ap.add_argument("--skip-genebe", action="store_true")
    # Keep the established production SpliceAI paths (not env-driven).
    ap.add_argument("--spliceai-snv", type=Path, default=home / "biotools/spliceai/spliceai_scores.raw.snv.hg38.vcf.gz")
    ap.add_argument("--spliceai-indel", type=Path, default=home / "biotools/spliceai/spliceai_scores.raw.indel.hg38.vcf.gz")
    ap.add_argument("--skip-spliceai", "--skip-extra-vep", dest="skip_spliceai", action="store_true")
    ap.add_argument("--giab-strat-dir", type=Path, default=None)
    ap.add_argument("--skip-giab", action="store_true")
    ap.add_argument("--inhouse-af-db", type=Path,
                    default=Path(os.environ.get("NGS_UI_INHOUSE_AF_DB") or home / "biotools/inhouse_af/inhouse_af.hg38.vcf.gz"))
    ap.add_argument("--skip-inhouse-af", action="store_true")
    ap.add_argument("--mane-summary", type=Path,
                    default=Path(os.environ.get("NGS_UI_MANE_SUMMARY") or home / "biotools/MANE.GRCh38.v1.5.summary.txt.gz"))
    ap.add_argument("--skip-mane-refseq", action="store_true")
    ap.add_argument("--litvar2-db", type=Path,
                    default=Path(os.environ.get("NGS_UI_LITVAR2_DB") or home / "biotools/litvar2/litvar2.sqlite"))
    ap.add_argument("--skip-litvar2", action="store_true")
    ap.add_argument("--skip-acmg-engine", action="store_true",
                    help="do not run the ACMG rule engine (plan §5.2)")
    ap.add_argument("--gpn-msa-db", type=Path,
                    default=Path(os.environ.get("NGS_UI_GPN_MSA_DB") or home / "biotools/gpn_msa/scores.tsv.bgz"))
    ap.add_argument("--gpn-star-dir", type=Path,
                    default=Path(os.environ.get("NGS_UI_GPN_STAR_DIR") or home / "biotools/gpn_star"))
    args = ap.parse_args(argv)
    if args.giab_strat_dir is None and os.environ.get("NGS_UI_GIAB_STRAT_DIR"):
        args.giab_strat_dir = Path(os.environ["NGS_UI_GIAB_STRAT_DIR"])
    return args


def load_samples(args: argparse.Namespace) -> list[dict]:
    if args.manifest:
        raw = json.loads(args.manifest.read_text(encoding="utf-8"))
        entries = raw.get("samples", raw) if isinstance(raw, dict) else raw
    else:
        if not args.work_tsv:
            raise SystemExit("ERROR: --tsv/--work-tsv or --manifest required")
        entries = [{
            "sample": args.sample, "work_tsv": str(args.work_tsv),
            "raw_tsv": str(args.raw_tsv or args.work_tsv),
            "post_dir": str(args.post_dir or args.work_tsv.parent), "seq_type": args.seq_type,
        }]
    samples = []
    for entry in entries:
        work = Path(entry["work_tsv"])
        raw = Path(entry.get("raw_tsv") or work)
        post = Path(entry.get("post_dir") or work.parent)
        for label, path in (("--tsv", work), ("--raw-tsv", raw)):
            if not path.is_file():
                raise SystemExit(f"ERROR: {label} not found: {path}")
        seq = str(entry.get("seq_type") or "WES").upper()
        if seq not in ("WES", "WGS"):
            raise SystemExit(f"ERROR: --seq-type must be WES or WGS (got: {seq})")
        sid = str(entry.get("sample") or "")
        if not sid:
            # Working TSVs live below tertiary_output/<SID>/08_postprocessing/;
            # legacy invocations point directly below tertiary_output/<SID>/.
            sid = (work.parent.parent.name if work.parent.name == "08_postprocessing" else work.parent.name)
        post.mkdir(parents=True, exist_ok=True)
        samples.append({"sample": sid, "work": work, "raw": raw, "post": post, "seq": seq})
    if len({s["sample"] for s in samples}) != len(samples):
        raise SystemExit("ERROR: duplicate sample in post-processing batch")
    return samples


def _tool_manifest(samples: list[dict], directory: Path) -> Path:
    path = directory / "tool-batch.json"
    path.write_text(json.dumps([
        {"tsv": str(s["work"]), "test_type": s["seq"], "sample": s["sample"]} for s in samples
    ], ensure_ascii=False), encoding="utf-8")
    return path


def _run_tool(cmd: list[str]) -> None:
    """Stream a batch tool's output; any failure aborts post-processing."""
    proc = subprocess.run(cmd, check=False)
    if proc.returncode != 0:
        raise SystemExit(proc.returncode or 1)


def run_annotation_steps(args: argparse.Namespace, samples: list[dict], workdir: Path) -> None:
    py = sys.executable
    tool_manifest = _tool_manifest(samples, workdir)

    # 1. Weekly ClinVar comparison: one shared lookup cache for every sample.
    if not args.skip_clinvar_latest:
        if args.clinvar_latest_db.is_file():
            clinvar = _load_script("annotate_clinvar_latest")
            with Step("clinvar-latest"):
                lookup = clinvar.CachedLookup(args.clinvar_latest_db)

                def clinvar_keys(s):
                    import csv
                    keys, ids = set(), set()
                    with s["work"].open("r", encoding="utf-8", newline="") as handle:
                        for row in csv.DictReader(handle, delimiter="\t"):
                            keys.add(clinvar._lookup_key(row))
                            ids.add(clinvar._variation_id(row))
                    return keys, ids

                collected = fork_map(clinvar_keys, samples, args.workers)
                # One lookup for every sample's keys/IDs; workers then merge
                # from the (inherited) prefilled cache without touching the DB.
                lookup(args.clinvar_latest_db,
                       set().union(*(k for k, _i in collected)), set().union(*(i for _k, i in collected)))

                def clinvar_merge(s):
                    try:
                        return clinvar.annotate_tsv(
                            s["work"], args.clinvar_latest_db,
                            s["post"] / f"{s['sample']}.clinvar_comparison.json",
                            baseline_release="2026-07-20", lookup=lookup,
                        )["stats"]
                    except (OSError, RuntimeError, ValueError) as exc:
                        return {"error": str(exc)}

                for s, stats in zip(samples, fork_map(clinvar_merge, samples, args.workers)):
                    if "error" in stats:
                        print(f"ERROR: {stats['error']}", file=sys.stderr, flush=True)
                        raise SystemExit(2)
                    emit(f"[clinvar-latest] {s['sample']}: {json.dumps(stats, ensure_ascii=False)}")
        else:
            emit(f"  - weekly ClinVar comparison skipped (DB not found: {args.clinvar_latest_db})")

    # 2. GeneBe: one local-DB lookup / API-cache pass / API run for the union.
    if not args.skip_genebe:
        with Step("genebe"):
            if not args.genebe_db.is_file():
                print(f"ERROR: GeneBe DB not found: {args.genebe_db}", file=sys.stderr, flush=True)
                print("       set NGS_UI_GENEBE_DB / --genebe-db, or pass --skip-genebe", file=sys.stderr, flush=True)
                raise SystemExit(2)
            _run_tool([py, str(SCRIPT_DIR / "annotate_acmg_genebe.py"),
                       "--batch-manifest", str(tool_manifest), "--genebe-db", str(args.genebe_db)])

    # 3. SpliceAI: one VEP run on the union of review sites.
    if not args.skip_spliceai:
        cmd = [py, str(SCRIPT_DIR / "annotate_spliceai.py"), "--batch-manifest", str(tool_manifest),
               "--spliceai-snv", str(args.spliceai_snv), "--spliceai-indel", str(args.spliceai_indel)]
        if not args.skip_candidate_bed and args.candidate_bed.is_file():
            cmd += ["--candidate-bed", str(args.candidate_bed)]
        with Step("spliceai"):
            _run_tool(cmd)
    else:
        emit("  - SpliceAI skipped (--skip-spliceai / Research-only unchecked)")

    # 3b. GIAB: load the stratification BEDs once.
    if not args.skip_giab:
        giab = _load_script("annotate_giab_strata")
        with Step("giab-strata"):
            strat_dir = args.giab_strat_dir or Path(giab.DEFAULT_STRAT_DIR)
            beds = giab.load_strata(strat_dir, None)
            if beds:
                fork_map(lambda s: giab.annotate_with(s["work"], beds, giab.DEFAULT_COLUMN), samples, args.workers)

    # 3b2. In-house AF: one sites-VCF scan per memory-bounded chunk.
    if not args.skip_inhouse_af:
        if args.inhouse_af_db.is_file():
            inhouse = _load_script("annotate_inhouse_af")
            import csv
            csv.field_size_limit(inhouse._CSV_MAX)
            with Step("inhouse-af"):
                inhouse.annotate_many([(s["work"], s["sample"]) for s in samples], str(args.inhouse_af_db),
                                      workers=args.workers)
        else:
            emit(f"  - in-house AF skipped (DB not found: {args.inhouse_af_db})")

    # 3c. MANE RefSeq: read the summary once.
    if not args.skip_mane_refseq:
        mane_mod = _load_script("annotate_mane_refseq")
        with Step("mane-refseq"):
            if not args.mane_summary.is_file():
                emit(f"[mane-refseq] MANE summary not found, skipped: {args.mane_summary}")
            else:
                mapping = mane_mod._load_mane(args.mane_summary)
                counts = fork_map(lambda s: mane_mod.annotate_tsv(s["work"], args.mane_summary, mapping),
                                  samples, args.workers)
                for s, (total, matched) in zip(samples, counts):
                    emit(f"[mane-refseq] {s['sample']}: rows={total} matched={matched} summary={args.mane_summary}")

    # 4. LitVar2: one read-only connection + memo for every sample.
    if not args.skip_litvar2:
        litvar = _load_script("annotate_litvar2")
        with Step("litvar2"):
            if not args.litvar2_db.is_file():
                emit(f"[litvar2] local database not found, skipped: {args.litvar2_db}")
            else:
                shared = None
                try:
                    shared = litvar.SharedLookup(args.litvar2_db)
                except Exception as exc:  # noqa: BLE001 — reviewer context only
                    print(f"[litvar2] WARNING: local annotation failed, skipped: {exc}", file=sys.stderr, flush=True)
                if shared is not None:
                    try:
                        def groups_of(s):
                            try:
                                return litvar._candidate_groups(s["work"], test_type=s["seq"])[0]
                            except Exception:  # noqa: BLE001 — reported by the merge pass
                                return {}

                        # One memoised lookup per distinct identifier group,
                        # then workers merge from the inherited memo.
                        for groups in fork_map(groups_of, samples, args.workers):
                            for identifiers in groups.values():
                                shared.lookup(identifiers)
                        shared.close()

                        def litvar_merge(s):
                            try:
                                return litvar.annotate_tsv(
                                    s["work"], args.litvar2_db, test_type=s["seq"],
                                    marker_path=s["post"] / f"{s['sample']}.litvar2_annotation.json",
                                    shared=shared,
                                )
                            except Exception as exc:  # noqa: BLE001 — never a release gate
                                return {"error": str(exc)}

                        for s, stats in zip(samples, fork_map(litvar_merge, samples, args.workers)):
                            if "error" in stats:
                                print(f"[litvar2] WARNING: {s['sample']}: local annotation failed, skipped: "
                                      f"{stats['error']}", file=sys.stderr, flush=True)
                                continue
                            emit(f"[litvar2] {s['sample']}: "
                                 + " ".join(f"{k}={v}" for k, v in stats.items()) + f" db={args.litvar2_db}")
                    finally:
                        shared.close()


def prefetch_gpn(args: argparse.Namespace, samples: list[dict], workdir: Path) -> Path:
    from app.services import gpn_prefetch

    path = workdir / "gpn_prefetch.sqlite"
    with Step("gpn-prefetch"):
        gpn_prefetch.build_for_tsvs(
            path, [(s["work"], s["seq"]) for s in samples],
            msa_db=args.gpn_msa_db, star_dir=args.gpn_star_dir,
            research_only=args.research_only, log=emit, workers=args.workers,
        )
    return path


def _capture(cmd: list[str], env: dict) -> tuple[int, str]:
    proc = subprocess.run(cmd, check=False, env=env, stdout=subprocess.PIPE,
                          stderr=subprocess.STDOUT, text=True, errors="replace")
    return proc.returncode, proc.stdout


def sample_tail(args: argparse.Namespace, s: dict, env: dict, raw_log_lock: threading.Lock) -> None:
    """snv-overlay → review-tsv → gene-index → acmg-engine for one sample."""
    sid, post = s["sample"], s["post"]
    overlay = post / f"{sid}.snv_annotations.sqlite"
    review_python = os.environ.get("NGS_UI_PYTHON_BIN") or sys.executable
    gpn_args = ["--gpn-msa-db", str(args.gpn_msa_db), "--gpn-star-dir", str(args.gpn_star_dir)]
    if args.research_only:
        gpn_args.append("--research-only")
    steps = [
        ("snv-overlay", [sys.executable, str(SCRIPT_DIR / "build_snv_annotation_overlay.py"),
                         "--raw", str(s["raw"]), "--annotated", str(s["work"]), "--out", str(overlay)]),
        ("review-tsv", [review_python, str(SCRIPT_DIR / "build_snv_review_tsv.py"),
                        "--tsv", str(s["raw"]), "--output-dir", str(post),
                        "--output-path", str(post / f"{sid}.snv_indel.review.tsv"),
                        "--manifest-path", str(post / f"{sid}.snv_indel.review.tsv.source.json"),
                        "--overlay", str(overlay), "--test-type", s["seq"], "--sample", sid, *gpn_args]),
        ("gene-index", [sys.executable, str(SCRIPT_DIR / "build_snv_gene_index.py"),
                        "--tsv", str(s["raw"]), "--out", str(post / f"{sid}.snv_gene_index.sqlite")]),
    ]
    if not args.skip_acmg_engine:
        steps.append(("acmg-engine", [sys.executable, str(SCRIPT_DIR / "acmg_engine" / "build_engine_results.py"),
                                      "--review", str(post / f"{sid}.snv_indel.review.tsv"),
                                      "--out", str(post / f"{sid}.snv_acmg_engine.sqlite"),
                                      "--sample", sid, "--test-type", s["seq"]]))
    for name, cmd in steps:
        if name == "acmg-engine":
            # New rule engine: never block publication of the case.
            with Step(name, marker="sample-step", sample=sid):
                code, output = _capture(cmd, env)
                if args.raw_log is not None:
                    with raw_log_lock, args.raw_log.open("a", encoding="utf-8") as handle:
                        handle.write(f"\n[{_stamp()}] {sid} {name}\n$ {' '.join(cmd)}\n{output}")
                for line in output.splitlines():
                    if line.startswith("[acmg-engine]"):
                        emit(line)
                if code != 0:
                    emit(f"[acmg-engine] WARNING sample={sid} failed (exit {code}); "
                         "results can be rebuilt with scripts/acmg_engine/build_engine_results.py")
                    emit(output.rstrip("\n")[-2000:])
            continue
        with Step(name, marker="sample-step", sample=sid):
            code, output = _capture(cmd, env)
            if args.raw_log is not None:
                with raw_log_lock, args.raw_log.open("a", encoding="utf-8") as handle:
                    handle.write(f"\n[{_stamp()}] {sid} {name}\n$ {' '.join(cmd)}\n{output}")
            elif name == "review-tsv" or code != 0:
                emit(output.rstrip("\n"))
            if code != 0:
                if args.raw_log is not None:
                    emit(output.rstrip("\n")[-4000:])
                raise RuntimeError(f"{name} failed for {sid} (exit {code})")
            if name == "review-tsv" and args.raw_log is not None:
                # Keep the compact predictor status visible in the job log.
                for line in output.splitlines():
                    if line.startswith(("[gpn-msa]", "[gpn-star]", "[gpn-star-v]")):
                        emit(f"{line.split(' manifest=')[0]} sample={sid}")


def run_tails(args: argparse.Namespace, samples: list[dict], prefetch: Path | None) -> None:
    env = dict(os.environ)
    if prefetch is not None and prefetch.is_file():
        env["NGS_UI_GPN_PREFETCH_DB"] = str(prefetch)
    raw_log_lock = threading.Lock()
    errors: list[BaseException] = []
    with ThreadPoolExecutor(max_workers=max(1, min(args.workers, len(samples)))) as pool:
        futures = [pool.submit(sample_tail, args, s, env, raw_log_lock) for s in samples]
        for future in futures:
            try:
                future.result()
            except BaseException as exc:  # noqa: BLE001 — report every sample, then fail
                errors.append(exc)
    if errors:
        for exc in errors:
            print(f"ERROR: {exc}", file=sys.stderr, flush=True)
        raise SystemExit(1)


def main(argv=None) -> int:
    args = parse_args(argv)
    samples = load_samples(args)
    # Every candidate-scoped annotator imports the same review predicate and
    # BED path; export the resolved value so GeneBe API, SpliceAI, LitVar2,
    # GPN prefetch and the review builder cannot use different regions.
    if args.skip_candidate_bed:
        os.environ["NGS_UI_CDS_CANDIDATE_BED"] = str(samples[0]["post"] / ".candidate-bed-disabled")
        bed_note = "skipped (--skip-candidate-bed)"
    else:
        os.environ["NGS_UI_CDS_CANDIDATE_BED"] = str(args.candidate_bed)
        bed_note = (str(args.candidate_bed) if args.candidate_bed.is_file()
                    else f"not found at {args.candidate_bed} (final-review rule uses AF-only candidates)")
    emit("================================================================")
    emit(f"  post-processing : {len(samples)} sample(s): {', '.join(s['sample'] for s in samples)}")
    emit(f"  candidate BED: {bed_note}")
    emit(f"  GPN-MSA DB: {args.gpn_msa_db} (best-effort)")
    emit(f"  GPN-Star DB: {args.gpn_star_dir} (research-only={int(args.research_only)}, best-effort)")
    emit("================================================================")
    with tempfile.TemporaryDirectory(prefix=".postprocess-", dir=str(samples[0]["post"])) as tmp:
        workdir = Path(tmp)
        run_annotation_steps(args, samples, workdir)
        prefetch = prefetch_gpn(args, samples, workdir)
        run_tails(args, samples, prefetch)
    emit("================================================================")
    for s in samples:
        emit(f"  done. {s['sample']}: raw TSV {s['raw']}; sparse overlay {s['post'] / (s['sample'] + '.snv_annotations.sqlite')}")
    emit("================================================================")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
