"""Version-pinned PharmCAT revisions. Only complete, verified reports are published.

The original report contains native GeneReports. Replay those losslessly; replace
only reviewed genes using the phenotyper, then run the reporter on the whole set.
No matcher output, phase, phenotype or activity score is inferred by this service.
"""
from __future__ import annotations

import copy
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import time
import uuid
from collections import Counter
from contextlib import contextmanager

from .. import config
from . import sample_layout

VERSION = "3.2.0"
DATA_VERSION = "2026-02-09-10-28"
JAR_SHA256 = "a1ceda4c4f8e2c53ed17ff4d104df82360c04d4b15226d09cb7f854bc585a0e7"
# Hospital operator reported `-version` and sha256sum on 2026-11-06.
# Registration identifies a build; per-case replay still proves compatibility.
HOSPITAL_JAR_SHA256 = "ff3e16f4c2819a5bf315f2c2852e420836944e82548f364758fab39cfe711835"
APPROVED_JARS = {JAR_SHA256: "official-release", HOSPITAL_JAR_SHA256: "hospital-3.2.0"}
SCHEMA = 1
ACTIVE = {"queued", "running", "cancelling"}


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, ensure_ascii=False).encode()).hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def atomic(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(mode="w", dir=path.parent, prefix=".pgx-", delete=False) as f:
        tmp = Path(f.name)
        try:
            json.dump(value, f, ensure_ascii=False, indent=2)
            f.flush()
            os.fsync(f.fileno())
            tmp.replace(path)
        finally:
            tmp.unlink(missing_ok=True)


def root():
    return config.DATA_ROOT / "jobs" / "pgx_reannotation"


def job_dir(jid):
    if not re.fullmatch(r"[a-f0-9]{32}", str(jid)):
        raise ValueError("無效的 PGx 工作")
    return root() / jid


@contextmanager
def review_lock(sid):
    from .pgx_concordance import REVIEW_NAME
    path = sample_layout.state_file(sid, REVIEW_NAME, for_write=True)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.with_suffix(".lock").open("a+") as handle:
        fcntl.flock(handle, fcntl.LOCK_EX)
        yield path


def key(choices):
    return digest({g: {k: c.get(k) for k in ("source", "genotype", "fingerprint")} for g, c in choices.items()})


def request_key(pgx, choices, engine_sha256=JAR_SHA256):
    return digest([SCHEMA, engine_sha256, pgx.get("source_identity"), key(choices),
                   {g: r.get("fingerprint") for g, r in pgx["concordance"]["genes"].items()}])


def _alive(job):
    directory = job_dir(job["job_id"])
    # The lease belongs to the worker, never a PID that might have been reused.
    with (directory / "lease").open("a+") as handle:
        try:
            fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            return True
    return job.get("status") == "queued" and time.time() - job.get("created", 0) < 30


def job_status(jid):
    try:
        job = read(job_dir(jid) / "state.json")
        if not isinstance(job, dict):
            raise ValueError("invalid PGx job")
        if job.get("status") in ACTIVE and not _alive(job):
            job = {**job, "status": "failed", "error": "更新中斷，請重試"}
        return job
    except (OSError, ValueError, TypeError):
        return {"status": "failed", "error": "更新紀錄無法讀取，請重試"}


def active_ids():
    ids = set()
    if root().exists():
        for path in root().glob("*/state.json"):
            try:
                job = job_status(path.parent.name)
                if job.get("status") in ACTIVE:
                    ids.update([job["sample_id"], job["source_sample_id"]])
            except (KeyError, OSError, ValueError):
                continue
    return ids


def revision_dir(sid, jid):
    job_dir(jid)  # validate before constructing a reviewer-owned path
    return sample_layout.state_dir(sid) / f"{sid}.pgx_annotations" / jid


def _validated_revision(sid, pointer, pgx):
    directory = revision_dir(sid, pointer["job_id"])
    manifest = read(directory / "manifest.json")
    report_path = directory / "report.json"
    if (manifest.get("schema") != SCHEMA or manifest.get("engine_sha256") not in APPROVED_JARS
            or manifest.get("baseline_replay_verified") is not True
            or manifest.get("sample_id") != sid
            or manifest.get("request_key") != pointer.get("request_key")
            or manifest.get("request_key") != request_key(pgx, manifest["choices"], manifest["engine_sha256"])
            or hashlib.sha256(report_path.read_bytes()).hexdigest() != manifest.get("report_sha256")):
        raise ValueError("PGx 更新結果與目前來源不符")
    return manifest, report_path


def project(sid, pgx, review):
    """Called after baseline concordance: never fingerprint a derived report."""
    from ..adapters.pgx_tsv import load_pgx, _is_actionable
    choices = review.get("choices") or {}
    pointer = review.get("published")
    status = {"status": "original", "busy": False, "error": ""}
    published_choices = {}
    if pointer:
        try:
            manifest, report_path = _validated_revision(sid, pointer, pgx)
            raw = read(report_path)
            clinical = load_pgx(sample_layout.pgx_tsv(sid), report_path)
            # Remove every old TSV advice row, including shared gene guidance.
            for gene, payload in clinical["genes"].items():
                details = payload.get("details") or {}
                source = ((raw.get("genes", {}).get(gene) or {}).get("sourceDiplotypes") or [{}])[0]
                if gene in manifest["choices"]:
                    payload.update(diplotype=details.get("label", ""),
                                   phenotype="; ".join(details.get("phenotypes") or []),
                                   activity_score=source.get("activityScore") or "")
                annotations = [a for a in clinical["guideline_annotations"] if gene in a.get("genes", [])]
                payload["drugs"] = [{"drug": a["drug"], "recommendations": [{"source": a["source"],
                    "recommendation": a["recommendation"], "evidence": a["classification"]}], "best_rank": 0}
                    for a in annotations]
                payload["actionable"] = _is_actionable(gene, "; ".join(details.get("phenotypes") or []), [])
            for field in ("genes", "guideline_annotations", "messages", "timestamp", "gene_order",
                          "additional_genes", "summary"):
                pgx[field] = clinical[field]
            pgx["actionable"] = [g for g, p in pgx["genes"].items() if p.get("actionable")]
            pgx["routine"] = [g for g, p in pgx["genes"].items() if not p.get("actionable")]
            pgx["summary"].update(actionable_genes=len(pgx["actionable"]),
                recommendations=len(pgx["guideline_annotations"]),
                strong_recommendations=sum(a.get("classification") == "Strong" for a in pgx["guideline_annotations"]))
            published_choices = manifest["choices"]
            status["status"] = "complete"
        except (OSError, ValueError, KeyError, TypeError):
            pgx["review_blocked"] = True
            status.update(status="failed", error="PGx 更新結果需重新確認")
    for gene, row in pgx["concordance"]["genes"].items():
        requested = choices.get(gene) or {}
        row.update(requested_source="dragen" if requested else "platform", requested_genotype=requested.get("genotype", ""))
        active = published_choices.get(gene)
        if active:
            if requested and key({gene: requested}) == key({gene: active}):
                active = requested
            row.update(selected_source="dragen", selected_genotype=active["genotype"],
                       reviewer=active.get("reviewer_username", ""), reviewed_at=active.get("updated_at", ""))
    job = job_status(review["annotation_job"]) if review.get("annotation_job") else {}
    published_key = pointer.get("request_key") if isinstance(pointer, dict) else None
    awaiting_job = bool(review.get("annotation_job")) and job.get("request_key") != published_key
    matches = key(choices) == key(published_choices) and not awaiting_job
    if not matches or status["status"] == "failed":
        pgx["review_blocked"] = True
        status.update(status=job.get("status", "needs_update"), error=job.get("error") or status["error"])
        if status["status"] == "complete":
            status["status"] = "needs_update"
        status["busy"] = status["status"] in ACTIVE
    elif not choices and not pointer:
        status["status"] = "original"
    pgx["reannotation"] = status
    pgx["review_fingerprint"] = digest([pgx.get("review_fingerprint"), pointer, status["status"], pgx["review_blocked"]])


def _launch(jid):
    directory = job_dir(jid)
    env = os.environ.copy()
    env["PYTHONPATH"] = str(config.REPO_ROOT / "backend") + os.pathsep + env.get("PYTHONPATH", "")
    with (directory / "worker.log").open("ab") as log:
        subprocess.Popen([sys.executable, "-m", "app.workers.pgx_reannotate", jid],
                         cwd=config.REPO_ROOT, env=env, stdin=subprocess.DEVNULL,
                         stdout=log, stderr=subprocess.STDOUT, start_new_session=True)


def schedule(sid, pgx, review):
    """Under sample/review locks; persist request before spawning its worker."""
    from .pgx_concordance import REVIEW_NAME
    path = sample_layout.state_file(sid, REVIEW_NAME, for_write=True)
    choices = review.get("choices") or {}
    if not choices:
        review.pop("published", None)
        review.pop("annotation_job", None)
        atomic(path, review)
        return
    engine_error = ""
    try:
        engine_sha256 = runtime().get("sha256", JAR_SHA256)
    except (OSError, ValueError, subprocess.SubprocessError) as exc:
        engine_sha256 = None
        engine_error = str(exc) if isinstance(exc, ValueError) else "PharmCAT 執行環境無法啟動，請管理者檢查"
    wanted = request_key(pgx, choices, engine_sha256)
    # Reuse only immutable revisions that still validate against current inputs.
    for cached in reversed(review.get("revisions") or []):
        if not engine_error and cached.get("request_key") == wanted:
            try:
                _validated_revision(sid, cached, pgx)
                review.update(published=cached)
                review.pop("annotation_job", None)
                review["audit"].append({"action": "annotation_reused", **cached, "time": time.time()})
                atomic(path, review)
                return
            except (OSError, ValueError, KeyError, TypeError):
                continue
    jid = uuid.uuid4().hex
    directory = job_dir(jid)
    directory.mkdir(parents=True)
    job = {"job_id": jid, "sample_id": sid, "source_sample_id": sample_layout.source_sample_id(sid),
           "status": "failed" if engine_error else "queued", "error": engine_error,
           "created": time.time(), "request_key": wanted, "engine_sha256": engine_sha256,
           "choices": copy.deepcopy(choices)}
    atomic(directory / "state.json", job)
    review["annotation_job"] = jid
    atomic(path, review)
    if engine_error:
        return
    try:
        _launch(jid)
    except OSError:
        atomic(directory / "state.json", {**job, "status": "failed", "error": "無法啟動更新，請重試"})


def action(sid, operation, user):
    """Retry or restore the original complete report; never accept executable input."""
    from . import pgx_concordance as pc
    from .sample_loader import load_sample_pgx
    from ..workers.dragen_run import _acquire_sample_locks, _release_sample_locks
    if operation not in {"retry", "restore"} or not re.fullmatch(r"[A-Za-z0-9_][A-Za-z0-9_.-]*", sid):
        raise ValueError("無效的 PGx 操作")
    source_id = sample_layout.source_sample_id(sid)
    if not re.fullmatch(r"[A-Za-z0-9_][A-Za-z0-9_.-]*", source_id):
        raise ValueError("無效的樣本 ID")
    with review_lock(sid):
        review = pc._review(sid)
        jid = review.get("annotation_job")
        if jid and job_status(jid).get("status") in ACTIVE:
            if operation != "restore":
                raise ValueError("PGx 更新中，請稍後")
            atomic(job_dir(jid) / "restore.json", {"username": user.get("username", ""), "time": time.time()})
            return load_sample_pgx(sid)
    handles = _acquire_sample_locks([sid, source_id])
    try:
        with review_lock(sid) as path:
            review = pc._review(sid)
            pgx = load_sample_pgx(sid)["pgx"]
            if pgx.get("reannotation", {}).get("busy"):
                raise ValueError("PGx 更新中，請稍後再試")
            if not sample_layout.state_file(sid, "sample_metadata.json").is_file():
                raise ValueError("個案尚未登錄")
            if operation == "restore":
                review["choices"] = {}
            elif not review.get("choices") or any(pgx["concordance"]["genes"].get(g, {}).get("stale", True) for g in review["choices"]):
                raise ValueError("比對來源已變更，請重新確認")
            review["audit"].append({"action": operation, "reviewer_username": user.get("username", ""), "time": time.time()})
            schedule(sid, pgx, review)
        return load_sample_pgx(sid)
    finally:
        _release_sample_locks(handles)


def runtime():
    """Only registered builds; identify the actual JAR, including inside SIF."""
    jar = Path(os.environ.get("NGS_UI_PGX_PHARMCAT_JAR", config.BIOTOOLS_DIR / "pharmcat" / f"pharmcat-{VERSION}-all.jar"))
    image = os.environ.get("NGS_UI_PGX_PHARMCAT_IMAGE", "")
    if not jar.is_file() and not image:
        for base in ("/home/datalake_Intermediate/pipeline/nextflow_containers", "/home/pipeline/nextflow_containers"):
            candidate = Path(base) / f"pharmcat_{VERSION}.sif"
            if candidate.is_file():
                image = str(candidate)
                break
    if jar.is_file():
        sha256 = hashlib.sha256(jar.read_bytes()).hexdigest()
        engine = {"jar": str(jar.resolve()), "image": image, "sha256": sha256, "jar_on_host": True}
    elif image:
        jar_inside = "/pharmcat/pharmcat.jar"
        result = subprocess.run(["apptainer", "exec", "--cleanenv", image, "sha256sum", jar_inside],
                                capture_output=True, text=True, timeout=30, check=True, env=execution_env())
        words = result.stdout.split()
        sha256 = words[0] if words else ""
        engine = {"jar": jar_inside, "image": image, "sha256": sha256, "jar_on_host": False}
    else:
        raise ValueError("尚未準備固定版本 PharmCAT，請管理者完成設定後重試")
    if sha256 not in APPROVED_JARS:
        raise ValueError(f"PharmCAT {'容器' if image else ''}程式校驗不符：未登錄的 JAR SHA256 {sha256 or '（無輸出）'}")
    with tempfile.TemporaryDirectory(prefix="pgx-version-") as name:
        directory = Path(name)
        (directory / "tmp").mkdir()
        result = subprocess.run(command(engine, directory, ["-version"]), cwd=directory,
                                capture_output=True, text=True, timeout=30, check=True, env=execution_env(directory))
    if f"PharmCAT {VERSION}" not in [line.strip() for line in (result.stdout + "\n" + result.stderr).splitlines()]:
        raise ValueError("PharmCAT 版本不符，必須為 3.2.0")
    return engine


def execution_env(directory=None):
    env = os.environ.copy()
    for name in list(env):
        if name.startswith(("APPTAINER_", "SINGULARITY_")) or name in {"JAVA_TOOL_OPTIONS", "_JAVA_OPTIONS", "JDK_JAVA_OPTIONS"}:
            env.pop(name, None)
    if directory is not None:
        env["TMPDIR"] = str(Path(directory) / "tmp")
    return env


def command(engine, directory, arguments):
    args = [config.JAVA_BIN if not engine["image"] else "java", "-Xmx1g",
            f"-Djava.io.tmpdir={directory / 'tmp'}", "-jar", engine["jar"], *arguments]
    if engine["image"]:
        binds = [str(directory)]
        if engine.get("jar_on_host"):
            binds.append(engine["jar"] + ":" + engine["jar"] + ":ro")
        args = ["apptainer", "exec", "--cleanenv", "--bind", ",".join(binds), engine["image"], *args]
    return args


def check_original(raw):
    if raw.get("pharmcatVersion") != VERSION or raw.get("dataVersion") != DATA_VERSION or not raw.get("genes"):
        raise ValueError("原 PharmCAT 版本不相容，暫不更新")
    return {"matcherMetadata": raw.get("matcherMetadata"), "geneReports": copy.deepcopy(raw["genes"]),
            "unannotatedGeneCalls": raw.get("unannotatedGeneCalls") or []}


def semantic(raw):
    from ..adapters.pgx_tsv import _compact_drug_annotations
    return {"genes": raw.get("genes"), "annotations": Counter(digest(a) for a in _compact_drug_annotations(raw)),
            "unannotated": raw.get("unannotatedGeneCalls") or []}


def reannotate(raw, choices, directory, run):
    """run(args) is the worker's bounded, cancellable process executor."""
    from .pgx_concordance import canonical
    native = check_original(raw)
    atomic(directory / "baseline.phenotype.json", native)
    reporter = lambda name: ["-reporter", "-ri", str(directory / f"{name}.phenotype.json"),
                            "-reporterJson", "-rt", str(raw.get("title") or "PGx"), "-o", str(directory), "-bf", name]
    run(reporter("baseline"))
    baseline = read(directory / "baseline.report.json")
    check_original(baseline)
    if semantic(baseline) != semantic(raw):
        raise ValueError("原用藥建議無法完整重現，暫不更新")
    with (directory / "override.tsv").open("w") as f:
        for gene, choice in sorted(choices.items()):
            gt = choice["genotype"]
            if gene not in raw["genes"] or any(c in gt for c in "\t\r\n"):
                raise ValueError("人工選擇無法轉換為 PharmCAT 輸入")
            f.write(f"{gene}\t{gt}\t\t\n")
    run(["-phenotyper", "-po", str(directory / "override.tsv"), "-o", str(directory), "-bf", "override"])
    overrides = read(directory / "override.phenotype.json")["geneReports"]
    for gene, choice in choices.items():
        new = overrides.get(gene) or {}
        sources = new.get("sourceDiplotypes") or []
        if (len(sources) != 1 or canonical(gene, sources[0].get("label")) != canonical(gene, choice["genotype"])
                or not sources[0].get("phenotypes")
                or any(str(p).lower() in {"unknown", "indeterminate", "no result"} for p in sources[0]["phenotypes"])
                or new.get("uncalledHaplotypes")):
            raise ValueError(f"{gene}：PharmCAT 無法辨識或註解此選擇，請恢復原結果")
        native["geneReports"][gene] = new
    atomic(directory / "revised.phenotype.json", native)
    run(reporter("revised"))
    result = read(directory / "revised.report.json")
    check_original(result)
    if set(result["genes"]) != set(raw["genes"]):
        raise ValueError("更新結果的基因範圍不符")
    for gene in raw["genes"]:
        expected = overrides[gene] if gene in choices else baseline["genes"][gene]
        # Reporter fills relatedDrugs; all call, phenotype, score and evidence
        # fields must survive exactly. Related drugs must be recomputed.
        before = {k: v for k, v in expected.items() if k != "relatedDrugs"}
        after = {k: v for k, v in result["genes"][gene].items() if k != "relatedDrugs"}
        if before != after:
            raise ValueError(f"{gene}：更新後的基因資料未通過核對")
    return result
