"""Detached PGx worker; only a final review pointer publishes a revision."""
import fcntl
import hashlib
import json
import os
import signal
import subprocess
import sys
import time
import traceback

from ..services import pgx_reannotation as jobs, pgx_concordance as pc, sample_layout
from ..services.sample_loader import load_sample_pgx
from .dragen_run import _acquire_sample_locks, _release_sample_locks


class Cancelled(Exception):
    pass


class Worker:
    def __init__(self, jid):
        self.directory = jobs.job_dir(jid)
        self.job = jobs.read(self.directory / "state.json")
        self.sid = self.job["sample_id"]
        self.process = None

    def update(self, **values):
        self.job.update(values, updated=time.time())
        jobs.atomic(self.directory / "state.json", self.job)

    def check(self):
        if (self.directory / "restore.json").exists():
            raise Cancelled("已恢復原結果")

    def stop(self):
        process, self.process = self.process, None
        if process is None:
            return
        try:
            os.killpg(process.pid, signal.SIGTERM)
            try:
                process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                pass
            # A container launcher can exit before a descendant. Reap the
            # remaining owned group even when the immediate child has exited.
            os.killpg(process.pid, signal.SIGKILL)
            process.wait(timeout=5)
        except ProcessLookupError:
            process.wait(timeout=5)

    def command(self, args):
        self.check()
        with (self.directory / "pharmcat.log").open("ab") as log:
            env = os.environ.copy()
            # No inherited site binds or Java options can substitute the pinned tool.
            for name in list(env):
                if name.startswith(("APPTAINER_", "SINGULARITY_")) or name in {"JAVA_TOOL_OPTIONS", "_JAVA_OPTIONS", "JDK_JAVA_OPTIONS"}:
                    env.pop(name, None)
            env["TMPDIR"] = str(self.directory / "tmp")
            self.process = subprocess.Popen(jobs.command(self.engine, self.directory, args),
                cwd=self.directory, env=env, stdout=log, stderr=subprocess.STDOUT,
                stdin=subprocess.DEVNULL, start_new_session=True)
            started = time.monotonic()
            try:
                while self.process.poll() is None:
                    self.check()
                    if time.monotonic() - started > 300:
                        raise ValueError("PharmCAT 更新逾時，請重試")
                    time.sleep(0.2)
                if self.process.returncode:
                    raise ValueError("PharmCAT 更新失敗，請管理者檢查執行紀錄")
            finally:
                self.stop()

    def verify(self):
        self.check()
        if not sample_layout.state_file(self.sid, "sample_metadata.json").is_file():
            raise ValueError("個案已取消登錄，無法更新 PGx")
        review = pc._review(self.sid)
        pgx = load_sample_pgx(self.sid)["pgx"]
        if (review.get("annotation_job") != self.job["job_id"]
                or jobs.request_key(pgx, review["choices"]) != self.job["request_key"]
                or any(pgx["concordance"]["genes"].get(g, {}).get("stale", True) for g in self.job["choices"])):
            raise ValueError("比對來源或人工選擇已更新，請重新確認")
        return review, pgx

    def run(self):
        handles = []
        def interrupted(*_):
            raise Cancelled("更新中斷，請重試")
        signal.signal(signal.SIGTERM, interrupted)
        signal.signal(signal.SIGINT, interrupted)
        with (self.directory / "lease").open("a+") as lease:
            fcntl.flock(lease, fcntl.LOCK_EX | fcntl.LOCK_NB)
            try:
                # save_choice may still hold its short admission lock at spawn.
                for attempt in range(50):
                    self.check()
                    try:
                        handles = _acquire_sample_locks([self.sid, self.job["source_sample_id"]])
                        break
                    except RuntimeError:
                        if attempt == 49:
                            raise
                        time.sleep(0.1)
                self.update(status="running")
                with (jobs.root() / "slot.lock").open("a+") as slot:
                    while True:
                        self.check()
                        try:
                            fcntl.flock(slot, fcntl.LOCK_EX | fcntl.LOCK_NB)
                            break
                        except BlockingIOError:
                            time.sleep(0.2)
                    _, pgx = self.verify()
                    original = sample_layout.pharmcat_json(self.sid).read_bytes()
                    if hashlib.sha256(original).hexdigest() != pgx.get("source_identity", {}).get("report_sha256"):
                        raise ValueError("PharmCAT 原報告已變更")
                    raw = json.loads(original)
                    self.engine = jobs.runtime()
                    (self.directory / "tmp").mkdir(exist_ok=True)
                    jobs.atomic(self.directory / "original.report.json", raw)
                    result = jobs.reannotate(raw, self.job["choices"], self.directory, self.command)
                    # Revalidate the actual program resource immediately before publishing.
                    if jobs.runtime() != self.engine:
                        raise ValueError("PharmCAT 執行設定已變更")
                    with jobs.review_lock(self.sid) as path:
                        review, _ = self.verify()
                        revision = jobs.revision_dir(self.sid, self.job["job_id"])
                        revision.mkdir(parents=True, exist_ok=False)
                        jobs.atomic(revision / "report.json", result)
                        jobs.atomic(revision / "manifest.json", {
                            "schema": jobs.SCHEMA, "sample_id": self.sid, "engine_sha256": jobs.JAR_SHA256,
                            "request_key": self.job["request_key"], "choices": self.job["choices"],
                            "baseline_replay_verified": True,
                            "report_sha256": hashlib.sha256((revision / "report.json").read_bytes()).hexdigest(),
                            "original_report_sha256": hashlib.sha256(original).hexdigest(), "created": time.time(),
                        })
                        pointer = {"job_id": self.job["job_id"], "request_key": self.job["request_key"]}
                        review["published"] = pointer
                        review.setdefault("revisions", []).append(pointer)
                        review["audit"].append({"action": "annotation_published", **pointer, "time": time.time()})
                        jobs.atomic(path, review)
                        self.update(status="complete", error="")
            except (Exception, SystemExit) as exc:
                self.stop()
                traceback.print_exc()
                restored = False
                if (self.directory / "restore.json").exists():
                    with jobs.review_lock(self.sid) as path:
                        review = pc._review(self.sid)
                        if review.get("annotation_job") == self.job["job_id"]:
                            review["choices"] = {}
                            review.pop("published", None)
                            review.pop("annotation_job", None)
                            review["audit"].append({"action": "restore", **jobs.read(self.directory / "restore.json")})
                            jobs.atomic(path, review)
                            restored = True
                error = str(exc) if isinstance(exc, (ValueError, Cancelled)) else "PGx 更新失敗，請管理者檢查執行紀錄"
                if isinstance(exc, RuntimeError):
                    error = "個案正由其他分析使用，請稍後重試"
                self.update(status="cancelled" if isinstance(exc, Cancelled) else "failed", error="" if restored else error)
            finally:
                self.stop()
                _release_sample_locks(handles)


if __name__ == "__main__":
    Worker(sys.argv[1]).run()
