#!/usr/bin/env python3
"""DGX duplicate secondary-output review. Python 3.8+, standard library only.

Default: report only. --review additionally asks on /dev/tty before deleting.
Never infer a directory's creation date from ctime, mtime or its batch name.
"""
import argparse
from collections import defaultdict
from contextlib import contextmanager
import csv
from datetime import datetime, timedelta, timezone
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import stat
import subprocess
import sys
import tempfile
import uuid

DEFAULT_ROOT = "/datalake_Intermediate/pipeline/nextflow_output"
TAIPEI = timezone(timedelta(hours=8))
SAMPLE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]{0,127}$")


def within(path, root):
    try:
        Path(path).relative_to(root)
        return True
    except ValueError:
        return False


def birth_times(paths):
    found, pending = {}, []
    for path in paths:
        value = getattr(path.stat(), "st_birthtime", 0)
        if value > 0:
            found[str(path)] = value
        else:
            pending.append(path)
    if sys.platform.startswith("linux"):
        for start in range(0, len(pending), 200):
            batch = pending[start:start+200]
            try:
                proc = subprocess.run(["stat", r"--printf=%n\0%W\0", "--", *map(str, batch)],
                                      capture_output=True, check=False, timeout=30)
                parts = proc.stdout.split(b"\0")
                for i in range(0, len(parts)-1, 2):
                    name = os.fsdecode(parts[i])
                    try:
                        value = float(parts[i+1])
                        if name in set(map(str, batch)) and value > 0:
                            found[name] = value
                    except ValueError:
                        pass
            except (OSError, subprocess.TimeoutExpired):
                pass
    return found


def signature(path):
    st = path.lstat()
    return [st.st_dev, st.st_ino, st.st_mode, st.st_size, st.st_mtime_ns, st.st_ctime_ns]


def tree_state(path):
    """Metadata digest, not large BAM content hashing. Never follow links."""
    digest = hashlib.sha256()
    newest = 0
    device = path.stat().st_dev
    stack = [path]
    while stack:
        p = stack.pop()
        sig = signature(p)
        if sig[0] != device or (p != path and os.path.ismount(p)):
            raise RuntimeError("目錄包含另一個掛載點：" + str(p))
        target = os.readlink(p) if p.is_symlink() else ""
        digest.update(json.dumps([str(p.relative_to(path)), sig, target], ensure_ascii=True).encode())
        newest = max(newest, sig[4]/1e9, sig[5]/1e9)
        if stat.S_ISDIR(sig[2]):
            stack.extend(sorted(p.iterdir(), reverse=True))
    return digest.hexdigest(), newest


def complete(path, sid):
    def present(p):
        return p.is_file() and p.stat().st_size > 0
    bam = path / "02_alignment" / (sid + ".aligned.sorted.bam")
    vcf = path / "04_snv_indel" / (sid + ".ensemble.fixed.vcf.gz")
    return (present(bam) and any(present(p) for p in [Path(str(bam)+".bai"), bam.with_suffix(".bai")])
            and present(vcf) and any(present(Path(str(vcf)+s)) for s in (".tbi", ".csi")))


def core_state(path, sid):
    bam = path / "02_alignment" / (sid + ".aligned.sorted.bam")
    vcf = path / "04_snv_indel" / (sid + ".ensemble.fixed.vcf.gz")
    paths = [bam, Path(str(bam)+".bai"), bam.with_suffix(".bai"), vcf, Path(str(vcf)+".tbi"), Path(str(vcf)+".csi")]
    return [[str(p), str(p.resolve()), signature(p.resolve())] for p in paths if p.is_file()]


def samples(root):
    rows = []
    for batch in sorted(root.iterdir()):
        if batch.name.startswith(".") or batch.is_symlink() or not batch.is_dir():
            continue
        for path in sorted(batch.iterdir()):
            if path.is_symlink() or not path.is_dir() or not SAMPLE.fullmatch(path.name):
                continue
            if not any((path / stage).is_dir() for stage in ("02_alignment", "04_snv_indel")):
                continue
            rows.append(dict(sample=path.name, batch=batch.name, path=str(path)))
    dates = birth_times([Path(row["path"]) for row in rows])
    for row in rows:
        path = Path(row["path"])
        row.update(birth=dates.get(str(path)), modified=path.stat().st_mtime,
                   identity=signature(path)[:2], complete=complete(path, path.name),
                   outputs=core_state(path, path.name),
                   seq_type=next((t for t in ("WES", "WGS") if re.search(r"(?:^|[_-])"+t+r"(?:$|[_-])", row["batch"], re.I)), ""))
    return rows


def referenced_candidates(root, reference_root, candidates):
    blocked = defaultdict(list)
    candidate_paths = {str(Path(p)): Path(p) for p in candidates}

    def record(target, origin):
        text = str(target)
        # The same shared volume is mounted with /home on DGM, without it on DGX.
        if text.startswith("/home/datalake_Intermediate/"):
            text = text[len("/home"):]
        target = Path(os.path.normpath(text))
        if within(target, root):
            parts = target.relative_to(root).parts
            if len(parts) >= 2:
                candidate = candidate_paths.get(str(root / parts[0] / parts[1]))
                relevant = [candidate] if candidate else []
            else:
                relevant = [p for p in candidate_paths.values() if within(p, target)]
        elif within(root, target):
            relevant = list(candidate_paths.values())
        else:
            return
        for candidate in relevant:
            if within(target, candidate) or within(candidate, target):
                if not within(origin, candidate):
                    blocked[str(candidate)].append(str(origin))

    def onerror(exc):
        raise exc
    for base in (root, reference_root):
        if not base.is_dir():
            raise RuntimeError("無法核對相依資料根：" + str(base))
        for parent, directories, files in os.walk(str(base), followlinks=False, onerror=onerror):
            directories[:] = [d for d in directories if not d.startswith(".")]
            for name in directories + files:
                p = Path(parent) / name
                if p.is_symlink():
                    record(p.resolve(strict=False), p)
                if base == reference_root and p.parent.name == "08_postprocessing" and name.endswith("pipeline_source.json"):
                    if p.stat().st_size > 2000000:
                        raise RuntimeError("來源紀錄過大：" + str(p))
                    data = json.loads(p.read_text(encoding="utf-8"))
                    for key in ("source_path", "source_vcf_path"):
                        if isinstance(data.get(key), str) and data[key].startswith("/"):
                            record(data[key], p)
    return blocked


def audit(root, reference_root, cutoff):
    grouped = defaultdict(list)
    for row in samples(root):
        grouped[row["sample"]].append(row)
    groups, candidates = {}, []
    for sid, rows in sorted(grouped.items()):
        if len(rows) < 2:
            continue
        groups[sid] = sorted(rows, key=lambda r: (r["birth"] or 0, r["path"]), reverse=True)
        latest = max((r["birth"] or 0 for r in rows))
        uncertain = any(r["birth"] is None for r in rows)
        mixed = len({r["seq_type"] for r in rows if r["seq_type"]}) > 1
        latest_ready = all(r["complete"] for r in rows if r["birth"] == latest)
        needs_tree = not uncertain and not mixed and latest_ready and any(r["birth"] < min(latest, cutoff) for r in rows)
        for row in rows:
            digest, changed = tree_state(Path(row["path"])) if needs_tree else (None, 0)
            row.update(action="KEEP", reason="", fingerprint=digest)
            if uncertain:
                row["reason"] = "有版本建立時間未知，整組暫不刪"
            elif mixed:
                row["reason"] = "同ID同時有WES與WGS，整組暫不刪"
            elif row["birth"] == latest:
                row["reason"] = "最新版本（同時間並列亦保留）"
            elif row["birth"] >= cutoff:
                row["reason"] = "2026/11/5起的保護期間"
            elif not latest_ready:
                row["reason"] = "最新版本的BAM/VCF及index未齊，整組暫不刪"
            else:
                if changed >= cutoff:
                    row["reason"] = "保護期間仍有內容或狀態變動"
                else:
                    row.update(action="DELETE", reason="舊版本且早於保護期間")
                    candidates.append(row["path"])
    dependency_error = None
    if candidates:
        try:
            references = referenced_candidates(root, reference_root, candidates)
        except (OSError, ValueError, RuntimeError) as exc:
            references = {}
            dependency_error = str(exc)
        for rows in groups.values():
            for row in rows:
                if row["action"] == "DELETE" and dependency_error:
                    row.update(action="KEEP", reason="相依資料無法核對，暫不刪")
                elif row["path"] in references:
                    row.update(action="KEEP", reason="仍被其他資料引用", references=references[row["path"]][:20])
    return dict(root=str(root), reference_root=str(reference_root), cutoff=cutoff, groups=groups, dependency_error=dependency_error)


def token(plan):
    return hashlib.sha256(json.dumps(plan, sort_keys=True, ensure_ascii=True).encode()).hexdigest()[:12]


def stamp(value):
    return datetime.fromtimestamp(value, TAIPEI).strftime("%Y-%m-%d %H:%M:%S") if value else "無法取得"


def save_report(plan, destination):
    destination.mkdir(mode=0o700, parents=True, exist_ok=True)
    (destination / "plan.json").write_text(json.dumps(plan, ensure_ascii=False, indent=2)+"\n", encoding="utf-8")
    with (destination / "duplicates.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample", "建立時間(台北)", "修改時間(不是建立時間)", "動作", "原因", "路徑"])
        for sid, rows in plan["groups"].items():
            print("\nSample: " + sid)
            for row in rows:
                line = [sid, stamp(row["birth"]), stamp(row["modified"]), row["action"], row["reason"], row["path"]]
                writer.writerow(line)
                print("  {}  建立={}  修改={}  {}\n    {}".format(row["action"], line[1], line[2], row["reason"], row["path"]))
    count = sum(r["action"] == "DELETE" for rows in plan["groups"].values() for r in rows)
    if plan.get("dependency_error"):
        print("\n相依檢查未完成，本次全部保留：" + plan["dependency_error"])
    print("\n重複sample={}；待刪個案目錄={}；報告={}".format(len(plan["groups"]), count, destination))
    return count


def no_nextflow():
    result = subprocess.run(["ps", "-e", "-ww", "-o", "pid=,comm=,args="], capture_output=True, text=True, check=True)
    for line in result.stdout.splitlines():
        parts = line.split(None, 2)
        if len(parts) != 3:
            continue
        _, comm, args = parts
        executable = Path(comm).name.lower()
        if executable == "nextflow" or (executable == "java" and re.search(r"nextflow(?:\.cli\.Launcher|[^\s]*\.jar)", args)) or (
            executable in {"bash", "sh"} and re.search(r"(?:^|[/\s])(?:nextflow|\.command\.(?:run|sh))(?:\s|$)", args)
        ):
            raise RuntimeError("仍有Nextflow／工具工作，這次不刪：" + line.strip())


@contextmanager
def idle_guard(work_root):
    work_root = work_root.resolve(strict=True)
    directory = work_root.parent / ("." + work_root.name + ".ngs-ui-locks")
    if directory.is_symlink() or not directory.is_dir():
        raise RuntimeError("找不到既有共用工作鎖；未刪除：" + str(directory))
    handles = []
    try:
        for path in [directory / "admission.lock"] + sorted(p for p in directory.glob("*.lock") if p.name != "admission.lock"):
            fd = os.open(str(path), os.O_RDWR | getattr(os, "O_NOFOLLOW", 0))
            handle = os.fdopen(fd, "a+")
            try:
                fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BaseException:
                handle.close()
                raise
            if path.name == "admission.lock":
                handles.append(handle)
            else:
                # Admission remains held, so new cooperating jobs cannot
                # acquire a lease. Don't exhaust FDs on thousands of old scopes.
                handle.close()
        no_nextflow()
        yield
    finally:
        for handle in reversed(handles):
            handle.close()


def nested_mounts(path):
    # ismount does not recognize every same-filesystem bind mount.
    text = Path("/proc/self/mountinfo").read_text()
    for line in text.splitlines():
        raw = line.split()[4]
        decoded = re.sub(r"\\([0-7]{3})", lambda m: chr(int(m.group(1), 8)), raw)
        if within(Path(decoded), path):
            raise RuntimeError("待刪目錄包含掛載點：" + decoded)


def delete_reviewed(plan, report_dir, work_root):
    if not getattr(shutil.rmtree, "avoids_symlink_attacks", False):
        raise RuntimeError("此Python不支援安全rmtree，未刪除")
    root = Path(plan["root"])
    with idle_guard(work_root):
        current = audit(root, Path(plan["reference_root"]), plan["cutoff"])
        if token(current) != token(plan):
            raise RuntimeError("確認前後的檔案或相依狀態已變動；未刪除，請重新盤點")
        for rows in plan["groups"].values():
            for row in rows:
                if row["action"] != "DELETE":
                    continue
                path = Path(row["path"])
                if path.is_symlink() or not within(path, root) or len(path.relative_to(root).parts) != 2:
                    raise RuntimeError("不是可刪除的batch/sample目錄：" + str(path))
                nested_mounts(path)
                no_nextflow()
                for keeper in rows:
                    if keeper["birth"] == max(r["birth"] for r in rows):
                        kept = Path(keeper["path"])
                        if core_state(kept, keeper["sample"]) != keeper["outputs"] or tree_state(kept)[0] != keeper["fingerprint"]:
                            raise RuntimeError("最新保留版本已變動，停止刪除：" + str(kept))
                digest, _ = tree_state(path)
                if digest != row["fingerprint"] or signature(path)[:2] != row["identity"]:
                    raise RuntimeError("目錄已變動，停止刪除：" + str(path))
                quarantine = path.parent / (".ngs-secondary-cleanup-" + uuid.uuid4().hex)
                moved = quarantine / path.name
                event = dict(path=str(path), temporary_path=str(moved), status="staging", at=datetime.now(TAIPEI).isoformat())
                def log():
                    with (report_dir / "actions.jsonl").open("a", encoding="utf-8") as handle:
                        handle.write(json.dumps(event, ensure_ascii=False)+"\n"); handle.flush(); os.fsync(handle.fileno())
                log()
                quarantine.mkdir(mode=0o700)
                os.rename(str(path), str(moved))
                if moved.is_symlink() or signature(moved)[:2] != row["identity"]:
                    if not path.exists(): os.rename(str(moved), str(path))
                    event["status"] = "identity_changed_no_delete"; log()
                    raise RuntimeError("目錄身分改變，停止刪除：" + str(path))
                event["status"] = "deleting"; log()
                try:
                    shutil.rmtree(str(moved))
                    quarantine.rmdir()
                except BaseException:
                    event["status"] = "failed_remaining_in_temporary_path"; log()
                    raise
                event["status"] = "deleted"; log()
                print("DELETED " + str(path), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", default=DEFAULT_ROOT)
    parser.add_argument("--reference-root", help="DGX tertiary_output root; default sibling of nextflow_output")
    parser.add_argument("--work-root", default="/raid/DGM/work")
    parser.add_argument("--report-dir")
    parser.add_argument("--review", action="store_true", help="After reporting, ask on /dev/tty before any deletion")
    args = parser.parse_args()
    root = Path(args.root).resolve(strict=True)
    reference_root = Path(args.reference_root).resolve() if args.reference_root else root.parent / "tertiary_output"
    cutoff = datetime(2026, 11, 5, tzinfo=TAIPEI).timestamp()
    report_dir = Path(args.report_dir) if args.report_dir else Path(tempfile.mkdtemp(prefix="ngs-secondary-review-"))
    if within(report_dir.resolve(), root):
        raise RuntimeError("報告必須放在二級資料根以外")
    print("搜尋根：{}\n保護：2026-11-05 00:00:00 +08:00起建立，或之後仍有變動的資料\n只刪batch/sample，不刪整個batch。".format(root))
    plan = audit(root, reference_root, cutoff)
    count = save_report(plan, report_dir)
    if not args.review or not count:
        print("本次只有盤點，沒有刪除。")
        return 0
    phrase = "DELETE " + token(plan)
    print("\n請確認DGX與DGM的分析工作均未執行，清理期間也不啟動新分析。")
    print("確認只永久刪除以上標記DELETE的{}個目錄。輸入 {} 才執行；Enter或其他輸入一律不刪。".format(count, phrase), flush=True)
    with open("/dev/tty", "r+") as terminal:
        terminal.write("確認："); terminal.flush()
        answer = terminal.readline().strip()
    if answer != phrase:
        print("未確認，沒有刪除。")
        return 0
    delete_reviewed(plan, report_dir, Path(args.work_root))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError, RuntimeError, subprocess.SubprocessError) as exc:
        print("STOP: " + str(exc), file=sys.stderr)
        raise SystemExit(2)
