"""One atomic discovery generation for FASTQs, pipeline inputs and file lists.

Only the configured FASTQ roots establish inventory samples. Pipeline views
retain VCF-only cases. Analysis eligibility still belongs to the existing
secondary/tertiary parsers, fed the same single directory walk.
"""
from __future__ import annotations

from collections import defaultdict
from contextlib import contextmanager
from datetime import datetime, timedelta, timezone
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import sqlite3
import stat
import subprocess
import sys
import tempfile
import time
import uuid

from .. import config

SCHEMA = 1
TAIPEI = timezone(timedelta(hours=8))
READ_PAIR = re.compile(r"^(.+?)(?:_S\d+)?(?:_L\d{3})?[_.]R([12])(?:_001|\.clean|_merged)?\.(?:fastq|fq)\.gz$", re.I)
OTHER = {
    ".cnv.vcf.gz": "cnv", ".gcnv.vcf.gz": "cnv", ".sv.vcf.gz": "sv",
    ".delly.vcf.gz": "sv", ".cnv_sv.vcf.gz": "cnv_sv",
    ".hard-filtered.gvcf.gz": "gvcf", ".gvcf.gz": "gvcf", ".g.vcf.gz": "gvcf",
    ".targeted.json": "targeted",
}
SUFFIXES = {**OTHER, ".hard-filtered.vcf.gz": "vcf", ".ensemble.fixed.vcf.gz": "vcf",
            ".aligned.sorted.bam.bai": "bam", ".aligned.sorted.bai": "bam",
            ".aligned.sorted.bam": "bam", ".bam.bai": "bam", ".bai": "bam", ".bam": "bam"}


def identity(path):
    return hashlib.sha256(str(path).encode()).hexdigest()


def data_roots():
    from . import secondary_analysis as secondary, dragen_jobs as tertiary
    paths = lambda key, fallback: [Path(p) for p in os.environ.get(key, fallback).split(":") if p]
    return {
        "wes": list(secondary.SECONDARY_WES_FASTQ_ROOTS),
        "wgs": list(secondary.SECONDARY_WGS_FASTQ_ROOTS),
        "dragen": list(tertiary.DRAGEN_VCF_ROOTS),
        "inhouse": list(tertiary.INHOUSE_VCF_ROOTS),
        "dragen_bam": paths("NGS_UI_DRAGEN_BAM_ROOT", "/home/datalake_Raw/Novaseq:/home/datalake_Intermediate/n102968"),
        "inhouse_bam": paths("NGS_UI_INHOUSE_BAM_ROOT", os.environ.get("NGS_UI_BAM_ROOT", "/home/datalake_Intermediate/pipeline/nextflow_output")),
    }


def tertiary_roots():
    # Inventory/browsing deliberately excludes both legacy roots. Other
    # clinical consumers keep their existing layout compatibility rules.
    return [config.PIPELINE_OUT_ROOT.resolve()]


@contextmanager
def connection():
    db = None
    try:
        db = sqlite3.connect(config.FILE_INDEX_PATH.resolve().as_uri() + "?mode=ro", uri=True)
        db.row_factory = sqlite3.Row
        if db.execute("PRAGMA user_version").fetchone()[0] != SCHEMA:
            raise sqlite3.DatabaseError("unsupported file index schema")
        yield db
    finally:
        if db is not None:
            db.close()


def load_view(name):
    try:
        with connection() as db:
            row = db.execute("SELECT payload FROM views WHERE name=?", (name,)).fetchone()
            return json.loads(row[0]) if row else None
    except (OSError, sqlite3.Error, ValueError):
        return None


def metadata():
    return load_view("meta")


def settings_signature():
    settings = {k: [str(p.resolve()) for p in paths] for k, paths in data_roots().items()}
    settings["tertiary"] = [str(p) for p in tertiary_roots()]
    return identity(json.dumps(settings, sort_keys=True))


def stale(meta):
    from . import secondary_analysis as secondary
    ttl = min(secondary.SECONDARY_FASTQ_INDEX_TTL_HOURS, config.PIPELINE_VCF_INDEX_TTL_HOURS)
    return (not meta or meta.get("settings_signature") != settings_signature()
            or time.time_ns() - meta.get("completed_at_ns", 0) >= ttl * 3600 * 1e9)


def _within(path, roots):
    return any(path.is_relative_to(root) for root in roots)


def discover(roots, previous):
    """Walk overlapping physical roots once; ignore irrelevant test FASTQs."""
    resolved = {key: [p.resolve() for p in paths] for key, paths in roots.items()}
    all_roots = sorted({p for paths in resolved.values() for p in paths})
    files, seen_dirs, errors, skipped, unknown = {}, set(), [], [], []
    old_roots = set((previous or {}).get("scanned_roots", []))
    optional = {Path("/datalake_Raw/Other/Reanalysis")}

    def visit(folder):
        try:
            st = folder.stat()
            key = (st.st_dev, st.st_ino)
            if key in seen_dirs:
                return
            seen_dirs.add(key)
            with os.scandir(folder) as entries:
                for entry in entries:
                    if entry.name.startswith("."):
                        continue
                    p = Path(entry.path)
                    try:
                        if entry.is_dir(follow_symlinks=False):
                            visit(p)
                            continue
                        # Directory links are not followed: explicit roots can
                        # point at their resolved target without traversing cycles.
                        if entry.is_symlink() and p.is_dir():
                            skipped.append(str(p))
                            continue
                        fastq = entry.name.lower().endswith((".fastq.gz", ".fq.gz"))
                        if fastq:
                            if not _within(p, resolved["wes"] + resolved["wgs"]):
                                continue
                        elif not any(entry.name.endswith(s) for s in SUFFIXES):
                            continue
                        if entry.name.endswith((".repeats.bam", ".repeats.bam.bai", ".repeats.bai")):
                            continue
                        st = p.stat()
                        if not stat.S_ISREG(st.st_mode):
                            continue
                        files[str(p)] = (p, st)
                        if fastq and not READ_PAIR.match(entry.name):
                            unknown.append(str(p))
                    except OSError as exc:
                        errors.append(f"{p}: {exc}")
        except OSError as exc:
            errors.append(f"{folder}: {exc}")

    scanned = []
    for root in all_roots:
        if not root.exists() and root in optional and str(root) not in old_roots:
            continue
        visit(root)
        scanned.append(str(root))
    if errors:
        raise OSError("檔案索引未更新，保留上次結果：" + "\n".join(errors[:20]))
    return files, resolved, scanned, skipped, unknown


def build(previous=None):
    from . import secondary_analysis as secondary, dragen_jobs as tertiary
    start = time.monotonic()
    found, roots, scanned, skipped, unknown = discover(data_roots(), previous)
    paths = [p for p, _ in found.values()]
    # Existing parser behavior remains authoritative for runnable inputs.
    wes = secondary.list_wes_fastqs(paths=paths)
    wgs = secondary.list_wgs_fastqs(paths=paths)
    dragen = tertiary.list_dragen_vcfs(paths=paths)
    inhouse = tertiary.list_inhouse_vcfs(paths=paths)
    samples, aliases, registry = {}, defaultdict(set), {}

    def sample(sid, source):
        row = samples.setdefault(sid, dict(sample_id=sid, fastq=[], bam=[], vcf=[], other=[], tertiary=[]))
        aliases[source].add(sid)
        aliases[sid].add(sid)
        return row

    # Include single-ended/raw-only cases; samplesheet view still requires pairs.
    for path, st in found.values():
        if not path.name.lower().endswith((".fastq.gz", ".fq.gz")):
            continue
        match = READ_PAIR.match(path.name)
        if not match:
            continue
        source = match.group(1)
        reanalysis = secondary._REANALYSIS_RE.match(path.name)
        sid = path.parent.name if reanalysis else source
        row = sample(sid, source)
        entry = dict(id=identity(path), name=path.name, path=str(path))
        row["fastq"].append(entry)
        registry[entry["id"]] = str(path)

    for path, st in found.values():
        suffix = next((s for s in SUFFIXES if path.name.endswith(s)), None)
        if not suffix:
            continue
        source, kind = path.name[:-len(suffix)], SUFFIXES[suffix]
        # Only requested files/layouts; in-house never offers gVCF/targeted/joint.
        is_dragen = _within(path, roots["dragen"] + roots["dragen_bam"])
        is_inhouse = _within(path, roots["inhouse"] + roots["inhouse_bam"])
        if not (is_dragen or is_inhouse):
            continue
        if kind in {"gvcf", "targeted", "cnv_sv"} and not is_dragen:
            continue
        if kind == "bam" and path.parent.name not in {"bam", "02_alignment"}:
            continue
        if kind == "vcf" and path.parent.name not in {"vcf.gz", "04_snv_indel"}:
            continue
        column = kind if kind in {"bam", "vcf"} else "other"
        for sid in aliases.get(source, ()):
            entry = dict(id=identity(path), name=path.name, path=str(path))
            samples[sid][column].append(entry)
            registry[entry["id"]] = str(path)

    # One-level tertiary enumeration only. Its arbitrary contents are browsed
    # lazily, never recursively indexed as part of the daily raw/input scan.
    for root in tertiary_roots():
        if not root.exists():
            if root == config.PIPELINE_OUT_ROOT.resolve() or str(root) in (previous or {}).get("tertiary_roots", []):
                raise OSError(f"無法讀取三級資料根：{root}")
            continue
        for folder in root.iterdir():
            if folder.name.startswith(".") or not folder.is_dir() or not folder.resolve().is_relative_to(root):
                continue
            sid = folder.name
            marker = folder / "08_postprocessing" / f"{sid}.pipeline_source.json"
            if not marker.exists():
                marker = folder / "08_postprocessing" / "pipeline_source.json"
            source = sid
            if marker.is_file():
                if not marker.resolve().is_relative_to(root) or marker.stat().st_size > 2_000_000:
                    continue
                value = json.loads(marker.read_text())
                if not isinstance(value, dict):
                    raise ValueError(f"三級來源紀錄格式錯誤：{marker}")
                source = value.get("source_sample_id") or sid
                if not isinstance(source, str):
                    raise ValueError(f"三級來源sample ID格式錯誤：{marker}")
            targets = set(aliases.get(source, ())) | set(aliases.get(sid, ()))
            if not marker.is_file():
                # Exact names emitted by our own tertiary worker, no broad
                # prefix matching (e.g. S1 must never match S10).
                for mode in ("dragen", "inhouse"):
                    original = sid.removesuffix("-dragen" if mode == "dragen" else "-nckuh")
                    if tertiary._default_ui_sample_id(original, mode) == sid:
                        targets.update(aliases.get(original, ()))
            for target in targets:
                samples[target]["tertiary"].append(dict(name=sid, root=identity(root), path=sid))
                annotation = _annotation_entry(root, folder, source)
                if annotation and not any(f["id"] == annotation["id"] for f in samples[target]["vcf"]):
                    samples[target]["vcf"].append(annotation)
                    registry[annotation["id"]] = annotation["path"]

    for row in samples.values():
        for column in ("fastq", "bam", "vcf", "other"):
            row[column].sort(key=lambda f: (str(Path(f["path"]).parent), f["name"]))
    meta = dict(generation=uuid.uuid4().hex, updated_at=datetime.now(TAIPEI).strftime("%Y-%m-%d %H:%M:%S"),
                completed_at_ns=time.time_ns(), settings_signature=settings_signature(),
                scan_duration_sec=round(time.monotonic()-start, 2),
                sample_count=len(samples), scanned_roots=scanned,
                tertiary_roots=[str(p) for p in tertiary_roots() if p.is_dir()],
                unrecognized_fastq_count=len(unknown), unrecognized_fastq_examples=unknown[:30],
                skipped_directory_links=skipped[:30])
    return dict(meta=meta, secondary=dict(**meta, schema_version=secondary._FASTQ_INDEX_SCHEMA_VERSION, wes=wes, wgs=wgs),
                tertiary=dict(**meta, dragen=dragen, inhouse=inhouse)), samples, registry


def refresh(*, force=False, requested_at_ns=None):
    """Coalesce overlapping processes; publish a complete DB via atomic replace."""
    requested = requested_at_ns or time.time_ns()
    observed = metadata()
    observed_generation = observed.get("generation") if observed else None
    path = config.FILE_INDEX_PATH
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.with_suffix(".refresh.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        old = metadata()
        if old and old.get("settings_signature") == settings_signature() and (
            old["generation"] != observed_generation or old["completed_at_ns"] >= requested
            or (not force and not stale(old))
        ):
            return {name: load_view(name) for name in ("meta", "secondary", "tertiary")}
        views, samples, registry = build(old)
        fd, temp = tempfile.mkstemp(prefix=".sample-files-", suffix=".sqlite", dir=path.parent)
        os.close(fd)
        try:
            db = sqlite3.connect(temp)
            try:
                db.executescript(f"""PRAGMA user_version={SCHEMA};
                    CREATE TABLE views(name TEXT PRIMARY KEY, payload TEXT NOT NULL);
                    CREATE TABLE samples(sample_id TEXT PRIMARY KEY, payload TEXT NOT NULL);
                    CREATE TABLE files(id TEXT PRIMARY KEY, path TEXT NOT NULL);""")
                db.executemany("INSERT INTO views VALUES (?,?)", [(n, json.dumps(v, ensure_ascii=False)) for n,v in views.items()])
                db.executemany("INSERT INTO samples VALUES (?,?)", [(s, json.dumps(v, ensure_ascii=False)) for s,v in samples.items()])
                db.executemany("INSERT INTO files VALUES (?,?)", registry.items())
                db.commit()
            finally:
                db.close()
            with open(temp, "rb") as complete:
                os.fsync(complete.fileno())
            os.replace(temp, path)
        finally:
            Path(temp).unlink(missing_ok=True)
        return views


def list_samples(query="", offset=0, limit=50):
    with connection() as db:
        pattern = "%" + query.replace("\\", "\\\\").replace("%", "\\%").replace("_", "\\_") + "%"
        count = db.execute("SELECT count(*) FROM samples WHERE sample_id LIKE ? ESCAPE '\\'", (pattern,)).fetchone()[0]
        rows = db.execute("SELECT payload FROM samples WHERE sample_id LIKE ? ESCAPE '\\' ORDER BY sample_id DESC LIMIT ? OFFSET ?",
                          (pattern, limit, offset)).fetchall()
        meta = json.loads(db.execute("SELECT payload FROM views WHERE name='meta'").fetchone()[0])
        samples = [json.loads(row[0]) for row in rows]
    _tertiary_details(samples)
    return dict(meta=meta, total=count, offset=offset, limit=limit, samples=samples)


def _native_birth_time(path):
    return getattr(path.stat(), "st_birthtime", None)


def _source_record(folder):
    post = folder / "08_postprocessing"
    path = post / f"{folder.name}.pipeline_source.json"
    if not path.exists():
        path = post / "pipeline_source.json"
    try:
        if not path.resolve().is_relative_to(folder.resolve()) or path.stat().st_size > 2_000_000:
            return {}
        value = json.loads(path.read_text(encoding="utf-8"))
        return value if isinstance(value, dict) else {}
    except (OSError, ValueError):
        return {}


def _annotation_entry(root, folder, source=None):
    # Exact source filename first. No wildcard or sample-prefix matching.
    candidates = [source] if source else []
    candidates += [folder.name.removesuffix(s) for s in ("-dragen", "-nckuh", "-inhouse")]
    for sid in dict.fromkeys(s for s in candidates if isinstance(s, str)):
        if sid in {"", ".", ".."} or Path(sid).name != sid:
            continue
        path = folder / "03_acmg" / f"{sid}.snv_indel.acmg.tsv"
        try:
            if path.is_file() and path.resolve().is_relative_to(root):
                return dict(id=identity(path), name=path.name, path=str(path),
                            tertiary_root=identity(root), relative_path=str(path.relative_to(root)))
        except OSError:
            continue
    return None


def _recorded_annotation_time(value):
    if not isinstance(value, str):
        return None
    try:
        stamp = datetime.fromisoformat(value.replace("Z", "+00:00"))
        if stamp.tzinfo is None:
            stamp = stamp.replace(tzinfo=TAIPEI)
        return stamp.astimezone(TAIPEI)
    except (ValueError, OverflowError):
        return None


def _birth_times(paths):
    """True filesystem birth time only; ctime/mtime are never substitutes.

    Python exposes st_birthtime on macOS. GNU stat uses statx on Linux;
    batch only the requested page's folders, so existing snapshots work
    immediately and no complete raw-data rescan is needed.
    """
    result, pending = {}, []
    for path in paths:
        try:
            stamp = _native_birth_time(path)
            if stamp is not None and stamp > 0:
                result[path] = stamp
            else:
                pending.append(path)
        except OSError:
            continue
    if pending and sys.platform.startswith("linux"):
        try:
            proc = subprocess.run(
                ["stat", "--dereference", r"--printf=%n\0%W\0", "--", *map(str, pending)],
                capture_output=True, check=False, timeout=5,
            )
            fields = proc.stdout.split(b"\0")
            allowed = set(pending)
            for i in range(0, len(fields)-1, 2):
                path = Path(os.fsdecode(fields[i]))
                try:
                    stamp = float(fields[i+1])
                    if path in allowed and stamp > 0:
                        result[path] = stamp
                except ValueError:
                    continue
        except (OSError, subprocess.TimeoutExpired):
            pass
    return result


def _tertiary_details(samples):
    roots = {identity(root): root for root in tertiary_roots()}
    folders = defaultdict(list)
    for sample in samples:
        # Drop old cached annotation links and re-add exact, currently present
        # files from the allowed primary folder. No full index rebuild needed.
        sample["vcf"] = [f for f in sample.get("vcf", []) if not f.get("tertiary_root")]
        valid_folders = []
        for folder in sample.get("tertiary", []):
            base = roots.get(folder.get("root"))
            relative = Path(folder.get("path", ""))
            if base is None or relative.is_absolute() or any(p.startswith(".") for p in relative.parts):
                continue
            absolute = base / relative
            try:
                resolved = absolute.resolve()
                if not resolved.is_relative_to(base):
                    continue
            except OSError:
                continue
            folder.update(absolute_path=str(absolute), created_date=None, created_at=None,
                          display_date=None, date_at=None, date_source=None)
            record = _source_record(absolute)
            annotated = _recorded_annotation_time(record.get("annotated_at"))
            if annotated:
                folder.update(display_date=annotated.strftime("%Y-%m-%d"),
                              date_at=annotated.isoformat(timespec="seconds"), date_source="tertiary_annotation")
            annotation = _annotation_entry(base, absolute, record.get("source_sample_id"))
            if annotation and not any(f["id"] == annotation["id"] for f in sample["vcf"]):
                sample["vcf"].append(annotation)
            folders[resolved].append(folder)
            valid_folders.append(folder)
        # This also removes legacy links already present in an older DB.
        sample["tertiary"] = valid_folders
    for path, stamp in _birth_times(folders).items():
        try:
            created = datetime.fromtimestamp(stamp, TAIPEI)
        except (ValueError, OverflowError, OSError):
            continue
        for folder in folders[path]:
            folder.update(created_date=created.strftime("%Y-%m-%d"), created_at=created.isoformat(timespec="seconds"),
                          display_date=created.strftime("%Y-%m-%d"), date_at=created.isoformat(timespec="seconds"),
                          date_source="filesystem_birth")


def download_path(file_id):
    with connection() as db:
        row = db.execute("SELECT path FROM files WHERE id=?", (file_id,)).fetchone()
    if not row:
        raise FileNotFoundError("找不到檔案")
    return Path(row[0])


def open_bounded(path, roots, *, directory=False):
    """Resolve allowed links, then open every component without following links.

    Returned descriptor remains bound even if the indexed path is replaced.
    Reject special files and parent-component symlink races before reading.
    """
    resolved = Path(path).resolve(strict=True)
    for base in roots:
        base = Path(base).resolve()
        if not resolved.is_relative_to(base):
            continue
        flags = os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK
        fd = os.open(base, flags | os.O_DIRECTORY)
        try:
            parts = resolved.relative_to(base).parts
            for i, part in enumerate(parts):
                new = os.open(part, flags | (os.O_DIRECTORY if directory or i < len(parts)-1 else 0), dir_fd=fd)
                os.close(fd)
                fd = new
            st = os.fstat(fd)
            if not (stat.S_ISDIR(st.st_mode) if directory else stat.S_ISREG(st.st_mode)):
                raise PermissionError("不支援此檔案類型")
            return fd, st
        except BaseException:
            os.close(fd)
            raise
    raise PermissionError("路徑超出可讀範圍")
