"""HPO database update: download → clean (phase 0) → knowledge base (phase 1) → activate.

Used by the home-page "更新 HPO" button (detached worker, same pattern as the
LitVar2 updater) and by ``scripts/build_hpo_kb.py`` /
``scripts/hpo/build_hpo_gene_table.py``.

The running web service picks the new table and hp.obo up by itself
(``phenotype_scorer`` / ``hpo_ontology`` check the files' signatures), so no
restart is needed.  Each case's pheno_score.tsv is recomputed when it is next
opened (table stamp, v10.51).
"""
from __future__ import annotations

import fcntl
import json
import os
import shutil
import subprocess
import sys
import time
import traceback
import urllib.request
import uuid
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path

from ..config import JOBS_DIR, PHENO_DATA_DIR, REPO_ROOT
from . import hpo_gene_table as H

RELEASE_URL = "https://github.com/obophenotype/human-phenotype-ontology/releases/{part}/{name}"
RELEASE_FILES = ("hp.obo", "phenotype.hpoa", "genes_to_disease.txt")
# orphadata.com first; Orphanet's own GitHub copy (same product6 file, may lag
# one release) when the server cannot reach orphadata.com.
ORPHANET_URLS = ("https://www.orphadata.com/data/xml/en_product6.xml",
                 "https://raw.githubusercontent.com/Orphanet/orphanet_nomenclature_binder/main/xml/en_product6.xml")
MONDO_URLS = ("https://github.com/monarch-initiative/mondo/releases/latest/download/mondo.obo",
              "https://purl.obolibrary.org/obo/mondo.obo")
OPTIONAL_MAX_AGE_DAYS = 30
PICKER_DATA = "hpo_data.json"           # the /phenotype/ tool's term list

STATE_DIR = JOBS_DIR / "hpo_update"
STATE_PATH = STATE_DIR / "state.json"
LOG_PATH = STATE_DIR / "update.log"


def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


# ── downloads ──────────────────────────────────────────────────────────────
def download(url: str, dest: Path, *, attempts: int = 3, timeout: int = 120) -> None:
    tmp = dest.with_name(dest.name + ".part")
    dest.parent.mkdir(parents=True, exist_ok=True)
    for attempt in range(attempts):
        try:
            request = urllib.request.Request(url, headers={"User-Agent": "ngs-ui-hpo/1.0"})
            with urllib.request.urlopen(request, timeout=timeout) as response, open(tmp, "wb") as out:
                shutil.copyfileobj(response, out)
            tmp.replace(dest)
            return
        except OSError:
            tmp.unlink(missing_ok=True)
            if attempt == attempts - 1:
                raise
            time.sleep(2 ** (attempt + 1))


def fetch_release(version: str, data_dir: Path = PHENO_DATA_DIR, *, download=download, progress=print) -> Path:
    """Download one HPO release's three files; returns the directory holding them."""
    version = (version or "latest").strip().lstrip("v")
    if version == "latest":
        staging = data_dir / "hpo_release_latest.part"
        staging.mkdir(parents=True, exist_ok=True)
        for name in RELEASE_FILES:
            progress(f"下載 HPO latest/{name} …")
            download(RELEASE_URL.format(part="latest/download", name=name), staging / name)
        version = H._date_version(H.parse_obo(staging / "hp.obo").version)
        final = data_dir / f"hpo_release_{version}"
        if final.exists():
            shutil.rmtree(staging)
        else:
            staging.rename(final)
        return final
    target = data_dir / f"hpo_release_{version}"
    target.mkdir(parents=True, exist_ok=True)
    for name in RELEASE_FILES:
        if not (target / name).is_file():
            progress(f"下載 HPO v{version}/{name} …")
            download(RELEASE_URL.format(part=f"download/v{version}", name=name), target / name)
    return target


def install_obo(source: Path, data_dir: Path = PHENO_DATA_DIR) -> Path | None:
    """Install a release's hp.obo as the live file (old one kept as hp.obo.bak-time)."""
    live = data_dir / "hp.obo"
    if live.is_file() and live.read_bytes() == source.read_bytes():
        return None
    if live.is_file():
        shutil.copy2(live, live.with_name(f"hp.obo.bak-{time.strftime('%Y%m%d%H%M%S')}"))
    tmp = live.with_name(".hp.obo.part")
    shutil.copy2(source, tmp)
    tmp.replace(live)
    return live


def fetch_optional(label: str, urls, dest: Path, *, download=download, progress=print,
                   max_age_days: int = OPTIONAL_MAX_AGE_DAYS) -> tuple[Path | None, str]:
    """Refresh an optional source when older than ``max_age_days``.

    Returns (file to use or None, warning).  A failed download keeps the
    existing copy."""
    if dest.is_file() and time.time() - dest.stat().st_mtime < max_age_days * 86400:
        return dest, ""
    errors = []
    for url in urls:
        try:
            progress(f"下載 {label}：{url} …")
            download(url, dest)
            return dest, ""
        except OSError as exc:
            errors.append(f"{url}：{getattr(exc, 'reason', None) or exc}")
    if dest.is_file():
        return dest, f"{label} 下載失敗，沿用本機舊檔（{time.strftime('%Y-%m-%d', time.localtime(dest.stat().st_mtime))}）"
    return None, f"{label} 下載失敗（{'；'.join(errors)}）"


# ── /phenotype/ picker data ────────────────────────────────────────────────
def write_picker_data(hp_obo: Path, gene_table: Path, out: Path) -> int:
    """The /phenotype/ tool's term list (id, name, synonyms, first parent, gene counts)."""
    from . import hpo_ontology
    terms = hpo_ontology._parse_obo(hp_obo)
    counts = {term: len(genes) for term, genes in H.read_gene_table(gene_table).items()}
    rows = []
    for hid in sorted(terms):
        t = terms[hid]
        row = {"id": hid, "name": t.name}
        if t.synonyms:
            row["syn"] = list(t.synonyms)
        parent = next((p for p in t.parent_ids if p in terms), "")
        if parent:
            row["par"] = parent
        row["g"] = counts.get(hid, 0)
        if parent:
            row["pn"] = terms[parent].name
            row["pg"] = counts.get(parent, 0)
        row["n"] = int(hid.split(":")[1])
        rows.append(row)
    tmp = out.with_name(out.name + ".part")
    tmp.write_text(json.dumps(rows, ensure_ascii=False, separators=(",", ":")), encoding="utf-8")
    tmp.replace(out)
    return len(rows)


# ── build ──────────────────────────────────────────────────────────────────
def build_snapshot(hp_obo: Path, hpoa: Path, genes_to_disease: Path, *, kb_root: Path | None = None,
                   orphanet: Path | None = None, mondo: Path | None = None, allow_missing=(),
                   with_kb: bool = True, omim=None, progress=print) -> tuple[Path, list[str], dict]:
    """Clean HPO (phase 0) and, when possible, add the knowledge base (phase 1).

    Returns (snapshot, warnings, manifest).  A knowledge-base QC failure aborts
    the whole snapshot; a missing optional source becomes a warning when it is
    listed in ``allow_missing``."""
    from .pheno_engine import build as B
    warnings: list[str] = []
    progress("清理 HPO 註記（階段 0）…")
    result = H.build(hp_obo, hpoa, genes_to_disease)
    extra = None
    if with_kb:
        try:
            import numpy  # noqa: F401
        except ImportError:
            warnings.append("伺服器沒有 numpy，這次只更新 HPO 查表、未建表型知識庫（.venv/bin/pip install numpy）")
        else:
            if omim is None:
                try:
                    omim = B.platform_omim_rows()
                except Exception as exc:  # noqa: BLE001
                    omim = ([], {})
                    warnings.append(f"讀不到平台 OMIM（{type(exc).__name__}），知識庫的疾病–基因不含 OMIM")
            omim_rows, omim_info = omim
            if not omim_rows and "omim" in allow_missing:
                warnings.append("知識庫的疾病–基因不含 OMIM")
            if orphanet is None and "orphanet" not in allow_missing:
                warnings.append("Orphanet product6 不在本機，疾病–基因改用 HPO genes_to_disease 內的 Orphanet 對應")
            if mondo is None and "mondo" in allow_missing:
                warnings.append("沒有 MONDO，疾病未分組")
            extra = B.make_extra(result, omim_rows=omim_rows, omim_info=omim_info, orphanet=orphanet, mondo=mondo,
                                 allow_missing=set(allow_missing), progress=progress)
    sources = {"hp_obo": hp_obo, "hpoa": hpoa, "genes_to_disease": genes_to_disease}
    progress("寫入快照 …")
    snap = H.write_snapshot(result, Path(kb_root or H.KB_ROOT), sources, extra=extra)
    manifest = json.loads((snap / "manifest.json").read_text(encoding="utf-8"))
    return snap, warnings, manifest


def summary_text(snap: Path, manifest: dict, warnings: list[str]) -> str:
    """Plain-language summary printed at the end of the command and shown in the UI."""
    qc = manifest.get("qc") or {}
    kb = manifest.get("kb") or {}
    lines = [f"HPO 版本：{manifest.get('hp_version', '')}（註記 {manifest.get('hpoa_version', '')}）",
             f"快照：{snap}"]
    pairs = qc.get("pairs") or qc.get("term_gene_pairs")
    if pairs:
        lines.append(f"HPO → 基因配對：{int(pairs):,}")
    if kb:
        kqc = kb.get("qc") or {}
        dg = kb.get("disease_gene") or {}
        lines.append(f"表型知識庫：{kqc.get('N', 0):,} 個疾病、{kqc.get('terms', 0):,} 個 term；"
                     f"QC {'全部通過' if kqc.get('passed') else '未通過'}")
        lines.append(f"疾病–基因：採用 {dg.get('used', 0):,} 筆，有基因的疾病 {dg.get('diseases_with_gene', 0):,} 個")
        srcs = kb.get("sources") or {}
        lines.append("來源：" + "、".join(
            f"{k} {('（經 genes_to_disease）' if v.get('fallback') else '') + str(v.get('version') or v.get('source') or '')}"
            for k, v in srcs.items()))
    lines += [f"注意：{w}" for w in warnings]
    return "\n".join(lines)


def run(release: str = "latest", *, data_dir: Path = PHENO_DATA_DIR, kb_root: Path | None = None,
        activate: bool = True, download=download, progress=print) -> dict:
    """The whole update; returns a summary dict."""
    folder = fetch_release(release, data_dir, download=download, progress=progress)
    warnings: list[str] = []
    orphanet, warn = fetch_optional("Orphanet", ORPHANET_URLS, data_dir / "orphanet" / "en_product6.xml",
                                    download=download, progress=progress)
    if warn:
        warnings.append(warn)
    mondo, warn = fetch_optional("MONDO", MONDO_URLS, data_dir / "mondo" / "mondo.obo",
                                 download=download, progress=progress)
    if warn:
        warnings.append(warn)
    snap, more, manifest = build_snapshot(folder / "hp.obo", folder / "phenotype.hpoa",
                                          folder / "genes_to_disease.txt", kb_root=kb_root, orphanet=orphanet,
                                          mondo=mondo, allow_missing=("omim", "mondo"), progress=progress)
    warnings += more
    installed = None
    if activate:
        H.activate(snap, kb_root)
        installed = install_obo(folder / "hp.obo", data_dir)
        try:
            n = write_picker_data(data_dir / "hp.obo", snap / H.GENE_TABLE, data_dir / PICKER_DATA)
            progress(f"更新「輸入臨床表徵」頁的 HPO 清單（{n:,} 個 term）")
        except Exception as exc:  # noqa: BLE001
            warnings.append(f"「輸入臨床表徵」頁的 HPO 清單沒更新（{type(exc).__name__}: {exc}）")
    return {"snapshot": snap.name, "path": str(snap), "hp_version": manifest.get("hp_version"),
            "hpoa_version": manifest.get("hpoa_version"), "activated": activate,
            "hp_obo_installed": bool(installed), "warnings": warnings,
            "summary": summary_text(snap, manifest, warnings)}


# ── background job (button) ───────────────────────────────────────────────
def _write_state(payload: dict) -> None:
    STATE_DIR.mkdir(parents=True, exist_ok=True)
    tmp = STATE_PATH.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(payload, ensure_ascii=False, indent=1), encoding="utf-8")
    os.replace(tmp, STATE_PATH)


def _read_state() -> dict:
    try:
        value = json.loads(STATE_PATH.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {}
    return value if isinstance(value, dict) else {}


def _pid_alive(pid) -> bool:
    try:
        os.kill(int(pid), 0)
        return True
    except (TypeError, ValueError, OSError):
        return False


@contextmanager
def _lock(path: Path, *, blocking: bool):
    path.parent.mkdir(parents=True, exist_ok=True)
    handle = path.open("a+")
    try:
        fcntl.flock(handle.fileno(), fcntl.LOCK_EX | (0 if blocking else fcntl.LOCK_NB))
        yield
    finally:
        fcntl.flock(handle.fileno(), fcntl.LOCK_UN)
        handle.close()


def status() -> dict:
    state = _read_state()
    running = state.get("state") in {"queued", "running"} and _pid_alive(state.get("pid"))
    if state.get("state") in {"queued", "running"} and not running and state.get("pid"):
        state = {**state, "state": "failed", "finished_at": _now(),
                 "error": state.get("error") or "更新程序已中止"}
        _write_state(state)
    return {**state, "running": running}


def start_background(release: str = "latest", user: str = "") -> dict:
    with _lock(STATE_DIR / ".spawn.lock", blocking=True):
        current = status()
        if current.get("running"):
            return {**current, "already_running": True}
        run_id = uuid.uuid4().hex
        _write_state({"run_id": run_id, "state": "queued", "step": "準備中", "release": release, "user": user,
                      "created_at": _now(), "started_at": None, "finished_at": None, "pid": None,
                      "error": None, "log": []})
        env = os.environ.copy()
        backend = str(REPO_ROOT / "backend")
        env["PYTHONPATH"] = backend + (os.pathsep + env["PYTHONPATH"] if env.get("PYTHONPATH") else "")
        with LOG_PATH.open("a", buffering=1, encoding="utf-8") as log:
            proc = subprocess.Popen([sys.executable, "-m", "app.workers.hpo_update", "--run-id", run_id,
                                     "--release", release], cwd=str(REPO_ROOT), env=env, stdout=log,
                                    stderr=subprocess.STDOUT, start_new_session=True)
        state = _read_state()
        if state.get("run_id") == run_id and state.get("state") == "queued":
            state["pid"] = proc.pid
            _write_state(state)
        return status()


def run_job(run_id: str, release: str = "latest") -> int:
    """Foreground worker body for the button."""
    try:
        with _lock(STATE_DIR / ".update.lock", blocking=False):
            state = {**_read_state(), "run_id": run_id, "state": "running", "step": "開始",
                     "started_at": _now(), "pid": os.getpid(), "error": None, "log": []}
            _write_state(state)

            def progress(message: str) -> None:
                print(f"[hpo-update] {message}", flush=True)
                current = _read_state()
                log = (current.get("log") or [])[-40:] + [message]
                current.update(step=message, log=log, pid=os.getpid(), state="running")
                _write_state(current)

            result = run(release, progress=progress)
            _write_state({**_read_state(), "state": "done", "step": "完成", "finished_at": _now(),
                          "result": result})
            print(result["summary"], flush=True)
            return 0
    except BlockingIOError:
        _write_state({**_read_state(), "state": "failed", "finished_at": _now(), "error": "已有更新在進行"})
        return 1
    except BaseException as exc:  # noqa: BLE001
        traceback.print_exc()
        _write_state({**_read_state(), "state": "failed", "finished_at": _now(),
                      "error": f"{type(exc).__name__}: {exc}"})
        return 1
