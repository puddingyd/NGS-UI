"""Targeted somatic jobs and additive, germline-first SNV results.

Only a completed run is published. Raw germline annotation is never modified.
Commands and resource paths come from an administrator-owned JSON file, never
from request payloads. Workers use the tertiary sample lock as well.
"""
from __future__ import annotations

import csv
import fcntl
import json
import os
import re
import subprocess
import shutil
import sys
import time
import uuid
from contextlib import contextmanager
from pathlib import Path

from .. import config
from . import panel_deadzone, sample_layout

ACTIVE = {"queued", "running", "cancelling"}
SID_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]{0,127}$")
RUN_RE = re.compile(r"^[0-9a-f]{32}$")


def validate_sid(sid: str) -> str:
    if not SID_RE.fullmatch(sid):
        raise ValueError("Invalid sample ID")
    return sid


def job_root() -> Path:
    return config.DATA_ROOT / "jobs" / "somatic"


def job_dir(run_id: str) -> Path:
    if not RUN_RE.fullmatch(run_id):
        raise ValueError("Invalid run ID")
    return job_root() / run_id


def result_dir(sid: str, run_id: str) -> Path:
    validate_sid(sid)
    job_dir(run_id)
    return config.PIPELINE_OUT_ROOT / sid / "09_somatic" / run_id


def index_path(sid: str) -> Path:
    validate_sid(sid)
    return sample_layout.state_file(sid, "somatic.json")


def read_json(path: Path, default=None):
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except FileNotFoundError:
        return {} if default is None else default


def atomic_json(path: Path, payload) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.{uuid.uuid4().hex}.tmp")
    try:
        tmp.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
        os.replace(tmp, path)
    finally:
        tmp.unlink(missing_ok=True)


def signature(path: Path) -> list:
    stat = path.stat()
    return [str(path.resolve()), stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns]


def settings() -> dict:
    path = Path(os.environ.get("NGS_UI_SOMATIC_CONFIG", config.DATA_ROOT / "somatic_config.json"))
    if not path.is_file():
        raise ValueError("尚未設定 Somatic 執行環境：請設定 NGS_UI_SOMATIC_CONFIG")
    cfg = read_json(path)
    if cfg.get("assembly") != "GRCh38":
        raise ValueError("Somatic reference assembly 必須為 GRCh38")
    for key in ("reference", "gene_regions", "vep_cache", "germline_resource", "clinvar_vcf"):
        if not cfg.get(key) or not Path(cfg[key]).exists():
            raise ValueError(f"Somatic 資源缺失：{key}")
    if not cfg.get("clinvar_release") or not cfg.get("gene_regions_release") or not str(cfg.get("vep_cache_version", "")).isdigit():
        raise ValueError("必須設定 clinvar_release、gene_regions_release 與數字 vep_cache_version")
    for suffix in (".fai",):
        if not Path(cfg["reference"] + suffix).is_file():
            raise ValueError(f"Reference 缺少 {suffix}")
    if not Path(cfg["reference"]).with_suffix(".dict").is_file():
        raise ValueError("Reference 缺少 sequence dictionary (.dict)")
    for key in ("pon", "contamination_sites"):
        if cfg.get(key) and not Path(cfg[key]).is_file():
            raise ValueError(f"Somatic 資源缺失：{key}")
    for tool in ("gatk", "bcftools", "samtools", "vep"):
        cmd = cfg.setdefault(tool + "_command", [tool])
        if not isinstance(cmd, list) or not cmd or not all(isinstance(x, str) and x for x in cmd):
            raise ValueError(f"{tool}_command 必須是非空 JSON 字串陣列")
        if not shutil.which(cmd[0]):
            raise ValueError(f"找不到執行工具：{cmd[0]}")
    cfg["config_signature"] = signature(path)
    return cfg


def canonical_gene(value: str) -> str:
    return panel_deadzone.canonical_panel_gene_symbol(value.upper())


def contigs(cfg: dict) -> dict[str, int]:
    with open(cfg["reference"] + ".fai") as handle:
        return {p[0]: int(p[1]) for line in handle if (p := line.split())}


def chromosome(value: str, lengths: dict) -> str:
    bare = re.sub(r"^chr", "", value, flags=re.I).upper()
    if bare not in {str(i) for i in range(1, 23)} | {"X", "Y"}:
        raise ValueError(f"不支援的核基因組染色體：{value}")
    chrom = "chr" + bare
    if chrom not in lengths:
        raise ValueError("Reference 必須使用 chr1–chr22/chrX/chrY 命名")
    return chrom


def merge_intervals(intervals: list[list]) -> list[list]:
    out = []
    for chrom, start, end in sorted(intervals):
        if out and out[-1][0] == chrom and start <= out[-1][2] + 1:
            out[-1][2] = max(end, out[-1][2])
        else:
            out.append([chrom, start, end])
    return out


def resolve_targets(payload: dict, cfg: dict) -> dict:
    """Public coordinates are 1-based inclusive; region reference is BED."""
    genes_text = payload.get("genes", "")
    positions_text = payload.get("positions", "")
    if not isinstance(genes_text, str) or not isinstance(positions_text, str):
        raise ValueError("基因及座標必須為文字")
    if len(genes_text) + len(positions_text) > 20000:
        raise ValueError("輸入內容過長")
    genes = sorted({canonical_gene(g) for g in re.split(r"[\s,;，；]+", genes_text.strip()) if g})
    if len(genes) > 200:
        raise ValueError("每次最多 200 個基因")
    mode = payload.get("region_mode", "exons")
    if mode not in {"exons", "gene"}:
        raise ValueError("分析範圍須為 exons 或 gene")
    lengths = contigs(cfg)
    intervals, found = [], set()
    padding = 20 if mode == "exons" else 0
    if genes:
        with open(cfg["gene_regions"], encoding="utf-8") as handle:
            for line in handle:
                if not line.strip() or line.startswith("#"):
                    continue
                p = line.rstrip().split("\t")
                if len(p) != 5:
                    raise ValueError("gene_regions 須為五欄 BED：chrom/start/end/gene/exon或gene")
                gene = canonical_gene(p[3])
                if gene not in genes or p[4] != ("exon" if mode == "exons" else "gene"):
                    continue
                chrom = chromosome(p[0], lengths)
                start, end = int(p[1]) + 1, int(p[2])
                if not 1 <= start <= end <= lengths[chrom]:
                    raise ValueError(f"基因 reference 座標不正確：{gene}")
                intervals.append([chrom, max(1, start - padding), min(lengths[chrom], end + padding)])
                found.add(gene)
        if set(genes) - found:
            raise ValueError("找不到基因範圍：" + ", ".join(sorted(set(genes) - found)))
    positions = []
    for token in re.split(r"[\s,;，；]+", positions_text.strip()):
        if not token:
            continue
        match = re.fullmatch(r"((?:chr)?(?:[0-9]+|X|Y)):(\d+)(?:-(\d+))?(?::([ACGT]+)>([ACGT]+))?", token, re.I)
        if not match:
            raise ValueError(f"座標格式錯誤：{token}；例如 chr1:100 或 chr1:100-200 或 chr1:100:A>G")
        chrom = chromosome(match[1], lengths)
        start, end = int(match[2]), int(match[3] or match[2])
        ref, alt = (match[4] or "").upper(), (match[5] or "").upper()
        if ref and (match[3] or ref == alt):
            raise ValueError("指定 allele 須為單一位置且 REF 不等於 ALT")
        if ref:
            end = start + len(ref) - 1
        if not 1 <= start <= end <= lengths[chrom]:
            raise ValueError(f"座標超出 reference 範圍：{token}")
        positions.append({"chrom": chrom, "start": start, "end": end, "ref": ref, "alt": alt})
        intervals.append([chrom, start, end])
    merged = merge_intervals(intervals)
    total = sum(end - start + 1 for _, start, end in merged)
    if not merged:
        raise ValueError("請輸入至少一個基因或座標")
    if len(positions) > 200 or total > int(cfg.get("max_target_bases", 20000000)):
        raise ValueError("指定範圍過大；請分批執行")
    return {"genes": genes, "positions": positions, "intervals": merged, "total_bases": total,
            "region_mode": mode, "exon_padding": padding, "assembly": "GRCh38",
            "gene_regions_release": cfg["gene_regions_release"]}


def sample_bams(sid: str) -> list[dict]:
    from ..routers import igv
    primary, source_sid, _ = igv._resolve_primary_bam(sid)
    hits = igv._bam_hits(source_sid, source=igv._pipeline_type_from_sidecar(sid))
    if primary and not hits:
        hits = [primary]
    return [hit for hit in hits if Path(hit["path"]).is_file() and igv._bam_index_for(Path(hit["path"]))]


def read_job(run_id: str) -> dict:
    job = read_json(job_dir(run_id) / "state.json")
    if job.get("status") in ACTIVE and time.time() - job.get("updated", 0) > 30:
        pid = job.get("pid")
        try:
            if not pid:
                raise ProcessLookupError()
            os.kill(pid, 0)
        except ProcessLookupError:
            job.update(status="failed", error="Somatic worker 已停止；既有結果未受影響")
            atomic_json(job_dir(run_id) / "state.json", job)
    return job


def jobs(sid: str | None = None) -> list[dict]:
    out = []
    if job_root().exists():
        for path in job_root().glob("*/state.json"):
            job = read_job(path.parent.name)
            if sid is None or job.get("sample_id") == sid:
                out.append(job)
    return sorted(out, key=lambda j: j.get("created", 0), reverse=True)


def active_ids() -> set[str]:
    return {j["sample_id"] for j in jobs() if j.get("status") in ACTIVE}


@contextmanager
def submission_lock():
    job_root().mkdir(parents=True, exist_ok=True)
    with (job_root() / ".submit.lock").open("a") as handle:
        fcntl.flock(handle, fcntl.LOCK_EX)
        yield


def start(sid: str, payload: dict, username: str = "") -> dict:
    from . import dragen_jobs
    validate_sid(sid)
    meta = read_json(sample_layout.state_file(sid, "sample_metadata.json"))
    if not meta:
        raise ValueError("請先載入個案")
    if (meta.get("genome_build") or "hg38").lower() not in {"hg38", "grch38"}:
        raise ValueError("Somatic 分析僅支援 GRCh38 個案")
    cfg = settings()
    targets = resolve_targets(payload, cfg)
    bams = sample_bams(sid)
    bam = payload.get("bam_path")
    if not isinstance(bam, str) or bam not in {x["path"] for x in bams}:
        raise ValueError("請選擇目前個案的有效 BAM 與 index")
    raw = sample_layout.snv_raw_tsv(sid)
    if not raw.is_file():
        raise ValueError("缺少完整 germline TSV，無法安全排除既有點位")
    replace = payload.get("replace_run_id")
    if replace is not None and not isinstance(replace, str):
        raise ValueError("Invalid replacement run ID")
    if replace:
        prior = read_job(replace)
        if prior.get("sample_id") != sid or prior.get("status") != "completed":
            raise ValueError("只能重新執行此個案已完成的分析")
    with submission_lock():
        if sid in dragen_jobs.active_sample_ids():
            raise ValueError("此個案已有分析工作執行中")
        if len(active_ids()) >= int(cfg.get("max_concurrent_jobs", 1)):
            raise ValueError("Somatic 執行名額已滿，請待目前工作完成後再試")
        run_id = uuid.uuid4().hex
        directory = job_dir(run_id)
        directory.mkdir()
        record = {"run_id": run_id, "sample_id": sid, "status": "queued", "step": "queued",
                  "created": time.time(), "updated": time.time(), "requested_by": username,
                  "targets": targets, "bam_path": bam, "bam_signature": signature(Path(bam)),
                  "raw_signature": signature(raw), "replace_run_id": replace,
                  "request": {k: payload.get(k, "") for k in ("genes", "positions", "region_mode")}}
        atomic_json(directory / "state.json", record)
        atomic_json(directory / "config.json", cfg)
        env = dict(os.environ, PYTHONPATH=str(config.REPO_ROOT / "backend"))
        try:
            with (directory / "log.txt").open("a") as log:
                subprocess.Popen([sys.executable, "-m", "app.workers.somatic_run", run_id],
                                 cwd=config.REPO_ROOT, env=env, stdout=log, stderr=subprocess.STDOUT,
                                 start_new_session=True)
        except Exception as exc:
            record.update(status="failed", error=str(exc))
            atomic_json(directory / "state.json", record)
            raise
    return record


def manifest(sid: str) -> dict:
    return read_json(index_path(sid), {"runs": [], "selected_filtered": {}})


def summary(sid: str) -> dict:
    data = manifest(sid)
    raw = sample_layout.snv_raw_tsv(sid)
    sig = signature(raw) if raw.exists() else []
    stale = any(r.get("raw_signature") != sig for r in data.get("runs", []) if not r.get("archived"))
    return {"completed": bool(data.get("runs")), "stale": stale,
            "run_count": len(data.get("runs", []))}


def missing_review_ids(sid: str, marked: set[str], available: set[str]) -> list[str]:
    """Do not silently lose a selected historical result after germline changes."""
    missing = marked - available
    if not missing:
        return []
    found = set()
    for run in manifest(sid).get("runs", []):
        path = result_dir(sid, run["run_id"]) / "candidates.json"
        for candidate in read_json(path, []):
            if candidate.get("id") in missing:
                found.add(candidate["id"])
    return sorted(found)


def load_variants(sid: str, *, wanted: set[str] | None = None, genes: set[str] | None = None,
                  keep_ids: set[str] | None = None) -> dict:
    """Load only additive results. Changed germline requires a new dedup run.

    No germline card is ever modified; exclusion was against the full raw TSV,
    before annotation, not against a filtered review or a visible card list.
    """
    from ..adapters.snv_tsv import _row_to_variant, merge_snv_variant_row
    data = manifest(sid)
    if not data.get("runs"):
        return {}
    raw = sample_layout.snv_raw_tsv(sid)
    sig = signature(raw) if raw.exists() else []
    out = {}
    keep_ids = (keep_ids or set()) | (wanted or set())
    for run in sorted(data["runs"], key=lambda r: r["created"]):
        if run.get("raw_signature") != sig:
            continue
        path = result_dir(sid, run["run_id"]) / "annotations.tsv"
        if not path.is_file():
            raise RuntimeError("已發布 Somatic annotation 遺失，請檢查伺服器檔案")
        current = {}
        with path.open(encoding="utf-8") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                vid = f'{row["CHROM"]}-{row["POS"]}-{row["REF"]}-{row["ALT"]}'
                if wanted is not None and vid not in wanted:
                    continue
                if run.get("archived") and vid not in keep_ids:
                    continue
                if genes is not None and canonical_gene(row.get("GENE", "")) not in genes:
                    continue
                filt = row.get("SOMATIC_FILTER", "")
                selected = vid in data.get("selected_filtered", {}).get(run["run_id"], [])
                if filt != "PASS" and not selected and vid not in keep_ids:
                    if not run.get("archived"):
                        out.pop(vid, None)
                    continue
                v = _row_to_variant(row)
                v.update(somatic=True, somatic_run_id=run["run_id"], somatic_filter=filt,
                         somatic_clinvar_release=run.get("clinvar_release", ""),
                         somatic_bam=run.get("bam_path", ""),
                         somatic_qc={key: row.get("SOMATIC_" + key, "") for key in
                                     ("TLOD", "MMQ", "MBQ", "MPOS", "F1R2", "F2R1")},
                         somatic_historical=bool(run.get("archived")),
                         somatic_validation="未驗證", low_depth=(v.get("depth") or 0) < 10)
                merge_snv_variant_row(current, v)
        # Latest observation wins without combining AD/DP from different runs.
        out.update(current)
    return out


def select_filtered(sid: str, run_id: str, vid: str) -> None:
    data = manifest(sid)
    if not any(r["run_id"] == run_id and not r.get("archived") for r in data.get("runs", [])):
        raise ValueError("找不到已完成的分析")
    rows = read_json(result_dir(sid, run_id) / "candidates.json", [])
    if not any(r["id"] == vid for r in rows):
        raise ValueError("找不到候選點位")
    selected = data.setdefault("selected_filtered", {}).setdefault(run_id, [])
    if vid not in selected:
        selected.append(vid)
    atomic_json(index_path(sid), data)
