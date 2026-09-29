"""Targeted somatic jobs and additive, germline-first SNV results.

Only a completed run is published. Raw germline annotation is never modified.
Commands and resource paths come from an administrator-owned JSON file, never
from request payloads. Workers use the tertiary sample lock as well.
"""
from __future__ import annotations

import csv
import fcntl
import gzip
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
SOMATIC_DBNSFP_FIELDS = (
    "PKNN_LLR", "AlphaMissense_score", "BayesDel_noAF_score", "ESM1b_score",
    "VARITY_R_score", "SIFT_score", "SIFT_pred", "DANN_score",
    "PHACTboost_score", "phyloP100way_vertebrate", "GERP++_RS", "REVEL_score",
    "MutPred2_score", "MutPred2_pred", "VEST4_score", "CADD_phred",
)
LOG_STEP_LABELS = {
    "preflight": "檢查輸入與工具",
    "mutect2": "Mutect2 變異偵測",
    "filtering:orientation": "校正方向性偏差",
    "filtering:pileup": "估計樣本污染",
    "filtering:contamination": "計算污染比例",
    "filtering:mutect-calls": "套用品質過濾",
    "subtract-germline": "排除 germline 已有點位",
    "annotation:vep": "VEP、dbNSFP 與 SpliceAI",
    "annotation:clinvar": "固定版 ClinVar",
    "annotation:clinvar-latest": "最新版 ClinVar 比對",
    "annotation:genebe": "GeneBe ACMG",
    "annotation:giab": "GIAB 困難區域",
    "annotation:inhouse-af": "本院族群頻率",
    "annotation:mane": "MANE RefSeq",
    "annotation:litvar2": "LitVar2 文獻",
    "annotation:gpn-msa": "GPN-MSA",
    "coverage": "檢查指定範圍覆蓋",
    "publishing": "發布結果",
    "completed": "完成",
    "failed": "失敗",
    "cancelled": "已取消",
}


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


def validate_dbnsfp_header(path: Path) -> None:
    """Fail before a run if the configured tertiary dbNSFP lacks UI fields."""
    try:
        with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
            header = next((line.rstrip("\n").lstrip("#").split("\t")
                           for line in handle if line.startswith("#chr\t")), None)
    except OSError as exc:
        raise ValueError(f"Somatic dbNSFP 無法讀取：{path}") from exc
    if not header:
        raise ValueError(f"Somatic dbNSFP 缺少 #chr header：{path}")
    missing = [field for field in SOMATIC_DBNSFP_FIELDS if field not in header]
    if missing:
        raise ValueError(f"Somatic dbNSFP 缺少欄位：{', '.join(missing)}")


def _command_label(args: list[str]) -> str:
    names = [Path(token).name for token in args]
    for action, label in (
        ("Mutect2", "Mutect2 calling"),
        ("LearnReadOrientationModel", "建立方向性偏差模型"),
        ("GetPileupSummaries", "整理污染估計位點"),
        ("CalculateContamination", "計算污染比例"),
        ("FilterMutectCalls", "Mutect2 品質過濾"),
    ):
        if action in args:
            return label
    for tool in ("bcftools", "samtools", "vep"):
        if tool in names:
            index = names.index(tool)
            action = args[index + 1] if index + 1 < len(args) else ""
            return f"{tool} {action}".strip()
    for token in names:
        if token.endswith(".py"):
            return token.removesuffix(".py").replace("_", " ")
    return names[0] if names else "command"


def format_log(raw: str) -> str:
    """Turn verbose third-party output into a reviewer-readable timeline."""
    output: list[str] = []
    for raw_line in raw.splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line in LOG_STEP_LABELS:
            output.extend(([""] if output else []) + [f"【{LOG_STEP_LABELS[line]}】"])
            continue
        if line.startswith("RUN "):
            try:
                args = json.loads(line[4:])
            except json.JSONDecodeError:
                args = []
            output.append(f"  執行：{_command_label([str(value) for value in args])}")
            continue
        lower = line.lower()
        keep = (
            "warning" in lower or "error" in lower or "traceback" in lower
            or "total reads filtered" in lower or line.startswith("Lines   total/")
            or (line.startswith("[clinvar]") and any(word in lower for word in ("matched", "scanned", "backfilled", "done")))
            or (line.startswith("[genebe]") and not any(word in lower for word in (" db:", "sqlite ready")))
            or line.startswith("[gpn-msa]") or line.startswith("[giab-strata]")
            or line.startswith("[inhouse-af]") or line.startswith("[mane-refseq]")
            or line.startswith("[litvar2]")
        )
        if keep:
            output.append("  " + line)
    if not output:
        return "目前沒有可顯示的執行摘要。"
    return "\n".join(output).strip()


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
    # Somatic always exposes the complete tertiary in-silico panel. The
    # academic dbNSFP file is the tertiary pipeline's 5.3a build with P-KNN
    # merged in, so it contains both the original and Research-only fields.
    cfg.setdefault(
        "dbnsfp_academic",
        str(Path(cfg["reference"]).parent / "tertiary/dbnsfp/dbNSFP5.3a_with_pknn_grch38.gz"),
    )
    cfg.setdefault("dbnsfp_academic_version", "5.3a")
    cfg.setdefault("spliceai_snv", str(config.BIOTOOLS_DIR / "spliceai/spliceai_scores.raw.snv.hg38.vcf.gz"))
    cfg.setdefault("spliceai_indel", str(config.BIOTOOLS_DIR / "spliceai/spliceai_scores.raw.indel.hg38.vcf.gz"))
    for key in ("dbnsfp_academic", "spliceai_snv", "spliceai_indel"):
        resource = Path(cfg[key])
        for path in (resource, Path(str(resource) + ".tbi")):
            if not path.is_file():
                raise ValueError(f"Somatic complete predictor 資源缺失：{path}")
    validate_dbnsfp_header(Path(cfg["dbnsfp_academic"]))
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
    usable = []
    broken = []
    for run in data.get("runs", []):
        annotation = result_dir(sid, run["run_id"]) / "annotations.tsv"
        (usable if annotation.is_file() else broken).append(run)
    stale = any(r.get("raw_signature") != sig for r in usable if not r.get("archived"))
    return {"completed": bool(usable), "stale": stale,
            "run_count": len(data.get("runs", [])), "published_run_count": len(usable),
            "broken_run_ids": [r["run_id"] for r in broken]}


def delete_run(sid: str, run_id: str) -> dict:
    """Delete one terminal job and its published result, including broken publications."""
    validate_sid(sid)
    directory = job_dir(run_id)
    job = read_job(run_id)
    data = manifest(sid)
    indexed = [run for run in data.get("runs", []) if run.get("run_id") == run_id]
    if job and job.get("sample_id") != sid:
        raise FileNotFoundError("找不到分析工作")
    if not job and not indexed:
        raise FileNotFoundError("找不到分析工作")
    if job.get("status") in ACTIVE:
        raise RuntimeError("分析仍在執行，請先終止後再刪除")

    if indexed:
        data["runs"] = [run for run in data.get("runs", []) if run.get("run_id") != run_id]
        selected = data.get("selected_filtered", {})
        if isinstance(selected, dict):
            selected.pop(run_id, None)
        atomic_json(index_path(sid), data)

    published = result_dir(sid, run_id)
    result_existed = published.exists()
    job_existed = directory.exists()
    shutil.rmtree(published, ignore_errors=True)
    shutil.rmtree(directory, ignore_errors=True)
    return {"deleted": True, "result_deleted": result_existed, "job_deleted": job_existed}


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
            # A manually removed or partially lost publication must not make the
            # entire case unloadable.  summary() exposes it for cleanup.
            continue
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
                         somatic_validation="未確認", low_depth=(v.get("depth") or 0) < 10)
                merge_snv_variant_row(current, v)
        # Latest observation wins without combining AD/DP from different runs.
        out.update(current)
    return out


def candidate_variants(sid: str, run_id: str) -> list[dict]:
    """Return every candidate in one run as a full SNV card payload."""
    from ..adapters.snv_tsv import _row_to_variant, merge_snv_variant_row

    data = manifest(sid)
    run = next((item for item in data.get("runs", []) if item.get("run_id") == run_id), None)
    if run is None:
        return []
    annotation = result_dir(sid, run_id) / "annotations.tsv"
    if not annotation.is_file():
        return []
    current: dict[str, dict] = {}
    with annotation.open(encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            vid = f'{row["CHROM"]}-{row["POS"]}-{row["REF"]}-{row["ALT"]}'
            variant = _row_to_variant(row)
            variant.update(
                somatic=True,
                somatic_run_id=run_id,
                somatic_filter=row.get("SOMATIC_FILTER", ""),
                somatic_clinvar_release=run.get("clinvar_release", ""),
                somatic_bam=run.get("bam_path", ""),
                somatic_qc={key: row.get("SOMATIC_" + key, "") for key in
                            ("TLOD", "MMQ", "MBQ", "MPOS", "F1R2", "F2R1")},
                somatic_validation="未確認",
                low_depth=(variant.get("depth") or 0) < 10,
            )
            merge_snv_variant_row(current, variant)
    selected = set(data.get("selected_filtered", {}).get(run_id, []))
    ordered = read_json(result_dir(sid, run_id) / "candidates.json", [])
    result = []
    for item in ordered:
        vid = item.get("id", "")
        variant = current.get(vid)
        if not variant:
            continue
        result.append({
            "id": vid,
            "variant": variant,
            "included": item.get("filter") == "PASS" or vid in selected,
        })
    return result


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
