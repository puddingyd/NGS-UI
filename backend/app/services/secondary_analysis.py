"""Secondary-analysis FASTQ indexing and DGX-2 samplesheet helpers."""
from __future__ import annotations

import csv
import fcntl
import json
import re
import shlex
import subprocess
import time
from datetime import datetime, timedelta, timezone
from pathlib import Path

from ..config import (
    SECONDARY_DGX_ENV_SCRIPT,
    SECONDARY_DGX_LAUNCH_ROOT,
    SECONDARY_DGX_OUTPUT_ROOT,
    SECONDARY_DGX_RAW_ROOT,
    SECONDARY_DGX_SAMPLESHEET_STAGING_ROOT,
    SECONDARY_DGX_WORK_ROOT,
    SECONDARY_FASTQ_INDEX_PATH,
    SECONDARY_FASTQ_INDEX_TTL_HOURS,
    SECONDARY_OUTPUT_ROOT,
    SECONDARY_SAMPLESHEET_STAGING_ROOT,
    SECONDARY_WES_FASTQ_ROOTS,
    SECONDARY_WGS_FASTQ_ROOTS,
)

TAIPEI_TZ = timezone(timedelta(hours=8))
_SAMPLE_ID_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]{0,127}$")
_NEXTSEQ_RE = re.compile(r"^(.+)_S\d+_R([12])_001\.fastq\.gz$", re.I)
_WGS_LANE_RE = re.compile(r"^(.+)_S\d+_L(\d{3})_R([12])_001\.fastq\.gz$", re.I)
_REANALYSIS_RE = re.compile(r"^(.+)\.R([12])\.clean\.fastq\.gz$", re.I)
_FASTQ_INDEX_SCHEMA_VERSION = 2


def _now() -> str:
    return datetime.now(TAIPEI_TZ).strftime("%Y-%m-%d %H:%M:%S")


def _validate_sample_id(sample_id: str) -> str:
    sid = (sample_id or "").strip()
    if not _SAMPLE_ID_RE.match(sid):
        raise ValueError("sample_id must match [A-Za-z0-9][A-Za-z0-9._-]{0,127}")
    return sid


def _find_fastqs(root: Path) -> list[Path]:
    if not root.is_dir():
        return []
    try:
        proc = subprocess.run(
            ["find", str(root), "-type", "f", "-name", "*.fastq.gz"],
            check=False,
            capture_output=True,
            text=True,
        )
        if proc.returncode == 0:
            return [Path(line) for line in proc.stdout.splitlines() if line]
    except OSError:
        pass
    return list(root.glob("**/*.fastq.gz"))


def _first_existing_parent(path: Path, roots: list[Path]) -> Path | None:
    try:
        rp = path.resolve()
    except OSError:
        return None
    for root in roots:
        if not root.exists():
            continue
        try:
            rr = root.resolve()
        except OSError:
            continue
        try:
            rp.relative_to(rr)
            return root
        except ValueError:
            continue
    return None


def _run_name_for_nextseq(path: Path) -> str:
    parts = path.parts
    if "Analysis" in parts:
        idx = parts.index("Analysis")
        if idx > 0:
            return parts[idx - 1]
    return path.parent.name


def _batch_date_from_run(run: str) -> str:
    m8 = re.search(r"(?<!\d)(20\d{6})(?!\d)", run or "")
    if m8:
        return m8.group(1)[2:]
    m6 = re.search(r"(?<!\d)(\d{6})(?!\d)", run or "")
    if m6:
        return m6.group(1)
    return datetime.now(TAIPEI_TZ).strftime("%y%m%d")


def _run_sort_date(run: str) -> int:
    try:
        return int(_batch_date_from_run(run))
    except ValueError:
        return 0


def _entry(
    *,
    seq_type: str,
    sample_id: str,
    fastq_1: Path,
    fastq_2: Path,
    run: str,
    input_dir: Path,
    lane: str = "",
    reanalysis: bool = False,
    source_sample_id: str = "",
) -> dict:
    try:
        st = fastq_1.stat()
    except OSError:
        st = None
    return {
        "seq_type": seq_type,
        "sample_id": sample_id,
        "source_sample_id": source_sample_id or sample_id,
        "fastq_1": str(fastq_1),
        "fastq_2": str(fastq_2),
        "run": run,
        "input_dir": str(input_dir),
        "lane": lane,
        "reanalysis": reanalysis,
        "size": (st.st_size if st else 0),
        "mtime": (st.st_mtime if st else 0),
    }


def list_wes_fastqs() -> list[dict]:
    grouped: dict[tuple[str, str, str], dict[str, Path]] = {}
    reanalysis_grouped: dict[Path, dict[str, Path]] = {}
    for root in SECONDARY_WES_FASTQ_ROOTS:
        for p in _find_fastqs(root):
            name = p.name
            m_re = _REANALYSIS_RE.match(name)
            if "Reanalysis" in str(root) or m_re:
                if m_re:
                    reanalysis_grouped.setdefault(p.parent, {})[m_re.group(2)] = p
                continue
            m = _NEXTSEQ_RE.match(name)
            if not m:
                continue
            sample, read = m.group(1), m.group(2)
            run = _run_name_for_nextseq(p)
            key = (sample, run, str(p.parent))
            grouped.setdefault(key, {})[read] = p

    out: list[dict] = []
    for (sample, run, input_dir), reads in grouped.items():
        if "1" in reads and "2" in reads:
            out.append(_entry(
                seq_type="WES",
                sample_id=sample,
                fastq_1=reads["1"],
                fastq_2=reads["2"],
                run=run,
                input_dir=Path(input_dir),
            ))
    for folder, reads in reanalysis_grouped.items():
        if "1" not in reads or "2" not in reads:
            continue
        sample = folder.name
        source = _REANALYSIS_RE.match(reads["1"].name).group(1)
        out.append(_entry(
            seq_type="WES",
            sample_id=sample,
            source_sample_id=source,
            fastq_1=reads["1"],
            fastq_2=reads["2"],
            run=folder.parent.name,
            input_dir=folder,
            reanalysis=True,
        ))
    out.sort(key=lambda r: (-_run_sort_date(r.get("run", "")), r["sample_id"]))
    return out


def list_wgs_fastqs() -> list[dict]:
    lanes: dict[tuple[str, str, str, str], dict[str, Path]] = {}
    for root in SECONDARY_WGS_FASTQ_ROOTS:
        for p in _find_fastqs(root):
            if p.name.endswith(".md5"):
                continue
            run = p.parent.parent.name if p.parent.name == "fastq.gz" else p.parent.name
            m = _WGS_LANE_RE.match(p.name)
            if m:
                sample, lane, read = m.group(1), f"L{m.group(2)}", m.group(3)
                lanes.setdefault((sample, run, str(p.parent), lane), {})[read] = p

    grouped_samples: dict[tuple[str, str, str], list[dict]] = {}
    for (sample, run, input_dir, lane), reads in lanes.items():
        if "1" in reads and "2" in reads:
            grouped_samples.setdefault((sample, run, input_dir), []).append({
                "lane": lane,
                "fastq_1": str(reads["1"]),
                "fastq_2": str(reads["2"]),
            })

    out: list[dict] = []
    for (sample, run, input_dir), sample_lanes in grouped_samples.items():
        sample_lanes.sort(key=lambda item: item["lane"])
        first = sample_lanes[0]
        total_size = 0
        latest_mtime = 0.0
        for lane_payload in sample_lanes:
            for path_key in ("fastq_1", "fastq_2"):
                try:
                    stat = Path(lane_payload[path_key]).stat()
                except OSError:
                    continue
                total_size += stat.st_size
                latest_mtime = max(latest_mtime, stat.st_mtime)
        row = _entry(
            seq_type="WGS",
            sample_id=sample,
            fastq_1=Path(first["fastq_1"]),
            fastq_2=Path(first["fastq_2"]),
            run=run,
            input_dir=Path(input_dir),
        )
        row.update({
            "lanes": sample_lanes,
            "lane_count": len(sample_lanes),
            "fastq_file_count": len(sample_lanes) * 2,
            "size": total_size,
            "mtime": latest_mtime,
        })
        out.append(row)
    out.sort(key=lambda r: (-_run_sort_date(r.get("run", "")), r["sample_id"]))
    return out


def load_index() -> dict | None:
    p = SECONDARY_FASTQ_INDEX_PATH
    if not p.is_file():
        return None
    try:
        idx = json.loads(p.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, OSError):
        return None
    if idx.get("schema_version") != _FASTQ_INDEX_SCHEMA_VERSION:
        return None
    return idx


def save_index(idx: dict) -> None:
    p = SECONDARY_FASTQ_INDEX_PATH
    p.parent.mkdir(parents=True, exist_ok=True)
    tmp = p.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(idx, ensure_ascii=False, indent=2), encoding="utf-8")
    tmp.replace(p)


def refresh_index() -> dict:
    """Rescan WES/WGS FASTQs and atomically persist the shared index.

    The daily systemd timer and the manual UI refresh can overlap. A
    cross-process lock prevents both callers from writing the same
    temporary file concurrently.
    """
    lock_path = SECONDARY_FASTQ_INDEX_PATH.with_suffix(".refresh.lock")
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    with lock_path.open("a", encoding="utf-8") as lock_file:
        fcntl.flock(lock_file.fileno(), fcntl.LOCK_EX)
        t0 = time.time()
        wes = list_wes_fastqs()
        wgs = list_wgs_fastqs()
        idx = {
            "schema_version": _FASTQ_INDEX_SCHEMA_VERSION,
            "updated_at": _now(),
            "scan_duration_sec": round(time.time() - t0, 2),
            "wes": wes,
            "wgs": wgs,
        }
        save_index(idx)
        return idx


def index_is_stale(idx: dict | None) -> bool:
    if not idx or not idx.get("updated_at"):
        return True
    try:
        ts = datetime.fromisoformat(str(idx["updated_at"]).replace("Z", "+00:00"))
    except ValueError:
        return True
    if ts.tzinfo is None:
        ts = ts.replace(tzinfo=TAIPEI_TZ)
    age_h = (datetime.now(timezone.utc) - ts).total_seconds() / 3600
    return age_h >= SECONDARY_FASTQ_INDEX_TTL_HOURS


def _server_path_to_dgx(path: str) -> str:
    raw = str(path)
    replacements = [
        (str(SECONDARY_SAMPLESHEET_STAGING_ROOT), str(SECONDARY_DGX_SAMPLESHEET_STAGING_ROOT)),
        (str(SECONDARY_OUTPUT_ROOT), str(SECONDARY_DGX_OUTPUT_ROOT)),
        ("/home/datalake_Raw", str(SECONDARY_DGX_RAW_ROOT)),
        ("/home/datalake_Intermediate", "/datalake_Intermediate"),
    ]
    for src, dst in replacements:
        if raw == src or raw.startswith(src + "/"):
            return dst + raw[len(src):]
    if raw == "/home" or raw.startswith("/home/"):
        return raw[5:] or "/"
    return raw


def _unique_batch_name(base: str) -> str:
    name = base
    i = 2
    while (
        (SECONDARY_OUTPUT_ROOT / name).exists()
        or (SECONDARY_SAMPLESHEET_STAGING_ROOT / name).exists()
    ):
        name = f"{base}_{i}"
        i += 1
    return name


def suggest_batch_name(seq_type: str, samples: list[dict]) -> str:
    seq = (seq_type or samples[0].get("seq_type") or "WES").upper()
    runs = [str(s.get("run") or "") for s in samples if s.get("run")]
    date = _batch_date_from_run(runs[0] if runs else "")
    return _unique_batch_name(f"{date}_{seq}")


def _validate_fastq_pair(f1: Path, f2: Path) -> None:
    roots = SECONDARY_WES_FASTQ_ROOTS + SECONDARY_WGS_FASTQ_ROOTS
    for p in (f1, f2):
        if not p.is_file():
            raise FileNotFoundError(f"FASTQ not found: {p}")
        if _first_existing_parent(p, roots) is None:
            raise ValueError(f"FASTQ is outside configured roots: {p}")


def _normalize_wgs_lanes(payload: dict) -> list[dict]:
    raw_lanes = payload.get("lanes")
    if not isinstance(raw_lanes, list) or not raw_lanes:
        raw_lanes = [payload]
    normalized: list[dict] = []
    seen: set[str] = set()
    source_sample_id = ""
    for raw in raw_lanes:
        lane = str(raw.get("lane") or "").strip().upper()
        f1 = Path(str(raw.get("fastq_1") or ""))
        f2 = Path(str(raw.get("fastq_2") or ""))
        m1 = _WGS_LANE_RE.match(f1.name)
        m2 = _WGS_LANE_RE.match(f2.name)
        if (
            not re.fullmatch(r"L\d{3}", lane)
            or not m1
            or not m2
            or m1.group(1) != m2.group(1)
            or m1.group(2) != m2.group(2)
            or m1.group(3) != "1"
            or m2.group(3) != "2"
            or lane != f"L{m1.group(2)}"
        ):
            raise ValueError("WGS samplesheet must use matching L00x R1/R2 lane FASTQs, not merged FASTQs")
        if lane in seen:
            raise ValueError(f"duplicate WGS lane: {lane}")
        if source_sample_id and m1.group(1) != source_sample_id:
            raise ValueError("all WGS lanes in one selection must belong to the same source sample")
        _validate_fastq_pair(f1, f2)
        source_sample_id = m1.group(1)
        seen.add(lane)
        normalized.append({"lane": lane, "fastq_1": str(f1), "fastq_2": str(f2)})
    normalized.sort(key=lambda item: item["lane"])
    return normalized


def _normalize_sample(payload: dict, seq_type: str) -> dict:
    sample_id = _validate_sample_id(payload.get("sample_id") or "")
    if seq_type == "WGS":
        wgs_lanes = _normalize_wgs_lanes(payload)
        f1 = Path(wgs_lanes[0]["fastq_1"])
        f2 = Path(wgs_lanes[0]["fastq_2"])
        lane = ""
    else:
        f1 = Path(str(payload.get("fastq_1") or ""))
        f2 = Path(str(payload.get("fastq_2") or ""))
        _validate_fastq_pair(f1, f2)
        lane = (payload.get("lane") or "").strip()
        wgs_lanes = []
    return {
        "sample_id": sample_id,
        "source_sample_id": payload.get("source_sample_id") or sample_id,
        "fastq_1": str(f1),
        "fastq_2": str(f2),
        "sex": (payload.get("sex") or "unknown").strip().lower() or "unknown",
        "lane": lane,
        "seq_type": seq_type,
        "run": payload.get("run") or "",
        "input_dir": payload.get("input_dir") or str(f1.parent),
        "reanalysis": bool(payload.get("reanalysis")),
        "lanes": wgs_lanes,
    }


def _launch_command(batch_name: str, seq_type: str) -> str:
    session = f"ngs2_{batch_name}"
    # DGX2 launches now always use the single-runner profile, including
    # samplesheets containing more than one primary sample.
    profile = "dgx_single"
    run_gcnv = " \\\n    --run_gcnv true" if seq_type == "WES" else ""
    extended_analysis = (
        " \\\n    --run_manta"
        " \\\n    --run_expansionhunter"
        " \\\n    --run_automap"
    )
    script_path = f"/tmp/{session}.sh"
    staged_sheet = SECONDARY_DGX_SAMPLESHEET_STAGING_ROOT / batch_name / "samplesheet.csv"
    return f"""cat > "{script_path}" <<'NGS2_EOF'
set -euo pipefail
BATCH_NAME="{batch_name}"
OUT_DIR="{SECONDARY_DGX_OUTPUT_ROOT}/${{BATCH_NAME}}"
LAUNCH_DIR="{SECONDARY_DGX_LAUNCH_ROOT}/${{BATCH_NAME}}"
WORK_DIR="{SECONDARY_DGX_WORK_ROOT}/${{BATCH_NAME}}"
STAGED_SAMPLESHEET="{staged_sheet}"
NEXTFLOW_LOG="${{OUT_DIR}}/nextflow.log"

finish() {{
    status=$?
    set +e
    trap - EXIT
    if [ -f "${{LAUNCH_DIR}}/.nextflow.log" ]; then
        cp -f "${{LAUNCH_DIR}}/.nextflow.log" "${{NEXTFLOW_LOG}}"
        echo
        echo "[NGS2] Nextflow log copied to: ${{NEXTFLOW_LOG}}"
    else
        echo
        echo "[NGS2] Nextflow log not found at: ${{LAUNCH_DIR}}/.nextflow.log"
    fi
    if [ "${{status}}" -eq 0 ]; then
        echo "[NGS2] DONE: Nextflow finished successfully."
    else
        echo "[NGS2] FAILED: exit status ${{status}}"
    fi
    echo "[NGS2] tmux pane kept open. Type 'exit' to close this shell."
    exec bash -i
}}
trap finish EXIT

source "{SECONDARY_DGX_ENV_SCRIPT}"
mkdir -p "${{OUT_DIR}}" "${{LAUNCH_DIR}}" "${{WORK_DIR}}"
cp "${{STAGED_SAMPLESHEET}}" "${{OUT_DIR}}/samplesheet.csv"
cd "${{LAUNCH_DIR}}"

nextflow -c "${{PIPELINE_CONFIG}}" run "${{PIPELINE_CODE}}/main.nf" \\
    -profile {profile} \\
    --input_csv "${{OUT_DIR}}/samplesheet.csv" \\
    --seq_type {seq_type}{run_gcnv} \\
    --out_dir "${{OUT_DIR}}"{extended_analysis} \\
    -w "${{WORK_DIR}}" \\
    -resume
NGS2_EOF
tmux new-session -d -s "{session}" "bash {script_path}"
tmux attach -t "{session}"
"""


def cleanup_nf_work_command() -> dict:
    """Return a guarded cleanup script to run manually on DGX2."""
    root = SECONDARY_DGX_WORK_ROOT
    if not root.is_absolute() or str(root) in {"", "/"} or root == root.parent:
        raise ValueError(f"不安全的二級分析 Nextflow work root：{root}")
    root_text = str(root)
    command = f'''SECONDARY_WORK_ROOT={shlex.quote(root_text)}

if [ ! -d "${{SECONDARY_WORK_ROOT}}" ]; then
    echo "找不到二級分析 Nextflow 暫存路徑：${{SECONDARY_WORK_ROOT}}"
    exit 1
fi

if pgrep -af '[n]extflow' >/dev/null; then
    echo "偵測到 Nextflow 程序，停止清理："
    pgrep -af '[n]extflow'
    exit 1
fi

echo "準備清理 ${{SECONDARY_WORK_ROOT}}；以下項目將被刪除："
find "${{SECONDARY_WORK_ROOT}}" -mindepth 1 -maxdepth 1 -print
read -r -p "確定刪除以上二級分析 Nextflow 暫存？[y/N] " answer
case "${{answer}}" in
    y|Y|yes|YES)
        find "${{SECONDARY_WORK_ROOT}}" -mindepth 1 -maxdepth 1 -exec rm -rf -- {{}} +
        echo "已清理：${{SECONDARY_WORK_ROOT}}"
        ;;
    *)
        echo "已取消。"
        ;;
esac'''
    return {"path": root_text, "command": command}


def create_samplesheet(seq_type: str, samples: list[dict], batch_name: str = "") -> dict:
    seq = (seq_type or "").strip().upper()
    if seq not in {"WES", "WGS"}:
        raise ValueError("seq_type must be WES or WGS")
    if not samples:
        raise ValueError("samples must be a non-empty list")
    normalized = [_normalize_sample(s, seq) for s in samples]
    batch = (batch_name or "").strip() or suggest_batch_name(seq, normalized)
    _validate_sample_id(batch)

    output_dir = SECONDARY_OUTPUT_ROOT / batch
    staging_dir = SECONDARY_SAMPLESHEET_STAGING_ROOT / batch
    staging_dir.mkdir(parents=True, exist_ok=True)
    try:
        staging_dir.chmod(0o755)
    except OSError:
        pass

    has_lane = seq == "WGS" or any(s.get("lane") for s in normalized)
    fields = ["sample", "fastq_1", "fastq_2", "sex"] + (["lane"] if has_lane else [])
    sheet_rows: list[tuple[dict, dict | None]] = []
    for sample in normalized:
        if seq == "WGS":
            sheet_rows.extend((sample, lane) for lane in sample["lanes"])
        else:
            sheet_rows.append((sample, None))
    sheet = staging_dir / "samplesheet.csv"
    with sheet.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for sample, lane_payload in sheet_rows:
            row = {
                "sample": sample["sample_id"],
                "fastq_1": _server_path_to_dgx((lane_payload or sample)["fastq_1"]),
                "fastq_2": _server_path_to_dgx((lane_payload or sample)["fastq_2"]),
                "sex": sample["sex"] if sample["sex"] in {"male", "female", "unknown"} else "unknown",
            }
            if has_lane:
                row["lane"] = (lane_payload or sample).get("lane", "")
            writer.writerow(row)
    try:
        sheet.chmod(0o644)
    except OSError:
        pass

    dgx_output_dir = str(SECONDARY_DGX_OUTPUT_ROOT / batch)
    return {
        "batch_name": batch,
        "seq_type": seq,
        "sample_count": len(normalized),
        "samplesheet_row_count": len(sheet_rows),
        "samplesheet_path": str(sheet),
        "staged_samplesheet_path": str(sheet),
        "output_dir": str(output_dir),
        "dgx_output_dir": dgx_output_dir,
        "dgx_staged_samplesheet_path": str(SECONDARY_DGX_SAMPLESHEET_STAGING_ROOT / batch / "samplesheet.csv"),
        "tmux_session": f"ngs2_{batch}",
        "command": _launch_command(batch, seq),
        "warnings": [],
    }
