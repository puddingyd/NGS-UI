#!/usr/bin/env python3
"""Prepare three WES rerun batches using the UI's current secondary workflow.

Run on the DGM host, where the WES FASTQ roots and DGX2 staging mount exist.
This script prepares files only; it never starts Nextflow or connects to DGX2.
"""

from __future__ import annotations

import argparse
import sys
from collections import defaultdict
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from backend.app.services import secondary_analysis as secondary  # noqa: E402


BATCHES = {
    "23WE": "WES_23WE_re",
    "24WE": "WES_24WE_re",
    "25WE": "WES_25WE_re",
}


def select_samples(discovered: list[dict]) -> dict[str, list[dict]]:
    """Choose complete WES pairs and reject ambiguous sample IDs."""
    grouped: dict[str, dict[str, dict]] = {prefix: {} for prefix in BATCHES}
    conflicts: dict[str, list[str]] = defaultdict(list)
    for entry in discovered:
        sample_id = str(entry.get("sample_id") or "")
        prefix = next((p for p in BATCHES if sample_id.startswith(p)), None)
        if prefix is None:
            continue
        f1 = Path(entry["fastq_1"])
        f2 = Path(entry["fastq_2"])
        if entry.get("reanalysis"):
            m1 = secondary._REANALYSIS_RE.fullmatch(f1.name)
            m2 = secondary._REANALYSIS_RE.fullmatch(f2.name)
        else:
            m1 = secondary._NEXTSEQ_RE.fullmatch(f1.name)
            m2 = secondary._NEXTSEQ_RE.fullmatch(f2.name)
        if not m1 or not m2 or m1.group(1) != m2.group(1) or m1.group(2) != "1" or m2.group(2) != "2":
            raise ValueError(f"FASTQ R1/R2 不屬於同一來源樣本：{sample_id}: {f1}, {f2}")
        previous = grouped[prefix].get(sample_id)
        if previous:
            old_pair = (Path(previous["fastq_1"]).resolve(), Path(previous["fastq_2"]).resolve())
            new_pair = (f1.resolve(), f2.resolve())
            if old_pair != new_pair:
                conflicts[sample_id].extend((str(previous["fastq_1"]), str(f1)))
            continue
        grouped[prefix][sample_id] = entry

    if conflicts:
        details = "; ".join(f"{sid}: {', '.join(paths)}" for sid, paths in sorted(conflicts.items()))
        raise ValueError(f"同一 Sample ID 找到不同 FASTQ pair，需先選定來源：{details}")
    empty = [prefix for prefix, entries in grouped.items() if not entries]
    if empty:
        raise ValueError(f"找不到完整配對的 WES FASTQ：{', '.join(empty)}；未建立任何 samplesheet")
    return {prefix: [entries[sid] for sid in sorted(entries)] for prefix, entries in grouped.items()}


def prepare(*, dry_run: bool = False, rescan: bool = False) -> list[dict]:
    existing_roots = [root for root in secondary.SECONDARY_WES_FASTQ_ROOTS if root.is_dir()]
    if not existing_roots:
        raise FileNotFoundError("找不到 WES FASTQ roots；請在掛載原始資料的 DGM 主機執行")

    if rescan:
        print("正在重新掃描 WES FASTQ roots；大型原始資料目錄可能需要數分鐘以上...", flush=True)
        discovered = secondary.list_wes_fastqs(strict_prefixes=tuple(BATCHES))
    else:
        index = secondary.load_index()
        if secondary.index_is_stale(index):
            raise ValueError("二級 FASTQ 索引不存在或已過期；請先更新 UI 索引，或加 --rescan 直接掃描")
        print(f"使用二級 FASTQ 索引（更新時間：{index['updated_at']}）；若需檢查新檔或缺配對檔，請加 --rescan", flush=True)
        discovered = index.get("wes") or []
    grouped = select_samples(discovered)
    for samples in grouped.values():
        for sample in samples:
            secondary._validate_fastq_pair(Path(sample["fastq_1"]), Path(sample["fastq_2"]))
    for batch in BATCHES.values():
        staging = secondary.SECONDARY_SAMPLESHEET_STAGING_ROOT / batch
        output = secondary.SECONDARY_OUTPUT_ROOT / batch
        if staging.exists() or output.exists():
            raise FileExistsError(f"批次已存在，停止以避免覆寫：{staging if staging.exists() else output}")

    results = []
    for prefix, batch in BATCHES.items():
        samples = grouped[prefix]
        if dry_run:
            results.append({"batch_name": batch, "sample_count": len(samples), "samples": samples})
            continue
        result = secondary.create_samplesheet("WES", samples, batch_name=batch)
        launch_file = Path(result["samplesheet_path"]).with_name("run_on_dgx2.sh")
        launch_file.write_text("#!/usr/bin/env bash\n" + result["command"], encoding="utf-8")
        launch_file.chmod(0o644)
        result["launch_file"] = str(launch_file)
        result["dgx_launch_file"] = str(secondary.SECONDARY_DGX_SAMPLESHEET_STAGING_ROOT / batch / launch_file.name)
        results.append(result)
    return results


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dry-run", action="store_true", help="只列出各批樣本，不寫入檔案")
    parser.add_argument("--rescan", action="store_true", help="不用每日索引，重新掃描全部 WES FASTQ（可能較慢）")
    args = parser.parse_args()
    try:
        results = prepare(dry_run=args.dry_run, rescan=args.rescan)
    except (FileNotFoundError, FileExistsError, ValueError, OSError) as exc:
        parser.exit(1, f"錯誤：{exc}\n")
    for result in results:
        print(f"{result['batch_name']}: {result['sample_count']} samples")
        if args.dry_run:
            for sample in result["samples"]:
                print(f"  {sample['sample_id']}  {sample['fastq_1']}  {sample['fastq_2']}")
        else:
            print(f"  samplesheet: {result['samplesheet_path']}")
            print(f"  DGX2 貼上: bash {result['dgx_launch_file']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
