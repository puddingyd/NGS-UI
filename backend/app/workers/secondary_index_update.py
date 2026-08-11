"""Refresh the secondary-analysis FASTQ discovery index for systemd."""
from __future__ import annotations

from ..services import secondary_analysis


def main() -> int:
    idx = secondary_analysis.refresh_index()
    lane_count = sum(int(row.get("lane_count") or 1) for row in idx.get("wgs", []))
    print(
        "secondary index refreshed: "
        f"wes={len(idx.get('wes', []))} "
        f"wgs={len(idx.get('wgs', []))} "
        f"wgs_lanes={lane_count} "
        f"duration={idx.get('scan_duration_sec')}s "
        f"updated_at={idx.get('updated_at')}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
