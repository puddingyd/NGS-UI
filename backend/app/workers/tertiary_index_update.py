"""Refresh the tertiary-analysis VCF discovery index for systemd."""
from __future__ import annotations

from ..services import dragen_jobs


def main() -> int:
    idx = dragen_jobs.refresh_index()
    print(
        "tertiary index refreshed: "
        f"dragen={len(idx.get('dragen', []))} "
        f"inhouse={len(idx.get('inhouse', []))} "
        f"duration={idx.get('scan_duration_sec')}s "
        f"updated_at={idx.get('updated_at')}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
