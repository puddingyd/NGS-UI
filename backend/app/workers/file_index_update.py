"""Daily shared FASTQ/BAM/VCF inventory update (also used by old workers)."""
from datetime import datetime, timedelta
from ..services import file_index


def main():
    now = datetime.now(file_index.TAIPEI)
    scheduled = now.replace(hour=2, minute=0, second=0, microsecond=0)
    if scheduled > now:
        scheduled -= timedelta(days=1)
    # Even when both legacy timers are still installed, each daily window
    # produces one generation, including when their processes don't overlap.
    views = file_index.refresh(force=True, requested_at_ns=int(scheduled.timestamp()*1e9))
    meta = views["meta"]
    print(f"shared file index refreshed: samples={meta['sample_count']} "
          f"generation={meta['generation']} updated_at={meta['updated_at']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
