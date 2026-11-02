#!/usr/bin/env python3
"""Link a completed Method2 summary to an existing NGS-UI case (operator only).

PYTHONPATH=backend python scripts/import_prs_summary.py --sample-id ID --summary /path/summary.json
"""
import argparse
import json
from pathlib import Path

from app.services import prs


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sample-id', required=True)
    parser.add_argument('--summary', required=True, type=Path)
    args = parser.parse_args()
    try:
        print(json.dumps(prs.import_score(args.sample_id, args.summary), ensure_ascii=False))
    except (ValueError, OSError) as exc:
        parser.exit(2, f'ERROR: {exc}\n')


if __name__ == '__main__':
    main()
