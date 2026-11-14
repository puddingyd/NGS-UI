#!/usr/bin/env python3
"""Hospital-side cohort preflight, isolated scoring and research reports."""
from pathlib import Path
import sys

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO/'backend'))

from research.prs.cohort_validation import main

if __name__ == '__main__':
    raise SystemExit(main())
