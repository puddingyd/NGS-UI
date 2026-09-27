#!/usr/bin/env python3
"""Backfill INHOUSE_AC/AN/AF into already-analysed samples.

Run this for samples analysed before the in-house AF step existed, and after
every refresh of the in-house AF DB (the cohort grows, so AN and every AF
change). New analyses don't need it — run_stopgaps.sh runs the same annotate
step inline.

For each sample it does the three things that must happen together:

  1. annotate_inhouse_af.py   fill/refresh INHOUSE_AC/AN/AF on the pipeline SNV
                              TSV (03_acmg in the unified layout)
  2. build_snv_review_tsv.py  rebuild the main-screen review TSV so the new
                              columns reach the UI
  3. build_snv_gene_index.py  REBUILD the gene index — step 1 rewrites the whole
                              TSV (atomic replace), shifting every byte offset,
                              and the index is only auto-rebuilt when missing.
                              Skipping this makes gene search seek stale bytes.

Paths come from backend/app/services/sample_layout.py, so this follows the
unified layout (<TERTIARY_ROOT>/<sample>/03_acmg + 08_postprocessing) and the
legacy UI tree identically to the backend. Derived filenames are sample-
prefixed under 08_postprocessing, matching run_stopgaps.sh.

--test-type and --gpn-msa-db are deliberately NOT passed through: the review
builder infers the test type from the sample metadata (run_stopgaps passes an
explicit WES default, which would be wrong for WGS samples here) and reads
NGS_UI_GPN_MSA_DB itself.

Usage:
    scripts/backfill_inhouse_af.py                  # every sample
    scripts/backfill_inhouse_af.py SID1 SID2 ...    # only these
    scripts/backfill_inhouse_af.py --dry-run        # show what would run
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "backend"))

from app import config  # noqa: E402
from app.services import sample_layout as layout  # noqa: E402

SCRIPTS = REPO / "scripts"


def derived_paths(sample_id: str):
    """(post_dir, review, manifest, gene_index, overlay) — same names as
    run_stopgaps.sh, i.e. sample-prefixed under 08_postprocessing and bare in
    the legacy tree."""
    post = layout.state_dir(sample_id, for_write=True)
    pre = f"{sample_id}." if post.name == layout.POSTPROCESSING_DIRNAME else ""
    return (
        post,
        post / f"{pre}snv_indel.review.tsv",
        post / f"{pre}snv_indel.review.tsv.source.json",
        post / f"{pre}snv_gene_index.sqlite",
        post / f"{pre}snv_annotations.sqlite",
    )


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("samples", nargs="*", help="sample IDs (default: all)")
    ap.add_argument("--db", default=str(config.INHOUSE_AF_DB),
                    help=f"in-house AF sites VCF (default {config.INHOUSE_AF_DB})")
    ap.add_argument("--dry-run", action="store_true",
                    help="list the samples and resolved paths, run nothing")
    ap.add_argument("--continue-on-error", action="store_true",
                    help="keep going when one sample fails (default: stop)")
    args = ap.parse_args()

    db = Path(args.db)
    if not db.is_file():
        print(f"ERROR: in-house AF DB not found: {db}\n"
              f"       deploy it first with scripts/inhouse_af/deploy_inhouse_af_db.sh",
              file=sys.stderr)
        return 1

    sids = args.samples or sorted(layout.iter_sample_ids())
    if not sids:
        print(f"no samples under {layout.unified_root()}")
        return 0

    print(f"backfilling INHOUSE_AF for {len(sids)} sample(s)")
    print(f"  DB   : {db}")
    print(f"  root : {layout.unified_root()}")

    n_ok = n_skip = n_fail = 0
    for sid in sids:
        raw = layout.snv_raw_tsv(sid)
        if not raw.is_file():
            print(f"  - {sid}: no SNV TSV ({raw}), skip", file=sys.stderr)
            n_skip += 1
            continue
        post, review, manifest, gene_index, overlay = derived_paths(sid)
        if args.dry_run:
            print(f"  • {sid}\n      tsv   {raw}\n      post  {post}")
            n_ok += 1
            continue
        post.mkdir(parents=True, exist_ok=True)
        steps = [
            [SCRIPTS / "annotate_inhouse_af.py", "--tsv", raw, "--db", db],
            [SCRIPTS / "build_snv_review_tsv.py", "--tsv", raw,
             "--output-dir", post, "--output-path", review,
             "--manifest-path", manifest, "--overlay", overlay],
            [SCRIPTS / "build_snv_gene_index.py", "--tsv", raw, "--out", gene_index],
        ]
        print(f"  • {sid}")
        try:
            for cmd in steps:
                subprocess.run([str(c) for c in cmd], check=True)
            n_ok += 1
        except subprocess.CalledProcessError as e:
            n_fail += 1
            print(f"  ! {sid}: {Path(e.cmd[0]).name} failed (exit {e.returncode})",
                  file=sys.stderr)
            if not args.continue_on_error:
                print("stopping (use --continue-on-error to keep going)", file=sys.stderr)
                break

    print(f"done. {n_ok} annotated, {n_skip} skipped, {n_fail} failed.")
    return 1 if n_fail else 0


if __name__ == "__main__":
    raise SystemExit(main())
