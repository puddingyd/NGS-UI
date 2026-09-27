#!/usr/bin/env python3
"""Backfill INHOUSE_AC/AN/AF into already-analysed samples.

Run this for samples analysed before the in-house AF step existed, and after
every refresh of the in-house AF DB (the cohort grows, so AN and every AF
change). New analyses don't need it — run_stopgaps.sh does the same work
inline.

It mirrors what the tertiary worker does. Which of the two modes applies is
decided by THE RESOLVED FILE, not by sample_layout.uses_unified_layout(): a
"legacy" sample can still resolve to a pipeline 03_acmg file, because
snv_raw_tsv() falls back to <legacy pipeline root>/<sample>/03_acmg/… when no
UI copy exists. Only a file literally named snv_indel.annotated.tsv outside
03_acmg is the UI's own copy and safe to rewrite; do_inplace() additionally
refuses any path containing 03_acmg.

OVERLAY mode — the resolved TSV is a pipeline source (03_acmg, under either
the unified or the legacy pipeline root):
    03_acmg/<source>.snv_indel.acmg.tsv   IMMUTABLE pipeline source of truth;
                                          exported reports read it. NEVER written.
    08_postprocessing/<sid>.snv_annotations.sqlite
                                          sparse overlay holding every field
                                          post-processing adds (GeneBe,
                                          SpliceAI, MANE, LitVar2, INHOUSE_*)
    08_postprocessing/<sid>.snv_indel.review.tsv
                                          main-screen rows, built from raw+overlay
  So the steps are: materialise raw+overlay -> a temp working TSV, add
  INHOUSE_* to THAT, re-diff it against raw into the overlay, rebuild the
  review TSV, delete the working file. The gene index is built from the raw
  TSV, which never changes, so it is NOT rebuilt here.

  Materialising first is essential: build_overlay() REPLACES the overlay, so
  diffing raw against a fresh copy carrying only INHOUSE_* would silently drop
  every other annotation.

IN-PLACE mode — the resolved TSV is the legacy UI copy
(<NGS_UI_HOME>/tertiary_output/<sample>/snv_indel.annotated.tsv):
    that file is itself the annotated copy, so it is rewritten in place and the
    gene index MUST be rebuilt with it — the atomic replace shifts every byte
    offset and the index is only auto-rebuilt when missing.

Derived artifacts (overlay, review TSV, gene index) always live where
state_dir() says, which for most samples is the unified 08_postprocessing even
when the source TSV still sits under the legacy pipeline root.

Paths come from backend/app/services/sample_layout.py, so this follows the
backend exactly.

Usage:
    scripts/backfill_inhouse_af.py                  # every sample
    scripts/backfill_inhouse_af.py SID1 SID2 ...    # only these
    scripts/backfill_inhouse_af.py --dry-run        # show plan, run nothing
"""
from __future__ import annotations

import argparse
import csv
import signal
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "backend"))

from app import config  # noqa: E402
from app.services import sample_layout as layout  # noqa: E402
from app.services.snv_overlay import OverlayReader  # noqa: E402

SCRIPTS = REPO / "scripts"
csv.field_size_limit(4 * 1024 * 1024)


def derived_paths(sample_id: str):
    """Same names run_stopgaps.sh uses: sample-prefixed under 08_postprocessing,
    bare in the legacy tree."""
    post = layout.state_dir(sample_id, for_write=True)
    pre = f"{sample_id}." if post.name == layout.POSTPROCESSING_DIRNAME else ""
    return {
        "post": post,
        "review": post / f"{pre}snv_indel.review.tsv",
        "manifest": post / f"{pre}snv_indel.review.tsv.source.json",
        "gene_index": post / f"{pre}snv_gene_index.sqlite",
        "overlay": post / f"{pre}snv_annotations.sqlite",
    }


def materialise(raw: Path, overlay: Path, out: Path) -> None:
    """raw + existing overlay -> full annotated TSV, so the re-diff preserves
    annotations this backfill does not touch."""
    have = overlay.is_file()
    with OverlayReader(raw, overlay if have else None) as rd, \
            open(raw, "r", encoding="utf-8", newline="") as fi, \
            open(out, "w", encoding="utf-8", newline="") as fo:
        rdr = csv.DictReader(fi, delimiter="\t")
        fields = list(rdr.fieldnames or [])
        fields += [f for f in (rd.fields if rd.active else []) if f not in fields]
        wtr = csv.DictWriter(fo, fieldnames=fields, delimiter="\t",
                             extrasaction="ignore", lineterminator="\n")
        wtr.writeheader()
        if rd.active:
            for row in rdr:
                wtr.writerow(rd.apply(row))
        else:
            wtr.writerows(rdr)


def is_ui_copy(raw: Path) -> bool:
    """True only for the legacy UI *annotated copy*, which is ours to rewrite.

    Do NOT decide this from sample_layout.uses_unified_layout(): a sample can
    be 'legacy' yet still resolve to a pipeline 03_acmg file, because
    snv_raw_tsv() falls back to <legacy pipeline root>/<sample>/03_acmg/… when
    no UI copy exists. Anything under 03_acmg is a pipeline source of truth
    (exported reports read it) and must never be modified."""
    return raw.name == "snv_indel.annotated.tsv" and "03_acmg" not in raw.parts


def run(cmd) -> None:
    subprocess.run([str(c) for c in cmd], check=True)


def do_overlay(sid: str, raw: Path, p: dict, db: Path) -> None:
    work = p["post"] / f".snv_indel.backfill.working.tsv"
    try:
        materialise(raw, p["overlay"], work)
        run([SCRIPTS / "annotate_inhouse_af.py", "--tsv", work, "--db", db])
        run([SCRIPTS / "build_snv_annotation_overlay.py",
             "--raw", raw, "--annotated", work, "--out", p["overlay"]])
        run([SCRIPTS / "build_snv_review_tsv.py", "--tsv", raw,
             "--output-dir", p["post"], "--output-path", p["review"],
             "--manifest-path", p["manifest"], "--overlay", p["overlay"]])
    finally:
        work.unlink(missing_ok=True)


def do_inplace(sid: str, raw: Path, p: dict, db: Path) -> None:
    # belt and braces: never rewrite a pipeline source of truth
    if "03_acmg" in raw.parts:
        raise RuntimeError(f"refusing to rewrite pipeline source in place: {raw}")
    run([SCRIPTS / "annotate_inhouse_af.py", "--tsv", raw, "--db", db])
    run([SCRIPTS / "build_snv_review_tsv.py", "--tsv", raw,
         "--output-dir", p["post"], "--output-path", p["review"],
         "--manifest-path", p["manifest"], "--overlay", p["overlay"]])
    # raw was rewritten in place -> every byte offset moved
    run([SCRIPTS / "build_snv_gene_index.py", "--tsv", raw, "--out", p["gene_index"]])


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("samples", nargs="*", help="sample IDs (default: all)")
    ap.add_argument("--db", default=str(config.INHOUSE_AF_DB),
                    help=f"in-house AF sites VCF (default {config.INHOUSE_AF_DB})")
    ap.add_argument("--dry-run", action="store_true",
                    help="list each sample, its layout and resolved paths; run nothing")
    ap.add_argument("--continue-on-error", action="store_true",
                    help="keep going when one sample fails (default: stop)")
    args = ap.parse_args()
    try:               # behave like a normal tool under `| head`
        signal.signal(signal.SIGPIPE, signal.SIG_DFL)
    except (AttributeError, ValueError):
        pass

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
        p = derived_paths(sid)
        inplace = is_ui_copy(raw)
        if args.dry_run:
            kind = "legacy UI copy" if inplace else "pipeline source + overlay"
            fate = "REWRITTEN IN PLACE" if inplace else "read-only"
            print(f"  • {sid}  [{kind}]")
            print(f"      raw      {raw}   ({fate})")
            print(f"      overlay  {p['overlay']}")
            print(f"      review   {p['review']}")
            if inplace:
                print(f"      index    {p['gene_index']}   (rebuilt)")
            n_ok += 1
            continue
        p["post"].mkdir(parents=True, exist_ok=True)
        print(f"  • {sid}  [{'in-place' if inplace else 'overlay'}]")
        try:
            (do_inplace if inplace else do_overlay)(sid, raw, p, db)
            n_ok += 1
        except subprocess.CalledProcessError as e:
            n_fail += 1
            print(f"  ! {sid}: {Path(e.cmd[0]).name} failed (exit {e.returncode})",
                  file=sys.stderr)
            if not args.continue_on_error:
                print("stopping (use --continue-on-error to keep going)", file=sys.stderr)
                break
        except Exception as e:  # noqa: BLE001
            n_fail += 1
            print(f"  ! {sid}: {e}", file=sys.stderr)
            if not args.continue_on_error:
                break

    print(f"done. {n_ok} annotated, {n_skip} skipped, {n_fail} failed.")
    return 1 if n_fail else 0


if __name__ == "__main__":
    raise SystemExit(main())
