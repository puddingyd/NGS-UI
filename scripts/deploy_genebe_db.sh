#!/usr/bin/env bash
# =========================================================
# deploy_genebe_db.sh — clean + index + atomically install the GeneBe DB
# =========================================================
# The GeneBe DB (genebe_hg38.tsv.gz) is rebuilt periodically and can
# contain malformed placeholder rows (e.g. "chr.\t.\t.\t…" where chr/pos
# are '.'). Those break `tabix` indexing (pos must be an integer), and a
# stale .tbi (older than the .gz) breaks tabix random access. This script
# makes deploying a fresh build a one-liner:
#
#   1. verify the source bgzip stream (bgzip -t)
#   2. drop rows whose pos (col 2) isn't an integer  (--sort to also sort)
#   3. bgzip the cleaned TSV
#   4. tabix -s 1 -b 2 -e 2 -c '#'  (fresh index, newer than the data)
#   5. sanity-probe a known region
#   6. atomically install to the destination (.gz then .tbi)
#
# NGS-UI itself streams the .gz (it does NOT use the .tbi), so a bad
# index can't break the platform — but a clean, indexed DB is still
# wanted for tabix-based tools and the GENEBE_USAGE.md workflow.
#
# Usage:
#   scripts/deploy_genebe_db.sh [SRC.tsv.gz] [--dest DEST] [--sort] [--keep-src]
#
#   SRC   freshly built/downloaded DB (default: the DEST itself = clean
#         it in place). DEST default = $NGS_UI_GENEBE_DB or
#         $HOME/NGS_UI/biotools/genebe/genebe_hg38.tsv.gz.
# =========================================================
set -euo pipefail

DEST="${NGS_UI_GENEBE_DB:-$HOME/NGS_UI/biotools/genebe/genebe_hg38.tsv.gz}"
SRC=""
SORT=0
KEEP_SRC=0
PROBE_REGION="chr7:140753336-140753336"

while [ $# -gt 0 ]; do
  case "$1" in
    --dest)     DEST="$2"; shift 2;;
    --sort)     SORT=1; shift;;
    --keep-src) KEEP_SRC=1; shift;;
    --probe)    PROBE_REGION="$2"; shift 2;;
    -h|--help)  sed -n '2,38p' "$0"; exit 0;;
    -*)         echo "unknown arg: $1" >&2; exit 2;;
    *)          SRC="$1"; shift;;
  esac
done
[ -n "$SRC" ] || SRC="$DEST"

for bin in bgzip tabix awk; do
  command -v "$bin" >/dev/null 2>&1 || { echo "ERROR: need '$bin' on PATH" >&2; exit 1; }
done
[ -s "$SRC" ] || { echo "ERROR: source not found / empty: $SRC" >&2; exit 1; }

echo "GeneBe DB deploy"
echo "  source: $SRC"
echo "  dest:   $DEST"

echo "[1/5] verifying source bgzip integrity…"
bgzip -t "$SRC"

DEST_DIR="$(dirname "$DEST")"
mkdir -p "$DEST_DIR"
# Work in the destination dir so the final mv is an atomic same-fs rename.
TMP="$(mktemp "$DEST_DIR/.genebe_deploy.XXXXXX")"
TMP_GZ="$TMP.tsv.gz"
cleanup() { rm -f "$TMP" "$TMP_GZ" "$TMP_GZ.tbi"; }
trap cleanup EXIT

echo "[2/5] cleaning malformed rows (pos must be integer)$( [ "$SORT" = 1 ] && echo ' + sorting')…"
if [ "$SORT" = 1 ]; then
  { zcat "$SRC" | head -1
    zcat "$SRC" | awk -F'\t' 'NR>1 && $2 ~ /^[0-9]+$/' | sort -k1,1 -k2,2n
  } | awk -F'\t' '
      NR==1 {print; next}
      {kept++; print}
      END {printf("      kept=%d\n", kept) > "/dev/stderr"}
    ' | bgzip -c > "$TMP_GZ"
else
  zcat "$SRC" | awk -F'\t' '
      NR==1 {print; header=1; next}
      $2 ~ /^[0-9]+$/ {kept++; print; next}
      {dropped++}
      END {printf("      kept=%d dropped=%d\n", kept, dropped) > "/dev/stderr"}
    ' | bgzip -c > "$TMP_GZ"
fi

echo "[3/5] indexing (tabix -s 1 -b 2 -e 2 -c '#')…"
if ! tabix -s 1 -b 2 -e 2 -c '#' "$TMP_GZ" 2> "$TMP.tabixerr"; then
  cat "$TMP.tabixerr" >&2
  rm -f "$TMP.tabixerr"
  echo "ERROR: tabix indexing failed — if it says 'not sorted', re-run with --sort" >&2
  exit 1
fi
rm -f "$TMP.tabixerr"

echo "[4/5] sanity probe: $PROBE_REGION"
if ! tabix "$TMP_GZ" "$PROBE_REGION" | head -3; then
  echo "WARNING: probe returned nothing (region may genuinely be absent)" >&2
fi

echo "[5/5] installing atomically…"
# .gz first, then .tbi. The brief window (new .gz + old .tbi) fails safe
# for any tabix consumer that checks index freshness.
mv -f "$TMP_GZ"     "$DEST"
mv -f "$TMP_GZ.tbi" "$DEST.tbi"
trap - EXIT
rm -f "$TMP"
[ "$KEEP_SRC" = 1 ] || [ "$SRC" = "$DEST" ] || true

echo "done → $DEST (+ .tbi)"
ls -la "$DEST" "$DEST.tbi"
