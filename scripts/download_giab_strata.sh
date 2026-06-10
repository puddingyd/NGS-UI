#!/usr/bin/env bash
# =========================================================
# download_giab_strata.sh — fetch the GIAB stratification BEDs
# the NGS-UI variant cards badge against, and write the manifest.
# =========================================================
# Works from the GIAB "@all" GRCh38 stratification tarball (default
# v3.6). It downloads (or reuses a local copy of) the tarball, extracts
# it, then picks a curated subset of BEDs by *substring pattern* — so it
# keeps working across GIAB versions even though the exact slop/threshold
# filenames change between releases. The selected BEDs are copied into
# NGS_UI_GIAB_STRAT_DIR and a strata_manifest.json (file -> badge label) is
# written for scripts/annotate_giab_strata.py.
#
# Run this ONCE on the dev machine (needs outbound HTTPS, tar, ~GB temp).
# The BEDs are NOT committed to git — they live under the reference tree
# like the MITOMAP / OMIM resources.
#
# Usage:
#   scripts/download_giab_strata.sh [--dir <out-dir>] \
#       [--tarball <local.tar.gz>] [--url <tarball-url>] [--keep] [--force]
#
# Env:
#   NGS_UI_GIAB_STRAT_DIR      output dir (default
#                              /home/pipeline/reference/hg38/tertiary/giab_stratification)
#   GIAB_STRAT_TARBALL_URL     tarball URL (default v3.6 GRCh38 @all)
# =========================================================
set -euo pipefail

OUT_DIR="${NGS_UI_GIAB_STRAT_DIR:-/home/pipeline/reference/hg38/tertiary/giab_stratification}"
TARBALL_URL="${GIAB_STRAT_TARBALL_URL:-https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.6/genome-stratifications-GRCh38@all.tar.gz}"
LOCAL_TARBALL=""
KEEP_EXTRACT=0
FORCE=0

while [ $# -gt 0 ]; do
  case "$1" in
    --dir)     OUT_DIR="$2"; shift 2;;
    --tarball) LOCAL_TARBALL="$2"; shift 2;;
    --url)     TARBALL_URL="$2"; shift 2;;
    --keep)    KEEP_EXTRACT=1; shift;;
    --force)   FORCE=1; shift;;
    -h|--help) sed -n '2,33p' "$0"; exit 0;;
    *) echo "unknown arg: $1" >&2; exit 2;;
  esac
done

# Each entry: <label>|<display>|<tooltip>|<pattern1>[ <pattern2> ...]
# Patterns are case-insensitive shell globs matched against extracted BED
# basenames; the first pattern that matches wins, and among matches the
# SHORTEST basename is taken (a good proxy for the broad "union" file over
# its sub-stratified variants). Stable substrings (AllHomopolymers,
# AllTandemRepeats, segdups, lowmappabilityall) have held since v2.0.
STRATA=(
  "homopolymer|Homopolymer|GIAB: homopolymer region — low-confidence indel calls|*AllHomopolymers*slop5.bed.gz *AllHomopolymers*.bed.gz"
  "tandem_repeat|Tandem repeat|GIAB: tandem / simple repeat region|*AllTandemRepeats.bed.gz *AllTandemRepeats*.bed.gz *allTandemRepeats*.bed.gz"
  "segdup|Segdup|GIAB: segmental duplication — possible mismapping / false positives|*segdups.bed.gz *SegDup*.bed.gz"
  "low_mappability|Low mappability|GIAB: low-mappability region — reads hard to place uniquely|*lowmappabilityall.bed.gz *lowmappability*.bed.gz"
  "gc_extreme|GC extreme|GIAB: extreme GC content — coverage / calling bias|*gclt*orgt*_slop50.bed.gz *gclt*orgt*.bed.gz"
  "other_difficult|Other difficult|GIAB: other difficult region (MHC/KIR/VDJ, false dup, gaps, ...)|*allOtherDifficult*.bed.gz *OtherDifficult*.bed.gz"
)

mkdir -p "$OUT_DIR"
echo "GIAB stratification → $OUT_DIR"

# --- obtain tarball -------------------------------------------------
WORK="$(mktemp -d "${TMPDIR:-/tmp}/giab-strat.XXXXXX")"
cleanup() { [ "$KEEP_EXTRACT" -eq 1 ] || rm -rf "$WORK"; }
trap cleanup EXIT

if [ -n "$LOCAL_TARBALL" ]; then
  TARBALL="$LOCAL_TARBALL"
  echo "tarball: $TARBALL (local)"
else
  TARBALL="$WORK/$(basename "$TARBALL_URL")"
  echo "tarball: $TARBALL_URL"
  if command -v wget >/dev/null 2>&1; then
    wget -q --show-progress -O "$TARBALL" "$TARBALL_URL"
  elif command -v curl >/dev/null 2>&1; then
    curl -fSL --progress-bar -o "$TARBALL" "$TARBALL_URL"
  else
    echo "ERROR: need wget or curl" >&2; exit 1
  fi
fi
[ -s "$TARBALL" ] || { echo "ERROR: tarball missing/empty: $TARBALL" >&2; exit 1; }

# --- extract --------------------------------------------------------
EXTRACT="$WORK/extract"
mkdir -p "$EXTRACT"
echo "extracting…"
tar -xzf "$TARBALL" -C "$EXTRACT"

# --- select BEDs by pattern + write manifest ------------------------
# pick_shortest <dir> <glob>  → echoes the shortest-basename match (or empty)
pick_shortest() {
  local dir="$1" glob="$2"
  find "$dir" -type f -iname "$glob" -printf '%f\t%p\n' 2>/dev/null \
    | awk -F'\t' '{ if (n=="" || length($1)<min) { min=length($1); n=$1; p=$2 } } END { if (p!="") print p }'
}

MANIFEST="$OUT_DIR/strata_manifest.json"
tmp_manifest="$WORK/manifest.json"
{
  echo "["
  first=1
  for entry in "${STRATA[@]}"; do
    IFS='|' read -r label display tooltip patterns <<<"$entry"
    match=""
    for glob in $patterns; do
      match="$(pick_shortest "$EXTRACT" "$glob")"
      [ -n "$match" ] && break
    done
    if [ -z "$match" ]; then
      echo "  ✗ $label: no BED matched ($patterns) — skipped" >&2
      continue
    fi
    fname="$(basename "$match")"
    dest="$OUT_DIR/$fname"
    if [ "$FORCE" -eq 1 ] || [ ! -s "$dest" ]; then
      cp -f "$match" "$dest"
    fi
    echo "  ✓ $label → $fname" >&2
    [ "$first" -eq 1 ] || echo ","
    first=0
    printf '  {"file": "%s", "label": "%s", "display": "%s", "tooltip": "%s"}' \
      "$fname" "$label" "$display" "$tooltip"
  done
  echo
  echo "]"
} >"$tmp_manifest"
mv "$tmp_manifest" "$MANIFEST"

echo
echo "manifest → $MANIFEST"
echo "selected BEDs:"
ls -la "$OUT_DIR"/*.bed.gz 2>/dev/null || true
echo
echo "done. Re-run scripts/run_stopgaps.sh (or re-analyse) to populate GIAB_STRATA."
