#!/usr/bin/env bash
# =========================================================
# download_giab_strata.sh — fetch the GIAB stratification BEDs
# the NGS-UI variant cards badge against, and write the manifest.
# =========================================================
# Downloads a curated subset of the GIAB / GA4GH genome
# stratifications (v3.1, GRCh38) into NGS_UI_GIAB_STRAT_DIR and
# writes strata_manifest.json so annotate_giab_strata.py knows which
# BED maps to which badge label.
#
# Run this ONCE on the dev machine (it needs outbound HTTPS). The BEDs
# are large-ish (~100s of MB total) and are NOT committed to git — they
# live under the reference tree like the MITOMAP / OMIM resources.
#
# Usage:
#   scripts/download_giab_strata.sh [--dir <out-dir>] [--force]
#
# Env:
#   NGS_UI_GIAB_STRAT_DIR   output dir (default
#                           /home/pipeline/reference/hg38/tertiary/giab_stratification)
#   GIAB_STRAT_BASE_URL     download base (default NCBI FTP v3.1/GRCh38)
# =========================================================
set -euo pipefail

OUT_DIR="${NGS_UI_GIAB_STRAT_DIR:-/home/pipeline/reference/hg38/tertiary/giab_stratification}"
BASE_URL="${GIAB_STRAT_BASE_URL:-https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.1/GRCh38}"
FORCE=0

while [ $# -gt 0 ]; do
  case "$1" in
    --dir)   OUT_DIR="$2"; shift 2;;
    --force) FORCE=1; shift;;
    -h|--help) sed -n '2,30p' "$0"; exit 0;;
    *) echo "unknown arg: $1" >&2; exit 2;;
  esac
done

# Each entry: <relative-path-under-BASE_URL>|<label>|<display>|<tooltip>
# label  → written into the TSV GIAB_STRATA column / drives the badge.
# display/tooltip → recorded in the manifest (frontend has its own map
#                   but keeping them here documents intent).
STRATA=(
  "LowComplexity/GRCh38_AllHomopolymers_gt6bp_imperfectgt10bp_slop5.bed.gz|homopolymer|Homopolymer|GIAB: homopolymer region — low-confidence indel calls"
  "LowComplexity/GRCh38_allTandemRepeats.bed.gz|tandem_repeat|Tandem repeat|GIAB: tandem / simple repeat region"
  "SegmentalDuplications/GRCh38_segdups.bed.gz|segdup|Segdup|GIAB: segmental duplication — possible false positives / mismapping"
  "mappability/GRCh38_lowmappabilityall.bed.gz|low_mappability|Low mappability|GIAB: low-mappability region — reads hard to place uniquely"
  "GCcontent/GRCh38_gclt30orgt55_slop50.bed.gz|gc_extreme|GC extreme|GIAB: extreme GC content (<30% or >55%) — coverage/calling bias"
  "OtherDifficult/GRCh38_allOtherDifficultregions.bed.gz|other_difficult|Other difficult|GIAB: other difficult region (MHC/KIR/VDJ, false dup, gaps, ...)"
)

mkdir -p "$OUT_DIR"
echo "GIAB stratification → $OUT_DIR"
echo "base URL: $BASE_URL"
echo

fetch() {
  local url="$1" dest="$2"
  if command -v wget >/dev/null 2>&1; then
    wget -q --show-progress -O "$dest" "$url"
  elif command -v curl >/dev/null 2>&1; then
    curl -fSL --progress-bar -o "$dest" "$url"
  else
    echo "ERROR: need wget or curl" >&2; exit 1
  fi
}

MANIFEST="$OUT_DIR/strata_manifest.json"
{
  echo "["
  first=1
  for entry in "${STRATA[@]}"; do
    IFS='|' read -r rel label display tooltip <<<"$entry"
    fname="$(basename "$rel")"
    dest="$OUT_DIR/$fname"
    if [ "$FORCE" -eq 1 ] || [ ! -s "$dest" ]; then
      echo "  ↓ $rel" >&2
      fetch "$BASE_URL/$rel" "$dest"
    else
      echo "  ✓ $fname (exists, skip; --force to re-download)" >&2
    fi
    [ "$first" -eq 1 ] || echo ","
    first=0
    printf '  {"file": "%s", "label": "%s", "display": "%s", "tooltip": "%s"}' \
      "$fname" "$label" "$display" "$tooltip"
  done
  echo
  echo "]"
} >"$MANIFEST"

echo
echo "manifest → $MANIFEST"
echo "done. Re-run scripts/run_stopgaps.sh (or re-analyse) to populate GIAB_STRATA."
