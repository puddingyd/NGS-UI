#!/usr/bin/env bash
# =========================================================
# update_inhouse_af.sh — fold a new sequencing batch into the in-house AF DB
# =========================================================
# Per-batch SOP (the "每批下機後更新" workflow). Maintains a persistent
# cohort under --cohort-dir:
#
#     <cohort-dir>/cohort_manifest.tsv   audit: every sample ever considered
#     <cohort-dir>/cohort_gvcfs.txt      included gVCF paths (the cohort)
#     <cohort-dir>/inhouse_af.hg38.vcf.gz[.tbi]   the AF DB (atomic-swapped)
#     <cohort-dir>/updates.log           one line per update (when/added/total)
#
# What it does:
#   1. select_cohort.py on the NEW batch (same exclusion rules) → candidates.
#   2. APPEND only samples whose sample_id is not already in the manifest
#      (dedup guard — re-running a batch never double-counts).
#   3. Re-genotype the FULL updated cohort with build_inhouse_af.sh
#      (always-correct baseline; see README "Incremental update design").
#   4. Atomic-swap the new AF VCF into place.
#
# Bootstrap (Phase 0 → first build): point --new-sizes at the whole cohort's
# size list with an empty cohort-dir; it becomes the initial DB.
#
# Usage:
#   scripts/inhouse_af/update_inhouse_af.sh \
#     --cohort-dir $NGS_UI_HOME/biotools/inhouse_af \
#     --new-sizes  new_batch_sizes.txt \
#     [--exclude-range 'VAL-:37-54'] [--exclude-id-file controls.txt] \
#     [--threads 32] [--mem-gbytes 192] \
#     [--manifest-only]      # update bookkeeping but skip the rebuild
#     [--strip-path-prefix /home]   # DGX
#   (build engine selection via the same env vars as build_inhouse_af.sh:
#    GLNEXUS_BIN/_SIF, BCFTOOLS_BIN/_SIF, REF/--ref …)
# =========================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

COHORT_DIR=""
NEW_SIZES=""
NEW_PATHS=""
THREADS=32
MEM_GB=192
MANIFEST_ONLY=0
STRIP_PREFIX=""
SEL_ARGS=()        # forwarded to select_cohort.py (exclusions, min-size…)
BUILD_ARGS=()      # forwarded to build_inhouse_af.sh (ref, config…)

while [ $# -gt 0 ]; do
  case "$1" in
    --cohort-dir)        COHORT_DIR="$2"; shift 2;;
    --new-sizes)         NEW_SIZES="$2"; shift 2;;
    --new-paths)         NEW_PATHS="$2"; shift 2;;
    --threads)           THREADS="$2"; shift 2;;
    --mem-gbytes)        MEM_GB="$2"; shift 2;;
    --manifest-only)     MANIFEST_ONLY=1; shift;;
    --strip-path-prefix) STRIP_PREFIX="$2"; shift 2;;
    --exclude-range)     SEL_ARGS+=(--exclude-range "$2"); shift 2;;
    --exclude-id-file)   SEL_ARGS+=(--exclude-id-file "$2"); shift 2;;
    --min-size-mb)       SEL_ARGS+=(--min-size-mb "$2"); shift 2;;
    --ref)               BUILD_ARGS+=(--ref "$2"); shift 2;;
    --config)            BUILD_ARGS+=(--config "$2"); shift 2;;
    -h|--help)           sed -n '2,55p' "$0"; exit 0;;
    *) echo "unknown arg: $1" >&2; exit 2;;
  esac
done

[ -n "$COHORT_DIR" ] || { echo "ERROR: --cohort-dir required" >&2; exit 2; }
[ -n "$NEW_SIZES$NEW_PATHS" ] || { echo "ERROR: --new-sizes or --new-paths required" >&2; exit 2; }
mkdir -p "$COHORT_DIR"

MANIFEST="$COHORT_DIR/cohort_manifest.tsv"
GVCFS="$COHORT_DIR/cohort_gvcfs.txt"
OUT="$COHORT_DIR/inhouse_af.hg38.vcf.gz"
LOG="$COHORT_DIR/updates.log"
TS="$(date '+%Y-%m-%d %H:%M:%S')"

# ---- 1. select_cohort on the new batch ----------------------------------
NEW_MANIFEST="$COHORT_DIR/.new_manifest.tsv"
if [ -n "$NEW_SIZES" ]; then SEL_INPUT=(--sizes "$NEW_SIZES"); else SEL_INPUT=(--paths "$NEW_PATHS"); fi
echo "[update] selecting new batch…"
python3 "$SCRIPT_DIR/select_cohort.py" "${SEL_INPUT[@]}" "${SEL_ARGS[@]}" \
  --out-manifest "$NEW_MANIFEST"

# ---- 2. append only samples not already in the cohort -------------------
BEFORE=0
[ -f "$MANIFEST" ] && BEFORE=$(awk -F'\t' 'NR>1 && $4=="include"' "$MANIFEST" | wc -l)

# header (first time only)
[ -f "$MANIFEST" ] || head -1 "$NEW_MANIFEST" > "$MANIFEST"

ADD_ROWS="$COHORT_DIR/.to_add.tsv"
awk -F'\t' '
  FNR==NR { if(FNR>1 && $4=="include") have[$1]=1; next }   # existing manifest
  FNR==1  { next }                                          # skip new header
  $4=="include" && !($1 in have) { print }                  # new & unseen
' "$MANIFEST" "$NEW_MANIFEST" > "$ADD_ROWS"

ADDED=$(wc -l < "$ADD_ROWS")
if [ "$ADDED" -gt 0 ]; then
  cat "$ADD_ROWS" >> "$MANIFEST"
  cut -f5 "$ADD_ROWS" >> "$GVCFS"            # col 5 = gvcf_path
fi
AFTER=$(awk -F'\t' 'NR>1 && $4=="include"' "$MANIFEST" | wc -l)
rm -f "$NEW_MANIFEST" "$ADD_ROWS"

echo "[update] cohort: $BEFORE → $AFTER  (+$ADDED new samples)"
if [ "$ADDED" -eq 0 ]; then
  echo "[update] nothing new to add — AF DB unchanged."
  echo "$TS	added=0	total=$AFTER	(no rebuild)" >> "$LOG"
  exit 0
fi

if [ "$MANIFEST_ONLY" -eq 1 ]; then
  echo "[update] --manifest-only: manifest + gvcf list updated, skipping rebuild."
  echo "$TS	added=$ADDED	total=$AFTER	(manifest-only)" >> "$LOG"
  exit 0
fi

# ---- 3. full re-genotype → temp output ----------------------------------
NEW_OUT="$OUT.new.vcf.gz"
echo "[update] rebuilding AF DB over $AFTER samples…"
BUILD_CMD=("$SCRIPT_DIR/build_inhouse_af.sh"
  --list "$GVCFS" --out "$NEW_OUT"
  --threads "$THREADS" --mem-gbytes "$MEM_GB"
  "${BUILD_ARGS[@]}")
[ -n "$STRIP_PREFIX" ] && BUILD_CMD+=(--strip-path-prefix "$STRIP_PREFIX")
"${BUILD_CMD[@]}"

# ---- 4. atomic swap -----------------------------------------------------
[ -f "$NEW_OUT" ] && [ -f "$NEW_OUT.tbi" ] || { echo "ERROR: build did not produce $NEW_OUT(.tbi)" >&2; exit 4; }
mv -f "$NEW_OUT"     "$OUT"
mv -f "$NEW_OUT.tbi" "$OUT.tbi"
echo "$TS	added=$ADDED	total=$AFTER	rebuilt=$OUT" >> "$LOG"

echo
echo "[update] DONE. AF DB now covers $AFTER samples → $OUT"
echo "[update] log: $LOG"
