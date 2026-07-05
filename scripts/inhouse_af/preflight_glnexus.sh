#!/usr/bin/env bash
# =========================================================
# preflight_glnexus.sh — find GLnexus-incompatible gVCFs BEFORE the long load
# =========================================================
# GLnexus aborts the whole joint-genotyping run (hours of bulk load) if any one
# input has a malformed record. This scans every gVCF in a cohort list IN
# PARALLEL (fast C decompression piped into the stdlib validator), so you learn
# which samples to drop up front instead of 11 h in.
#
# It writes:
#   <out-prefix>.clean.txt   gVCF paths that passed  (feed this to GLnexus)
#   <out-prefix>.bad.txt      "sample_id<TAB>gvcf<TAB>locus reason" for failures
# and, unless --no-append, appends any NEW bad sample ids to known_bad_gvcfs.txt
# next to this script (so the next quarterly rebuild skips them automatically).
#
# Usage:
#   scripts/inhouse_af/preflight_glnexus.sh \
#     --list  cohort_gvcfs.txt \
#     --out-prefix $DB/glnexus_preflight \
#     [--jobs 32] [--no-append]
#
# Decompression: bgzip -dc if on PATH (fastest), else zcat, else python gzip.
# Air-gap friendly: no bcftools needed.
# =========================================================
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VALIDATOR="$SCRIPT_DIR/validate_gvcf_glnexus.py"
KNOWN_BAD="$SCRIPT_DIR/known_bad_gvcfs.txt"

LIST=""; OUT_PREFIX=""; JOBS=32; APPEND=1
while [ $# -gt 0 ]; do
  case "$1" in
    --list)       LIST="$2"; shift 2;;
    --out-prefix) OUT_PREFIX="$2"; shift 2;;
    --jobs)       JOBS="$2"; shift 2;;
    --no-append)  APPEND=0; shift;;
    -h|--help)    sed -n '2,30p' "$0"; exit 0;;
    *) echo "unknown arg: $1" >&2; exit 2;;
  esac
done
[ -n "$LIST" ] && [ -f "$LIST" ] || { echo "ERROR: --list <cohort_gvcfs.txt> required" >&2; exit 2; }
[ -n "$OUT_PREFIX" ] || { echo "ERROR: --out-prefix required" >&2; exit 2; }

if command -v bgzip >/dev/null 2>&1; then DECOMP="bgzip -dc";
elif command -v zcat  >/dev/null 2>&1; then DECOMP="zcat";
else DECOMP=""; fi   # empty -> let the validator open the .gz itself (python gzip)

N=$(grep -cvE '^\s*(#|$)' "$LIST")
echo "[preflight] scanning $N gVCFs, jobs=$JOBS, decomp=${DECOMP:-python-gzip}"

BAD_DIR="$(mktemp -d)"
trap 'rm -rf "$BAD_DIR"' EXIT
export VALIDATOR DECOMP BAD_DIR

one() {
  gvcf="$1"
  id=$(basename "$gvcf" | sed 's/\..*//')
  if [ -n "$DECOMP" ]; then
    out=$($DECOMP "$gvcf" 2>/dev/null | python3 "$VALIDATOR" --stdin --sample-id "$id" 2>&1) && st=0 || st=$?
  else
    out=$(python3 "$VALIDATOR" --gvcf "$gvcf" --sample-id "$id" 2>&1) && st=0 || st=$?
  fi
  if [ "$st" -eq 2 ]; then
    # record every [BAD] line, tab-joined with the path
    echo "$out" | awk -v g="$gvcf" '/^\[BAD\]/{sub(/^\[BAD\] /,""); print $0"\t"g}' >> "$BAD_DIR/bad.txt"
    echo "BAD  $id"
  elif [ "$st" -ne 0 ]; then
    echo "ERR  $id (validator exit $st): $out" >&2
    echo -e "$id\t?\tvalidator error" >> "$BAD_DIR/bad.txt"
  else
    echo "$gvcf" >> "$BAD_DIR/clean.txt"
  fi
}
export -f one

grep -vE '^\s*(#|$)' "$LIST" | xargs -r -P "$JOBS" -I{} bash -c 'one "$@"' _ {}

# stable, sorted outputs
sort -u "$BAD_DIR/clean.txt" 2>/dev/null > "${OUT_PREFIX}.clean.txt" || : > "${OUT_PREFIX}.clean.txt"
sort -u "$BAD_DIR/bad.txt"   2>/dev/null > "${OUT_PREFIX}.bad.txt"   || : > "${OUT_PREFIX}.bad.txt"
n_clean=$(grep -c . "${OUT_PREFIX}.clean.txt" || true)
n_bad=$(cut -f1 "${OUT_PREFIX}.bad.txt" 2>/dev/null | sort -u | grep -c . || true)
echo "[preflight] clean=$n_clean  bad_samples=$n_bad"
echo "[preflight]   clean list -> ${OUT_PREFIX}.clean.txt   (feed to build_inhouse_af.sh)"
[ "$n_bad" -gt 0 ] && echo "[preflight]   bad detail -> ${OUT_PREFIX}.bad.txt"

if [ "$APPEND" -eq 1 ] && [ "$n_bad" -gt 0 ] && [ -f "$KNOWN_BAD" ]; then
  ts=$(date '+%Y-%m-%d')
  while read -r sid; do
    [ -n "$sid" ] || continue
    grep -qE "^\s*${sid}\b" "$KNOWN_BAD" && continue
    echo "${sid}        # preflight ${ts}: GLnexus GT-entry defect" >> "$KNOWN_BAD"
    echo "[preflight]   + recorded new known-bad: $sid"
  done < <(cut -f1 "${OUT_PREFIX}.bad.txt" | sort -u)
  echo "[preflight] known_bad_gvcfs.txt updated — commit it so the next rebuild skips them."
fi
