#!/usr/bin/env bash
# =========================================================
# make_glnexus_list.sh — cohort list -> GLnexus-safe list (drop known-bad)
# =========================================================
# The incremental cohort (cohort_gvcfs.txt) intentionally KEEPS every sample,
# including ones with a malformed record that GLnexus refuses to load. This
# helper writes a filtered copy for the GLnexus audit baseline, excluding the
# ids in known_bad_gvcfs.txt (next to this script, unless --known-bad given).
#
# Usage:
#   scripts/inhouse_af/make_glnexus_list.sh \
#     --in  $DB/cohort_gvcfs.txt \
#     --out $DB/cohort_gvcfs.glnexus.txt \
#     [--known-bad <file>]
# =========================================================
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
KNOWN_BAD="$SCRIPT_DIR/known_bad_gvcfs.txt"
IN=""; OUT=""
while [ $# -gt 0 ]; do
  case "$1" in
    --in)        IN="$2"; shift 2;;
    --out)       OUT="$2"; shift 2;;
    --known-bad) KNOWN_BAD="$2"; shift 2;;
    -h|--help)   sed -n '2,18p' "$0"; exit 0;;
    *) echo "unknown arg: $1" >&2; exit 2;;
  esac
done
[ -n "$IN" ]  && [ -f "$IN" ]  || { echo "ERROR: --in <cohort_gvcfs.txt> required" >&2; exit 2; }
[ -n "$OUT" ] || { echo "ERROR: --out required" >&2; exit 2; }

# collect bare sample ids (strip comments/blank) from the known-bad file
ids=""
if [ -f "$KNOWN_BAD" ]; then
  ids=$(sed 's/#.*//' "$KNOWN_BAD" | awk 'NF{print $1}')
fi
if [ -z "$ids" ]; then
  cp "$IN" "$OUT"
  echo "[glnexus-list] no known-bad ids; copied $(grep -c . "$OUT") paths -> $OUT"
  exit 0
fi

# drop any path whose gVCF basename is "<id>.hard-filtered.gvcf.gz"
pat=$(echo "$ids" | paste -sd'|' -)
grep -vE "/(${pat})\.hard-filtered\.gvcf\.gz$" "$IN" > "$OUT"
n_in=$(grep -c . "$IN"); n_out=$(grep -c . "$OUT")
echo "[glnexus-list] dropped $((n_in - n_out)) known-bad; $n_out -> $OUT"
echo "[glnexus-list] excluded ids: $(echo "$ids" | paste -sd' ' -)"
