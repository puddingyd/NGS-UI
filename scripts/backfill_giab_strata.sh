#!/usr/bin/env bash
# =========================================================
# backfill_giab_strata.sh — add GIAB_STRATA to already-analysed samples
# =========================================================
# For samples that were analysed before the GIAB stratification step
# existed. For each tertiary_output/<SID>/snv_indel.annotated.tsv it:
#
#   1. annotate_giab_strata.py   — add/refresh the GIAB_STRATA column
#   2. build_snv_review_tsv.py   — rebuild the main-screen review TSV so
#                                  the new column reaches the UI
#   3. build_snv_gene_index.py   — REBUILD the gene index: step 1 rewrites
#                                  the whole TSV (atomic replace), shifting
#                                  every byte offset, so the old index
#                                  would seek stale bytes for gene search
#
# New analyses don't need this — run_stopgaps.sh already runs step 1 inline.
# Idempotent: safe to re-run. No-op for GIAB when the BED dir is missing.
#
# Usage:
#   scripts/backfill_giab_strata.sh                 # all samples under root
#   scripts/backfill_giab_strata.sh SID1 SID2 ...   # only these samples
#
# Env:
#   TERTIARY_OUTPUT_ROOT   sample root (default $HOME/NGS_UI/tertiary_output)
#   NGS_UI_GIAB_STRAT_DIR  BED dir (default $HOME/NGS_UI/biotools/giab_stratification)
# =========================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="${TERTIARY_OUTPUT_ROOT:-$HOME/NGS_UI/tertiary_output}"

[ -d "$ROOT" ] || { echo "ERROR: sample root not found: $ROOT" >&2; exit 1; }

# Build the sample list: explicit args, or every dir with an annotated TSV.
SAMPLES=()
if [ "$#" -gt 0 ]; then
  SAMPLES=("$@")
else
  for d in "$ROOT"/*/; do
    [ -f "$d/snv_indel.annotated.tsv" ] && SAMPLES+=("$(basename "$d")")
  done
fi

[ "${#SAMPLES[@]}" -gt 0 ] || { echo "no samples to process under $ROOT"; exit 0; }

echo "backfilling GIAB_STRATA for ${#SAMPLES[@]} sample(s) under $ROOT"
n_ok=0
n_skip=0
for sid in "${SAMPLES[@]}"; do
  tsv="$ROOT/$sid/snv_indel.annotated.tsv"
  if [ ! -f "$tsv" ]; then
    echo "  - $sid: no snv_indel.annotated.tsv, skip" >&2
    n_skip=$((n_skip + 1))
    continue
  fi
  echo "  • $sid"
  "$SCRIPT_DIR/annotate_giab_strata.py" --tsv "$tsv"
  "$SCRIPT_DIR/build_snv_review_tsv.py" --tsv "$tsv" >/dev/null
  "$SCRIPT_DIR/build_snv_gene_index.py" --tsv "$tsv" >/dev/null
  n_ok=$((n_ok + 1))
done

echo "done. $n_ok annotated, $n_skip skipped."
