#!/usr/bin/env bash
# =========================================================
# backfill_inhouse_af.sh — add INHOUSE_AC/AN/AF to already-analysed samples
# =========================================================
# For samples analysed before the in-house AF step existed, or after the
# in-house AF DB is refreshed with a new batch. For each
# tertiary_output/<SID>/snv_indel.annotated.tsv it:
#
#   1. annotate_inhouse_af.py    — add/refresh INHOUSE_AC/AN/AF columns
#   2. build_snv_review_tsv.py   — rebuild the main-screen review TSV so the
#                                  new columns reach the UI
#   3. build_snv_gene_index.py   — REBUILD the gene index: step 1 rewrites the
#                                  whole TSV (atomic replace), shifting every
#                                  byte offset, so the old index would seek
#                                  stale bytes for gene search
#
# New analyses don't need this — run_stopgaps.sh runs step 1 inline.
# Idempotent; no-op for in-house AF when the DB is missing.
#
# Usage:
#   scripts/backfill_inhouse_af.sh                 # all samples under root
#   scripts/backfill_inhouse_af.sh SID1 SID2 ...   # only these samples
#
# Env:
#   TERTIARY_OUTPUT_ROOT   sample root (default $HOME/NGS_UI/tertiary_output)
#   NGS_UI_INHOUSE_AF_DB   sites VCF (default $HOME/NGS_UI/biotools/inhouse_af/
#                          inhouse_af.hg38.vcf.gz)
# =========================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="${TERTIARY_OUTPUT_ROOT:-$HOME/NGS_UI/tertiary_output}"
DB="${NGS_UI_INHOUSE_AF_DB:-$HOME/NGS_UI/biotools/inhouse_af/inhouse_af.hg38.vcf.gz}"

[ -d "$ROOT" ] || { echo "ERROR: sample root not found: $ROOT" >&2; exit 1; }
if [ ! -f "$DB" ]; then
  echo "ERROR: in-house AF DB not found: $DB" >&2
  echo "       deploy it first with scripts/inhouse_af/deploy_inhouse_af_db.sh" >&2
  exit 1
fi

SAMPLES=()
if [ "$#" -gt 0 ]; then
  SAMPLES=("$@")
else
  for d in "$ROOT"/*/; do
    [ -f "$d/snv_indel.annotated.tsv" ] && SAMPLES+=("$(basename "$d")")
  done
fi
[ "${#SAMPLES[@]}" -gt 0 ] || { echo "no samples to process under $ROOT"; exit 0; }

echo "backfilling INHOUSE_AF for ${#SAMPLES[@]} sample(s) under $ROOT"
echo "  DB: $DB"
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
  "$SCRIPT_DIR/annotate_inhouse_af.py" --tsv "$tsv" --db "$DB"
  "$SCRIPT_DIR/build_snv_review_tsv.py" --tsv "$tsv" >/dev/null
  "$SCRIPT_DIR/build_snv_gene_index.py" --tsv "$tsv" >/dev/null
  n_ok=$((n_ok + 1))
done

echo "done. $n_ok annotated, $n_skip skipped."
