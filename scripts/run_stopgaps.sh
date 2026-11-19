#!/usr/bin/env bash
# =========================================================
# run_stopgaps.sh — one-sample post-processing (thin wrapper)
# =========================================================
# Runs the tertiary post-processing chain on a disposable working TSV. The
# immutable pipeline TSV remains under 03_acmg; only sparse/derived
# artifacts persist. All logic lives in scripts/postprocess_batch.py, which
# the tertiary worker also uses for whole batches (one DB pass per resource
# for every sample); this wrapper is the same code path with one sample.
#
#   1. clinvar-latest  annotate_clinvar_latest.py — weekly ClinVar comparison
#                      (UI only; the fixed pipeline ClinVar stays untouched)
#   2. genebe          annotate_acmg_genebe.py — SECOND-opinion ACMG to
#                      GENEBE_* (local DB → API cache → review-filtered API)
#   3. spliceai        annotate_spliceai.py — SpliceAI only (Research-only;
#                      dbNSFP stays in Nextflow)
#   3b. giab-strata    annotate_giab_strata.py
#   3b2. inhouse-af    annotate_inhouse_af.py
#   3c. mane-refseq    annotate_mane_refseq.py
#   4. litvar2         annotate_litvar2.py (review candidates only)
#   5. gpn-prefetch    GPN-MSA (+ Research-only GPN-Star V/M) review keys
#   6. snv-overlay → review-tsv → gene-index (sample steps)
#
# CNV/SV AnnotSV output is owned by the Nextflow pipeline (06_cnv_sv/); the
# worker refuses to post-process a sample whose pipeline CNV/SV TSVs are
# missing, so there is no local AnnotSV fallback here.
#
# Usage:
#   scripts/run_stopgaps.sh \
#       --work-tsv <post_dir>/.snv_indel.<job>.working.tsv \
#       --raw-tsv  <sample>/03_acmg/<SRC>.snv_indel.acmg.tsv \
#       --post-dir <sample>/08_postprocessing \
#       [--sample SID] [--seq-type WES|WGS] [--research-only]
#
# Resource flags (same names/defaults as before; env in parentheses):
#   --genebe-db (NGS_UI_GENEBE_DB) / --skip-genebe
#   --clinvar-latest-db (NGS_UI_CLINVAR_LATEST_DB) / --skip-clinvar-latest
#   --spliceai-snv / --spliceai-indel / --skip-spliceai (--skip-extra-vep)
#   --candidate-bed (NGS_UI_CDS_CANDIDATE_BED) / --skip-candidate-bed
#   --giab-strat-dir (NGS_UI_GIAB_STRAT_DIR) / --skip-giab
#   --inhouse-af-db (NGS_UI_INHOUSE_AF_DB) / --skip-inhouse-af
#   --mane-summary (NGS_UI_MANE_SUMMARY) / --skip-mane-refseq
#   --litvar2-db (NGS_UI_LITVAR2_DB) / --skip-litvar2
#   --gpn-msa-db (NGS_UI_GPN_MSA_DB) / --gpn-star-dir (NGS_UI_GPN_STAR_DIR)
#   NGS_UI_PYTHON_BIN — interpreter with the worker's dependencies (Polars)
#   NGS_UI_POSTPROC_WORKERS — parallel sample steps (batch use)
# =========================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [ "${1:-}" = "-h" ] || [ "${1:-}" = "--help" ]; then
  sed -n '2,46p' "$0"
  exit 0
fi
exec "${NGS_UI_PYTHON_BIN:-python3}" "$SCRIPT_DIR/postprocess_batch.py" "$@"
