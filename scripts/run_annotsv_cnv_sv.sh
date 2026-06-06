#!/usr/bin/env bash
# =========================================================
# run_annotsv_cnv_sv.sh — NGS-UI CNV/SV AnnotSV entry point
# =========================================================
# Runs AnnotSV for the CNV/SV files that NGS-UI can discover from a
# tertiary job, then writes the GUI-expected files:
#
#   <out-dir>/cnv.annotated.tsv
#   <out-dir>/sv.annotated.tsv
#
# DRAGEN mode receives the hard-filtered SNV VCF path and looks for
# sibling <sample>.cnv.vcf.gz and <sample>.sv.vcf.gz in the same vcf.gz
# directory. In-house mode receives explicit gCNV and Delly VCF paths.
#
# AnnotSV install defaults:
#   ANNOTSV_BIN=$HOME/NGS_UI/biotools/AnnotSV/bin/AnnotSV
#   ANNOTSV_ANNOTATIONS=$HOME/NGS_UI/biotools/AnnotSV/share/AnnotSV
#
# Usage:
#   scripts/run_annotsv_cnv_sv.sh \
#     --sample 26WE0001 \
#     --out-dir $HOME/NGS_UI/tertiary_output/26WE0001 \
#     --dragen-cnv-source /path/to/26WE0001.hard-filtered.vcf.gz
#
#   scripts/run_annotsv_cnv_sv.sh \
#     --sample 26WE0001 \
#     --out-dir $HOME/NGS_UI/tertiary_output/26WE0001 \
#     --inhouse-cnv-vcf /path/to/26WE0001.gcnv.vcf.gz \
#     --inhouse-sv-vcf  /path/to/26WE0001.delly.vcf.gz
# =========================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

SID=""
OUT_DIR=""
DRAGEN_VCF=""
INHOUSE_CNV_VCF=""
INHOUSE_SV_VCF=""

while [ $# -gt 0 ]; do
  case "$1" in
    --sample)             SID="$2"; shift 2;;
    --out-dir)            OUT_DIR="$2"; shift 2;;
    --dragen-cnv-source)  DRAGEN_VCF="$2"; shift 2;;
    --inhouse-cnv-vcf)    INHOUSE_CNV_VCF="$2"; shift 2;;
    --inhouse-sv-vcf)     INHOUSE_SV_VCF="$2"; shift 2;;
    -h|--help)            sed -n '2,38p' "$0"; exit 0;;
    *) echo "unknown arg: $1" >&2; exit 2;;
  esac
done

[ -n "$SID" ] || { echo "ERROR: --sample required" >&2; exit 2; }
[ -n "$OUT_DIR" ] || { echo "ERROR: --out-dir required" >&2; exit 2; }

if [ -n "$DRAGEN_VCF" ] && { [ -n "$INHOUSE_CNV_VCF" ] || [ -n "$INHOUSE_SV_VCF" ]; }; then
  echo "ERROR: use either --dragen-cnv-source or --inhouse-*-vcf, not both" >&2
  exit 2
fi

if [ -n "$DRAGEN_VCF" ]; then
  if [ ! -f "$DRAGEN_VCF" ]; then
    echo "[annotsv] skipped (--dragen-cnv-source not found: $DRAGEN_VCF)"
    exit 0
  fi
  "$SCRIPT_DIR/annotate_dragen_cnv_annotsv.sh" \
    --dragen-vcf "$DRAGEN_VCF" \
    --sample "$SID" \
    --out-dir "$OUT_DIR"
elif [ -n "$INHOUSE_CNV_VCF" ] || [ -n "$INHOUSE_SV_VCF" ]; then
  "$SCRIPT_DIR/annotate_inhouse_cnv_sv_annotsv.sh" \
    --cnv-vcf "$INHOUSE_CNV_VCF" \
    --sv-vcf  "$INHOUSE_SV_VCF" \
    --sample  "$SID" \
    --out-dir "$OUT_DIR"
else
  echo "[annotsv] skipped (no --dragen-cnv-source / --inhouse-*-vcf)"
fi
