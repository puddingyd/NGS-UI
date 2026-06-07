#!/usr/bin/env bash
# =========================================================
# annotate_dragen_mito_vcf.sh — DRAGEN hard-filtered VCF → mito.annotated.tsv
# =========================================================
# Wrapper for scripts/parse_mito_vcf.py when mtDNA calls are embedded in
# a DRAGEN *.hard-filtered.vcf.gz. The parser scans the VCF, keeps only
# chrM/MT/chrMT rows, derives m. HGVS + rCRS gene/locus locally, and joins
# MITOMAP tables. FORMAT/AF is used as heteroplasmy; FORMAT/AD and DP are
# copied through; INFO/TLOD is used when present, otherwise FORMAT/SQ is
# used as the best available DRAGEN quality value in the TLOD output field.
#
# Usage on the dev machine:
#   scripts/annotate_dragen_mito_vcf.sh \
#     --in /home/datalake_Raw/.../vcf.gz/VAL-33.hard-filtered.vcf.gz \
#     --sample VAL-33 \
#     --outdir /home/n102968/NGS_UI/tertiary_output/VAL-33_legacy_mito
#
# Output:
#   <outdir>/mito.annotated.tsv
# =========================================================
set -euo pipefail

REF_DIR="${REF_DIR:-/home/pipeline/reference/hg38}"
MITOMAP_DIR="${MITOMAP_DIR:-${REF_DIR}/tertiary/mitomap}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

IN=""
SAMPLE=""
OUTDIR="."

while [ $# -gt 0 ]; do
  case "$1" in
    --in|--dragen-vcf) IN="$2"; shift 2;;
    --sample)          SAMPLE="$2"; shift 2;;
    --outdir)          OUTDIR="$2"; shift 2;;
    -h|--help)         sed -n '2,32p' "$0"; exit 0;;
    *) echo "unknown arg: $1" >&2; exit 2;;
  esac
done

[ -n "$IN" ] || { echo "ERROR: --in <VAL-33.hard-filtered.vcf.gz> required" >&2; exit 2; }
[ -f "$IN" ] || { echo "ERROR: input not found: $IN" >&2; exit 2; }

if [ -z "$SAMPLE" ]; then
  SAMPLE="$(basename "$IN")"
  SAMPLE="${SAMPLE%.hard-filtered.vcf.gz}"
  SAMPLE="${SAMPLE%%.*}"
fi

MM_CC="${MITOMAP_DIR}/mitomap_mutations_coding_control.tsv"
MM_RNA="${MITOMAP_DIR}/mitomap_mutations_rna.tsv"
[ -f "$MM_CC" ] || { echo "ERROR: MITOMAP coding/control table not found: $MM_CC" >&2; exit 2; }
[ -f "$MM_RNA" ] || { echo "ERROR: MITOMAP RNA table not found: $MM_RNA" >&2; exit 2; }

mkdir -p "$OUTDIR"
TSV="${OUTDIR}/mito.annotated.tsv"

echo "[dragen-mito] sample       : $SAMPLE" >&2
echo "[dragen-mito] input        : $IN" >&2
echo "[dragen-mito] mitomap_dir  : $MITOMAP_DIR" >&2
echo "[dragen-mito] output       : $TSV" >&2

python3 "${SCRIPT_DIR}/parse_mito_vcf.py" \
  --vcf "$IN" \
  --mitomap_cc "$MM_CC" \
  --mitomap_rna "$MM_RNA" \
  --sample_id "$SAMPLE" \
  --output "$TSV"

echo "[dragen-mito] done         : $TSV" >&2
