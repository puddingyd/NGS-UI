#!/usr/bin/env bash
# =========================================================
# annotate_mito_vcf.sh — mitochondrial VCF → mito.annotated.tsv
# =========================================================
# MITOMAP-only annotation (no VEP). Thin wrapper around
# scripts/parse_mito_vcf.py: reads the GATK Mutect2-mito VCF directly,
# derives gene/locus from a hard-coded rCRS map, joins the local
# MITOMAP tables (disease / status / MitoTIP / GenBank freq / amino-
# acid change). Pure Python — needs only the VCF + the MITOMAP TSVs.
#
# Input : a GATK Mutect2 --mitochondria-mode VCF (chrM / rCRS coords,
#         FILTER applied by FilterMutectCalls; FORMAT has AF=heteroplasmy,
#         AD, DP).
# Output: <outdir>/mito.annotated.tsv   ← the file the UI reads; point
#         --outdir at the sample's tertiary dir
#         (tertiary_output/<LIS_ID>/).
#
# $MITOMAP_DIR (default <ref_dir>/tertiary/mitomap, matching the
# tertiary pipeline) must hold mitomap_mutations_coding_control.tsv
# and mitomap_mutations_rna.tsv.
#
# Usage:
#   scripts/annotate_mito_vcf.sh --in 26WE0040.mito.vcf.gz [--sample 26WE0040] [--outdir .]
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
    --in)     IN="$2"; shift 2;;
    --sample) SAMPLE="$2"; shift 2;;
    --outdir) OUTDIR="$2"; shift 2;;
    -h|--help) sed -n '2,30p' "$0"; exit 0;;
    *) echo "unknown arg: $1" >&2; exit 2;;
  esac
done
[ -n "$IN" ] || { echo "ERROR: --in <vcf> required" >&2; exit 2; }
[ -f "$IN" ] || { echo "ERROR: input not found: $IN" >&2; exit 2; }
if [ -z "$SAMPLE" ]; then SAMPLE="$(basename "$IN")"; SAMPLE="${SAMPLE%%.*}"; fi
mkdir -p "$OUTDIR"

MM_CC="${MITOMAP_DIR}/mitomap_mutations_coding_control.tsv"
MM_RNA="${MITOMAP_DIR}/mitomap_mutations_rna.tsv"
TSV="${OUTDIR}/mito.annotated.tsv"

echo "[mito] $SAMPLE  parse → $TSV" >&2
python3 "${SCRIPT_DIR}/parse_mito_vcf.py" \
  --vcf "$IN" \
  --mitomap_cc  "$MM_CC" \
  --mitomap_rna "$MM_RNA" \
  --sample_id "$SAMPLE" \
  --output "$TSV"

# ClinVar annotation — fills CLINVAR_SIG / STARS / DN / SIGCONF on the
# TSV by joining (CHROM, POS, REF, ALT) against a chrM-only ClinVar
# VCF. Built once by scripts/build_mito_clinvar_vcf.sh; reused for
# every sample. Silently skipped when the file doesn't exist so this
# wrapper still works on hosts without ClinVar.
CLINVAR_MITO_VCF="${CLINVAR_MITO_VCF:-$HOME/NGS_UI/biotools/clinvar/clinvar_mito.vcf.gz}"
if [ -f "$CLINVAR_MITO_VCF" ]; then
  echo "[mito] $SAMPLE  ClinVar annotate ($CLINVAR_MITO_VCF)" >&2
  python3 "${SCRIPT_DIR}/annotate_clinvar.py" \
    --tsv "$TSV" \
    --clinvar "$CLINVAR_MITO_VCF"
else
  echo "[mito] $SAMPLE  ClinVar annotation skipped (no $CLINVAR_MITO_VCF)" >&2
fi

echo "[mito] $SAMPLE  done → $TSV" >&2
