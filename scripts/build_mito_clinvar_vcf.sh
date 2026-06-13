#!/usr/bin/env bash
# =========================================================
# build_mito_clinvar_vcf.sh — extract chrM records from a full
# ClinVar VCF, normalising the chromosome name to chrM (UCSC
# convention) so it joins cleanly against our mito.annotated.tsv.
# =========================================================
# Why a separate small VCF: annotate_clinvar.py (the existing SNV
# post-processing) reads its whole input VCF into memory keyed by
# (CHROM, POS, REF, ALT). The full ClinVar GRCh38 file is ~1 GB and
# loading it just to annotate the ~200 chrM rows is wasteful — extract
# the mtDNA subset once, hand the small file to the mito annotation
# step.
#
# Usage:
#   scripts/build_mito_clinvar_vcf.sh \\
#       --in  /home/pipeline/reference/hg38/tertiary/clinvar/clinvar_20260510.vcf.gz \\
#       --out $HOME/NGS_UI/biotools/clinvar/clinvar_mito_20260510.vcf.gz
#
# The output VCF always uses `chrM` (rename MT→chrM if the input uses
# the Ensembl name). Indexed with tabix.
#
# Requires apptainer + the bcftools_1.23.1 container the rest of the
# scripts use.
# =========================================================
set -euo pipefail

BCF_SIF="${BCFTOOLS_SIF:-/home/pipeline/nextflow_containers/bcftools_1.23.1.sif}"
IN=""
OUT=""
while [ $# -gt 0 ]; do
  case "$1" in
    --in)  IN="$2"; shift 2;;
    --out) OUT="$2"; shift 2;;
    -h|--help) sed -n '2,25p' "$0"; exit 0;;
    *) echo "unknown arg: $1" >&2; exit 2;;
  esac
done
[ -n "$IN"  ] && [ -f "$IN"  ] || { echo "ERROR: --in not found: $IN"  >&2; exit 2; }
[ -n "$OUT" ] || { echo "ERROR: --out required" >&2; exit 2; }

mkdir -p "$(dirname "$OUT")"
tmp_dir=$(mktemp -d)
trap 'rm -rf "$tmp_dir"' EXIT

# Detect what name the input uses for mtDNA: chrM (UCSC) or MT
# (Ensembl). bcftools view -r returns exit 0 even on no hits, so we
# grep the resulting count.
echo "[mito-clinvar] detecting chrM contig name in $IN"
N_CHRM=$(apptainer exec --bind /home,"$tmp_dir" "$BCF_SIF" \
  bcftools view -r chrM "$IN" 2>/dev/null | grep -vc '^#' || true)
N_MT=$(apptainer exec --bind /home,"$tmp_dir" "$BCF_SIF" \
  bcftools view -r MT "$IN" 2>/dev/null | grep -vc '^#' || true)
echo "[mito-clinvar]   chrM rows: $N_CHRM"
echo "[mito-clinvar]   MT   rows: $N_MT"

REGION="chrM"
[ "$N_CHRM" -gt 0 ] || REGION="MT"

# rename-chrs map — used unconditionally; harmless if the input already
# uses chrM.
RENAME="$tmp_dir/rename.tsv"
printf "MT\tchrM\n" > "$RENAME"

echo "[mito-clinvar] extracting region $REGION → $OUT"
apptainer exec --bind /home,"$tmp_dir","$(dirname "$OUT")" "$BCF_SIF" bash -c "
  set -e
  bcftools view -r '$REGION' '$IN' \
    | bcftools annotate --rename-chrs '$RENAME' \
    -Oz -o '$OUT'
  bcftools index -t '$OUT'
"
N_OUT=$(apptainer exec --bind "$(dirname "$OUT")" "$BCF_SIF" \
  bcftools view "$OUT" | grep -vc '^#')
echo "[mito-clinvar] done — $N_OUT records → $OUT"
