#!/usr/bin/env bash
# =========================================================
# build_gnomad_af_vcf.sh — ANNOVAR gnomAD txt → AF-only VCF
# =========================================================
# Two-stage wrapper:
#   1. scripts/annovar_gnomad_to_vcf.py runs on HOST python3
#      (stdlib only — uses FASTA .fai index for indel padding, no
#      pysam needed) to produce an unsorted VCF.
#   2. Pipe through `bcftools sort` (in bcftools sif) → bgzipped +
#      tabix-indexed final VCF for `bcftools annotate -a`.
#
# Usage:
#   scripts/build_gnomad_af_vcf.sh
#     [--txt    $HOME/NGS_UI/biotools/hg38_gnomad41_genome.txt]
#     [--ref    /home/pipeline/reference/hg38/Homo_sapiens_assembly38.fasta]
#     [--out    $HOME/NGS_UI/biotools/gnomad/gnomad_af.hg38.vcf.gz]
#     [--af-col gnomad41_genome_AF]
#     [--min-af 0.01]
#     [--snv-only]
#
# Conversion is single-threaded; gnomAD v4.1 genome ANNOVAR txt
# (~870M rows) takes 30-60 min with --min-af 0.01, or hours without
# a cutoff.
# =========================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

BCF_SIF="${BCFTOOLS_SIF:-/home/pipeline/nextflow_containers/bcftools_1.23.1.sif}"
TXT="${TXT:-$HOME/NGS_UI/biotools/hg38_gnomad41_genome.txt}"
REF="${REF:-/home/pipeline/reference/hg38/Homo_sapiens_assembly38.fasta}"
OUT="${OUT:-$HOME/NGS_UI/biotools/gnomad/gnomad_af.hg38.vcf.gz}"
EXTRA=()
while [ $# -gt 0 ]; do
  case "$1" in
    --txt)       TXT="$2"; shift 2;;
    --ref)       REF="$2"; shift 2;;
    --out)       OUT="$2"; shift 2;;
    --snv-only)  EXTRA+=(--snv-only); shift;;
    --af-col)    EXTRA+=(--af-col "$2"); shift 2;;
    --min-af)    EXTRA+=(--min-af "$2"); shift 2;;
    -h|--help)   sed -n '2,30p' "$0"; exit 0;;
    *) echo "unknown arg: $1" >&2; exit 2;;
  esac
done

[ -f "$TXT" ] || { echo "ERROR: --txt not found: $TXT" >&2; exit 2; }
[ -f "$REF" ] || { echo "ERROR: --ref not found: $REF" >&2; exit 2; }
[ -f "${REF}.fai" ] || { echo "ERROR: ${REF}.fai not found" >&2; exit 2; }
[ -f "$BCF_SIF" ] || { echo "ERROR: container not found: $BCF_SIF" >&2; exit 2; }

mkdir -p "$(dirname "$OUT")"
UNSORTED="${OUT}.unsorted.vcf"

echo "[build] txt        : $TXT"
echo "[build] ref        : $REF"
echo "[build] out        : $OUT"
echo "[build] extra args : ${EXTRA[*]:-(none)}"
echo "[build] bcftools   : $BCF_SIF"
echo

# Stage 1 — convert on host (stdlib, no pysam)
echo "[build] stage 1: ANNOVAR txt → unsorted VCF (host python)"
python3 "$SCRIPT_DIR/annovar_gnomad_to_vcf.py" \
    --txt "$TXT" \
    --ref "$REF" \
    --out "$UNSORTED" \
    "${EXTRA[@]}"

echo
echo "[build] stage 2: bcftools sort + bgzip + tabix"
apptainer exec --bind /home,"$(dirname "$OUT")" "$BCF_SIF" bash -c "
  set -e
  bcftools sort --max-mem 4G '$UNSORTED' -Oz -o '$OUT'
  bcftools index -t -f '$OUT'
"
rm -f "$UNSORTED"

echo
echo "[build] verify:"
apptainer exec --bind /home "$BCF_SIF" bash -c "
  bcftools view -h '$OUT' | tail -5
  echo '---'
  bcftools view '$OUT' | head -3
"
ls -la "$OUT" "${OUT}.tbi"
