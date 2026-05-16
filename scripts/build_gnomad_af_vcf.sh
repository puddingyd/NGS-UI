#!/usr/bin/env bash
# =========================================================
# build_gnomad_af_vcf.sh — ANNOVAR gnomAD txt → AF-only VCF
# =========================================================
# One-time wrapper: runs scripts/annovar_gnomad_to_vcf.py inside the
# pipeline's tertiary_python container (where pysam + bcftools live),
# producing a sorted + tabix-indexed sites VCF that stage_dragen_for_
# tertiary.sh can pass to `bcftools annotate -a` for pre-filtering.
#
# Usage:
#   scripts/build_gnomad_af_vcf.sh
#     [--txt /home/n102968/NGS_UI/biotools/hg38_gnomad41_genome.txt]
#     [--ref /home/pipeline/reference/hg38/Homo_sapiens_assembly38.fasta]
#     [--out $HOME/NGS_UI/biotools/gnomad/gnomad_af.hg38.vcf.gz]
#     [--snv-only]
#
# Conversion is single-threaded; gnomAD v4.1 genome (~100M rows)
# typically takes 30-90 min on a workstation.
# =========================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

PY_SIF="${PY_SIF:-/home/pipeline/nextflow_containers/tertiary_python_1.0.0.sif}"
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
    -h|--help)   sed -n '2,25p' "$0"; exit 0;;
    *) echo "unknown arg: $1" >&2; exit 2;;
  esac
done

[ -f "$TXT" ] || { echo "ERROR: --txt not found: $TXT" >&2; exit 2; }
[ -f "$REF" ] || { echo "ERROR: --ref not found: $REF" >&2; exit 2; }
[ -f "$PY_SIF" ] || { echo "ERROR: container not found: $PY_SIF" >&2; exit 2; }

mkdir -p "$(dirname "$OUT")"

echo "[build] txt : $TXT"
echo "[build] ref : $REF"
echo "[build] out : $OUT"
echo "[build] container: $PY_SIF"
echo "[build] extra args: ${EXTRA[*]:-(none)}"
echo

apptainer exec --bind /home,"$(dirname "$OUT")","$(dirname "$REF")","$(dirname "$TXT")" \
  "$PY_SIF" \
  python3 "$SCRIPT_DIR/annovar_gnomad_to_vcf.py" \
    --txt "$TXT" \
    --ref "$REF" \
    --out "$OUT" \
    "${EXTRA[@]}"

echo
echo "[build] verify:"
apptainer exec --bind /home "$PY_SIF" \
  bcftools view -h "$OUT" | tail -5
apptainer exec --bind /home "$PY_SIF" \
  bcftools view "$OUT" | head -3
ls -la "$OUT" "${OUT}.tbi"
