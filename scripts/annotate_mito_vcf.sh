#!/usr/bin/env bash
# =========================================================
# annotate_mito_vcf.sh — annotate a mitochondrial VCF end-to-end
# =========================================================
# Input : a GATK Mutect2 --mitochondria-mode VCF (chrM / rCRS coords,
#         FILTER applied by FilterMutectCalls). FORMAT carries AF
#         (heteroplasmy fraction), AD, DP.
# Steps : 1. bcftools norm (left-align indels, split multiallelics)
#         2. VEP (gene / consequence / HGVS; vertebrate-mito codon table)
#         3. parse_mito_vep.py — join the local MITOMAP tables (disease /
#            status / MitoTIP / GenBank freq) → mito.annotated.tsv
# Output: <outdir>/mito.annotated.tsv       ← the file the UI reads;
#                                              point --outdir at the
#                                              sample's tertiary dir.
#         <outdir>/<sample>.mito.vep.vcf.gz (+ .tbi)  ← intermediate
#
# Paths are env-overridable; defaults follow the tertiary pipeline's
# `dgm` (production) profile. If $VEP_SIF / $BCFTOOLS_SIF are set the
# tools run inside those apptainer images; otherwise `vep` / `bcftools`
# are expected on PATH. python3 is always run on the host.
#
# Usage:
#   scripts/annotate_mito_vcf.sh --in 26WE0040.mito.vcf.gz [--sample 26WE0040] [--outdir .]
# =========================================================
set -euo pipefail

# ---- config (override via env) ----
# Defaults follow the tertiary pipeline's `dgm` (production) profile;
# the `local` profile uses /scratch/.../hg38 + /data/pylin1991/nf-containers.
REF_DIR="${REF_DIR:-/home/pipeline/reference/hg38}"
REF_FASTA="${REF_FASTA:-${REF_DIR}/Homo_sapiens_assembly38.fasta}"
VEP_CACHE="${VEP_CACHE:-${REF_DIR}/tertiary/vep_cache}"
VEP_SIF="${VEP_SIF:-/home/pipeline/nextflow_containers/vep_115.sif}"
BCFTOOLS_SIF="${BCFTOOLS_SIF:-/home/pipeline/nextflow_containers/bcftools_1.23.1.sif}"
MITOMAP_DIR="${MITOMAP_DIR:-${REF_DIR}/tertiary/mitomap}"
APPTAINER_BIND="${APPTAINER_BIND:---bind /home}"
THREADS="${THREADS:-4}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ---- args ----
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
if [ -z "$SAMPLE" ]; then
  SAMPLE="$(basename "$IN")"; SAMPLE="${SAMPLE%%.*}"
fi
mkdir -p "$OUTDIR"

# ---- tool wrappers (bare or apptainer) ----
run_vep() {      # run_vep <args...>
  if [ -n "$VEP_SIF" ]; then apptainer exec $APPTAINER_BIND "$VEP_SIF" vep "$@"
  else vep "$@"; fi
}
run_bcftools() { # run_bcftools <args...>
  if [ -n "$BCFTOOLS_SIF" ]; then apptainer exec $APPTAINER_BIND "$BCFTOOLS_SIF" bcftools "$@"
  else bcftools "$@"; fi
}
run_tabix() {    # tabix lives in the bcftools image too
  if [ -n "$BCFTOOLS_SIF" ]; then apptainer exec $APPTAINER_BIND "$BCFTOOLS_SIF" tabix "$@"
  else tabix "$@"; fi
}

NORM_VCF="${OUTDIR}/${SAMPLE}.mito.norm.vcf.gz"
VEP_VCF="${OUTDIR}/${SAMPLE}.mito.vep.vcf.gz"

# ---- Step 1: normalise (left-align indels, split any multiallelics) ----
# Mutect2-mito is run with --max-alt-allele-count 1 so multiallelics
# are rare, but left-alignment against the reference is still good
# hygiene before VEP. FILTER values are preserved — downstream code
# decides whether to keep e.g. weak_evidence / possible_numt calls.
echo "[mito] $SAMPLE  norm…" >&2
run_bcftools norm -m -any -f "$REF_FASTA" "$IN" -Oz -o "$NORM_VCF"
run_tabix -p vcf "$NORM_VCF"

# ---- Step 2: VEP ----
# No plugins (dbNSFP / LOFTEE / Pangolin don't apply to mtDNA).
# --hgvs / --hgvsg give HGVSc / HGVSp / genomic HGVS; for chrM VEP
# emits the m. notation. --mane / --flag_pick are mostly no-ops here
# (mito genes are single-transcript) but harmless and keep the CSQ
# format aligned with the nuclear pipeline.
echo "[mito] $SAMPLE  VEP…" >&2
run_vep \
  --input_file "$NORM_VCF" \
  --output_file "$VEP_VCF" \
  --vcf --compress_output bgzip \
  --offline --cache --dir_cache "$VEP_CACHE" \
  --species homo_sapiens --assembly GRCh38 \
  --fasta "$REF_FASTA" \
  --fork "$THREADS" \
  --hgvs --hgvsg \
  --symbol --canonical --mane --numbers --biotype --domains --uniprot \
  --flag_pick \
  --pick_order mane_select,mane_plus_clinical,canonical,appris,tsl,biotype,ccds,rank,length \
  --force_overwrite --no_stats --safe

run_tabix -p vcf "$VEP_VCF"

# ---- Step 3: parse to TSV (join the local MITOMAP tables) ----
# Final TSV gets the canonical, sample-free name so it can be dropped
# straight into tertiary_output/<LIS_ID>/mito.annotated.tsv where the
# UI looks for it (mirrors snv_indel.annotated.tsv). The intermediate
# VEP / norm VCFs keep the <sample>. prefix.
TSV="${OUTDIR}/mito.annotated.tsv"
MM_CC="${MITOMAP_DIR}/mitomap_mutations_coding_control.tsv"
MM_RNA="${MITOMAP_DIR}/mitomap_mutations_rna.tsv"
echo "[mito] $SAMPLE  parse → TSV…" >&2
python3 "${SCRIPT_DIR}/parse_mito_vep.py" \
  --vep_vcf "$VEP_VCF" \
  --mitomap_cc  "$MM_CC" \
  --mitomap_rna "$MM_RNA" \
  --sample_id "$SAMPLE" \
  --output "$TSV"

echo "[mito] $SAMPLE  done → $TSV  (+ $VEP_VCF)" >&2
run_bcftools stats "$VEP_VCF" | grep "^SN" >&2 || true
