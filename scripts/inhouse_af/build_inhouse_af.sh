#!/usr/bin/env bash
# =========================================================
# build_inhouse_af.sh — DRAGEN gVCFs → in-house AF sites VCF
# =========================================================
# Joint-genotype a cohort of DRAGEN per-sample gVCFs and emit a small,
# genotype-stripped "AF-only" sites VCF that `bcftools annotate -a` can
# merge into snv_indel.annotated.tsv (exactly like the gnomAD VCF built
# by scripts/build_gnomad_af_vcf.sh):
#
#     INFO/INHOUSE_AC   alt allele count across the cohort
#     INFO/INHOUSE_AN   total called alleles at the site (the denominator;
#                       comes from the gVCF reference blocks, so hom-ref
#                       and no-call samples are counted correctly)
#     INFO/INHOUSE_AF   INHOUSE_AC / INHOUSE_AN
#     INFO/INHOUSE_NHOM number of homozygous-alt samples
#
# WHY joint genotyping (not just stacking per-sample VCFs): a per-sample
# hard-filtered VCF only lists sites where THAT sample had a variant. It
# cannot tell "everyone else is hom-ref" from "everyone else was no-call",
# so naive stacking inflates rare-variant AF (denominator too small). The
# gVCF reference blocks give an accurate AN, which is exactly what we need
# for rare-variant AF as an annotation.
#
# ENGINE: GLnexus (open-source joint genotyper, container-friendly). DRAGEN
# gVCFs are GATK-style (<NON_REF>, GQ/DP/MIN_DP/PL), so the `gatk` config is
# the starting point — Phase 0 validates this is right for our data. The
# native DRAGEN *iterative* gVCF Genotyper is the true-incremental
# alternative (see scripts/inhouse_af/README.md); this script is the
# always-correct full-rebuild baseline.
#
# PASS-only: GLnexus replaces low-confidence genotypes with ./. (so AN
# already excludes no-calls); we additionally keep only site FILTER=PASS|'.'.
#
# Usage:
#   scripts/inhouse_af/build_inhouse_af.sh \
#     --list   cohort_gvcfs.txt \
#     --out    $HOME/NGS_UI/biotools/inhouse_af/inhouse_af.hg38.vcf.gz \
#     [--ref     /home/pipeline/reference/hg38/Homo_sapiens_assembly38.fasta] \
#     [--config  gatk] \
#     [--threads 16] [--mem-gbytes 96] \
#     [--scratch <dir>] [--max-samples N] [--keep-cohort]
#
# Containers (override via env):
#   GLNEXUS_SIF  apptainer image with glnexus_cli   (or GLNEXUS_BIN=host binary)
#   BCFTOOLS_SIF /home/pipeline/nextflow_containers/bcftools_1.23.1.sif
# =========================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

REF="${REF:-/home/pipeline/reference/hg38/Homo_sapiens_assembly38.fasta}"
BCF_SIF="${BCFTOOLS_SIF:-/home/pipeline/nextflow_containers/bcftools_1.23.1.sif}"
GLNEXUS_SIF="${GLNEXUS_SIF:-}"          # apptainer image w/ glnexus_cli
GLNEXUS_BIN="${GLNEXUS_BIN:-glnexus_cli}"  # host binary fallback
CONFIG="gatk"
THREADS=16
MEM_GB=96
LIST=""
OUT=""
SCRATCH=""
MAX_SAMPLES=0
KEEP_COHORT=0

while [ $# -gt 0 ]; do
  case "$1" in
    --list)        LIST="$2"; shift 2;;
    --out)         OUT="$2"; shift 2;;
    --ref)         REF="$2"; shift 2;;
    --config)      CONFIG="$2"; shift 2;;
    --threads)     THREADS="$2"; shift 2;;
    --mem-gbytes)  MEM_GB="$2"; shift 2;;
    --scratch)     SCRATCH="$2"; shift 2;;
    --max-samples) MAX_SAMPLES="$2"; shift 2;;
    --keep-cohort) KEEP_COHORT=1; shift;;
    -h|--help)     sed -n '2,60p' "$0"; exit 0;;
    *) echo "unknown arg: $1" >&2; exit 2;;
  esac
done

[ -n "$LIST" ] || { echo "ERROR: --list <cohort_gvcfs.txt> required" >&2; exit 2; }
[ -n "$OUT" ]  || { echo "ERROR: --out <inhouse_af.hg38.vcf.gz> required" >&2; exit 2; }
[ -f "$LIST" ] || { echo "ERROR: --list not found: $LIST" >&2; exit 2; }
[ -f "$REF" ]  || { echo "ERROR: --ref not found: $REF" >&2; exit 2; }
[ -f "${REF}.fai" ] || { echo "ERROR: ${REF}.fai not found (samtools faidx)" >&2; exit 2; }
[ -f "$BCF_SIF" ]   || { echo "ERROR: bcftools container not found: $BCF_SIF" >&2; exit 2; }

OUT_DIR="$(dirname "$OUT")"
mkdir -p "$OUT_DIR"
[ -n "$SCRATCH" ] || SCRATCH="$OUT_DIR/.glnexus_scratch"
mkdir -p "$SCRATCH"

# Build the effective gVCF list (drop blanks/comments, optional head for a
# quick Phase-0 smoke test).
EFF_LIST="$SCRATCH/cohort_gvcfs.effective.txt"
grep -vE '^\s*(#|$)' "$LIST" > "$EFF_LIST"
if [ "$MAX_SAMPLES" -gt 0 ]; then
  head -n "$MAX_SAMPLES" "$EFF_LIST" > "$EFF_LIST.tmp" && mv "$EFF_LIST.tmp" "$EFF_LIST"
fi
N_SAMPLES=$(wc -l < "$EFF_LIST")
[ "$N_SAMPLES" -gt 0 ] || { echo "ERROR: no gVCFs in list" >&2; exit 2; }

echo "[inhouse-af] samples    : $N_SAMPLES"
echo "[inhouse-af] config     : $CONFIG"
echo "[inhouse-af] ref        : $REF"
echo "[inhouse-af] out        : $OUT"
echo "[inhouse-af] scratch    : $SCRATCH"
echo "[inhouse-af] threads/mem: $THREADS / ${MEM_GB}G"
echo

# Resolve how to call GLnexus: container image, else host binary on PATH.
glnexus() {
  if [ -n "$GLNEXUS_SIF" ]; then
    apptainer exec --bind /home,"$SCRATCH" "$GLNEXUS_SIF" glnexus_cli "$@"
  elif command -v "$GLNEXUS_BIN" >/dev/null 2>&1; then
    "$GLNEXUS_BIN" "$@"
  else
    echo "ERROR: need GLnexus — set GLNEXUS_SIF=<image.sif> or put glnexus_cli on PATH (GLNEXUS_BIN)." >&2
    echo "       e.g. apptainer pull glnexus.sif docker://ghcr.io/dnanexus-rnd/glnexus:v1.4.1" >&2
    exit 3
  fi
}

# ---- Stage 1: joint genotyping → raw cohort BCF -------------------------
# GLnexus needs its scratch DB dir to NOT pre-exist.
GLDB="$SCRATCH/GLnexus.DB"
rm -rf "$GLDB"
COHORT_BCF="$SCRATCH/cohort.raw.bcf"
echo "[inhouse-af] stage 1: GLnexus joint genotyping ($N_SAMPLES gVCFs)…"
glnexus \
  --config "$CONFIG" \
  --dir "$GLDB" \
  --threads "$THREADS" \
  --mem-gbytes "$MEM_GB" \
  --list "$EFF_LIST" \
  > "$COHORT_BCF"
rm -rf "$GLDB"   # the DB is large; the BCF is the only thing we keep
echo "[inhouse-af] stage 1 done: $COHORT_BCF ($(du -h "$COHORT_BCF" | cut -f1))"
echo

# ---- Stage 2: normalize → PASS → AF tags → strip genotypes → rename -----
# rename map (created in scratch so it is bind-mounted into the container)
RENAME="$SCRATCH/rename_annots.txt"
cat > "$RENAME" <<'EOF'
INFO/AC	INFO/INHOUSE_AC
INFO/AN	INFO/INHOUSE_AN
INFO/AF	INFO/INHOUSE_AF
INFO/nhomalt	INFO/INHOUSE_NHOM
EOF

echo "[inhouse-af] stage 2: norm + PASS + fill-tags + sites-only → $OUT"
apptainer exec --bind /home,"$OUT_DIR","$SCRATCH" "$BCF_SIF" bash -c "
  set -euo pipefail
  bcftools norm -m- -f '$REF' '$COHORT_BCF' -Ou \
  | bcftools view -f 'PASS,.' -Ou \
  | bcftools +fill-tags -Ou -- -t AN,AC,AF,nhomalt \
  | bcftools view -G -Ou \
  | bcftools annotate -x '^INFO/AC,INFO/AN,INFO/AF,INFO/nhomalt' -Ou \
  | bcftools annotate --rename-annots '$RENAME' -Oz -o '$OUT'
  bcftools index -t -f '$OUT'
"

if [ "$KEEP_COHORT" -eq 0 ]; then
  rm -f "$COHORT_BCF"
fi

echo
echo "[inhouse-af] verify:"
apptainer exec --bind /home,"$OUT_DIR" "$BCF_SIF" bash -c "
  echo -n '  sites: '; bcftools index -n '$OUT'
  echo '  --- header INFO ---'
  bcftools view -h '$OUT' | grep -E '^##INFO=<ID=INHOUSE' || true
  echo '  --- first records ---'
  bcftools view '$OUT' | grep -vE '^#' | head -3 | cut -f1-8
  echo '  --- AN distribution (should peak near 2 x $N_SAMPLES) ---'
  bcftools query -f '%INFO/INHOUSE_AN\n' '$OUT' \
    | sort -n | awk '{a[NR]=\$1} END{if(NR){print \"    min=\"a[1], \"median=\"a[int(NR/2)], \"max=\"a[NR]}}'
"
ls -la "$OUT" "${OUT}.tbi"
echo
echo "[inhouse-af] DONE. Annotate per-sample TSVs with:"
echo "  bcftools annotate -a '$OUT' \\"
echo "    -c INFO/INHOUSE_AC,INFO/INHOUSE_AN,INFO/INHOUSE_AF,INFO/INHOUSE_NHOM <sample.vcf>"
