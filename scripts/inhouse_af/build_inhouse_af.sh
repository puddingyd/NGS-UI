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
# ENGINE: GLnexus (open-source joint genotyper). DRAGEN gVCFs are GATK-style
# (<NON_REF>, GQ/DP/MIN_DP/PL), so the `gatk` config is the starting point —
# Phase 0 validates this is right for our data. The native DRAGEN *iterative*
# gVCF Genotyper is the true-incremental alternative (see README.md); this
# script is the always-correct full-rebuild baseline.
#
# PORTABILITY (DGM → air-gapped DGX): every external tool is pluggable so the
# whole pipeline can run with ZERO containers — just two static binaries:
#   * GLnexus : GLNEXUS_SIF=<image.sif>  OR  GLNEXUS_BIN=<glnexus_cli> (static)
#   * bcftools: BCFTOOLS_SIF=<image.sif> OR  BCFTOOLS_BIN=<bcftools>
# When a *_SIF is set we apptainer-exec it; otherwise we run the host binary.
# On the DGX the datalake paths drop the leading /home — feed the same cohort
# list but add `--strip-path-prefix /home`, and pass `--ref` without /home.
#
# PASS-only: GLnexus replaces low-confidence genotypes with ./. (so AN
# already excludes no-calls); we additionally keep only site FILTER=PASS|'.'.
#
# Usage:
#   scripts/inhouse_af/build_inhouse_af.sh \
#     --list   cohort_gvcfs.txt \
#     --out    $NGS_UI_HOME/biotools/inhouse_af/inhouse_af.hg38.vcf.gz \
#     [--ref     /home/datalake_Intermediate/pipeline/reference/hg38/Homo_sapiens_assembly38.fasta] \
#     [--config  gatk] [--threads 16] [--mem-gbytes 96] \
#     [--scratch <dir>] [--max-samples N] [--keep-cohort] \
#     [--strip-path-prefix /home]    # DGX: cohort list paths drop /home
#
# Containers / binaries (override via env):
#   GLNEXUS_SIF / GLNEXUS_BIN     joint genotyper (BIN default: glnexus_cli)
#   BCFTOOLS_SIF / BCFTOOLS_BIN   bcftools        (BIN default: bcftools)
#   APPTAINER_BIND                bind spec when using *_SIF (default: /home)
# =========================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

REF="${REF:-/home/datalake_Intermediate/pipeline/reference/hg38/Homo_sapiens_assembly38.fasta}"
GLNEXUS_SIF="${GLNEXUS_SIF:-}"
GLNEXUS_BIN="${GLNEXUS_BIN:-glnexus_cli}"
BCFTOOLS_SIF="${BCFTOOLS_SIF:-}"
BCFTOOLS_BIN="${BCFTOOLS_BIN:-bcftools}"
APPTAINER_BIND="${APPTAINER_BIND:-/home}"
CONFIG="gatk"
THREADS=16
MEM_GB=96
LIST=""
OUT=""
SCRATCH=""
MAX_SAMPLES=0
KEEP_COHORT=0
STRIP_PREFIX=""

while [ $# -gt 0 ]; do
  case "$1" in
    --list)              LIST="$2"; shift 2;;
    --out)               OUT="$2"; shift 2;;
    --ref)               REF="$2"; shift 2;;
    --config)            CONFIG="$2"; shift 2;;
    --threads)           THREADS="$2"; shift 2;;
    --mem-gbytes)        MEM_GB="$2"; shift 2;;
    --scratch)           SCRATCH="$2"; shift 2;;
    --max-samples)       MAX_SAMPLES="$2"; shift 2;;
    --keep-cohort)       KEEP_COHORT=1; shift;;
    --strip-path-prefix) STRIP_PREFIX="$2"; shift 2;;
    -h|--help)           sed -n '2,60p' "$0"; exit 0;;
    *) echo "unknown arg: $1" >&2; exit 2;;
  esac
done

[ -n "$LIST" ] || { echo "ERROR: --list <cohort_gvcfs.txt> required" >&2; exit 2; }
[ -n "$OUT" ]  || { echo "ERROR: --out <inhouse_af.hg38.vcf.gz> required" >&2; exit 2; }
[ -f "$LIST" ] || { echo "ERROR: --list not found: $LIST" >&2; exit 2; }
[ -f "$REF" ]  || { echo "ERROR: --ref not found: $REF" >&2; exit 2; }
[ -f "${REF}.fai" ] || { echo "ERROR: ${REF}.fai not found (samtools faidx)" >&2; exit 2; }

# ---- resolve engines: container image preferred, else host binary --------
if [ -n "$GLNEXUS_SIF" ]; then
  [ -f "$GLNEXUS_SIF" ] || { echo "ERROR: GLNEXUS_SIF not found: $GLNEXUS_SIF" >&2; exit 2; }
  GLNEXUS_MODE="sif"
elif command -v "$GLNEXUS_BIN" >/dev/null 2>&1; then
  GLNEXUS_MODE="bin"
else
  echo "ERROR: need GLnexus — set GLNEXUS_SIF=<image.sif> or put glnexus_cli on PATH (GLNEXUS_BIN)." >&2
  echo "       static binary:  https://github.com/dnanexus-rnd/GLnexus/releases (glnexus_cli)" >&2
  echo "       container:      apptainer pull glnexus.sif docker://ghcr.io/dnanexus-rnd/glnexus:v1.4.1" >&2
  exit 3
fi
if [ -n "$BCFTOOLS_SIF" ]; then
  [ -f "$BCFTOOLS_SIF" ] || { echo "ERROR: BCFTOOLS_SIF not found: $BCFTOOLS_SIF" >&2; exit 2; }
  BCF_MODE="sif"; BCF="bcftools"
elif command -v "$BCFTOOLS_BIN" >/dev/null 2>&1; then
  BCF_MODE="bin"; BCF="$BCFTOOLS_BIN"
else
  echo "ERROR: need bcftools — set BCFTOOLS_SIF=<image.sif> or put bcftools on PATH (BCFTOOLS_BIN)." >&2
  exit 3
fi

OUT_DIR="$(dirname "$OUT")"
mkdir -p "$OUT_DIR"
[ -n "$SCRATCH" ] || SCRATCH="$OUT_DIR/.glnexus_scratch"
mkdir -p "$SCRATCH"

# run a bcftools snippet in the chosen mode (container exec or host shell)
bcf_run() {  # $1 = bash snippet referencing $BCF
  if [ "$BCF_MODE" = "sif" ]; then
    apptainer exec --bind "$APPTAINER_BIND","$OUT_DIR","$SCRATCH" "$BCFTOOLS_SIF" bash -c "$1"
  else
    bash -c "$1"
  fi
}

# ---- effective gVCF list: drop blanks/comments, optional prefix strip ----
EFF_LIST="$SCRATCH/cohort_gvcfs.effective.txt"
grep -vE '^\s*(#|$)' "$LIST" > "$EFF_LIST"
if [ -n "$STRIP_PREFIX" ]; then
  sed -i "s#^${STRIP_PREFIX}##" "$EFF_LIST"
fi
if [ "$MAX_SAMPLES" -gt 0 ]; then
  head -n "$MAX_SAMPLES" "$EFF_LIST" > "$EFF_LIST.tmp" && mv "$EFF_LIST.tmp" "$EFF_LIST"
fi
N_SAMPLES=$(wc -l < "$EFF_LIST")
[ "$N_SAMPLES" -gt 0 ] || { echo "ERROR: no gVCFs in list" >&2; exit 2; }
# sanity: first gVCF must be readable from here
FIRST=$(head -1 "$EFF_LIST")
[ -f "$FIRST" ] || { echo "ERROR: first gVCF not readable: $FIRST  (wrong --strip-path-prefix?)" >&2; exit 2; }

echo "[inhouse-af] samples    : $N_SAMPLES"
echo "[inhouse-af] config     : $CONFIG"
echo "[inhouse-af] ref        : $REF"
echo "[inhouse-af] out        : $OUT"
echo "[inhouse-af] scratch    : $SCRATCH"
echo "[inhouse-af] glnexus    : $GLNEXUS_MODE ($([ "$GLNEXUS_MODE" = sif ] && echo "$GLNEXUS_SIF" || command -v "$GLNEXUS_BIN"))"
echo "[inhouse-af] bcftools   : $BCF_MODE ($([ "$BCF_MODE" = sif ] && echo "$BCFTOOLS_SIF" || command -v "$BCFTOOLS_BIN"))"
echo "[inhouse-af] threads/mem: $THREADS / ${MEM_GB}G"
echo

# ---- Stage 1: joint genotyping → raw cohort BCF -------------------------
# GLnexus needs its scratch DB dir to NOT pre-exist.
GLDB="$SCRATCH/GLnexus.DB"
rm -rf "$GLDB"
COHORT_BCF="$SCRATCH/cohort.raw.bcf"
echo "[inhouse-af] stage 1: GLnexus joint genotyping ($N_SAMPLES gVCFs)…"
if [ "$GLNEXUS_MODE" = "sif" ]; then
  apptainer exec --bind "$APPTAINER_BIND","$SCRATCH" "$GLNEXUS_SIF" glnexus_cli \
    --config "$CONFIG" --dir "$GLDB" --threads "$THREADS" --mem-gbytes "$MEM_GB" \
    --list "$EFF_LIST" > "$COHORT_BCF"
else
  "$GLNEXUS_BIN" \
    --config "$CONFIG" --dir "$GLDB" --threads "$THREADS" --mem-gbytes "$MEM_GB" \
    --list "$EFF_LIST" > "$COHORT_BCF"
fi
rm -rf "$GLDB"   # the DB is large; the BCF is the only thing we keep
echo "[inhouse-af] stage 1 done: $COHORT_BCF ($(du -h "$COHORT_BCF" | cut -f1))"
echo

# ---- Stage 2: normalize → PASS → AF tags → strip genotypes → rename -----
RENAME="$SCRATCH/rename_annots.txt"
cat > "$RENAME" <<'EOF'
INFO/AC	INFO/INHOUSE_AC
INFO/AN	INFO/INHOUSE_AN
INFO/AF	INFO/INHOUSE_AF
INFO/nhomalt	INFO/INHOUSE_NHOM
EOF

echo "[inhouse-af] stage 2: norm + PASS + fill-tags + sites-only → $OUT"
bcf_run "
  set -euo pipefail
  $BCF norm -m- -f '$REF' '$COHORT_BCF' -Ou \
  | $BCF view -f 'PASS,.' -Ou \
  | $BCF +fill-tags -Ou -- -t AN,AC,AF,nhomalt \
  | $BCF view -G -Ou \
  | $BCF annotate -x '^INFO/AC,INFO/AN,INFO/AF,INFO/nhomalt' -Ou \
  | $BCF annotate --rename-annots '$RENAME' -Oz -o '$OUT'
  $BCF index -t -f '$OUT'
"

[ "$KEEP_COHORT" -eq 1 ] || rm -f "$COHORT_BCF"

echo
echo "[inhouse-af] verify:"
bcf_run "
  echo -n '  sites: '; $BCF index -n '$OUT'
  echo '  --- header INFO ---'
  $BCF view -h '$OUT' | grep -E '^##INFO=<ID=INHOUSE' || true
  echo '  --- first records ---'
  $BCF view '$OUT' | grep -vE '^#' | head -3 | cut -f1-8
  echo '  --- AN distribution (should peak near 2 x $N_SAMPLES) ---'
  $BCF query -f '%INFO/INHOUSE_AN\n' '$OUT' \
    | sort -n | awk '{a[NR]=\$1} END{if(NR){print \"    min=\"a[1], \"median=\"a[int(NR/2)], \"max=\"a[NR]}}'
"
ls -la "$OUT" "${OUT}.tbi"
echo
echo "[inhouse-af] DONE. Annotate per-sample TSVs with:"
echo "  bcftools annotate -a '$OUT' \\"
echo "    -c INFO/INHOUSE_AC,INFO/INHOUSE_AN,INFO/INHOUSE_AF,INFO/INHOUSE_NHOM <sample.vcf>"
