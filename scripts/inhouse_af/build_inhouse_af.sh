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
#     [--bed     primary_contigs.bed]   # restrict to chr1-22,X,Y,M — skip the
#                                       # ~3340 decoy/alt/HLA contigs (much faster)
#     [--config  glnexus_config_dragen.yml] [--threads 16] [--mem-gbytes 96] \
#     [--scratch <dir>] [--max-samples N] [--keep-cohort] \
#     [--cohort-bcf <cohort.raw.bcf>]   # re-derive AF from an existing genotyped
#                                       # BCF, skipping GLnexus stage 1
#     [--keep-cohort]                   # keep cohort.raw.bcf for --cohort-bcf reuse
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
# DRAGEN gVCFs need revise_genotypes:false (the built-in `gatk` preset aborts
# with "couldn't find genotype likelihoods"). Default to our DRAGEN config file;
# pass --config <preset|file> to override.
DEFAULT_CONFIG="$SCRIPT_DIR/glnexus_config_dragen.yml"
CONFIG="$DEFAULT_CONFIG"
THREADS=16
MEM_GB=96
LIST=""
OUT=""
SCRATCH=""
MAX_SAMPLES=0
KEEP_COHORT=0
STRIP_PREFIX=""
BED=""
EXT_COHORT_BCF=""

while [ $# -gt 0 ]; do
  case "$1" in
    --list)              LIST="$2"; shift 2;;
    --out)               OUT="$2"; shift 2;;
    --ref)               REF="$2"; shift 2;;
    --bed)               BED="$2"; shift 2;;
    --cohort-bcf)        EXT_COHORT_BCF="$2"; shift 2;;
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

[ -n "$OUT" ]  || { echo "ERROR: --out <inhouse_af.hg38.vcf.gz> required" >&2; exit 2; }
if [ -n "$EXT_COHORT_BCF" ]; then
  # Re-derive AF from an existing genotyped cohort BCF (skip GLnexus stage 1).
  [ -f "$EXT_COHORT_BCF" ] || { echo "ERROR: --cohort-bcf not found: $EXT_COHORT_BCF" >&2; exit 2; }
else
  [ -n "$LIST" ] || { echo "ERROR: --list <cohort_gvcfs.txt> required (or --cohort-bcf)" >&2; exit 2; }
  [ -f "$LIST" ] || { echo "ERROR: --list not found: $LIST" >&2; exit 2; }
fi
[ -f "$REF" ]  || { echo "ERROR: --ref not found: $REF" >&2; exit 2; }
[ -f "${REF}.fai" ] || { echo "ERROR: ${REF}.fai not found (samtools faidx)" >&2; exit 2; }
[ -z "$BED" ] || [ -f "$BED" ] || { echo "ERROR: --bed not found: $BED" >&2; exit 2; }
# CONFIG may be a built-in preset name (e.g. gatk_unfiltered) or a YAML file.
# If it looks like a path (has a / or .yml) it must exist.
case "$CONFIG" in
  */*|*.yml|*.yaml)
    [ -f "$CONFIG" ] || { echo "ERROR: --config file not found: $CONFIG" >&2; exit 2; } ;;
esac

# ---- resolve engines: container image preferred, else host binary --------
# GLnexus is only needed for stage 1; skip the check when reusing --cohort-bcf.
GLNEXUS_MODE="skip"
if [ -z "$EXT_COHORT_BCF" ]; then
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
if [ -z "$EXT_COHORT_BCF" ]; then
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
else
  N_SAMPLES="?"
fi

echo "[inhouse-af] samples    : $N_SAMPLES"
echo "[inhouse-af] config     : $CONFIG"
echo "[inhouse-af] ref        : $REF"
echo "[inhouse-af] bed        : ${BED:-(none — all contigs)}"
echo "[inhouse-af] out        : $OUT"
echo "[inhouse-af] scratch    : $SCRATCH"
if [ "$GLNEXUS_MODE" = "skip" ]; then
  echo "[inhouse-af] glnexus    : skip (reuse --cohort-bcf $EXT_COHORT_BCF)"
else
  echo "[inhouse-af] glnexus    : $GLNEXUS_MODE ($([ "$GLNEXUS_MODE" = sif ] && echo "$GLNEXUS_SIF" || command -v "$GLNEXUS_BIN"))"
fi
echo "[inhouse-af] bcftools   : $BCF_MODE ($([ "$BCF_MODE" = sif ] && echo "$BCFTOOLS_SIF" || command -v "$BCFTOOLS_BIN"))"
echo "[inhouse-af] threads/mem: $THREADS / ${MEM_GB}G"
echo

# ---- Stage 1: joint genotyping → raw cohort BCF -------------------------
if [ -n "$EXT_COHORT_BCF" ]; then
  COHORT_BCF="$EXT_COHORT_BCF"
  KEEP_COHORT=1          # never delete a user-supplied BCF
  echo "[inhouse-af] stage 1: skipped — reusing $COHORT_BCF ($(du -h "$COHORT_BCF" | cut -f1))"
else
  # GLnexus needs its scratch DB dir to NOT pre-exist.
  GLDB="$SCRATCH/GLnexus.DB"
  rm -rf "$GLDB"
  COHORT_BCF="$SCRATCH/cohort.raw.bcf"
  echo "[inhouse-af] stage 1: GLnexus joint genotyping ($N_SAMPLES gVCFs)…"
  # Build the GLnexus invocation (container or host binary).
  if [ "$GLNEXUS_MODE" = "sif" ]; then
    GLNEXUS_CMD=(apptainer exec --bind "$APPTAINER_BIND","$SCRATCH" "$GLNEXUS_SIF" glnexus_cli)
  else
    GLNEXUS_CMD=("$GLNEXUS_BIN")
  fi
  GLNEXUS_CMD+=(--config "$CONFIG" --dir "$GLDB" --threads "$THREADS" --mem-gbytes "$MEM_GB" --list "$EFF_LIST")
  [ -n "$BED" ] && GLNEXUS_CMD+=(--bed "$BED")
  # Pipe GLnexus stdout through bcftools so the cohort BCF is a properly
  # finalized, BGZF-EOF-terminated file — GLnexus's raw stdout can omit the EOF
  # marker, which later reads report as "EOF marker absent / truncated".
  if [ "$BCF_MODE" = "sif" ]; then
    "${GLNEXUS_CMD[@]}" | apptainer exec --bind "$APPTAINER_BIND","$SCRATCH" "$BCFTOOLS_SIF" bcftools view -Ob -o "$COHORT_BCF"
  else
    "${GLNEXUS_CMD[@]}" | "$BCF" view -Ob -o "$COHORT_BCF"
  fi
  rm -rf "$GLDB"   # the DB is large; the BCF is the only thing we keep
  echo "[inhouse-af] stage 1 done: $COHORT_BCF ($(du -h "$COHORT_BCF" | cut -f1))"
fi
echo

# ---- Stage 2: normalize → PASS → AF tags → strip genotypes → rename -----
RENAME="$SCRATCH/rename_annots.txt"
cat > "$RENAME" <<'EOF'
INFO/AC	INFO/INHOUSE_AC
INFO/AN	INFO/INHOUSE_AN
INFO/AF	INFO/INHOUSE_AF
INFO/AC_Hom	INFO/INHOUSE_AC_HOM
EOF

echo "[inhouse-af] stage 2: norm + PASS + fill-tags + sites-only → $OUT"
# --check-ref x drops the (few) sites whose REF disagrees with the FASTA instead
# of aborting; those can't be matched against the primary-assembly TSV anyway.
# norm's stderr (incl. the skipped-sites summary) is captured for reporting.
NORM_LOG="$SCRATCH/norm.log"
bcf_run "
  set -euo pipefail
  $BCF norm -m- -f '$REF' --check-ref x '$COHORT_BCF' -Ou 2> '$NORM_LOG' \
  | $BCF view -f 'PASS,.' -Ou \
  | $BCF +fill-tags -Ou -- -t AN,AC,AF,AC_Hom \
  | $BCF view -G -Ou \
  | $BCF annotate -x '^INFO/AC,INFO/AN,INFO/AF,INFO/AC_Hom' -Ou \
  | $BCF annotate --rename-annots '$RENAME' -Oz -o '$OUT'
  $BCF index -t -f '$OUT'
"
echo "[inhouse-af] norm summary (ref-mismatch sites are skipped, not fatal):"
grep -iE 'mismatch|skipped|total|lines' "$NORM_LOG" | tail -5 | sed 's/^/    /' || true

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
echo "    -c INFO/INHOUSE_AC,INFO/INHOUSE_AN,INFO/INHOUSE_AF,INFO/INHOUSE_AC_HOM <sample.vcf>"
echo "  (INHOUSE_AC_HOM = alleles in hom-alt genotypes = 2 x #homozygous individuals)"
