#!/usr/bin/env bash
# =========================================================
# stage_dragen_for_tertiary.sh
# =========================================================
# Stage a single-sample DRAGEN germline hard-filter VCF or an in-house
# 2-sample ensemble VCF so the tertiary pipeline (which expects
# ensemble DV+HC input) accepts it. Does:
#
#   1. Drop chrM rows (mtDNA is handled separately by
#      scripts/annotate_mito_vcf.sh — DRAGEN's chrM calls have
#      heteroplasmy in FORMAT/AF, no TLOD; the mito script's
#      auto-detect handles that).
#   2. bcftools norm — left-align indels, split multi-allelics. DRAGEN
#      raw VCFs ship un-normalised; downstream tools (Pangolin
#      especially) segfault on the corner cases (IUPAC bases, un-split
#      multi-allelics, etc.).
#   3. Reheader sample columns to ${SID}_DV (+ ${SID}_HC):
#      - 1-sample input (DRAGEN): rename the one column to ${SID}_DV,
#        then append a synthetic ${SID}_HC column populated with "./."
#        for every row so PREPARE_VCF:ADD_CALLERS_TAG sees the expected
#        (DV, HC) shape. CALLERS = DV everywhere.
#      - 2-sample input (in-house ensemble): both columns already
#        carry the (DV, HC) shape; just rename to the new SID prefix.
#   5. bgzip + tabix the result.
#   6. Plant it at $STAGE_HOME/$SID/04_snv_indel/${SID}.ensemble.fixed.vcf.gz
#      so the tertiary pipeline's `--input_dir $STAGE_HOME/$SID` finds
#      it without any pipeline-side modification.
#
# Usage:
#   scripts/stage_dragen_for_tertiary.sh \
#     --in /path/to/dragen.hard-filtered.vcf.gz \
#     --sample VAL-58-dragen \
#     [--stage-home $HOME/NGS_UI/nf_stage] \
#     [--ref-fasta /home/pipeline/reference/hg38/Homo_sapiens_assembly38.fasta] \
#     [--skip-norm]
#
# After this:
#   nextflow ... --sample_id $SID --input_dir $STAGE_HOME/$SID ...
#
# Requires: bcftools (any container with it; default uses the
# pipeline's bcftools_1.23.1.sif).
# =========================================================
set -euo pipefail

BCF_SIF="${BCFTOOLS_SIF:-/home/pipeline/nextflow_containers/bcftools_1.23.1.sif}"
REF_FASTA="${REF_FASTA:-/home/pipeline/reference/hg38/Homo_sapiens_assembly38.fasta}"
GENE_BED="${GENE_BED:-$HOME/NGS_UI/biotools/gene_body.protein_coding.bed}"
GNOMAD_AF_VCF="${GNOMAD_AF_VCF:-$HOME/NGS_UI/biotools/gnomad/gnomad_af.hg38.vcf.gz}"
GNOMAD_AF_CUTOFF="${GNOMAD_AF_CUTOFF:-0.01}"

IN=""
SID=""
STAGE_HOME="${STAGE_HOME:-$HOME/NGS_UI/nf_stage}"
SKIP_NORM=0
SKIP_BED=0
SKIP_GNOMAD=0
while [ $# -gt 0 ]; do
  case "$1" in
    --in)              IN="$2"; shift 2;;
    --sample)          SID="$2"; shift 2;;
    --stage-home)      STAGE_HOME="$2"; shift 2;;
    --ref-fasta)       REF_FASTA="$2"; shift 2;;
    --gene-bed)        GENE_BED="$2"; shift 2;;
    --gnomad-af-vcf)   GNOMAD_AF_VCF="$2"; shift 2;;
    --gnomad-af-cutoff) GNOMAD_AF_CUTOFF="$2"; shift 2;;
    --skip-norm)       SKIP_NORM=1; shift;;
    --skip-bed)        SKIP_BED=1; shift;;
    --skip-gnomad)     SKIP_GNOMAD=1; shift;;
    -h|--help)         sed -n '2,40p' "$0"; exit 0;;
    *) echo "unknown arg: $1" >&2; exit 2;;
  esac
done
[ -n "$IN" ]  || { echo "ERROR: --in required"  >&2; exit 2; }
[ -n "$SID" ] || { echo "ERROR: --sample required" >&2; exit 2; }
[ -f "$IN" ]  || { echo "ERROR: input not found: $IN" >&2; exit 2; }
[ "$SKIP_NORM" -eq 1 ] || [ -f "$REF_FASTA" ] || \
  { echo "ERROR: --ref-fasta not found: $REF_FASTA (use --skip-norm to bypass)" >&2; exit 2; }
[ "$SKIP_BED" -eq 1 ] || [ -f "$GENE_BED" ] || \
  { echo "ERROR: --gene-bed not found: $GENE_BED" >&2
    echo "       build one first:" >&2
    echo "         scripts/build_gene_body_bed.sh" >&2
    echo "       or pass --skip-bed to bypass" >&2; exit 2; }
# gnomAD AF VCF is optional — silently disabled if the default path
# doesn't exist (the user may not have built it yet). Only error out
# when the user explicitly pointed --gnomad-af-vcf at a missing file.
if [ "$SKIP_GNOMAD" -eq 0 ] && [ ! -f "$GNOMAD_AF_VCF" ]; then
  if [ "${GNOMAD_AF_VCF:-}" != "$HOME/NGS_UI/biotools/gnomad/gnomad_af.hg38.vcf.gz" ]; then
    echo "ERROR: --gnomad-af-vcf not found: $GNOMAD_AF_VCF" >&2
    exit 2
  fi
  SKIP_GNOMAD=1
fi

STAGE_DIR="$STAGE_HOME/$SID/04_snv_indel"
mkdir -p "$STAGE_DIR"
OUT="$STAGE_DIR/${SID}.ensemble.fixed.vcf.gz"
TMP="$STAGE_DIR/.${SID}.stage.tmp.vcf"
MID="$STAGE_DIR/.${SID}.stage.mid.vcf.gz"

echo "[stage] in            : $IN"
echo "[stage] out           : $OUT"
echo "[stage] sample        : $SID  (→ ${SID}_DV + empty ${SID}_HC)"
[ "$SKIP_NORM" -eq 0 ]   && echo "[stage] norm          : on  ($REF_FASTA)" \
                         || echo "[stage] norm          : OFF (--skip-norm)"
[ "$SKIP_BED" -eq 0 ]    && echo "[stage] gene BED      : on  ($GENE_BED)" \
                         || echo "[stage] gene BED      : OFF (--skip-bed)"
[ "$SKIP_GNOMAD" -eq 0 ] && echo "[stage] gnomAD filter : on  (drop AF > $GNOMAD_AF_CUTOFF; $GNOMAD_AF_VCF)" \
                         || echo "[stage] gnomAD filter : OFF (no $GNOMAD_AF_VCF — build with scripts/build_gnomad_af_vcf.sh)"

# Inspect the input sample names
SAMPLES=$(apptainer exec --bind /home,"$STAGE_HOME" "$BCF_SIF" \
  bcftools query -l "$IN")
N_SAMPLES=$(echo "$SAMPLES" | wc -l)
echo "[stage] input has $N_SAMPLES sample(s): $(echo "$SAMPLES" | tr '\n' ',' | sed 's/,$//')"

RENAME_TSV="$STAGE_DIR/.rename.tsv"
: > "$RENAME_TSV"
SYNTH_HC=0   # whether to append a synthetic ${SID}_HC column

if [ "$N_SAMPLES" -eq 1 ]; then
  # DRAGEN germline path — single sample becomes the DV column, HC
  # synthesised below as missing rows.
  ORIG_SAMPLE=$(echo "$SAMPLES" | head -1)
  printf "%s\t%s_DV\n" "$ORIG_SAMPLE" "$SID" > "$RENAME_TSV"
  SYNTH_HC=1
elif [ "$N_SAMPLES" -eq 2 ]; then
  # In-house ensemble path — VCF already carries (_DV, _HC) shape.
  # Map by suffix so we don't depend on column order. If neither
  # sample carries the suffix, fall back to position (1=DV, 2=HC).
  S1=$(echo "$SAMPLES" | sed -n '1p')
  S2=$(echo "$SAMPLES" | sed -n '2p')
  case "$S1" in
    *_DV) printf "%s\t%s_DV\n" "$S1" "$SID" >> "$RENAME_TSV" ;;
    *_HC) printf "%s\t%s_HC\n" "$S1" "$SID" >> "$RENAME_TSV" ;;
    *)    printf "%s\t%s_DV\n" "$S1" "$SID" >> "$RENAME_TSV" ;;
  esac
  case "$S2" in
    *_DV) printf "%s\t%s_DV\n" "$S2" "$SID" >> "$RENAME_TSV" ;;
    *_HC) printf "%s\t%s_HC\n" "$S2" "$SID" >> "$RENAME_TSV" ;;
    *)    printf "%s\t%s_HC\n" "$S2" "$SID" >> "$RENAME_TSV" ;;
  esac
  echo "[stage] 2-sample input → rename map:"
  sed 's/^/[stage]   /' "$RENAME_TSV"
else
  echo "ERROR: this stager supports 1- or 2-sample VCFs only (got $N_SAMPLES)" >&2
  exit 2
fi

BIND_DIRS="/home,$STAGE_HOME,$(dirname "$REF_FASTA"),$(dirname "$GENE_BED")"
[ "$SKIP_GNOMAD" -eq 0 ] && BIND_DIRS="$BIND_DIRS,$(dirname "$GNOMAD_AF_VCF")"

# Phase 1 — drop chrM + (optional) gene-body BED + (optional) norm.
# Writes bgzipped intermediate. Splitting the pipeline at a real file
# avoids the long-pipe failure modes htslib runs into ("not bgzip" /
# "unknown file type") when many bcftools processes share a single
# stdin chain in some container builds.
echo "[stage] phase 1: drop chrM → gene-body BED → norm → bgzipped intermediate"
BED_STEP=""
[ "$SKIP_BED" -eq 0 ]  && BED_STEP="| bcftools view -T '$GENE_BED' -"
NORM_STEP=""
[ "$SKIP_NORM" -eq 0 ] && NORM_STEP="| bcftools norm -f '$REF_FASTA' -m -any --check-ref w --keep-sum AD"
apptainer exec --bind "$BIND_DIRS" "$BCF_SIF" bash -c "
  set -e
  bcftools view --regions-overlap pos -t ^chrM,^MT '$IN' \
    $BED_STEP \
    $NORM_STEP \
    -Oz -o '$MID'
  bcftools index -t '$MID'
"

# Phase 2 — (optional) gnomAD AF annotate + drop common + reheader sample
# columns, → final $TMP (uncompressed VCF; the synthesise-HC step
# below works on text).
echo "[stage] phase 2: gnomAD AF annotate + drop common + reheader → $TMP"
if [ "$SKIP_GNOMAD" -eq 0 ]; then
  apptainer exec --bind "$BIND_DIRS" "$BCF_SIF" bash -c "
    set -e
    bcftools annotate -a '$GNOMAD_AF_VCF' -c INFO/gnomAD_AF '$MID' \
      | bcftools view -e 'INFO/gnomAD_AF > $GNOMAD_AF_CUTOFF' \
      | bcftools reheader -s '$RENAME_TSV' \
      > '$TMP'
  "
else
  apptainer exec --bind "$BIND_DIRS" "$BCF_SIF" bash -c "
    set -e
    bcftools reheader -s '$RENAME_TSV' '$MID' \
      | bcftools view \
      > '$TMP'
  "
fi
rm -f "$MID" "$MID.tbi"

# 3. Append a synthetic ${SID}_HC sample column (DRAGEN single-sample
#    case only). For in-house 2-sample input the column is already
#    present; skip the awk step entirely.
if [ "$SYNTH_HC" -eq 1 ]; then
  echo "[stage] appending synthetic ${SID}_HC column (all ./.) …"
  awk -v sid="$SID" '
    BEGIN { OFS="\t" }
    /^##/ { print; next }
    /^#CHROM/ {
      print $0 "\t" sid "_HC"; next
    }
    {
      # Build a "missing" sample value matching the FORMAT field count.
      # FORMAT is column 9, sample data starts at 10. Synthesise one
      # extra column with "./." plus "." for each FORMAT key beyond GT.
      n_fmt = split($9, fmt_arr, ":")
      miss = "./."
      for (i = 2; i <= n_fmt; i++) miss = miss ":."
      print $0 "\t" miss
    }
  ' "$TMP" > "$TMP.synth"
else
  echo "[stage] 2-sample input — keeping ${SID}_DV + ${SID}_HC columns as-is"
  mv "$TMP" "$TMP.synth"
fi

# 4. bgzip + tabix.
echo "[stage] bgzip + tabix …"
apptainer exec --bind /home,"$STAGE_HOME" "$BCF_SIF" bash -c "
  set -e
  bgzip -f '$TMP.synth'
  mv '$TMP.synth.gz' '$OUT'
  tabix -p vcf -f '$OUT'
"
rm -f "$TMP" "$RENAME_TSV"

echo
echo "[stage] verify (#CHROM line) ====="
apptainer exec --bind /home,"$STAGE_HOME" "$BCF_SIF" \
  bcftools view -h "$OUT" | tail -1

echo
echo "[stage] ready. next step — run tertiary pipeline:"
cat <<EOF

  cd ~/NGS_UI/NGS-UI
  source /home/pipeline/pipeline_code/NGS2ndAnalysis_env.sh
  export NXF_TEMP=\$HOME/NGS_UI/nf_tmp
  export APPTAINER_TMPDIR=\$HOME/NGS_UI/apptainer_tmp
  export APPTAINER_CACHEDIR=\$HOME/NGS_UI/apptainer_cache
  mkdir -p \$NXF_TEMP \$APPTAINER_TMPDIR \$APPTAINER_CACHEDIR

  nextflow -c /home/pipeline/tertiary_code/nextflow_tertiary.config \\
    run /home/pipeline/tertiary_code/main_tertiary.nf \\
    -profile dgm \\
    -work-dir \$HOME/NGS_UI/nf_work/$SID \\
    --sample_id $SID \\
    --input_dir $STAGE_HOME/$SID \\
    --seq_type WGS \\
    --out_dir \$HOME/NGS_UI/tertiary_output \\
    2>&1 | tee \$HOME/NGS_UI/tertiary_output/$SID.nf.log

EOF
