#!/usr/bin/env bash
# =========================================================
# stage_dragen_for_tertiary.sh
# =========================================================
# Stage a DRAGEN germline hard-filter VCF so the existing tertiary
# pipeline (which expects ensemble DV+HC input) accepts it. Does:
#
#   1. Drop chrM rows (mtDNA is handled separately by
#      scripts/annotate_mito_vcf.sh — DRAGEN's chrM calls have
#      heteroplasmy in FORMAT/AF, no TLOD; the mito script's
#      auto-detect handles that).
#   2. bcftools norm — left-align indels, split multi-allelics. DRAGEN
#      raw VCFs ship un-normalised; downstream tools (Pangolin
#      especially) segfault on the corner cases (IUPAC bases, un-split
#      multi-allelics, etc.).
#   3. Reheader the single DRAGEN sample column to `${SID}_DV`.
#   4. Append a synthetic `${SID}_HC` sample column populated with
#      "./." for every row so PREPARE_VCF:ADD_CALLERS_TAG sees the
#      expected (DV, HC) shape. CALLERS = DV everywhere (DRAGEN-only).
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

IN=""
SID=""
STAGE_HOME="${STAGE_HOME:-$HOME/NGS_UI/nf_stage}"
SKIP_NORM=0
while [ $# -gt 0 ]; do
  case "$1" in
    --in)         IN="$2"; shift 2;;
    --sample)     SID="$2"; shift 2;;
    --stage-home) STAGE_HOME="$2"; shift 2;;
    --ref-fasta)  REF_FASTA="$2"; shift 2;;
    --skip-norm)  SKIP_NORM=1; shift;;
    -h|--help)    sed -n '2,40p' "$0"; exit 0;;
    *) echo "unknown arg: $1" >&2; exit 2;;
  esac
done
[ -n "$IN" ]  || { echo "ERROR: --in required"  >&2; exit 2; }
[ -n "$SID" ] || { echo "ERROR: --sample required" >&2; exit 2; }
[ -f "$IN" ]  || { echo "ERROR: input not found: $IN" >&2; exit 2; }
[ "$SKIP_NORM" -eq 1 ] || [ -f "$REF_FASTA" ] || \
  { echo "ERROR: --ref-fasta not found: $REF_FASTA (use --skip-norm to bypass)" >&2; exit 2; }

STAGE_DIR="$STAGE_HOME/$SID/04_snv_indel"
mkdir -p "$STAGE_DIR"
OUT="$STAGE_DIR/${SID}.ensemble.fixed.vcf.gz"
TMP="$STAGE_DIR/.${SID}.stage.tmp.vcf"

echo "[stage] in     : $IN"
echo "[stage] out    : $OUT"
echo "[stage] sample : $SID  (will end up as ${SID}_DV + empty ${SID}_HC)"
[ "$SKIP_NORM" -eq 0 ] && echo "[stage] norm   : on  (bcftools norm -f $REF_FASTA -m -any)" \
                       || echo "[stage] norm   : OFF (--skip-norm)"

# Inspect the input sample names
SAMPLES=$(apptainer exec --bind /home,"$STAGE_HOME" "$BCF_SIF" \
  bcftools query -l "$IN")
echo "[stage] input has $(echo "$SAMPLES" | wc -l) sample(s): $(echo "$SAMPLES" | tr '\n' ',' | sed 's/,$//')"
N_SAMPLES=$(echo "$SAMPLES" | wc -l)
if [ "$N_SAMPLES" -ne 1 ]; then
  echo "ERROR: this stager only supports single-sample DRAGEN VCFs (got $N_SAMPLES)" >&2
  exit 2
fi
ORIG_SAMPLE=$(echo "$SAMPLES" | head -1)

# 1+2+3. Drop chrM rows + normalise + rename the single sample to ${SID}_DV.
echo "[stage] dropping chrM + normalising + renaming sample → ${SID}_DV …"
RENAME_TSV="$STAGE_DIR/.rename.tsv"
printf "%s\t%s_DV\n" "$ORIG_SAMPLE" "$SID" > "$RENAME_TSV"

# Build the bcftools pipeline. The norm step is what unblocks Pangolin
# (DRAGEN un-normalised input causes segfaults on IUPAC bases / multi-
# allelics). --check-ref w warns instead of erroring on ref mismatches
# (DRAGEN occasionally calls on alt contigs we've already excluded,
# but belt-and-braces). --keep-sum AD preserves AD sums when splitting
# multi-allelics so DP_DV / AD_DV stay consistent.
NORM_STEP=""
if [ "$SKIP_NORM" -eq 0 ]; then
  NORM_STEP="| bcftools norm -f '$REF_FASTA' -m -any --check-ref w --keep-sum AD"
fi
apptainer exec --bind /home,"$STAGE_HOME","$(dirname "$REF_FASTA")" "$BCF_SIF" bash -c "
  set -e
  bcftools view --regions-overlap pos -t ^chrM,^MT '$IN' \
    | bcftools reheader -s '$RENAME_TSV' \
    $NORM_STEP \
    > '$TMP'
"

# 3. Append a synthetic ${SID}_HC sample column. We do this in pure
#    awk — much cleaner than trying to splice a per-row FORMAT through
#    bcftools merge. For every data row append a tab + './.' + ':.' for
#    each remaining FORMAT key (so column count stays valid). For the
#    #CHROM header line append the new sample name.
echo "[stage] appending synthetic ${SID}_HC column (all ./.) …"
awk -v sid="$SID" '
  BEGIN { OFS="\t" }
  /^##/ { print; next }
  /^#CHROM/ {
    print $0 "\t" sid "_HC"; next
  }
  {
    # Build a "missing" sample value matching the FORMAT field count.
    # FORMAT is column 9, sample data starts at 10. We synthesise one
    # extra column with "./." plus "." for each FORMAT key beyond GT.
    n_fmt = split($9, fmt_arr, ":")
    miss = "./."
    for (i = 2; i <= n_fmt; i++) miss = miss ":."
    print $0 "\t" miss
  }
' "$TMP" > "$TMP.synth"

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
