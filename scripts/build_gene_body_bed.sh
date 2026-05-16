#!/usr/bin/env bash
# =========================================================
# build_gene_body_bed.sh — GENCODE GTF → BED for pre-filtering
# =========================================================
# One-time script: extract `feature == "gene"` rows from the GENCODE
# GTF and emit a sorted BED of gene-body intervals. Default keeps
# only protein-coding genes (the only category clinically reported
# in this pipeline). Use --include-lncRNA / --include-all to widen.
#
# Output is consumed by stage_dragen_for_tertiary.sh's `bcftools
# view -T` step — drops intergenic variants (~50% of WGS calls) so
# downstream VEP / Pangolin processes fewer records.
#
# Usage:
#   scripts/build_gene_body_bed.sh                       # → $HOME/NGS_UI/biotools/gene_body.protein_coding.bed
#   scripts/build_gene_body_bed.sh --include-lncRNA      # also keep lincRNA / lncRNA
#   scripts/build_gene_body_bed.sh --include-all         # every gene biotype
#   scripts/build_gene_body_bed.sh --gtf /path/to/file.gtf.gz --out /path/to/output.bed
# =========================================================
set -euo pipefail

GTF="${GTF:-/home/pipeline/reference/hg38/tertiary/pangolin/gencode.v47.annotation.gtf.gz}"
OUT="${OUT:-$HOME/NGS_UI/biotools/gene_body.protein_coding.bed}"
INCLUDE="protein_coding"          # default
while [ $# -gt 0 ]; do
  case "$1" in
    --gtf)              GTF="$2"; shift 2;;
    --out)              OUT="$2"; shift 2;;
    --include-lncRNA)   INCLUDE="protein_coding_lncRNA"; shift;;
    --include-all)      INCLUDE="all"; shift;;
    -h|--help)          sed -n '2,20p' "$0"; exit 0;;
    *) echo "unknown arg: $1" >&2; exit 2;;
  esac
done
[ -f "$GTF" ] || { echo "ERROR: GTF not found: $GTF" >&2; exit 2; }
mkdir -p "$(dirname "$OUT")"

echo "[bed] GTF     : $GTF"
echo "[bed] include : $INCLUDE"
echo "[bed] out     : $OUT"

# GTF is 1-based inclusive; BED is 0-based half-open → subtract 1 from start.
# We pick `feature == "gene"` rows so each gene contributes one interval
# (whole body, introns + exons + UTRs included — important for splice
# analysis downstream).
case "$INCLUDE" in
  protein_coding)
    FILT='/gene_type "protein_coding"/' ;;
  protein_coding_lncRNA)
    FILT='/gene_type "(protein_coding|lncRNA|lincRNA)"/' ;;
  all)
    FILT='1' ;;
esac

zcat "$GTF" \
  | awk -v filt="$FILT" -F'\t' '
      $1 ~ /^#/ { next }
      $3 == "gene" {
        if (eval(filt)) {
          print $1 "\t" ($4 - 1) "\t" $5
        }
      }
      function eval(_) { return 1 }   # awk has no eval; pattern below
    ' 2>/dev/null \
  | sort -k1,1 -k2,2n \
  > "$OUT.unfilt"

# awk's pattern matching doesn't let us splice the regex from a variable
# cleanly, so apply the include-filter as a second awk pass.
case "$INCLUDE" in
  protein_coding)
    zcat "$GTF" \
      | awk -F'\t' '$1 !~ /^#/ && $3 == "gene" && $9 ~ /gene_type "protein_coding"/ \
                    { print $1 "\t" ($4 - 1) "\t" $5 }' \
      | sort -k1,1 -k2,2n > "$OUT" ;;
  protein_coding_lncRNA)
    zcat "$GTF" \
      | awk -F'\t' '$1 !~ /^#/ && $3 == "gene" && $9 ~ /gene_type "(protein_coding|lncRNA|lincRNA)"/ \
                    { print $1 "\t" ($4 - 1) "\t" $5 }' \
      | sort -k1,1 -k2,2n > "$OUT" ;;
  all)
    zcat "$GTF" \
      | awk -F'\t' '$1 !~ /^#/ && $3 == "gene" { print $1 "\t" ($4 - 1) "\t" $5 }' \
      | sort -k1,1 -k2,2n > "$OUT" ;;
esac
rm -f "$OUT.unfilt"

LINES=$(wc -l < "$OUT")
COV=$(awk '{s += $3 - $2} END{printf "%.1f Mb", s/1e6}' "$OUT")
echo "[bed] done : $LINES intervals, total coverage ~$COV"
echo "[bed] head -3 :"
head -3 "$OUT" | sed 's/^/  /'
