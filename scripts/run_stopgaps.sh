#!/usr/bin/env bash
# =========================================================
# run_stopgaps.sh — one-shot stop-gap annotation chain
# =========================================================
# Runs every stop-gap on a single snv_indel.annotated.tsv in the
# order they make sense:
#
#   1. filter_snv_tsv.py        — filter alt-contig / * only
#   2. annotate_acmg_genebe.py  — write SECOND-opinion ACMG to GENEBE_*
#                                  columns (pipeline's ACMG_* untouched)
#   3. annotate_extra_vep.py    — add MetaRNN + SpliceAI as new columns
#   4. run_annotsv_cnv_sv.sh — DRAGEN sibling CNV/SV VCFs or in-house
#                              gCNV + Delly VCFs. Skipped if none supplied.
#   5. build_snv_review_tsv.py  — pre-build the compact main-screen TSV
#   6. build_snv_gene_index.py  — pre-build complete-TSV gene search index
#
# (ClinVar annotation was a stop-gap before the new pipeline shipped
#  CLINVAR_SIG / STARS / DN / SIGCONF / VARIATION_ID natively. The
#  pipeline now owns it; the legacy annotate_clinvar.py is left on
#  disk for emergencies.)
#
# All steps are idempotent fill-or-augment. Step 4 produces fresh
# cnv.annotated.tsv / sv.annotated.tsv whenever CNV VCFs are
# pointed to; pass --skip-cnv to bypass.
#
# Usage:
#   scripts/run_stopgaps.sh \\
#       --tsv tertiary_output/<SID>/<SID>.snv_indel.annotated.tsv \\
#       [--dragen-cnv-source /path/to/<sample>.hard-filtered.vcf.gz] \\
#       [--sample SID]
#
# Env / flags:
#   GENEBE_USER / GENEBE_API_KEY       — required (step 3)
#   NGS_UI_CDS_CANDIDATE_BED           — optional, default
#                                         $HOME/NGS_UI/biotools/cds_combined.bed
#   --spliceai-snv / --spliceai-indel  — optional, default
#                                         $HOME/NGS_UI/biotools/spliceai/...
#   --candidate-bed / --skip-candidate-bed
#                                      — restrict GeneBe/Extra VEP candidates
#   --skip-spliceai / --skip-extra-vep — disable MetaRNN/SpliceAI step 4
#   --skip-cnv                         — disable AnnotSV step 5 even
#                                         when --dragen-cnv-source is set
# =========================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

timestamp() { date -u +"%Y-%m-%dT%H:%M:%SZ"; }
step_start() {
  STOPGAP_STEP="$1"
  STOPGAP_STEP_STARTED="$(date +%s)"
  echo "[$(timestamp)] [stopgaps-step] $STOPGAP_STEP start"
}
step_done() {
  local ended
  ended="$(date +%s)"
  echo "[$(timestamp)] [stopgaps-step] $STOPGAP_STEP done elapsed=$((ended - STOPGAP_STEP_STARTED))s"
}

TSV=""
SID=""
DRAGEN_VCF=""
INHOUSE_CNV_VCF=""
INHOUSE_SV_VCF=""
SKIP_SPLICEAI=0
SKIP_EXTRA_VEP=0
SKIP_CANDIDATE_BED=0
SKIP_CNV=0
SPLICEAI_SNV="$HOME/NGS_UI/biotools/spliceai/spliceai_scores.raw.snv.hg38.vcf.gz"
SPLICEAI_INDEL="$HOME/NGS_UI/biotools/spliceai/spliceai_scores.raw.indel.hg38.vcf.gz"
CANDIDATE_BED="${NGS_UI_CDS_CANDIDATE_BED:-$HOME/NGS_UI/biotools/cds_combined.bed}"
while [ $# -gt 0 ]; do
  case "$1" in
    --tsv)                TSV="$2"; shift 2;;
    --sample)             SID="$2"; shift 2;;
    --dragen-cnv-source)  DRAGEN_VCF="$2"; shift 2;;
    --inhouse-cnv-vcf)    INHOUSE_CNV_VCF="$2"; shift 2;;
    --inhouse-sv-vcf)     INHOUSE_SV_VCF="$2"; shift 2;;
    --spliceai-snv)       SPLICEAI_SNV="$2"; shift 2;;
    --spliceai-indel)     SPLICEAI_INDEL="$2"; shift 2;;
    --candidate-bed)      CANDIDATE_BED="$2"; shift 2;;
    --skip-candidate-bed) SKIP_CANDIDATE_BED=1; shift;;
    --skip-spliceai)      SKIP_SPLICEAI=1; shift;;
    --skip-extra-vep)     SKIP_EXTRA_VEP=1; shift;;
    --skip-cnv)           SKIP_CNV=1; shift;;
    -h|--help) sed -n '2,40p' "$0"; exit 0;;
    *) echo "unknown arg: $1" >&2; exit 2;;
  esac
done
[ -n "$TSV" ] || { echo "ERROR: --tsv required" >&2; exit 2; }
[ -f "$TSV" ] || { echo "ERROR: --tsv not found: $TSV" >&2; exit 2; }
# Derive sample id from path if not supplied: tertiary_output/<SID>/...
if [ -z "$SID" ]; then SID="$(basename "$(dirname "$TSV")")"; fi

echo "================================================================"
echo "  run_stopgaps : $TSV"
echo "================================================================"

CANDIDATE_BED_ARGS=()
if [ "$SKIP_CANDIDATE_BED" -eq 1 ]; then
  echo "  candidate BED: skipped (--skip-candidate-bed)"
elif [ -f "$CANDIDATE_BED" ]; then
  CANDIDATE_BED_ARGS=(--candidate-bed "$CANDIDATE_BED")
  echo "  candidate BED: $CANDIDATE_BED"
else
  echo "  candidate BED: not found at $CANDIDATE_BED (GeneBe/Extra VEP use AF-only candidates)"
fi

# 1. ClinVar fallback. The new pipeline's CLINVAR_SIG column has been
# flaky (often blank); annotate_clinvar.py is fill-empty-only, so
# pipeline-supplied values win when present. Once the pipeline's
# ClinVar step is reliable this step turns into a no-op.
echo
echo "[stopgaps] 1/5  annotate_clinvar.py (fill-empty fallback)"
step_start "clinvar"
"$SCRIPT_DIR/annotate_clinvar.py" --tsv "$TSV"
step_done

# 2. Filter.
echo
echo "[stopgaps] 2/5  filter_snv_tsv.py"
step_start "filter-snv"
"$SCRIPT_DIR/filter_snv_tsv.py" --tsv "$TSV"
step_done

# 3. GeneBe ACMG (writes to GENEBE_* columns; pipeline ACMG_* preserved).
echo
echo "[stopgaps] 3/5  annotate_acmg_genebe.py"
step_start "genebe"
if [ -z "${GENEBE_USER:-}" ] || [ -z "${GENEBE_API_KEY:-}" ]; then
  echo "ERROR: GENEBE_USER + GENEBE_API_KEY must be exported" >&2
  exit 2
fi
"$SCRIPT_DIR/annotate_acmg_genebe.py" --tsv "$TSV" "${CANDIDATE_BED_ARGS[@]}"
step_done

# 4. Extra VEP (MetaRNN + optional SpliceAI). Skippable.
echo
echo "[stopgaps] 4/5  annotate_extra_vep.py"
step_start "extra-vep"
if [ "$SKIP_EXTRA_VEP" -eq 1 ]; then
  echo "  - skipped (--skip-extra-vep)"
else
  EXTRA_VEP_ARGS=(--tsv "$TSV")
  if [ "$SKIP_SPLICEAI" -eq 0 ] && [ -f "$SPLICEAI_SNV" ] && [ -f "$SPLICEAI_INDEL" ]; then
    EXTRA_VEP_ARGS+=(--spliceai-snv "$SPLICEAI_SNV" --spliceai-indel "$SPLICEAI_INDEL")
    echo "  + SpliceAI enabled ($SPLICEAI_SNV)"
  else
    if [ "$SKIP_SPLICEAI" -eq 1 ]; then
      echo "  - SpliceAI skipped (--skip-spliceai)"
    else
      echo "  - SpliceAI VCFs not found at $SPLICEAI_SNV — MetaRNN only"
    fi
  fi
  "$SCRIPT_DIR/annotate_extra_vep.py" "${EXTRA_VEP_ARGS[@]}" "${CANDIDATE_BED_ARGS[@]}"
fi
step_done

# 5. CNV/SV via AnnotSV.
SAMPLE_DIR="$(dirname "$TSV")"
echo
echo "[stopgaps] 5/5  AnnotSV CNV/SV"
step_start "annotsv"
if [ "$SKIP_CNV" -eq 1 ]; then
  echo "  - skipped (--skip-cnv)"
else
  ANNOTSV_ARGS=(--sample "$SID" --out-dir "$SAMPLE_DIR")
  if [ -n "$DRAGEN_VCF" ]; then
    ANNOTSV_ARGS+=(--dragen-cnv-source "$DRAGEN_VCF")
  fi
  if [ -n "$INHOUSE_CNV_VCF" ]; then
    ANNOTSV_ARGS+=(--inhouse-cnv-vcf "$INHOUSE_CNV_VCF")
  fi
  if [ -n "$INHOUSE_SV_VCF" ]; then
    ANNOTSV_ARGS+=(--inhouse-sv-vcf "$INHOUSE_SV_VCF")
  fi
  "$SCRIPT_DIR/run_annotsv_cnv_sv.sh" "${ANNOTSV_ARGS[@]}"
fi
step_done

# Drop a copy at the GUI-expected path (no sample prefix).
GUI_TSV="$SAMPLE_DIR/snv_indel.annotated.tsv"
if [ "$TSV" != "$GUI_TSV" ]; then
  echo
  echo "[stopgaps] copy → $GUI_TSV  (GUI-expected path)"
  cp -v "$TSV" "$GUI_TSV"
fi

echo
echo "[stopgaps] review TSV  build_snv_review_tsv.py"
step_start "review-tsv"
"$SCRIPT_DIR/build_snv_review_tsv.py" --tsv "$GUI_TSV"
step_done

echo
echo "[stopgaps] gene index  build_snv_gene_index.py"
step_start "gene-index"
"$SCRIPT_DIR/build_snv_gene_index.py" --tsv "$GUI_TSV"
step_done

echo
echo "================================================================"
echo "  done. final TSV: $GUI_TSV"
echo "================================================================"
wc -l "$GUI_TSV"
ls -la "$SAMPLE_DIR"/{cnv,sv,mito}.annotated.tsv 2>/dev/null || true
