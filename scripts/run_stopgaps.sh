#!/usr/bin/env bash
# =========================================================
# run_stopgaps.sh — one-shot post-processing annotation chain
# =========================================================
# Runs every post-processing step on a disposable working TSV. The immutable
# pipeline TSV remains under 03_acmg; only sparse/derived artifacts persist.
#
#   1. annotate_acmg_genebe.py  — write SECOND-opinion ACMG to GENEBE_*
#                                  columns via local DB first, then query
#                                  review-filtered DB misses through the live
#                                  API when credentials/SIF are configured;
#                                  pipeline's ACMG_* stay untouched
#   2. annotate_extra_vep.py    — add MetaRNN + REVEL + SpliceAI columns
#   3. annotate_mane_refseq.py — map Ensembl transcript IDs to MANE RefSeq
#   4. run_annotsv_cnv_sv.sh — DRAGEN sibling CNV/SV VCFs or in-house
#                              gCNV + Delly VCFs. Skipped if none supplied.
#   5. build_snv_annotation_overlay.py — persist sparse field differences
#   6. build_snv_review_tsv.py / build_snv_gene_index.py
#
# (ClinVar annotation was a compatibility step before the new pipeline shipped
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
#       [--sample SID] [--seq-type WES|WGS]
#
# Env / flags:
#   NGS_UI_GENEBE_DB / --genebe-db     — local GeneBe DB for step 1
#                                         (default $HOME/NGS_UI/biotools/
#                                         genebe/genebe_hg38.tsv.gz); builds
#                                         genebe_hg38.sqlite lazily.
#                                         --skip-genebe to disable.
#   GENEBE_USER / GENEBE_API_KEY / GENEBE_SIF
#                                      — optional live fallback; only DB misses
#                                        retained by the review TSV filter are
#                                        submitted. Results are cached and
#                                        saved as import-ready seven-column TSV.
#   NGS_UI_CDS_CANDIDATE_BED           — optional, default
#                                         $HOME/NGS_UI/biotools/cds_combined.bed
#   --spliceai-snv / --spliceai-indel  — optional, default
#                                         $HOME/NGS_UI/biotools/spliceai/...
#   NGS_UI_EXTRA_VEP_DBNSFP / --extra-vep-dbnsfp
#                                      — dbNSFP VEP-ready BGZF used by Extra
#                                        VEP; default
#                                        $NGS_UI_HOME/biotools/dbnsfp/
#                                        dbNSFP5.3.1a_grch38.gz
#   --candidate-bed / --skip-candidate-bed
#                                      — restrict GeneBe/Extra VEP candidates
#   --skip-spliceai / --skip-extra-vep — disable SpliceAI or all of step 2
#   NGS_UI_MANE_SUMMARY / --mane-summary
#                                      — MANE summary for RefSeq mapping
#   --skip-mane-refseq                 — disable Ensembl→RefSeq mapping
#   --skip-cnv                         — disable AnnotSV step 4 even
#                                         when --dragen-cnv-source is set
# =========================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

timestamp() { TZ=Asia/Taipei date +"%Y-%m-%dT%H:%M:%S+08:00"; }
step_start() {
  STOPGAP_STEP="$1"
  STOPGAP_MARKER="${2:-post-processing-step}"
  STOPGAP_STEP_STARTED="$(date +%s)"
  echo "[$(timestamp)] [$STOPGAP_MARKER] $STOPGAP_STEP start"
}
step_done() {
  local ended
  ended="$(date +%s)"
  echo "[$(timestamp)] [$STOPGAP_MARKER] $STOPGAP_STEP done elapsed=$((ended - STOPGAP_STEP_STARTED))s"
}
run_silent_step() {
  local label="$1"
  shift
  local tmp_log
  tmp_log="$(mktemp "${TMPDIR:-/tmp}/ngs-${label}.XXXXXX.log")"
  if ! "$@" >"$tmp_log" 2>&1; then
    cat "$tmp_log" >&2
    rm -f "$tmp_log"
    exit 1
  fi
  rm -f "$tmp_log"
}

TSV=""
RAW_TSV=""
POST_DIR=""
SID=""
DRAGEN_VCF=""
INHOUSE_CNV_VCF=""
INHOUSE_SV_VCF=""
SKIP_SPLICEAI=0
SKIP_EXTRA_VEP=0
SKIP_CANDIDATE_BED=0
SKIP_CNV=0
SKIP_GIAB=0
SKIP_GENEBE=0
SKIP_MANE_REFSEQ=0
SKIP_INHOUSE_AF=0
SEQ_TYPE="${SEQ_TYPE:-WES}"
NGS_HOME_DEFAULT="${NGS_UI_HOME:-$HOME/NGS_UI}"
GENEBE_DB="${NGS_UI_GENEBE_DB:-$HOME/NGS_UI/biotools/genebe/genebe_hg38.tsv.gz}"
MANE_SUMMARY="${NGS_UI_MANE_SUMMARY:-$NGS_HOME_DEFAULT/biotools/MANE.GRCh38.v1.5.summary.txt.gz}"
GIAB_STRAT_DIR="${NGS_UI_GIAB_STRAT_DIR:-}"
INHOUSE_AF_DB="${NGS_UI_INHOUSE_AF_DB:-$NGS_HOME_DEFAULT/biotools/inhouse_af/inhouse_af.hg38.vcf.gz}"
EXTRA_VEP_DBNSFP="${NGS_UI_EXTRA_VEP_DBNSFP:-$NGS_HOME_DEFAULT/biotools/dbnsfp/dbNSFP5.3.1a_grch38.gz}"
# Keep the established production SpliceAI paths. On n102968 these resolve to
# /home/n102968/NGS_UI/biotools/spliceai/...; do not point runtime at a desktop
# /Volumes mount.
SPLICEAI_SNV="$HOME/NGS_UI/biotools/spliceai/spliceai_scores.raw.snv.hg38.vcf.gz"
SPLICEAI_INDEL="$HOME/NGS_UI/biotools/spliceai/spliceai_scores.raw.indel.hg38.vcf.gz"
CANDIDATE_BED="${NGS_UI_CDS_CANDIDATE_BED:-$HOME/NGS_UI/biotools/cds_combined.bed}"
while [ $# -gt 0 ]; do
  case "$1" in
    --tsv|--work-tsv)     TSV="$2"; shift 2;;
    --raw-tsv)            RAW_TSV="$2"; shift 2;;
    --post-dir)           POST_DIR="$2"; shift 2;;
    --sample)             SID="$2"; shift 2;;
    --seq-type)           SEQ_TYPE="$2"; shift 2;;
    --dragen-cnv-source)  DRAGEN_VCF="$2"; shift 2;;
    --inhouse-cnv-vcf)    INHOUSE_CNV_VCF="$2"; shift 2;;
    --inhouse-sv-vcf)     INHOUSE_SV_VCF="$2"; shift 2;;
    --spliceai-snv)       SPLICEAI_SNV="$2"; shift 2;;
    --spliceai-indel)     SPLICEAI_INDEL="$2"; shift 2;;
    --extra-vep-dbnsfp)   EXTRA_VEP_DBNSFP="$2"; shift 2;;
    --candidate-bed)      CANDIDATE_BED="$2"; shift 2;;
    --skip-candidate-bed) SKIP_CANDIDATE_BED=1; shift;;
    --skip-spliceai)      SKIP_SPLICEAI=1; shift;;
    --skip-extra-vep)     SKIP_EXTRA_VEP=1; shift;;
    --skip-cnv)           SKIP_CNV=1; shift;;
    --genebe-db)          GENEBE_DB="$2"; shift 2;;
    --skip-genebe)        SKIP_GENEBE=1; shift;;
    --mane-summary)       MANE_SUMMARY="$2"; shift 2;;
    --skip-mane-refseq)   SKIP_MANE_REFSEQ=1; shift;;
    --giab-strat-dir)     GIAB_STRAT_DIR="$2"; shift 2;;
    --skip-giab)          SKIP_GIAB=1; shift;;
    --inhouse-af-db)      INHOUSE_AF_DB="$2"; shift 2;;
    --skip-inhouse-af)    SKIP_INHOUSE_AF=1; shift;;
    -h|--help) sed -n '2,40p' "$0"; exit 0;;
    *) echo "unknown arg: $1" >&2; exit 2;;
  esac
done
[ -n "$TSV" ] || { echo "ERROR: --tsv required" >&2; exit 2; }
[ -f "$TSV" ] || { echo "ERROR: --tsv not found: $TSV" >&2; exit 2; }
if [ -z "$RAW_TSV" ]; then RAW_TSV="$TSV"; fi
[ -f "$RAW_TSV" ] || { echo "ERROR: --raw-tsv not found: $RAW_TSV" >&2; exit 2; }
if [ -z "$POST_DIR" ]; then POST_DIR="$(dirname "$TSV")"; fi
mkdir -p "$POST_DIR"
SEQ_TYPE="$(printf '%s' "$SEQ_TYPE" | tr '[:lower:]' '[:upper:]')"
case "$SEQ_TYPE" in
  WES|WGS) ;;
  *) echo "ERROR: --seq-type must be WES or WGS (got: $SEQ_TYPE)" >&2; exit 2;;
esac
# Derive sample id from path if not supplied. Working TSVs now live below
# tertiary_output/<SID>/08_postprocessing/, while legacy invocations may
# still point directly below tertiary_output/<SID>/.
if [ -z "$SID" ]; then
  TSV_DIR="$(dirname "$TSV")"
  if [ "$(basename "$TSV_DIR")" = "08_postprocessing" ]; then
    SID="$(basename "$(dirname "$TSV_DIR")")"
  else
    SID="$(basename "$TSV_DIR")"
  fi
fi
OVERLAY_PATH="$POST_DIR/$SID.snv_annotations.sqlite"
REVIEW_PATH="$POST_DIR/$SID.snv_indel.review.tsv"
REVIEW_MANIFEST_PATH="$POST_DIR/$SID.snv_indel.review.tsv.source.json"
GENE_INDEX_PATH="$POST_DIR/$SID.snv_gene_index.sqlite"

echo "================================================================"
echo "  post-processing : $TSV"
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

# 1. GeneBe ACMG second opinion. The local DB is always queried first across
#    the complete TSV. Live API fallback is best-effort and only sees DB misses
#    that meet the same WES/WGS filter as review.tsv.
if [ "$SKIP_GENEBE" -eq 0 ]; then
  echo
  echo "[post-processing] annotate_acmg_genebe.py"
  step_start "genebe"
  if [ ! -f "$GENEBE_DB" ]; then
    echo "ERROR: GeneBe DB not found: $GENEBE_DB" >&2
    echo "       set NGS_UI_GENEBE_DB / --genebe-db, or pass --skip-genebe" >&2
    exit 2
  fi
  # Whole-TSV lookup (no candidate gate) — the DB read is one streaming
  # pass regardless of how many variants are queried.
  "$SCRIPT_DIR/annotate_acmg_genebe.py" \
    --tsv "$TSV" --genebe-db "$GENEBE_DB" --test-type "$SEQ_TYPE"
  step_done
fi

# 3. Extra VEP (MetaRNN + REVEL when available + optional SpliceAI). Skippable.
if [ "$SKIP_EXTRA_VEP" -eq 0 ]; then
  echo
  echo "[post-processing] annotate_extra_vep.py"
  step_start "extra-vep"
  EXTRA_VEP_ARGS=(--tsv "$TSV" --dbnsfp "$EXTRA_VEP_DBNSFP")
  echo "  + dbNSFP: $EXTRA_VEP_DBNSFP"
  if [ "$SKIP_SPLICEAI" -eq 0 ] && [ -f "$SPLICEAI_SNV" ] && [ -f "$SPLICEAI_INDEL" ]; then
    EXTRA_VEP_ARGS+=(--spliceai-snv "$SPLICEAI_SNV" --spliceai-indel "$SPLICEAI_INDEL")
    echo "  + SpliceAI enabled ($SPLICEAI_SNV)"
  else
    if [ "$SKIP_SPLICEAI" -eq 1 ]; then
      echo "  - SpliceAI skipped (--skip-spliceai)"
    else
      echo "  - SpliceAI VCFs not found at $SPLICEAI_SNV — dbNSFP predictors only"
    fi
  fi
  "$SCRIPT_DIR/annotate_extra_vep.py" "${EXTRA_VEP_ARGS[@]}" "${CANDIDATE_BED_ARGS[@]}"
  step_done
fi

# 3b. GIAB stratification labels (difficult-region badges). Runs on the
#     whole TSV; cheap interval lookup. No-op when the BED dir is absent.
if [ "$SKIP_GIAB" -eq 0 ]; then
  echo
  echo "[post-processing] annotate_giab_strata.py"
  step_start "giab-strata"
  GIAB_ARGS=(--tsv "$TSV")
  if [ -n "$GIAB_STRAT_DIR" ]; then
    GIAB_ARGS+=(--strat-dir "$GIAB_STRAT_DIR")
  fi
  "$SCRIPT_DIR/annotate_giab_strata.py" "${GIAB_ARGS[@]}"
  step_done
fi

# 3b2. In-house allele frequency. Joins INHOUSE_AC/AN/AF from the local cohort
#      sites VCF (single streaming pass). Runs on the whole TSV, before review
#      TSV / gene index so the columns reach the UI. No-op when the DB is absent.
if [ "$SKIP_INHOUSE_AF" -eq 0 ]; then
  if [ -f "$INHOUSE_AF_DB" ]; then
    echo
    echo "[post-processing] annotate_inhouse_af.py"
    step_start "inhouse-af"
    "$SCRIPT_DIR/annotate_inhouse_af.py" --tsv "$TSV" --db "$INHOUSE_AF_DB"
    step_done
  else
    echo "  - in-house AF skipped (DB not found: $INHOUSE_AF_DB)"
  fi
fi

# 3c. MANE RefSeq mapping. Runs before review TSV / gene index so every
#     derived artifact carries the display transcript IDs.
if [ "$SKIP_MANE_REFSEQ" -eq 0 ]; then
  echo
  echo "[post-processing] annotate_mane_refseq.py"
  step_start "mane-refseq"
  "$SCRIPT_DIR/annotate_mane_refseq.py" --tsv "$TSV" --mane-summary "$MANE_SUMMARY"
  step_done
fi

# 4. CNV/SV via AnnotSV.
SAMPLE_DIR="$POST_DIR"
if [ "$SKIP_CNV" -eq 0 ] && { [ -n "$DRAGEN_VCF" ] || [ -n "$INHOUSE_CNV_VCF" ] || [ -n "$INHOUSE_SV_VCF" ]; }; then
  echo
  echo "[post-processing] run_annotsv_cnv_sv.sh"
  step_start "annotsv"
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
  step_done
fi

echo
echo "[sample] sparse overlay  build_snv_annotation_overlay.py"
step_start "snv-overlay" "sample-step"
run_silent_step "snv-overlay" "$SCRIPT_DIR/build_snv_annotation_overlay.py" \
  --raw "$RAW_TSV" --annotated "$TSV" --out "$OVERLAY_PATH"
step_done

echo
echo "[sample] review TSV  build_snv_review_tsv.py"
step_start "review-tsv" "sample-step"
run_silent_step "review-tsv" "$SCRIPT_DIR/build_snv_review_tsv.py" \
  --tsv "$RAW_TSV" --output-dir "$POST_DIR" \
  --output-path "$REVIEW_PATH" --manifest-path "$REVIEW_MANIFEST_PATH" \
  --overlay "$OVERLAY_PATH" --test-type "$SEQ_TYPE"
step_done

echo
echo "[sample] gene index  build_snv_gene_index.py"
step_start "gene-index" "sample-step"
run_silent_step "gene-index" "$SCRIPT_DIR/build_snv_gene_index.py" \
  --tsv "$RAW_TSV" --out "$GENE_INDEX_PATH"
step_done

echo
echo "================================================================"
echo "  done. raw TSV: $RAW_TSV"
echo "  sparse overlay: $OVERLAY_PATH"
echo "================================================================"
wc -l "$RAW_TSV"
ls -la "$SAMPLE_DIR"/"$SID".{cnv,sv,mito}.annotated.tsv 2>/dev/null || true
