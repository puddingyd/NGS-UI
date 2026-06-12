#!/usr/bin/env bash
# =========================================================
# run_stopgaps.sh — one-shot stop-gap annotation chain
# =========================================================
# Runs every stop-gap on a single snv_indel.annotated.tsv in the
# order they make sense:
#
#   1. filter_snv_tsv.py        — filter alt-contig / * only (silent)
#   2. annotate_acmg_genebe.py  — write SECOND-opinion ACMG to GENEBE_*
#                                  columns via the local GeneBe DB
#                                  (offline SQLite lookup, no API/creds;
#                                  pipeline's ACMG_* untouched)
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
#   NGS_UI_GENEBE_DB / --genebe-db     — local GeneBe DB for step 2
#                                         (default $HOME/NGS_UI/biotools/
#                                         genebe/genebe_hg38.tsv.gz); builds
#                                         genebe_hg38.sqlite lazily.
#                                         --skip-genebe to disable.
#   NGS_UI_CDS_CANDIDATE_BED           — optional, default
#                                         $HOME/NGS_UI/biotools/cds_combined.bed
#   --spliceai-snv / --spliceai-indel  — optional, default
#                                         $HOME/NGS_UI/biotools/spliceai/...
#   --candidate-bed / --skip-candidate-bed
#                                      — restrict GeneBe/Extra VEP candidates
#   --skip-spliceai / --skip-extra-vep — disable MetaRNN/SpliceAI step 3
#   --skip-cnv                         — disable AnnotSV step 4 even
#                                         when --dragen-cnv-source is set
# =========================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

timestamp() { TZ=Asia/Taipei date +"%Y-%m-%dT%H:%M:%S+08:00"; }
step_start() {
  STOPGAP_STEP="$1"
  STOPGAP_MARKER="${2:-stopgaps-step}"
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
GENEBE_DB="${NGS_UI_GENEBE_DB:-$HOME/NGS_UI/biotools/genebe/genebe_hg38.tsv.gz}"
GIAB_STRAT_DIR="${NGS_UI_GIAB_STRAT_DIR:-}"
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
    --genebe-db)          GENEBE_DB="$2"; shift 2;;
    --skip-genebe)        SKIP_GENEBE=1; shift;;
    --giab-strat-dir)     GIAB_STRAT_DIR="$2"; shift 2;;
    --skip-giab)          SKIP_GIAB=1; shift;;
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

# 1. Filter. Keep the cleanup, but keep it out of the user-facing log.
run_silent_step "filter-snv" "$SCRIPT_DIR/filter_snv_tsv.py" --tsv "$TSV"

# 2. GeneBe ACMG second opinion via the local DB (writes GENEBE_* columns;
#    pipeline ACMG_* preserved). Offline SQLite lookup — no API / creds.
if [ "$SKIP_GENEBE" -eq 0 ]; then
  echo
  echo "[stopgaps] annotate_acmg_genebe.py"
  step_start "genebe"
  if [ ! -f "$GENEBE_DB" ]; then
    echo "ERROR: GeneBe DB not found: $GENEBE_DB" >&2
    echo "       set NGS_UI_GENEBE_DB / --genebe-db, or pass --skip-genebe" >&2
    exit 2
  fi
  # Whole-TSV lookup (no candidate gate) — the DB read is one streaming
  # pass regardless of how many variants are queried.
  "$SCRIPT_DIR/annotate_acmg_genebe.py" --tsv "$TSV" --genebe-db "$GENEBE_DB"
  step_done
fi

# 3. Extra VEP (MetaRNN + optional SpliceAI). Skippable.
if [ "$SKIP_EXTRA_VEP" -eq 0 ]; then
  echo
  echo "[stopgaps] annotate_extra_vep.py"
  step_start "extra-vep"
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
  step_done
fi

# 3b. GIAB stratification labels (difficult-region badges). Runs on the
#     whole TSV; cheap interval lookup. No-op when the BED dir is absent.
if [ "$SKIP_GIAB" -eq 0 ]; then
  echo
  echo "[stopgaps] annotate_giab_strata.py"
  step_start "giab-strata"
  GIAB_ARGS=(--tsv "$TSV")
  if [ -n "$GIAB_STRAT_DIR" ]; then
    GIAB_ARGS+=(--strat-dir "$GIAB_STRAT_DIR")
  fi
  "$SCRIPT_DIR/annotate_giab_strata.py" "${GIAB_ARGS[@]}"
  step_done
fi

# 4. CNV/SV via AnnotSV.
SAMPLE_DIR="$(dirname "$TSV")"
if [ "$SKIP_CNV" -eq 0 ] && { [ -n "$DRAGEN_VCF" ] || [ -n "$INHOUSE_CNV_VCF" ] || [ -n "$INHOUSE_SV_VCF" ]; }; then
  echo
  echo "[stopgaps] run_annotsv_cnv_sv.sh"
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

# Drop a copy at the GUI-expected path (no sample prefix).
GUI_TSV="$SAMPLE_DIR/snv_indel.annotated.tsv"
if [ "$TSV" != "$GUI_TSV" ]; then
  echo
  echo "[stopgaps] copy → $GUI_TSV  (GUI-expected path)"
  cp -v "$TSV" "$GUI_TSV"
fi

echo
echo "[sample] review TSV  build_snv_review_tsv.py"
step_start "review-tsv" "sample-step"
run_silent_step "review-tsv" "$SCRIPT_DIR/build_snv_review_tsv.py" --tsv "$GUI_TSV"
step_done

echo
echo "[sample] gene index  build_snv_gene_index.py"
step_start "gene-index" "sample-step"
run_silent_step "gene-index" "$SCRIPT_DIR/build_snv_gene_index.py" --tsv "$GUI_TSV"
step_done

echo
echo "================================================================"
echo "  done. final TSV: $GUI_TSV"
echo "================================================================"
wc -l "$GUI_TSV"
ls -la "$SAMPLE_DIR"/{cnv,sv,mito}.annotated.tsv 2>/dev/null || true
