#!/usr/bin/env bash
# =========================================================
# run_stopgaps.sh — one-shot stop-gap annotation chain
# =========================================================
# Runs every stop-gap on a single snv_indel.annotated.tsv in the
# order they make sense:
#
#   1. annotate_clinvar.py      — backfill CLINVAR_SIG/STARS/DN/SIGCONF
#   2. filter_snv_tsv.py        — filter MODIFIER / common / alt-contig / *
#   3. annotate_acmg_genebe.py  — backfill ACMG_* via GeneBe REST
#   4. annotate_extra_vep.py    — re-run VEP for MetaRNN (+ SpliceAI)
#   5. annotate_dragen_cnv_annotsv.sh — DRAGEN <sample>.cnv.vcf.gz +
#                                  cnv_sv.vcf.gz → AnnotSV TSV (only
#                                  when --dragen-cnv-source supplied)
#
# Steps 1-4 are idempotent fill-empty-only. Step 5 produces fresh
# cnv.annotated.tsv / sv.annotated.tsv whenever DRAGEN CNV VCFs are
# pointed to; pass --skip-cnv to bypass.
#
# Before step 1 we snapshot the TSV to <tsv>.raw (only the first time)
# so you can always restart from scratch with `cp <tsv>.raw <tsv>`.
#
# Usage:
#   scripts/run_stopgaps.sh \\
#       --tsv tertiary_output/<SID>/<SID>.snv_indel.annotated.tsv \\
#       [--dragen-cnv-source /path/to/<sample>.hard-filtered.vcf.gz] \\
#       [--sample SID]
#
# Env / flags:
#   GENEBE_USER / GENEBE_API_KEY       — required (step 3)
#   --spliceai-snv / --spliceai-indel  — optional, default
#                                         $HOME/NGS_UI/biotools/spliceai/...
#   --skip-spliceai / --skip-extra-vep — disable MetaRNN/SpliceAI step 4
#   --skip-cnv                         — disable AnnotSV step 5 even
#                                         when --dragen-cnv-source is set
#   --no-backup                        — don't write .raw snapshot
# =========================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

TSV=""
SID=""
DRAGEN_VCF=""
NO_BACKUP=0
SKIP_SPLICEAI=0
SKIP_EXTRA_VEP=0
SKIP_CNV=0
SPLICEAI_SNV="$HOME/NGS_UI/biotools/spliceai/spliceai_scores.raw.snv.hg38.vcf.gz"
SPLICEAI_INDEL="$HOME/NGS_UI/biotools/spliceai/spliceai_scores.raw.indel.hg38.vcf.gz"
while [ $# -gt 0 ]; do
  case "$1" in
    --tsv)                TSV="$2"; shift 2;;
    --sample)             SID="$2"; shift 2;;
    --dragen-cnv-source)  DRAGEN_VCF="$2"; shift 2;;
    --spliceai-snv)       SPLICEAI_SNV="$2"; shift 2;;
    --spliceai-indel)     SPLICEAI_INDEL="$2"; shift 2;;
    --skip-spliceai)      SKIP_SPLICEAI=1; shift;;
    --skip-extra-vep)     SKIP_EXTRA_VEP=1; shift;;
    --skip-cnv)           SKIP_CNV=1; shift;;
    --no-backup)          NO_BACKUP=1; shift;;
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

# 0. Snapshot.
if [ "$NO_BACKUP" -eq 0 ] && [ ! -f "$TSV.raw" ]; then
  echo "[stopgaps] 0/5  snapshot → ${TSV}.raw"
  cp -v "$TSV" "$TSV.raw"
else
  echo "[stopgaps] 0/5  snapshot skipped (already exists or --no-backup)"
fi

# 1. ClinVar.
echo
echo "[stopgaps] 1/5  annotate_clinvar.py"
"$SCRIPT_DIR/annotate_clinvar.py" --tsv "$TSV"

# 2. Filter.
echo
echo "[stopgaps] 2/5  filter_snv_tsv.py"
"$SCRIPT_DIR/filter_snv_tsv.py" --tsv "$TSV"

# 3. GeneBe ACMG.
echo
echo "[stopgaps] 3/5  annotate_acmg_genebe.py"
if [ -z "${GENEBE_USER:-}" ] || [ -z "${GENEBE_API_KEY:-}" ]; then
  echo "ERROR: GENEBE_USER + GENEBE_API_KEY must be exported" >&2
  exit 2
fi
"$SCRIPT_DIR/annotate_acmg_genebe.py" --tsv "$TSV"

# 4. Extra VEP (MetaRNN + optional SpliceAI). Skippable.
echo
echo "[stopgaps] 4/5  annotate_extra_vep.py"
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
  "$SCRIPT_DIR/annotate_extra_vep.py" "${EXTRA_VEP_ARGS[@]}"
fi

# 5. DRAGEN CNV/SV via AnnotSV. Only when --dragen-cnv-source supplied
#    (i.e. we're processing a DRAGEN VCF whose <sample>.cnv.vcf.gz +
#    cnv_sv.vcf.gz live in the same directory).
SAMPLE_DIR="$(dirname "$TSV")"
echo
echo "[stopgaps] 5/5  annotate_dragen_cnv_annotsv.sh"
if [ "$SKIP_CNV" -eq 1 ]; then
  echo "  - skipped (--skip-cnv)"
elif [ -z "$DRAGEN_VCF" ]; then
  echo "  - skipped (no --dragen-cnv-source)"
elif [ ! -f "$DRAGEN_VCF" ]; then
  echo "  - skipped (--dragen-cnv-source not found: $DRAGEN_VCF)"
else
  "$SCRIPT_DIR/annotate_dragen_cnv_annotsv.sh" \
    --dragen-vcf "$DRAGEN_VCF" \
    --sample "$SID" \
    --out-dir "$SAMPLE_DIR"
fi

# Drop a copy at the GUI-expected path (no sample prefix).
GUI_TSV="$SAMPLE_DIR/snv_indel.annotated.tsv"
if [ "$TSV" != "$GUI_TSV" ]; then
  echo
  echo "[stopgaps] copy → $GUI_TSV  (GUI-expected path)"
  cp -v "$TSV" "$GUI_TSV"
fi

echo
echo "================================================================"
echo "  done. final TSV: $GUI_TSV"
echo "================================================================"
wc -l "$GUI_TSV"
ls -la "$SAMPLE_DIR"/{cnv,sv,mito}.annotated.tsv 2>/dev/null || true
