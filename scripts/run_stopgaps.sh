#!/usr/bin/env bash
# =========================================================
# run_stopgaps.sh — one-shot stop-gap annotation chain
# =========================================================
# Runs every stop-gap on a single snv_indel.annotated.tsv in the
# order they make sense:
#
#   1. annotate_clinvar.py      — backfill CLINVAR_SIG/STARS/DN/SIGCONF
#                                  (pipeline emits all '.')
#   2. filter_snv_tsv.py        — filter MODIFIER / common / alt-contig
#                                  / '*'-allele, keep ClinVar P/LP
#   3. annotate_acmg_genebe.py  — backfill ACMG_POINTS/EVIDENCE/CLASS
#                                  via GeneBe REST API
#   4. annotate_extra_vep.py    — re-run VEP with MetaRNN (and SpliceAI
#                                  if scores VCF available)
#
# All four are idempotent fill-empty-only — re-running a step does
# nothing harmful, and once the pipeline grows the matching step
# upstream, the corresponding stop-gap silently becomes a no-op.
#
# Before step 1 we snapshot the TSV to <tsv>.raw (only the first time)
# so you can always restart from scratch with `cp <tsv>.raw <tsv>`.
#
# Usage:
#   scripts/run_stopgaps.sh --tsv tertiary_output/<SID>/<SID>.snv_indel.annotated.tsv
#
# Env / flags:
#   GENEBE_USER / GENEBE_API_KEY            — required (step 3)
#   --spliceai-snv / --spliceai-indel       — optional, default
#                                              $HOME/NGS_UI/biotools/spliceai/...
#   --skip-spliceai                         — disable SpliceAI even if
#                                              defaults exist
#   --no-backup                             — don't write .raw snapshot
# =========================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

TSV=""
NO_BACKUP=0
SKIP_SPLICEAI=0
SPLICEAI_SNV="$HOME/NGS_UI/biotools/spliceai/spliceai_scores.raw.snv.hg38.vcf.gz"
SPLICEAI_INDEL="$HOME/NGS_UI/biotools/spliceai/spliceai_scores.raw.indel.hg38.vcf.gz"
while [ $# -gt 0 ]; do
  case "$1" in
    --tsv)            TSV="$2"; shift 2;;
    --spliceai-snv)   SPLICEAI_SNV="$2"; shift 2;;
    --spliceai-indel) SPLICEAI_INDEL="$2"; shift 2;;
    --skip-spliceai)  SKIP_SPLICEAI=1; shift;;
    --no-backup)      NO_BACKUP=1; shift;;
    -h|--help) sed -n '2,30p' "$0"; exit 0;;
    *) echo "unknown arg: $1" >&2; exit 2;;
  esac
done
[ -n "$TSV" ] || { echo "ERROR: --tsv required" >&2; exit 2; }
[ -f "$TSV" ] || { echo "ERROR: --tsv not found: $TSV" >&2; exit 2; }

echo "================================================================"
echo "  run_stopgaps : $TSV"
echo "================================================================"

# 0. Snapshot.
if [ "$NO_BACKUP" -eq 0 ] && [ ! -f "$TSV.raw" ]; then
  echo "[stopgaps] 0/4  snapshot → ${TSV}.raw"
  cp -v "$TSV" "$TSV.raw"
else
  echo "[stopgaps] 0/4  snapshot skipped (already exists or --no-backup)"
fi

# 1. ClinVar.
echo
echo "[stopgaps] 1/4  annotate_clinvar.py"
"$SCRIPT_DIR/annotate_clinvar.py" --tsv "$TSV"

# 2. Filter.
echo
echo "[stopgaps] 2/4  filter_snv_tsv.py"
"$SCRIPT_DIR/filter_snv_tsv.py" --tsv "$TSV"

# 3. GeneBe ACMG.
echo
echo "[stopgaps] 3/4  annotate_acmg_genebe.py"
if [ -z "${GENEBE_USER:-}" ] || [ -z "${GENEBE_API_KEY:-}" ]; then
  echo "ERROR: GENEBE_USER + GENEBE_API_KEY must be exported" >&2
  exit 2
fi
"$SCRIPT_DIR/annotate_acmg_genebe.py" --tsv "$TSV"

# 4. Extra VEP (MetaRNN + optional SpliceAI).
echo
echo "[stopgaps] 4/4  annotate_extra_vep.py"
EXTRA_VEP_ARGS=(--tsv "$TSV")
if [ "$SKIP_SPLICEAI" -eq 0 ] && [ -f "$SPLICEAI_SNV" ] && [ -f "$SPLICEAI_INDEL" ]; then
  EXTRA_VEP_ARGS+=(
    --spliceai-snv   "$SPLICEAI_SNV"
    --spliceai-indel "$SPLICEAI_INDEL"
  )
  echo "  + SpliceAI enabled ($SPLICEAI_SNV)"
else
  if [ "$SKIP_SPLICEAI" -eq 1 ]; then
    echo "  - SpliceAI skipped (--skip-spliceai)"
  else
    echo "  - SpliceAI VCFs not found at $SPLICEAI_SNV — running MetaRNN only"
  fi
fi
"$SCRIPT_DIR/annotate_extra_vep.py" "${EXTRA_VEP_ARGS[@]}"

# 5. Drop a copy at the GUI-expected path (no sample prefix) so the
#    sample loader finds it without manual rename.
SAMPLE_DIR="$(dirname "$TSV")"
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
