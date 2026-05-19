#!/usr/bin/env bash
# =========================================================
# annotate_inhouse_cnv_sv_annotsv.sh — in-house CNV/SV → AnnotSV TSV
# =========================================================
# In-house ensemble Nextflow emits CNV and SV with two different
# callers (GATK gCNV → <SID>.gcnv.vcf.gz and Delly → <SID>.delly.vcf.gz)
# whose event sets don't overlap, so unlike DRAGEN there is no
# cnv-vs-sv subtraction step here — just run AnnotSV on each and drop
# the result into the GUI's expected files.
#
# Usage:
#   scripts/annotate_inhouse_cnv_sv_annotsv.sh \\
#     --cnv-vcf /path/to/<SID>.gcnv.vcf.gz \\
#     --sv-vcf  /path/to/<SID>.delly.vcf.gz \\
#     --sample  VAL-31-inhouse \\
#     --out-dir $HOME/NGS_UI/tertiary_output/VAL-31-inhouse/
#
# Either --cnv-vcf or --sv-vcf may be empty/missing — that step is then
# skipped. AnnotSV env (ANNOTSV_BIN, ANNOTSV_ANNOTATIONS) and output
# layout match annotate_dragen_cnv_annotsv.sh exactly so the GUI's
# annotsv_tsv.py adapter consumes both pipelines uniformly.
# =========================================================
set -euo pipefail

ANNOTSV_BIN="${ANNOTSV_BIN:-$HOME/NGS_UI/biotools/AnnotSV/bin/AnnotSV}"
ANNOTSV_ANNOTATIONS="${ANNOTSV_ANNOTATIONS:-$HOME/NGS_UI/biotools/AnnotSV/share/AnnotSV}"

CNV_VCF=""
SV_VCF=""
SID=""
OUT_DIR=""
while [ $# -gt 0 ]; do
  case "$1" in
    --cnv-vcf) CNV_VCF="$2"; shift 2;;
    --sv-vcf)  SV_VCF="$2";  shift 2;;
    --sample)  SID="$2";     shift 2;;
    --out-dir) OUT_DIR="$2"; shift 2;;
    -h|--help) sed -n '2,30p' "$0"; exit 0;;
    *) echo "unknown arg: $1" >&2; exit 2;;
  esac
done
[ -n "$SID" ] && [ -n "$OUT_DIR" ] || {
  echo "usage: $0 --cnv-vcf <path> --sv-vcf <path> --sample <SID> --out-dir <dir>" >&2
  exit 2; }
[ -x "$ANNOTSV_BIN" ] || {
  echo "ERROR: AnnotSV not found / not executable: $ANNOTSV_BIN" >&2; exit 2; }
[ -d "$ANNOTSV_ANNOTATIONS" ] || {
  echo "ERROR: AnnotSV annotations dir not found: $ANNOTSV_ANNOTATIONS" >&2; exit 2; }

mkdir -p "$OUT_DIR"

run_annotsv() {
  local input_vcf="$1" kind="$2"   # kind: cnv | sv
  if [ -z "$input_vcf" ]; then
    echo "[annotsv] $kind: (no input) — skipping"
    return 0
  fi
  [ -f "$input_vcf" ] || { echo "[annotsv] $kind: not found $input_vcf — skipping"; return 0; }

  echo "[annotsv] $kind: $input_vcf"
  local tmp_out="$OUT_DIR/_annotsv_${kind}"
  rm -rf "$tmp_out"; mkdir -p "$tmp_out"
  local out_name="${SID}.${kind}.annotated.tsv"

  "$ANNOTSV_BIN" \
    -SVinputFile "$input_vcf" \
    -outputDir "$tmp_out" \
    -outputFile "$out_name" \
    -genomeBuild GRCh38 \
    -annotationsDir "$ANNOTSV_ANNOTATIONS" \
    -SVinputInfo 1 \
    >&2 || { echo "[annotsv] $kind: AnnotSV failed (see stderr)" >&2; return 1; }

  local produced=""
  for cand in "$tmp_out/$out_name" "$tmp_out/${out_name%.tsv}.tsv" "$tmp_out"/*.annotated.tsv; do
    [ -f "$cand" ] && { produced="$cand"; break; }
  done
  if [ -z "$produced" ]; then
    echo "[annotsv] $kind: no output TSV in $tmp_out" >&2
    ls -la "$tmp_out" >&2
    return 1
  fi

  local dst="$OUT_DIR/${kind}.annotated.tsv"
  mv "$produced" "$dst"
  rm -rf "$tmp_out"
  local n
  n=$(awk 'END{print NR-1}' "$dst")
  echo "[annotsv] $kind → $dst  ($n rows)"
}

run_annotsv "$CNV_VCF" cnv
run_annotsv "$SV_VCF"  sv

ls -la "$OUT_DIR"/{cnv,sv}.annotated.tsv 2>/dev/null || true
