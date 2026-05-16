#!/usr/bin/env bash
# =========================================================
# annotate_dragen_cnv_annotsv.sh — DRAGEN CNV/SV → AnnotSV TSV
# =========================================================
# DRAGEN writes two SV-style VCFs alongside the SNV/Indel VCF:
#   <sample>.cnv.vcf.gz         copy-number events (CNV; DEL/DUP)
#   <sample>.cnv_sv.vcf.gz      breakpoint-level structural variants
#
# Run AnnotSV on each, drop the result into the same directory the
# GUI reads (`tertiary_output/<SID>/cnv.annotated.tsv` and
# `sv.annotated.tsv`). GUI's existing annotsv_tsv.py adapter consumes
# both — no UI changes needed.
#
# AnnotSV install (user's home):
#   $HOME/NGS_UI/biotools/AnnotSV/bin/AnnotSV
#   $HOME/NGS_UI/biotools/AnnotSV/share/AnnotSV/
# Override via env:
#   ANNOTSV_BIN          path to AnnotSV binary
#   ANNOTSV_ANNOTATIONS  path to annotation root (Annotations_Human/, ...)
#
# Usage:
#   scripts/annotate_dragen_cnv_annotsv.sh \\
#     --dragen-vcf /path/to/<sample>.hard-filtered.vcf.gz \\
#     --sample VAL-58-dragen \\
#     --out-dir $HOME/NGS_UI/tertiary_output/VAL-58-dragen/
#
# Source CNV/SV VCFs are looked up by replacing the `.hard-filtered`
# suffix on the same directory.
# =========================================================
set -euo pipefail

ANNOTSV_BIN="${ANNOTSV_BIN:-$HOME/NGS_UI/biotools/AnnotSV/bin/AnnotSV}"
ANNOTSV_ANNOTATIONS="${ANNOTSV_ANNOTATIONS:-$HOME/NGS_UI/biotools/AnnotSV/share/AnnotSV}"

DRAGEN_VCF=""
SID=""
OUT_DIR=""
while [ $# -gt 0 ]; do
  case "$1" in
    --dragen-vcf) DRAGEN_VCF="$2"; shift 2;;
    --sample)     SID="$2"; shift 2;;
    --out-dir)    OUT_DIR="$2"; shift 2;;
    -h|--help)    sed -n '2,30p' "$0"; exit 0;;
    *) echo "unknown arg: $1" >&2; exit 2;;
  esac
done
[ -n "$DRAGEN_VCF" ] && [ -n "$SID" ] && [ -n "$OUT_DIR" ] || {
  echo "usage: $0 --dragen-vcf <path> --sample <SID> --out-dir <dir>" >&2
  exit 2; }
[ -f "$DRAGEN_VCF" ] || { echo "ERROR: --dragen-vcf not found: $DRAGEN_VCF" >&2; exit 2; }
[ -x "$ANNOTSV_BIN" ] || {
  echo "ERROR: AnnotSV not found / not executable: $ANNOTSV_BIN" >&2
  echo "       install AnnotSV under \$HOME/NGS_UI/biotools/AnnotSV or set ANNOTSV_BIN" >&2
  exit 2; }
[ -d "$ANNOTSV_ANNOTATIONS" ] || {
  echo "ERROR: AnnotSV annotations dir not found: $ANNOTSV_ANNOTATIONS" >&2
  echo "       set ANNOTSV_ANNOTATIONS to the install's share/AnnotSV/ root" >&2
  exit 2; }

mkdir -p "$OUT_DIR"

# DRAGEN convention: <sample>.{hard-filtered,cnv,cnv_sv}.vcf.gz all
# live in the same directory; sample base is the prefix before
# `.hard-filtered`.
DRAGEN_DIR=$(dirname "$DRAGEN_VCF")
BASE=$(basename "$DRAGEN_VCF" .hard-filtered.vcf.gz)
CNV_VCF="$DRAGEN_DIR/$BASE.cnv.vcf.gz"
SV_VCF="$DRAGEN_DIR/$BASE.cnv_sv.vcf.gz"

run_annotsv() {
  local input_vcf="$1" kind="$2"   # kind: cnv | sv
  [ -f "$input_vcf" ] || { echo "[annotsv] $kind: no $input_vcf — skipping"; return 0; }

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

  # AnnotSV may name the output as supplied or with subtle variations
  # (.annotated.tsv suffix is standard but check defensively).
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
  local n=$(awk 'END{print NR-1}' "$dst")  # subtract header
  echo "[annotsv] $kind → $dst  ($n rows)"
}

run_annotsv "$CNV_VCF" cnv
run_annotsv "$SV_VCF"  sv

ls -la "$OUT_DIR"/{cnv,sv}.annotated.tsv 2>/dev/null || true
