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

# DRAGEN's cnv_sv.vcf.gz from the Manta-style SV caller also re-calls
# the large DEL/DUP that the dedicated CNV caller emits. Annotating
# both as-is duplicates events between the GUI's CNV and SV cards.
# Pre-filter cnv_sv to drop any record whose (CHROM, POS, INFO/END,
# INFO/SVTYPE) already appears in cnv.vcf.gz; BND breakpoints have
# SVTYPE=BND so they never collide and stay in SV. The filtered VCF
# is what AnnotSV actually sees, so ranking_score isn't polluted.
subtract_cnv_from_sv() {
  local cnv="$1" sv="$2" out="$3"
  [ -f "$cnv" ] && [ -f "$sv" ] || return 1
  local keys
  keys=$(mktemp --suffix=.keys.tsv)
  zcat "$cnv" | awk -F'\t' '
    /^#/ { next }
    {
      end=""; svt=""
      n=split($8, kv, ";")
      for (i=1;i<=n;i++) {
        if (kv[i] ~ /^END=/)    { end=kv[i];  sub(/^END=/,    "", end)  }
        if (kv[i] ~ /^SVTYPE=/) { svt=kv[i];  sub(/^SVTYPE=/, "", svt)  }
      }
      print $1 "\t" $2 "\t" end "\t" svt
    }
  ' > "$keys"
  local n_cnv n_sv n_kept
  n_cnv=$(wc -l < "$keys")
  n_sv=$(zcat "$sv" | awk '!/^#/' | wc -l)
  zcat "$sv" | awk -F'\t' -v kf="$keys" '
    BEGIN {
      while ((getline line < kf) > 0) keys[line] = 1
      close(kf)
    }
    /^#/ { print; next }
    {
      end=""; svt=""
      n=split($8, kv, ";")
      for (i=1;i<=n;i++) {
        if (kv[i] ~ /^END=/)    { end=kv[i];  sub(/^END=/,    "", end)  }
        if (kv[i] ~ /^SVTYPE=/) { svt=kv[i];  sub(/^SVTYPE=/, "", svt)  }
      }
      key = $1 "\t" $2 "\t" end "\t" svt
      if (!(key in keys)) print
    }
  ' > "$out"
  n_kept=$(awk '!/^#/' "$out" | wc -l)
  rm -f "$keys"
  echo "[annotsv] sv-pre-filter: cnv keys=$n_cnv, cnv_sv rows=$n_sv → kept=$n_kept (dropped $((n_sv - n_kept)) duplicates)"
}

# Process CNV first (it's the keys source). Then pre-filter cnv_sv
# against those keys and annotate the SV-only remainder.
run_annotsv "$CNV_VCF" cnv

if [ -f "$CNV_VCF" ] && [ -f "$SV_VCF" ]; then
  SV_ONLY="$OUT_DIR/.${SID}.sv_only.vcf"
  subtract_cnv_from_sv "$CNV_VCF" "$SV_VCF" "$SV_ONLY"
  run_annotsv "$SV_ONLY" sv
  rm -f "$SV_ONLY"
else
  # No CNV VCF to subtract against — annotate cnv_sv as-is.
  run_annotsv "$SV_VCF"  sv
fi

ls -la "$OUT_DIR"/{cnv,sv}.annotated.tsv 2>/dev/null || true
