#!/usr/bin/env bash
# Download the official static GRCh38 GPN-MSA pre-computed score table.
#
# The destination is intentionally a single fixed BGZF + TBI pair. There is no
# release/current hierarchy and no update timer because this is a pre-computed
# immutable reference used only for review-TSV annotation.
set -euo pipefail

NGS_HOME_DEFAULT="${NGS_UI_HOME:-$HOME/NGS_UI}"
DEST="${NGS_UI_GPN_MSA_DB:-$NGS_HOME_DEFAULT/biotools/gpn_msa/scores.tsv.bgz}"
BASE_URL="https://huggingface.co/datasets/songlab/gpn-msa-hg38-scores/resolve/main"
FORCE=0

while [ $# -gt 0 ]; do
  case "$1" in
    --dest) DEST="$2"; shift 2;;
    --url) BASE_URL="${2%/}"; shift 2;;
    --force) FORCE=1; shift;;
    -h|--help)
      sed -n '2,11p' "$0"
      exit 0
      ;;
    *) echo "unknown arg: $1" >&2; exit 2;;
  esac
done

if [ "$FORCE" -ne 1 ] && { [ -e "$DEST" ] || [ -e "$DEST.tbi" ]; }; then
  echo "ERROR: destination already exists; use --force to replace: $DEST" >&2
  exit 2
fi
command -v curl >/dev/null 2>&1 || { echo "ERROR: curl is required" >&2; exit 2; }
command -v tabix >/dev/null 2>&1 || { echo "ERROR: tabix is required" >&2; exit 2; }

mkdir -p "$(dirname "$DEST")"
DATA_PART="$DEST.download"
INDEX_PART="$DEST.tbi.download"

echo "Downloading GPN-MSA scores (large file; curl will resume a partial download)..."
curl --fail --location --continue-at - \
  --output "$DATA_PART" "$BASE_URL/scores.tsv.bgz"
curl --fail --location --continue-at - \
  --output "$INDEX_PART" "$BASE_URL/scores.tsv.bgz.tbi"

# Probe a documented BRCA1 interval before switching either destination file.
cp "$INDEX_PART" "$DATA_PART.tbi"
if ! tabix "$DATA_PART" "17:43044295-43044295" | awk -F '\t' '
  $1 == "17" && $2 == "43044295" && $3 == "T" && NF >= 5 { found=1 }
  END { exit(found ? 0 : 1) }
'; then
  echo "ERROR: downloaded GPN-MSA BGZF/TBI pair failed the sentinel query" >&2
  exit 1
fi
rm -f "$DATA_PART.tbi"
mv -f "$DATA_PART" "$DEST"
mv -f "$INDEX_PART" "$DEST.tbi"

echo "Installed: $DEST"
echo "Source/license: songlab/gpn-msa-hg38-scores (MIT)"
