#!/usr/bin/env bash
# =========================================================
# ingest_batch.sh — Phase A in parallel over many DRAGEN samples
# =========================================================
# Runs ingest_sample.py for each sample's other/{id}/ dir, in parallel.
#
# Usage:
#   scripts/inhouse_af/ingest_batch.sh \
#     --out-dir $NGS_UI_HOME/biotools/inhouse_af/per_sample \
#     [--ref <fasta>] [--jobs 8] \
#     ( --run /home/datalake_Raw/Novaseq/<RUN>   # all other/* under a run
#     | --dirs-file <file of sample dirs>         # one other/{id} dir per line
#     | <dir> [<dir> ...] )                        # explicit sample dirs
#
# Skips a sample whose per_sample/{id}/qc.json already exists (resumable).
# bcftools via BCFTOOLS_SIF / BCFTOOLS_BIN (same as ingest_sample.py).
# =========================================================
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

REF="${REF:-/home/datalake_Intermediate/pipeline/reference/hg38/Homo_sapiens_assembly38.fasta}"
OUT_DIR=""
JOBS=8
RUN=""
DIRS_FILE=""
DIRS=()

while [ $# -gt 0 ]; do
  case "$1" in
    --out-dir)   OUT_DIR="$2"; shift 2;;
    --ref)       REF="$2"; shift 2;;
    --jobs)      JOBS="$2"; shift 2;;
    --run)       RUN="$2"; shift 2;;
    --dirs-file) DIRS_FILE="$2"; shift 2;;
    -h|--help)   sed -n '2,22p' "$0"; exit 0;;
    -*)          echo "unknown arg: $1" >&2; exit 2;;
    *)           DIRS+=("$1"); shift;;
  esac
done

[ -n "$OUT_DIR" ] || { echo "ERROR: --out-dir required" >&2; exit 2; }
mkdir -p "$OUT_DIR"

# Build the sample-dir list
LIST="$(mktemp)"
trap 'rm -f "$LIST"' EXIT
if [ -n "$RUN" ]; then
  # other/<id>/ dirs that actually contain a gVCF
  find "$RUN/other" -maxdepth 2 -name '*.hard-filtered.gvcf.gz' 2>/dev/null \
    | xargs -r -n1 dirname | sort -u >> "$LIST"
fi
[ -n "$DIRS_FILE" ] && cat "$DIRS_FILE" >> "$LIST"
for d in "${DIRS[@]:-}"; do [ -n "$d" ] && echo "$d" >> "$LIST"; done

N=$(grep -c . "$LIST" || true)
[ "$N" -gt 0 ] || { echo "ERROR: no sample dirs to ingest" >&2; exit 2; }
echo "[ingest-batch] $N sample dirs, jobs=$JOBS, out=$OUT_DIR"

export REF OUT_DIR SCRIPT_DIR
one() {
  d="$1"
  id=$(basename "$(find "$d" -maxdepth 2 -name '*.hard-filtered.gvcf.gz' | head -1)" | sed 's/\..*//')
  [ -n "$id" ] || { echo "[skip] no gVCF in $d"; return 0; }
  if [ -f "$OUT_DIR/$id/qc.json" ]; then echo "[skip] $id (already ingested)"; return 0; fi
  python3 "$SCRIPT_DIR/ingest_sample.py" --sample-dir "$d" --ref "$REF" --out-dir "$OUT_DIR" \
    >/dev/null 2>"$OUT_DIR/.$id.log" \
    && echo "[ok] $id" || { echo "[FAIL] $id (see $OUT_DIR/.$id.log)"; }
}
export -f one

grep . "$LIST" | xargs -r -P "$JOBS" -I{} bash -c 'one "$@"' _ {}
echo "[ingest-batch] done."
