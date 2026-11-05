#!/usr/bin/env bash
set -euo pipefail

PROJECT=/home/n102968/2026_INTERVENE_PRS
PIPE="$PROJECT/pipeline/PGS004879_integrated_v2"

RAW_BASE=/home/datalake_Raw/Novaseq/20260923_LH00873_0036_B255WGJLT4/other

RUN_ID=method2_v2_fresh_validation_01

LOG_DIR="$PROJECT/integrated_runs/$RUN_ID/batch_logs"
mkdir -p "$LOG_DIR"

cd "$PIPE"

SAMPLES=(
  26T00970
  26T00969
  26T00968
)

echo "=== PGS004879 FRESH VALIDATION ==="
echo "Pipeline: $PIPE"
echo "Run ID:   $RUN_ID"
echo "Samples:  ${SAMPLES[*]}"

# Verify all inputs before running anything.
for SAMPLE in "${SAMPLES[@]}"; do

  GVCF="$RAW_BASE/$SAMPLE/germline_seq/$SAMPLE.hard-filtered.gvcf.gz"

  echo
  echo "Checking $SAMPLE"

  [[ -s "$GVCF" ]] || {
    echo "ERROR: Missing gVCF: $GVCF"
    exit 1
  }

  [[ -s "$GVCF.tbi" || -s "$GVCF.csi" ]] || {
    echo "ERROR: Missing gVCF index: $GVCF"
    exit 1
  }

  ACTUAL_SAMPLE="$(bcftools query -l "$GVCF")"

  [[ "$ACTUAL_SAMPLE" == "$SAMPLE" ]] || {
    echo "ERROR: Sample ID mismatch: $ACTUAL_SAMPLE"
    exit 1
  }

  DEST="$PROJECT/integrated_runs/$RUN_ID/$SAMPLE"

  [[ ! -e "$DEST" ]] || {
    echo "ERROR: Output already exists: $DEST"
    echo "Use a new RUN_ID; do not overwrite."
    exit 1
  }

  echo "INPUT CHECK: PASS"

done

# Execute each sample independently, in sequence.
for SAMPLE in "${SAMPLES[@]}"; do

  GVCF="$RAW_BASE/$SAMPLE/germline_seq/$SAMPLE.hard-filtered.gvcf.gz"

  echo
  echo "=========================================="
  echo "START FRESH RUN: $SAMPLE"
  echo "=========================================="

  python3 pipeline_v2.py run \
    --config config.json \
    --mode fresh \
    --sample "$SAMPLE" \
    --gvcf "$GVCF" \
    --run-id "$RUN_ID" \
    2>&1 | tee "$LOG_DIR/${SAMPLE}.console.log"

  STATUS="$PROJECT/integrated_runs/$RUN_ID/$SAMPLE/method2/status.json"

  python3 - "$STATUS" <<'PY'
import json, sys
from pathlib import Path

p = Path(sys.argv[1])
assert p.is_file(), f"Missing status: {p}"
d = json.loads(p.read_text())
assert d.get("status") == "PASS", d
print("METHOD2 STATUS: PASS")
PY

  echo "FRESH RUN COMPLETED: $SAMPLE"

done

echo
echo "=========================================="
echo "ALL THREE FRESH RUNS COMPLETED"
echo "=========================================="
