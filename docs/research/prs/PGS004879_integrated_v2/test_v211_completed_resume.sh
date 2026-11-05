#!/usr/bin/env bash
set -euo pipefail

PROJECT=/home/n102968/2026_INTERVENE_PRS
PIPE="$PROJECT/pipeline/PGS004879_integrated_v2"

RUN_ID=method2_v211_fresh_validation_01
RUN="$PROJECT/integrated_runs/$RUN_ID"

cd "$PIPE"

BEFORE="$(mktemp)"
AFTER="$(mktemp)"
trap 'rm -f "$BEFORE" "$AFTER"' EXIT

snapshot() {
    for SAMPLE in 26T00970 26T00969 26T00968; do
        for FILE in scores.tsv summary.json inputs.lock.json; do
            TARGET="$RUN/$SAMPLE/method2/$FILE"
            sha256sum "$TARGET"
        done
    done
}

echo "=== BEFORE RESUME ==="
snapshot > "$BEFORE"

echo
echo "=== TEST COMPLETED-RUN RESUME ==="

python3 pipeline_v21.py batch \
    --config config.json \
    --manifest manifest.v211_fresh_3cases.tsv \
    --run-id "$RUN_ID" \
    --jobs 1 \
    --resume

echo
echo "=== AFTER RESUME ==="
snapshot > "$AFTER"

echo
echo "=== FILE INTEGRITY CHECK ==="

diff -u "$BEFORE" "$AFTER"

echo
echo "COMPLETED-RUN RESUME: PASS"
echo "ALL 9 CORE OUTPUT FILES: IDENTICAL"
