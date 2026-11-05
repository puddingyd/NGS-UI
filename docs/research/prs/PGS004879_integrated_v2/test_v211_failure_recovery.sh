#!/usr/bin/env bash
set -euo pipefail

PROJECT=/home/n102968/2026_INTERVENE_PRS
PIPE="$PROJECT/pipeline/PGS004879_integrated_v2"

SAMPLE=26T00970
RUN_ID=method2_v211_recovery_test_01

GVCF="/home/datalake_Raw/Novaseq/20260923_LH00873_0036_B255WGJLT4/other/$SAMPLE/germline_seq/$SAMPLE.hard-filtered.gvcf.gz"

DEST="$PROJECT/integrated_runs/$RUN_ID/$SAMPLE"
RAW="$DEST/raw"

cd "$PIPE"

# Never overwrite any previous test run.
if [[ -e "$PROJECT/integrated_runs/$RUN_ID" ]]; then
    echo "ERROR: Test run already exists."
    echo "Use a new RUN_ID; do not remove existing outputs."
    exit 1
fi

[[ -s "$GVCF" ]] || {
    echo "ERROR: Missing gVCF"
    exit 1
}

REAL_NEXTFLOW="$(command -v nextflow)"
[[ -x "$REAL_NEXTFLOW" ]] || {
    echo "ERROR: Nextflow executable not found"
    exit 1
}

TESTDIR="$(mktemp -d)"
SHIM="$TESTDIR/shim"
mkdir -p "$SHIM"

MANIFEST="$TESTDIR/manifest.tsv"
MARKER="$TESTDIR/fault_injected.marker"

printf 'sample\tmode\tgvcf\traw_source\n' > "$MANIFEST"
printf '%s\tfresh\t%s\t\n' "$SAMPLE" "$GVCF" >> "$MANIFEST"

# Intercept only "nextflow run" on the first attempt.
cat > "$SHIM/nextflow" <<'SH'
#!/usr/bin/env bash
set -euo pipefail

if [[ "${1:-}" == "run" ]]; then
    echo "CONTROLLED_NEXTFLOW_FAILURE_INJECTED" >&2
    printf 'injected\n' > "$PGS_V21_FAULT_MARKER"
    exit 86
fi

exec "$PGS_V21_REAL_NEXTFLOW" "$@"
SH

chmod +x "$SHIM/nextflow"

echo "=== TEST CONFIGURATION ==="
echo "Sample: $SAMPLE"
echo "Run: $RUN_ID"
echo "Real Nextflow: $REAL_NEXTFLOW"
echo "Test logs: $TESTDIR"

echo
echo "=== PHASE 1: EXPECTED FAILURE ==="

if env \
    PATH="$SHIM:$PATH" \
    PGS_V21_REAL_NEXTFLOW="$REAL_NEXTFLOW" \
    PGS_V21_FAULT_MARKER="$MARKER" \
    python3 pipeline_v21.py batch \
        --config config.json \
        --manifest "$MANIFEST" \
        --run-id "$RUN_ID" \
        --jobs 1 \
        2>&1 | tee "$TESTDIR/first_attempt.log"
then
    echo "ERROR: Expected failure did not occur."
    exit 1
fi

[[ -s "$MARKER" ]] || {
    echo "ERROR: Failure was not injected at Nextflow."
    exit 1
}

echo
echo "=== PHASE 2: VERIFY UPSTREAM CHECKPOINT ==="

for FILE in \
    "$RAW/prep/score_ready.vcf.gz" \
    "$RAW/samplesheet.csv" \
    "$RAW/roundtrip_22chr/pgen/chr1.pgen"
do
    [[ -s "$FILE" ]] || {
        echo "ERROR: Required upstream checkpoint missing: $FILE"
        exit 1
    }
done

[[ ! -f "$RAW/audit/final_reconciliation.json" ]] || {
    echo "ERROR: Raw reconciliation unexpectedly completed."
    exit 1
}

sha256sum \
    "$RAW/prep/score_ready.vcf.gz" \
    "$RAW/samplesheet.csv" \
    "$RAW/roundtrip_22chr/pgen/chr1.pgen" \
    > "$TESTDIR/upstream_before.sha256"

echo "UPSTREAM CHECKPOINT: PRESENT"
echo "EXPECTED NEXTFLOW FAILURE: CONFIRMED"

echo
echo "=== PHASE 3: RESUME USING REAL NEXTFLOW ==="

# The injected PATH is deliberately absent here.
python3 pipeline_v21.py batch \
    --config config.json \
    --manifest "$MANIFEST" \
    --run-id "$RUN_ID" \
    --jobs 1 \
    --resume \
    2>&1 | tee "$TESTDIR/resume.log"

echo
echo "=== PHASE 4: CHECK UPSTREAM INTEGRITY ==="

sha256sum \
    "$RAW/prep/score_ready.vcf.gz" \
    "$RAW/samplesheet.csv" \
    "$RAW/roundtrip_22chr/pgen/chr1.pgen" \
    > "$TESTDIR/upstream_after.sha256"

diff -u \
    "$TESTDIR/upstream_before.sha256" \
    "$TESTDIR/upstream_after.sha256"

grep -Fq "RAW_PGSC_CALC_RESUME" \
    "$TESTDIR/resume.log" || {
        echo "ERROR: Expected resume path was not observed."
        exit 1
    }

echo "UPSTREAM SHA256: UNCHANGED"

echo
echo "=== PHASE 5: VALIDATE RECOVERED PRS ==="

python3 - "$DEST" <<'PY'
import csv
import sys
from pathlib import Path

dest = Path(sys.argv[1])
scorefile = dest / "method2/scores.tsv"

expected = {
    "A": {
        "sum": 0.366424148908984,
        "percentile": 35.10204081632653,
        "observed": 609032,
        "missing": 1273,
    },
    "B": {
        "sum": 0.363472617655232,
        "percentile": 35.374149659863946,
        "observed": 602628,
        "missing": 810,
    },
}

assert scorefile.is_file()

with scorefile.open() as f:
    rows = list(csv.DictReader(f, delimiter="\t"))

assert len(rows) == 2

seen = set()

for r in rows:
    model = r["model"]
    assert model in expected
    assert model not in seen
    seen.add(model)

    x = expected[model]

    assert r["sample"] == "26T00970"
    assert abs(float(r["total_sum"]) - x["sum"]) <= 2e-6
    assert abs(float(r["eas_percentile"]) - x["percentile"]) <= 1e-6
    assert int(r["observed_loci"]) == x["observed"]
    assert int(r["missing_loci"]) == x["missing"]

    print(
        f"MODEL {model}: "
        f"SUM={r['total_sum']}, "
        f"percentile={r['eas_percentile']} PASS"
    )

assert seen == {"A", "B"}

print("RECOVERED PRS BENCHMARK: PASS (2/2)")
PY

echo
echo "POST-ROUNDTRIP FAILURE RECOVERY: PASS"
echo "Test logs retained at: $TESTDIR"
