#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 3 ]]; then
  echo "usage: $0 HAPLOTYPECALLER_VCF WORK_DIR SAMPLE_ID" >&2
  exit 2
fi

VCF=$1
WORK_DIR=$2
SAMPLE_ID=$3
DEFAULT_SIF=/home/datalake_Intermediate/pipeline/nextflow_containers/automap_1.3.sif
FALLBACK_SIF=/home/pipeline/nextflow_containers/automap_1.3.sif
AUTOMAP_SIF=${NGS_UI_AUTOMAP_SIF:-$DEFAULT_SIF}

if [[ ! -f "$AUTOMAP_SIF" && -f "$FALLBACK_SIF" && -z "${NGS_UI_AUTOMAP_SIF:-}" ]]; then
  AUTOMAP_SIF=$FALLBACK_SIF
fi
[[ -f "$VCF" ]] || { echo "HaplotypeCaller VCF not found: $VCF" >&2; exit 2; }
[[ -f "$AUTOMAP_SIF" ]] || { echo "AutoMap image not found: $AUTOMAP_SIF" >&2; exit 2; }
command -v apptainer >/dev/null || { echo "apptainer not found" >&2; exit 2; }

mkdir -p "$WORK_DIR"
BIND_ARGS=(--bind /home)
if [[ "$WORK_DIR" != /home/* ]]; then
  BIND_ARGS+=(--bind "$WORK_DIR")
fi
VCF_DIR=$(dirname "$VCF")
if [[ "$VCF_DIR" != /home/* && "$VCF_DIR" != "$WORK_DIR" ]]; then
  BIND_ARGS+=(--bind "$VCF_DIR")
fi

# AutoMap writes into its installation directory while unpacking resources, so
# use a writable copy inside this job-private staging directory.  Arguments are
# passed positionally to the container shell rather than interpolated into it.
apptainer exec \
  "${BIND_ARGS[@]}" \
  "$AUTOMAP_SIF" \
  bash -c '
    set -euo pipefail
    work=$1
    vcf=$2
    sample=$3
    cd "$work"
    cp -r /opt/AutoMap ./AutoMap_local
    bcftools view "$vcf" -O v -o input.vcf
    bash ./AutoMap_local/AutoMap_v1.3.sh \
      --vcf input.vcf \
      --out . \
      --genome hg38 \
      --id "$sample" \
      --chrX \
      --minsize 1.0
    test -f "$sample/$sample.HomRegions.tsv"
    mv "$sample/$sample.HomRegions.tsv" "$work/"
    if test -f "$sample/$sample.HomRegions.pdf"; then
      mv "$sample/$sample.HomRegions.pdf" "$work/"
    fi
  ' bash "$WORK_DIR" "$VCF" "$SAMPLE_ID"
