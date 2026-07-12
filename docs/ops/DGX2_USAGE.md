# DGX-2 NGS Secondary Analysis Pipeline — Usage Guide

## System Overview

| Item | Details |
|------|---------|
| Host | DGX-2 (`10.11.33.75`) — accessible from internal hospital network only |
| GPUs | 6 × NVIDIA V100 32GB (GPU indices 10–15) |
| RAM | 1.5 TB |
| Analysis output | `/datalake_Intermediate/pipeline/nextflow_output/` |
| Raw FASTQ | `/datalake_Raw/` |

---

## Step 1 — Prepare the Samplesheet

```bash
# Set the batch name (recommended format: date_type_batch)
BATCH_NAME="20260425_WES_BATCH1"
OUT_DIR="/datalake_Intermediate/pipeline/nextflow_output/${BATCH_NAME}"
LAUNCH_DIR="/datalake_Intermediate/pipeline/nextflow_launch/${BATCH_NAME}"
WORK_DIR="/raid/DGM/work/${BATCH_NAME}"

mkdir -p "${LAUNCH_DIR}" "${WORK_DIR}"
nano ${OUT_DIR}/samplesheet.csv
```

Samplesheet format (`samplesheet.csv`):

```csv
sample,fastq_1,fastq_2,sex
SAMPLE001,/datalake_Raw/.../SAMPLE001_R1.fastq.gz,/datalake_Raw/.../SAMPLE001_R2.fastq.gz,female
SAMPLE002,/datalake_Raw/.../SAMPLE002_R1.fastq.gz,/datalake_Raw/.../SAMPLE002_R2.fastq.gz,male
```

**Notes:**
- For multi-lane samples, add a `lane` column (e.g., `L001`, `L002`). FASTP will QC each lane separately; FQ2BAM will merge them automatically.
- Without a `lane` column, every `sample` value must be unique.
- `sex`: use `male`, `female`, or `unknown`.
- All paths must be absolute paths as seen from DGX-2.

### WGS FASTQ rule

- WGS must use the original lane pairs named like `SAMPLE_S1_L001_R1_001.fastq.gz` / `R2`; do not use `SAMPLE_R1_merged.fastq.gz` / `R2_merged.fastq.gz` even when both forms exist in the run folder.
- The NGS-UI picker shows one option per main sample and reports its lane count. Creating the samplesheet expands that option into one row per lane, repeats the same `sample` value, and writes the `lane` column.
- For example, a sample with `L001` through `L008` appears once in the picker but produces eight samplesheet rows (16 FASTQ paths).

---

## Step 2 — Load the Environment

```bash
source /datalake_Intermediate/pipeline/pipeline_code/NGS2ndAnalysis_env.sh
```

On success, the script prints environment info and a usage example. It also automatically clears any stale GPU lock files from previous runs.

---

## Step 3 — Run the Pipeline

```bash
# Always run from the work directory
cd "${LAUNCH_DIR}"
```

### WES (with gCNV CNV calling)

```bash
nextflow -c ${PIPELINE_CONFIG} run ${PIPELINE_CODE}/main.nf \
    -profile dgx \
    --input_csv "${OUT_DIR}/samplesheet.csv" \
    --seq_type WES \
    --run_gcnv true \
    --out_dir "${OUT_DIR}" \
    -w "${WORK_DIR}" \
    -resume
```

### WGS

```bash
nextflow -c ${PIPELINE_CONFIG} run ${PIPELINE_CODE}/main.nf \
    -profile dgx \
    --input_csv ${OUT_DIR}/samplesheet.csv \
    --seq_type WGS \
    --out_dir ${OUT_DIR} \
    -w "${WORK_DIR}" \
    -resume
```

### Single-sample accelerated mode (`dgx_single`)

When running only one sample, use `-profile dgx_single` to allocate all 6 GPUs to that sample. Alignment completes in ~40 minutes instead of ~3 hours.

> **Note:** If multiple samples are being processed in the same batch, the standard `-profile dgx` (one GPU per sample, up to 6 in parallel) gives better overall throughput.

```bash
nextflow -c ${PIPELINE_CONFIG} run ${PIPELINE_CODE}/main.nf \
    -profile dgx_single \
    --input_csv ${OUT_DIR}/samplesheet.csv \
    --seq_type WGS \
    --out_dir ${OUT_DIR} \
    -resume
```

---

## Step 4 — Check Output

After the pipeline completes, each sample produces the following directory structure:

```
/datalake_Intermediate/pipeline/nextflow_output/<BATCH_NAME>/<SAMPLE>/
├── 01_preprocessing/     ← FASTP QC reports
├── 02_alignment/         ← BAM files
├── 03_alignment_qc/      ← Mosdepth coverage reports
├── 04_snv_indel/         ← SNV/Indel VCFs (DeepVariant + HaplotypeCaller + Ensemble)
├── 05_cnv_sv/            ← CNV (gCNV/CNVkit) + SV (Manta)
├── 06_repeat/            ← STR (ExpansionHunter)
├── 07_mitochondria/      ← mtDNA variants
├── 08_roh/               ← ROH analysis (AutoMap)
├── 09_postprocessing/    ← bcftools stats
└── pipeline_info/        ← Execution reports (HTML + timeline)
```

---

## GPU Management

The pipeline includes an automatic GPU lock mechanism (`gpu_lock.sh` / `gpu_unlock.sh`) that uses `flock` to atomically assign GPUs 10–15 to processes. **You do not need to specify GPUs manually.**

Lock files are stored in `/tmp/nxf_gpu_locks/`. Stale locks from crashed runs are automatically cleared when you `source NGS2ndAnalysis_env.sh`.

---

## Disk Cleanup

Intermediate files accumulate in `/raid/DGM/work/`. After confirming outputs are correct, clean up with:

```bash
cd /raid/DGM/work
nextflow clean -f
```

---

## Resuming an Interrupted Run

If the pipeline is interrupted for any reason, re-run the original command with `-resume`. Nextflow will skip all already-completed steps as long as `/raid/DGM/work/` is intact.

---

## License and Acknowledgements

This pipeline was developed for the Department of Genomic Medicine and Neurology, National Cheng Kung University Hospital.

### Pipeline Code

Copyright © 2026 pylin1991. All rights reserved.  
This code is proprietary and intended for internal hospital use only.  
Redistribution or use outside of NCKUH without explicit written permission is prohibited.

### Third-party Tools

This pipeline orchestrates the following third-party tools, each governed by their respective licenses:

| Tool | Version | License |
|------|---------|---------|
| [Nextflow](https://github.com/nextflow-io/nextflow) | ≥ 23.x | Apache License 2.0 |
| [Apptainer](https://github.com/apptainer/apptainer) | ≥ 1.x | BSD 3-Clause License |
| [NVIDIA Clara Parabricks](https://docs.nvidia.com/clara/parabricks/latest/index.html) | 4.4.0 / 4.7.0 | [NVIDIA AI Product Agreement](https://docs.nvidia.com/clara/parabricks/latest/license.html) — free to use; enterprise support available via NVIDIA AI Enterprise |
| [GATK](https://github.com/broadinstitute/gatk) | 4.6.2.0 | Apache License 2.0 |
| [fastp](https://github.com/OpenGene/fastp) | 1.3.0 | MIT License |
| [SAMtools](https://github.com/samtools/samtools) | 1.23.1 | MIT License |
| [BCFtools](https://github.com/samtools/bcftools) | 1.23.1 | MIT License (default; GPL v3 if compiled with GNU Scientific Library) |
| [Mosdepth](https://github.com/brentp/mosdepth) | 0.3.13 | MIT License |
| [Manta](https://github.com/Illumina/manta) | 1.6.0 | GNU GPL v3 |
| [CNVkit](https://github.com/etal/cnvkit) | 0.9.12 | Apache License 2.0 |
| [ExpansionHunter](https://github.com/Illumina/ExpansionHunter) | 5.0.0 | GNU GPL v3 |
| [BWA](https://github.com/lh3/bwa) | 0.7.19 | GNU GPL v3 |
| [MultiQC](https://github.com/MultiQC/MultiQC) | 1.33 | GNU GPL v3 |
| [AutoMap](https://github.com/mquinodo/AutoMap) | 1.3 | MIT License |

### Reference Data

Reference files used by this pipeline are sourced from the [Broad Institute Google Cloud Storage bucket](https://console.cloud.google.com/storage/browser/gcp-public-data--broad-references) and are subject to the [Broad Institute data use terms](https://software.broadinstitute.org/gatk/download/bundle).
