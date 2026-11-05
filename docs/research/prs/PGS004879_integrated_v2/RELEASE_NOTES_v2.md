# PGS004879 integrated v2.0 release notes — 2026-10-31

## What changed

- **Added `pipeline_v2.py`**: a single CLI for new raw gVCF → unchanged legacy raw reconstruction → Method 2 Model A/B. Optional `--with-legacy-ancestry` also invokes unchanged v1 ancestry/QC; default mode skips its target-dependent reference matching.
- **Added `method2_score.py`**: frozen 610,305-locus Model A and 603,438-locus Model B scoring, allele/GT crosswalk, separate observed/imputed contributions, unrelated 3,330 AF imputation, 735 EAS ECDF, boundary flags, high-weight missing variants and reports.
- **Added model freeze**: source PGS Catalog content hash, existing 17-file audit inventory, Model A/B scoring-file SHA256, source-to-canonical provenance issues, per-variant master-vs-scorefile check, `.afreq` effect allele frequency consistency and 3,330/735 IID checks. New server lock protects model files on every execution.
- **Added three-target regression**: directly consume the prevalidated `pilot_fresh_02` raw audits and compare six target SUM/percentile benchmarks; does not write to original three pilot runs.
- **Added safeguards**: no QC-failed genotype as 0/0; checks actual VCF genotype/alleles vs canonical effect dosage; checks accepted noncanonical raw contributions so they reconcile to legacy raw SUM; sample-lock, input+script hashes, fail on stale resources, warning on unknown ancestry.
- **Added synthetic tests** and new integration SOP.

## Files explicitly unchanged from supplied integrated_v1

`pipeline.py`, `vendor/PGS004879_FINAL_gVCF_resume.sh`, `vendor/PGS004879_FINAL_gVCF_ORIGINAL.sh`, `qc_report.py`, `compare_benchmarks.py`, `benchmarks.json` and the v1 test module. They are retained to preserve previously reviewed scientific logic and golden regression baselines.

## Validation boundary

The package has passed author-side Python compilation, `bash -n`, and **12 local synthetic unit tests**; all three real-case raw audits and reference model assets are on the user's server and are **not included in the uploaded ZIP**. Thus the new v2 entrypoint is **not yet server regression validated**. Run `freeze-model`, then `regression` (expected `6/6 PASS`), then a new indexed gVCF fresh run. Do not replace the prevalidated v1 production/golden package before server validation.

**Not included:** 50,000-sample batch execution orchestration, Taiwan-specific percentile/CHD outcome calibration, age-/sex-specific absolute-risk conversion. These remain separate downstream milestones.
