# In-house allele frequency (in-house AF) — Phase 0

Build a cohort allele-frequency database from our own NovaSeq + DRAGEN WGS
gVCFs, and annotate `snv_indel.annotated.tsv` with `INHOUSE_AF` the same way
we annotate `GNOMAD_G_AF`. Goal: an in-house AF (incl. **rare** variants) for
annotation, and a filter for locally common polymorphisms that gnomAD
under-represents.

> **Scope (Phase 0/1):** SNV/indel only. SV in-house frequency (breakpoint /
> reciprocal-overlap matching) is a later phase. Mito is separate.

## Cohort

613 WGS samples (one DRAGEN gVCF each), selected from the raw datalake by
`select_cohort.py`:

- de-duplicated to one path per sample (same sample is staged under both
  `.../<run>/other/<sample>/` and `.../<run>/vcf.gz/`; we prefer the
  `vcf.gz` copy),
- excluded 1 broken/empty gVCF (`< 50 MB`),
- excluded 18 standard-reference controls (`VAL-37`..`VAL-54`).

No MRN/family de-duplication yet — a patient sequenced twice under different
ids, or related individuals, still count more than once. Acceptable for now
(small effect on *common* variants; revisit before using in-house AF as
population-frequency evidence). The manifest carries `sample_id` + `run` so
MRN/family columns can be bolted on later.

## Files

| file | what |
|---|---|
| `select_cohort.py` | gVCF path list → deduplicated, filtered cohort (`--out-list`) + audit manifest. Stdlib only. |
| `build_inhouse_af.sh` | cohort gVCFs → joint genotyping (GLnexus) → `inhouse_af.hg38.vcf.gz` (sites-only, `INHOUSE_{AC,AN,AF,NHOM}`). |
| `update_inhouse_af.sh` | per-batch SOP: append new samples (dedup), re-genotype the full cohort, atomic-swap the AF DB. Maintains `cohort_manifest.tsv` + `updates.log`. |

The cohort list / manifest / AF DB all contain patient sample ids and
datalake paths — **keep them out of git** (put them under
`$NGS_UI_HOME/biotools/inhouse_af/`, like the gnomAD/GeneBe DBs). Only these
scripts are committed.

## Phase 0 runbook (validate the toolchain first)

Run on the cluster head node where the datalake + reference live.

### 1. Select the cohort

```bash
# the size list is `du -h .../*.hard-filtered.gvcf.gz` over the run dirs
scripts/inhouse_af/select_cohort.py \
  --sizes hard_filtered_gvcf_sizes.txt \
  --exclude-range 'VAL-:37-54' \
  --out-manifest $NGS_UI_HOME/biotools/inhouse_af/cohort_manifest.tsv \
  --out-list     $NGS_UI_HOME/biotools/inhouse_af/cohort_gvcfs.txt
# → 613 included across 10 runs
```

### 2. Smoke test on ONE run (~64 samples) before the full cohort

```bash
# GLnexus — prefer the static binary (single file, zero deps, air-gap friendly):
#   download glnexus_cli from https://github.com/dnanexus-rnd/GLnexus/releases
export GLNEXUS_BIN=$NGS_UI_HOME/biotools/glnexus/glnexus_cli   # OR GLNEXUS_SIF=...
# bcftools: host binary on PATH is fine; or BCFTOOLS_SIF=<bcftools.sif> on DGM.

# restrict to one run by grepping the list, OR use --max-samples for a quick run
grep 20251118_LH00873_0004 $NGS_UI_HOME/biotools/inhouse_af/cohort_gvcfs.txt \
  > /tmp/run1_gvcfs.txt

scripts/inhouse_af/build_inhouse_af.sh \
  --list /tmp/run1_gvcfs.txt \
  --out  /tmp/inhouse_af.run1.vcf.gz \
  --threads 16 --mem-gbytes 96
```

The reference defaults to
`/home/datalake_Intermediate/pipeline/reference/hg38/Homo_sapiens_assembly38.fasta`
(DGM). Override with `--ref` when it lives elsewhere.

**What to check (this is the Phase 0 validation):**

1. **GLnexus accepts DRAGEN gVCFs with `--config gatk`** — no crash, plausible
   site count. If `gatk` over/under-filters DRAGEN calls, try `gatk_unfiltered`
   (then we rely on our own `view -f PASS` / thresholds). This config choice is
   the main open question Phase 0 resolves.
2. **AN denominator is right** — the verify block prints the AN distribution;
   at common sites it should peak near `2 × N_samples` (≈128 for one run).
   If AN is way below that, the gVCF reference blocks aren't being read
   (wrong config / gVCFs lack `<NON_REF>` blocks).
3. **Known common variant sanity** — pick a few well-known common SNPs and
   confirm `INHOUSE_AF` is in the right ballpark vs gnomAD EAS.

### 3. Demonstrate accumulation (the "every batch" requirement)

```bash
# run1 alone, then run1+run2 — AC should grow, AN ≈ 2×N, AF stays sane
grep -E '20251118_LH00873_0004|20251121_LH00873_0005' \
  $NGS_UI_HOME/biotools/inhouse_af/cohort_gvcfs.txt > /tmp/run12_gvcfs.txt
scripts/inhouse_af/build_inhouse_af.sh --list /tmp/run12_gvcfs.txt \
  --out /tmp/inhouse_af.run12.vcf.gz --threads 16 --mem-gbytes 96
```

### 4. Full cohort

**Restrict to primary contigs.** DRAGEN gVCFs carry ~3340 decoy/alt/HLA
contigs; variants there are useless for annotating `snv_indel.annotated.tsv`
(primary assembly only) and processing them is slow. Build a one-time BED of
chr1-22, X, Y, M from the reference `.fai` and pass it with `--bed`:

```bash
REF=/home/datalake_Intermediate/pipeline/reference/hg38/Homo_sapiens_assembly38.fasta
awk 'BEGIN{OFS="\t"} $1 ~ /^chr([0-9]+|X|Y|M)$/ {print $1,0,$2}' "$REF.fai" \
  > $NGS_UI_HOME/biotools/inhouse_af/primary_contigs.bed

scripts/inhouse_af/build_inhouse_af.sh \
  --list $NGS_UI_HOME/biotools/inhouse_af/cohort_gvcfs.txt \
  --out  $NGS_UI_HOME/biotools/inhouse_af/inhouse_af.hg38.vcf.gz \
  --bed  $NGS_UI_HOME/biotools/inhouse_af/primary_contigs.bed \
  --threads 32 --mem-gbytes 192
```

**Timing note.** GLnexus bulk-loads every gVCF into a scratch DB before
genotyping, and that load is data-volume bound (~2.5 h / 64 WGS observed →
roughly a day for 600+). `--bed` speeds the discover/genotype tail but not the
bulk load. If the full-rebuild-per-batch time becomes painful at scale, that is
the signal to switch to the DRAGEN iterative gVCF Genotyper (see below).

## Running on DGM now, air-gapped DGX later

Both machines can read the reference and the gVCFs; the only difference is
that **on the DGX the datalake paths drop the leading `/home`** (reference at
`/datalake_Intermediate/.../hg38/...`, gVCFs at `/datalake_Raw/Novaseq/...`).

Every external tool is pluggable, so the same script runs in both places:

| | DGM (now) | air-gapped DGX (later) |
|---|---|---|
| GLnexus | `GLNEXUS_BIN` (static) or `GLNEXUS_SIF` | `GLNEXUS_BIN` static binary (no network) |
| bcftools | host `bcftools` or `BCFTOOLS_SIF` | host `bcftools` binary |
| reference | `--ref /home/datalake_Intermediate/.../hg38/Homo_sapiens_assembly38.fasta` | `--ref /datalake_Intermediate/.../hg38/Homo_sapiens_assembly38.fasta` |
| cohort list | absolute `/home/...` paths | reuse same list + `--strip-path-prefix /home` |

To move to the DGX, stage offline once: the `glnexus_cli` binary, a `bcftools`
binary, and (already present) the reference + gVCFs. GLnexus is **CPU-only**
(the DGX GPUs are irrelevant) — it just needs many cores, RAM, and scratch
disk. Example DGX invocation:

```bash
GLNEXUS_BIN=/opt/glnexus/glnexus_cli BCFTOOLS_BIN=bcftools \
scripts/inhouse_af/build_inhouse_af.sh \
  --list cohort_gvcfs.txt --strip-path-prefix /home \
  --ref /datalake_Intermediate/pipeline/reference/hg38/Homo_sapiens_assembly38.fasta \
  --out /path/inhouse_af.hg38.vcf.gz --threads 32 --mem-gbytes 192
```

When a `*_SIF` is set the script apptainer-execs it (binding `APPTAINER_BIND`,
default `/home`); otherwise it runs the host binary with no container at all.

## Incremental update design (each new run)

Two options; Phase 0 uses **A** as the always-correct baseline.

**A. Full re-genotype every batch (baseline, what these scripts do).**
Append the new run's gVCFs to `cohort_gvcfs.txt`, re-run `build_inhouse_af.sh`
over the whole cohort, atomic-swap the output. Simple and always correct.
GLnexus scales to thousands of WGS, so for ~600→growing this is acceptable
(hours, run off-peak). The pitfall it avoids: you cannot just genotype the new
batch alone and add counts, because a site that is variant only in old batches
has no record in the new batch's joint VCF — its AN contribution from the new
(hom-ref) samples would be lost.

**B. DRAGEN iterative gVCF Genotyper (true incremental, later).**
DRAGEN's native iterative genotyper maintains a cohort "census" and folds in
new gVCFs without re-processing the old ones, back-filling AN at known sites
from each gVCF's reference blocks. This is the right long-term engine if
re-genotyping time becomes a problem, but it needs DRAGEN hardware/license
time and its own validation. Keep the same `select_cohort.py` front-end and
the same `INHOUSE_*` sites-VCF output contract so the UI side doesn't change.

Either way, maintain `cohort_manifest.tsv` as the audit record of which
samples (and runs) are in the current DB, and **dedup on append** so a re-run
batch is never counted twice.

## Phase 1: each new batch (after Phase 0 validates the config)

Once `--config` is confirmed, you never call `build_inhouse_af.sh` by hand
again — `update_inhouse_af.sh` is the single per-batch entry point:

```bash
COHORT=$NGS_UI_HOME/biotools/inhouse_af

# new run lands → list its gVCFs (exact bytes, tab-separated)
find /home/datalake_Raw/Novaseq/<NEW_RUN> -name '*.hard-filtered.gvcf.gz' \
  -printf '%s\t%p\n' > /tmp/new_batch_sizes.txt

scripts/inhouse_af/update_inhouse_af.sh \
  --cohort-dir "$COHORT" \
  --new-sizes  /tmp/new_batch_sizes.txt \
  --exclude-range 'VAL-:37-54' \
  --threads 32 --mem-gbytes 192
# → appends only unseen samples, re-genotypes the whole cohort, atomic-swaps
#   $COHORT/inhouse_af.hg38.vcf.gz, logs to $COHORT/updates.log
```

Re-running the same batch is safe (dedup by `sample_id` → "added=0, no
rebuild"). Use `--manifest-only` to update bookkeeping without rebuilding, and
`--strip-path-prefix /home` on the DGX.

## Integrating into NGS-UI (Phase 2, not done yet)

Mirror the gnomAD pattern:

1. `config.py`: add `NGS_UI_INHOUSE_AF_DB` (default
   `$NGS_UI_HOME/biotools/inhouse_af/inhouse_af.hg38.vcf.gz`).
2. `run_stopgaps.sh`: after the gnomAD annotate step, add
   `bcftools annotate -a $INHOUSE_AF_DB -c INFO/INHOUSE_AC,INFO/INHOUSE_AN,INFO/INHOUSE_AF,INFO/INHOUSE_NHOM`
   onto `snv_indel.annotated.tsv`. Re-annotating rewrites the whole TSV, so
   rebuild `snv_gene_index.sqlite` too (same gotcha as `backfill_giab_strata.sh`).
3. `adapters/snv_tsv.py`: surface `inhouse_af` / `inhouse_ac` / `inhouse_an`
   next to the gnomAD fields.
4. Front-end: an in-house badge on the SNV card showing **AC/AN** (not just %,
   so N=613 sample-size is visible) and an "In-house AF < x" display filter
   next to the existing `gnomAD_G_AF < 0.01` toggle.
5. A `backfill_inhouse_af.sh` (cf. `backfill_giab_strata.sh`) for existing
   samples.

## Caveats to keep in mind

- **Disease-referral cohort, not population controls.** Use in-house AF to
  flag recurrent artifacts / locally common variants; do **not** use it as
  ACMG BA1/BS1 population-frequency evidence (that stays gnomAD).
- **N=613** → minimum resolvable allele frequency ≈ 1/1226 ≈ 0.08%; show AC/AN.
- **No relatedness/MRN dedup yet** (see Cohort note).
