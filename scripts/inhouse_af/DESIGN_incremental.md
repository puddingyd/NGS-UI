# In-house AF — incremental architecture design

Status: **design only, not implemented.** This document specifies the
incremental ("每批下機只算新批") production design. The current working
pipeline (`build_inhouse_af.sh`, full GLnexus joint genotyping) stays as the
periodic ground-truth / reconciliation engine.

---

## 1. Locked decisions

| # | Decision | Value |
|---|---|---|
| 1 | Callable depth threshold | **DP ≥ 10** |
| 2 | Chromosomes | **autosomes + X + Y + M, all from v1** (sex-aware ploidy) |
| 3 | AC (numerator) source | **per-sample carrier counting** (truly per-sample incremental) |
| 4 | Sex per sample | **inferred from gVCF coverage** (portable; covers unregistered samples) |
| 5 | Ground-truth check | **quarterly full GLnexus rebuild**, compared against the incremental DB |
| 6 | Scope | SNV/indel (incl. chrM). SV is a separate later phase. |

Why per-sample counting is correct here: the denominator (AN) no longer comes
from the variant calls — it comes from a coverage-derived **AN track**. So the
classic "naive stacking under-counts AN" problem is gone. Per-sample counting
then differs from full joint genotyping only by a few borderline calls, which
the quarterly GLnexus reconciliation catches. In return we get true per-sample
incrementality, the simplest moving parts, and CPU-only portability (runs on the
air-gapped DGX, unlike DRAGEN's FPGA-bound iterative genotyper).

This is the same shape as gnomAD: allele **counts** from genotypes + an
**allele-number / coverage track** for the denominator.

---

## 2. Data stores

```
$NGS_UI_HOME/biotools/inhouse_af/
  counts.sqlite              # per normalized site: n_hom / n_het / n_hemi
  an_track.bed.gz (+ .tbi)   # genome-wide AN = Σ ploidy-weight of callable samples
  samples_manifest.tsv       # ingested samples (dedup, sex call, QC)
  inhouse_af.hg38.vcf.gz      # PUBLISHED sites VCF (the annotation DB; atomic-swapped)
  per_sample/{sample_id}/     # audit + re-do material
      callable.weighted.bed.gz   # DP>=10 intervals, ploidy-weighted
      counts.tsv                 # this sample's variant contributions
      qc.json                    # sex call, x/y ratios, callable bp, n_variants
  reconcile/                  # quarterly GLnexus snapshots + comparison reports
```

Everything except this repo's scripts is patient-derived → stays out of git,
under `biotools/inhouse_af/` like the gnomAD/GeneBe DBs.

---

## 3. Output contract (what the UI eventually consumes)

A sites-only VCF, identical in spirit to the gnomAD AF VCF, with INFO:

| INFO | meaning |
|---|---|
| `INHOUSE_AC` | alt allele count = `2·n_hom + n_het + n_hemi` |
| `INHOUSE_AN` | total called alleles at the site (from the AN track) |
| `INHOUSE_AF` | `INHOUSE_AC / INHOUSE_AN` |
| `INHOUSE_NHOM` | number of homozygous-alt **individuals** (`n_hom`) |
| `INHOUSE_HEMI` | number of hemizygous-alt individuals (`n_hemi`; males on non-PAR X/Y, chrM) |

> Contract-alignment note: the GLnexus path currently emits `INHOUSE_AC_HOM`
> (= alleles in hom genotypes = `2·n_hom`). We standardize the canonical contract
> on `INHOUSE_NHOM` (individuals, clinician-friendly) and will have the GLnexus
> path divide by 2 so both engines emit identical fields. The UI integration
> (separate branch, later) only ever sees this contract, so it is agnostic to
> which engine produced the DB.

---

## 3a. Per-sample data sources (confirmed on real DRAGEN output)

Anchor each sample on its **`other/{sample}/` directory** and drive everything
from the **indexed gVCF there** — it contains both the reference blocks (for AN)
and the PASS variant records (for AC), so one file covers both:

| need | file | notes |
|---|---|---|
| sex / karyotype | `{sample}.ploidy_estimation_metrics.csv` | `Ploidy estimation,XX/XY/…` (4.1) |
| **AC + callable/AN** | `{sample}.hard-filtered.gvcf.gz` (+ `.tbi`) | **indexed** here. Ref blocks carry `MIN_DP` (→ callable, 4.2); variant rows carry `FILTER` + `GT` (→ AC, 4.4). chrM variants are present too (4.6). |
| (future) ROH / STR | `{sample}.roh.bed`, `{sample}.repeats.vcf.gz` | unrelated to AF; material for the UI's ROH/STR cards |

The separate `{sample}.hard-filtered.vcf.gz` (variant-only) lives in the
`vcf.gz/` delivery dir and is **not indexed there**; its PASS calls are the same
ones already in the gVCF, so we don't need it. Datalake files are read-only →
ingest never writes indexes; it uses the existing `.tbi` in `other/` and streams
the gVCF for the genome-wide callable pass.

## 4. Components

### 4.1 Sex / karyotype (from DRAGEN, not inferred)

**Primary source: DRAGEN's own ploidy estimation.** Each sample's
`other/{sample}/{sample}.ploidy_estimation_metrics.csv` ends with:

```
PLOIDY ESTIMATION,,Ploidy estimation,XX
```

This is DRAGEN's final karyotype call (XX / XY / X0 / XXY / …), computed from
**median** per-contig coverage — robust, aneuploidy-aware, and free. We parse
that one field. No need to roll our own inference.

> Why not infer it ourselves: median beats mean badly here. On a confirmed XX
> sample, DRAGEN's `Y median / Autosomal = 0.00` cleanly says female, but a
> *mean* MIN_DP over non-PAR Y gives ~4× (Ry≈0.14) because of Y/X homology
> hotspots — right on a naive 0.1 threshold. DRAGEN already solved this.

**Fallback (older samples missing the CSV):** mean MIN_DP coverage ratios
(Rx/Ry) from the gVCF, with conservative "ambiguous" handling (4.2). Flag in
`qc.json` which source was used.

Karyotype → ploidy weighting in 4.2 (`XX`→female column, `XY`→male column,
anything else → ambiguous/excluded for X/Y, flagged for review).

### 4.2 Callable BED (DP ≥ 10) + ploidy weighting

From the gVCF: emit an interval wherever the sample is callable:
- reference blocks with `MIN_DP ≥ 10` → whole block callable
- variant records with `DP ≥ 10` → that position callable

Then assign a **ploidy weight** per interval by region and sex:

| Region | XX (female) | XY (male) | ambiguous |
|---|---|---|---|
| autosomes (chr1–22) | 2 | 2 | 2 |
| chrX PAR1 / PAR2 | 2 | 2 | 0 (X excluded) |
| chrX non-PAR | 2 | 1 | 0 |
| chrY non-PAR | 0 (excluded) | 1 | 0 |
| chrM | 1 | 1 | 1 |

hg38 PAR coordinates (constants):
- PAR1: `chrX:10,001–2,781,479`, `chrY:10,001–2,781,479`
- PAR2: `chrX:155,701,383–156,030,895`, `chrY:56,887,903–57,217,415`

Output `per_sample/{id}/callable.weighted.bed.gz`: `chrom start end weight`.

### 4.3 AN track (the denominator, accumulated)

`an_track` is a genome-wide bedGraph: for every interval, the value = **sum of
ploidy weights of all callable samples** there. That value **is** the AN at any
position in the interval.

- **Accumulate** (per new sample): `an_track ⊕ sample.callable.weighted.bed`
  (interval add: where they overlap, AN += sample weight). Implemented as a
  sort + sweep-line merge of breakpoints; the track is run-length so size stays
  bounded (it is a coverage track, long flat runs).
- **Store** bgzipped + tabix-indexed → O(log n) AN lookup at any site on publish.
- Periodic recompaction keeps breakpoints tidy as the cohort grows.

Because AN comes from coverage (not from who-carries-a-variant), a site that is
polymorphic in only one batch still gets the full denominator from every other
callable sample — this is the whole reason the design is correct **and**
incremental.

### 4.4 Variant counts (the numerator, accumulated) — autosomes + X/Y

**Source: PASS variant records in the indexed gVCF** (`FILTER=PASS`, real ALT,
drop the `<NON_REF>` symbolic allele). We only need `FILTER` + `GT`, both always
present, so the gVCF's variable variant FORMAT (some rows are `GT:GQ:PS`) is a
non-issue. chrM is handled separately in 4.6.

Per new sample: normalize its PASS variant calls and classify each ALT-bearing
record by genotype, then UPSERT-add into `counts.sqlite`:

- diploid `1/1` → `n_hom += 1`
- diploid `0/1` / `1/2`-split het → `n_het += 1`
- haploid alt (male non-PAR X, male Y, chrM) → `n_hemi += 1`

Normalization: `bcftools norm -m- -f ref --check-ref x` (same convention as the
GLnexus path, so keys match the annotation TSV and the two engines agree).
Skip a sample's X/Y records when sex is ambiguous (matches 4.2).

`counts.sqlite` schema (additive, idempotent per sample via the manifest):

```
variant_counts(
  chrom TEXT, pos INTEGER, ref TEXT, alt TEXT,   -- normalized; PRIMARY KEY
  n_hom INTEGER, n_het INTEGER, n_hemi INTEGER
)
```

### 4.5 Publish (counts + AN track → sites VCF)

For each row in `counts.sqlite`:
- `AC = 2·n_hom + n_het + n_hemi`
- `AN = an_track lookup at (chrom,pos)`
- `AF = AC / AN` (guard `AN > 0`)
- `NHOM = n_hom`, `HEMI = n_hemi`

Emit sites VCF → bgzip + tabix → **atomic swap** into
`inhouse_af.hg38.vcf.gz`. Same atomic-install discipline as `deploy_genebe_db.sh`.

### 4.6 chrM (mitochondrial) — haploid carrier frequency

chrM needs its own counting rule, confirmed against real DRAGEN output:

- chrM is in the same gVCF, with PASS variant rows and high depth (DP ~1000–3000).
- **GT is written diploid-style** (`1/1` ≈ homoplasmic, `0/1` ≈ heteroplasmic) —
  do **not** apply the autosome `2·hom + het` arithmetic.
- **`FORMAT/AF` is the heteroplasmy fraction** (e.g. `1`, `0.9994`, `0.0541`,
  `0.1963`).

Counting (one mito genome per sample):
- `AN_chrM` at a site = number of **callable** samples on chrM (MIN_DP ≥ 10,
  weight **1**) — chrM is effectively always callable given the coverage.
- a sample is a **carrier** of a chrM ALT if it has a **PASS** record there with
  a real ALT in GT; count **1** per carrier (never 2, even for `1/1`).
- `AF_chrM = carriers / AN_chrM`.
- store the **heteroplasmy** alongside: e.g. `n_carrier`, plus
  `n_homoplasmic` (AF ≥ 0.95) vs `n_heteroplasmic` (AF < 0.95), so the UI can
  later distinguish "homoplasmic in N" from "low-level heteroplasmy in M".

**v1 carrier threshold = PASS only** (DRAGEN's filter already gates quality; we
keep even low-heteroplasmy PASS calls but record their level). A minimum
heteroplasmy cutoff can be added later from the stored stats without re-ingest.

Output fields for chrM rows reuse the contract: `INHOUSE_AC` = carriers,
`INHOUSE_AN` = callable mito genomes, `INHOUSE_AF` = carrier frequency,
`INHOUSE_NHOM` = `n_homoplasmic`, plus an extra `INHOUSE_HET_MT` = `n_heteroplasmic`.

---

## 5. Incremental flow (each new run)

```
new run gVCFs
   └─ select_cohort.py (existing) → new, non-duplicate samples
        └─ for each sample (parallel, ~minutes each):
             1. infer sex (4.1)            → qc.json
             2. callable weighted BED (4.2)→ per_sample/{id}/
             3. add to AN track (4.3)
             4. count variants (4.4)       → counts.sqlite
             5. append manifest (dedup)
        └─ publish (4.5) → atomic-swap inhouse_af.hg38.vcf.gz
```

- **Per-sample, embarrassingly parallel.** Cost is O(new samples), independent of
  cohort size — 600 or 6000, each batch is the same few-hours-for-~64.
- **Idempotent.** A sample already in the manifest is skipped, so re-running a
  batch never double-counts (AN-track and counts both gated by the manifest).
- **Crash-safe.** AN-track and counts updates are staged then committed per
  sample; a half-ingested sample is rolled back and retried.

---

## 5a. Stratified AF by disease category (cancer / neuro / healthy / …)

Fully supported by this architecture — it is exactly the gnomAD "AF per
population" model, and a **late-binding** one: the per-sample artifacts
(`per_sample/{id}/` callable BED + counts) are **category-agnostic**, so disease
labels can be added/changed **without re-processing any gVCF**. Stratification
lives entirely in the accumulate/publish layer.

To add it:
1. A `sample_category(sample_id, category)` table (or TSV) — many-to-one (a
   sample can also belong to several groups, e.g. "cancer" + "all").
2. **Counts per category:** add `category` to the variant_counts key
   (`PRIMARY KEY(chrom,pos,ref,alt,category)`), accumulating each sample into its
   category bucket (plus an `ALL` bucket).
3. **One AN track per category:** the AN track is just "sum of weights of the
   callable samples in that group", so build `an_track.<category>.bg.gz` from
   only that category's BEDs — same event/delta machinery, partitioned by label.
4. **Publish:** emit overall `INHOUSE_AF` plus per-category
   `INHOUSE_AF_cancer`, `INHOUSE_AN_cancer`, … (gnomAD-style suffixes).

Because labels bind at accumulate time, you can (re)derive a fully stratified DB
from the existing `per_sample/` cache by re-running accumulate with the
`sample_category` map — no gVCF re-processing, no re-ingest. Caveats: small
strata → noisy AF (show AC/AN); a "healthy" stratum is the most useful control
but hardest to define in a referral cohort. Recommend keeping `ALL` as the
default and adding strata once labels exist. (Not built yet — flagged as the
next extension after Phase C validates the un-stratified path.)

## 6. Reconciliation with full GLnexus (quarterly)

1. Full `build_inhouse_af.sh` over the whole cohort → a GLnexus snapshot.
2. Compare incremental `INHOUSE_AF` vs GLnexus `INHOUSE_AF` at shared sites:
   correlation, AF scatter, and a list of largest discrepancies.
3. Expected, benign differences: borderline genotype rescue (joint vs
   per-sample), complex-indel representation, AN coverage-proxy vs exact joint
   call. Large/systematic drift → investigate.

Production DB = the incremental one. GLnexus = audit/ground-truth. (Optionally,
adopt the GLnexus snapshot as the canonical DB each quarter and continue
incrementing from there — decide after we see the first comparison.)

### 6a. Phase C validation result (64-sample run, on the DGX)

The incremental DB and a full GLnexus joint-genotyping of the **same 64 WGS**
were compared by `INHOUSE_AF` at shared, normalized sites. On sites where both
methods well-call the locus (`AN ≥ 120` of 128):

| variant class | n | Pearson r |
|---|---|---|
| SNP | 9.38 M | **0.9998** |
| indel | 1.86 M | **0.9962** |

→ the incremental AC/AN/AF machinery is **correct** (near-perfect agreement
with the joint-genotyping ground truth). The *overall* correlation looks low
(~0.5) only because ~14% of sites fall in a low-`AN` tail where GLnexus
**no-calls** low-confidence samples (smaller denominator → higher GLnexus AF),
while the incremental counts every DP≥10-callable sample. That is a
**methodological difference, not a bug** — and "count all callable samples" is
the more complete choice for an in-house AF. Repeat/STR indels additionally
differ because per-sample DRAGEN calls stutter there (incremental may
over-count); those are difficult regions flagged elsewhere anyway.

This validation also caught and fixed a real bug: indel callable intervals must
be a **single anchor base** (§4.2) and the AN-track lookup must use `pos-1`
(0-based) — otherwise a deletion's multi-base REF span overlaps the same
sample's adjacent ref block and double-counts AN.

### 6b. Phase C validation result (full cohort: 677 incremental vs 675 GLnexus)

Repeated at production scale (`compare_inhouse_af.py`). Incremental = 677
samples (43.5 M sites); GLnexus baseline = 675 (2 malformed gVCFs dropped, see
§7 / `known_bad_gvcfs.txt`), 25.59 M sites. Pearson r of `INHOUSE_AF` at shared
normalized sites, by AN floor (2×675 = 1350):

| class | AN≥1300 (≥96%) | AN≥1200 | AN≥1000 | all (AN≥0) |
|---|---|---|---|---|
| SNP   | **0.9999** | 0.9998 | 0.9994 | 0.1497 |
| indel | **0.9960** | 0.9945 | 0.9892 | 0.1294 |

Two structural checks reinforce the result:
- **Incremental is a superset of GLnexus.** 99.4 % of GLnexus sites (25.44 M /
  25.59 M) are present in the incremental DB; the incremental adds ~18 M more —
  mostly rare / low-quality sites GLnexus's unifier is more conservative about
  emitting. Exactly what we want for a *complete* in-house AF.
- **The low overall r (~0.15) is the rare/low-AN tail, not a regression.** Only
  ~6 % of sites have AN<1000, but at those sites the two denominators diverge
  most (GLnexus no-calls low-confidence samples → smaller AN → inflated AF),
  and rare variants there swing AF hardest. It reads *lower* than the 64-sample
  ~0.5 purely because at N=675 the rare/singleton tail dominates the site count.
  The metric that matters — well-called sites (AN≥1200/1300) — is near-perfect.

Conclusion: the incremental AC/AN/AF machinery matches joint genotyping on
well-called sites at production scale, and is strictly more complete. The
incremental DB is the production source of truth; GLnexus stays as a quarterly
audit.

---

## 7. Pitfalls / edge cases

- **Indel AN double-count (fixed).** A variant's callable interval is a SINGLE
  anchor base `[pos-1, pos)`; using the full REF span makes a deletion overlap
  the same sample's adjacent ref block → AN counted twice. Publish looks up AN
  at `pos-1` (0-based) to match the bedGraph. Caught in Phase C (indel `incAN`
  reached 226 for 64 samples, impossible >128).
- **gVCF depth field.** DRAGEN ref blocks expose `MIN_DP` (we already rely on
  this via GLnexus `ref_dp_format: MIN_DP`); variant records use `DP`. Confirm on
  real data before building.
- **Normalization must match** the annotation TSV and the GLnexus path
  (`norm -m- -f ref`), or keys won't join.
- **Sex ambiguity / aneuploidy** (XXY, X0, XYY, contamination, mosaics):
  flagged and excluded from X/Y (autosome+M still counted). Surfaced in
  `manifest`/`qc.json` for review.
- **Female chrY:** excluded entirely (AN weight 0 and Y variants dropped) to
  avoid spurious Y calls from X/Y homology inflating counts.
- **chrM:** treated as haploid **carrier frequency** — `AN += 1` per callable
  sample, `AC += 1` per carrier (PASS, DP≥10). Heteroplasmy fraction is a
  per-sample property shown elsewhere in the UI, **not** folded into cohort AF.
- **`AN=0` / `AC=0` sites** are not emitted (same cleanup as the GLnexus path).
- **AN-track growth.** Breakpoints accumulate with cohort size; mitigate with
  RLE + periodic recompaction; track is bgzipped+tabixed.
- **AN as a coverage proxy.** "callable (DP≥10)" ≈ "has a confident genotype."
  This is an approximation vs joint genotyping's actual GT; quantified by the
  quarterly reconciliation. Acceptable for an annotation/filtering DB.
- **No MRN/family dedup yet** (manifest carries the hooks): related individuals
  and re-sequenced patients still count more than once. Revisit before using
  in-house AF as anything beyond annotation/artifact-filtering.

---

## 8. Portability (DGM now → air-gapped DGX later)

Every component is stdlib Python + `bcftools` + `bedtools`/sweep-line — **no
FPGA, no network, no container required**. Same `--strip-path-prefix` / env-path
conventions as the existing scripts. The only piece that does *not* move to the
DGX is the optional GLnexus reconciliation (CPU, fine) — and DRAGEN's iterative
genotyper (FPGA, excluded by this design on purpose).

---

## 9. Proposed build order (when we implement)

| Phase | Deliverable | Validation |
|---|---|---|
| A | per-sample tools: sex inference, callable weighted BED, variant counts | on the 8-sample set: sex calls sane; per-sample AC vs GLnexus 8-sample |
| B | AN-track accumulator + publish | 8-sample incremental `INHOUSE_AF` vs GLnexus 8-sample AF (scatter ≈ y=x) |
| C | full 633 ingest | compare against the GLnexus full-633 snapshot (ground truth) |
| D | per-batch automation + reconciliation harness | re-run a batch → idempotent; drift report |
| E | UI integration (separate branch, later) | unchanged output contract → adapter + badge + filter |

Components A/B reuse `select_cohort.py` and the normalization/atomic-install
conventions already in this folder. Nothing here touches `run_stopgaps.sh` or any
UI script — that is Phase E, on the appropriate branch, later.

---

## 10. Open items

1. ~~gVCF fields~~ — **confirmed**: ref blocks carry `MIN_DP`; variant counts
   come from `hard-filtered.vcf.gz` so the gVCF's variable variant FORMAT is moot.
2. ~~Sex-call thresholds~~ — **resolved**: use DRAGEN's `Ploidy estimation`
   field; our coverage inference is fallback only.
3. **AN-track tooling** — `bedtools v2.31.1` is available on DGM; leaning toward
   using it (one static binary to stage on the DGX) vs a pure-Python sweep-line.
   Decide at build time.
4. ~~chrM carrier definition~~ — **resolved (4.6)**: GT is diploid-style, `AF`
   is heteroplasmy; count haploid carriers (PASS, weight 1), AN = callable mito
   genomes, store homoplasmic vs heteroplasmic. v1 threshold = PASS only.
   *To confirm: is "PASS-only" the right carrier threshold, or do you want a
   minimum heteroplasmy (e.g. AF ≥ 0.10)?*
5. **Confirm the `other/{sample}/` layout is uniform** across runs (some newer
   runs nest under `other/{sample}/germline_seq/`) so `select_cohort.py` can find
   the ploidy CSV + indexed gVCF for every sample. (`ploidy CSV` coverage already
   confirmed: 0 cohort samples missing.)
