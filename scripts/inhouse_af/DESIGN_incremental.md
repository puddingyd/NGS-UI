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

Anchor each sample on its **`other/{sample}/` directory** (richer than the
`vcf.gz/` delivery copy):

| need | file | notes |
|---|---|---|
| sex / karyotype | `{sample}.ploidy_estimation_metrics.csv` | `Ploidy estimation,XX/XY/…` (4.1) |
| callable / AN (MIN_DP) | `{sample}.hard-filtered.gvcf.gz` (+ `.tbi`) | ref blocks carry `MIN_DP`; **indexed** here (the `vcf.gz/` copy is not) |
| variant counts (AC) | `{sample}.hard-filtered.vcf.gz` (+ `.tbi`) | clean PASS + `GT:AD:AF:DP` (4.4) |
| (future) ROH / STR | `{sample}.roh.bed`, `{sample}.repeats.vcf.gz` | unrelated to AF; material for the UI's ROH/STR cards |

`select_cohort.py` should record the `other/{sample}/` dir per sample so ingest
finds all four. (Datalake files are read-only and the `vcf.gz/` copies are
un-indexed → ingest must not assume it can write indexes; it streams the gVCF
and uses the existing `.tbi` in `other/` when present.)

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

### 4.4 Variant counts (the numerator, accumulated)

**Source: the per-sample `{sample}.hard-filtered.vcf.gz`** (the regular VCF,
not the gVCF). It carries clean `FILTER=PASS` and a consistent
`GT:AD:AF:DP:…` FORMAT, so we sidestep the gVCF's variable variant-record
FORMAT (some gVCF variant rows are only `GT:GQ:PS`). The hard-filtered VCF is
small + indexed.

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

## 6. Reconciliation with full GLnexus (quarterly)

1. Full `build_inhouse_af.sh` over the whole cohort → a GLnexus snapshot.
2. Compare incremental `INHOUSE_AF` vs GLnexus `INHOUSE_AF` at shared sites:
   correlation, AF scatter, and a list of largest discrepancies.
3. Expected, benign differences: borderline genotype rescue (joint vs
   per-sample), complex-indel representation, AN coverage-proxy vs exact joint
   call. Large/systematic drift → investigate.

Production DB = the incremental one. GLnexus = audit/ground-truth. (Optionally,
adopt the GLnexus snapshot as the canonical DB each quarter and continue
incrementing from there — decide after we see the first comparison.) The very
first ground truth is the GLnexus full-633 run currently in progress.

---

## 7. Pitfalls / edge cases

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
4. **chrM carrier definition** — chrM has PASS variants with `AF` (heteroplasmy)
   + `DP`. Pending: confirm GT ploidy (`1` vs `0/1`) and whether to set a
   heteroplasmy threshold (e.g. count any PASS ALT as a carrier vs AF≥x).
5. **Confirm the `other/{sample}/` layout is uniform** across runs (some newer
   runs nest under `other/{sample}/germline_seq/`) so `select_cohort.py` can find
   the ploidy CSV + indexed gVCF + hard-filtered VCF for every sample.
