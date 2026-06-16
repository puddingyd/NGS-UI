# Dead-zone integration guide (for the clinical reporting platform)

This document describes the panel dead-zone outputs and how to ingest them
when generating per-sample clinical reports. It is intended for the
downstream analysis / reporting platform engineer, not for the validation
analyst (who should read `README.md` first).

---

## What "dead-zone" means here

For every protein-coding CDS exon in the validated gene panel
(`panel/panel_loose_plus_clinical.hgnc_canonical.txt`, **7,363** HGNC-current
gene symbols = curated disease list 6,240 ∪ clinician specialty panels +1,123),
we ran `mosdepth` across a cohort of clinical-grade samples and recorded, per
exon:

> The cohort coverage (`cohort.exon_detail.tsv`) is genome-wide (all MANE
> genes), so the panel is applied as a report-time filter — no re-sequencing is
> needed when the panel changes.

- **median mean coverage** across the cohort
- **fraction of bp at ≥ N×** for thresholds N ∈ {10, 15, 20, 30}
- a binary `is_dead_<N>x` flag: 1 if cohort-median fraction-at-N× < 95 %

A "dead exon" at threshold N× means the cohort cannot reliably detect
variants in that exon at depth ≥ N. For an individual sample, this is the
prior probability that the sample will also fail there — we recommend the
report adds a limitation note for variants reported in (or absent from)
these exons.

The dead-zone is computed per pipeline because alignment quality, error
profiles, and reference indexing differ:

| Pipeline | Production status | Dataset | Thresholds available |
| --- | --- | --- | --- |
| WGS DRAGEN | **★ primary clinical pipeline** | 60 NA-validation samples | 10, 15, 20, 30 |
| WGS in_house | fallback (DV + HaplotypeCaller+VQSR union) | 60 samples | 10, 20, 30 |
| WES in_house | WES capture-based | 36 samples | 10, 20, 30 |

---

## Threshold to use

Use the current reporting threshold selected by the clinical reporting
platform. DRAGEN WGS ships multiple threshold columns in the same reference
package, and the platform currently reports DRAGEN WGS dead-zone at 10×.

| Pipeline | Threshold | Rationale |
| --- | ---: | --- |
| **WGS DRAGEN** | **10×** | Current NGS-UI reporting threshold; 15×/20×/30× columns remain available for comparison and audit. |
| WGS in_house | 20× | Both SNV and INDEL LoD = 20×; sensitivity ≥ 96 % at 20-25× bin. |
| WES in_house | 20× | SNV LoD = 20× (90.5 % sens); INDEL nominally 30× but 20× is the common clinical floor with a documented caveat (84.3 % INDEL sensitivity). |

Numbers above the reporting threshold (e.g. dead at 30×) are informational
only and overstate the clinical limitation. Stay at the platform-selected
threshold for clinical reporting unless the reporting policy changes.

---

## Files to ingest

For each pipeline, the relevant directory is
`results/wgs/<pipeline>/dead_zone_cds_canonical/` (or
`results/wes/expand_cds_canonical/` for WES). Same file naming scheme in
all three. Examples below use the DRAGEN WGS variant.

### Long format — primary ingest target ★

`wgs_dragen_panel_dead_exons.tsv`

One row per (panel gene, transcript, exon) that is dead at *any* threshold.
Each threshold has its own boolean column. Universal CSV/TSV — no list-encoded
columns to parse.

| Column | Type | Notes |
| --- | --- | --- |
| `gene` | str | HGNC current symbol |
| `transcript` | str | Ensembl ENST or RefSeq NM_ |
| `exon` | int | exon number per the MANE Select transcript (matches clinical exon numbering — BRCA1 exon 11 = our row exon=11) |
| `coord` | str | `chr17:43090933-43091042` |
| `exon_len` | int | bp |
| `median_cov` | float | cohort median of per-sample mean coverage |
| `is_dead_10x` | 0/1 | |
| `is_dead_15x` | 0/1 | DRAGEN only |
| `is_dead_20x` | 0/1 | |
| `is_dead_30x` | 0/1 | |
| `gene_cds_len` | int | total CDS length (bp) of this gene's transcript |
| `gene_dead_10x_pct` | float | **% of this gene's CDS that is dead at 10x** (gene-level; same value on every exon row of the gene) |
| `gene_dead_15x_pct` | float | gene-level CDS dead % at 15x (DRAGEN only) |
| `gene_dead_20x_pct` | float | gene-level CDS dead % at 20x |
| `gene_dead_30x_pct` | float | gene-level CDS dead % at 30x |

The `gene_*` columns are **gene-level** (identical across all exon rows of a
gene) and let a per-exon view show "this exon is dead, and the gene is X% dead
CDS overall".  Appended after the original columns — existing parsers unaffected.

```python
import pandas as pd
dz = pd.read_csv("dead_zone/dragen_wgs/wgs_dragen_panel_dead_exons.tsv",
                 sep="\t")
# Dead exons at the clinical threshold (DRAGEN 10x):
dead_10 = dz[dz["is_dead_10x"] == 1]
# Lookup: is BRCA1 exon 11 dead?
brca1_11 = dz[(dz.gene == "BRCA1") & (dz.exon == 11)]
```

### Per-gene wide format — convenient for gene-level UI ★

`wgs_dragen_panel_dead_exon_summary.tsv`

One row per (gene, transcript). Dead-exon numbers are **expanded** comma
lists (no run-length encoding), so they're trivially split into integers.

Columns: `gene  transcript  n_cds_exons  dead_10x_n  dead_10x_exons  dead_15x_n  dead_15x_exons  dead_20x_n  dead_20x_exons  dead_30x_n  dead_30x_exons  cds_len  dead_10x_len  dead_10x_pct  dead_15x_len  dead_15x_pct  dead_20x_len  dead_20x_pct  dead_30x_len  dead_30x_pct`

`dead_<N>x_exons` is `"5,7,8,9,15"` (each entry a positive integer; empty
string if no dead exons at that threshold).

**Per-gene CDS dead fraction (appended columns):**
- `cds_len` — total CDS length (bp) of the transcript (sum of CDS-exon lengths).
- `dead_<N>x_len` — bp of CDS lying in dead exons at threshold N.
- `dead_<N>x_pct` — `dead_<N>x_len / cds_len × 100`, 1 d.p. (e.g. `33.2`).
  This is the share of the gene's coding sequence that falls in dead exons —
  display it next to the dead-exon marks for a gene-level "how much is affected".

These columns were **appended** after the original ones, so existing
column positions / parsers are unaffected.

```python
import pandas as pd
summary = pd.read_csv("…/wgs_dragen_panel_dead_exon_summary.tsv", sep="\t")

# Get DRC4's clinical dead exons (10x):
row = summary[summary.gene == "DRC4"].iloc[0]
dead_exons = [int(x) for x in row["dead_10x_exons"].split(",") if x] if row["dead_10x_exons"] else []
```

### JSON Lines — easiest universal ingest

`wgs_dragen_panel_dead_exons.jsonl`

One JSON object per gene-transcript, exon lists as JSON arrays.

```json
{"gene": "BRCA1", "transcript": "ENST00000357654", "n_cds_exons": 22, "cds_len": 5592,
 "dead_10x": [], "dead_10x_len": 0, "dead_10x_pct": 0.0, "dead_15x": [], "dead_15x_len": 0,
 "dead_15x_pct": 0.0, "dead_20x": [], "dead_20x_len": 0, "dead_20x_pct": 0.0,
 "dead_30x": [], "dead_30x_len": 0, "dead_30x_pct": 0.0}
```

Each `dead_<N>x` array (the dead exon numbers) now has a companion
`dead_<N>x_len` (bp in dead exons) and `dead_<N>x_pct` (% of `cds_len` that is
dead) — the per-gene "fraction of CDS that is dead" to show alongside the marks.

Use this when your platform can ingest JSON Lines natively (pandas
`read_json(lines=True)`, DuckDB `read_json_auto`, R `jsonlite::stream_in`,
…).

### RLE summary — humans only, do NOT parse

`wgs_dragen_panel_dead_exon_summary_rle.tsv`

Same shape as the expanded summary, but dead-exon lists are run-length-encoded
(e.g. `"1,5-12,18"`). For human review / docx tables only. Parsers should
use the expanded `.tsv` or `.jsonl` instead — RLE is fragile to changes in
delimiter convention.

### Standalone markdown report

`wgs_dragen_panel_dead_exons.md`

Header counts + the top RLE table per pipeline. Used to drop into the
formal validation docx (Section 2.4). Not for downstream parsing.

---

## Recommended ingest pattern per clinical report

```python
import pandas as pd

DZ_PATH = "dead_zone/dragen_wgs/wgs_dragen_panel_dead_exons.tsv"
THR_COL = "is_dead_10x"      # DRAGEN clinical threshold

dz = pd.read_csv(DZ_PATH, sep="\t")

def dead_exons_for(gene):
    """Return list of dead exon numbers for `gene` at the clinical threshold."""
    hits = dz[(dz.gene == gene) & (dz[THR_COL] == 1)]
    return sorted(hits["exon"].astype(int).tolist())

def limitation_note(gene):
    exs = dead_exons_for(gene)
    if not exs:
        return None
    # Compact for human display: "1, 5-12, 18"
    runs = []
    s = e = exs[0]
    for x in exs[1:]:
        if x == e + 1: e = x
        else:
            runs.append(str(s) if s == e else f"{s}-{e}")
            s = e = x
    runs.append(str(s) if s == e else f"{s}-{e}")
    return (f"{gene} exon(s) {', '.join(runs)} 之 cohort 中位涵蓋深度 < 15×, "
            "本檢測對該區段之變異判讀敏感度不足；陰性結果不代表該區段無變異。")

# Use during report generation:
note = limitation_note("USP9Y")
if note:
    report.add_limitation(note)
```

If the per-sample VCF has the affected gene's variants but the variant is in
a dead exon, also flag the **call** with reduced confidence (positive
findings in a low-coverage region need orthogonal confirmation per ACMG
guidelines).

---

## What's NOT in dead-zone

- **lncRNA / miRNA / pseudogene / mtDNA tRNA / locus-only entries** — 73
  panel genes have no MANE Select / RefSeq protein-coding CDS, so they
  don't appear in the dead-zone tables. They're listed in
  `results/panel_alignment/panel_unmapped.tsv`. If the per-sample report
  needs to comment on these, treat them as "outside the validated
  dead-zone analysis"; they're still in the panel for SNV/INDEL reporting.
- **UTR regions** — dead-zone is CDS-only by design. UTR variants in panel
  genes still get annotated by VEP but aren't covered by this analysis.
- **Structural variants / CNV** — coverage tells you whether short-read
  small-variant calling will work; CNV calling has its own performance
  envelope (gCNV / CNVkit / Delly / Manta / DRAGEN-CNV).
- **Per-sample dead-zones** — the analysis is *cohort-level*. An
  individual sample with abnormally low overall coverage can have many
  more dead exons than the cohort median; rely on the per-sample QC
  metrics (`scripts/wes_qc_metrics.py` / `wgs_qc_metrics*.py`) to gate
  whether the cohort dead-zone is even a fair prior.

---

## Regenerating these files

On DGM:

```bash
cd ~/NGS-validation
git pull origin ngs-validation

# Reruns mosdepth + aggregator + panel report; one entry point each.
bash scripts/dead_zone/rerun_with_canonical_panel.sh    # WES + in_house WGS
bash scripts/dead_zone/dragen_wgs_dead_zone.sh          # DRAGEN WGS (10x/15x/20x/30x)

git add results/{wes,wgs}/**/dead_zone_cds_canonical
git commit -m '...'
git push
```

Each entry point auto-detects the threshold list from `mosdepth`'s output
header — you can change the threshold set by exporting `THRESHOLDS=…`
before invoking the wrapper (e.g. `THRESHOLDS=10,15,20,25` to add 25×).
The output schema adapts automatically.

---

## Versioning

When this file's contents change (different threshold, different panel,
re-run on a different cohort), update the docx (Section 2.4 of
`WGS_validation_report_sections.docx`) by re-running
`python3 scripts/wgs_validation/build_wgs_report.py` and check the
cohort + threshold quoted in the platform code matches.

The current state (2026-05-26):

| | DRAGEN WGS | in_house WGS | in_house WES |
| --- | --- | --- | --- |
| panel | `panel_loose_plus_clinical.hgnc_canonical.txt` (7,363 genes) | same | same |
| BAMs cohort | 60 samples NovaSeq 2026-04-28 | 60 samples (in-house pipeline) | 36 samples (WES_VAL_BATCH1) |
| clinical threshold | **15×** | 20× | 20× |
| panel-gene coverage at clinical threshold | 99.53 % (29 genes with ≥ 1 dead exon) | 89.88 % (624 genes) | 89.88 % (effectively similar; capture-bias-driven) |
