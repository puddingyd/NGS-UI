# Dead-zone integration guide (for the clinical reporting platform)

This document describes the panel dead-zone outputs and how to ingest them
when generating per-sample clinical reports. It is intended for the
downstream analysis / reporting platform engineer, not for the validation
analyst (who should read `README.md` first).

---

## What "dead-zone" means here

For every protein-coding CDS exon in the validated gene panel
(`results/panel_alignment/panel_loose.hgnc_canonical.txt`, 6,240 HGNC-current
gene symbols), we ran `mosdepth` across a cohort of clinical-grade samples
and recorded, per exon:

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

Anchored to the validated LoD floor in Section 2.3.1 of the WGS validation
report (`results/wgs/validation/{in_house,dragen}/WGS_validation_report_sections.docx`):

| Pipeline | Threshold | Rationale |
| --- | ---: | --- |
| **WGS DRAGEN** | **15×** | DRAGEN INDEL LoD = 15× (SNV LoD ≥ 5× from data; 15× is the conservative of the two). At 15× SNV sensitivity = 98.6 %, INDEL = 90.2 %. |
| WGS in_house | 20× | Both SNV and INDEL LoD = 20×; sensitivity ≥ 96 % at 20-25× bin. |
| WES in_house | 20× | SNV LoD = 20× (90.5 % sens); INDEL nominally 30× but 20× is the common clinical floor with a documented caveat (84.3 % INDEL sensitivity). |

Numbers below the threshold (e.g. dead at 10×) are **informational only** —
they overstate the limitation. Numbers above (e.g. dead at 30×) overstate
the limitation in the other direction (because at WGS ~30× mean coverage,
~50 % of bp drop below 30× by definition). Stay at the LoD-anchored
threshold for clinical reporting.

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

```python
import pandas as pd
dz = pd.read_csv("dead_zone/dragen_wgs/wgs_dragen_panel_dead_exons.tsv",
                 sep="\t")
# Dead exons at the clinical threshold (DRAGEN 15x):
dead_15 = dz[dz["is_dead_15x"] == 1]
# Lookup: is BRCA1 exon 11 dead?
brca1_11 = dz[(dz.gene == "BRCA1") & (dz.exon == 11)]
```

### Per-gene wide format — convenient for gene-level UI ★

`wgs_dragen_panel_dead_exon_summary.tsv`

One row per (gene, transcript). Dead-exon numbers are **expanded** comma
lists (no run-length encoding), so they're trivially split into integers.

Columns: `gene  transcript  n_cds_exons  dead_10x_n  dead_10x_exons  dead_15x_n  dead_15x_exons  dead_20x_n  dead_20x_exons  dead_30x_n  dead_30x_exons`

`dead_<N>x_exons` is `"5,7,8,9,15"` (each entry a positive integer; empty
string if no dead exons at that threshold).

```python
import pandas as pd
summary = pd.read_csv("…/wgs_dragen_panel_dead_exon_summary.tsv", sep="\t")

# Get DRC4's clinical dead exons (15x):
row = summary[summary.gene == "DRC4"].iloc[0]
dead_exons = [int(x) for x in row["dead_15x_exons"].split(",") if x] if row["dead_15x_exons"] else []
```

### JSON Lines — easiest universal ingest

`wgs_dragen_panel_dead_exons.jsonl`

One JSON object per gene-transcript, exon lists as JSON arrays.

```json
{"gene": "BRCA1", "transcript": "ENST00000357654", "n_cds_exons": 22,
 "dead_10x": [], "dead_15x": [], "dead_20x": [], "dead_30x": []}
{"gene": "USP9Y", "transcript": "ENST00000338981", "n_cds_exons": 44,
 "dead_10x": [3,4,…,46], "dead_15x": [3,4,…,46], "dead_20x": […], "dead_30x": […]}
```

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
THR_COL = "is_dead_15x"      # DRAGEN clinical threshold

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
bash scripts/dead_zone/dragen_wgs_dead_zone.sh          # DRAGEN WGS (15x added)

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
| panel | `panel_loose.hgnc_canonical.txt` (6,240 genes) | same | same |
| BAMs cohort | 60 samples NovaSeq 2026-04-28 | 60 samples (in-house pipeline) | 36 samples (WES_VAL_BATCH1) |
| clinical threshold | **15×** | 20× | 20× |
| panel-gene coverage at clinical threshold | 99.53 % (29 genes with ≥ 1 dead exon) | 89.88 % (624 genes) | 89.88 % (effectively similar; capture-bias-driven) |
