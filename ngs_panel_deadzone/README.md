# NGS panel + dead-zone delivery package

Self-contained data + docs for the clinical reporting UI / backend. Drop
this whole folder into the UI repo. Everything the backend needs to
(a) translate VEP gene names to HGNC-current and (b) annotate per-exon
coverage limitations is here — no access to the validation repo required.

Generated from the NGS-validation repo. Re-delivered whenever the panel,
HGNC release, or VEP cache changes (see "Versioning" in each guide).

---

## Start here

1. **`panel_annotation.md`** — how to translate VEP's per-variant `SYMBOL`
   into the HGNC-current names used in reports. Read this first; it has the
   ingest pseudo-code.
2. **`dead_zone.md`** — how to annotate per-exon coverage ("dead-zone")
   limitations for a sample's reported genes.

Both guides' file paths have been rewritten to be **relative to this
package** (e.g. `panel/vep_alias_map.tsv`, `dead_zone/dragen_wgs/...`).

---

## Layout

```
ngs_panel_deadzone/
├── README.md                    ← you are here (index)
├── panel_annotation.md          ← VEP→HGNC name translation guide
├── dead_zone.md                 ← per-exon coverage limitation guide
│
├── panel/
│   ├── panel_loose.hgnc_canonical.txt   panel scope: 6,240 HGNC-current symbols (one/line)
│   ├── vep_alias_map.tsv                31 panel genes whose VEP output needs translation
│   └── hgnc_id_to_symbol.tsv            HGNC_ID → current symbol, all 44,989 Approved genes
│
└── dead_zone/
    ├── dragen_wgs/      ★ PRIMARY clinical pipeline (threshold 15×)
    ├── inhouse_wgs/        in-house WGS fallback (threshold 20×)
    └── wes/                WES capture (threshold 20×)
```

Each `dead_zone/<pipeline>/` folder contains the same five files (prefix
differs by pipeline — `wgs_dragen_`, `wgs_`, `wes_`):

| File | Use |
| --- | --- |
| `*_panel_dead_exons.tsv` | ★ long format, one row per dead exon — **primary ingest** |
| `*_panel_dead_exons.jsonl` | same data, JSON Lines (exon arrays per gene) |
| `*_panel_dead_exon_summary.tsv` | one row per gene, exon numbers expanded as comma list |
| `*_panel_dead_exon_summary_rle.tsv` | one row per gene, run-length-encoded (human review only) |
| `*_panel_dead_exons.md` | human-readable summary |

---

## Which pipeline?

`dragen_wgs/` is the production clinical pipeline. Use it unless the
sample was run on a different pipeline. Threshold rationale and the full
per-pipeline comparison are in `dead_zone.md`.

---

## Two-minute sanity check

```bash
cd ngs_panel_deadzone

# panel scope size
wc -l panel/panel_loose.hgnc_canonical.txt          # 6240

# alias-map kinds
cut -f5 panel/vep_alias_map.tsv | tail -n +2 | sort | uniq -c
#   1 ensembl_id_drift
#  11 no_vep_record
#   4 vep_missing_gene
#  15 vep_uses_old_symbol

# dead-exon rows for the primary pipeline
wc -l dead_zone/dragen_wgs/wgs_dragen_panel_dead_exons.tsv
```

---

## What's intentionally NOT here

- Raw VCFs, BAMs, validation truth sets — not needed for reporting.
- The full HGNC complete set (54 cols, 16 MB) — replaced by the 2-col
  `hgnc_id_to_symbol.tsv` the guides actually use.
- The validation scripts themselves — they live in the NGS-validation
  repo; this package is data + integration docs only.
