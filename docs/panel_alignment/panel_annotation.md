# Gene-name annotation guide (for the backend analysis platform)

How the backend should translate VEP's per-variant SYMBOL output into the
HGNC-current names used in our clinical reports.

This is the companion to `dead_zone.md` (CDS coverage limitations) — both
are intended for whoever wires VEP CSQ output into the clinical report
generator.

---

## The problem

The reporting layer wants to emit, e.g., **GBA1** or **DRC4** or **SACK1H**
(HGNC current symbols, panel scope). VEP 115's cache, however, sometimes
outputs older symbols (e.g. `GAS8` for what HGNC now calls **DRC4**), and
in a few cases doesn't have a transcript for the panel gene at all
(reports a different overlapping gene instead).

We did a per-gene alignment check (synthetic-VCF probe of all 6,240 panel
genes against the production VEP cache). The check produced the alias map
described below — apply it once at the variant-annotation step and the
report will always use HGNC-current names.

Per-gene status:

| Outcome | Count | What the backend should do |
| --- | ---: | --- |
| `match` | 6,250 | nothing — VEP's SYMBOL already matches HGNC current |
| `vep_uses_old_symbol` | 15 unique genes | translate VEP SYMBOL → HGNC current (alias map) |
| `vep_missing_gene` | 4 | translate by HGNC_ID (VEP labels these regions with a *different* gene's HGNC_ID; use cross-reference) |
| `ensembl_id_drift` | 1 (CAST) | nothing — same SYMBOL + HGNC_ID, Ensembl gene_id drifted |
| `no_vep_record` | 11 | accept that pre-verification wasn't possible; trust VEP for real-sample variants and use HGNC_ID join (see "Untested edge cases" below) |

---

## The two files the backend needs

| File | Purpose | Rows |
| --- | --- | ---: |
| `results/panel_alignment/panel_loose.hgnc_canonical.txt` | the panel scope (HGNC-current symbols) | 6,240 |
| `results/panel_alignment/vep_alias_map.tsv` | every panel gene whose VEP output needs translation; rest of the panel matches VEP exactly | 30 |

`vep_alias_map.tsv` columns:

| Column | Meaning |
| --- | --- |
| `panel_hgnc_id` | HGNC ID to use in the report (the *canonical* key) |
| `panel_hgnc_symbol` | HGNC-current symbol to display in the report |
| `vep_symbol` | what VEP will put in the CSQ `SYMBOL` field |
| `vep_hgnc_id` | what VEP will put in the CSQ `HGNC_ID` field |
| `kind` | one of: `vep_uses_old_symbol`, `vep_missing_gene`, `ensembl_id_drift`, `no_vep_record` |
| `notes` | short human-readable description |

For `no_vep_record` rows the `vep_symbol` / `vep_hgnc_id` columns are
empty (we couldn't pre-validate). For all other kinds they are populated
with whatever VEP will actually emit.

---

## Recommended translation algorithm (HGNC_ID first)

VEP's `HGNC_ID` field is more reliable than its `SYMBOL` field — even when
VEP shows an older symbol, the HGNC_ID typically matches. **Use HGNC_ID as
the primary join key.**

```python
import csv

# --- one-time load ---
panel_hgnc_ids = set()        # the panel scope
with open("panel_loose.hgnc_canonical.txt") as fh:
    panel_symbols = [line.strip() for line in fh if line.strip()]

hgnc_to_display = {}          # HGNC_ID -> HGNC-current symbol
# (load from data/hgnc/hgnc_complete_set.txt: status=Approved, hgnc_id, symbol)
with open("hgnc_complete_set.txt", newline="") as fh:
    r = csv.DictReader(fh, delimiter="\t")
    for row in r:
        if row["status"] == "Approved":
            hgnc_to_display[row["hgnc_id"]] = row["symbol"]

# Special-case map for vep_missing_gene (Ensembl-split / read-through cases).
# Key = VEP's HGNC_ID, value = (panel HGNC_ID, panel display symbol, kind).
# In real-sample VEP output these variants will carry VEP's HGNC_ID (e.g.
# HGNC:14862 for POLR2M). The reporter may want to ALSO surface the panel
# gene (e.g. GCOM1) on those variants.
vep_to_panel_alias = {}       # VEP HGNC_ID -> (panel HGNC_ID, panel symbol)
with open("vep_alias_map.tsv", newline="") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        if r["kind"] in ("vep_missing_gene",) and r["vep_hgnc_id"]:
            vep_to_panel_alias[r["vep_hgnc_id"]] = (
                r["panel_hgnc_id"], r["panel_hgnc_symbol"]
            )

# --- per variant, per CSQ entry ---
def report_gene_for_csq(csq):
    """csq is a dict-like with keys SYMBOL, HGNC_ID, ... from VEP CSQ.
    Returns the (HGNC_ID, display_symbol) to put in the report row, or None
    if this CSQ entry isn't in the panel."""
    hid = csq.get("HGNC_ID", "")
    if hid:
        # 1. Direct: panel gene reachable via VEP HGNC_ID
        if hid in panel_hgnc_ids:
            return hid, hgnc_to_display.get(hid, csq["SYMBOL"])
        # 2. vep_missing_gene cross-reference (POLR2M -> GCOM1, etc.)
        if hid in vep_to_panel_alias:
            panel_hid, panel_sym = vep_to_panel_alias[hid]
            return panel_hid, panel_sym
    # 3. Fallback by SYMBOL (covers the no_vep_record edge cases where
    #    HGNC_ID is empty but SYMBOL is set, and pseudogene cases where
    #    VEP returns a SYMBOL not in HGNC)
    sym = csq.get("SYMBOL", "")
    return None  # not in panel
```

The 15 `vep_uses_old_symbol` cases are handled automatically by **step
1** — VEP's `HGNC_ID` matches the panel's HGNC_ID, and we look up the
display symbol from the HGNC table (not from VEP). VEP's older SYMBOL is
ignored.

The 4 `vep_missing_gene` cases are handled by **step 2** — VEP labels
those regions with a different gene's HGNC_ID; the lookup translates
back. This expects the report layer to be OK with showing **GCOM1** as
the panel context for variants VEP labelled as POLR2M, etc.

---

## How each of the four `vep_missing_gene` cases works in practice

| Panel gene (clinical report) | VEP will actually label as | What the backend should do |
| --- | --- | --- |
| **GCOM1** (HGNC:26424) | POLR2M (HGNC:14862) and others | When the report needs to discuss GCOM1: also surface any variants VEP annotated as POLR2M / MYZAP / ARMH3 in the chr15:57.7 Mb locus. The GCOM1 designation is a read-through — Ensembl 115 split it. |
| **LRTOMT** (HGNC:25033) | TOMT (HGNC:55527), LRRC51 | When reporting LRTOMT: surface TOMT + LRRC51 variants in the chr11:72 Mb locus. |
| **DISC2** (HGNC:2889) | DISC1 (HGNC:2888) | DISC2 is an antisense lncRNA of DISC1; Ensembl 115 doesn't have a separate DISC2 transcript. Variants reported under DISC2 will be labeled as DISC1 by VEP (intronic). Recommend not reporting DISC2 as a primary gene without orthogonal review. |
| **TRU-TCA1-1** (HGNC:12348) | FOSB (HGNC:3797) (downstream artifact) | TRU-TCA1-1 is a single-exon tRNA; VEP doesn't annotate it as a gene. The "FOSB" label is the nearest downstream gene VEP found at the position we probed. **Do not use this alias** — TRU-TCA1-1 variants will not be discoverable through VEP output. Flag as a known limitation. |

---

## How each of the 15 `vep_uses_old_symbol` cases works

These are the simple cases — same gene, VEP just uses an older name.
The backend's HGNC_ID lookup handles them automatically. They're listed
here for the validation report and for human review:

| HGNC current (report) | VEP cache (output) | HGNC_ID |
| --- | --- | --- |
| ASPNAT | NAT8L | HGNC:26742 |
| CATSPERT | C2CD6 | HGNC:14438 |
| COXFA4 | NDUFA4 | HGNC:7687 |
| DRC2 | CCDC65 | HGNC:29937 |
| DRC4 | GAS8 | HGNC:4166 |
| FAM194C | C3orf20 | HGNC:25320 |
| IFT38 | CLUAP1 | HGNC:19009 |
| IFT54 | TRAF3IP1 | HGNC:17861 |
| QTMAN | GTDC1 | HGNC:20887 |
| RMP64 | NEPRO | HGNC:24496 |
| SACK1G | FAM83G | HGNC:32554 |
| SACK1H | FAM83H | HGNC:24797 |
| SLC38A12 | TMEM104 | HGNC:25984 |
| SLC9D1 | TMCO3 | HGNC:20329 |
| TIMCC | FAM136A | HGNC:25911 |

When VEP 115 is upgraded to a future cache that uses HGNC-current names,
these rows will simply move into the `match` category — no backend code
changes needed.

---

## Untested edge cases (`no_vep_record`)

These 11 panel genes failed our pre-validation. The reasons are physical
limitations, not pipeline bugs:

| Panel gene | HGNC_ID | Reason | Real-sample handling |
| --- | --- | --- | --- |
| FEB1, TCL4 | (none) | OMIM-only locus with no HGNC entry | Real samples will never have VEP output labeled "FEB1" / "TCL4" — these aren't VEP-known genes. If the report needs to discuss them, treat as locus-level (cytoband only). |
| ABO, ORAI1, CCL3L1, ATXN8, FCGR2C | HGNC:79, HGNC:25896, HGNC:10628, HGNC:32925, HGNC:15626 | alt-locus Ensembl gene IDs; not in MANE Select or VEP primary-assembly cache at the expected positions | Real-sample variants in these may or may not be annotated; cross-check with a UCSC/Ensembl direct lookup before reporting. |
| CASP12, GULOP, UOX, GGT2P | HGNC:19004, HGNC:4695, HGNC:12575, HGNC:4251 | pseudogenes; VEP cache may or may not have an entry | Pseudogene variants are usually not clinically reportable; flag if anything shows up. |

For all 11: **trust VEP's HGNC_ID if it appears in real-sample output**.
The alignment check just couldn't pre-verify the naming because we
couldn't synthesize a variant at the right position.

---

## Verification

The alias map is regenerated from `gene_alignment_report.tsv` whenever
the panel alignment is re-run:

```bash
python3 scripts/panel_alignment/build_alias_map.py
```

A complete rebuild of the alignment (after a panel change, HGNC update,
or VEP cache upgrade) goes:

```bash
# on DGM
bash scripts/panel_alignment/rebuild_alignment.sh
# back on laptop
python3 scripts/panel_alignment/parse_vep_alignment.py \
    --vep-vcf      results/panel_alignment/synthetic.vep.vcf.gz \
    --panel-hgnc   results/panel_alignment/panel_loose.hgnc.tsv \
    --out-report   results/panel_alignment/gene_alignment_report.tsv \
    --out-disagree results/panel_alignment/vep_disagreements.tsv
python3 scripts/panel_alignment/build_alias_map.py
```

After every rebuild, **re-send `vep_alias_map.tsv` and
`panel_loose.hgnc_canonical.txt` to the backend** and rerun their ingest.

---

## Versioning

| Item | Current value |
| --- | --- |
| Panel | `panel_loose.hgnc_canonical.txt`, 6,240 unique HGNC symbols |
| Panel source databases | OMIM (2026-05-10) + ClinGen (2026-05-10) + GenCC (2026-05-10), HGNC-aligned |
| VEP version | Ensembl VEP 115, `vep_115.sif`, cache at `/home/pipeline/reference/hg38/tertiary/vep_cache/` |
| Alignment date | as of `gene_alignment_report.tsv` commit (see `git log results/panel_alignment/gene_alignment_report.tsv`) |
| Alias map rows | 30 (15 `vep_uses_old_symbol` + 4 `vep_missing_gene` + 1 `ensembl_id_drift` + 10 `no_vep_record`) |
