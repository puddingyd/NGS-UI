# Gene-name annotation guide (for the backend analysis platform)

How the backend should translate VEP's per-variant SYMBOL output into the
HGNC-current names used in our clinical reports.

This is the companion to `dead_zone.md` (CDS coverage limitations) — both
are intended for whoever wires VEP CSQ output into the clinical report
generator.

> **Scope (2026-06-12):** the VEP-alignment probe described below was re-run
> against the **expanded reportable panel** — `panel_loose_plus_clinical.hgnc_canonical.txt`,
> **7,363 genes** (curated disease list 6,240 ∪ clinician specialty panels +1,123).
> All 7,363 genes have been through the synthetic-VCF alignment check; the
> alias map and dead-zone files are both on the same 7,363 scope.

---

## The problem

The reporting layer wants to emit, e.g., **GBA1** or **DRC4** or **SACK1H**
(HGNC current symbols, panel scope). VEP 115's cache, however, sometimes
outputs older symbols (e.g. `GAS8` for what HGNC now calls **DRC4**), and
in a few cases doesn't have a transcript for the panel gene at all
(reports a different overlapping gene instead).

We did a per-gene alignment check (synthetic-VCF probe of all 7,363 panel
genes against the production VEP cache). The check produced the alias map
described below — apply it once at the variant-annotation step and the
report will always use HGNC-current names.

Per-gene status:

| Outcome | Count | What the backend should do |
| --- | ---: | --- |
| `match` | 7,325 | nothing — VEP's SYMBOL already matches HGNC current |
| `vep_uses_old_symbol` | 16 | translate VEP SYMBOL → HGNC current (alias map) |
| `vep_missing_gene` | 12 | translate by HGNC_ID (VEP labels these regions with a *different* gene's HGNC_ID; use cross-reference) |
| `ensembl_id_drift` | 1 (CAST) | nothing — same SYMBOL + HGNC_ID, Ensembl gene_id drifted |
| `disagree_other` | 1 (MAFIP→MIP) | translate by HGNC_ID (MAFIP pseudogene overlaps MIP) |
| `no_vep_record` | 8 | accept that pre-verification wasn't possible; trust VEP for real-sample variants and use HGNC_ID join (see "Untested edge cases" below) |

---

## The two files the backend needs

| File | Purpose | Rows |
| --- | --- | ---: |
| `panel/panel_loose_plus_clinical.hgnc_canonical.txt` | the panel scope (HGNC-current symbols) | 7,363 |
| `panel/vep_alias_map.tsv` | every panel gene whose VEP output needs translation; rest of the panel matches VEP exactly | 38 |
| `panel/hgnc_id_to_symbol.tsv` | HGNC_ID → current symbol lookup (all Approved genes) | 44,989 |

`vep_alias_map.tsv` columns:

| Column | Meaning |
| --- | --- |
| `panel_hgnc_id` | HGNC ID to use in the report (the *canonical* key) |
| `panel_hgnc_symbol` | HGNC-current symbol to display in the report |
| `vep_symbol` | what VEP will put in the CSQ `SYMBOL` field |
| `vep_hgnc_id` | what VEP will put in the CSQ `HGNC_ID` field |
| `kind` | one of: `vep_uses_old_symbol`, `vep_missing_gene`, `ensembl_id_drift`, `disagree_other`, `no_vep_record` |
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

# --- one-time load (paths are relative to this package) ---

# HGNC_ID <-> current symbol (2-col file shipped in this package)
hgnc_to_display = {}          # HGNC_ID -> HGNC-current symbol
symbol_to_hgnc  = {}          # HGNC-current symbol -> HGNC_ID
with open("panel/hgnc_id_to_symbol.tsv", newline="") as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        hgnc_to_display[row["hgnc_id"]] = row["symbol"]
        symbol_to_hgnc[row["symbol"]]   = row["hgnc_id"]

# Panel scope is a list of HGNC-current SYMBOLS; turn it into HGNC IDs
# (HGNC_ID is the join key we use against VEP output).
panel_symbols = set()
with open("panel/panel_loose_plus_clinical.hgnc_canonical.txt") as fh:
    panel_symbols = {line.strip() for line in fh if line.strip()}
panel_hgnc_ids = {symbol_to_hgnc[s] for s in panel_symbols if s in symbol_to_hgnc}

# Special-case map for vep_missing_gene / disagree_other (Ensembl-split,
# read-through, or overlapping-pseudogene cases).
# Key = VEP's HGNC_ID, value = (panel HGNC_ID, panel display symbol).
# In real-sample VEP output these variants will carry VEP's HGNC_ID (e.g.
# HGNC:14862 for POLR2M). The reporter may want to ALSO surface the panel
# gene (e.g. GCOM1) on those variants.
vep_to_panel_alias = {}       # VEP HGNC_ID -> (panel HGNC_ID, panel symbol)
with open("panel/vep_alias_map.tsv", newline="") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        if r["kind"] in ("vep_missing_gene", "disagree_other") and r["vep_hgnc_id"]:
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
        # 2. vep_missing_gene / disagree_other cross-reference
        #    (POLR2M -> GCOM1, MIP -> MAFIP, etc.)
        if hid in vep_to_panel_alias:
            panel_hid, panel_sym = vep_to_panel_alias[hid]
            return panel_hid, panel_sym
    # 3. Fallback by SYMBOL (covers the no_vep_record edge cases where
    #    HGNC_ID is empty but SYMBOL is set, and pseudogene cases where
    #    VEP returns a SYMBOL not in HGNC)
    sym = csq.get("SYMBOL", "")
    return None  # not in panel
```

The 16 `vep_uses_old_symbol` cases are handled automatically by **step
1** — VEP's `HGNC_ID` matches the panel's HGNC_ID, and we look up the
display symbol from the HGNC table (not from VEP). VEP's older SYMBOL is
ignored.

The 12 `vep_missing_gene` cases (plus the 1 `disagree_other` MAFIP case)
are handled by **step 2** — VEP labels those regions with a different
gene's HGNC_ID; the lookup translates back. This expects the report layer
to be OK with showing **GCOM1** as the panel context for variants VEP
labelled as POLR2M, etc.

---

## How each of the 12 `vep_missing_gene` cases works in practice

For each, VEP does not carry the panel gene's own transcript at the probed
position; it annotates the listed overlapping/adjacent gene instead. The
HGNC_ID cross-reference (step 2) maps VEP's label back to the panel gene.

| Panel gene (clinical report) | HGNC_ID | VEP will actually label as | VEP's HGNC_ID |
| --- | --- | --- | --- |
| **ATXN8**     | HGNC:32925 | ATXN8OS | HGNC:10561 |
| **CCL3L1**    | HGNC:10628 | CCL3L3  | HGNC:30554 |
| **DISC2**     | HGNC:2889  | DISC1   | HGNC:2888  |
| **FCGR2C**    | HGNC:15626 | FCGR3B  | HGNC:3620  |
| **GCOM1**     | HGNC:26424 | POLR2M  | HGNC:14862 |
| **GGT2P**     | HGNC:4251  | GNAZ    | HGNC:4395  |
| **GULOP**     | HGNC:4695  | PTK2B   | HGNC:9612  |
| **LRTOMT**    | HGNC:25033 | TOMT    | HGNC:55527 |
| **PGBD3**     | HGNC:19400 | ERCC6   | HGNC:3438  |
| **SNORA93**   | HGNC:50397 | FBLN2   | HGNC:3601  |
| **TRU-TCA1-1**| HGNC:12348 | FOSB    | HGNC:3797  |
| **UOX**       | HGNC:12575 | SSX2IP  | HGNC:16509 |

Notes on the trickier ones:

- **GCOM1** is a read-through that Ensembl 115 split — variants in the
  chr15:57.7 Mb locus get labelled POLR2M / MYZAP / ARMH3. When the report
  needs to discuss GCOM1, also surface those.
- **LRTOMT** — surface TOMT + LRRC51 variants in the chr11:72 Mb locus.
- **DISC2** is an antisense lncRNA of DISC1; Ensembl 115 has no separate
  DISC2 transcript. Variants land under DISC1 (intronic). Don't report
  DISC2 as a primary gene without orthogonal review.
- **TRU-TCA1-1** is a single-exon tRNA; VEP doesn't annotate it as a gene,
  so the "FOSB" label is just the nearest downstream gene at the probed
  position. **Do not use this alias for discovery** — TRU-TCA1-1 variants
  are not discoverable through VEP output. Flag as a known limitation.
- **PGBD3 / SNORA93 / ATXN8 / GGT2P / GULOP / UOX** are pseudogenes,
  snoRNA, or non-coding loci nested inside a host gene (ERCC6, FBLN2,
  ATXN8OS, GNAZ, PTK2B, SSX2IP); they are rarely clinically reportable on
  their own. Flag if anything shows up.

---

## The 1 `disagree_other` case (MAFIP → MIP)

| Panel gene | HGNC_ID | VEP labels as | VEP's HGNC_ID |
| --- | --- | --- | --- |
| **MAFIP** | HGNC:31102 | MIP | HGNC:7103 |

MAFIP is a pseudogene whose probed position overlaps **MIP**; VEP returns
MIP's SYMBOL, HGNC_ID and Ensembl gene_id (all three differ from MAFIP).
Treated like a `vep_missing_gene` cross-reference in the algorithm above
(step 2 picks it up via VEP's HGNC_ID). MAFIP is not independently
reportable.

---

## How each of the 16 `vep_uses_old_symbol` cases works

These are the simple cases — same gene, VEP just uses an older name.
The backend's HGNC_ID lookup handles them automatically. They're listed
here for the validation report and for human review:

| HGNC current (report) | VEP cache (output) | HGNC_ID |
| --- | --- | --- |
| ASPNAT | NAT8L | HGNC:26742 |
| CATSPERT | C2CD6 | HGNC:14438 |
| COXFA4 | NDUFA4 | HGNC:7687 |
| COXFA4L2 | NDUFA4L2 | HGNC:29836 |
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

These 8 panel genes returned no CSQ entry at the synthetic-probe position,
so we couldn't pre-verify their naming. The reasons are physical
limitations (pseudogenes, copy-number-variable loci, non-coding RNAs), not
pipeline bugs:

| Panel gene | HGNC_ID | Likely reason |
| --- | --- | --- |
| GSTT1     | HGNC:4641  | common germline copy-number-deletion locus; often absent from the primary-assembly cache at the probed position |
| HELLPAR   | HGNC:43984 | long non-coding RNA (macro-lncRNA); no protein transcript to probe |
| LRP5L     | HGNC:25323 | partial-duplication paralog; alt-locus annotation |
| NXF5      | HGNC:8075  | X-linked, paralog-rich family; cache coverage gap |
| POM121L9P | HGNC:30080 | pseudogene |
| SLC22A20P | HGNC:29867 | pseudogene |
| SSPOP     | HGNC:21998 | pseudogene of the SSPO locus |
| WTAPP1    | HGNC:44115 | pseudogene |

For all 8: **trust VEP's HGNC_ID if it appears in real-sample output**.
The alignment check just couldn't pre-verify the naming because we
couldn't synthesize a variant that VEP would annotate at the right
position. Most are pseudogenes / non-coding loci that are rarely
clinically reportable on their own.

---

## Verification

The alias map is regenerated from `gene_alignment_report.tsv` whenever
the panel alignment is re-run:

```bash
python3 scripts/panel_alignment/build_alias_map.py \
    --report  results/panel_alignment/gene_alignment_report_plus.tsv \
    --vep-vcf results/panel_alignment/synthetic_plus.vep.vcf.gz \
    --out     results/panel_alignment/vep_alias_map_plus.tsv
```

A complete rebuild of the alignment (after a panel change, HGNC update,
or VEP cache upgrade) goes:

```bash
# on DGM
bash scripts/panel_alignment/rebuild_alignment_expanded.sh
# back on laptop
python3 scripts/panel_alignment/parse_vep_alignment.py \
    --vep-vcf      results/panel_alignment/synthetic_plus.vep.vcf.gz \
    --panel-hgnc   results/panel_alignment/panel_loose_plus_clinical.hgnc.tsv \
    --out-report   results/panel_alignment/gene_alignment_report_plus.tsv \
    --out-disagree results/panel_alignment/vep_disagreements_plus.tsv
python3 scripts/panel_alignment/build_alias_map.py \
    --report  results/panel_alignment/gene_alignment_report_plus.tsv \
    --vep-vcf results/panel_alignment/synthetic_plus.vep.vcf.gz \
    --out     results/panel_alignment/vep_alias_map_plus.tsv
```

After every rebuild, **re-send `vep_alias_map.tsv` and
`panel_loose_plus_clinical.hgnc_canonical.txt` to the backend** and rerun
their ingest.

---

## Versioning

| Item | Current value |
| --- | --- |
| Panel | `panel_loose_plus_clinical.hgnc_canonical.txt`, 7,363 unique HGNC symbols |
| Panel source databases | OMIM (2026-05-10) + ClinGen (2026-05-10) + GenCC (2026-05-10), HGNC-aligned, ∪ clinician specialty panels (2026-06-12) |
| VEP version | Ensembl VEP 115, `vep_115.sif`, cache at `/home/pipeline/reference/hg38/tertiary/vep_cache/` |
| Alignment date | 2026-06-12 (see `git log results/panel_alignment/gene_alignment_report_plus.tsv`) |
| Alias map rows | 38 (16 `vep_uses_old_symbol` + 12 `vep_missing_gene` + 1 `ensembl_id_drift` + 1 `disagree_other` + 8 `no_vep_record`) |
