# Using `genebe_hg38.tsv.gz` for annotation (backend integration)

This guide is for the **downstream analysis platform** that needs to look up
pre-computed ACMG classifications for hg38 variants. It does **not** cover
building the database — see `README.md` for that.

The database is a single **bgzip-compressed, tabix-indexed TSV**:

```
output/genebe_hg38.tsv.gz       # the data
output/genebe_hg38.tsv.gz.tbi   # the index (must sit next to it)
```

---

## 1. File format

- **Coordinates**: hg38 / GRCh38, 1-based, `chr`-prefixed (`chr1`, `chrX`, …).
- **One row per variant**, keyed by `(chr, pos, ref, alt)`. Multi-allelic sites
  are already **split** — every ALT is its own row.
- **Header line** starts with `#` (tabix skips it via `-c '#'`).
- **Missing values** are `.` (a single dot), not empty string.
- Tabix index: sequence = col 1, begin = col 2, end = col 2.

### Columns

| # | Column | Type | Description |
|---|--------|------|-------------|
| 1 | `chr` | string | Chromosome (`chr1` … `chr22`, `chrX`, `chrY`) |
| 2 | `pos` | integer | 1-based position |
| 3 | `ref` | string | Reference allele |
| 4 | `alt` | string | Alternate allele |
| 5 | `acmg_classification` | string | e.g. `Pathogenic`, `Likely_pathogenic`, `VUS`, `Likely_benign`, `Benign` |
| 6 | `acmg_score` | integer | GeneBe ACMG point score (positive = pathogenic direction) |
| 7 | `acmg_criteria` | string | Comma-separated criteria, e.g. `PVS1,PM2,PP3` |

> Full per-variant annotations (ClinVar, gnomAD AF, REVEL, SpliceAI, AlphaMissense,
> etc.) remain in the individual `genebe_results/` chunk TSVs. See `02_annotate.py`
> (`CANONICAL_COLUMNS`) for the complete 55-column schema.

---

## 2. Prerequisites on the backend

`htslib` / `tabix` must be installed:

```bash
# macOS
brew install htslib
# Debian/Ubuntu
apt-get install tabix
# conda
conda install -c bioconda htslib
```

Python (`pysam`) is the most ergonomic for programmatic use:

```bash
pip install pysam
```

---

## 3. Querying

### A. Single-variant / region lookup (CLI)

```bash
DB=output/genebe_hg38.tsv.gz

# A genomic region
tabix "$DB" chr1:925000-926000

# A single position (then filter ref/alt yourself)
tabix "$DB" chr7:140753336-140753336

# Include the header line
tabix -h "$DB" chr17:43044295-43125483
```

### B. Point lookup of an exact variant (Python / pysam)

This is the core function the platform should wrap. Given a normalized
`(chr, pos, ref, alt)`, return the annotation dict or `None`.

```python
import pysam

class GeneBeDB:
    def __init__(self, path="output/genebe_hg38.tsv.gz"):
        self.tbx = pysam.TabixFile(path)
        # Read column names from the '#'-prefixed header line.
        header = next(iter(self.tbx.header))   # "#chr\tpos\t..."
        self.cols = header.lstrip("#").split("\t")

    def get(self, chrom, pos, ref, alt):
        """Exact match on (chr, pos, ref, alt). Returns dict or None."""
        if not chrom.startswith("chr"):
            chrom = "chr" + chrom
        try:
            rows = self.tbx.fetch(chrom, pos - 1, pos)   # 0-based half-open
        except ValueError:
            return None   # contig not in index
        for line in rows:
            f = line.split("\t")
            if f[1] == str(pos) and f[2] == ref and f[3] == alt:
                return dict(zip(self.cols, f))
        return None

db  = GeneBeDB()
hit = db.get("chr7", 140753336, "A", "T")   # BRAF V600E
if hit:
    print(hit["acmg_classification"])   # e.g. "Pathogenic"
    print(hit["acmg_score"])            # e.g. "10"
    print(hit["acmg_criteria"])         # e.g. "PVS1,PM2,PP3"
```

### C. Annotating a whole VCF / variant list (batch)

Sort your input by coordinate before querying — tabix random-access is fast,
but sequential locality helps the OS page cache for large inputs.

```python
def annotate_vcf(in_vcf, out_tsv, db):
    out_cols = ["acmg_classification", "acmg_score", "acmg_criteria"]
    with open(in_vcf) as fi, open(out_tsv, "w") as fo:
        fo.write("chr\tpos\tref\talt\t" + "\t".join(out_cols) + "\n")
        for line in fi:
            if line.startswith("#"):
                continue
            c, p, _id, ref, alt = line.split("\t")[:5]
            for a in alt.split(","):           # split multi-allelic
                hit = db.get(c, int(p), ref, a)
                vals = [hit[k] for k in out_cols] if hit else ["."] * len(out_cols)
                fo.write(f"{c}\t{p}\t{ref}\t{a}\t" + "\t".join(vals) + "\n")
```

---

## 4. Variant normalization — **read this before matching**

A lookup only hits if your `(chr, pos, ref, alt)` is byte-identical to what's
stored. Normalize your query variants **before** calling `get()`:

1. **`chr` prefix** — DB rows are `chr1`…`chrX`. Add the prefix if your source
   omits it. (`get()` above does this for you.)
2. **Left-align + trim indels** — use `bcftools norm -f ref.fa` (or equivalent).
   This is the **#1 cause** of "variant exists but I get a miss".
3. **Split multi-allelic** — one ALT per query (`bcftools norm -m -`). DB rows
   are already split.
4. **Build** = GRCh38/hg38. Liftover hg19 → hg38 first if needed.
5. **Case / symbols** — ref/alt are uppercase `A/C/G/T`. Symbolic ALTs (`<DEL>`)
   are never in the DB.

Recommended one-liner to canonicalize an input VCF:

```bash
bcftools norm -m - -f GRCh38.fa input.vcf.gz -Oz -o input.norm.vcf.gz
```

---

## 5. What a "miss" means (and the fallback)

This DB is **pre-computed**, not exhaustive. A variant is present only if it
came from one of the build phases:

- **dbNSFP** — essentially all coding non-synonymous + canonical splice SNVs
- **ClinVar** — clinically reported variants (incl. indels, some non-coding)
- **gnomAD v4.1 joint** — rare observed variants (global AF < 0.01)

So expect misses for: novel variants, most deep-intronic / intergenic positions,
and indels never seen in ClinVar/gnomAD. A miss is **not** an error.

**Recommended backend policy:**

1. Look up in `genebe_hg38.tsv.gz` (fast, offline, free).
2. On a miss, fall back to the **live GeneBe API** and cache the result:

   ```python
   import genebe as gnb
   df = gnb.annotate(["7-140753336-A-T"], genome="hg38",
                     username="...", api_key="...", use_netrc=False,
                     output_format="dataframe")
   hit = df.iloc[0] if df is not None and len(df) else None
   ```

---

## 6. Operational notes

- **Grows daily**: the file is rebuilt every night as more chunks finish. It is
  usable while partial — coverage increases over time.
- **Atomic swap**: `03_merge.sh` writes to a `.tmp` file and `mv`s it in place,
  so consumers never read a half-written file. If your platform keeps a
  long-lived `pysam.TabixFile` handle open, **reopen it** after detecting a
  new mtime.
- **`.tbi` must travel with the `.gz`**. If you copy the DB elsewhere, copy both
  files, or re-run `tabix -s 1 -b 2 -e 2 -c '#' file.tsv.gz`.
- **Version pinning**: record the `.gz` mtime alongside your annotation results
  so you can reproduce which DB snapshot was used.
