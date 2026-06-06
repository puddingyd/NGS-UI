# Using `genebe_hg38.tsv.gz` for annotation (backend integration)

This guide is for the **downstream analysis platform** that needs to attach
GeneBe ACMG + ClinVar + gnomAD + in-silico annotations to variants. It does
**not** cover building the database — see `README.md` for that.

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
- **Header line** starts with `#` (tabix is told to skip it with `-c '#'`).
- **Missing values** are `.` (a single dot), not empty string.
- Two columns are **JSON-encoded** (compact, no spaces): `acmg_by_gene` and
  `custom_annotations`. All others are plain scalars.
- Tabix index: sequence = col 1, begin = col 2, end = col 2.

### Columns (in order)

| # | Column | # | Column |
|---|--------|---|--------|
| 1 | `chr` | 29 | `frequency_reference_population` |
| 2 | `pos` | 30 | `allele_count_reference_population` |
| 3 | `ref` | 31 | `hom_count_reference_population` |
| 4 | `alt` | 32 | `dbsnp` |
| 5 | `gene_symbol` | 33 | `computational_prediction_selected` |
| 6 | `gene_hgnc_id` | 34 | `computational_score_selected` |
| 7 | `transcript` | 35 | `computational_source_selected` |
| 8 | `hgvs_c` | 36 | `revel_prediction` |
| 9 | `consequences` | 37 | `revel_score` |
| 10 | `effect` | 38 | `alphamissense_prediction` |
| 11 | `acmg_classification` | 39 | `alphamissense_score` |
| 12 | `acmg_score` | 40 | `bayesdelnoaf_prediction` |
| 13 | `acmg_criteria` | 41 | `bayesdelnoaf_score` |
| 14 | `acmg_by_gene` (JSON) | 42 | `phylop100way_prediction` |
| 15 | `clinvar_classification` | 43 | `phylop100way_score` |
| 16 | `clinvar_disease` | 44 | `splice_prediction_selected` |
| 17 | `clinvar_review_status` | 45 | `splice_score_selected` |
| 18 | `clinvar_submissions_summary` | 46 | `splice_source_selected` |
| 19 | `pathogenicity_classification_combined` | 47 | `spliceai_max_prediction` |
| 20 | `phenotype_combined` | 48 | `spliceai_max_score` |
| 21 | `gnomad_exomes_af` | 49 | `dbscsnv_ada_prediction` |
| 22 | `gnomad_exomes_ac` | 50 | `dbscsnv_ada_score` |
| 23 | `gnomad_exomes_homalt` | 51 | `apogee2_prediction` |
| 24 | `gnomad_genomes_af` | 52 | `apogee2_score` |
| 25 | `gnomad_genomes_ac` | 53 | `mitotip_prediction` |
| 26 | `gnomad_genomes_homalt` | 54 | `mitotip_score` |
| 27 | `gnomad_mito_heteroplasmic` | 55 | `custom_annotations` (JSON) |
| 28 | `gnomad_mito_homoplasmic` | | |

> The canonical, authoritative order lives in `02_annotate.py`
> (`CANONICAL_COLUMNS`). Don't hard-code column numbers in long-lived code —
> read the `#` header row once and build a name→index map.

---

## 2. Prerequisites on the backend

`htslib` / `tabix` must be installed.

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

# A genomic region (returns every annotated variant in it)
tabix "$DB" chr1:925000-926000

# A single position (then filter ref/alt yourself)
tabix "$DB" chr7:140753336-140753336

# Keep the header for readability
tabix -h "$DB" chr17:43044295-43125483
```

### B. Point lookup of an exact variant (Python / pysam)

This is the core function the platform should wrap. Given a normalized
`(chr, pos, ref, alt)`, return the annotation row or `None`.

```python
import pysam

class GeneBeDB:
    def __init__(self, path="output/genebe_hg38.tsv.gz"):
        self.tbx = pysam.TabixFile(path)
        # Build a name -> index map from the '#'-prefixed header line.
        header = next(l for l in self.tbx.header)  # e.g. "#chr\tpos\t..."
        self.cols = header.lstrip("#").split("\t")
        self.idx  = {name: i for i, name in enumerate(self.cols)}

    def get(self, chrom, pos, ref, alt):
        """Exact match on (chr, pos, ref, alt). Returns dict or None."""
        if not chrom.startswith("chr"):
            chrom = "chr" + chrom
        try:
            rows = self.tbx.fetch(chrom, pos - 1, pos)  # 0-based half-open
        except ValueError:
            return None  # contig not in index
        for line in rows:
            f = line.split("\t")
            if f[1] == str(pos) and f[2] == ref and f[3] == alt:
                return dict(zip(self.cols, f))
        return None

db = GeneBeDB()
hit = db.get("chr7", 140753336, "A", "T")   # BRAF V600E
if hit:
    print(hit["acmg_classification"], hit["clinvar_classification"])
```

### C. Annotating a whole VCF / variant list

For batch jobs, fetch each variant's region once and match in memory. To
annotate millions of variants efficiently, sort your input by coordinate and
stream — tabix random-access is fast, but sequential locality helps the OS
page cache.

```python
import json

def annotate_vcf(in_vcf, out_tsv, db):
    want = ["acmg_classification", "acmg_score",
            "clinvar_classification", "gnomad_genomes_af",
            "revel_score", "spliceai_max_score"]
    with open(in_vcf) as fi, open(out_tsv, "w") as fo:
        fo.write("chr\tpos\tref\talt\t" + "\t".join(want) + "\n")
        for line in fi:
            if line.startswith("#"):
                continue
            c, p, _id, ref, alt = line.split("\t")[:5]
            for a in alt.split(","):                 # split multi-allelic
                hit = db.get(c, int(p), ref, a)
                vals = [hit[k] for k in want] if hit else ["."] * len(want)
                fo.write(f"{c}\t{p}\t{ref}\t{a}\t" + "\t".join(vals) + "\n")

# JSON columns, when you need them:
#   by_gene = json.loads(hit["acmg_by_gene"])   # {} / "." when absent
```

---

## 4. Variant normalization — **read this before matching**

A lookup only hits if your `(chr, pos, ref, alt)` is byte-identical to what's
stored. The database was built from dbNSFP / ClinVar / gnomAD coordinates, so
normalize your query variants the same way **before** calling `get()`:

1. **`chr` prefix** — DB rows are `chr1`…`chrX`. Add the prefix if your source
   (e.g. plain Ensembl VCF) omits it. (`get()` above does this for you.)
2. **Left-align + trim indels** — use `bcftools norm -f ref.fa` (or equivalent)
   so indels are in canonical minimal representation. This is the #1 cause of
   "variant exists but I get a miss".
3. **Split multi-allelic** — one ALT per query (`bcftools norm -m -`). DB rows
   are already split.
4. **Build** = GRCh38/hg38. Liftover hg19 → hg38 first if needed
   (`genebe_hg19.R` in this repo is for the hg19 flow).
5. **Case / symbols** — SNV ref/alt are uppercase `A/C/G/T`. Don't pass `N`,
   IUPAC codes, or symbolic ALTs (`<DEL>`); those are never in the DB.

Recommended one-liner to canonicalize an input VCF first:

```bash
bcftools norm -m - -f GRCh38.fa input.vcf.gz -Oz -o input.norm.vcf.gz
```

---

## 5. What a "miss" means (and the fallback)

This DB is **pre-computed**, not exhaustive. A variant is present only if it came
from one of the build phases:

- **dbNSFP** — essentially all coding non-synonymous + canonical splice SNVs
- **ClinVar** — clinically reported variants (incl. indels, some non-coding)
- **gnomAD v4.1 joint** — rare observed variants (global AF < 0.01)

So expect **misses** for: novel variants, most deep-intronic / intergenic
positions, and indels never seen in ClinVar/gnomAD. A miss is **not** an error —
it just means "not pre-annotated."

**Recommended backend policy:**

1. Look up in `genebe_hg38.tsv.gz` (fast, offline, free).
2. On a miss, optionally fall back to the **live GeneBe API** for that one
   variant and cache the result:

   ```python
   import genebe as gnb
   df = gnb.annotate(["7-140753336-A-T"], genome="hg38",
                     username="...", api_key="...", use_netrc=False,
                     output_format="dataframe")
   ```

   This keeps coverage at 100% while serving the common case from the local file.

---

## 6. Operational notes

- **The file is rebuilt daily** (after each `02_annotate.py` run) and grows as
  more chunks finish. It is **usable while partial** — coverage simply increases
  over time. Until the build completes, expect more misses, especially for
  genome/ClinVar/gnomAD-only variants.
- **Atomic swap**: `03_merge.sh` writes then `bgzip`s in place. If the platform
  keeps a long-lived `pysam.TabixFile` handle open, **reopen it** after a
  rebuild (e.g. detect a new mtime) so you don't read a half-written file.
- **`.tbi` must travel with the `.gz`** and be newer than it. If you copy the DB
  elsewhere, copy both files, or re-run `tabix -s 1 -b 2 -e 2 -c '#' file.tsv.gz`.
- **Version pinning**: record the `.gz` mtime / a checksum alongside your
  annotation results so you can reproduce exactly which DB snapshot was used.
