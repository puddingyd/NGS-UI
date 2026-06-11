# GeneBe DB (`genebe_hg38.tsv.gz`) — requirements for the downstream platform

The NGS-UI tertiary platform consumes `genebe_hg38.tsv.gz` as an **offline
replacement for the live GeneBe API**. To keep it consumable on every rebuild,
the build/publish step must satisfy the following. (Most are also required for
any `tabix`-based consumer, e.g. the `GENEBE_USAGE.md` workflow.)

## 1. No malformed / placeholder rows

Every data row must be a real variant with an **integer `pos`** and a concrete
`chr/ref/alt`. The current build emits rows like:

```
chr.	.	.	.	.	.	.
```

i.e. `chr = "chr."`, `pos = "."`, everything `.`. These appear to come from
`"chr" + str(chrom)` where `chrom` was missing/NaN.

Problems they cause:

- **`tabix` indexing fails**: `tabix -s 1 -b 2 -e 2 -c '#'` parses column 2 as
  the begin coordinate; a `.` there aborts with
  `Failed to parse TBX_GENERIC ... offending line was: "chr.\t.\t..."`.
- They carry no usable information.

**Fix at the source**: drop any record with a missing/non-integer `pos` (and a
missing `chr`) before writing — don't `"chr"`-prefix a missing chromosome.

## 2. Coordinate-sorted

The body must be sorted by `chr` then numeric `pos` (`sort -k1,1 -k2,2n`) so
`tabix` can index it. If it isn't, `tabix` errors with "chromosome blocks not
continuous" / "not sorted".

## 3. `#`-prefixed header, stable column **names**

Keep the header line starting with `#`, with these column **names** present
(order may vary; downstream reads them by name, not position):

```
#chr  pos  ref  alt  acmg_classification  acmg_score  acmg_criteria
```

- `chr` is `chr`-prefixed (`chr1`…`chr22`, `chrX`, `chrY`, `chrM`).
- `pos` is 1-based integer; `ref`/`alt` are uppercase `A/C/G/T…`, multi-allelic
  already split (one ALT per row), indels left-aligned + trimmed.
- `acmg_score` integer; `acmg_criteria` comma-separated (e.g. `PVS1,PM2,PP3`);
  missing values are a single `.`.

The slim 7-column layout above is fine; a full 55-column layout is fine too.

## 4. bgzip + fresh tabix index, published atomically

- Compress with **`bgzip`** (not plain `gzip`).
- Re-create the index every rebuild: `tabix -s 1 -b 2 -e 2 -c '#' file.tsv.gz`,
  and make sure the resulting `.tbi` is **newer than** the `.gz`. A `.tbi`
  older than its `.gz` makes `tabix` read stale byte offsets and fail with
  `Invalid BGZF header at offset …` (this is *not* file corruption — it's a
  stale index).
- **Publish atomically**: write to a temp path, `bgzip`, `tabix`, then `mv` the
  `.gz` and `.tbi` into place (same filesystem) so a consumer never sees a
  half-written file or a `.gz`/`.tbi` mismatch.

## 5. Self-check before publishing

```bash
bgzip -t genebe_hg38.tsv.gz                      # 0 = intact
zcat genebe_hg38.tsv.gz | awk -F'\t' 'NR>1 && $2 !~ /^[0-9]+$/' | head   # must be empty
tabix genebe_hg38.tsv.gz chr7:140753336-140753336   # one+ clean rows, no W/E
```

> NGS-UI streams the `.gz` directly and skips any non-integer-`pos` row
> defensively, so it won't crash on a dirty DB — but cleaning at the source
> keeps the DB correct for `tabix` consumers and avoids needless misses.
