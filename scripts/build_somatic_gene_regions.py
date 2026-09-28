#!/usr/bin/env python3
"""Build the five-column gene/exon BED from a matching GRCh38 GTF (all transcripts)."""
import argparse
import gzip
import re
from pathlib import Path


def build(gtf: Path, output: Path, release: str):
    opener = gzip.open if gtf.suffix == ".gz" else open
    rows, genes = set(), {}
    with opener(gtf, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            p = line.rstrip().split("\t")
            if len(p) != 9 or p[2] not in {"gene", "exon"}:
                continue
            attrs = dict(re.findall(r'(\w+) "([^"]+)"', p[8]))
            name = attrs.get("gene_name")
            chrom = p[0] if p[0].startswith("chr") else "chr" + p[0]
            if not name or chrom not in {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}:
                continue
            start, end = int(p[3]) - 1, int(p[4])
            if start < 0 or end <= start:
                raise ValueError("Invalid GTF interval")
            if p[2] == "exon":
                rows.add((chrom, start, end, name, "exon"))
            key = chrom, name
            old = genes.get(key, (start, end))
            genes[key] = min(start, old[0]), max(end, old[1])
    rows.update((chrom, start, end, name, "gene") for (chrom, name), (start, end) in genes.items())
    if not rows:
        raise ValueError("No primary nuclear gene/exon records in GTF")
    with output.open("w") as out:
        out.write(f"#assembly=GRCh38 release={release}\n")
        for row in sorted(rows):
            out.write("\t".join(map(str, row)) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gtf", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--release", required=True)
    args = parser.parse_args()
    build(args.gtf, args.out, args.release)
