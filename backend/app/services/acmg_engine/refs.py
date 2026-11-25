"""Reference data the rule engine consults (plan section 5.2.11).

The engine only talks to a :class:`References` object.  Every lookup returns
``None`` when the underlying table is not loaded, which the engine reports as
"not evaluated" rather than "not met".  :class:`StaticReferences` holds plain
in-memory tables; file loaders (step 2) will build it from the versioned
reference files.
"""
from __future__ import annotations

from dataclasses import dataclass, field


class References:
    """Interface with every table missing."""

    version: dict = {}
    available: frozenset = frozenset()

    def lof_mechanism(self, gene: str) -> tuple[str, str] | None:
        """("yes" | "no", reason) or None when unknown/unloaded."""
        return None

    def truncating_downstream(self, gene: str, protein_pos: int) -> bool | None:
        """ClinVar P/LP truncating variant downstream of this protein position."""
        return None

    def clinvar_missense(self, gene: str, protein_pos: int) -> list[dict] | None:
        """ClinVar P/LP missense at this residue: [{alt_aa, variant_key, stars}]."""
        return None

    def clinvar_missense_window(self, gene: str, start: int, end: int) -> list[dict] | None:
        """ClinVar P/LP missense between two residues: [{protein_pos, alt_aa, variant_key, stars, category}]."""
        return None

    def benign_missense_window(self, gene: str, start: int, end: int) -> list[dict] | None:
        """ClinVar B/LB missense between two residues: [{protein_pos, variant_key, stars}]."""
        return None

    def missense_gene_counts(self, gene: str) -> tuple[int, int] | None:
        """(P/LP, B/LB) ≥1★ missense variants of the gene on its MANE transcript."""
        return None

    def exon_skip(self, transcript: str, exon_number: int) -> dict | None:
        """{"frameshift": bool, "fraction": float, "last_exon": bool} for skipping an exon."""
        return None

    def cds_exons(self, transcript: str) -> list[int] | None:
        """CDS length of every exon of the MANE transcript, 5'→3' (UTR-only exons are 0)."""
        return None

    def ba1_exception(self, variant_key: str) -> bool:
        return False

    def wes_recurrence(self, variant_key: str) -> tuple[int, int] | None:
        """(patients carrying the variant, patients in the WES table)."""
        return None

    def vcep(self, gene: str) -> dict | None:
        """Gene-specific ClinGen VCEP specification (cspec.build_vcep_specs entry):
        ``criteria`` {code: allowed strengths}, ``thresholds`` {ba1, bs1, pm2},
        ``pm1`` [{strength, ranges}], ``pm1_unparsed``; older ``bs1`` /
        ``pm1_regions`` keys are still honoured."""
        return None

    def collagen_gly(self, gene: str) -> dict | None:
        """{"accession", "length", "positions"} of Gly-X-Y triple-helix glycines."""
        return None

    def inhouse_probands(self, variant_key: str) -> int | None:
        """Other unrelated patients whose case marked this variant Causative (PS4 A)."""
        return None

    def clinvar_affected(self, variation_id: str) -> dict | None:
        """ClinVar affected-case evidence for a Variation ID (PS4 B):
        {"submitters": n, "individuals": m}; ``{}`` when ClinVar has none."""
        return None

    def clinvar_biallelic(self, variation_id: str) -> dict | None:
        """ClinVar submitters describing homozygous / compound-heterozygous patients (PM3):
        {"submitters", "hom", "comphet", "detail": [{"scv", "snippet"}]}; ``{}`` when none."""
        return None


@dataclass
class StaticReferences(References):
    """In-memory tables; any table left as ``None`` behaves as not loaded."""

    lof_genes: dict | None = None            # gene -> ("yes"|"no", reason)
    truncating_positions: dict | None = None  # gene -> sorted protein positions
    missense_by_residue: dict | None = None   # (gene, pos) -> [ {alt_aa, variant_key, stars, category?} ]
    benign_missense: dict | None = None       # gene -> [ {protein_pos, variant_key, stars} ]
    missense_genes: dict | None = None        # gene -> (P/LP count, B/LB count)
    exons: dict | None = None                 # (transcript, exon) -> {frameshift, fraction, last_exon}
    ba1_exceptions: set = field(default_factory=set)
    recurrence: dict | None = None            # variant_key -> patients
    recurrence_patients: int = 0
    vcep_specs: dict = field(default_factory=dict)
    collagen: dict = field(default_factory=dict)     # gene -> {"length", "positions"}
    probands: dict | None = None             # variant_key -> other unrelated Causative patients
    clinvar_cases: dict | None = None        # variation_id -> {"submitters", "individuals"}
    biallelic: dict | None = None            # variation_id -> {"submitters", "hom", "comphet", "detail"}
    version: dict = field(default_factory=dict)

    @property
    def available(self) -> frozenset:  # type: ignore[override]
        tables = {
            "lof_genes": self.lof_genes,
            "truncating_positions": self.truncating_positions,
            "missense_by_residue": self.missense_by_residue,
            "exons": self.exons,
            "recurrence": self.recurrence,
            "missense_window": self.benign_missense if self.missense_by_residue is not None else None,
            "missense_genes": self.missense_genes,
        }
        return frozenset(name for name, value in tables.items() if value is not None)

    def lof_mechanism(self, gene):
        if self.lof_genes is None:
            return None
        return self.lof_genes.get(gene)

    def truncating_downstream(self, gene, protein_pos):
        if self.truncating_positions is None:
            return None
        return any(pos > protein_pos for pos in self.truncating_positions.get(gene, ()))

    def clinvar_missense(self, gene, protein_pos):
        if self.missense_by_residue is None:
            return None
        return list(self.missense_by_residue.get((gene, protein_pos), ()))

    def clinvar_missense_window(self, gene, start, end):
        if self.missense_by_residue is None or self.benign_missense is None:
            return None
        return [dict(e, protein_pos=pos) for (g, pos), entries in self.missense_by_residue.items()
                if g == gene and start <= pos <= end for e in entries]

    def benign_missense_window(self, gene, start, end):
        if self.benign_missense is None:
            return None
        return [e for e in self.benign_missense.get(gene, ()) if start <= e["protein_pos"] <= end]

    def missense_gene_counts(self, gene):
        if self.missense_genes is None:
            return None
        return tuple(self.missense_genes.get(gene, (0, 0)))

    def exon_skip(self, transcript, exon_number):
        if self.exons is None:
            return None
        return self.exons.get((transcript.split(".")[0], exon_number))

    def cds_exons(self, transcript):
        if self.exons is None or not transcript:
            return None
        tx = str(transcript).split(".")[0]
        hits = sorted((n, v.get("cds_length", 0)) for (t, n), v in self.exons.items() if t == tx)
        return [int(length) for _n, length in hits] or None

    def ba1_exception(self, variant_key):
        return variant_key in self.ba1_exceptions

    def wes_recurrence(self, variant_key):
        if self.recurrence is None:
            return None
        return self.recurrence.get(variant_key, 0), self.recurrence_patients

    def vcep(self, gene):
        return self.vcep_specs.get(gene)

    def collagen_gly(self, gene):
        return self.collagen.get(gene)

    def inhouse_probands(self, variant_key):
        if self.probands is None:
            return None
        return int(self.probands.get(variant_key, 0))

    def clinvar_affected(self, variation_id):
        if self.clinvar_cases is None:
            return None
        if not variation_id:
            return {}
        return dict(self.clinvar_cases.get(str(variation_id), {}))

    def clinvar_biallelic(self, variation_id):
        if self.biallelic is None:
            return None
        if not variation_id:
            return {}
        return dict(self.biallelic.get(str(variation_id), {}))
