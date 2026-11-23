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

    def exon_skip(self, transcript: str, exon_number: int) -> dict | None:
        """{"frameshift": bool, "fraction": float, "last_exon": bool} for skipping an exon."""
        return None

    def ba1_exception(self, variant_key: str) -> bool:
        return False

    def wes_recurrence(self, variant_key: str) -> tuple[int, int] | None:
        """(patients carrying the variant, patients in the WES table)."""
        return None

    def vcep(self, gene: str) -> dict | None:
        """Gene-specific VCEP specification (thresholds, PM1 regions)."""
        return None


@dataclass
class StaticReferences(References):
    """In-memory tables; any table left as ``None`` behaves as not loaded."""

    lof_genes: dict | None = None            # gene -> ("yes"|"no", reason)
    truncating_positions: dict | None = None  # gene -> sorted protein positions
    missense_by_residue: dict | None = None   # (gene, pos) -> [ {alt_aa, variant_key, stars} ]
    exons: dict | None = None                 # (transcript, exon) -> {frameshift, fraction, last_exon}
    ba1_exceptions: set = field(default_factory=set)
    recurrence: dict | None = None            # variant_key -> patients
    recurrence_patients: int = 0
    vcep_specs: dict = field(default_factory=dict)
    version: dict = field(default_factory=dict)

    @property
    def available(self) -> frozenset:  # type: ignore[override]
        tables = {
            "lof_genes": self.lof_genes,
            "truncating_positions": self.truncating_positions,
            "missense_by_residue": self.missense_by_residue,
            "exons": self.exons,
            "recurrence": self.recurrence,
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

    def exon_skip(self, transcript, exon_number):
        if self.exons is None:
            return None
        return self.exons.get((transcript.split(".")[0], exon_number))

    def ba1_exception(self, variant_key):
        return variant_key in self.ba1_exceptions

    def wes_recurrence(self, variant_key):
        if self.recurrence is None:
            return None
        return self.recurrence.get(variant_key, 0), self.recurrence_patients

    def vcep(self, gene):
        return self.vcep_specs.get(gene)
