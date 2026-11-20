"""Gene-list exclusions for WES diagnostic exports only.

These rules describe the lab's report scope, not variant filtering or a
change to the underlying HPO/panel definitions. Inputs use HGNC symbols.
"""
from __future__ import annotations

from . import test_types


WES_EXCLUDED_GENES = frozenset({
    "SMN1", "SMN2", "STRC", "CYP21A2", "NCF1", "IKBKG",
    "OPN1LW", "OPN1MW", "MUC1",
})


def include_in_gene_list(gene: str, test_type: str) -> bool:
    if test_types.normalize_test_type(test_type) != test_types.WES:
        return True
    symbol = str(gene or "").strip().upper()
    return not (symbol.startswith("MT-") or symbol in WES_EXCLUDED_GENES)
