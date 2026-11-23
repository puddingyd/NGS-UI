"""NGS-UI ACMG rule engine (plan docs/research/NGS_Integrated_Prioritization_Plan_v3.md §5.2).

Replaces the pipeline ``ACMG_*`` columns as the platform's third ACMG source
(manual → GeneBe → rule engine).  Pure functions; reference tables are
injected through :class:`References`.
"""
from .engine import Context, evaluate_variant, load_rulebook
from .refs import References, StaticReferences

__all__ = ["Context", "References", "StaticReferences", "evaluate_variant", "load_rulebook"]
