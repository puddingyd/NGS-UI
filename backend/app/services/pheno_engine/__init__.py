"""HPO → disease → gene phenotype engine (plan docs/research/NGS_Integrated_Prioritization_Plan_v3.md ch. 4).

Phase 1 (this package so far): the versioned knowledge base — information
content, the pairwise match kernel, per-disease baseline scores, disease–gene
links (platform OMIM, Orphanet product6) and MONDO groups — written next to the
phase-0 cleaned HPO → gene table under ``phenotype_data/hpo_kb/{hp}_{hpoa}/``.
"""
