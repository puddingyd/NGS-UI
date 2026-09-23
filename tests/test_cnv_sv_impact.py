"""Display impact must not create false coding hits or remove source variants."""
import csv

import pytest

from app.services.cnv_sv_impact import attach, gene_impact, summarize
from app.adapters.annotsv_tsv import load_annotsv_tsv


def gene(**kwargs):
    return {"gene": "TEST", "location": "intron2-intron2", "location2": "CDS",
            "overlap_cds_len": 0, "splice_distance": 100, "splice_type": "donor",
            "pheno_score": 50, "hpo_score": 30, **kwargs}


@pytest.mark.parametrize("kind", ["DEL", "DUP"])
@pytest.mark.parametrize("changes, category", [
    ({}, "noncoding"),
    ({"overlap_cds_len": 1}, "functional"),
    ({"splice_distance": 2}, "functional"),
    ({"splice_distance": 2, "splice_type": "NA"}, "unknown"),
    ({"splice_distance": 3}, "noncoding"),
    ({"splice_distance": None}, "unknown"),
    ({"overlap_cds_len": None}, "unknown"),
    ({"location": "exon1-exon1", "location2": "5'UTR"}, "noncoding"),
    ({"location": "intron2-intron4", "overlap_cds_len": None}, "unknown"),
    ({"location": "", "location2": ""}, "unknown"),
])
def test_location_classification(kind, changes, category):
    assert gene_impact(gene(**changes), kind)[0] == category


@pytest.mark.parametrize("kind", ["INV", "BND", "TRA", "INS", "CPX"])
def test_balanced_complex_spans_do_not_imply_exon_disruption(kind):
    assert gene_impact(gene(overlap_cds_len=100), kind)[0] == "unknown"


def test_only_matched_gene_can_drive_clinical_impact():
    variant = {"sv_type": "DEL", "genes": [gene(), gene(gene="OTHER", pheno_score=0,
                                                         overlap_cds_len=50, hpo_score=0)]}
    attach(variant)
    assert variant["impact_clinical"]["category"] == "noncoding"
    assert variant["impact_all"]["category"] == "functional"
    assert variant["impact_all"]["hpo_score"] == 0


def test_unknown_relevant_gene_prevents_hiding_alongside_noncoding_hit():
    assert summarize([gene(), gene(overlap_cds_len=None)], "DEL")["category"] == "unknown"


def test_directional_dosage_is_not_a_generic_pathogenicity_score():
    g = gene(overlap_cds_len=100, hi=3, ts=0)
    assert summarize([g], "DEL")["mechanism"] == 1
    assert summarize([g], "DUP")["mechanism"] == 0
    g["ts"] = 3
    assert summarize([g], "DUP")["mechanism"] == 0  # partial DUP
    g["location"] = "txStart-txEnd"
    assert summarize([g], "DUP")["mechanism"] == 1
    assert summarize([gene(hi=30)], "DEL")["mechanism"] == 0


@pytest.mark.parametrize("source", ["cnv", "sv"])
def test_adapter_computes_before_trim_and_preserves_out_of_order_split_rows(tmp_path, monkeypatch, source):
    from app.services import gene_disease_store, panel_deadzone
    monkeypatch.setattr(gene_disease_store, "ensure_loaded", lambda: None)
    monkeypatch.setattr(gene_disease_store, "merged_associations", lambda *a, **kw: [])
    monkeypatch.setattr(panel_deadzone, "canonical_gene_symbol", lambda g: (g, None))
    headers = ["AnnotSV_ID", "Annotation_mode", "SV_chrom", "SV_start", "SV_end", "SV_type",
               "Gene_name", "ACMG_class", "Location", "Location2", "Overlapped_CDS_length",
               "Dist_nearest_SS", "Nearest_SS_type", "HI", "TS", "AnnotSV_ranking_score"]
    path = tmp_path / "annotated.tsv"
    base = dict(AnnotSV_ID="test", SV_chrom="1", SV_start="100", SV_end="200", SV_type="DEL")
    genes = ["G" + str(i) for i in range(12)]
    rows = [{**base, "Annotation_mode": "split", "Gene_name": name,
             "Location": "intron2-intron2", "Location2": "CDS",
             "Overlapped_CDS_length": "10" if name == "G11" else "0",
             "Dist_nearest_SS": "100", "Nearest_SS_type": "donor", "HI": "3"} for name in genes]
    rows.append({**base, "Annotation_mode": "full", "ACMG_class": "3", "AnnotSV_ranking_score": "0.1"})
    with path.open("w") as f:
        writer = csv.DictWriter(f, fieldnames=headers, delimiter="\t")
        writer.writeheader(); writer.writerows(rows)
    variants, categories = load_annotsv_tsv(path, source=source,
        pheno_by_gene={g: 20 for g in genes}, pheno_matched={g: 1 for g in genes},
        hpo_by_gene={"G11": 70})
    v = variants["test"]
    assert v["genes_total"] == 12
    assert v["impact_clinical"]["category"] == "functional"
    assert v["impact_clinical"]["hpo_score"] == 70
    assert v["impact_clinical"]["reasons"][0]["gene"] == "G11"
    assert v["acmg_class"] == 3
    assert categories["CNV-1A" if source == "cnv" else "SV-2A"] == ["test"]
    assert categories["CNV-1B" if source == "cnv" else "SV-2B"] == []


def test_backend_parent_preserves_segment_scope():
    from app.services.cnv_sv_merge import build_parent
    first = {"id": "a", "source": "cnv", "CHROM": "1", "POS": 100, "END": 200,
             "sv_type": "DEL", "in_panel": True, "cnv_sv_sort_score": 1, "genes": [gene()]}
    second = {**first, "id": "b", "POS": 300, "END": 400, "in_panel": False,
              "genes": [gene(pheno_score=0, overlap_cds_len=100)]}
    attach(first); attach(second)
    parent = build_parent({"member_ids": ["a", "b"]}, {"a": first, "b": second})
    assert parent["impact_clinical"]["category"] == "noncoding"
    assert parent["impact_all"]["category"] == "functional"


@pytest.mark.parametrize("loader_name", ["load_sample_cnv", "load_sample_sv", "load_sample_cnv_sv"])
def test_staged_loaders_pass_hpo_only_scores_without_panel_inflation(tmp_path, monkeypatch, loader_name):
    from app.services import sample_loader, phenotype_scorer, sample_layout
    import app.adapters.annotsv_tsv as adapter
    path = tmp_path / "annotation.tsv"
    path.write_text("fixture")
    monkeypatch.setattr(phenotype_scorer, "_LOADED", True)
    monkeypatch.setattr(phenotype_scorer, "_HPO_TO_GENES", {"HP:1": {"HPOGENE"}, "PANEL": {"PANELGENE"}})
    monkeypatch.setattr(phenotype_scorer, "_PANEL_TO_GENES", {"PANEL": {"PANELGENE"}})
    monkeypatch.setattr(sample_loader, "_load_pheno_context", lambda *a: (
        tmp_path, tmp_path, [{"hpo_id": "HP:1", "weight": 1}], ["PANEL"],
        {"HPOGENE": 50, "PANELGENE": 50}))
    monkeypatch.setattr(sample_layout, "cnv_tsv", lambda _: path)
    monkeypatch.setattr(sample_layout, "sv_tsv", lambda _: path)
    observed = []
    def load(path, **kwargs):
        observed.append(kwargs)
        return {}, {}
    monkeypatch.setattr(adapter, "load_annotsv_tsv", load)
    getattr(sample_loader, loader_name)("SAMPLE")
    assert observed
    for args in observed:
        assert args["hpo_by_gene"] == {"HPOGENE": 100}
        assert args["pheno_by_gene"]["PANELGENE"] == 50
