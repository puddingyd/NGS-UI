from app.services import (
    docx_export,
    gene_disease_store,
    litvar2_on_demand,
    omim_store,
    sample_loader,
)


def _omim_row_with_disease16() -> dict:
    row = omim_store._empty_row()
    row.update({
        "OMIM_id": "120140",
        "OMIM_disease": "Late disease (654321)(AD)",
        "Inheritance": "AD",
        "Disease16": "Late disease (654321)(AD)\n\nCurator detail.",
    })
    return row


def test_old_five_slot_workbook_row_remains_compatible():
    headers = (
        "OMIM_id", "gene_symbol", "OMIM_disease", "Inheritance",
        "Disease1", "Disease2", "Disease3", "Disease4", "Disease5", "Done",
    )
    raw = (
        120140, "COL2A1", "Czech dysplasia (609162)(AD)", "AD",
        "Czech dysplasia (609162)(AD)", None, None, None, None, None,
    )

    row = omim_store._row_to_dict(headers, raw)

    assert row["Disease1"] == "Czech dysplasia (609162)(AD)"
    assert all(row[field] == "" for field in omim_store.DISEASE_FIELDS[5:])


def test_new_workbook_row_keeps_disease16():
    headers = (
        "OMIM_id", "gene_symbol", "OMIM_disease", "Inheritance",
        *omim_store.DISEASE_FIELDS, "Done",
    )
    raw = [120140, "COL2A1", "Late disease (654321)(AD)", "AD"]
    raw.extend([None] * 15)
    raw.extend(["Late disease (654321)(AD)\n\nCurator detail.", None])

    row = omim_store._row_to_dict(headers, tuple(raw))

    assert row["Disease16"].endswith("Curator detail.")


def test_synthesis_uses_all_16_slots_only_when_every_curated_slot_is_empty():
    row = omim_store._empty_row()
    row["OMIM_disease"] = "\n".join(
        f"Disease {idx} ({600000 + idx})(AD)" for idx in range(1, 17)
    )

    omim_store._synthesize_diseases(row)

    assert row["Disease1"] == "Disease 1 (600001)(AD)"
    assert row["Disease16"] == "Disease 16 (600016)(AD)"

    curated = omim_store._empty_row()
    curated["Disease1"] = "Curated disease (600001)(AD)\n\nRich detail."
    curated["OMIM_disease"] = row["OMIM_disease"]

    omim_store._synthesize_diseases(curated)

    assert curated["Disease1"].endswith("Rich detail.")
    assert curated["Disease2"] == ""
    assert curated["Disease16"] == ""


def test_disease16_becomes_an_omim_association_and_receives_evidence(monkeypatch):
    supplemental = gene_disease_store._row_to_item({
        "gene": "COL2A1",
        "disease_name": "Late disease",
        "source": "GenCC",
        "classification": "Definitive",
        "phenotype_mim": "654321",
        "source_disease_id": "GenCC:654321",
        "inheritance": "AD",
    })
    monkeypatch.setattr(
        gene_disease_store,
        "lookup_cached",
        lambda gene: [supplemental],
    )

    associations = gene_disease_store.merged_associations(
        "COL2A1", _omim_row_with_disease16(), refresh=False
    )

    assert [item["id"] for item in associations] == ["omim-slot:16"]
    assert associations[0]["omim_slot"] == 16
    assert associations[0]["evidence"] == ["OMIM", "GenCC Definitive"]


def test_sample_payload_case_summary_and_docx_accept_disease16(
    tmp_path, monkeypatch
):
    omim_row = _omim_row_with_disease16()
    monkeypatch.setattr(omim_store, "ensure_loaded", lambda: None)
    monkeypatch.setattr(omim_store, "lookup_cached", lambda **kwargs: omim_row)
    monkeypatch.setattr(gene_disease_store, "ensure_loaded", lambda: None)
    monkeypatch.setattr(
        gene_disease_store,
        "merged_associations",
        lambda gene, row, refresh=False: gene_disease_store._omim_associations(row),
    )
    monkeypatch.setattr(litvar2_on_demand, "apply_cached", lambda *args: None)
    variants = {"v1": {"gene_symbol": "COL2A1"}}

    sample_loader._enrich_snv_variants(variants, "S1", tmp_path)

    assert variants["v1"]["Disease16"].endswith("Curator detail.")
    assert variants["v1"]["disease_associations"][0]["omim_slot"] == 16

    edits = {"report_diseases": {"16": True}}
    assert sample_loader._case_selected_diseases(variants["v1"], edits) == [
        "Late disease (654321)(AD)"
    ]
    assert docx_export._picked_disease_for_snv(variants["v1"], edits) == (
        "Late disease (654321)(AD)\n\nCurator detail."
    )
