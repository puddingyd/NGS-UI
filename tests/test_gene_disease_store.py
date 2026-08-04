from app.services import gene_disease_store


def _supplemental(
    disease_name: str,
    phenotype_mim: str,
    classification: str,
    *,
    mondo_id: str,
) -> dict:
    return gene_disease_store._row_to_item(
        {
            "gene": "OTX2",
            "disease_name": disease_name,
            "source": "GenCC",
            "classification": classification,
            "mondo_id": mondo_id,
            "phenotype_mim": phenotype_mim,
            "source_disease_id": f"GenCC:{phenotype_mim}:{classification}",
            "inheritance": "AD",
        }
    )


def _otx2_supplemental_rows() -> list[dict]:
    return [
        _supplemental(
            "syndromic microphthalmia type 5",
            "610125",
            "Definitive",
            mondo_id="MONDO:0012413",
        ),
        _supplemental(
            "syndromic microphthalmia type 5",
            "610125",
            "Strong",
            mondo_id="MONDO:0012413",
        ),
        _supplemental(
            "syndromic microphthalmia type 5",
            "610125",
            "Moderate",
            mondo_id="MONDO:0012413",
        ),
        _supplemental(
            "pituitary hormone deficiency, combined, 6",
            "613986",
            "Strong",
            mondo_id="MONDO:0013518",
        ),
        _supplemental(
            "pituitary hormone deficiency, combined, 6",
            "613986",
            "Moderate",
            mondo_id="MONDO:0013518",
        ),
        _supplemental(
            "agnathia-otocephaly complex",
            "202650",
            "Strong",
            mondo_id="MONDO:0008740",
        ),
    ]


def test_duplicate_phenotype_mim_preserves_every_omim_slot(monkeypatch):
    omim_row = {
        "OMIM_id": "600037",
        "Inheritance": "AD",
        "Disease1": "Microphthalmia, syndromic 5 (610125)(AD)\n\nFirst OMIM detail.",
        "Disease2": (
            "Retinal dystrophy, early-onset, with or without pituitary "
            "dysfunction (610125)(AD)\n\nSecond OMIM detail."
        ),
        "Disease3": (
            "Pituitary hormone deficiency, combined, 6 (613986)(AD)"
            "\n\nThird OMIM detail."
        ),
        "Disease4": "",
        "Disease5": "",
    }
    monkeypatch.setattr(
        gene_disease_store,
        "lookup_cached",
        lambda gene: _otx2_supplemental_rows(),
    )

    associations = gene_disease_store.merged_associations(
        "OTX2", omim_row, refresh=False
    )

    assert [item["id"] for item in associations] == [
        "omim-slot:1",
        "omim-slot:2",
        "omim-slot:3",
        "mim:202650",
    ]
    assert [item["display_name"] for item in associations[:3]] == [
        "Microphthalmia, syndromic 5",
        "Retinal dystrophy, early-onset, with or without pituitary dysfunction",
        "Pituitary hormone deficiency, combined, 6",
    ]
    assert [item["omim_slot"] for item in associations[:3]] == [1, 2, 3]
    assert associations[0]["detail"] == omim_row["Disease1"]
    assert associations[1]["detail"] == omim_row["Disease2"]
    assert associations[2]["detail"] == omim_row["Disease3"]

    # MIM 610125 applies to two curator-owned OMIM slots.  Evidence enriches
    # both while their distinct OMIM labels/details remain authoritative.
    assert associations[0]["evidence"] == [
        "OMIM",
        "GenCC Definitive",
        "GenCC Strong",
        "GenCC Moderate",
    ]
    assert associations[1]["evidence"] == associations[0]["evidence"]

    # MIM 613986 maps to one OMIM slot, so evidence enriches that row while
    # its OMIM label and detail remain authoritative.
    assert associations[2]["evidence"] == [
        "OMIM",
        "GenCC Strong",
        "GenCC Moderate",
    ]
    assert associations[3]["display_name"] == "agnathia-otocephaly complex"


def test_supplemental_only_associations_still_group_by_relationship(monkeypatch):
    monkeypatch.setattr(
        gene_disease_store,
        "lookup_cached",
        lambda gene: _otx2_supplemental_rows(),
    )

    associations = gene_disease_store.merged_associations(
        "OTX2", None, refresh=False
    )

    assert [item["id"] for item in associations] == [
        "mim:610125",
        "mim:613986",
        "mim:202650",
    ]
    assert associations[0]["evidence"] == [
        "GenCC Definitive",
        "GenCC Strong",
        "GenCC Moderate",
    ]
    assert associations[1]["evidence"] == [
        "GenCC Strong",
        "GenCC Moderate",
    ]
