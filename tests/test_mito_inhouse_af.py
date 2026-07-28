from pathlib import Path

from app.adapters import mito_tsv
from app.services import inhouse_af_mito


REPO_ROOT = Path(__file__).resolve().parents[1]
APP_JS = (REPO_ROOT / "frontend" / "app.js").read_text(encoding="utf-8")


def test_chrM_vcf_parser_preserves_carrier_frequency_contract():
    table = inhouse_af_mito._table_from_lines([
        (
            "chrM\t73\t.\tA\tG,C\t.\t.\t"
            "INHOUSE_AC=12,3;INHOUSE_AN=677;"
            "INHOUSE_AF=0.0177253,0.00443131;"
            "INHOUSE_NHOM=8,1;INHOUSE_HET_MT=4,2"
        ),
    ])

    assert table[(73, "A", "G")] == {
        "inhouse_ac": 12,
        "inhouse_an": 677,
        "inhouse_af": 0.0177253,
        "inhouse_nhom": 8,
        "inhouse_het_mt": 4,
    }
    assert table[(73, "A", "C")]["inhouse_ac"] == 3
    assert table[(73, "A", "C")]["inhouse_an"] == 677


def test_chrM_lookup_uses_same_minimal_representation_as_db(monkeypatch):
    monkeypatch.setattr(
        inhouse_af_mito,
        "_TABLE",
        {
            (101, "T", "G"): {
                "inhouse_ac": 2,
                "inhouse_an": 600,
                "inhouse_af": 2 / 600,
                "inhouse_nhom": 0,
                "inhouse_het_mt": 2,
            },
        },
    )

    hit = inhouse_af_mito.lookup(100, "AT", "AG")

    assert hit["inhouse_ac"] == 2
    assert hit["inhouse_het_mt"] == 2


def test_mito_adapter_surfaces_inhouse_af_without_changing_tier(tmp_path, monkeypatch):
    path = tmp_path / "sample.mito.tsv"
    path.write_text(
        "CHROM\tPOS\tREF\tALT\tGENE\tFILTER\tGNOMAD_MITO_AF\n"
        "chrM\t73\tA\tG\tMT-TF\tPASS\t0.02\n",
        encoding="utf-8",
    )
    monkeypatch.setattr(mito_tsv.clinvar_mito, "lookup", lambda *_args: {})
    monkeypatch.setattr(
        mito_tsv.inhouse_af_mito,
        "lookup",
        lambda *_args: {
            "inhouse_ac": 12,
            "inhouse_an": 677,
            "inhouse_af": 12 / 677,
            "inhouse_nhom": 8,
            "inhouse_het_mt": 4,
        },
    )

    variants, categories = mito_tsv.load_mito_tsv(path)
    variant = variants["chrM-73-A-G"]

    assert variant["inhouse_af"] == 12 / 677
    assert variant["inhouse_ac"] == 12
    assert variant["inhouse_an"] == 677
    assert variant["inhouse_nhom"] == 8
    assert variant["inhouse_het_mt"] == 4
    # gnomAD AF=0.02 keeps this otherwise unreported variant in MITO-3;
    # the in-house AF is contextual and must not change classification.
    assert categories["MITO-3"] == ["chrM-73-A-G"]


def test_mito_card_labels_inhouse_af_as_chrM_carrier_frequency():
    assert "function _formatMitoInhouseAf(v)" in APP_JS
    assert "<strong>AF_nckuh:</strong>" in APP_JS
    assert "mtDNA carrier frequency" in APP_JS
    assert "AC 是帶有此 ALT 的樣本數" in APP_JS
    assert 'carrierKinds.push(`hom ${v.inhouse_nhom}`)' in APP_JS
    assert 'carrierKinds.push(`het ${v.inhouse_het_mt}`)' in APP_JS
