import json
import subprocess
from pathlib import Path

import pytest

from app.services import disease_picks, docx_export, omim_store, sample_loader

ROOT = Path(__file__).resolve().parents[1]

OTX2 = {
    "Disease1": "Microphthalmia, syndromic 5 (610125)(AD)\n\nDetail A.",
    "Disease2": (
        "Retinal dystrophy, early-onset, with or without pituitary "
        "dysfunction (610125)(AD)\n\nDetail B."
    ),
    "Disease3": "Pituitary hormone deficiency, combined, 6 (613986)(AD)",
}


def _picks(record, *slots):
    return {disease_picks.PICKS_FIELD: disease_picks.picks_from_slots(record, slots)}


def test_key_strips_metadata_and_normalizes_name():
    assert disease_picks.pick_key(
        "?Mitochondrial DNA depletion syndrome 16 (hepatic type) (618528)(AR)\nx"
    ) == "618528|mitochondrial dna depletion syndrome 16"
    assert disease_picks.pick_key("Gaucher disease (AR)") == "|gaucher disease"


def test_legacy_slots_still_work_without_new_field():
    texts = disease_picks.selected_texts(OTX2, {"report_diseases": {"2": True}})
    assert texts == [OTX2["Disease2"]]


def test_same_mim_different_names_stay_distinct():
    edits = _picks(OTX2, 2)
    assert disease_picks.selected_texts(OTX2, edits) == [OTX2["Disease2"]]


def test_reordered_workbook_keeps_the_ticked_disease():
    edits = _picks(OTX2, 3)
    reordered = {
        "Disease1": OTX2["Disease3"],
        "Disease2": OTX2["Disease1"],
        "Disease3": OTX2["Disease2"],
    }
    # Legacy slot 3 would now point at the retinal dystrophy label.
    assert disease_picks.selected_texts(reordered, {"report_diseases": {"3": True}}) == [
        OTX2["Disease2"]
    ]
    assert disease_picks.selected_texts(reordered, edits) == [OTX2["Disease3"]]


def test_rename_with_same_mim_follows_current_text():
    before = {"Disease1": "Mental retardation, autosomal dominant 39 (616521)(AD)"}
    after = {
        "Disease1": "Intellectual developmental disorder, autosomal dominant 39 (616521)(AD)",
    }
    edits = _picks(before, 1)
    assert disease_picks.selected_texts(after, edits) == [after["Disease1"]]


def test_mim_added_later_matches_by_name():
    before = {"Disease1": "Deafness 5 (XLR)"}
    after = {"Disease1": "Other (300001)(AD)", "Disease2": "Deafness 5 (300614)(XLR)"}
    assert disease_picks.selected_texts(after, _picks(before, 1)) == [after["Disease2"]]


def test_disease_removed_from_workbook_keeps_snapshot_label():
    edits = _picks(OTX2, 3)
    shrunk = {"Disease1": OTX2["Disease1"]}
    assert disease_picks.selected_texts(shrunk, edits) == [
        "Pituitary hormone deficiency, combined, 6 (613986)(AD)"
    ]


def test_empty_picks_fall_back_to_first_disease_only_when_asked():
    edits = {disease_picks.PICKS_FIELD: {}}
    assert disease_picks.selected_texts(OTX2, edits) == [OTX2["Disease1"]]
    assert disease_picks.selected_texts(OTX2, edits, fallback_first=False) == []


def test_new_picks_override_stale_legacy_slots():
    edits = {**_picks(OTX2, 1), "report_diseases": {"3": True}}
    assert disease_picks.selected_texts(OTX2, edits) == [OTX2["Disease1"]]


def test_docx_and_case_summary_use_resolved_picks(monkeypatch):
    edits = _picks(OTX2, 3)
    reordered = {
        "gene_symbol": "OTX2",
        "Disease1": OTX2["Disease3"],
        "Disease2": OTX2["Disease1"],
    }
    assert docx_export._picked_diseases_for_snv(reordered, edits) == [OTX2["Disease3"]]
    assert "613986" in docx_export._omim_block_for_snv(reordered, edits)

    monkeypatch.setattr(omim_store, "lookup_cached", lambda **_: dict(reordered))
    assert sample_loader._case_selected_diseases({"gene_symbol": "OTX2"}, edits) == [
        "Pituitary hormone deficiency, combined, 6 (613986)(AD)"
    ]
    assert sample_loader._case_selected_diseases(
        {"gene_symbol": "OTX2"}, {disease_picks.PICKS_FIELD: {}}
    ) == []


@pytest.mark.parametrize("text", [
    "Microphthalmia, syndromic 5 (610125)(AD)\n\nDetail",
    "?Nephrotic syndrome, type 19 (618178)(AR)",
    "{Preeclampsia, susceptibility to}",
    "Acute myeloid leukemia, somatic (601626)()",
    "Hypophosphatasia, adult (146300)(AR, AD)",
    "Ectodermal dysplasia 16 (odontoonychodermal dysplasia) (257980)(AR)",
    "Alpha-thalassemia/impaired intellectual development syndrome (301040)(XLD)",
    "Nephrotic syndrome, type 20, (301028)(XL)",
    "Lysosomal & storage (123456)(AR)",
    "45832",
])
def test_frontend_key_matches_backend(text):
    source = (ROOT / "frontend" / "app.js").read_text(encoding="utf-8")
    start = source.index("// disease-pick-key:start")
    end = source.index("// disease-pick-key:end")
    script = source[start:end] + f"\nprocess.stdout.write(diseasePickKey({json.dumps(text)}));"
    result = subprocess.run(
        ["node", "-e", script], capture_output=True, text=True, check=True
    )
    assert result.stdout == disease_picks.pick_key(text)


def test_migration_script_freezes_legacy_slots(tmp_path, monkeypatch):
    import importlib.util

    spec = importlib.util.spec_from_file_location(
        "migrate_report_disease_picks", ROOT / "scripts" / "migrate_report_disease_picks.py"
    )
    mig = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mig)

    meta_path = tmp_path / "sample_metadata.json"
    meta_path.write_text(json.dumps({
        "edits": {
            "v1": {"report_diseases": {"3": True}, "comment": "keep"},
            "v2": {"report_diseases": {"1": True}, disease_picks.PICKS_FIELD: {}},
            "v3": {"report_diseases": {"1": True}},
        },
    }))
    monkeypatch.setattr(mig.report_store, "_meta_path", lambda _sid: meta_path)
    monkeypatch.setattr(
        mig.sample_loader, "load_sample",
        lambda *_a, **_k: {"variants": {"v1": dict(OTX2), "v2": dict(OTX2)}},
    )

    conv, skip, _ = mig.migrate_sample("S1", apply=False)
    assert (conv, skip) == (1, 1)
    assert disease_picks.PICKS_FIELD not in json.loads(meta_path.read_text())["edits"]["v1"]

    mig.migrate_sample("S1", apply=True)
    edits = json.loads(meta_path.read_text())["edits"]
    assert edits["v1"]["comment"] == "keep"
    assert edits["v1"]["report_diseases"] == {"3": True}
    assert disease_picks.selected_texts(OTX2, edits["v1"]) == [OTX2["Disease3"]]
    assert edits["v2"][disease_picks.PICKS_FIELD] == {}
    assert disease_picks.PICKS_FIELD not in edits["v3"]
    assert mig.migrate_sample("S1", apply=True)[:2] == (0, 1)
