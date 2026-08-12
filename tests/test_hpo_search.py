from pathlib import Path

import pytest

from app.services import hpo_ontology


@pytest.fixture
def isolated_hpo_index(monkeypatch):
    monkeypatch.setattr(hpo_ontology, "_TERMS", {})
    monkeypatch.setattr(hpo_ontology, "_NAME_INDEX", [])
    monkeypatch.setattr(hpo_ontology, "_SYN_INDEX", [])
    monkeypatch.setattr(hpo_ontology, "_SEARCH_TEXTS", {})
    monkeypatch.setattr(hpo_ontology, "_GRAM_INDEX", {})


def _write_test_obo(path: Path) -> None:
    path.write_text(
        """format-version: 1.2

[Term]
id: HP:0000001
name: Abnormal nervous system physiology

[Term]
id: HP:0001250
name: Seizure
synonym: \"Epileptic convulsion\" EXACT []
is_a: HP:0000001 ! Abnormal nervous system physiology

[Term]
id: HP:0001324
name: Muscle weakness
is_a: HP:0000001 ! Abnormal nervous system physiology
""",
        encoding="utf-8",
    )


def test_search_handles_typos_and_returns_immediate_parent(tmp_path, isolated_hpo_index):
    obo = tmp_path / "hp.obo"
    _write_test_obo(obo)
    assert hpo_ontology.load(obo) == 3

    result = hpo_ontology.search("seizuer", limit=20)

    assert result[0]["hpo_id"] == "HP:0001250"
    assert result[0]["parents"] == [
        {"hpo_id": "HP:0000001", "name": "Abnormal nervous system physiology"}
    ]


def test_fuzzy_search_also_matches_misspelled_synonym(tmp_path, isolated_hpo_index):
    obo = tmp_path / "hp.obo"
    _write_test_obo(obo)
    hpo_ontology.load(obo)

    assert hpo_ontology.search("convulsin", limit=20)[0]["hpo_id"] == "HP:0001250"


def test_main_and_new_case_hpo_pickers_share_parent_renderer():
    app_js = (Path(__file__).resolve().parents[1] / "frontend" / "app.js").read_text(
        encoding="utf-8"
    )

    assert "function _renderHpoSearchResults" in app_js
    assert "_renderHpoSearchResults(dropdown, list);" in app_js
    assert "_renderHpoSearchResults(drop, rows, { newCase: true });" in app_js
    assert "hpo-parent-option" in app_js
