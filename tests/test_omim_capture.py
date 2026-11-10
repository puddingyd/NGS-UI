"""OMIM capture parser / storage / API (synthetic page text only).

The fixtures imitate the *structure* of a rendered omim.org entry page
(what Ctrl+A / Ctrl+C or the extension's innerText yields) with made-up
disease content, so no OMIM text is committed.
"""
import gzip
import json

import pytest
from fastapi import FastAPI
from fastapi.testclient import TestClient

from app import auth
from app.routers import omim_capture as router_mod
from app.services import capture_tokens, omim_capture, omim_store, users

NAV = "OMIM\nAbout\nStatistics \nDownloads \nDisplay:   Feature IDs\n"
TABLE_HEADER = (
    "Location\tPhenotype\tPhenotype\nMIM number\tInheritance\tPhenotype\n"
    "mapping key\tGene/Locus\tGene/Locus\nMIM number\n"
)
FOOTER = (
    "\n\nContributors: Curator A - updated : 08/25/2025\n"
    "Curator B - updated : 01/02/2024\n"
    "Creation Date: Founder : 6/3/1986\n"
    "Edit History: carol : 09/04/2025\n"
    "alopez : 11/12/2024\n"
    "\n\nNOTE: OMIM is intended for use primarily by physicians.\n"
)


def phenotype_page(*, synopsis=True, description=True, rows=None, toc_synopsis=True):
    rows = rows or ["1p36\tTestosis, infantile-onset\t600001\tAR\t3\tGENE1\t100001"]
    toc = ["Title", "Phenotype-Gene Relationships"]
    if toc_synopsis:
        toc.append("Clinical Synopsis")
    toc += ["Text", "Description" if description else "Clinical Features", "References"]
    syn = (
        "   \nINHERITANCE\n- Autosomal recessive\nHEAD & NECK\nFace\n- Long face\n"
        "MOLECULAR BASIS\n- Caused by mutation in GENE1 (GENE1, 100001.0001)\n▲ Close\n"
    ) if synopsis else ""
    desc = (
        "▼ Description\n\nTestosis is a made-up disorder (summary by Doe et al., 2020).\n\n"
        "Genetic Heterogeneity of Testosis\n\nSee also TST2.\n\n\n"
    ) if description else ""
    return (
        NAV + "#600001\nTable of Contents\n\n" + "\n".join(toc)
        + "\n▼  \tExternal Links\n \nICD+\n# 600001\nTESTOSIS, INFANTILE-ONSET; TSTI\n\n"
        "Alternative titles; symbols\n\nTESTOSIS TYPE I\nTST1\n\n\n"
        "Other entities represented in this entry:\n\nMILD TESTOSIS, INCLUDED\n\n"
        "Phenotype-Gene Relationships\n\n" + TABLE_HEADER + "\n".join(rows) + "\n"
        + syn + "\n▼ TEXT\n\nA number sign (#) is used with this entry.\n\n\n" + desc
        + "▼ Clinical Features\n\nPatients had long faces.\n\nClinical Variability\n\nSome were mild.\n\n\n"
        "▼ REFERENCES\n\nDoe, J. Made up. J. Test. 1: 1-2, 2020.\n" + FOOTER
    )


GENE_PAGE = (
    NAV + "*100001\nTable of Contents\n\nTitle\nGene-Phenotype Relationships\nText\n"
    "▼  \tExternal Links\n \n* 100001\nTEST GENE 1; GENE1\n\n"
    "Alternative titles; symbols\n\nTG1\n\n\n"
    "Gene-Phenotype Relationships\n\n"
    "Location\tPhenotype\tPhenotype\nMIM number\tInheritance\tPhenotype\nmapping key\n"
    "1p36\tTestosis, infantile-onset\t600001\tAR\t3\n"
    "Testosis, adult-onset, new name\t600002\tAD\t3\n"
    "{Susceptibility to testing}\t\t\t3\n"
    "\n▼ TEXT\n\n▼ Description\n\nGENE1 encodes a test protein.\n" + FOOTER
)


def test_parses_phenotype_entry_structure():
    p = omim_capture.parse_entry_text(phenotype_page(), url_mim="600001")
    assert (p["mim"], p["prefix"], p["entry_type"]) == ("600001", "#", "phenotype")
    assert (p["title"], p["symbol"]) == ("TESTOSIS, INFANTILE-ONSET", "TSTI")
    assert p["alternative_titles"] == ["TESTOSIS TYPE I", "TST1"]
    assert p["included_entities"] == ["MILD TESTOSIS, INCLUDED"]
    assert p["relationships"] == [{
        "location": "1p36", "phenotype": "Testosis, infantile-onset",
        "phenotype_mim": "600001", "inheritance": "AR", "mapping_key": "3",
        "gene": "GENE1", "gene_mim": "100001",
    }]
    assert p["description_source"] == "Description"
    assert "Genetic Heterogeneity of Testosis" in p["description"]
    assert p["synopsis"].splitlines()[0] == "INHERITANCE"
    assert p["synopsis"].splitlines()[-1].startswith("- Caused by mutation")
    assert p["synopsis_status"] == "present"
    assert p["omim_edit_date"] == "2025-09-04"
    assert p["omim_contributor_date"] == "2025-08-25"
    assert p["omim_creation_date"] == "1986-06-03"
    assert p["warnings"] == []
    cell = omim_capture.disease_cell(p, p["relationships"][0])
    assert cell.startswith("Testosis, infantile-onset (600001)(AR)\n\nTestosis is a made-up")
    assert cell.endswith("(GENE1, 100001.0001)")


def test_falls_back_to_clinical_features_and_flags_collapsed_synopsis():
    p = omim_capture.parse_entry_text(
        phenotype_page(synopsis=False, description=False), url_mim="600001"
    )
    assert p["description_source"] == "Clinical Features"
    assert p["description"].startswith("Patients had long faces.")
    assert "Clinical Variability" in p["description"]
    assert p["synopsis_status"] == "not_loaded"
    assert any("Clinical Features" in w for w in p["warnings"])
    assert any("Clinical Synopsis" in w for w in p["warnings"])


def test_synopsis_absent_when_entry_has_none():
    p = omim_capture.parse_entry_text(
        phenotype_page(synopsis=False, toc_synopsis=False), url_mim="600001"
    )
    assert p["synopsis_status"] == "absent"


def test_multi_gene_phenotype_keeps_every_row():
    p = omim_capture.parse_entry_text(phenotype_page(rows=[
        "1p36\tTestosis, infantile-onset\t600001\tAR\t3\tGENE1\t100001",
        "2q11\tTestosis, infantile-onset\t600001\tAR\t3\tGENE2\t100002",
    ]), url_mim="600001")
    assert [r["gene"] for r in p["relationships"]] == ["GENE1", "GENE2"]
    assert any("2 個基因" in w for w in p["warnings"])


def test_extension_synopsis_page_fallback():
    extracted = {"synopsis_text": (
        "OMIM\n#600001\nTESTOSIS\nClinical Synopsis\nINHERITANCE\n- Autosomal recessive\n"
        "MISCELLANEOUS\n- Onset in infancy\nContributors: Someone\n"
    )}
    p = omim_capture.parse_entry_text(
        phenotype_page(synopsis=False), url_mim="600001", extracted=extracted
    )
    assert p["synopsis_source"] == "clinical_synopsis_page"
    assert p["synopsis"] == "INHERITANCE\n- Autosomal recessive\nMISCELLANEOUS\n- Onset in infancy"


def test_mim_mismatch_and_dom_table_fallback():
    p = omim_capture.parse_entry_text(phenotype_page(), url_mim="600009")
    assert p["mim"] == "600009"
    assert any("找不到 entry" in w for w in p["warnings"])
    text = phenotype_page().replace(
        "1p36\tTestosis, infantile-onset\t600001\tAR\t3\tGENE1\t100001\n", ""
    )
    p = omim_capture.parse_entry_text(text, url_mim="600001", extracted={
        "relationship_rows": [["1p36", "Testosis, infantile-onset", "600001", "AR", "3", "GENE1", "100001"]],
    })
    assert p["relationships"][0]["gene"] == "GENE1"
    assert any("表格結構" in w for w in p["warnings"])


def test_gene_page_rows_without_gene_columns():
    p = omim_capture.parse_entry_text(GENE_PAGE, url_mim="100001")
    assert (p["prefix"], p["entry_type"], p["symbol"]) == ("*", "gene", "GENE1")
    assert [(r["phenotype"], r["phenotype_mim"], r["inheritance"], r["gene"]) for r in p["relationships"]] == [
        ("Testosis, infantile-onset", "600001", "AR", "GENE1"),
        ("Testosis, adult-onset, new name", "600002", "AD", "GENE1"),
        ("{Susceptibility to testing}", "", "", "GENE1"),
    ]
    assert p["synopsis_status"] == "not_applicable"


def test_content_hash_ignores_volatile_fields():
    a = omim_capture.parse_entry_text(phenotype_page(), url_mim="600001")
    b = omim_capture.parse_entry_text(
        phenotype_page().replace("09/04/2025", "10/04/2025"), url_mim="600001"
    )
    assert a["content_hash"] == b["content_hash"]


@pytest.fixture
def workbook(monkeypatch):
    rows = {
        "GENE1": {
            **omim_store._empty_row(),
            "OMIM_id": "100001", "gene_symbol": "GENE1", "Done": "2025-01-01",
            "OMIM_disease": (
                "Testosis, infantile-onset (600001)(AR)\n"
                "Testosis, adult-onset (600002)(AD)\n"
                "Old removed disease (600003)(AR)\n"
                "Label without mim (AD)"
            ),
            "Disease1": "Testosis, infantile-onset (600001)(AR)\n\nOld detail.",
            "Disease2": "Testosis, adult-onset (600002)(AD)",
        },
    }
    monkeypatch.setattr(omim_store, "_state", {
        "mtime": 1, "by_gene": rows, "by_omim_id": {100001: rows["GENE1"]},
    })
    monkeypatch.setattr(omim_store, "ensure_loaded", lambda: None)
    monkeypatch.setattr(omim_store, "cache_signature", lambda: ("test", id(rows)))
    return rows


def test_compare_phenotype_capture_with_workbook(workbook):
    newer = omim_capture.compare_with_workbook("600001", edit_date="2025-09-04")
    assert newer["status"] == "outdated"
    assert newer["phenotype_links"][0]["slots"] == [1]
    older = omim_capture.compare_with_workbook("600001", edit_date="2024-01-01")
    assert older["status"] == "up_to_date"
    assert omim_capture.compare_with_workbook("600002")["status"] == "missing_detail"
    assert omim_capture.compare_with_workbook("699999")["status"] == "not_in_workbook"


def test_compare_gene_capture_lists_skeleton_differences(workbook):
    p = omim_capture.parse_entry_text(GENE_PAGE, url_mim="100001")
    gene = omim_capture.compare_with_workbook(
        "100001", entry_type="gene", relationships=p["relationships"]
    )["gene"]
    assert gene["gene"] == "GENE1"
    assert gene["added_in_omim"] == []
    assert gene["missing_from_omim"] == ["Old removed disease (600003)(AR)"]
    assert gene["workbook_lines_without_mim"] == ["Label without mim (AD)"]
    assert [r["phenotype_mim"] for r in gene["renamed"]] == ["600002"]


@pytest.fixture
def client(tmp_path, monkeypatch, workbook):
    monkeypatch.setattr(users, "_DB_PATH", tmp_path / "users.db")
    monkeypatch.setattr(users, "DATA_ROOT", tmp_path)
    monkeypatch.setattr(omim_capture, "OMIM_CAPTURE_DIR", tmp_path / "captures")
    user = users.create_user("reviewer", "secret-pass")
    app = FastAPI()
    app.include_router(router_mod.router)
    app.dependency_overrides[auth.current_user] = lambda: user
    return TestClient(app), user, tmp_path


def test_capture_api_requires_valid_token(client):
    c, user, _ = client
    body = {"mim": "600001", "url": "https://omim.org/entry/600001", "text": phenotype_page()}
    assert c.post("/api/omim/captures", json=body, headers={
        "Authorization": "Bearer ngsc_wrong"}).status_code == 401
    token = capture_tokens.create(user, "laptop")["token"]
    ping = c.get("/api/omim/ping", headers={"Authorization": f"Bearer {token}"})
    assert ping.json() == {"ok": True, "username": "reviewer"}
    assert capture_tokens.list_for_user(user)[0]["last_used_at"]


def test_capture_api_stores_parses_and_dedupes(client):
    c, user, tmp_path = client
    headers = {"Authorization": f"Bearer {capture_tokens.create(user)['token']}"}
    body = {"mim": "600001", "url": "https://omim.org/entry/600001", "text": phenotype_page(),
            "extension_version": "0.1.0"}
    first = c.post("/api/omim/captures", json=body, headers=headers).json()
    assert first["capture"]["duplicate"] is False
    summary = first["summary"]
    assert summary["synopsis_status"] == "present"
    assert summary["comparison"]["status"] == "outdated"
    assert summary["disease_cells"][0]["gene"] == "GENE1"
    files = list((tmp_path / "captures" / "600001").glob("*.json.gz"))
    assert len(files) == 1
    with gzip.open(files[0], "rt", encoding="utf-8") as fh:
        stored = json.load(fh)
    assert stored["meta"]["user"] == "reviewer"
    assert stored["payload"]["text"] == phenotype_page()

    again = c.post("/api/omim/captures", json=body, headers=headers).json()
    assert again["capture"]["duplicate"] is True
    assert len(list((tmp_path / "captures" / "600001").glob("*.json.gz"))) == 1

    listing = c.get("/api/omim/captures").json()["captures"]
    assert len(listing) == 2 and listing[0]["mim"] == "600001"
    mim, name = first["capture"]["id"].split("/")
    detail = c.get(f"/api/omim/captures/{mim}/{name}").json()
    assert detail["summary"]["title"] == "TESTOSIS, INFANTILE-ONSET"
    status = c.get("/api/omim/status?mim=600001", headers=headers).json()
    assert status["latest_capture"]["omim_edit_date"] == "2025-09-04"
    assert status["comparison"]["status"] == "outdated"


@pytest.mark.parametrize("body,code", [
    ({"mim": "60000", "text": "x"}, 400),
    ({"mim": "600001", "text": " "}, 400),
    ({"mim": "600001", "text": "no omim entry here"}, 200),
])
def test_capture_api_validation(client, body, code):
    c, user, _ = client
    headers = {"Authorization": f"Bearer {capture_tokens.create(user)['token']}"}
    response = c.post("/api/omim/captures", json=body, headers=headers)
    assert response.status_code == code
    if code == 200:
        assert any("找不到 entry" in w for w in response.json()["summary"]["warnings"])


def test_token_management_and_extension_zip(client):
    c, user, _ = client
    created = c.post("/api/omim/capture-tokens", json={"label": "office PC"}).json()
    assert created["token"].startswith("ngsc_")
    tokens = c.get("/api/omim/capture-tokens").json()["tokens"]
    assert tokens[0]["label"] == "office PC" and "token" not in tokens[0]
    assert c.delete(f"/api/omim/capture-tokens/{created['id']}").json() == {"ok": True}
    assert capture_tokens.verify(created["token"]) is None
    zipped = c.get("/api/omim/extension.zip")
    assert zipped.status_code == 200
    assert zipped.content[:2] == b"PK"


def test_status_uses_page_edit_date(client):
    c, user, _ = client
    headers = {"Authorization": f"Bearer {capture_tokens.create(user)['token']}"}
    assert c.get("/api/omim/status?mim=600001&edit_date=2024-01-01", headers=headers).json()[
        "comparison"]["status"] == "up_to_date"
    assert c.get("/api/omim/status?mim=600001&edit_date=2026-01-01", headers=headers).json()[
        "comparison"]["status"] == "outdated"
    assert c.get("/api/omim/status?mim=600001&edit_date=bad", headers=headers).status_code == 400


def test_synopsis_drops_expander_label_line():
    text = phenotype_page().replace(
        "   \nINHERITANCE\n", "► Clinical Synopsis\n   \nINHERITANCE\n"
    )
    p = omim_capture.parse_entry_text(text, url_mim="600001")
    assert p["synopsis"].splitlines()[0] == "INHERITANCE"
    assert p["synopsis_status"] == "present"


def test_synopsis_drops_page_widgets_and_blank_lines():
    noisy = (
        "   \nToggle Dropdown\n\nPhenotypic Series\nToggle Dropdown\n\nPheneGene Graphics\n\n"
        "INHERITANCE\n- Autosomal recessive \n\nHEAD & NECK\nFace\n- Long face\n\n"
        "MOLECULAR BASIS\n- Caused by mutation in GENE1 (GENE1, 100001.0001)\n▲ Close\n"
    )
    clean = (
        "   \nINHERITANCE\n- Autosomal recessive\nHEAD & NECK\nFace\n- Long face\n"
        "MOLECULAR BASIS\n- Caused by mutation in GENE1 (GENE1, 100001.0001)\n▲ Close\n"
    )
    base = phenotype_page()
    assert clean in base
    p = omim_capture.parse_entry_text(base.replace(clean, noisy), url_mim="600001")
    assert p["synopsis"] == (
        "INHERITANCE\n- Autosomal recessive\nHEAD & NECK\nFace\n- Long face\n"
        "MOLECULAR BASIS\n- Caused by mutation in GENE1 (GENE1, 100001.0001)"
    )
    cell = omim_capture.disease_cell(p, p["relationships"][0])
    assert "\n\n\n" not in cell and "Toggle" not in cell
    assert "\n\n" not in p["synopsis"]
    assert cell.endswith("\n\n" + p["synopsis"])


def test_synopsis_page_fallback_drops_widgets():
    extracted = {"synopsis_text": (
        "Clinical Synopsis\nToggle Dropdown\nPhenotypic Series\nINHERITANCE\n\n- Autosomal recessive\n"
        "\nMISCELLANEOUS\n- Onset in infancy\nNOTE: OMIM is intended\n"
    )}
    p = omim_capture.parse_entry_text(
        phenotype_page(synopsis=False), url_mim="600001", extracted=extracted
    )
    assert p["synopsis"] == "INHERITANCE\n- Autosomal recessive\nMISCELLANEOUS\n- Onset in infancy"


def test_old_captures_are_reparsed_with_current_rules(client):
    c, user, tmp_path = client
    headers = {"Authorization": f"Bearer {capture_tokens.create(user)['token']}"}
    noisy = phenotype_page().replace("   \nINHERITANCE\n", "Toggle Dropdown\n\nPheneGene Graphics\n\nINHERITANCE\n")
    created = c.post("/api/omim/captures", json={"mim": "600001", "text": noisy}, headers=headers).json()
    path = next((tmp_path / "captures" / "600001").glob("*.json.gz"))
    with gzip.open(path, "rt", encoding="utf-8") as fh:
        data = json.load(fh)
    data["parsed"]["parser_version"] = 1
    data["parsed"]["synopsis"] = "Toggle Dropdown\nINHERITANCE"
    with gzip.open(path, "wt", encoding="utf-8") as fh:
        json.dump(data, fh)
    mim, name = created["capture"]["id"].split("/")
    cell = c.get(f"/api/omim/captures/{mim}/{name}").json()["summary"]["disease_cells"][0]["cell"]
    assert "Toggle" not in cell and "\nINHERITANCE\n- Autosomal recessive\n" in cell
