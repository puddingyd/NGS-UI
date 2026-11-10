"""Platform OMIM database: migration from OMIM.xlsx, rendering, adoption."""
import io
import json

import openpyxl
import pytest
from fastapi import FastAPI
from fastapi.testclient import TestClient

from app import auth
from app.routers import omim_capture as router_mod
from app.services import capture_tokens, disease_picks, omim_capture, omim_db, omim_store, users

HEADERS = ["OMIM_id", "gene_symbol", "OMIM_disease", "Inheritance"] + [
    f"Disease{i}" for i in range(1, 17)
] + ["Done"]


def _row(mim, symbol, lines, inh="", slots=None, done=None):
    values = [mim, symbol, "\n".join(lines) if lines else None, inh]
    cells = [None] * 16
    for slot, text in (slots or {}).items():
        cells[slot - 1] = text
    return values + cells + [done]


ALIGNED = "Disease A (600001)(AD)\n\nDetail A.\n\nINHERITANCE\n- Autosomal dominant"
WORKBOOK_ROWS = [
    # aligned: two lines, one rich cell, one hidden line (only slot 1 filled)
    _row(100001, "GENE1", ["Disease A (600001)(AD)", "Disease B (600002)(AR)"], "AD, AR",
         {1: ALIGNED}),
    # OTX2-style duplicate MIM with two names, cells in a different order
    _row(100002, "OTX9", ["Micro 5 (610125)(AD)", "Retinal (610125)(AD)"], "AD",
         {1: "Retinal (610125)(AD)\n\nR detail", 2: "Micro 5 (610125)(AD)\n\nM detail"}),
    # renamed: workbook line has no MIM, Disease cell has the new name + MIM
    _row(100003, "CERTX", ["Mental retardation, autosomal dominant 34 (AD)"], "AD",
         {1: "Intellectual developmental disorder, autosomal dominant 34 (616351)(AD)\n\nX"}),
    # different MIMs must not be merged (SCA1 vs SCA10)
    _row(100004, "ATXNX", ["Spinocerebellar ataxia 1 (164400)(AD)"], "AD",
         {1: "Spinocerebellar ataxia 10 (603516)(AD)\n\nwrong gene text"}),
    # alias rows: same gene MIM, old/new symbol, only the alias row has text
    _row(100005, "BLMX", ["Bloom (210900)(AR)"], "AR"),
    _row(100005, "RECQLX", ["Bloom (210900)(AR)"], "AR", {1: "Bloom (210900)(AR)\n\nBloom detail"}),
    # placeholder
    _row(100006, "NILX", ["Nil"]),
    # no diseases at all
    _row(100007, "QUIET", []),
]


def _workbook(path):
    wb = openpyxl.Workbook()
    ws = wb.active
    ws.append(HEADERS)
    for row in WORKBOOK_ROWS:
        ws.append(row)
    wb.save(path)
    return path


@pytest.fixture
def db(tmp_path, monkeypatch):
    xlsx = _workbook(tmp_path / "OMIM.xlsx")
    target = tmp_path / "omim.sqlite"
    built = omim_db.build_from_xlsx(xlsx, target)
    omim_db.install(__import__("pathlib").Path(built["path"]), target)
    monkeypatch.setattr(omim_db, "OMIM_DB", target)
    monkeypatch.setattr(omim_store, "OMIM_DB", target)
    monkeypatch.setattr(omim_store, "OMIM_XLSX", xlsx)
    monkeypatch.setattr(omim_store, "_state", {"mtime": None, "by_omim_id": {}, "by_gene": {}})
    return {"path": target, "xlsx": xlsx, "stats": built["stats"]}


def cells(row):
    return [row[f] for f in omim_store.DISEASE_FIELDS if row.get(f)]


def test_migration_stats_and_reviews(db):
    stats = db["stats"]
    assert stats["genes"] == 7
    assert stats["duplicate_gene_mim"] == 1
    assert stats["placeholder_lines"] == 1
    assert stats["newly_visible"] == 2  # Disease B, SCA1
    assert stats["slot_only"] == 1      # SCA10 cell kept as its own row
    assert stats["review_items"] == {
        "missing_mim": 1, "renamed_match": 1, "slot_not_in_omim_disease": 1,
    }


def test_rendering_keeps_curated_text_and_shows_hidden_lines(db):
    by_id, by_gene = omim_db.load_rows(db["path"])
    assert cells(by_gene["GENE1"]) == [ALIGNED, "Disease B (600002)(AR)"]
    # skeleton order wins; texts follow their own disease, legacy slots map
    otx = by_gene["OTX9"]
    assert cells(otx) == ["Micro 5 (610125)(AD)\n\nM detail", "Retinal (610125)(AD)\n\nR detail"]
    assert otx["Disease_legacy_slots"] == {"1": 2, "2": 1}
    assert cells(by_gene["CERTX"])[0].startswith("Intellectual developmental disorder")
    assert cells(by_gene["ATXNX"]) == [
        "Spinocerebellar ataxia 1 (164400)(AD)",
        "Spinocerebellar ataxia 10 (603516)(AD)\n\nwrong gene text",
    ]
    assert by_gene["BLMX"] is by_gene["RECQLX"] is by_id[100005]
    assert cells(by_gene["BLMX"]) == ["Bloom (210900)(AR)\n\nBloom detail"]
    assert cells(by_gene["NILX"]) == [] and by_gene["NILX"]["OMIM_disease"] == ""


def test_legacy_slot_ticks_follow_the_reordered_disease(db):
    _by_id, by_gene = omim_db.load_rows(db["path"])
    otx = by_gene["OTX9"]
    # old workbook slot 1 was "Retinal"
    assert disease_picks.selected_texts(otx, {"report_diseases": {"1": True}}) == [
        "Retinal (610125)(AD)\n\nR detail"
    ]
    picks = disease_picks.picks_from_slots(otx, [1])
    assert list(picks) == [disease_picks.pick_key("Retinal (610125)(AD)")]


def test_omim_store_switches_to_database(db, monkeypatch):
    assert omim_store.source() == "db"
    rec = omim_store.lookup(gene="RECQLX")
    assert rec["Disease1"].endswith("Bloom detail")
    db["path"].rename(db["path"].with_suffix(".off"))
    assert omim_store.source() == "xlsx"
    assert omim_store.lookup(gene="BLMX")["Disease1"] == "Bloom (210900)(AR)"


def _phenotype_capture(mim="600002", gene_mim="100001", gene="GENE1", name="Disease B renamed",
                       inh="AR", desc="New description.", synopsis="INHERITANCE\n- Autosomal recessive",
                       edit="2026-08-14"):
    return {
        "mim": mim, "entry_type": "phenotype", "prefix": "#", "title": name.upper(),
        "relationships": [{"phenotype": name, "phenotype_mim": mim, "inheritance": inh,
                           "mapping_key": "3", "gene": gene, "gene_mim": gene_mim, "location": ""}],
        "description": desc, "description_source": "Description", "synopsis": synopsis,
        "synopsis_status": "present", "omim_edit_date": edit, "content_hash": f"h-{mim}-{desc}",
    }


def test_adopt_phenotype_capture_updates_label_text_and_status(db):
    result = omim_db.adopt_capture(_phenotype_capture(), user="rev", capture_id="600002/x")
    assert result["text"] == "added"
    assert result["genes"][0]["changes"][0]["action"] == "updated"
    _by_id, by_gene = omim_db.load_rows(db["path"])
    row = by_gene["GENE1"]
    assert cells(row)[1] == (
        "Disease B renamed (600002)(AR)\n\nNew description.\n\nINHERITANCE\n- Autosomal recessive"
    )
    assert row["Done"] == omim_db._today()
    assert "Disease B renamed (600002)(AR)" in row["OMIM_disease"]
    # status follows the adopted capture's OMIM date
    omim_store._state.update(mtime=None)
    assert omim_capture.compare_with_workbook("600002", edit_date="2026-08-14")["status"] == "up_to_date"
    assert omim_capture.compare_with_workbook("600002", edit_date="2026-09-01")["status"] == "outdated"
    # re-adopting replaces
    again = omim_db.adopt_capture(_phenotype_capture(desc="Newer."), user="rev", capture_id="600002/y")
    assert again["text"] == "replaced"
    with omim_db.connect(db["path"]) as con:
        actions = [r["action"] for r in con.execute("SELECT action FROM change_log ORDER BY id")]
    assert actions.count("adopt_text") == 2


def test_adopt_resolves_migration_review_for_that_mim(db):
    omim_db.adopt_capture(
        _phenotype_capture(mim="616351", gene_mim="100003", gene="CERTX",
                           name="Intellectual developmental disorder, autosomal dominant 34", inh="AD"),
        user="rev", capture_id="616351/x",
    )
    counts = omim_db.list_reviews()["open_counts"]
    assert "renamed_match" not in counts
    _by_id, by_gene = omim_db.load_rows(db["path"])
    assert by_gene["CERTX"]["OMIM_disease"] == (
        "Intellectual developmental disorder, autosomal dominant 34 (616351)(AD)"
    )


def test_adopt_new_gene_and_new_phenotype(db):
    omim_db.adopt_capture(
        _phenotype_capture(mim="699999", gene_mim="199999", gene="NEWG", name="Brand new"),
        user="rev", capture_id="699999/x",
    )
    _by_id, by_gene = omim_db.load_rows(db["path"])
    assert cells(by_gene["NEWG"])[0].startswith("Brand new (699999)(AR)\n\nNew description.")


def test_adopt_gene_page_syncs_skeleton_without_deleting(db):
    parsed = {
        "mim": "100004", "entry_type": "gene", "prefix": "*", "symbol": "ATXNX", "title": "ATAXIN X",
        "relationships": [
            {"phenotype": "Spinocerebellar ataxia 1", "phenotype_mim": "164400", "inheritance": "AD",
             "mapping_key": "3", "gene": "ATXNX", "gene_mim": "100004", "location": ""},
            {"phenotype": "New thing", "phenotype_mim": "620000", "inheritance": "AR",
             "mapping_key": "3", "gene": "ATXNX", "gene_mim": "100004", "location": ""},
        ],
        "description": "", "synopsis": "", "content_hash": "g",
    }
    result = omim_db.adopt_capture(parsed, user="rev", capture_id="100004/x")
    actions = [c["action"] for c in result["genes"][0]["changes"]]
    assert actions == ["unchanged", "added"]
    assert result["reviews"] == 1  # SCA10 row is not listed by OMIM
    items = omim_db.list_reviews(kind="not_in_omim")["items"]
    assert items[0]["current_label"].startswith("Spinocerebellar ataxia 10")
    _by_id, by_gene = omim_db.load_rows(db["path"])
    assert len(cells(by_gene["ATXNX"])) == 3  # nothing deleted yet
    omim_db.resolve_review(items[0]["id"], "remove", user="rev")
    _by_id, by_gene = omim_db.load_rows(db["path"])
    assert [c.splitlines()[0] for c in cells(by_gene["ATXNX"])] == [
        "Spinocerebellar ataxia 1 (164400)(AD)", "New thing (620000)(AR)",
    ]
    assert by_gene["ATXNX"]["Inheritance"] == "AD, AR"


def test_review_actions_validate(db):
    item = omim_db.list_reviews(kind="missing_mim")["items"][0]
    with pytest.raises(ValueError):
        omim_db.resolve_review(item["id"], "delete-everything", user="rev")
    assert omim_db.resolve_review(item["id"], "keep", user="rev")["resolution"] == "keep"
    assert omim_db.resolve_review(item["id"], "keep", user="rev")["status"] == "resolved"
    with pytest.raises(KeyError):
        omim_db.resolve_review(99999, "keep", user="rev")


def test_moved_entries_cannot_be_adopted(db):
    with pytest.raises(ValueError):
        omim_db.adopt_capture({"mim": "600009", "entry_type": "moved_removed"}, user="r", capture_id="x")


def test_export_roundtrip(db, tmp_path):
    omim_db.adopt_capture(_phenotype_capture(), user="rev", capture_id="600002/x")
    data = omim_db.export_xlsx(db["path"])
    wb = openpyxl.load_workbook(io.BytesIO(data))
    ws = wb["工作表1"]
    rows = list(ws.iter_rows(values_only=True))
    assert list(rows[0]) == HEADERS
    gene1 = next(r for r in rows if r[1] == "GENE1")
    assert gene1[0] == 100001 and gene1[4] == ALIGNED and gene1[5].startswith("Disease B renamed")
    assert "Overflow" in wb.sheetnames
    # Rebuilding from the export renders the same disease cells.
    exported = tmp_path / "exported.xlsx"
    exported.write_bytes(data)
    rebuilt = omim_db.build_from_xlsx(exported, tmp_path / "rebuilt.sqlite")
    _a, before = omim_db.load_rows(db["path"])
    _b, after = omim_db.load_rows(__import__("pathlib").Path(rebuilt["path"]))
    assert cells(after["GENE1"]) == cells(before["GENE1"])


@pytest.fixture
def client(db, tmp_path, monkeypatch):
    monkeypatch.setattr(users, "_DB_PATH", tmp_path / "users.db")
    monkeypatch.setattr(users, "DATA_ROOT", tmp_path)
    monkeypatch.setattr(omim_capture, "OMIM_CAPTURE_DIR", tmp_path / "captures")
    user = users.create_user("reviewer", "secret-pass")
    app = FastAPI()
    app.include_router(router_mod.router)
    app.dependency_overrides[auth.current_user] = lambda: user
    return TestClient(app), user


PAGE = (
    "OMIM\n#600002\nTable of Contents\nTitle\nPhenotype-Gene Relationships\nClinical Synopsis\n"
    "▼ External Links\n# 600002\nDISEASE B RENAMED; DBR\n\nPhenotype-Gene Relationships\n\n"
    "Location\tPhenotype\tPhenotype\nMIM number\tInheritance\tPhenotype\nmapping key\tGene/Locus\tGene/Locus\nMIM number\n"
    "1p36\tDisease B renamed\t600002\tAR\t3\tGENE1\t100001\n"
    "Toggle Dropdown\n\nINHERITANCE\n- Autosomal recessive\n▲ Close\n\n▼ TEXT\n\nNote.\n\n"
    "▼ Description\n\nDisease B is made up.\n\n▼ REFERENCES\n\nRef.\n\n"
    "Edit History: carol : 08/14/2026\n\nNOTE: OMIM is intended\n"
)


def test_api_capture_preview_adopt_reviews_export(client):
    c, user = client
    headers = {"Authorization": f"Bearer {capture_tokens.create(user)['token']}"}
    created = c.post("/api/omim/captures", json={"mim": "600002", "text": PAGE}, headers=headers).json()
    summary = created["summary"]
    assert summary["db_enabled"] is True
    assert summary["adopt_preview"]["genes"][0]["current_labels"] == ["Disease B (600002)(AR)"]
    mim, name = created["capture"]["id"].split("/")
    adopted = c.post(f"/api/omim/captures/{mim}/{name}/adopt", headers=headers).json()
    assert adopted["result"]["text"] == "added"
    assert adopted["summary"]["adopt_preview"]["same_content"] is True
    rec = omim_store.lookup(gene="GENE1")
    assert rec["Disease2"] == (
        "Disease B renamed (600002)(AR)\n\nDisease B is made up.\n\nINHERITANCE\n- Autosomal recessive"
    )
    # duplicate capture returns the stored id, which can be adopted again
    dup = c.post("/api/omim/captures", json={"mim": "600002", "text": PAGE}, headers=headers).json()
    assert dup["capture"]["duplicate"] is True and dup["capture"]["id"] == created["capture"]["id"]
    status = c.get("/api/omim/db/status").json()
    assert status["enabled"] and status["source"] == "db" and status["counts"]["phenotype_texts"] == 1
    reviews = c.get("/api/omim/reviews").json()
    assert reviews["db_enabled"] and reviews["items"]
    item = reviews["items"][0]
    assert c.post(f"/api/omim/reviews/{item['id']}", json={"action": "keep"},
                  headers=headers).json()["resolution"] == "keep"
    assert c.post(f"/api/omim/reviews/{item['id']}", json={"action": "bogus"},
                  headers=headers).status_code == 400
    export = c.get("/api/omim/export.xlsx")
    assert export.status_code == 200 and export.content[:2] == b"PK"


def test_api_adopt_requires_database(client, monkeypatch, tmp_path):
    c, user = client
    headers = {"Authorization": f"Bearer {capture_tokens.create(user)['token']}"}
    created = c.post("/api/omim/captures", json={"mim": "600002", "text": PAGE}, headers=headers).json()
    monkeypatch.setattr(omim_db, "OMIM_DB", tmp_path / "missing.sqlite")
    mim, name = created["capture"]["id"].split("/")
    assert c.post(f"/api/omim/captures/{mim}/{name}/adopt", headers=headers).status_code == 409
    assert c.get("/api/omim/reviews").json() == {"items": [], "open_counts": {}, "db_enabled": False}


def test_pending_list_hides_adopted_and_bulk_adopt_skips_blocked(client):
    c, user = client
    headers = {"Authorization": f"Bearer {capture_tokens.create(user)['token']}"}
    ok = c.post("/api/omim/captures", json={"mim": "600002", "text": PAGE}, headers=headers).json()
    collapsed = PAGE.replace("600002", "600001").replace("Disease B renamed", "Disease A").replace(
        "Toggle Dropdown\n\nINHERITANCE\n- Autosomal recessive\n▲ Close\n", "")
    blocked = c.post("/api/omim/captures", json={"mim": "600001", "text": collapsed}, headers=headers).json()
    pending = c.get("/api/omim/captures").json()
    assert pending["total"] == 2
    by_id = {e["id"]: e for e in pending["captures"]}
    assert by_id[blocked["capture"]["id"]]["blockers"] == ["Clinical Synopsis 沒有擷取到"]
    assert by_id[ok["capture"]["id"]]["blockers"] == []

    result = c.post("/api/omim/captures/adopt-pending", json={"ids": list(by_id)}).json()
    assert [a["id"] for a in result["adopted"]] == [ok["capture"]["id"]]
    assert result["skipped"][0]["reasons"] == ["Clinical Synopsis 沒有擷取到"]

    pending = c.get("/api/omim/captures").json()
    assert [e["id"] for e in pending["captures"]] == [blocked["capture"]["id"]]
    adoptions = c.get("/api/omim/adoptions").json()["adoptions"]
    assert adoptions[0]["capture_id"] == ok["capture"]["id"] and adoptions[0]["genes"] == ["GENE1"]
    again = c.post("/api/omim/captures/adopt-pending", json={"ids": [ok["capture"]["id"]]}).json()
    assert again["adopted"] == [] and again["skipped"][0]["reasons"] == ["已採用或已有較新的擷取"]


def test_single_adopt_keeps_synopsis_when_new_capture_missed_it(db):
    omim_db.adopt_capture(_phenotype_capture(), user="rev", capture_id="600002/a")
    collapsed = _phenotype_capture(desc="Updated.", synopsis="")
    collapsed["synopsis_status"] = "not_loaded"
    result = omim_db.adopt_capture(collapsed, user="rev", capture_id="600002/b")
    assert result["kept_previous_synopsis"] is True
    _by_id, by_gene = omim_db.load_rows(db["path"])
    assert cells(by_gene["GENE1"])[1].endswith("Updated.\n\nINHERITANCE\n- Autosomal recessive")


def test_old_database_gets_adoptions_table_backfilled(db):
    import sqlite3
    con = sqlite3.connect(db["path"])
    con.execute("DROP TABLE adoptions")
    con.execute("INSERT INTO change_log(at, user, action, phenotype_mim, capture_id) "
                "VALUES ('2026-11-10T00:00:00+00:00', 'rev', 'adopt_text', '600002', '600002/old')")
    con.commit()
    con.close()
    omim_db._schema_checked.clear()
    assert omim_db.adopted_capture_ids(db["path"]) == {"600002/old"}


def _gene_page(mim, symbol, rows):
    return {
        "mim": mim, "entry_type": "gene", "prefix": "*", "symbol": symbol, "title": symbol,
        "relationships": [
            {"phenotype": name, "phenotype_mim": pm, "inheritance": inh, "mapping_key": "3",
             "gene": symbol, "gene_mim": mim, "location": ""}
            for name, pm, inh in rows
        ],
        "description": "", "synopsis": "", "content_hash": f"gene-{mim}",
    }


@pytest.fixture
def queue_env(client, monkeypatch):
    from app.services import omim_queue
    monkeypatch.setattr(omim_queue.panel_deadzone, "disease_associated_genes",
                        lambda: frozenset({"ATXNX", "OTX9"}))
    monkeypatch.setattr(omim_queue, "reported_usage", lambda ttl=0: {"600002": 3})
    omim_queue._cache.update(sig=None)
    c, user = client
    headers = {"Authorization": f"Bearer {capture_tokens.create(user)['token']}"}
    return c, headers


def test_queue_is_gene_centred_and_ordered(queue_env):
    c, _headers = queue_env
    data = c.get("/api/omim/queue").json()
    genes = [i["gene"] for i in data["items"]]
    # nothing synced yet: panel genes needing a gene-page sync come first
    assert genes[:2] == ["ATXNX", "OTX9"]
    atxn = data["items"][0]
    assert atxn["status"] == "gene_sync"
    # gene page, then the label-only SCA1 (missing) before SCA10 (curated, no Done → stale)
    assert [(s["type"], s["mim"], s["kinds"]) for s in atxn["steps"]] == [
        ("gene", "100004", ["gene_page"]), ("phenotype", "164400", ["missing"]),
        ("phenotype", "603516", ["stale"]),
    ]
    assert [r["kind"] for r in atxn["reviews"]] == ["slot_not_in_omim_disease"]
    # non-panel gene with a disease to fill still gets its gene page first
    gene1 = next(i for i in data["items"] if i["gene"] == "GENE1")
    assert gene1["reported"] == 3
    assert [(s["type"], s["kinds"]) for s in gene1["steps"]] == [
        ("gene", ["gene_page"]), ("phenotype", ["missing"]), ("phenotype", ["stale"])]
    assert data["counts"]["panel_gene_sync"] == 2
    assert c.get("/api/omim/queue?status=review").json()["total"] >= 1
    assert c.get("/api/omim/queue?scope=bogus").status_code == 400


def test_next_stays_on_gene_and_gene_sync_puts_it_in_progress(queue_env):
    c, headers = queue_env
    first = c.get("/api/omim/queue/next", headers=headers).json()
    assert first["step"]["type"] == "gene" and first["gene"]["gene"] == "ATXNX"
    assert first["gene"]["total_steps"] == 3
    # adopt the gene page: OMIM lists SCA1 and a disease the platform lacked
    result = omim_db.adopt_capture(
        _gene_page("100004", "ATXNX", [("Spinocerebellar ataxia 1", "164400", "AD"),
                                       ("Brand new ataxia", "620001", "AR")]),
        user="rev", capture_id="100004/g",
    )
    # migration review settled; SCA10 row flagged as not listed by OMIM
    assert result["reviews"] == 1
    omim_store.ensure_loaded()
    nxt = c.get("/api/omim/queue/next?current=100004&gene=100004", headers=headers).json()
    assert nxt["continuing"] is True and nxt["gene"]["gene"] == "ATXNX"
    assert nxt["step"]["mim"] in {"164400", "620001"}
    item = c.get("/api/omim/queue").json()["items"][0]
    assert item["gene"] == "ATXNX" and item["status"] == "in_progress"
    assert {s["mim"] for s in item["steps"]} == {"164400", "620001", "603516"}
    assert [r["kind"] for r in item["reviews"]] == ["not_in_omim"]
    # skipping both diseases leaves the gene; next gene comes up
    other = c.get("/api/omim/queue/next?current=164400&gene=100004&skip=620001,603516",
                  headers=headers).json()
    assert other["gene"]["gene"] != "ATXNX" and other["continuing"] is False
    # skipping a whole gene with g<MIM>
    skipped = c.get("/api/omim/queue/next?skip=g100004,g100002", headers=headers).json()
    assert skipped["gene"]["gene"] not in {"ATXNX", "OTX9"}


def test_gene_sync_resolves_reviews_for_listed_diseases(queue_env):
    c, headers = queue_env
    # CERTX: stale label without MIM, curated cell already new name + MIM
    omim_db.adopt_capture(_gene_page("100003", "CERTX", [
        ("Intellectual developmental disorder, autosomal dominant 34", "616351", "AD")]),
        user="rev", capture_id="100003/g")
    reviews = omim_db.open_reviews_by_gene().get("100003", [])
    assert reviews == []
    # GENE1: OMIM lists A and an unknown C, but not B (renamed elsewhere)
    result = omim_db.adopt_capture(_gene_page("100001", "GENE1", [
        ("Disease A", "600001", "AD"), ("Disease B renamed", "600002", "AR")]),
        user="rev", capture_id="100001/g")
    assert result["reviews"] == 0


def test_adopt_response_lists_gene_reviews_and_token_can_resolve(queue_env, tmp_path):
    c, headers = queue_env
    page = (
        "OMIM\n*100004\nTable of Contents\nTitle\nGene-Phenotype Relationships\n"
        "▼ External Links\n* 100004\nATAXIN X; ATXNX\n\nGene-Phenotype Relationships\n\n"
        "Location\tPhenotype\tPhenotype\nMIM number\tInheritance\tPhenotype\nmapping key\n"
        "6p22\tSpinocerebellar ataxia 1\t164400\tAD\t3\n\n▼ TEXT\n\nText.\n\n"
        "Edit History: carol : 08/14/2026\n"
    )
    created = c.post("/api/omim/captures", json={"mim": "100004", "text": page}, headers=headers).json()
    mim, name = created["capture"]["id"].split("/")
    adopted = c.post(f"/api/omim/captures/{mim}/{name}/adopt", headers=headers).json()
    reviews = adopted["result"]["gene_reviews"]
    assert [r["kind"] for r in reviews] == ["not_in_omim"]
    assert reviews[0]["current_label"].startswith("Spinocerebellar ataxia 10")
    done = c.post(f"/api/omim/reviews/{reviews[0]['id']}", json={"action": "remove"}, headers=headers)
    assert done.json()["resolution"] == "remove"


def test_old_database_gets_legacy_done_column(db):
    import sqlite3
    con = sqlite3.connect(db["path"])
    con.execute("CREATE TABLE genes_old AS SELECT gene_mim, gene_symbol, alias_symbols, inheritance, "
                "done, row_order, updated_at FROM genes")
    con.execute("DROP TABLE genes")
    con.execute("ALTER TABLE genes_old RENAME TO genes")
    con.execute("UPDATE genes SET done = '2020-01-01' WHERE gene_symbol = 'GENE1'")
    con.commit()
    con.close()
    omim_db._schema_checked.clear()
    _by_id, by_gene = omim_db.load_rows(db["path"])
    assert by_gene["GENE1"]["Legacy_done"] == "2020-01-01"


def test_preview_reports_mim_filled_for_old_label(db):
    omim_db.adopt_capture(_gene_page("100004", "ATXNX", [("Mystery", "", "")]),
                          user="rev", capture_id="100004/a")
    preview = omim_db.preview_adoption(_gene_page("100004", "ATXNX", [("Mystery", "620002", "AR")]))
    row = preview["genes"][0]
    assert row["current_labels"] == [] and row["fills_mim_for"] == "Mystery"


def test_shared_mim_adoption_does_not_hide_other_stale_text(db):
    """Adopting a disease shared by two genes moves both genes' Done date,
    but the other gene's untouched curated text is still judged by the
    workbook Done (here: none → no date, not "up to date")."""
    import sqlite3
    con = sqlite3.connect(db["path"])
    con.execute("UPDATE genes SET legacy_done = '2020-01-01', done = '2020-01-01' WHERE gene_symbol = 'GENE1'")
    con.commit()
    con.close()
    # 600002 is linked to GENE1; adopting it moves GENE1's Done to today.
    omim_db.adopt_capture(_phenotype_capture(), user="rev", capture_id="600002/x")
    omim_store._state.update(mtime=None)
    status = omim_capture.compare_with_workbook("600001", edit_date="2025-06-01")
    link = status["phenotype_links"][0]
    assert link["done"] == omim_db._today()
    assert link["status"] == "outdated"
