import json
import threading

import pytest
from fastapi import HTTPException

from app import config
from app.services import category_store as categories, report_store, sample_layout, sample_loader


@pytest.fixture
def cases(tmp_path, monkeypatch):
    from app.workers import dragen_run
    monkeypatch.setattr(dragen_run, "TERTIARY_JOBS_DIR", tmp_path / "jobs")
    monkeypatch.setattr(config, "DATA_ROOT", tmp_path / "data")
    monkeypatch.setattr(config, "PIPELINE_OUT_ROOT", tmp_path / "unified")
    monkeypatch.setattr(config, "LEGACY_TERTIARY_OUTPUT_ROOT", tmp_path / "legacy")
    monkeypatch.setattr(config, "LEGACY_PIPELINE_OUT_ROOT", tmp_path / "old_pipeline")
    monkeypatch.setattr(sample_loader, "update_case_table_row", lambda sid: None)

    def create(sid, category, *, unified=False):
        directory = config.PIPELINE_OUT_ROOT / sid / "08_postprocessing" if unified else config.LEGACY_TERTIARY_OUTPUT_ROOT / sid
        directory.mkdir(parents=True)
        if unified:
            (directory / f"{sid}.layout.json").write_text('{"layout_version": 3}')
        filename = f"{sid}.sample_metadata.json" if unified else "sample_metadata.json"
        path = directory / filename
        path.write_text(json.dumps({"sample_id": sid, "category": category, "comment": "keep", "edits": {"v1": {"comment": "keep too"}}}))
        return path
    return create


def edit(snapshot, mapping=None, add=()):
    mapping = mapping or {}
    return {
        "category_revision": snapshot["category_revision"],
        "categories": [
            {"original": name, "name": mapping.get(name, name)}
            for name in snapshot["category_options"] if mapping.get(name, name) is not None
        ] + [{"original": None, "name": name} for name in add],
    }


def test_usage_counts_registered_cases_once_and_includes_legacy_categories(cases):
    cases("A", "GI", unified=True)
    cases("B", "GI")
    path = cases("C", "Legacy name")
    # Same ID in a fallback root and extra analysis versions do not add cases.
    cases("A", "Neurology")
    (path.parent / "analyses" / "v2").mkdir(parents=True)
    result = categories.options()
    assert result["category_options"][:2] == ["GI", "Legacy name"]
    assert result["category_counts"]["GI"] == 2
    assert result["category_counts"]["Neurology"] == 0
    assert result["category_options"][2:4] == ["Neurology", "Endocrinology"]


def test_rename_delete_add_persist_and_preserve_other_report_fields(cases):
    a, b, c = cases("A", "GI", unified=True), cases("B", "MCA"), cases("C", "Neurology")
    originals = [json.loads(p.read_text()) for p in (a, b, c)]
    cache = sample_layout.global_cache_path("_case_table.json")
    cache.parent.mkdir(parents=True, exist_ok=True)
    cache.write_text("{}")
    result = categories.save(edit(categories.options(), {"GI": "Gastro", "MCA": None}, ["New"]))
    assert result["updated_samples"] == 2
    assert result["category_options"][0:2] == ["Neurology", "Gastro"]
    assert result["category_counts"]["Gastro"] == 1
    assert result["category_counts"]["New"] == 0
    for path, original, expected in zip((a, b, c), originals, ("Gastro", "", "Neurology")):
        original["category"] = expected
        assert json.loads(path.read_text()) == original
    assert not cache.exists()
    assert categories.options() == {k: result[k] for k in ("category_options", "category_counts", "category_revision")}
    assert not (config.DATA_ROOT / "categories.pending.json").exists()


def test_simultaneous_swaps_are_not_sequential_replacements(cases):
    a, b = cases("A", "GI"), cases("B", "MCA")
    categories.save(edit(categories.options(), {"GI": "MCA", "MCA": "GI"}))
    assert json.loads(a.read_text())["category"] == "MCA"
    assert json.loads(b.read_text())["category"] == "GI"


@pytest.mark.parametrize("names", [[""], ["GI", "gi"], ["a\nb"], ["a" * 101]])
def test_invalid_names_do_not_write_any_files(cases, names):
    path = cases("A", "GI")
    before = path.read_bytes()
    with pytest.raises(HTTPException) as exc:
        categories.save({"category_revision": categories.options()["category_revision"], "categories": [{"original": None, "name": n} for n in names]})
    assert exc.value.status_code == 400
    assert path.read_bytes() == before
    assert not (config.DATA_ROOT / "categories.json").exists()


def test_stale_manager_cannot_overwrite_another_edit(cases):
    cases("A", "GI")
    snapshot = categories.options()
    categories.save(edit(snapshot, {"GI": "Gastro"}))
    with pytest.raises(HTTPException) as exc:
        categories.save(edit(snapshot, {"GI": None}))
    assert exc.value.status_code == 409
    assert "Gastro" in categories.options()["category_options"]


def test_failed_write_rolls_back_every_case_and_settings(cases, monkeypatch):
    a, b = cases("A", "GI"), cases("B", "GI")
    originals = [p.read_bytes() for p in (a, b)]
    original_write = categories._write
    failed = False

    def write(path, value):
        nonlocal failed
        if path == b and not failed:
            failed = True
            raise OSError("synthetic disk error")
        original_write(path, value)
    monkeypatch.setattr(categories, "_write", write)
    with pytest.raises(OSError):
        categories.save(edit(categories.options(), {"GI": "Gastro"}))
    assert [json.loads(p.read_text()) for p in (a, b)] == [json.loads(data) for data in originals]
    assert not (config.DATA_ROOT / "categories.json").exists()
    assert "GI" in categories.options()["category_options"]


def test_interrupted_transaction_recovers_before_next_read(cases, monkeypatch):
    a, b = cases("A", "GI"), cases("B", "GI")
    original_write = categories._write

    def interrupted(path, value):
        if path == b:
            raise RuntimeError("synthetic process interruption")
        original_write(path, value)
    monkeypatch.setattr(categories, "_write", interrupted)
    with pytest.raises(RuntimeError):
        categories.save(edit(categories.options(), {"GI": "Gastro"}))
    assert json.loads(a.read_text())["category"] == "Gastro"
    monkeypatch.setattr(categories, "_write", original_write)
    assert categories.options()["category_counts"]["Gastro"] == 2
    assert json.loads(b.read_text())["category"] == "Gastro"


def test_old_report_cannot_resurrect_removed_category_or_undo_swap(cases, monkeypatch):
    cases("A", "GI")
    report = report_store.load("A")
    categories.save(edit(categories.options(), {"GI": "MCA", "MCA": "GI"}))
    with pytest.raises(HTTPException) as exc:
        report_store.save("A", report)
    assert exc.value.status_code == 409
    report = report_store.load("A")
    categories.save(edit(categories.options(), {"MCA": None}))
    with pytest.raises(HTTPException):
        report_store.save("A", {"category": "MCA"})
    assert report_store.load("A")["category"] == ""


def test_unreadable_registered_case_prevents_partial_update(cases):
    a, b = cases("A", "GI"), cases("B", "MCA")
    payload = edit(categories.options(), {"GI": "Gastro"})
    b.write_text("{bad json")
    with pytest.raises(ValueError):
        categories.save(payload)
    assert json.loads(a.read_text())["category"] == "GI"


def test_empty_list_clears_all_registered_categories(cases):
    a = cases("A", "GI")
    categories.save({"category_revision": categories.options()["category_revision"], "categories": []})
    assert categories.options()["category_options"] == []
    assert json.loads(a.read_text())["category"] == ""


def test_metadata_writers_wait_for_category_transaction(cases):
    cases("A", "GI")
    started, done = threading.Event(), threading.Event()

    @categories.serialized
    def writer():
        started.set()
        done.set()
    with categories.locked():
        thread = threading.Thread(target=writer)
        thread.start()
        assert not started.wait(0.05)
    thread.join(1)
    assert done.is_set()


def test_authenticated_api_updates_old_cases_and_retains_tag_options(cases, monkeypatch):
    from fastapi import FastAPI
    from fastapi.testclient import TestClient
    from starlette.middleware.sessions import SessionMiddleware
    from app.auth import current_user
    from app.routers import samples

    path = cases("A", "GI")
    sample_layout.global_cache_path("_options.json").parent.mkdir(parents=True, exist_ok=True)
    sample_layout.global_cache_path("_options.json").write_text(json.dumps({"tag_suggestions": ["keep tag"]}))
    app = FastAPI()
    app.add_middleware(SessionMiddleware, secret_key="synthetic-category-test")
    app.include_router(samples.router)
    client = TestClient(app)
    assert client.put("/api/options/categories", json={}).status_code == 401
    app.dependency_overrides[current_user] = lambda: {"id": 1, "username": "test"}
    snapshot = client.get("/api/options").json()
    assert snapshot["tag_suggestions"] == ["keep tag"]
    response = client.put("/api/options/categories", json=edit(snapshot, {"GI": "Gastro"}))
    assert response.status_code == 200
    assert client.get("/api/samples").json()[0]["category"] == "Gastro"
    assert client.get("/api/samples/A/report").json()["category"] == "Gastro"
    # Metadata edits share validation with report edits.
    stale = client.put("/api/samples/A/metadata", json={"category": "GI", "category_revision": snapshot["category_revision"]})
    assert stale.status_code == 409
    assert json.loads(path.read_text())["category"] == "Gastro"
    assert client.get("/api/options").json()["tag_suggestions"] == ["keep tag"]


def test_active_analysis_allows_category_changes(cases):
    from app.workers.dragen_run import _acquire_sample_locks, _release_sample_locks
    path = cases("A", "GI")
    payload = edit(categories.options(), {"GI": "Gastro"})
    handles = _acquire_sample_locks(["A"])
    try:
        result = categories.save(payload)
        assert result["updated_samples"] == 1
        assert json.loads(path.read_text())["category"] == "Gastro"
        assert not (config.DATA_ROOT / "categories.pending.json").exists()
    finally:
        _release_sample_locks(handles)


@pytest.mark.parametrize("layout_version", [2, 3])
def test_active_analysis_publication_and_rollback_preserve_category_changes(cases, tmp_path, layout_version):
    from app.workers import dragen_run
    path = cases("A", "GI", unified=True)
    post = path.parent
    marker = post / "A.layout.json"
    marker.write_text(json.dumps({"layout_version": layout_version}))
    if layout_version == 2:
        path.rename(post / "sample_metadata.json")
    staged = tmp_path / "staged" / "A"
    staged_post = staged / "08_postprocessing"
    staged_post.mkdir(parents=True)
    (staged_post / "A.layout.json").write_text('{"layout_version": 3}')
    handles = dragen_run._acquire_sample_locks(["A"])
    try:
        categories.save(edit(categories.options(), {"GI": "Gastro"}))
        operations = dragen_run._promote_staged_generation(
            sample_id="A", staged_sample_dir=staged,
            live_sample_dir=config.PIPELINE_OUT_ROOT / "A",
            rollback_sample_dir=tmp_path / "rollback" / "A",
        )
        assert report_store.load("A")["category"] == "Gastro"
        categories.save(edit(categories.options(), {"Gastro": None}))
        dragen_run._rollback_promotion_operations(operations)
        assert report_store.load("A")["category"] == ""
        assert report_store.load("A")["comment"] == "keep"
    finally:
        dragen_run._release_sample_locks(handles)


def test_publication_waits_for_category_commit_instead_of_whole_analysis(cases, tmp_path):
    from app.workers import dragen_run
    cases("A", "GI", unified=True)
    staged = tmp_path / "staged" / "A"
    staged_post = staged / "08_postprocessing"
    staged_post.mkdir(parents=True)
    staged_marker = staged_post / "A.layout.json"
    staged_marker.write_text('{"layout_version": 3}')
    started, done = threading.Event(), threading.Event()
    failures = []

    def publish():
        started.set()
        try:
            dragen_run._promote_staged_generation(
                sample_id="A", staged_sample_dir=staged,
                live_sample_dir=config.PIPELINE_OUT_ROOT / "A",
                rollback_sample_dir=tmp_path / "rollback" / "A",
            )
        except Exception as exc:
            failures.append(exc)
        finally:
            done.set()

    with categories.locked():
        thread = threading.Thread(target=publish)
        thread.start()
        assert started.wait(1)
        assert not done.wait(0.05)
        assert staged_marker.exists()
        categories.save(edit(categories.options(), {"GI": "Gastro"}))
    thread.join(2)
    assert done.is_set()
    assert not failures
    assert not staged_marker.exists()
    assert report_store.load("A")["category"] == "Gastro"
