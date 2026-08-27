import ast
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_standalone_phenotype_has_matching_patient_action_buttons_and_gc_modal():
    html = (ROOT / "frontend" / "phenotype" / "index.html").read_text(encoding="utf-8")
    js = (ROOT / "frontend" / "phenotype" / "app.js").read_text(encoding="utf-8")
    css = (ROOT / "frontend" / "phenotype" / "style.css").read_text(encoding="utf-8")

    assert 'id="btn-load"' in html
    assert 'id="btn-emr-link"' in html
    assert 'id="btn-gc-records"' in html
    assert html.count('class="btn btn-secondary"') >= 3
    assert 'id="gc-record-modal"' in html
    assert "/consultation" in js
    assert 'credentials: "same-origin"' in js
    assert ".patient-action-buttons .btn" in css
    assert ".gc-modal[hidden]" in css


def test_new_case_form_distinguishes_untouched_from_explicitly_empty_phenotype():
    js = (ROOT / "frontend" / "app.js").read_text(encoding="utf-8")

    assert "phenotype_explicit" in js
    assert "newCaseEdit.edited" in js


def test_main_save_and_gc_endpoint_use_shared_backend_services():
    phenotype_router = (ROOT / "backend" / "app" / "routers" / "phenotype.py").read_text(encoding="utf-8")
    emr_router = (ROOT / "backend" / "app" / "routers" / "emr.py").read_text(encoding="utf-8")

    assert "patient_phenotype_store.save(" in phenotype_router
    assert '"patient_phenotype": patient_snapshot' in phenotype_router
    assert 'if target_version == "default"' in phenotype_router
    assert '@router.get("/emr/{mrn}/consultation")' in emr_router
    assert "emr_client.fetch_consultation(safe_mrn)" in emr_router


def test_only_default_main_analysis_syncs_patient_snapshot():
    source = (ROOT / "backend" / "app" / "routers" / "phenotype.py").read_text(
        encoding="utf-8",
    )
    tree = ast.parse(source)
    update = next(
        node for node in tree.body
        if isinstance(node, ast.FunctionDef) and node.name == "update_phenotype"
    )
    default_guard = next(
        node for node in ast.walk(update)
        if isinstance(node, ast.If)
        and isinstance(node.test, ast.Compare)
        and isinstance(node.test.left, ast.Name)
        and node.test.left.id == "target_version"
        and any(
            isinstance(comparator, ast.Constant) and comparator.value == "default"
            for comparator in node.test.comparators
        )
    )

    def is_patient_save(node):
        return (
            isinstance(node, ast.Call)
            and isinstance(node.func, ast.Attribute)
            and node.func.attr == "save"
            and isinstance(node.func.value, ast.Name)
            and node.func.value.id == "patient_phenotype_store"
        )

    all_saves = [node for node in ast.walk(update) if is_patient_save(node)]
    guarded_saves = [node for node in ast.walk(default_guard) if is_patient_save(node)]
    assert len(all_saves) == len(guarded_saves) == 1
    assert '"patient_phenotype_synced": patient_snapshot is not None' in source

    frontend = (ROOT / "frontend" / "app.js").read_text(encoding="utf-8")
    assert "result.patient_phenotype_synced" in frontend
    assert "此組合不影響病人 default phenotype" in frontend
