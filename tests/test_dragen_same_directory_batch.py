from pathlib import Path

from app.services import dragen_jobs


REPO_ROOT = Path(__file__).resolve().parents[1]
APP_JS = (REPO_ROOT / "frontend" / "app.js").read_text(encoding="utf-8")
INDEX_HTML = (REPO_ROOT / "frontend" / "index.html").read_text(encoding="utf-8")


def _function_body(name: str) -> str:
    marker = f"function {name}("
    start = APP_JS.index(marker)
    brace = APP_JS.index(") {", start) + 2
    depth = 0
    for index in range(brace, len(APP_JS)):
        char = APP_JS[index]
        if char == "{":
            depth += 1
        elif char == "}":
            depth -= 1
            if depth == 0:
                return APP_JS[brace + 1:index]
    raise AssertionError(f"unterminated JavaScript function: {name}")


def test_same_directory_button_is_left_of_add_batch_and_wired():
    same_dir = INDEX_HTML.index('id="dragen-add-folder-btn"')
    add_batch = INDEX_HTML.index('id="dragen-add-batch-btn"')

    assert same_dir < add_batch
    assert "加入同目錄全部檢體" in INDEX_HTML[same_dir:add_batch]
    assert (
        'document.getElementById("dragen-add-folder-btn")?.addEventListener('
        '"click", _dragenAddDirectoryToBatch);'
    ) in APP_JS


def test_same_directory_action_uses_directory_key_defaults_and_deduplication():
    directory_body = _function_body("_dragenAddDirectoryToBatch")
    add_body = _function_body("_dragenAddSampleToBatch")

    assert "_dragenDirectoryKey(mode, current)" in directory_body
    assert "_dragenDirectoryKey(mode, row) === directory" in directory_body
    assert '_dragenSuggestSid(row.sample_id || "", mode)' in directory_body
    assert "_dragenAddSampleToBatch(_dragenBuildSample" in directory_body
    assert "row.source_sample_id === sample.source_sample_id" in add_body
    assert "row.vcf_path === sample.vcf_path" in add_body


def test_dragen_index_rows_share_the_run_input_directory(tmp_path, monkeypatch):
    root = tmp_path / "dragen"
    run_dir = root / "RUN-1"
    vcf_dir = run_dir / "vcf.gz"
    vcf_dir.mkdir(parents=True)
    (vcf_dir / "S1.hard-filtered.vcf.gz").touch()
    (vcf_dir / "S2.hard-filtered.vcf.gz").touch()
    monkeypatch.setattr(dragen_jobs, "DRAGEN_VCF_ROOTS", [root])

    rows = dragen_jobs.list_dragen_vcfs()

    assert {row["sample_id"] for row in rows} == {"S1", "S2"}
    assert {row["batch_dir"] for row in rows} == {str(run_dir)}


def test_inhouse_index_rows_share_the_parent_batch_directory(tmp_path, monkeypatch):
    root = tmp_path / "inhouse"
    batch_dir = root / "RUN-2"
    for sample_id in ("N1", "N2"):
        vcf_dir = batch_dir / sample_id / "04_snv_indel"
        vcf_dir.mkdir(parents=True)
        (vcf_dir / f"{sample_id}.ensemble.fixed.vcf.gz").touch()
    monkeypatch.setattr(dragen_jobs, "INHOUSE_VCF_ROOTS", [root])

    rows = dragen_jobs.list_inhouse_vcfs()

    assert {row["sample_id"] for row in rows} == {"N1", "N2"}
    assert {row["batch_dir"] for row in rows} == {str(batch_dir)}
