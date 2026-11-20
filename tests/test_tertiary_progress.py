from pathlib import Path

import pytest

from app.services import tertiary_progress as tp


def test_simple_process_name_handles_full_and_truncated_names():
    assert tp.simple_process_name("ANNOTATE_SNV:SNV_ANNOTATE:VEP_ANNOTATE (23WE0128)") == "VEP_ANNOTATE"
    # Nextflow's ANSI view truncates from the left; the old token match missed these.
    assert tp.simple_process_name("ANN…NNOTATE:VEP_ANNOTATE (80)") == "VEP_ANNOTATE"
    assert tp.simple_process_name("ANNOTATE_SNV:ACMG_CLASSIFY") == "ACMG_CLASSIFY"


def test_parse_stdout_task_lines_from_plain_log():
    assert tp.parse_stdout_task(
        "[6c/2b4d1e] Submitted process > ANNOTATE_SNV:SNV_ANNOTATE:VEP_ANNOTATE (12)"
    ) == ("6c/2b4d1e", "submitted", "VEP_ANNOTATE")
    assert tp.parse_stdout_task(
        "[aa/bbbbbb] Cached process > MITO_ANNOTATE:MITO_VEP (3)"
    ) == ("aa/bbbbbb", "cached", "MITO_VEP")
    assert tp.parse_stdout_task("executor >  local (6)") is None


@pytest.mark.parametrize("text,seconds", [
    ("1h 2m 3s", 3723.0), ("4.5s", 4.5), ("350ms", 0.35), ("2d 1h", 176400.0),
    ("-", None), ("", None), ("12000", 12.0),
])
def test_parse_duration(text, seconds):
    assert tp.parse_duration(text) == seconds


TRACE_HEADER = "task_id\thash\tnative_id\tname\tstatus\texit\tsubmit\tduration\trealtime\n"


def _trace_row(task_id, task_hash, name, status, realtime="10s"):
    return f"{task_id}\t{task_hash}\t1\t{name}\t{status}\t0\t-\t{realtime}\t{realtime}\n"


def test_read_trace_is_incremental_and_waits_for_complete_lines(tmp_path):
    trace = tmp_path / "trace.txt"
    trace.write_text(TRACE_HEADER + _trace_row(1, "aa/000001", "P:A (1)", "COMPLETED") + "2\tbb/00")
    rows, offset = tp.read_trace(trace, 0)
    assert [r["hash"] for r in rows] == ["aa/000001"]
    with trace.open("a") as handle:
        handle.write("0002\t1\tP:A (2)\tCOMPLETED\t0\t-\t5s\t5s\n")
    rows, _ = tp.read_trace(trace, offset)
    assert [r["hash"] for r in rows] == ["bb/000002"]
    assert rows[0]["realtime"] == "5s"


def _profile(seconds):
    return {"version": 1, "kinds": {"nckuh:WES": {"processes": {
        name: {"seconds": secs, "tasks_per_sample": 1.0} for name, secs in seconds.items()
    }}}}


def test_vep_dominates_nextflow_progress_and_running_tasks_earn_partial_credit():
    profile = _profile({"ADD_CALLERS_TAG": 10, "MITO_VEP": 10, "VEP_ANNOTATE": 800, "ACMG_CLASSIFY": 180})
    nf = tp.NextflowProgress(kinds=["nckuh:WES"] * 2, profile=profile, with_pgx=False)
    now = 1000.0
    for i, name in enumerate(["ADD_CALLERS_TAG", "MITO_VEP"]):
        for sample in (1, 2):
            h = f"{i}{sample}/000000"
            nf.on_stdout(f"[{h}] Submitted process > X:{name} ({sample})", now)
            nf.on_trace_row({"hash": h, "name": f"X:{name} ({sample})", "status": "COMPLETED"}, now)
    small_done = nf.fraction(now)
    assert small_done == pytest.approx(40 / 2000)
    nf.on_stdout("[ee/000001] Submitted process > ANNOTATE_SNV:SNV_ANNOTATE:VEP_ANNOTATE (1)", now)
    halfway = nf.fraction(now + 400)
    assert halfway == pytest.approx((40 + 400) / 2000)
    nf.on_trace_row({"hash": "ee/000001", "name": "ANNOTATE_SNV:SNV_ANNOTATE:VEP_ANNOTATE (1)",
                     "status": "COMPLETED"}, now + 800)
    assert nf.fraction(now + 800) == pytest.approx(840 / 2000)
    assert nf.current()["process"] == "VEP_ANNOTATE"
    assert nf.current()["done"] == 1 and nf.current()["total"] == 2


def test_cached_tasks_leave_the_denominator_and_done_event_fires_once():
    nf = tp.NextflowProgress(kinds=["nckuh:WES"] * 2, profile=_profile({"VEP_ANNOTATE": 100}), with_pgx=False)
    first = nf.on_stdout("[aa/000001] Cached process > S:VEP_ANNOTATE (1)", 0)
    assert [e["event"] for e in first] == ["start"]
    nf.on_stdout("[bb/000002] Submitted process > S:VEP_ANNOTATE (2)", 0)
    assert nf.fraction(0) == 0.0
    events = nf.on_trace_row({"hash": "bb/000002", "name": "S:VEP_ANNOTATE (2)", "status": "COMPLETED"}, 60)
    assert events == [{"process": "VEP_ANNOTATE", "event": "done", "done": 2, "total": 2, "elapsed": 60}]
    assert nf.fraction(60) == pytest.approx(0.99)  # capped until the phase ends
    assert nf.on_trace_row({"hash": "bb/000002", "name": "S:VEP_ANNOTATE (2)", "status": "COMPLETED"}, 61) == []


def test_default_table_matches_pipeline_and_seq_type():
    nf = tp.NextflowProgress(kinds=["nckuh:WES"], profile={"version": 1, "kinds": {}}, with_pgx=False)
    names = set(nf._processes)
    assert {"VEP_ANNOTATE", "PANGOLIN_SCORE", "ANNOTSV_CNV_NCKUH_WES", "ADD_CALLERS_TAG"} <= names
    assert not names & {"ANNOTSV_CNV_DRAGEN", "ANNOTSV_CNV_NCKUH_WGS", "CNVKIT_TO_BED", "PGX_PHARMCAT"}


def test_fit_linear_separates_fixed_and_per_sample_cost():
    default = (100.0, 10.0)
    assert tp.fit_linear([], default) == default
    # one batch size: keep the default shape, scale to the observation
    a, b = tp.fit_linear([[10, 400]], default)
    assert a + 10 * b == pytest.approx(400) and a / b == pytest.approx(10)
    # two sizes: least squares (GeneBe: 330 s for 1 sample, 843 s for 36)
    a, b = tp.fit_linear([[1, 330], [36, 843]], default)
    assert (round(a), round(b, 2)) == (315, 14.66)
    # never negative
    assert tp.fit_linear([[1, 500], [10, 100]], default) == (300.0, 0.0)


def test_phase_bands_depend_on_batch_size():
    empty = {"version": 1, "kinds": {}}
    one = tp.phase_bands(empty, ["nckuh:WES"], with_prs=False)
    many = tp.phase_bands(empty, ["nckuh:WES"] * 36, with_prs=False)
    assert one["pre"] == (0.0, 2.0) and one["promote"][1] == 99.0
    # 1 sample: fixed post-processing costs dominate; 36: Nextflow dominates
    assert 40 < one["nextflow"][1] < 50
    assert 80 < many["nextflow"][1] < 86
    with_prs = tp.phase_bands(empty, ["nckuh:WES"] * 36, with_prs=True)
    assert with_prs["promote"][1] == 95.0 and with_prs["prs"] == (95.0, 99.0)
    learned = {"version": 1, "kinds": {"nckuh:WES": {"observations": {
        "phase:nextflow": [[1, 600], [36, 17000]], "phase:postprocessing": [[1, 100], [36, 400]],
        "phase:promote": [[1, 5], [36, 100]]}}}}
    bands = tp.phase_bands(learned, ["nckuh:WES"] * 36, with_prs=False)
    assert bands["nextflow"][1] == pytest.approx(2 + 97 * 17000 / 17500, abs=0.2)
    assert tp.band_pct(bands, "postprocessing", 0.5) == pytest.approx(
        round((bands["postprocessing"][0] + bands["postprocessing"][1]) / 2, 1))


def test_postprocessing_progress_uses_chunk_cost_and_parallel_sample_steps():
    pp = tp.PostprocessingProgress(chunks=[["nckuh:WES"] * 4], profile={"version": 1, "kinds": {}}, workers=2)
    plan = pp._chunk_plan(0)
    assert "spliceai" not in plan and "cnv-rescue" not in plan
    genebe = tp.estimate_seconds({}, "nckuh:WES", "step:genebe", 4)
    assert plan["genebe"] == (genebe, genebe)
    assert plan["review-tsv"] == (2.0 * 4 / 2, 2.0)  # 4 samples, 2 at a time
    total = sum(w for w, _p in plan.values())
    pp.start_chunk(0)
    pp.on_marker("genebe", "start", now=0)
    assert pp.fraction(now=genebe / 2) == pytest.approx(genebe / 2 / total)
    pp.on_marker("genebe", "done", elapsed=300, now=300)
    assert pp.step_observations == [("genebe", 4.0, 300.0)]
    pp.on_marker("review-tsv", "start", sample="A", now=300)
    pp.on_marker("review-tsv", "done", sample="A", elapsed=3, now=303)
    assert pp.sample_seconds["review-tsv"] == [3.0]
    assert pp.fraction(now=303) == pytest.approx((genebe + 4.0 / 4) / total)
    pp.finish_chunk(0)
    assert pp.fraction() == pytest.approx(1.0)


def _row(name, status, realtime="10m"):
    return {"name": f"A:{name} (1)", "status": status, "realtime": realtime}


def test_cached_share_is_weighted_by_time_not_task_count():
    profile = {"version": 1, "kinds": {}}
    rows = [_row(p, "CACHED") for p in ("ADD_CALLERS_TAG", "MITO_VEP", "MITO_PARSE", "STR_PARSE_NCKUH",
                                        "ANNOTSV_SV_NCKUH", "FILTER_FOR_ANNOTATION", "ANNOTSV_CNV_NCKUH_WES",
                                        "STR_PREPARE_NCKUH")]
    rows += [_row(p, "COMPLETED") for p in ("VEP_ANNOTATE", "PANGOLIN_SCORE", "PARSE_CSQ", "ACMG_CLASSIFY")]
    # 8 of 12 tasks were cached, but almost all of the expected time really ran.
    assert tp.noncached_share(profile, "nckuh:WES", rows) > 0.8
    tp.update_profile(profile, "nckuh:WES", trace_rows=rows, sample_count=1, nextflow_seconds=587)
    obs = profile["kinds"]["nckuh:WES"]["observations"]["phase:nextflow"]
    assert len(obs) == 1 and obs[0][1] == 587.0 and 0.8 < obs[0][0] <= 1
    all_cached = [_row("VEP_ANNOTATE", "CACHED")]
    tp.update_profile(profile, "nckuh:WES", trace_rows=all_cached, sample_count=1, nextflow_seconds=6)
    assert len(profile["kinds"]["nckuh:WES"]["observations"]["phase:nextflow"]) == 1


def test_update_profile_learns_from_trace_and_round_trips(tmp_path):
    rows = [
        {"name": "A:VEP_ANNOTATE (1)", "status": "COMPLETED", "realtime": "10m"},
        {"name": "A:VEP_ANNOTATE (2)", "status": "COMPLETED", "realtime": "20m"},
        {"name": "A:MITO_VEP (1)", "status": "CACHED", "realtime": "1m"},
    ]
    old = {"version": 1, "kinds": {"nckuh:WES": {"phases": {"nextflow": 1}, "postprocessing_batch_steps": {}}}}
    profile = tp.update_profile(old, "nckuh:WES", trace_rows=rows, sample_count=2,
                                nextflow_seconds=2000, postprocessing_seconds=400, promote_seconds=20,
                                step_observations=[("genebe", 2, 330.0), ("bogus", 2, 1.0)],
                                sample_step_seconds_seen={"review-tsv": [2, 4]})
    entry = profile["kinds"]["nckuh:WES"]
    assert entry["processes"]["VEP_ANNOTATE"] == {"seconds": 900.0, "tasks_per_sample": 1.0}
    assert "seconds" not in entry["processes"]["MITO_VEP"]  # cached tasks carry no timing
    assert entry["observations"]["phase:postprocessing"] == [[2.0, 400.0]]
    assert entry["observations"]["phase:promote"] == [[2.0, 20.0]]
    assert entry["observations"]["step:genebe"] == [[2.0, 330.0]]
    assert "step:bogus" not in entry["observations"]
    assert entry["sample_steps"]["review-tsv"] == 3.0
    assert "phases" not in entry and "postprocessing_batch_steps" not in entry  # superseded keys dropped
    path = tmp_path / "_timing_profile.json"
    tp.save_profile(path, profile)
    loaded = tp.load_profile(path)
    for _ in range(20):
        tp.update_profile(loaded, "nckuh:WES", sample_count=1, postprocessing_seconds=700)
    assert len(loaded["kinds"]["nckuh:WES"]["observations"]["phase:postprocessing"]) == tp.MAX_OBSERVATIONS
    assert tp.load_profile(tmp_path / "missing.json") == {"version": 1, "kinds": {}}


def test_eta_needs_meaningful_progress():
    assert tp.eta_seconds(0, 3, now=100) is None
    assert tp.eta_seconds(0, 50, now=3600) == 3600


def test_worker_progress_is_monotonic_and_prs_maps_onto_overall_pct(monkeypatch):
    from app.workers import dragen_run
    store = {}
    monkeypatch.setattr(dragen_run.dragen_jobs, "load_state", lambda jid: dict(store))
    monkeypatch.setattr(dragen_run.dragen_jobs, "save_state", lambda jid, st: store.update(st))
    dragen_run._update("job", progress_pct=40.0)
    dragen_run._update("job", progress_pct=35.0)  # late ticker update must not move backwards
    assert store["progress_pct"] == 40.0
    dragen_run._update("job", prs_progress_pct=96.5)
    assert store["progress_pct"] == 96.5


def test_postprocessing_chunks_respect_scratch_budget(tmp_path, monkeypatch):
    from app.workers import dragen_run
    sizes = {"A": 10, "B": 10, "C": 30, "D": 5}
    raw = {}
    for sid, size in sizes.items():
        raw[sid] = tmp_path / f"{sid}.tsv"
        raw[sid].write_bytes(b"x" * size)
    samples = [{"sample_id": sid} for sid in sizes]
    monkeypatch.setenv("NGS_UI_POSTPROC_CHUNK_GB", str(25 / 1024 ** 3))
    chunks = dragen_run._postprocessing_chunks(samples, raw)
    assert [[s["sample_id"] for s in c] for c in chunks] == [["A", "B"], ["C"], ["D"]]
    monkeypatch.delenv("NGS_UI_POSTPROC_CHUNK_GB")
    assert len(dragen_run._postprocessing_chunks(samples, raw)) == 1


def test_marker_regex_reads_sample_and_elapsed():
    from app.workers import dragen_run
    m = dragen_run._POSTPROC_MARKER_RE.search("[t] [sample-step] review-tsv done sample=S1-nckuh elapsed=12s")
    assert m.groups() == ("sample-step", "review-tsv", "done", "S1-nckuh", "12")
    m = dragen_run._POSTPROC_MARKER_RE.search("[t] [post-processing-step] genebe done elapsed=4s")
    assert m.groups() == ("post-processing-step", "genebe", "done", None, "4")


def test_log_summary_attributes_sample_steps_inside_a_batch():
    from app.services import postprocessing_log
    events = []
    summary = postprocessing_log.Summary("batch 1/1 (2 samples)", "postprocessing", events.append)
    summary.feed("[t] [post-processing-step] genebe start")
    summary.feed("[t] [post-processing-step] genebe done elapsed=4s")
    summary.feed("[t] [sample-step] review-tsv start sample=S1")
    summary.feed("[t] [sample-step] review-tsv start sample=S2")
    summary.feed("[t] [sample-step] review-tsv done sample=S1 elapsed=3s")
    summary.finish(0)
    assert "sample=batch 1/1 (2 samples) step=genebe status=done elapsed=4.0s" in events[1]
    assert any("sample=S1 step=review-tsv status=done elapsed=3.0s" in e for e in events)
    assert any("sample=S2 step=review-tsv status=done" in e for e in events)


def test_snapshot_lists_processes_in_pipeline_order_like_the_console_block():
    profile = _profile({"ADD_CALLERS_TAG": 10, "VEP_ANNOTATE": 800, "ACMG_CLASSIFY": 100})
    nf = tp.NextflowProgress(kinds=["nckuh:WES"] * 2, profile=profile, with_pgx=False)
    nf.on_stdout("[aa/000001] Cached process > PREPARE_VCF:ADD_CALLERS_TAG (1)", 0)
    nf.on_stdout("[aa/000002] Cached process > PREPARE_VCF:ADD_CALLERS_TAG (2)", 0)
    nf.on_stdout("[bb/000001] Submitted process > ANNOTATE_SNV:SNV_ANNOTATE:VEP_ANNOTATE (1)", 0)
    rows = nf.snapshot()
    assert [r["process"] for r in rows] == ["ADD_CALLERS_TAG", "VEP_ANNOTATE", "ACMG_CLASSIFY"]
    assert rows[0] | {} == {"process": "ADD_CALLERS_TAG", "name": "PREPARE_VCF:ADD_CALLERS_TAG", "done": 2,
                            "total": 2, "cached": 2, "running": 0, "failed": 0, "status": "done"}
    assert (rows[1]["name"], rows[1]["running"], rows[1]["status"]) == (
        "ANNOTATE_SNV:SNV_ANNOTATE:VEP_ANNOTATE", 1, "running")
    assert rows[2]["status"] == "pending" and rows[2]["total"] == 2
    nf.on_trace_row({"hash": "bb/000001", "name": "ANNOTATE_SNV:SNV_ANNOTATE:VEP_ANNOTATE (1)",
                     "status": "FAILED"}, 5)
    assert nf.snapshot()[1]["status"] == "failed"
    # after the run, processes that never started are dropped
    assert [r["process"] for r in nf.snapshot(finished=True)] == ["ADD_CALLERS_TAG", "VEP_ANNOTATE"]
