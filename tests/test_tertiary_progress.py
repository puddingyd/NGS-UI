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


def test_postprocessing_progress_is_weighted_by_step_time():
    pp = tp.PostprocessingProgress(sample_kinds=["nckuh:WES"] * 2, profile={"version": 1, "kinds": {}})
    weights = pp._weights("nckuh:WES")
    assert "spliceai" not in weights and "cnv-rescue" not in weights
    total = sum(weights.values())
    pp.start_sample(0)
    pp.on_marker("clinvar-latest", "start", now=0)
    pp.on_marker("clinvar-latest", "done", elapsed=10, now=10)
    pp.on_marker("genebe", "start", now=10)
    expected = (weights["clinvar-latest"] + weights["genebe"] * 0.5) / total / 2
    assert pp.fraction(now=10 + weights["genebe"] * 0.5) == pytest.approx(expected)
    pp.start_sample(1)
    assert pp.fraction(now=0) == pytest.approx(0.5)
    assert pp.step_seconds["clinvar-latest"] == [10.0]


def test_phase_bands_follow_learned_wall_time():
    default = tp.phase_bands({"version": 1, "kinds": {}}, ["nckuh:WES"] * 10, with_prs=False)
    assert default["pre"] == (0.0, 2.0)
    assert default["promote"] == (97.0, 99.0)
    assert 45 < default["nextflow"][1] < 52  # ~even split from the measured batch
    learned = {"version": 1, "kinds": {"nckuh:WES": {"phases": {"nextflow": 900, "postprocessing": 100}}}}
    bands = tp.phase_bands(learned, ["nckuh:WES"] * 10, with_prs=True)
    assert bands["nextflow"][1] == pytest.approx(2 + 92 * 0.9)
    assert bands["prs"] == (95.0, 99.0)
    assert tp.band_pct(bands, "postprocessing", 0.5) == pytest.approx(round((bands["postprocessing"][0] + 94) / 2, 1))


def test_update_profile_learns_from_trace_and_round_trips(tmp_path):
    rows = [
        {"name": "A:VEP_ANNOTATE (1)", "status": "COMPLETED", "realtime": "10m"},
        {"name": "A:VEP_ANNOTATE (2)", "status": "COMPLETED", "realtime": "20m"},
        {"name": "A:MITO_VEP (1)", "status": "CACHED", "realtime": "1m"},
    ]
    profile = tp.update_profile({"version": 1, "kinds": {}}, "nckuh:WES", trace_rows=rows, sample_count=2,
                                nextflow_seconds=2000, postprocessing_seconds=400,
                                step_seconds={"genebe": [100, 140]})
    entry = profile["kinds"]["nckuh:WES"]
    assert entry["processes"]["VEP_ANNOTATE"] == {"seconds": 900.0, "tasks_per_sample": 1.0}
    assert "seconds" not in entry["processes"]["MITO_VEP"]  # cached tasks carry no timing
    assert entry["phases"] == {"nextflow": 1000.0, "postprocessing": 200.0}
    assert entry["postprocessing_steps"]["genebe"] == 120.0
    path = tmp_path / "_timing_profile.json"
    tp.save_profile(path, profile)
    again = tp.update_profile(tp.load_profile(path), "nckuh:WES", sample_count=2, nextflow_seconds=1000)
    assert again["kinds"]["nckuh:WES"]["phases"]["nextflow"] == pytest.approx(750.0)  # EMA
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
