"""Synthetic CNV fixtures; no patient VCF data is checked into the repository."""
import csv
import gzip
import json
from pathlib import Path

import pytest

from app.services import dragen_cnv_rescue as rescue


def record(pos=1000, end=8000, *, identifier="cnv", kind="DEL", filt="cnvLength",
           qual="8", gt="0/1", chrom="chr1", svtype="CNV", **info):
    fields = [chrom, str(pos), identifier, "N", f"<{kind}>", qual, filt, "",
              "GT:CN:PR:SR", f"{gt}:1:30,14:23,10"]
    data = {"SVTYPE": svtype, "END": str(end), "SVLEN": str(end - pos), **info}
    fields[7] = ";".join(f"{k}={v}" for k, v in data.items())
    return rescue.Record(fields, data)


def paired(*, start=1000, end=8000, sv_start=4500, sv_end=8000, **kwargs):
    raw = record(start, end)
    cnv = record(start, end, SVCLAIM="DJ", RIGHT_BND="original-sv.end")
    sv = record(sv_start, sv_end, identifier="integrated-sv", filt="PASS", qual="150",
                svtype="DEL", SVCLAIM="J", END_RIGHT_BND_OF="cnv", **kwargs)
    return raw, cnv, sv


def write_vcf(path, rows, sample="SRC"):
    text = ("##fileformat=VCFv4.2\n"
            '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="type">\n'
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + sample + "\n")
    text += "".join("\t".join(row.fields) + "\n" for row in rows)
    if path.suffix == ".gz":
        with gzip.open(path, "wt") as handle:
            handle.write(text)
    else:
        path.write_text(text)


HEADER = ["AnnotSV_ID", "SV_chrom", "SV_start", "SV_end", "SV_type", "SV_length",
          "ID", "QUAL", "FILTER", "FORMAT", "SRC", "Annotation_mode", "Gene_name", "ACMG_class"]


def write_tsv(path, rows, headers=HEADER):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=headers, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def baseline_row():
    return dict(zip(HEADER, ["existing-id", "1", "100000", "200000", "DEL", "-100000",
                            "baseline-vcf-id", "50", "PASS", "GT:CN", "0/1:1", "full", "", "3"]))


def fake_annotate(vcf, output):
    _, records = rescue.read_vcf(vcf, "SRC")
    rows = []
    for r in records:
        # Simulate AnnotSV's full/split topology and arbitrary per-run AnnotSV_IDs.
        row = dict(zip(HEADER, ["unstable-annotsv-id", r.key[0], r.key[1], r.key[2], r.kind,
                               r.key[1] - r.key[2], r.id, r.fields[5], r.fields[6],
                               r.fields[8], r.fields[9], "full", "", "3"]))
        rows.extend([row, {**row, "Annotation_mode": "split", "Gene_name": "TESTGENE"}])
    write_tsv(output, rows)


@pytest.mark.parametrize("filt", ["cnvLength", "cnvQual", "cnvLength;cnvQual"])
def test_a_requires_original_length_quality_filters_and_preserves_evidence(filt):
    raw = record(1000, 8540, filt=filt, qual="19")
    cnv = record(1399, 8938, filt="PASS", qual="150", SVCLAIM="DJ", MatchSv="sv-full",
                 OrigCnvPos="1000", OrigCnvEnd="8540")
    selected, counts = rescue.select_rescues([raw], [cnv])
    assert counts == {"rescued_A": 1}
    e = selected[0]["evidence"]
    assert e["original"]["qual"] == "19"
    assert e["original"]["filter"] == filt
    assert e["integrated"]["qual"] == "150"
    assert e["original"]["pos"] == 1000 and e["integrated"]["pos"] == 1399
    assert e["sv_support"] == [{"original_sv_id": "sv-full"}]


@pytest.mark.parametrize("filt", ["PASS", ".", "cnvLength;dinucQual", "cnvQual;highCN",
                                   "cnvBinSupportRatio", "cnvCopyRatio"])
def test_integrated_pass_cannot_rescue_other_original_filters(filt):
    cnv = record(filt="PASS", SVCLAIM="DJ", MatchSv="sv-full")
    assert not rescue.select_rescues([record(filt=filt)], [cnv])[0]


def test_known_low_quality_geometry_uses_integrated_end_and_has_no_qual_floor():
    # Same lengths as the discussed borderline example, on synthetic coordinates.
    raw = record(1000, 8071, filt="cnvLength;cnvQual", qual="5")
    cnv = record(1000, 7748, filt="cnvLength;cnvQual", qual="5", SVCLAIM="DJ",
                 RIGHT_BND="original-sv.end", OrigCnvEnd="8071")
    sv = record(4344, 7748, identifier="sv", svtype="DEL", filt="PASS", SVCLAIM="J",
                END_RIGHT_BND_OF="cnv")
    selected, _ = rescue.select_rescues([raw], [cnv, sv])
    evidence = selected[0]["evidence"]
    assert evidence["rule"] == "B"
    assert evidence["sv_support"][0]["cnv_overlap"] == pytest.approx(3404 / 6748)
    assert evidence["sv_support"][0]["sv_overlap"] == 1
    assert evidence["original"]["end"] == 8071
    assert rescue._supported(raw, sv) is None  # original coordinates would be <50%


@pytest.mark.parametrize("mutation", ["opposite", "chrom", "sv_filter", "sample_filter",
    "ref_gt", "cnv_ref", "gt_conflict", "missing_link", "missing_forward", "too_short", "too_long", "adjacent"])
def test_b_rejects_insufficient_or_conflicting_support(mutation):
    raw, cnv, sv = paired()
    if mutation == "opposite":
        sv.fields[4] = "<DUP>"
    elif mutation == "chrom":
        sv.fields[0] = "chr2"
    elif mutation == "sv_filter":
        sv.fields[6] = "MaxDepth"
    elif mutation == "sample_filter":
        sv.fields[8] += ":FT"
        sv.fields[9] += ":MinGQ"
    elif mutation == "ref_gt":
        sv.fields[9] = "0/0:2:10,0:10,0"
    elif mutation == "cnv_ref":
        cnv.fields[9] = "0/0:2:10,0:10,0"
    elif mutation == "gt_conflict":
        sv.fields[9] = "1/1:0:0,10:0,10"
    elif mutation == "missing_link":
        sv.info.pop("END_RIGHT_BND_OF")
    elif mutation == "missing_forward":
        cnv.info.pop("RIGHT_BND")
    elif mutation == "too_short":
        sv.fields[1] = "4501"
    elif mutation == "too_long":
        sv.fields[1], sv.info["END"] = "1", "16000"
    elif mutation == "adjacent":
        sv.fields[1], sv.info["END"] = "8000", "12000"
    assert not rescue.select_rescues([raw], [cnv, sv])[0]


def test_b_includes_exact_50_percent_and_phased_genotypes():
    raw, cnv, sv = paired(gt="1|0")
    assert rescue.select_rescues([raw], [cnv, sv])[0][0]["evidence"]["rule"] == "B"


def test_b_dup_partial_missing_genotype_and_haploid_are_supported():
    raw, cnv, sv = paired()
    for r in (raw, cnv, sv):
        r.fields[4] = "<DUP>"
    cnv.fields[9] = "./1:3:30,10:20,10"
    assert rescue.select_rescues([raw], [cnv, sv])[0]
    for r in (raw, cnv, sv):
        r.fields[0] = "chrX"
        r.fields[9] = "1:0:0,10:0,10"
    assert rescue.select_rescues([raw], [cnv, sv])[0]


def test_multimatch_is_not_union_of_unrelated_sv_intervals():
    raw, cnv, sv = paired(sv_start=6500)
    other = record(1000, 2500, identifier="other", filt="PASS", svtype="DEL", SVCLAIM="J",
                   LEFT_BND_OF="cnv")
    cnv.info["LEFT_BND"] = "other-original"
    assert not rescue.select_rescues([raw], [cnv, sv, other])[0]


def test_ambiguous_original_and_nonprimary_contig_are_rejected():
    raw, cnv, sv = paired()
    assert rescue.select_rescues([raw, raw], [cnv, sv])[1]["unmatched_or_ambiguous_original"] == 1
    for r in (raw, cnv, sv):
        r.fields[0] = "chr1_alt"
    assert rescue.select_rescues([raw], [cnv, sv])[1]["non_primary_contig"] == 1


def test_vcf_requires_exact_single_sample(tmp_path):
    path = tmp_path / "source.vcf.gz"
    write_vcf(path, [record()], sample="WRONG")
    with pytest.raises(ValueError, match="sample mismatch"):
        rescue.read_vcf(path, "SRC")
    path = tmp_path / "broken.vcf"
    path.write_text("##fileformat=VCFv4.2\n")
    with pytest.raises(ValueError, match="sample header"):
        rescue.read_vcf(path, "SRC")


def setup_build(tmp_path):
    raw, cnv, sv = paired()
    orig, joint, base = tmp_path / "raw.vcf.gz", tmp_path / "joint.vcf.gz", tmp_path / "base.tsv"
    write_vcf(orig, [raw])
    write_vcf(joint, [cnv, sv])
    write_tsv(base, [baseline_row()])
    return dict(raw_cnv=orig, joint_cnv=joint, base_tsv=base, post_dir=tmp_path / "08",
                sample_id="LIS", source_sample="SRC", annotate=fake_annotate)


def test_annotation_merge_loaders_and_rerun_preserve_sources_and_stable_ids(tmp_path, monkeypatch):
    args = setup_build(tmp_path)
    before = {key: args[key].read_bytes() for key in ("raw_cnv", "joint_cnv", "base_tsv")}
    manifest = rescue.build_rescue(**args)
    assert manifest["added_events"] == 1
    stable_id = manifest["records"][0]["id"]
    review = args["post_dir"] / "LIS.cnv.review.tsv"
    _, rows = rescue._read_tsv(review)
    assert [r["AnnotSV_ID"] for r in rows] == ["existing-id", stable_id, stable_id]
    assert rows[1]["FILTER"] == "cnvLength"  # no fabricated PASS
    from app.adapters import annotsv_tsv
    monkeypatch.setattr(annotsv_tsv.gene_disease_store, "ensure_loaded", lambda: None)
    monkeypatch.setattr(annotsv_tsv.gene_disease_store, "merged_associations", lambda *a, **kw: [])
    variants, _ = annotsv_tsv.load_annotsv_tsv(review, source="cnv")
    assert variants[stable_id]["cnv_rescue"]["rule"] == "B"
    assert variants[stable_id]["acmg_class"] == 3
    assert variants[stable_id]["GT"] == "0/1"
    by_id = annotsv_tsv.load_annotsv_variants_by_ids(review, source="cnv", ids={stable_id})
    assert by_id[stable_id]["cnv_rescue"] == variants[stable_id]["cnv_rescue"]
    first = review.read_bytes()
    rescue.build_rescue(**args)
    assert review.read_bytes() == first
    assert not list(args["post_dir"].glob(".cnv-rescue-*"))
    for key, data in before.items():
        assert args[key].read_bytes() == data


@pytest.mark.parametrize("missing", [True, False])
def test_no_rescues_or_missing_input_removes_old_generation(tmp_path, missing):
    args = setup_build(tmp_path)
    rescue.build_rescue(**args)
    if missing:
        args["joint_cnv"].unlink()
    else:
        write_vcf(args["raw_cnv"], [record(filt="PASS")])
    result = rescue.build_rescue(**args)
    assert result["status"] == ("skipped" if missing else "complete")
    assert not (args["post_dir"] / "LIS.cnv.review.tsv").exists()
    assert not (args["post_dir"] / "LIS.cnv.rescued.annotated.tsv").exists()


@pytest.mark.parametrize("failure", ["tool", "missing_full", "lost_id"])
def test_annotation_failure_never_publishes_partial_review(tmp_path, failure):
    args = setup_build(tmp_path)
    def annotate(vcf, out):
        if failure == "tool":
            raise RuntimeError("tool failed")
        fake_annotate(vcf, out)
        headers, rows = rescue._read_tsv(out)
        if failure == "missing_full":
            rows = [r for r in rows if r["Annotation_mode"] != "full"]
        else:
            for r in rows:
                r["ID"] = "."
        write_tsv(out, rows, headers)
    args["annotate"] = annotate
    with pytest.raises((RuntimeError, ValueError)):
        rescue.build_rescue(**args)
    assert not (args["post_dir"] / "LIS.cnv.review.tsv").exists()
    assert not (args["post_dir"] / "LIS.cnv_rescue.json").exists()
    assert not list(args["post_dir"].glob(".cnv-rescue-*"))


def test_duplicate_event_keeps_baseline_id_and_does_not_duplicate_gene_rows(tmp_path):
    args = setup_build(tmp_path)
    row = baseline_row()
    row.update(SV_start="1000", SV_end="8000")
    write_tsv(args["base_tsv"], [row])
    result = rescue.build_rescue(**args)
    _, rows = rescue._read_tsv(args["post_dir"] / "LIS.cnv.review.tsv")
    assert result["added_events"] == 0
    assert len(result["already_annotated"]) == 1
    assert [r["AnnotSV_ID"] for r in rows] == ["existing-id"]


def test_worker_dispatch_is_dragen_only_and_passes_source_sample(tmp_path, monkeypatch):
    from app.workers import dragen_run
    calls = []
    monkeypatch.setattr(dragen_run, "_run", lambda command, **kwargs: calls.append(command))
    args = dict(sample={"sample_id": "LIS", "source_sample_id": "SRC",
                        "vcf_path": "/input/SRC.hard-filtered.vcf.gz"},
                base_tsv=tmp_path / "base.tsv", post_dir=tmp_path / "08", scripts=tmp_path)
    dragen_run._run_cnv_rescue_for_sample(mode="inhouse", **args)
    assert not calls
    dragen_run._run_cnv_rescue_for_sample(mode="dragen", **args)
    assert calls[0][-2:] == ["--source-sample", "SRC"]
    assert "--skip-cnv" not in calls[0]
    assert rescue.MANAGED_NAMES <= dragen_run.MANAGED_POSTPROCESSING_NAMES


def test_merged_parent_keeps_rescue_provenance_from_all_members():
    from app.services.cnv_sv_merge import build_parent
    common = dict(source="cnv", CHROM="1", sv_type="DEL", POS=1000, END=8000, genes=[])
    v = {"a": {**common, "id": "a", "cnv_sv_sort_score": 90},
         "b": {**common, "id": "b", "cnv_rescue": {"rule": "B"}, "cnv_sv_sort_score": 20}}
    parent = build_parent({"member_ids": ["a", "b"]}, v)
    assert parent["cnv_rescue_events"] == [{"rule": "B"}]


def test_rescue_artifacts_promote_clear_and_rollback_with_generation(tmp_path):
    from app.workers import dragen_run
    stage, live, backup = [tmp_path / name for name in ("stage", "live", "backup")]
    stage_post, live_post = stage / "08_postprocessing", live / "08_postprocessing"
    stage_post.mkdir(parents=True)
    live_post.mkdir(parents=True)
    for name in rescue.MANAGED_NAMES:
        (live_post / f"LIS.{name}").write_text("old")
    (live_post / "LIS.sample_metadata.json").write_text("reviewer state")
    (stage_post / "LIS.layout.json").write_text("new marker")
    (stage_post / "LIS.cnv_rescue.json").write_text('{"status":"skipped"}')
    operations = dragen_run._promote_staged_generation(
        sample_id="LIS", staged_sample_dir=stage, live_sample_dir=live, rollback_sample_dir=backup)
    assert not (live_post / "LIS.cnv.review.tsv").exists()
    assert not (live_post / "LIS.cnv.rescued.annotated.tsv").exists()
    assert (live_post / "LIS.sample_metadata.json").read_text() == "reviewer state"
    dragen_run._rollback_promotion_operations(operations)
    assert (live_post / "LIS.cnv.review.tsv").read_text() == "old"


def test_rescue_manifest_paths_follow_live_generation(tmp_path):
    from app.workers import dragen_run
    post = tmp_path / "stage" / "08_postprocessing"
    post.mkdir(parents=True)
    manifest = post / "LIS.cnv_rescue.json"
    manifest.write_text(json.dumps({"inputs": {"base_tsv": {
        "path": str(tmp_path / "stage/06_cnv_sv/SRC.cnv.annotated.tsv"), "sha256": "unchanged"}}}))
    dragen_run._rebase_staged_derived_paths(
        sample_id="LIS", stage_post_dir=post, final_raw_tsv=tmp_path / "live/03_acmg/SRC.tsv",
        final_post_dir=tmp_path / "live/08_postprocessing")
    baseline = json.loads(manifest.read_text())["inputs"]["base_tsv"]
    assert baseline["path"] == str(tmp_path / "live/06_cnv_sv/SRC.cnv.annotated.tsv")
    assert baseline["sha256"] == "unchanged"


def test_annotsv_container_and_override_configuration(tmp_path, monkeypatch):
    for name in ("ANNOTSV_BIN", "ANNOTSV_ANNOTATIONS", "NGS_UI_ANNOTSV_SIF"):
        monkeypatch.delenv(name, raising=False)
    annotations = tmp_path / "annotations"
    annotations.mkdir()
    monkeypatch.setenv("ANNOTSV_ANNOTATIONS", str(annotations))
    image = tmp_path / "annotsv.sif"
    image.touch()
    monkeypatch.setenv("NGS_UI_ANNOTSV_SIF", str(image))
    monkeypatch.setattr(rescue.shutil, "which", lambda name: "/bin/apptainer")
    command = rescue.annotsv_command(tmp_path / "input.vcf", tmp_path / "out.tsv", tmp_path)
    assert command[:2] == ["/bin/apptainer", "exec"]
    assert str(image) in command and "/tmp" in command
    assert command[-2:] == ["-includeCI", "0"]
    monkeypatch.setenv("ANNOTSV_BIN", "/custom/AnnotSV")
    assert rescue.annotsv_command(tmp_path / "input.vcf", tmp_path / "out.tsv", tmp_path)[0] == "/custom/AnnotSV"
