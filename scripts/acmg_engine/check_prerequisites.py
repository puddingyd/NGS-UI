#!/usr/bin/env python3
"""Read-only prerequisite check for the ACMG rule engine (plan section 5.2).

Run on the DGM dev machine; nothing is written anywhere:

    cd ~/NGS_UI/NGS-UI && git fetch -q origin claude/plan-ngs-ui-RQW8J && \
      git show FETCH_HEAD:scripts/acmg_engine/check_prerequisites.py | python3 -

Options: --home DIR (default $NGS_UI_HOME or ~/NGS_UI), --cases N (recent cases
per test type whose review TSV header is inspected, default 10), --no-network.

Reports only paths, sizes, dates, column names and counts. Sample IDs of the
inspected recent cases are printed; patient names / MRNs are never read.
Only stdlib; works on Python 3.8+.
"""
from __future__ import annotations

import argparse
import collections
import gzip
import hashlib
import json
import os
import re
import subprocess
import sys
import time
import urllib.request

REVIEWED_CLASSIFIER_SHA256 = "ccaa944736c6b7f48af552fa45125a424fa8e3d5cfc57afa46ec7a238f014b2b"
CLASSIFIER_MARKERS = {
    "HI=30(AR_only)": "PVS1 只給 ClinGen HI=3（隱性基因不給）",
    "def min_eas_af": "PM2 取東亞來源的最小值",
    "cascade": "PP3／BP4 依序換工具",
    "PVS1_REASON": "v3.6 的 PVS1 決策樹輸出",
    "CLINGEN_AGREEMENT": "v3.6 的 ClinGen 專家判讀對照",
    "HAPLOID_HET": "v3.8 欄位",
    "GNOMAD41": "使用 gnomAD v4.1",
    "INHOUSE": "使用院內頻率",
}
REVIEW_COLUMNS = [
    "GIAB_STRATA", "INHOUSE_AF", "INHOUSE_AN", "SPLICEAI_MAX", "PANGOLIN_SCORE",
    "PKNN_EVIDENCE", "BAYESDEL_NOAF", "NMD", "PROTEIN_POSITION", "HGVS_C", "HGVS_P",
    "TRANSCRIPT_TYPE", "LOFTEE", "LOFTEE_FILTER", "GNOMAD41_JOINT_AF",
    "GNOMAD41_JOINT_EAS_AF", "CLINVAR_SIG", "CLINVAR_STARS", "CLINGEN_VCEP_CLASS",
    "CLINGEN_VCEP_CRITERIA", "PHYLOP100", "HAPLOID_HET", "STRAND_BIAS",
]
PATH_KEY_RE = re.compile(r"^(NGS_UI_[A-Z0-9_]*(DIR|DB|ROOT|HOME|SUMMARY|VCF|CONFIG|PATH)|TERTIARY_OUTPUT_ROOT)$")
CODE_DIRS = [
    "/home/pipeline/tertiary_code",
    "/home/datalake_Intermediate/pipeline/tertiary_code",
    "/home/datalake_Intermediate/pipeline/pipeline_code",
]
REF_DIRS = ["/home/pipeline/reference/hg38/tertiary", "/home/pipeline/reference/hg38"]
URLS = [
    ("ClinVar variant_summary", "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"),
    ("ClinVar VCF 歷史版本目錄", "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/archive_2.0/"),
    ("ClinGen dosage", "https://ftp.clinicalgenome.org/ClinGen_gene_curation_list_GRCh38.tsv"),
    ("MANE 目前版本目錄", "https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/"),
]


def head(title: str) -> None:
    print(f"\n== {title}")


def show(label: str, value: object) -> None:
    print(f"   {label}：{value}")


def fmt_size(n: int) -> str:
    for unit in ("B", "KB", "MB", "GB"):
        if n < 1024 or unit == "GB":
            return f"{n:.0f} {unit}" if unit == "B" else f"{n:.1f} {unit}"
        n /= 1024.0
    return str(n)


def describe_file(path: str) -> str:
    try:
        st = os.stat(path)
    except OSError:
        return "不存在"
    stamp = time.strftime("%Y-%m-%d", time.localtime(st.st_mtime))
    return f"存在，{fmt_size(st.st_size)}，修改日期 {stamp}"


def open_text(path: str):
    if path.endswith((".gz", ".bgz")):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, encoding="utf-8", errors="replace")


def first_lines(path: str, n: int = 3) -> list[str]:
    out = []
    try:
        with open_text(path) as handle:
            for line in handle:
                out.append(line.rstrip("\n")[:200])
                if len(out) >= n:
                    break
    except OSError as exc:
        out.append(f"(無法讀取：{exc})")
    return out


def header_columns(path: str) -> list:
    """Full first line split on tabs (never truncated)."""
    try:
        with open_text(path) as handle:
            return handle.readline().rstrip("\n").split("\t")
    except OSError:
        return []


def walk_limited(root: str, depth: int):
    root = root.rstrip("/")
    base = root.count(os.sep)
    for dirpath, dirnames, filenames in os.walk(root):
        if dirpath.count(os.sep) - base >= depth:
            dirnames[:] = []
        dirnames[:] = [d for d in dirnames if not d.startswith(".")]
        yield dirpath, filenames


def service_env() -> dict:
    """Path-like NGS_UI_* settings of the running ngs-ui service (never secrets)."""
    env: dict = {}
    try:
        raw = subprocess.run(["systemctl", "show", "ngs-ui", "-p", "Environment", "--value"],
                             capture_output=True, text=True, timeout=10).stdout
        for token in raw.split():
            key, _, value = token.partition("=")
            if PATH_KEY_RE.match(key):
                env[key] = value
        files = subprocess.run(["systemctl", "show", "ngs-ui", "-p", "EnvironmentFiles", "--value"],
                               capture_output=True, text=True, timeout=10).stdout
        for token in files.split():
            path = token.lstrip("-")
            if not os.path.isfile(path):
                continue
            try:
                with open(path, encoding="utf-8", errors="replace") as handle:
                    for line in handle:
                        key, _, value = line.strip().partition("=")
                        key = key.replace("export ", "").strip()
                        if PATH_KEY_RE.match(key):
                            env[key] = value.strip().strip('"').strip("'")
            except OSError:
                pass
    except (OSError, subprocess.SubprocessError):
        pass
    return env


def resolve(env: dict, key: str, default: str) -> str:
    return os.path.expanduser(env.get(key) or os.environ.get(key) or default)


def sha256(path: str) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def check_environment(args) -> dict:
    head("1. 執行環境與平台設定")
    show("Python", sys.version.split()[0])
    env = service_env()
    if env:
        show("ngs-ui 服務的路徑設定", "")
        for key in sorted(env):
            print(f"      {key}={env[key]}")
    else:
        show("ngs-ui 服務的路徑設定", "沒有讀到（使用預設值）")
    home = os.path.expanduser(args.home or env.get("NGS_UI_HOME") or os.environ.get("NGS_UI_HOME") or "~/NGS_UI")
    biotools = resolve(env, "NGS_UI_BIOTOOLS_DIR", os.path.join(home, "biotools"))
    paths = {
        "home": home,
        "biotools": biotools,
        "tertiary_root": resolve(env, "NGS_UI_TERTIARY_ROOT",
                                 env.get("NGS_UI_PIPELINE_OUT_ROOT") or "/home/datalake_Intermediate/pipeline/tertiary_output"),
        "clinvar_dir": resolve(env, "NGS_UI_CLINVAR_LATEST_DIR", os.path.join(biotools, "clinvar_latest")),
        "mane_summary": resolve(env, "NGS_UI_MANE_SUMMARY", os.path.join(biotools, "MANE.GRCh38.v1.5.summary.txt.gz")),
        "inhouse_af": resolve(env, "NGS_UI_INHOUSE_AF_DB", os.path.join(biotools, "inhouse_af", "inhouse_af.hg38.vcf.gz")),
        "giab": resolve(env, "NGS_UI_GIAB_STRAT_DIR", os.path.join(biotools, "giab_stratification")),
        "genebe": resolve(env, "NGS_UI_GENEBE_DB", os.path.join(biotools, "genebe", "genebe_hg38.tsv.gz")),
        "tertiary_config": resolve(env, "NGS_UI_TERTIARY_CONFIG", "/home/pipeline/tertiary_code/nextflow_tertiary.config"),
    }
    show("NGS_UI_HOME", home)
    show("三級輸出 root", f"{paths['tertiary_root']}（{'存在' if os.path.isdir(paths['tertiary_root']) else '不存在'}）")
    return paths


def check_classifier(paths: dict) -> None:
    head("2. pipeline 的 ACMG 分類程式")
    dirs = list(dict.fromkeys([os.path.dirname(paths["tertiary_config"])] + CODE_DIRS))
    found = []
    for d in dirs:
        if not os.path.isdir(d):
            continue
        for dirpath, filenames in walk_limited(d, 3):
            if "acmg_classifier.py" in filenames:
                found.append(os.path.join(dirpath, "acmg_classifier.py"))
    if not found:
        show("acmg_classifier.py", f"在 {', '.join(dirs)} 下找不到")
        return
    for path in dict.fromkeys(found):
        print(f"   {path}：{describe_file(path)}")
        digest = sha256(path)
        same = "是" if digest == REVIEWED_CLASSIFIER_SHA256 else "否"
        print(f"      sha256 {digest[:16]}…（與審查的版本相同：{same}）")
        with open(path, encoding="utf-8", errors="replace") as handle:
            text = handle.read()
        desc = re.search(r'description="([^"]+)"', text)
        print(f"      命令列說明：{desc.group(1) if desc else '（未找到）'}")
        for marker, meaning in CLASSIFIER_MARKERS.items():
            print(f"      {'有' if marker in text else '無'}  {marker:<20} {meaning}")


def check_pipeline_refs(paths: dict) -> None:
    head("3. pipeline 用的 ClinGen HI／遺傳模式表")
    cfg = paths["tertiary_config"]
    show("Nextflow config", f"{cfg}（{describe_file(cfg)}）")
    if os.path.isfile(cfg):
        with open(cfg, encoding="utf-8", errors="replace") as handle:
            for line in handle:
                low = line.lower()
                if any(k in low for k in ("clingen", "gene_moi", "dosage", "ref_dir", "reference")):
                    print(f"      {line.strip()[:180]}")
    hits = []
    for d in REF_DIRS:
        if not os.path.isdir(d):
            continue
        for dirpath, filenames in walk_limited(d, 3):
            for name in filenames:
                low = name.lower()
                if any(k in low for k in ("clingen", "dosage", "gene_moi", "mane")) and not low.endswith((".tbi", ".csi")):
                    hits.append(os.path.join(dirpath, name))
    for path in sorted(set(hits))[:20]:
        print(f"   {path}：{describe_file(path)}")
        if path.endswith((".tsv", ".tsv.gz", ".txt", ".txt.gz")):
            for line in first_lines(path, 2):
                print(f"      | {line}")


def check_resources(paths: dict) -> None:
    head("4. 平台現有的參考資料")
    cv = paths["clinvar_dir"]
    show("每週 ClinVar VCF", describe_file(os.path.join(cv, "clinvar.vcf.gz")))
    show("每週 ClinVar SQLite", describe_file(os.path.join(cv, "clinvar_latest.sqlite")))
    manifest = os.path.join(cv, "clinvar_latest_manifest.json")
    if os.path.isfile(manifest):
        try:
            data = json.load(open(manifest, encoding="utf-8"))
            show("ClinVar 版本日期", data.get("release_date") or data.get("release") or "（manifest 沒有日期欄位）")
        except (OSError, ValueError):
            show("ClinVar manifest", "無法解析")
    show("MANE summary", f"{paths['mane_summary']}：{describe_file(paths['mane_summary'])}")
    gtfs = []
    for d in [paths["biotools"]] + REF_DIRS:
        if os.path.isdir(d):
            for dirpath, filenames in walk_limited(d, 3):
                gtfs += [os.path.join(dirpath, f) for f in filenames if "mane" in f.lower() and ".gtf" in f.lower()]
    show("MANE GTF", "、".join(sorted(set(gtfs))) if gtfs else "找不到（步驟 2 需要下載）")
    show("GeneBe DB", describe_file(paths["genebe"]))
    giab = paths["giab"]
    manifest = os.path.join(giab, "strata_manifest.json")
    if os.path.isfile(manifest):
        try:
            labels = [item.get("label") for item in json.load(open(manifest, encoding="utf-8"))]
            show("GIAB strata", "、".join(str(x) for x in labels))
        except (OSError, ValueError, AttributeError):
            show("GIAB strata", "manifest 無法解析")
    else:
        show("GIAB strata", f"{giab}：沒有 strata_manifest.json")
    for name in ("snv", "indel"):
        p = os.path.join(paths["biotools"], "spliceai", f"spliceai_scores.raw.{name}.hg38.vcf.gz")
        show(f"SpliceAI {name}", describe_file(p))
    check_inhouse_af(paths["inhouse_af"])


def check_inhouse_af(path: str) -> None:
    show("院內 AF（nckuh-af）", f"{path}：{describe_file(path)}")
    if not os.path.isfile(path):
        return
    info_ids, meta_keys, an_values = [], [], []
    try:
        with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                if line.startswith("##INFO=<ID="):
                    info_ids.append(line[11:].split(",", 1)[0])
                elif line.startswith("##"):
                    key, _, value = line[2:].strip().partition("=")
                    if key not in ("contig", "FILTER", "FORMAT", "ALT") and not key.startswith("bcftools"):
                        numeric = value if value.replace(".", "", 1).isdigit() else ""
                        meta_keys.append(f"{key}={numeric}" if numeric else key)
                elif line.startswith("#"):
                    continue
                else:
                    fields = line.split("\t", 8)
                    if len(fields) > 7:
                        for item in fields[7].split(";"):
                            k, _, v = item.partition("=")
                            if k in ("AN", "INHOUSE_AN"):
                                try:
                                    an_values.append(int(v.split(",")[0]))
                                except ValueError:
                                    pass
                    if len(an_values) >= 5000:
                        break
    except (OSError, EOFError) as exc:
        show("院內 AF 讀取", f"失敗：{exc}")
        return
    print(f"      INFO 欄位：{', '.join(info_ids)}")
    print(f"      其他 header：{', '.join(sorted(set(meta_keys)))[:300]}")
    if an_values:
        an_values.sort()
        print(f"      前 {len(an_values)} 筆的 AN：最小 {an_values[0]}、中位 {an_values[len(an_values)//2]}、最大 {an_values[-1]}"
              f"（最大值約等於 2 × 院內人數）")


def scoped(directory: str, sid: str, name: str) -> str:
    prefixed = os.path.join(directory, f"{sid}.{name}")
    return prefixed if os.path.exists(prefixed) else os.path.join(directory, name)


def load_json(path: str) -> dict:
    try:
        with open(path, encoding="utf-8") as handle:
            data = json.load(handle)
        return data if isinstance(data, dict) else {}
    except (OSError, ValueError):
        return {}


def test_type_of(sid: str, post: str) -> str:
    manifest = load_json(scoped(post, sid, "snv_indel.review.tsv.source.json"))
    value = str(manifest.get("test_type") or "").upper()
    if not value:
        value = str(load_json(scoped(post, sid, "pipeline_source.json")).get("seq_type") or "").upper()
    if not value and sid.endswith("-dragen"):
        value = "WGS"
    if re.match(r"^\d{2}T", sid):
        value = "TITAN-WGS"
    return value or "未知"


def check_cases(paths: dict, n_recent: int) -> None:
    head("5. 三級分析個案")
    root = paths["tertiary_root"]
    if not os.path.isdir(root):
        show("三級輸出 root", "不存在，略過")
        return
    by_type = collections.Counter()
    with_review = collections.Counter()
    raw_cols = collections.defaultdict(collections.Counter)
    recent = collections.defaultdict(list)
    wes_patients: set = set()
    wes_test = 0
    for sid in sorted(os.listdir(root)):
        sample_dir = os.path.join(root, sid)
        if sid.startswith((".", "_")) or not os.path.isdir(sample_dir):
            continue
        post = os.path.join(sample_dir, "08_postprocessing")
        ttype = test_type_of(sid, post)
        group = "WES" if ttype == "WES" else ("WGS" if "WGS" in ttype else ttype)
        by_type[group] += 1
        acmg_dir = os.path.join(sample_dir, "03_acmg")
        raws = sorted(f for f in os.listdir(acmg_dir) if f.endswith(".snv_indel.acmg.tsv")) if os.path.isdir(acmg_dir) else []
        if raws:
            raw_cols[group][len(header_columns(os.path.join(acmg_dir, raws[0])))] += 1
        else:
            raw_value = str(load_json(scoped(post, sid, "layout.json")).get("raw_tsv")
                            or load_json(scoped(post, sid, "pipeline_source.json")).get("source_path") or "")
            if raw_value and not os.path.isabs(raw_value):
                raw_value = os.path.join(sample_dir, raw_value)
            if raw_value and os.path.isfile(raw_value):
                raw_cols[group][f"沒有 03_acmg（原始檔在別處，{len(header_columns(raw_value))} 欄）"] += 1
            else:
                raw_cols[group]["沒有 03_acmg（找不到原始檔）"] += 1
        if group == "WES":
            base = re.sub(r"-(nckuh|dragen|inhouse)$", "", sid)
            if re.search(r"\d+T$", base):
                wes_test += 1
            else:
                wes_patients.add(re.sub(r"_\d+$", "", base))
        review = scoped(post, sid, "snv_indel.review.tsv")
        if os.path.isfile(review):
            with_review[group] += 1
            recent[group].append((os.path.getmtime(review), sid, review))
    for group in sorted(by_type):
        cols = "、".join(f"{k} 欄 {v} 例" if isinstance(k, int) else f"{k} {v} 例" for k, v in sorted(raw_cols[group].items(), key=str))
        show(group, f"{by_type[group]} 例；有 review TSV {with_review[group]} 例；03_acmg 欄位數：{cols}")
    if by_type.get("WES"):
        show("WES 出現次數表的基數", f"排除 T（test）{wes_test} 例、合併 _1～_N 後約 {len(wes_patients)} 位病人（尚未用病歷號合併）")
    for group in sorted(recent):
        picks = sorted(recent[group], reverse=True)[:n_recent]
        if not picks:
            continue
        print(f"\n   {group} 最近 {len(picks)} 個個案的 review TSV 欄位（有幾例含該欄）：")
        print(f"      個案：{', '.join(sid for _, sid, _ in picks)}")
        counts = collections.Counter()
        for _, _, review in picks:
            header = set(header_columns(review))
            for col in REVIEW_COLUMNS:
                counts[col] += col in header
        for col in REVIEW_COLUMNS:
            print(f"      {counts[col]:>3}/{len(picks)}  {col}")


def check_network(enabled: bool) -> None:
    head("6. 下載來源是否連得到")
    if not enabled:
        show("略過", "--no-network")
        return
    for label, url in URLS:
        req = urllib.request.Request(url, method="HEAD" if not url.endswith("/") else "GET",
                                     headers={"User-Agent": "ngs-ui-prereq-check"})
        try:
            with urllib.request.urlopen(req, timeout=20) as resp:
                size = resp.headers.get("Content-Length")
                extra = f"，{fmt_size(int(size))}" if size and size.isdigit() else ""
                body = resp.read(200000).decode("utf-8", "replace") if url.endswith("/") else ""
            gtf = sorted(set(re.findall(r"MANE\.GRCh38\.v[\d.]+\.ensembl_genomic\.gtf\.gz", body)))
            note = f"（{', '.join(gtf)}）" if gtf else ""
            show(label, f"可連線 HTTP {resp.status}{extra}{note}")
        except Exception as exc:  # noqa: BLE001 - report every failure kind
            show(label, f"連不到：{exc}")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="ACMG rule engine prerequisite check (read-only)")
    ap.add_argument("--home", help="NGS_UI_HOME (default $NGS_UI_HOME or ~/NGS_UI)")
    ap.add_argument("--cases", type=int, default=10, help="recent cases per test type to inspect")
    ap.add_argument("--no-network", action="store_true", help="skip download reachability checks")
    ap.add_argument("--only-cases", action="store_true", help="run only the case section (section 5)")
    args = ap.parse_args(argv)
    print("ACMG 規則引擎前提檢查（只讀，不寫入任何檔案）")
    paths = check_environment(args)
    if args.only_cases:
        check_cases(paths, max(1, args.cases))
        print("\n完成。請把以上輸出整段貼回。")
        return 0
    check_classifier(paths)
    check_pipeline_refs(paths)
    check_resources(paths)
    check_cases(paths, max(1, args.cases))
    check_network(not args.no_network)
    print("\n完成。請把以上輸出整段貼回。")
    return 0


if __name__ == "__main__":
    sys.exit(main())
