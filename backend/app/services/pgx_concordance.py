"""DRAGEN comparison and reviewer choices. Original calls are never rewritten.

The comparison uses source calls, not recommendation diplotypes. Choices bind to
both inputs and the rule version; stale choices cannot silently enter a report.
"""
from __future__ import annotations

import copy
import fcntl
import hashlib
import json
import re
import tempfile
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

from . import sample_layout, specialized

RULE_VERSION = "2"
REVIEW_NAME = "pgx_review.json"  # reviewer-owned; deliberately not pipeline-managed
GENES = (
    "ABCG2", "CACNA1S", "CFTR", "CYP2B6", "CYP2C19", "CYP2C9", "CYP2D6",
    "CYP3A5", "CYP4F2", "DPYD", "G6PD", "HLA-A", "HLA-B", "MT-RNR1",
    "NAT2", "NUDT15", "RYR1", "SLCO1B1", "TPMT", "UGT1A1", "VKORC1",
)
_STAR_REF = {"CYP2B6", "CYP2C19", "CYP2C9", "CYP2D6", "CYP3A5", "CYP4F2", "NUDT15", "TPMT", "UGT1A1"}
_NO_CALL = {"", ".", "./.", "none", "unknown", "noresult", "n/a", "na"}
# PharmVar NAT2 change log v1.7 / GeneFocus doi:10.1002/cpt.70168.
# Only the explicitly identified legacy DRAGEN definitions use these aliases.
_NAT2_LEGACY = {"*13": "*4", "*13a": "*4", "*7a": "*7", "*7b": "*7"}


def _hash(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, ensure_ascii=False).encode()).hexdigest()


def _token(gene, allele, legacy=False):
    value = re.sub(r"\s+", "", str(allele)).lower()
    value = re.sub(r"^" + re.escape(gene.lower()) + r"(?=\*)", "", value)
    # An HGVS variant's star-name annotation does not change its identity.
    if gene in {"DPYD", "RYR1"} and re.search(r"[cmng]\.", value):
        value = re.sub(r"\(\*\d+[a-z0-9.]*\)", "", value)
    if gene == "NAT2" and legacy:
        value = _NAT2_LEGACY.get(value, value)
    if value in {"reference", "ref", "wildtype", "wild-type"}:
        return "*1" if gene in _STAR_REF else "reference"
    if gene == "G6PD" and value in {"b(reference)", "b(wildtype)", "b(wild-type)"}:
        return "reference"
    if gene == "CFTR" and value == "ivacaftornon-responsivecftrsequence":
        return "reference"
    return value


def _explicit_diploid(label):
    """Expand an explicit het/hom suffix, never assume a bare variant is het."""
    match = re.fullmatch(r"(.+)\((heterozygous|het|homozygous|hom)\)", label)
    if match and "/" not in match[1]:
        variant = match[1]
        return variant + "/" + (variant if match[2] in {"homozygous", "hom"} else "reference")
    return label


def canonical(gene, label, *, legacy=False):
    """Order-insensitive within each chromosome, preserving dosage and grouping.

    Per reviewer policy + order is ignored. A bracketed single-copy combination
    remains distinct from a duplication; phase and ploidy are not flattened.
    """
    label = re.sub(r"\s+", "", str(label or "")).lower()
    label = re.sub(re.escape(gene.lower()) + r"(?=\*)", "", label)
    if gene in {"DPYD", "RYR1"}:
        label = _explicit_diploid(_token(gene, label, legacy))
    if label in _NO_CALL:
        return None
    chromosomes = []
    for part in label.split("/"):
        if not part or part in _NO_CALL:
            return None
        if part.startswith("[") and part.endswith("]"):
            chromosomes.append(("combination", tuple(sorted(_token(gene, x, legacy) for x in part[1:-1].split("+")))))
            continue
        # '+' in an HGVS coordinate is not a gene duplication separator.
        copies = part.split("+") if re.fullmatch(r"\*[^/]+", part) else [part]
        expanded = []
        for allele in copies:
            match = re.fullmatch(r"(\*\d+(?:\.\d+)?)[x×](\d+)", allele)
            if match and 1 <= int(match[2]) <= 99:
                expanded.extend([_token(gene, match[1], legacy)] * int(match[2]))
            else:
                expanded.append(_token(gene, allele, legacy))
        chromosomes.append(("copies", tuple(sorted(expanded))))
    # Phase-2 validation: G6PD B/B versus B is a representational difference.
    # Apply only to one-copy identical homozygotes; never collapse a het or CNV.
    if gene == "G6PD" and len(chromosomes) == 2 and chromosomes[0] == chromosomes[1] and len(chromosomes[0][1]) == 1:
        chromosomes = chromosomes[:1]
    return tuple(sorted(chromosomes))


def _allele_observations(gene, allele):
    value = _token(gene, allele).strip("[]")
    if value == "reference":
        return Counter()
    parts = re.split(r"[+,](?=[cmng]\.)", value)
    if not all(re.fullmatch(r"[cmng]\.[a-z0-9_*+>.-]+", part) for part in parts):
        return None
    return Counter(parts)


def _diplotype_observations(gene, label):
    label = _explicit_diploid(_token(gene, label))
    if len(label.split("/")) != 2:
        return None
    observations = Counter()
    for allele in label.split("/"):
        values = _allele_observations(gene, allele)
        if values is None:
            return None
        observations.update(values)
    return observations


def _source_observations(gene, labels, details):
    """Separate DPYD/RYR1 observations require explicit zygosity or JSON calls.

    Returning a dosage multiset does not assign the variants to chromosomes.
    recommendationDiplotypes may encode assumptions and are not a phase source.
    """
    if not labels or any("/" in label for label in labels):
        return None
    variants = details.get("comparison_variants", details.get("variants", []))
    observations = Counter()
    seen = set()
    for label in labels:
        normalized = _token(gene, label)
        explicit = _diplotype_observations(gene, normalized)
        if explicit is not None:
            observations.update(explicit)
            continue
        values = _allele_observations(gene, normalized)
        if values is None or len(values) != 1:
            return None
        name = next(iter(values))
        if name in seen:
            continue  # the position-level genotype supplies the total dosage
        matches = [v for v in variants if name in {_token(gene, a) for a in v.get("allele_names", [])}]
        if len(matches) != 1 or matches[0].get("warnings"):
            return None
        v = matches[0]
        copies = re.split(r"[/|]", v.get("call") or "")
        ref = v.get("reference")
        alt = {x for x in copies if x and x != ref and x not in {".", "?"}}
        if len(copies) != 2 or not ref or len(alt) != 1 or any(x in {"", ".", "?"} for x in copies):
            return None
        observations[name] = sum(x != ref for x in copies)
        seen.add(name)
    return observations


def _nudt15_phase_compatible(labels, candidates, details, dragen):
    """Bounded, evidence-backed *1/*3 vs legacy *1/*2;*3/*6 ambiguity.

    CPIC/ClinPGx: *2 is now a *3 suballele. Different phase solutions can have
    different phenotypes; compatibility never rewrites a genotype/phenotype.
    """
    if details.get("definition_source") != "CLINPGX" or details.get("definition_version") != "2026-02-09-10-28":
        return False
    if dragen.get("starAlleleDatabaseSource") != "PharmGKB Database Version: Snapshot-2022.01.01":
        return False
    if details.get("phased") or details.get("effectively_phased"):
        return False
    ours = {canonical("NUDT15", label) for label in labels}
    expected_ours = {canonical("NUDT15", x) for x in ("*1/*3", "*3/*6")}
    theirs = {canonical("NUDT15", x) for x in candidates}
    expected_theirs = {canonical("NUDT15", x) for x in ("*1/*2", "*3/*6")}
    if not ours or not ours <= expected_ours or theirs != expected_theirs:
        return False
    variants = details.get("comparison_variants", [])
    for rsid, position, ref, alleles in (
        ("rs746071566", 48037782, "AGGAGTC", {"AGGAGTC", "AGGAGTCGGAGTC"}),
        ("rs116855232", 48045719, "C", {"C", "T"}),
    ):
        rows = [v for v in variants if v.get("rsid") == rsid]
        if len(rows) != 1:
            return False
        v = rows[0]
        if (str(v.get("chr", "")).removeprefix("chr") != "13" or v.get("pos") != position
                or v.get("reference") != ref or v.get("warnings") or v.get("phased") or v.get("phase_set")
                or "|" in v.get("call", "") or set(v.get("call", "").split("/")) != alleles
                or len(v.get("call", "").split("/")) != 2):
            return False
    return True


def _labels(details, fallback=""):
    labels = [x.get("label") for x in details.get("source_diplotypes", []) if x.get("label")]
    return labels or ([details.get("label") or fallback] if details.get("label") or fallback else [])


def _schema(label):
    value = str(label).lower()
    if re.search(r"rs\d+", value):
        return "rs"
    if re.search(r"[cmng]\.\d|[ACGT]>[ACGT]", value, re.I):
        return "hgvs"
    if "*" in value:
        return "star"
    return "named"


def compare_gene(gene, payload, dragen):
    details = payload.get("details") or {}
    labels = _labels(details)
    if gene == "MT-RNR1" and not any(canonical(gene, x) for x in labels):
        labels = [payload.get("diplotype") or ""]
    raw = str(dragen.get("genotype") or "")
    candidates = [s.strip() for s in raw.split(";") if canonical(gene, s)]
    legacy = gene == "NAT2" and dragen.get("starAlleleDatabaseSource") == "PharmGKB Database Version: Snapshot-2022.01.01"
    # Older reports cannot be assumed to use the post-2024 NAT2 definitions.
    if legacy and str(details.get("definition_version") or "") < "2024-03":
        legacy = False
    ours = {v for label in labels if (v := canonical(gene, label)) is not None}
    theirs = {canonical(gene, x, legacy=legacy) for x in candidates}
    status = "unavailable"
    observation_phase = False
    if ours and theirs:
        component_sources = gene in {"DPYD", "RYR1"} and all("/" not in label for label in labels)
        if component_sources:
            observed = _source_observations(gene, labels, details)
            target_sets = [_diplotype_observations(gene, label) for label in candidates]
            if observed is not None and target_sets and all(x is not None for x in target_sets):
                if observed in target_sets:
                    if len(observed) > 1:
                        if details.get("effectively_phased") or details.get("phased"):
                            # Position dosage alone cannot verify a known phase
                            # that the component-style labels fail to preserve.
                            status = "uncomparable"
                        else:
                            observation_phase = True
                            status = "phase_compatible"
                    else:
                        status = "same"
                else:
                    status = "different"
            elif ours & theirs:
                status = "same"
            else:
                status = "uncomparable"
        elif ours & theirs:
            status = "same"
        elif gene == "NUDT15" and _nudt15_phase_compatible(labels, candidates, details, dragen):
            status = "phase_compatible"
            observation_phase = True
        elif {_schema(s) for s in labels} != {_schema(s) for s in candidates}:
            status = "uncomparable"
        else:
            status = "different"
    low_gq = "lowgq" in json.dumps([dragen.get("missingVariantSites") or [], dragen.get("genotypeFilter") or ""]).lower()
    phase = observation_phase or (status == "same" and (
        len(theirs) > 1 or (len(ours) > 1 and not details.get("phased") and not details.get("effectively_phased"))
    ))
    badges = (["結果不一致"] if status == "different" else ["無法比對"] if status == "uncomparable" else [])
    if low_gq:
        badges.append("LowGQ")
    if phase:
        badges.append("相位未定")
    return {"status": status, "badges": badges, "platform_genotypes": labels,
            "platform_phenotypes": details.get("phenotypes") or [],
            "dragen_genotype": raw, "candidates": candidates,
            "dragen_phenotype": dragen.get("phenotypeDatabaseAnnotation") or "",
            "missing_sites": dragen.get("missingVariantSites") or [],
            "genotype_quality": dragen.get("genotypeQuality"), "genotype_filter": dragen.get("genotypeFilter") or ""}


def _read_json(path):
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError("PGx 檔案格式錯誤")
    return value


def _targeted(sample_id):
    """Prefer the validated local snapshot. Never bypass a damaged bundle."""
    source = sample_layout._pipeline_source(sample_id)
    source_id = sample_layout.source_sample_id(sample_id)
    manifest_path = sample_layout.state_file(sample_id, specialized.NAMES["manifest"])
    if manifest_path.is_file():
        specialized.load(sample_id)  # validates all checksums and generation
        manifest = _read_json(manifest_path)
        if manifest.get("source_sample_id") != source_id:
            raise ValueError("DRAGEN 快照來源樣本不符")
        item = manifest["files"]["source"]
        content = (manifest_path.parent / item["name"]).read_bytes()
        if hashlib.sha256(content).hexdigest() != item["sha256"]:
            raise ValueError("DRAGEN 快照校驗失敗")
        data = json.loads(content)
        if not data:  # an explicit not_available bundle is not an alternate run
            return {}, ""
    else:
        if source.get("pipeline_type") != "dragen" or not source.get("source_vcf_path"):
            return {}, ""
        from .gba import discover_json
        path = discover_json(Path(source["source_vcf_path"]), source_id)
        if path is None:
            return {}, ""
        data = _read_json(path)
    if not isinstance(data, dict) or data.get("sampleId") != source_id:
        raise ValueError("DRAGEN targeted.json 樣本不符")
    if str(data.get("genomeBuild")) not in {"hg38", "GRCh38", "hg19", "GRCh37"}:
        raise ValueError("DRAGEN genome build 無法確認")
    return data, _hash(data)


def _review(sample_id):
    path = sample_layout.state_file(sample_id, REVIEW_NAME)
    data = _read_json(path) if path.is_file() else {"choices": {}, "audit": []}
    if not isinstance(data.get("choices"), dict) or not isinstance(data.get("audit"), list) or any(not isinstance(c, dict) for c in data["choices"].values()):
        raise ValueError("PGx 人工選擇紀錄格式錯誤")
    if any(gene not in GENES or choice.get("source") != "dragen" or not isinstance(choice.get("genotype"), str) or not choice.get("fingerprint") for gene, choice in data["choices"].items()):
        raise ValueError("PGx 人工選擇紀錄格式錯誤")
    return data


def attach(sample_id, pgx):
    """Attach source comparison and apply only source-bound reviewer choices."""
    pgx["concordance"] = {"genes": {}, "message": "", "rule_version": RULE_VERSION}
    pgx["review_blocked"] = False
    pgx["recommendations_pending"] = []
    choices = {}
    review_exists = sample_layout.state_file(sample_id, REVIEW_NAME).is_file()
    review_valid = False
    try:
        review = _review(sample_id)
        review_valid = True
        choices = review.get("choices") or {}
        data, digest = _targeted(sample_id)
        if not data:
            pgx["concordance"]["message"] = "無 DRAGEN 比對資料"
        identity = pgx.get("source_identity") or {}
        source_id = sample_layout.source_sample_id(sample_id)
        if data and identity.get("sample_id") and identity["sample_id"] != source_id:
            raise ValueError("PharmCAT 來源樣本不符")
        build = str(identity.get("genome_build") or "")
        target_build = str(data.get("genomeBuild") or "")
        build_group = lambda x: "38" if x in {"hg38", "GRCh38"} or x.startswith("GRCh38.") else "37" if x in {"hg19", "GRCh37"} or x.startswith("GRCh37.") else ""
        if data and build and (not build_group(build) or build_group(build) != build_group(target_build)):
            raise ValueError("PharmCAT／DRAGEN genome build 不符")
        mapping = {x["gene"]: x for x in data.get("locusAnnotations", []) if isinstance(x, dict) and x.get("gene")}
        for key in ("cyp2d6", "cyp2b6"):
            if isinstance(data.get(key), dict):
                mapping[key.upper()] = data[key]
        blocked = []
        platform_genes = copy.deepcopy(pgx.get("genes") or {})
        for gene in GENES:
            payload = platform_genes.get(gene) or {}
            dr = mapping.get(gene) or {}
            row = compare_gene(gene, payload, dr)
            fingerprint = _hash([RULE_VERSION, source_id, sample_layout._pipeline_source(sample_id),
                                 identity, pgx.get("data_version"), pgx.get("pharmcat_version"),
                                 payload, digest])
            row.update(fingerprint=fingerprint, selected_source="platform", selected_genotype="",
                       reviewer="", reviewed_at="", stale=False)
            choice = choices.get(gene) or {}
            if choice:
                row.update(reviewer=choice.get("reviewer_username", ""), reviewed_at=choice.get("updated_at", ""))
                if choice.get("fingerprint") != fingerprint or (choice.get("source") == "dragen" and choice.get("genotype") not in row["candidates"]):
                    row["stale"] = True
                    row["badges"].append("來源已更新，請重新確認")
                    blocked.append(gene)
                elif choice.get("source") == "dragen":
                    row.update(selected_source="dragen", selected_genotype=choice["genotype"])
                    # An aggregate phenotype for several candidates cannot be
                    # attributed to the one candidate selected by a reviewer.
                    phenotype = row["dragen_phenotype"] if len(row["candidates"]) == 1 else ""
                    _apply(pgx, gene, choice["genotype"], phenotype)
            pgx["concordance"]["genes"][gene] = row
        pgx["review_blocked"] = bool(blocked)
        pgx["review_fingerprint"] = _hash([digest, identity, choices, [r["fingerprint"] for r in pgx["concordance"]["genes"].values()]])
    except (OSError, ValueError, KeyError, TypeError) as exc:
        pgx["concordance"]["message"] = str(exc)
        # A lost source must not silently undo a previously selected call.
        pgx["review_blocked"] = bool(choices) or (review_exists and not review_valid)
        pgx["review_fingerprint"] = _hash(["error", str(exc)])
        for gene in choices:
            if gene in GENES:
                row = compare_gene(gene, pgx.get("genes", {}).get(gene) or {}, {})
                row.update(stale=True, selected_source="platform", fingerprint="",
                           badges=["來源已更新，請重新確認"])
                pgx["concordance"]["genes"][gene] = row
    return pgx


def _apply(pgx, gene, genotype, phenotype):
    payload = pgx["genes"][gene]
    payload["original_details"] = copy.deepcopy(payload.get("details") or {})
    alleles = genotype.split("/")
    details = {"label": genotype, "phenotypes": [phenotype] if phenotype else [],
               "allele1_name": alleles[0], "allele2_name": alleles[1] if len(alleles) == 2 else "",
               "allele1_function": "", "allele2_function": "", "source_diplotypes": [],
               "variants": [], "messages": [], "uncalled": []}
    details["source_diplotypes"] = [{"label": genotype, "allele1_name": details["allele1_name"],
                                     "allele2_name": details["allele2_name"], "phenotypes": details["phenotypes"]}]
    payload.update(details=details, diplotype=genotype, phenotype=phenotype, activity_score="", drugs=[], actionable=False)
    # These recommendations were computed for the original genotype. Suppress
    # the whole joint annotation, including other genes, rather than relabel it.
    removed = [x for x in pgx.get("guideline_annotations", []) if gene in (x.get("genes") or [])]
    pending = set(pgx["recommendations_pending"]) | {gene}
    for annotation in removed:
        pending.update(g for g in annotation.get("genes", []) if g in GENES)
    pgx["recommendations_pending"] = sorted(pending)
    pgx["guideline_annotations"] = [x for x in pgx.get("guideline_annotations", []) if x not in removed]
    # Legacy PDF consumers use TSV drug lists rather than the presenter.
    for affected in pending:
        if affected in pgx["genes"]:
            pgx["genes"][affected]["drugs"] = []


def save_choice(sample_id, gene, source, genotype, fingerprint, user):
    """The browser submits a choice, never an arbitrary genotype or phenotype."""
    if gene not in GENES or source not in {"platform", "dragen"}:
        raise ValueError("無效的 PGx 選擇")
    if not re.fullmatch(r"[A-Za-z0-9_][A-Za-z0-9_.-]*", sample_id):
        raise ValueError("無效的樣本 ID")
    source_id = sample_layout.source_sample_id(sample_id)
    if not re.fullmatch(r"[A-Za-z0-9_][A-Za-z0-9_.-]*", source_id):
        raise ValueError("無效的樣本 ID")
    from ..workers.dragen_run import _acquire_sample_locks, _release_sample_locks
    from .sample_loader import load_sample_pgx
    handles = _acquire_sample_locks([sample_id, source_id])
    try:
        path = sample_layout.state_file(sample_id, REVIEW_NAME, for_write=True)
        if not sample_layout.state_file(sample_id, "sample_metadata.json").is_file():
            raise ValueError("個案尚未登錄")
        lock_path = path.with_suffix(".lock")
        with lock_path.open("a+") as lock:
            fcntl.flock(lock.fileno(), fcntl.LOCK_EX)
            current = load_sample_pgx(sample_id)["pgx"]
            row = current["concordance"]["genes"].get(gene)
            if source == "dragen":
                if not row or row["status"] != "different" or row["fingerprint"] != fingerprint or genotype not in row["candidates"]:
                    raise ValueError("比對來源已變更，請重新載入後確認")
            review = _review(sample_id)
            review.setdefault("choices", {})
            old = review["choices"].get(gene)
            choice = {"source": source, "genotype": genotype if source == "dragen" else "",
                      "fingerprint": row["fingerprint"] if row else "", "reviewer_user_id": user.get("id"),
                      "reviewer_username": user.get("username") or "",
                      "updated_at": datetime.now(timezone.utc).isoformat()}
            choice["evidence"] = {key: row.get(key) for key in (
                "platform_genotypes", "platform_phenotypes", "dragen_genotype", "dragen_phenotype",
                "missing_sites", "genotype_quality", "genotype_filter",
            )} if row else {}
            choice.update(source_sample_id=source_id, rule_version=RULE_VERSION)
            review.setdefault("audit", []).append({"gene": gene, "previous": old, **choice})
            if source == "platform":
                review["choices"].pop(gene, None)
            else:
                review["choices"][gene] = choice
            path.parent.mkdir(parents=True, exist_ok=True)
            with tempfile.NamedTemporaryFile(mode="w", dir=path.parent, prefix=".pgx-review-", delete=False) as tmp:
                tmp_path = Path(tmp.name)
                try:
                    json.dump(review, tmp, ensure_ascii=False, indent=2)
                    tmp.flush()
                    tmp_path.replace(path)
                finally:
                    tmp_path.unlink(missing_ok=True)
            return load_sample_pgx(sample_id)
    finally:
        _release_sample_locks(handles)


def assert_exportable(pgx):
    if pgx.get("review_blocked"):
        raise ValueError("PGx 比對來源或人工選擇需要重新確認，請重新載入個案")
