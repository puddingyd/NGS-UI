"""Parse copied DRAGEN/NCKUH ploidy VCF sidecars for review UI use."""
from __future__ import annotations

import gzip
import re
from pathlib import Path

from . import sample_layout


_HEADER_FIELDS = {
    "estimatedSexKaryotype": "karyotype",
    "referenceSexKaryotype": "reference_karyotype",
    "source": "pipeline_source",
    "seqType": "seq_type",
    "autosomeDepthOfCoverage": "autosome_depth",
}


def _clean_karyotype(value: str) -> str:
    return re.sub(r"[^A-Z0-9+-]", "", str(value or "").upper())


def _number(value: str):
    text = str(value or "").strip()
    if not text or text.upper() in {".", "NA", "N/A"}:
        return None
    try:
        number = float(text)
        return int(number) if number.is_integer() else number
    except ValueError:
        return text


def _is_nuclear_target(chrom: str) -> bool:
    value = str(chrom or "").strip()
    if value.lower().startswith("chr"):
        value = value[3:]
    value = value.upper()
    return value in {"X", "Y"} or (
        value.isdigit() and 1 <= int(value) <= 22
    )


def _chrom_name(chrom: str) -> str:
    value = str(chrom or "").strip()
    return value[3:] if value.lower().startswith("chr") else value


def _as_float(value) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number


def _pipeline_kind(source: str) -> str:
    value = str(source or "").upper()
    if "DRAGEN" in value:
        return "dragen"
    if "NCKUH" in value or "MOSDEPTH" in value:
        return "nckuh"
    return "unknown"


def _expected_ratio(chrom: str, reference: str, estimated: str) -> float | None:
    name = _chrom_name(chrom).upper()
    if name.isdigit() and 1 <= int(name) <= 22:
        return 1.0
    baseline = reference if reference in {"XX", "XY"} else estimated
    if baseline not in {"XX", "XY"}:
        return None
    if name == "X":
        return 1.0 if baseline == "XX" else 0.5
    if name == "Y":
        return 0.0 if baseline == "XX" else 0.5
    return None


def _karyotype_interpretation(karyotype: str) -> str:
    value = str(karyotype or "").upper()
    labels = {
        "X": "possible 45,X",
        "X0": "possible 45,X",
        "XXX": "possible 47,XXX",
        "XXY": "possible 47,XXY",
        "XYY": "possible 47,XYY",
        "XXXY": "possible 48,XXXY",
    }
    if value in {"", "XX", "XY"}:
        return ""
    return labels.get(value, f"possible sex-chromosome dosage abnormality ({value})")


def _row_interpretation(chrom: str, dosage_call: str) -> str:
    name = _chrom_name(chrom).upper()
    if dosage_call not in {"gain", "loss"}:
        return "possible chromosome dosage abnormality"
    if name.isdigit() and 1 <= int(name) <= 22:
        label = "trisomy" if dosage_call == "gain" else "monosomy"
        return f"possible {label} {name}"
    if name in {"X", "Y"}:
        return "possible sex-chromosome dosage abnormality"
    return "possible chromosome dosage abnormality"


def _annotate_row(
    row: dict,
    *,
    pipeline_kind: str,
    autosome_depth,
    reference_karyotype: str,
    estimated_karyotype: str,
) -> dict:
    nuclear = _is_nuclear_target(row["chrom"])
    filter_value = str(row.get("filter") or ".")
    filter_upper = filter_value.upper()
    alt_upper = str(row.get("alt") or ".").upper()
    svtype_upper = str(row.get("svtype") or "").upper()
    native_ratio = _as_float(row.get("RATIO"))
    dc = _as_float(row.get("DC"))
    autosome = _as_float(autosome_depth)
    observed_ratio = native_ratio
    ratio_source = "native" if native_ratio is not None else ""
    if observed_ratio is None and dc is not None and autosome and autosome > 0:
        observed_ratio = dc / autosome
        ratio_source = "derived"

    confidence = (
        "pass"
        if filter_upper == "PASS"
        else "low"
        if filter_upper == "LOWQUAL"
        else "suspect"
        if filter_upper == "SUSPECT"
        else "warning"
    )
    dosage_call = "not_assessed" if not nuclear else "normal"
    explicit_gain = "<DUP>" in alt_upper or svtype_upper == "DUP"
    explicit_loss = "<DEL>" in alt_upper or svtype_upper == "DEL"
    if nuclear and explicit_gain:
        dosage_call = "gain"
    elif nuclear and explicit_loss:
        dosage_call = "loss"
    elif nuclear and pipeline_kind == "nckuh" and filter_upper != "PASS":
        ndc = _as_float(row.get("NDC"))
        signal = ndc if ndc is not None else observed_ratio
        dosage_call = "gain" if signal is not None and signal > 1 else (
            "loss" if signal is not None and signal < 1 else "suspect"
        )
    elif nuclear and pipeline_kind == "unknown" and filter_upper != "PASS":
        ndc = _as_float(row.get("NDC"))
        dosage_call = "gain" if ndc is not None and ndc > 1 else (
            "loss" if ndc is not None and ndc < 1 else "suspect"
        )

    is_abnormal = nuclear and dosage_call in {"gain", "loss", "suspect"}
    if dosage_call == "gain":
        call_label = "Gain" if confidence == "pass" else "Gain signal"
    elif dosage_call == "loss":
        call_label = "Loss" if confidence == "pass" else "Loss signal"
    elif dosage_call == "suspect":
        call_label = "Dosage signal"
    elif dosage_call == "not_assessed":
        call_label = "Not assessed"
    else:
        call_label = "Normal"

    row.update({
        "pipeline_kind": pipeline_kind,
        "expected_ratio": _expected_ratio(
            row["chrom"], reference_karyotype, estimated_karyotype
        ),
        "observed_ratio": observed_ratio,
        "ratio_source": ratio_source,
        "confidence": confidence,
        "dosage_call": dosage_call,
        "call_label": call_label,
        "is_abnormal": is_abnormal,
        "interpretation": (
            _row_interpretation(row["chrom"], dosage_call)
            if is_abnormal
            else ""
        ),
    })
    return row


def parse_ploidy_vcf(path: Path) -> dict:
    """Return VCF header metadata plus every chromosome dosage row."""
    path = Path(path)
    result = {
        "exists": path.is_file(),
        "karyotype": "",
        "reference_karyotype": "",
        "pipeline_source": "",
        "pipeline_kind": "unknown",
        "seq_type": "",
        "autosome_depth": None,
        "source": path.name if path.is_file() else "",
        "chromosomes": [],
        "warnings": [],
        "abnormal_chromosomes": [],
        "qc_warnings": [],
        "aneuploidy_suspected": False,
        "dosage_status": "unavailable",
        "karyotype_interpretation": "",
    }
    if not path.is_file():
        return result
    try:
        opener = gzip.open if path.suffix == ".gz" else open
        sample_columns: list[str] = []
        with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                line = line.rstrip("\n")
                if line.startswith("##"):
                    for header_key, result_key in _HEADER_FIELDS.items():
                        prefix = f"##{header_key}="
                        if line.startswith(prefix):
                            value = line[len(prefix):].strip()
                            result[result_key] = (
                                _clean_karyotype(value)
                                if result_key in {"karyotype", "reference_karyotype"}
                                else _number(value)
                                if result_key == "autosome_depth"
                                else value
                            )
                            break
                    continue
                if line.startswith("#CHROM"):
                    sample_columns = line.lstrip("#").split("\t")
                    continue
                if not line or line.startswith("#"):
                    continue
                fields = line.split("\t")
                if len(fields) < 8:
                    continue
                info = {}
                for item in fields[7].split(";"):
                    key, sep, value = item.partition("=")
                    if sep:
                        info[key] = value
                format_keys = fields[8].split(":") if len(fields) > 8 else []
                sample_values = fields[9].split(":") if len(fields) > 9 else []
                format_values = dict(zip(format_keys, sample_values))
                result["chromosomes"].append({
                    "chrom": fields[0],
                    "alt": fields[4] or ".",
                    "qual": _number(fields[5]),
                    "filter": fields[6] or ".",
                    "end": _number(info.get("END", "")),
                    "svtype": info.get("SVTYPE", ""),
                    "DC": _number(format_values.get("DC", "")),
                    "NDC": _number(format_values.get("NDC", "")),
                    "RATIO": _number(format_values.get("RATIO", "")),
                    "sample": sample_columns[9] if len(sample_columns) > 9 else "",
                })
    except OSError:
        result["exists"] = False
        result["source"] = ""
        result["chromosomes"] = []
        return result

    result["pipeline_kind"] = _pipeline_kind(result["pipeline_source"])
    result["chromosomes"] = [
        _annotate_row(
            row,
            pipeline_kind=result["pipeline_kind"],
            autosome_depth=result["autosome_depth"],
            reference_karyotype=result["reference_karyotype"],
            estimated_karyotype=result["karyotype"],
        )
        for row in result["chromosomes"]
    ]
    result["abnormal_chromosomes"] = [
        row for row in result["chromosomes"] if row["is_abnormal"]
    ]
    result["qc_warnings"] = [
        row
        for row in result["chromosomes"]
        if _is_nuclear_target(row["chrom"])
        and str(row["filter"]).upper() != "PASS"
    ]
    # ``warnings`` is retained as a payload compatibility alias. It now means
    # an actual dosage signal rather than every non-PASS quality record.
    result["warnings"] = result["abnormal_chromosomes"]
    karyotype = str(result["karyotype"] or "")
    result["karyotype_interpretation"] = _karyotype_interpretation(karyotype)
    result["aneuploidy_suspected"] = (
        bool(karyotype and karyotype not in {"XX", "XY"})
        or bool(result["abnormal_chromosomes"])
    )
    result["dosage_status"] = (
        "aneuploidy_signal" if result["aneuploidy_suspected"] else "no_signal"
    )
    return result


def read_karyotype(path: Path) -> str:
    """Read estimatedSexKaryotype, preserving calls such as XXY or X."""
    return str(parse_ploidy_vcf(path).get("karyotype") or "")


def _path_for_sample(sample: str | Path) -> Path:
    if isinstance(sample, str):
        return sample_layout.state_file(sample, "ploidy.vcf.gz")
    directory = Path(sample)
    sample_id = (
        directory.parent.name
        if directory.name == sample_layout.POSTPROCESSING_DIRNAME
        else directory.name
    )
    prefixed = directory / sample_layout.prefixed_filename(sample_id, "ploidy.vcf.gz")
    if prefixed.is_file():
        return prefixed
    return directory / "ploidy.vcf.gz"


def load_sample_ploidy(sample: str | Path) -> dict:
    """Load a sample's copied VCF only; ploidy_qc.txt is intentionally ignored."""
    return parse_ploidy_vcf(_path_for_sample(sample))
