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


def parse_ploidy_vcf(path: Path) -> dict:
    """Return VCF header metadata plus every chromosome dosage row."""
    path = Path(path)
    result = {
        "exists": path.is_file(),
        "karyotype": "",
        "reference_karyotype": "",
        "pipeline_source": "",
        "seq_type": "",
        "source": path.name if path.is_file() else "",
        "chromosomes": [],
        "warnings": [],
        "aneuploidy_suspected": False,
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
                    "filter": fields[6] or ".",
                    "end": _number(info.get("END", "")),
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

    result["warnings"] = [
        row
        for row in result["chromosomes"]
        if _is_nuclear_target(row["chrom"])
        and str(row["filter"]).upper() != "PASS"
    ]
    karyotype = str(result["karyotype"] or "")
    result["aneuploidy_suspected"] = (
        bool(karyotype and karyotype not in {"XX", "XY"})
        or bool(result["warnings"])
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
