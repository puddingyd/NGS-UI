"""Parse the user's phenotype.txt + write helpers.

Input format (tab-separated, first row is header):

    phenotype<TAB>hpo_name<TAB>weight
    HP:0001508<TAB>Failure to thrive<TAB>1
    HP:0002902<TAB>Hyponatremia<TAB>2
    growth_panel<TAB><TAB>1

Rows where the second column has a name are HPO terms. Rows where the
second column is empty (column 1 is anything that isn't an HP: id) are
treated as panel names. Weight defaults to 1 when missing or unparseable.

Filenames on disk land as <LIS_ID>_<MRN>_phenotype.txt inside the
analysis version directory so the source-of-truth is auditable.
"""
from __future__ import annotations

from pathlib import Path


HEADER_TOKENS = {"phenotype", "hpo_name", "weight"}


def parse(text: str) -> tuple[list[dict], list[dict]]:
    """Return (hpo_list, panels_list) parsed from the txt body.

    Each HPO entry: {"phenotype": "HP:...", "label": "...", "weight": float}
    Each panel entry: {"name": "...", "weight": float}
    """
    hpo: list[dict] = []
    panels: list[dict] = []
    if not text:
        return hpo, panels

    lines = text.splitlines()
    if not lines:
        return hpo, panels

    # Detect and drop a header row: any line whose tokens overlap with
    # the canonical header keywords. Keeps body rows that happen to
    # start without an HP: prefix (e.g., a panel) parseable.
    first = [c.strip().lower() for c in lines[0].split("\t")]
    if any(tok in HEADER_TOKENS for tok in first):
        lines = lines[1:]

    for raw in lines:
        if not raw.strip():
            continue
        parts = raw.split("\t")
        col1 = parts[0].strip() if len(parts) > 0 else ""
        col2 = parts[1].strip() if len(parts) > 1 else ""
        col3 = parts[2].strip() if len(parts) > 2 else ""
        if not col1:
            continue
        try:
            weight = float(col3) if col3 else 1.0
        except ValueError:
            weight = 1.0
        is_hpo = col1.upper().startswith("HP:")
        if is_hpo:
            hpo.append({
                "phenotype": col1,
                "label":     col2 or col1,
                "weight":    weight,
            })
        else:
            # Panel name. Second column should be empty per the spec
            # but we don't enforce it — some hand-edited files leave a
            # display name there.
            panels.append({"name": col1, "weight": weight})
    return hpo, panels


def write(hpo: list[dict], panels: list[dict], path: Path) -> Path:
    """Write a phenotype.txt-formatted file, with the canonical header.

    Used when load-new-case re-emits the parsed file into the analysis
    version directory so the on-disk artefact matches what the backend
    actually used (handy for audit + manual re-runs).
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["phenotype\thpo_name\tweight"]
    for h in hpo or []:
        hid = h.get("phenotype") or h.get("hpo_id") or ""
        name = h.get("label") or h.get("hpo_name") or ""
        w = h.get("weight", 1)
        lines.append(f"{hid}\t{name}\t{_fmt_weight(w)}")
    for p in panels or []:
        name = p.get("name") if isinstance(p, dict) else str(p)
        w = p.get("weight", 1) if isinstance(p, dict) else 1
        lines.append(f"{name}\t\t{_fmt_weight(w)}")
    body = "\n".join(lines) + "\n"
    path.write_text(body, encoding="utf-8")
    return path


def _fmt_weight(w) -> str:
    try:
        n = float(w)
    except (TypeError, ValueError):
        return "1"
    return str(int(n)) if n.is_integer() else f"{n:g}"
