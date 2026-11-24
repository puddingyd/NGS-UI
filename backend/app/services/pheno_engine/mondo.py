"""MONDO equivalence groups for OMIM / Orphanet diseases (plan §4.2.4)."""
from __future__ import annotations

import json
import re
from pathlib import Path

_XREF = re.compile(r"^xref:\s*(OMIM|Orphanet|ORPHA):(\d+)(.*)$")


def _norm(prefix: str, number: str) -> str:
    return f"{'OMIM' if prefix == 'OMIM' else 'ORPHA'}:{number}"


def parse(path: Path) -> tuple[dict[str, list[str]], dict[str, str], dict]:
    """({mondo_id: [disease ids]}, {mondo_id: label}, info).

    OBO: only ``xref`` lines annotated ``MONDO:equivalentTo``.  JSON (OBO
    graphs) carries no such qualifier in older releases, so every OMIM /
    Orphanet xref counts there; ``info["mode"]`` records which was used.
    """
    path = Path(path)
    groups: dict[str, list[str]] = {}
    labels: dict[str, str] = {}
    version = ""
    if path.suffix == ".json":
        data = json.loads(path.read_text(encoding="utf-8"))
        graph = (data.get("graphs") or [{}])[0]
        version = str((graph.get("meta") or {}).get("version") or "")
        for node in graph.get("nodes") or []:
            nid = str(node.get("id") or "")
            if "MONDO_" not in nid:
                continue
            meta = node.get("meta") or {}
            if meta.get("deprecated"):
                continue
            mondo = "MONDO:" + nid.rsplit("MONDO_", 1)[1]
            members = []
            for xref in meta.get("xrefs") or []:
                val = xref.get("val") if isinstance(xref, dict) else str(xref)
                m = re.match(r"^(OMIM|Orphanet|ORPHA):(\d+)$", str(val or ""))
                if m:
                    members.append(_norm(m.group(1), m.group(2)))
            if members:
                groups[mondo] = sorted(set(members))
                labels[mondo] = str(node.get("lbl") or "")
        return groups, labels, {"mode": "json-all-xrefs", "version": version}
    cur = None
    obsolete = False
    label = ""
    members: list[str] = []

    def flush():
        if cur and not obsolete and members:
            groups[cur] = sorted(set(members))
            labels[cur] = label

    with open(path, encoding="utf-8") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line.startswith("data-version:"):
                version = line.split(":", 1)[1].strip()
            if line.startswith("["):
                flush()
                cur, obsolete, label, members = None, False, "", []
                continue
            if line.startswith("id: MONDO:"):
                cur = line[4:].strip()
            elif cur is None:
                continue
            elif line.startswith("name: "):
                label = line[6:].strip()
            elif line.startswith("is_obsolete:") and "true" in line:
                obsolete = True
            else:
                m = _XREF.match(line)
                if m and "MONDO:equivalentTo" in m.group(3):
                    members.append(_norm(m.group(1), m.group(2)))
        flush()
    return groups, labels, {"mode": "obo-equivalentTo", "version": version}


def disease_groups(groups: dict[str, list[str]], universe: set[str]) -> tuple[list[tuple], dict]:
    """Rows (mondo_id, disease_id) for diseases in ``universe`` claimed by exactly one MONDO term."""
    owners: dict[str, list[str]] = {}
    for mondo, members in groups.items():
        for d in members:
            owners.setdefault(d, []).append(mondo)
    rows = []
    ambiguous = 0
    for d in sorted(universe):
        claim = owners.get(d) or []
        if len(claim) == 1:
            rows.append((claim[0], d))
        elif len(claim) > 1:
            ambiguous += 1
    by_group: dict[str, set] = {}
    for mondo, d in rows:
        by_group.setdefault(mondo, set()).add(d)
    mixed = sum(1 for ds in by_group.values()
                if any(x.startswith("OMIM:") for x in ds) and any(x.startswith("ORPHA:") for x in ds))
    return sorted(rows), {"grouped_diseases": len(rows), "groups": len(by_group),
                          "groups_omim_and_orpha": mixed, "ambiguous_diseases": ambiguous}
