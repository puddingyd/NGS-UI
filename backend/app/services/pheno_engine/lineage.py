"""HPO lineage for explaining a phenotype match (direct parents, shortest paths, organ systems).

The knowledge base stores only the ancestor closure of each term (inside the
phenotype subtree); direct parents are recovered by transitive reduction of
that closure, so the answer always matches the snapshot in use.  Everything
is computed lazily for the handful of terms a picture needs.
"""
from __future__ import annotations

from collections import deque

OTHER_ANNOTATION_LIMIT = 80


class Lineage:
    def __init__(self, engine):
        self.engine = engine
        self.kb = engine.kb
        self.root = engine.root                 # HP:0000118 Phenotypic abnormality
        self._parents: dict[int, tuple] = {}

    def anc(self, t: int) -> frozenset:
        return self.engine._anc[t]

    def parents(self, t: int) -> tuple:
        """Direct is_a parents inside the phenotype subtree (transitive reduction of the closure)."""
        hit = self._parents.get(t)
        if hit is None:
            up = [a for a in self.anc(t) if a != t]
            hit = tuple(sorted(p for p in up if not any(a != p and p in self.anc(a) for a in up)))
            self._parents[t] = hit
        return hit

    def is_organ(self, t: int) -> bool:
        return self.root is not None and t != self.root and self.root in self.parents(t)

    def path_up(self, start: int, target: int) -> list[int] | None:
        """Shortest chain start → … → target following parents (target must be an ancestor)."""
        if start == target:
            return [start]
        if target not in self.anc(start):
            return None
        prev = {start: None}
        queue = deque([start])
        while queue:
            x = queue.popleft()
            for p in self.parents(x):
                if p in prev or target not in self.anc(p):
                    continue
                prev[p] = x
                if p == target:
                    path = [p]
                    while prev[path[-1]] is not None:
                        path.append(prev[path[-1]])
                    return path[::-1]
                queue.append(p)
        return None

    def path_to_organ(self, start: int) -> list[int]:
        """Shortest chain from a term up to its nearest organ-system term (a child of the root)."""
        if self.root is None or start == self.root:
            return [start]
        prev = {start: None}
        queue = deque([start])
        while queue:
            x = queue.popleft()
            if self.is_organ(x):
                path = [x]
                while prev[path[-1]] is not None:
                    path.append(prev[path[-1]])
                return path[::-1]
            for p in self.parents(x):
                if p not in prev and p != self.root:
                    prev[p] = x
                    queue.append(p)
        return [start]


def lineage(engine) -> Lineage:
    hit = getattr(engine, "_lineage", None)
    if hit is None:
        hit = Lineage(engine)
        engine._lineage = hit
    return hit


def pair_steps(lin: Lineage, q: int, d: int, relation: str, mica: int) -> dict:
    """Steps and paths of one patient term ↔ disease annotation pairing (patient side first)."""
    if relation == "same":
        return {"steps_patient": 0, "steps_annotation": 0, "path_patient": [q], "path_annotation": [d]}
    if relation == "patient_specific":                 # annotation is an ancestor of the patient term
        path = lin.path_up(q, d) or [q, d]
        return {"steps_patient": len(path) - 1, "steps_annotation": 0, "path_patient": path, "path_annotation": [d]}
    if relation == "patient_broader":                  # annotation is a descendant of the patient term
        path = lin.path_up(d, q) or [d, q]
        return {"steps_patient": 0, "steps_annotation": len(path) - 1, "path_patient": [q], "path_annotation": path}
    if relation == "collateral":
        pq = lin.path_up(q, mica) or [q, mica]
        pd = lin.path_up(d, mica) or [d, mica]
        return {"steps_patient": len(pq) - 1, "steps_annotation": len(pd) - 1, "path_patient": pq, "path_annotation": pd}
    return {"steps_patient": None, "steps_annotation": None, "path_patient": [q], "path_annotation": []}


def disease_graph(engine, result, disease: int, *, all_annotations: bool = False) -> dict:
    """Sub-graph of the HPO tree linking the patient terms to one disease's annotations.

    Nodes: patient terms, their paired annotations, the meeting term of each
    pair, the shortest chains between them and up to each organ system;
    optionally the disease's other annotations.  Edges: is_a child → parent."""
    kb = engine.kb
    lin = lineage(engine)
    roles: dict[int, set] = {}

    def add(path, role=None):
        for t in path:
            roles.setdefault(int(t), set())
        if role and path:
            roles[int(path[0])].add(role)

    pairs = []
    paired = set()
    for k, q in enumerate(result.q):
        q = int(q)
        d = engine.best_annotation(q, result.evecs[k], disease)
        p = engine.pair(q, d)
        rel = p["relation"]
        steps = pair_steps(lin, q, d, rel, int(p["mica"]))
        add(steps["path_patient"], "patient")
        if rel == "none":
            add(lin.path_to_organ(q))
        else:
            add(steps["path_annotation"], "annotation")
            paired.add(d)
            mica = int(p["mica"])
            if rel == "collateral":          # only side branches meet at a third term
                roles.setdefault(mica, set()).add("meet")
            add(lin.path_to_organ(mica))
        pairs.append({"hpo_id": kb.terms[q], "annotation": kb.terms[d] if rel != "none" else "",
                      "relation": rel, "mica": kb.terms[int(p["mica"])] if rel != "none" else "",
                      "path_patient": [kb.terms[t] for t in steps["path_patient"]],
                      "path_annotation": [kb.terms[t] for t in steps["path_annotation"]],
                      "steps_patient": steps["steps_patient"], "steps_annotation": steps["steps_annotation"]})
    others = 0
    if all_annotations:
        for a in kb.dir_idx[kb.dir_ptr[disease]:kb.dir_ptr[disease + 1]]:
            a = int(a)
            if a in paired:
                continue
            if others >= OTHER_ANNOTATION_LIMIT:
                break
            others += 1
            add(lin.path_to_organ(a), "other_annotation")
    for t in list(roles):
        if lin.is_organ(t):
            roles[t].add("organ")
    nodes = [{"id": kb.terms[t], "name": kb.names[t], "ic": round(float(kb.ic_star[t]), 3),
              "n_diseases": int(kb.n[t]), "roles": sorted(roles[t])} for t in sorted(roles)]
    edges = [[kb.terms[t], kb.terms[p]] for t in sorted(roles) for p in lin.parents(t) if p in roles]
    n_ann = int(kb.dir_ptr[disease + 1] - kb.dir_ptr[disease])
    return {"disease_id": kb.diseases[disease], "disease_name": kb.disease_names[disease], "n_annotations": n_ann,
            "nodes": nodes, "edges": edges, "pairs": pairs, "others_shown": others,
            "others_total": max(n_ann - len(paired), 0)}
