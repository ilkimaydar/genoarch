"""Human Phenotype Ontology (HPO) loading.

Reads the official release files (``hp.json`` in OBO Graphs JSON format, or
``hp.obo``) and returns an :class:`Ontology` with current terms, ``is_a``
parents, synonyms and a map from alternative / obsolete identifiers to the
current term.

The ontology content is used unchanged: no terms or relationships are added
or altered (HPO license, item 2).
"""
from __future__ import annotations

import json
import re
from dataclasses import dataclass, field
from pathlib import Path

ROOT = "HP:0000001"
PHENOTYPIC_ABNORMALITY = "HP:0000118"

_PURL = "http://purl.obolibrary.org/obo/"
_ALT_ID = "http://www.geneontology.org/formats/oboInOwl#hasAlternativeId"
_REPLACED_BY = "http://purl.obolibrary.org/obo/IAO_0100001"
_VERSION_INFO = "http://www.w3.org/2002/07/owl#versionInfo"
_HP_RE = re.compile(r"^HP:\d{7}$")
_HP_LOOSE = re.compile(r"^(?:HP)[:_]?(\d{1,7})$")


def curie(iri: str) -> str:
    """Convert an OBO PURL (``.../obo/HP_0000118``) to a CURIE (``HP:0000118``)."""
    if iri.startswith(_PURL):
        local = iri[len(_PURL):]
        prefix, sep, rest = local.partition("_")
        if sep:
            return f"{prefix}:{rest}"
    return iri


def normalize_hpo_id(raw: object) -> str | None:
    """Return a canonical ``HP:NNNNNNN`` identifier or ``None`` if not parseable.

    Accepts ``HP:0001250``, ``hp:1250``, ``HP_0001250`` and OBO PURLs.
    """
    if raw is None:
        return None
    s = str(raw).strip()
    if not s:
        return None
    if s.lower().startswith("http"):
        s = curie(s)
    m = _HP_LOOSE.match(s.upper())
    if not m:
        return None
    return "HP:" + m.group(1).zfill(7)


@dataclass
class Term:
    id: str
    name: str
    parents: list[str] = field(default_factory=list)
    synonyms: list[str] = field(default_factory=list)
    alt_ids: list[str] = field(default_factory=list)
    obsolete: bool = False
    replaced_by: str | None = None


@dataclass
class Ontology:
    """Current (non-obsolete) HPO terms plus identifier mapping."""

    terms: dict[str, Term]
    version: str | None
    alt_map: dict[str, str]
    n_obsolete: int = 0

    def resolve(self, term_id: str) -> str | None:
        """Map any HPO id (current, alternative or obsolete) to a current id."""
        tid = normalize_hpo_id(term_id)
        if tid is None:
            return None
        if tid in self.terms:
            return tid
        return self.alt_map.get(tid)

    def children(self) -> dict[str, list[str]]:
        out: dict[str, list[str]] = {t: [] for t in self.terms}
        for t in self.terms.values():
            for p in t.parents:
                out.setdefault(p, []).append(t.id)
        return out

    def subtree(self, root: str = PHENOTYPIC_ABNORMALITY) -> set[str]:
        """Return ``root`` and all of its descendants."""
        kids = self.children()
        seen = {root}
        stack = [root]
        while stack:
            cur = stack.pop()
            for c in kids.get(cur, ()):  # noqa: B020
                if c not in seen:
                    seen.add(c)
                    stack.append(c)
        return seen


def _finalize(terms: dict[str, Term], obsolete: dict[str, Term], version: str | None) -> Ontology:
    alt_map: dict[str, str] = {}
    for t in terms.values():
        t.parents = [p for p in dict.fromkeys(t.parents) if p in terms and p != t.id]
        for a in t.alt_ids:
            a = normalize_hpo_id(a)
            if a and a not in terms:
                alt_map[a] = t.id
    for o in obsolete.values():
        target = o.replaced_by
        hops = 0
        while target and target not in terms and hops < 5:
            nxt = obsolete.get(target)
            target = nxt.replaced_by if nxt else alt_map.get(target)
            hops += 1
        if target and target in terms:
            alt_map.setdefault(o.id, target)
    return Ontology(terms=terms, version=version, alt_map=alt_map, n_obsolete=len(obsolete))


def load_hpo_json(path: str | Path) -> Ontology:
    """Load ``hp.json`` (OBO Graphs JSON)."""
    with open(path, encoding="utf-8") as fh:
        doc = json.load(fh)
    graph = doc["graphs"][0]
    gmeta = graph.get("meta", {}) or {}
    version = None
    for bpv in gmeta.get("basicPropertyValues", []) or []:
        if bpv.get("pred") == _VERSION_INFO:
            version = bpv.get("val")
    if version is None:
        version = gmeta.get("version")
    terms: dict[str, Term] = {}
    obsolete: dict[str, Term] = {}
    for node in graph.get("nodes", []):
        if node.get("type") != "CLASS":
            continue
        tid = curie(node.get("id", ""))
        if not _HP_RE.match(tid):
            continue
        meta = node.get("meta", {}) or {}
        exact: list[str] = []
        other: list[str] = []
        for syn in meta.get("synonyms", []) or []:
            if str(syn.get("synonymType", "")).endswith("obsolete_synonym"):
                continue
            val = str(syn.get("val", "")).strip()
            if not val:
                continue
            (exact if syn.get("pred") == "hasExactSynonym" else other).append(val)
        alt_ids: list[str] = []
        replaced_by = None
        for bpv in meta.get("basicPropertyValues", []) or []:
            if bpv.get("pred") == _ALT_ID:
                alt_ids.append(str(bpv.get("val")))
            elif bpv.get("pred") == _REPLACED_BY:
                replaced_by = normalize_hpo_id(curie(str(bpv.get("val"))))
        term = Term(
            id=tid,
            name=str(node.get("lbl") or tid),
            synonyms=list(dict.fromkeys(exact + other)),
            alt_ids=alt_ids,
            obsolete=bool(meta.get("deprecated")),
            replaced_by=replaced_by,
        )
        (obsolete if term.obsolete else terms)[tid] = term
    for edge in graph.get("edges", []):
        if edge.get("pred") != "is_a":
            continue
        sub, obj = curie(edge.get("sub", "")), curie(edge.get("obj", ""))
        if sub in terms and _HP_RE.match(obj):
            terms[sub].parents.append(obj)
    return _finalize(terms, obsolete, version)


_SYN_RE = re.compile(r'^"((?:[^"\\]|\\.)*)"\s+(EXACT|RELATED|BROAD|NARROW)(.*)$')


def load_hpo_obo(path: str | Path) -> Ontology:
    """Load ``hp.obo`` (OBO 1.2 flat file)."""
    terms: dict[str, Term] = {}
    obsolete: dict[str, Term] = {}
    version = None
    cur: dict | None = None
    in_term = False

    def flush(rec: dict | None) -> None:
        if not rec or not rec.get("id") or not _HP_RE.match(rec["id"]):
            return
        exact = [s for s, kind in rec["syn"] if kind == "EXACT"]
        other = [s for s, kind in rec["syn"] if kind != "EXACT"]
        term = Term(
            id=rec["id"],
            name=rec.get("name") or rec["id"],
            parents=rec["is_a"],
            synonyms=list(dict.fromkeys(exact + other)),
            alt_ids=rec["alt"],
            obsolete=rec["obsolete"],
            replaced_by=rec.get("replaced_by"),
        )
        (obsolete if term.obsolete else terms)[term.id] = term

    with open(path, encoding="utf-8") as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if line.startswith("["):
                if in_term:
                    flush(cur)
                in_term = line.strip() == "[Term]"
                cur = {"syn": [], "is_a": [], "alt": [], "obsolete": False} if in_term else None
                continue
            if not in_term:
                if line.startswith("data-version:") and version is None:
                    version = line.split(":", 1)[1].strip().split("/")[-1]
                continue
            if not line or ":" not in line:
                continue
            key, _, val = line.partition(":")
            val = val.strip()
            if key == "id":
                cur["id"] = val
            elif key == "name":
                cur["name"] = val
            elif key == "is_a":
                # "is_a: HP:0030991 {xref=\"PMID:1\"} ! name" -> first token only
                target = val.split("!")[0].split("{")[0].strip().split()
                if target:
                    cur["is_a"].append(target[0])
            elif key == "alt_id":
                cur["alt"].append(val)
            elif key == "is_obsolete":
                cur["obsolete"] = val.lower() == "true"
            elif key == "replaced_by":
                cur["replaced_by"] = normalize_hpo_id(val)
            elif key == "synonym":
                m = _SYN_RE.match(val)
                if m and "obsolete_synonym" not in m.group(3):
                    cur["syn"].append((m.group(1).replace('\\"', '"'), m.group(2)))
    if in_term:
        flush(cur)
    return _finalize(terms, obsolete, version)


def load_ontology(path: str | Path) -> Ontology:
    """Load ``hp.json`` or ``hp.obo`` depending on the file extension."""
    p = Path(path)
    if p.suffix.lower() == ".obo":
        return load_hpo_obo(p)
    return load_hpo_json(p)
