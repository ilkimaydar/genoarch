"""Compact browser index (gzip-compressed JSON) shared by the Python and JS engines.

Layout (format ``genoarch-index`` v1). All term / disease references are
integer positions; sorted integer lists are delta-encoded.

* ``terms.hp``    HP numbers (term i is ``HP:%07d``), Phenotypic abnormality subtree
* ``terms.name``, ``terms.syn``  labels and synonyms (search)
* ``terms.par``   ``is_a`` parents (positions)
* ``terms.anc``   ancestor closure incl. self (delta-encoded positions)
* ``terms.ic6``   information content x 1e6 (integers), ``terms.cnt`` annotated-disease counts
* ``alt``         ``[alt HP number, position]`` for alternative / obsolete ids
* ``other``       HP numbers outside the Phenotypic abnormality subtree
* ``diseases``    ``id``, ``name``, ``ann`` (delta-encoded term positions), ``f4`` (frequency x 1e4, -1 unknown)
* ``genes``       ``sym``, ``ncbi``, ``dis`` (delta-encoded disease positions), ``assoc`` (0 Mendelian, 1 polygenic, 2 unknown)

Gene phenotype profiles are not stored: they are the union of the propagated
annotations of the gene's linked diseases (see docs/DATA_CARD.md).
"""
from __future__ import annotations

import gzip
import json
from pathlib import Path

import numpy as np

from .knowledge import ASSOCIATION_CODES, FREQ_SCALE, IC_SCALE, KnowledgeBase

FORMAT = "genoarch-index"
FORMAT_VERSION = 1


def _delta(values) -> list[int]:
    out, prev = [], 0
    for v in values:
        v = int(v)
        out.append(v - prev)
        prev = v
    return out


def _undelta(values) -> list[int]:
    out, acc = [], 0
    for v in values:
        acc += int(v)
        out.append(acc)
    return out


def _hp_num(term_id: str) -> int:
    return int(term_id.split(":", 1)[1])


def index_dict(kb: KnowledgeBase, *, include_ancestors: bool = True) -> dict:
    ic6 = np.rint(kb.ic * IC_SCALE).astype(np.int64)
    terms = {
        "hp": [_hp_num(t) for t in kb.term_ids],
        "name": kb.names,
        "syn": kb.synonyms,
        "par": [[int(x) for x in p] for p in kb.parents],
        "ic6": [int(x) for x in ic6],
        "cnt": [int(x) for x in kb.term_count],
    }
    if include_ancestors:
        terms["anc"] = [_delta(a) for a in kb.ancestors]
    f4 = []
    for f in kb.ann_freqs:
        f4.append([-1 if np.isnan(x) else int(np.rint(x * FREQ_SCALE)) for x in f])
    meta = {k: v for k, v in kb.meta.items() if k not in ("hpo_dir",)}
    return {
        "format": FORMAT,
        "version": FORMAT_VERSION,
        "meta": meta,
        "terms": terms,
        "alt": sorted([_hp_num(a), int(i)] for a, i in kb.alt_map.items()),
        "other": sorted(_hp_num(t) for t in kb.other_terms),
        "other_names": [kb.other_terms[t] for t in sorted(kb.other_terms, key=_hp_num)],
        "diseases": {
            "id": kb.disease_ids,
            "name": kb.disease_names,
            "ann": [_delta(t) for t in kb.ann_terms],
            "f4": f4,
        },
        "genes": {
            "sym": kb.gene_symbols,
            "ncbi": kb.gene_ncbi,
            "dis": [_delta(d) for d in kb.gene_diseases],
            "assoc": [[ASSOCIATION_CODES.index(a) if a in ASSOCIATION_CODES else 2 for a in assoc] for assoc in kb.gene_assoc],
        },
    }


def write_index(kb: KnowledgeBase, path: str | Path, *, include_ancestors: bool = True) -> dict:
    """Write the gzip-compressed index. Returns size information."""
    doc = index_dict(kb, include_ancestors=include_ancestors)
    raw = json.dumps(doc, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    # mtime=0 keeps the archive byte-identical for identical content
    with open(path, "wb") as fh:
        with gzip.GzipFile(filename="", mode="wb", fileobj=fh, compresslevel=9, mtime=0) as gz:
            gz.write(raw)
    return {"json_bytes": len(raw), "gzip_bytes": path.stat().st_size}


def load_index_doc(path: str | Path) -> dict:
    data = Path(path).read_bytes()
    if data[:2] == b"\x1f\x8b":
        data = gzip.decompress(data)
    doc = json.loads(data.decode("utf-8"))
    if doc.get("format") != FORMAT:
        raise ValueError(f"not a {FORMAT} file: {path}")
    return doc


def kb_from_doc(doc: dict) -> KnowledgeBase:
    t = doc["terms"]
    term_ids = [f"HP:{n:07d}" for n in t["hp"]]
    ic = np.asarray(t["ic6"], dtype=np.float64) / IC_SCALE
    alt_map = {f"HP:{a:07d}": int(i) for a, i in doc["alt"]}
    other_terms = {f"HP:{n:07d}": name for n, name in zip(doc["other"], doc.get("other_names", [""] * len(doc["other"])))}
    d = doc["diseases"]
    ann_terms = [_undelta(x) for x in d["ann"]]
    ann_freqs = [[np.nan if k < 0 else k / FREQ_SCALE for k in f] for f in d["f4"]]
    g = doc["genes"]
    kb = KnowledgeBase(
        term_ids=term_ids,
        names=t["name"],
        synonyms=t["syn"],
        parents=t["par"],
        alt_map=alt_map,
        disease_ids=d["id"],
        disease_names=d["name"],
        ann_terms=ann_terms,
        ann_freqs=ann_freqs,
        gene_symbols=g["sym"],
        gene_ncbi=g["ncbi"],
        gene_diseases=[_undelta(x) for x in g["dis"]],
        gene_assoc=[[ASSOCIATION_CODES[c] for c in a] for a in g["assoc"]],
        ic=ic,
        term_count=np.asarray(t["cnt"], dtype=np.int64),
        other_terms=other_terms,
        meta=dict(doc.get("meta", {}), source="index"),
    )
    if "anc" in t:
        stored = [_undelta(a) for a in t["anc"]]
        for i, a in enumerate(stored):
            if a != [int(x) for x in kb.ancestors[i]]:
                raise ValueError(f"stored ancestor closure differs from is_a closure at {term_ids[i]}")
    return kb


def read_index(path: str | Path) -> KnowledgeBase:
    return kb_from_doc(load_index_doc(path))
