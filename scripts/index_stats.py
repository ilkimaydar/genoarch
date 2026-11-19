#!/usr/bin/env python3
"""Descriptive statistics of the knowledge base built from an HPO release (for docs/DATA_CARD.md).

Usage: python scripts/index_stats.py --hpo-dir data/hpo --out stats.json
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import Counter
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "core"))

from genoarch_core.annotations import combine_frequencies, read_phenotype_hpoa  # noqa: E402
from genoarch_core.knowledge import KnowledgeBase  # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--hpo-dir", default=str(ROOT / "data" / "hpo"))
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    hpo = Path(args.hpo_dir)
    kb = KnowledgeBase.from_hpo_files(hpo)
    src = kb.disease_source
    alen = np.asarray([len(a) for a in kb.ann_terms])
    per_source = {}
    for s in sorted(set(src)):
        m = np.asarray([x == s for x in src])
        per_source[s] = {"diseases": int(m.sum()), "annotations_median": float(np.median(alen[m])), "annotations_mean": round(float(alen[m].mean()), 3),
                         "diseases_with_gene": int(sum(1 for d in np.nonzero(m)[0] if kb.disease_genes[int(d)]))}
    gene_src = Counter()
    for g in range(kb.n_genes):
        s = {src[int(d)] for d in kb.gene_diseases[g]}
        gene_src["with_omim_link" if "OMIM" in s else ("orphanet_only" if s == {"ORPHA"} else "other")] += 1

    # unmapped disease-term pairs by aspect
    hpoa = read_phenotype_hpoa(hpo / "phenotype.hpoa")
    aspect: dict[tuple[str, str], str] = {}
    with open(hpo / "phenotype.hpoa", encoding="utf-8") as fh:
        rows = [line for line in fh if not line.startswith("#")]
    for row in csv.DictReader(rows, delimiter="\t"):
        if (row.get("qualifier") or "").upper() != "NOT":
            aspect.setdefault((row["database_id"], row["hpo_id"]), row.get("aspect") or "")
    unmapped = Counter()
    zero_pairs = set()
    for did, tmap in hpoa.positive.items():
        for tid, obs in tmap.items():
            i, status = kb.resolve(tid)
            if i is None:
                unmapped[f"aspect_{aspect.get((did, tid), '?')}_{status}"] += 1
            f = combine_frequencies(obs)
            if f is not None and f <= 0:
                zero_pairs.add((did, tid))

    # genes_to_phenotype profile vs union of disease annotations: explain the differences
    union = {}
    for g in range(kb.n_genes):
        s: set[int] = set()
        for d in kb.gene_diseases[g]:
            s.update(int(x) for x in kb.prop[int(d)])
        union[kb.gene_symbols[g]] = s
    g2p_rows: dict[str, list[tuple[str, str]]] = {}
    with open(hpo / "genes_to_phenotype.txt", encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            g2p_rows.setdefault(row["gene_symbol"], []).append((row["disease_id"], row["hpo_id"]))
    equal = only_profile_genes = only_union_genes = 0
    reasons = Counter()
    for g in range(kb.n_genes):
        sym = kb.gene_symbols[g]
        prof = kb.gene_prop[g] if kb.gene_prop is not None else []
        prof_set = set(int(x) for x in prof)
        if prof_set == union[sym]:
            equal += 1
            continue
        if prof_set - union[sym]:
            only_profile_genes += 1
        if union[sym] - prof_set:
            only_union_genes += 1
        linked = {kb.disease_ids[int(d)] for d in kb.gene_diseases[g]}
        for did, tid in g2p_rows.get(sym, []):
            i, _ = kb.resolve(tid)
            if i is None or i in union[sym]:
                continue
            if (did, tid) in zero_pairs:
                reasons["row_pair_has_pooled_frequency_0"] += 1
            elif did not in linked:
                reasons["row_disease_not_linked_to_gene_in_genes_to_disease"] += 1
            elif (did, tid) in hpoa.negative:
                reasons["row_pair_is_NOT_qualified"] += 1
            else:
                reasons["other"] += 1
    out = {
        "hpo_dir": str(hpo), "ontology_version": kb.meta.get("ontology_version"),
        "annotation_version": kb.meta.get("annotation_header", {}).get("version"),
        "counts": {"terms": kb.n_terms, "diseases": kb.n_diseases, "genes": kb.n_genes, "annotations": int(alen.sum())},
        "meta": {k: v for k, v in kb.meta.items() if k.startswith("n_")},
        "per_source": per_source, "genes_by_link_source": dict(gene_src),
        "unmapped_pairs_by_aspect_and_status": dict(unmapped),
        "gene_profile_check": {"genes": kb.n_genes, "equal": equal, "genes_with_extra_profile_terms": only_profile_genes,
                               "genes_with_union_terms_missing_from_profile": only_union_genes,
                               "extra_profile_rows_by_reason": dict(reasons)},
        "ic": {"max": float(kb.ic.max()), "terms_with_annotations": int((kb.term_count > 0).sum())},
    }
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(out, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
