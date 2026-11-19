#!/usr/bin/env python3
"""Build the compact Genoarch index (web/data/index.json.gz) from an HPO release directory.

Usage: python scripts/build_index.py --hpo-dir data/hpo --out web/data/index.json.gz
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "core"))

from genoarch_core.index_io import read_index, write_index  # noqa: E402
from genoarch_core.knowledge import KnowledgeBase  # noqa: E402

HPO_FILES = ["hp.json", "phenotype.hpoa", "genes_to_disease.txt", "genes_to_phenotype.txt"]
HPO_RELEASE_URL = "https://github.com/obophenotype/human-phenotype-ontology/releases/download/{tag}/{name}"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--hpo-dir", default=str(ROOT / "data" / "hpo"))
    ap.add_argument("--out", default=str(ROOT / "web" / "data" / "index.json.gz"))
    ap.add_argument("--no-ancestors", action="store_true", help="omit the stored ancestor closure")
    ap.add_argument("--report", default=None, help="optional JSON build report path")
    args = ap.parse_args()

    hpo_dir = Path(args.hpo_dir)
    t0 = time.time()
    kb = KnowledgeBase.from_hpo_files(hpo_dir)
    tag = "v" + str(kb.meta.get("ontology_version") or "unknown")
    kb.meta["inputs"] = {
        name: {"sha256": sha256(hpo_dir / name), "bytes": (hpo_dir / name).stat().st_size,
               "url": HPO_RELEASE_URL.format(tag=tag, name=name)}
        for name in HPO_FILES if (hpo_dir / name).exists()
    }
    kb.meta["counts"] = {"terms": kb.n_terms, "diseases": kb.n_diseases, "genes": kb.n_genes,
                         "annotations": int(sum(len(a) for a in kb.ann_terms))}
    kb.meta["profile_union_check"] = {k: v for k, v in kb.gene_profile_union_check().items() if k != "examples"}
    sizes = write_index(kb, args.out, include_ancestors=not args.no_ancestors)
    kb2 = read_index(args.out)  # round-trip check (raises on mismatch)
    same = (kb2.term_ids == kb.term_ids and kb2.disease_ids == kb.disease_ids
            and bool((kb2.ic == kb.ic).all()) and kb2.gene_symbols == kb.gene_symbols
            and np.array_equal(kb2.flat_ann_freq_eff, kb.flat_ann_freq_eff)
            and np.array_equal(kb2.flat_ann_term, kb.flat_ann_term))
    report = {"out": str(args.out), **sizes, "seconds": round(time.time() - t0, 1), "roundtrip_identical": same,
              "counts": kb.meta["counts"], "ontology_version": kb.meta.get("ontology_version"),
              "annotation_version": kb.meta.get("annotation_header", {}).get("version")}
    print(json.dumps(report, indent=1))
    if args.report:
        Path(args.report).write_text(json.dumps(report, indent=1) + "\n")
    return 0 if same else 1


if __name__ == "__main__":
    raise SystemExit(main())
