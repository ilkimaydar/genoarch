"""Integer-indexed knowledge base used by the analysis engine.

The same structure is built either from the raw HPO release files
(:meth:`KnowledgeBase.from_hpo_files`) or from the compact browser index
written by ``scripts/build_index.py`` (:func:`genoarch_core.index_io.read_index`).

Definitions (D = disease, t = HPO term, N = number of annotated diseases):

* ``A(D)``  direct phenotype annotations of D (Phenotypic abnormality subtree)
* ``P(D)``  propagated annotations: every ancestor (including itself) of every t in A(D)
* ``IC(t) = ln(N / |{D : t in P(D)}|)``  information content from propagated annotations
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from .annotations import (
    combine_frequencies,
    read_genes_to_disease,
    read_genes_to_phenotype,
    read_phenotype_hpoa,
)
from .ontology import PHENOTYPIC_ABNORMALITY, load_ontology, normalize_hpo_id

UNKNOWN_FREQUENCY = 0.5
ASSOCIATION_PRIORITY = {"MENDELIAN": 0, "POLYGENIC": 1, "UNKNOWN": 2}
ASSOCIATION_CODES = ["MENDELIAN", "POLYGENIC", "UNKNOWN"]
# Fixed-point scales shared with the browser index and the JavaScript engine.
# Values are stored as integers k and decoded as k / SCALE in both languages, so the
# Python and JavaScript engines operate on bit-identical float64 inputs.
IC_SCALE = 1_000_000
FREQ_SCALE = 10_000


def quantize(values, scale: int) -> np.ndarray:
    """Round to a fixed-point grid: ``rint(v * scale) / scale`` (NaN preserved)."""
    v = np.asarray(values, dtype=np.float64)
    return np.rint(v * scale) / scale


class KnowledgeBase:
    """Terms, diseases and genes with precomputed propagation and IC."""

    def __init__(
        self,
        *,
        term_ids,
        names,
        synonyms,
        parents,
        alt_map,
        disease_ids,
        disease_names,
        ann_terms,
        ann_freqs,
        gene_symbols,
        gene_ncbi,
        gene_diseases,
        gene_assoc,
        gene_profiles=None,
        ic=None,
        term_count=None,
        other_terms=None,
        meta=None,
        calibration=None,
    ):
        self.term_ids: list[str] = list(term_ids)
        self.term_index: dict[str, int] = {t: i for i, t in enumerate(self.term_ids)}
        self.names: list[str] = list(names)
        self.synonyms: list[list[str]] = [list(s) for s in synonyms]
        self.parents = [np.asarray(sorted(set(int(x) for x in p)), dtype=np.int32) for p in parents]
        self.alt_map: dict[str, int] = dict(alt_map)
        self.other_terms: dict[str, str] = dict(other_terms or {})
        self.meta: dict = dict(meta or {})
        self.calibration: dict | None = calibration

        self.disease_ids: list[str] = list(disease_ids)
        self.disease_index = {d: i for i, d in enumerate(self.disease_ids)}
        self.disease_names: list[str] = list(disease_names)
        self.ann_terms: list[np.ndarray] = []
        self.ann_freqs: list[np.ndarray] = []
        for t, f in zip(ann_terms, ann_freqs):
            t = np.asarray(t, dtype=np.int32)
            f = quantize(f, FREQ_SCALE)
            order = np.argsort(t, kind="stable")
            self.ann_terms.append(t[order])
            self.ann_freqs.append(f[order])

        self.gene_symbols: list[str] = list(gene_symbols)
        self.gene_ncbi: list[int] = [int(x) for x in gene_ncbi]
        self.gene_diseases: list[np.ndarray] = []
        self.gene_assoc: list[list[str]] = []
        for dis, assoc in zip(gene_diseases, gene_assoc):
            pairs = sorted(zip((int(d) for d in dis), assoc))
            self.gene_diseases.append(np.asarray([p[0] for p in pairs], dtype=np.int32))
            self.gene_assoc.append([p[1] for p in pairs])
        self.gene_profiles = None if gene_profiles is None else [np.asarray(sorted(set(p)), dtype=np.int32) for p in gene_profiles]
        self._derive(ic, term_count)

    # ------------------------------------------------------------------ build
    def _derive(self, ic, term_count) -> None:
        n = self.n_terms = len(self.term_ids)
        self.n_diseases = len(self.disease_ids)
        self.n_genes = len(self.gene_symbols)
        self.root = self.term_index[PHENOTYPIC_ABNORMALITY]

        children: list[list[int]] = [[] for _ in range(n)]
        for i, ps in enumerate(self.parents):
            for p in ps:
                children[int(p)].append(i)
        self.children = [np.asarray(sorted(c), dtype=np.int32) for c in children]

        anc_sets: list[frozenset | None] = [None] * n

        def ancestors_of(i: int) -> frozenset:
            stack = [i]
            while stack:
                cur = stack[-1]
                if anc_sets[cur] is not None:
                    stack.pop()
                    continue
                pending = [int(p) for p in self.parents[cur] if anc_sets[int(p)] is None]
                if pending:
                    stack.extend(pending)
                    continue
                s = {cur}
                for p in self.parents[cur]:
                    s |= anc_sets[int(p)]
                anc_sets[cur] = frozenset(s)
                stack.pop()
            return anc_sets[i]

        for i in range(n):
            ancestors_of(i)
        self.anc_sets: list[frozenset] = anc_sets  # type: ignore[assignment]
        self.ancestors = [np.asarray(sorted(s), dtype=np.int32) for s in self.anc_sets]
        lens = np.asarray([len(a) for a in self.ancestors], dtype=np.int64)
        self.anc_off = np.concatenate([[0], np.cumsum(lens)]).astype(np.int64)
        self.anc_flat = np.concatenate(self.ancestors).astype(np.int32)

        top = set(int(c) for c in self.children[self.root])
        self.categories = [sorted(int(a) for a in s if int(a) in top) for s in self.anc_sets]

        # propagated disease annotation sets
        self.prop = [self._propagate(t) for t in self.ann_terms]
        nd = self.n_diseases
        flat = np.concatenate(self.prop) if nd else np.zeros(0, dtype=np.int32)
        if term_count is None:
            term_count = np.bincount(flat, minlength=n)
        self.term_count = np.asarray(term_count, dtype=np.int64)
        if ic is None:
            ic = np.zeros(n, dtype=np.float64)
            mask = self.term_count > 0
            ic[mask] = np.log(nd / self.term_count[mask])
            ic = quantize(ic, IC_SCALE)
        self.ic = np.asarray(ic, dtype=np.float64)
        self.postings = self._postings(self.prop, n)

        alen = np.asarray([len(t) for t in self.ann_terms], dtype=np.int64)
        self.ann_len = alen
        self.ann_off = np.concatenate([[0], np.cumsum(alen)[:-1]]).astype(np.int64) if nd else np.zeros(0, dtype=np.int64)
        self.disease_source = [d.split(":", 1)[0] for d in self.disease_ids]

        # flattened annotation table (disease, term, effective frequency)
        self.flat_ann_disease = np.repeat(np.arange(nd, dtype=np.int32), alen)
        self.flat_ann_term = np.concatenate(self.ann_terms) if nd else np.zeros(0, dtype=np.int32)
        f = np.concatenate(self.ann_freqs) if nd else np.zeros(0)
        self.flat_ann_freq_eff = np.where(np.isnan(f), UNKNOWN_FREQUENCY, f)

        # genes -> diseases (flattened for reduceat)
        glen = np.asarray([len(d) for d in self.gene_diseases], dtype=np.int64)
        self.gene_dis_flat = np.concatenate(self.gene_diseases).astype(np.int32) if self.n_genes else np.zeros(0, dtype=np.int32)
        self.gene_dis_off = np.concatenate([[0], np.cumsum(glen)[:-1]]).astype(np.int64) if self.n_genes else np.zeros(0, dtype=np.int64)
        self.disease_genes: list[list[int]] = [[] for _ in range(nd)]
        for g, dis in enumerate(self.gene_diseases):
            for d in dis:
                self.disease_genes[int(d)].append(g)
        if self.gene_profiles is not None:
            self.gene_prop = [self._propagate(p) for p in self.gene_profiles]
            self.gene_postings = self._postings(self.gene_prop, n)
        else:
            self.gene_prop = None
            self.gene_postings = None
        self.annotated_pool = np.unique(self.flat_ann_term)
        self.prop_sets: list[frozenset] | None = None
        self._desc_cache: dict[int, np.ndarray] = {}
        self._anc_by_ic_cache: dict[int, np.ndarray] = {}
        self._anc_by_ic_desc_cache: dict[int, np.ndarray] = {}

    def _propagate(self, terms) -> np.ndarray:
        s: set[int] = set()
        for t in terms:
            s |= self.anc_sets[int(t)]
        return np.asarray(sorted(s), dtype=np.int32)

    @staticmethod
    def _postings(prop: list[np.ndarray], n_terms: int) -> list[np.ndarray]:
        if not prop:
            return [np.zeros(0, dtype=np.int32) for _ in range(n_terms)]
        owner = np.repeat(np.arange(len(prop), dtype=np.int32), [len(p) for p in prop])
        flat = np.concatenate(prop)
        order = np.argsort(flat, kind="stable")
        flat_s, owner_s = flat[order], owner[order]
        bounds = np.searchsorted(flat_s, np.arange(n_terms + 1))
        return [owner_s[bounds[i]: bounds[i + 1]] for i in range(n_terms)]

    # --------------------------------------------------------------- helpers
    def resolve(self, raw: str) -> tuple[int | None, str]:
        """Resolve an HPO id. Returns ``(index, status)``.

        status: ``ok`` | ``alt`` (mapped from alternative/obsolete id) |
        ``not_phenotype`` (outside Phenotypic abnormality) | ``unknown`` | ``invalid``
        """
        tid = normalize_hpo_id(raw)
        if tid is None:
            return None, "invalid"
        if tid in self.term_index:
            return self.term_index[tid], "ok"
        if tid in self.alt_map:
            return self.alt_map[tid], "alt"
        if tid in self.other_terms:
            return None, "not_phenotype"
        return None, "unknown"

    def descendants(self, i: int) -> np.ndarray:
        """All descendants of term ``i`` including itself (sorted)."""
        hit = self._desc_cache.get(i)
        if hit is not None:
            return hit
        seen = {i}
        stack = [i]
        while stack:
            cur = stack.pop()
            for c in self.children[cur]:
                c = int(c)
                if c not in seen:
                    seen.add(c)
                    stack.append(c)
        arr = np.asarray(sorted(seen), dtype=np.int32)
        if len(self._desc_cache) < 4096:
            self._desc_cache[i] = arr
        return arr

    def ancestors_by_ic(self, i: int) -> np.ndarray:
        """Ancestors of ``i`` (incl. itself) ordered by IC ascending, ties by index."""
        hit = self._anc_by_ic_cache.get(i)
        if hit is not None:
            return hit
        anc = self.ancestors[i]
        order = np.lexsort((anc, self.ic[anc]))
        arr = anc[order]
        if len(self._anc_by_ic_cache) < 8192:
            self._anc_by_ic_cache[i] = arr
        return arr

    def ancestors_by_ic_desc(self, i: int) -> np.ndarray:
        """Ancestors of ``i`` (incl. itself) ordered by IC descending, ties by index ascending."""
        hit = self._anc_by_ic_desc_cache.get(i)
        if hit is not None:
            return hit
        anc = self.ancestors[i]
        arr = anc[np.lexsort((anc, -self.ic[anc]))]
        if len(self._anc_by_ic_desc_cache) < 8192:
            self._anc_by_ic_desc_cache[i] = arr
        return arr

    def disease_prop_set(self, d: int) -> frozenset:
        if self.prop_sets is None:
            self.prop_sets = [None] * self.n_diseases  # type: ignore[list-item]
        s = self.prop_sets[d]
        if s is None:
            s = frozenset(int(x) for x in self.prop[d])
            self.prop_sets[d] = s
        return s

    def term_id_str(self, i: int) -> str:
        return self.term_ids[i]

    def term_label(self, i: int) -> dict:
        return {"id": self.term_ids[i], "name": self.names[i]}

    def disease_label(self, d: int) -> dict:
        return {"id": self.disease_ids[d], "name": self.disease_names[d]}

    # ----------------------------------------------------------- constructors
    @classmethod
    def from_hpo_files(cls, hpo_dir: str | Path, *, ontology_file: str | None = None) -> "KnowledgeBase":
        """Build from an HPO release directory containing ``hp.json`` (or ``hp.obo``),
        ``phenotype.hpoa``, ``genes_to_disease.txt`` and ``genes_to_phenotype.txt``."""
        hpo_dir = Path(hpo_dir)
        if ontology_file is None:
            ontology_file = "hp.json" if (hpo_dir / "hp.json").exists() else "hp.obo"
        onto = load_ontology(hpo_dir / ontology_file)
        sub = onto.subtree(PHENOTYPIC_ABNORMALITY)
        term_ids = sorted(sub, key=lambda x: int(x[3:]))
        tindex = {t: i for i, t in enumerate(term_ids)}
        names = [onto.terms[t].name for t in term_ids]
        synonyms = [onto.terms[t].synonyms for t in term_ids]
        parents = [[tindex[p] for p in onto.terms[t].parents if p in tindex] for t in term_ids]
        alt_map = {a: tindex[cur] for a, cur in onto.alt_map.items() if cur in tindex}
        other_terms = {t: onto.terms[t].name for t in onto.terms if t not in tindex}
        for a, cur in onto.alt_map.items():
            if cur not in tindex and cur in onto.terms:
                other_terms.setdefault(a, onto.terms[cur].name)

        def map_term(raw: str) -> int | None:
            tid = normalize_hpo_id(raw)
            if tid in tindex:
                return tindex[tid]
            if tid in alt_map:
                return alt_map[tid]
            return None

        hpoa = read_phenotype_hpoa(hpo_dir / "phenotype.hpoa")
        diseases: dict[str, dict[int, float | None]] = {}
        unmapped = 0
        dropped_zero = 0
        for did, tmap in hpoa.positive.items():
            merged: dict[int, list] = {}
            for raw_t, obs in tmap.items():
                ti = map_term(raw_t)
                if ti is None:
                    unmapped += 1
                    continue
                merged.setdefault(ti, []).extend(obs)
            ann: dict[int, float | None] = {}
            for ti, obs in merged.items():
                f = combine_frequencies(obs)
                if f is not None and f <= 0.0:
                    dropped_zero += 1
                    continue
                ann[ti] = f
            if ann:
                diseases[did] = ann
        disease_ids = sorted(diseases)
        dindex = {d: i for i, d in enumerate(disease_ids)}
        disease_names = [hpoa.names.get(d, d) for d in disease_ids]
        ann_terms = [list(diseases[d].keys()) for d in disease_ids]
        ann_freqs = [[np.nan if v is None else v for v in diseases[d].values()] for d in disease_ids]

        links = read_genes_to_disease(hpo_dir / "genes_to_disease.txt")
        genes: dict[str, dict] = {}
        n_links_unannotated = 0
        for link in links:
            if not link.symbol:
                continue
            g = genes.setdefault(link.symbol, {"ncbi": link.ncbi_gene_id, "dis": {}})
            d = dindex.get(link.disease_id)
            if d is None:
                n_links_unannotated += 1
                continue
            prev = g["dis"].get(d)
            if prev is None or ASSOCIATION_PRIORITY.get(link.association_type, 9) < ASSOCIATION_PRIORITY.get(prev, 9):
                g["dis"][d] = link.association_type
        all_gene_symbols = set(genes)
        gene_symbols = sorted(s for s, g in genes.items() if g["dis"])

        # gene phenotype profiles from genes_to_phenotype.txt; rows that correspond to
        # NOT-qualified annotations in phenotype.hpoa are removed.
        g2p = read_genes_to_phenotype(hpo_dir / "genes_to_phenotype.txt")
        positive_pairs = {(d, t) for d, tmap in hpoa.positive.items() for t in tmap}
        profiles: dict[str, set[int]] = {}
        n_not_rows = 0
        for row in g2p:
            pair = (row.disease_id, row.hpo_id)
            if pair not in positive_pairs and pair in hpoa.negative:
                n_not_rows += 1
                continue
            ti = map_term(row.hpo_id)
            if ti is None:
                continue
            profiles.setdefault(row.symbol, set()).add(ti)
        gene_profiles = [sorted(profiles.get(s, set())) for s in gene_symbols]

        meta = {
            "source": "raw",
            "hpo_dir": str(hpo_dir),
            "ontology_version": onto.version,
            "ontology_file": ontology_file,
            "annotation_header": hpoa.header,
            "n_obsolete_terms": onto.n_obsolete,
            "n_terms_ontology": len(onto.terms),
            "n_hpoa_rows": hpoa.n_rows,
            "n_annotations_unmapped": unmapped,
            "n_annotations_zero_frequency_dropped": dropped_zero,
            "n_genes_total": len(all_gene_symbols),
            "n_gene_links_to_unannotated_diseases": n_links_unannotated,
            "n_g2p_rows_not_qualified_removed": n_not_rows,
        }
        return cls(
            term_ids=term_ids,
            names=names,
            synonyms=synonyms,
            parents=parents,
            alt_map=alt_map,
            disease_ids=disease_ids,
            disease_names=disease_names,
            ann_terms=ann_terms,
            ann_freqs=ann_freqs,
            gene_symbols=gene_symbols,
            gene_ncbi=[genes[s]["ncbi"] for s in gene_symbols],
            gene_diseases=[list(genes[s]["dis"].keys()) for s in gene_symbols],
            gene_assoc=[list(genes[s]["dis"].values()) for s in gene_symbols],
            gene_profiles=gene_profiles,
            other_terms=other_terms,
            meta=meta,
        )

    def gene_profile_union_check(self) -> dict:
        """Compare genes_to_phenotype profiles with the union of the gene's disease
        annotations. Equality justifies the compact browser index, which stores only
        gene -> disease links and derives gene profiles as unions."""
        if self.gene_prop is None:
            return {"available": False}
        equal = 0
        diffs = []
        for g in range(self.n_genes):
            union: set[int] = set()
            for d in self.gene_diseases[g]:
                union.update(int(x) for x in self.prop[int(d)])
            mine = set(int(x) for x in self.gene_prop[g])
            if mine == union:
                equal += 1
            elif len(diffs) < 10:
                diffs.append({"gene": self.gene_symbols[g], "only_profile": len(mine - union), "only_union": len(union - mine)})
        return {"available": True, "genes": self.n_genes, "equal": equal, "examples": diffs}
