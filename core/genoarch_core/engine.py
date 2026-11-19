"""Genoarch analysis engine (numpy only).

Disease ranking by information-content (IC) weighted Resnik best-match average
(BMA), a temperature-scaled posterior, uncertainty outputs, an open-world
(abstain) test against an IC-matched permutation null, next-best-phenotype
suggestions by expected information gain, and novel gene-association
hypotheses. ``web/js/engine.js`` is a line-by-line port: every floating-point
reduction that feeds a ranking runs in the same order in both languages, so
scores and null statistics are bit-identical (verified by scripts/parity_check.py).

Definitions (patient present terms P, absent terms E, disease D with direct
annotations A(D) and propagated annotations P(D)):

* sim(p, D)  = max IC(a) over a in anc(p) and P(D)               (Resnik, best match)
* S_PD(D)    = mean_p sim(p, D)                                  (asymmetric, patient -> disease)
* S_DP(D)    = mean_{d in A(D)} max IC(a), a in anc(d) and anc(P)  (disease -> patient)
* S_sym(D)   = (S_PD + S_DP) / 2                                 (symmetric BMA)
* S_w(D)     = (1 - w) S_PD + w S_DP                             (blended; w fitted)
* pen(D)     = sum_e IC(e) * f_D(e), f_D(e) = max frequency of annotations of D at or below e
* score(D)   = S(D) - beta * pen(D) / |P|
* posterior  = softmax(score / T) over the selected disease universe
"""
from __future__ import annotations

import math
import re

import numpy as np

from .knowledge import KnowledgeBase
from .prng import Mulberry32, fnv1a32

ENGINE_VERSION = "0.1.0"
MODES = ("asymmetric", "blended", "symmetric")
UNIVERSES = {"omim": ("OMIM", "DECIPHER"), "all": ("OMIM", "DECIPHER", "ORPHA")}
RESEARCH_USE_ONLY = (
    "Research use only. Genoarch is not a medical device and must not be used for diagnosis, "
    "treatment or any clinical decision."
)
DEFAULTS = {
    "mode": "blended",
    "universe": "omim",
    "T": 1.0,
    "beta": 0.0,
    "w": 0.1,
    "alpha": 0.05,
    "null_samples": 100,
    "ic_bins": 10,
    "eps": 0.05,
    "q_max": 0.99,
    "eig_top_k": 50,
    "eig_min_ic": 1.0,
    "hyp_top_diseases": 10,
    "hyp_min_ic": 3.0,
    "hyp_min_terms": 2,
}
_ID_QUERY = re.compile(r"^(?:hp[:_]?\s*)?(\d{1,7})$")
_TOKEN = re.compile(r"[a-z0-9]+")


def _seq_sum(x: np.ndarray) -> float:
    """Left-to-right float64 sum (matches a JavaScript for-loop)."""
    return float(np.add.accumulate(x)[-1]) if len(x) else 0.0


def _qkey(x: float, scale: float) -> float:
    return math.floor(x * scale + 0.5)


def _h2(q: np.ndarray) -> np.ndarray:
    q = np.asarray(q, dtype=np.float64)
    out = np.zeros_like(q)
    m = (q > 0) & (q < 1)
    qm = q[m]
    out[m] = -(qm * np.log2(qm) + (1.0 - qm) * np.log2(1.0 - qm))
    return out


class QueryError(ValueError):
    pass


class Engine:
    def __init__(self, kb: KnowledgeBase, calibration: dict | None = None) -> None:
        self.kb = kb
        self.calibration = calibration or {}
        self.version = ENGINE_VERSION
        self.k_genes = np.asarray([len(g) for g in kb.disease_genes], dtype=np.int64)
        self._universes: dict[str, dict] = {}
        self._search_index = None

    # ------------------------------------------------------------------ settings
    def settings(self, mode=None, universe=None, beta=None, T=None, null_samples=None, alpha=None, w=None) -> dict:
        cal = self.calibration
        mode = mode or cal.get("default_mode") or DEFAULTS["mode"]
        universe = universe or cal.get("default_universe") or DEFAULTS["universe"]
        if mode not in MODES:
            raise QueryError(f"mode must be one of {MODES}")
        if universe not in UNIVERSES:
            raise QueryError(f"universe must be one of {tuple(UNIVERSES)}")
        fitted = (cal.get("params", {}).get(mode, {}) or {}).get(universe, {}) or {}
        s = dict(DEFAULTS)
        s.update({k: v for k, v in (cal.get("constants") or {}).items() if k in DEFAULTS})
        s["mode"], s["universe"] = mode, universe
        s["T"] = float(T if T is not None else fitted.get("T", DEFAULTS["T"]))
        s["beta"] = float(beta if beta is not None else fitted.get("beta", DEFAULTS["beta"]))
        if mode == "blended":
            s["w"] = float(w if w is not None else fitted.get("w", DEFAULTS["w"]))
        else:
            s["w"] = 0.5 if mode == "symmetric" else 0.0
        if not 0.0 <= s["w"] <= 1.0:
            raise QueryError("w must be in [0, 1]")
        if s["beta"] < 0:
            raise QueryError("beta must be >= 0")
        s["null_samples"] = int(null_samples if null_samples is not None else s["null_samples"])
        s["alpha"] = float(alpha if alpha is not None else s["alpha"])
        s["calibrated"] = bool(fitted) and T is None and beta is None and w is None
        if s["T"] <= 0:
            raise QueryError("T must be positive")
        if not 0 <= s["null_samples"] <= 1000:
            raise QueryError("null_samples must be in [0, 1000]")
        return s

    def universe(self, name: str, exclude=()) -> dict:
        """Disease universe (source prefixes), optionally minus excluded disease positions."""
        exclude = sorted(set(int(d) for d in exclude))
        key = name if not exclude else None
        hit = self._universes.get(name) if key else None
        if hit is not None:
            return hit
        kb = self.kb
        prefixes = set(UNIVERSES[name])
        mask = np.asarray([s in prefixes for s in kb.disease_source], dtype=bool)
        mask[exclude] = False
        uidx = np.nonzero(mask)[0].astype(np.int64)
        # IC-stratified pool of directly annotated terms for the permutation null
        pool = np.unique(kb.flat_ann_term[mask[kb.flat_ann_disease]]).astype(np.int64)
        pool = pool[np.lexsort((pool, kb.ic[pool]))]
        nb = max(1, min(DEFAULTS["ic_bins"], len(pool)))
        bounds = [(b * len(pool)) // nb for b in range(nb + 1)]
        bins = [pool[bounds[b]: bounds[b + 1]] for b in range(nb)]
        edges = [float(kb.ic[b[-1]]) for b in bins]
        # gene <- universe diseases (gene-major) and disease -> gene (disease-major)
        g_dis, g_off, g_ids = [], [], []
        for g in range(kb.n_genes):
            ds = [int(d) for d in kb.gene_diseases[g] if mask[int(d)]]
            if ds:
                g_off.append(len(g_dis))
                g_dis.extend(ds)
                g_ids.append(g)
        dg_dis, dg_gene = [], []
        for d in uidx:
            for g in sorted(kb.disease_genes[int(d)]):
                dg_dis.append(int(d))
                dg_gene.append(int(g))
        u = {
            "name": name, "mask": mask, "uidx": uidx, "bins": bins, "edges": edges,
            "g_dis": np.asarray(g_dis, dtype=np.int64), "g_off": np.asarray(g_off, dtype=np.int64),
            "g_ids": np.asarray(g_ids, dtype=np.int64),
            "dg_dis": np.asarray(dg_dis, dtype=np.int64), "dg_gene": np.asarray(dg_gene, dtype=np.int64),
            "dis_with_genes": np.asarray(sorted(set(dg_dis)), dtype=np.int64),
            "excluded": exclude,
        }
        if key:
            self._universes[name] = u
        return u

    # ------------------------------------------------------------------- query
    def prepare(self, present, absent=()) -> dict:
        kb = self.kb
        warnings: list[str] = []
        unresolved: list[dict] = []

        def resolve_all(items, role):
            out = []
            for raw in items or []:
                i, status = kb.resolve(str(raw))
                if i is None:
                    unresolved.append({"input": str(raw), "status": status, "role": role})
                    continue
                if status == "alt":
                    warnings.append(f"{raw} was mapped to its current id {kb.term_ids[i]}.")
                out.append(int(i))
            return sorted(set(out))

        pres = resolve_all(present, "present")
        absn = resolve_all(absent, "absent")
        implied_present: set[int] = set()
        for p in pres:
            implied_present |= kb.anc_sets[p]
        conflict = [e for e in absn if e in implied_present]
        for e in conflict:
            warnings.append(f"{kb.term_ids[e]} was marked absent but is implied by a present term; ignored as absent.")
        absn = [e for e in absn if e not in implied_present]
        pres_used = [p for p in pres if not any(q != p and p in kb.anc_sets[q] for q in pres)]
        abs_used = [e for e in absn if not any(f != e and f in kb.anc_sets[e] for f in absn)]
        if len(pres_used) < len(pres):
            warnings.append("Present terms that are ancestors of other present terms were folded into the more specific term.")
        return {"present": pres, "absent": absn, "present_used": pres_used, "absent_used": abs_used,
                "unresolved": unresolved, "warnings": warnings}

    # ---------------------------------------------------------------- scoring
    def bm_matrix(self, terms) -> np.ndarray:
        """Row r: best-match Resnik similarity of term r against every disease."""
        kb = self.kb
        out = np.zeros((len(terms), kb.n_diseases), dtype=np.float64)
        for r, p in enumerate(terms):
            row = out[r]
            for a in kb.ancestors_by_ic(int(p)):
                ica = kb.ic[a]
                if ica > 0:
                    row[kb.postings[a]] = ica
        return out

    def s_dp(self, terms) -> np.ndarray:
        kb = self.kb
        in_q = np.zeros(kb.n_terms, dtype=bool)
        for p in terms:
            in_q[kb.ancestors[int(p)]] = True
        vals = np.where(in_q[kb.anc_flat], kb.ic[kb.anc_flat], 0.0)
        m = np.maximum.reduceat(vals, kb.anc_off[:-1])
        sums = np.bincount(kb.flat_ann_disease, weights=m[kb.flat_ann_term], minlength=kb.n_diseases)
        return sums / kb.ann_len

    def absent_frequencies(self, e: int) -> np.ndarray:
        """f_D(e): max effective frequency of annotations of each disease at or below e."""
        kb = self.kb
        in_desc = np.zeros(kb.n_terms, dtype=bool)
        in_desc[kb.descendants(int(e))] = True
        vals = np.where(in_desc[kb.flat_ann_term], kb.flat_ann_freq_eff, 0.0)
        return np.maximum.reduceat(vals, kb.ann_off)

    def penalty(self, absent) -> np.ndarray:
        kb = self.kb
        pen = np.zeros(kb.n_diseases, dtype=np.float64)
        for e in absent:
            pen += kb.ic[int(e)] * self.absent_frequencies(e)
        return pen

    def components(self, present, absent=(), need_dp=True) -> dict:
        bm = self.bm_matrix(present)
        acc = np.zeros(self.kb.n_diseases, dtype=np.float64)
        for r in range(len(present)):
            acc += bm[r]
        n = max(len(present), 1)
        return {"bm": bm, "s_pd": acc / n, "s_dp": self.s_dp(present) if need_dp else None,
                "pen": self.penalty(absent) if len(absent) else None, "n": n}

    @staticmethod
    def combine(comp: dict, mode: str, beta: float, w: float = 0.0) -> np.ndarray:
        if mode == "symmetric":
            base = 0.5 * (comp["s_pd"] + comp["s_dp"])
        elif mode == "blended":
            base = (1.0 - w) * comp["s_pd"] + w * comp["s_dp"]
        else:
            base = comp["s_pd"]
        if beta and comp["pen"] is not None:
            return base - beta * comp["pen"] / comp["n"]
        return base

    def order(self, score: np.ndarray, uidx: np.ndarray) -> np.ndarray:
        """Universe diseases sorted by score descending, ties by index ascending."""
        s = score[uidx]
        return uidx[np.lexsort((uidx, -s))]

    def posterior(self, score: np.ndarray, uidx: np.ndarray, T: float) -> np.ndarray:
        s = score[uidx]
        w = np.exp((s - s.max()) / T)
        post = np.zeros(self.kb.n_diseases, dtype=np.float64)
        post[uidx] = w / _seq_sum(w)
        return post

    # ------------------------------------------------------------ null model
    def null_test(self, q: dict, s: dict, u: dict, comp: dict, score: np.ndarray) -> dict:
        """Open-world test (max-T permutation null).

        Each null patient replaces every present term by a random directly-annotated
        term from the same IC decile (seeded, reproducible). Per disease the observed
        score is standardised by the null mean and variance of that disease plus a
        pooled variance floor v0 (median of the per-disease null variances); each null
        patient is standardised leave-one-out, so it is exchangeable with the
        observation. Statistic: maximum standardised score over the universe;
        p = (1 + #{null max >= observed max}) / (R + 1). Abstain when p > alpha.
        """
        kb = self.kb
        R = s["null_samples"]
        uidx = u["uidx"]
        if R < 20:
            return {"samples": R, "p_value": None, "abstain": False, "alpha": s["alpha"], "z": None}
        seed_text = "genoarch-null|{}|{}|{:.4f}|{:.4f}|{}|{}".format(
            s["mode"], s["universe"], s["beta"], s["w"],
            ",".join(kb.term_ids[p] for p in q["present_used"]),
            ",".join(kb.term_ids[e] for e in q["absent_used"]),
        )
        rng = Mulberry32(fnv1a32(seed_text))
        edges = u["edges"]
        term_bins = []
        for p in q["present_used"]:
            icp = float(kb.ic[kb.ancestors[p]].max())  # best IC this term can reach
            term_bins.append(next((i for i, e in enumerate(edges) if icp <= e), len(edges) - 1))
        n = comp["n"]
        S = np.zeros((R, len(uidx)), dtype=np.float64)
        for r in range(R):
            terms = []
            for b in term_bins:
                pool = u["bins"][b]
                terms.append(int(pool[int(rng.next() * len(pool))]))
            bm = self.bm_matrix(terms)
            acc = np.zeros(kb.n_diseases, dtype=np.float64)
            for i in range(len(terms)):
                acc += bm[i]
            nc = {"s_pd": acc / n, "s_dp": self.s_dp(terms) if s["mode"] != "asymmetric" else None,
                  "pen": comp["pen"], "n": n}
            S[r] = self.combine(nc, s["mode"], s["beta"], s["w"])[uidx]
        s1 = np.add.accumulate(S, axis=0)[-1]
        s2 = np.add.accumulate(S * S, axis=0)[-1]
        mean = s1 / R
        var = np.maximum(s2 / R - mean * mean, 0.0)
        srt = np.sort(var)
        k = len(srt)
        v0 = float(srt[k // 2]) if k % 2 else float((srt[k // 2 - 1] + srt[k // 2]) / 2)
        z_obs = (score[uidx] - mean) / np.sqrt(var + v0 + 1e-12)
        m1 = (s1[None, :] - S) / (R - 1)
        v1 = np.maximum((s2[None, :] - S * S) / (R - 1) - m1 * m1, 0.0)
        z_null_max = ((S - m1) / np.sqrt(v1 + v0 + 1e-12)).max(axis=1)
        j = int(np.argmax(z_obs))
        z_max = float(z_obs[j])
        p_value = (1 + int((z_null_max >= z_max).sum())) / (R + 1)
        obs_raw = float(score[uidx].max())
        raw_max = S.max(axis=1)
        p_raw = (1 + int((raw_max >= obs_raw).sum())) / (R + 1)
        z_full = np.full(kb.n_diseases, np.nan)
        z_full[uidx] = z_obs
        return {
            "samples": R, "alpha": s["alpha"], "p_value": p_value, "abstain": p_value > s["alpha"],
            "statistic": "max standardized score (max-T, leave-one-out, pooled variance floor)", "z_max": z_max,
            "variance_floor": v0,
            "z_max_disease": kb.disease_label(int(uidx[j])),
            "null_z_max_mean": _seq_sum(z_null_max) / R, "null_z_max_max": float(z_null_max.max()),
            "p_value_raw_max": p_raw, "observed_max_score": obs_raw, "null_max_score_mean": _seq_sum(raw_max) / R,
            "_z": z_full,
        }

    # -------------------------------------------------------------- evidence
    def evidence(self, d: int, q: dict) -> dict:
        kb = self.kb
        prop = kb.disease_prop_set(d)
        ann = [int(t) for t in kb.ann_terms[d]]
        freqs = kb.ann_freqs[d]
        rows = []
        for p in q["present_used"]:
            mica = next(int(a) for a in kb.ancestors_by_ic_desc(p) if int(a) in prop)
            best = -1
            for j, t in enumerate(ann):
                if mica in kb.anc_sets[t] and (best < 0 or kb.ic[t] > kb.ic[ann[best]]):
                    best = j
            f = None if best < 0 or np.isnan(freqs[best]) else float(freqs[best])
            rows.append({"term": kb.term_label(p), "match": kb.term_label(ann[best]) if best >= 0 else None,
                         "mica": kb.term_label(mica), "ic": float(kb.ic[mica]), "frequency": f})
        absent_rows = []
        for e in q["absent_used"]:
            desc = set(int(x) for x in kb.descendants(e))
            best, bf = -1, 0.0
            for j, t in enumerate(ann):
                if t in desc:
                    fe = 0.5 if np.isnan(freqs[j]) else float(freqs[j])
                    if fe > bf:
                        best, bf = j, fe
            if best >= 0:
                absent_rows.append({"term": kb.term_label(e), "annotated": kb.term_label(ann[best]),
                                    "frequency": bf, "ic": float(kb.ic[e])})
        return {"disease": kb.disease_label(d), "present": rows, "absent": absent_rows}

    # -------------------------------------------------- next-best phenotypes
    def next_phenotypes(self, q: dict, s: dict, post: np.ndarray, order: np.ndarray, n_out: int = 10) -> list[dict]:
        kb = self.kb
        top = [int(d) for d in order[: s["eig_top_k"]]]
        if not top:
            return []
        pi = post[top]
        pi = pi / _seq_sum(pi)
        qmaps = []
        cand: set[int] = set()
        for d in top:
            qm: dict[int, float] = {}
            for t, f in zip(kb.ann_terms[d], kb.ann_freqs[d]):
                fe = 0.5 if np.isnan(f) else float(f)
                for a in kb.ancestors[int(t)]:
                    a = int(a)
                    if fe > qm.get(a, -1.0):
                        qm[a] = fe
            qmaps.append(qm)
            cand.update(int(t) for t in kb.ann_terms[d])  # curated (direct) annotations only
        blocked: set[int] = set()
        for p in q["present"]:
            blocked |= kb.anc_sets[p]
        for e in q["absent"]:
            blocked.update(int(x) for x in kb.descendants(e))
        eps, qmax = s["eps"], s["q_max"]
        terms = [t for t in sorted(cand) if t not in blocked and kb.ic[t] >= s["eig_min_ic"]]
        if not terms:
            return []
        col = {t: j for j, t in enumerate(terms)}
        Q = np.full((len(top), len(terms)), eps, dtype=np.float64)
        for i, qm in enumerate(qmaps):
            for t, v in qm.items():
                j = col.get(t)
                if j is not None:
                    Q[i, j] = min(max(v, eps), qmax)
        weighted = pi[:, None] * Q
        p_yes = np.add.accumulate(weighted, axis=0)[-1]
        cond = np.add.accumulate(pi[:, None] * _h2(Q), axis=0)[-1]
        eig = _h2(p_yes) - cond
        keys = [(-_qkey(float(eig[j]), 1e9), terms[j], j) for j in range(len(terms))]
        keys.sort()
        out = []
        for _, t, j in keys[:n_out]:
            col = Q[:, j]
            sup = [i for i in range(len(top)) if col[i] >= 0.5][:3]
            against = [i for i in range(len(top)) if col[i] == eps][:3]
            out.append({
                **kb.term_label(t), "eig_bits": float(eig[j]), "p_present": float(p_yes[j]), "ic": float(kb.ic[t]),
                "supports": [kb.disease_label(top[i]) for i in sup],
                "against": [kb.disease_label(top[i]) for i in against],
            })
        return out

    # ------------------------------------------------------------ hypotheses
    def hypotheses(self, q: dict, s: dict, u: dict, comp: dict, post: np.ndarray, order: np.ndarray, n_out: int = 10) -> list[dict]:
        kb = self.kb
        if len(u["g_ids"]) == 0 or not q["present_used"]:
            return []
        top = [int(d) for d in order[: s["hyp_top_diseases"]]]
        excluded: set[int] = set()
        for d in top:
            excluded.update(int(g) for g in kb.disease_genes[d])
        bm = comp["bm"]
        gbm = np.maximum.reduceat(bm[:, u["g_dis"]], u["g_off"], axis=1)
        hit = gbm >= s["hyp_min_ic"]
        masked = np.where(hit, gbm, 0.0)
        acc = np.zeros(gbm.shape[1], dtype=np.float64)
        for r in range(gbm.shape[0]):
            acc += masked[r]
        score = acc / comp["n"]
        n_ev = hit.sum(axis=0)
        need = min(s["hyp_min_terms"], len(q["present_used"]))
        gids = u["g_ids"]
        ok = np.asarray([int(g) not in excluded for g in gids]) & (n_ev >= need)
        cols = np.nonzero(ok)[0]
        cols = cols[np.lexsort((gids[cols], -n_ev[cols], -score[cols]))][:n_out]
        mask = u["mask"]
        out = []
        for rank, c in enumerate(cols, 1):
            g = int(gids[c])
            udis = [int(d) for d in kb.gene_diseases[g] if mask[int(d)]]
            ev = []
            for r, p in enumerate(q["present_used"]):
                if not hit[r, c]:
                    continue
                mica, via = None, None
                for a in kb.ancestors_by_ic_desc(p):
                    a = int(a)
                    via = next((d for d in udis if a in kb.disease_prop_set(d)), None)
                    if via is not None:
                        mica = a
                        break
                ev.append({"term": kb.term_label(p), "mica": kb.term_label(mica), "ic": float(kb.ic[mica]),
                           "via_disease": kb.disease_label(via)})
            known = sorted(udis, key=lambda d: (-post[d], d))[:3]
            out.append({
                "rank": rank, "symbol": kb.gene_symbols[g], "ncbi": kb.gene_ncbi[g], "score": float(score[c]),
                "n_evidence": int(n_ev[c]), "evidence": ev,
                "known_diseases": [kb.disease_label(d) for d in known],
                "label": "Hypothesis: shares specific phenotypes with this patient but is not linked to any of the "
                         f"top {len(top)} ranked diseases. Not an established association; research follow-up only.",
            })
        return out

    # --------------------------------------------------------------- analyze
    def resolve_diseases(self, ids) -> tuple[list[int], list[str]]:
        found, missing = [], []
        for raw in ids or []:
            d = self.kb.disease_index.get(str(raw).strip())
            (found.append(d) if d is not None else missing.append(str(raw)))
        return sorted(set(found)), missing

    def analyze(self, present, absent=(), *, mode=None, universe=None, beta=None, T=None, w=None, null_samples=None,
                alpha=None, exclude_diseases=(), top_k: int = 20, n_next: int = 10, n_hyp: int = 10,
                evidence_k: int = 10, parts=("rank", "next", "hypotheses")) -> dict:
        kb = self.kb
        s = self.settings(mode, universe, beta, T, null_samples, alpha, w)
        q = self.prepare(present, absent)
        if not q["present_used"]:
            raise QueryError("At least one valid present HPO term (Phenotypic abnormality subtree) is required.")
        excl, missing = self.resolve_diseases(exclude_diseases)
        for m in missing:
            q["warnings"].append(f"Excluded disease {m} is not in the knowledge base.")
        u = self.universe(s["universe"], excl)
        if len(u["uidx"]) < 2:
            raise QueryError("The disease universe is empty after exclusions.")
        comp = self.components(q["present_used"], q["absent_used"], need_dp=s["mode"] != "asymmetric")
        score = self.combine(comp, s["mode"], s["beta"], s["w"])
        order = self.order(score, u["uidx"])
        post = self.posterior(score, u["uidx"], s["T"])
        uidx = u["uidx"]
        pu = post[uidx]
        logs = np.zeros_like(pu)
        nz = pu > 0
        logs[nz] = pu[nz] * np.log2(pu[nz])
        entropy = -_seq_sum(logs)
        # gene posterior: P(g) = sum_D P(D) / k(D) over universe diseases linked to g
        w = post[u["dg_dis"]] / self.k_genes[u["dg_dis"]]
        gpost = np.bincount(u["dg_gene"], weights=w, minlength=kb.n_genes)
        gene_cand = np.unique(u["dg_gene"])
        gkeys = np.floor(gpost[gene_cand] * 1e12 + 0.5)
        gorder = gene_cand[np.lexsort((gene_cand, -gkeys))]
        no_gene = 1.0 - _seq_sum(post[u["dis_with_genes"]])

        result: dict = {
            "engine": f"genoarch-engine {self.version}",
            "query": {
                "present": [kb.term_label(p) for p in q["present"]],
                "absent": [kb.term_label(e) for e in q["absent"]],
                "present_used": [kb.term_label(p) for p in q["present_used"]],
                "absent_used": [kb.term_label(e) for e in q["absent_used"]],
                "unresolved": q["unresolved"], "warnings": q["warnings"],
            },
            "settings": {k: s[k] for k in ("mode", "universe", "T", "beta", "w", "alpha", "null_samples", "calibrated")}
            | {"n_universe": int(len(uidx)), "excluded_diseases": [kb.disease_ids[d] for d in excl]},
            "research_use_only": RESEARCH_USE_ONLY,
        }
        top = order[:top_k]
        if "rank" in parts:
            diseases = []
            for r, d in enumerate(top, 1):
                d = int(d)
                diseases.append({"rank": r, **kb.disease_label(d), "score": float(score[d]), "posterior": float(post[d]),
                                 "genes": [kb.gene_symbols[g] for g in sorted(kb.disease_genes[d])]})
            genes = []
            for r, g in enumerate(gorder[:top_k], 1):
                g = int(g)
                linked = sorted((int(d) for d in kb.gene_diseases[g] if u["mask"][int(d)]), key=lambda d: (-post[d], d))
                genes.append({"rank": r, "symbol": kb.gene_symbols[g], "ncbi": kb.gene_ncbi[g], "posterior": float(gpost[g]),
                              "diseases": [{**kb.disease_label(d), "posterior": float(post[d])} for d in linked[:3]],
                              "n_diseases": len(linked)})
            p1 = float(post[top[0]])
            p2 = float(post[top[1]]) if len(top) > 1 else 0.0
            null = self.null_test(q, s, u, comp, score)
            zf = null.pop("_z", None)
            for row, d in zip(diseases, top):
                zv = None if zf is None or np.isnan(zf[int(d)]) else float(zf[int(d)])
                row["z_null"] = zv
            result["diseases"] = diseases
            result["genes"] = genes
            result["uncertainty"] = {
                "top1_posterior": p1, "top1_margin": p1 - p2,
                "entropy_bits": entropy, "entropy_norm": entropy / math.log2(len(uidx)) if len(uidx) > 1 else 0.0,
                "gene_top1_posterior": float(gpost[gorder[0]]) if len(gorder) else 0.0,
                "no_known_gene_mass": no_gene, "open_world": null,
            }
            result["evidence"] = [self.evidence(int(d), q) for d in top[:evidence_k]]
        if "next" in parts:
            result["next_phenotypes"] = self.next_phenotypes(q, s, post, order, n_next)
        if "hypotheses" in parts:
            result["hypotheses"] = self.hypotheses(q, s, u, comp, post, order, n_hyp)
        return result

    # ----------------------------------------------------------------- search
    def _build_search(self) -> list:
        kb = self.kb
        idx = []
        for i in range(kb.n_terms):
            name = kb.names[i].lower()
            syns = [x.lower() for x in kb.synonyms[i]]
            idx.append((name, syns, _TOKEN.findall(name), [_TOKEN.findall(x) for x in syns]))
        return idx

    def search(self, query: str, limit: int = 20) -> list[dict]:
        kb = self.kb
        qs = (query or "").strip().lower()
        if not qs:
            return []
        out: list[tuple] = []
        m = _ID_QUERY.match(qs)
        if m and (qs.startswith("hp") or len(m.group(1)) == 7):
            i, status = kb.resolve("HP:" + m.group(1).zfill(7))
            if i is not None:
                out.append((-1, 0, 0, i, None))
        if self._search_index is None:
            self._search_index = self._build_search()
        qtok = _TOKEN.findall(qs)
        seen = {o[3] for o in out}
        for i, (name, syns, ntok, stoks) in enumerate(self._search_index):
            if i in seen:
                continue
            tier, syn = None, None
            if name == qs:
                tier = 0
            elif qs in syns:
                tier, syn = 1, kb.synonyms[i][syns.index(qs)]
            elif name.startswith(qs):
                tier = 2
            else:
                j = next((k for k, x in enumerate(syns) if x.startswith(qs)), -1)
                if j >= 0:
                    tier, syn = 3, kb.synonyms[i][j]
                elif qtok and all(any(w.startswith(t) for w in ntok) for t in qtok):
                    tier = 4
                else:
                    j = next((k for k, tk in enumerate(stoks) if qtok and all(any(w.startswith(t) for w in tk) for t in qtok)), -1)
                    if j >= 0:
                        tier, syn = 5, kb.synonyms[i][j]
                    elif qs in name:
                        tier = 6
                    else:
                        j = next((k for k, x in enumerate(syns) if qs in x), -1)
                        if j >= 0:
                            tier, syn = 7, kb.synonyms[i][j]
            if tier is not None:
                out.append((tier, -int(kb.term_count[i]), len(kb.names[i]), i, syn))
        out.sort()
        res = []
        for tier, _, _, i, syn in out[:limit]:
            res.append({**kb.term_label(i), "synonym": syn, "ic": float(kb.ic[i]), "n_diseases": int(kb.term_count[i]),
                        "match": "id" if tier < 0 else ("synonym" if syn else "name")})
        return res
