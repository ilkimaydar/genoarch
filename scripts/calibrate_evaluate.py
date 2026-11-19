#!/usr/bin/env python3
"""Fit (beta, T) and evaluate the Genoarch engine on Phenopacket Store.

Split: disease-disjoint. A case belongs to the calibration slice when
int(sha256(disease_id)[:8], 16) % 10 < 3, otherwise to the test slice.

Calibration slice: beta in {0, 0.02, 0.05, 0.1, 0.25, 0.5, 1} (and, for the blended
mode, w in {0.02, 0.05, 0.1, 0.2, 0.3}) chosen by disease MRR per (mode, universe); temperature T fitted by minimum NLL of the true disease
(48-point log grid 0.02..2 plus parabolic refinement in log T); default
(mode, universe) chosen by gene MRR, then disease MRR.
Test slice: ranking, calibration and uncertainty metrics for every
(mode, universe); open-world (abstain) and next-best-phenotype simulations on a
fixed 200-case subsample with the default configuration.

Outputs (all raw): <out>/summary.json, <out>/test_cases.jsonl,
<out>/open_world_cases.jsonl, <out>/next_phenotype_sim.jsonl,
<out>/parity_patients.json, and <calibration> (web/data/calibration.json).
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "core"))

from genoarch_core.engine import DEFAULTS, ENGINE_VERSION, MODES, UNIVERSES, Engine  # noqa: E402
from genoarch_core.index_io import read_index  # noqa: E402
from genoarch_core.phenopackets import read_store  # noqa: E402
from genoarch_core.prng import Mulberry32, fnv1a32  # noqa: E402

BETAS = [0.0, 0.02, 0.05, 0.1, 0.25, 0.5, 1.0]
WS = [0.02, 0.05, 0.1, 0.2, 0.3]


def w_grid(mode: str) -> list[float]:
    return WS if mode == "blended" else [0.5 if mode == "symmetric" else 0.0]
T_GRID = np.exp(np.linspace(math.log(0.02), math.log(2.0), 48))
N_SUBSAMPLE = 200


def h8(text: str) -> int:
    return int(hashlib.sha256(text.encode("utf-8")).hexdigest()[:8], 16)


def split_of(disease_id: str) -> str:
    return "calibration" if h8(disease_id) % 10 < 3 else "test"


def rank_of(score: np.ndarray, uidx: np.ndarray, truth: int) -> int:
    s, st = score[uidx], score[truth]
    return int((s > st).sum() + ((s == st) & (uidx < truth)).sum()) + 1


def log_softmax_truth(score, uidx, truth, T) -> float:
    s = score[uidx] / T
    m = s.max()
    return float(score[truth] / T - m - math.log(np.exp(s - m).sum()))


def gene_rank_and_post(eng: Engine, u: dict, post: np.ndarray, genes_truth: list[int]):
    kb = eng.kb
    w = post[u["dg_dis"]] / eng.k_genes[u["dg_dis"]]
    gpost = np.bincount(u["dg_gene"], weights=w, minlength=kb.n_genes)
    cand = np.unique(u["dg_gene"])
    keys = np.floor(gpost[cand] * 1e12 + 0.5)
    order = cand[np.lexsort((cand, -keys))]
    pos = {int(g): i + 1 for i, g in enumerate(order)}
    ranks = [pos[g] for g in genes_truth if g in pos]
    top = int(order[0])
    return (min(ranks) if ranks else None), float(gpost[top]), top


def ece(conf, correct, bins=15) -> float:
    conf, correct = np.asarray(conf), np.asarray(correct, dtype=float)
    edges = np.linspace(0, 1, bins + 1)
    total = 0.0
    for i in range(bins):
        m = (conf > edges[i]) & (conf <= edges[i + 1]) if i else (conf >= 0) & (conf <= edges[1])
        if m.any():
            total += m.sum() / len(conf) * abs(conf[m].mean() - correct[m].mean())
    return float(total)


def reliability(conf, correct, bins=10) -> list[dict]:
    conf, correct = np.asarray(conf), np.asarray(correct, dtype=float)
    edges = np.linspace(0, 1, bins + 1)
    out = []
    for i in range(bins):
        m = (conf > edges[i]) & (conf <= edges[i + 1]) if i else (conf >= 0) & (conf <= edges[1])
        out.append({"bin": [round(edges[i], 2), round(edges[i + 1], 2)], "n": int(m.sum()),
                    "mean_conf": float(conf[m].mean()) if m.any() else None,
                    "accuracy": float(correct[m].mean()) if m.any() else None})
    return out


def auroc(score, label) -> float | None:
    score, label = np.asarray(score, dtype=float), np.asarray(label, dtype=bool)
    pos, neg = score[label], score[~label]
    if not len(pos) or not len(neg):
        return None
    allv = np.concatenate([pos, neg])
    ranks = np.empty(len(allv))
    order = np.argsort(allv, kind="mergesort")
    sv = allv[order]
    i = 0
    while i < len(sv):
        j = i
        while j + 1 < len(sv) and sv[j + 1] == sv[i]:
            j += 1
        ranks[order[i:j + 1]] = (i + j) / 2 + 1
        i = j + 1
    return float((ranks[: len(pos)].sum() - len(pos) * (len(pos) + 1) / 2) / (len(pos) * len(neg)))


def topk(ranks, k) -> float:
    r = [x for x in ranks if x is not None]
    return float(np.mean([x <= k for x in r])) if r else float("nan")


def mrr(ranks) -> float:
    r = [x for x in ranks if x is not None]
    return float(np.mean([1.0 / x for x in r])) if r else float("nan")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--index", default=str(ROOT / "web" / "data" / "index.json.gz"))
    ap.add_argument("--store", required=True, help="Phenopacket Store release directory (extracted)")
    ap.add_argument("--store-version", default="0.1.27")
    ap.add_argument("--out", required=True, help="log directory, e.g. ../genoarch_proj/exp/platform_eval/run_<ts>")
    ap.add_argument("--calibration", default=str(ROOT / "web" / "data" / "calibration.json"))
    ap.add_argument("--limit", type=int, default=0, help="debug: limit number of cases")
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    t_start = time.time()
    log_lines: list[str] = []

    def log(msg: str) -> None:
        line = f"[{time.strftime('%H:%M:%S')}] {msg}"
        print(line, flush=True)
        log_lines.append(line)

    kb = read_index(args.index)
    eng = Engine(kb)
    cases = read_store(args.store)
    if args.limit:
        cases = sorted(cases, key=lambda c: h8(c.id))[: args.limit]
    log(f"loaded index ({kb.n_terms} terms, {kb.n_diseases} diseases, {kb.n_genes} genes) and {len(cases)} phenopackets")

    gene_index = {g: i for i, g in enumerate(kb.gene_symbols)}
    usable, skipped = [], {"no_disease_id": 0, "disease_not_in_index": 0, "no_valid_present_term": 0}
    for c in cases:
        if not c.disease_id:
            skipped["no_disease_id"] += 1
            continue
        d = kb.disease_index.get(c.disease_id)
        if d is None:
            skipped["disease_not_in_index"] += 1
            continue
        q = eng.prepare(c.present, c.absent)
        if not q["present_used"]:
            skipped["no_valid_present_term"] += 1
            continue
        usable.append((c, d, q, [gene_index[g] for g in c.genes if g in gene_index]))
    cal = [x for x in usable if split_of(x[0].disease_id) == "calibration"]
    test = [x for x in usable if split_of(x[0].disease_id) == "test"]
    log(f"usable {len(usable)} (calibration {len(cal)}, test {len(test)}); skipped {skipped}")
    universes = {name: eng.universe(name) for name in UNIVERSES}

    # ---------------------------------------------------------------- calibration
    cache = []
    for i, (c, d, q, g) in enumerate(cal):
        comp = eng.components(q["present_used"], q["absent_used"], need_dp=True)
        cache.append({"s_pd": comp["s_pd"], "s_dp": comp["s_dp"], "pen": comp["pen"], "n": comp["n"]})
        if i % 500 == 0:
            log(f"calibration components {i}/{len(cal)}")
    selection = {}
    for mode in MODES:
        for uname, u in universes.items():
            rows = []
            for w in w_grid(mode):
                for beta in BETAS:
                    ranks = []
                    for (c, d, q, g), comp in zip(cal, cache):
                        if not u["mask"][d]:
                            continue
                        ranks.append(rank_of(Engine.combine(comp, mode, beta, w), u["uidx"], d))
                    rows.append({"w": w, "beta": beta, "mrr": mrr(ranks), "top1": topk(ranks, 1), "top10": topk(ranks, 10),
                                 "n": len(ranks)})
            best = max(rows, key=lambda r: (r["mrr"], -r["beta"], -r["w"]))
            selection.setdefault(mode, {})[uname] = {"grid": rows, "beta": best["beta"], "w": best["w"]}
            log(f"selection {mode}/{uname}: best w={best['w']} beta={best['beta']} mrr={best['mrr']:.4f}; "
                f"beta=0 by w: {[(r['w'], round(r['mrr'], 4)) for r in rows if r['beta'] == 0.0]}")
    for mode in MODES:
        for uname, u in universes.items():
            beta, w = selection[mode][uname]["beta"], selection[mode][uname]["w"]
            lse = np.zeros(len(T_GRID))
            n = 0
            for (c, d, q, g), comp in zip(cal, cache):
                if not u["mask"][d]:
                    continue
                sc = Engine.combine(comp, mode, beta, w)
                s = sc[u["uidx"]]
                rel = s - s.max()
                st = sc[d] - s.max()
                z = np.outer(1.0 / T_GRID, rel)
                lse += np.log(np.exp(z).sum(axis=1)) - st / T_GRID
                n += 1
            nll = lse / n
            j = int(np.argmin(nll))
            T = float(T_GRID[j])
            if 0 < j < len(T_GRID) - 1:
                x = np.log(T_GRID[j - 1: j + 2])
                y = nll[j - 1: j + 2]
                den = (x[0] - x[1]) * (x[0] - x[2]) * (x[1] - x[2])
                a = (x[2] * (y[1] - y[0]) + x[1] * (y[0] - y[2]) + x[0] * (y[2] - y[1])) / den
                b = (x[2] ** 2 * (y[0] - y[1]) + x[1] ** 2 * (y[2] - y[0]) + x[0] ** 2 * (y[1] - y[2])) / den
                if a > 0:
                    T = float(math.exp(-b / (2 * a)))
            T = float(f"{T:.6g}")
            nll_fit = 0.0
            for (c, d, q, g), comp in zip(cal, cache):
                if u["mask"][d]:
                    nll_fit -= log_softmax_truth(Engine.combine(comp, mode, beta, w), u["uidx"], d, T)
            selection[mode][uname].update({"T": T, "nll_cal": nll_fit / n, "nll_grid_min": float(nll[j]), "n_cal": n})
            log(f"temperature {mode}/{uname}: T={T} NLL={nll_fit / n:.4f} (grid min {nll[j]:.4f} at {T_GRID[j]:.4f})")
    # default (mode, universe): gene MRR then disease MRR on the calibration slice at fitted (beta, T)
    for mode in MODES:
        for uname, u in universes.items():
            sel = selection[mode][uname]
            gr, dr = [], []
            for (c, d, q, g), comp in zip(cal, cache):
                if not u["mask"][d]:
                    continue
                sc = Engine.combine(comp, mode, sel["beta"], sel["w"])
                dr.append(rank_of(sc, u["uidx"], d))
                post = eng.posterior(sc, u["uidx"], sel["T"])
                gr.append(gene_rank_and_post(eng, u, post, g)[0] if g else None)
            sel.update({"cal_disease_mrr": mrr(dr), "cal_gene_mrr": mrr(gr), "cal_gene_top1": topk(gr, 1)})
    default = max(((m, un) for m in MODES for un in UNIVERSES),
                  key=lambda k: (round(selection[k[0]][k[1]]["cal_gene_mrr"], 6), round(selection[k[0]][k[1]]["cal_disease_mrr"], 6)))
    log(f"default configuration: mode={default[0]} universe={default[1]}")
    del cache

    # ----------------------------------------------------------------------- test
    per_case = []
    metrics_raw = {f"{m}/{un}": {"d_rank": [], "g_rank": [], "conf": [], "correct": [], "g_conf": [], "g_correct": [],
                                 "nll": [], "conf_T1": [], "entropy": [], "margin": [], "disease": []}
                   for m in MODES for un in UNIVERSES}
    for i, (c, d, q, g) in enumerate(test):
        comp = eng.components(q["present_used"], q["absent_used"], need_dp=True)
        row = {"id": c.id, "cohort": c.cohort, "disease": c.disease_id, "genes": c.genes,
               "n_present": len(q["present_used"]), "n_absent": len(q["absent_used"]), "configs": {}}
        for mode in MODES:
            for uname, u in universes.items():
                if not u["mask"][d]:
                    continue
                sel = selection[mode][uname]
                sc = Engine.combine(comp, mode, sel["beta"], sel["w"])
                order = eng.order(sc, u["uidx"])
                post = eng.posterior(sc, u["uidx"], sel["T"])
                post1 = eng.posterior(sc, u["uidx"], 1.0)
                dr = rank_of(sc, u["uidx"], d)
                pu = post[u["uidx"]]
                nz = pu > 0
                ent = -float((pu[nz] * np.log2(pu[nz])).sum())
                top1, top2 = int(order[0]), int(order[1])
                grank, gconf, gtop = gene_rank_and_post(eng, u, post, g) if g else (None, None, None)
                key = f"{mode}/{uname}"
                m = metrics_raw[key]
                m["d_rank"].append(dr)
                m["g_rank"].append(grank)
                m["conf"].append(float(post[top1]))
                m["correct"].append(top1 == d)
                m["conf_T1"].append(float(post1[top1]))
                m["nll"].append(-math.log(max(float(post[d]), 1e-300)))
                m["entropy"].append(ent)
                m["margin"].append(float(post[top1] - post[top2]))
                m["disease"].append(c.disease_id)
                if g:
                    m["g_conf"].append(gconf)
                    m["g_correct"].append(gtop in g)
                row["configs"][key] = {"disease_rank": dr, "gene_rank": grank, "top1_posterior": float(post[top1]),
                                       "truth_posterior": float(post[d]), "entropy_bits": ent,
                                       "margin": float(post[top1] - post[top2]), "top1": kb.disease_ids[top1]}
        per_case.append(row)
        if i % 1000 == 0:
            log(f"test {i}/{len(test)}")
    results = {}
    for key, m in metrics_raw.items():
        corr = np.asarray(m["correct"], dtype=bool)
        conf = np.asarray(m["conf"])
        by_dis: dict[str, list[bool]] = {}
        for dis, ok in zip(m["disease"], m["correct"]):
            by_dis.setdefault(dis, []).append(ok)
        sel_mask = conf >= 0.5
        results[key] = {
            "n_cases": len(m["d_rank"]), "n_cases_with_gene_in_index": len([x for x in m["g_rank"] if x is not None]),
            "n_diseases": len(by_dis),
            "disease_top1": topk(m["d_rank"], 1), "disease_top3": topk(m["d_rank"], 3), "disease_top10": topk(m["d_rank"], 10),
            "disease_mrr": mrr(m["d_rank"]), "disease_top1_macro": float(np.mean([np.mean(v) for v in by_dis.values()])),
            "gene_top1": topk(m["g_rank"], 1), "gene_top3": topk(m["g_rank"], 3), "gene_top10": topk(m["g_rank"], 10),
            "gene_mrr": mrr(m["g_rank"]),
            "nll": float(np.mean(m["nll"])), "ece_15": ece(conf, corr), "ece_15_T1": ece(m["conf_T1"], corr),
            "brier_top1": float(np.mean((conf - corr) ** 2)), "mean_top1_confidence": float(conf.mean()),
            "auroc_conf_vs_correct": auroc(conf, corr),
            "auroc_neg_entropy_vs_correct": auroc(-np.asarray(m["entropy"]), corr),
            "auroc_margin_vs_correct": auroc(m["margin"], corr),
            "coverage_conf_ge_0.5": float(sel_mask.mean()),
            "accuracy_conf_ge_0.5": float(corr[sel_mask].mean()) if sel_mask.any() else None,
            "gene_ece_15": ece(m["g_conf"], m["g_correct"]) if m["g_conf"] else None,
            "reliability_10": reliability(conf, corr),
        }
        log(f"test {key}: d@1={results[key]['disease_top1']:.4f} d@10={results[key]['disease_top10']:.4f} "
            f"g@1={results[key]['gene_top1']:.4f} g@10={results[key]['gene_top10']:.4f} ECE={results[key]['ece_15']:.4f} "
            f"(T=1: {results[key]['ece_15_T1']:.4f})")
    with open(out / "test_cases.jsonl", "w") as fh:
        for row in per_case:
            fh.write(json.dumps(row) + "\n")

    # ------------------------------------------- open-world and next-phenotype simulations
    calib = {
        "format": "genoarch-calibration", "version": 1, "engine": ENGINE_VERSION,
        "created_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "default_mode": default[0], "default_universe": default[1],
        "params": {m: {un: {"T": selection[m][un]["T"], "beta": selection[m][un]["beta"], "w": selection[m][un]["w"]}
                       for un in UNIVERSES} for m in MODES},
        "constants": {k: DEFAULTS[k] for k in ("alpha", "null_samples", "eps", "q_max", "eig_top_k", "eig_min_ic",
                                               "hyp_top_diseases", "hyp_min_ic", "hyp_min_terms", "ic_bins")},
    }
    eng = Engine(kb, calib)
    sub = sorted(test, key=lambda x: h8("subsample|" + x[0].id))[:N_SUBSAMPLE]
    ow_rows, sim_rows = [], []
    for i, (c, d, q, g) in enumerate(sub):
        r_in = eng.analyze(c.present, c.absent, parts=("rank",), top_k=1, evidence_k=0)
        gene_dis = sorted({int(x) for gi in g for x in kb.gene_diseases[gi]} | {d})
        r_out = eng.analyze(c.present, c.absent, exclude_diseases=[kb.disease_ids[x] for x in gene_dis],
                            parts=("rank",), top_k=1, evidence_k=0)
        u = eng.universe(default[1])
        rng = Mulberry32(fnv1a32("random-patient|" + c.id))
        rand_terms = []
        for p in q["present_used"]:
            icp = float(kb.ic[kb.ancestors[p]].max())
            b = next((k for k, e in enumerate(u["edges"]) if icp <= e), len(u["edges"]) - 1)
            pool = u["bins"][b]
            rand_terms.append(kb.term_ids[int(pool[int(rng.next() * len(pool))])])
        r_rand = eng.analyze(rand_terms, [], parts=("rank",), top_k=1, evidence_k=0)
        ow_rows.append({"id": c.id, "disease": c.disease_id, "n_present": len(q["present_used"]),
                        "in_kb": {k: r_in["uncertainty"]["open_world"][k] for k in ("p_value", "abstain", "z_max")},
                        "gene_removed": {k: r_out["uncertainty"]["open_world"][k] for k in ("p_value", "abstain", "z_max")},
                        "n_removed": len(gene_dis),
                        "random_patient": {k: r_rand["uncertainty"]["open_world"][k] for k in ("p_value", "abstain", "z_max")}})
        # next-best-phenotype simulation: start from half of the present terms, ask one question
        pres = q["present_used"]
        if len(pres) >= 2:
            half = [kb.term_ids[p] for p in pres[: max(1, len(pres) // 2)]]
            hidden_present = set()
            for p in q["present"]:
                hidden_present |= kb.anc_sets[p]
            excluded = set()
            for e in q["absent"]:
                excluded.update(int(x) for x in kb.descendants(e))
            base = eng.analyze(half, [], parts=("rank", "next"), top_k=1, evidence_k=0, null_samples=0, n_next=200)
            sc_base = rank_after(eng, half, [], d)

            def answer(term_id):
                t = kb.term_index[term_id]
                if t in hidden_present:
                    return "present"
                if t in excluded:
                    return "absent"
                return "unknown"

            def after(term_id, ans):
                if ans == "present":
                    return rank_after(eng, half + [term_id], [], d)
                if ans == "absent":
                    return rank_after(eng, half, [term_id], d)
                return sc_base

            cands = base["next_phenotypes"]
            if cands:
                top_t = cands[0]["id"]
                a_top = answer(top_t)
                rng2 = Mulberry32(fnv1a32("random-question|" + c.id))
                rnd_t = cands[int(rng2.next() * len(cands))]["id"]
                a_rnd = answer(rnd_t)
                sim_rows.append({"id": c.id, "rank_start": sc_base, "eig_term": top_t, "eig_answer": a_top,
                                 "eig_rank_after": after(top_t, a_top), "random_term": rnd_t, "random_answer": a_rnd,
                                 "random_rank_after": after(rnd_t, a_rnd), "n_candidates": len(cands)})
        if i % 25 == 0:
            log(f"open-world / next-phenotype simulation {i}/{len(sub)}")

    def rate(key):
        return float(np.mean([r[key]["abstain"] for r in ow_rows]))
    open_world = {"n_cases": len(ow_rows), "null_samples": DEFAULTS["null_samples"], "alpha": DEFAULTS["alpha"],
                  "abstain_rate_in_kb": rate("in_kb"), "abstain_rate_gene_removed": rate("gene_removed"),
                  "abstain_rate_random_patient": rate("random_patient"),
                  "median_p_in_kb": float(np.median([r["in_kb"]["p_value"] for r in ow_rows])),
                  "median_p_gene_removed": float(np.median([r["gene_removed"]["p_value"] for r in ow_rows])),
                  "auroc_p_gene_removed_vs_in_kb": auroc([r["gene_removed"]["p_value"] for r in ow_rows] + [r["in_kb"]["p_value"] for r in ow_rows],
                                                         [True] * len(ow_rows) + [False] * len(ow_rows))}
    log(f"open-world: {open_world}")

    def sim_summary(prefix):
        ans = [r[f"{prefix}_answer"] for r in sim_rows]
        known = [r for r in sim_rows if r[f"{prefix}_answer"] != "unknown"]
        delta = [math.log2(r["rank_start"]) - math.log2(r[f"{prefix}_rank_after"]) for r in known]
        return {"n": len(sim_rows), "answer_present": ans.count("present"), "answer_absent": ans.count("absent"),
                "answer_unknown": ans.count("unknown"),
                "mean_log2_rank_gain_when_answered": float(np.mean(delta)) if delta else None,
                "fraction_rank_improved_when_answered": float(np.mean([x > 0 for x in delta])) if delta else None}
    next_sim = {"protocol": "start from the first half (by term position) of the used present terms, ask the top suggestion; "
                            "oracle answers present if implied by the full present set, absent if implied by an excluded term, else unknown",
                "eig_top1": sim_summary("eig"), "random_candidate": sim_summary("random")}
    log(f"next-phenotype simulation: {next_sim}")
    with open(out / "open_world_cases.jsonl", "w") as fh:
        for r in ow_rows:
            fh.write(json.dumps(r) + "\n")
    with open(out / "next_phenotype_sim.jsonl", "w") as fh:
        for r in sim_rows:
            fh.write(json.dumps(r) + "\n")

    # ------------------------------------------------------------- parity patients
    parity = []
    for c, d, q, g in sorted(test, key=lambda x: h8("parity|" + x[0].id)):
        if len(q["present_used"]) >= 3 and len(parity) < 16:
            parity.append({"name": f"phenopacket:{c.id}", "present": c.present, "absent": c.absent, "options": {}})
    base_case = parity[0]
    alt_pairs = [(a, i) for a, i in kb.alt_map.items() if kb.ic[i] > 3]
    alt_id = sorted(alt_pairs)[0][0] if alt_pairs else None
    ex = next(c for c, d, q, g in test if len(q["present_used"]) >= 3)
    pid = kb.term_index[ex.present[0]] if ex.present[0] in kb.term_index else None
    parent = kb.term_ids[int(kb.parents[pid][0])] if pid is not None and len(kb.parents[pid]) else None
    parity.append({"name": "edge:alternative-id", "present": base_case["present"][:3] + ([alt_id] if alt_id else []),
                   "absent": [], "options": {}})
    parity.append({"name": "edge:present-absent-conflict", "present": ex.present[:4], "absent": ([parent] if parent else []) + ex.absent[:2],
                   "options": {}})
    single = max((t for t in base_case["present"] if t in kb.term_index), key=lambda t: kb.ic[kb.term_index[t]])
    parity.append({"name": "edge:single-term", "present": [single], "absent": [], "options": {}})
    many_abs = next(c for c, d, q, g in test if len(q["absent_used"]) >= 10 and len(q["present_used"]) >= 3)
    parity.append({"name": "edge:asymmetric-all-sources-many-absent", "present": many_abs.present, "absent": many_abs.absent,
                   "options": {"mode": "asymmetric", "universe": "all"}})
    (out / "parity_patients.json").write_text(json.dumps(parity, indent=1) + "\n")

    calib["fitted_on"] = {"dataset": f"Phenopacket Store {args.store_version}", "slice": "calibration",
                          "split_rule": "disease-disjoint: int(sha256(disease_id)[:8],16) % 10 < 3",
                          "n_cases": len(cal), "n_diseases": len({x[0].disease_id for x in cal})}
    calib["selection"] = selection
    dflt = results[f"{default[0]}/{default[1]}"]
    calib["evaluation"] = {
        "slice": "test", "n_cases": dflt["n_cases"], "n_diseases": dflt["n_diseases"],
        "default_config": {k: v for k, v in dflt.items() if k != "reliability_10"},
        "all_configs": {k: {kk: v[kk] for kk in ("disease_top1", "disease_top10", "disease_mrr", "gene_top1", "gene_top10",
                                                  "gene_mrr", "ece_15", "nll")} for k, v in results.items()},
        "open_world": open_world, "next_phenotype_simulation": next_sim,
        "log": "exp/platform_eval/" + out.name,
    }
    Path(args.calibration).write_text(json.dumps(calib, indent=1) + "\n")
    commit = subprocess.run(["git", "-C", str(ROOT), "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    summary = {"created_utc": calib["created_utc"], "seconds": round(time.time() - t_start, 1), "engine": ENGINE_VERSION,
               "git_commit_at_run": commit, "python": platform.python_version(), "numpy": np.__version__,
               "index": {"path": args.index, "sha256": hashlib.sha256(Path(args.index).read_bytes()).hexdigest(), "meta_counts": kb.meta.get("counts")},
               "store": {"path": args.store, "version": args.store_version, "n_phenopackets": len(cases)},
               "skipped": skipped, "n_usable": len(usable), "n_calibration": len(cal), "n_test": len(test),
               "selection": selection, "default": {"mode": default[0], "universe": default[1]},
               "test_results": results, "open_world": open_world, "next_phenotype_simulation": next_sim,
               "parity_patients": len(parity)}
    (out / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    (out / "run.log").write_text("\n".join(log_lines) + "\n")
    shutil.copy2(__file__, out / Path(__file__).name)
    log(f"done in {time.time() - t_start:.0f}s -> {out}")
    return 0


def rank_after(eng: Engine, present, absent, truth: int) -> int:
    s = eng.settings()
    q = eng.prepare(present, absent)
    u = eng.universe(s["universe"])
    comp = eng.components(q["present_used"], q["absent_used"], need_dp=s["mode"] != "asymmetric")
    return rank_of(Engine.combine(comp, s["mode"], s["beta"], s["w"]), u["uidx"], truth)


if __name__ == "__main__":
    raise SystemExit(main())
