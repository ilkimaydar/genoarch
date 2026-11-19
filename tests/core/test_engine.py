import math

import numpy as np
import pytest
from genoarch_core.engine import Engine, QueryError


def resnik(kb, a, b):
    return max(kb.ic[x] for x in kb.anc_sets[a] & kb.anc_sets[b])


def brute_scores(kb, present):
    s_pd = np.zeros(kb.n_diseases)
    s_dp = np.zeros(kb.n_diseases)
    for d in range(kb.n_diseases):
        ann = [int(t) for t in kb.ann_terms[d]]
        s_pd[d] = np.mean([max(resnik(kb, p, t) for t in ann) for p in present])
        s_dp[d] = np.mean([max(resnik(kb, t, p) for p in present) for t in ann])
    return s_pd, s_dp


@pytest.mark.parametrize("terms", [["HP:0001166"], ["HP:0001166", "HP:0001083", "HP:0000545"],
                                   ["HP:0002373", "HP:0001249"], ["HP:0007359", "HP:0004970", "HP:0001519"]])
def test_bma_matches_bruteforce(mini_engine, terms):
    kb = mini_engine.kb
    idx = sorted(kb.term_index[t] for t in terms)
    comp = mini_engine.components(idx, [], need_dp=True)
    s_pd, s_dp = brute_scores(kb, idx)
    assert np.allclose(comp["s_pd"], s_pd, atol=1e-12)
    assert np.allclose(comp["s_dp"], s_dp, atol=1e-12)
    sym = Engine.combine(comp, "symmetric", 0.0)
    assert np.allclose(sym, 0.5 * (s_pd + s_dp), atol=1e-12)
    bl = Engine.combine(comp, "blended", 0.0, 0.25)
    assert np.allclose(bl, 0.75 * s_pd + 0.25 * s_dp, atol=1e-12)


def test_prepare_conflicts_and_folding(mini_engine):
    q = mini_engine.prepare(["HP:0002373", "HP:0001250", "HP:0009999", "nonsense", "HP:0000006"], ["HP:0001250", "HP:0000545"])
    kb = mini_engine.kb
    ids = lambda xs: [kb.term_ids[i] for i in xs]  # noqa: E731
    assert ids(q["present_used"]) == ["HP:0001083", "HP:0002373"]       # Seizure folded into Febrile seizure
    assert ids(q["absent_used"]) == ["HP:0000545"]                      # Seizure absent conflicts with present
    assert {u["status"] for u in q["unresolved"]} == {"invalid", "not_phenotype"}
    assert any("mapped to its current id HP:0001083" in w for w in q["warnings"])
    assert any("implied by a present term" in w for w in q["warnings"])


def test_posterior_gene_mass_and_outputs(mini_engine):
    r = mini_engine.analyze(["HP:0001166", "HP:0001083"], [], top_k=50)
    post = [d["posterior"] for d in r["diseases"]]
    assert r["settings"]["n_universe"] == 6 and math.isclose(sum(post), 1.0, abs_tol=1e-12)
    assert post == sorted(post, reverse=True)
    assert r["diseases"][0]["id"] in ("OMIM:100001", "ORPHA:200001")
    assert math.isclose(sum(g["posterior"] for g in r["genes"]) + r["uncertainty"]["no_known_gene_mass"], 1.0, abs_tol=1e-12)
    assert r["genes"][0]["symbol"] == "FBN1"
    u = r["uncertainty"]
    assert 0 <= u["entropy_norm"] <= 1 and u["top1_margin"] >= 0 and u["top1_posterior"] == post[0]
    ev = r["evidence"][0]
    assert {row["term"]["id"] for row in ev["present"]} == {"HP:0001166", "HP:0001083"}
    assert all(row["ic"] >= 0 for row in ev["present"])
    assert "Research use only" in r["research_use_only"]


def test_absent_penalty_lowers_annotated_disease(mini_engine):
    base = mini_engine.analyze(["HP:0000545"], [], beta=0.0, null_samples=0, top_k=50)
    pen = mini_engine.analyze(["HP:0000545"], ["HP:0001083"], beta=1.0, null_samples=0, top_k=50)
    score = lambda r, d: next(x["score"] for x in r["diseases"] if x["id"] == d)  # noqa: E731
    assert score(pen, "OMIM:100001") < score(base, "OMIM:100001")
    assert score(pen, "OMIM:100003") == score(base, "OMIM:100003")  # not annotated with Ectopia lentis


def test_next_phenotypes_match_bruteforce_information_gain(mini_engine):
    r = mini_engine.analyze(["HP:0000545"], [], null_samples=0, top_k=50, n_next=20)
    s = mini_engine.settings()
    kb = mini_engine.kb
    post = {d["id"]: d["posterior"] for d in r["diseases"]}
    top = [d["id"] for d in r["diseases"]][: s["eig_top_k"]]
    z = sum(post[d] for d in top)

    def q(d, t):
        di = kb.disease_index[d]
        vals = [0.5 if np.isnan(f) else f for a, f in zip(kb.ann_terms[di], kb.ann_freqs[di]) if t in kb.anc_sets[int(a)]]
        return min(max(max(vals), s["eps"]), s["q_max"]) if vals else s["eps"]

    def h(p):
        return 0.0 if p <= 0 or p >= 1 else -(p * math.log2(p) + (1 - p) * math.log2(1 - p))

    assert r["next_phenotypes"]
    for item in r["next_phenotypes"]:
        t = kb.term_index[item["id"]]
        pyes = sum(post[d] / z * q(d, t) for d in top)
        eig = h(pyes) - sum(post[d] / z * h(q(d, t)) for d in top)
        assert item["eig_bits"] >= -1e-12 and math.isclose(item["eig_bits"], eig, abs_tol=1e-9)
        assert item["id"] not in ("HP:0000545", "HP:0000478", "HP:0000118")  # answered or implied
    eigs = [x["eig_bits"] for x in r["next_phenotypes"]]
    assert eigs == sorted(eigs, reverse=True)


def test_hypotheses_exclude_top_disease_genes(mini_engine):
    r = mini_engine.analyze(["HP:0007359", "HP:0000545"], [], top_k=50)
    top_genes = set(r["diseases"][0]["genes"])
    for h in r["hypotheses"]:
        assert h["symbol"] not in top_genes
        assert h["n_evidence"] >= 1 and all(e["ic"] >= 1.0 for e in h["evidence"])
        assert h["label"].startswith("Hypothesis")


def test_open_world_null_is_deterministic(mini_engine):
    a = mini_engine.analyze(["HP:0001166", "HP:0001083"], [], null_samples=40)
    b = mini_engine.analyze(["HP:0001166", "HP:0001083"], [], null_samples=40)
    ow = a["uncertainty"]["open_world"]
    assert ow == b["uncertainty"]["open_world"]
    assert 0 < ow["p_value"] <= 1 and ow["abstain"] == (ow["p_value"] > ow["alpha"]) and ow["samples"] == 40
    off = mini_engine.analyze(["HP:0001166"], [], null_samples=0)["uncertainty"]["open_world"]
    assert off["p_value"] is None and off["abstain"] is False


def test_exclusions_and_errors(mini_engine):
    r = mini_engine.analyze(["HP:0001166"], [], exclude_diseases=["OMIM:100001", "OMIM:404"], top_k=50)
    assert "OMIM:100001" not in [d["id"] for d in r["diseases"]] and r["settings"]["n_universe"] == 5
    assert any("OMIM:404" in w for w in r["query"]["warnings"])
    with pytest.raises(QueryError):
        mini_engine.analyze(["HP:0000006"], [])
    with pytest.raises(QueryError):
        mini_engine.analyze(["HP:0001166"], [], mode="bogus")
    with pytest.raises(QueryError):
        mini_engine.analyze(["HP:0001166"], [], T=0)


def test_search(mini_engine):
    assert mini_engine.search("seizure")[0]["id"] == "HP:0001250"
    hit = mini_engine.search("long slender")[0]
    assert hit["id"] == "HP:0001166" and hit["match"] == "synonym"
    assert mini_engine.search("HP:0009999")[0] == {**mini_engine.search("HP:0009999")[0], "id": "HP:0001083", "match": "id"}
    assert mini_engine.search("dislocated")[0]["id"] == "HP:0001083"
    assert mini_engine.search("") == [] and mini_engine.search("zzzz") == []


def test_real_index_marfan_like_patient(real_engine):
    r = real_engine.analyze(["HP:0001166", "HP:0001083", "HP:0002616", "HP:0001519", "HP:0000545"], [])
    assert "FBN1" in [g["symbol"] for g in r["genes"][:3]]
    assert r["uncertainty"]["open_world"]["abstain"] is False
    assert r["next_phenotypes"] and r["evidence"] and r["settings"]["calibrated"] in (True, False)
