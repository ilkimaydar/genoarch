"""Method card (served by GET /method and rendered on the web Method page)."""
from __future__ import annotations

from .engine import DEFAULTS, ENGINE_VERSION, RESEARCH_USE_ONLY, UNIVERSES, Engine

AI_USE = ("Genoarch code, documentation and evaluation scripts were written with assistance from an AI model "
          "(Anthropic Claude) under the direction of the project owner; all reported numbers come from logged runs.")


def method_card(engine: Engine) -> dict:
    kb, cal = engine.kb, engine.calibration or {}
    meta = kb.meta or {}
    header = meta.get("annotation_header") or {}
    return {
        "name": "Genoarch phenotype-driven disease and gene ranking",
        "engine_version": ENGINE_VERSION,
        "research_use_only": RESEARCH_USE_ONLY,
        "intended_use": "Research prioritisation of candidate diseases and genes for undiagnosed cases from a set of HPO phenotype terms, "
                        "with calibrated uncertainty, suggested next phenotypes to examine and gene-association hypotheses.",
        "not_intended_for": ["clinical diagnosis or treatment decisions", "screening of individuals", "use as a medical device",
                             "variant interpretation (no genotype input in this version)"],
        "inputs": "HPO term identifiers marked present or absent (Phenotypic abnormality subtree, HP:0000118). "
                  "Alternative and obsolete identifiers are mapped to current terms.",
        "method": {
            "information_content": "IC(t) = ln(N / n_t), N = annotated diseases in the index, n_t = diseases annotated to t or a descendant.",
            "similarity": "Resnik: sim(p, D) = max IC of a common ancestor of p and any annotation of D (best match).",
            "modes": {
                "asymmetric": "S_PD = mean over patient terms of the best match in the disease (patient -> disease BMA).",
                "symmetric": "(S_PD + S_DP) / 2, where S_DP = mean over disease annotations of the best match in the patient.",
                "blended": "(1 - w) S_PD + w S_DP with w fitted on the calibration slice; breaks ties of the asymmetric score.",
            },
            "absent_terms": "Penalty beta * sum_e IC(e) f_D(e) / |P| where f_D(e) is the highest annotated frequency of e or a "
                            "descendant in D (unknown frequency = 0.5). beta is fitted on the calibration slice.",
            "posterior": "softmax(score / T) over the selected disease universe; T fitted by minimum negative log-likelihood "
                         "of the true disease on the calibration slice (temperature scaling).",
            "gene_posterior": "P(g) = sum over linked diseases D of P(D) / k(D), k(D) = number of genes linked to D. "
                              "Mass on diseases without a linked gene is reported as no_known_gene_mass.",
            "uncertainty": "Top-1 posterior (calibrated confidence), top-1 margin, Shannon entropy (bits and normalised).",
            "open_world": "Permutation null: each present term is replaced by a random directly-annotated term from the same "
                          "IC decile (R null patients, seeded). Per-disease scores are standardised by the null mean and "
                          "variance plus a pooled variance floor (leave-one-out for null patients); statistic = maximum "
                          "standardised score; abstain when p > alpha.",
            "next_phenotypes": "Expected information gain (bits) about the disease among the top-K posterior diseases for each "
                               "curated annotation term not yet answered, under a Bernoulli model P(term | D) = annotated "
                               "frequency clamped to [eps, q_max], eps for diseases without the term.",
            "hypotheses": "Genes not linked to any of the top diseases whose known phenotype profile (union of their diseases' "
                          "annotations) matches at least hyp_min_terms patient terms with IC >= hyp_min_ic. Listed with the "
                          "matching terms as evidence and labelled as hypotheses.",
            "universes": {k: list(v) for k, v in UNIVERSES.items()},
        },
        "parameters": {
            "default_mode": cal.get("default_mode", DEFAULTS["mode"]),
            "default_universe": cal.get("default_universe", DEFAULTS["universe"]),
            "fitted": cal.get("params", {}),
            "constants": cal.get("constants", {k: DEFAULTS[k] for k in ("alpha", "null_samples", "eps", "q_max", "eig_top_k",
                                                                        "eig_min_ic", "hyp_top_diseases", "hyp_min_ic", "hyp_min_terms")}),
            "fitted_on": cal.get("fitted_on"),
        },
        "data": {
            "hpo_ontology_version": meta.get("ontology_version"),
            "hpo_annotation_version": header.get("version"),
            "annotation_description": header.get("description"),
            "counts": meta.get("counts") or {"terms": kb.n_terms, "diseases": kb.n_diseases, "genes": kb.n_genes},
            "inputs": meta.get("inputs"),
            "evaluation_dataset": "Phenopacket Store 0.1.27 (used only to fit and evaluate; not part of the index).",
        },
        "evaluation": cal.get("evaluation"),
        "limitations": [
            "Phenotype-only: no variant, segregation or inheritance input; ranks reflect HPO annotation similarity only.",
            "Knowledge is limited to diseases and gene links annotated in the HPO release; incomplete or biased "
            "annotations propagate to the ranking.",
            "The evaluation uses solved phenopacket cases (Phenopacket Store 0.1.27) as a retrospective proxy for undiagnosed cases; "
            "no open cohort of truly undiagnosed patients is available to us, so the reported metrics describe solved, published cases "
            "and performance on undiagnosed cohorts is expected to be lower.",
            "Calibration was fitted on one disease-disjoint slice of Phenopacket Store; confidence values may be miscalibrated "
            "for other populations or phenotyping styles.",
            "The open-world test detects phenotype sets not better matched than random sets of similar specificity; "
            "it is not a reliable detector of novel diseases when similar known diseases exist.",
            "Gene hypotheses are generated from phenotype overlap only and are not evidence of causality.",
            "OMIM and Orphanet entries for the same condition are not merged in the 'all' universe.",
        ],
        "license": "Code: MIT (see LICENSE). Data: HPO release files (hpo.jax.org) used unchanged; see docs/DATA_CARD.md.",
        "ai_use_disclosure": AI_USE,
    }
