# Genoarch architecture

## Modules

| Path | Role |
| --- | --- |
| `core/genoarch_core/ontology.py` | HPO loader (`hp.json` OBO Graphs or `hp.obo`), alternative / obsolete id map |
| `core/genoarch_core/annotations.py` | Parsers for `phenotype.hpoa`, `genes_to_disease.txt`, `genes_to_phenotype.txt`; frequency parsing and pooling |
| `core/genoarch_core/knowledge.py` | Integer-indexed knowledge base: ancestors, propagated annotations, postings, IC, gene links |
| `core/genoarch_core/index_io.py` | Compact browser index (gzip JSON) writer / reader |
| `core/genoarch_core/engine.py` | Ranking, posterior, uncertainty, open-world test, next phenotypes, hypotheses, search |
| `core/genoarch_core/prng.py` | FNV-1a and mulberry32, bit-identical with the JavaScript port |
| `core/genoarch_core/method_card.py` | Method card served by `GET /method` and rendered on the Method page |
| `core/genoarch_core/phenopackets.py` | Phenopacket v2 reader for evaluation |
| `api/genoarch_api/` | FastAPI app and Pydantic models; mounts `web/` at `/` |
| `web/js/engine.js` | JavaScript port of `engine.py` (browser and Node) |
| `web/js/app.js`, `web/index.html`, `web/css/app.css` | Static web app (Home, Analyze, Method, Pricing, About, Legal) |
| `scripts/` | Index build, calibration and evaluation, parity check, static build, data download |

## Data flow

```
data/hpo/*.{json,hpoa,txt}
   │  KnowledgeBase.from_hpo_files: Phenotypic abnormality subtree, frequency pooling,
   │  zero-frequency rows dropped, NOT rows excluded, IC = ln(N / n_t) on a 1e-6 grid
   ▼
web/data/index.json.gz  ── read_index ──> KnowledgeBase ──> Engine (Python, API)
   │                                                        
   └── decodeIndexBytes ──> buildKB ──> Engine (JavaScript, browser)

Phenopacket Store ──> calibrate_evaluate.py ──> web/data/calibration.json
                                          └──> exp/platform_eval/run_*/ (raw logs)
```

## Index format (`genoarch-index` v1)

Gzip-compressed JSON (1,453,542 bytes for HPO 2026-09-01; 5,916,188 bytes uncompressed). Terms, diseases and
genes are referenced by position; sorted integer lists are delta-encoded. IC is stored as integers x 1e6 and
frequencies as integers x 1e4 (-1 = unknown), so both engines decode identical float64 values. The stored
ancestor closure is checked against the `is_a` closure on load. Gene phenotype profiles are not stored; they are
the union of the propagated annotations of each gene's diseases.

## Scoring pipeline (per query)

1. Resolve ids (alternative and obsolete ids map to current terms), drop absent terms implied by present
   terms, fold present ancestors into their most specific present descendant.
2. For each present term p, a best-match vector over diseases: scatter IC(a) over the postings of the
   ancestors a of p in increasing IC order (the last write is the maximum).
3. S_PD = sequential sum of best-match vectors / |P|; S_DP from the max IC of shared ancestors of each
   disease annotation; blended score (1 - w) S_PD + w S_DP; optional absent penalty.
4. Posterior = softmax(score / T) over the selected universe; gene posterior by splitting disease mass
   equally over linked genes.
5. Open-world test: R seeded null patients, per-disease standardisation with a pooled variance floor,
   leave-one-out for null patients, max-T p-value.
6. Next phenotypes: Bernoulli information gain over the top 50 diseases for curated annotation terms.
7. Hypotheses: per-gene best matches (max over the gene's diseases), genes outside the top 10 diseases with
   at least two matches of IC >= 3.

## Python / JavaScript parity

Every reduction that feeds a ranking runs in the same order in both languages (sequential sums, the same
max scans, the same tie-breaks by position). Sorting keys for posteriors and information gain are quantised
(1e-12, 1e-9) to make near-ties deterministic across `exp` implementations. The PRNG for the null model is
32-bit FNV-1a + mulberry32 in both languages, seeded from the canonical query string. `scripts/parity_check.py`
runs both engines and compares full outputs; `tests/js/test_js_parity.py` repeats the check on the mini fixture
in CI.

## Web app

Static HTML, CSS and JavaScript with hash routes (`#home`, `#analyze`, `#method`, `#pricing`, `#about`,
`#legal`). On load it reads the index embedded in the page (single-file build) or `data/index.json.gz`, the
calibration and the method card. When served by the FastAPI app it offers an optional server mode
(`POST /analyze`); otherwise all computation stays in the browser.
