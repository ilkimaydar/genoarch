# Genoarch

Genoarch ranks candidate diseases and genes for undiagnosed cases from a patient's phenotype terms (Human Phenotype
Ontology, HPO), states how confident the ranking is, flags phenotype sets that no known disease explains, suggests the next
phenotype worth examining, and lists gene-association hypotheses for follow-up.

**Research use only. Genoarch is not a medical device and must not be used for diagnosis, treatment or any clinical decision.**

## Live

* **Website: [ilkimaydar.github.io/genoarch](https://ilkimaydar.github.io/genoarch/)** (GitHub Pages, redeployed from
  `main` by `.github/workflows/pages.yml` on every push).
* In-browser build (Claude artifact, private to the owner until shared): https://claude.ai/artifact/5r1wD8Pn2ezzAmgU5BNbsS
* Hugging Face Space: [huggingface.co/spaces/ilkimaydar/genoarch](https://huggingface.co/spaces/ilkimaydar/genoarch)
  (static Space; the files of `dist/hf-space/` were uploaded by hand, see docs/DEPLOY.md; the workflow "Publish static
  build to the Hugging Face Space" stays optional and needs the repository secret `HF_TOKEN`).

## What it does

| Output | How |
| --- | --- |
| Ranked diseases | Resnik best-match average (BMA) between patient terms and HPO disease annotations; asymmetric, symmetric and a blended mode |
| Ranked genes | Gene probability = sum of linked disease probabilities, split across a disease's genes |
| Calibrated confidence | Softmax over diseases with a temperature fitted on a held-out slice of Phenopacket Store |
| Uncertainty | Top-1 probability, top-1 margin, entropy (bits, normalised), probability on diseases without a linked gene |
| Open-world flag | Permutation null (IC-matched random phenotype sets), max-T statistic; abstain when p > 0.05 |
| Next phenotypes | Expected information gain over the current disease posterior |
| Gene hypotheses | Genes outside the top diseases that share specific (high-IC) phenotypes with the patient, with the evidence terms |
| Evidence | Per disease: patient term, closest annotation, most informative shared ancestor, IC, annotated frequency |

The same engine exists twice: `core/genoarch_core/engine.py` (Python, numpy only) and `web/js/engine.js`
(browser). They produce the same ranking for the same input; see Parity below.

## Quickstart

Open the single-file build: `make static`, then open `dist/genoarch.html` in a browser. Nothing is uploaded; the
analysis runs in the page.

Python and API:

```sh
make install          # pip install -e ".[api,dev]"
make test             # pytest: core, API, JS parity on the mini fixture
make serve            # API + web app at http://127.0.0.1:8000 (OpenAPI docs at /docs)
```

```python
import json, sys
sys.path.insert(0, "core")
from genoarch_core.engine import Engine
from genoarch_core.index_io import read_index

engine = Engine(read_index("web/data/index.json.gz"), json.load(open("web/data/calibration.json")))
result = engine.analyze(["HP:0001166", "HP:0001083", "HP:0002616"], absent=[])
print(result["genes"][0]["symbol"], result["uncertainty"]["open_world"]["abstain"])
```

Docker (one image: FastAPI serving the API and `web/`):

```sh
docker compose up --build   # http://127.0.0.1:8000
```

Rebuild data, index and calibration from the public releases:

```sh
make data       # HPO v2026-09-01 and Phenopacket Store 0.1.27 from GitHub releases, sha256-verified
make index      # web/data/index.json.gz
make calibrate  # web/data/calibration.json + evaluation logs
make parity     # JS vs Python on the 20 fixed patients
make static     # dist/site, dist/hf-space, dist/genoarch.html
```

## API

| Method | Path | Purpose |
| --- | --- | --- |
| GET | `/health` | status, data versions, defaults |
| GET | `/hpo/search?q=` | term search by label, synonym or id |
| POST | `/rank` | diseases, genes, uncertainty, open-world test, evidence |
| POST | `/next-phenotypes` | next-best phenotypes by expected information gain |
| POST | `/hypotheses` | gene-association hypotheses |
| POST | `/analyze` | all of the above in one call |
| GET | `/method` | method card |

## Architecture

```
HPO release files ──> scripts/build_index.py ──> web/data/index.json.gz (1.45 MB gzip)
                                                     │
Phenopacket Store ──> scripts/calibrate_evaluate.py ─┼─> web/data/calibration.json (beta, w, T)
                                                     │
                 core/genoarch_core (Python) <───────┼───────> web/js/engine.js (browser)
                          │                                         │
                 api/genoarch_api (FastAPI)                 web/ (static app, works with no backend)
```

Details: [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md). Deployment, including Cloudflare in front of genoarch.ai:
[docs/DEPLOY.md](docs/DEPLOY.md). Data: [docs/DATA_CARD.md](docs/DATA_CARD.md). Model and evaluation:
[docs/MODEL_CARD.md](docs/MODEL_CARD.md).

## Data versions

* HPO ontology release 2026-09-01 and annotation release 2026-09-02 (hpo.jax.org; files from the
  `obophenotype/human-phenotype-ontology` GitHub release `v2026-09-01`, sha256 in docs/DATA_CARD.md)
* Index content: 19,178 phenotype terms (Phenotypic abnormality subtree), 12,865 diseases (8,463 OMIM,
  4,355 Orphanet, 47 DECIPHER), 5,268 genes, 264,562 disease-phenotype annotations
* Phenopacket Store 0.1.27 (10,377 phenopackets) for calibration and evaluation only; no case data is in the index

## Evaluation (test slice)

Disease-disjoint split of Phenopacket Store 0.1.27: calibration slice 3,619 cases; test slice 6,726 cases from
537 diseases. Default configuration (blended mode, OMIM + Orphanet + DECIPHER):

| Metric | Value |
| --- | --- |
| Gene top-1 / top-3 / top-10 | 0.658 / 0.737 / 0.818 |
| Disease top-1 / top-10 | 0.613 / 0.764 |
| Expected calibration error (15 bins), after / before temperature scaling | 0.119 / 0.609 |
| Open-world test: abstain rate on in-index cases / random phenotype sets (200 cases) | 0.245 / 0.960 |

All numbers come from the logged run `exp/platform_eval/run_20261118T224839Z` (copy of the summary in
[docs/evaluation/](docs/evaluation/)). The evaluation uses solved phenopacket cases (Phenopacket Store 0.1.27) as a
retrospective proxy for undiagnosed cases; no open cohort of truly undiagnosed patients is available to us, so the
metrics above describe solved, published cases. Expect lower accuracy on undiagnosed cohorts.

## Parity (Python vs JavaScript)

On 20 fixed patients (16 Phenopacket Store test cases and 4 edge cases), all 20 full analysis outputs match:
10,575 of 10,852 floating-point values are bit-identical, the largest difference is 8.9e-16 (posterior
exponentials), and every ranking-defining value (scores, null statistics) is bit-identical. Log:
`exp/platform_eval/parity_20261118T225909Z/parity_report.json`.

## Limits

* Phenotype-only ranking; no variant, inheritance or segregation input in this version.
* Accuracy and calibration are bounded by HPO annotation coverage and by how cases were phenotyped.
* Absent phenotypes are accepted and used for suggestions, but their ranking weight was fitted to 0 on the
  calibration slice (penalising them lowered accuracy there); users can override it.
* The open-world test rarely detects a missing disease when similar diseases exist (abstain 0.285 after removing
  all diseases of the causal gene vs 0.245 with them).
* OMIM and Orphanet entries for the same condition are not merged.

## License

Code: MIT ([LICENSE](LICENSE)). Data: this product uses the Human Phenotype Ontology (hpo.jax.org). HPO files
are used as downloaded (the index derived from them is described in docs/DATA_CARD.md) and their versions are displayed; cite the HPO as
the HPO project asks. OMIM and Orphanet content reaches Genoarch through the HPO annotation files; check their
own terms before any commercial use (not verified here).

## AI-use disclosure

The code, documentation and evaluation scripts were written with assistance from an AI model (Anthropic Claude)
under the direction of the project owner. Every reported number comes from a logged run.
