# Genoarch developer tasks. Raw data goes to data/ (git-ignored); built artifacts to dist/.
PY ?= python3
STORE_TAG ?= 0.1.27
STORE_DIR ?= data/phenopacket-store/$(STORE_TAG)/extracted/$(STORE_TAG)
EVAL_OUT ?= dist/eval
PORT ?= 8000

.PHONY: install data index calibrate parity static test lint serve docker-build docker-run clean

install:          ## install the package with API and dev extras
	$(PY) -m pip install -e ".[api,dev]"

data:             ## download pinned HPO release files and Phenopacket Store, verify sha256
	sh scripts/fetch_data.sh --with-store

index:            ## build web/data/index.json.gz from data/hpo
	$(PY) scripts/build_index.py --hpo-dir data/hpo --out web/data/index.json.gz

calibrate:        ## fit (beta, w, T) on the calibration slice and evaluate on the test slice
	$(PY) scripts/calibrate_evaluate.py --store $(STORE_DIR) --out $(EVAL_OUT) --calibration web/data/calibration.json

parity:           ## JS-vs-Python parity on the fixed patients written by `make calibrate`
	$(PY) scripts/parity_check.py --patients $(EVAL_OUT)/parity_patients.json --out $(EVAL_OUT)/parity

static:           ## method card, dist/site, dist/hf-space and the single-file dist/genoarch.html
	$(PY) scripts/build_static.py

test:
	$(PY) -m pytest -q

lint:
	ruff check core api scripts tests
	node --check web/js/engine.js
	node --check web/js/app.js

serve:            ## API + web app on http://127.0.0.1:$(PORT)
	PYTHONPATH=core:api uvicorn genoarch_api.main:app --reload --port $(PORT)

docker-build:
	docker build -t genoarch:0.1.0 .

docker-run:
	docker run --rm -p $(PORT):7860 genoarch:0.1.0

clean:
	rm -rf dist .pytest_cache .ruff_cache
