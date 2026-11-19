"""Genoarch FastAPI service: JSON API plus the static web app.

Environment variables:
  GENOARCH_INDEX        compact index (default: <repo>/web/data/index.json.gz)
  GENOARCH_CALIBRATION  calibration JSON (default: <repo>/web/data/calibration.json)
  GENOARCH_WEB_DIR      static web app directory served at / (default: <repo>/web)
"""
from __future__ import annotations

import json
import os
from contextlib import asynccontextmanager
from functools import lru_cache
from pathlib import Path

from fastapi import FastAPI, HTTPException, Query
from fastapi.middleware.cors import CORSMiddleware
from fastapi.staticfiles import StaticFiles
from genoarch_core.engine import ENGINE_VERSION, RESEARCH_USE_ONLY, Engine, QueryError
from genoarch_core.index_io import read_index
from genoarch_core.method_card import method_card

from .models import (
    AnalyzeRequest,
    AnalyzeResponse,
    Health,
    HypothesesRequest,
    HypothesesResponse,
    NextRequest,
    NextResponse,
    RankResponse,
    SearchResponse,
)

REPO = Path(__file__).resolve().parents[2]


def _path(var: str, default: Path) -> Path:
    return Path(os.environ.get(var) or default)


@lru_cache(maxsize=1)
def get_engine() -> Engine:
    index = _path("GENOARCH_INDEX", REPO / "web" / "data" / "index.json.gz")
    cal_path = _path("GENOARCH_CALIBRATION", REPO / "web" / "data" / "calibration.json")
    calibration = json.loads(cal_path.read_text()) if cal_path.exists() else {}
    return Engine(read_index(index), calibration)


@asynccontextmanager
async def lifespan(app: FastAPI):
    get_engine()
    yield


app = FastAPI(
    title="Genoarch API",
    version=ENGINE_VERSION,
    description="HPO-based disease and gene ranking with calibrated uncertainty, next-best-phenotype suggestions and "
                "gene-association hypotheses. " + RESEARCH_USE_ONLY,
    lifespan=lifespan,
    docs_url="/docs",
    redoc_url="/redoc",
    openapi_url="/openapi.json",
)
app.add_middleware(CORSMiddleware, allow_origins=["*"], allow_methods=["GET", "POST"], allow_headers=["*"])


def _run(req: AnalyzeRequest, parts, **extra) -> dict:
    try:
        return get_engine().analyze(
            req.present, req.absent, mode=req.mode, universe=req.universe, beta=req.beta, w=req.w, T=req.T,
            null_samples=req.null_samples, alpha=req.alpha, exclude_diseases=req.exclude_diseases, top_k=req.top_k,
            parts=parts, **extra)
    except QueryError as exc:
        raise HTTPException(status_code=422, detail=str(exc)) from exc


@app.get("/health", response_model=Health, tags=["meta"])
def health() -> dict:
    eng = get_engine()
    meta = eng.kb.meta
    s = eng.settings()
    return {"status": "ok", "engine": f"genoarch-engine {eng.version}", "ontology_version": meta.get("ontology_version"),
            "annotation_version": (meta.get("annotation_header") or {}).get("version"),
            "counts": meta.get("counts") or {"terms": eng.kb.n_terms, "diseases": eng.kb.n_diseases, "genes": eng.kb.n_genes},
            "default_mode": s["mode"], "default_universe": s["universe"], "calibrated": s["calibrated"],
            "research_use_only": RESEARCH_USE_ONLY}


@app.get("/hpo/search", response_model=SearchResponse, tags=["hpo"])
def hpo_search(q: str = Query(..., min_length=1, max_length=200, description="label, synonym or HPO id"),
               limit: int = Query(20, ge=1, le=100)) -> dict:
    return {"query": q, "results": get_engine().search(q, limit)}


@app.post("/rank", response_model=RankResponse, tags=["analysis"])
def rank(req: AnalyzeRequest) -> dict:
    """Ranked diseases and genes with calibrated posteriors, uncertainty, open-world test and evidence."""
    return _run(req, ("rank",))


@app.post("/next-phenotypes", response_model=NextResponse, tags=["analysis"])
def next_phenotypes(req: NextRequest) -> dict:
    """Phenotypes to examine next, ranked by expected information gain over the current posterior."""
    return _run(req, ("next",), n_next=req.n)


@app.post("/hypotheses", response_model=HypothesesResponse, tags=["analysis"])
def hypotheses(req: HypothesesRequest) -> dict:
    """Genes outside the top-ranked diseases that share specific phenotypes with the patient (hypotheses only)."""
    return _run(req, ("hypotheses",), n_hyp=req.n)


@app.post("/analyze", response_model=AnalyzeResponse, tags=["analysis"])
def analyze(req: AnalyzeRequest) -> dict:
    """Everything above in one call (used by the web app in server mode)."""
    return _run(req, ("rank", "next", "hypotheses"))


@app.get("/method", tags=["meta"])
def method() -> dict:
    """Method card: definitions, fitted parameters, data versions, evaluation summary and limits."""
    return method_card(get_engine())


_web = _path("GENOARCH_WEB_DIR", REPO / "web")
if _web.is_dir():
    app.mount("/", StaticFiles(directory=str(_web), html=True), name="web")
