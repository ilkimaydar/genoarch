"""Pydantic request / response models for the Genoarch API."""
from __future__ import annotations

from typing import Literal

from pydantic import BaseModel, ConfigDict, Field


class AnalyzeRequest(BaseModel):
    present: list[str] = Field(..., min_length=1, max_length=200, description="HPO ids observed in the patient",
                               examples=[["HP:0001166", "HP:0001083", "HP:0002616", "HP:0001519"]])
    absent: list[str] = Field(default_factory=list, max_length=500, description="HPO ids explicitly excluded")
    mode: Literal["asymmetric", "blended", "symmetric"] | None = Field(None, description="BMA variant; default from calibration")
    universe: Literal["omim", "all"] | None = Field(None, description="omim = OMIM + DECIPHER; all = + Orphanet")
    beta: float | None = Field(None, ge=0, le=10, description="absent-term weight (override of the fitted value)")
    w: float | None = Field(None, ge=0, le=1, description="disease->patient weight in blended mode (override)")
    T: float | None = Field(None, gt=0, le=100, description="softmax temperature (override of the fitted value)")
    null_samples: int | None = Field(None, ge=0, le=500, description="permutation null size (0 disables the open-world test)")
    alpha: float | None = Field(None, gt=0, lt=1, description="abstain threshold on the open-world p-value")
    exclude_diseases: list[str] = Field(default_factory=list, max_length=500, description="disease ids ruled out by the user")
    top_k: int = Field(20, ge=1, le=200)


class NextRequest(AnalyzeRequest):
    n: int = Field(10, ge=1, le=50)


class HypothesesRequest(AnalyzeRequest):
    n: int = Field(10, ge=1, le=50)


class Loose(BaseModel):
    model_config = ConfigDict(extra="allow")


class TermLabel(BaseModel):
    id: str
    name: str


class DiseaseLabel(BaseModel):
    id: str
    name: str


class SearchHit(TermLabel):
    synonym: str | None = None
    ic: float
    n_diseases: int
    match: Literal["id", "name", "synonym"]


class SearchResponse(BaseModel):
    query: str
    results: list[SearchHit]


class Unresolved(BaseModel):
    input: str
    status: str
    role: str


class QueryEcho(BaseModel):
    present: list[TermLabel]
    absent: list[TermLabel]
    present_used: list[TermLabel]
    absent_used: list[TermLabel]
    unresolved: list[Unresolved]
    warnings: list[str]


class Settings(BaseModel):
    mode: str
    universe: str
    T: float
    beta: float
    w: float
    alpha: float
    null_samples: int
    calibrated: bool
    n_universe: int
    excluded_diseases: list[str]


class DiseaseRow(DiseaseLabel):
    rank: int
    score: float
    posterior: float
    genes: list[str]
    z_null: float | None = None


class LinkedDisease(DiseaseLabel):
    posterior: float


class GeneRow(BaseModel):
    rank: int
    symbol: str
    ncbi: int
    posterior: float
    diseases: list[LinkedDisease]
    n_diseases: int


class OpenWorld(Loose):
    samples: int
    alpha: float
    p_value: float | None = None
    abstain: bool


class Uncertainty(BaseModel):
    top1_posterior: float
    top1_margin: float
    entropy_bits: float
    entropy_norm: float
    gene_top1_posterior: float
    no_known_gene_mass: float
    open_world: OpenWorld


class EvidencePresent(BaseModel):
    term: TermLabel
    match: TermLabel | None
    mica: TermLabel
    ic: float
    frequency: float | None


class EvidenceAbsent(BaseModel):
    term: TermLabel
    annotated: TermLabel
    frequency: float
    ic: float


class Evidence(BaseModel):
    disease: DiseaseLabel
    present: list[EvidencePresent]
    absent: list[EvidenceAbsent]


class NextPhenotype(TermLabel):
    eig_bits: float
    p_present: float
    ic: float
    supports: list[DiseaseLabel]
    against: list[DiseaseLabel]


class HypothesisEvidence(BaseModel):
    term: TermLabel
    mica: TermLabel
    ic: float
    via_disease: DiseaseLabel


class Hypothesis(BaseModel):
    rank: int
    symbol: str
    ncbi: int
    score: float
    n_evidence: int
    evidence: list[HypothesisEvidence]
    known_diseases: list[DiseaseLabel]
    label: str


class Base(BaseModel):
    engine: str
    query: QueryEcho
    settings: Settings
    research_use_only: str


class RankResponse(Base):
    diseases: list[DiseaseRow]
    genes: list[GeneRow]
    uncertainty: Uncertainty
    evidence: list[Evidence]


class NextResponse(Base):
    next_phenotypes: list[NextPhenotype]


class HypothesesResponse(Base):
    hypotheses: list[Hypothesis]


class AnalyzeResponse(RankResponse):
    next_phenotypes: list[NextPhenotype]
    hypotheses: list[Hypothesis]


class Health(BaseModel):
    status: str
    engine: str
    ontology_version: str | None
    annotation_version: str | None
    counts: dict
    default_mode: str
    default_universe: str
    calibrated: bool
    research_use_only: str
