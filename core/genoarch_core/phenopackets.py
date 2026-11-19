"""Minimal GA4GH Phenopacket (v2 JSON) reader for evaluation cases."""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path


@dataclass
class Case:
    id: str
    cohort: str
    present: list[str]
    absent: list[str]
    disease_id: str | None
    disease_label: str | None
    genes: list[str] = field(default_factory=list)


def read_phenopacket(path: str | Path) -> Case:
    path = Path(path)
    doc = json.loads(path.read_text(encoding="utf-8"))
    present, absent = [], []
    for pf in doc.get("phenotypicFeatures", []) or []:
        tid = ((pf.get("type") or {}).get("id") or "").strip()
        if not tid.startswith("HP:"):
            continue
        (absent if pf.get("excluded") else present).append(tid)
    disease_id = disease_label = None
    genes: list[str] = []
    for interp in doc.get("interpretations", []) or []:
        diag = interp.get("diagnosis") or {}
        dis = diag.get("disease") or {}
        if disease_id is None and dis.get("id"):
            disease_id, disease_label = dis["id"], dis.get("label")
        for gi in diag.get("genomicInterpretations", []) or []:
            vi = gi.get("variantInterpretation") or {}
            sym = (((vi.get("variationDescriptor") or {}).get("geneContext") or {}).get("symbol"))
            if not sym:
                sym = (gi.get("gene") or {}).get("symbol")
            if sym and sym not in genes:
                genes.append(sym)
    if disease_id is None:
        for dis in doc.get("diseases", []) or []:
            term = dis.get("term") or {}
            if term.get("id"):
                disease_id, disease_label = term["id"], term.get("label")
                break
    cohort = path.parent.name if path.parent.name != "phenopackets" else path.parent.parent.name
    return Case(id=str(doc.get("id") or path.stem), cohort=cohort, present=list(dict.fromkeys(present)),
                absent=list(dict.fromkeys(absent)), disease_id=disease_id, disease_label=disease_label, genes=genes)


def read_store(release_dir: str | Path) -> list[Case]:
    """Read every phenopacket under a Phenopacket Store release directory (sorted by path)."""
    root = Path(release_dir)
    files = sorted(set(root.glob("*/*.json")) | set(root.glob("*/phenopackets/*.json")))
    return [read_phenopacket(f) for f in files]
