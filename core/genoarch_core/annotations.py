"""Parsers for HPO annotation release files.

* ``phenotype.hpoa``          disease -> phenotype annotations with frequencies
* ``genes_to_disease.txt``    gene -> disease links (NCBI Gene, association type)
* ``genes_to_phenotype.txt``  gene -> phenotype profile (derived by HPO from the two above)
"""
from __future__ import annotations

import csv
import re
from dataclasses import dataclass, field
from pathlib import Path

# Midpoints of the HPO frequency classes (HP:0040279 subtree).
FREQUENCY_TERMS: dict[str, float] = {
    "HP:0040280": 1.0,    # Obligate (100%)
    "HP:0040281": 0.895,  # Very frequent (80-99%)
    "HP:0040282": 0.545,  # Frequent (30-79%)
    "HP:0040283": 0.17,   # Occasional (5-29%)
    "HP:0040284": 0.025,  # Very rare (1-4%)
    "HP:0040285": 0.0,    # Excluded (0%)
}
FREQ_DECIMALS = 4
_FRACTION = re.compile(r"^\s*(\d+)\s*/\s*(\d+)\s*$")
_PERCENT = re.compile(r"^\s*(\d+(?:\.\d+)?)\s*%\s*$")


def parse_frequency(raw: str | None) -> tuple[str, float, int, int] | None:
    """Parse one frequency cell.

    Returns ``(kind, value, numerator, denominator)`` where kind is
    ``fraction``, ``percent`` or ``term``; ``None`` when the cell is empty
    or not interpretable.
    """
    if raw is None:
        return None
    s = raw.strip()
    if not s or s == "-":
        return None
    m = _FRACTION.match(s)
    if m:
        num, den = int(m.group(1)), int(m.group(2))
        if den <= 0:
            return None
        return ("fraction", min(1.0, num / den), num, den)
    m = _PERCENT.match(s)
    if m:
        return ("percent", min(1.0, float(m.group(1)) / 100.0), 0, 0)
    if s in FREQUENCY_TERMS:
        return ("term", FREQUENCY_TERMS[s], 0, 0)
    return None


def combine_frequencies(observations: list[tuple[str, float, int, int] | None]) -> float | None:
    """Combine several frequency observations for one (disease, term) pair.

    * all observations are fractions -> pooled ``sum(num) / sum(den)``
    * otherwise                      -> maximum of the available values
    * no usable observation          -> ``None`` (unknown; 0.5 is used downstream)

    The result is rounded to 4 decimals so that every consumer (Python engine,
    browser index, JS engine) sees exactly the same number.
    """
    obs = [o for o in observations if o is not None]
    if not obs:
        return None
    if all(o[0] == "fraction" for o in obs):
        num = sum(o[2] for o in obs)
        den = sum(o[3] for o in obs)
        value = num / den if den else 0.0
    else:
        value = max(o[1] for o in obs)
    return round(min(1.0, max(0.0, value)), FREQ_DECIMALS)


@dataclass
class HpoaData:
    header: dict[str, str]
    names: dict[str, str]
    # disease -> term -> list of parsed frequency observations (positive rows only)
    positive: dict[str, dict[str, list]] = field(default_factory=dict)
    negative: set[tuple[str, str]] = field(default_factory=set)
    n_rows: int = 0


def read_phenotype_hpoa(path: str | Path) -> HpoaData:
    header: dict[str, str] = {}
    lines: list[str] = []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#"):
                key, _, val = line[1:].partition(":")
                header[key.strip()] = val.strip().strip('"')
            else:
                lines.append(line)
    data = HpoaData(header=header, names={})
    for row in csv.DictReader(lines, delimiter="\t"):
        data.n_rows += 1
        did = (row.get("database_id") or "").strip()
        hid = (row.get("hpo_id") or "").strip()
        if not did or not hid:
            continue
        data.names.setdefault(did, (row.get("disease_name") or did).strip())
        if (row.get("qualifier") or "").strip().upper() == "NOT":
            data.negative.add((did, hid))
            continue
        data.positive.setdefault(did, {}).setdefault(hid, []).append(parse_frequency(row.get("frequency")))
    return data


@dataclass
class GeneDiseaseLink:
    ncbi_gene_id: int
    symbol: str
    association_type: str
    disease_id: str
    source: str


def read_genes_to_disease(path: str | Path) -> list[GeneDiseaseLink]:
    out: list[GeneDiseaseLink] = []
    with open(path, encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            raw_id = (row.get("ncbi_gene_id") or "").strip()
            try:
                gid = int(raw_id.split(":")[-1])
            except ValueError:
                continue
            out.append(
                GeneDiseaseLink(
                    ncbi_gene_id=gid,
                    symbol=(row.get("gene_symbol") or "").strip(),
                    association_type=(row.get("association_type") or "UNKNOWN").strip() or "UNKNOWN",
                    disease_id=(row.get("disease_id") or "").strip(),
                    source=(row.get("source") or "").strip(),
                )
            )
    return out


@dataclass
class GenePhenotypeRow:
    ncbi_gene_id: int
    symbol: str
    hpo_id: str
    disease_id: str


def read_genes_to_phenotype(path: str | Path) -> list[GenePhenotypeRow]:
    out: list[GenePhenotypeRow] = []
    with open(path, encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            try:
                gid = int(str(row.get("ncbi_gene_id", "")).split(":")[-1])
            except ValueError:
                continue
            out.append(
                GenePhenotypeRow(
                    ncbi_gene_id=gid,
                    symbol=(row.get("gene_symbol") or "").strip(),
                    hpo_id=(row.get("hpo_id") or "").strip(),
                    disease_id=(row.get("disease_id") or "").strip(),
                )
            )
    return out
