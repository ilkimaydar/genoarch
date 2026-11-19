#!/usr/bin/env python3
"""Write the synthetic mini HPO release used by the unit tests (hp.obo, hp.json, annotations).

Term ids are real HPO ids where possible but the content is a toy subset; do not use for analysis.
"""
import json
from pathlib import Path

HERE = Path(__file__).parent
VERSION = "2026-01-01"
# id, name, parents, exact synonyms, alt_ids
TERMS = [
    ("HP:0000001", "All", [], [], []),
    ("HP:0000005", "Mode of inheritance", ["HP:0000001"], [], []),
    ("HP:0000006", "Autosomal dominant inheritance", ["HP:0000005"], [], []),
    ("HP:0000118", "Phenotypic abnormality", ["HP:0000001"], [], []),
    ("HP:0000707", "Abnormality of the nervous system", ["HP:0000118"], [], []),
    ("HP:0001250", "Seizure", ["HP:0000707"], ["Seizures"], []),
    ("HP:0002373", "Febrile seizure", ["HP:0001250"], [], []),
    ("HP:0007359", "Focal-onset seizure", ["HP:0001250"], ["Focal seizure"], []),
    ("HP:0001249", "Intellectual disability", ["HP:0000707"], ["Mental retardation"], []),
    ("HP:0000478", "Abnormality of the eye", ["HP:0000118"], [], []),
    ("HP:0001083", "Ectopia lentis", ["HP:0000478"], ["Dislocated lens"], ["HP:0009999"]),
    ("HP:0000545", "Myopia", ["HP:0000478"], ["Nearsightedness"], []),
    ("HP:0033127", "Abnormality of the musculoskeletal system", ["HP:0000118"], [], []),
    ("HP:0001166", "Arachnodactyly", ["HP:0033127"], ["Long slender fingers"], []),
    ("HP:0001519", "Disproportionate tall stature", ["HP:0033127"], [], []),
    ("HP:0002616", "Aortic root aneurysm", ["HP:0000118"], [], []),
    ("HP:0004970", "Ascending tubular aorta aneurysm", ["HP:0002616", "HP:0033127"], [], []),
]
OBSOLETE = [("HP:0000999", "obsolete Old seizure term", "HP:0001250")]
QUALIFIED_IS_A = {"HP:0007359"}  # written with a trailing {xref=...} qualifier in hp.obo

# database_id, disease_name, qualifier, hpo_id, frequency, aspect
HPOA = [
    ("OMIM:100001", "Marfan-like A", "", "HP:0001166", "HP:0040281", "P"),
    ("OMIM:100001", "Marfan-like A", "", "HP:0001083", "12/20", "P"),
    ("OMIM:100001", "Marfan-like A", "", "HP:0000545", "", "P"),
    ("OMIM:100001", "Marfan-like A", "", "HP:0001519", "HP:0040280", "P"),
    ("OMIM:100001", "Marfan-like A", "", "HP:0004970", "1/2", "P"),
    ("OMIM:100001", "Marfan-like A", "", "HP:0000006", "", "I"),
    ("OMIM:100002", "Epilepsy B", "", "HP:0001250", "HP:0040282", "P"),
    ("OMIM:100002", "Epilepsy B", "", "HP:0002373", "3/4", "P"),
    ("OMIM:100002", "Epilepsy B", "", "HP:0001249", "", "P"),
    ("OMIM:100002", "Epilepsy B", "NOT", "HP:0001083", "", "P"),
    ("OMIM:100003", "Focal epilepsy C", "", "HP:0007359", "HP:0040283", "P"),
    ("OMIM:100003", "Focal epilepsy C", "", "HP:0000545", "", "P"),
    ("OMIM:100004", "Zero frequency E", "", "HP:0001250", "0/5", "P"),
    ("OMIM:100004", "Zero frequency E", "", "HP:0000545", "40%", "P"),
    ("ORPHA:200001", "Orphan Marfan", "", "HP:0001166", "", "P"),
    ("ORPHA:200001", "Orphan Marfan", "", "HP:0001083", "", "P"),
    ("DECIPHER:1", "CNV D", "", "HP:0001249", "", "P"),
    ("DECIPHER:1", "CNV D", "", "HP:0001519", "", "P"),
]
G2D = [
    (2200, "FBN1", "MENDELIAN", "OMIM:100001", "mim2gene"),
    (2200, "FBN1", "UNKNOWN", "ORPHA:200001", "orphadata"),
    (6323, "SCN1A", "MENDELIAN", "OMIM:100002", "mim2gene"),
    (1, "GENEC", "MENDELIAN", "OMIM:100003", "mim2gene"),
    (2, "GENED", "POLYGENIC", "OMIM:100003", "mim2gene"),
    (3, "NOPE", "MENDELIAN", "OMIM:999999", "mim2gene"),
]


def purl(t: str) -> str:
    return "http://purl.obolibrary.org/obo/" + t.replace(":", "_")


def main() -> None:
    names = {t[0]: t[1] for t in TERMS}
    lines = ["format-version: 1.2", f"data-version: hp/releases/{VERSION}", "ontology: hp", ""]
    for tid, name, parents, syns, alts in TERMS:
        lines += ["[Term]", f"id: {tid}", f"name: {name}"]
        lines += [f"alt_id: {a}" for a in alts]
        lines += [f'synonym: "{s}" EXACT []' for s in syns]
        for p in parents:
            q = ' {xref="PMID:1"}' if tid in QUALIFIED_IS_A else ""
            lines.append(f"is_a: {p}{q} ! {names[p]}")
        lines.append("")
    for tid, name, rep in OBSOLETE:
        lines += ["[Term]", f"id: {tid}", f"name: {name}", "is_obsolete: true", f"replaced_by: {rep}", ""]
    lines += ["[Typedef]", "id: part_of", "name: part of", ""]
    (HERE / "hp.obo").write_text("\n".join(lines))

    nodes, edges = [], []
    for tid, name, parents, syns, alts in TERMS:
        meta = {"synonyms": [{"pred": "hasExactSynonym", "val": s} for s in syns],
                "basicPropertyValues": [{"pred": "http://www.geneontology.org/formats/oboInOwl#hasAlternativeId", "val": a} for a in alts]}
        nodes.append({"id": purl(tid), "lbl": name, "type": "CLASS", "meta": meta})
        edges += [{"sub": purl(tid), "pred": "is_a", "obj": purl(p)} for p in parents]
    for tid, name, rep in OBSOLETE:
        nodes.append({"id": purl(tid), "lbl": name, "type": "CLASS", "meta": {
            "deprecated": True, "basicPropertyValues": [{"pred": "http://purl.obolibrary.org/obo/IAO_0100001", "val": purl(rep)}]}})
    graph = {"id": "http://purl.obolibrary.org/obo/hp.json",
             "meta": {"basicPropertyValues": [{"pred": "http://www.w3.org/2002/07/owl#versionInfo", "val": VERSION}]},
             "nodes": nodes, "edges": edges}
    (HERE / "hp.json").write_text(json.dumps({"graphs": [graph]}, indent=1))

    hpoa = ["#description: \"mini test annotations\"", f"#version: {VERSION}", f"#hpo-version: http://purl.obolibrary.org/obo/hp/releases/{VERSION}/hp.json",
            "database_id\tdisease_name\tqualifier\thpo_id\treference\tevidence\tonset\tfrequency\tsex\tmodifier\taspect\tbiocuration"]
    for did, dname, qual, hid, freq, aspect in HPOA:
        hpoa.append("\t".join([did, dname, qual, hid, "PMID:1", "PCS", "", freq, "", "", aspect, "HPO:test[2026-01-01]"]))
    (HERE / "phenotype.hpoa").write_text("\n".join(hpoa) + "\n")

    g2d = ["ncbi_gene_id\tgene_symbol\tassociation_type\tdisease_id\tsource"]
    g2d += [f"NCBIGene:{g}\t{s}\t{a}\t{d}\t{src}" for g, s, a, d, src in G2D]
    (HERE / "genes_to_disease.txt").write_text("\n".join(g2d) + "\n")

    g2p = ["ncbi_gene_id\tgene_symbol\thpo_id\thpo_name\tfrequency\tdisease_id"]
    for g, s, a, d, src in G2D:
        for did, dname, qual, hid, freq, aspect in HPOA:
            if did == d and qual != "NOT":
                g2p.append(f"{g}\t{s}\t{hid}\t{names.get(hid, '')}\t{freq or '-'}\t{d}")
    (HERE / "genes_to_phenotype.txt").write_text("\n".join(g2p) + "\n")


if __name__ == "__main__":
    main()
