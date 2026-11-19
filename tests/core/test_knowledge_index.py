import math

import numpy as np
from genoarch_core.engine import Engine
from genoarch_core.index_io import load_index_doc, read_index, write_index


def test_knowledge_base_counts(mini_kb):
    kb = mini_kb
    assert kb.n_terms == 14
    assert kb.disease_ids == ["DECIPHER:1", "OMIM:100001", "OMIM:100002", "OMIM:100003", "OMIM:100004", "ORPHA:200001"]
    assert kb.gene_symbols == ["FBN1", "GENEC", "GENED", "SCN1A"]
    assert kb.meta["n_annotations_unmapped"] == 1          # inheritance row (aspect I)
    assert kb.meta["n_annotations_zero_frequency_dropped"] == 1  # Seizure 0/5
    assert kb.meta["n_gene_links_to_unannotated_diseases"] == 1


def test_information_content(mini_kb):
    kb = mini_kb
    n = kb.n_diseases
    root = kb.term_index["HP:0000118"]
    assert kb.ic[root] == 0.0
    arach = kb.term_index["HP:0001166"]
    assert kb.term_count[arach] == 2
    assert kb.ic[arach] == round(math.log(n / 2), 6)
    eye = kb.term_index["HP:0000478"]  # OMIM:100001, 100003, 100004, ORPHA:200001
    assert kb.term_count[eye] == 4


def test_frequencies_and_resolution(mini_kb):
    kb = mini_kb
    d = kb.disease_index["OMIM:100001"]
    freqs = dict(zip((kb.term_ids[t] for t in kb.ann_terms[d]), kb.ann_freqs[d]))
    assert freqs["HP:0001083"] == 0.6 and freqs["HP:0001519"] == 1.0 and np.isnan(freqs["HP:0000545"])
    assert kb.resolve("HP:0001083") == (kb.term_index["HP:0001083"], "ok")
    assert kb.resolve("HP:0009999") == (kb.term_index["HP:0001083"], "alt")
    assert kb.resolve("HP:0000999") == (kb.term_index["HP:0001250"], "alt")
    assert kb.resolve("HP:0000006") == (None, "not_phenotype")
    assert kb.resolve("HP:7777777") == (None, "unknown")
    assert kb.resolve("banana") == (None, "invalid")


def test_gene_profiles_equal_disease_union(mini_kb):
    check = mini_kb.gene_profile_union_check()
    assert check["available"] and check["equal"] == check["genes"] == 4


def test_index_roundtrip_is_exact(mini_kb, tmp_path):
    path = tmp_path / "i.json.gz"
    info = write_index(mini_kb, path)
    assert info["gzip_bytes"] > 0
    doc = load_index_doc(path)
    assert doc["format"] == "genoarch-index" and doc["version"] == 1
    kb2 = read_index(path)
    assert kb2.term_ids == mini_kb.term_ids and kb2.disease_ids == mini_kb.disease_ids
    assert np.array_equal(kb2.ic, mini_kb.ic)
    assert np.array_equal(kb2.flat_ann_freq_eff, mini_kb.flat_ann_freq_eff)
    assert [list(a) for a in kb2.gene_diseases] == [list(a) for a in mini_kb.gene_diseases]
    a = Engine(mini_kb).analyze(["HP:0001166", "HP:0001083"], ["HP:0001250"], null_samples=30)
    b = Engine(kb2).analyze(["HP:0001166", "HP:0001083"], ["HP:0001250"], null_samples=30)
    assert a == b
