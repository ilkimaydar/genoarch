import math

import pytest
from genoarch_core.annotations import combine_frequencies, parse_frequency
from genoarch_core.ontology import load_hpo_json, load_hpo_obo, normalize_hpo_id
from genoarch_core.prng import Mulberry32, fnv1a32


def test_obo_and_json_loaders_agree(mini_dir):
    a, b = load_hpo_obo(mini_dir / "hp.obo"), load_hpo_json(mini_dir / "hp.json")
    assert a.version == b.version == "2026-01-01"
    assert set(a.terms) == set(b.terms)
    for t in a.terms:
        assert sorted(a.terms[t].parents) == sorted(b.terms[t].parents)
        assert a.terms[t].synonyms == b.terms[t].synonyms
        assert a.terms[t].name == b.terms[t].name
    assert a.alt_map == b.alt_map == {"HP:0009999": "HP:0001083", "HP:0000999": "HP:0001250"}


def test_obo_is_a_with_qualifier_is_parsed(mini_dir):
    onto = load_hpo_obo(mini_dir / "hp.obo")
    assert onto.terms["HP:0007359"].parents == ["HP:0001250"]


def test_multi_parent_and_subtree(mini_dir):
    onto = load_hpo_obo(mini_dir / "hp.obo")
    assert sorted(onto.terms["HP:0004970"].parents) == ["HP:0002616", "HP:0033127"]
    sub = onto.subtree()
    assert "HP:0000006" not in sub and "HP:0004970" in sub and len(sub) == 14


@pytest.mark.parametrize("raw,expected", [
    ("HP:0001250", "HP:0001250"), ("hp:1250", "HP:0001250"), ("HP_0001250", "HP:0001250"),
    ("http://purl.obolibrary.org/obo/HP_0001250", "HP:0001250"), (" HP0001250 ", "HP:0001250"),
    ("OMIM:154700", None), ("", None), (None, None), ("HP:12345678", None),
])
def test_normalize_hpo_id(raw, expected):
    assert normalize_hpo_id(raw) == expected


def test_parse_frequency():
    assert parse_frequency("3/4") == ("fraction", 0.75, 3, 4)
    assert parse_frequency("40%") == ("percent", 0.4, 0, 0)
    assert parse_frequency("HP:0040281") == ("term", 0.895, 0, 0)
    assert parse_frequency("") is None and parse_frequency("-") is None and parse_frequency("1/0") is None


def test_combine_frequencies():
    assert combine_frequencies([("fraction", 0.5, 1, 2), ("fraction", 1.0, 3, 3)]) == 0.8
    assert combine_frequencies([("fraction", 0.5, 1, 2), ("term", 0.895, 0, 0)]) == 0.895
    assert combine_frequencies([None, None]) is None
    assert combine_frequencies([("fraction", 1 / 3, 1, 3)]) == 0.3333


def test_prng_matches_javascript_vectors():
    # vectors produced by web/js/engine.js under Node 22
    assert fnv1a32("") == 2166136261
    seed = fnv1a32("genoarch")
    assert seed == 3252211830
    r = Mulberry32(seed)
    assert [r.next() for _ in range(3)] == [0.7283493769355118, 0.3494753597769886, 0.9045961715746671]
    q = Mulberry32(1)
    assert [q.next(), q.next()] == [0.6270739405881613, 0.002735721180215478]
    assert all(0 <= Mulberry32(s).next() < 1 for s in range(100))
    assert not math.isnan(Mulberry32(0).next())
