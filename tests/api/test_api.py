import importlib

import pytest

fastapi = pytest.importorskip("fastapi")
from fastapi.testclient import TestClient  # noqa: E402


@pytest.fixture(scope="module")
def client(mini_index, tmp_path_factory):
    import os

    web = tmp_path_factory.mktemp("web")
    (web / "index.html").write_text("<!doctype html><title>Genoarch</title><p>Research use only</p>")
    os.environ["GENOARCH_INDEX"] = str(mini_index)
    os.environ["GENOARCH_CALIBRATION"] = str(mini_index.parent / "calibration.json")
    os.environ["GENOARCH_WEB_DIR"] = str(web)
    import genoarch_api.main as main

    main = importlib.reload(main)
    main.get_engine.cache_clear()
    with TestClient(main.app) as c:
        yield c


def test_health(client):
    r = client.get("/health")
    assert r.status_code == 200
    body = r.json()
    assert body["status"] == "ok" and body["counts"]["diseases"] == 6 and body["default_universe"] == "all"


def test_search(client):
    r = client.get("/hpo/search", params={"q": "arachno"})
    assert r.status_code == 200 and r.json()["results"][0]["id"] == "HP:0001166"
    assert client.get("/hpo/search").status_code == 422


def test_rank_next_hypotheses_analyze(client):
    body = {"present": ["HP:0001166", "HP:0001083"], "absent": ["HP:0001250"], "null_samples": 30}
    r = client.post("/rank", json=body)
    assert r.status_code == 200
    j = r.json()
    assert j["genes"][0]["symbol"] == "FBN1" and "open_world" in j["uncertainty"] and j["evidence"]
    n = client.post("/next-phenotypes", json={**body, "n": 3}).json()
    assert len(n["next_phenotypes"]) <= 3
    h = client.post("/hypotheses", json=body).json()
    assert "hypotheses" in h and h["research_use_only"].startswith("Research use only")
    a = client.post("/analyze", json=body).json()
    assert {"diseases", "genes", "next_phenotypes", "hypotheses"} <= set(a)


def test_validation_errors(client):
    assert client.post("/rank", json={"present": []}).status_code == 422
    assert client.post("/rank", json={"present": ["HP:0000006"]}).status_code == 422
    assert client.post("/rank", json={"present": ["HP:0001166"], "mode": "x"}).status_code == 422


def test_method_card_and_static(client):
    m = client.get("/method").json()
    assert m["engine_version"] and m["limitations"] and "Research use only" in m["research_use_only"]
    assert client.get("/openapi.json").status_code == 200
    page = client.get("/")
    assert page.status_code == 200 and "Research use only" in page.text
