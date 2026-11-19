"""JS engine (web/js/engine.js) vs Python engine on the mini index (requires Node)."""
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))

from parity_check import compare  # noqa: E402

PATIENTS = [
    {"name": "marfan", "present": ["HP:0001166", "HP:0001083"], "absent": ["HP:0001250"], "options": {}},
    {"name": "epilepsy", "present": ["HP:0002373", "HP:0001249"], "absent": [], "options": {"mode": "symmetric"}},
    {"name": "alt-id+conflict", "present": ["HP:0009999", "HP:0007359"], "absent": ["HP:0001250", "HP:0000545"],
     "options": {"universe": "omim", "beta": 1.0}},
    {"name": "excluded", "present": ["HP:0000545"], "absent": [], "options": {"exclude_diseases": ["OMIM:100003"], "mode": "asymmetric"}},
]


def test_js_matches_python_on_mini_index(mini_engine, mini_index, node, tmp_path):
    spec = tmp_path / "spec.json"
    spec.write_text(json.dumps({"patients": PATIENTS, "search": ["seizure", "long slender", "HP:0009999"]}))
    out = tmp_path / "js.json"
    proc = subprocess.run([node, str(ROOT / "scripts" / "parity_check.mjs"), str(mini_index), str(mini_index.parent / "calibration.json"),
                           str(spec), str(out)], capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr
    js = json.loads(out.read_text())
    for p, jr in zip(PATIENTS, js["results"]):
        py = mini_engine.analyze(p["present"], p["absent"], **p["options"])
        stats = {"floats": 0, "identical": 0, "max_abs_diff": 0.0, "ranking_float_not_identical": 0}
        mism = []
        compare(py, jr["result"], p["name"], 1e-9, stats, mism)
        assert not mism, mism[:5]
        assert stats["ranking_float_not_identical"] == 0
    for q, jr in zip(["seizure", "long slender", "HP:0009999"], js["search"]):
        assert mini_engine.search(q, 10) == jr["result"]
