#!/usr/bin/env python3
"""JS-vs-Python parity check on fixed patients.

Runs core/genoarch_core/engine.py and web/js/engine.js (via Node) on the same
index, calibration and patients, then compares the full analysis outputs:
strings, integers and booleans must be equal; floats must be equal within
|a - b| <= tol * max(1, |a|, |b|) (tol 1e-9). Ranking-defining floats (disease
scores, null statistics, hypothesis scores) are additionally reported as
bit-identical or not.

Usage: python scripts/parity_check.py --patients P.json --out DIR
"""
from __future__ import annotations

import argparse
import json
import math
import platform
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "core"))

from genoarch_core.engine import Engine, QueryError  # noqa: E402
from genoarch_core.index_io import read_index  # noqa: E402

SEARCH_QUERIES = ["seizure", "long fingers", "HP:0001166", "ectopia", "short stature", "hypotonia", "cleft pal", "0000252"]
EXACT_KEYS = {"score", "z_max", "observed_max_score", "variance_floor", "z_null", "p_value", "p_value_raw_max", "ic"}


def compare(a, b, path, tol, stats, mism):
    if isinstance(a, dict) and isinstance(b, dict):
        if set(a) != set(b):
            mism.append(f"{path}: keys differ {sorted(set(a) ^ set(b))}")
        for k in sorted(set(a) & set(b)):
            compare(a[k], b[k], f"{path}.{k}", tol, stats, mism)
    elif isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b):
            mism.append(f"{path}: length {len(a)} vs {len(b)}")
        for i, (x, y) in enumerate(zip(a, b)):
            compare(x, y, f"{path}[{i}]", tol, stats, mism)
    elif isinstance(a, bool) or isinstance(b, bool) or a is None or b is None or isinstance(a, str) or isinstance(b, str):
        if a != b:
            mism.append(f"{path}: {a!r} vs {b!r}")
    elif isinstance(a, (int, float)) and isinstance(b, (int, float)):
        stats["floats"] += 1
        if float(a) == float(b):
            stats["identical"] += 1
        else:
            d = abs(float(a) - float(b))
            stats["max_abs_diff"] = max(stats["max_abs_diff"], d)
            key = path.rsplit(".", 1)[-1]
            if key in EXACT_KEYS:
                stats["ranking_float_not_identical"] += 1
            if d > tol * max(1.0, abs(float(a)), abs(float(b))) or (math.isnan(float(a)) != math.isnan(float(b))):
                mism.append(f"{path}: {a!r} vs {b!r}")
    else:
        mism.append(f"{path}: type {type(a).__name__} vs {type(b).__name__}")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--index", default=str(ROOT / "web" / "data" / "index.json.gz"))
    ap.add_argument("--calibration", default=str(ROOT / "web" / "data" / "calibration.json"))
    ap.add_argument("--patients", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--tol", type=float, default=1e-9)
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    patients = json.loads(Path(args.patients).read_text())
    cal = json.loads(Path(args.calibration).read_text()) if Path(args.calibration).exists() else {}
    kb = read_index(args.index)
    eng = Engine(kb, cal)
    py_results, timings = [], []
    for p in patients:
        t = time.time()
        try:
            res = eng.analyze(p["present"], p["absent"], **(p.get("options") or {}))
        except QueryError as exc:
            res = {"error": str(exc)}
        timings.append(round((time.time() - t) * 1000))
        py_results.append({"name": p["name"], "result": res})
    py_search = [{"query": q, "result": eng.search(q, 10)} for q in SEARCH_QUERIES]
    spec = out / "parity_spec.json"
    spec.write_text(json.dumps({"patients": patients, "search": SEARCH_QUERIES}))
    (out / "parity_python.json").write_text(json.dumps({"engine": "python", "python": platform.python_version(),
                                                       "analyze_ms": timings, "results": py_results, "search": py_search}))
    js_out = out / "parity_js.json"
    proc = subprocess.run(["node", str(ROOT / "scripts" / "parity_check.mjs"), args.index, args.calibration, str(spec), str(js_out)],
                          capture_output=True, text=True)
    if proc.returncode != 0:
        print(proc.stdout, proc.stderr)
        return 2
    js = json.loads(js_out.read_text())
    per_patient, all_ok = [], True
    totals = {"floats": 0, "identical": 0, "max_abs_diff": 0.0, "ranking_float_not_identical": 0}
    for pr, jr in zip(py_results, js["results"]):
        stats = {"floats": 0, "identical": 0, "max_abs_diff": 0.0, "ranking_float_not_identical": 0}
        mism: list[str] = []
        compare(pr["result"], jr["result"], pr["name"], args.tol, stats, mism)
        ok = not mism
        all_ok &= ok
        top_py = [d["id"] for d in pr["result"].get("diseases", [])]
        top_js = [d["id"] for d in jr["result"].get("diseases", [])]
        genes_py = [g["symbol"] for g in pr["result"].get("genes", [])]
        genes_js = [g["symbol"] for g in jr["result"].get("genes", [])]
        per_patient.append({"name": pr["name"], "pass": ok, "mismatches": mism[:20], **stats,
                            "top20_disease_ids_identical": top_py == top_js,
                            "top20_gene_symbols_identical": genes_py == genes_js,
                            "abstain_identical": (pr["result"].get("uncertainty", {}).get("open_world", {}).get("abstain")
                                                  == jr["result"].get("uncertainty", {}).get("open_world", {}).get("abstain"))})
        for k in ("floats", "identical", "ranking_float_not_identical"):
            totals[k] += stats[k]
        totals["max_abs_diff"] = max(totals["max_abs_diff"], stats["max_abs_diff"])
    search_mism: list[str] = []
    sstats = {"floats": 0, "identical": 0, "max_abs_diff": 0.0, "ranking_float_not_identical": 0}
    compare(py_search, js["search"], "search", args.tol, sstats, search_mism)
    report = {
        "created_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "n_patients": len(patients), "n_pass": sum(p["pass"] for p in per_patient), "all_pass": all_ok and not search_mism,
        "tolerance": args.tol, "totals": totals, "search_queries": len(SEARCH_QUERIES), "search_pass": not search_mism,
        "search_mismatches": search_mism[:20], "python_ms": timings, "js_ms": js["analyze_ms"], "js_load_ms": js["load_ms"],
        "node": js.get("node"), "python": platform.python_version(), "patients": per_patient,
    }
    (out / "parity_report.json").write_text(json.dumps(report, indent=1) + "\n")
    print(json.dumps({k: v for k, v in report.items() if k != "patients"}, indent=1))
    return 0 if report["all_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
