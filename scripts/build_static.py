#!/usr/bin/env python3
"""Build the static deliverables.

* web/data/method_card.json  method card for the Method page (same content as GET /method)
* dist/site/                 multi-file static site (GitHub Pages, any static host)
* dist/hf-space/             dist/site plus the Hugging Face Space README (sdk: static)
* dist/genoarch.html         single self-contained page (index, calibration and method card embedded)

Usage: python scripts/build_static.py [--eval-summary ../genoarch_proj/exp/platform_eval/run_X/summary.json]
"""
from __future__ import annotations

import argparse
import base64
import json
import re
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "core"))

from genoarch_core.engine import Engine  # noqa: E402
from genoarch_core.index_io import read_index  # noqa: E402
from genoarch_core.method_card import method_card  # noqa: E402

HF_README = """---
title: Genoarch
colorFrom: indigo
colorTo: pink
sdk: static
app_file: index.html
pinned: false
license: mit
short_description: HPO gene and disease ranking with calibrated uncertainty
---

# Genoarch

Rank candidate diseases and genes for undiagnosed cases from HPO phenotype terms, with calibrated confidence, an
open-world (abstain) test, next-phenotype suggestions and gene-association hypotheses. The analysis runs entirely
in the browser on a compact index built from the HPO release {hpo}.

Research use only. Not a medical device. Not for diagnosis or treatment decisions.

Source: https://github.com/ilkimaydar/genoarch (MIT).
"""


def inline_json(obj) -> str:
    return json.dumps(obj, separators=(",", ":"), ensure_ascii=False).replace("</", "<\\/")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--web", default=str(ROOT / "web"))
    ap.add_argument("--out", default=str(ROOT / "dist"))
    ap.add_argument("--eval-summary", default=None, help="summary.json of the evaluation run (adds reliability bins)")
    args = ap.parse_args()
    web, out = Path(args.web), Path(args.out)
    data = web / "data"
    cal = json.loads((data / "calibration.json").read_text())
    kb = read_index(data / "index.json.gz")
    card = method_card(Engine(kb, cal))
    if args.eval_summary and card.get("evaluation"):
        summary = json.loads(Path(args.eval_summary).read_text())
        key = f"{cal['default_mode']}/{cal['default_universe']}"
        card["evaluation"]["reliability_10"] = summary["test_results"][key]["reliability_10"]
    (data / "method_card.json").write_text(json.dumps(card, indent=1, ensure_ascii=False) + "\n")

    site = out / "site"
    if site.exists():
        shutil.rmtree(site)
    shutil.copytree(web, site)
    hf = out / "hf-space"
    if hf.exists():
        shutil.rmtree(hf)
    shutil.copytree(site, hf)
    (hf / "README.md").write_text(HF_README.format(hpo=kb.meta.get("ontology_version")))

    html = (web / "index.html").read_text()
    css = (web / "css" / "app.css").read_text()
    eng_js = (web / "js" / "engine.js").read_text().replace("</script", "<\\/script")
    app_js = (web / "js" / "app.js").read_text().replace("</script", "<\\/script")
    b64 = base64.b64encode((data / "index.json.gz").read_bytes()).decode("ascii")
    b64 = "\n".join(b64[i:i + 4096] for i in range(0, len(b64), 4096))
    html = html.replace('<link rel="stylesheet" href="css/app.css">', f"<style>\n{css}</style>")
    payload = (f'<script type="application/octet-stream" id="ga-index">\n{b64}\n</script>\n'
               f'<script type="application/json" id="ga-calibration">{inline_json(cal)}</script>\n'
               f'<script type="application/json" id="ga-method">{inline_json(card)}</script>\n')
    html, n1 = re.subn(r'<script src="js/engine\.js"></script>', lambda m: payload + f"<script>\n{eng_js}</script>", html)
    html, n2 = re.subn(r'<script src="js/app\.js"></script>', lambda m: f"<script>\n{app_js}</script>", html)
    if n1 != 1 or n2 != 1:
        raise SystemExit("script tags not found in index.html")
    single = out / "genoarch.html"
    single.write_text(html)
    report = {"single_file": str(single), "single_file_bytes": single.stat().st_size,
              "index_gzip_bytes": (data / "index.json.gz").stat().st_size, "site": str(site), "hf_space": str(hf),
              "method_card": str(data / "method_card.json")}
    print(json.dumps(report, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
