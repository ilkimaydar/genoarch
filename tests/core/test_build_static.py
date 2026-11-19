"""The static build must reproduce the committed method card, including the reliability bins.

CI builds without --eval-summary; the bins then come from web/data/eval_reliability.json.
"""
import json
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def test_build_without_eval_summary_keeps_reliability_bins(tmp_path):
    web = tmp_path / "web"
    shutil.copytree(ROOT / "web", web)
    out = tmp_path / "out"
    subprocess.run([sys.executable, str(ROOT / "scripts" / "build_static.py"), "--web", str(web), "--out", str(out)],
                   check=True, capture_output=True)
    built = json.loads((out / "site" / "data" / "method_card.json").read_text())
    committed = json.loads((ROOT / "web" / "data" / "method_card.json").read_text())
    assert len(built["evaluation"]["reliability_10"]) == 10
    assert built == committed
