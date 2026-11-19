import json
import shutil
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "core"))
sys.path.insert(0, str(ROOT / "api"))

from genoarch_core.engine import Engine  # noqa: E402
from genoarch_core.index_io import read_index, write_index  # noqa: E402
from genoarch_core.knowledge import KnowledgeBase  # noqa: E402

MINI = ROOT / "tests" / "fixtures" / "mini"
REAL_INDEX = ROOT / "web" / "data" / "index.json.gz"
REAL_CAL = ROOT / "web" / "data" / "calibration.json"
MINI_CAL = {"default_mode": "blended", "default_universe": "all",
            "params": {"blended": {"all": {"T": 0.5, "beta": 0.5, "w": 0.1}, "omim": {"T": 0.5, "beta": 0.5, "w": 0.1}}},
            "constants": {"null_samples": 40, "hyp_min_ic": 1.0, "eig_min_ic": 0.0, "hyp_top_diseases": 1}}


@pytest.fixture(scope="session")
def mini_dir():
    return MINI


@pytest.fixture(scope="session")
def mini_kb():
    return KnowledgeBase.from_hpo_files(MINI)


@pytest.fixture(scope="session")
def mini_index(tmp_path_factory, mini_kb):
    path = tmp_path_factory.mktemp("index") / "index.json.gz"
    write_index(mini_kb, path)
    cal = path.parent / "calibration.json"
    cal.write_text(json.dumps(MINI_CAL))
    return path


@pytest.fixture(scope="session")
def mini_engine(mini_index):
    return Engine(read_index(mini_index), MINI_CAL)


@pytest.fixture(scope="session")
def real_engine():
    if not REAL_INDEX.exists():
        pytest.skip("real index not built (run make index)")
    cal = json.loads(REAL_CAL.read_text()) if REAL_CAL.exists() else {}
    return Engine(read_index(REAL_INDEX), cal)


@pytest.fixture(scope="session")
def node():
    exe = shutil.which("node")
    if not exe:
        pytest.skip("node not installed")
    return exe
