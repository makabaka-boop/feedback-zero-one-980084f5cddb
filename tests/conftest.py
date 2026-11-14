"""pytest 共享设施：夹具生成、CLI 运行、清单加载。"""

import json
import os
import subprocess
import sys

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

FIX = os.path.join(ROOT, "tests", "fixtures")


@pytest.fixture(scope="session", autouse=True)
def fixtures():
    """会话开始时（重新）生成全部夹具（确定性，幂等）。"""
    sys.path.insert(0, os.path.join(ROOT, "tests"))
    import make_fixtures
    make_fixtures.main()
    return FIX


@pytest.fixture(scope="session")
def manifest(fixtures):
    with open(os.path.join(fixtures, "manifest.json"), "rb") as fh:
        return json.load(fh)


def run_cli(in_path, out_path):
    """以子进程运行 CLI 入口；cwd 为仓库根（python -m 据此找到 bamproj）。"""
    return subprocess.run(
        [sys.executable, "-m", "bamproj", str(in_path), str(out_path)],
        cwd=ROOT, capture_output=True, text=True, timeout=120)
