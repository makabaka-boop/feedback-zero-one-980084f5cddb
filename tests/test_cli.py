"""CLI 行为：原子替换、失败保留旧文件、退出码、库 API。"""

import json
import os
import subprocess
import sys

import pytest

from conftest import ROOT, run_cli

import bamproj
from bamproj import bgzf as bgzf_mod
from bamproj import bam as bam_mod


def test_success_replaces_output_atomically(manifest, fixtures, tmp_path):
    out = tmp_path / "out.json"
    out.write_text("STALE")
    res = run_cli(os.path.join(fixtures, "valid_main.bam"), out)
    assert res.returncode == 0, res.stderr
    assert res.stdout == "" and res.stderr == ""
    got = json.loads(out.read_text())
    assert got[0]["columns"] == [[0, 0, "A"], [1, 2, "G"]]
    assert not list(tmp_path.glob("*.tmp.*"))


def test_failure_keeps_old_file(manifest, fixtures, tmp_path):
    out = tmp_path / "out.json"
    out.write_text("STALE")
    res = run_cli(os.path.join(fixtures, "rec_dup_aux.bam"), out)
    assert res.returncode == 1
    assert out.read_text() == "STALE"


def test_no_output_created_on_failure(manifest, fixtures, tmp_path):
    out = tmp_path / "out.json"
    res = run_cli(os.path.join(fixtures, "bgzf_bad_crc.bam"), out)
    assert res.returncode == 1
    assert not out.exists()


def test_missing_input(manifest, fixtures, tmp_path):
    out = tmp_path / "out.json"
    res = run_cli(os.path.join(fixtures, "does_not_exist.bam"), out)
    assert res.returncode == 1
    assert res.stderr.strip() == "0"
    assert not out.exists()


def test_usage_error():
    res = subprocess.run([sys.executable, "-m", "bamproj"], cwd=ROOT,
                         capture_output=True, text=True)
    assert res.returncode == 2
    assert "usage" in res.stderr


def test_library_api_offsets(manifest, fixtures, tmp_path):
    """库 API：异常携带 .offset；头错为 0，块错为压缩起点，记录错为 vo。"""
    with pytest.raises(bgzf_mod.BgzfError) as ei:
        bamproj.convert(os.path.join(fixtures, "bgzf_second_block_crc.bam"),
                        str(tmp_path / "o1.json"))
    assert ei.value.offset == manifest["bgzf_second_block_crc"]["offset"] > 0

    with pytest.raises(bam_mod.HeaderError) as ei:
        bamproj.convert(os.path.join(fixtures, "hdr_ref_dup.bam"),
                        str(tmp_path / "o2.json"))
    assert ei.value.offset == 0

    with pytest.raises(bam_mod.RecordError) as ei:
        bamproj.convert(os.path.join(fixtures, "rec_qual_mixed.bam"),
                        str(tmp_path / "o3.json"))
    assert ei.value.offset == manifest["rec_qual_mixed"]["offset"]

    # 成功路径：返回 None，输出已写入
    out = tmp_path / "ok.json"
    assert bamproj.convert(os.path.join(fixtures, "valid_main.bam"),
                           str(out)) is None
    assert json.loads(out.read_text())[1]["vo"].isdigit()


def test_output_to_missing_directory(manifest, fixtures, tmp_path):
    out = tmp_path / "no_such_dir" / "out.json"
    res = run_cli(os.path.join(fixtures, "valid_main.bam"), out)
    assert res.returncode == 1
    assert not out.exists()
