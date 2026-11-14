"""合法输入验收：投影对照 pysam，vo 另由测试侧独立走查核对。"""

import json
import os
import struct
import zlib

import pytest

from conftest import run_cli

pysam = pytest.importorskip("pysam", reason="pysam 未安装，跳过对照")


def pysam_projection(path):
    """用 pysam 独立推导期望输出（含经 tell() 取得的虚拟偏移）。"""
    af = pysam.AlignmentFile(str(path), "rb")
    refs = af.references
    out = []
    it = iter(af)
    while True:
        vo = af.tell()
        try:
            r = next(it)
        except StopIteration:
            break
        if r.is_unmapped or r.query_length == 0:
            continue
        cols = []
        pos = r.reference_start
        q = 0
        seq = r.query_sequence
        for op, ln in r.cigartuples:
            if op in (0, 7, 8):  # M = X
                for _ in range(ln):
                    cols.append([pos, q, seq[q]])
                    pos += 1
                    q += 1
            elif op in (1, 4):  # I S
                q += ln
            elif op in (2, 3):  # D N
                pos += ln
            # H P 不推进
        out.append({"name": r.query_name, "ref": refs[r.reference_id],
                    "vo": str(vo), "columns": cols})
    af.close()
    return out


def walk_collected_vos(comp):
    """测试侧独立 BGZF/BAM 走查：返回收录记录（未比对/空 SEQ 除外）的
    [(vo, name)]，仅依赖 zlib，与被测实现相互独立。"""
    blocks = []
    out = bytearray()
    off = 0
    while off < len(comp):
        cstart = off
        xlen = struct.unpack_from("<H", comp, off + 10)[0]
        p = off + 12
        bsize = None
        while p < off + 12 + xlen:
            slen = struct.unpack_from("<H", comp, p + 2)[0]
            if comp[p:p + 2] == b"BC":
                bsize = struct.unpack_from("<H", comp, p + 4)[0]
            p += 4 + slen
        total = bsize + 1
        payload = comp[off + 12 + xlen: off + total - 8]
        d = zlib.decompressobj(-15)
        plain = d.decompress(payload) + d.flush()
        if plain:
            blocks.append((len(out), len(plain), cstart))
            out += plain
        off += total
    data = bytes(out)

    def vo_of(u):
        for us, sz, cs in blocks:
            if us <= u < us + sz:
                return (cs << 16) | (u - us)
        raise AssertionError("offset outside blocks")

    l_text = struct.unpack_from("<i", data, 4)[0]
    pos = 8 + l_text
    n_ref = struct.unpack_from("<i", data, pos)[0]
    pos += 4
    for _ in range(n_ref):
        l_name = struct.unpack_from("<i", data, pos)[0]
        pos += 4 + l_name + 4
    got = []
    while pos < len(data):
        u = pos
        bs = struct.unpack_from("<i", data, pos)[0]
        flag = struct.unpack_from("<H", data, pos + 4 + 14)[0]
        l_seq = struct.unpack_from("<i", data, pos + 4 + 16)[0]
        l_read_name = data[pos + 4 + 8]
        name = data[pos + 4 + 32: pos + 4 + 32 + l_read_name - 1].decode("utf-8")
        if not (flag & 0x4) and l_seq > 0:
            got.append((vo_of(u), name))
        pos += 4 + bs
    return got


def valid_entries(manifest):
    return sorted(k for k, v in manifest.items() if v["status"] == "ok")


def test_task_example(manifest, fixtures, tmp_path):
    """任务书示例：pos=0、SEQ=ACG、CIGAR=1M1I1M → [[0,0,"A"],[1,2,"G"]]。"""
    out = tmp_path / "out.json"
    res = run_cli(os.path.join(fixtures, "valid_main.bam"), out)
    assert res.returncode == 0, res.stderr
    got = json.loads(out.read_text())
    r001 = got[0]
    assert r001["name"] == "r001"
    assert r001["columns"] == [[0, 0, "A"], [1, 2, "G"]]


@pytest.mark.parametrize("name", ["valid_main", "valid_multiblock",
                                  "valid_random", "valid_cg", "valid_empty",
                                  "valid_all_unmapped", "valid_no_text"])
def test_valid(name, manifest, fixtures, tmp_path):
    entry = manifest[name]
    bam = os.path.join(fixtures, entry["bam"])
    out = tmp_path / "out.json"
    res = run_cli(bam, out)
    assert res.returncode == 0, res.stderr
    got = json.loads(out.read_text())
    # 1) 与夹具生成器（独立写出端）的期望一致
    expected = json.loads(
        open(os.path.join(fixtures, entry["expected"]), "rb").read())
    assert got == expected
    # 2) 投影（含 vo）与 pysam 对照一致
    assert got == pysam_projection(bam)
    # 3) vo 与测试侧独立走查一致（重名按文件序对齐）
    with open(bam, "rb") as fh:
        comp = fh.read()
    assert [(int(g["vo"]), g["name"]) for g in got] == walk_collected_vos(comp)
    # 4) 输出形状：文件序 {name, ref, vo, columns}
    for item in got:
        assert list(item.keys()) == ["name", "ref", "vo", "columns"]
        assert isinstance(item["vo"], str) and item["vo"].isdigit()
        for r, q, b in item["columns"]:
            assert isinstance(r, int) and isinstance(q, int)
            assert isinstance(b, str) and len(b) == 1


def test_empty_output_is_empty_array(manifest, fixtures, tmp_path):
    out = tmp_path / "out.json"
    res = run_cli(os.path.join(fixtures, "valid_empty.bam"), out)
    assert res.returncode == 0, res.stderr
    assert json.loads(out.read_text()) == []
